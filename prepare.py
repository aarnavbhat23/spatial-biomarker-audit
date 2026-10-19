"""Build a small, patient-disjoint H&E/marker pilot with complete local graphs."""
import concurrent.futures
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path

import fsspec
import numpy as np
import pandas as pd
import tifffile
import zarr
from scipy.ndimage import uniform_filter
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent
SPLITS = {'CRC01':'train', 'CRC03':'train', 'CRC19':'validation', 'CRC02':'test', 'CRC11':'test'}
SPLITS.update({'CRC04':'confirmation','CRC12':'confirmation'})
TILE_COUNTS = {'train':24, 'validation':12, 'test':32, 'confirmation':32}
CORE = 2048
RADIUS = 80
LEVEL = 2
DOWNSAMPLE = 4
FEATURE_MARGIN = 32
MARKERS = ['Pan-CK_pos', 'CD3e_pos']


def image_features(crop, x, y):
    """Fixed image-only color/optical-density descriptors at three local scales."""
    rgb = crop.astype(np.float32) / 255
    od = -np.log((crop.astype(np.float32) + 1) / 256)
    features = [rgb[y, x], od[y, x]]
    for image in (rgb, od):
        for width in (3, 9, 17):
            mean = uniform_filter(image, size=(width, width, 1), mode='reflect')
            second = uniform_filter(image * image, size=(width, width, 1), mode='reflect')
            features.extend([mean[y, x], np.sqrt(np.maximum(second[y, x] - mean[y, x] ** 2, 0))])
    return np.concatenate(features, axis=1)


def prepare_patient(patient):
    source = json.loads((ROOT/'data/manifest.json').read_text())[patient]
    destination = ROOT/'data'/f'{patient}-prepared.npz'
    record_path = ROOT/'data'/f'{patient}-preparation.json'
    if destination.exists() and record_path.exists():
        print(patient, 'using cached prepared data', flush=True)
        return patient, json.loads(record_path.read_text())
    table = pd.read_csv(ROOT/'data'/f'{patient}-cells.csv',
                        usecols=['label','x','y'] + MARKERS,
                        dtype={'label':'int32','x':'int32','y':'int32', **{m:'bool' for m in MARKERS}})
    split = SPLITS[patient]
    rng = np.random.default_rng(1729 + int(patient[3:]))
    xy = table[['x','y']].to_numpy()
    tile_xy = xy // CORE
    candidates, counts = np.unique(tile_xy, axis=0, return_counts=True)
    feature_rows, cell_rows = [], []
    with fsspec.open(source['image_url'], 'rb', block_size=262144, cache_type='blockcache') as file:
        with tifffile.TiffFile(file) as image:
            native_shape = image.series[0].shape
            if native_shape[-1] != 3 or xy[:,0].max() >= native_shape[1] or xy[:,1].max() >= native_shape[0]:
                raise ValueError('Cell coordinates do not fit the registered H&E image')
            metadata = ET.fromstring(image.ome_metadata)
            pixels = metadata.find('.//{*}Pixels')
            physical = {key:pixels.attrib.get(key) for key in ('PhysicalSizeX','PhysicalSizeXUnit','PhysicalSizeY','PhysicalSizeYUnit')}
            low = image.series[0].levels[LEVEL]
            if low.shape[:2] != tuple(math.ceil(a/DOWNSAMPLE) for a in native_shape[:2]):
                raise ValueError('Unexpected image pyramid downsample')
            margin = RADIUS + FEATURE_MARGIN
            valid = ((counts >= 100) & (candidates[:,0]*CORE >= margin) &
                     (candidates[:,1]*CORE >= margin) &
                     ((candidates[:,0]+1)*CORE + margin < native_shape[1]) &
                     ((candidates[:,1]+1)*CORE + margin < native_shape[0]))
            candidates = candidates[valid]
            tiles = candidates[rng.choice(len(candidates), TILE_COUNTS[split], replace=False)]
            tiles = tiles[np.lexsort((tiles[:,1],tiles[:,0]))]
            store = low.aszarr()
            array = zarr.open(store, mode='r')
            for number, (tx,ty) in enumerate(tiles):
                x0,y0 = int(tx)*CORE, int(ty)*CORE
                include = ((xy[:,0]>=x0-RADIUS) & (xy[:,0]<x0+CORE+RADIUS) &
                           (xy[:,1]>=y0-RADIUS) & (xy[:,1]<y0+CORE+RADIUS))
                rows = np.flatnonzero(include)
                lx,ly = (x0-margin)//DOWNSAMPLE, (y0-margin)//DOWNSAMPLE
                rx,ry = math.ceil((x0+CORE+margin)/DOWNSAMPLE), math.ceil((y0+CORE+margin)/DOWNSAMPLE)
                crop = np.array(array[ly:ry,lx:rx,:])
                local_x = xy[rows,0]//DOWNSAMPLE - lx
                local_y = xy[rows,1]//DOWNSAMPLE - ly
                features = image_features(crop,local_x,local_y)
                feature_rows.append(features)
                cell_rows.append(rows)
                if number%8==0:
                    print(patient,'image regions',number+1,'/',len(tiles),flush=True)
            store.close()
    all_rows = np.concatenate(cell_rows)
    selected, first = np.unique(all_rows, return_index=True)
    features = np.concatenate(feature_rows)[first]
    coords = xy[selected]
    core_tiles = coords//CORE
    tile_map = {tuple(tile):i for i,tile in enumerate(tiles)}
    owner = np.array([tile_map.get(tuple(tile),-1) for tile in core_tiles], dtype=np.int32)
    anchors = np.flatnonzero(owner>=0)
    neighborhoods = cKDTree(coords).query_ball_point(coords[anchors],RADIUS)
    left = np.repeat(anchors,[len(row) for row in neighborhoods])
    right = np.concatenate(neighborhoods)
    distinct = left!=right
    edges = np.stack([left[distinct],right[distinct]],axis=1).astype(np.int32)
    edge_tiles = owner[edges[:,0]]
    rectangles = np.array([[int(x)*CORE-RADIUS,int(y)*CORE-RADIUS,
                            (int(x)+1)*CORE+RADIUS,(int(y)+1)*CORE+RADIUS] for x,y in tiles])
    np.savez_compressed(destination,features=features,labels=table[MARKERS].to_numpy()[selected],
                        coords=coords,cell_ids=table['label'].to_numpy()[selected],
                        edges=edges,edge_tiles=edge_tiles,rectangles=rectangles,core_tiles=tiles)
    record = {**source,'split':split,'n_source_cells':len(table),'n_cells':len(selected),
              'n_anchor_cells':len(anchors),'n_directed_pairs':len(edges),'n_tiles':len(tiles),
              'core_tiles_xy':tiles.tolist(),'native_image_shape':native_shape,
              'pyramid_level':LEVEL,'downsample':DOWNSAMPLE,'physical_pixel_size':physical,
              'core_pixels':CORE,'measurement_halo_pixels':RADIUS,
              'candidate_occupied_tiles':len(candidates),'features':int(features.shape[1]),
              'markers':MARKERS,'graph_dependency':'Source measured-image cell segmentation',
              'scope':'Selected tile panel, with all source cells in its local halos. Not a whole-slide census.'}
    record_path.write_text(json.dumps(record,indent=2)+'\n')
    print(patient,'prepared',len(selected),'cells and',len(edges),'pairs',flush=True)
    return patient,record


def main():
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        records = dict(pool.map(prepare_patient,SPLITS))
    (ROOT/'results/data_manifest.json').write_text(json.dumps(records,indent=2)+'\n')


if __name__ == '__main__':
    main()
