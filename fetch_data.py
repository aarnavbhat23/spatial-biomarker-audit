"""Fetch selected public cell tables by ZIP range, never the full image archive."""
import concurrent.futures
import hashlib
import json
import struct
import urllib.request
import xml.etree.ElementTree as ET
import zlib
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCES = json.loads((ROOT / 'data_sources.json').read_text())


def fetch_range(url, start, size):
    request = urllib.request.Request(url, headers={'Range':f'bytes={start}-{start + size - 1}'})
    with urllib.request.urlopen(request, timeout=120) as response:
        if response.status != 206 or not response.headers.get('Content-Range','').startswith(f'bytes {start}-'):
            raise ValueError('Server did not return the requested byte range')
        content = response.read(size + 1)
    if len(content) != size:
        raise ValueError('Incorrect range length')
    return content


def fetch_table(item):
    patient, entry = item
    output = ROOT / 'data' / f'{patient}-cells.csv'
    if not output.exists() or output.stat().st_size != entry['size']:
        header = struct.unpack('<4s5H3I2H', fetch_range(SOURCES['archive'], entry['offset'], 30))
        if header[0] != b'PK\x03\x04' or header[3] != 8:
            raise ValueError('Unexpected local ZIP header')
        start = entry['offset'] + 30 + header[-2] + header[-1]
        content = zlib.decompress(fetch_range(SOURCES['archive'], start, entry['compressed']), -15)
        if len(content) != entry['size'] or zlib.crc32(content) != header[6]:
            raise ValueError('Cell table checksum or length mismatch')
        output.write_bytes(content)
    checksum = hashlib.file_digest(output.open('rb'), 'sha256').hexdigest()
    if entry.get('sha256') and checksum != entry['sha256']:
        raise ValueError('Cell table differs from the pinned source checksum')
    # Discover exact registered-image keys from the authors' bucket listing.
    url = SOURCES['raw_images'] + '?list-type=2&prefix=data/' + patient + '/'
    with urllib.request.urlopen(url, timeout=30) as response:
        tree = ET.fromstring(response.read())
    namespace = {'s':'http://s3.amazonaws.com/doc/2006-03-01/'}
    keys = [node.findtext('s:Key', namespaces=namespace) for node in tree.findall('s:Contents', namespace)]
    image = next(key for key in keys if '-registered.ome.tif' in key)
    record = {'patient':patient, 'cell_table_bytes':entry['size'], 'sha256':checksum,
              'archive_member':entry['name'], 'image_url':SOURCES['raw_images'] + image}
    print(patient, 'cell table verified', flush=True)
    return patient, record


def main():
    (ROOT / 'data').mkdir(exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        records = dict(pool.map(fetch_table, SOURCES['cell_tables'].items()))
    (ROOT / 'data' / 'manifest.json').write_text(json.dumps(records, indent=2) + '\n')


if __name__ == '__main__':
    main()
