"""Verified byte-range fetch of selected public IMC tables; no image archive."""
import hashlib,json,struct,zlib
from pathlib import Path
from fetch_data import fetch_range

ROOT=Path(__file__).resolve().parent
URL='https://zenodo.org/records/7324285/files/MBTMEStrPublicProcessedDataCode.zip?download=1'


def entries():
    central=fetch_range(URL,910656804,9346);result={};offset=0
    while central[offset:offset+4]==b'PK\x01\x02':
        fixed=struct.unpack_from('<4s6H3L5H2L',central,offset);nl,el,cl=fixed[10:13]
        name=central[offset+46:offset+46+nl].decode()
        result[name]={'name':name,'compressed':fixed[8],'size':fixed[9],'offset':fixed[-1],'crc32':fixed[7]}
        offset+=46+nl+el+cl
    return result


def fetch_member(entry,folder):
    output=folder/Path(entry['name']).name
    if not output.exists():
        header=struct.unpack('<4s5H3I2H',fetch_range(URL,entry['offset'],30))
        if header[0]!=b'PK\x03\x04' or header[3]!=8:raise ValueError('Unexpected ZIP member')
        start=entry['offset']+30+header[-2]+header[-1];decoder=zlib.decompressobj(-15)
        temporary=output.with_suffix(output.suffix+'.partial');crc=0;size=0
        with temporary.open('wb') as target:
            for offset in range(0,entry['compressed'],8*1024**2):
                data=decoder.decompress(fetch_range(URL,start+offset,min(8*1024**2,entry['compressed']-offset)))
                target.write(data);crc=zlib.crc32(data,crc);size+=len(data)
                print(Path(entry['name']).name,offset+min(8*1024**2,entry['compressed']-offset),'/',entry['compressed'],flush=True)
            tail=decoder.flush();target.write(tail);crc=zlib.crc32(tail,crc);size+=len(tail)
        if not decoder.eof or size!=entry['size'] or crc!=entry['crc32']:raise ValueError('ZIP size or CRC mismatch')
        temporary.rename(output)
    raw_hash=hashlib.sha256();crc=0;size=0
    with output.open('rb') as source:
        for chunk in iter(lambda:source.read(8*1024**2),b''):raw_hash.update(chunk);crc=zlib.crc32(chunk,crc);size+=len(chunk)
    if size!=entry['size'] or crc!=entry['crc32']:raise ValueError('Cached file checksum mismatch')
    return {**entry,'sha256':raw_hash.hexdigest(),'path':str(output.relative_to(ROOT))}


def main():
    folder=ROOT/'data/danenberg';folder.mkdir(exist_ok=True);members=entries()
    report={'record_doi':'10.5281/zenodo.7324285','archive_url':URL,'license':'CC-BY-4.0','files':{}}
    for name in ('SingleCells.csv','SingleCellsAnnotation.xlsx','IMCClinical.csv','SI_metaclusters.csv'):
        entry=members['MBTMEIMCPublic/'+name];report['files'][name]=fetch_member(entry,folder)
        (folder/'fetch_manifest.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
