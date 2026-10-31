"""Fetch small public metadata and the ZIP index without downloading the archive."""
import hashlib
import io
import json
from pathlib import Path
import struct
import urllib.request
import zipfile

from acquire_spatial import ROOT, SIZE, response

CHECKSUMS = {'codebook.csv': '80753788495b287711dbbed533c13e96',
             'dataorganization.csv': 'fce81a8f522d59f04b21be9fd2daebe2',
             'microscope.json': '3819047ac6b973a2050193fb4b92ac60',
             'positions.csv': '6f11a3b825fbf1d6ecf0886a08ac4619'}


def fetch(url, limit=2_000_000):
    request = urllib.request.Request(url, headers={'User-Agent': 'measurement-dataset-curation/0.1'})
    with urllib.request.urlopen(request, timeout=60) as response:
        data = response.read(limit+1)
    if len(data) > limit:
        raise ValueError('Source exceeds bounded metadata size')
    return data


def parse_index(tail):
    start = SIZE-len(tail)
    position = tail.rfind(b'PK\x06\x06')
    if position < 0:
        raise ValueError('ZIP64 end record absent')
    record = struct.unpack_from('<4sQ2H2L4Q', tail, position)
    _, _, _, _, disk, central_disk, disk_entries, count, length, offset = record
    if (disk, central_disk, disk_entries, count) != (0, 0, 488, 488):
        raise ValueError('Unexpected archive topology')
    if offset < start or offset+length > SIZE or offset+length != start+position:
        raise ValueError('Central directory lies outside verified tail bounds')
    directory = tail[offset-start:offset-start+length]
    # This temporary in-memory ZIP contains just the central directory. Its
    # entries retain the original ZIP64 local-header offsets used for fetching.
    eocd = struct.pack('<4s4H2LH', b'PK\x05\x06', 0, 0, count, count, length, 0, 0)
    with zipfile.ZipFile(io.BytesIO(directory+eocd)) as archive:
        return [dict(name=x.filename, compressed=x.compress_size, size=x.file_size,
                     offset=x.header_offset, compression=x.compress_type, crc32=x.CRC)
                for x in archive.infolist()]


def main():
    source, metadata = ROOT/'source', ROOT/'spatial/metadata'
    source.mkdir(exist_ok=True); metadata.mkdir(parents=True, exist_ok=True)
    record_path = metadata/'zenodo_record.json'
    if not record_path.exists():
        record_path.write_bytes(fetch('https://zenodo.org/api/records/15678603'))
    record = json.loads(record_path.read_text())
    if record['metadata']['doi'] != '10.5281/zenodo.15678603':
        raise ValueError('Wrong source record')
    for name, checksum in CHECKSUMS.items():
        path = metadata/name
        data = path.read_bytes() if path.exists() else fetch(f'https://zenodo.org/records/15678603/files/{name}?download=1')
        if hashlib.md5(data).hexdigest() != checksum:
            raise ValueError(f'Publisher MD5 mismatch: {name}')
        path.write_bytes(data)
    tail_path = source/'spatial_zip_tail.bin'
    if not tail_path.exists():
        with response(SIZE-131072, SIZE-1) as stream:
            tail = stream.read(131073)
        if len(tail) != 131072:
            raise ValueError('Truncated ZIP tail')
        tail_path.write_bytes(tail)
    entries = parse_index(tail_path.read_bytes())
    index_path = source/'spatial_archive_index.json'
    if index_path.exists() and json.loads(index_path.read_text()) != entries:
        raise ValueError('Archive index differs from pinned release')
    index_path.write_text(json.dumps(entries, indent=2)+'\n')
    print('Verified source metadata and 488 ZIP index entries')


if __name__ == '__main__':
    main()
