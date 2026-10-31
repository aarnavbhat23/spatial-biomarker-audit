"""Retrieve one complete field from a public ZIP via bounded byte ranges."""
import concurrent.futures
import hashlib
import json
from pathlib import Path
import re
import struct
import time
import urllib.request
import zlib

ROOT = Path(__file__).resolve().parent
URL = 'https://zenodo.org/records/15678603/files/MERFISH_Tissue_raw_unaligned_images.zip?download=1'
SIZE = 170983659800
CHUNK = 8*1024*1024


def response(start, stop):
    req = urllib.request.Request(URL, headers={'Range': f'bytes={start}-{stop}'})
    result = urllib.request.urlopen(req, timeout=60)
    if result.status != 206 or result.headers.get('Content-Range') != f'bytes {start}-{stop}/{SIZE}':
        result.close()
        raise ValueError('Server did not return the exact requested range')
    return result


def get_chunk(bounds):
    start, stop = bounds
    cache = ROOT/'spatial/chunks'/f'{start}-{stop}.bin'
    cache.parent.mkdir(parents=True, exist_ok=True)
    if cache.exists() and cache.stat().st_size == stop-start+1:
        return cache
    for attempt in range(3):
        try:
            with response(start, stop) as incoming:
                blob = incoming.read(stop-start+2)
            if len(blob) != stop-start+1:
                raise ValueError(f'Truncated range {start}: {len(blob)} bytes')
            cache.write_bytes(blob)
            return cache
        except Exception:
            if attempt == 2:
                raise
            time.sleep(2**attempt)


def acquire(entry):
    name = Path(entry['name']).name
    if not re.fullmatch(r'merFISH_merged_0[1-9]_000.tiff', name) or entry['compression'] != 8:
        raise ValueError('Entry is outside the selected raw-field scope')
    destination = ROOT/'spatial/raw'/name
    receipt = destination.with_suffix('.receipt.json')
    if destination.exists() and receipt.exists():
        value = json.loads(receipt.read_text())
        with destination.open('rb') as handle:
            digest = hashlib.file_digest(handle, 'sha256').hexdigest()
        if (digest == value['sha256'] and destination.stat().st_size == entry['size']
            and value['zip_crc32'] == entry['crc32'] and value['archive_entry'] == entry['name']):
            print('Verified existing', name, flush=True)
            return value
        raise ValueError('Existing source differs from its receipt')
    destination.parent.mkdir(parents=True, exist_ok=True)
    with response(entry['offset'], entry['offset']+1023) as incoming:
        header = incoming.read(1024)
    if header[:4] != b'PK\x03\x04':
        raise ValueError('Invalid ZIP local header')
    name_length, extra_length = struct.unpack_from('<HH', header, 26)
    if header[30:30+name_length].decode() != entry['name']:
        raise ValueError('Local-header filename differs from central directory')
    start = entry['offset'] + 30 + name_length + extra_length
    temporary = destination.with_suffix('.partial')
    inflater = zlib.decompressobj(-15)
    crc = total = compressed = 0
    sha = hashlib.sha256()
    bounds = [(x, min(x+CHUNK-1, start+entry['compressed']-1)) for x in range(start, start+entry['compressed'], CHUNK)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as executor:
        chunks = list(executor.map(get_chunk, bounds))
    with temporary.open('wb') as outgoing:
        for chunk in chunks:
            block = chunk.read_bytes()
            compressed += len(block)
            expanded = inflater.decompress(block)
            total += len(expanded)
            if compressed > entry['compressed'] or total > entry['size']:
                raise ValueError('Source exceeds declared size')
            outgoing.write(expanded); sha.update(expanded); crc = zlib.crc32(expanded, crc)
        tail = inflater.flush()
        outgoing.write(tail); sha.update(tail); crc = zlib.crc32(tail, crc); total += len(tail)
    if not inflater.eof or inflater.unused_data:
        raise ValueError('Incomplete or trailing compressed stream')
    if (compressed, total, crc) != (entry['compressed'], entry['size'], entry['crc32']):
        raise ValueError('ZIP member size or CRC32 mismatch')
    temporary.replace(destination)
    value = dict(source_url=URL, archive_bytes=SIZE, archive_entry=entry['name'],
                 byte_range=[start, start+entry['compressed']-1], compressed_bytes=compressed,
                 bytes=total, zip_crc32=crc, sha256=sha.hexdigest(),
                 license='CC-BY-4.0', source_doi='10.5281/zenodo.15678603',
                 split='development', independent_biological_unit='unknown', field_id='000')
    receipt.write_text(json.dumps(value, indent=2)+'\n')
    for chunk in chunks:
        chunk.unlink()
    print('Acquired', name, total, sha.hexdigest(), flush=True)
    return value


def main():
    entries = json.loads((ROOT/'source/spatial_archive_index.json').read_text())
    selected = [e for e in entries if re.fullmatch(r'merlin_data/merFISH_merged_0[1-9]_000.tiff', e['name'])]
    if len(selected) != 9 or sum(e['size'] for e in selected) >= 4_000_000_000:
        raise ValueError('Selected source count or storage cap exceeded')
    receipts = [acquire(entry) for entry in selected]
    (ROOT/'spatial/source_manifest.json').write_text(json.dumps(receipts, indent=2)+'\n')


if __name__ == '__main__':
    main()
