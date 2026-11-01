"""Fetch exact public reference assets, retaining source terms and checksums."""
import base64
import hashlib
import json
from pathlib import Path
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parent/'reference_calibration_v01'
FILES = {'cross_channel_matches.csv': '024fe5233ece13c91db1116c9aa43bbbef791cf0',
         'SOURCE_LICENSE.md': 'b7159508ea58b59bc2db25ebd2373654a34ba826'}


def acquire():
    ROOT.mkdir(exist_ok=True)
    manifest = {}
    for name, blob in FILES.items():
        destination = ROOT/name
        if destination.exists():
            data = destination.read_bytes()
        else:
            url = f'https://api.github.com/repos/CaiGroup/UntanglingBarcodes/git/blobs/{blob}'
            request = Request(url, headers={'User-Agent': 'reference-calibration-research'})
            with urlopen(request, timeout=60) as response:
                source = json.load(response)
            data = base64.b64decode(source['content'])
        actual = hashlib.sha1(f'blob {len(data)}\0'.encode()+data).hexdigest()
        if actual != blob:
            raise ValueError(f'Source Git blob mismatch: {name}')
        destination.write_bytes(data)
        manifest[name] = dict(git_blob=blob, sha256=hashlib.sha256(data).hexdigest(), bytes=len(data))
    (ROOT/'source_checksums.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    acquire()
