"""Record the rebuilt small inputs consumed by the measured barcode screen."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
reference = json.loads((ROOT/'reference_checksums.json').read_text())['files']
records = []
for item in reference:
    path = ROOT/item['path']
    with path.open('rb') as handle:
        digest = hashlib.file_digest(handle, 'sha256').hexdigest()
    if item['path'] in ['spatial/middle_z14_oriented_unaligned.npy', 'spatial/metadata/codebook.csv']:
        if digest != item['sha256']:
            raise ValueError(f"Immutable source-derived input differs: {item['path']}")
    records.append(dict(path=item['path'], bytes=path.stat().st_size, sha256=digest))
(ROOT/'release_checksums.json').write_text(json.dumps(dict(files=records), indent=2)+'\n')
print('Recorded four barcode-screen input hashes; source pixels/codebook match the reference')
