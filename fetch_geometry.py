"""Fetch unchanged CC0 BBBC039 archives and verify the frozen checksums."""
import hashlib
import json
import ssl
import urllib.request
from pathlib import Path

SOURCES = {
    'images': '6f30a5d4fe38c928ded972704f085975f8dc0d65d9aa366df00e5a9d449fddd7',
    'masks': 'f9e6043d8ca56344a4886f96a700d804d6ee982f31e2b2cd3194af2a053c2710',
    'metadata': 'a2c1f900bed9ba92a99553efd4c2ae98598433691c7401d818653ab61110deb2'}


def fetch(root=Path('data/bbbc039')):
    root=Path(root);root.mkdir(parents=True,exist_ok=True);manifest={}
    for name,digest in SOURCES.items():
        url=f'https://data.broadinstitute.org/bbbc/BBBC039/{name}.zip';path=root/(name+'.zip')
        if not path.exists():
            ca=Path('/etc/ssl/cert.pem');context=ssl.create_default_context(cafile=str(ca) if ca.exists() else None)
            path.write_bytes(urllib.request.urlopen(url,context=context,timeout=120).read())
        actual=hashlib.sha256(path.read_bytes()).hexdigest()
        if actual!=digest:raise ValueError(f'{name}: checksum mismatch')
        manifest[name]={'url':url,'sha256':digest,'bytes':path.stat().st_size}
    (root/'fetch_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return root,manifest

if __name__=='__main__':fetch()
