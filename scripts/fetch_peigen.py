#!/usr/bin/env python3
"""Fetch and verify the exact, unmodified PEIGEN headers used by this artifact."""
import hashlib
import io
import json
import tarfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    metadata = json.loads((ROOT / 'third_party/peigen.json').read_text())
    target = ROOT / 'third_party/PEIGEN'
    expected = metadata['header_sha256']
    if all((target / name).is_file() and
           hashlib.sha256((target / name).read_bytes()).hexdigest() == digest
           for name, digest in expected.items()):
        (target / '.verified').write_text(metadata['commit'] + '\n')
        print('PEIGEN headers verified.')
        return
    url = metadata['repository'] + '/archive/' + metadata['commit'] + '.tar.gz'
    with urllib.request.urlopen(url, timeout=120) as response:
        archive = response.read()
    with tarfile.open(fileobj=io.BytesIO(archive), mode='r:gz') as tar:
        prefix = 'PEIGEN-' + metadata['commit'] + '/'
        # Only named, hash-checked headers are extracted; archive paths are not trusted.
        for name, digest in expected.items():
            member = tar.getmember(prefix + name)
            if not member.isfile() or '..' in Path(name).parts:
                raise ValueError('Invalid upstream file: ' + name)
            content = tar.extractfile(member).read()
            if hashlib.sha256(content).hexdigest() != digest:
                raise ValueError('PEIGEN checksum mismatch: ' + name)
            destination = target / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(content)
    (target / 'LICENSE').write_bytes((ROOT / 'third_party/PEIGEN-LICENSE').read_bytes())
    (target / '.verified').write_text(metadata['commit'] + '\n')
    print(f'Fetched and verified {len(expected)} PEIGEN headers at {metadata["commit"]}.')


if __name__ == '__main__':
    main()
