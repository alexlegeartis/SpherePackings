"""Decode the ASCII chunks of a deterministic gzip-compressed JSON receipt."""
import json,gzip,base64
from pathlib import Path
p=Path(__file__).resolve().parent
paths=sorted(p.glob('receipt.[0-9][0-9]'))
assert paths
payload=b''.join(b''.join(a.read_bytes().split()) for a in paths)
read=json.loads(gzip.decompress(base64.b85decode(payload)))
