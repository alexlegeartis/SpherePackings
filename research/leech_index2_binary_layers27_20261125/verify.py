"""Run from a checkout. Requires Python >=3.10 and a C++17 compiler."""
from pathlib import Path
from tempfile import TemporaryDirectory
from subprocess import run
import sys
here=Path(__file__).resolve().parent
with TemporaryDirectory() as td:
 def call(*a):run(a,cwd=td,check=True)
 for src in ('marked','dodecad','tri'):call('c++','-std=c++17','-O3',str(here/(src+'.cpp')),'-o',src)
 for m in ['3','15']:
  call('./marked',m);call('./tri','far'+m+'.bin')
 call('./dodecad');call('./tri','fardod.bin')
 call(sys.executable,str(here/'tail_exact.py'))
print('OK: all three far sets are sum-free; squared vertical covolume >= 400.')
