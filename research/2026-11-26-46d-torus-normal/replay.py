from mixed import exact
from two import go
import subprocess,sys
for a in (6,8):
 for b in range(a,13,2):
  for g in range(a//2+1):
   if a*b-g*g<=68:
    exact(a,b,g);go(f'mixed_{a}_{b}_{g}.json')
subprocess.run([sys.executable,'-B','audit.py'],check=True)
