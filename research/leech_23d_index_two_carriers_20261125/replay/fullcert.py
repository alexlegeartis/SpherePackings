exec(open('cert.py').read().split('s=linprog(-ff')[0])
mom={}
for l in open('moments.txt'):
 a=list(map(int,l.split()));c,o,e,n=a[:4];mom[c,o,e]=[n]+[F(v,n) for v in a[4:]]
O=[0,1,5,6,7,9,11,23];eps=lambda i:(-1 if i==23 else 1) if i in O else 0
ijs=[(i,j) for i in range(24) for j in range(i+1)]
rows2=[];names2=[]
for name,row in zip(names,rows):
 if name=='pure-octad':opts=[(0,64,0)];x=y=0;i=j=0
 elif name=='pure-exterior':opts=[(0,0,64)];x=y=0;i=j=0
 else:
  i,j,u,v,o,e=name;c=A[A[C[j]][mul(C[i],-1)]][A[mul(4,u)][mul(33,v)]]
  opts=[(c,o,e)];x=F(3*u+j-i,3);y=F(3*v+j-i,3)
 for key in opts:
  m=mom[key];base=list(row[:7])+[row[-1]];mean=m[1:25];cov=m[25:];ext=cov+[x*h for h in mean]+[y*h for h in mean]+[int((j==k)-(i==k))*h for k in (1,2) for h in mean]
  rows2.append(base+ext);names2.append((name,m[0]))
target=list(-F(t) for t in f[:7])+[F(-f[-1])]+[64*(F(int(i==j and i not in [2,3]))-F(eps(i)*eps(j),8)) for i,j in ijs]+[0]*96
M=np.array(rows2,dtype=float).T;Y=np.array(target,dtype=float)
out=linprog(np.zeros(len(rows2)),A_eq=M,b_eq=Y,bounds=[(0,None)]*len(rows2),method='highs')
print('full moments',len(rows2),'classes, target dimension',len(target),'LP',out.message)
if out.success:
 I=np.nonzero(out.x>1e-8)[0];E=sp.Matrix([[sp.Rational(rows2[k][j]) for k in I] for j in range(len(target))]);b=sp.Matrix([sp.Rational(v) for v in target]);w=E.gauss_jordan_solve(b)[0];assert E*w==b and all(t>0 for t in w);print('exact total weights grouped by shell:',list(zip((names2[k] for k in I),w)))
