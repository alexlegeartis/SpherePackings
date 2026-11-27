
const V=Array.from({length:256},(_,i)=>i),ones=n=>{let s=0;for(;n;n&=n-1)s++;return s;};
const b=(x,y)=>(ones((x&15)&(y>>4))+ones((y&15)&(x>>4)))%2;
const q=x=>ones((x&15)&(x>>4))%2;
const U=V.filter(x=>q(x)),D=new Int32Array(2**16);
for(let x of U)for(let y of U)if(q(x^y))D[x+256*y]=1;
for(let h=1;h<D.length;h*=2)for(let i=0;i<D.length;i+=2*h)for(let j=0;j<h;j++){
 let x=D[i+j],y=D[i+j+h];D[i+j]=x+y;D[i+j+h]=x-y;
}
const f=x=>((x&15)<<4)|(x>>4);
const key=(u,v)=>!u&&!v?"same":(!u||!v||u===v)?"double"+q(u||v):"distinct"+(q(u)+q(v)+q(u^v));
const values={};
for(let u of V)for(let v of V){const k=key(u,v),z=D[f(u)+256*f(v)];if(k in values&&values[k]!==z)throw Error("not invariant");values[k]=z;}


const vs=values;
const R=[];function rec(a,p){if(a.length===5){R.push(a);return;}for(let i=0;i<=p+1;i++)rec([...a,i],Math.max(p,i));};rec([0],0);
const lines=[];
for(const p of R){
 const n=Math.max(...p)+1,ind={};let c=0;for(let i=0;i<n;i++)for(let j=0;j<i;j++)ind[j+","+i]=c++;
 for(let m=0;m<(1<<c);m++){
   const q=(a,b)=>a===b?0:((m>>ind[Math.min(a,b)+","+Math.max(a,b)])&1);
   let E2=0,E3=0;
   for(let i=0;i<5;i++)for(let j=0;j<i;j++){
    E2+=p[i]===p[j]?120:q(p[i],p[j])?8:-8;
    for(let k=0;k<j;k++){
      const s=q(p[i],p[j])+q(p[i],p[k])+q(p[j],p[k]);
      const d=new Set([p[i],p[j],p[k]]).size;
      E3+=vs[d===1?"same":d===2?"double"+(s/2):"distinct"+s];
    }
   }
   lines.push({p,m,E2,E3});
 }
}

const dirs=Array.from({length:8},(_,i)=>i);
const E=[dirs.map(()=>1),...dirs.slice(1,7).map(i=>dirs.map(j=>j===i?2:j===i+1?-2:0)),dirs.map(j=>j>=6?2:0)];
function residue(y){
 let a=0;
 for(let i=0;i<8;i++){const d=E[i].reduce((s,x,j)=>s+x*y[j],0);if(d%4)throw Error("nonintegral");
 if((d/4)&1)a|=1<<i;}
 return a;
}
function choose(n,k){
 const O=[];
 function run(i,x){if(x.length===k){O.push(x);return;}for(let t=i;t<=n-k+x.length;t++)run(t+1,[...x,t]);}
 run(0,[]);return O;
}
const check=(x,m)=>{if(!x)throw Error(m);};
function listIntegral(m,q){
 const out=[];
 for(const a of choose(8,m))for(let s=0;s<(1<<m);s++){
  let w=dirs.map(()=>0);for(let i=0;i<m;i++)w[a[i]]=s>>i&1?q:-q;out.push(w);
 }
 return out;
}
const root=listIntegral(2,2),norm4=[...listIntegral(1,4),...listIntegral(4,2)];
for(let z=0;z<256;z++){
 const w=dirs.map(i=>z>>i&1?1:-1);
 if(w.reduce((a,b)=>a+b,0)%4===0)root.push(w);
 for(let i=0;i<8;i++){const v=w.map((x,j)=>i===j?3*x:x);if(v.reduce((a,b)=>a+b,0)%4===0)norm4.push(v);}
}
const cl2=new Set(root.map(residue)),cl4=new Set(norm4.map(residue));
check(root.length===240&&cl2.size===120&&!cl2.has(0),"E8 roots");
check(norm4.length===2160&&cl4.size===135&&!cl4.has(0)&&[...cl2].every(x=>!cl4.has(x)),"E8 singular nonzero cosets");
check(values.same===6720&&values.double0===-448&&values.double1===448&&values.distinct0===64&&values.distinct1===0&&values.distinct2===-64&&values.distinct3===128,"triangle eigenvalues");
const eig=lines.map(x=>12*x.E2+x.E3);
check(lines.length===1895,"partition support");
check(Math.min(...eig)===-320&&Math.max(...eig)===81600,"weighted Cayley spectral bound");
check((2n**32n)*320n===(2n**24n)*(81600n+320n),"Hoffman bound");
check(640n*2n**16n<3n**16n&& (2n**25n)<3n**16n,"no gain");
check(40n**3n<=16n**4n&&41n**3n>16n**4n,"Loomis–Whitney");
check((2n**32n)*3n**20n<3n**16n*4n**20n,"minimum 6 case");
check(2n**20n<3n**16n,"minimum at most 4 case");
check(9n*416n**47n>4n*423n**47n,"47D three-phase gate");
check(17n*2n**46n*13n**46n<3n**139n && 3n**139n<18n*2n**46n*13n**46n,"46D determinant gate");
console.log("PASS: E8 depth alphabet, all "+lines.length+" fourier profile states, and exact density arithmetic");

