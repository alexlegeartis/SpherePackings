const m=[[21,9,9,12,12],[9,29,13,20,16],[9,13,29,16,20],[12,20,16,32,16],[12,16,20,16,32]];
const rowSpec=[[0,0,0,0,2,0,0,0,2,0,0,0],[2,0,0,0,0,2,0,0,0,0,0,0],[0,2,0,0,0,0,0,0,0,2,0,0],[-1,-1,-1,-1,-1,-1,1,1,0,0,0,0],[-1,-1,1,1,0,0,0,0,-1,-1,1,1]];
const B=rowSpec.map(r=>r.concat(Array(12).fill(0)));
const sum=a=>a.reduce((s,x)=>s+x,0);
const dot=(a,b)=>sum(a.map((x,i)=>x*b[i]));
const W=x=>{x>>>=0;let n=0;while(x){x&=x-1;n++}return n};
const mul=z=>m.map(r=>dot(r,z));
const key=y=>y.reduce((k,c)=>k*60+((c%60+60)%60),0);
const C0=[];for(let u=0;u<4096;u++){let v=0;for(let i=0;i<12;i++)if(u>>i&1)v^=0xae3<<i; C0.push(v|((W(v)%2)<<23))}
const oct=C0.filter(x=>W(x)==8),o=oct[0],p=oct.find(x=>W(x&o)==4);
const coords=[];for(const [a,b] of [[1,1],[1,0],[0,1],[0,0]])for(let k=0;k<24;k++)if(((o>>k)&1)==a&&((p>>k)&1)==b)coords.push(k);
const code=C0.map(c=>coords.reduce((s,k,i)=>s|(((c>>k)&1)<<i),0));
if(code.filter(x=>W(x)==8).length!=759||code.some(c=>c&&W(c)<8))throw Error("Golay");
const K=B.map(r=>B.map(s=>dot(r,s)/2));if(K.some((r,i)=>mul(r).some((w,j)=>w!=(i==j?60:0))))throw Error("wrong K");
const Cs=new Map;
for(let a=-4;a<=4;a++)for(let b=-4;b<=4;b++)for(let c=-4;c<=4;c++)for(let d=-4;d<=4;d++)for(let e=-4;e<=4;e++){
const z=[a,b,c,d,e],y=mul(z),k=key(y),q=dot(y,z);let t=Cs.get(k);if(!t||q<t.normal)Cs.set(k,{y:y.map(c=>(c%60+60)%60),z,normal:q,shell:0,act:null});
}
if(Cs.size!=300)throw Error("class count "+Cs.size);

let echelon=Array(24).fill(0);for(let val of code){let w=val;for(let i=0;i<24;i++)if((w>>i)&1){if(echelon[i])w^=echelon[i];else{echelon[i]=w;break}}}
for(let i=23;i>=0;i--)if(echelon[i])for(let j=0;j<i;j++)if((echelon[j]>>i)&1)echelon[j]^=echelon[i];
const piv=Array.from({length:24},(_,i)=>i).filter(i=>echelon[i]),free=Array.from({length:24},(_,i)=>i).filter(i=>!echelon[i]),last=free[11],u=Array.from({length:24},(_,i)=>i==0?-3:1),p0=piv.find(i=>u[i]==1);
const bas=[u,...piv.filter(i=>i!=p0).map(i=>u.map((_,j)=>2*((echelon[i]>>j)&1))),...free.slice(0,11).map(i=>u.map((_,j)=>4*((i==j)-(last==j)))),u.map((_,j)=>8*(last==j))];
const fastHints=bas.map(r=>r.map((x,i)=>[i,x]).filter(e=>e[1]));
function repr(x){
let w=x.slice(),mask=0,bit=0;
function sub(a){if(a%1)throw Error("nonintegral basis");if(a&1)mask|=1<<bit;if(a)for(let [i,c] of fastHints[bit])w[i]-=a*c;bit++}
sub(w[p0]);for(let i of piv)if(i!=p0)sub(w[i]/2);for(let i of free)if(i!=last)sub(w[i]/4);sub(w[last]/8);if(w.some(x=>x!=0))throw Error("bad basis");return mask;
}
function put(C,w){for(let i=23;i>=0;i--){if(w>>i&1){if(C[i])w^=C[i];else {C[i]=w;return true}}}return false}
const samples=[]; const ext=[];
let nv=0;const x=Array(24).fill(0);
function emit(){
nv++;let z=B.map(r=>dot(x,r)/4);if(z.some(v=>v%1))throw Error("nonintegral");
let y=mul(z),q=dot(z,y),k=key(y),t=Cs.get(k);
if(k!=0&&q!=t.normal)throw Error("norm nonminimal "+k+","+q+","+t.normal);

if(k!=0&&z.every((a,i)=>a==t.z[i])){if(!t.rs)t.rs=[];t.rs.push(repr(x))}
if(nv%43==0)samples.push(x.slice());t.shell++; if(k!=0)t.act=240-q;
}
for(let i=0;i<24;i++)for(let j=0;j<i;j++)for(let a of [-4,4])for(let b of [-4,4]){x[i]=a;x[j]=b;emit();x[i]=x[j]=0}
for(let o of code.filter(c=>W(c)==8)){let I=[];for(let i=0;i<24;i++)if(o>>i&1)I.push(i);for(let s=0;s<256;s++)if(W(s)%2==0){I.forEach((i,b)=>x[i]=(s>>b&1)?-2:2);emit()}x.fill(0)}
for(let v of code){for(let i=0;i<24;i++) x[i]=(v>>i&1)?-1:1;for(let i=0;i<24;i++){let was=x[i];x[i]*=-3;emit();x[i]=was}}
const hist={};for(let [k,t] of Cs) {let w=t.normal+":"+(t.shell>0?"hit":"MISS");hist[w]=(hist[w]||0)+1};
if(nv!=196560)throw Error("shell");

let stats=[];for(const [k,t] of Cs){let R=Array(24).fill(0),rs=t.rs||[];for(let a of rs)put(R,a^rs[0]);stats.push([k,t.normal,rs.length,R.filter(x=>x).length]);t.ranks=R.filter(x=>x);t.first=rs[0]??null;delete t.rs}

const dict=new Map();for(let x of samples){let z=B.map(r=>dot(x,r)/4),k=z.join(",");if(!dict.has(k))dict.set(k,[]);dict.get(k).push(x)}
const cert={};const filedata=[];
for(let [k,t] of Cs)if(!t.shell){let rank=Array(24).fill(0),first=null,n=0,v=[];
for(let [s,A] of dict){let z=s.split(",").map(Number),BB=dict.get(t.z.map((c,i)=>c-z[i]).join(","));if(!BB)continue;
for(let x of A)for(let b of BB){if(dot(x,b)!=-8)continue;let w=x.map((v,i)=>v+b[i]),c=repr(w); n++;if(first==null){first=c;v.push(w)}else if(put(rank,c^first))v.push(w);}
if(v.length==20)break;}
cert[k]={r:rank.filter(x=>x),first,n,v};filedata.push([k,t.normal,n,rank.filter(x=>x).length]);}

if(!Object.values(cert).every(t=>t.r.length==19)||Object.keys(cert).length!=18)throw Error("norm6 ranks");
const rows=[...Cs.values()],kv=B.map(b=>repr(b.map(x=>2*x)));
const pivL=rows.find(t=>t.ranks.length==19).ranks.map(x=>31-Math.clz32(x));
function redBasis(ww){let a=Array(19).fill(0);for(let z of ww){let w=pivL.reduce((s,c,i)=>s|((z>>c&1)<<i),0);for(let i=18;i>=0;i--)if(w>>i&1){if(a[i])w^=a[i];else{a[i]=w;break}}}for(let i=0;i<19;i++)for(let j=i+1;j<19;j++)if(a[j]>>i&1)a[j]^=a[i];return a}
const groupMap=new Map;for(let [i,t] of rows.entries()){if(!t.shell||!t.normal||t.ranks.length==19)continue;let r=redBasis(t.ranks),s=r.join(",");if(!groupMap.has(s))groupMap.set(s,{r,idx:[]});groupMap.get(s).idx.push(i)}
const groups=[...groupMap.values()];if(groups.length!=19)throw Error("groups");

function masksFor(groups){const sig=new Uint32Array(1<<19);for(let [j,t] of groups.entries()){let a=t.r,free=a.map((v,i)=>v?-1:i).filter(i=>i>=0);for(let s=0;s<1<<free.length;s++){let f=0;for(let i=0;i<free.length;i++)if(s>>i&1)f|=1<<free[i];for(let i=0;i<19;i++)if(a[i]){let z=(a[i]&f);z^=z>>>16;z^=z>>>8;z^=z>>>4;if((0x6996>>>(z&15))&1)f|=1<<i;}sig[f]|=1<<j;}}return sig;}

const sigsByChar=masksFor(groups),counts=new Map;for(let f=1;f<sigsByChar.length;f++){let s=sigsByChar[f];counts.set(s,(counts.get(s)||0)+1)}
const sigs=[...counts].map(([sig,n])=>({sig,n}));if(sigs.length!=31)throw Error("sigs");


const byKey=new Map(rows.map((t,i)=>[key(t.y),i]));
const ix=(u,v)=>byKey.get(key(rows[u].y.map((a,k)=>a-rows[v].y[k])));
const sep=rows.map(t=>t.normal==0?240:(t.act||360-t.normal));
const L=rows.map((t,i)=>rows.map((s,j)=>ix(i,j)));
const zero=byKey.get(0);
function calc(sig,th,goal,lim=5e6){let level=sep.map(s=>s>=th?2:0);for(let [j,t] of groups.entries())if((sig>>j)&1)for(let i of t.idx)level[i]=Math.max(level[i],1);
let D=L.map(r=>r.map(k=>level[k])),steps=0,best=goal-1,Wit=null,stop=false,finished=true;
function rec(P,S,score){if(++steps>lim){finished=false;stop=true;return}
if(score>best){best=score;Wit=S}if(score>=goal){stop=true;return}
let sum=0,U=P.slice(),oo=[],bb=[];while(U.length){let rem=[],col=[],c=0;for(let e of U){if(col.every(j=>!edge(e,j))){c=Math.max(c,e%2+1);col.push(e)}else rem.push(e)}sum+=c;for(let e of col){oo.push(e);bb.push(sum)}U=rem}
for(let k=oo.length-1;k>=0;k--){if(score+bb[k]<=best)break;let e=oo[k];rec(oo.slice(0,k).filter(i=>edge(e,i)),S.concat(e),score+1+(e&1));if(stop)return}}
function edge(a,b){let i=a>>1,j=b>>1;if(i==j)return false;let l=D[i][j];return l==2||l==1&&!(a&1)&&!(b&1)}
for(let w of [1,0]){let p=[];for(let i=0;i<rows.length;i++){for(let b of [0,1])if(edge(2*zero+w,2*i+b))p.push(2*i+b)}rec(p,[2*zero+w],1+w);if(stop)break}
return {best,steps,finished,Wit};}

const par=x=>{x^=x>>>16;x^=x>>>8;x^=x>>>4;return (0x6996>>>(x&15))&1};
const signed=L.map((r,i)=>r.map((k,j)=>{let z=rows[k].z,dt=z.map((x,h)=>rows[i].z[h]-rows[j].z[h]-x),t=m.map(v=>v.reduce((s,c,h)=>s+c*dt[h],0)/60);
if(t.some(x=>x%1))throw Error("sign");return t.reduce((v,c,h)=>v^((c&1)?kv[h]:0),rows[k].first||0)}));
const ids=Array.from({length:600},(_,i)=>i);
function light(th,goal,sig,f){let c=pivL.reduce((s,j,i)=>s|(((f>>i)&1)<<j),0),spl=new Set;
for(let [j,g] of groups.entries())if(sig>>j&1)for(let i of g.idx)spl.add(i);
function edge(a,b){if(a==b)return false;let i=a>>1,j=b>>1,k=L[i][j];return sep[k]>=th||spl.has(k)&&((a^b^par(signed[i][j]&c))&1)}
let best=[],steps=0;
function rec(P,S){steps++;if(S.length>best.length)best=S;if(S.length>=goal)return true;
let U=P.slice(),ord=[],bounds=[],n=0;while(U.length){n++;let rem=[],g=[];for(let x of U){if(g.every(y=>!edge(x,y))){g.push(x);ord.push(x);bounds.push(n)}else rem.push(x)}U=rem}
for(let k=ord.length-1;k>=0;k--){if(S.length+bounds[k]<=best.length)break;let x=ord[k];if(rec(ord.slice(0,k).filter(y=>edge(y,x)),S.concat(x)))return true}return false}
let root=2*zero;rec(ids.filter(x=>edge(root,x)),[root]);return [best.length,steps];}
function goal(t){let v=0;while(BigInt(v+1)**2n*BigInt(t)**19n<=144n*208n**19n)v++;return v+1}
let table=[];for(let th of [...new Set(sep)].sort((a,b)=>b-a)){let n=goal(th);if(n>600)continue;for(let s of sigs){let q=calc(s.sig,th,n,2e5);table.push([th,n,s.sig,q.best>=n,q.finished,q.steps]);}}
const ex=table.filter(x=>x[3]);if(table.some(x=>!x[4])||ex.length!=8)throw Error("bounds");
let totals=[];for(let [t,g,s] of ex){let mx=0,n=0;for(let f=1;f<sigsByChar.length;f++)if(sigsByChar[f]==s){mx=Math.max(mx,light(t,g,s,f)[0]);n++}if(mx>=g)throw Error("open");totals.push([t,g,n,mx]);}
console.log(JSON.stringify({n_classes:rows.length,shell:nv,norm6_full_rank:Object.keys(cert).length,characters:sigsByChar.length-1,signatures:sigs.length,upper_tests:table.length,exact_signed:totals}));

const six=Array.from({length:6},(_,i)=>Array.from({length:5},(_,j)=>+(i===j+1)));
const tside=[0,1,1,1,0,0], cross=[];
for(let i=0;i<6;i++) for(let j=0;j<i;j++) {
 const z=six[i].map((v,b)=>v-six[j][b]),k=byKey.get(key(mul(z)));
 if(tSideSafe(i,j)!=(sep[k]===208)) throw Error("wrong contact triangles");
 if(!tSideSafe(i,j)) cross.push(sep[k]);
}
function tSideSafe(i,j){return tside[i]===tside[j]}
if(!cross.every(x=>x===211||x===219)||Math.min(...cross)!==211)throw Error("wrong slide margin");
if(!rows.every(x=>x.normal<=180))throw Error("normal search threshold");
