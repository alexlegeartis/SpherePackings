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
const target=[[0,0,0,1,0],[0,0,0,0,1],[0,0,0,1,-1],[1,-1,0,0,0],[1,0,-1,0,0],[0,1,-1,0,0]];
const sums=target.map(()=>Array(24).fill(0)),cov=target.map(()=>Array.from({length:24},()=>Array(24).fill(0))),cnt=target.map(()=>0),hs=target.map(()=>[]),P=x=>x;let nv=0;
const x=Array(24).fill(0);
function emit(){
nv++;let z=B.map(r=>dot(x,r)/4);if(z.some(v=>v%1))throw Error("nonintegral");
let y=mul(z),q=dot(z,y),k=key(y),t=Cs.get(k);
if(k!=0&&q!=t.normal)throw Error("norm nonminimal "+k+","+q+","+t.normal);
let idx=target.findIndex(a=>a.every((w,j)=>w==z[j]));if(idx>=0){cnt[idx]++;let H=x.map((v,j)=>60*v-2*B.reduce((s,r,i)=>s+r[j]*y[i],0));hs[idx].push(H);for(let i=0;i<24;i++){sums[idx][i]+=H[i];for(let j=0;j<=i;j++)cov[idx][i][j]+=H[i]*H[j];}}t.shell++; if(k!=0)t.act=240-q;
}
for(let i=0;i<24;i++)for(let j=0;j<i;j++)for(let a of [-4,4])for(let b of [-4,4]){x[i]=a;x[j]=b;emit();x[i]=x[j]=0}
for(let o of code.filter(c=>W(c)==8)){let I=[];for(let i=0;i<24;i++)if(o>>i&1)I.push(i);for(let s=0;s<256;s++)if(W(s)%2==0){I.forEach((i,b)=>x[i]=(s>>b&1)?-2:2);emit()}x.fill(0)}
for(let v of code){for(let i=0;i<24;i++) x[i]=(v>>i&1)?-1:1;for(let i=0;i<24;i++){let was=x[i];x[i]*=-3;emit();x[i]=was}}
const hist={};for(let [k,t] of Cs) {let w=t.normal+":"+(t.shell>0?"hit":"MISS");hist[w]=(hist[w]||0)+1};
if(nv!=196560)throw Error("shell");

const labels='063367777111177770000111111110000222233333333222211110011888811111111118811118888222288331111818118181111111111111111584885848485485855555151515112011021102112014114144114414114327777327723237715155151414114140303030311616116111623555523553232551541154114511451120110211021120111111111327777237732237722221415415151411514611161655665564114411403303003115555111155551188222288557777555577775588222288327777237723327723233232555555555555555523555532552332553277773277232377222211111541145115411451222222225454454554544545584885848485485841141441151515155555555532777732772323773663222236632322232232223222222322232232223223555532553223553277773277232377545445452222222212011021102112011415514141511514151514144141515135535225255253351144441116161155666655355352252552533555555555555555551111111111111111118888111111828228283362266226515151515555548888458854458882825757575782822828757575752828566556656556655614144141151551511145541145111154011212012110102154544545151415144151415135525325255352351154451145111154111101121201211010210112102121101201151515151441411454888845885445886226444444446226622644444444622651747011223224757355323325747244222224155115156154518184661455545123717255033221051585815114554155551450356522356451515868114455546175701122233502318681414555461545514181811514554555555142222304184846614555115145652325717244022321515115517172550332210235118184661455544515565111515757355323325747244222224747011223224511515455623265503514555585815114554154130222214555558181151455454151451518681414555462375701122233508681144555461511155717244022321522365455184846614555111515455623265503514555585815114554155171725503322102351181846614555445155651115157573553233257472442222247470112232245151115571724402232152236545518484661455511413022221455555818115145545415145151868141455546237570112223350868114455546158581511455415555145035652235645151551747011223224757355323325747244222224155115156154518184661455545123717255033221051848466145551151456523257172440223215151158681144555461757011222335023186814145554615455141818115145545555551422223041';

const J=[3,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23],H=hs[0],N=H.length;
const den=27360000n,ws=[20046672n,24149415n,24842880n,24971472n,26486400n,26610336n,28252648n,41926776n,43561224n];
if(labels.length!=2055||cnt.some(n=>n!=2055)||sums.flat().some(x=>x!=0))throw Error("contacts");
for(let i=0;i<19;i++){let sm=0n;for(let k=0;k<N;k++)sm+=ws[+labels[k]]*BigInt(H[k][J[i]]);if(sm!=0n)throw Error("mean");
for(let j=0;j<=i;j++){let s=0n;for(let k=0;k<N;k++)s+=ws[+labels[k]]*BigInt(2280*H[k][J[i]]*H[k][J[j]]);
let nn=0;for(let a=0;a<5;a++)for(let b=0;b<5;b++)nn+=B[a][J[i]]*B[b][J[j]]*m[a][b];if(s!=den*BigInt((+(i==j)*120-nn)*N*99840))throw Error("cov");}}
let inv=[];for(let i=1;i<1009;i++){let a=1,v=i,b=1007;while(b){if(b&1)a=a*v%1009;v=v*v%1009;b>>=1}inv[i]=a}
function rank(V,n){let bs=Array(n).fill(null),r=0;
for(let z0 of V){let z=Int32Array.from(z0,x=>(x%1009+1009)%1009);for(let i=0;i<n;i++){if(!z[i])continue;if(bs[i]){let t=z[i],u=bs[i];for(let j=i;j<n;j++)z[j]=(z[j]-t*u[j]+1018081)%1009;}else{let a=inv[z[i]];for(let j=i;j<n;j++)z[j]=z[j]*a%1009;bs[i]=z;r++;break}}if(r==n)break}return r}
function feature(w){let x=[];for(let i=0;i<19;i++){x.push(w[J[i]]);for(let j=0;j<=i;j++)x.push(w[J[i]]*w[J[j]]);}return x}
let rankQ=rank(H.map(feature),209),rankL=hs.map(r=>rank(r.map(w=>J.map(j=>w[j])),19));
if(rankQ!=209||rankL.some(r=>r!=19))throw Error("rigidity rank");

console.log(JSON.stringify({contacts:cnt,full_quadratic_linear_rank:rankQ,other_linear_ranks:rankL,weight_classes:ws.length}));
