'use strict';
const assert=require('node:assert/strict');
const p2=[1n,2n,3n,4n,4n,4n,3n,2n,1n];
const add=(a,b)=>[a[0]*b[1]+b[0]*a[1],a[1]*b[1]];
const mul=(a,b)=>[a[0]*b[0],a[1]*b[1]];
const lt=(a,b)=>a[0]*b[1]<b[0]*a[1];
const eq=(a,b)=>a[0]*b[1]===b[0]*a[1];
function h6(j) {
 const m=BigInt(j);if(j<=8)return [p2[j]*3n**m,1n];
 const c=(j===10)?[2500n,9n]:[625n,4n];
 return mul([3n**m,2n**m],c);
}
const T={
  39:[21870n,1n],40:[6561n,1n],41:[4374n,1n],42:[2187n,1n],
  43:[2n**180n*3n**85n,137n**43n],
  44:[3n**44n*605n*110n**44n,2n**44n*25n*157n**44n],
  45:[2n**138n*3n**45n,23n**45n],
  46:[3n*27n**46n,26n**46n],
  47:[2n**189n*3n**47n,47n**47n],48:[1n,1n],
  49:[2n,1n],50:[3n,1n],51:[4n,1n],52:[4n,1n],53:[4n,1n],
};
for(const [d,U] of [[43,728n],[44,216n],[45,55n],[46,18n],[47,6n]])
 assert(lt(T[d],[U,1n]));
assert(eq(h6(10),[4100625n,256n]));
const cert=[];
for(let d=39;d<=48;d++){
 let low=null;
 for(let s=1;s<=5;s++){
  const j=48-d+s, v=mul(h6(j),[p2[s],1n]);assert(lt(T[d],v),[d,s].join(','));
  if(!low||lt(v,low))low=v;
 }
 cert.push(Number(low[0]/low[1]));
 if(d>=40)assert(!lt(h6(48-d),T[d]));
}
assert.deepEqual(cert,[32036,12013,13122,8748,4374,1944,648,216,54,12]);
for(let d=49;d<=53;d++)for(let s=d-48;s<=5;s++){
 const j=48-d+s, v=mul(h6(j),[p2[s],1n]);
 assert(j>0?lt(T[d],v):eq(T[d],v),[d,s].join(','));
}
console.log('PASS: all 50 positive-tail inequalities, the zero-tail comparisons, and all tight pure-vertical targets.');
