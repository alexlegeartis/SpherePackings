'use strict';
const assert = require('node:assert/strict');
const Q = (n,d=1n) => [BigInt(n),BigInt(d)];
const add = ([n,d],[m,e]) => [n*e+m*d,d*e];
const neg = ([n,d]) => [-n,d];
const sub = (a,b) => add(a,neg(b));
const mul = ([n,d],[m,e]) => [n*m,d*e];
const div = ([n,d],[m,e]) => !m ? (()=>{throw Error('divide by zero')})() :
    m>0n ? [n*e,d*m] : [-n*e,-d*m];
const pow = ([n,d],p) => [n**BigInt(p),d**BigInt(p)];
const cmp = ([n,d],[m,e]) => n*e-m*d;
const zero=Q(0), one=Q(1);
const lt=(a,b)=>cmp(a,b)<0n, le=(a,b)=>cmp(a,b)<=0n;
const gt=(a,b)=>lt(b,a), eq=(a,b)=>cmp(a,b)===0n;
function str([n,d]) {const f=(a,b)=>b?f(b,a%b):a; const g=f(n,d); return (n/g)+'/'+(d/g)}

// Complete vertex check of the polar in sorted gaps (five gaps and t).
const d1=[1,2,3,4,5,6].map(j=>Q(j,6)), d2=[2,4,3,2,1,0].map(j=>Q(j,3));
let polar=zero, pivotCount=1;
function accept(v) {
    assert(v.every(x=>le(zero,x)));
    const dot=C=>C.reduce((a,c,j)=>add(a,mul(c,v[j])),zero);
    assert(le(dot(d1),one)&&le(dot(d2),one));
    const w=Array.from({length:6},(_,k)=>v.slice(0,5).reduce(
        (a,g,i)=>add(a,mul(g,Q((k>i?6:0)-(5-i),6))),zero));
    const s=add(mul(Q(1,2),w.reduce((a,x)=>add(a,pow(x,2)),zero)),pow(v[5],2));
    if(gt(s,polar))polar=s;
}
for(let i=0;i<6;i++){
    let u=div(one,d1[i]);
    if(gt(d2[i],zero)&&gt(mul(u,d2[i]),one))u=div(one,d2[i]);
    const v=new Array(6).fill(zero); v[i]=u;accept(v);pivotCount++;
    for(let j=i+1;j<6;j++){
        const det=sub(mul(d1[i],d2[j]),mul(d2[i],d1[j]));
        if(eq(det,zero))continue;
        const a=div(sub(d2[j],d1[j]),det), b=div(sub(d1[i],d2[i]),det);
        if(lt(a,zero)||lt(b,zero))continue;
        const u=new Array(6).fill(zero);u[i]=a;u[j]=b;accept(u);pivotCount++;
    }
}
assert(eq(polar,Q(3,2)));

// Raw A^(3) contact congruences. An Eisenstein pair denotes a+b*omega.
const zmod=(x,m)=>(x%m+m)%m;
const rot=([a,b])=>[-b,a-b];
const qn=([a,b])=>a*a-a*b+b*b;
const subE=([a,b],[c,d])=>[a-c,b-d];
function divTheta([a,b]) {
    // (a+b omega)/(1+2 omega) = ((2b-a)+(b-2a) omega)/3.
    assert((2*b-a)%3===0);assert((b-2*a)%3===0);
    return [(2*b-a)/3,(b-2*a)/3];
}
function isA(v) {
    const C=v.map(([a,b])=>zmod(a+b,3));
    return C.every(c=>c===C[0])&&[0,1].every(j=>v.reduce((s,w)=>s+w[j],0)%3===0);
}
const S=[];
for(let i=0;i<6;i++){
    S.push(Array.from({length:6},(_,j)=>[j===i?3:0,1]));
    S.push(Array.from({length:6},(_,j)=>[(j===i?3:0)-1,-1]));
    for(let j=i+1;j<6;j++) S.push(Array.from({length:6},(_,k)=>[k===i||k===j?-2:1,0]));
}
const all=S.flatMap(v=>[v,v.map(rot),v.map(rot).map(rot)]);
assert.equal(new Set(all.map(JSON.stringify)).size,81);
for(const v of all){
    assert(isA(v));assert.equal(v.reduce((s,w)=>s+qn(w),0),12);
    assert(isA(v.map((w,i)=>divTheta(subE(w,all[0][i])))));
}
// Each continuum estimate is checked by p < sqrt(q), using p^2 < q.
function belowSqrt(p,q){ return lt(p,zero)||lt(pow(p,2),q); }
const R2=Q(8,3), U2=Q(8,9);
assert(belowSqrt(sub(add(Q(2),U2),Q(5,2)),Q(8,27)));
for(const [h,r] of [[Q(5,2),Q(1,5)],[Q(21,8),Q(1,25)]]){
    assert(belowSqrt(sub(add(R2,pow(r,2)),h),mul(Q(4,3),pow(r,2))));
    assert(belowSqrt(sub(add(R2,U2),h),Q(32,27)));
}
assert(lt(pow(Q(3,5),2),Q(4,3)));
assert(lt(pow(Q(4,25),2),Q(4,3)));
// Last assertion gives sqrt(2)+2/25 < sqrt(23/10) (strict value checked).
assert(gt(sub(Q(23,10),add(Q(2),pow(Q(2,25),2))),zero));
assert(lt(Q(2),pow(div(sub(Q(23,10),add(Q(2),pow(Q(2,25),2))),Q(4,25)),2)));

// P10c from Best's cyclic quaternary code, automorphism, and literal clique.
const gray=[[0,0],[0,1],[1,1],[1,0]];
const words=new Set();
for(const b of [-1,1])for(const c of [-1,1])for(const d of [-1,1]){
    let w=[c-d,b,c,d,b+c].map(x=>zmod(x,4)).join('');
    for(let i=0;i<5;i++){words.add(w);w=w.slice(1)+w[0];}
}
assert.equal(words.size,40);
const rotWord=w=>w.slice(1)+w[0];
const genWord=w=>[zmod(-Number(w[0]),4),zmod(Number(w[1])+2,4),w[2],w[3],zmod(Number(w[4])+2,4)].join('');
for(const w of words){assert(words.has(rotWord(w)));assert(words.has(genWord(w)));}
const orbit=new Set(['00133']), queue=[...orbit];
for(let i=0;i<queue.length;i++)for(const f of [rotWord,genWord]){
    const w=f(queue[i]); if(!orbit.has(w)){orbit.add(w);queue.push(w)}
}
assert.equal(orbit.size,40);
const bits=w=>Array.from(w).flatMap(c=>gray[Number(c)]);
const C=new Set([...words].map(w=>bits(w).join('')));
const D2=(a,b)=>a.reduce((s,x,i)=>s+(x-b[i])**2,0);
for(const a of C)for(const b of C)if(a!==b)
    assert(D2(Array.from(a,Number),Array.from(b,Number))>=4);
const offsets=[
 [0,0,0,0,0,0,0,0,0,0],[-2,0,0,0,0,0,0,0,0,0],
 [-1,0,0,0,-1,0,0,0,-1,-1],[-1,0,0,0,-1,0,0,0,-1,1],
 [-1,0,0,-1,-1,0,0,-1,0,0],[-1,0,0,-1,-1,0,0,1,0,0],
 [-1,0,-1,0,-1,0,-1,0,0,0],[-1,0,-1,0,-1,0,1,0,0,0],
 [-1,-1,0,0,-1,-1,0,0,0,0],[-1,-1,0,0,-1,1,0,0,0,0]
];
for(let i=0;i<offsets.length;i++){
    assert(C.has(offsets[i].map((x,j)=>zmod(x+bits('00133')[j],2)).join('')));
    for(let j=0;j<i;j++)assert.equal(D2(offsets[i],offsets[j]),4);
}

const comp22Sq=mul(Q(3,4),pow(Q(11,12),22));
const full22=s=>div(Q(5,3456),pow(s,5));
assert(eq(full22(Q(1,3)),Q(45,128)));
assert(lt(Q(225,2048),comp22Sq)); // (15 sqrt(2)/64)^2
assert(lt(pow(Q(63,256),2),Q(225,2048)));
assert(lt(pow(mul(Q(1,10),full22(Q(2,9))),2),Q(225,2048)));
assert(eq(mul(Q(17,20),full22(Q(1,3))),Q(153,512)));
assert(lt(pow(Q(153,512),2),Q(225,2048)));
assert(lt(pow(mul(Q(1,2),full22(Q(5,17))),2),Q(225,2048)));

// The D4 cross-polytope is an 8-point translate window with norms 2 and 4.
const X4=[...Array(4).keys()].flatMap(i=>[-1,1].map(s=>Array.from({length:4},(_,j)=>(j===0?1:0)+(i===j?s:0))));
let roots=0, antip=0;
for(let i=0;i<X4.length;i++){assert(X4[i].reduce((a,b)=>a+b,0)%2===0);
 for(let j=0;j<i;j++){
    const n=D2(X4[i],X4[j]);
    if(n===2)roots++;else if(n===4)antip++;else throw Error('4-window norm');
 }}
assert.deepEqual([roots,antip],[24,4]);

const DD=(m,u)=>div(pow(div(m,Q(4)),8),mul(Q(54),pow(u,2)));
const cover=m=>div(sub(m,R2),Q(2)), period=m=>div(m,Q(18));
for(const m of [Q(24,7),Q(4)]) assert(lt(mul(Q(4,3),DD(m,cover(m))),Q(1,16)));
for(const m of [Q(10,3),Q(4)]) assert(lt(mul(Q(4,3),DD(m,cover(m))),Q(1,16)));
assert(lt(DD(Q(59,20),period(Q(59,20))),Q(1,16)));
assert(lt(DD(Q(121,40),cover(Q(121,40))),Q(1,16)));
assert(lt(DD(Q(4),cover(Q(4))),Q(1,16)));
assert(lt(mul(Q(7,8),DD(Q(3),Q(1,6))),Q(1,16)));
const scaleSq=m=>mul(Q(8,27),pow(div(m,Q(4)),8));
for(const m of [Q(59,20),Q(121,40)]){
    for(const [d,h] of [[2,Q(21,8)],[4,Q(23,10)]]){
        assert(gt(m,h));
        assert(lt(scaleSq(m),pow(div(sub(m,h),Q(d)),2)));
    }
}
console.log('PASS: deep polar',str(polar),pivotCount,'81 contacts; rounded thresholds; P10c 10-window; 22D and 16D exact comparisons');
