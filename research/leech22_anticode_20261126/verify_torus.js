const forms=(function main(){
let out=[];
function solve(a,b,c){
 let D=a*c-b*b; const norm=(x,y)=>c*x*x-2*b*x*y+a*y*y;
 let E=new Map();
 for(let x=-2*(a+b+c);x<=2*(a+b+c);x++) for(let y=-2*(a+b+c);y<=2*(a+b+c);y++){
   let s=(p)=>((p%D)+D)%D; let key=s(c*x-b*y)+D*s(a*y-b*x);
   if(!E.has(key)||E.get(key)>norm(x,y))E.set(key,norm(x,y))
 }
 if(E.size!==D)throw Error('group order');
 let labels=[...E.keys()].sort((a,b)=>a-b),F=labels.map(l=>E.get(l)), lev=[...new Set(F)].filter(q=>0<q&&9*q<8*D);
 let s=x=>((x%D)+D)%D,sub=(x,y)=>s(x%D-y%D)+D*s(Math.floor(x/D)-Math.floor(y/D));
 let best=0,win=null,levels=[];
 for(let q of lev){
  let adj=labels.map(x=>labels.reduce((m,y,j)=>x!==y&&E.get(sub(x,y))<=q?m|1n<<BigInt(j):m,0n));
  let max=1,idx=[];
  function dfs(P,ls){
    if(!P){if(ls.length>max){max=ls.length;idx=ls;}return}
    while(P){
      let v=P&-P,j=v.toString(2).length-1;P^=v;
      dfs(P&adj[j],[...ls,j]);
    }
  }
  dfs(adj[0],[0]);
  let mu=4-q/D, delta=max/Math.sqrt(D)*(mu/4)**11;
  let lhs=12n*BigInt(max)**2n*(3n*BigInt(4*D-q))**22n,rhs=9n*BigInt(D)*(11n*BigInt(D))**22n;
  levels.push([q,max,Math.sign(Number(lhs-rhs))]);
  if(delta>best){best=delta;win=[q,max,idx]}
 }
 return {a,b,c,D,best,win,levels};
}
for(let a of [4,6])for(let b=0;b<=a/2;b++)for(let c=a;a*c-b*b<48;c+=2)out.push(solve(a,b,c));
return out;})();
if(forms.length!==21)throw Error('forms');
const equal=[];
for(const f of forms)for(const [q,s,relation] of f.levels){
 if(relation>0)throw Error('counterexample');
 if(relation===0)equal.push([f.a,f.b,f.c,q,s]);
}
if(JSON.stringify(equal)!=='[[4,2,4,4,3]]')throw Error('equality cases');
console.log(forms.length,equal);
