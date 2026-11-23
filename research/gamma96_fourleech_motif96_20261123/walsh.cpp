#include <bits/stdc++.h>
using namespace std;
vector<int> readW() {
  vector<int> xs; for(int x;cin>>x;){ assert(0<=x && x<(1<<24)); xs.push_back(x);}
  sort(xs.begin(),xs.end());
  assert(!xs.empty()&&xs.front()!=0&&adjacent_find(xs.begin(),xs.end())==xs.end());
  return xs;
}
void walsh(vector<int>& a) {
  assert(a.size()==(1<<24));
  for(int m=1;m<(1<<24);m*=2)
    for(int j=0;j<int(a.size());j+=2*m)
      for(int k=j;k<j+m;k++){
        int x=a[k],y=a[k+m]; a[k]=x+y;a[k+m]=x-y;
      }
}
int main(int argc,char**argv){
  assert(argc==2&&(string(argv[1])=="source"||string(argv[1])=="target"));
  vector<int> xs=readW();
  vector<int> a(1<<24);
  for(int x:xs)a[x]=1;
  walsh(a);
  if(string(argv[1])=="source"){
    const map<int,int> target{{-24,16678935},{4072,98280},{98280,1}};
    map<int,int> hist;for(int v:a)++hist[v];
    assert(xs.size()==98280&&a[0]==98280&&hist==target);
    cout<<"PASS Leech Q4 least eigenvalue -24, degree 98280"<<endl;
  }else{
    assert(xs.size()==5850&&a[0]==5850);
    assert(*min_element(a.begin(),a.end())==-390);
    int n=0;
    for(int& z:a){z=(z==-390);n+=z;}
    assert(n==20280&&a[0]==0);
    walsh(a);
    for(int x:xs)assert(a[x]==-1352);
    assert(a[0]==20280);
    cout<<"PASS target 20280-character embedding with inner product -1/15 on every allowed difference"<<endl;
  }
}
