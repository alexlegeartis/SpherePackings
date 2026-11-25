#include <vector>
#include <fstream>
#include <cassert>
#include <cstdint>
#include <iostream>
using namespace std;int main(int ac,char**av){const int N=1<<24;vector<int>w(N);ifstream f(av[1],ios::binary);int a,m=0;while(f.read((char*)&a,4)){assert(0<a&&a<N&&!w[a]);w[a]=1;m++;}for(int s=1;s<N;s*=2)for(int i=0;i<N;i+=2*s)for(int j=i;j<i+s;j++){int x=w[j],y=w[j+s];w[j]=x+y;w[j+s]=x-y;}__int128 t=0;for(int x:w)t+=(__int128)x*x*x;assert(w[0]==m&&t==0);cout<<av[1]<<" "<<m<<" exact Walsh cubic moment zero\n";}
