// Local research driver. No production-code or mechanical-state writes.
// Coordinates are x,y,z in micrometres. Obstacles are convex xy capsules
// extruded across each physical slab; swept-sphere collision is tested against
// their continuous surfaces, not against endpoint occupancy alone.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using V = std::array<double,3>;
V add(V a,V b){for(int j=0;j<3;++j)a[j]+=b[j];return a;}
V mul(V a,double s){for(auto &v:a)v*=s;return a;}
V sub(V a,V b){return add(a,mul(b,-1));}
double dot(V a,V b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
struct Primitive { double ax,ay,bx,by,r,lo,hi; V low,high; };
struct Node {V low,high; int begin,end,left=-1,right=-1;};
struct Scene {
 std::vector<Primitive> p; std::vector<int> order; std::vector<Node> tree;
 Scene(const std::string& path){std::ifstream f(path); if(!f)throw std::runtime_error("scene open");
   Primitive q; while(f>>q.ax>>q.ay>>q.bx>>q.by>>q.r>>q.lo>>q.hi){
    q.low={std::min(q.ax,q.bx)-q.r,std::min(q.ay,q.by)-q.r,q.lo};
    q.high={std::max(q.ax,q.bx)+q.r,std::max(q.ay,q.by)+q.r,q.hi};p.push_back(q);}
   order.resize(p.size());std::iota(order.begin(),order.end(),0);if(!p.empty())build(0,p.size());}
 int build(int a,int b){Node n; n.begin=a;n.end=b;n.low={INFINITY,INFINITY,INFINITY};n.high={-INFINITY,-INFINITY,-INFINITY};
  for(int i=a;i<b;++i)for(int j=0;j<3;++j){n.low[j]=std::min(n.low[j],p[order[i]].low[j]);n.high[j]=std::max(n.high[j],p[order[i]].high[j]);}
  int id=tree.size();tree.push_back(n);if(b-a>6){int axis=0;for(int j=1;j<3;++j)if(n.high[j]-n.low[j]>n.high[axis]-n.low[axis])axis=j;
   int mid=(a+b)/2;std::nth_element(order.begin()+a,order.begin()+mid,order.begin()+b,[&](int u,int v){return p[u].low[axis]+p[u].high[axis]<p[v].low[axis]+p[v].high[axis];});
   int l=build(a,mid),r=build(mid,b);tree[id].left=l;tree[id].right=r;}return id;}
 // squared distance and half-gradient projected on trajectory direction
 static std::pair<double,double> eval(const Primitive&q,V x,V direction){
  double vx=q.bx-q.ax,vy=q.by-q.ay,den=vx*vx+vy*vy;
  double u=den?std::clamp(((x[0]-q.ax)*vx+(x[1]-q.ay)*vy)/den,0.,1.):0.;
  double ex=x[0]-(q.ax+u*vx),ey=x[1]-(q.ay+u*vy),rho=std::hypot(ex,ey);
  double fac=rho>q.r?(rho-q.r)/rho:0.;ex*=fac;ey*=fac;
  double ez=x[2]-std::clamp(x[2],q.lo,q.hi);
  return {ex*ex+ey*ey+ez*ez,ex*direction[0]+ey*direction[1]+ez*direction[2]};}
 static double minimum(const Primitive&q,V a,V b){V v=sub(b,a);auto fa=eval(q,a,v);if(fa.second>=0)return fa.first;
  auto fb=eval(q,b,v);if(fb.second<=0)return fb.first;
  double l=0.,h=1.;for(int i=0;i<46;++i){double m=(l+h)*.5;auto f=eval(q,add(a,mul(v,m)),v);if(f.second>0)h=m;else l=m;}
  return eval(q,add(a,mul(v,(l+h)*.5)),v).first;}
 static bool bbox(V a,V b,V low,V high,double r){double enter=0,leave=1;V v=sub(b,a);
  for(int j=0;j<3;++j){double l=low[j]-r,h=high[j]+r;if(v[j]==0){if(a[j]<l||a[j]>h)return false;}
   else {double t=(l-a[j])/v[j],s=(h-a[j])/v[j];if(t>s)std::swap(t,s);enter=std::max(enter,t);leave=std::min(leave,s);if(enter>leave)return false;}}
  return true;}
 bool collision(V a,V b,double r)const {if(tree.empty())return false;std::array<int,128> stack;int top=0;stack[top++]=0;
  double guard=r+1e-12;while(top){int id=stack[--top];const Node& n=tree[id];if(!bbox(a,b,n.low,n.high,guard))continue;
   if(n.left<0){for(int i=n.begin;i<n.end;++i){const auto&q=p[order[i]];if(bbox(a,b,q.low,q.high,guard)&&minimum(q,a,b)<=guard*guard)return true;}}
   else {if(top+2>=128)throw std::runtime_error("BVH stack");stack[top++]=n.left;stack[top++]=n.right;}}return false;}
 double clearance(V a)const {double best=INFINITY;for(auto&q:p)best=std::min(best,eval(q,a,{0,0,0}).first);return std::sqrt(best);}
};
struct Mask {
 int side,nz;double dx,hz,half;std::vector<uint8_t> bits;
 Mask(const std::string& path,int n,int m):side(n/8),nz(50*m),dx(47.6/n),hz(.1/m),half(2.975){
  std::ifstream f(path,std::ios::binary);if(!f)throw std::runtime_error("mask open");char magic[8];f.read(magic,8);
  if(std::string(magic,6)!=std::string("\x93NUMPY",6))throw std::runtime_error("npy signature");
  uint32_t h=0;unsigned char bytes[4];int c=magic[6]==1?2:4;f.read((char*)bytes,c);for(int i=0;i<c;++i)h|=uint32_t(bytes[i])<<(8*i);f.seekg(h,std::ios::cur);
  bits.resize(size_t(nz)*side*side/8);f.read((char*)bits.data(),bits.size());if(size_t(f.gcount())!=bits.size())throw std::runtime_error("mask length");}
 bool valid(int x,int y,int z)const{if(x<0||x>=side||y<0||y>=side||z<0||z>=nz)return false;
  size_t q=(size_t(z)*side+y)*(side/8)+x/8;return (bits[q]>>(7-x%8))&1;}
 bool point(V a)const{return valid(int(std::floor((a[0]+half)/dx)),int(std::floor((a[1]+half)/dx)),int(std::floor(a[2]/hz)));}
 int trace(V a,V b)const{
  V s={(a[0]+half)/dx,(a[1]+half)/dx,a[2]/hz},t={(b[0]+half)/dx,(b[1]+half)/dx,b[2]/hz};V d=sub(t,s);
  std::array<int,3> q,step;V next,inc;int stationary=0;
  for(int j=0;j<3;++j){q[j]=int(std::floor(s[j]));step[j]=(d[j]>0)-(d[j]<0);inc[j]=step[j]?1/std::abs(d[j]):INFINITY;
   if(d[j]==0&&s[j]==q[j])stationary|=1<<j;
   next[j]=step[j]>0?(q[j]+1-s[j])/d[j]:(step[j]<0?(q[j]-s[j])/d[j]:INFINITY);}
  auto cell=[&](std::array<int,3> c){if(c[0]<0||c[0]>=side||c[1]<0||c[1]>=side)return 4;if(c[2]<0)return 1;if(c[2]>=nz)return 2;return valid(c[0],c[1],c[2])?0:3;};
  auto covered=[&](std::array<int,3> c){int e=cell(c);for(int code=1;code<8;++code)if((stationary&code)==code){auto v=c;for(int j=0;j<3;++j)if(code&(1<<j))--v[j];e=std::max(e,cell(v));}return e;};
  int first=covered(q);if(first)return first;const double eps=32*std::numeric_limits<double>::epsilon();
  for(int guard=0;guard<100000;++guard){double hit=std::min({next[0],next[1],next[2]});if(hit>1.+eps||!std::isfinite(hit))return 0;
   int axes[3],k=0;for(int j=0;j<3;++j)if(std::abs(next[j]-hit)<=eps*std::max(1.,std::abs(hit)))axes[k++]=j;
   int status=0;
   for(int set=1;set<(1<<k);++set){auto v=q;for(int j=0;j<k;++j)if(set&(1<<j))v[axes[j]]+=step[axes[j]];
    status=std::max(status,covered(v));}if(status)return status;
   for(int j=0;j<k;++j){int ax=axes[j];q[ax]+=step[ax];next[ax]+=inc[ax];}}
  throw std::runtime_error("DDA guard");}
};
// 0 accept, 1 bottom, 2 top, 3 solid rejection, 4 lateral rejection.
int event(const Scene& scene,const Mask* mask,V a,V proposal,double radius,double half,double L){
 if(mask)return mask->trace(a,proposal);
 V d=sub(proposal,a);double zt=2,xt=2;int face=0;
 if(proposal[2]<0){zt=-a[2]/d[2];face=1;}else if(proposal[2]>=L){zt=(L-a[2])/d[2];face=2;}
 for(int j=0;j<2;++j){if(proposal[j]<-half)xt=std::min(xt,(-half-a[j])/d[j]);else if(proposal[j]>=half)xt=std::min(xt,(half-a[j])/d[j]);}
 double stop=std::min({1.,zt,xt});V b=add(a,mul(d,stop));
 if(scene.collision(a,b,radius))return 3;
 if(xt<=1&&xt<=zt+1e-12)return 4;
 if(zt<=1)return face;return 0;
}
uint64_t seed_hash(uint64_t x){x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);}
int main(int argc,char**argv){try{
 std::string mode=argv[1];Scene scene(argv[2]);
 if(mode=="query"){std::ifstream f(argv[3]);std::ofstream o(argv[4]);o<<std::setprecision(17);V a,b;double r;
  while(f>>a[0]>>a[1]>>a[2]>>b[0]>>b[1]>>b[2]>>r)o<<scene.collision(a,b,r)<<","<<scene.clearance(a)<<"\n";return 0;}
 if(mode=="events"){int n=std::stoi(argv[4]),m=std::stoi(argv[5]);Mask mask(argv[3],n,m);std::ifstream f(argv[6]);std::ofstream o(argv[7]);V a,b;double r;
  while(f>>a[0]>>a[1]>>a[2]>>b[0]>>b[1]>>b[2]>>r)o<<event(scene,&mask,a,b,r,2.975,5.)<<","<<event(scene,nullptr,a,b,r,2.975,5.)<<"\n";return 0;}
 if(mode=="control"){
  // Synthetic free interval / plane with a reflecting solid face. Same
  // full-proposal rejection as the real local experiment, at two paired dt.
  std::ofstream o(argv[3]);int count=std::stoi(argv[4]);double half=std::stod(argv[5]);uint64_t seed=std::stoull(argv[6]);
  o<<std::setprecision(17)<<"particle,method,time_s,outcome,x0,y0,z0,proposals,solid_reject,lateral_reject,xf,yf,zf\n";
  for(int i=0;i<count;++i){std::mt19937_64 rng(seed_hash(seed)^seed_hash(i));std::uniform_real_distribution<double> uniform(2.,3.);V a={0.,0.,uniform(rng)};
   std::normal_distribution<double> normal(0,std::sqrt(.000125));std::array<V,2> pos={a,a};std::array<double,2> times={100.,100.};std::array<int,2> outs{};std::array<uint64_t,2> pr{},so{},la{};
   auto take=[&](int k,V delta,double clock){if(outs[k])return;V q=add(pos[k],delta);int e=event(scene,nullptr,pos[k],q,.1,half,5.);++pr[k];if(e==0)pos[k]=q;else if(e<3){outs[k]=e;times[k]=clock;}else if(e==3)++so[k];else ++la[k];};
   for(int s=0;s<800000;++s){if(outs[0]&&outs[1])break;V d1,d2;for(int j=0;j<3;++j){d1[j]=normal(rng);d2[j]=normal(rng);}take(0,add(d1,d2),(s+1)*.000125);take(1,d1,s*.000125+.0000625);take(1,d2,(s+1)*.000125);}
   for(int k=0;k<2;++k)o<<i<<","<<k<<","<<times[k]<<","<<outs[k]<<",0,0,"<<a[2]<<","<<pr[k]<<","<<so[k]<<","<<la[k]<<","<<pos[k][0]<<","<<pos[k][1]<<","<<pos[k][2]<<"\n";}
  return 0;}
 // run scene initial.tsv prefix masks_dir radius D horizon count seed
 if(mode!="run")throw std::runtime_error("unknown mode");
 std::ifstream in(argv[3]);std::string prefix=argv[4],folder=argv[5];double radius=std::stod(argv[6]),D=std::stod(argv[7]),T=std::stod(argv[8]);int count=std::stoi(argv[9]);uint64_t seed=std::stoull(argv[10]);
 std::string rs=radius<.2?"0.10":"0.25";
 std::vector<Mask> masks;masks.emplace_back(folder+"/packed_n4096_m32_r"+rs+".npy",4096,32);masks.emplace_back(folder+"/packed_n8192_m16_r"+rs+".npy",8192,16);masks.emplace_back(folder+"/packed_n8192_m32_r"+rs+".npy",8192,32);
 std::ofstream out(prefix+".csv");out<<std::setprecision(17)<<"particle,method,time_s,outcome,x0,y0,z0,proposals,solid_reject,lateral_reject,xf,yf,zf\n";
 const double dt=.000125,fine=dt/2;const int ns=int(std::llround(T/dt));if(std::abs(ns*dt-T)>1e-9)throw std::runtime_error("horizon multiple");
 auto tick=std::chrono::steady_clock::now();
 for(int i=0;i<count;++i){V a;if(!(in>>a[0]>>a[1]>>a[2]))throw std::runtime_error("initial EOF");for(auto& mask:masks)if(!mask.point(a))throw std::runtime_error("invalid common start");if(scene.collision(a,a,radius))throw std::runtime_error("invalid continuous start");
  std::array<V,5> pos;pos.fill(a);std::array<double,5> time;time.fill(T);std::array<int,5> outcome{};std::array<uint64_t,5> proposed{},solid{},lateral{};
  std::mt19937_64 rng(seed_hash(seed)^seed_hash(i));std::normal_distribution<double> normal(0,std::sqrt(2*D*fine));
  auto take=[&](int k,V delta,double clock){if(outcome[k])return;V q=add(pos[k],delta);int e=event(scene,k<3?&masks[k]:nullptr,pos[k],q,radius,2.975,5.);++proposed[k];
   if(e==0)pos[k]=q;else if(e==1||e==2){outcome[k]=e;time[k]=std::min(clock,T);}else if(e==3)++solid[k];else ++lateral[k];};
  for(int s=0;s<ns;++s){if(std::all_of(outcome.begin(),outcome.end(),[](int e){return e>0;}))break;V d1,d2;for(int j=0;j<3;++j){d1[j]=normal(rng);d2[j]=normal(rng);}V total=add(d1,d2);
   for(int k=0;k<4;++k)take(k,total,(s+1)*dt);take(4,d1,s*dt+fine);take(4,d2,(s+1)*dt);}
  for(int k=0;k<5;++k)out<<i<<","<<k<<","<<time[k]<<","<<outcome[k]<<","<<a[0]<<","<<a[1]<<","<<a[2]<<","<<proposed[k]<<","<<solid[k]<<","<<lateral[k]<<","<<pos[k][0]<<","<<pos[k][1]<<","<<pos[k][2]<<"\n";
  if((i+1)%8==0){out.flush();std::cout<<"particles "<<i+1<<"/"<<count<<" elapsed "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-tick).count()<<" s"<<std::endl;}}
 return 0;
 }catch(std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}}
