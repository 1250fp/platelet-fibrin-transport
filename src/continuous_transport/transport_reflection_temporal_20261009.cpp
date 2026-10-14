// Isolated research prototype: force-free, constant-D reflected increments.
// Not the general Behringer drift propagator. No production-code modifications.
// Continuous swept-sphere obstacles; raster-defined common lateral domain.
#define main historical_bridge_main
#include "transport_geometry_bridge.cpp"
#undef main

constexpr double TOI_TOL=2e-11, START_TOL=2e-10;
int reflection_limit=256;
bool diagnostic=false;
double norm(V x){return std::sqrt(dot(x,x));}
V unit(V x){double n=norm(x);if(n==0)throw std::runtime_error("zero normal");return mul(x,1/n);}
struct Hit {double t=2;std::vector<V> normals;int count=0;};
void merge(Hit& h,double t,V n){
 if(t<h.t-TOI_TOL){h.t=t;h.normals.clear();h.count=0;}
 if(std::abs(t-h.t)<=TOI_TOL){h.normals.push_back(n);++h.count;}
}
struct ReflectScene : Scene {
 using Scene::Scene;
 static V separation(const Primitive&q,V x){
  double vx=q.bx-q.ax,vy=q.by-q.ay,den=vx*vx+vy*vy;
  double u=den?std::clamp(((x[0]-q.ax)*vx+(x[1]-q.ay)*vy)/den,0.,1.):0.;
  double ex=x[0]-q.ax-u*vx,ey=x[1]-q.ay-u*vy,rho=std::hypot(ex,ey);
  double scale=rho>q.r?(rho-q.r)/rho:0.;
  return {ex*scale,ey*scale,x[2]-std::clamp(x[2],q.lo,q.hi)};
 }
 static double first(const Primitive&q,V a,V d,double r,double cap){
  const double rr=r*r;auto fa=eval(q,a,d);
  if(std::sqrt(fa.first)<r-START_TOL)throw std::runtime_error("penetrating start");
  // Squared distance to a convex primitive is convex along a segment.
  if(fa.second>=0)return 2;
  if(fa.first<=rr)return 0;
  auto fb=eval(q,add(a,mul(d,cap)),d);
  double len=norm(d)*cap;
  if(std::sqrt(fa.first)+std::sqrt(fb.first)>2*r+len+1e-12)return 2;
  double tm=cap;
  if(fb.second>0){
   double lo=0,hi=cap,dl=std::sqrt(fa.first),dh=std::sqrt(fb.first);
   for(int it=0;it<48;++it){
    if(dl+dh>2*r+(hi-lo)*norm(d)+1e-12)return 2;
    double m=(lo+hi)*.5;auto f=eval(q,add(a,mul(d,m)),d);
    if(f.second>0){hi=m;dh=std::sqrt(f.first);}else{lo=m;dl=std::sqrt(f.first);}
   }tm=(lo+hi)*.5;
  }
  if(eval(q,add(a,mul(d,tm)),d).first>rr)return 2;
  double lo=0,hi=tm;
  for(int it=0;it<48;++it){double m=(lo+hi)*.5;
   if(eval(q,add(a,mul(d,m)),d).first>rr)lo=m;else hi=m;}
  return lo; // approach from the free side; no artificial positional push
 }
 Hit hit(V a,V d,double r,double cap=1)const{
  Hit result;if(tree.empty())return result;
  std::array<int,128> stack;int top=0;stack[top++]=0;
  while(top){auto& n=tree[stack[--top]];double lim=std::min(cap,result.t+TOI_TOL);
   V b=add(a,mul(d,lim));if(!bbox(a,b,n.low,n.high,r))continue;
   if(n.left<0){for(int i=n.begin;i<n.end;++i){auto&q=p[order[i]];
    lim=std::min(cap,result.t+TOI_TOL);b=add(a,mul(d,lim));
    if(!bbox(a,b,q.low,q.high,r))continue;
    double t=first(q,a,d,r,lim);if(t>lim)continue;
    V normal=unit(separation(q,add(a,mul(d,t))));
    // A grazing contact does not change direction; do not create a zero loop.
    if(dot(normal,d)>=-1e-13*std::max(1.,norm(d)))continue;
    merge(result,t,normal);
   }}else{if(top+2>=128)throw std::runtime_error("BVH stack");stack[top++]=n.left;stack[top++]=n.right;}
  }return result;
 }
 double clearance_fast(V a,double cap=1)const{
  double best=cap*cap;if(tree.empty())return cap;std::array<int,128> st;int top=0;st[top++]=0;
  while(top){auto&n=tree[st[--top]];double b2=0;
   for(int j=0;j<3;++j){double q=std::max({n.low[j]-a[j],a[j]-n.high[j],0.});b2+=q*q;}
   if(b2>best)continue;
   if(n.left<0){for(int i=n.begin;i<n.end;++i)best=std::min(best,eval(p[order[i]],a,{0,0,0}).first);}
   else{if(top+2>=128)throw std::runtime_error("BVH stack");st[top++]=n.left;st[top++]=n.right;}
  }return std::sqrt(best);
 }
};

struct Region {
 int nx,ny;double xmin,ymin,dx,L;std::vector<uint8_t> inside;
 Region(const std::string& path){std::ifstream in(path);if(!(in>>nx>>ny>>xmin>>ymin>>dx>>L))throw std::runtime_error("domain header");
  if(nx<=0||ny<=0||dx<=0||L<=0)throw std::runtime_error("domain dimensions");
  inside.resize(size_t(nx)*ny);int b;for(auto&v:inside){if(!(in>>b))throw std::runtime_error("domain EOF");v=b!=0;}}
 bool valid(int x,int y)const{return x>=0&&y>=0&&x<nx&&y<ny&&inside[size_t(y)*nx+x];}
 bool point(V a)const{return valid(int(std::floor((a[0]-xmin)/dx)),int(std::floor((a[1]-ymin)/dx)));}
 bool closed_point(V a)const{
  double x=(a[0]-xmin)/dx,y=(a[1]-ymin)/dx;
  int ix=int(std::floor(x)),iy=int(std::floor(y));if(valid(ix,iy))return true;
  const double tol=2e-9;
  if(std::abs(x-std::round(x))<tol){ix=int(std::round(x));if(valid(ix-1,iy))return true;}
  if(std::abs(y-std::round(y))<tol){iy=int(std::round(y));if(valid(ix,iy-1))return true;
   if(std::abs(x-std::round(x))<tol&&valid(ix-1,iy-1))return true;}
  return false;
 }
 Hit hit(V a,V d,double cap=1)const{
  Hit h;double s[2]={(a[0]-xmin)/dx,(a[1]-ymin)/dx};
  double v[2]={d[0]/dx,d[1]/dx},next[2],inc[2];int q[2],step[2];
  for(int j=0;j<2;++j){
   double nearest=std::round(s[j]);if(std::abs(s[j]-nearest)<2e-12)s[j]=nearest;
   q[j]=int(std::floor(s[j]));step[j]=(v[j]>0)-(v[j]<0);
   if(v[j]<0&&s[j]==std::floor(s[j]))--q[j];
   inc[j]=step[j]?1/std::abs(v[j]):INFINITY;
   next[j]=step[j]>0?(q[j]+1-s[j])/v[j]:(step[j]<0?(q[j]-s[j])/v[j]:INFINITY);
  }
  if(!valid(q[0],q[1])){
   // Exact boundary starts directed outward are collisions at t=0.
   if(!closed_point(a))throw std::runtime_error("lateral invalid start");
   bool found=false;
   for(int j=0;j<2;++j)if(step[j]&&s[j]==std::round(s[j])){
    int c[2]={q[0],q[1]};c[j]-=step[j];
    if(valid(c[0],c[1])){V n={0,0,0};n[j]=-step[j];merge(h,0,n);found=true;}}
   if(!found){for(int j=0;j<2;++j)if(step[j]&&s[j]==std::round(s[j])){V n={0,0,0};n[j]=-step[j];merge(h,0,n);}}
   if(h.normals.empty())throw std::runtime_error("unresolved boundary start");return h;
  }
  for(int guard=0;guard<100000;++guard){
   double t=std::min(next[0],next[1]);if(t>cap+1e-13||!std::isfinite(t))return h;
   bool ax[2]={std::abs(next[0]-t)<=1e-13,std::abs(next[1]-t)<=1e-13};
   for(int j=0;j<2;++j)if(ax[j]){
    int c[2]={q[0],q[1]};c[j]+=step[j];
    if(!valid(c[0],c[1])){V n={0,0,0};n[j]=-step[j];merge(h,std::max(0.,t),n);}}
   if(h.normals.empty()&&ax[0]&&ax[1]&&!valid(q[0]+step[0],q[1]+step[1])){
    // Measure-zero re-entrant corner: both active coordinate normals.
    merge(h,t,{-double(step[0]),0,0});merge(h,t,{0,-double(step[1]),0});}
   if(!h.normals.empty())return h;
   for(int j=0;j<2;++j)if(ax[j]){q[j]+=step[j];next[j]+=inc[j];}
  }throw std::runtime_error("lateral DDA guard");
 }
};

V reflected(V d,const std::vector<V>& normals){
 std::vector<V> basis;
 for(auto n:normals){for(auto b:basis)n=sub(n,mul(b,dot(n,b)));
  if(norm(n)>1e-9)basis.push_back(unit(n));}
 if(basis.empty())throw std::runtime_error("empty reflection basis");
 V result=d;for(auto b:basis)result=sub(result,mul(b,2*dot(d,b)));return result;
}
struct ReflectionBudget : std::runtime_error {ReflectionBudget():std::runtime_error("local reflection budget"){} };
struct Result {V position,final_direction;int outcome=0,solid_hits=0,lateral_hits=0,multi_hits=0,subdivisions=0;
 double elapsed_fraction=1;};
Result advance_raw(const ReflectScene&s,const Region&reg,V a,V d,double r,bool reject=false){
 Result out;out.position=a;out.final_direction=d;V pos=a,rem=d;double elapsed=0,left=1;
 for(int bounces=0;bounces<reflection_limit;++bounces){
  double zt=2;int face=0;
  if(rem[2]<0&&pos[2]+rem[2]<=0){zt=-pos[2]/rem[2];face=1;}
  else if(rem[2]>0&&pos[2]+rem[2]>=reg.L){zt=(reg.L-pos[2])/rem[2];face=2;}
  double cap=std::min(1.,zt);Hit lat=reg.hit(pos,rem,cap);
  Hit solid=s.hit(pos,rem,r,std::min(cap,lat.t+TOI_TOL));double t=std::min(lat.t,solid.t);
  if(t>cap+TOI_TOL){
   out.position=add(pos,mul(rem,cap));
   if(zt<=1){out.outcome=face;out.elapsed_fraction=elapsed+left*zt;out.position[2]=face==1?0:reg.L;}
   return out;
  }
  bool sh=std::abs(solid.t-t)<=TOI_TOL,lh=std::abs(lat.t-t)<=TOI_TOL;
  if(diagnostic&&(bounces<4||bounces%10000==0)){std::cerr<<std::setprecision(17)<<"debug bounce "<<bounces<<" t "<<t<<" pos "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<" rem "<<rem[0]<<" "<<rem[1]<<" "<<rem[2]<<" normals";
   for(auto n:solid.normals)std::cerr<<" ["<<n[0]<<","<<n[1]<<","<<n[2]<<"]";std::cerr<<"\n";}
  out.solid_hits+=sh;out.lateral_hits+=lh;
  if(reject){out.position=a;return out;}
  std::vector<V> normals;if(sh)normals=solid.normals;
  if(lh)normals.insert(normals.end(),lat.normals.begin(),lat.normals.end());
  if(normals.size()>1)++out.multi_hits;
  t=std::clamp(t,0.,1.);pos=add(pos,mul(rem,t));elapsed+=left*t;left*=1-t;
  out.final_direction=reflected(out.final_direction,normals);
  rem=reflected(mul(rem,1-t),normals);
  if(norm(rem)<1e-15){out.position=pos;return out;}
 }
 throw ReflectionBudget();
}

// Resolve rare long billiard paths without changing the Gaussian increment.
// This is geometric subdivision, NOT Brownian time-step refinement/resampling.
// The second half continues in the direction resulting from the first half.
Result advance(const ReflectScene&s,const Region&reg,V a,V d,double r,bool reject=false,int depth=0){
 try{return advance_raw(s,reg,a,d,r,reject);}
 catch(const ReflectionBudget&){
  if(depth>=16)throw std::runtime_error("geometric subdivision depth: numerical failure, NOT censoring");
  Result first=advance(s,reg,a,mul(d,.5),r,reject,depth+1);
  if(first.outcome){first.elapsed_fraction*=.5;first.final_direction=mul(first.final_direction,2);++first.subdivisions;return first;}
  Result last=advance(s,reg,first.position,first.final_direction,r,reject,depth+1);
  last.elapsed_fraction=.5+.5*last.elapsed_fraction;last.final_direction=mul(last.final_direction,2);
  last.solid_hits+=first.solid_hits;last.lateral_hits+=first.lateral_hits;last.multi_hits+=first.multi_hits;
  last.subdivisions+=first.subdivisions+1;return last;
 }
}

int main(int argc,char**argv){try{
 if(argc<2)throw std::runtime_error("mode required");std::string mode=argv[1];
 if(mode=="contacts"){
  if(argc!=5)throw std::runtime_error("contacts scene queries output");ReflectScene scene(argv[2]);
  std::ifstream in(argv[3]);std::ofstream out(argv[4]);out<<std::setprecision(17);V a,d;double r;
  while(in>>a[0]>>a[1]>>a[2]>>d[0]>>d[1]>>d[2]>>r){auto h=scene.hit(a,d,r);V n=h.normals.empty()?V{0,0,0}:h.normals[0];
   out<<h.t<<" "<<n[0]<<" "<<n[1]<<" "<<n[2]<<" "<<h.count<<"\n";}return 0;
 }
 if(mode=="lateral"){
  if(argc!=5)throw std::runtime_error("lateral domain queries output");Region reg(argv[2]);std::ifstream in(argv[3]);std::ofstream out(argv[4]);out<<std::setprecision(17);V a,d;
  while(in>>a[0]>>a[1]>>a[2]>>d[0]>>d[1]>>d[2]){auto h=reg.hit(a,d);out<<h.t<<" "<<h.count<<"\n";}return 0;
 }
 if(mode=="steps"||mode=="diagnostic_step"){
  if(mode=="diagnostic_step"){reflection_limit=200000;diagnostic=true;}
  if(argc!=6)throw std::runtime_error("steps scene domain queries output");ReflectScene s(argv[2]);Region reg(argv[3]);
  std::ifstream in(argv[4]);std::ofstream out(argv[5]);out<<std::setprecision(17);V a,d;double r;
  while(in>>a[0]>>a[1]>>a[2]>>d[0]>>d[1]>>d[2]>>r){auto q=advance(s,reg,a,d,r);
   out<<q.position[0]<<" "<<q.position[1]<<" "<<q.position[2]<<" "<<q.outcome<<" "<<q.solid_hits<<" "<<q.lateral_hits<<" "<<q.elapsed_fraction<<" "
    <<q.final_direction[0]<<" "<<q.final_direction[1]<<" "<<q.final_direction[2]<<" "<<q.multi_hits<<" "<<s.clearance_fast(q.position,r+.01)<<" "<<q.subdivisions<<"\n";}
  return 0;
 }
 if(mode!="run"||(argc!=13&&argc!=14))throw std::runtime_error("run scene domain initials output N T dt seed radius D reject [particle_offset]");
 ReflectScene s(argv[2]);Region reg(argv[3]);std::ifstream in(argv[4]);std::ofstream out(argv[5]);
 int N=std::stoi(argv[6]);double T=std::stod(argv[7]),dt=std::stod(argv[8]);uint64_t seed=std::stoull(argv[9]);
 double r=std::stod(argv[10]),D=std::stod(argv[11]);bool rejection=std::stoi(argv[12])!=0;
 int offset=argc==14?std::stoi(argv[13]):0;if(offset<0)throw std::runtime_error("negative particle offset");
 if(dt<=0||T<=0||r<=0||D<=0)throw std::runtime_error("run parameters");
 if(std::abs(T/dt-std::round(T/dt))>1e-7)throw std::runtime_error("horizon not divisible");
 if(!in||!out)throw std::runtime_error("run IO");
 out<<std::setprecision(17)<<"particle,time_s,outcome,x0,y0,z0,xf,yf,zf,steps,solid_hits,lateral_hits,multi_hits,min_clearance_um,linear_exit_time_s,geometric_subdivisions\n";
 auto started=std::chrono::steady_clock::now();
 // Fixed fine Gaussian stream; coarser paths sum consecutive increments.
 // Only the finest Gaussian interval changes in this verification executable.
 // A legacy build with TEMPORAL_NOISE_DT=0.000015625 must replay old paths exactly.
#ifndef TEMPORAL_NOISE_DT
#define TEMPORAL_NOISE_DT 0.0000078125
#endif
 const double noise_dt=TEMPORAL_NOISE_DT;int stride=int(std::llround(dt/noise_dt));
 if(stride<1||std::abs(stride*noise_dt-dt)>1e-14)throw std::runtime_error("dt must be multiple of fine noise interval");
 for(int i=0;i<N;++i){V initial;if(!(in>>initial[0]>>initial[1]>>initial[2]))throw std::runtime_error("initial EOF");
  if(!reg.point(initial)||initial[2]<=0||initial[2]>=reg.L||s.clearance_fast(initial,r+.01)<r-START_TOL)throw std::runtime_error("invalid initial center");
  std::mt19937_64 rng(seed_hash(seed)^seed_hash(i+offset));std::normal_distribution<double> normal(0,std::sqrt(2*D*noise_dt));
  V pos=initial;uint64_t steps=0,solid=0,lateral=0,multi=0,subdivisions=0;int outcome=0;double time=T,linear=T,minclear=r+.01;
  for(uint64_t k=0;k<uint64_t(std::llround(T/dt));++k){
   V d={0,0,0};for(int c=0;c<stride;++c)for(int j=0;j<3;++j)d[j]+=normal(rng);
   Result a;
   try{a=advance(s,reg,pos,d,r,rejection);}catch(std::exception&e){std::cerr<<std::setprecision(17)<<"failed particle "<<i<<" step "<<k<<" pos "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<" delta "<<d[0]<<" "<<d[1]<<" "<<d[2]<<"\n";throw;}
   ++steps;solid+=a.solid_hits;lateral+=a.lateral_hits;multi+=a.multi_hits;subdivisions+=a.subdivisions;pos=a.position;
   if(a.solid_hits||a.lateral_hits||a.outcome){
    double clearance=s.clearance_fast(pos,r+.01);minclear=std::min(minclear,clearance);
    if(clearance<r-START_TOL||!reg.closed_point(pos))throw std::runtime_error("invalid accepted endpoint");}
   if(a.outcome){outcome=a.outcome;time=(k+1)*dt;linear=(k+a.elapsed_fraction)*dt;break;}
  }
  out<<i+offset<<","<<time<<","<<outcome<<","<<initial[0]<<","<<initial[1]<<","<<initial[2]<<","<<pos[0]<<","<<pos[1]<<","<<pos[2]<<","<<steps<<","<<solid<<","<<lateral<<","<<multi<<","<<minclear<<","<<linear<<","<<subdivisions<<"\n";
  if((i+1)%32==0||i+1==N){out.flush();std::cout<<"particles "<<i+1<<"/"<<N<<" seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<std::endl;}
 }
 return 0;
 }catch(std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}}
