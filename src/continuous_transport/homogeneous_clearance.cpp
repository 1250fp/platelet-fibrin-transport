// Geometry queries only. Reuse the existing continuous primitive definition.
// This executable does not integrate trajectories or alter mechanical data.
#define main historical_bridge_main
#include "transport_geometry_bridge.cpp"
#undef main

struct DistanceScene : Scene {
 using Scene::Scene;
 static double box_distance2(V x,const Node& n){
  double d=0;for(int j=0;j<3;++j){double e=std::max({n.low[j]-x[j],x[j]-n.high[j],0.});d+=e*e;}return d;
 }
 double capped_clearance(V x,double cap)const {
  double best=cap*cap;if(tree.empty())return cap;
  std::array<int,128> stack;int top=0;stack[top++]=0;
  while(top){const auto& n=tree[stack[--top]];
   if(box_distance2(x,n)>best)continue;
   if(n.left<0){for(int i=n.begin;i<n.end;++i){best=std::min(best,eval(p[order[i]],x,{0,0,0}).first);if(best==0)return 0;}}
   else {if(top+2>=128)throw std::runtime_error("BVH stack");stack[top++]=n.left;stack[top++]=n.right;}
  }return std::sqrt(best);
 }
};
int main(int argc,char**argv){try{
 if(argc<2)throw std::runtime_error("points scene input_xyz.tsv output.tsv cap | grid scene output.bin side nz half L cap");
 const std::string mode=argv[1];DistanceScene scene(argv[2]);
 if(mode=="points" || mode=="check"){
  if(argc!=6)throw std::runtime_error("points/check argument count");
  std::ifstream in(argv[3]);std::ofstream out(argv[4]);double cap=std::stod(argv[5]);
  if(!in||!out)throw std::runtime_error("query IO");out<<std::setprecision(17);V x;
  while(in>>x[0]>>x[1]>>x[2]){
   const double d=scene.capped_clearance(x,cap);out<<d;
   if(mode=="check")out<<" "<<std::min(cap,scene.clearance(x));out<<"\n";
  }return 0;
 }
 if(mode=="grid"){
  if(argc!=9)throw std::runtime_error("grid argument count");
  int side=std::stoi(argv[4]),nz=std::stoi(argv[5]);double half=std::stod(argv[6]),L=std::stod(argv[7]),cap=std::stod(argv[8]);
  if(side<=0||nz<=0||half<=0||L<=0||cap<=0)throw std::runtime_error("grid dimensions");
  std::ofstream out(argv[3],std::ios::binary);if(!out)throw std::runtime_error("grid output");
  for(int z=0;z<nz;++z)for(int y=0;y<side;++y)for(int x=0;x<side;++x){
   V q={-half+(x+.5)*2*half/side,-half+(y+.5)*2*half/side,(z+.5)*L/nz};
   double d=scene.capped_clearance(q,cap);out.write(reinterpret_cast<const char*>(&d),sizeof(d));
  }if(!out)throw std::runtime_error("grid write");return 0;
 }
 throw std::runtime_error("unknown mode");
}catch(std::exception& e){std::cerr<<e.what()<<"\n";return 1;}}
