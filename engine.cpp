#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>
#include <cstdio>
#include <cassert>
#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#define API extern "C" EMSCRIPTEN_KEEPALIVE
#else
#define API extern "C"
#endif
struct V { double x=0,y=0; V operator+(V b)const{return{x+b.x,y+b.y};} V operator-(V b)const{return{x-b.x,y-b.y};} V operator*(double a)const{return{x*a,y*a};} };
double dot(V a,V b){return a.x*b.x+a.y*b.y;} double len(V a){return std::sqrt(dot(a,a));}
struct Body{V p,v;double radius,inv;};
struct Joint{int a,b;double rest;};
struct Contact{int a,b;V n;double bias,rebound,normal=0,tangent=0;};
static std::vector<Body> bodies;
static std::vector<Joint> joints;
static std::vector<Contact> contacts;
static constexpr double H=1.0/120,WIDTH=12,HEIGHT=8;
static int iterations=16,steps=0,scene=0,selected=-1;
static double gravity=9.81,restitution=.25,friction=.4,peakError=0;
static uint32_t seed=1234567;
static double rnd(){seed=seed*1664525u+1013904223u;return (seed>>8)/16777216.0;}
static int add(V p,double r,double mass=1){if(bodies.size()>=160)return -1;bodies.push_back({p,{},r,mass>0?1/mass:0});return (int)bodies.size()-1;}
static void joint(int a,int b,double distance){joints.push_back({a,b,distance});}
API void reset(int mode){
 bodies.clear();joints.clear();contacts.clear();scene=std::clamp(mode,0,3);steps=0;peakError=0;selected=-1;seed=1234567;
 if(scene==0){for(int row=0;row<7;row++)for(int col=0;col<7-row;col++)add({6+(col-(6-row)*.5)*.7,7.5-row*.64},.34,1);}
 if(scene==1){int prev=add({2,1},.16,0);for(int i=1;i<15;i++){int b=add({2+i*.55,1.15+i*.04},.16,.7);joint(prev,b,len(bodies[b].p-bodies[prev].p));prev=b;}}
 if(scene==2){for(int i=0;i<90;i++)add({.6+rnd()*10.8,.6+rnd()*5.8},.09+rnd()*.09,.4);}
 if(scene==3){int a=add({6,1},.14,0),b=add({7.3,2.5},.28,1),c=add({8.4,4},.28,1);joint(a,b,len(bodies[b].p-bodies[a].p));joint(b,c,len(bodies[c].p-bodies[b].p));}
}
static void impulse(int a,int b,V n,double lambda){bodies[a].v=bodies[a].v-n*(lambda*bodies[a].inv);if(b>=0)bodies[b].v=bodies[b].v+n*(lambda*bodies[b].inv);}
static void makeContact(int a,int b,V n,double penetration){
 if(penetration < -.015)return;
 V relative=(b<0?V{}:bodies[b].v)-bodies[a].v;
 double closing=dot(relative,n);
 contacts.push_back({a,b,n,-.2*std::max(penetration-.003,0.0)/H,closing < -1 ? restitution*closing:0,0,0});
}
static void tick(){
 for(auto &b:bodies)if(b.inv){b.v.y+=gravity*H;b.v=b.v*.9999;}
 contacts.clear();
 for(int i=0;i<(int)bodies.size();i++){
 auto&a=bodies[i];if(a.inv){makeContact(i,-1,{-1,0},a.radius-a.p.x);makeContact(i,-1,{1,0},a.p.x+a.radius-WIDTH);makeContact(i,-1,{0,-1},a.radius-a.p.y);makeContact(i,-1,{0,1},a.p.y+a.radius-HEIGHT);}
 for(int j=i+1;j<(int)bodies.size();j++){
 auto&b=bodies[j];if(a.inv+b.inv==0)continue;V delta=b.p-a.p;double d=len(delta);if(d<a.radius+b.radius+.015)makeContact(i,j,d>1e-10?delta*(1/d):V{1,0},a.radius+b.radius-d);
 }}
 for(int k=0;k<iterations;k++){
  for(auto &j:joints){auto&a=bodies[j.a];auto&b=bodies[j.b];V d=b.p-a.p;double distance=len(d),eff=a.inv+b.inv;if(distance<1e-10||eff==0)continue;V n=d*(1/distance);double C=distance-j.rest;double lambda=-(dot(b.v-a.v,n)+.15*C/H)/eff;impulse(j.a,j.b,n,lambda);}
  for(auto &c:contacts){auto&a=bodies[c.a];double eff=a.inv+(c.b<0?0:bodies[c.b].inv);if(eff==0)continue;V relative=(c.b<0?V{}:bodies[c.b].v)-a.v;double lambda=-(dot(relative,c.n)+c.bias+c.rebound)/eff;double old=c.normal;c.normal=std::max(0.0,old+lambda);impulse(c.a,c.b,c.n,c.normal-old);
   V t={-c.n.y,c.n.x};relative=(c.b<0?V{}:bodies[c.b].v)-a.v;old=c.tangent;c.tangent=std::clamp(old-dot(relative,t)/eff,-friction*c.normal,friction*c.normal);impulse(c.a,c.b,t,c.tangent-old);
  }
 }
 for(auto &b:bodies)if(b.inv){double speed=len(b.v);if(speed>100)b.v=b.v*(100/speed);b.p=b.p+b.v*H;}
 steps++;peakError=0;for(auto&j:joints)peakError=std::max(peakError,std::abs(len(bodies[j.b].p-bodies[j.a].p)-j.rest));
}
API void advance(int count){for(int i=0;i<std::clamp(count,0,24);i++)tick();}
API int body_count(){return (int)bodies.size();}
API int joint_count(){return (int)joints.size();}
API double body_value(int i,int field){if(i<0||i>=(int)bodies.size())return 0;auto&b=bodies[i];switch(field){case 0:return b.p.x;case 1:return b.p.y;case 2:return b.radius;case 3:return b.inv;case 4:return b.v.x;case 5:return b.v.y;default:return 0;}}
API int joint_value(int i,int field){if(i<0||i>=(int)joints.size())return -1;return field==0?joints[i].a:joints[i].b;}
API double metric(int field){if(field==0)return steps*H;if(field==1)return peakError;if(field==2){double e=0;for(auto&b:bodies)if(b.inv)e+=.5*dot(b.v,b.v)/b.inv;return e;}if(field==3)return contacts.size();return 0;}
API void configure(int it,double g,double e,double mu){if(!std::isfinite(g)||!std::isfinite(e)||!std::isfinite(mu))return;iterations=std::clamp(it,2,40);gravity=std::clamp(g,0.0,20.0);restitution=std::clamp(e,0.0,.9);friction=std::clamp(mu,0.0,1.0);}
API void spawn(double x,double y){if(!std::isfinite(x)||!std::isfinite(y))return;add({std::clamp(x,.4,11.6),std::clamp(y,.4,7.6)},.18+rnd()*.12,1);}
API int grab(double x,double y){if(!std::isfinite(x)||!std::isfinite(y))return -1;selected=-1;double best=.7;for(int i=0;i<(int)bodies.size();i++)if(bodies[i].inv){double d=len(bodies[i].p-V{x,y});if(d<best){selected=i;best=d;}}return selected;}
API void drag(double x,double y){if(selected<0||!std::isfinite(x)||!std::isfinite(y))return;auto&b=bodies[selected];V target{std::clamp(x,b.radius,WIDTH-b.radius),std::clamp(y,b.radius,HEIGHT-b.radius)};b.v=(target-b.p)*12;}
API void release(){selected=-1;}
#ifdef TEST_ENGINE
static bool close(double a,double b,double tol){return std::abs(a-b)<tol;}
int main(){
 // Semi-implicit Euler free fall before any contact.
 bodies.clear();joints.clear();add({6,1},.2);configure(16,9.81,0,.4);advance(12);assert(close(bodies[0].v.y,9.81*H*12,.002));assert(bodies[0].p.y>1);
 // Gravity-free equal masses exchange normal velocity for near-elastic contact.
 bodies.clear();joints.clear();int a=add({5,4},.5),b=add({6,4},.5);bodies[a].v={2,0};configure(20,0,.9,0);advance(1);assert(bodies[b].v.x>1.8);assert(bodies[a].v.x<.11);
 // Distance joint residual converges, static anchor stays fixed.
 reset(3);configure(32,9.81,.1,.4);for(int i=0;i<1200;i++)advance(1);assert(metric(1)<.03);assert(close(bodies[0].p.x,6,1e-10));assert(close(bodies[0].p.y,1,1e-10));
 // Long runs for all scenes: finite state, bounded world, bounded joints.
 for(int sc=0;sc<4;sc++){reset(sc);configure(24,9.81,.25,.4);for(int i=0;i<6000;i++)advance(1);for(auto&x:bodies){assert(std::isfinite(x.p.x)&&std::isfinite(x.p.y));assert(x.p.x>-.1&&x.p.x<12.1&&x.p.y>-.1&&x.p.y<8.1);}if(sc==1||sc==3)assert(metric(1)<.05);}
 // Determinism and fixed step replay.
 reset(2);advance(24);double pos=body_value(7,1);reset(2);for(int i=0;i<24;i++)advance(1);assert(pos==body_value(7,1));
 // Public API rejects non-finite inputs and caps work/storage.
 reset(0);int before=body_count();spawn(NAN,3);assert(body_count()==before);for(int i=0;i<300;i++)spawn(5,3);assert(body_count()==160);assert(body_value(-1,0)==0);configure(500,NAN,0,0);advance(-4);assert(std::isfinite(metric(0)));reset(0);grab(6,7.5);drag(7,7);advance(1);release();
 puts("PASS: freefall, impulse transfer, anchored joints, 4x50s stress, deterministic replay, API bounds and drag");
}
#endif
