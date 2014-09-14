#include <math.h>
#include <iostream>
#include <map>
#include <vector>
#include <list>
#include <string>
#include <assert.h>

// Need this to avoid a glm warning
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "utils.h"

using std::vector;

using glm::vec3;
using glm::vec4;
using glm::ivec3;
using glm::ivec4;
using glm::dot;
using glm::cross;
using glm::length;

// Global variables

// Drawing options
int drawtype = 0;
bool snubify = false;

int style = 0;
bool facecolor = true;
bool doinvert = true;

bool dolog = false;

int nstyles = 2;
static const double pi = 3.14159265358979324;
static const double twopi = 2*pi;
double eps = 1e-4;
double meps = 1e4;

// Trilinears
vec3 afact; // Use to translate between trilinear and barycentric coords
vec3 snuba; // Tri coords of the centre of a snub face

// Sphere radius to perform inversion with
float midsphere2;
float snubsphere2;

// Set this when we need to recalculate
bool needparams = true;

// Geometry

int compound = 0;

// Constant data per Schwarz triangle
std::vector<vec3> points;       // Vertex positions for regions
std::vector<ivec4> regions;     // Fundamental regions in Wythoff construction
std::map<int,ivec3> adjacent;   // Adjacency relation for regions

// Construct the Wythoff faces here. Have 3 maps
// for the 3 different kinds of faces. Each face is
// represented by a list of fundamental regions.
std::map<int,std::vector<int> > facemap[3];

// Data varying by point in Schwarz triangle
std::vector<vec3> regionpoints; // Region point in Wythoff construction
// Distance of face centre from origin.
vec3 facecentres;
std::vector<vec3> snubcentres;

// Color and lighting
vec4 *facecolors[] = { &red, &green, &blue, &yellow, &white, &magenta, &cyan, &orange, &purple, &grey };
int nfacecolors = 10;

// Whether to display faces
bool omit[4] = { false, false, false, false };
float explode = 0;

bool approxeq(const vec3 &p, const vec3 &q)
{
  return fabs(p.x - q.x) < eps &&
    fabs(p.y - q.y) < eps &&
    fabs(p.z - q.z) < eps;
}

// A little dubious (no w-scaling), but we only use it to check for duplicates
// in point lists.
bool approxeq(const vec4 &p, const vec4 &q)
{
  return fabs(p.x - q.x) < eps && fabs(p.y - q.y) < eps &&
    fabs(p.z - q.z) < eps && fabs(p.w - q.w) < eps;
}

vec3 normalize(vec3 p)
{
  assert(length(p) > eps);   // Never normalize a short vector!
  return glm::normalize(p);
}

// Geometry stuff
// TBD: use getbary for this
// TBD: use vec4 for points
vec3 tri2bary(const vec3 &p, const vec3 &q, const vec3 &r,
              const vec3 &tri)
{
   vec3 np = normalize(cross(q,r));
   vec3 nq = normalize(cross(r,p));
   vec3 nr = normalize(cross(p,q));
   vec3 t;
   if (snubify) {
     // TBD: This can assert...
     // TBD: Why the negation?
     t = -normalize(p*tri[0]*snuba[0] + q*tri[1]*snuba[1] + r*tri[2]*snuba[2]);
   } else {
     // nq x nr = norm(r x p) x norm(p x q) = kp
     t = -normalize(tri[0]*cross(nq,nr) + tri[1]*cross(nr,np) + tri[2]*cross(np,nq));
   }
   float P = dot(t,np)/dot(p,np);
   float Q = dot(t,nq)/dot(q,nq);
   float R = dot(t,nr)/dot(r,nr);
   return vec3(P,Q,R);
}

void cosine (double A, double B, double C, 
             double &a, double &b, double &c)
{
  a = acos((cos(A) + cos(B)*cos(C))/(sin(B)*sin(C)));
  b = acos((cos(B) + cos(C)*cos(A))/(sin(C)*sin(A)));
  c = acos((cos(C) + cos(A)*cos(B))/(sin(A)*sin(B)));
  assert(!isnan(a));
  assert(!isnan(b));
  assert(!isnan(c));
}

// Construct a basic Schwarz triangle with given base angles
void triangle1(double a, double b, double c,
               vec3 &p, vec3 &q, vec3 &r)
{
  p = vec3(0,0,1);
  q = vec3(0,sin(c), cos(c));
  double rz = cos(b);
  double ry = (cos(a)-q.z*rz)/q.y;
  double t = 1-rz*rz-ry*ry;
  r = vec3(sqrt(std::max(t,0.0)),ry,rz);
}

vec3 reflect(vec3 p, vec3 q, vec3 r)
{
  // Reflect r in the plane of p and q (and the origin)
  vec3 n = normalize(cross(p,q));
  return r - 2*dot(n,r)*n;
}

// Drawing functions

// Interface to OpenGL
void setColor(const vec4 &color);
void setNormal(const vec3 &norm);
void drawvertex(const vec4 &p, const vec3 &offset);

// Initialization functions

// Point and vertex store
int getnode(const vec3 &p)
{
  // Linear scan, but the most points we have is 120
  int nnodes = points.size();
  for (int i = 0; i < nnodes; i++) {
    if (approxeq(p,points[i])) return i;
  }
  points.push_back(p);
  return nnodes;
}

bool addtriangle(const ivec4 &t, int &index)
{
  int ntriangles = regions.size();
  for (int i = 0; i < ntriangles; i++) {
    if (regions[i] == t) {
      index = i;
      return false;
    }
  }
  regions.push_back(t);
  index = ntriangles;
  return true;
}

// Solving snub equation.
// Reflect s in 3 sides to make triangle, find departure from
// equilateral as difference between longest and shortest side.
float eval(const vec3 &p, const vec3 &q, const vec3 &r, const vec3 &s)
{
  vec3 s0 = reflect(p,q,s);
  vec3 s1 = reflect(q,r,s);
  vec3 s2 = reflect(r,p,s);
  float l1 = length(s0-s1);
  float l2 = length(s1-s2);
  float l3 = length(s2-s0);
  float max = std::max(l1,std::max(l2,l3));
  float min = std::min(l1,std::min(l2,l3));
  float res = (max-min)/(max+min);
  return res;
}

// Compounds

vec3 ip,iq,ir;
int zrotarg = 0;
double zrottot = 0;

void makeico(ivec3 iargs)
{
  double A,B,C,a,b,c;
  A = pi/iargs[0];
  B = pi/iargs[1];
  C = pi/iargs[2];
  cosine(A,B,C,a,b,c);
  triangle1(a,b,c,ip,iq,ir);
}

vec3 zrot(vec3 p, float theta)
{
  return vec3(cos(theta)*p.x + sin(theta)*p.y,
              -sin(theta)*p.x + cos(theta)*p.y,
              p.z);
}

void makeschwarz(int args[6])
{
  double A,B,C,a,b,c;
  vec3 p0,q0,r0;
  A = pi*args[1]/args[0];
  B = pi*args[3]/args[2];
  C = pi*args[5]/args[4];
  cosine(A,B,C,a,b,c);
  triangle1(a,b,c,p0,q0,r0);
  if (zrotarg != 0) {
    float theta = pi/zrotarg;
    p0 = zrot(p0,theta);
    q0 = zrot(q0,theta);
    r0 = zrot(r0,theta);
    zrottot -= theta;
  }
  snuba = vec3(1,1,1);
  // Compute snub coordinates
  find(snuba[0],snuba[1],snuba[2],0.1, 60,
       [&](float a0, float b0, float c0) {
         vec3 s = normalize(a0*p0 + b0*q0 + c0*r0);
         float t = eval(p0,q0,r0,s);
         return t;
       });
  vec3 s = normalize(snuba[0]*p0 + snuba[1]*q0 + snuba[2]*r0);
  float res = eval(p0,q0,r0,s);
  if (res > eps || snuba[0] < eps || snuba[1] < eps || snuba[2] < eps) {
    fprintf(stderr,"Snubification failure: %g %g %g %g\n", snuba[0], snuba[1], snuba[2], res);
  }

  struct Entry {
    Entry(ivec4 &t, int i) : triangle(t), index(i) {}
    ivec4 triangle;
    int index;
  };
  std::list<Entry> qq;
  getnode(p0); getnode(q0); getnode(r0);
  ivec4 t0(0,1,2,1);
  int index;
  addtriangle(t0,index);
  assert(index == 0);
  qq.push_back(Entry(t0,index));
  for ( ; !qq.empty(); qq.pop_front()) {
    Entry &entry = qq.front();
    ivec4 &t = entry.triangle;
    vec3 p = points[t[0]];
    vec3 q = points[t[1]];
    vec3 r = points[t[2]];
    // We could compute p1, q1 and r1 just with barycentric coordinates.
    int p1 = getnode(reflect(q,r,p));
    int q1 = getnode(reflect(r,p,q));
    int r1 = getnode(reflect(p,q,r));
    ivec4 t1(p1,t[1],t[2],-t[3]);
    ivec4 t2(t[0],q1,t[2],-t[3]);
    ivec4 t3(t[0],t[1],r1,-t[3]);
    int a1,a2,a3;
    if (addtriangle(t1,a1)) qq.push_back(Entry(t1,a1));
    if (addtriangle(t2,a2)) qq.push_back(Entry(t2,a2));
    if (addtriangle(t3,a3)) qq.push_back(Entry(t3,a3));
    adjacent[entry.index] = ivec3(a1,a2,a3);
  }
}

// Sort out the faces for the given Schwarz triangles
// v is the point index
void makefaces(int v, int p, int q, int r) 
{
  std::vector<int> facelist;
  // Make list of all faces that contain this point as p'th element.
  int n = regions.size();
  for (int i = 0; i < n; i++) {
    if (regions[i][p] == v) {
      facelist.push_back(i);
    }
  }
  if (!facelist.empty()) {
    int fsize = facelist.size();
    for (int i = 0; i < fsize; i++) {
      if (regions[facelist[i]][3] > 0) {
        std::swap(facelist[0],facelist[i]);
        break;
      }
    }
    assert(regions[facelist[0]][3] == 1);
    // Match up successive regions.
    for (int i = 0; i < fsize; i++) {
      for (int j = i+1; j < fsize; j++) {
        if (((i%2 == 0) && regions[facelist[i]][r] == regions[facelist[j]][r]) ||
            ((i%2 == 1) && regions[facelist[i]][q] == regions[facelist[j]][q])) {
          std::swap(facelist[j],facelist[i+1]);
          break;
        }
      }
    }
    // Insert into the p'th facemap.
    facemap[p][v].swap(facelist);
  }
}

// For snubification, we want, for every -ve triangle, to draw a triangle between
// the Wythoff points of the surrounding +ve trianges.

// Compute the distance of each face type from origin
void setpointface(int k, int v, std::vector<int> &facelist)
{
  assert(!facelist.empty());
  vec3 centre;
  int npoints = 0;
  for (int i = 0; i < (int)facelist.size(); i++) {
    int n = facelist[i];
    vec3 &p = regionpoints[n];
    centre += p;
    npoints++;
  }
  centre /= npoints;
  float len = dot(centre,points[v]);
  if (fabs(len) < eps) len = 0;
  facecentres[k] = len;
}

vec3 getnormal0(const vec4 &p0, const vec4 &p1, const vec4 &p2)
{
  // Avoid multiplying by infinity (it's the same as dividing by zero).
  vec3 d1;
  if (p1[3] > meps) d1 = vec3(p0/p0[3]);
  else d1 = vec3(p1) - vec3(p0*p1[3]/p0[3]);
  vec3 d2;
  if (p2[3] > meps) d2 = vec3(p0/p0[3]);
  else d2 = vec3(p2) - vec3(p0*p2[3]/p0[3]);
  return cross(d1,d2);
}

vec3 getnormal(const vec4 &p0, const vec4 &p1, const vec4 &p2)
{
#if 1
  if (fabs(p0[3] > eps)) return getnormal0(p0,p1,p2);
  else if (fabs(p1[3] > eps)) return getnormal0(p1,p2,p0);
  else if (fabs(p2[3] > eps)) return getnormal0(p2,p0,p1);
  else assert(0);
#else
  vec3 q0 = vec3(p0)/p0[3];
  vec3 q1 = vec3(p1)/p1[3];
  vec3 q2 = vec3(p2)/p2[3];
  return cross(q1-q0,q2-q0);
#endif
}

void drawtriangle0(const vec4 &p0, const vec4 &p1, const vec4 &p2,
                   const vec4 &c0, const vec4 &c1, const vec4 &c2,
                   const vec3 &offset)
{
  vec3 norm = getnormal(p0,p1,p2);
  if (dolog) {
    std::cout << "Normal:" << norm << p0 << p1 << p2 << "\n";
  }
  if (length(norm) > eps) {
    setNormal(normalize(norm));
    setColor(c0);
    drawvertex(p0,offset);
    setColor(c1);
    drawvertex(p1,offset);
    setColor(c2);
    drawvertex(p2,offset);
  }
}

void drawtriangle0(const vec4 &p0, const vec4 &p1, const vec4 &p2,
                   const vec4 &c,  const vec3 &offset)
{
  drawtriangle0(p0,p1,p2,c,c,c,offset);
}

void drawpointface(int k, int v, std::vector<int> &facelist, int snubface)
{
  assert(!facelist.empty());
  std::vector<vec4> plist;
  for (int i = 0; i < (int)facelist.size(); i++) {
    int n = facelist[i];
    ivec4 &f = regions[n];
    if (!snubify || f[3] != snubface) {
      vec4 p = vec4(regionpoints[n],1);
      // Suppress equal points so triangle code gets a chance to run
      if (plist.empty() || !approxeq(p,plist.back())) {
        plist.push_back(p);
      }
    }
  }
  if (approxeq(plist[0],plist.back())) plist.pop_back();
  int npoints = plist.size();
  if (npoints <= 2) return;
  vec3 centre = facecentres[k]*points[v];
  vec3 offset = centre;
  int color = facecolor ? k : compound;
  if (npoints == 3) {
    drawtriangle0(plist[0],plist[1],plist[2],
                  *facecolors[color%nfacecolors],
                  offset);
  } else {
    for (int i = 0; i < npoints; i++) {
      int p1 = i;
      int p2 = (i+1)%npoints;
      drawtriangle0(vec4(centre,1),plist[p1],plist[p2],
                    *facecolors[color%nfacecolors],
                    offset);
    }
  }
}

int imax = 300;
int ifact = 0;
int iinc = -1;

vec4 invert(float r2, const vec4 &p4)
{
  vec3 p = vec3(p4);
  float w = p4[3];
  assert(w >= 0);
  assert(r2 >= 0);
#if 0
  // Cute animation. Needs fixing for points at infinity.
  float k = dot(p,p)/(r2*w);
  float t = (w+(k-w)*float(ifact)/imax);
  return vec4(p,t); 
#else
  float k = dot(p,p)/(r2*w);
  return vec4(p,k);
#endif
}

vec4 makevec4(const vec3 &p, float k)
{
  if (k < 0) return vec4(-p,-k);
  else return vec4(p,k);
}

void drawregion(int r) {
  std::vector <vec4> plist;
  for (int i = 0; i < 3; i++) {
    int v = regions[r][i];
    float len = facecentres[i];
    vec4 p = makevec4(points[v],1/len);
    if (dolog) std::cout << "# A: " << p << "\n";
    if (doinvert) p = invert(midsphere2,p);
    if (dolog) std::cout << "# B: " << p << "\n";
    plist.push_back(p);
  }
  if (regions[r][3] < 0) std::swap(plist[1],plist[2]);
  assert(length(regionpoints[r]) > eps);
  vec4 centre = makevec4(regionpoints[r],1);
  if (doinvert) centre = invert(midsphere2,centre);
  if (dolog) std::cout << "Centre: " << centre << "\n";
  vec3 offset = vec3(centre);
  if (style == 0) {
    // TBD: use compound color here too.
    int color = (regions[r][3] < 0) ? 3 : 4;
    drawtriangle0(plist[0],plist[1],plist[2],*facecolors[color],offset);
  } else {
    for (int i = 0; i < 3; i++) {
      int p1 = i, p2 = (i+1)%3;
      drawtriangle0(centre,plist[p1],plist[p2],
                    white, *facecolors[p1], *facecolors[p2],
                    offset);
    }
  }
}

void setsnubface(int i)
{
  // TBD: This should compute barycentric coords for the centre
  // of the snub triangle.
  vec3 triangle[3];
  for (int j = 0; j < 3; j++) {
    triangle[j] = regionpoints[adjacent[i][2-j]];
  }
  if (regions[i][3] > 0) std::swap(triangle[1],triangle[2]);
  snubsphere2 = mid3(square(length(triangle[0]+triangle[1])/2),
                     square(length(triangle[1]+triangle[2])/2),
                     square(length(triangle[2]+triangle[0])/2));
  vec3 norm = cross(triangle[0]-triangle[1],triangle[1]-triangle[2]);
  if (length(norm) <= eps) {
    snubcentres[i] = vec3(0,0,0);
  } else {
    norm = normalize(norm);
    float a = dot(triangle[0],norm); // Can probably simplify this
    snubcentres[i] = a*norm;
  }
}

void drawsnubface(int i)
{
  vec4 triangle[3];
  vec4 centre;
  for (int j = 0; j < 3; j++) {
    triangle[j] = vec4(regionpoints[adjacent[i][2-j]],1);
    centre += triangle[j];
  }
  centre /= 3.0f;
  if (regions[i][3] > 0) std::swap(triangle[1],triangle[2]);
  vec3 offset = vec3(centre);
  int color;
  if (facecolor) {
    color = 3;
  } else {
    color = compound;
  }
  for (int j = 0; j < 3; j++) {
    int p1 = j;
    int p2 = (j+1)%3;
    drawtriangle0(centre,triangle[p1],triangle[p2],
                  *facecolors[color%nfacecolors],
                  offset);
  }
}

// Find A,B,C such that Ap + Bq + Cr = t
// t.(q x r) = Ap.(q x r) ie. A = t.(q x r)/k (k as above)
vec4 getbary(const vec3 &p, const vec3 &q, const vec3 &r, const vec4 &t)
{
  vec3 P = cross(q,r);
  vec3 Q = cross(r,p);
  vec3 R = cross(p,q);
  float k = dot(p,P);
  vec3 t3 = vec3(t);
  float A = dot(t3,P);
  float B = dot(t3,Q);
  float C = dot(t3,R);
  vec4 s = vec4(vec3(A,B,C),k*t[3]);
  return s;
}

void drawinit()
{
  regionpoints.resize(0);
  for (ivec4 &f : regions) {
    vec3 &p0 = points[f[0]];
    vec3 &p1 = points[f[1]];
    vec3 &p2 = points[f[2]];
    vec3 p = afact[0]*p0 + afact[1]*p1 + afact[2]*p2;
    regionpoints.push_back(p);
  }
  snubcentres.resize(regions.size());
  for (int i = 0; i < 3; i++) {
    for (auto &j : facemap[i]) {
      setpointface(i,j.first,j.second);
      break;
    }
  }
  for (int i = 0; i < (int)regions.size(); i++) {
    setsnubface(i);
  }
}

void doparams(const vec3 &abc)
{
  // Compute trilinear->barycentric factors if needed
  float A = abc[0];
  float B = abc[1];
  float C = abc[2];

  vec3 &p = points[regions[0][0]];
  vec3 &q = points[regions[0][1]];
  vec3 &r = points[regions[0][2]];
  vec3 tri(A,B,C);
  afact = tri2bary(p,q,r,tri);
  vec3 s = afact[0]*p + afact[1]*q + afact[2]*r;
  // Centres of edges. Take min to avoid zero sized edges.
  midsphere2 = square(min3(length(s+reflect(p,q,s)),
                           length(s+reflect(q,r,s)),
                           length(s+reflect(r,p,s)))/2);
  drawinit();
}

void drawDual()
{
  if (!snubify) {
    for (int i = 0; i < (int)regions.size(); i++) {
      drawregion(i);
    }
  } else {
    int snubface = 1;
    for (int i = 0; i < (int)regions.size(); i++) {
      if (regions[i][3] == snubface) {
        std::vector<vec4> plist;
        for (int j = 0; j < 3; j++) {
          int v = regions[i][j];
          vec4 p = vec4(points[v]*facecentres[j],1);
          if (doinvert) p = invert(snubsphere2,p);
          plist.push_back(p);
          vec4 tmp = vec4(snubcentres[adjacent[i][(j+2)%3]],1);
          if (length(vec3(tmp)) > eps) {
            if (doinvert) tmp = invert(snubsphere2,tmp);
            plist.push_back(tmp);
          }
        }
        int N = plist.size();
        assert(length(regionpoints[i]) > eps);
        vec4 centre = vec4(regionpoints[i],1);
        if (doinvert) centre = invert(snubsphere2,centre);
        vec3 offset = vec3(centre);
        for (int k = 0; k < N; k++) {
          int p1 = k;
          int p2 = (k+1)%N;
          drawtriangle0(centre,plist[p1],plist[p2],
                        white,
                        *facecolors[p1%nfacecolors],
                        *facecolors[p2%nfacecolors],
                        offset);
        }
      }
    }
  }
}

void drawMain()
{
  int snubface = -1;
  for (int i = 0; i < 3; i++) {
    if (!omit[i]) {
      for (auto &j : facemap[i]) {
        if (snubify) drawpointface(i,j.first,j.second,snubface);
        else drawpointface(i,j.first,j.second,0);
      }
    }
  }
  if (snubify && !omit[3]) {
    for (int i = 0; i < (int)regions.size(); i++) {
      if (regions[i][3] == snubface) drawsnubface(i);
    }
  }
}


#include <GL/glut.h>
#include <GL/gl.h>

// Drawing stuff

int args[6];
ivec3 iargs;
bool docompound = false;
std::string trans;
int dozrot = 0;
bool rotonly = false;
bool docull = false;
bool doboth = true;

Interpolator<vec3> interpolator;
bool pause = true;
int screen_width=800;
int screen_height=800;

bool filling = true; // When false, draw a wire frame of each face
// Animation
double rotx = 0;
double roty = 0;
double rotz = 0;
double rotxinc = 0;
double rotyinc = 0;
double rotzinc = 0;
double zpos = -30;
double zinc = 0.0;
float size = 8;

double tt = 0;
float ttinc = 0.0015;
vec3 abc;

// Lights settings
float amb = 0.2;
float diff = 0.2;
float spec = 0.3;
GLfloat light_ambient[]= { amb, amb, amb, 0.0f };
GLfloat light_diffuse[]= { diff, diff, diff, 0.0f };
GLfloat light_specular[]= { spec, spec, spec, 0.0f };

// Lights are at infinity
GLfloat light_position1[]= { 100.0f, 100.0f, 100.0f, 0.0f };
GLfloat light_position2[]= { -100.0f, 100.0f, 100.0f, 0.0f };

// Materials settings
GLfloat mat_specular[]= { 0, 0, 0, 0 };
GLfloat mat_shininess = 100;
GLfloat specular = 1;

void setColor(const vec4 &color)
{
  if (dolog) {
    printf("color %g %g %g %g\n", color[0], color[1], color[2], color[3]);
  }
  const float *p = glm::value_ptr(color);
  glMaterialfv (GL_FRONT_AND_BACK, GL_DIFFUSE, p);
  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT, p);
}

vec3 P(vec3 p) { return reflect(iq,ir,p); }
vec3 Q(vec3 p) { return reflect(ir,ip,p); }
vec3 R(vec3 p) { return reflect(ip,iq,p); }

vec3 apply(const std::string &s, vec3 p)
{
  for (int i = 0; i < (int)s.length(); i++) {
    if (s[i] == 'P') p = P(p);
    else if (s[i] == 'Q') p = Q(p);
    else if (s[i] == 'R') p = R(p);
  }
  return p;
}

void setNormal(const vec3 &norm)
{
  vec3 n = apply(trans,norm);
  if (dolog) printf("normal %g %g %g\n", norm[0], norm[1], norm[2]);
  glNormal3fv(glm::value_ptr(n));
}

void drawvertex(const vec4 &p, const vec3 &offset)
{
  assert(p[3] >= 0);
  vec4 p0 = vec4(apply(trans,vec3(p)),p[3]);
  if (p[3] > meps) p0 = vec4(0,0,0,1);
  if (explode != 0) {
    p0 += vec4(explode*offset*p0[3],0);
  }
  if (dolog) printf("vertex %g %g %g %g\n", p0[0], p0[1], p0[2], p0[3]);
  glVertex4f(size*p0[0],size*p0[1],size*p0[2],p0[3]);
}

void drawone()
{
  glBegin(GL_TRIANGLES);
  if (drawtype <= 1) drawMain();
  if (drawtype >= 1) drawDual();
  glEnd();
}

void init(void)
{
  glLightfv (GL_LIGHT1, GL_AMBIENT, light_ambient);
  glLightfv (GL_LIGHT1, GL_DIFFUSE, light_diffuse);
  glLightfv (GL_LIGHT1, GL_SPECULAR, light_specular);
  glLightfv (GL_LIGHT1, GL_POSITION, light_position1);
  glEnable (GL_LIGHT1);

  glLightfv (GL_LIGHT2, GL_AMBIENT, light_ambient);
  glLightfv (GL_LIGHT2, GL_DIFFUSE, light_diffuse);
  glLightfv (GL_LIGHT2, GL_SPECULAR, light_specular);
  glLightfv (GL_LIGHT2, GL_POSITION, light_position2);
  glEnable (GL_LIGHT2);

  glEnable (GL_LIGHTING);

  glShadeModel(GL_SMOOTH);
  glHint (GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
  glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
  glEnable(GL_DEPTH_TEST);
}

void resize (int p_width, int p_height)
{
  assert(screen_width != 0 && screen_height != 0);
  screen_width = p_width;
  screen_height = p_height;
  glutPostRedisplay ();
}

bool contains(const std::vector<vec3> &s, const vec3 &p)
{
  for (auto q: s) {
    if (approxeq(p,q)) return true;
  }
  return false;
}

void addpoint(std::vector<vec3> &s, const vec3 &p)
{
  if (!contains(s,p)) s.push_back(p);
}

// We could do this more efficiently
bool eqset(const std::vector<vec3> &s1, const std::vector<vec3> &s2)
{
  if (s1.size() != s2.size()) return false;
  for (auto p: s1) {
    if (!contains(s2,p)) return false;
  }
  return true;
}

void display(void)
{
  if (dolog) {
    printf ("start\n");
    printf (";; %d/%d %d/%d %d/%d\n", 
            args[0],args[1],args[2],args[3],args[4],args[5]);
    printf (";; %g %g %g%s: %d\n",
            abc[0],abc[1],abc[2],
            snubify?" snub":"",
            drawtype);
  }

  ifact += iinc;
  if (ifact < 0) ifact = 0;
  if (ifact > imax) ifact = imax;

  rotx += rotxinc;
  roty += rotyinc;
  rotz += rotzinc;

  if (rotx >= 360) rotx -= 360;
  if (roty >= 360) roty -= 360;
  if (rotz >= 360) rotz -= 360;

  if (dozrot) {
    zrottot += 0.001;
    ip = zrot(ip,0.001);
    iq = zrot(iq,0.001);
    ir = zrot(ir,0.001);
  }

  if (docull) glEnable(GL_CULL_FACE); // Enable the back face culling
  else glDisable(GL_CULL_FACE); 
  if (doboth) glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  else glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

  glClearColor(0, 0, 0.2, 0.0);
  glViewport(0,0,screen_width,screen_height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0f,(GLfloat)screen_width/(GLfloat)screen_height,1.0f,1000.0f);

  glTranslatef(0.0,0.0,zpos);
  glMatrixMode(GL_MODELVIEW);
  glRotatef(rotxinc,1.0,0.0,0.0);
  glRotatef(rotyinc,0.0,1.0,0.0);
  glRotatef(rotzinc,0.0,0.0,1.0);
  zpos += zinc;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // This clear the background color to dark blue

  mat_specular[0] = mat_specular[1] = mat_specular[2] = specular;
  glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
  glMaterialfv (GL_FRONT_AND_BACK, GL_SHININESS, &mat_shininess);    
  if (needparams) {
    abc = interpolator(tt);
    doparams(abc);
    needparams = false;
  }
  if (!pause) {
    needparams = true;
    tt += ttinc;
  }
  std::list<std::string> qq;
  typedef std::vector<vec3> pointset;
  std::vector<pointset> seen;
  int rsize = regions.size();
  qq.push_back("");
  compound = 0;
  while (!qq.empty()) {
    trans = qq.front();
    qq.pop_front();
    pointset newpoints;
    for (int i = 0; i < rsize; i++) {
      if (!snubify || regions[i][3] == 1) {
        vec3 pp = apply(trans,regionpoints[i]);
        addpoint(newpoints,pp);
      }
    }
    bool done = false;
    for (auto s: seen) {
      if (eqset(newpoints,s)) {
        done = true;
        break;
      }
    }
    if (!done) {
      seen.push_back(pointset());
      newpoints.swap(seen.back()); // swap in new data
      if (rotonly) {
        qq.push_back(trans+"PQ");
        qq.push_back(trans+"QR");
        qq.push_back(trans+"RP");
      } else {
        qq.push_back(trans+"P");
        qq.push_back(trans+"Q");
        qq.push_back(trans+"R");
      }
      glFrontFace(trans.length()%2 == 0 ? GL_CCW : GL_CW);
      drawone();
      compound++;
      if (!docompound) break;
    }
  }
  static int lastcompound = -1;
  if (compound != lastcompound){
    lastcompound = compound;
    printf("Compound: %d (%d points)\n",
           compound, (int)seen[0].size());
  }
  if (dolog) {
    printf ("end\n");
    dolog = false;
  }
  glFlush();
  glutSwapBuffers();
}


void keyboard(unsigned char p_key, int p_x, int p_y)
{  
  switch (p_key) {
  case 'i':
    doinvert = !doinvert;
    iinc = -iinc;
    needparams = true;
    break;
  case 'b':
    doboth = !doboth;
    printf("doboth = %d\n", doboth);
    break;
  case 'f':
    facecolor = !facecolor;
    break;
  case 'z':
    dozrot = !dozrot;
    if (!dozrot) printf("zrot = %f\n", pi/zrottot);
    break;
  case 'x':
    rotonly = !rotonly;
    break;
  case 'c':
    if (length(ip) > eps) docompound = !docompound;
    break;
  case 'l':
    dolog = true;
    break;
  case 'k':
    style = (style+1)%nstyles;
    break;
  case ']':
    tt = floor(tt+1);
    needparams = true;
    break;
  case '[':
    tt = ceil(tt-1.1);
    needparams = true;
    break;
  case ',':
    zpos += 0.25;
    break;
  case '.':
    zpos -= 0.25;
    break;
  case '-':
    ttinc -= 0.0001;
    break;
  case '=':
    ttinc += 0.0001;
    break;
  case 'p':
    rotxinc = rotyinc = rotzinc = 0;
    break;
  case ' ':
    pause = !pause;
    break;
  case 'v':
    for (int i = 1; i < nfacecolors; i++) {
      std::swap(facecolors[i-1],facecolors[i]);
    }
    break;
  case 's':
    snubify = !snubify;
    needparams = true; // We use different coord transform for snub figures
    break;
  case 'r': case 'R':
    filling = !filling;
    glPolygonMode (GL_FRONT_AND_BACK, filling?GL_FILL:GL_LINE);
    break;
  case 'd':
    drawtype = (drawtype+1)%3;
    break;
  case 'e':
    explode += 0.1;
    break;
  case 'w':
    explode -= 0.1;
    break;
  case '1':
    omit[0] = !omit[0];
    break;
  case '2':
    omit[1] = !omit[1];
    break;
  case '3':
    omit[2] = !omit[2];
    break;
  case '4':
    omit[3] = !omit[3];
    break;
  case 27:
    exit(0);
    break;
  }
}

void keyboard_s (int p_key, int p_x, int py)
{
  switch (p_key) {
  case GLUT_KEY_UP:
    rotxinc += 0.005;
    break;
  case GLUT_KEY_DOWN:
    rotxinc -= 0.005;
    break;
  case GLUT_KEY_LEFT:
    rotyinc += 0.005;
    break;
  case GLUT_KEY_RIGHT:
    rotyinc -= 0.005;
    break;
  }
}

int main(int argc, char **argv)
{
  //char *progname = argv[0];
  argc--; argv++;
  while (true) {
    if (strcmp(argv[0],"--zrot") == 0) {
      argc--; argv++;
      zrotarg = atoi(argv[0]);
      argc--; argv++;
    } else {
      break;
    }
  }
  for (int i = 0; i<3; i++) {
    args[2*i] = 1;
    args[2*i+1] = 1;
    int d = sscanf(argv[i],"%d/%d",&args[2*i],&args[2*i+1]);
    assert(d == 1 || d == 2);
  }
  makeschwarz(args);
  if (argc > 3) {
    iargs[0] = atoi(argv[3]);
    iargs[1] = atoi(argv[4]);
    iargs[2] = atoi(argv[5]);
    makeico(iargs);
  }
  {
    int npoints = points.size();
    for (int i = 0; i < npoints; i++) {
      makefaces(i,0,1,2);
      makefaces(i,1,2,0);
      makefaces(i,2,0,1);
    }
  }
  interpolator.add(vec3(0,1,0));
  interpolator.add(vec3(0,1,1));
  interpolator.add(vec3(0,0,1));
  interpolator.add(vec3(1,1,1));
  interpolator.add(vec3(1,0,0));
  interpolator.add(vec3(1,1,0));
  interpolator.add(vec3(0,1,0));
  interpolator.add(vec3(0,1,1));
  interpolator.add(vec3(0,0,1));
  interpolator.add(vec3(1,0,1));
  interpolator.add(vec3(1,0,0));
  interpolator.add(vec3(1,1,1));
  interpolator.add(vec3(1,1,0));
  interpolator.add(vec3(0,1,0));
  interpolator.add(vec3(1,1,1));
  interpolator.add(vec3(0,1,1));
  interpolator.add(vec3(0,0,1));

  glutInit(&argc, argv);   
 
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  glutInitWindowSize(screen_width,screen_height);
  glutInitWindowPosition(0,0);
  glutCreateWindow("OpenGL demo - To exit press ESC");    
  glutDisplayFunc(display);
  glutIdleFunc(display);
  glutReshapeFunc (resize);
  glutKeyboardFunc (keyboard);
  glutSpecialFunc (keyboard_s);
  init();
  glutMainLoop();
  return(0);    
}
