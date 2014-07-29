#include <math.h>
#include <iostream>
#include <map>
#include <vector>
#include <list>
#include <assert.h>

// Need this to avoid a glm warning
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "utils.h"

// The classics
// ./poly 2 3 3
// ./poly 2 3 4
// ./poly 2 3 5
// Nice snub forms
// ./poly 3 5/2 5/3  // All 4 colours in snub
// ./poly 5/2 3/2 3/2 // Excellent dual animation
// ./poly 2 3/2 5/3 // Strange incomplete dual momentary form.
// ./poly 3 2 5/3    // Nice duals
// ./poly 2 2 37/19  // Nice in wireframe
// ./poly 3/2 3/2 3/2 // "partially degenerate" dual - concentric faces in main

// 3 5/3 (3/2 5/2) | is interesting intersection polyhedron

// Check
//./poly 5 2 3/2
//Snubification failure: 1.05876 1.80118 8.79375e-05 1.1451e-05

using std::max;
using std::vector;
using std::list;
using std::ostream;
using std::cerr;
using std::cout;

using glm::vec3;
using glm::vec4;
using glm::ivec3;
using glm::ivec4;
using glm::dot;
using glm::cross;
using glm::length;

static const double pi = 3.14159265358979324;
static const double twopi = 2*pi;
double eps = 1e-4;

typedef vec3 Point;
typedef vec4 Point4;

////////////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////////////

// Drawing options
bool dual = false;
bool snubify = false;

int style = 0;
int nstyles = 2;
bool doclamp = false;

// Trilinears
vec3 afact; // Use to translate between trilinear and barycentric coords
vec3 snuba; // Ditto in snub mode - produces uniform snub at (1,1,1)
vec3 dfact; // Perform dual transform on Schwarz triangle
vec3 efact; // Perform dual transform on Schwarz triangle

float midsphere2;
float snubsphere2;

bool needparams = true;

// Geometry

// Constant data per Schwarz triangle
std::vector<ivec4> regions;     // Fundamental regions in Wythoff construction
std::vector<ivec3> adjacent;    // Adjacency relation for regions
std::vector<vec3> points;       // Vertex positions for regions

// Construct the Wythoff faces here. Have 3 maps
// for the 3 different kinds of faces. Each face is
// represented by a list of fundamental regions.
std::map<int,std::vector<int> > facemap[3];

// Data varying by point in Schwarz triangle
std::vector<vec3> regionpoints; // Region point in Wythoff construction
//std::map<int,std::vector<vec3>> facecentres;
std::map<int,std::vector<float>> facecentres;
std::vector<vec3> snubcentres;

// Colour and lighting
vec4 *facecolours[] = { &red, &green, &blue, &yellow, &white, &white, &white, &white };
int nfacecolours = 8;

// Whether to display faces
bool omit[4] = { false, false, false, false };
float explode = 0;

// Never normalize a short vector!
vec3 normalize(vec3 p)
{
  assert(length(p) > eps);
  return glm::normalize(p);
}

// Geometry stuff
vec3 tri2bary(const vec3 &p, const vec3 &q, const vec3 &r,
              const vec3 &tri)
{
   vec3 np = normalize(cross(q,r));
   vec3 nq = normalize(cross(r,p));
   vec3 nr = normalize(cross(p,q));
   vec3 t;
   if (snubify) {
     t = normalize(p*tri[0]*snuba[0] + q*tri[1]*snuba[1] + r*tri[2]*snuba[2]);
   } else {
     t = normalize(tri[0]*cross(nq,nr) + tri[1]*cross(nr,np) + tri[2]*cross(np,nq));
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
  q = vec3(sin(c), 0, cos(c));
  double rz = cos(b);
  double rx = (cos(a)-q.z*rz)/q.x;
  double t = 1-rz*rz-rx*rx;
  r = vec3(rx,sqrt(max(t,0.0)),rz);
}

vec3 reflect(vec3 p, vec3 q, vec3 r)
{
  // Reflect r in the plane of p and q (and the origin)
  vec3 n = normalize(cross(p,q));
  return r - 2*dot(n,r)*n;
}

// Drawing functions

// Interface to OpenGL
void setColour(const vec4 &colour);
void setColour(const vec4 *colour);
void setNormal(const vec3 &norm);
void drawvertex(const vec3 &p);
void drawvertex(const vec4 &p); // Homogeneous

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
    for (int i = 0; i < (int)facelist.size(); i++) {
      if (regions[facelist[i]][3] > 0) {
        std::swap(facelist[0],facelist[i]);
        break;
      }
    }
    assert(regions[facelist[0]][3] == 1);
    // Match up successive regions.
    for (int i = 0; i < (int)facelist.size(); i++) {
      for (int j = i+1; j < (int)facelist.size(); j++) {
        if (((i%2 == 0) && regions[facelist[i]][r] == regions[facelist[j]][r]) ||
            ((i%2 == 1) && regions[facelist[i]][q] == regions[facelist[j]][q])) {
          std::swap(facelist[j],facelist[i+1]);
          break;
        }
      }
    }

#if 0
    cout << v << ": " << facelist.size();
    for (int i = 0; i < (int)facelist.size(); i++) {
      cout << " " << facelist[i];
    }
    cout << "\n";
#endif
    // Insert into the p'th facemap.
    facemap[p][v].swap(facelist);
  }
}

// For snubification, we want, for every -ve triangle, to draw a triangle between
// the Wythoff points of the surrounding +ve trianges.
// nextface: face -> 3 -> face
// wythoff: face -> point
// draw(wythoff(nextface(f,0)),wythoff(nextface(f,1)),wythoff(nextface(f,2)));
// facelist: vertex -> face list ;; The regions incident at that vertex

// Compute the distance of each face from origin
// For {2n} this is just |p[0]+p[n]|/2
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
  facecentres[k][v] = len;
}

float sgn(float x) {
  if (x < 0) return -1; else return 1;
}

void drawpointface(int k, int v, std::vector<int> &facelist, int snubface)
{
  assert(!facelist.empty());
  vector<vec3> pointlist;
  for (int i = 0; i < (int)facelist.size(); i++) {
    int n = facelist[i];
    ivec4 &f = regions[n];
    if (!snubify || f[3] != snubface) {
      vec3 &p = regionpoints[n];
      pointlist.push_back(p);
    }
  }
  int npoints = pointlist.size();
  vec3 centre = facecentres[k][v]*points[v];
  
  if (npoints <= 2) return;
  // Faces are always square on so the centre is the normal
  //setNormal(points[v]*sgn(facecentres[k][v]));
  setNormal(points[v]*sgn(facecentres[k][v]));
  vec3 offset = explode*centre;
  for (int i = 0; i < npoints; i++) {
    if (npoints == 3) {
      drawvertex(pointlist[i]+offset);
    } else {
      drawvertex(centre+offset);
      drawvertex(pointlist[i]+offset);
      drawvertex(pointlist[(i+1)%npoints]+offset);
    }
  }
}

// Point and vertex store
int getnode(const vec3 &p)
{
  // Linear scan, but the most points we have is 120
  int nnodes = points.size();
  for (int i = 0; i < nnodes; i++) {
    if (fabs(p.x - points[i].x) < eps &&
        fabs(p.y - points[i].y) < eps &&
        fabs(p.z - points[i].z) < eps) {
      return i;
    }
  }
  points.push_back(p);
  return nnodes;
}

bool addtriangle(const ivec4 &t)
{
  int ntriangles = regions.size();
  for (int i = 0; i < ntriangles; i++) {
    if (regions[i] == t) return false;
  }
  regions.push_back(t);
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

void makeschwarz(int args[6])
{
  double A,B,C,a,b,c;
  vec3 p0,q0,r0;
  A = pi*args[1]/args[0];
  B = pi*args[3]/args[2];
  C = pi*args[5]/args[4];
  cosine(A,B,C,a,b,c);
  triangle1(a,b,c,p0,q0,r0);
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

  list<ivec4> qq;
  getnode(p0);
  getnode(q0);
  getnode(r0);
  ivec4 t0(0,1,2,1);
  qq.push_back(t0);
  addtriangle(t0);
  while (!qq.empty()) {
    ivec4 t = qq.front();
    qq.pop_front();
    vec3 p = points[t[0]];
    vec3 q = points[t[1]];
    vec3 r = points[t[2]];
    int p1 = getnode(reflect(q,r,p));
    int q1 = getnode(reflect(r,p,q));
    int r1 = getnode(reflect(p,q,r));
    ivec4 t1(p1,t[1],t[2],-t[3]);
    ivec4 t2(t[0],q1,t[2],-t[3]);
    ivec4 t3(t[0],t[1],r1,-t[3]);
    if (addtriangle(t1)) qq.push_back(t1);
    if (addtriangle(t2)) qq.push_back(t2);
    if (addtriangle(t3)) qq.push_back(t3);
  }
}

vec4 invert(float r2, const vec3 &p)
{
  return vec4(p*r2,dot(p,p));
}

void drawregion(const int r) {
  ivec4 &f = regions[r];
  vec3 n = normalize(regionpoints[r]);
  setNormal(n);

  vector <vec4> plist;
  for (int i = 0; i < 3; i++) {
    plist.push_back(vec4(points[f[i]]*dfact[i],efact[i]));
  }
  //if (doclamp) clamp(plist);
  //setColour(mist);
  for (int i = 0; i < 3; i++) {
    if (f.w > 0) {
      setColour(facecolours[i]);
      drawvertex(plist[i]);
    } else {
      setColour(facecolours[2-i]);
      drawvertex(plist[2-i]);
    }
  }
}

void setsnubface(int i)
{
  vec3 triangle[3];
  for (int j = 0; j < 3; j++) {
    triangle[j] = regionpoints[adjacent[i][2-j]];
  }
  if (regions[i][3] > 0) std::swap(triangle[1],triangle[2]);
  snubsphere2 = square(length(triangle[0]+triangle[1])/2);
  vec3 norm1 = cross(triangle[0]-triangle[1],triangle[1]-triangle[2]);
  if (length(norm1) <= eps) {
    snubcentres[i] = vec3(0,0,0);
    return;
  }
  vec3 norm = normalize(norm1);
  float a = dot(triangle[0],norm); // Can probably simplify this
  snubcentres[i] = a*norm;
#if 0
  cout << triangle[0] << "\n";
  cout << triangle[1] << "\n";
  cout << triangle[2] << "\n";
  cout << norm1 << "\n";
  cout << a << "\n";
  cout << snubsphere2 << snubcentres[i] << "\n";
#endif
}

void drawsnubface(int i)
{
  vec3 triangle[3];
  vec3 centre;
  for (int j = 0; j < 3; j++) {
    triangle[j] = regionpoints[adjacent[i][2-j]];
    centre += triangle[j];
  }
  centre /= 3.0f;
  if (regions[i][3] > 0) std::swap(triangle[1],triangle[2]);

  vec3 norm = cross(triangle[0]-triangle[1],triangle[1]-triangle[2]);
  if (length(norm) < eps) return;
  norm = normalize(norm);
  // Orientation of snub faces is uncertain, I guess some faces can be
  // retrograde.
  if (dot(norm,centre) < 0) norm = -norm;
#if 0
  vec3 norm1 = normalize(centre);
  assert(fabs(dot(norm,norm1)-1) < eps);
  assert(fabs(length(triangle[0]-triangle[1]) -
              length(triangle[1]-triangle[2])) < eps);
#endif
  setNormal(norm);
  for (int j = 0; j < 3; j++) {
    drawvertex(centre+explode*centre);
    drawvertex(triangle[j]+explode*centre);
    drawvertex(triangle[(j+1)%3]+explode*centre);
  }
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
  for (int i = 0; i < 3; i++) {
    facecentres[i].resize(points.size());
  }
  snubcentres.resize(regions.size());
  for (int i = 0; i < 3; i++) {
    for (auto &j : facemap[i]) {
      setpointface(i,j.first,j.second);
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
  // Want (kp-s).s  = 0 
  // k(s.p) = s.s ie. k = (s.s)/(p.s)
  // Also scale by mm*|p| - why?
  dfact[0] = midsphere2*length(p)*length(s);
  dfact[1] = midsphere2*length(q)*length(s);
  dfact[2] = midsphere2*length(r)*length(s);
  efact[0] = dot(s,p);
  efact[1] = dot(s,q);
  efact[2] = dot(s,r);

  if (efact[0] < 0) { dfact[0] *= -1; efact[0] *= -1; }
  if (efact[1] < 0) { dfact[1] *= -1; efact[1] *= -1; }
  if (efact[2] < 0) { dfact[2] *= -1; efact[2] *= -1; }
  //printf("%g\n", midsphere2);
  drawinit();
}

vec4 hmul(float k, vec4 p)
{
  return vec4(k*p.x, k*p.y, k*p.z, p.w);
}

vec4 hadd(vec4 p, vec4 q)
{
  if (p.w > q.w) std::swap(p,q);
  q *= p.w/q.w;
  p[0] += q[0];
  p[1] += q[1];
  p[2] += q[2];
  return p;
}

void drawDual()
{
  if (!snubify) {
    if (style == 0) {
      for (int i = 0; i < (int)regions.size(); i++) {
        drawregion(i);
        //break;
      }
    } else {
      for (int i = 0; i < (int)regions.size(); i++) {
        setColour(facecolours[3%nfacecolours]);
        vector<vec4> plist;
        //float len = length(regionpoints[i]);
        int N = 3;
        for (int j = 0; j < N; j++) {
          // This should be reciprocation relative to the edge centre
          // ie. the maximum of dot(s,p[i]) - set up in initialization.
          // |r]|kr] = k|r||r| = mm
          int v = regions[i][j];
          float len = facecentres[j][v];
          plist.push_back(vec4(midsphere2*points[v],len));
        }
        //if (doclamp) clamp(plist);
        vec3 norm = normalize(regionpoints[i]);
        setNormal(norm);
        vec4 centre;
        // This will work since all homogeneous coords have same w
        for (int k = 0; k < N; k+=1) {
          centre += plist[k];
        }
        centre /= N;
        //cout << centre << "\n";
        for (int k = 0; k < N; k++) {
          //drawvertex(centre+explode*centre);
          drawvertex(plist[k]/*+explode*centre*/);
          //drawvertex(plist[(k+1)%N]+explode*centre);
        }
        //break;
      }
    }
  } else {
    int snubface = -1;
    for (int i = 0; i < (int)regions.size(); i++) {
      if (regions[i][3] == snubface) {
        vector<vec4> plist;
        for (int j = 0; j < 3; j++) {
          plist.push_back(vec4(snubsphere2*points[regions[i][j]],
                               facecentres[j][regions[i][j]]));
          vec3 &tmp = snubcentres[adjacent[i][(j+2)%3]];
          if (length(tmp) > eps) {
            plist.push_back(invert(snubsphere2,tmp));
          } else {
            assert(0);
          }
        }
        int N = plist.size();
        //if (doclamp) clamp(plist);
        vec3 norm = normalize(regionpoints[i]);
        setNormal(norm);
        assert(length(regionpoints[i]) > eps);
        vec4 centre = invert(snubsphere2,regionpoints[i]);
        vec4 offset = hmul(explode,centre);
        for (int k = 0; k < N; k++) {
          setColour(black); 
          drawvertex(hadd(centre,offset));
          setColour(facecolours[k%nfacecolours]);
          drawvertex(hadd(plist[k],offset));
          setColour(facecolours[((k+1)%N)%nfacecolours]);
          drawvertex(hadd(plist[(k+1)%N],offset));
        }
        //break;
      }
    }
  }
}

void drawMain()
{
  int snubface = 1;
  for (int i = 0; i < 3; i++) {
    if (!omit[i]) {
      setColour(facecolours[i%nfacecolours]);
      for (auto &j : facemap[i]) {
        if (snubify) {
          drawpointface(i,j.first,j.second,snubface);
          //drawpointface(i,j.first,j.second,-snubface);
        } else {
          drawpointface(i,j.first,j.second,0);
        }
      }
    }
  }
  if (snubify && !omit[3]) {
    setColour(facecolours[3%nfacecolours]);
    for (int i = 0; i < (int)regions.size(); i++) {
      if (regions[i][3] == snubface) {
        drawsnubface(i);
      }
    }
  }
}


#include <GL/glut.h>
#include <GL/gl.h>

// Drawing stuff

int args[6];

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

static double tt = 0;

// Lights settings
float amb = 0.2;
float diff = 0.2;
float spec = 0.3;
GLfloat light_ambient[]= { amb, amb, amb, 0.0f };
GLfloat light_diffuse[]= { diff, diff, diff, 0.0f };
GLfloat light_specular[]= { spec, spec, spec, 0.0f };

// Lights are at infinity
GLfloat light_position1[]= { 10.0f, 10.0f, 10.0f, 0.0f };
GLfloat light_position2[]= { -10.0f, 10.0f, 10.0f, 0.0f };

//Materials settings
GLfloat mat_specular[]= { 0, 0, 0, 0 };
GLfloat mat_shininess = 100;
GLfloat specular = 1;

bool dolog = false;

void setColour(const vec4 *colour)
{
  setColour(*colour);
}
void setColour(const vec4 &colour)
{
  if (dolog) printf("colour %g %g %g %g\n", 
                    colour[0], colour[1], colour[2], colour[3]);
  const float *p = glm::value_ptr(colour);
  glMaterialfv (GL_FRONT_AND_BACK, GL_DIFFUSE, p);
  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT, p);
}

void setNormal(const vec3 &norm)
{
  if (dolog) printf("normal %g %g %g\n", norm[0], norm[1], norm[2]);
  glNormal3fv(glm::value_ptr(norm));
}

void drawvertex(const vec4 &p)
{
  if (p[3] < 0) {
    if (dolog) printf("vertex %g %g %g %g\n", -p[0], -p[1], -p[2], -p[3]);
    glVertex4f(-size*p[0],-size*p[1],-size*p[2],-p[3]);
  } else {
    if (dolog) printf("vertex %g %g %g %g\n", p[0], p[1], p[2], p[3]);
    glVertex4f(size*p[0],size*p[1],size*p[2],p[3]);
  }
}

void drawvertex(const vec3 &p)
{
  drawvertex(vec4(p,1));
}

void drawone()
{
  glBegin(GL_TRIANGLES);
  if (dual) {
    //drawMain();
    drawDual();
  } else {
    drawMain();
  }
  glEnd();
}
void init(void)
{
  //glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

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

#if 1
  glEnable (GL_BLEND); 
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
#endif

  glShadeModel(GL_SMOOTH);
  glHint (GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
  glPolygonMode (GL_FRONT_AND_BACK, GL_FILL); // Polygon rasterization mode (polygon filled)
  glEnable(GL_DEPTH_TEST); // Enable the depth test (also called z buffer)
  //glEnable(GL_CULL_FACE); // Enable the back face culling
}

void resize (int p_width, int p_height)
{
  if (screen_width==0 && screen_height==0) exit(0);
  screen_width = p_width;
  screen_height = p_height;
  glutPostRedisplay ();
}

float ttinc = 0.0015;
vec3 abc;

void display(void)
{
  if (dolog) {
    printf ("start\n");
    printf (";; %d/%d %d/%d %d/%d\n", 
            args[0],args[1],args[2],args[3],args[4],args[5]);
    printf (";; %g %g %g%s%s\n",
            abc[0],abc[1],abc[2],
            snubify?" snub":"",
            dual?" dual":"");
  }
  rotx += rotxinc;
  roty += rotyinc;
  rotz += rotzinc;

  if (rotx >= 360) rotx -= 360;
  if (roty >= 360) roty -= 360;
  if (rotz >= 360) rotz -= 360;

  glClearColor(0, 0, 0.2, 0.0);
  glViewport(0,0,screen_width,screen_height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0f,(GLfloat)screen_width/(GLfloat)screen_height,1.0f,10000.0f);

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
  drawone();
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
  case 'l':
    dolog = true;
    break;
  case 'k':
    style = (style+1)%nstyles;
    break;
  case ']':
    // eps should be optional probably
    // avoids apparent rounding errors at singular positions
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
    rotxinc = 0;
    rotyinc = 0;
    rotzinc = 0;
    break;
  case ' ':
    pause = !pause;
    break;
  case 'v': {
    for (int i = 1; i < nfacecolours; i++) {
      std::swap(facecolours[i-1],facecolours[i]);
    }
    break;
  }
  case 's':
    snubify = !snubify;
    needparams = true; // We use different coord transform for snub figures
    break;
  case 'r': case 'R':
    filling = !filling;
    glPolygonMode (GL_FRONT_AND_BACK, filling?GL_FILL:GL_LINE);
    break;
  case 'd':
    dual = !dual;
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
  for (int i = 0; i<3; i++) {
    args[2*i] = 1;
    args[2*i+1] = 1;
    int d = sscanf(argv[1+i],"%d/%d",&args[2*i],&args[2*i+1]);
    assert(d == 1 || d == 2);
  }
  makeschwarz(args);
  // Make triangle adjacency map
  {
    int n = regions.size();
    adjacent.resize(n);
    for (int i = 0; i < n; i++) {
      int total = 0;
      for (int j = 0; j < n; j++) {
        // Find triangle that differs from face[i] only in jth vertex
        int c = 0;
        int d = 0;
        for (int k = 0; k < 3; k++) {
          if (regions[i][k] == regions[j][k]) {
            c++;
          } else {
            d = k;
          }
        }
        if (c == 2) {
          adjacent[i][d] = j;
          total++;
        }
      }
      assert(total == 3);
    }
  }
  {
    int npoints = points.size();
    for (int i = 0; i < npoints; i++) {
      makefaces(i,0,1,2);
      makefaces(i,1,2,0);
      makefaces(i,2,0,1);
    }
  }
#if 0
  interpolator.add(vec3(-0.5,1,1));
  interpolator.add(vec3(1,-0.5,1));
  interpolator.add(vec3(1,1,-0.5));
#elif 0
  interpolator.add(vec3(1,0,0));
  interpolator.add(vec3(1,1,1));
  interpolator.add(vec3(0,1,1));
  interpolator.add(vec3(-1,1,1));
  interpolator.add(vec3(1,-1,1));
  interpolator.add(vec3(1,0,1));
  interpolator.add(vec3(1,1,1));
  interpolator.add(vec3(0,1,0));
#elif 0
  interpolator.add(vec3(1,0,0));
  interpolator.add(vec3(-1,1,0));
  interpolator.add(vec3(0,0,1));
  interpolator.add(vec3(1,-1,1));
  interpolator.add(vec3(1,0,0));
  interpolator.add(vec3(1,0,-1));
  interpolator.add(vec3(0,0,1));
  interpolator.add(vec3(-1,1,1));
  interpolator.add(vec3(1,-1,1));
#else
  interpolator.add(vec3(1,1,0));
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
#endif
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