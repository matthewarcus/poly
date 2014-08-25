#if !defined UTILS_H
#define UTILS_H
template <typename T>
T max3(T t1, T t2, T t3) {
  return std::max(t1,std::max(t2,t3));
}

template <typename T>
T min3(T t1, T t2, T t3) {
  return std::min(t1,std::min(t2,t3));
}
template <typename T>
T mid3(T t1, T t2, T t3) {
  if (t1 > t2) std::swap(t1,t2);
  if (t2 > t3) std::swap(t2,t3);
  if (t1 > t2) std::swap(t1,t2);
  return t2;
}
template <typename T>
T square(T t) { return t * t; }

std::ostream &operator<<(std::ostream &str, const glm::vec3 &a)
{
   str << "[ " << a.x << ", " << a.y << ", " << a.z << " ]";
   return str;
}

glm::vec4 black( 0, 0, 0, 1 );
glm::vec4 white( 1, 1, 1, 1 );
glm::vec4 red( 1, 0, 0, 1 );
glm::vec4 pink( 1, 0.5, 0.5, 1 );
glm::vec4 grey( 0.5, 0.5, 0.5, 1 );
glm::vec4 orange( 1, 0.5, 0, 1 );
glm::vec4 green( 0, 1, 0, 1 );
glm::vec4 blue( 0, 0, 1, 1 );
glm::vec4 yellow( 1, 1, 0, 1 );
glm::vec4 cyan( 0, 1, 1, 1 );
glm::vec4 magenta( 1, 0, 1, 1 );
glm::vec4 purple( 1, 0.5, 1, 1 );
glm::vec4 mist( 1, 1, 1, 0.3 );

// Generate series of values in a range.
template<typename T>
class Interpolator
{
public:
  size_t size() {
    return a.size();
  }
  void add(const T &t) {
    a.push_back(t);
  }
  T operator()(double t) {
    int k = floor(t);
    float r = t - k;
    return (1-r)*a[k%a.size()] + r*a[(k+1)%a.size()];
  }
private:
  std::vector<T> a;
};

template <typename F>
bool adjust(float &a, float delta, F f)
{
   float t = f(a);
   bool found = false;
   while (true) {
     float a1 = a+delta;
     if (a > 100 || a == a1) return false;
     float t1 = f(a1);
     if (t1 >= t) break;
     a = a1; t = t1; found = true;
   }
   return found;
}

// Scan 2 dimensions, but using all 3 trilinear coords gives better
// results (sometimes gets stuck in local minimum with just 2).
template <typename G>
void find(float &a, float &b, float &c, float delta, int n, G g)
{
   for (int i = 0; i < n && delta != 0; i++, delta/=10) {
     while(true) {
       int changed = 0;
       auto gbc = [&](float a1){ return g(a1,b,c); };
       changed += adjust(a, delta, gbc) || adjust(a, -delta, gbc);
       auto gca = [&](float b1){ return g(a,b1,c); };
       changed += adjust(b, delta, gca) || adjust(b, -delta, gca);
       auto gab = [&](float c1){ return g(a,b,c1); };
       changed += adjust(c, delta, gab) || adjust(c, -delta, gab);
       if (!changed) break;
     }
   }
}

void clamp(std::vector<glm::vec3> &plist)
{
  int N = plist.size();
  float maxlen = 0;
  for (int k = 0; k < N; k+=1) {
    float len = glm::length(plist[k]);
    if (len > maxlen) maxlen = len;
  }
  if (maxlen > 100) maxlen = 100;
  for (int k = 0; k < N; k+=1) {
    plist[k] /= maxlen;
  }
}

#endif
