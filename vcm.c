#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define CAST(t, x) (t)(x)
#define ALLOC(t, n) malloc(sizeof(t) * (n))
#define PI 3.141592653589793
#define EPS 1.0e-6
#define TIME_TRI 1.0
#define TIME_BOX 1.0
#define SEED 0
#define DEST_PATH "o.pfm"

struct bvh_node;
typedef struct bvh_node BVH_NODE;

typedef struct {
  double x, y, z;
} IO_VEC;
typedef struct {
  double xx, xy, xz, yx, yy, yz, zx, zy, zz;
} IO_MAT;
typedef struct {
  int t;
  IO_VEC a;
} IO_EDF;
typedef struct {
  int t;
  IO_VEC a;
} IO_BSDF;
typedef struct {
  double eta;
  IO_VEC a;
} IO_MEDIUM;
typedef struct {
  int edf, bsdf, mp, mn, v0, v1, v2;
} IO_FACE;
typedef struct {
  IO_VEC t;
  IO_MAT m;
  int num_edfs, num_bsdfs, num_media, num_verts, num_faces;
  IO_EDF *edfs;
  IO_BSDF *bsdfs;
  IO_MEDIUM *media;
  IO_VEC *verts;
  IO_FACE *faces;
} IO_SCENE;

typedef struct {
  uint64_t s0, s1;
} PRNG_STATE;
typedef struct {
  double x, y, z;
} VEC;
typedef struct {
  double xx, xy, xz, yx, yy, yz, zx, zy, zz;
} MAT;
typedef struct {
  int n, w, h;
  double *data;
} IMAGE;
typedef struct {
  int t;
  VEC a;
} EDF;
typedef struct {
  int t;
  VEC a;
} BSDF;
typedef struct {
  double eta;
  VEC a;
} MEDIUM;
typedef struct {
  EDF *edf;
  BSDF *bsdf;
  MEDIUM *mp, *mn;
  VEC *v0, *v1, *v2;
  double p;
} FACE;
typedef struct {
  FACE *f;
  bool c;
  double pf, pb;
  VEC a, x, n, wi;
} VERTEX;
typedef struct {
  int n;
  VERTEX *vs;
} PATH;
typedef struct {
  PATH *p;
  int j;
} KDT_NODE;
typedef struct {
  VEC min, max;
} BOX;
struct bvh_node {
  BOX box;
  BVH_NODE *l, *r;
  int num_faces;
  FACE *faces;
};
typedef struct {
  int num_threads, num_samples, num_verts, num_iters, num_bvh_nodes, num_phots, num_kdt_nodes, wi, hi;
  double ef, wf, hf, df, r, f, ir;
} CONF;
typedef struct {
  int num_edfs, num_bsdfs, num_media, num_verts, num_faces;
  EDF *edfs;
  BSDF *bsdfs;
  MEDIUM *media;
  VEC *verts;
  FACE *faces;
  int num_lights;
  FACE **lights;
  double *cdf;
  BVH_NODE *bvh_nodes;
} SCENE;
typedef struct {
	int i;
	double x, y;
	VEC c;
} IMG_SPL;

void io_free_scene(IO_SCENE ioscene) {
  free(ioscene.edfs);
  free(ioscene.bsdfs);
  free(ioscene.media);
  free(ioscene.verts);
  free(ioscene.faces);
}

IO_SCENE io_read_scene(char *path) {
  FILE *f = fopen(path, "r");
  if (f == NULL) {
    printf("cannot open file: %s\n", path);
    exit(0);
  }

  int buf_size = 256, id_size = 256;
  char buf[buf_size], id[id_size];

  IO_SCENE scene;

  while (true) {
    if (fgets(buf, buf_size, f) == NULL) {
      break;
    }
    sscanf(buf, "%s ", id);
    if (strcmp(id, "translate") == 0) {
      IO_VEC t;
      sscanf(buf, "%*s %lf %lf %lf ", &t.x, &t.y, &t.z);
      scene.t = t;
    } else if (strcmp(id, "transform") == 0) {
      IO_MAT m;
      sscanf(buf, "%*s %lf %lf %lf %lf %lf %lf %lf %lf %lf ", &m.xx, &m.xy, &m.xz, &m.yx, &m.yy, &m.yz, &m.zx, &m.zy, &m.zz);
      scene.m = m;
    } else if (strcmp(id, "edfs") == 0) {
      sscanf(buf, "%*s %i ", &scene.num_edfs);
      scene.edfs = ALLOC(IO_EDF, scene.num_edfs);
      for (int i = 0; i < scene.num_edfs; ++i) {
        IO_EDF edf;
        fgets(buf, buf_size, f);
        sscanf(buf, "%i %lf %lf %lf ", &edf.t, &edf.a.x, &edf.a.y, &edf.a.z);
        scene.edfs[i] = edf;
      }
    } else if (strcmp(id, "bsdfs") == 0) {
      sscanf(buf, "%*s %i ", &scene.num_bsdfs);
      scene.bsdfs = ALLOC(IO_BSDF, scene.num_bsdfs);
      for (int i = 0; i < scene.num_bsdfs; ++i) {
        IO_BSDF bsdf;
        fgets(buf, buf_size, f);
        sscanf(buf, "%i %lf %lf %lf ", &bsdf.t, &bsdf.a.x, &bsdf.a.y, &bsdf.a.z);
        scene.bsdfs[i] = bsdf;
      }
    } else if (strcmp(id, "media") == 0) {
      sscanf(buf, "%*s %i ", &scene.num_media);
      scene.media = ALLOC(IO_MEDIUM, scene.num_media);
      for (int i = 0; i < scene.num_media; ++i) {
        IO_MEDIUM medium;
        fgets(buf, buf_size, f);
        sscanf(buf, "%lf %lf %lf %lf ", &medium.eta, &medium.a.x, &medium.a.y, &medium.a.z);
        scene.media[i] = medium;
      }
    } else if (strcmp(id, "verticies") == 0) {
      sscanf(buf, "%*s %i ", &scene.num_verts);
      scene.verts = ALLOC(IO_VEC, scene.num_verts);
      for (int i = 0; i < scene.num_verts; ++i) {
        IO_VEC vert;
        fgets(buf, buf_size, f);
        sscanf(buf, "%lf %lf %lf ", &vert.x, &vert.y, &vert.z);
        scene.verts[i] = vert;
      }
    } else if (strcmp(id, "faces") == 0) {
      sscanf(buf, "%*s %i ", &scene.num_faces);
      scene.faces = ALLOC(IO_FACE, scene.num_faces);
      for (int i = 0; i < scene.num_faces; ++i) {
        IO_FACE face;
        fgets(buf, buf_size, f);
        sscanf(buf, "%i %i %i %i %i %i %i ", &face.edf, &face.bsdf, &face.mp, &face.mn, &face.v0, &face.v1, &face.v2);
        scene.faces[i] = face;
      }
    }
  }

  fclose(f);
  return scene;
}

double sign(double x) {
  return x < 0.0 ? -1.0 : x > 0.0 ? 1.0 : 0.0;
}
double sqr(double x) {
  return x * x;
}
double cube(double x) {
  return x * x * x;
}
double quad(double x) {
  return x * x * x * x;
}
double opst(double x) {
  return sqrt(fmax(0.0, 1.0 - sqr(x)));
}

double time_get() {
  return clock() / CAST(double, CLOCKS_PER_SEC);
}

//splitmix64
uint64_t sm64_next(uint64_t *s) {
  uint64_t x = *s += 0x9E3779B97F4A7C15;
  x = 0xBF58476D1CE4E5B9 * (x ^ (x >> 30));
  x = 0x94D049BB133111EB * (x ^ (x >> 27));
  x = x ^ (x >> 31);
  return x;
}
//xoroshiro128+
uint64_t rotl(uint64_t a, uint64_t b) {
  return (a << b) | (a >> (64 - b));
}
uint64_t prng_next(PRNG_STATE *s) {
  PRNG_STATE a = *s;
  uint64_t ret = a.s0 + a.s1;
  a.s1 ^= a.s0;
  s->s0 = rotl(a.s0, 55) ^ a.s1 ^ (a.s1 << 14);
  s->s1 = rotl(a.s1, 36);
  return ret;
}
double prng_db(PRNG_STATE *s) {
  return 0x1.0p-53 * (prng_next(s) >> 11);
}
PRNG_STATE prng_seed(uint64_t seed) {
  PRNG_STATE s = {sm64_next(&seed), sm64_next(&seed)};
  return s;
}
PRNG_STATE prng_jump(PRNG_STATE s) {
  PRNG_STATE jump = {0xBEAC0467EBA5FACB, 0xD86B048B86AA9922};
  PRNG_STATE next = {0, 0};
  for (int i = 0; i < 64; ++i) {
    if (jump.s0 & 1 << i) {
      next.s0 ^= s.s0;
      next.s1 ^= s.s1;
    }
    prng_next(&s);
  }
  for (int i = 0; i < 64; ++i) {
    if (jump.s1 & 1 << i) {
      next.s0 ^= s.s0;
      next.s1 ^= s.s1;
    }
    prng_next(&s);
  }
  return next;
}

VEC vec_init(double x, double y, double z) {
  VEC v = {x, y, z};
  return v;
}
VEC vec_inita(double a) {
  return vec_init(a, a, a);
}
VEC vec_add(VEC a, VEC b) {
  return vec_init(a.x + b.x, a.y + b.y, a.z + b.z);
}
VEC vec_sub(VEC a, VEC b) {
  return vec_init(a.x - b.x, a.y - b.y, a.z - b.z);
}
VEC vec_mul(VEC a, VEC b) {
  return vec_init(a.x * b.x, a.y * b.y, a.z * b.z);
}
VEC vec_scale(double a, VEC b) {
  return vec_init(a * b.x, a * b.y, a * b.z);
}
VEC vec_cross(VEC a, VEC b) {
  return vec_init(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}
double vec_dot(VEC a, VEC b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}
double vec_len(VEC a) {
  return sqrt(vec_dot(a, a));
}
VEC vec_norm(VEC a) {
  return vec_scale(1.0 / vec_len(a), a);
}
VEC vec_exp(VEC a) {
  return vec_init(exp(a.x), exp(a.y), exp(a.z));
}
VEC vec_min(VEC a, VEC b) {
  return vec_init(fmin(a.x, b.x), fmin(a.y, b.y), fmin(a.z, b.z));
}
VEC vec_max(VEC a, VEC b) {
  return vec_init(fmax(a.x, b.x), fmax(a.y, b.y), fmax(a.z, b.z));
}
bool vec_eq(VEC a, VEC b) {
  double eps = EPS;
  VEC d = vec_sub(a, b);
  return fabs(d.x) <= eps && fabs(d.y) <= eps && fabs(d.z) <= eps;
}
void vec_pr(VEC a, FILE *o) {
  fprintf(o, "%f %f %f\n", a.x, a.y, a.z);
}
MAT mat_init(double xx, double xy, double xz, double yx, double yy, double yz, double zx, double zy, double zz) {
  MAT a = {xx, xy, xz, yx, yy, yz, zx, zy, zz};
  return a;
}
MAT mat_initb(VEC x, VEC y, VEC z) {
  return mat_init(x.x, x.y, x.z, y.x, y.y, y.z, z.x, z.y, z.z);
}
MAT mat_scale(double a, MAT b) {
  return mat_init(a * b.xx, a * b.xy, a * b.xz, a * b.yx, a * b.yy, a * b.yz, a * b.zx, a * b.zy, a * b.zz);
}
VEC mat_trafo(MAT a, VEC b) {
  return vec_init(
  a.xx * b.x + a.yx * b.y + a.zx * b.z,
  a.xy * b.x + a.yy * b.y + a.zy * b.z,
  a.xz * b.x + a.yz * b.y + a.zz * b.z);
}
MAT mat_inv(MAT a) {
  MAT b = mat_init(
  a.yy * a.zz - a.yz * a.zy, a.xz * a.zy - a.xy * a.zz, a.xy * a.yz - a.xz * a.yy,
  a.yz * a.zx - a.yx * a.zz, a.xx * a.zz - a.xz * a.zx, a.xz * a.yx - a.xx * a.yz,
  a.yx * a.zy - a.yy * a.zx, a.xy * a.zx - a.xx * a.zy, a.xx * a.yy - a.xy * a.yx);
  double det = a.xx * b.xx + a.xy * b.yx + a.xz * b.zx;
  return mat_scale(1.0 / det, b);
}
MAT mat_trapo(MAT a) {
  return mat_init(a.xx, a.yx, a.zx, a.xy, a.yy, a.zy, a.xz, a.yz, a.zz);
}
MAT tan_mat(VEC z) {
  VEC x, y;
  if (fabs(z.x) > fabs(z.y)) {
    x = vec_scale(1.0 / sqrt(sqr(z.x) + sqr(z.z)), vec_init(-z.z, 0.0, z.x));
  } else {
    x = vec_scale(1.0 / sqrt(sqr(z.y) + sqr(z.z)), vec_init(0.0, z.z, -z.y));
  }
  y = vec_cross(z, x);
  return mat_initb(x, y, z);
}
VEC spherical_coord(double cost, double sint, double cosp, double sinp) {
  return vec_init(sint * cosp, sint * sinp, cost);
}
VEC refl(VEC n, VEC d, double cos) {
  return vec_sub(d, vec_scale(2.0 * cos, n));
}
VEC refr(VEC n, VEC d, double eta, double cosi, double cost) {
  return vec_sub(vec_scale(eta, d), vec_scale(eta * cosi - cost, n));
}

double luminance(VEC a) {
  return vec_dot(vec_init(0.2126, 0.7152, 0.0722), a);
}

BOX box_default() {
  BOX box = {vec_inita(INFINITY), vec_inita(-INFINITY)};
  return box;
}
double box_area(BOX box) {
  VEC d = vec_sub(box.max, box.min);
  return 2.0 * (d.x * d.y + d.y * d.z + d.z * d.x);
}
VEC box_center(BOX box) {
  return vec_scale(0.5, vec_add(box.min, box.max));
}
BOX box_merge(BOX a, BOX b) {
  BOX c = {vec_min(a.min, b.min), vec_max(a.max, b.max)};
  return c;
}
bool box_isect(BOX box, VEC o, VEC d) {
  double eps = EPS;
  double tmin = -INFINITY, tmax = INFINITY;
  double tx0 = (box.min.x - o.x) / d.x;
  double tx1 = (box.max.x - o.x) / d.x;
  tmin = fmax(tmin, fmin(tx0, tx1));
  tmax = fmin(tmax, fmax(tx0, tx1));
  if (tmax < tmin) {
	  return false;
  }
  double ty0 = (box.min.y - o.y) / d.y;
  double ty1 = (box.max.y - o.y) / d.y;
  tmin = fmax(tmin, fmin(ty0, ty1));
  tmax = fmin(tmax, fmax(ty0, ty1));
  if (tmax < tmin) {
	  return false;
  }
  double tz0 = (box.min.z - o.z) / d.z;
  double tz1 = (box.max.z - o.z) / d.z;
  tmin = fmax(tmin, fmin(tz0, tz1));
  tmax = fmin(tmax, fmax(tz0, tz1));
  if (tmax < tmin) {
	  return false;
  }
  return tmax > eps;
}
BOX tri_box(VEC v0, VEC v1, VEC v2) {
  BOX box = box_default();
  box.min = vec_min(box.min, v0);
  box.min = vec_min(box.min, v1);
  box.min = vec_min(box.min, v2);
  box.max = vec_max(box.max, v0);
  box.max = vec_max(box.max, v1);
  box.max = vec_max(box.max, v2);
  return box;
}
double tri_area(VEC v0, VEC v1, VEC v2) {
  return 0.5 * vec_len(vec_cross(vec_sub(v1, v0), vec_sub(v2, v0)));
}
double tri_pdf(VEC v0, VEC v1, VEC v2) {
  return 1.0 / tri_area(v0, v1, v2);
}
VEC tri_sample(VEC v0, VEC v1, VEC v2, PRNG_STATE *s) {
  double u0, u1;
  u0 = prng_db(s);
  u1 = prng_db(s);
  double a = sqrt(u0);
  double u = 1.0 - a;
  double v = u1 * a;
  return vec_add(vec_add(vec_scale(1.0 - u - v, v0), vec_scale(u, v1)), vec_scale(v, v2));
}
VEC tri_norm(VEC v0, VEC v1, VEC v2) {
  return vec_norm(vec_cross(vec_sub(v1, v0), vec_sub(v2, v0)));
}
bool tri_isect(VEC v0, VEC v1, VEC v2, VEC o, VEC d, double *t) {
  double eps = EPS;
  VEC e0 = vec_sub(v1, v0);
  VEC e1 = vec_sub(v2, v0);
  VEC pv = vec_cross(d, e1);
  double det = vec_dot(e0, pv);
  if (fabs(det) < eps) {
    return false;
  }
  VEC tv = vec_sub(o, v0);
  double u = vec_dot(tv, pv) / det;
  if (u < 0.0 || u > 1.0) {
    return false;
  }
  VEC qv = vec_cross(tv, e0);
  double v = vec_dot(d, qv) / det;
  if (v < 0.0 || u + v > 1.0) {
    return false;
  }
  *t = vec_dot(e1, qv) / det;
  return *t > eps;
}

double disc_area(double r) {
  return PI * sqr(r);
}
double disc_pdf(double r) {
  return 1.0 / disc_area(r);
}
VEC disc_sample(double r, PRNG_STATE *s) {
  double u = r * sqrt(prng_db(s));
  double v = 2.0 * PI * prng_db(s);
  return vec_scale(u, vec_init(cos(v), sin(v), 0.0));
}

IMAGE image_init(int n, int w, int h) {
  IMAGE img = {n, w, h, ALLOC(double, n * w * h)};
  return img;
}
void image_del(IMAGE img) {
  free(img.data);
}
void image_fill(IMAGE img, double x) {
  for (int i = 0; i < img.n * img.w * img.h; ++i) {
    img.data[i] = x;
  }
}
void image_save_pfm(IMAGE img, char *path) {
  FILE *f = fopen(path, "wb");
  if (f == NULL) {
    printf("cannot open file: %s\n", path);
    exit(0);
  }
  fprintf(f, "PF %d %d -1 ", img.w, img.h);
  float *data = ALLOC(float, img.n * img.w * img.h);
  for (int iy = 0; iy < img.h; ++iy) {
    for (int ix = 0; ix < img.w; ++ix) {
      for (int i = 0; i < img.n; ++i) {
        data[img.n * (ix + img.w * iy) + i] = CAST(float, img.data[img.n * (ix + img.w * (img.h - iy - 1)) + i]);
      }
    }
  }
  fwrite(data, sizeof(float), img.n * img.w * img.h, f);
  free(data);
}

double fresnel_dielectric(double eta, double cosi, double cost) {
  double perp = (eta * cosi - cost) / (eta * cosi + cost);
  double paral = (eta * cost - cosi) / (eta * cost + cosi);
  return (sqr(perp) + sqr(paral)) / 2.0;
}
double fresnel_conductor(double eta, double k, double cosi) {
  cosi = fabs(cosi);
  double s = 1.0 - sqr(cosi);
  double i = sqr(eta) - sqr(k) - s;
  double w = sqrt(sqr(i) + sqr(2.0 * eta * k));
  double a = sqrt(fmax(0.0, (w + i) / 2.0));
  double perp = (w + sqr(cosi) - 2.0 * a * cosi) / (w + sqr(cosi) + 2.0 * a * cosi);
  double paral = (sqr(cosi) * w + sqr(s) - 2.0 * a * cosi * s) / (sqr(cosi) * w + sqr(s) + 2.0 * a * cosi * s);
  return (perp + perp * paral) / 2.0;
}

double hemisphere_cosine_weighted_pdf(VEC n, VEC wo) {
  double cos = vec_dot(n, wo);
  return cos > 0.0 ? cos / PI : 0.0;
}
VEC hemisphere_cosine_weighted_sample(VEC n, PRNG_STATE *s) {
  double u0, u1;
  u0 = prng_db(s);
  u1 = prng_db(s);
  double phi = 2.0 * PI * u1;
  double cost, sint, cosp, sinp;
  cost = sqrt(u0);
  sint = opst(cost);
  cosp = cos(phi);
  sinp = sin(phi);
  MAT m = tan_mat(n);
  return mat_trafo(m, spherical_coord(cost, sint, cosp, sinp));
}
double cos_pow_pdf(double a, VEC n, VEC wo) {
  double cos = vec_dot(n, wo);
  return cos > 0.0 ? (2.0 + a) * pow(cos, 1.0 + a) / (2.0 * PI) : 0.0;
}
VEC cos_pow_sample(double a, VEC n, PRNG_STATE *s) {
  double u0, u1;
  u0 = prng_db(s);
  u1 = prng_db(s);
  double phi = 2.0 * PI * u1;
  double cost, sint, cosp, sinp;
  cost = pow(1.0 - u0, 1.0 / (2.0 + a));
  sint = opst(cost);
  cosp = cos(phi);
  sinp = sin(phi);
  MAT m = tan_mat(n);
  return mat_trafo(m, spherical_coord(cost, sint, cosp, sinp));
}
int discrete_sample(int n, double *cdf, PRNG_STATE *s) {
  double x = prng_db(s);
  for(int i = 0; i < n; ++i) {
    if (x < cdf[i]) {
      return i;
    }
  }
  return -1;
}

VEC edf_eval(EDF edf, VEC n, VEC wo) {
  switch (edf.t) {
    case 0 : {
      return vec_inita(0.0);
    }
    case 1 : {
      double cos = vec_dot(n, wo);
      return vec_scale(cos > 0.0 ? 1.0 / PI : 0.0, edf.a);
    }
    case 2 : {
      double a = 65536.0;
      double cos = vec_dot(n, wo);
      return vec_scale(cos > 0.0 ? (2.0 + a) * pow(cos, 1.0 + a) / (2.0 * PI) : 0.0, edf.a);
    }
    case 3 : {
      return vec_scale(vec_eq(wo, n) ? 1.0 : 0.0, edf.a);
    }
    default : {
      return vec_inita(0.0);
    }
  }
}
double edf_pdf(EDF edf, VEC n, VEC wo) {
  switch (edf.t) {
    case 1 : {
      return hemisphere_cosine_weighted_pdf(n, wo);
    }
    case 2 : {
      double a = 65536.0;
      return cos_pow_pdf(a, n, wo);
    }
    case 3 : {
      return vec_eq(wo, n) ? 1.0 : 0.0;
    }
    default : {
      return 0.0;
    }
  }
}
VEC edf_sample(EDF edf, VEC n, PRNG_STATE *s) {
  switch (edf.t) {
    case 1 : {
      return hemisphere_cosine_weighted_sample(n, s);
    }
    case 2 : {
      double a = 65536.0;
      return cos_pow_sample(a, n, s);
    }
    case 3 : {
      return n;
    }
    default : {
      return vec_inita(0.0);
    }
  }
}
bool edf_connect(EDF edf) {
  switch (edf.t) {
    case 0 : {
      return false;
    }
    case 1 : {
      return true;
    }
    case 2 : {
      return true;
    }
    case 3 : {
      return false;
    }
    default : {
      return false;
    }
  }
}

VEC bsdf_eval(BSDF bsdf, double etai, double etat, VEC n, VEC wi, VEC wo) {
  switch (bsdf.t) {
    case 0 : {
      double cosi = vec_dot(n, wi);
      return vec_scale(vec_eq(wi, wo) ? 1.0 / fabs(cosi) : 0.0, bsdf.a);
    }
    case 1 : {
      double cosi = vec_dot(n, wi), coso = vec_dot(n, wo);
      return vec_scale(cosi * coso < 0.0 ? 1.0 / PI : 0.0, bsdf.a);
    }
    case 2 : {
      double m = 2.7732, k = 2.9278;
      double cosi = vec_dot(n, wi);
      VEC r = refl(n, wi, cosi);
      return vec_scale(vec_eq(wo, r) ? fresnel_conductor(etai / m, k, cosi) / fabs(cosi) : 0.0, bsdf.a);
    }
    case 3 : {
      double cosi = vec_dot(n, wi);
      double eta = etai / etat;
      double eta2 = sqr(eta);
      double det = 1.0 - eta2 * (1.0 - sqr(cosi));
      if (det > 0.0) {
        double cost = sign(cosi) * sqrt(det);
        double fr = fresnel_dielectric(eta, cosi, cost);
        VEC r = refl(n, wi, cosi);
        if (vec_eq(r, wo)) {
          return vec_scale(fr / fabs(cosi), bsdf.a);
        }
        VEC t = refr(n, wi, eta, cosi, cost);
        if (vec_eq(t, wo)) {
          return vec_scale((1.0 - fr) / (eta2 * fabs(cosi)), bsdf.a);
        }
      }
      VEC r = refl(n, wi, cosi);
      return vec_scale(vec_eq(r, wo) ? 1.0 / fabs(cosi) : 0.0, bsdf.a);
    }
    default : {
      return vec_inita(0.0);
    }
  }
}
double bsdf_pdf(BSDF bsdf, double etai, double etat, VEC n, VEC wi, VEC wo) {
  switch (bsdf.t) {
    case 0 : {
      return vec_eq(wi, wo) ? 1.0 : 0.0;
    }
    case 1 : {
      double cosi = vec_dot(n, wi);
      return hemisphere_cosine_weighted_pdf(vec_scale(-sign(cosi), n), wo);
    }
    case 2 : {
      double cosi = vec_dot(n, wi);
      VEC r = refl(n, wi, cosi);
      return vec_eq(wo, r) ? 1.0 : 0.0;
    }
    case 3 : {
      double cosi = vec_dot(n, wi);
      double eta = etai / etat;
      double eta2 = sqr(eta);
      double det = 1.0 - eta2 * (1.0 - sqr(cosi));
      if (det > 0.0) {
        double cost = sign(cosi) * sqrt(det);
        double fr = fresnel_dielectric(eta, cosi, cost);
        VEC r = refl(n, wi, cosi);
        if (vec_eq(r, wo)) {
          return fr;
        }
        VEC t = refr(n, wi, eta, cosi, cost);
        if (vec_eq(t, wo)) {
          return 1.0 - fr;
        }
      }
      VEC r = refl(n, wi, cosi);
      return vec_eq(r, wo) ? 1.0 : 0.0;
    }
    default : {
      return 0.0;
    }
  }
}
VEC bsdf_sample(BSDF bsdf, double etai, double etat, VEC n, VEC wi, PRNG_STATE *s) {
  switch (bsdf.t) {
    case 0 : {
      return wi;
    }
    case 1 : {
      double cosi = vec_dot(n, wi);
      return hemisphere_cosine_weighted_sample(vec_scale(-sign(cosi), n), s);
    }
    case 2 : {
      double cosi = vec_dot(n, wi);
      return refl(n, wi, cosi);
    }
    case 3 : {
      double cosi = vec_dot(n, wi);
      double eta = etai / etat;
      double eta2 = sqr(eta);
      double det = 1.0 - eta2 * (1.0 - sqr(cosi));
      if (det > 0.0) {
        double cost = sign(cosi) * sqrt(det);
        double fr = fresnel_dielectric(eta, cosi, cost);
        if (prng_db(s) >= fr) {
          return refr(n, wi, eta, cosi, cost);
        }
      }
      return refl(n, wi, cosi);
    }
    default : {
      return vec_inita(0.0);
    }
  }
}
bool bsdf_connect(BSDF bsdf) {
  switch (bsdf.t) {
    case 0 : {
      return false;
    }
    case 1 : {
      return true;
    }
    case 2 : {
      return false;
    }
    case 3 : {
      return false;
    }
    default : {
      return false;
    }
  }
}

MEDIUM medium_init(double eta, VEC a) {
  MEDIUM m = {eta, a};
  return m;
}
VEC transmittance(VEC tr, double t) {
  return vec_exp(vec_scale(-t, tr));
}
void media_sides(double cos, MEDIUM mp, MEDIUM mn, MEDIUM *mr, MEDIUM *mt) {
  bool into = cos < 0.0;
  *mr = into ? mp : mn;
  *mt = into ? mn : mp;
}
void scene_prep(SCENE *scene) {
  scene->num_lights = 0;
  for (int i = 0; i < scene->num_faces; ++i) {
    if (scene->faces[i].edf->t != 0) {
      ++scene->num_lights;
    }
  }
  scene->lights = ALLOC(FACE *, scene->num_lights);
  scene->cdf = ALLOC(double, scene->num_lights);
  double *pdf = ALLOC(double, scene->num_lights);
  int num_lights = 0;
  for (int i = 0; i < scene->num_faces; ++i) {
    FACE *face = scene->faces + i;
    if (face->edf->t != 0) {
      scene->lights[num_lights] = face;
      ++num_lights;
    }
  }
  double sum = 0.0;
  for (int i = 0; i < scene->num_lights; ++i) {
    FACE face = *scene->lights[i];
    double a = luminance(face.edf->a) * tri_area(*face.v0, *face.v1, *face.v2);
    sum += a;
    pdf[i] = a;
    scene->cdf[i] = sum;
  }
  for (int i = 0; i < scene->num_lights; ++i) {
    pdf[i] /= sum;
    scene->cdf[i] /= sum;
  }
  for (int i = 0; i < scene->num_lights; ++i) {
    scene->lights[i]->p = pdf[i];
  }
  //for (int i = 0; i < scene->num_lights; ++i) {
  //  printf("%i/%i: %f %f\n", i, scene->num_lights, pdf[i], scene->cdf[i]);
  //}
  free(pdf);
}
int bvh_cmp_x(void const *A, void const *B) {
  FACE const *fa = A;
  FACE const *fb = B;
  BOX ba = tri_box(*fa->v0, *fa->v1, *fa->v2);
  BOX bb = tri_box(*fb->v0, *fb->v1, *fb->v2);
  VEC ca = box_center(ba);
  VEC cb = box_center(bb);
  return ca.x < cb.x ? -1 : ca.x > cb.x ? 1 : 0;
}
int bvh_cmp_y(void const *A, void const *B) {
  FACE const *fa = A;
  FACE const *fb = B;
  BOX ba = tri_box(*fa->v0, *fa->v1, *fa->v2);
  BOX bb = tri_box(*fb->v0, *fb->v1, *fb->v2);
  VEC ca = box_center(ba);
  VEC cb = box_center(bb);
  return ca.y < cb.y ? -1 : ca.y > cb.y ? 1 : 0;
}
int bvh_cmp_z(void const *A, void const *B) {
  FACE const *fa = A;
  FACE const *fb = B;
  BOX ba = tri_box(*fa->v0, *fa->v1, *fa->v2);
  BOX bb = tri_box(*fb->v0, *fb->v1, *fb->v2);
  VEC ca = box_center(ba);
  VEC cb = box_center(bb);
  return ca.z < cb.z ? -1 : ca.z > cb.z ? 1 : 0;
}
BVH_NODE *bvh_make_recur(int num_faces, FACE *faces, double *al, double *ar, BVH_NODE *nodes, int *num_nodes) {
  //if (num_faces == 0) {
  //  return NULL;
  //}
  double time_tri = TIME_TRI, time_box = TIME_BOX;
  BOX box = box_default();
  for (int i = 0; i < num_faces; ++i) {
    FACE face = faces[i];
    box = box_merge(box, tri_box(*face.v0, *face.v1, *face.v2));
  }
  double area = box_area(box);
  double cost = time_tri * num_faces;
  int axis = -1, idx = -1;
  for (int i = 0; i < 3; ++i) {
    switch (i) {
      case 0 : {
        qsort(faces, num_faces, sizeof(FACE), bvh_cmp_x);
      } break;
      case 1 : {
        qsort(faces, num_faces, sizeof(FACE), bvh_cmp_y);
      } break;
      case 2 : {
        qsort(faces, num_faces, sizeof(FACE), bvh_cmp_z);
      } break;
    }
    BOX bl = box_default();
    al[0] = 0.0;
    for (int j = 0; j < num_faces; ++j) {
      FACE face = faces[j];
      bl = box_merge(bl, tri_box(*face.v0, *face.v1, *face.v2));
      al[j + 1] = box_area(bl);
    }
    BOX br = box_default();
    ar[num_faces] = 0.0;
    for (int j = num_faces - 1; j >= 0; --j) {
      FACE face = faces[j];
      br = box_merge(br, tri_box(*face.v0, *face.v1, *face.v2));
      ar[j] = box_area(br);
    }
    for (int j = 0; j < num_faces + 1; ++j) {
      double cost_tmp = 2.0 * time_box + (al[j] * j + ar[j] * (num_faces - j)) * time_tri / area;
      if (cost_tmp < cost) {
        cost = cost_tmp;
        axis = i;
        idx = j;
      }
    }
  }
  if (axis == -1) {
    BVH_NODE *node = nodes + *num_nodes;
    ++*num_nodes;
    node->box = box;
    node->l = NULL;
    node->r = NULL;
    node->num_faces = num_faces;
    node->faces = faces;
    return node;
  } else {
    switch (axis) {
      case 0 : {
        qsort(faces, num_faces, sizeof(FACE), bvh_cmp_x);
      } break;
      case 1 : {
        qsort(faces, num_faces, sizeof(FACE), bvh_cmp_y);
      } break;
      case 2 : {
        qsort(faces, num_faces, sizeof(FACE), bvh_cmp_z);
      } break;
    }
    BVH_NODE *node = nodes + *num_nodes;
    ++*num_nodes;
    node->box = box;
    node->l = bvh_make_recur(idx, faces, al, ar, nodes, num_nodes);
    node->r = bvh_make_recur(num_faces - idx, faces + idx, al, ar, nodes, num_nodes);
    node->num_faces = 0;
    node->faces = NULL;
    return node;
  }
}
void bvh_make(CONF conf, SCENE *scene) {
  scene->bvh_nodes = ALLOC(BVH_NODE, conf.num_bvh_nodes);
  double *al = ALLOC(double, scene->num_faces + 1);
  double *ar = ALLOC(double, scene->num_faces + 1);
  int num_nodes = 0;
  scene->bvh_nodes = bvh_make_recur(scene->num_faces, scene->faces, al, ar, scene->bvh_nodes, &num_nodes);
  free(al);
  free(ar);
  printf("bvh usage: %i/%i\n", num_nodes, conf.num_bvh_nodes);
}
FACE *isect_faces(int num_faces, FACE *faces, VEC o, VEC d, double *t) {
  FACE *face = NULL;
  *t = INFINITY;
  for (int i = 0; i < num_faces; ++i) {
    FACE *facec = faces + i;
    double tc;
    if (tri_isect(*facec->v0, *facec->v1, *facec->v2, o, d, &tc) && tc < *t) {
      face = facec;
      *t = tc;
    }
  }
  return face;
}
bool visible_faces(int num_faces, FACE *faces, VEC o, VEC d, double t) {
  double eps = EPS;
  for (int i = 0; i < num_faces; ++i) {
    FACE *facec = faces + i;
    double tc;
    if (tri_isect(*facec->v0, *facec->v1, *facec->v2, o, d, &tc) && tc + eps < t) {
      return false;
    }
  }
  return true;
}
FACE *bvh_isect(BVH_NODE *node, VEC o, VEC d, double *t) {
  //if (node == NULL) {
  //  return NULL;
  //}
  if (!box_isect(node->box, o, d)) {
    return NULL;
  }
  if (node->faces == NULL) {
    FACE *face = NULL;
    *t = INFINITY;
    double tl, tr;
    FACE *fl = bvh_isect(node->l, o, d, &tl);
    FACE *fr = bvh_isect(node->r, o, d, &tr);
    if (fl != NULL && tl < *t) {
      face = fl;
      *t = tl;
    }
    if (fr != NULL && tr < *t) {
      face = fr;
      *t = tr;
    }
    return face;
  }
  return isect_faces(node->num_faces, node->faces, o, d, t);
}
bool bvh_visible(BVH_NODE *node, VEC o, VEC d, double t) {
  //if (node == NULL) {
  //  return NULL;
  //}
  if (!box_isect(node->box, o, d)) {
    return true;
  }
  if (node->faces == NULL) {
    if (!bvh_visible(node->l, o, d, t)) {
      return false;
    }
    if (!bvh_visible(node->r, o, d, t)) {
      return false;
    }
    return true;
  }
  return visible_faces(node->num_faces, node->faces, o, d, t);
}

SCENE scene_read(int argc, char **argv) {
  int num_ioscenes = argc - 3;
  IO_SCENE *ioscenes = ALLOC(IO_SCENE, num_ioscenes);
  for (int i = 0; i < num_ioscenes; ++i) {
    ioscenes[i] = io_read_scene(argv[i + 3]);
  }

  SCENE scene;

  scene.num_edfs = 0;
  scene.num_bsdfs = 0;
  scene.num_media = 0;
  scene.num_verts = 0;
  scene.num_faces = 0;
  for (int i = 0; i < num_ioscenes; ++i) {
    scene.num_edfs += ioscenes[i].num_edfs;
    scene.num_bsdfs += ioscenes[i].num_bsdfs;
    scene.num_media += ioscenes[i].num_media;
    scene.num_verts += ioscenes[i].num_verts;
    scene.num_faces += ioscenes[i].num_faces;
  }

  scene.edfs = ALLOC(EDF, scene.num_edfs);
  scene.bsdfs = ALLOC(BSDF, scene.num_bsdfs);
  scene.media = ALLOC(MEDIUM, scene.num_media);
  scene.verts = ALLOC(VEC, scene.num_verts);
  scene.faces = ALLOC(FACE, scene.num_faces);
  int iedf = 0, ibsdf = 0, imedium = 0, ivert = 0, iface = 0;
  for (int i = 0; i < num_ioscenes; ++i) {
    IO_SCENE ioscene = ioscenes[i];
    VEC t = {ioscene.t.x, ioscene.t.y, ioscene.t.z};
    MAT mat = {
      ioscene.m.xx, ioscene.m.xy, ioscene.m.xz,
      ioscene.m.yx, ioscene.m.yy, ioscene.m.yz,
      ioscene.m.zx, ioscene.m.zy, ioscene.m.zz};
    for (int j = 0; j < ioscene.num_edfs; ++j) {
      IO_EDF s = ioscene.edfs[j];
      EDF d = {s.t, {s.a.x, s.a.y, s.a.z}};
      scene.edfs[iedf + j] = d;
    }
    for (int j = 0; j < ioscene.num_bsdfs; ++j) {
      IO_BSDF s = ioscene.bsdfs[j];
      BSDF d = {s.t, {s.a.x, s.a.y, s.a.z}};
      scene.bsdfs[ibsdf + j] = d;
    }
    for (int j = 0; j < ioscene.num_media; ++j) {
      IO_MEDIUM s = ioscene.media[j];
      MEDIUM d = {s.eta, {s.a.x, s.a.y, s.a.z}};
      scene.media[imedium + j] = d;
    }
    for (int j = 0; j < ioscene.num_verts; ++j) {
      IO_VEC s = ioscene.verts[j];
      VEC d = {s.x, s.y, s.z};
      d = vec_add(t, mat_trafo(mat, d));
      scene.verts[ivert + j] = d;
    }
    for (int j = 0; j < ioscene.num_faces; ++j) {
      IO_FACE s = ioscene.faces[j];
      FACE d = {
        scene.edfs + iedf + s.edf,
        scene.bsdfs + ibsdf + s.bsdf,
        scene.media + imedium + s.mp,
        scene.media + imedium + s.mn,
        scene.verts + ivert + s.v0,
        scene.verts + ivert + s.v1,
        scene.verts + ivert + s.v2,
        0.0
      };
      scene.faces[iface + j] = d;
    }
    iedf += ioscene.num_edfs;
    ibsdf += ioscene.num_bsdfs;
    imedium += ioscene.num_media;
    ivert += ioscene.num_verts;
    iface += ioscene.num_faces;
  }

  for (int i = 0; i < num_ioscenes; ++i) {
    io_free_scene(ioscenes[i]);
  }
  free(ioscenes);

  return scene;
}
CONF conf_read(char *path) {
  FILE *f = fopen(path, "r");
  if (f == NULL) {
    printf("cannot open file: %s\n", path);
    exit(0);
  }

  int buflen = 256, idlen = 256;
  char buf[buflen], id[idlen];

  CONF conf;

  while (true) {
    if (fgets(buf, buflen, f) == NULL) {
      break;
    }
    sscanf(buf, "%s ", id);
    if (strcmp(id, "num-threads") == 0) {
      sscanf(buf, "%*s %i", &conf.num_threads);
    } else if (strcmp(id, "num-samples") == 0) {
      sscanf(buf, "%*s %i", &conf.num_samples);
    } else if (strcmp(id, "num-verticies") == 0) {
      sscanf(buf, "%*s %i", &conf.num_verts);
    } else if (strcmp(id, "num-iterations") == 0) {
      sscanf(buf, "%*s %i", &conf.num_iters);
    } else if (strcmp(id, "num-bvh-nodes") == 0) {
      sscanf(buf, "%*s %i", &conf.num_bvh_nodes);
    } else if (strcmp(id, "num-photons") == 0) {
      sscanf(buf, "%*s %i", &conf.num_phots);
    } else if (strcmp(id, "num-kdt-nodes") == 0) {
      sscanf(buf, "%*s %i", &conf.num_kdt_nodes);
    } else if (strcmp(id, "image-width") == 0) {
      sscanf(buf, "%*s %d", &conf.wi);
    } else if (strcmp(id, "image-height") == 0) {
      sscanf(buf, "%*s %i", &conf.hi);
    } else if (strcmp(id, "film-exposure") == 0) {
      sscanf(buf, "%*s %lf", &conf.ef);
    } else if (strcmp(id, "film-width") == 0) {
      sscanf(buf, "%*s %lf", &conf.wf);
    } else if (strcmp(id, "film-height") == 0) {
      sscanf(buf, "%*s %lf", &conf.hf);
    } else if (strcmp(id, "film-depth") == 0) {
      sscanf(buf, "%*s %lf", &conf.df);
    } else if (strcmp(id, "lens-radius") == 0) {
      sscanf(buf, "%*s %lf", &conf.r);
    } else if (strcmp(id, "lens-focal-length") == 0) {
      sscanf(buf, "%*s %lf", &conf.f);
    } else if (strcmp(id, "initial-radius") == 0) {
      sscanf(buf, "%*s %lf", &conf.ir);
    }
  }

  fclose(f);
  return conf;
}
void sample_light(SCENE scene, PRNG_STATE *s,
FACE **f, VEC *n, VEC *o, VEC *d, VEC *fs, double *pa, double *ps, bool *c) {
  *f = scene.lights[discrete_sample(scene.num_lights, scene.cdf, s)];
  FACE face = **f;
  *n = tri_norm(*face.v0, *face.v1, *face.v2);
  *o = tri_sample(*face.v0, *face.v1, *face.v2, s);
  *d = edf_sample(*face.edf, *n, s);
  *fs = edf_eval(*face.edf, *n, *d);
  *pa = face.p * tri_pdf(*face.v0, *face.v1, *face.v2);
  *ps = edf_pdf(*face.edf, *n, *d);
  *c = edf_connect(*face.edf);
}
void sample_eye(CONF conf, PRNG_STATE *s, double x, double y,
FACE **f, VEC *n, VEC *o, VEC *d, VEC *fs, double *pa, double *ps, bool *c) {
  *f = NULL;
  *n = vec_init(0.0, 0.0, -1.0);
  *o = disc_sample(conf.r, s);
  *d = vec_init(conf.wf * (0.5 - x / conf.wi), conf.hf * (y / conf.hi - 0.5), conf.df);
  *d = vec_norm(vec_sub(vec_scale(conf.f / (conf.f - conf.df), *d), *o));
  *fs = vec_inita(
  (conf.ef * conf.wi * conf.hi) *
  (1.0 / (PI * sqr(conf.r))) *
  (sqr(conf.df) / (quad(fabs(d->z)) * conf.wf * conf.hf))
  );
  *pa = disc_pdf(conf.r);
  *ps = sqr(conf.df) / (cube(fabs(d->z)) * conf.wf * conf.hf);
  *c = true;
}
void debug(CONF conf, SCENE scene, IMAGE img, PRNG_STATE *s) {
  for (int i = omp_get_thread_num(); i < img.w * img.h; i += omp_get_max_threads()) {
    int ix = i % img.w, iy = i / img.w;
    VEC sum0 = vec_inita(0.0);
    for (int i = 0; i < conf.num_samples; ++i) {
      FACE *f;
      VEC n, o, d, fs;
      double pa, ps;
      bool c;
      sample_eye(conf, s, ix, iy, &f, &n, &o, &d, &fs, &pa, &ps, &c);
      VEC a = vec_scale(fabs(vec_dot(n, d)) / (pa * ps), fs);
      double t;
      FACE *face = bvh_isect(scene.bvh_nodes, o, d, &t);
      if (face != NULL) {
        VEC normal = tri_norm(*face->v0, *face->v1, *face->v2);
        normal.x = fabs(normal.x);
        normal.y = fabs(normal.y);
        normal.z = fabs(normal.z);
        sum0 = vec_add(sum0, vec_mul(a, vec_scale(0.1, normal)));
      }
    }
    sum0 = vec_scale(1.0 / conf.num_samples, sum0);
    sum0 = vec_scale(1.0 / (img.w * img.h), sum0);
    img.data[3 * i + 0] = sum0.x;
    img.data[3 * i + 1] = sum0.y;
    img.data[3 * i + 2] = sum0.z;
  }
}
IMAGE render_debug(CONF conf, SCENE scene) {
  IMAGE img = image_init(3, conf.wi, conf.hi);
  omp_set_num_threads(conf.num_threads);
  #pragma omp parallel
  {
    PRNG_STATE s = prng_seed(SEED);
    for (int i = 0; i < omp_get_thread_num(); ++i) {
      s = prng_jump(s);
    }
    debug(conf, scene, img, &s);
  }
  return img;
}
void trace(CONF conf, SCENE scene, PRNG_STATE *s, PATH *path, FACE *f, VEC n, VEC o, VEC d, VEC fs,
double pa, double ps, bool c, bool light) {
  VEC a = vec_inita(1.0 / pa);
  VERTEX v;
  v.f = f;
  v.c = c;
  v.pf = pa;
  v.a = a;
  v.x = o;
  v.n = n;
  path->vs[0] = v;
  a = vec_scale(fabs(vec_dot(n, d)) / ps, vec_mul(a, fs));
  double pr = 1.0;
  for (int i = 1; i < conf.num_verts; ++i) {
    if (prng_db(s) >= pr) {
      path->n = i;
      return;
    }
    double t;
    f = bvh_isect(scene.bvh_nodes, o, d, &t);
    if (f == NULL) {
      path->n = i;
      return;
    }

    o = vec_add(o, vec_scale(t, d));
    n = tri_norm(*f->v0, *f->v1, *f->v2);

    double cosi = vec_dot(n, d);
    MEDIUM mir, mit;
    media_sides(cosi, *f->mp, *f->mn, &mir, &mit);

    VEC wo = bsdf_sample(*f->bsdf, mir.eta, mit.eta, n, d, s);

    double coso = vec_dot(n, wo);
    MEDIUM mor, mot;
    media_sides(coso, *f->mp, *f->mn, &mor, &mot);

    VEC tr = transmittance(mir.a, t);

    a = vec_scale(1.0 / pr, vec_mul(tr, a));

    if (i >= 2) {
      VERTEX vm1 = path->vs[i - 1], *vm2 = path->vs + i - 2;
      double pbrm2 = luminance(vec_mul(tr, vm1.f->bsdf->a));
      double tm1 = vec_len(vec_sub(vm1.x, vm2->x));
      vm2->pb = pbrm2 * ps * fabs(vec_dot(vm2->n, vm1.wi)) / sqr(tm1);
    }

    v.f = f;
    v.c = bsdf_connect(*f->bsdf);
    v.pf = pr * ps * fabs(cosi) / sqr(t);
    v.a = a;
    v.x = o;
    v.n = n;
    v.wi = d;
    path->vs[i] = v;

    fs = light ?
    bsdf_eval(*f->bsdf, mir.eta, mit.eta, n, vec_scale(1.0, d), vec_scale(1.0, wo)) :
    bsdf_eval(*f->bsdf, mot.eta, mor.eta, n, vec_scale(-1.0, wo), vec_scale(-1.0, d));
    ps = bsdf_pdf(*f->bsdf, mir.eta, mit.eta, n, d, wo);
    a = vec_scale(fabs(coso) / ps, vec_mul(a, fs));
    pr = luminance(vec_mul(tr, f->bsdf->a));
    d = wo;
  }
  path->n = conf.num_verts;
}
void pt(CONF conf, SCENE scene, IMAGE img, PRNG_STATE *s) {
  PATH path;
  path.vs = ALLOC(VERTEX, conf.num_verts);
  for (int i = omp_get_thread_num(); i < img.w * img.h; i += omp_get_max_threads()) {
    int ix = i % img.w, iy = i / img.w;
    VEC sum0 = vec_inita(0.0);
    for (int i = 0; i < conf.num_samples; ++i) {
      FACE *f;
      VEC n, o, d, fs;
      double pa, ps;
      bool c;
      sample_eye(conf, s, ix, iy, &f, &n, &o, &d, &fs, &pa, &ps, &c);
      trace(conf, scene, s, &path, f, n, o, d, fs, pa, ps, c, false);
      VEC sum1 = vec_inita(0.0);
      for (int i = 1; i < path.n; ++i) {
        VERTEX v = path.vs[i];
        FACE *f = v.f;
        if (edf_connect(*f->edf)) {
          sum1 = vec_add(sum1, vec_mul(v.a, edf_eval(*f->edf, v.n, vec_scale(-1.0, v.wi))));
        }
      }
      sum0 = vec_add(sum0, sum1);
    }
    sum0 = vec_scale(1.0 / conf.num_samples, sum0);
    sum0 = vec_scale(1.0 / (img.w * img.h), sum0);
    img.data[3 * i + 0] = sum0.x;
    img.data[3 * i + 1] = sum0.y;
    img.data[3 * i + 2] = sum0.z;
  }
}
IMAGE render_pt(CONF conf, SCENE scene) {
  IMAGE img = image_init(3, conf.wi, conf.hi);
  omp_set_num_threads(conf.num_threads);
  #pragma omp parallel
  {
    PRNG_STATE s = prng_seed(SEED);
    for (int i = 0; i < omp_get_thread_num(); ++i) {
      s = prng_jump(s);
    }
    pt(conf, scene, img, &s);
  }
  return img;
}
int kdt_cmp_x(void const *A, void const *B) {
  KDT_NODE const *a = A, *b = B;
  double d = a->p->vs[a->j].x.x - b->p->vs[b->j].x.x;
  return d < 0.0 ? -1 : d > 0.0 ? 1 : 0;
}
int kdt_cmp_y(void const *A, void const *B) {
  KDT_NODE const *a = A, *b = B;
  double d = a->p->vs[a->j].x.y - b->p->vs[b->j].x.y;
  return d < 0.0 ? -1 : d > 0.0 ? 1 : 0;
}
int kdt_cmp_z(void const *A, void const *B) {
  KDT_NODE const *a = A, *b = B;
  double d = a->p->vs[a->j].x.z - b->p->vs[b->j].x.z;
  return d < 0.0 ? -1 : d > 0.0 ? 1 : 0;
}
void kdt_make(int a, int b, KDT_NODE *nodes, int depth) {
  if (b - a <= 1) {
    return;
  }
  int axis = depth % 3;
  switch (axis) {
    case 0 : {
      qsort(nodes + a, b - a, sizeof(KDT_NODE), kdt_cmp_x);
    } break;
    case 1 : {
      qsort(nodes + a, b - a, sizeof(KDT_NODE), kdt_cmp_y);
    } break;
    case 2 : {
      qsort(nodes + a, b - a, sizeof(KDT_NODE), kdt_cmp_z);
    } break;
    default : {
    } break;
  }
  int im = (a + b) / 2;
  kdt_make(a, im, nodes, depth + 1);
  kdt_make(im + 1, b, nodes, depth + 1);
}
void kdt_acc(int a, int b, KDT_NODE *nodes, int depth, int N, double r, PATH pe, VERTEX ve, int je, VEC *sum) {
  if (b - a <= 0) {
    return;
  }
  int im = (a + b) / 2;
  KDT_NODE node = nodes[im];
  int jl = node.j;
  PATH pl = *node.p;
  VERTEX vl = pl.vs[jl];
  double d = vec_len(vec_sub(vl.x, ve.x));
  if (d < r) {
    VERTEX vlm1 = pl.vs[jl - 1];
    VERTEX vem1 = pe.vs[je - 1];

    double tl = vec_len(vec_sub(vl.x, vlm1.x));
    double te = vec_len(vec_sub(ve.x, vem1.x));

    double gl = fabs(vec_dot(vlm1.n, vl.wi)) / sqr(tl);
    double ge = fabs(vec_dot(vem1.n, ve.wi)) / sqr(te);

    MEDIUM mlr, mlt, mer, met;
    media_sides(vec_dot(vl.n, vl.wi), *ve.f->mp, *ve.f->mn, &mlr, &mlt);
    media_sides(vec_dot(ve.n, ve.wi), *ve.f->mp, *ve.f->mn, &mer, &met);

    double psl = bsdf_pdf(*ve.f->bsdf, mlr.eta, mlt.eta, vl.n, vl.wi, vec_scale(-1.0, ve.wi));
    double pse = bsdf_pdf(*ve.f->bsdf, mer.eta, met.eta, ve.n, ve.wi, vec_scale(-1.0, vl.wi));

    double pblr = luminance(vec_mul(transmittance(mer.a, te), ve.f->bsdf->a));
    double pbla = gl * pse;
    double pbl = pblr * pbla;

    double pber = luminance(vec_mul(transmittance(mlr.a, tl), ve.f->bsdf->a));
    double pbea = ge * psl;
    double pbe = pber * pbea;

    double w = 1.0;

    if (jl >= 2) {
      double p = 1.0;
      p *= pbl / vl.pf;
      if (vlm1.c) {
        w += sqr(p);
      }
      for (int kl = jl - 2; kl >= 1; --kl) {
        VERTEX vlp1 = pl.vs[kl + 1], vl = pl.vs[kl];
        p *= vl.pb / vlp1.pf;
        if (vl.c) {
          w += sqr(p);
        }
      }
    }

    if (jl == 1) {
      double p = 1.0 / (N * disc_area(r));
      p *= 1.0 / vl.pf;
      if (vlm1.c) {
        w += sqr(p);
      }
      p *= pbl / vlm1.pf;
      if (vlm1.c) {
        w += sqr(p);
      }
    } else if (jl == 2) {
      VERTEX vlm2 = pl.vs[0];
      double p = 1.0 / (N * disc_area(r));
      p *= 1.0 / vl.pf;
      if (vlm1.c) {
        w += sqr(p);
      }
      p *= pbl / vlm1.pf;
      if (vlm1.c && vlm2.c) {
        w += sqr(p);
      }
      p *= vlm2.pb / vlm2.pf;
      if (vlm2.c) {
        w += sqr(p);
      }
    } else {
      VERTEX vlm2 = pl.vs[jl - 2];
      double p = 1.0 / (N * disc_area(r));
      p *= 1.0 / vl.pf;
      if (vlm1.c) {
        w += sqr(p);
      }
      p *= pbl / vlm1.pf;
      if (vlm1.c && vlm2.c) {
        w += sqr(p);
      }
      for (int kl = jl - 2; kl >= 2; --kl) {
        VERTEX vl = pl.vs[kl], vlm1 = pl.vs[kl - 1];
        p *= vl.pb / vl.pf;
        if (vl.c && vlm1.c) {
          w += sqr(p);
        }
      }
      VERTEX vl1 = pl.vs[1], vl0 = pl.vs[0];
      p *= vl1.pb / vl1.pf;
      if (vl1.c && vl0.c) {
        w += sqr(p);
      }
      p *= vl0.pb / vl0.pf;
      if (vl0.c) {
        w += sqr(p);
      }
    }

    if (je >= 2) {
      double p = 1.0;
      p *= pbe / ve.pf;
      if (vem1.c) {
        w += sqr(p);
      }
      for (int ke = je - 2; ke >= 1; --ke) {
        VERTEX vep1 = pe.vs[ke + 1], ve = pe.vs[ke];
        p *= ve.pb / vep1.pf;
        if (ve.c) {
          w += sqr(p);
        }
      }
    }

    if (je == 1) {
    } else if (je == 2) {
      double p = 1.0 / (N * disc_area(r));
      p *= 1.0 / ve.pf;
      if (vem1.c) {
        w += sqr(p);
      }
    } else {
      VERTEX vem2 = pe.vs[je - 2];
      double p = 1.0 / (N * disc_area(r));
      p *= 1.0 / ve.pf;
      if (vem1.c) {
        w += sqr(p);
      }
      p *= pbe / vem1.pf;
      if (vem1.c && vem2.c) {
        w += sqr(p);
      }
      for (int ke = je - 2; ke >= 2; --ke) {
        VERTEX ve = pe.vs[ke], vem1 = pe.vs[ke - 1];
        p *= ve.pb / ve.pf;
        if (ve.c && vem1.c) {
          w += sqr(p);
        }
      }
    }

    w = 1.0 / w;

    double k = opst(vec_dot(ve.n, vec_norm(vec_sub(vl.x, ve.x)))) / disc_area(r);
    VEC fs = bsdf_eval(*ve.f->bsdf, mlr.eta, mlt.eta, ve.n, vl.wi, vec_scale(-1.0, ve.wi));
    VEC I = vec_scale(w * k, vec_mul(vec_mul(ve.a, vl.a), fs));
    *sum = vec_add(*sum, I);
  }
  double t;
  int axis = depth % 3;
  switch (axis) {
    case 0 : {
      t = ve.x.x - vl.x.x;
    } break;
    case 1 : {
      t = ve.x.y - vl.x.y;
    } break;
    case 2 : {
      t = ve.x.z - vl.x.z;
    } break;
    default : {
      t = 0.0;
    } break;
  }
  if (t < r) {
    kdt_acc(a, im, nodes, depth + 1, N, r, pe, ve, je, sum);
  }
  if (t > -r) {
    kdt_acc(im + 1, b, nodes, depth + 1, N, r, pe, ve, je, sum);
  }
}
VEC vm(int N, double r, int num_kdt_nodes, KDT_NODE *nodes, PATH pe) {
  VEC sum0 = vec_inita(0.0);
  for (int je = 1; je < pe.n; ++je) {
    VERTEX ve = pe.vs[je];
    if (ve.c) {
      VEC sum1 = vec_inita(0.0);
      kdt_acc(0, num_kdt_nodes, nodes, 0, N, r, pe, ve, je, &sum1);
      sum0 = vec_add(sum0, sum1);
    }
  }
  sum0 = vec_scale(1.0 / N, sum0);
  return sum0;
}
VEC vc(SCENE scene, int N, double r, PATH pl, PATH pe) {
  VEC sum0 = vec_inita(0.0);
  //s = 0, t >= 2
  for (int je = 1; je < pe.n; ++je) {
    VERTEX ve = pe.vs[je];
    bool ce = edf_connect(*ve.f->edf);
    if (ce) {
      VERTEX vem1 = pe.vs[je - 1];

      double te = vec_len(vec_sub(ve.x, vem1.x));

      double ge1 = fabs(vec_dot(vem1.n, ve.wi)) / sqr(te);

      double pse = edf_pdf(*ve.f->edf, ve.n, vec_scale(-1.0, ve.wi));

      double pbea0 = ve.f->p * tri_pdf(*ve.f->v0, *ve.f->v1, *ve.f->v2);
      double pbe0 = pbea0;

      double pbea1 = ge1 * pse;
      double pbe1 = pbea1;

      double w = 1.0;

      if (je == 1) {
      } else if (je == 2) {
        double p = 1.0;
        p *= pbe0 / ve.pf;
        if (vem1.c) {
          w += sqr(p);
        }
      } else {
        VERTEX vem2 = pe.vs[je - 2];
        double p = 1.0;
        p *= pbe0 / ve.pf;
        if (vem1.c) {
          w += sqr(p);
        }
        p *= pbe1 / vem1.pf;
        if (vem1.c && vem2.c) {
          w += sqr(p);
        }
        for (int ke = je - 2; ke >= 2; --ke) {
          VERTEX ve = pe.vs[ke], vem1 = pe.vs[ke - 1];
          p *= ve.pb / ve.pf;
          if (ve.c && vem1.c) {
            w += sqr(p);
          }
        }
      }

      if (je == 1) {
      } else {
        double p = N * disc_area(r);
        p *= pbe0;
        p *= pbe1 / ve.pf;
        if (vem1.c) {
          w += sqr(p);
        }
        for (int ke = je - 2; ke >= 1; --ke) {
          VERTEX vep1 = pe.vs[ke + 1], ve = pe.vs[ke];
          p *= ve.pb / vep1.pf;
          if (ve.c) {
            w += sqr(p);
          }
        }
      }

      w = 1.0 / w;

      VEC fsl = edf_eval(*ve.f->edf, ve.n, vec_scale(-1.0, ve.wi));
      VEC I = vec_scale(w, vec_mul(ve.a, fsl));
      sum0 = vec_add(sum0, I);
    }
  }
  //s = 1, t >= 2
  {
    VERTEX vl = pl.vs[0];
    bool cl = vl.c;
    if (cl) {
      for (int je = 1; je < pe.n; ++je) {
        VERTEX ve = pe.vs[je];
        bool ce = ve.c;
        if (ce) {
          VEC d = vec_sub(vl.x, ve.x);
          double t = vec_len(d);
          VEC dl = vec_scale(-1.0 / t, d);
          VEC de = vec_scale(1.0 / t, d);
          if (bvh_visible(scene.bvh_nodes, ve.x, de, t)) {
            VERTEX vem1 = pe.vs[je - 1];

            double te = vec_len(vec_sub(ve.x, vem1.x));

            double g = fabs(vec_dot(vl.n, dl) * vec_dot(ve.n, de)) / sqr(t);
            double gl0 = fabs(vec_dot(vl.n, dl)) / sqr(t);
            double ge0 = fabs(vec_dot(ve.n, de)) / sqr(t);
            double ge1 = fabs(vec_dot(vem1.n, ve.wi)) / sqr(te);

            MEDIUM meir, meit, meor, meot;
            media_sides(vec_dot(ve.n, ve.wi), *ve.f->mp, *ve.f->mn, &meir, &meit);
            media_sides(vec_dot(ve.n, de), *ve.f->mp, *ve.f->mn, &meor, &meot);

            VEC tr = transmittance(meot.a, t);

            double psl = edf_pdf(*vl.f->edf, vl.n, dl);
            double pse = bsdf_pdf(*ve.f->bsdf, meir.eta, meit.eta, ve.n, ve.wi, de);

            double pblr0 = luminance(vec_mul(transmittance(meir.a, te), ve.f->bsdf->a));
            double pbla0 = gl0 * pse;
            double pbl0 = pblr0 * pbla0;

            double pbea0 = ge0 * psl;
            double pbe0 = pbea0;

            double pber1 = luminance(vec_mul(tr, ve.f->bsdf->a));
            double pbea1 = ge1 * pse;
            double pbe1 = pber1 * pbea1;

            double w = 1.0;

            {
              double p = 1.0;
              p *= pbl0 / vl.pf;
              w += sqr(p);
            }

            if (je == 1) {
            } else if (je == 2) {
              double p = 1.0;
              p *= pbe0 / ve.pf;
              if (vem1.c) {
                w += sqr(p);
              }
            } else {
              VERTEX vem2 = pe.vs[je - 2];
              double p = 1.0;
              p *= pbe0 / ve.pf;
              if (vem1.c) {
                w += sqr(p);
              }
              p *= pbe1 / vem1.pf;
              if (vem1.c && vem2.c) {
                w += sqr(p);
              }
              for (int ke = je - 2; ke >= 2; --ke) {
                VERTEX ve = pe.vs[ke], vem1 = pe.vs[ke - 1];
                p *= ve.pb / ve.pf;
                if (ve.c && vem1.c) {
                  w += sqr(p);
                }
              }
            }

            if (je == 1) {
              double p = N * disc_area(r);
              p *= pbe0;
              w += sqr(p);
            } else {
              double p = N * disc_area(r);
              p *= pbe0;
              w += sqr(p);
              p *= pbe1 / ve.pf;
              if (vem1.c) {
                w += sqr(p);
              }
              for (int ke = je - 2; ke >= 1; --ke) {
                VERTEX vep1 = pe.vs[ke + 1], ve = pe.vs[ke];
                p *= ve.pb / vep1.pf;
                if (ve.c) {
                  w += sqr(p);
                }
              }
            }

            w = 1.0 / w;

            VEC fsl = edf_eval(*vl.f->edf, vl.n, dl);
            VEC fse = bsdf_eval(*ve.f->bsdf, meot.eta, meor.eta, ve.n, dl, vec_scale(-1.0, ve.wi));
            VEC I = vec_scale(w * g, vec_mul(vec_mul(vec_mul(vec_mul(vl.a, ve.a), fsl), fse), tr));
            sum0 = vec_add(sum0, I);
          }
        }
      }
    }
  }
  //s >= 2, t >= 2
  for (int jl = 1; jl < pl.n; ++jl) {
    VERTEX vl = pl.vs[jl];
    bool cl = vl.c;
    if (cl) {
      for (int je = 1; je < pe.n; ++je) {
        VERTEX ve = pe.vs[je];
        bool ce = ve.c;
        if (ce) {
          VEC d = vec_sub(vl.x, ve.x);
          double t = vec_len(d);
          VEC dl = vec_scale(-1.0 / t, d);
          VEC de = vec_scale(1.0 / t, d);
          if (bvh_visible(scene.bvh_nodes, ve.x, de, t)) {
            VERTEX vlm1 = pl.vs[jl - 1];
            VERTEX vem1 = pe.vs[je - 1];

            double tl = vec_len(vec_sub(vl.x, vlm1.x));
            double te = vec_len(vec_sub(ve.x, vem1.x));

            double g = fabs(vec_dot(vl.n, dl) * vec_dot(ve.n, de)) / sqr(t);
            double gl0 = fabs(vec_dot(vl.n, dl)) / sqr(t);
            double gl1 = fabs(vec_dot(vlm1.n, vl.wi)) / sqr(tl);
            double ge0 = fabs(vec_dot(ve.n, de)) / sqr(t);
            double ge1 = fabs(vec_dot(vem1.n, ve.wi)) / sqr(te);

            MEDIUM mlir, mlit, mlor, mlot, meir, meit, meor, meot;
            media_sides(vec_dot(vl.n, vl.wi), *vl.f->mp, *vl.f->mn, &mlir, &mlit);
            media_sides(vec_dot(vl.n, dl), *vl.f->mp, *vl.f->mn, &mlor, &mlot);
            media_sides(vec_dot(ve.n, ve.wi), *ve.f->mp, *ve.f->mn, &meir, &meit);
            media_sides(vec_dot(ve.n, de), *ve.f->mp, *ve.f->mn, &meor, &meot);

            VEC tr = transmittance(meot.a, t);

            double psl = bsdf_pdf(*vl.f->bsdf, mlir.eta, mlit.eta, vl.n, vl.wi, dl);
            double pse = bsdf_pdf(*ve.f->bsdf, meir.eta, meit.eta, ve.n, ve.wi, de);

            double pblr0 = luminance(vec_mul(transmittance(meir.a, te), ve.f->bsdf->a));
            double pbla0 = gl0 * pse;
            double pbl0 = pblr0 * pbla0;

            double pblr1 = luminance(vec_mul(tr, vl.f->bsdf->a));
            double pbla1 = gl1 * psl;
            double pbl1 = pblr1 * pbla1;

            double pber0 = luminance(vec_mul(transmittance(mlir.a, tl), vl.f->bsdf->a));
            double pbea0 = ge0 * psl;
            double pbe0 = pber0 * pbea0;

            double pber1 = luminance(vec_mul(tr, ve.f->bsdf->a));
            double pbea1 = ge1 * pse;
            double pbe1 = pber1 * pbea1;

            double w = 1.0;

            if (jl == 1) {
              double p = 1.0;
              p *= pbl0 / vl.pf;
              if (vlm1.c) {
                w += sqr(p);
              }
              p *= pbl1 / vlm1.pf;
              if (vlm1.c) {
                w += sqr(p);
              }
            } else if (jl == 2) {
              VERTEX vlm2 = pl.vs[jl - 2];
              double p = 1.0;
              p *= pbl0 / vl.pf;
              if (vlm1.c) {
                w += sqr(p);
              }
              p *= pbl1 / vlm1.pf;
              if (vlm1.c && vlm2.c) {
                w += sqr(p);
              }
              p *= vlm2.pb / vlm2.pf;
              if (vlm2.c) {
                w += sqr(p);
              }
            } else {
              VERTEX vlm2 = pl.vs[jl - 2];
              double p = 1.0;
              p *= pbl0 / vl.pf;
              if (vlm1.c) {
                w += sqr(p);
              }
              p *= pbl1 / vlm1.pf;
              if (vlm1.c && vlm2.c) {
                w += sqr(p);
              }
              for (int kl = jl - 2; kl >= 2; --kl) {
                VERTEX vl = pl.vs[kl], vlm1 = pl.vs[kl - 1];
                p *= vl.pb / vl.pf;
                if (vl.c && vlm1.c) {
                  w += sqr(p);
                }
              }
              VERTEX vl1 = pl.vs[1], vl0 = pl.vs[0];
              p *= vl1.pb / vl1.pf;
              if (vl1.c && vl0.c) {
                w += sqr(p);
              }
              p *= vl0.pb / vl0.pf;
              if (vl0.c) {
                w += sqr(p);
              }
            }

            if (jl == 1) {
              double p = N * disc_area(r);
              p *= pbl0;
              w += sqr(p);
            } else {
              double p = N * disc_area(r);
              p *= pbl0;
              w += sqr(p);
              p *= pbl1 / vl.pf;
              if (vlm1.c) {
                w += sqr(p);
              }
              for (int kl = jl - 2; kl >= 1; --kl) {
                VERTEX vlp1 = pl.vs[kl + 1], vl = pl.vs[kl];
                p *= vl.pb / vlp1.pf;
                if (vl.c) {
                  w += sqr(p);
                }
              }
            }

            if (je == 1) {
            } else if (je == 2) {
              double p = 1.0;
              p *= pbe0 / ve.pf;
              if (vem1.c) {
                w += sqr(p);
              }
            } else {
              VERTEX vem2 = pe.vs[je - 2];
              double p = 1.0;
              p *= pbe0 / ve.pf;
              if (vem1.c) {
                w += sqr(p);
              }
              p *= pbe1 / vem1.pf;
              if (vem1.c && vem2.c) {
                w += sqr(p);
              }
              for (int ke = je - 2; ke >= 2; --ke) {
                VERTEX ve = pe.vs[ke], vem1 = pe.vs[ke - 1];
                p *= ve.pb / ve.pf;
                if (ve.c && vem1.c) {
                  w += sqr(p);
                }
              }
            }

            if (je == 1) {
              double p = N * disc_area(r);
              p *= pbe0;
              w += sqr(p);
            } else {
              double p = N * disc_area(r);
              p *= pbe0;
              w += sqr(p);
              p *= pbe1 / ve.pf;
              if (vem1.c) {
                w += sqr(p);
              }
              for (int ke = je - 2; ke >= 1; --ke) {
                VERTEX vep1 = pe.vs[ke + 1], ve = pe.vs[ke];
                p *= ve.pb / vep1.pf;
                if (ve.c) {
                  w += sqr(p);
                }
              }
            }

            w = 1.0 / w;

            VEC fsl = bsdf_eval(*vl.f->bsdf, mlir.eta, mlit.eta, vl.n, vl.wi, dl);
            VEC fse = bsdf_eval(*ve.f->bsdf, meot.eta, meor.eta, ve.n, dl, vec_scale(-1.0, ve.wi));
            VEC I = vec_scale(w * g, vec_mul(vec_mul(vec_mul(vec_mul(vl.a, ve.a), fsl), fse), tr));
            sum0 = vec_add(sum0, I);
          }
        }
      }
    }
  }
  return sum0;
}
void vcm(CONF conf, SCENE scene, IMAGE img, PRNG_STATE *states) {
  image_fill(img, 0.0);
  int N = conf.num_samples * conf.num_samples * conf.wi * conf.hi;

  IMG_SPL *imgspls = ALLOC(IMG_SPL, N);
  PRNG_STATE s1 = prng_seed(SEED);

  VERTEX *vss[omp_get_max_threads()];
  int num_phots_per_thread = conf.num_phots / omp_get_max_threads();
  for (int i = 0; i < omp_get_max_threads(); ++i) {
    vss[i] = ALLOC(VERTEX, num_phots_per_thread);
  }
  PATH *pls = ALLOC(PATH, N);
  PATH pes[omp_get_max_threads()];
  for (int i = 0; i < omp_get_max_threads(); ++i) {
    pes[i].vs = ALLOC(VERTEX, conf.num_verts);
  }
  KDT_NODE *nodes = ALLOC(KDT_NODE, conf.num_kdt_nodes);
  for (int iter_num = 0; iter_num < conf.num_iters; ++iter_num) {
    double t;
    printf("number of iterations: %d\n", iter_num);
    t = time_get();
    #pragma omp parallel
    {
      PRNG_STATE *s = states + omp_get_thread_num();
      int offset = 0;
      for (int i = omp_get_thread_num(); i < N; i += omp_get_max_threads()) {
        PATH *pl = pls + i;
        pl->vs = vss[omp_get_thread_num()] + offset;
        FACE *f;
        VEC n, o, d, fs;
        double pa, ps;
        bool c;
        sample_light(scene, s, &f, &n, &o, &d, &fs, &pa, &ps, &c);
        trace(conf, scene, s, pl, f, n, o, d, fs, pa, ps, c, true);
        offset += pl->n;
      }
    }
    t = time_get() - t;
    printf("light paths trace time: %f\n", t);
    int num_phots = 0;
    for (int i = 0; i < N; ++i) {
      num_phots += pls[i].n;
    }
    printf("number of photons: %d/%d\n", num_phots, conf.num_phots);
    int num_kdt_nodes = 0;
    for (int i = 0; i < N; ++i) {
      PATH *p = pls + i;
      for (int j = 1; j < p->n; ++j) {
        if (bsdf_connect(*p->vs[j].f->bsdf)) {
          KDT_NODE node = {p, j};
          nodes[num_kdt_nodes] = node;
          ++num_kdt_nodes;
        }
      }
    }
    printf("number of nodes: %d/%d\n", num_kdt_nodes, conf.num_kdt_nodes);
    printf("average number of verticies: %f\n", num_phots / CAST(double, N));
    t = time_get();
    kdt_make(0, num_kdt_nodes, nodes, 0);
    t = time_get() - t;
    printf("kd tree build time: %f\n", t);
    t = time_get();

	for (int iy = 0; iy < conf.hi; ++iy) {
		for (int ix = 0; ix < conf.wi; ++ix) {
			for (int idy = 0; idy < conf.num_samples; ++idy) {
				for (int idx = 0; idx < conf.num_samples; ++idx) {
					int n = conf.num_samples;

					int i = ix + conf.wi * iy;

					double u = prng_db(&s1);
					double v = prng_db(&s1);

					double x = (u + idx) / n + ix;
					double y = (v + idy) / n + iy;

					VEC c = vec_inita(0.0);

					IMG_SPL spl = {i, x, y, c};

					imgspls[idx + n * idy + n * n * ix + n * n * conf.wi * iy] = spl;
				}
			}
		}
	}

	for (int i = N - 1; i > 0; --i) {
		int j = prng_db(&s1) * i;
		IMG_SPL tmp = imgspls[i];
		imgspls[i] = imgspls[j];
		imgspls[j] = tmp;
	}

    #pragma omp parallel
    {
      PRNG_STATE *s = states + omp_get_thread_num();
      for (int i = omp_get_thread_num(); i < N; i += omp_get_max_threads()) {
		IMG_SPL *spl = imgspls + i;

        PATH *pl = pls + i;
        PATH *pe = pes + omp_get_thread_num();

        FACE *f;
        VEC n, o, d, fs;
        double pa, ps;
        bool c;
		double r = conf.ir * pow(i + 1, -1.0 / 12.0);
        sample_eye(conf, s, spl->x, spl->y, &f, &n, &o, &d, &fs, &pa, &ps, &c);
        trace(conf, scene, s, pe, f, n, o, d, fs, pa, ps, c, false);
        spl->c = vec_scale(1.0 / N, vec_add(vm(N, r, num_kdt_nodes, nodes, *pe), vc(scene, N, r, *pl, *pe)));
      }
    }
    t = time_get() - t;
    printf("evaluation time: %f\n", t);

	for (int i = 0; i < N; ++i) {
		IMG_SPL spl = imgspls[i];
		img.data[3 * spl.i + 0] += spl.c.x;
		img.data[3 * spl.i + 1] += spl.c.y;
		img.data[3 * spl.i + 2] += spl.c.z;
	}
  }
  for (int i = 0; i < img.n * img.w * img.h; ++i) {
    img.data[i] /= conf.num_iters;
  }
}
IMAGE render_vcm(CONF conf, SCENE scene) {
  IMAGE img = image_init(3, conf.wi, conf.hi);
  omp_set_num_threads(conf.num_threads);
  PRNG_STATE states[omp_get_max_threads()];
  for (int i = 0; i < omp_get_max_threads(); ++i) {
    states[i] = prng_seed(SEED);
    for (int j = 0; j < i + 1; ++j) {
      states[i] = prng_jump(states[i]);
    }
  }
  vcm(conf, scene, img, states);
  return img;
}
VEC bpt_vc(SCENE scene, PATH pl, PATH pe) {
  VEC sum0 = vec_inita(0.0);
  //s = 0, t >= 2
  for (int je = 1; je < pe.n; ++je) {
    VERTEX ve = pe.vs[je];
    bool ce = edf_connect(*ve.f->edf);
    if (ce) {
      VERTEX vem1 = pe.vs[je - 1];

      double te = vec_len(vec_sub(ve.x, vem1.x));

      double ge1 = fabs(vec_dot(vem1.n, ve.wi)) / sqr(te);

      double pse = edf_pdf(*ve.f->edf, ve.n, vec_scale(-1.0, ve.wi));

      double pbea0 = ve.f->p * tri_pdf(*ve.f->v0, *ve.f->v1, *ve.f->v2);
      double pbe0 = pbea0;

      double pbea1 = ge1 * pse;
      double pbe1 = pbea1;

      double w = 1.0;

      if (je == 1) {
      } else if (je == 2) {
        double p = 1.0;
        p *= pbe0 / ve.pf;
        if (vem1.c) {
          w += sqr(p);
        }
      } else {
        VERTEX vem2 = pe.vs[je - 2];
        double p = 1.0;
        p *= pbe0 / ve.pf;
        if (vem1.c) {
          w += sqr(p);
        }
        p *= pbe1 / vem1.pf;
        if (vem1.c && vem2.c) {
          w += sqr(p);
        }
        for (int ke = je - 2; ke >= 2; --ke) {
          VERTEX ve = pe.vs[ke], vem1 = pe.vs[ke - 1];
          p *= ve.pb / ve.pf;
          if (ve.c && vem1.c) {
            w += sqr(p);
          }
        }
      }

      w = 1.0 / w;

      VEC fsl = edf_eval(*ve.f->edf, ve.n, vec_scale(-1.0, ve.wi));
      VEC I = vec_scale(w, vec_mul(ve.a, fsl));
      sum0 = vec_add(sum0, I);
    }
  }
  //s = 1, t >= 2
  {
    VERTEX vl = pl.vs[0];
    bool cl = vl.c;
    if (cl) {
      for (int je = 1; je < pe.n; ++je) {
        VERTEX ve = pe.vs[je];
        bool ce = ve.c;
        if (ce) {
          VEC d = vec_sub(vl.x, ve.x);
          double t = vec_len(d);
          VEC dl = vec_scale(-1.0 / t, d);
          VEC de = vec_scale(1.0 / t, d);
          if (bvh_visible(scene.bvh_nodes, ve.x, de, t)) {
            VERTEX vem1 = pe.vs[je - 1];

            double te = vec_len(vec_sub(ve.x, vem1.x));

            double g = fabs(vec_dot(vl.n, dl) * vec_dot(ve.n, de)) / sqr(t);
            double gl0 = fabs(vec_dot(vl.n, dl)) / sqr(t);
            double ge0 = fabs(vec_dot(ve.n, de)) / sqr(t);
            double ge1 = fabs(vec_dot(vem1.n, ve.wi)) / sqr(te);

            MEDIUM meir, meit, meor, meot;
            media_sides(vec_dot(ve.n, ve.wi), *ve.f->mp, *ve.f->mn, &meir, &meit);
            media_sides(vec_dot(ve.n, de), *ve.f->mp, *ve.f->mn, &meor, &meot);

            VEC tr = transmittance(meot.a, t);

            double psl = edf_pdf(*vl.f->edf, vl.n, dl);
            double pse = bsdf_pdf(*ve.f->bsdf, meir.eta, meit.eta, ve.n, ve.wi, de);

            double pblr0 = luminance(vec_mul(transmittance(meir.a, te), ve.f->bsdf->a));
            double pbla0 = gl0 * pse;
            double pbl0 = pblr0 * pbla0;

            double pbea0 = ge0 * psl;
            double pbe0 = pbea0;

            double pber1 = luminance(vec_mul(tr, ve.f->bsdf->a));
            double pbea1 = ge1 * pse;
            double pbe1 = pber1 * pbea1;

            double w = 1.0;

            {
              double p = 1.0;
              p *= pbl0 / vl.pf;
              w += sqr(p);
            }

            if (je == 1) {
            } else if (je == 2) {
              double p = 1.0;
              p *= pbe0 / ve.pf;
              if (vem1.c) {
                w += sqr(p);
              }
            } else {
              VERTEX vem2 = pe.vs[je - 2];
              double p = 1.0;
              p *= pbe0 / ve.pf;
              if (vem1.c) {
                w += sqr(p);
              }
              p *= pbe1 / vem1.pf;
              if (vem1.c && vem2.c) {
                w += sqr(p);
              }
              for (int ke = je - 2; ke >= 2; --ke) {
                VERTEX ve = pe.vs[ke], vem1 = pe.vs[ke - 1];
                p *= ve.pb / ve.pf;
                if (ve.c && vem1.c) {
                  w += sqr(p);
                }
              }
            }

            w = 1.0 / w;

            VEC fsl = edf_eval(*vl.f->edf, vl.n, dl);
            VEC fse = bsdf_eval(*ve.f->bsdf, meot.eta, meor.eta, ve.n, dl, vec_scale(-1.0, ve.wi));
            VEC I = vec_scale(w * g, vec_mul(vec_mul(vec_mul(vec_mul(vl.a, ve.a), fsl), fse), tr));
            sum0 = vec_add(sum0, I);
          }
        }
      }
    }
  }
  //s >= 2, t >= 2
  for (int jl = 1; jl < pl.n; ++jl) {
    VERTEX vl = pl.vs[jl];
    bool cl = vl.c;
    if (cl) {
      for (int je = 1; je < pe.n; ++je) {
        VERTEX ve = pe.vs[je];
        bool ce = ve.c;
        if (ce) {
          VEC d = vec_sub(vl.x, ve.x);
          double t = vec_len(d);
          VEC dl = vec_scale(-1.0 / t, d);
          VEC de = vec_scale(1.0 / t, d);
          if (bvh_visible(scene.bvh_nodes, ve.x, de, t)) {
            VERTEX vlm1 = pl.vs[jl - 1];
            VERTEX vem1 = pe.vs[je - 1];

            double tl = vec_len(vec_sub(vl.x, vlm1.x));
            double te = vec_len(vec_sub(ve.x, vem1.x));

            double g = fabs(vec_dot(vl.n, dl) * vec_dot(ve.n, de)) / sqr(t);
            double gl0 = fabs(vec_dot(vl.n, dl)) / sqr(t);
            double gl1 = fabs(vec_dot(vlm1.n, vl.wi)) / sqr(tl);
            double ge0 = fabs(vec_dot(ve.n, de)) / sqr(t);
            double ge1 = fabs(vec_dot(vem1.n, ve.wi)) / sqr(te);

            MEDIUM mlir, mlit, mlor, mlot, meir, meit, meor, meot;
            media_sides(vec_dot(vl.n, vl.wi), *vl.f->mp, *vl.f->mn, &mlir, &mlit);
            media_sides(vec_dot(vl.n, dl), *vl.f->mp, *vl.f->mn, &mlor, &mlot);
            media_sides(vec_dot(ve.n, ve.wi), *ve.f->mp, *ve.f->mn, &meir, &meit);
            media_sides(vec_dot(ve.n, de), *ve.f->mp, *ve.f->mn, &meor, &meot);

            VEC tr = transmittance(meot.a, t);

            double psl = bsdf_pdf(*vl.f->bsdf, mlir.eta, mlit.eta, vl.n, vl.wi, dl);
            double pse = bsdf_pdf(*ve.f->bsdf, meir.eta, meit.eta, ve.n, ve.wi, de);

            double pblr0 = luminance(vec_mul(transmittance(meir.a, te), ve.f->bsdf->a));
            double pbla0 = gl0 * pse;
            double pbl0 = pblr0 * pbla0;

            double pblr1 = luminance(vec_mul(tr, vl.f->bsdf->a));
            double pbla1 = gl1 * psl;
            double pbl1 = pblr1 * pbla1;

            double pber0 = luminance(vec_mul(transmittance(mlir.a, tl), vl.f->bsdf->a));
            double pbea0 = ge0 * psl;
            double pbe0 = pber0 * pbea0;

            double pber1 = luminance(vec_mul(tr, ve.f->bsdf->a));
            double pbea1 = ge1 * pse;
            double pbe1 = pber1 * pbea1;

            double w = 1.0;

            if (jl == 1) {
              double p = 1.0;
              p *= pbl0 / vl.pf;
              if (vlm1.c) {
                w += sqr(p);
              }
              p *= pbl1 / vlm1.pf;
              if (vlm1.c) {
                w += sqr(p);
              }
            } else if (jl == 2) {
              VERTEX vlm2 = pl.vs[jl - 2];
              double p = 1.0;
              p *= pbl0 / vl.pf;
              if (vlm1.c) {
                w += sqr(p);
              }
              p *= pbl1 / vlm1.pf;
              if (vlm1.c && vlm2.c) {
                w += sqr(p);
              }
              p *= vlm2.pb / vlm2.pf;
              if (vlm2.c) {
                w += sqr(p);
              }
            } else {
              VERTEX vlm2 = pl.vs[jl - 2];
              double p = 1.0;
              p *= pbl0 / vl.pf;
              if (vlm1.c) {
                w += sqr(p);
              }
              p *= pbl1 / vlm1.pf;
              if (vlm1.c && vlm2.c) {
                w += sqr(p);
              }
              for (int kl = jl - 2; kl >= 2; --kl) {
                VERTEX vl = pl.vs[kl], vlm1 = pl.vs[kl - 1];
                p *= vl.pb / vl.pf;
                if (vl.c && vlm1.c) {
                  w += sqr(p);
                }
              }
              VERTEX vl1 = pl.vs[1], vl0 = pl.vs[0];
              p *= vl1.pb / vl1.pf;
              if (vl1.c && vl0.c) {
                w += sqr(p);
              }
              p *= vl0.pb / vl0.pf;
              if (vl0.c) {
                w += sqr(p);
              }
            }

            if (je == 1) {
            } else if (je == 2) {
              double p = 1.0;
              p *= pbe0 / ve.pf;
              if (vem1.c) {
                w += sqr(p);
              }
            } else {
              VERTEX vem2 = pe.vs[je - 2];
              double p = 1.0;
              p *= pbe0 / ve.pf;
              if (vem1.c) {
                w += sqr(p);
              }
              p *= pbe1 / vem1.pf;
              if (vem1.c && vem2.c) {
                w += sqr(p);
              }
              for (int ke = je - 2; ke >= 2; --ke) {
                VERTEX ve = pe.vs[ke], vem1 = pe.vs[ke - 1];
                p *= ve.pb / ve.pf;
                if (ve.c && vem1.c) {
                  w += sqr(p);
                }
              }
            }

            w = 1.0 / w;

            VEC fsl = bsdf_eval(*vl.f->bsdf, mlir.eta, mlit.eta, vl.n, vl.wi, dl);
            VEC fse = bsdf_eval(*ve.f->bsdf, meot.eta, meor.eta, ve.n, dl, vec_scale(-1.0, ve.wi));
            VEC I = vec_scale(w * g, vec_mul(vec_mul(vec_mul(vec_mul(vl.a, ve.a), fsl), fse), tr));
            sum0 = vec_add(sum0, I);
          }
        }
      }
    }
  }
  return sum0;
}
void bpt(CONF conf, SCENE scene, IMAGE img, PRNG_STATE *ss) {
  #pragma omp parallel
  {
    PATH pl, pe;
    pl.vs = ALLOC(VERTEX, conf.num_verts);
    pe.vs = ALLOC(VERTEX, conf.num_verts);
    PRNG_STATE *s = ss + omp_get_thread_num();
    for (int i = omp_get_thread_num(); i < img.w * img.h; i += omp_get_max_threads()) {
      int ix = i % img.w, iy = i / img.w;
      VEC sum0 = vec_inita(0.0);
      for (int j = 0; j < conf.num_samples; ++j) {
        {
          FACE *f;
          VEC n, o, d, fs;
          double pa, ps;
          bool c;
          sample_light(scene, s, &f, &n, &o, &d, &fs, &pa, &ps, &c);
          trace(conf, scene, s, &pl, f, n, o, d, fs, pa, ps, c, true);
        }
        {
          FACE *f;
          VEC n, o, d, fs;
          double pa, ps;
          bool c;
          sample_eye(conf, s, ix, iy, &f, &n, &o, &d, &fs, &pa, &ps, &c);
          trace(conf, scene, s, &pe, f, n, o, d, fs, pa, ps, c, false);
        }
        sum0 = vec_add(sum0, bpt_vc(scene, pl, pe));
      }
      sum0 = vec_scale(1.0 / conf.num_samples, sum0);
      sum0 = vec_scale(1.0 / (conf.wi * conf.hi), sum0);
      img.data[3 * i + 0] += sum0.x;
      img.data[3 * i + 1] += sum0.y;
      img.data[3 * i + 2] += sum0.z;
    }
  }
}
IMAGE render_bpt(CONF conf, SCENE scene) {
  IMAGE img = image_init(3, conf.wi, conf.hi);
  omp_set_num_threads(conf.num_threads);
  PRNG_STATE ss[omp_get_max_threads()];
  for (int i = 0; i < omp_get_max_threads(); ++i) {
    ss[i] = prng_seed(SEED);
    for (int j = 0; j < i; ++j) {
      ss[i] = prng_jump(ss[i]);
    }
  }
  bpt(conf, scene, img, ss);
  return img;
}
int main(int argc, char **argv) {
  double t;

  t = time_get();
  CONF conf = conf_read(argv[2]);
  t = time_get() - t;
  printf("read render configuration file: %f\n", t);

  t = time_get();
  SCENE scene = scene_read(argc, argv);
  t = time_get() - t;
  printf("read scene files: %f\n", t);

  t = time_get();
  bvh_make(conf, &scene);
  t = time_get() - t;
  printf("made bvh: %f\n", t);

  t = time_get();
  scene_prep(&scene);
  t = time_get() - t;
  printf("prepared light sampler: %f\n", t);

  t = time_get();
  IMAGE img = render_vcm(conf, scene);
  t = time_get() - t;
  printf("image rendered: %f\n", t);

  t = time_get();
  image_save_pfm(img, argv[1]);
  t = time_get() - t;
  printf("image saved: %f\n", t);

  puts("\a");
}
