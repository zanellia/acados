/* This file was automatically generated by CasADi.
   The CasADi copyright holders make no ownership claim of its contents. */
#ifdef __cplusplus
extern "C" {
#endif

#ifdef CODEGEN_PREFIX
  #define NAMESPACE_CONCAT(NS, ID) _NAMESPACE_CONCAT(NS, ID)
  #define _NAMESPACE_CONCAT(NS, ID) NS ## ID
  #define CASADI_PREFIX(ID) NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)
#else /* CODEGEN_PREFIX */
  #define CASADI_PREFIX(ID) ls_res_end_quadcopter_ ## ID
#endif /* CODEGEN_PREFIX */

#include <math.h>

#ifndef real_t
#define real_t double
#endif /* real_t */

#define to_double(x) (double) x
#define to_int(x) (int) x
#define CASADI_CAST(x,y) (x) y
/* Pre-c99 compatibility */
#if __STDC_VERSION__ < 199901L
real_t CASADI_PREFIX(fmin)(real_t x, real_t y) { return x<y ? x : y;}
#define fmin(x,y) CASADI_PREFIX(fmin)(x,y)
real_t CASADI_PREFIX(fmax)(real_t x, real_t y) { return x>y ? x : y;}
#define fmax(x,y) CASADI_PREFIX(fmax)(x,y)
#endif

#define PRINTF printf
real_t CASADI_PREFIX(sq)(real_t x) { return x*x;}
#define sq(x) CASADI_PREFIX(sq)(x)

real_t CASADI_PREFIX(sign)(real_t x) { return x<0 ? -1 : x>0 ? 1 : x;}
#define sign(x) CASADI_PREFIX(sign)(x)

static const int CASADI_PREFIX(s0)[15] = {11, 1, 0, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
#define s0 CASADI_PREFIX(s0)
static const int CASADI_PREFIX(s1)[8] = {4, 1, 0, 4, 0, 1, 2, 3};
#define s1 CASADI_PREFIX(s1)
static const int CASADI_PREFIX(s2)[135] = {11, 11, 0, 11, 22, 33, 44, 55, 66, 77, 88, 99, 110, 121, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
#define s2 CASADI_PREFIX(s2)
/* ls_res_end_Fun */
int ls_res_end_Fun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem) {
  real_t a0=3.1622776601683795e+00;
  real_t a1=arg[0] ? arg[0][0] : 0;
  a1=(a0*a1);
  if (res[0]!=0) res[0][0]=a1;
  a1=arg[0] ? arg[0][1] : 0;
  a1=(a0*a1);
  if (res[0]!=0) res[0][1]=a1;
  a1=arg[0] ? arg[0][2] : 0;
  a1=(a0*a1);
  if (res[0]!=0) res[0][2]=a1;
  a1=arg[0] ? arg[0][3] : 0;
  a1=(a0*a1);
  if (res[0]!=0) res[0][3]=a1;
  a1=1.0000000000000001e-01;
  real_t a2=arg[0] ? arg[0][4] : 0;
  a2=(a1*a2);
  if (res[0]!=0) res[0][4]=a2;
  a2=arg[0] ? arg[0][5] : 0;
  a2=(a1*a2);
  if (res[0]!=0) res[0][5]=a2;
  a2=arg[0] ? arg[0][6] : 0;
  a2=(a1*a2);
  if (res[0]!=0) res[0][6]=a2;
  a2=arg[0] ? arg[0][7] : 0;
  a2=(a1*a2);
  if (res[0]!=0) res[0][7]=a2;
  a2=arg[0] ? arg[0][8] : 0;
  a2=(a1*a2);
  if (res[0]!=0) res[0][8]=a2;
  a2=arg[0] ? arg[0][9] : 0;
  a2=(a1*a2);
  if (res[0]!=0) res[0][9]=a2;
  a2=arg[0] ? arg[0][10] : 0;
  a2=(a1*a2);
  if (res[0]!=0) res[0][10]=a2;
  if (res[1]!=0) res[1][0]=a0;
  a2=0.;
  if (res[1]!=0) res[1][1]=a2;
  if (res[1]!=0) res[1][2]=a2;
  if (res[1]!=0) res[1][3]=a2;
  if (res[1]!=0) res[1][4]=a2;
  if (res[1]!=0) res[1][5]=a2;
  if (res[1]!=0) res[1][6]=a2;
  if (res[1]!=0) res[1][7]=a2;
  if (res[1]!=0) res[1][8]=a2;
  if (res[1]!=0) res[1][9]=a2;
  if (res[1]!=0) res[1][10]=a2;
  if (res[1]!=0) res[1][11]=a2;
  if (res[1]!=0) res[1][12]=a0;
  if (res[1]!=0) res[1][13]=a2;
  if (res[1]!=0) res[1][14]=a2;
  if (res[1]!=0) res[1][15]=a2;
  if (res[1]!=0) res[1][16]=a2;
  if (res[1]!=0) res[1][17]=a2;
  if (res[1]!=0) res[1][18]=a2;
  if (res[1]!=0) res[1][19]=a2;
  if (res[1]!=0) res[1][20]=a2;
  if (res[1]!=0) res[1][21]=a2;
  if (res[1]!=0) res[1][22]=a2;
  if (res[1]!=0) res[1][23]=a2;
  if (res[1]!=0) res[1][24]=a0;
  if (res[1]!=0) res[1][25]=a2;
  if (res[1]!=0) res[1][26]=a2;
  if (res[1]!=0) res[1][27]=a2;
  if (res[1]!=0) res[1][28]=a2;
  if (res[1]!=0) res[1][29]=a2;
  if (res[1]!=0) res[1][30]=a2;
  if (res[1]!=0) res[1][31]=a2;
  if (res[1]!=0) res[1][32]=a2;
  if (res[1]!=0) res[1][33]=a2;
  if (res[1]!=0) res[1][34]=a2;
  if (res[1]!=0) res[1][35]=a2;
  if (res[1]!=0) res[1][36]=a0;
  if (res[1]!=0) res[1][37]=a2;
  if (res[1]!=0) res[1][38]=a2;
  if (res[1]!=0) res[1][39]=a2;
  if (res[1]!=0) res[1][40]=a2;
  if (res[1]!=0) res[1][41]=a2;
  if (res[1]!=0) res[1][42]=a2;
  if (res[1]!=0) res[1][43]=a2;
  if (res[1]!=0) res[1][44]=a2;
  if (res[1]!=0) res[1][45]=a2;
  if (res[1]!=0) res[1][46]=a2;
  if (res[1]!=0) res[1][47]=a2;
  if (res[1]!=0) res[1][48]=a1;
  if (res[1]!=0) res[1][49]=a2;
  if (res[1]!=0) res[1][50]=a2;
  if (res[1]!=0) res[1][51]=a2;
  if (res[1]!=0) res[1][52]=a2;
  if (res[1]!=0) res[1][53]=a2;
  if (res[1]!=0) res[1][54]=a2;
  if (res[1]!=0) res[1][55]=a2;
  if (res[1]!=0) res[1][56]=a2;
  if (res[1]!=0) res[1][57]=a2;
  if (res[1]!=0) res[1][58]=a2;
  if (res[1]!=0) res[1][59]=a2;
  if (res[1]!=0) res[1][60]=a1;
  if (res[1]!=0) res[1][61]=a2;
  if (res[1]!=0) res[1][62]=a2;
  if (res[1]!=0) res[1][63]=a2;
  if (res[1]!=0) res[1][64]=a2;
  if (res[1]!=0) res[1][65]=a2;
  if (res[1]!=0) res[1][66]=a2;
  if (res[1]!=0) res[1][67]=a2;
  if (res[1]!=0) res[1][68]=a2;
  if (res[1]!=0) res[1][69]=a2;
  if (res[1]!=0) res[1][70]=a2;
  if (res[1]!=0) res[1][71]=a2;
  if (res[1]!=0) res[1][72]=a1;
  if (res[1]!=0) res[1][73]=a2;
  if (res[1]!=0) res[1][74]=a2;
  if (res[1]!=0) res[1][75]=a2;
  if (res[1]!=0) res[1][76]=a2;
  if (res[1]!=0) res[1][77]=a2;
  if (res[1]!=0) res[1][78]=a2;
  if (res[1]!=0) res[1][79]=a2;
  if (res[1]!=0) res[1][80]=a2;
  if (res[1]!=0) res[1][81]=a2;
  if (res[1]!=0) res[1][82]=a2;
  if (res[1]!=0) res[1][83]=a2;
  if (res[1]!=0) res[1][84]=a1;
  if (res[1]!=0) res[1][85]=a2;
  if (res[1]!=0) res[1][86]=a2;
  if (res[1]!=0) res[1][87]=a2;
  if (res[1]!=0) res[1][88]=a2;
  if (res[1]!=0) res[1][89]=a2;
  if (res[1]!=0) res[1][90]=a2;
  if (res[1]!=0) res[1][91]=a2;
  if (res[1]!=0) res[1][92]=a2;
  if (res[1]!=0) res[1][93]=a2;
  if (res[1]!=0) res[1][94]=a2;
  if (res[1]!=0) res[1][95]=a2;
  if (res[1]!=0) res[1][96]=a1;
  if (res[1]!=0) res[1][97]=a2;
  if (res[1]!=0) res[1][98]=a2;
  if (res[1]!=0) res[1][99]=a2;
  if (res[1]!=0) res[1][100]=a2;
  if (res[1]!=0) res[1][101]=a2;
  if (res[1]!=0) res[1][102]=a2;
  if (res[1]!=0) res[1][103]=a2;
  if (res[1]!=0) res[1][104]=a2;
  if (res[1]!=0) res[1][105]=a2;
  if (res[1]!=0) res[1][106]=a2;
  if (res[1]!=0) res[1][107]=a2;
  if (res[1]!=0) res[1][108]=a1;
  if (res[1]!=0) res[1][109]=a2;
  if (res[1]!=0) res[1][110]=a2;
  if (res[1]!=0) res[1][111]=a2;
  if (res[1]!=0) res[1][112]=a2;
  if (res[1]!=0) res[1][113]=a2;
  if (res[1]!=0) res[1][114]=a2;
  if (res[1]!=0) res[1][115]=a2;
  if (res[1]!=0) res[1][116]=a2;
  if (res[1]!=0) res[1][117]=a2;
  if (res[1]!=0) res[1][118]=a2;
  if (res[1]!=0) res[1][119]=a2;
  if (res[1]!=0) res[1][120]=a1;
  return 0;
}

void ls_res_end_Fun_incref(void) {
}

void ls_res_end_Fun_decref(void) {
}

int ls_res_end_Fun_n_in(void) { return 2;}

int ls_res_end_Fun_n_out(void) { return 2;}

const char* ls_res_end_Fun_name_in(int i){
  switch (i) {
    case 0: return "i0";
    case 1: return "i1";
    default: return 0;
  }
}

const char* ls_res_end_Fun_name_out(int i){
  switch (i) {
    case 0: return "o0";
    case 1: return "o1";
    default: return 0;
  }
}

const int* ls_res_end_Fun_sparsity_in(int i) {
  switch (i) {
    case 0: return s0;
    case 1: return s1;
    default: return 0;
  }
}

const int* ls_res_end_Fun_sparsity_out(int i) {
  switch (i) {
    case 0: return s0;
    case 1: return s2;
    default: return 0;
  }
}

int ls_res_end_Fun_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w) {
  if (sz_arg) *sz_arg = 2;
  if (sz_res) *sz_res = 2;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 3;
  return 0;
}


#ifdef __cplusplus
} /* extern "C" */
#endif
