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
static const int CASADI_PREFIX(s2)[18] = {14, 1, 0, 14, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
#define s2 CASADI_PREFIX(s2)
static const int CASADI_PREFIX(s3)[168] = {14, 11, 0, 14, 28, 42, 56, 70, 84, 98, 112, 126, 140, 154, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
#define s3 CASADI_PREFIX(s3)
/* ls_res_end_Fun */
int ls_res_end_Fun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem) {
  real_t a0=arg[0] ? arg[0][0] : 0;
  real_t a1=arg[0] ? arg[0][1] : 0;
  real_t a2=(a0*a1);
  real_t a3=arg[0] ? arg[0][2] : 0;
  real_t a4=arg[0] ? arg[0][3] : 0;
  real_t a5=(a3*a4);
  a2=(a2+a5);
  a5=2.;
  a2=(a5*a2);
  real_t a6=sq(a1);
  real_t a7=sq(a3);
  a6=(a6+a7);
  a6=(a5*a6);
  a7=1.;
  a6=(a7-a6);
  real_t a8=atan2(a2,a6);
  real_t a9=2.2360679774997898e+01;
  a8=(a9*a8);
  if (res[0]!=0) res[0][0]=a8;
  a8=(a0*a3);
  real_t a10=(a4*a1);
  a8=(a8-a10);
  a8=(a5*a8);
  a10=asin(a8);
  a10=(a9*a10);
  if (res[0]!=0) res[0][1]=a10;
  a10=(a0*a4);
  real_t a11=(a1*a3);
  a10=(a10+a11);
  a10=(a5*a10);
  a11=sq(a3);
  real_t a12=sq(a4);
  a11=(a11+a12);
  a11=(a5*a11);
  a11=(a7-a11);
  a12=atan2(a10,a11);
  a12=(a9*a12);
  if (res[0]!=0) res[0][2]=a12;
  a12=3.1622776601683791e-02;
  real_t a13=(a12*a0);
  if (res[0]!=0) res[0][3]=a13;
  a13=(a12*a1);
  if (res[0]!=0) res[0][4]=a13;
  a13=(a12*a3);
  if (res[0]!=0) res[0][5]=a13;
  a13=(a12*a4);
  if (res[0]!=0) res[0][6]=a13;
  a13=arg[0] ? arg[0][4] : 0;
  a13=(a12*a13);
  if (res[0]!=0) res[0][7]=a13;
  a13=arg[0] ? arg[0][5] : 0;
  a13=(a12*a13);
  if (res[0]!=0) res[0][8]=a13;
  a13=arg[0] ? arg[0][6] : 0;
  a13=(a12*a13);
  if (res[0]!=0) res[0][9]=a13;
  a13=arg[0] ? arg[0][7] : 0;
  a13=(a12*a13);
  if (res[0]!=0) res[0][10]=a13;
  a13=arg[0] ? arg[0][8] : 0;
  a13=(a12*a13);
  if (res[0]!=0) res[0][11]=a13;
  a13=arg[0] ? arg[0][9] : 0;
  a13=(a12*a13);
  if (res[0]!=0) res[0][12]=a13;
  a13=arg[0] ? arg[0][10] : 0;
  a13=(a12*a13);
  if (res[0]!=0) res[0][13]=a13;
  a13=sq(a2);
  real_t a14=sq(a6);
  a13=(a13+a14);
  a6=(a6/a13);
  a14=(a5*a1);
  a14=(a6*a14);
  a14=(a9*a14);
  if (res[1]!=0) res[1][0]=a14;
  a14=(a5*a3);
  a8=sq(a8);
  a7=(a7-a8);
  a7=sqrt(a7);
  a14=(a14/a7);
  a14=(a9*a14);
  if (res[1]!=0) res[1][1]=a14;
  a14=sq(a10);
  a8=sq(a11);
  a14=(a14+a8);
  a11=(a11/a14);
  a8=(a5*a4);
  a8=(a11*a8);
  a8=(a9*a8);
  if (res[1]!=0) res[1][2]=a8;
  if (res[1]!=0) res[1][3]=a12;
  a8=0.;
  if (res[1]!=0) res[1][4]=a8;
  if (res[1]!=0) res[1][5]=a8;
  if (res[1]!=0) res[1][6]=a8;
  if (res[1]!=0) res[1][7]=a8;
  if (res[1]!=0) res[1][8]=a8;
  if (res[1]!=0) res[1][9]=a8;
  if (res[1]!=0) res[1][10]=a8;
  if (res[1]!=0) res[1][11]=a8;
  if (res[1]!=0) res[1][12]=a8;
  if (res[1]!=0) res[1][13]=a8;
  real_t a15=(a5*a0);
  a15=(a6*a15);
  a2=(a2/a13);
  a13=(a1+a1);
  a13=(a5*a13);
  a13=(a2*a13);
  a15=(a15+a13);
  a15=(a9*a15);
  if (res[1]!=0) res[1][14]=a15;
  a15=(a5*a4);
  a15=(a15/a7);
  a15=(a9*a15);
  a15=(-a15);
  if (res[1]!=0) res[1][15]=a15;
  a15=(a5*a3);
  a15=(a11*a15);
  a15=(a9*a15);
  if (res[1]!=0) res[1][16]=a15;
  if (res[1]!=0) res[1][17]=a8;
  if (res[1]!=0) res[1][18]=a12;
  if (res[1]!=0) res[1][19]=a8;
  if (res[1]!=0) res[1][20]=a8;
  if (res[1]!=0) res[1][21]=a8;
  if (res[1]!=0) res[1][22]=a8;
  if (res[1]!=0) res[1][23]=a8;
  if (res[1]!=0) res[1][24]=a8;
  if (res[1]!=0) res[1][25]=a8;
  if (res[1]!=0) res[1][26]=a8;
  if (res[1]!=0) res[1][27]=a8;
  a15=(a5*a4);
  a15=(a6*a15);
  a13=(a3+a3);
  a13=(a5*a13);
  a2=(a2*a13);
  a15=(a15+a2);
  a15=(a9*a15);
  if (res[1]!=0) res[1][28]=a15;
  a15=(a5*a0);
  a15=(a15/a7);
  a15=(a9*a15);
  if (res[1]!=0) res[1][29]=a15;
  a15=(a5*a1);
  a15=(a11*a15);
  a10=(a10/a14);
  a14=(a3+a3);
  a14=(a5*a14);
  a14=(a10*a14);
  a15=(a15+a14);
  a15=(a9*a15);
  if (res[1]!=0) res[1][30]=a15;
  if (res[1]!=0) res[1][31]=a8;
  if (res[1]!=0) res[1][32]=a8;
  if (res[1]!=0) res[1][33]=a12;
  if (res[1]!=0) res[1][34]=a8;
  if (res[1]!=0) res[1][35]=a8;
  if (res[1]!=0) res[1][36]=a8;
  if (res[1]!=0) res[1][37]=a8;
  if (res[1]!=0) res[1][38]=a8;
  if (res[1]!=0) res[1][39]=a8;
  if (res[1]!=0) res[1][40]=a8;
  if (res[1]!=0) res[1][41]=a8;
  a3=(a5*a3);
  a6=(a6*a3);
  a6=(a9*a6);
  if (res[1]!=0) res[1][42]=a6;
  a1=(a5*a1);
  a1=(a1/a7);
  a1=(a9*a1);
  a1=(-a1);
  if (res[1]!=0) res[1][43]=a1;
  a0=(a5*a0);
  a11=(a11*a0);
  a4=(a4+a4);
  a5=(a5*a4);
  a10=(a10*a5);
  a11=(a11+a10);
  a11=(a9*a11);
  if (res[1]!=0) res[1][44]=a11;
  if (res[1]!=0) res[1][45]=a8;
  if (res[1]!=0) res[1][46]=a8;
  if (res[1]!=0) res[1][47]=a8;
  if (res[1]!=0) res[1][48]=a12;
  if (res[1]!=0) res[1][49]=a8;
  if (res[1]!=0) res[1][50]=a8;
  if (res[1]!=0) res[1][51]=a8;
  if (res[1]!=0) res[1][52]=a8;
  if (res[1]!=0) res[1][53]=a8;
  if (res[1]!=0) res[1][54]=a8;
  if (res[1]!=0) res[1][55]=a8;
  if (res[1]!=0) res[1][56]=a8;
  if (res[1]!=0) res[1][57]=a8;
  if (res[1]!=0) res[1][58]=a8;
  if (res[1]!=0) res[1][59]=a8;
  if (res[1]!=0) res[1][60]=a8;
  if (res[1]!=0) res[1][61]=a8;
  if (res[1]!=0) res[1][62]=a8;
  if (res[1]!=0) res[1][63]=a12;
  if (res[1]!=0) res[1][64]=a8;
  if (res[1]!=0) res[1][65]=a8;
  if (res[1]!=0) res[1][66]=a8;
  if (res[1]!=0) res[1][67]=a8;
  if (res[1]!=0) res[1][68]=a8;
  if (res[1]!=0) res[1][69]=a8;
  if (res[1]!=0) res[1][70]=a8;
  if (res[1]!=0) res[1][71]=a8;
  if (res[1]!=0) res[1][72]=a8;
  if (res[1]!=0) res[1][73]=a8;
  if (res[1]!=0) res[1][74]=a8;
  if (res[1]!=0) res[1][75]=a8;
  if (res[1]!=0) res[1][76]=a8;
  if (res[1]!=0) res[1][77]=a8;
  if (res[1]!=0) res[1][78]=a12;
  if (res[1]!=0) res[1][79]=a8;
  if (res[1]!=0) res[1][80]=a8;
  if (res[1]!=0) res[1][81]=a8;
  if (res[1]!=0) res[1][82]=a8;
  if (res[1]!=0) res[1][83]=a8;
  if (res[1]!=0) res[1][84]=a8;
  if (res[1]!=0) res[1][85]=a8;
  if (res[1]!=0) res[1][86]=a8;
  if (res[1]!=0) res[1][87]=a8;
  if (res[1]!=0) res[1][88]=a8;
  if (res[1]!=0) res[1][89]=a8;
  if (res[1]!=0) res[1][90]=a8;
  if (res[1]!=0) res[1][91]=a8;
  if (res[1]!=0) res[1][92]=a8;
  if (res[1]!=0) res[1][93]=a12;
  if (res[1]!=0) res[1][94]=a8;
  if (res[1]!=0) res[1][95]=a8;
  if (res[1]!=0) res[1][96]=a8;
  if (res[1]!=0) res[1][97]=a8;
  if (res[1]!=0) res[1][98]=a8;
  if (res[1]!=0) res[1][99]=a8;
  if (res[1]!=0) res[1][100]=a8;
  if (res[1]!=0) res[1][101]=a8;
  if (res[1]!=0) res[1][102]=a8;
  if (res[1]!=0) res[1][103]=a8;
  if (res[1]!=0) res[1][104]=a8;
  if (res[1]!=0) res[1][105]=a8;
  if (res[1]!=0) res[1][106]=a8;
  if (res[1]!=0) res[1][107]=a8;
  if (res[1]!=0) res[1][108]=a12;
  if (res[1]!=0) res[1][109]=a8;
  if (res[1]!=0) res[1][110]=a8;
  if (res[1]!=0) res[1][111]=a8;
  if (res[1]!=0) res[1][112]=a8;
  if (res[1]!=0) res[1][113]=a8;
  if (res[1]!=0) res[1][114]=a8;
  if (res[1]!=0) res[1][115]=a8;
  if (res[1]!=0) res[1][116]=a8;
  if (res[1]!=0) res[1][117]=a8;
  if (res[1]!=0) res[1][118]=a8;
  if (res[1]!=0) res[1][119]=a8;
  if (res[1]!=0) res[1][120]=a8;
  if (res[1]!=0) res[1][121]=a8;
  if (res[1]!=0) res[1][122]=a8;
  if (res[1]!=0) res[1][123]=a12;
  if (res[1]!=0) res[1][124]=a8;
  if (res[1]!=0) res[1][125]=a8;
  if (res[1]!=0) res[1][126]=a8;
  if (res[1]!=0) res[1][127]=a8;
  if (res[1]!=0) res[1][128]=a8;
  if (res[1]!=0) res[1][129]=a8;
  if (res[1]!=0) res[1][130]=a8;
  if (res[1]!=0) res[1][131]=a8;
  if (res[1]!=0) res[1][132]=a8;
  if (res[1]!=0) res[1][133]=a8;
  if (res[1]!=0) res[1][134]=a8;
  if (res[1]!=0) res[1][135]=a8;
  if (res[1]!=0) res[1][136]=a8;
  if (res[1]!=0) res[1][137]=a8;
  if (res[1]!=0) res[1][138]=a12;
  if (res[1]!=0) res[1][139]=a8;
  if (res[1]!=0) res[1][140]=a8;
  if (res[1]!=0) res[1][141]=a8;
  if (res[1]!=0) res[1][142]=a8;
  if (res[1]!=0) res[1][143]=a8;
  if (res[1]!=0) res[1][144]=a8;
  if (res[1]!=0) res[1][145]=a8;
  if (res[1]!=0) res[1][146]=a8;
  if (res[1]!=0) res[1][147]=a8;
  if (res[1]!=0) res[1][148]=a8;
  if (res[1]!=0) res[1][149]=a8;
  if (res[1]!=0) res[1][150]=a8;
  if (res[1]!=0) res[1][151]=a8;
  if (res[1]!=0) res[1][152]=a8;
  if (res[1]!=0) res[1][153]=a12;
  a8=arg[2] ? arg[2][0] : 0;
  a8=(a9*a8);
  if (res[2]!=0) res[2][0]=a8;
  a8=arg[2] ? arg[2][1] : 0;
  a8=(a9*a8);
  if (res[2]!=0) res[2][1]=a8;
  a8=arg[2] ? arg[2][2] : 0;
  a9=(a9*a8);
  if (res[2]!=0) res[2][2]=a9;
  a9=arg[2] ? arg[2][3] : 0;
  a9=(a12*a9);
  if (res[2]!=0) res[2][3]=a9;
  a9=arg[2] ? arg[2][4] : 0;
  a9=(a12*a9);
  if (res[2]!=0) res[2][4]=a9;
  a9=arg[2] ? arg[2][5] : 0;
  a9=(a12*a9);
  if (res[2]!=0) res[2][5]=a9;
  a9=arg[2] ? arg[2][6] : 0;
  a9=(a12*a9);
  if (res[2]!=0) res[2][6]=a9;
  a9=arg[2] ? arg[2][7] : 0;
  a9=(a12*a9);
  if (res[2]!=0) res[2][7]=a9;
  a9=arg[2] ? arg[2][8] : 0;
  a9=(a12*a9);
  if (res[2]!=0) res[2][8]=a9;
  a9=arg[2] ? arg[2][9] : 0;
  a9=(a12*a9);
  if (res[2]!=0) res[2][9]=a9;
  a9=arg[2] ? arg[2][10] : 0;
  a9=(a12*a9);
  if (res[2]!=0) res[2][10]=a9;
  a9=arg[2] ? arg[2][11] : 0;
  a9=(a12*a9);
  if (res[2]!=0) res[2][11]=a9;
  a9=arg[2] ? arg[2][12] : 0;
  a9=(a12*a9);
  if (res[2]!=0) res[2][12]=a9;
  a9=arg[2] ? arg[2][13] : 0;
  a12=(a12*a9);
  if (res[2]!=0) res[2][13]=a12;
  return 0;
}

void ls_res_end_Fun_incref(void) {
}

void ls_res_end_Fun_decref(void) {
}

int ls_res_end_Fun_n_in(void) { return 3;}

int ls_res_end_Fun_n_out(void) { return 3;}

const char* ls_res_end_Fun_name_in(int i){
  switch (i) {
    case 0: return "i0";
    case 1: return "i1";
    case 2: return "i2";
    default: return 0;
  }
}

const char* ls_res_end_Fun_name_out(int i){
  switch (i) {
    case 0: return "o0";
    case 1: return "o1";
    case 2: return "o2";
    default: return 0;
  }
}

const int* ls_res_end_Fun_sparsity_in(int i) {
  switch (i) {
    case 0: return s0;
    case 1: return s1;
    case 2: return s2;
    default: return 0;
  }
}

const int* ls_res_end_Fun_sparsity_out(int i) {
  switch (i) {
    case 0: return s2;
    case 1: return s3;
    case 2: return s2;
    default: return 0;
  }
}

int ls_res_end_Fun_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w) {
  if (sz_arg) *sz_arg = 3;
  if (sz_res) *sz_res = 3;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 16;
  return 0;
}


#ifdef __cplusplus
} /* extern "C" */
#endif
