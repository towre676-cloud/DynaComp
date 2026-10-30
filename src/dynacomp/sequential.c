#include <float.h>
#include <fenv.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#pragma STDC FENV_ACCESS ON
#pragma STDC FP_CONTRACT OFF
#if FLT_RADIX != 2 || DBL_MANT_DIG != 53 || DBL_MAX_EXP != 1024 || FLT_EVAL_METHOD != 0
#error Requires directly evaluated IEEE binary64
#endif
/* Volatile stores enforce each product and addition, even on compilers that
   ignore STDC FP_CONTRACT. State and input products share one ordered dot. */
static double dot(const double *a, const double *x, size_t k, size_t q, size_t j) {
    volatile double sum=0.0;
    for(size_t i=0;i<k;i++) {
        volatile double product=a[i]*x[i*q+j];
        sum=sum+product;
    }
    return sum;
}
#ifdef _WIN32
__declspec(dllexport)
#endif
int run(size_t n,size_t r,size_t m,size_t p,size_t steps,size_t q,
        const double *e,const double *g,const double *h,const double *d,
        const double *x,const double *u,double *out) {
    if(sizeof(double)!=8 || fegetround()!=FE_TONEAREST) return 1;
    volatile double tiny=0x1p-1074, one=1.0, normal=0x1p-1022, half=0.5;
    volatile double transported=tiny*one, below=normal*half;
    double transported_copy=transported, below_copy=below;
    uint64_t transported_bits, below_bits;
    memcpy(&transported_bits,&transported_copy,8);
    memcpy(&below_bits,&below_copy,8);
    /* Compare encodings: an FP comparison under DAZ could regard both
       the expected subnormal and a flushed zero as equal. */
    if(transported_bits!=1 || below_bits!=UINT64_C(0x0008000000000000)) return 2;
    if(!q || !r || r>SIZE_MAX/q/sizeof(double)) return 3;
    double *state=malloc(r*q*sizeof(double)), *next=malloc(r*q*sizeof(double));
    if(!state || !next) {free(state);free(next);return 4;}
    for(size_t i=0;i<r;i++) for(size_t j=0;j<q;j++)
        state[i*q+j]=dot(e+i*n,x,n,q,j);
    for(size_t t=0;t<=steps;t++) {
        for(size_t i=0;i<p;i++) for(size_t j=0;j<q;j++)
            out[(t*p+i)*q+j]=dot(d+i*r,state,r,q,j);
        if(t==steps) break;
        for(size_t i=0;i<r;i++) for(size_t j=0;j<q;j++) {
            volatile double sum=dot(g+i*r,state,r,q,j);
            for(size_t k=0;k<m;k++) {
                volatile double product=h[i*m+k]*u[(t*m+k)*q+j];
                sum=sum+product;
            }
            next[i*q+j]=sum;
        }
        double *swap=state;state=next;next=swap;
    }
    free(state);free(next);return 0;
}
