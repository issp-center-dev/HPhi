#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "DefCommon.h"
struct BindStruct;
#include "struct.h"
#include "symmetry_kondo.h"
#include "symmetry_terms.h"
#include "symmetry_kondo_terms.h"
#include "symmetry_diagonal.h"

FILE *stdoutMPI = NULL;
int nproc = 1, myrank = 0;
static void require(int ok, const char *label)
{
  if (!ok) { fprintf(stderr, "FAIL: %s\n", label); exit(1); }
}
static void close_value(double complex a, double complex b, const char *label)
{
  require(isfinite(creal(a)) && isfinite(cimag(a)) && cabs(a-b) < 1e-12, label);
}
static struct DefineList definition(int model, int *local, unsigned int n)
{
  struct DefineList d = {0};
  unsigned int i;
  d.iCalcModel = model == KondoNConserved ? Kondo : model; d.Nsite = n; d.LocSpn = local;
  for (i=0; i<n; ++i) d.NLocSpn += local[i] == LOCSPIN;
  d.NCond = 1; d.Total2Sz = (int)(d.NLocSpn + 1) % 2;
  require(NormalizeSymmetryKondoQuantumNumbers(&d, model != KondoGC,
      model == Kondo, 0, 0) == 0, "normalized fixture");
  return d;
}
/* Independent tensor CAR oracle. Local basis: |0>, |up>, |down>, |up down>.
 * These literal creation matrices, their transposes, and parity strings are
 * the only fermion machinery used to build the expected matrices. */
static const double create4[2][4][4] = {
  {{0,0,0,0},{1,0,0,0},{0,0,0,0},{0,0,1,0}},
  {{0,0,0,0},{0,0,0,0},{1,0,0,0},{0,-1,0,0}}
};
static int physical(const struct DefineList *d, unsigned int state)
{
  unsigned int i;
  for (i=0;i<d->Nsite;++i)
    if (d->LocSpn[i] == LOCSPIN && ((state>>(2*i))&3) != 1 &&
        ((state>>(2*i))&3) != 2) return 0;
  return 1;
}
static void car(const struct DefineList *d, int orbital, int creation,
                double complex v[64])
{
  double complex next[64] = {0};
  unsigned int in, out, s, dim=1U<<(2*d->Nsite);
  unsigned int site=(unsigned int)orbital/2, spin=(unsigned int)orbital%2;
  for(in=0;in<dim;++in) {
    unsigned int digit=(in>>(2*site))&3;
    int parity=1;
    for(s=0;s<site;++s) {
      unsigned int x=(in>>(2*s))&3;
      if(x==1 || x==2) parity=-parity;
    }
    for(out=0;out<4;++out) {
      double x=creation ? create4[spin][out][digit] : create4[spin][digit][out];
      next[(in & ~(3U<<(2*site))) | (out<<(2*site))] += parity*x*v[in];
    }
  }
  memcpy(v,next,sizeof(next));
}
static void oracle(const struct DefineList *d, const struct SymmetryTerm *t,
                   unsigned int in, double complex v[64])
{
  unsigned int f, out, dim=1U<<(2*d->Nsite);
  memset(v,0,64*sizeof(*v)); v[in]=t->value;
  for(f=t->factors;f-->0;) {
    const int *x=t->index+4*f;
    car(d,2*x[2]+x[3],0,v); car(d,2*x[0]+x[1],1,v);
  }
  for(out=0;out<dim;++out) if(!physical(d,out)) v[out]=0;
}
struct CanonicalMatrix {
  const struct DefineList *def;
  unsigned int input;
  double complex v[64];
  unsigned int count;
};
static int canonical_matrix(const struct SymmetryKondoMonomial *m, void *context)
{
  struct CanonicalMatrix *a=context;
  double complex v[64]={0}, next[64];
  unsigned int f, in, dim=1U<<(2*a->def->Nsite);
  require(m->nlocal<=2 && m->ncreate<=2 && m->nannihilate<=2,"bounded monomial");
  v[a->input]=m->value;
  for(f=m->nannihilate;f-->0;) car(a->def,m->annihilate[f],0,v);
  for(f=m->ncreate;f-->0;) car(a->def,m->create[f],1,v);
  for(f=0;f<m->nlocal;++f) {
    unsigned int site=(unsigned int)m->local_site[f];
    require(m->local_out[f] || m->local_in[f],"no E00 in independent basis");
    require(f==0 || m->local_site[f-1]<m->local_site[f],"sorted local keys");
    memset(next,0,sizeof(next));
    for(in=0;in<dim;++in) if(((in>>(2*site))&3)==(1U<<m->local_in[f]))
      next[(in&~(3U<<(2*site))) | (1U<<(2*site+m->local_out[f]))]+=v[in];
    memcpy(v,next,sizeof(v));
  }
  for(in=0;in<dim;++in) a->v[in]+=v[in];
  ++a->count;
  return 0;
}
static void check_term(const struct DefineList *d, const struct SymmetryTerm *t)
{
  unsigned int in,out,dim=1U<<(2*d->Nsite);
  for(in=0;in<dim;++in) if(physical(d,in)) {
    double complex v[64],value=0;
    unsigned long target=999;
    int status;
    struct CanonicalMatrix canonical={0};
    oracle(d,t,in,v);
    canonical.def=d; canonical.input=in;
    require(CanonicalizeSymmetryKondoTerm(d,t,NULL,canonical_matrix,&canonical)==0,
            "canonicalize Kondo term");
    for(out=0;out<dim;++out)
      close_value(canonical.v[out],v[out],"canonical polynomial CAR matrix");
    status=ApplySymmetryTerm(d,t,in,&target,&value);
    require(status>=0,"Kondo term application accepted");
    for(out=0;out<dim;++out)
      close_value(status && out==target ? value : 0,v[out],"CAR full matrix");
  }
}
static void test_empty_and_application(void)
{
  int local[2]={LOCSPIN,ITINERANT}, identity[2]={0,1}, *group[]={identity};
  int models[]={Kondo,KondoNConserved,KondoGC};
  unsigned int m,a,b,c,e;
  for(m=0;m<3;++m) {
    struct DefineList d=definition(models[m],local,2);
    double diagonal=9;
    d.NSymTrans=1; d.SymTrans=group;
    require(SymmetryUsesExtendedTerms(&d)==1,"empty Kondo always uses extended terms");
    require(ValidateSymmetryTerms(&d)==0,"empty Kondo polynomial accepted");
    require(ValidateSymmetryKondoTerms(&d)==0,"empty dedicated validator accepted");
    require(EvaluateSymmetryStateDiagonal(&d,1,&diagonal)==0 && diagonal==0,
            "empty Kondo diagonal");
  }
  {
    struct DefineList d=definition(KondoGC,local,2);
    struct SymmetryTerm t={2,{0},0.3+0.2*I};
    unsigned long out;
    double complex value;
    double diagonal;
    for(a=0;a<4;++a) for(b=0;b<4;++b) for(c=0;c<4;++c) for(e=0;e<4;++e) {
      int x[]={a/2,a%2,b/2,b%2,c/2,c%2,e/2,e%2};
      memcpy(t.index,x,sizeof(x)); check_term(&d,&t);
    }
    t.factors=1;
    for(a=0;a<4;++a) for(b=0;b<4;++b) {
      t.index[0]=a/2; t.index[1]=a%2; t.index[2]=b/2; t.index[3]=b%2;
      check_term(&d,&t);
    }
    require(ApplySymmetryTerm(&d,&t,0,&out,&value)==-1,"unphysical input rejected");
    require(EvaluateSymmetryStateDiagonal(&d,0,&diagonal)==-1,"unphysical empty diagonal rejected");
  }
}
static int reject_monomial(const struct SymmetryKondoMonomial *m, void *context)
{
  (void)m; (void)context; return 7;
}
static void test_validation(void)
{
  int local[]={LOCSPIN,ITINERANT,LOCSPIN};
  int identity[]={0,1,2}, swap[]={2,1,0}, *group[]={identity,swap};
  int rows[4][4]={{0,0,0,0},{0,1,0,1},{0,0,0,1},{0,0,0,1}};
  int *ptrs[]={rows[0],rows[1],rows[2],rows[3]};
  double complex values[]={0.4,0.4,0.3,-0.3};
  struct DefineList d=definition(KondoGC,local,3);
  struct SymmetryTerm t={1,{0,0,0,1},1};
  struct CanonicalMatrix a={0};
  d.NSymTrans=2; d.SymTrans=group;
  d.EDNTransfer=4; d.EDGeneralTransfer=ptrs; d.EDParaGeneralTransfer=values;
  require(ValidateSymmetryKondoTerms(&d)==0,"local density identity and duplicate cancellation");
  values[1]=0.5;
  require(ValidateSymmetryTerms(&d)==-1,"local density asymmetry rejected");
  values[1]=0.4; values[2]=NAN;
  require(ValidateSymmetryTerms(&d)==-1,"NaN rejected before cancellation");
  values[2]=DBL_MAX; values[3]=DBL_MAX;
  require(ValidateSymmetryTerms(&d)==-1,"nonfinite coefficient sum rejected");
  values[2]=0.3; values[3]=-0.3; d.NSymTrans=1;
  d.EDGeneralTransfer=ptrs+2; d.EDParaGeneralTransfer=values+2; d.EDNTransfer=1;
  require(ValidateSymmetryTerms(&d)==0,"GC transverse term allowed");
  d.iCalcModel=Kondo; d.NCond=1;
  require(NormalizeSymmetryKondoQuantumNumbers(&d,1,0,0,0)==0,"normalize NConserved");
  require(ValidateSymmetryTerms(&d)==0,"NConserved transverse term allowed");
  d.iCalcModel=Kondo; d.Total2Sz=1;
  require(NormalizeSymmetryKondoQuantumNumbers(&d,1,1,0,0)==0,"normalize canonical");
  require(ValidateSymmetryTerms(&d)==-1,"canonical transverse term rejected");
  d.EDNTransfer=2;
  require(ValidateSymmetryTerms(&d)==0,"cancelled transverse terms conserve Sz");
  require(CanonicalizeSymmetryKondoTerm(&d,&t,NULL,reject_monomial,NULL)==-1,
          "callback errors propagated");
  a.def=&d; a.input=17;
  require(CanonicalizeSymmetryKondoTerm(&d,&t,swap,canonical_matrix,&a)==0,
          "permutation before normalization");
  close_value(a.v[17],0,"mapped local flip does not act on up");
  a.input=33; memset(a.v,0,sizeof(a.v));
  require(CanonicalizeSymmetryKondoTerm(&d,&t,swap,canonical_matrix,&a)==0,
          "mapped local flip acts on down");
  close_value(a.v[17],1,"mapped local flip amplitude");
  t.value=NAN;
  require(CanonicalizeSymmetryKondoTerm(&d,&t,NULL,canonical_matrix,&a)==-1,
          "direct canonicalizer rejects NaN");
  t.value=INFINITY;
  require(CanonicalizeSymmetryKondoTerm(&d,&t,NULL,canonical_matrix,&a)==-1,
          "direct canonicalizer rejects infinity");
  t.value=1;
  t.index[1]=2;
  require(CanonicalizeSymmetryKondoTerm(&d,&t,NULL,canonical_matrix,&a)==-1,
          "canonicalizer rejects invalid spin");
  t.index[1]=0; swap[0]=1;
  require(CanonicalizeSymmetryKondoTerm(&d,&t,swap,canonical_matrix,&a)==-1,
          "canonicalizer rejects type-changing permutation");
}
struct HamiltonianMatrix {
  const struct DefineList *def;
  double complex h[64][64];
};
static int accumulate(const struct SymmetryTerm *t, void *context)
{
  struct HamiltonianMatrix *h=context;
  unsigned int in,dim=1U<<(2*h->def->Nsite);
  check_term(h->def,t);
  for(in=0;in<dim;++in) if(physical(h->def,in)) {
    unsigned long out=0;
    double complex value=0;
    int status=ApplySymmetryTerm(h->def,t,in,&out,&value);
    require(status>=0,"family application");
    if(status) h->h[out][in]+=value;
  }
  return 0;
}
static void compare_family(const struct DefineList *d, double complex expected[64][64],
                           const char *label)
{
  struct HamiltonianMatrix h={0};
  unsigned int in,out,dim=1U<<(2*d->Nsite);
  h.def=d;
  require(EnumerateSymmetryTerms(d,-1,accumulate,&h)==0,"enumerate family");
  for(in=0;in<dim;++in) if(physical(d,in)) {
    double diagonal=999;
    require(EvaluateSymmetryStateDiagonal(d,in,&diagonal)==0,"family diagonal");
    close_value(diagonal,expected[in][in],label);
    for(out=0;out<dim;++out) close_value(h.h[out][in],expected[out][in],label);
  }
  require(ValidateSymmetryTerms(d)==0,"family invariant in identity group");
}
static void test_families(void)
{
  int local[]={LOCSPIN,ITINERANT,ITINERANT};
  int pair[]={0,1},*pairs[]={pair,pair};
  double coupling[]={0.37,0.37};
  unsigned int in;
  double complex expected[64][64];
  struct DefineList base=definition(KondoGC,local,2),d;
  int identity[]={0,1},*group[]={identity};
  base.NSymTrans=1; base.SymTrans=group;
  /* Individual U rows: local doublon is zero; conduction doublon is one. */
  for(pair[0]=0;pair[0]<2;++pair[0]) {
    d=base; d.NCoulombIntra=1; d.CoulombIntra=pairs; d.ParaCoulombIntra=coupling;
    memset(expected,0,sizeof(expected));
    for(in=0;in<16;++in) if(physical(&d,in) && pair[0]==1 && in/4==3)
      expected[in][in]=0.37;
    compare_family(&d,expected,"CoulombIntra physical doublon");
  }
  pair[0]=0;
  d=base; d.NCoulombInter=1; d.CoulombInter=pairs; d.ParaCoulombInter=coupling;
  memset(expected,0,sizeof(expected));
  for(in=0;in<16;++in) if(physical(&d,in))
    expected[in][in]=0.37*((in/4%2)+(in/8%2));
  compare_family(&d,expected,"CoulombInter local density is one");
  d=base; d.NHundCoupling=1; d.HundCoupling=pairs; d.ParaHundCoupling=coupling;
  memset(expected,0,sizeof(expected));
  for(in=0;in<16;++in) if(physical(&d,in))
    expected[in][in]=-0.37*((in%4==1 ? in/4 : in/8)%2);
  compare_family(&d,expected,"Hund negative parallel density");
  d=base; d.NExchangeCoupling=1; d.ExchangeCoupling=pairs; d.ParaExchangeCoupling=coupling;
  memset(expected,0,sizeof(expected));
  expected[6][9]=expected[9][6]=-0.37;
  compare_family(&d,expected,"Exchange negative physical spin flip and late projection");
  d.NExchangeCoupling=2;
  expected[6][9]=expected[9][6]=-0.74;
  compare_family(&d,expected,"duplicate Exchange accumulates");
  d.NExchangeCoupling=1; pair[0]=1; pair[1]=0;
  expected[6][9]=expected[9][6]=-0.37;
  compare_family(&d,expected,"reverse Exchange same physical sign");
  pair[0]=pair[1]=0;
  memset(expected,0,sizeof(expected));
  compare_family(&d,expected,"same-site Exchange zero");
  pair[1]=1;
  {
    double v=-0.25,h=-0.5;
    d=base; d.NIsingCoupling=1;
    d.NCoulombInter=1; d.CoulombInter=pairs; d.ParaCoulombInter=&v;
    d.NHundCoupling=1; d.HundCoupling=pairs; d.ParaHundCoupling=&h;
    memset(expected,0,sizeof(expected));
    for(in=0;in<16;++in) if(physical(&d,in))
      expected[in][in]=0.25*(in%4==1 ? 1 : -1)*((int)(in/4%2)-(int)(in/8%2));
    compare_family(&d,expected,"Ising reader expansion counted exactly once");
  }
  {
    int chemi=1,spin=1,transfer[]={1,0,1,1},*transfers[]={transfer};
    double mu=0.23;
    double complex t=0.19+0.13*I;
    d=base; d.EDNChemi=1; d.EDChemi=&chemi; d.EDSpinChemi=&spin; d.EDParaChemi=&mu;
    d.EDNTransfer=1; d.EDGeneralTransfer=transfers; d.EDParaGeneralTransfer=&t;
    memset(expected,0,sizeof(expected));
    for(in=0;in<16;++in) if(physical(&d,in)) {
      expected[in][in]=-0.23*(in/8%2);
      if(in/4==2) expected[in-4][in]=-t;
    }
    compare_family(&d,expected,"Transfer and EDChemi negative coefficients");
  }
  {
    int diag[]={0,0,1,1},off[]={0,0,0,1,1,1,1,0};
    int *diags[]={diag},*offs[]={off};
    double dv=0.37;
    double complex ov=0.2+0.1*I;
    d=base; d.NInterAll=17; /* Original storage is deliberately absent. */
    d.NInterAll_Diagonal=1; d.InterAll_Diagonal=diags; d.ParaInterAll_Diagonal=&dv;
    d.NInterAll_OffDiagonal=1; d.InterAll_OffDiagonal=offs; d.ParaInterAll_OffDiagonal=&ov;
    memset(expected,0,sizeof(expected));
    for(in=0;in<16;++in) if(physical(&d,in) && in%4==1 && in/8%2) expected[in][in]=dv;
    expected[9][6]=ov;
    compare_family(&d,expected,"InterAll physical onsite spin product and split storage");
  }
  {
    int identity3[]={0,1,2},*group3[]={identity3};
    d=definition(KondoGC,local,3); d.NSymTrans=1; d.SymTrans=group3;
    pair[0]=1;pair[1]=2;
    d.NPairHopping=1; d.PairHopping=pairs; d.ParaPairHopping=coupling;
    memset(expected,0,sizeof(expected));
    expected[13][49]=expected[14][50]=0.37;
    compare_family(&d,expected,"conduction PairHop positive amplitude");
    pair[0]=2;pair[1]=1;
    memset(expected,0,sizeof(expected));
    expected[49][13]=expected[50][14]=0.37;
    compare_family(&d,expected,"reverse conduction PairHop");
    pair[0]=pair[1]=1;
    memset(expected,0,sizeof(expected));
    for(in=0;in<64;++in) if(physical(&d,in) && (in/4)%4==3) expected[in][in]=0.37;
    compare_family(&d,expected,"onsite conduction PairHop doublon");
  }
}
static void test_mixed_limits_and_permutation(void)
{
  int layouts[4][3]={{LOCSPIN,LOCSPIN,ITINERANT},
                    {LOCSPIN,ITINERANT,LOCSPIN},
                    {ITINERANT,LOCSPIN,ITINERANT},
                    {ITINERANT,ITINERANT,ITINERANT}};
  unsigned int layout,a,b,c,e;
  for(layout=0;layout<4;++layout) {
    struct DefineList d=definition(KondoGC,layouts[layout],3);
    struct SymmetryTerm t={2,{0},0.31-0.17*I};
    for(a=0;a<6;++a) for(b=0;b<6;++b) for(c=0;c<6;++c) for(e=0;e<6;++e) {
      int ix[]={a/2,a%2,b/2,b%2,c/2,c%2,e/2,e%2};
      memcpy(t.index,ix,sizeof(ix)); check_term(&d,&t);
    }
  }
  {
    int local[]={LOCSPIN,LOCSPIN};
    struct DefineList d=definition(KondoGC,local,2);
    struct SymmetryTerm t={2,{0,0,0,1,1,1,1,0},1};
    check_term(&d,&t); /* Pure spin limit. */
  }
  {
    int local[]={LOCSPIN,ITINERANT,LOCSPIN},swap[]={2,1,0};
    struct DefineList d=definition(KondoGC,local,3);
    struct SymmetryTerm t={2,{0,0,2,0,2,1,0,1},0.3+0.2*I},mapped=t;
    unsigned int i,in,out;
    for(i=0;i<8;i+=2) mapped.index[i]=swap[t.index[i]];
    for(in=0;in<64;++in) if(physical(&d,in)) {
      struct CanonicalMatrix a={0};
      double complex expected[64];
      a.def=&d;a.input=in;
      oracle(&d,&mapped,in,expected);
      require(CanonicalizeSymmetryKondoTerm(&d,&t,swap,canonical_matrix,&a)==0,"mapped canonical term");
      for(out=0;out<64;++out) close_value(a.v[out],expected[out],"mapped crossed Exchange CAR");
    }
  }
}

int main(void)
{
  stdoutMPI=stdout;
  test_empty_and_application();
  test_validation();
  test_families();
  test_mixed_limits_and_permutation();
  puts("unittest_symmetry_kondo_terms: PASS");
  return 0;
}
