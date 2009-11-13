// -*- mode:C++ ; compile-command: "g++-3.4 -I.. -I../include -g -c gen.cc -Wall" -*-
#include "first.h"
/*
 *  Copyright (C) 2001,7 B. Parisse, Institut Fourier, 38402 St Martin d'Heres
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
using namespace std;
#define __USE_ISOC9X 1
#include <time.h>
#include <stdexcept>
#include <ctype.h>
#include <math.h>
#include <cstdlib>
#include <list>
#include <errno.h>
#include <string.h>
#include <iomanip>
#include "gen.h"
#include "gausspol.h"
#include "identificateur.h"
#include "poly.h"
#include "usual.h"
#include "input_lexer.h"
#include "sym2poly.h"
#include "vecteur.h"
#include "modpoly.h"
#include "alg_ext.h"
#include "prog.h"
#include "rpn.h"
#include "plot.h"
#include "intg.h"
#include "subst.h"
#include "derive.h"
#include "threaded.h"
#include "maple.h"
#include "solve.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  unsigned control_c_counter=0;
  unsigned control_c_counter_mask=0xff;

#ifdef HAVE_LIBPTHREAD
  pthread_mutex_t mpfr_mutex = PTHREAD_MUTEX_INITIALIZER;
  pthread_mutex_t locale_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif

  // bool is_inevalf=false;

  // FIXME GIAC_CONTEXT
  string last_evaled_function(GIAC_CONTEXT){
    std::vector<const std::string *> & last =last_evaled_function_name(contextptr);
    if (last.empty() || !last.back())
      return "";
    return *last.back()+" ";
  }

  void settypeerr(GIAC_CONTEXT){
    throw(std::runtime_error(last_evaled_function(contextptr)+"Bad Argument Type"));
  }

  void setsizeerr(GIAC_CONTEXT){
    throw(std::runtime_error(last_evaled_function(contextptr)+"Bad Argument Value"));
  }

  void setdimerr(GIAC_CONTEXT){
    throw(std::runtime_error(last_evaled_function(contextptr)+"Invalid dimension"));
  }

  void settypeerr(const string & s){
    throw(std::runtime_error(s+" Error: Bad Argument Type"));
  }

  void setsizeerr(const string & s){
    throw(std::runtime_error(s+" Error: Bad Argument Value"));
  }

  void setdimerr(const string & s){
    throw(std::runtime_error(s+" Error: Invalid dimension"));
  }

  void divisionby0err(const gen & e,GIAC_CONTEXT){
    throw(std::runtime_error(last_evaled_function(contextptr)+"Division of " + e.print(contextptr)+ " by 0"));
  }

  void cksignerr(const gen & e,GIAC_CONTEXT){
    throw(std::runtime_error(last_evaled_function(contextptr)+"Unable to check sign: "+e.print(contextptr)));
  }

  void invalidserieserr(const string & s,GIAC_CONTEXT){
    throw(std::runtime_error(last_evaled_function(contextptr)+"Invalid series expansion: "+s));
  }

  void toofewargs(const string & s,GIAC_CONTEXT){
    throw(std::runtime_error(last_evaled_function(contextptr)+"Too few arguments: "+s));
  }

  void toomanyargs(const string & s,GIAC_CONTEXT){
    throw(std::runtime_error(last_evaled_function(contextptr)+"Too many arguments: "+s));
  }

  void maxordererr(GIAC_CONTEXT){
    throw(std::runtime_error(last_evaled_function(contextptr)+"Max order ("+gen(max_series_expansion_order).print(contextptr)+") exceeded or non unidirectional series"));
  }

  void setstabilityerr(GIAC_CONTEXT){
    throw(std::runtime_error(last_evaled_function(contextptr)+"calculation size limit exceeded"));
  }

  void swap(int & a,int & b){
    int t=a;
    a=b;
    b=t;
  }

  gen vector2vecteur(const vecteur & v){
    gen g=v.back()-v.front();
    if (g.type!=_VECT)
      return makenewvecteur(re(g,context0),im(g,context0));
    return g;
  }

  // void parseerror(){
  //  throw(std::runtime_error("Parse error"));
  // }

  gen zero(0);
  gen minus_one(-1);
  gen plus_one(1);
  gen plus_two(2);
  gen plus_three(3);
  gen cst_i(0,1);
  double rad2deg_d(180/M_PI);
  double deg2rad_d(M_PI/180);
  gen rad2deg_g(rad2deg_d);
  gen deg2rad_g(deg2rad_d);
  
  bool is_inf(const gen & e){
    switch (e.type){
    case _IDNT:
      return *e._IDNTptr==_IDNT_infinity();
    case _SYMB:
      return is_inf(e._SYMBptr->feuille);
    default:
      return false;
    }
  }
  bool is_undef(const gen & e){
    if (e.type==_IDNT)
      return *e._IDNTptr==_IDNT_undef();
    if ( (e.type==_POLY) && (!e._POLYptr->coord.empty()) )
      return is_undef(e._POLYptr->coord.front().value);
    else
      return false;
  }
  
  const vecteur null_vector;
  
  enum { debugtype=_CPLX };
#define debugtypeptr _CPLXptr

  int absint(int a){
    if (a<0)
      return -a;
    else
      return a;
  }

  int min(int a, int b){
    if (a<b)
      return a;
    else
      return b;
  }

  int max(int a, int b){
    if (a<b)
      return b;
    else
      return a;
  }

  int invmod(int a,int b){
    if (a==1 || a==-1 || a==1-b)
      return a;
    int aa(1),ab(0),ar(0);
    div_t qr;
    while (b){
      qr=div(a,b);
      ar=aa-qr.quot*ab;
      a=b;
      b=qr.rem;
      aa=ab;
      ab=ar;
    }
    if (a==1)
      return aa;
    if (a!=-1)
      setsizeerr("Not invertible");
    return -aa;
  }

  unsigned invmod(unsigned a,int b){
    int i=invmod(int(a),int(b));
    if (i<0)
      i+=b;
    return unsigned(i);
  }

  int invmod(longlong a,int b){
    return invmod(int(a%b),b);
  }

  int powmod(int a,unsigned long n,int m){
    if (!n)
      return 1;
    int b=powmod(a,n/2,m);
    longlong tmp=b;
    b=(tmp*b) % m;
    tmp=b;
    if (n % 2)
      return (tmp*a) % m;
    else
      return b;
  }

  int smod(int r,int m){
    if (!m)
      return r;
    if (m<0)
      m=-m;
    r = r % m;
    register int tmp= r+r;
    if (tmp>m)
      return r-m;
    if (tmp<=-m)
      return r+m;
    return r;
  }

  int gcd(int a,int b){
    int r;
    while (b){
      r=a%b;
      a=b;
      b=r;
    }
    return absint(a);
  }

  int simplify(int & a,int & b){
    int d=gcd(a,b);
    a=a/d;
    b=b/d;
    return d;
  }

  // list for allocated mpz_t
  list<__mpz_struct> zint_list;

  inline mpz_t * new_mpz_t(){
    __mpz_struct tmp={0,0,0};
    zint_list.push_back(tmp);
    return (mpz_t *) &zint_list.back();
  }

  gen factorial(unsigned long int i){
    if (i>FACTORIAL_SIZE_LIMIT)
      setstabilityerr();
    mpz_t * e = (mpz_t *) malloc(sizeof(mpz_t)); // new mpz_t[1];
    //mpz_t * e=new_mpz_t();
    mpz_init(*e);
    mpz_fac_ui(*e,i);
    return e;
  }

  gen comb(unsigned long int i,unsigned long j){
    if (i>FACTORIAL_SIZE_LIMIT || j>FACTORIAL_SIZE_LIMIT)
      setstabilityerr();
    mpz_t * e =(mpz_t *) malloc(sizeof(mpz_t));
    mpz_init(*e);
    if (i<j)
      return e;
    mpz_set_ui(*e,1);
    for (unsigned long int k=i;k>i-j;--k)
      mpz_mul_ui(*e,*e,k);
    mpz_t tmp;
    mpz_init(tmp);
    mpz_fac_ui(tmp,j);
    mpz_fdiv_q(*e,*e,tmp);
    return e;
  }

  gen perm(unsigned long int i,unsigned long j){
    if (i>FACTORIAL_SIZE_LIMIT || j>FACTORIAL_SIZE_LIMIT)
      setstabilityerr();
    mpz_t * e =(mpz_t *) malloc(sizeof(mpz_t));//new mpz_t[1];
    mpz_init(*e);
    if (i<j)
      return e;
    mpz_set_ui(*e,1);
    for (unsigned long int k=i;k>i-j;--k)
      mpz_mul_ui(*e,*e,k);
    return e;
  }

  gen gen::change_subtype(int newsubtype){ 
    subtype=newsubtype; return *this; 
  }

  gen::gen(longlong i) { 
#ifdef COMPILE_FOR_STABILITY
    control_c();
#endif
    val=i;
    //    longlong temp=val;
    if (val==i && val!=1<<31){
      type=_INT_;
      subtype=0;
    }
    else {
      type =_ZINT;
      subtype=0;
      ptr_val.ref_count=new int(1);
      _ZINTptr = (mpz_t *) malloc(sizeof(mpz_t)); // new mpz_t[1];
      // convert longlong to mpz_t
      bool signe=(i<0);
      if (signe)
	i=-i;
      unsigned int i1=i>>32;
      unsigned int i2=i;
      mpz_init_set_ui(*_ZINTptr,i1);
      mpz_mul_2exp(*_ZINTptr,*_ZINTptr,32);
      mpz_add_ui(*_ZINTptr,*_ZINTptr,i2);
      if (signe)
	mpz_neg(*_ZINTptr,*_ZINTptr);
      /*
      longlong lbase=65536;
      long base=65536;
      longlong i1=i/lbase;
      long i2=i1/lbase; // i2=i/2^32
      //cout << "Initialization of " << _ZINTptr << endl ;
      mpz_init_set_si(*_ZINTptr,i2);
      mpz_mul_ui(*_ZINTptr,*_ZINTptr,base); // i/2^32 * 2^16
      long i2mod=i1 % lbase;
      if (i2mod>0)
	mpz_add_ui(*_ZINTptr,*_ZINTptr,i2mod);
      else
	mpz_sub_ui(*_ZINTptr,*_ZINTptr,-i2mod); // i/2^16
      mpz_mul_ui(*_ZINTptr,*_ZINTptr,base); // i/2^16 * 2^16
      long i1mod = i % lbase;
      if (i1mod>0)
	mpz_add_ui(*_ZINTptr,*_ZINTptr,i1mod);
      else
	mpz_sub_ui(*_ZINTptr,*_ZINTptr,-i1mod); // i
      */
    }
  }

  gen::gen(const mpz_t & m) { 
#ifdef COMPILE_FOR_STABILITY
    control_c();
#endif
    type =_ZINT;
    subtype=0;
    ptr_val.ref_count=new int(1);
    _ZINTptr= (mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1] ;
    mpz_init_set(*_ZINTptr,m);
  }

#ifdef HAVE_GMPXX_H
  gen::gen(const mpz_class & m){
    int l=mpz_sizeinbase(m.get_mpz_t(),2);
    subtype=0;
    if (l<32){
      type = _INT_;
      val = mpz_get_si(m.get_mpz_t());
    }
    else {
      type =_ZINT;
      ptr_val.ref_count=new int(1);
      _ZINTptr= (mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1] ;
      mpz_init_set(*_ZINTptr,m.get_mpz_t());
    }
  }
#endif

  gen::gen(const identificateur & s){
#ifdef COMPILE_FOR_STABILITY
    control_c();
#endif
    type=_IDNT;
    subtype=0;
    ptr_val.ref_count = new int(1);
    _IDNTptr= new identificateur(s);
  }

  gen::gen(const vecteur & v,short int s){
    type=_VECT;
    subtype=s;
    ptr_val.ref_count = new int(1);
    _VECTptr= new vecteur(v);
  }

  gen::gen(vecteur * vptr,short int s){
    type=_VECT;
    subtype=s;
    ptr_val.ref_count = new int(1);
    _VECTptr= vptr;
  }

  gen::gen(const symbolic & s){
#ifdef COMPILE_FOR_STABILITY
    control_c();
#endif
    type = _SYMB;
    subtype = 0;
    ptr_val.ref_count = new int(1);
    _SYMBptr = new symbolic(s) ;
  }

  gen::gen(symbolic * sptr){
#ifdef COMPILE_FOR_STABILITY
    control_c();
#endif
    type = _SYMB;
    subtype = 0;
    _SYMBptr = sptr;
    ptr_val.ref_count = new int(1);
  }

  gen::gen(const gen_user & g){
    type = _USER;
    subtype=0;
    ptr_val.ref_count = new int(1);
    _USERptr = g.memory_alloc() ;
  }

  gen::gen(const eqwdata & g){
    type = _EQW;
    subtype=0;
    ptr_val.ref_count = new int(1);
    _EQWptr = new eqwdata(g);
  }

  gen::gen(const grob & g){
    type = _GROB;
    subtype=0;
    ptr_val.ref_count = new int(1);
    _GROBptr = new grob(g);
  }

  gen makemap(){ 
    gen g;
    g.type=_MAP;
    g.subtype=0;
    g.ptr_val.ref_count = new int(1);
    g._MAPptr = new gen_map(ptr_fun(islesscomplexthanf));
    return g;
  }

  gen::gen(const gen_map & s){
    type = _MAP;
    subtype = 0;
    ptr_val.ref_count = new int(1);
    _MAPptr = new gen_map(s) ;
  }

  real_object::real_object() { 
#ifdef HAVE_LIBMPFR
    mpfr_init(inf); 
#else
    mpf_init(inf); 
#endif
  }

  real_object::real_object(const real_object & g){ 
#ifdef HAVE_LIBMPFR
    mpfr_init2(inf,mpfr_get_prec(g.inf));
    mpfr_set(inf,g.inf,GMP_RNDN);
#else
    mpf_init_set(inf,g.inf);
#endif
  }

  real_object::real_object(double d) { 
#ifdef HAVE_LIBMPFR
    mpfr_init_set_d(inf,d,GMP_RNDN); 
#else
    mpf_init_set_d(inf,d); 
#endif
  }

#ifdef HAVE_LIBMPFR
  real_object::real_object(const mpfr_t & d) { 
    mpfr_init2(inf,mpfr_get_prec(d));
    mpfr_set(inf,d,GMP_RNDN);
  }
#endif

  real_object::real_object(const mpf_t & d) { 
#ifdef HAVE_LIBMPFR
    mpfr_init(inf);
    mpfr_set_f(inf,d,GMP_RNDN);
#else
    mpf_init_set(inf,d); 
#endif
  }

  real_object::real_object(const gen & g){
    switch (g.type){
    case _INT_:
#ifdef HAVE_LIBMPFR
      mpfr_init_set_si(inf,g.val,GMP_RNDN);
#else
      mpf_init_set_si(inf,g.val);
#endif
      return;
    case _DOUBLE_:
#ifdef HAVE_LIBMPFR
      mpfr_init_set_d(inf,g._DOUBLE_val,GMP_RNDN);
#else
      mpf_init_set_d(inf,g._DOUBLE_val);
#endif
      return;
    case _ZINT:
#ifdef HAVE_LIBMPFR
      mpfr_init(inf);
      mpfr_set_z(inf,*g._ZINTptr,GMP_RNDN);
#else
      mpf_init(inf);
      mpf_set_z(inf,*g._ZINTptr);
#endif
      return;
    case _REAL:
#ifdef HAVE_LIBMPFR
      mpfr_init2(inf,mpfr_get_prec(g._REALptr->inf));
      mpfr_set(inf,g._REALptr->inf,GMP_RNDN);
#else
      mpf_init_set(inf,g._REALptr->inf);
#endif
      return;
    }
    if (g.type==_FRAC){
      gen tmp=real_object(g._FRACptr->num)/real_object(g._FRACptr->den);
      if (tmp.type==_REAL){
#ifdef HAVE_LIBMPFR
	mpfr_init2(inf,mpfr_get_prec(tmp._REALptr->inf));
	mpfr_set(inf,tmp._REALptr->inf,GMP_RNDN);
#else
	mpf_init_set(inf,tmp._REALptr->inf);
#endif
	return;
      }
    }
    setsizeerr("Unable to convert to real "+g.print(context0));
  }

  real_object::real_object(const gen & g,unsigned int precision){
    switch (g.type){
    case _INT_:
#ifdef HAVE_LIBMPFR
      mpfr_init2(inf,precision);
      mpfr_set_si(inf,g.val,GMP_RNDN);
#else
      mpf_init_set_si(inf,g.val);
#endif
      return;
    case _DOUBLE_:
#ifdef HAVE_LIBMPFR
      mpfr_init2(inf,precision);
      mpfr_set_d(inf,g._DOUBLE_val,GMP_RNDN);
#else
      mpf_init_set_d(inf,g._DOUBLE_val);
#endif
      return;
    case _ZINT:
#ifdef HAVE_LIBMPFR
      mpfr_init2(inf,precision);
      mpfr_set_z(inf,*g._ZINTptr,GMP_RNDN);
#else
      mpf_init(inf);
      mpf_set_z(inf,*g._ZINTptr);
#endif
      return;
    case _REAL:
#ifdef HAVE_LIBMPFR
      mpfr_init2(inf,precision);
      mpfr_set(inf,g._REALptr->inf,GMP_RNDN);
#else
      mpf_init_set(inf,g._REALptr->inf);
#endif
      return;
    }
    if (g.type==_FRAC){
      gen tmp=real_object(g._FRACptr->num,precision)/real_object(g._FRACptr->den,precision);
      if (tmp.type==_REAL){
#ifdef HAVE_LIBMPFR
	mpfr_init2(inf,mpfr_get_prec(tmp._REALptr->inf));
	mpfr_set(inf,tmp._REALptr->inf,GMP_RNDN);
#else
	mpf_init_set(inf,tmp._REALptr->inf);
#endif
	return;
      }
    }
    setsizeerr("Unable to convert to real "+g.print(context0));
  }

  gen set_precision(const gen & g,int nbits){
#ifdef HAVE_LIBMPFR
    return real_object(g,nbits);
#else
    return real_object(g);
#endif
  }

  gen accurate_evalf(const gen & g,int nbits){
    gen r(g.re(context0)),i(g.im(context0)); // only called for numeric values
    if (is_zero(i))
      return set_precision(r,nbits);
    else
      return gen(set_precision(r,nbits),set_precision(i,nbits));
  }

  vecteur accurate_evalf(const vecteur & v,int nbits){
    vecteur res(v);
    iterateur it=res.begin(),itend=res.end();
    for (;it!=itend;++it)
      *it = accurate_evalf(*it,nbits);
    return res;
  }

  gen::gen(const real_object & g){
    type = _REAL;
    subtype=0;
    ptr_val.ref_count = new int(1);
    _REALptr = new real_object;
#ifdef HAVE_LIBMPFR
    mpfr_set_prec(_REALptr->inf,mpfr_get_prec(g.inf));
    mpfr_set(_REALptr->inf,g.inf,GMP_RNDN);
#else
    mpf_set(_REALptr->inf,g.inf);
#endif
  }

  gen::gen(const polynome & p){
    subtype=0;
    if (p.coord.empty()){
      type = _INT_;
      val = 0;
    }
    else {
      if (Tis_constant<gen>(p) && is_atomic(p.coord.front().value) ){
	type = _INT_;
	* this = p.coord.front().value;
      }
      else {
	type = _POLY;
	ptr_val.ref_count = new int(1);
	_POLYptr = new polynome(p) ;
      }
    }
    /*
    type = poly;
    polyptr = new polynome(p) ;
    ptr_val.ref_count = new int(1);
    */
  }

  gen::gen(const real_complex_rootof & r){
    subtype=0;
    type = _ROOT;
    ptr_val.ref_count = new int(1);
    _ROOTptr = new real_complex_rootof(r) ;
  }

  gen::gen(const fraction & p){
    subtype=0;
    if (is_exactly_zero(p.num)){
      type=_INT_;
      val=0;
    }
    else {
      if (is_one(p.den)){
	type=_INT_;
	*this = p.num;
      }
      else {
          if (is_minus_one(p.den)){
              type=_INT_;
              *this = -p.num;
          }
          else {
              
              type = _FRAC;
              ptr_val.ref_count = new int(1);
              _FRACptr = new fraction(p) ;
          }
          
      }
    }
  }

  gen::gen(polynome * pptr){
    subtype=0;
    type = _POLY;
    _POLYptr = pptr ;
    ptr_val.ref_count = new int(1);
  }

  // WARNING coerce *mptr to an int if possible, in this case delete mptr
  // Pls do not use this constructor unless you know exactly what you do!!
  gen::gen(mpz_t * mptr){
    subtype=0;
    int l=mpz_sizeinbase(*mptr,2);
    // if (l<17){
    if (l<32){
      type = _INT_;
      val = mpz_get_si(*mptr);
      // cout << "Destruction by mpz_t * " << *mptr << endl;
      mpz_clear(*mptr);
      free(mptr);
    }
    else {
      type =_ZINT;
      _ZINTptr = mptr;
      ptr_val.ref_count=new int(1);
    }
    // cout << *this << endl;
  }

  gen::gen(const my_mpz& z){
    subtype=0;
    int l=mpz_sizeinbase(z.ptr,2);
    if (l<32){
      type = _INT_;
      val = mpz_get_si(z.ptr);
    }
    else {
      type =_ZINT;
      ptr_val.ref_count=new int(1);
      _ZINTptr = (mpz_t *) malloc(sizeof(mpz_t));
      mpz_init_set(*_ZINTptr,z.ptr);
    }    
  }

  gen::gen(const gen & e) { 
    // if (e.type==debugtype)
    // cout << e << "Construct @ " << *e.debugtypeptr << "[" << *e.ptr_val.ref_count <<"+1]" <<endl;
    type=e.type;
    subtype=e.subtype;
    if (!type){
      val=e.val;
      return;
    }
    switch (type ) {
    case _DOUBLE_:
      _DOUBLE_val=e._DOUBLE_val;
      break;
    case _POINTER_:
      ptr_val.ref_count=e.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _POINTER_val=e._POINTER_val;
      break;
    case _ZINT: 
      ptr_val.ref_count=e.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _ZINTptr=e._ZINTptr;
      break; 
    case _REAL: 
      ptr_val.ref_count=e.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _REALptr=e._REALptr;
      break; 
    case _CPLX: 
      ptr_val.ref_count=e.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _CPLXptr=e._CPLXptr;
      break; 
    case _IDNT:
      ptr_val.ref_count=e.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _IDNTptr=e._IDNTptr;
      break;
    case _POLY: 
      ptr_val.ref_count=e.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _POLYptr=e._POLYptr;
      break; 
    case _FRAC: 
      ptr_val.ref_count=e.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _FRACptr=e._FRACptr;
      break; 
    case _VECT: 
      ptr_val.ref_count=e.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _VECTptr=e._VECTptr;
      break; 
    case _SPOL1: 
      ptr_val.ref_count=e.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _SPOL1ptr=e._SPOL1ptr;
      break; 
    case _SYMB:
      ptr_val.ref_count=e.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _SYMBptr=e._SYMBptr;
      break; 
    case _USER:
      ptr_val.ref_count=e.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _USERptr=e._USERptr;
      break; 
    case _EXT:
      ptr_val.ref_count=e.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _EXTptr=e._EXTptr;
      break; 
    case _MOD:
      ptr_val.ref_count=e.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _MODptr=e._MODptr;
      break; 
    case _ROOT:
      ptr_val.ref_count=e.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _ROOTptr=e._ROOTptr;
      break; 
    case _STRNG:
      ptr_val.ref_count=e.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _STRNGptr=e._STRNGptr;
      break; 
    case _FUNC:
      ptr_val.ref_count=e.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _FUNCptr=e._FUNCptr;
      break; 
    case _MAP:
      ptr_val.ref_count=e.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _MAPptr=e._MAPptr;
      break; 
    case _EQW:
      ptr_val.ref_count=e.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _EQWptr=e._EQWptr;
      break; 
    default: 
      settypeerr("Gen constructor "+e.print(context0));
    }
  }

  gen::gen(int a,int b) {
    subtype=0;
    if (!b){
      type=_INT_;
      val=a;
    }
    else {
      type =_CPLX;
      subtype=0;
      _CPLXptr = new gen[2];
      _CPLXptr->type=_INT_;
      _CPLXptr->val=a;
      ++_CPLXptr;
      _CPLXptr->type=_INT_;
      _CPLXptr->val =b;
      --_CPLXptr;
      ptr_val.ref_count = new int(1);
    }
  }

  gen::gen(double a,double b){
    subtype=0;
    // cout << a << " " << b << " " << epsilon << endl;
    if (fabs(b)<1e-12*fabs(a)){ 
      type=_DOUBLE_;
      _DOUBLE_val=a;
    }
    else {
      type =_CPLX;
      subtype=3;
      _CPLXptr = new gen[2];
      _CPLXptr->type=_DOUBLE_;
      _CPLXptr->_DOUBLE_val=a;
      ++_CPLXptr;
      _CPLXptr->type=_DOUBLE_;
      _CPLXptr->_DOUBLE_val =b;
      --_CPLXptr;
      ptr_val.ref_count = new int(1);
    }
  }
  
  gen::gen(const gen & a,const gen & b) { // a and b must be type <2!
    subtype=0;
    if ( (a.type>=_CPLX) || (b.type>=_CPLX) )
      settypeerr("complex constructor");
    if (is_zero(b)){
      type=a.type;
      switch (type ) {
      case _INT_: 
	val=a.val;
	break; 
      case _DOUBLE_: 
	_DOUBLE_val=a._DOUBLE_val;
	break; 
      case _ZINT: 
	ptr_val.ref_count=a.ptr_val.ref_count;
	(*ptr_val.ref_count)++;
	_ZINTptr=a._ZINTptr; // a is a _ZINT
	break; 
      case _REAL: 
	ptr_val.ref_count=a.ptr_val.ref_count;
	(*ptr_val.ref_count)++;
	_REALptr=a._REALptr; // a is a _ZINT
	break; 
      default: 
	settypeerr("complex constructor");
      }
    }
    else {
      type =_CPLX;
      subtype= (a.type==_DOUBLE_) + (b.type==_DOUBLE_)*2;
      _CPLXptr = new gen[2];
      * _CPLXptr = a;
      ++_CPLXptr;
      * _CPLXptr = b;
      --_CPLXptr;
      ptr_val.ref_count = new int(1);
    }
  }
  
  gen::gen(const complex<double> & c) {
    subtype=0;
    type=_CPLX;
    subtype=3;
    _CPLXptr = new gen[2];
    * _CPLXptr = real(c);
    ++_CPLXptr;
    * _CPLXptr = imag(c);
    --_CPLXptr;
    ptr_val.ref_count = new int(1);    
  }

  gen gen::makegen(int i) const {
    switch (type){
    case _INT_: case _ZINT: case _CPLX:
      return gen(i);
    case _VECT:
      return vecteur(1,i);
    case _USER:
      return _USERptr->makegen(i);
    default:
      setsizeerr("makegen of type "+print(context0));
      return 0;
    }
  }

  complex<double> gen2complex_d(const gen & e){
    if (e.type==_CPLX){
      if (e.subtype==3)
	return complex<double>((*e._CPLXptr)._DOUBLE_val,(*(e._CPLXptr+1))._DOUBLE_val);
      gen ee=e.evalf_double(1,context0); // ok
      if (ee.type==_DOUBLE_) return complex<double>(ee._DOUBLE_val,0);
      if (ee.type!=_CPLX) setsizeerr();
      return complex<double>((*ee._CPLXptr)._DOUBLE_val,(*(ee._CPLXptr+1))._DOUBLE_val);
    }
    if (e.type==_DOUBLE_)
      return complex<double>(e._DOUBLE_val,0);
    if (e.type==_INT_) 
      return complex<double>(e.val,0);
    if (e.type==_ZINT)
      return complex<double>(e.evalf(1,context0)._DOUBLE_val,0); // ok
    setsizeerr("gen -> complex");
    return 0;
  }

  gen chartab2gen(char * & s,GIAC_CONTEXT){
    gen res;
    // subtype=0;
    // initialize as a null _INT_
    // type = _INT_;
    // val = 0;
    if (!*s)
      return res;
    int base=0;
    if (s[0]=='#' || s[0]=='0') {
      if (s[1]=='x' || s[1]=='X'){
	s[0]='0';
	base=16;
      }
      if (s[1]=='o' || s[1]=='O'){
	s[0]='0';
	s[1]='0';
	base=8;
      }
    }
    if (s[1]=='b' || s[1]=='B'){
      s[1]='0';
      base=2;
    }
#ifdef _LIB_CE_ERRNO_H
	__set_errno(0);
#else
     errno = 0;
#endif
    char * endchar;
    long ll=strtol(s,&endchar,base);
    int l =strlen(s);
    if (*endchar) {// non integer
      if (decimal_digits(contextptr)>14){
#ifdef HAVE_LIBMPFR
	int nbits=digits2bits(decimal_digits(contextptr));
#ifdef HAVE_LIBPTHREAD
	int locked=pthread_mutex_trylock(&mpfr_mutex);
	if (!locked)
	  mpfr_set_default_prec(nbits);
	// mpf_set_default_prec (decimal_digits);
	real_object r;
	int res=mpfr_set_str(r.inf,s,10,GMP_RNDN);
	if (!locked)
	  pthread_mutex_unlock(&mpfr_mutex);
#else
	real_object r;
	mpfr_set_default_prec(nbits);
	int res=mpfr_set_str(r.inf,s,10,GMP_RNDN);
#endif // HAVE_LIBPTHREAD
#else // LIBMPFR
	real_object r;
	int res=mpf_set_str(r.inf,s,10);
#endif // LIBMPFR
	gen rg(r);
	// rg.dbgprint();
	if (!res)
	  return rg;
      }
      double d;
#ifdef HAVE_LIBPTHREAD
      int locked=pthread_mutex_trylock(&locale_mutex);
      if (!locked){
	char * lc=setlocale(LC_NUMERIC,0);
	setlocale(LC_NUMERIC,"POSIX");
	d=strtod(s,&endchar);
	setlocale(LC_NUMERIC,lc);
	pthread_mutex_unlock(&locale_mutex);
      }
      else
	d=strtod(s,&endchar);	
#else
      char * lc=setlocale(LC_NUMERIC,0);
      setlocale(LC_NUMERIC,"POSIX");
      d=strtod(s,&endchar);
      setlocale(LC_NUMERIC,lc);
#endif
      if (*endchar)
	return gen(string(s),contextptr);
      return gen(d);
    }
    if (!errno ){
      if (ll==int(ll) || 
	  ( (base == 2 ||base == 8 || base == 16) && 
	    ll == (unsigned int)(ll) )
	  )
	return gen ( int(ll));
      else
	return gen(longlong(ll));
    }
    // check if a non 0-9 char is there
    if (!base){
      base=10;
      for (int i=0;i<l;++i){
	if ((s[i]<'0') || (s[i]>'9'))
	  base=16;
      }
    }
    int maxsize = 5 + (s[0]=='-');
    if (base==10 && l<maxsize){
      res.type=_INT_;
      res.val = atoi(s);
      return res;
    }
    else {
      mpz_t * ptr= (mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1] ;
      mpz_init_set_str(*ptr,s,base);
      res= gen(ptr);
      return res;
    }
  }

  gen string2gen(const string & ss,bool remove_ss_quotes){
    gen res;
    res.type=_STRNG;
    res.ptr_val.ref_count = new int(1);
    if (remove_ss_quotes)
      res._STRNGptr = new string(ss.substr(1,ss.size()-2));
    else
      res._STRNGptr = new string(ss);
    return res;
  }

  int giac_yyparse(void * scanner);

  int try_parse(const string & s,GIAC_CONTEXT){
    int res;
    try {
      void * scanner;
      YY_BUFFER_STATE state=set_lexer_string(s,scanner,contextptr);
      res=giac_yyparse(scanner);
      delete_lexer_string(state,scanner);
    }
    catch (std::runtime_error & error){
      cerr << "Error " << error.what() << endl;
      messages_to_print += string(error.what()) + '\n';
      return 1;
    }
    return res;
  }
  
  gen aplatir_plus(const gen & g){
    // Quick check for embedded + at the left coming from parser
    if (g.is_symb_of_sommet(at_plus) && g._SYMBptr->feuille.type==_VECT){
      iterateur it=g._SYMBptr->feuille._VECTptr->begin(),itend=g._SYMBptr->feuille._VECTptr->end();
      if (it==itend)
	return 0;
      vecteur v;
      v.reserve(itend-it+1);
      gen f;
      for (;it!=itend;){
	for (--itend;itend!=it;--itend)
	  v.push_back(aplatir_fois_plus(*itend));
	if (!it->is_symb_of_sommet(at_plus)){
	  v.push_back(aplatir_fois_plus(*it));
	  break;
	}
	f=it->_SYMBptr->feuille;
	if (f.type!=_VECT){
	  v.push_back(*it);
	  break;
	}
	it=f._VECTptr->begin();
	itend=f._VECTptr->end();
      }
      reverse(v.begin(),v.end());
      return new symbolic(at_plus,v);
    }
    return g;
  }

  gen aplatir_fois_plus(const gen & g){
    if (g.type==_VECT){
      vecteur v(*g._VECTptr);
      iterateur it=v.begin(),itend=v.end();
      for (;it!=itend;++it)
	*it=aplatir_fois_plus(*it);
      return gen(v,g.subtype);
    }
    if (g.type!=_SYMB)
      return g;
    if (g._SYMBptr->sommet==at_plus)
      return aplatir_plus(g);
    gen & f=g._SYMBptr->feuille;
    if (g._SYMBptr->sommet==at_prod && f.type==_VECT && f._VECTptr->size()==2)
      return sym_mult(aplatir_fois_plus(f._VECTptr->front()),aplatir_fois_plus(f._VECTptr->back()),context0);
    return new symbolic(g._SYMBptr->sommet,aplatir_fois_plus(f));
  }
  
  gen aplatir_plus_only(const gen & g){
    if (g.type==_VECT){
      vecteur v(*g._VECTptr);
      iterateur it=v.begin(),itend=v.end();
      for (;it!=itend;++it)
	*it=aplatir_plus_only(*it);
      return gen(v,g.subtype);
    }
    if (g.type!=_SYMB)
      return g;
    // Quick check for embedded + at the left coming from parser
    if (g.is_symb_of_sommet(at_plus) && g._SYMBptr->feuille.type==_VECT){
      iterateur it=g._SYMBptr->feuille._VECTptr->begin(),itend=g._SYMBptr->feuille._VECTptr->end();
      if (it==itend)
	return 0;
      vecteur v;
      v.reserve(itend-it+1);
      register gen * f;
      for (;it!=itend;){
	for (--itend;itend!=it;--itend)
	  v.push_back(*itend);
	if (it->type!=_SYMB || it->_SYMBptr->sommet!=at_plus || (f=&it->_SYMBptr->feuille,f->type!=_VECT) ){
	  v.push_back(*it);
	  break;
	}
	itend=f->_VECTptr->end();
	it=f->_VECTptr->begin();
      }
      reverse(v.begin(),v.end());
      return new symbolic(at_plus,v);
    }
    return new symbolic(g._SYMBptr->sommet,aplatir_plus_only(g._SYMBptr->feuille));
  }

  int protected_giac_yyparse(const string & chaine,gen & parse_result,GIAC_CONTEXT){
    int s;
    s=chaine.size();
    if (!s)
      return 1;
    int res=try_parse(chaine,contextptr);
    gen g=parsed_gen(contextptr);
    if (g.type<=20){
      parse_result=aplatir_plus_only(g);
      // parse_result=aplatir_fois_plus(g);
      if (g.type==_SYMB && parse_result.type==_SYMB)
	parse_result.subtype=g.subtype;
      return res;
    }
    parsed_gen(0,contextptr);
    parse_result.type=0;
    parse_result=0;
    cerr << "Incomplete parse" << endl;
    return res;
  }

  gen::gen(const string & s,GIAC_CONTEXT){
    subtype=0;
    string ss(s);
    /*
      string::iterator it=ss.begin(),itend=ss.end();
      for (;it!=itend;++it)
      if (*it=='\\')
      *it=' ';
      */
    type=_INT_;
    if (s==string(s.size(),' ')){
      *this=undef;
      return;
    }
    if (protected_giac_yyparse(s,*this,contextptr)){
      type=_STRNG;
      ptr_val.ref_count = new int(1);
      if (ss.empty())
	ss="""""";
      if (ss[0]!='"')
	ss = '"'+ss;
      if ((ss.size()==1) || (ss[ss.size()-1]!='"'))
	ss += '"';
      _STRNGptr = new string(ss.substr(1,ss.size()-2));
    }
  }

  gen::gen(const string & s,const vecteur & l,GIAC_CONTEXT){
    subtype=0;
    type=_INT_;
    if (protected_giac_yyparse(s,*this,contextptr)){
      type=_STRNG;
      ptr_val.ref_count = new int(1);
      string ss(s);
      if (ss.empty())
	ss="""""";
      if (ss[0]!='"')
	ss = '"'+ss;
      if ((ss.size()==1) || (ss[ss.size()-1]!='"'))
	ss += '"';
      _STRNGptr = new string(ss.substr(1,ss.size()-2));
    }
  }

  gen::gen(const sparse_poly1 & p){
    subtype=0;
    if (p.empty()){
      type=0;
      val=0;
    }
    else {
      type=_SPOL1;
      ptr_val.ref_count = new int(1);
      _SPOL1ptr= new sparse_poly1(p);
    }
  }

  gen::gen(const unary_function_ptr & f,int nargs){
    subtype=0;
    type=_FUNC;
    subtype=nargs;
    ptr_val.ref_count = new int(1);
    _FUNCptr= new unary_function_ptr(f);
  }

  gen::~gen() {  
    if ( type>_DOUBLE_ ){
      // if (type==debugtype)
      // cout << *this << " ~ Delete @ " << *debugtypeptr << "[" << *ptr_val.ref_count <<"-1]" << endl;
      (*ptr_val.ref_count)--;
      if (!*ptr_val.ref_count){
	delete ptr_val.ref_count;
	switch (type) {
	case _ZINT: 
	  mpz_clear(*_ZINTptr); 
	  free(_ZINTptr);
	  break; 
	case _REAL:  
	  delete _REALptr;
	  break; 
	case _CPLX: 
	  delete [] _CPLXptr;
	  break; 
	case _IDNT: 
	  delete _IDNTptr;
	  break;
	case _VECT: 
	  //_VECTptr->clear();
	  delete _VECTptr;
	  break;
	case _SYMB: 
	  delete _SYMBptr;
	  break;
	case _USER:
	  _USERptr->memory_free(_USERptr);
	  break;
	case _EXT: 
	  delete [] _EXTptr;
	  break;
	case _MOD: 
	  delete [] _MODptr;
	  break;
	case _ROOT: 
	  delete _ROOTptr;
	  break;
	case _POLY:
	  _POLYptr->coord.clear();
	  delete _POLYptr;
	  break;
	case _FRAC:
	  delete _FRACptr;
	  break;
	case _SPOL1:
	  delete _SPOL1ptr;
	  break;
	case _STRNG:
	  delete _STRNGptr;
	  break;
	case _FUNC:
	  delete _FUNCptr;
	  break;
	case _MAP:
	  delete _MAPptr;
	  break;
	case _EQW:
	  delete _EQWptr;
	  break;
	case _GROB:
	  delete _GROBptr;
	  break;
	case _POINTER_:
	  if (subtype==_FL_WIDGET_POINTER && fl_widget_delete_function)
	    fl_widget_delete_function(_POINTER_val);
	  break;
	default: 
	  settypeerr("Gen Destructor");
	}
      }
    }
  }

  gen evalf_VECT(const vecteur & v,int subtype,int level,const context * contextptr){
    // bool save_is_inevalf=is_inevalf;
    // is_inevalf=true;
    vecteur w;
    vecteur::const_iterator it=v.begin(), itend=v.end();
    w.reserve(itend-it);
    for (;it!=itend;++it){
      gen tmp=it->evalf(level,contextptr);
      if (subtype){
	if ((subtype==_SEQ__VECT)&&(tmp.type==_VECT) && (tmp.subtype==_SEQ__VECT)){
	  const_iterateur jt=tmp._VECTptr->begin(),jtend=tmp._VECTptr->end();
	  for (;jt!=jtend;++jt)
	    w.push_back(*jt);
	}
	else {
	  if ((subtype!=_SET__VECT) || (!equalposcomp(w,tmp)))
	    w.push_back(tmp);
	}
      }
      else
	w.push_back(tmp);
    }
    // is_inevalf=save_is_inevalf;
    return gen(w,subtype);
  }


  bool eval_VECT(const gen & g,gen & evaled,int subtype,int level,const context * contextptr){
    const vecteur & v = *g._VECTptr;
    const_iterateur it=v.begin(),itend=v.end();
    if (subtype!=_SET__VECT && subtype!=_SEQ__VECT){
      for (;it!=itend;++it){
	if (it->in_eval(level,evaled,contextptr))
	  break;
      }
      if (it==itend)
	return false;
    }
    vecteur * vptr = new vecteur;
    vptr->reserve(itend-v.begin());
    if (subtype!=_SET__VECT && subtype!=_SEQ__VECT){
      for (const_iterateur it1=v.begin();it1!=it;++it1)
	vptr->push_back(*it1);
      if (evaled.type==_VECT && evaled.subtype==_SEQ__VECT){
	const_iterateur jt=evaled._VECTptr->begin(),jtend=evaled._VECTptr->end();
	for (;jt!=jtend;++jt){
	  if ((subtype!=_SET__VECT) || (!equalposcomp(*vptr,*jt)))
	    vptr->push_back(*jt);
	}
      }
      else {
	if ( subtype!=_SET__VECT || (!equalposcomp(*vptr,evaled)))
	  vptr->push_back(evaled);
      }
    }
    else --it;
    evaled=gen(vptr,subtype);
    gen tmp;
    const gen * ansptr;
    for (++it;it!=itend;++it){
      if ((it->type==_SYMB) && (it->_SYMBptr->sommet==at_comment))
	continue;
      ansptr=(it->in_eval(level,tmp,contextptr))?&tmp:&*it;
      if (ansptr->type==_VECT && ansptr->subtype==_SEQ__VECT){
	const_iterateur jt=ansptr->_VECTptr->begin(),jtend=ansptr->_VECTptr->end();
	for (;jt!=jtend;++jt){
	  if ((subtype!=_SET__VECT) || (!equalposcomp(*vptr,*jt)))
	    vptr->push_back(*jt);
	}
      }
      else {
	if ( subtype!=_SET__VECT || (!equalposcomp(*vptr,*ansptr)))
	  vptr->push_back(*ansptr);
      }
    }
    // cerr << "End " << v << " " << w << endl;
    return true;
  }

      // evaluate _FUNCndary in RPN mode, f must be of type _FUNC
  gen rpneval_FUNC(const gen & f,GIAC_CONTEXT){
    // int s=history_out(contextptr).size();
      int nargs=max(f.subtype,0);
      if (!nargs){
	gen res=(*f._FUNCptr)(gen(history_out(contextptr),_RPN_STACK__VECT),contextptr);
	if ( (res.type!=_VECT) || (res.subtype!=_RPN_STACK__VECT))
	  res=gen(makenewvecteur(res),_RPN_STACK__VECT);
	history_out(contextptr)=*res._VECTptr;
	history_in(contextptr)=history_out(contextptr);
	return res;
      }
      vecteur v(nargs);
      for (int i=nargs-1;i>=0;--i){
          v[i]=history_out(contextptr).back();
          history_out(contextptr).pop_back();
          history_in(contextptr).pop_back();   
      }
      if (nargs==1)
          return (*f._FUNCptr)(v.front(),contextptr);
      else
          return (*f._FUNCptr)(v,contextptr);
  }

  gen evalcomment(const vecteur & v,int level,const context * contextptr){
    const_iterateur it=v.begin(),itend=v.end();
    for (;it!=itend;++it){
      if ( (it->type!=_SYMB) || (it->_SYMBptr->sommet!=at_comment) )
	break;
    }
    if (it+1==itend)
      return it->eval(level,contextptr);
    if (it!=itend){
      gen partial=vecteur(it,itend),evaled;
      if (eval_VECT(partial,evaled,_SEQ__VECT,level,contextptr))
	return evaled;
      else
	return partial;
    }
    return zero;
  }

  bool check_not_assume(const gen & not_evaled,gen & evaled, bool evalf_after,const context * contextptr){
    if ( evaled.type==_VECT && evaled.subtype==_ASSUME__VECT ){
      if ( evalf_after && evaled._VECTptr->size()==2 && (evaled._VECTptr->back().type<=_CPLX || evaled._VECTptr->back().type==_FRAC) ){
	evaled=evaled._VECTptr->back().evalf(1,contextptr);
	return true;
      }
      if (not_evaled.type==_IDNT && evaled._VECTptr->size()==1 && evaled._VECTptr->front().type==_INT_){
	gen tmp=not_evaled;
	tmp.subtype=evaled._VECTptr->front().val;
	evaled=tmp;
	return true;
      }
      return false;
    }
    else {
      if (evalf_after && evaled.type!=_IDNT){
	gen res;
	if (has_evalf(evaled,res,0,contextptr)){
	  evaled=res;
	  return true;
	}
      }
      return &evaled!=&not_evaled;
    }
    return false;
  }

  gen gen::eval(int level,const context * contextptr) const{
    // cerr << "eval " << *this << " " << level << endl;
    control_c();
    if (level==0)
      return *this;
    gen res;
    if (in_eval(level,res,contextptr))
      return res;
    else
      return *this;
  }

  bool gen::in_eval(int level,gen & evaled,const context * contextptr) const{
    if (!level)
      return false;
    switch (type) {
    case _INT_: case _DOUBLE_: case _ZINT: case _REAL: case _CPLX: case _POLY: case _FRAC: case _SPOL1: case _EXT: case _STRNG: case _MAP: case _EQW: case _GROB: case _POINTER_:
      return false;
    case _IDNT:
      if (_IDNTptr->name==_IDNT_pi.name || _IDNTptr->name==_IDNT_euler_gamma.name )
	return false;
      if (!contextptr && subtype==_GLOBAL__EVAL)
	evaled=global_eval(*_IDNTptr,level);
      else {
	if (!_IDNTptr->in_eval(level-1,*this,evaled,contextptr))
	  return false;
      }
      if ( evaled.type!=_VECT || evaled.subtype!=_ASSUME__VECT )
	return true;
      return check_not_assume(*this,evaled,false,contextptr);
    case _VECT:
      if (subtype==_SPREAD__VECT){
	makespreadsheetmatrice(*_VECTptr,contextptr);
	spread_eval(*_VECTptr,contextptr);
	return false;
      }
      if (subtype==_FOLDER__VECT || subtype==_RGBA__VECT)
	return false;
      if ( (subtype==_SEQ__VECT) && (!_VECTptr->empty()) && (_VECTptr->front().type==_SYMB) && (_VECTptr->front()._SYMBptr->sommet==at_comment) && (_VECTptr->back().type==_SYMB) && (_VECTptr->back()._SYMBptr->feuille==at_return ) ){
	evaled=evalcomment(*_VECTptr,level,contextptr);
	return true;
      }
      return eval_VECT(*this,evaled,subtype,level,contextptr);
    case _SYMB:
      if (subtype==_SPREAD__SYMB)
	return false;
      if (_SYMBptr->sommet==at_plus || _SYMBptr->sommet==at_prod || _SYMBptr->sommet==at_pow){
	if (_SYMBptr->feuille.in_eval(level,evaled,contextptr))
	  evaled=(*_SYMBptr->sommet.ptr)(evaled,contextptr);
	else // ?? return false; ??
	  evaled=(*_SYMBptr->sommet.ptr)(_SYMBptr->feuille,contextptr);
      }
      else
	evaled=_SYMBptr->eval(level,contextptr);
      return true;
    case _USER:
      evaled=_USERptr->eval(level,contextptr);
      return true;
    case _MOD:
      evaled=makemod(_MODptr->eval(level,contextptr),(_MODptr+1)->eval(level,contextptr));
      return true;
    case _FUNC:
      if (rpn_mode && (history_out(contextptr).size()>=unsigned(subtype)))
	evaled=rpneval_FUNC(*this,contextptr);
      else {
	if (subtype)
	  return false; 
	else
	  evaled=(*_FUNCptr)(gen(vecteur(0),_SEQ__VECT),contextptr);
      }
      return true;
    default: 
      settypeerr("Eval") ;
    }
    return false;
  }

  symbolic _FRAC2_SYMB(const fraction & f){
    if (is_one(f.num))
      return symb_inv(f.den);
    if (is_minus_one(f.num))
      return symb_inv(-f.den);      
    return symbolic(at_prod,makenewvecteur(f.num,symb_inv(f.den)));
  }

  symbolic _FRAC2_SYMB(const gen & e){
#ifdef DEBUG_SUPPORT
    if (e.type!=_FRAC) setsizeerr("gen.cc/_FRAC2_SYMB");
#endif
    return _FRAC2_SYMB(*e._FRACptr);
  }

  symbolic _FRAC2_SYMB(const gen & n,const gen & d){
    return symbolic(at_prod,makenewvecteur(n,symb_inv(d)));
  }

  polynome apply( const polynome & p, const gen_op & f){
    polynome res(p.dim);
    std::vector< monomial<gen> > :: const_iterator it=p.coord.begin(),itend=p.coord.end();
    res.coord.reserve(itend-it);
    for (;it!=itend;++it){
      gen tmp(f(it->value));
      if (!is_zero(tmp))
	res.coord.push_back(monomial<gen>(tmp,it->index));
    }
    return res;
  }

  polynome apply( const polynome & p, const context * contextptr, gen (* f) (const gen &, const context *)){
    polynome res(p.dim);
    std::vector< monomial<gen> > :: const_iterator it=p.coord.begin(),itend=p.coord.end();
    res.coord.reserve(itend-it);
    for (;it!=itend;++it){
      gen tmp(f(it->value,contextptr));
      if (!is_zero(tmp))
	res.coord.push_back(monomial<gen>(tmp,it->index));
    }
    return res;
  }

  gen evalf(const gen & e,int level,const context * contextptr){ 
    return e.evalf(level,contextptr); 
  }

  gen no_context_evalf(const gen & e){
    gen tmp;
    if (has_evalf(e,tmp,1,context0))
      return tmp;
    else
      return e;
  }

  double double0_15[]={0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0};
  double _int2double(unsigned i){
    if (i<16)
      return double0_15[i];
    else
      return _int2double(i/16)*16.0+double0_15[i%16];
  }
  // double(int) does not seem to work under GNUWINCE
  double int2double(int i){
    if (i<0)
      return -_int2double(-i);
    else
      return _int2double(i);
  }


  gen gen::evalf(int level,const context * contextptr) const{
    // cerr << "evalf " << *this << " " << level << endl;
    control_c();
    if (level==0)
      return *this;
    gen evaled;
    if (in_evalf(level,evaled,contextptr))
      return evaled;
    else
      return *this;
  }

  gen m_pi(GIAC_CONTEXT){
    int nbits=digits2bits(decimal_digits(contextptr));
    return m_pi(nbits);
  }

  gen m_pi(int nbits){
#ifdef HAVE_LIBMPFR
    if (nbits>52){
#ifdef HAVE_LIBPTHREAD
      int locked=pthread_mutex_lock(&mpfr_mutex);
#else // HAVE_LIBPTHREAD
      int locked=0;
#endif
      if (!locked)
	mpfr_set_default_prec(nbits);
      mpfr_t pi;
      mpfr_init(pi);
      mpfr_const_pi(pi,GMP_RNDN);
#ifdef HAVE_LIBPTHREAD
      if (!locked)
	pthread_mutex_unlock(&mpfr_mutex);
#endif
      gen res=real_object(pi);
      mpfr_clear(pi);
      return res;
    }
#endif
    return M_PI;
  }

  gen m_gamma(int nbits){
#ifdef HAVE_LIBMPFR
    if (nbits>15){
#ifdef HAVE_LIBPTHREAD
      int locked=pthread_mutex_lock(&mpfr_mutex);
#else // HAVE_LIBPTHREAD
      int locked=0;
#endif
      if (!locked)
	mpfr_set_default_prec(nbits);
      mpfr_t euler_gamma;
      mpfr_init(euler_gamma);
      mpfr_const_euler(euler_gamma,GMP_RNDN);
#ifdef HAVE_LIBPTHREAD
      if (!locked)
	pthread_mutex_unlock(&mpfr_mutex);
#endif
      gen res=real_object(euler_gamma);
      mpfr_clear(euler_gamma);
      return res;
    }
#endif
    return .577215664901533;
  }

  gen m_gamma(GIAC_CONTEXT){
    int nbits=digits2bits(decimal_digits(contextptr));
    return m_gamma(nbits);
  }

  bool gen::in_evalf(int level,gen & evaled,const context * contextptr) const{
    if (!level)
      return false;
    switch (type) {
    case _DOUBLE_: case _REAL: case _STRNG: case _MAP: case _EQW: case _GROB:
      return false;
    case _INT_:
      if (subtype)
	return false;
      if (decimal_digits(contextptr)>14){
	evaled=real_object(*this,digits2bits(decimal_digits(contextptr)));
	return true;
      }
      evaled=int2double(val);
      return true;
    case _ZINT:
      if (decimal_digits(contextptr)>14)
	evaled=real_object(*this,digits2bits(decimal_digits(contextptr)));
      else
	evaled=mpz_get_d(*_ZINTptr);
      return true;
    case _CPLX: 
      evaled=gen(_CPLXptr->evalf(level,contextptr),(_CPLXptr+1)->evalf(level,contextptr));
      return true;
    case _USER:
      evaled=_USERptr->evalf(level,contextptr);
      return true;
    case _IDNT:
      if (_IDNTptr->name==_IDNT_pi.name){
	evaled=m_pi(contextptr);
	return true;
      }
      if (_IDNTptr->name==_IDNT_euler_gamma.name){
	evaled=m_gamma(contextptr);
	return true;
      }
      if (!contextptr && subtype==_GLOBAL__EVAL)
	evaled=global_evalf(*_IDNTptr,level-1);
      else {
	if (!_IDNTptr->in_eval(level-1,*this,evaled,contextptr))
	  return false;
      }
      return check_not_assume(*this,evaled,true,contextptr);
    case _VECT:
      evaled=evalf_VECT(*_VECTptr,subtype,level,contextptr);
      return true;
    case _SYMB:
      if (_SYMBptr->sommet==at_pow && _SYMBptr->feuille._VECTptr->back().type==_INT_){
	evaled=pow(_SYMBptr->feuille._VECTptr->front().evalf(level,contextptr),_SYMBptr->feuille._VECTptr->back(),contextptr);
	return true;
      }
      if (_SYMBptr->sommet==at_integrate){
	evaled=_romberg(_SYMBptr->feuille,contextptr);
	return true;
      }
      if (_SYMBptr->sommet==at_rootof){
	evaled=approx_rootof(_SYMBptr->feuille.evalf(level,contextptr),contextptr);
	return true;
      }
      if (_SYMBptr->sommet==at_cell)
	return false;
      evaled=_SYMBptr->evalf(level,contextptr);
      return true;
    case _FRAC:
      evaled=rdiv(_FRACptr->num.evalf(level,contextptr),_FRACptr->den.evalf(level,contextptr));
      return true;
    case _FUNC: case _MOD: case _ROOT:
      return false; // replace in RPN mode
    case _EXT:
      evaled=alg_evalf(_EXTptr->eval(level,contextptr),(_EXTptr+1)->eval(level,contextptr),contextptr);
      return true;
    case _POLY:
      evaled=apply(*_POLYptr,giac::no_context_evalf);
      return true;
    default: 
      settypeerr("Evalf") ;
    }
    return false;
  }

  double real_object::evalf_double() const{
#ifdef HAVE_LIBMPFR
    return mpfr_get_d(inf,GMP_RNDN);
#else
    return mpf_get_d(inf);
#endif
  }

  gen real2int(const gen & g,GIAC_CONTEXT){
    if (g.type==_REAL){
      if (is_strictly_positive(-g,contextptr))
	return -real2int(-g,contextptr);
#ifdef HAVE_LIBMPFR
      mpz_t * m=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init(*m);
      int n=int(mpfr_get_z_exp(*m,g._REALptr->inf));
      gen res(*m);
      if (n>=0)
	return res*pow(plus_two,gen(n),contextptr);
      return _iquo(makenewvecteur(res,pow(plus_two,gen(-n),contextptr)));
#else
      return g;
#endif
    }
    if (g.type!=_VECT)
      return g;
    return apply(g,real2int,contextptr);
  }

  gen real2double(const gen & g){
    if (g.type==_REAL)
      return g._REALptr->evalf_double();
    if (g.type!=_VECT)
      return g;
    return apply(g,real2double);
  }

  gen gen::evalf_double(int level,const context * contextptr) const{
    gen g;
    if (has_evalf(*this,g,level,contextptr)){
      if (g.type==_CPLX)
	return gen(real2double(*g._CPLXptr),real2double(*(g._CPLXptr+1)));
      else
	return real2double(g);
    }
    else
      return *this;
  }

  gen gen::evalf2double(int level,const context * contextptr) const{
    gen g=evalf(level,contextptr);
    return g.evalf_double(level,contextptr);
  }
  
  real_object & real_object::operator = (const real_object & g) { 
#ifdef HAVE_LIBMPFR
    mpfr_clear(inf);
    mpfr_init2(inf,mpfr_get_prec(g.inf));
    mpfr_set(inf,g.inf,GMP_RNDN);
#else
    mpf_clear(inf);
    mpf_init_set(inf,g.inf);
#endif
    return *this;
  }

  real_object & real_interval::operator = (const real_interval & g) { 
#ifdef HAVE_LIBMPFR
    mpfr_clear(inf);
#else
    mpf_clear(inf);
#endif
#ifdef HAVE_LIBMPFI
    mpfi_clear(infsup); 
#else
#ifdef HAVE_LIBMPFR
    mpfr_clear(sup);
#else
    mpf_clear(sup);
#endif
#endif
#ifdef HAVE_LIBMPFR
    mpfr_init2(inf,mpfr_get_prec(g.inf));
    mpfr_set(inf,g.inf,GMP_RNDN);
#else
    mpf_init_set(inf,g.inf);
#endif
#ifdef HAVE_LIBMPFI
    mpfi_init_set(infsup,g.infsup);
#else
#ifdef HAVE_LIBMPFR
    mpfr_init2(sup,mpfr_get_prec(g.sup));
    mpfr_set(sup,g.sup,GMP_RNDN);
#else
    mpf_init_set(sup,g.sup);
#endif
#endif
    return *this;
  }

  real_object & real_interval::operator = (const real_object & g) { 
    const real_interval * ptr=dynamic_cast<const real_interval *> (&g);
    if (ptr)
      return *this=*ptr;
#ifdef HAVE_LIBMPFR
    mpfr_clear(inf);
#ifdef HAVE_LIBMPFI
    mpfi_clear(infsup); 
#else
    mpfr_clear(sup); 
#endif
    mpfr_init2(inf,mpfr_get_prec(g.inf));
    mpfr_set(inf,g.inf,GMP_RNDN);
#ifdef HAVE_LIBMPFI
    mpfi_init_set_fr(infsup,g.inf);
#else
    mpfr_init2(sup,mpfr_get_prec(g.inf));
    mpfr_set(sup,g.inf,GMP_RNDN);
#endif
#else // HAVE_LIBMPFR
    mpf_clear(inf);
#ifdef HAVE_LIBMPFI
    mpfi_clear(infsup); 
#else
    mpf_clear(sup); 
#endif
    mpf_init_set(inf,g.inf);
#ifdef HAVE_LIBMPFI
    mpfi_init_set_fr(infsup,g.inf);
#else
    mpf_init_set(sup,g.inf);
#endif
#endif // HAVE_LIBMPFR
    return *this;
  }

  gen & gen::operator = (const gen & a) { 
    register unsigned t=(type << _DECALAGE) | a.type;
    if (!t){
      subtype=a.subtype;
      val=a.val;
      return *this;
    }
    // Copy before deleting because the target might be embedded in a
    // with a ptr_val.ref_count of a equals to 1
    int * ref_count_save=ptr_val.ref_count;
    short int type_save=type; // short int subtype_save=subtype; 
    void * ptr_save = _ZINTptr;
    subtype=a.subtype;
    type=a.type;
    // if (type==debugtype)
    // cout << a << "Overwriting by = " << *a.debugtypeptr << "[" << *(a.ptr_val.ref_count) << "+1]"<< endl ;
    switch (type) {
    case _INT_:
      val=a.val;
      break;
    case _DOUBLE_:
      _DOUBLE_val=a._DOUBLE_val;
      break;
    case _POINTER_:
      ptr_val.ref_count=a.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _POINTER_val=a._POINTER_val;
      break;
    case _ZINT: 
      ptr_val.ref_count=a.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _ZINTptr=a._ZINTptr;
      break; 
    case _REAL: 
      ptr_val.ref_count=a.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _REALptr=a._REALptr;
      break; 
    case _CPLX: 
      ptr_val.ref_count=a.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _CPLXptr=a._CPLXptr;
      break;
    case _IDNT: 
      ptr_val.ref_count=a.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _IDNTptr=a._IDNTptr;
      break;
    case _VECT: 
      ptr_val.ref_count=a.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _VECTptr=a._VECTptr;
      break;
    case _SYMB: 
      ptr_val.ref_count=a.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _SYMBptr=a._SYMBptr;
      break;
    case _USER: 
      ptr_val.ref_count=a.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _USERptr=a._USERptr;
      break;
    case _EXT: 
      ptr_val.ref_count=a.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _EXTptr=a._EXTptr;
      break;
    case _MOD: 
      ptr_val.ref_count=a.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _MODptr=a._MODptr;
      break;
    case _ROOT: 
      ptr_val.ref_count=a.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _ROOTptr=a._ROOTptr;
      break;
    case _POLY: 
      ptr_val.ref_count=a.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _POLYptr=a._POLYptr;
      break;
    case _FRAC: 
      ptr_val.ref_count=a.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _FRACptr=a._FRACptr;
      break;
    case _SPOL1: 
      ptr_val.ref_count=a.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _SPOL1ptr=a._SPOL1ptr;
      break;
    case _STRNG: 
      ptr_val.ref_count=a.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _STRNGptr=a._STRNGptr;
      break;
    case _FUNC: 
      ptr_val.ref_count=a.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _FUNCptr=a._FUNCptr;
      break;
    case _MAP: 
      ptr_val.ref_count=a.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _MAPptr=a._MAPptr;
      break;
    case _EQW: 
      ptr_val.ref_count=a.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _EQWptr=a._EQWptr;
      break;
    case _GROB: 
      ptr_val.ref_count=a.ptr_val.ref_count;
      (*ptr_val.ref_count)++;
      _GROBptr=a._GROBptr;
      break;
    default: 
      break ;
    }
    // Now we delete the copied object
    if ( type_save>_DOUBLE_ ){
      control_c();
      // if (type==debugtype)
      // cout << *this << " Destruction by = of " << *debugtypeptr << "[" << *ptr_val.ref_count << "-1]" << endl;
      (*ref_count_save)--;
      if (!*ref_count_save){
	delete ref_count_save;
	switch (type_save) {
	case _ZINT: 
	  mpz_clear( * (mpz_t *) ptr_save);
	  free(ptr_save); // delete (mpz_t *) ptr_save;
	  break; 
	case _REAL: 
	  delete (real_object *) ptr_save;
	  break; 
	case _CPLX: 
	  delete [] (gen *) ptr_save ;
	  break; 
	case _IDNT: 
	  delete (identificateur *) ptr_save ;
	  break;
	case _SYMB: 
	  delete (symbolic *) ptr_save;
	  break;
	case _USER: 
	  ((gen_user *) ptr_save)->memory_free((gen_user *) ptr_save);
	  break;
	case _EXT: 
	  delete [] (gen * ) ptr_save;
	  break;
	case _MOD: 
	  delete [] (gen * ) ptr_save;
	  break;
	case _ROOT:
	  delete (real_complex_rootof *) ptr_save;
	  break;
	case _VECT: 
	  //((vecteur *) ptr_save)->clear();
	  delete (vecteur *) ptr_save;
	  break;
	case _POLY:
	  //((polynome *) ptr_save)->coord.clear();
	  delete (polynome *) ptr_save;
	  break;
	case _FRAC:
	  delete (fraction *) ptr_save;
	  break;
	case _SPOL1:
	  delete (sparse_poly1 *) ptr_save;
	  break;
	case _STRNG:
	  delete (string *) ptr_save;
	  break;
	case _FUNC:
	  delete (unary_function_ptr *) ptr_save;
	  break;
	case _MAP:
	  delete (gen_map *) ptr_save;
	  break;
	case _EQW:
	  delete (eqwdata *) ptr_save;
	  break;
	case _GROB:
	  delete (grob *) ptr_save;
	  break;
	case _POINTER_:
	  if (subtype==_FL_WIDGET_POINTER && fl_widget_delete_function)
	    fl_widget_delete_function(ptr_save);
	  break;
	default: 
	  settypeerr("Gen Operator =");
	}
      }
    }
    return *this;
  }
  
  int gen::to_int() const {
    switch (type ) {
    case _INT_: 
      return val;
    case _ZINT: 
      return mpz_get_si(*_ZINTptr);
    case _CPLX: 
      return _CPLXptr->to_int() ;
    default:
      settypeerr("To_int");
    }
    return 0;
  }
  
  bool poly_is_real(const polynome & p){
    vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    for (;it!=itend;++it){
      if (!it->value.is_real(0)) // context is 0 since coeff do not depend on
	return false;
    }
    return true;
  }

  bool gen::is_real(GIAC_CONTEXT) const {
    switch (type) {
    case _INT_: case _DOUBLE_: case _ZINT: case _REAL:
      return true; 
    case _CPLX: 
      return (is_zero(*(_CPLXptr+1)));
    case _POLY:
      return poly_is_real(*_POLYptr);
    default: 
      return is_zero(im(contextptr));
    }
  }

  bool gen::is_approx() const {
    switch(type){
    case _DOUBLE_: case _REAL:
      return true;
    case _CPLX:
      return subtype==3 || (_CPLXptr->is_approx() && (_CPLXptr+1)->is_approx());
    default:
      return false;
    }
  }
  
  bool gen::is_cinteger() const {
    switch (type ) {
    case _INT_: case _ZINT: 
      return true; 
    case _CPLX:
      return _CPLXptr->is_integer() && (_CPLXptr+1)->is_integer();
    default: 
      return false;
    }
  }
   
  bool gen::is_integer() const {
    switch (type ) {
    case _INT_: case _ZINT:
      return true;
    case _CPLX:
      return is_exactly_zero(*(_CPLXptr+1)) && _CPLXptr->is_integer();
    default: 
      return false;
    }
  }

  bool _VECT_is_constant(const vecteur & v){
      const_iterateur it=v.begin(),itend=v.end();
      for (;it!=itend;++it)
          if (!(it->is_constant()))
              return false;
      return true;
  }
  
  bool gen::is_constant() const {
    switch (type ) {
    case _INT_: case _DOUBLE_: case _REAL: case _ZINT: case _CPLX:
      return true;
    case _VECT:
      return _VECT_is_constant(*this->_VECTptr);
    case _EXT:
      return _EXTptr->is_constant() && (_EXTptr+1)->is_constant();
    case _POLY:
      return Tis_constant<gen>(*_POLYptr) && _POLYptr->coord.front().value.is_constant();
    default: 
      return false;
    }
  }

  bool is_atomic(const gen & e){
    return (e.type<_POLY);
  }

  void gen::uncoerce() {
    if (type==_INT_){
      int tmp =val;
      _ZINTptr = (mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init_set_si(*_ZINTptr,tmp); 
      ptr_val.ref_count=new int(1);
      type=_ZINT;
    }
  }

  vecteur _VECTconj(const vecteur & a,GIAC_CONTEXT){
    vecteur res;
    vecteur::const_iterator it=a.begin(),itend=a.end();
    for (;it!=itend;++it)
      res.push_back(it->conj(contextptr));
    return res;
  }

  gen gen::conj(GIAC_CONTEXT) const {
    switch (type ) {
    case _INT_: case _DOUBLE_: case _ZINT: case _REAL: case _STRNG:
      return *this;
    case _CPLX:
      return gen(*_CPLXptr,-(*(_CPLXptr+1)));
    case _VECT:
      return gen(_VECTconj(*_VECTptr,contextptr),subtype);
    case _USER:
      return _USERptr->conj(contextptr);
    case _IDNT: 
      if (!complex_variables(contextptr)) 
	return *this;
      /* if ( (_IDNTptr->value) && (is_zero(_IDNTptr->value->im())) )
	 return *this; */
      return new symbolic(at_conj,*this);
    case _SYMB:
      if (equalposcomp(plot_sommets,_SYMBptr->sommet) || equalposcomp(analytic_sommets,_SYMBptr->sommet))
	return new symbolic(_SYMBptr->sommet,_SYMBptr->feuille.conj(contextptr));
      else
	return new symbolic(at_conj,*this);
    case _FRAC:
      return fraction(_FRACptr->num.conj(contextptr),_FRACptr->den.conj(contextptr));
    case _MOD:
      return makemod(_MODptr->conj(contextptr),*(_MODptr+1));
    case _EXT:
      return algebraic_EXTension(_EXTptr->conj(contextptr),*(_EXTptr+1));
    default: 
      settypeerr("Conj");
    }
    return 0;
  }
   
  vecteur _VECTre(const vecteur & a,GIAC_CONTEXT){
    vecteur res;
    vecteur::const_iterator it=a.begin(),itend=a.end();
    for (;it!=itend;++it)
      res.push_back(it->re(contextptr));
    return res;
  }

  vecteur pascal_next_line(const vecteur & v){
    if (v.empty())
      return vecteur(1,plus_one);
    const_iterateur it=v.begin(),itend=v.end();
    gen current(*it);
    vecteur w;
    w.reserve(itend-it+1);
    w.push_back(current);
    for (++it;it!=itend;++it){
      w.push_back(*it+current);
      current=*it;
    }
    w.push_back(plus_one);
    return w;
  }

  vecteur pascal_nth_line(int n){
    n=absint(n);
    vecteur v(1,plus_one);
    for (int i=0;i<n;++i)
      v=pascal_next_line(v);
    return v;
  }

  gen algtrim(const gen & g){
    if (g.type!=_VECT)
      return g;
    vecteur tmp=trim(*g._VECTptr,0);
    if (tmp.empty())
      return zero;
    if (tmp.size()==1)
      return tmp.front();
    return tmp;
  }

  void symb_reim(const symbolic & s,gen & r,gen & i,GIAC_CONTEXT){
    unary_function_ptr u=s.sommet;
    gen f=s.feuille;
    if ( (u==at_re) || (u==at_im) || (u==at_abs) ){
      r=s;
      i=0;
      return;
    }
    if (u==at_conj){
      reim(f,r,i,contextptr);
      i=-i;
      return;
    }
    if (u==at_plus){
      reim(f,r,i,contextptr);
      r=_plus(r,contextptr);
      i=_plus(i,contextptr);
      return;
    }
    if (u==at_neg){
      reim(f,r,i,contextptr);
      r=-r;
      i=-i;
      return;
    }
    if (u==at_division){
      reim(f[0]*inv(f[1],contextptr),r,i,contextptr);
      return ;
    }
    if (u==at_sqrt){
      reim(pow(f,plus_one_half,contextptr),r,i,contextptr);
      return;
    }
    if (u==at_prod){
      if (f.type!=_VECT){
	reim(f,r,i,contextptr);
	return;
      }
      vecteur v(*f._VECTptr);
      if (v.empty()){
	r=plus_one;
	i=0;
	return;
      }
      if (v.size()==1){
	reim(v.front(),r,i,contextptr);
	return;
      }
      // cut v in 2 parts and recursive call
      // re(a*b)=re(a)*re(b)-im(a)*im(b)
      const_iterateur it=v.begin(),itend=v.end();
      const_iterateur itm=it+(itend-it+1)/2;
      gen a(new symbolic(u,vecteur(it,itm)));
      gen b(new symbolic(u,vecteur(itm,itend)));
      gen ra,rb,ia,ib;
      reim(a,ra,ia,contextptr);
      reim(b,rb,ib,contextptr);
      r=ra*rb-ia*ib;
      i=ra*ib+rb*ia;
      return;
    }
    if (u==at_pow){
      gen e=f._VECTptr->front(),expo=f._VECTptr->back();
      if (expo.type==_INT_){
	int n=expo.val;
	if (n==0){
	  r=1;
	  i=0;
	  return;
	}
	reim(e,r,i,contextptr);
	if (n==1)
	  return ;
	if (is_zero(i)){
	  r=pow(r,n);
	  return;
	}
	if (is_zero(r)){
	  if (n%2){
	    r=zero;
	    i=pow(i,n);
	    if (n%4==3)
	      i=-i;
	    return;
	  }
	  r=pow(i,n);
	  i=0;
	  if (n%4==2)
	    r=-r;
	  return;
	}
	bool n_pos=(n>0);
	if (!n_pos){
	  reim(inv(pow(e,-n),contextptr),r,i,contextptr);
	  return;
	}
	vecteur v=pascal_nth_line(n);
	vecteur sommer,sommei;
	gen signer=1; 
	const_iterateur it=v.begin(); 
	for (int j=0;j<=n;j+=2){
	  sommer.push_back(signer*(*it)*pow(r,n-j)*pow(i,j));
	  ++it;
	  ++it;
	  signer=-signer;
	}
	it=v.begin(); 
	gen signei=1; 
	++it;
	for (int j=1;j<=n;j+=2){
	  sommei.push_back(signei*(*it)*pow(r,n-j)*pow(i,j));
	  ++it;
	  ++it;
	  signei=-signei;
	}
	r=new symbolic(at_plus,sommer);
	i=new symbolic(at_plus,sommei);
	return ;
      } // end integer exponent
      if ( is_zero(im(expo,contextptr)) && is_zero(im(e,contextptr)) ){
	if (!is_integer(expo) && is_positive(-e,contextptr)){
	  r=pow(-e,expo,contextptr)*cos(cst_pi*expo,contextptr);
	  i=pow(-e,expo,contextptr)*sin(cst_pi*expo,contextptr);
	}
	else {
	  r=s;
	  i=0;
	}
	return;
      }
    }
    // FIXME, should check that the rootof is really real, see also for im
    if (u==at_rootof && f.type==_VECT && f._VECTptr->size()==2){
      vecteur tmp=*f._VECTptr;
      reim(tmp[0],r,i,contextptr);
      r=algtrim(r);
      if (r.type==_VECT){
	tmp[0]=r;
	r=new symbolic(u,gen(tmp,f.subtype));
      }
      i=algtrim(i);
      if (i.type==_VECT){
	tmp[0]=i;
	i=new symbolic(u,gen(tmp,f.subtype));
      }
      return;
    }
    gen ref,imf;
    reim(f,ref,imf,contextptr);
    if (is_zero(imf) && equalposcomp(reim_op,u)){
      r=s; i=0; return;
    }
    if (u==at_ln){ // FIXME?? might recurse
      r=ln(abs(f,contextptr),contextptr);
      i=arg(f,contextptr);
      return ;
    }
    if (u==at_tan){
      reim(rdiv(sin(f,contextptr),cos(f,contextptr)),r,i,contextptr);
      return;
    }
    if (u==at_tanh){
      reim(rdiv(sinh(f,contextptr),cosh(f,contextptr)),r,i,contextptr);
      return;
    }
    if ((u==at_asin || u==at_acos) && is_zero(imf) && is_greater(1,f,contextptr) && is_greater(f,-1,contextptr)){
      r=s; i=0; return; 
    }
    if (u==at_inv){
      gen tmp=inv(pow(ref,2)+pow(imf,2),contextptr);
      r=ref*tmp;
      i=-imf*tmp;
      return;
    }
    if (u==at_exp) {
      // FIXME?? exp might recurse
      r=exp(ref,contextptr)*cos(imf,contextptr);
      i=exp(ref,contextptr)*sin(imf,contextptr);
      return;
    }
    if (u==at_cos){
      r=cosh(imf,contextptr)*cos(ref,contextptr);
      i=-sinh(imf,contextptr)*sin(ref,contextptr);
      return;
    }
    if (u==at_sin){
      r=cosh(imf,contextptr)*sin(ref,contextptr);
      i=sinh(imf,contextptr)*cos(ref,contextptr);
      return;
    }
    if (u==at_cosh){
      r=cos(imf,contextptr)*cosh(ref,contextptr);
      i=sin(imf,contextptr)*sinh(ref,contextptr);
      return;
    }
    if (u==at_sinh){
      r=cos(imf,contextptr)*sinh(ref,contextptr);
      i=sin(imf,contextptr)*cosh(ref,contextptr);
      return;
    }
    r=new symbolic(at_re,gen(s));
    i=new symbolic(at_im,gen(s));
  }

  void reim_poly(const polynome & p,gen & r,gen & i,GIAC_CONTEXT){
    polynome R(p.dim),I(p.dim);
    vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    for (;it!=itend;++it){
      reim(it->value,r,i,contextptr);
      if (!is_zero(r))
	R.coord.push_back(monomial<gen>(r,it->index));
      if (!is_zero(i))
	I.coord.push_back(monomial<gen>(i,it->index));
    }
    r=R;
    i=I;
  }

  void reim_spol(const sparse_poly1 & p,gen & r,gen & i,GIAC_CONTEXT){
    sparse_poly1 R,I;
    sparse_poly1::const_iterator it=p.begin(),itend=p.end();
    for (;it!=itend;++it){
      reim(it->coeff,r,i,contextptr);
      if (!is_zero(r))
	R.push_back(monome(r,it->exponent));
      if (!is_zero(i))
	I.push_back(monome(i,it->exponent));
    }
    r=R;
    i=I;
  }

  void reim_vect(const vecteur & v,gen & r,gen & i,int subtype,GIAC_CONTEXT){
    const_iterateur it=v.begin(),itend=v.end();
    vecteur R,I;
    R.reserve(itend-it);
    I.reserve(itend-it);
    for (;it!=itend;++it){
      reim(*it,r,i,contextptr);
      R.push_back(r);
      I.push_back(i);
    }
    if (subtype==_POLY1__VECT){
      R=trim(R,0);
      I=trim(I,0);
    }
    r=gen(R,subtype);
    i=gen(I,subtype);
  }
  
  gen frac_reim(const gen & n,const gen & d,bool findre,GIAC_CONTEXT){
    gen dbar(conj(d,contextptr)),tmp(n*dbar);
    tmp=findre?re(tmp,contextptr):im(tmp,contextptr);
    return tmp/(d*dbar);
  }

  // compute simultaneously real and imaginary part
  void reim(const gen & g,gen & r,gen & i,GIAC_CONTEXT){
    switch (g.type ) {
    case _INT_: case _DOUBLE_: case _ZINT: case _REAL: case _STRNG:
      r=g;
      i=0;
      break;
    case _CPLX:
      r=*g._CPLXptr;
      i=*(g._CPLXptr+1);
      break;
    case _VECT:
      reim_vect(*g._VECTptr,r,i,g.subtype,contextptr);
      break;
    case _IDNT: 
      if (!complex_variables(contextptr) || g==cst_euler_gamma || g==cst_pi){
	r=g;
	i=0;
      }
      else {
	r=new symbolic(at_re,g);
	i=new symbolic(at_im,g);
      }
      break;
    case _SYMB:
      if (equalposcomp(plot_sommets,g._SYMBptr->sommet)){
	reim(g._SYMBptr->feuille,r,i,contextptr);
	r=new symbolic(g._SYMBptr->sommet,r);
	i=new symbolic(g._SYMBptr->sommet,i);
      }
      else {
	if (expand_re_im(contextptr))
	  symb_reim(*g._SYMBptr,r,i,contextptr);
	else {
	  r=new symbolic(at_re,g);
	  i=new symbolic(at_im,g);
	}
      }
      break;
    case _USER:
      r=g._USERptr->re(contextptr);
      i=g._USERptr->im(contextptr);
      break;
    case _FRAC:
      r=frac_reim(g._FRACptr->num,g._FRACptr->den,true,contextptr);
      i=frac_reim(g._FRACptr->num,g._FRACptr->den,false,contextptr);
      break;
    case _MOD:
      reim(*g._MODptr,r,i,contextptr);
      r=makemod(r,*(g._MODptr+1));
      i=makemod(i,*(g._MODptr+1));
      break;
    case _EXT:
      reim(*g._EXTptr,r,i,contextptr);
      r=algebraic_EXTension(r,*(g._EXTptr+1));
      i=algebraic_EXTension(i,*(g._EXTptr+1));
      break;
    case _POLY:
      reim_poly(*g._POLYptr,r,i,contextptr);
      break;
    case _SPOL1:
      reim_spol(*g._SPOL1ptr,r,i,contextptr);
      break;
    default: 
      settypeerr("reim");
    }
  }

  gen symb_re(const symbolic & s,GIAC_CONTEXT){
    unary_function_ptr u=s.sommet;
    gen f=s.feuille;
    if ( (u==at_re) || (u==at_im) || (u==at_abs) )// re(re), re(im), re(abs)
      return s;
    if (u==at_conj)
      return re(f,contextptr);
    if (u==at_plus)
      return _plus(re(f,contextptr),contextptr);
    if (u==at_neg)
      return -re(f,contextptr);
    if (u==at_pow){
      gen e=f._VECTptr->front(),expo=f._VECTptr->back();
      if (expo.type==_INT_){
	int n=expo.val;
	if (n==0)
	  return plus_one;
	// ? compute conj and use 1/2*(z+-zbar)?
	gen r=re(e,contextptr);
	if (n==1)
	  return r;
	gen i=im(e,contextptr);
	if (n==2)
	  return pow(r,2)-pow(i,2);
	if (is_zero(i))
	  return pow(r,n);
	if (is_zero(r)){
	  if (n%2)
	    return zero;
	  if (n%4==2)
	    return -pow(i,n);
	  else
	    return pow(i,n);
	}
	bool n_pos=(n>0);
	if (!n_pos)
	  return re(inv(pow(e,-n),contextptr),contextptr);
	vecteur v=pascal_nth_line(n);
	vecteur somme;
	gen signe=plus_one; 
	const_iterateur it=v.begin(); //,itend=v.end();
	for (int j=0;j<=n;j+=2){
	  somme.push_back(signe*(*it)*pow(r,n-j)*pow(i,j));
	  ++it;
	  ++it;
	  signe=-signe;
	}
	gen res=new symbolic(at_plus,somme);
	return res;
      } // end integer exponent
      if ( is_zero(im(expo,contextptr)) && is_zero(im(e,contextptr)) ){
	if (!is_integer(expo) && is_positive(-e,contextptr))
	  return pow(-e,expo,contextptr)*cos(cst_pi*expo,contextptr);
	return s;
      }
    }
    // FIXME, should check that the rootof is really real, see also for im
    if (u==at_rootof && f.type==_VECT && f._VECTptr->size()==2){
      vecteur tmp=*f._VECTptr;
      tmp[0]=algtrim(re(tmp[0],contextptr));
      if (tmp[0].type!=_VECT)
	return tmp[0];
      return new symbolic(u,gen(tmp,f.subtype));
    }
    if (u==at_ln) // FIXME?? might recurse
      return ln(abs(f,contextptr),contextptr); 
    gen r,i;
    symb_reim(s,r,i,contextptr);
    return r;
  }

  gen no_context_re(const gen & a){ 
    return re(a,context0); 
  }

  gen no_context_im(const gen & a){ 
    return im(a,context0); 
  }

  gen no_context_conj(const gen & a){ 
    return conj(a,context0); 
  }

  gen gen::re(GIAC_CONTEXT) const {
    switch (type ) {
    case _INT_: case _DOUBLE_: case _ZINT: case _REAL: case _STRNG:
      return *this;
    case _CPLX:
      return *_CPLXptr;
    case _VECT:
      return gen(subtype==_POLY1__VECT?trim(_VECTre(*_VECTptr,contextptr),0):_VECTre(*_VECTptr,contextptr),subtype);
    case _IDNT: 
      if (!complex_variables(contextptr))
	return *this;
      if ( (_IDNTptr->value) && (is_zero(_IDNTptr->value->im(contextptr))) )
	return *this;
      return new symbolic(at_re,*this);
    case _SYMB:
      if (equalposcomp(plot_sommets,_SYMBptr->sommet))
	return new symbolic(_SYMBptr->sommet,_SYMBptr->feuille.re(contextptr));
      if (expand_re_im(contextptr))
	return symb_re(*_SYMBptr,contextptr);
      else
	return new symbolic(at_re,*this);
    case _USER:
      return _USERptr->re(contextptr);
    case _FRAC:
      return frac_reim(_FRACptr->num,_FRACptr->den,true,contextptr);
    case _MOD:
      return makemod(_MODptr->re(contextptr),*(_MODptr+1));
    case _EXT:
      return algebraic_EXTension(_EXTptr->re(contextptr),*(_EXTptr+1));
    case _POLY:
      return apply(*_POLYptr,contextptr,giac::re);
    default: 
      settypeerr("Re");
    }
    return 0;
  }
   
  vecteur _VECTim(const vecteur & a,GIAC_CONTEXT){
    vecteur res;
    vecteur::const_iterator it=a.begin(),itend=a.end();
    for (;it!=itend;++it)
      res.push_back(it->im(contextptr));
    return res;
  }

  gen symb_im(const symbolic & s,GIAC_CONTEXT){
    unary_function_ptr u=s.sommet;
    gen f=s.feuille;
    if ( (u==at_re) || (u==at_im) || (u==at_abs) )// im of a real
      return zero;
    if (u==at_conj)
      return -im(f,contextptr);
    if (u==at_plus)
      return _plus(im(f,contextptr),contextptr);
    if (u==at_neg)
      return -im(f,contextptr);
    if (u==at_pow){
      gen e=f._VECTptr->front(),expo=f._VECTptr->back();
      if (expo.type==_INT_) {
	// ? compute conj and use 1/2*(z+-zbar)?
	gen r=re(e,contextptr);
	gen i=im(e,contextptr);
	int n=f._VECTptr->back().val;
	if (n==0)
	  return zero;
	if (is_zero(i))
	  return zero;
	if (is_zero(r)){
	  if (n%2==0)
	    return zero;
	  if (n%4==1)
	    return pow(i,n);
	  else
	    return -pow(i,n);
	}
	bool n_pos=(n>0);
	if (!n_pos)
	  return im(inv(pow(f,-n),contextptr),contextptr);
	vecteur v=pascal_nth_line(n);
	vecteur somme;
	gen signe=plus_one; 
	const_iterateur it=v.begin(); // ,itend=v.end();
	++it;
	for (int j=1;j<=n;j+=2){
	  somme.push_back(signe*(*it)*pow(r,n-j)*pow(i,j));
	  ++it;
	  ++it;
	  signe=-signe;
	}
	gen res=new symbolic(at_plus,somme);
	return res;
      } // end integer exponent
      if ( is_zero(im(expo,contextptr)) && is_zero(im(e,contextptr)) ){
	// e must also be positive for non-integral power
	if (!is_integer(expo) && is_positive(-e,contextptr))
	  return pow(-e,expo,contextptr)*sin(cst_pi*expo,contextptr);
	return zero;
      }
    }
    if (u==at_rootof && f.type==_VECT && f._VECTptr->size()==2){
      vecteur tmp=*f._VECTptr;
      tmp[0]=algtrim(im(tmp[0],contextptr));
      if (tmp[0].type!=_VECT)
	return tmp[0];
      return new symbolic(u,gen(tmp,f.subtype));
    }
    if (u==at_ln)
      return arg(f,contextptr);
    gen r,i;
    symb_reim(s,r,i,contextptr);
    return i;
  }

  gen gen::im(GIAC_CONTEXT) const {
    switch (type) {
    case _INT_: case _DOUBLE_: case _ZINT: case _REAL: case _STRNG:
      return 0;
    case _CPLX:
      return *(_CPLXptr+1);
    case _VECT:
      return gen(subtype==_POLY1__VECT?trim(_VECTim(*_VECTptr,contextptr),0):_VECTim(*_VECTptr,contextptr),subtype);
    case _IDNT: 
      if (!complex_variables(contextptr))
	return zero;
      if ( (_IDNTptr->value) && (is_zero(_IDNTptr->value->im(contextptr))) )
	return zero;
      return new symbolic(at_im,*this);
    case _SYMB:      
      if (equalposcomp(plot_sommets,_SYMBptr->sommet))
	return new symbolic(_SYMBptr->sommet,_SYMBptr->feuille.im(contextptr));
      if (expand_re_im(contextptr))
	return symb_im(*_SYMBptr,contextptr); 
      else
	return new symbolic(at_im,*this);
    case _USER:
      return _USERptr->im(contextptr);
    case _FRAC:
      return frac_reim(_FRACptr->num,_FRACptr->den,false,contextptr);
    case _MOD:
      return makemod(_MODptr->im(contextptr),*(_MODptr+1));
    case _EXT:
      return algebraic_EXTension(_EXTptr->im(contextptr),*(_EXTptr+1));
    case _POLY:
      return apply(*_POLYptr,contextptr,giac::im);
    default: 
      settypeerr("Im");
    }
    return 0;
  }

  gen _VECTabs(const vecteur & a,GIAC_CONTEXT){
    gen res(0);
    vecteur::const_iterator it=a.begin(), itend=a.end();
    for (;it!=itend;++it){
      res=max(res,abs(*it,contextptr),contextptr);
    }
    return res;
  }

  gen sq(const gen & a){
    return a*a;
  }

  gen real_abs(const gen & s,GIAC_CONTEXT){
    gen tmp=evalf_double(s,1,contextptr);
    if (tmp.type==_DOUBLE_){
      if (tmp._DOUBLE_val>epsilon(contextptr))
	return s;
      if (tmp._DOUBLE_val<-epsilon(contextptr))
	return -s;
      return 0.0;
    }
    int j=sturmsign(s,contextptr);
    if (!j)
      return new symbolic(at_abs,gen(s));
    return j*s;
  }

  gen idnt_abs(const gen & s,GIAC_CONTEXT){
    if (is_inf(s))
      return plus_inf;
    if (is_undef(s))
      return s;
    if (!eval_abs(contextptr))
      return new symbolic(at_abs,s);
    gen i=im(s,contextptr);
    if (is_zero(i))
      return real_abs(s,contextptr);
    else
      return sqrt(pow(re(s,contextptr),2)+pow(i,2),contextptr);
  }

  gen symb_abs(const symbolic & s,GIAC_CONTEXT){
    unary_function_ptr u=s.sommet;
    gen f=s.feuille;
    if (u==at_abs) // abs(abs)
      return s;
    if (u==at_neg)
      return abs(f,contextptr);
    if (!complex_mode(contextptr)){ 
      if (u==at_ln)
	return real_abs(s,contextptr);
      if (!has_i(s) && (u==at_exp || u==at_sqrt))
	return s;
    }
    else {
      if (do_lnabs(contextptr) && u==at_ln)
	return new symbolic(at_abs,s);
      if (u==at_exp)
	return exp(re(f,contextptr),contextptr);
    }
    if ( (u==at_pow) && (is_zero(im(f._VECTptr->back(),contextptr))))
      return new symbolic(u,makenewvecteur(abs(f._VECTptr->front(),contextptr),f._VECTptr->back()));
    if (u==at_inv)
      return inv(abs(f,contextptr),contextptr);
    if (u==at_prod)
      return new symbolic(u,apply(f,contextptr,giac::abs));
    return idnt_abs(s,contextptr);
  }

  gen abs(const gen & a,GIAC_CONTEXT){ 
    if (is_equal(a))
      return apply_to_equal(a,abs,contextptr);
    switch (a.type ) {
    case _INT_: 
      return(absint(a.val));
    case _ZINT: 
      if (mpz_sgn(*a._ZINTptr)<0)
	return(-a);
      else
	return(a);
    case _REAL:
      return a._REALptr->abs();
    case _CPLX: 
      return sqrt(sq(*a._CPLXptr)+sq(*(a._CPLXptr+1)),contextptr) ;
    case _DOUBLE_:
      return fabs(a._DOUBLE_val);
    case _VECT:
      return _VECTabs(*a._VECTptr,contextptr);
    case _IDNT:
      return idnt_abs(a,contextptr);
    case _SYMB:
      return symb_abs(*a._SYMBptr,contextptr);
    case _USER:
      return a._USERptr->abs(contextptr);
    case _FRAC:
      return fraction(abs(a._FRACptr->num,contextptr),abs(a._FRACptr->den,contextptr));
    default:
      settypeerr("Abs");
    }
    return 0;
  }

  gen linfnorm(const gen & a,GIAC_CONTEXT){ // L^inf norm is |re|+|im| for a complex
    switch (a.type ) {
    case _INT_: 
      return(absint(a.val));
    case _ZINT: 
      if (mpz_sgn(*a._ZINTptr)<0)
	return(-a);
      else
	return(a);
    case _CPLX: 
      return(abs(*a._CPLXptr,contextptr)+abs(*(a._CPLXptr+1),contextptr)) ;
    case _DOUBLE_:
      return fabs(a._DOUBLE_val);
    case _VECT:
      return _VECTabs(*a._VECTptr,contextptr);
    case _USER:
      return a._USERptr->abs(contextptr);
    case _IDNT: case _SYMB:
      return new symbolic(at_abs,a);
    default:
      settypeerr("Linfnorm");
    }
    return 0;
  }

  gen arg_CPLX(const gen & a,GIAC_CONTEXT){
    gen realpart=normal(a.re(contextptr),contextptr),
      imaginaire=normal(a.im(contextptr),contextptr);
    if (is_zero(realpart)){
      if (is_zero(imaginaire))
	return undef;
      return cst_pi_over_2*sign(imaginaire,contextptr);
    }
    if (is_zero(imaginaire))
      return (1-sign(realpart,contextptr))*cst_pi_over_2;
    if ( (realpart.type==_DOUBLE_) || (imaginaire.type==_DOUBLE_) )
      return eval(atan(rdiv(imaginaire,realpart),contextptr)+(1-sign(realpart,contextptr))*sign(imaginaire,contextptr)*evalf_double(cst_pi_over_2,1,contextptr),1,contextptr);
    else
      return atan(rdiv(imaginaire,realpart),contextptr)+(1-sign(realpart,contextptr))*sign(imaginaire,contextptr)*cst_pi_over_2;
  }
  
  gen _VECTarg(const vecteur & a,GIAC_CONTEXT){
    vecteur res;
    vecteur::const_iterator it=a.begin(), itend=a.end();
    for (;it!=itend;++it){
      res.push_back(arg(*it,contextptr));
    }
    return res;
  }

  gen arg(const gen & a,GIAC_CONTEXT){ 
    if (is_equal(a))
      return apply_to_equal(a,arg,contextptr);
    checkanglemode(contextptr);
    switch (a.type ) {
    case _INT_: case _ZINT: case _DOUBLE_: case _REAL:
      if (is_positive(a,contextptr))
	return 0;
      else
	return cst_pi;
    case _CPLX:
      return arg_CPLX(a,contextptr);
    case _VECT:
      return _VECTarg(*a._VECTptr,contextptr);
    case _IDNT: 
    case _SYMB:
      // if ( is_zero(im(a,contextptr)) || (evalf(a,eval_level(contextptr),contextptr).type==_CPLX) )
	return arg_CPLX(a,contextptr);
	// return new symbolic(at_arg,a);
    case _USER:
      return a._USERptr->arg();
    case _FRAC:
      return arg(a._FRACptr->num,contextptr)-arg(a._FRACptr->den,contextptr);
    default:
      settypeerr("Arg");
    }
    return 0;
  }

  gen gen::squarenorm(GIAC_CONTEXT) const {
    switch (type ) {
    case _INT_: case _DOUBLE_: case _ZINT: case _REAL:
      return (*this) * (*this);
    case _CPLX:
      return ( (*_CPLXptr)*(*_CPLXptr)+(*(_CPLXptr+1)*(*(_CPLXptr+1))) );      
    default: 
      return( (*this) * this->conj(contextptr));
    }    
  }

  int gen::bindigits() const{
    int res,valeur;
    switch (type ) {
    case _INT_: 
      res=0;
      valeur=val;
      for (;valeur;res++)
	valeur = valeur >> 1;
      return res; 
    case _ZINT:
      return mpz_sizeinbase(*_ZINTptr,2)+1;
    case _CPLX: 
      return max(_CPLXptr->bindigits(),(_CPLXptr+1)->bindigits() ) ;
    default:
      settypeerr("Bindigits");
    }
    return 0;
  }

  gen gen::operator [] (int i) const{
    return operator_at(i,context0);
  }

  gen gen::operator_at(int i,GIAC_CONTEXT) const{
    if (type==_SYMB){
      if (!i)
	return _SYMBptr->sommet;
      if (_SYMBptr->feuille.type!=_VECT){
	if (i==1)
	  return _SYMBptr->feuille;
	else
	  setdimerr(contextptr);
      }
      if (unsigned(i)>_SYMBptr->feuille._VECTptr->size())
	setdimerr(contextptr);
      return (*(_SYMBptr->feuille._VECTptr))[i-1];  
    }
    if (type==_IDNT)
      return symb_at(makenewvecteur(*this,i));
    if (type==_FUNC){
      if (*this==at_ln){
	i=i+(xcas_mode(contextptr)!=0);
	return inv(ln(i,contextptr),contextptr)*(*this);
      }
      if (*this==at_maple_root){
	identificateur tmp(" x");
	gen g=symb_program(tmp,zero,new symbolic(at_makesuite,i,tmp),contextptr);
	g=makenewvecteur(at_maple_root,g);
	return symb_compose(g);
      }
    }
    if (this->type!=_VECT)
      settypeerr("Gen [int]");
    if (unsigned(i)>=_VECTptr->size()){
      i=i+(xcas_mode(contextptr)!=0);
      setdimerr("Index outside range : "+ print_INT_(i)+", vector size is "+print_INT_(_VECTptr->size())+", syntax compatibility mode "+print_program_syntax(xcas_mode(contextptr))+"\n");
    }
    return (*(this->_VECTptr))[i];
  }

  gen gen::operator [] (const gen & i) const {
    return operator_at(i,context0);
  }
  
  gen gen::operator_at(const gen & i,GIAC_CONTEXT) const {
    if (i.type==_DOUBLE_){
      double id=i._DOUBLE_val;
      if (int(id)==id)
	return (*this)[int(id)];
    }
    if (i.type==_REAL){
      double id=i.evalf_double(1,contextptr)._DOUBLE_val;
      if (int(id)==id)
	return (*this)[int(id)];
    }
    if ((type==_STRNG) && (i.type==_INT_)){
      int s=_STRNGptr->size();
      if ( (i.val<s) && (i.val>=0))
	return string2gen(string()+'"'+(*_STRNGptr)[i.val]+'"');
    }
    if (type==_IDNT)
      return new symbolic(at_at,gen(makenewvecteur(*this,i),_SEQ__VECT));
    if (type==_USER)
      return (*_USERptr)[i];
    if (type==_MAP){
      gen_map::const_iterator it=_MAPptr->find(i),itend=_MAPptr->end();
      if (it!=itend)
	return it->second;
    }
    if (is_symb_of_sommet(at_at)){ // add i at the end of the index
      if (_SYMBptr->feuille.type==_VECT && _SYMBptr->feuille._VECTptr->size()==2){
	gen operand=_SYMBptr->feuille._VECTptr->front();
	vecteur indice=makevecteur(_SYMBptr->feuille._VECTptr->back());
	indice.push_back(i);
	return symb_at(makenewvecteur(operand,gen(indice,_SEQ__VECT)));
      }
    }
    if (i.type==_DOUBLE_)
      return (*this)[(int) i._DOUBLE_val];
    if (i.type==_SYMB){
      if (i._SYMBptr->sommet==at_interval) {
	if ((i._SYMBptr->feuille._VECTptr->front().type==_INT_) && (i._SYMBptr->feuille._VECTptr->back().type==_INT_) ){
	  int debut=i._SYMBptr->feuille._VECTptr->front().val,fin=i._SYMBptr->feuille._VECTptr->back().val;
	  if (fin<debut)
	    return (type==_STRNG)?string2gen("",false):gen(vecteur(0),subtype); // swap(debut,fin);
	  if (type==_STRNG)
	    return string2gen('"'+_STRNGptr->substr(debut,fin-debut+1)+'"');
	  if (type==_VECT){
	    debut=max(debut,0);
	    fin=min(fin,_VECTptr->size()-1);
	    return gen(vecteur(_VECTptr->begin()+debut,_VECTptr->begin()+fin+1),subtype);
	  }
	}
      }
    }
    if (i.type==_VECT){
      const_iterateur it=i._VECTptr->begin(),itend=i._VECTptr->end();
      gen res (*this);
      for (;it!=itend;++it){
	if (it->type==_VECT){
	  vecteur tmp;
	  const_iterateur jt=it->_VECTptr->begin(),jtend=it->_VECTptr->end();
	  for (;jt!=jtend;++jt){
	    tmp.push_back(res[*jt]);
	  }
	  return gen(tmp,it->subtype);
	}
	if ( (it->type==_SYMB) && (it->_SYMBptr->sommet==at_interval) && (it+1!=itend) ){
	  // submatrix extraction
	  if ((it->_SYMBptr->feuille._VECTptr->front().type==_INT_) && (it->_SYMBptr->feuille._VECTptr->back().type==_INT_) ){
	    int debut=it->_SYMBptr->feuille._VECTptr->front().val,fin=it->_SYMBptr->feuille._VECTptr->back().val;
	    if (fin<debut)
	      swap(debut,fin);
	    if (res.type==_VECT){
	      debut=max(debut,0);
	      fin=min(fin,res._VECTptr->size()-1);
	      iterateur jt=res._VECTptr->begin()+debut,jtend=_VECTptr->begin()+fin+1;
	      gen fin_it(vecteur(it+1,itend),_SEQ__VECT);
	      vecteur v;
	      v.reserve(jtend-jt);
	      for (;jt!=jtend;++jt)
		v.push_back((*jt)[fin_it]);
	      return gen(v,res.subtype);
	    }
	  }
	}
	res = res[*it];
      }
      return res;
    }
    if (i.type!=_INT_)
      return symb_at(makenewvecteur(*this,i));
    return (*this)[i.val];
  }

  /*
  gen & gen::operator [](int i){
    if (this->type!=_VECT)
      settypeerr("Gen [int]");
    if (i>=_VECTptr->size())
      setdimerr(contextptr);
    return (*(this->_VECTptr))[i];
  }
  
  gen & gen::operator [] (const gen & i) {
    if (i.type==_DOUBLE_)
      return (*this)[(int) i._DOUBLE_val];
    if (i.type!=_INT_)
      settypeerr("Gen [gen]");
    if (this->type!=_VECT)
      settypeerr("Gen [gen]");
    if (i.val>=_VECTptr->size())
      setdimerr("Gen [_VECT]");
    return (*(this->_VECTptr))[i.val];
  }
  */

  bool has_inf_or_undef(const gen & g){
    if (g.type!=_VECT)
      return is_inf(g) || is_undef(g);
    const_iterateur it=g._VECTptr->begin(),itend=g._VECTptr->end();
    for (;it!=itend;++it){
      if (has_inf_or_undef(*it))
	return true;
    }
    return false;
  }

  gen gen::operator () (const gen & i,const context * contextptr) const{
    if (type==_SYMB){
      // Functional case for sommet
      if (_SYMBptr->sommet==at_program) {
	gen tmp=_SYMBptr->feuille;
	if (tmp.type!=_VECT)
	  setsizeerr(contextptr);
	(*tmp._VECTptr)[1]=i;
	return _program(tmp,undef,contextptr);
      }
      if (_SYMBptr->sommet==at_rpn_prog){
	vecteur pile;
	if (rpn_mode)
	  pile=history_out(contextptr);
	if ( (i.type!=_VECT) || (i.subtype!=_SEQ__VECT))
	  pile.push_back(i);
	else 
	  pile=mergevecteur(pile,*i._VECTptr);
	vecteur prog;
	if (_SYMBptr->feuille.type==_VECT)
	  prog=*_SYMBptr->feuille._VECTptr;
	else
	  prog=vecteur(1,_SYMBptr->feuille);
	return gen(rpn_eval(prog,pile,contextptr),_RPN_STACK__VECT);
      }
      if (_SYMBptr->sommet==at_compose){
	gen tmp=_SYMBptr->feuille;
	if (tmp.type!=_VECT)
	  return tmp(i,contextptr);
	gen res=i;
	const_iterateur it=tmp._VECTptr->begin(),itend=tmp._VECTptr->end();
	for (;itend!=it;){
	  --itend;
	  res=(*itend)(res,contextptr);
	}
	return res;
      }
      if (_SYMBptr->sommet==at_composepow){
	gen tmp=_SYMBptr->feuille;
	if (tmp.type!=_VECT || tmp._VECTptr->size()!=2 || tmp._VECTptr->back().type!=_INT_)
	  return symb_of(tmp,i);
	gen res=i;
	int n=tmp._VECTptr->back().val;
	if (n<0)
	  setsizeerr(contextptr);
	if (!n)
	  return i;
	tmp=tmp._VECTptr->front();
	for (int j=0;j<n;++j)
	  res=tmp(res,contextptr);
	return res;
      }
      if (_SYMBptr->sommet==at_function_diff || _SYMBptr->sommet==at_of || _SYMBptr->sommet==at_at)
	return new symbolic(at_of,makenewvecteur(*this,i));
      gen & f=_SYMBptr->feuille;
      if (_SYMBptr->sommet.ptr->s=="pari"){
	vecteur argv(gen2vecteur(f));
	if (i.type==_VECT && i.subtype!=_SEQ__VECT)
	  argv=mergevecteur(argv,vecteur(1,i));
	else
	  argv=mergevecteur(argv,gen2vecteur(i));
	return _SYMBptr->sommet(gen(argv,_SEQ__VECT),contextptr);
      }
      if (f==makenewvecteur(zero)){
	return _SYMBptr->sommet(i,contextptr);
      }
      // other case, apply feuille to i then apply sommet
      if (f.type!=_VECT)
	return _SYMBptr->sommet(f(i,contextptr),contextptr);
      vecteur res(*f._VECTptr);
      iterateur it=res.begin(),itend=res.end();
      for (;it!=itend;++it)
	*it=(*it)(i,contextptr);
      return _SYMBptr->sommet(res,contextptr);
    }
    if (type==_FUNC){
      if ( (i.type==_VECT) && (i.subtype==_SEQ__VECT) && (i._VECTptr->size()==1))
	return (*_FUNCptr)(i._VECTptr->front(),contextptr);
      else
	return (*_FUNCptr)(i,contextptr);
    }
    if (i.type==_DOUBLE_ && giac_floor(i._DOUBLE_val)==i._DOUBLE_val )
      return (*this)((int) i._DOUBLE_val,contextptr);
    if (type==_INT_ && subtype==_INT_TYPE && i.type==_VECT){
      return gen(*i._VECTptr,type);
    }
    if (type<_IDNT || type==_STRNG)
      return *this;
    if (type==_USER)
      return (*_USERptr)(i,contextptr);
    if (type==_VECT){
      // Old code for _VECT type was just return (*this)[i];
      vecteur w(*_VECTptr);
      iterateur it=w.begin(),itend=w.end();
      for (;it!=itend;++it)
	*it=(*it)(i,contextptr);
      return gen(w,subtype);
    }
    else {
      if (has_inf_or_undef(i))
	return undef;
      return symb_of(*this,i);
    }
  }

  bool compare_VECT(const vecteur & v,const vecteur & w){
    int s1=v.size(),s2=w.size();
    if (s1!=s2)
      return s1<s2;
    const_iterateur it=v.begin(),itend=v.end(),jt=w.begin();
    for (;it!=itend;++it,++jt){
      if (*it!=*jt)
	return it->islesscomplexthan(*jt);
    }
    // setsizeerr(); should not happen... commented because it happens!
    return false;
  }

  // return true if *this is "strictly less complex" than other
  bool gen::islesscomplexthan (const gen & other ) const {
    // FIXME it is not the natural order, but used for pivot selection
    if (type<_IDNT && is_zero(*this))
      return false;
    if (other.type<_IDNT && is_zero(other))
      return true;
    if (type != other.type)
      return type < other.type;
    if (*this==other)
      return false;
    if (type<_POLY && *this==-other)
      return is_strictly_positive(*this,context0);
    switch ( type) {
    case _INT_:
      return absint(val)<absint(other.val);
    case _ZINT:
      return is_greater(abs(other,context0),abs(*this,context0),context0);
    case _DOUBLE_: case _CPLX: case _REAL:
      return is_greater(abs(*this,context0),abs(other,context0),context0);
    case _IDNT:
      return *_IDNTptr->name < *other._IDNTptr->name;
    case _POLY:
      if (_POLYptr->coord.size()!=other._POLYptr->coord.size())
	return _POLYptr->coord.size()<other._POLYptr->coord.size();
      return _POLYptr->coord.front().value.islesscomplexthan(other._POLYptr->coord.front().value);
    case _MOD:
      if (*(_MODptr+1)!=*(other._MODptr+1))
	setsizeerr();
      return _MODptr->islesscomplexthan(*other._MODptr);
    case _SYMB:
      if (_SYMBptr->sommet !=other._SYMBptr->sommet )
	return (unsigned long) _SYMBptr->sommet.ptr <(unsigned long) other._SYMBptr->sommet.ptr;
      return _SYMBptr->feuille.islesscomplexthan(other._SYMBptr->feuille);
      // return false;
    case _VECT:
      return compare_VECT(*_VECTptr,*other._VECTptr);
    case _EXT:
      if (*(_EXTptr+1)!=*(other._EXTptr+1))
	return (_EXTptr+1)->islesscomplexthan(*(other._EXTptr+1));
      return _EXTptr->islesscomplexthan(*(other._EXTptr));
    case _STRNG:
      return *_STRNGptr<*other._STRNGptr;
    default:
      return this->print(context0)< other.print(context0); 
    }
  }

  bool islesscomplexthanf(const gen & a,const gen & b){
    return a.islesscomplexthan(b);
  }

  bool islesscomplexthanf2(const gen & a,const gen & b){
    if (a.type==_VECT && b.type==_VECT && a._VECTptr->size()==2 && b._VECTptr->size()==2){
      gen & a2=a._VECTptr->back();
      gen & b2=b._VECTptr->back();
      if (a2!=b2)
	return !a2.islesscomplexthan(b2);
    }
    return !a.islesscomplexthan(b);
  }

  int gen::symb_size () const {
    if (type==_SYMB)
      return _SYMBptr->size();
    else
      return 1;
  }

  bool symb_size_less(const gen & a,const gen & b){
    return a.symb_size() < b.symb_size();
  }

  bool gen::is_symb_of_sommet(const unary_function_ptr & u) const {
    return type==_SYMB && _SYMBptr->sommet==u;
  }

  gen operator && (const gen & a,const gen & b){
    if (is_zero(a)){
      if (b.type==_DOUBLE_)
	return 0.0;
      return a;
    }
    if (is_zero(b)){
      if (a.type==_DOUBLE_)
	return 0.0;
      return b;
    }
    if (a.is_symb_of_sommet(at_and)){
      if (b.is_symb_of_sommet(at_and))
	return new symbolic(at_and,mergevecteur(*a._SYMBptr->feuille._VECTptr,*b._SYMBptr->feuille._VECTptr));
      vecteur v=*a._SYMBptr->feuille._VECTptr;
      v.push_back(b);
      return new symbolic(at_and,v);
    }
    if (b.is_symb_of_sommet(at_and)){
      vecteur v=*b._SYMBptr->feuille._VECTptr;
      v.push_back(a);
      return new symbolic(at_and,v);
    }
    if ( ((a.type==_IDNT) || (a.type==_SYMB)) || ((a.type==_IDNT) || (a.type==_SYMB)) )
      return symb_and(a,b);
    if ( (a.type==_DOUBLE_) || (b.type==_DOUBLE_) )
      return 1.0;
    return plus_one;
  }

  gen operator || (const gen & a,const gen & b){
    if (is_zero(a))
      return b;
    if (is_zero(b))
      return a;
    if (a.is_symb_of_sommet(at_ou)){
      if (b.is_symb_of_sommet(at_ou))
	return new symbolic(at_ou,mergevecteur(*a._SYMBptr->feuille._VECTptr,*b._SYMBptr->feuille._VECTptr));
      vecteur v=*a._SYMBptr->feuille._VECTptr;
      v.push_back(b);
      return new symbolic(at_ou,v);
    }
    if (b.is_symb_of_sommet(at_ou)){
      vecteur v=*b._SYMBptr->feuille._VECTptr;
      v.push_back(a);
      return new symbolic(at_ou,v);
    }
    if ( ((a.type==_IDNT) || (a.type==_SYMB)) || ((a.type==_IDNT) || (a.type==_SYMB)) )
      return symb_ou(a,b);
    if ( (a.type==_DOUBLE_) || (b.type==_DOUBLE_) )
      return 1.0;
    return plus_one;
  }

  gen addpoly(const gen & th, const gen & other){
    if ((th.type!=_POLY) || (other.type!=_POLY))
      settypeerr("addpoly");
    // Tensor addition
    vector< monomial<gen> >::const_iterator a=th._POLYptr->coord.begin();
    vector< monomial<gen> >::const_iterator a_end=th._POLYptr->coord.end();
    if (a == a_end) {
      return other;
    }
    vector< monomial<gen> >::const_iterator b=other._POLYptr->coord.begin();
    vector< monomial<gen> >::const_iterator b_end=other._POLYptr->coord.end();
    if (b==b_end){
      return th;
    }
    polynome * resptr=new polynome(th._POLYptr->dim);
    Add<gen>(a,a_end,b,b_end,resptr->coord,th._POLYptr->is_strictly_greater);
    return resptr;
  }

  polynome addpoly(const polynome & p,const gen & c){
    if (is_exactly_zero(c))
      return p;
    polynome pcopy(p);
    if ( (!p.coord.empty()) && is_zero(*(p.coord.back().index.iptr)) ) {
      pcopy.coord.back().value = pcopy.coord.back().value + c;
      if (is_exactly_zero(pcopy.coord.back().value))
	pcopy.coord.pop_back();
    }
    else
      pcopy.coord.push_back(monomial<gen>(c,pcopy.dim));
    return pcopy;
  }

  gen chkmod(const gen& a,const gen & b){
    if  ( (b.type!=_MOD) || ((a.type==_MOD) && (*(a._MODptr+1)==*(b._MODptr+1))) )
      return a;
    return makemodquoted(a,*(b._MODptr+1));
  }
  gen makemod(const gen & a,const gen & b){
    if (a.type==_VECT)
      return apply1st(a,b,makemod);
    if (a.type==_POLY){
      polynome res(a._POLYptr->dim);
      vector< monomial<gen> >::const_iterator it=a._POLYptr->coord.begin(),itend=a._POLYptr->coord.end();
      res.coord.reserve(itend-it);
      for (;it!=itend;++it){
	gen tmp=makemod(it->value,b);
	if (!is_zero(tmp))
	  res.coord.push_back(monomial<gen>(tmp,it->index));
      }
      return res;
    }
    if (a.type==_MOD){
      if (is_zero(b)) // unmodularize
	return *a._MODptr;
      if (*(a._MODptr+1)==b) // avoid e.g. 7 % 5 % 5
	return a;
    }
    if (is_zero(b)) 
      return a;
    gen res;
    res.type=_MOD;
    res._MODptr=new gen[2];
    if ( (b.type==_INT_) || (b.type==_ZINT) )
      *res._MODptr=smod(a,b);
    else {
      if (b.type!=_VECT){
	res.type=0;
	delete [] res._MODptr;
	setsizeerr("Bad mod:"+b.print(context0));
      }
      if (a.type==_VECT)
	*res._MODptr=(*a._VECTptr)%(*b._VECTptr);
      else
	*res._MODptr=a;
    }
    *(res._MODptr+1)=b;
    res.ptr_val.ref_count=new int(1);
    return res;
  }

  gen makemodquoted(const gen & a,const gen & b){
    gen res;
    res.type=_MOD;
    res._MODptr=new gen[2];
    *res._MODptr=a;
    *(res._MODptr+1)=b;
    res.ptr_val.ref_count=new int(1);
    return res;
  }

  gen modadd(const gen * a,const gen *b){
    if (*(a+1)!=*(b+1))
      setsizeerr("Mod are different");
    return makemod(*a+*b,*(a+1));
  }

  gen modsub(const gen * a,const gen *b){
    if (*(a+1)!=*(b+1))
      setsizeerr("Mod are different");
    return makemod(*a-*b,*(a+1));
  }

  gen modmul(const gen * a,const gen *b){
    if (*(a+1)!=*(b+1))
      setsizeerr("Mod are different");
    return makemod(*a*(*b),*(a+1));
  }

  gen modinv(const gen & a){
    gen modu=*(a._MODptr+1);
    if ( ( (modu.type==_INT_) || (modu.type==_ZINT) ) && 
	 a._MODptr->is_cinteger() )
      return makemod(invmod(*a._MODptr,modu),modu);
    if (modu.type==_VECT){
      modpoly polya,u,v,d;
      if (a._MODptr->type!=_VECT)
	polya.push_back(*a._MODptr);
      else
	polya=*a._MODptr->_VECTptr;
      egcd(polya,*modu._VECTptr,0,u,v,d);
      if (d.size()!=1)
	setsizeerr("Non invertible");
      return makemod(u/d.front(),modu);
    }
    return fraction(makemod(plus_one,*(a._MODptr+1)),a);
  }

  gen real_object::operator + (const gen & g) const{
    switch (g.type){
    case _REAL:
      return *this+*g._REALptr;
    case _INT_: case _DOUBLE_: case _ZINT: case _FRAC:
#ifdef HAVE_LIBMPFR
      return *this+real_object(g,mpfr_get_prec(inf));      
#else
      return *this+real_object(g);
#endif
    default:
      return sym_add(*this,g,context0);
    }
    setsizeerr("real_object + gen"+this->print(context0)+","+g.print(context0));
  }
  
  real_interval add(const real_interval & i,const real_interval & g){
    real_interval res(i);
#ifdef HAVE_LIBMPFR
    mpfr_add(res.inf,i.inf,g.inf,GMP_RNDD);
#ifdef HAVE_LIBMPFI
    mpfi_add(res.infsup,i.infsup,g.infsup);
#else
    mpfr_add(res.sup,i.sup,g.sup,GMP_RNDU);
#endif
#else // HAVE_LIBMPFR
    mpf_add(res.inf,i.inf,g.inf);
#ifdef HAVE_LIBMPFI
    mpfi_add(res.infsup,i.infsup,g.infsup);
#else
    mpf_add(res.sup,i.sup,g.sup);
#endif
#endif // HAVE_LIBMPFR
    return res;
  }

  real_interval real_interval::operator + (const real_interval & g) const{
    return add(*this,g);
  }

  real_interval add(const real_interval & i,const real_object & g){
    const real_interval * ptr=dynamic_cast<const real_interval *>(&g);
    if (ptr)
      return add(i,*ptr);
    real_interval res(i);
#ifdef HAVE_LIBMPFR
    mpfr_add(res.inf,i.inf,g.inf,GMP_RNDD);
#ifdef HAVE_LIBMPFI
    mpfi_add_fr(res.infsup,i.infsup,g.inf);
#else
    mpfr_add(res.sup,i.sup,g.inf,GMP_RNDU);
#endif
#else // HAVE_LIBMPFR
    mpf_add(res.inf,i.inf,g.inf);
#ifdef HAVE_LIBMPFI
    mpfi_add_fr(res.infsup,i.infsup,g.inf);
#else
    mpf_add(res.sup,i.sup,g.inf);
#endif
#endif // HAVE_LIBMPFR
    return res;
  }

  real_object real_interval::operator + (const real_object & g) const{
    return add(*this,g);
  }

  real_object real_object::operator + (const real_object & g) const{
    const real_interval * ptr=dynamic_cast<const real_interval *>(&g);
    if (ptr)
      return add(*ptr,*this);
#ifdef HAVE_LIBMPFR
    mpfr_t sum;
    mpfr_init2(sum,min(mpfr_get_prec(this->inf),mpfr_get_prec(g.inf)));
    mpfr_add(sum,this->inf,g.inf,GMP_RNDN);
    real_object res(sum);
    mpfr_clear(sum);
#else
    mpf_t sum;
    mpf_init(sum);
    mpf_add(sum,this->inf,g.inf);
    real_object res(sum);
    mpf_clear(sum);
#endif
    return res;
  }

  // a and b must be dense univariate polynomials
  // WARNING: may modify a in place (suitable inside a += operator)
  gen addgen_poly(const gen & a,const gen & b,bool inplace=false){
    vecteur & av=*a._VECTptr;
    vecteur & bv=*b._VECTptr;
    if (inplace){
      /*
      int as=av.size(),bs=bv.size();
      if (as<bs)
	av.insert(av.begin(),bs-as,0);
      */
      gen res(a);
      Addmodpoly(av.begin(),av.end(),bv.begin(),bv.end(),0,av);
      return res;
    }
    gen res(vecteur(0), _POLY1__VECT);
    addmodpoly(av,bv,0,*res._VECTptr);
    return res._VECTptr->empty()?0:res;
  }

  gen operator_plus_eq(gen &a,const gen & b,GIAC_CONTEXT){
    register unsigned t=(a.type<< _DECALAGE) | b.type;
    if (t==_DOUBLE___DOUBLE_)
      return a._DOUBLE_val += b._DOUBLE_val;
    if (!t){
      longlong tmp=((longlong) a.val+b.val);
      a.val=tmp;
      if (a.val==tmp)
	return a;
      return a=tmp;
    }
    // if (!( (++control_c_counter) & control_c_counter_mask))
      control_c();
    // FIXME: move _POINTER type below _DOUBLE
    if (a.type>_DOUBLE_ && *a.ptr_val.ref_count>1)
      return a=operator_plus(a,b,contextptr);
    switch ( t ) {
    case _ZINT__ZINT:
      mpz_add(*a._ZINTptr,*a._ZINTptr,*b._ZINTptr);
      if (mpz_sizeinbase(*a._ZINTptr,2)<32){
	return a=gen(*a._ZINTptr);
      }
      return a;
    case _ZINT__INT_:
      if (b.val<0)
	mpz_sub_ui(*a._ZINTptr,*a._ZINTptr,-b.val);
      else
	mpz_add_ui(*a._ZINTptr,*a._ZINTptr,b.val);
      if (mpz_sizeinbase(*a._ZINTptr,2)<32){
	return a=gen(*a._ZINTptr);
      }
      return a;
    case _VECT__VECT:
      if (a.subtype==_POLY1__VECT){
	if (addgen_poly(a,b,true)._VECTptr->empty())
	  a=0;
	return a;
      }
    default:
      return a=operator_plus(a,b,contextptr);
    }    
  }

  gen operator_plus(const gen & a,const gen & b,unsigned t,GIAC_CONTEXT){
    // if (!( (++control_c_counter) & control_c_counter_mask))
      control_c();
    register mpz_t * e;
    switch ( t ) {
    case _ZINT__ZINT:
      e =(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init(*e);
      mpz_add(*e,*a._ZINTptr,*b._ZINTptr);
      return(e);
    case _DOUBLE___DOUBLE_:
      return a._DOUBLE_val+b._DOUBLE_val;
    case _VECT__VECT:
      if (a.subtype==_POLY1__VECT)
	return addgen_poly(a,b);
      if (a.subtype==_PNT__VECT)
	return gen(makenewvecteur(a._VECTptr->front()+b,a._VECTptr->back()),a.subtype);
      if (a.subtype!=_POINT__VECT && equalposcomp(_GROUP__VECT_subtype,a.subtype))
	return sym_add(a,b,contextptr);
      if (b.subtype!=_POINT__VECT && equalposcomp(_GROUP__VECT_subtype,b.subtype))
	return sym_add(b,a,contextptr);
      if (a.subtype==_POINT__VECT && b.subtype==_POINT__VECT)
	return gen(addvecteur(*a._VECTptr,*b._VECTptr),0);	
      return gen(addvecteur(*a._VECTptr,*b._VECTptr),a.subtype?a.subtype:b.subtype);
    case _INT___ZINT: 
      e =(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init(*e);
      if (a.val<0)
	mpz_sub_ui(*e,*b._ZINTptr,-a.val);
      else
	mpz_add_ui(*e,*b._ZINTptr,a.val);
      return(e);
    case _ZINT__INT_:
      e =(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init(*e);
      if (b.val<0)
	mpz_sub_ui(*e,*a._ZINTptr,-b.val);
      else
	mpz_add_ui(*e,*a._ZINTptr,b.val);
      return(e);
    case _DOUBLE___INT_:
      return a._DOUBLE_val+b.val;
    case _INT___DOUBLE_:
      return a.val+b._DOUBLE_val;
    case _DOUBLE___ZINT:
      return a._DOUBLE_val+mpz_get_d(*b._ZINTptr);
    case _DOUBLE___REAL:
      return a._DOUBLE_val+real2double(*b._REALptr);
    case _REAL__DOUBLE_:
      return b._DOUBLE_val+real2double(*a._REALptr);
    case _ZINT__DOUBLE_:
      return b._DOUBLE_val+mpz_get_d(*a._ZINTptr);
    case _CPLX__INT_: case _CPLX__ZINT: case _CPLX__DOUBLE_: case _CPLX__REAL:
      return gen(*a._CPLXptr+b,*(a._CPLXptr+1));
    case _INT___CPLX: case _ZINT__CPLX: case _DOUBLE___CPLX: case _REAL__CPLX:
      return gen(a+*b._CPLXptr,*(b._CPLXptr+1));
    case _CPLX__CPLX:
      return gen(*a._CPLXptr + *b._CPLXptr, *(a._CPLXptr+1) + *(b._CPLXptr+1));
    case _POLY__POLY:
      return addpoly(a,b);
    case _FRAC__FRAC:
      return (*a._FRACptr)+(*b._FRACptr);
    case _SPOL1__SPOL1:
      return spadd(*a._SPOL1ptr,*b._SPOL1ptr,contextptr);
    case _EXT__EXT:
      return ext_add(a,b,contextptr);
    case _STRNG__STRNG:
      return string2gen('"'+(*a._STRNGptr)+(*b._STRNGptr)+'"');
    case _POLY__INT_: case _POLY__ZINT: case _POLY__DOUBLE_: case _POLY__CPLX: case _POLY__MOD: case _POLY__USER: case _POLY__REAL: 
      return addpoly(*a._POLYptr,b);
    case _INT___POLY: case _ZINT__POLY: case _DOUBLE___POLY: case _CPLX__POLY: case _MOD__POLY: case _USER__POLY: case _REAL__POLY: 
      return addpoly(*b._POLYptr,a);
    case _MOD__MOD:
      return modadd(a._MODptr,b._MODptr);
    case _REAL__REAL:
      return (*a._REALptr)+(*b._REALptr);
    case _IDNT__IDNT:
      if (a==undef || a==unsigned_inf)
	return a;
      if (b==undef || b==unsigned_inf)
	return b;
      return new symbolic(at_plus,makenewvecteur(a,b));
    default:
      if (a.type==_STRNG)
	return string2gen(*a._STRNGptr+b.print(contextptr),false);
      if (b.type==_STRNG)
	return string2gen(a.print(contextptr)+*b._STRNGptr,false);
      if (a.type==_USER)
	return (*a._USERptr)+b;
      if (b.type==_USER)
	return (*b._USERptr)+a;
      if (a.type==_REAL)
	return (*a._REALptr)+b;
      if (b.type==_REAL){
	return (*b._REALptr)+a;
      }
      return sym_add(a,b,contextptr);
    }
  }

  gen operator_plus (const gen & a,const gen & b,GIAC_CONTEXT){
    register unsigned t=(a.type<< _DECALAGE) | b.type;
    if (!t)
      return((longlong) a.val+b.val);
    return operator_plus(a,b,t,contextptr);
  }

  gen operator + (const gen & a,const gen & b){
    register unsigned t=(a.type<< _DECALAGE) | b.type;
    if (!t)
      return ((longlong) a.val+b.val);
    return operator_plus(a,b,t,context0);
  }
  
  // specialization of Tfraction<gen> operator +
  Tfraction<gen> operator + (const Tfraction<gen> & a,const Tfraction<gen> &b){
    if (is_one(a.den))
      return(Tfraction<gen> (a.num+b));
    if (is_one(b.den))
      return(Tfraction<gen> (b.num+a));
    gen da(a.den),db(b.den);
    gen den=simplify3(da,db),num;
    if (a.num.type==_POLY && b.num.type==_POLY && db.type==_POLY && da.type==_POLY)
      num=foisplus(*a.num._POLYptr,*db._POLYptr,*b.num._POLYptr,*da._POLYptr);
    else
      num=(a.num*db+b.num*da);
    if (is_exactly_zero(num))
      return Tfraction<gen>(num,1);
    simplify3(num,den);
    den=den*da*db;
    return Tfraction<gen> (num,den);
  }
  

  symbolic symbolic_plot_makevecteur(const unary_function_ptr & u,const gen & e,bool project,GIAC_CONTEXT){
    if ( (u!=at_pnt) || (e.type!=_VECT) || (e.subtype!=_PNT__VECT) )
      return symbolic(u,e);
    // e is a curve or a pnt
    vecteur w(*e._VECTptr);
    if ( (w.size()!=2) && (w.size()!=3))
      return symbolic(u,e);
    gen a0(w[0]);
    gen a1(w[1]);
    if ( a1.type==_VECT && a1._VECTptr->size()==2 ){
      if (project){
	// we must project a0
	gen param=a1._VECTptr->back(); // v= [ pnt() t ]
	if (param.type!=_VECT) setsizeerr("gen.cc/symbolic_plot_makevecteur");
	vecteur v=*param._VECTptr;
	v[1]=projection(v[0],a0,contextptr); 
	a0=remove_at_pnt(parameter2point(v,contextptr)); // same
	a1=makenewvecteur(a1._VECTptr->front(),v);
      }
      else
	a1=a1._VECTptr->front();
    }
    return symbolic(u,gen(makenewvecteur(a0,a1),_PNT__VECT));
  }

  gen collect(const gen & g,GIAC_CONTEXT){
    if (g.type==_VECT)
      return apply(g,collect,contextptr);
    return liste2symbolique(symbolique2liste(g,contextptr));
  }

  bool modified_islesscomplexthanf(const gen& a,const gen& b){
    if (a.is_symb_of_sommet(at_neg))
      return modified_islesscomplexthanf(a._SYMBptr->feuille,b);
    if (b.is_symb_of_sommet(at_neg))
      return modified_islesscomplexthanf(a,b._SYMBptr->feuille);
    if (a.is_symb_of_sommet(at_pow))
      return modified_islesscomplexthanf(a._SYMBptr->feuille[0],b);
    if (b.is_symb_of_sommet(at_pow))
      return modified_islesscomplexthanf(a,b._SYMBptr->feuille[0]);
    return islesscomplexthanf(a,b);
  }

  // return true if a is -basis^exp
  bool power_basis_exp(const gen& a,gen & basis,gen & expa){
    if (a.is_symb_of_sommet(at_neg))
      return !power_basis_exp(a._SYMBptr->feuille,basis,expa);
    if (a.is_symb_of_sommet(at_pow)){
      gen & tmp=a._SYMBptr->feuille;
      if (tmp.type!=_VECT || tmp._VECTptr->size()!=2)
	setsizeerr();
      expa=tmp._VECTptr->back();
      basis=tmp._VECTptr->front();
    }
    else {
      basis=a;
      expa=plus_one;
    }
    return false;
  }

  // Helpers for symbolic addition
  // from a product returns a list with the numeric coeff and the monomial
  vecteur terme2unitaire(const gen & x,bool sorted,GIAC_CONTEXT){
    if (x.type<_POLY)
      return makevecteur(x,plus_one);
    gen tmp;
    if (x.type!=_SYMB)
      return makevecteur(1,x);
    if (x._SYMBptr->sommet==at_neg){
      vecteur v=terme2unitaire(x._SYMBptr->feuille,sorted,contextptr);
      v[0]=-v[0];
      return v;
    }
    if (x._SYMBptr->sommet==at_prod && (tmp=x._SYMBptr->feuille).type==_VECT && !tmp._VECTptr->empty() ){
      vecteur & v = *tmp._VECTptr;
      int s=v.size();
      if (s==2 && (sorted || tmp.subtype==_SORTED__VECT))
	return makevecteur(v[0],v[1]);
      vecteur vtmp(v.begin(),v.end());
      sort(vtmp.begin(),vtmp.end(),modified_islesscomplexthanf);
      // collect term with the same power
      const_iterateur it=vtmp.begin(),itend=vtmp.end();
      vecteur vsorted;
      vsorted.reserve(itend-it);
      gen precbasis,precexpo,basis,expo,constcoeff(plus_one);
      for (;it!=itend;++it){
	if (it->type<=_CPLX)
	  constcoeff=constcoeff*(*it);
	else {
	  if (it->is_symb_of_sommet(at_inv) && it->_SYMBptr->feuille.type<=_CPLX)
	    constcoeff=constcoeff/it->_SYMBptr->feuille;
	  else
	    break;
	}
      }
      if (!is_one(constcoeff))
	vsorted.push_back(constcoeff);
      bool isneg(false),hasneg(false);
      if (it!=itend){
	power_basis_exp(*it,precbasis,precexpo);
	precexpo=zero;
	for (;it!=itend;++it){
	  if (power_basis_exp(*it,basis,expo)){
	    isneg=!isneg;
	    hasneg=true;
	  }
	  if (basis==precbasis)
	    precexpo=precexpo+expo;
	  else {
	    vsorted.push_back(pow(precbasis,precexpo,contextptr));
	    precbasis=basis;
	    precexpo=expo;
	  }
	}
	vsorted.push_back(pow(precbasis,precexpo,contextptr));
      }
      vecteur res;
      if (hasneg){
	res=terme2unitaire(_prod(vsorted,contextptr),sorted,contextptr);
	if (isneg)
	  res[0]=-res[0];
	return res;
      }
      if (vsorted.empty())
	vsorted.push_back(1);
      if (vsorted.front().type<_POLY){
	vtmp=vecteur(vsorted.begin()+1,vsorted.end());
	gen tt(1);
	if (!vtmp.empty()){
	  if (vtmp.size()==1)
	    tt=vtmp.front();
	  else
	    tt=new symbolic(at_prod,gen(vtmp,_SORTED__VECT));
	}
	res=makevecteur(vsorted.front(),tt);
      }
      else
	res=makevecteur(plus_one,new symbolic(at_prod,gen(vsorted,_SORTED__VECT)));
      return res;
    }
    // recurse
    return makevecteur(plus_one,new symbolic(x._SYMBptr->sommet,collect(x._SYMBptr->feuille,contextptr)));
  }

  // assumes v is a sorted list, shrink it
  // should be written to a gen of type _VECT and subtype _SORTED__VECT
  vecteur fusionliste(const vecteur & v){
    const_iterateur it=v.begin(),itend=v.end();
    if (itend-it<2)
      return v;
    vecteur res;
    gen current=(*it)[1];
    gen current_coeff=(*it)[0];
    ++it;
    for (;it!=itend;++it){
      if ((*it)[1]!=current){
	res.push_back(makenewvecteur(current_coeff,current));
	current_coeff=(*it)[0];
	current=(*it)[1];
      }
      else
	current_coeff=current_coeff+(*it)[0];
    }
    res.push_back(makenewvecteur(current_coeff,current));
    return res;
  }

  // from a sum in x returns a list of [coeff monomial]
  // e.g. 5+2x+3*x*y -> [ [5 1] [2 x] [ 3 x*y] ]
  vecteur symbolique2liste(const gen & x,GIAC_CONTEXT){
    if (!x.is_symb_of_sommet(at_plus))
      return vecteur(1,terme2unitaire(x,false,contextptr));
    bool sorted=x._SYMBptr->feuille.subtype==_SORTED__VECT;
    gen number;
    vecteur varg=gen2vecteur(x._SYMBptr->feuille);
    vecteur vres;
    const_iterateur it=varg.begin(),itend=varg.end();    
    for (;it!=itend;++it){
      if (it->type<_POLY)
	number=number+(*it);
      else
	vres.push_back(terme2unitaire(*it,sorted,contextptr));
    }
    if (!is_exactly_zero(number))
      vres.push_back(makenewvecteur(1,number));
    if (x._SYMBptr->feuille.subtype==_SORTED__VECT)
      return vres;
    sort(vres.begin(),vres.end(),islesscomplexthanf2);
    return fusionliste(vres);
  }

  // assumes v1, v2 are sorted and shrinked, merge them -> sorted and shrinked
  vecteur fusion2liste(const vecteur & v1,const vecteur & v2){
    const_iterateur it=v1.begin(),itend=v1.end(),jt=v2.begin(),jtend=v2.end();
    vecteur res;
    gen tmp;
    for (;it!=itend;){
      if (jt==jtend){
	for (;it!=itend;++it)
	  res.push_back(*it);
	return res;
      }
      // both iterator are valid
      vecteur & vi=*it->_VECTptr;
      vecteur & vj=*jt->_VECTptr;      
      if (vi[1]==vj[1]){
	tmp=vi[0]+vj[0];
	if (!is_exactly_zero(tmp))
	  res.push_back(makenewvecteur(tmp,vi[1]));
	++it;
	++jt;
      }
      else {
	if (vi[1].islesscomplexthan(vj[1])){
	  res.push_back(*it);
	  ++it;
	}
	else {
	  res.push_back(*jt);
	  ++jt;
	}
      } // end tests
    } // end for loop
    // finish jt
    for (;jt!=jtend;++jt)
      res.push_back(*jt);
    return res;
  }

  // v should be sorted and shrinked
  gen liste2symbolique(const vecteur & v){
    vecteur res;
    const_iterateur it=v.begin(),itend=v.end();
    res.reserve(itend-it);
    for (;it!=itend;++it){
      vecteur & vtmp(*it->_VECTptr);
      gen & tmp = vtmp.back();
      if (tmp.is_symb_of_sommet(at_prod) && tmp._SYMBptr->feuille.type==_VECT && tmp._SYMBptr->feuille._VECTptr->size()==1)
	res.push_back(vtmp.front()*tmp._SYMBptr->feuille._VECTptr->front());
      else
	res.push_back(vtmp.front()*tmp);
    }
    int s=res.size();
    if (!s)
      return zero;
    if (s==1)
      return res.front();
    return new symbolic(at_plus,gen(res,_SORTED__VECT));
  }

  gen sym_add(const gen & a,const gen & b,GIAC_CONTEXT){
    control_c();
    if (a.is_symb_of_sommet(at_unit)){
      if (b.is_symb_of_sommet(at_unit)){
	vecteur & va=*a._SYMBptr->feuille._VECTptr;
	vecteur & vb=*b._SYMBptr->feuille._VECTptr;
	if (va[1]==vb[1])
	  return new symbolic(at_unit,makenewvecteur(va[0]+vb[0],va[1]));
	gen g=mksa_reduce(vb[1]/va[1],contextptr);
	chk_not_unit(g);
	return new symbolic(at_unit,makenewvecteur(va[0]+g*vb[0],va[1]));
      }
      gen g=mksa_reduce(a,contextptr);
      chk_not_unit(g);
      return g+b;
    }
    if (b.is_symb_of_sommet(at_unit)){
      gen g=mksa_reduce(b,contextptr);
      chk_not_unit(g);
      return a+g;
    }
    if (a.is_approx()){
      gen b1;
      if (has_evalf(b,b1,1,contextptr) && b!=b1)
	return a+b1;
    }
    if (b.is_approx()){
      gen a1;
      if (has_evalf(a,a1,1,contextptr) && a!=a1)
	return a1+b;
    }
    if ( (a.type==_SYMB) && equalposcomp(plot_sommets,a._SYMBptr->sommet) ){
      if ( (b.type==_SYMB) && equalposcomp(plot_sommets,b._SYMBptr->sommet) )
	return a._SYMBptr->feuille._VECTptr->front()+b._SYMBptr->feuille._VECTptr->front();
      else {
	if (b.type==_VECT)
	  return translation(b,a,contextptr);
	return symbolic_plot_makevecteur( a._SYMBptr->sommet,a._SYMBptr->feuille+b,true,contextptr);
      }
    }
    if ( (b.type==_SYMB) && equalposcomp(plot_sommets,b._SYMBptr->sommet) ){
      if (a.type==_VECT)
	return translation(a,b,contextptr);
      return symbolic_plot_makevecteur(b._SYMBptr->sommet,b._SYMBptr->feuille+a,true,contextptr);
    }
    if (a.type==_VECT){
      if (is_zero(b))
	return a;
      vecteur res=*a._VECTptr;
      if (res.empty())
	return b;
      if (a.subtype==_VECTOR__VECT && a._VECTptr->size()==2){ 
	if (b.type==_VECT && b._VECTptr->size()==2){
	  vecteur & bv=*b._VECTptr;
	  if (b.subtype==_VECTOR__VECT && res.front()==bv.back())
	    return _vector(gen(makenewvecteur(bv.front(),bv.back()+res.back()-res.front()),_SEQ__VECT),contextptr);
	  return _vector(gen(makenewvecteur(res.front(),res.back()+bv.back()-bv.front()),_SEQ__VECT),contextptr);
	}
	return _point(b+res.back()-res.front(),contextptr);
      }
      if (b.type==_VECT && b.subtype==_VECTOR__VECT && b._VECTptr->size()==2)
	return a+vector2vecteur(*b._VECTptr);
      if (equalposcomp(_GROUP__VECT_subtype,a.subtype)){ // add to each element
	iterateur it=res.begin(),itend=res.end();
	for (;it!=itend;++it)
	  *it=*it+b;
	return gen(res,a.subtype);
      }
      if (a.subtype==_PNT__VECT){
	res.front()=res.front()+b;
	return gen(res,_PNT__VECT);
      }
      if (a.subtype!=_POLY1__VECT && ckmatrix(a)){ // matrix+cst
	int s=res.size();
	if (unsigned(s)==res.front()._VECTptr->size()){
	  for (int i=0;i<s;i++)
	    (*(res[i]._VECTptr))[i] = (*(res[i]._VECTptr))[i] + b;
	  return res;
	}
      }
      // polynomial+cst
      res.back()=res.back()+b;
      if ( (res.size()==1) && is_exactly_zero(res.back()))
          return zero;
      else
          return gen(res,_POLY1__VECT);
    }
    if (b.type==_VECT)
      return sym_add(b,a,contextptr);
    if ((a==undef) || (b==undef))
      return undef;
    if (is_inf(a)){
      if (is_inf(b)){
	if ((a==b) && (a!=unsigned_inf))
	  return a;
	else
	  return undef;
      }
      else
	return a;
    }
    if (is_inf(b))
      return b;
    if (b.is_symb_of_sommet(at_neg) && a==b._SYMBptr->feuille)
      return chkmod(zero,a);
    if (a.is_symb_of_sommet(at_neg) && b==a._SYMBptr->feuille)
      return chkmod(zero,b);
    if (is_exactly_zero(a))
      return b;
    if (is_exactly_zero(b))
      return a;
    if (a.type==_STRNG)
      return string2gen(*a._STRNGptr+b.print(context0),false);
    if (b.type==_STRNG)
      return string2gen(a.print(context0)+*b._STRNGptr,false);
    if (a.type==_FRAC){
      if ( (b.type!=_SYMB) && (b.type!=_IDNT) )
	return (*a._FRACptr)+b;
      return sym_add(_FRAC2_SYMB(a),b,contextptr);
    }
    if (b.type==_FRAC){
      if ( (a.type!=_SYMB) && (a.type!=_IDNT) )
	return a+(*b._FRACptr);
      return sym_add(a,_FRAC2_SYMB(b),contextptr);
    }
    if (a.type==_EXT){
      if (a.is_constant() && (b.type==_POLY))
	return addpoly(*b._POLYptr,a);
      /*
      if (b.type==_POLY && b.is_constant())
	return a+b._POLYptr->coord.front().value;
      */
      else
	return algebraic_EXTension(*a._EXTptr+b,*(a._EXTptr+1));
    }
    if (b.type==_EXT){
      if (b.is_constant() && (a.type==_POLY))
	return addpoly(*a._POLYptr,b);
      /*
      if (a.type==_POLY && a.is_constant())
	return a._POLYptr->coord.front().value+b;
      */
      else
	return algebraic_EXTension(a+*b._EXTptr,*(b._EXTptr+1));
    }
    int ia=is_inequality(a),ib=is_inequality(b);
    if (ia){
      vecteur & va=*a._SYMBptr->feuille._VECTptr;
      if (ia==ib || (ia==1 && ib)){
	if (ia==4) // <> + <>
	  return undef;
	vecteur & vb=*b._SYMBptr->feuille._VECTptr;
	return new symbolic(b._SYMBptr->sommet,makenewvecteur(va.front()+vb.front(),va.back()+vb.back()));
      }
      if (ia==1 || !ib) // = + 
	return new symbolic(a._SYMBptr->sommet,makenewvecteur(va.front()+b,va.back()+b));
      if ( (ia==5 && ib==6) || (ia==6 && ib==5)){
	vecteur & vb=*b._SYMBptr->feuille._VECTptr;
	return new symbolic(at_superieur_strict,makenewvecteur(va.front()+vb.front(),va.back()+vb.back()));
      }
    }
    if (ib)
      return b+a;
    if (a.is_symb_of_sommet(at_interval)){
      gen & f=a._SYMBptr->feuille;
      if (f.type==_VECT && f._VECTptr->size()==2){
	vecteur & v=*f._VECTptr;
	if (b.is_symb_of_sommet(at_interval)){
	  gen & g=b._SYMBptr->feuille;
	  if (g.type==_VECT && g._VECTptr->size()==2){
	    vecteur & w=*g._VECTptr;
	    return new symbolic(at_interval,gen(makenewvecteur(w[0]+v[0],w[1]+v[1]),_SEQ__VECT));
	  }
	}
	return new symbolic(at_interval,gen(makenewvecteur(b+v[0],b+v[1]),_SEQ__VECT));
      }
    }
    if (b.is_symb_of_sommet(at_interval))
      return b+a;
    /* if (xcas_mode(contextptr) && (a.type==_SYMB|| b.type==_SYMB) )
       return liste2symbolique(fusion2liste(symbolique2liste(a),symbolique2liste(b))); */
    if ((a.type==_SYMB) && (b.type==_SYMB)){
      if (a._SYMBptr->sommet==at_plus) {
	if (b._SYMBptr->sommet==at_plus)
	  return new symbolic(at_plus,mergevecteur(*(a._SYMBptr->feuille._VECTptr),*(b._SYMBptr->feuille._VECTptr)));
	else
	  return new symbolic(*a._SYMBptr,b);
      }
      else { 
	if (b._SYMBptr->sommet==at_plus)
	  return new symbolic(*(b._SYMBptr),a);
	else
	  return new symbolic(at_plus,gen(makenewvecteur(a,b)));
      }
    }
    if (b.type==_SYMB){
      if (b._SYMBptr->sommet==at_plus)
	return new symbolic(a,b._SYMBptr->sommet,b._SYMBptr->feuille); // new symbolic(*b._SYMBptr,a);
      else
	return new symbolic(at_plus,gen(makenewvecteur(a,b)));
    }
    if (a.type==_SYMB){
      if (b==plus_one && a._SYMBptr->sommet==at_plus && a._SYMBptr->feuille.type==_VECT && a._SYMBptr->feuille._VECTptr->size()>1 && a._SYMBptr->feuille._VECTptr->back()==minus_one){
	vecteur v=*a._SYMBptr->feuille._VECTptr;
	v.pop_back();
	if (v.size()==1)
	  return v.front();
	else
	  return new symbolic(at_plus,gen(v,a._SYMBptr->feuille.subtype));
      }
      if (a._SYMBptr->sommet==at_plus)
	return new symbolic(*a._SYMBptr,b);
      else
	return new symbolic(at_plus,gen(makenewvecteur(a,b)));
    }
    if ( (a.type==_IDNT) || (b.type==_IDNT))
      return new symbolic(at_plus,gen(makenewvecteur(a,b)));
    if (a.type==_MOD)
      return a+makemod(b,*(a._MODptr+1));
    if (b.type==_MOD)
      return makemod(a,*(b._MODptr+1))+b;
    return new symbolic(at_plus,makenewvecteur(a,b));
    // settypeerr("sym_add");
  }

  gen subpoly(const gen & th, const gen & other){
    if ((th.type!=_POLY) || (other.type!=_POLY))
      settypeerr("subpoly");
    vector< monomial<gen> >::const_iterator a=th._POLYptr->coord.begin();
    vector< monomial<gen> >::const_iterator a_end=th._POLYptr->coord.end();
    vector< monomial<gen> >::const_iterator b=other._POLYptr->coord.begin();
    vector< monomial<gen> >::const_iterator b_end=other._POLYptr->coord.end();
    if (b==b_end){
      return th;
    }
    polynome * resptr=new polynome(th._POLYptr->dim);
    Sub<gen>(a,a_end,b,b_end,resptr->coord,th._POLYptr->is_strictly_greater);
    return resptr;
  }

  polynome subpoly(const polynome & p,const gen & c){
    if (is_exactly_zero(c))
      return p;
    polynome pcopy(p);
    if ( (!p.coord.empty()) && is_zero(*(p.coord.back().index.iptr)) ) {
      pcopy.coord.back().value = pcopy.coord.back().value - c;
      if (is_exactly_zero(pcopy.coord.back().value))
	pcopy.coord.pop_back();
    }
    else
      pcopy.coord.push_back(monomial<gen>(-c,pcopy.dim));
    return pcopy;
  }

  polynome subpoly(const gen & c,const polynome & p){
    if (is_exactly_zero(c))
      return -p;
    polynome pcopy(-p);
    if ( (!p.coord.empty()) && is_zero(*(p.coord.back().index.iptr)) ) {
      pcopy.coord.back().value = pcopy.coord.back().value + c;
      if (is_exactly_zero(pcopy.coord.back().value))
	pcopy.coord.pop_back();
    }
    else
      pcopy.coord.push_back(monomial<gen>(c,pcopy.dim));
    return pcopy;
  }

  gen real_object::operator - (const gen & g) const{
    switch (g.type){
    case _REAL:
      return *this-*g._REALptr;
    case _INT_: case _DOUBLE_: case _ZINT: case _FRAC:
#ifdef HAVE_LIBMPFR
      return *this - real_object(g,mpfr_get_prec(inf));      
#else
      return *this - real_object(g);
#endif
    default:
      return sym_sub(*this,g,context0);
    }
    setsizeerr("real_object + gen"+this->print(context0)+","+g.print(context0));
  }
  
  real_interval sub(const real_interval & i,const real_interval & g){
    real_interval res(i);
#ifdef HAVE_LIBMPFI
    mpfi_sub(res.infsup,i.infsup,g.infsup);
    mpfr_sub(res.inf,i.inf,g.inf,GMP_RNDD);
#else
#ifdef HAVE_LIBMPFR
    mpfr_sub(res.inf,i.sup,g.inf,GMP_RNDD);
    mpfr_sub(res.sup,i.inf,g.sup,GMP_RNDU);    
#else
    mpf_sub(res.inf,i.sup,g.inf);
    mpf_sub(res.sup,i.inf,g.sup);    
#endif
#endif
    return res;
  }

  real_interval real_interval::operator - (const real_interval & g) const{
    return sub(*this,g);
  }

  real_interval sub(const real_interval & i,const real_object & g){
    const real_interval * ptr=dynamic_cast<const real_interval *>(&g);
    if (ptr)
      return sub(i,*ptr);
    real_interval res(i);
#ifdef HAVE_LIBMPFI
    mpfi_sub_fr(res.infsup,i.infsup,g.inf);
    mpfr_sub(res.inf,i.inf,g.inf,GMP_RNDD);
#else
#ifdef HAVE_LIBMPFR
    mpfr_sub(res.inf,i.sup,g.inf,GMP_RNDD);
    mpfr_sub(res.sup,i.inf,g.inf,GMP_RNDU);    
#else
    mpf_sub(res.inf,i.sup,g.inf);
    mpf_sub(res.sup,i.inf,g.inf);    
#endif
#endif
    return res;
  }

  real_object real_interval::operator - (const real_object & g) const{
    return sub(*this,g);
  }

  real_object real_object::operator - (const real_object & g) const{
    const real_interval * ptr=dynamic_cast<const real_interval *>(&g);
    if (ptr)
      return add(-*ptr,*this);
#ifdef HAVE_LIBMPFR
    mpfr_t sum;
    mpfr_init2(sum,min(mpfr_get_prec(this->inf),mpfr_get_prec(g.inf)));
    mpfr_sub(sum,this->inf,g.inf,GMP_RNDN);
    real_object res(sum);
    mpfr_clear(sum);
#else
    mpf_t sum;
    mpf_init(sum);
    mpf_sub(sum,this->inf,g.inf);
    real_object res(sum);
    mpf_clear(sum);
#endif
    return res;
  }

  gen subgen_poly(const gen & a,const gen & b,bool inplace=false){
    vecteur & av=*a._VECTptr;
    vecteur & bv=*b._VECTptr;
    if (inplace){
      /*
      int as=av.size(),bs=bv.size();
      if (as<bs)
	av.insert(av.begin(),bs-as,0);
      */
      gen res(a);
      Submodpoly(av.begin(),av.end(),bv.begin(),bv.end(),0,av);
      return res;
    }
    gen res(vecteur(0), _POLY1__VECT);
    submodpoly(av,bv,0,*res._VECTptr);
    return res._VECTptr->empty()?0:res;
  }

  gen operator_minus_eq (gen & a,const gen & b,GIAC_CONTEXT){
    register unsigned t=(a.type<< _DECALAGE) | b.type;
    if (t==_DOUBLE___DOUBLE_)
      return a._DOUBLE_val -= b._DOUBLE_val;
    if (!t){
      longlong tmp=((longlong) a.val-b.val);
      a.val=tmp;
      if (a.val==tmp)
	return a;
      return a=tmp;
    }
    // if (!( (++control_c_counter) & control_c_counter_mask))
      control_c();
    // FIXME: move _POINTER type below _DOUBLE
    if (a.type>_DOUBLE_ && *a.ptr_val.ref_count>1)
      return a=operator_minus(a,b,contextptr);
    switch ( t ) {
    case _ZINT__ZINT:
      mpz_sub(*a._ZINTptr,*a._ZINTptr,*b._ZINTptr);
      if (mpz_sizeinbase(*a._ZINTptr,2)<32){
	return a=gen(*a._ZINTptr);
      }
      return a;
    case _ZINT__INT_:
      if (b.val>0)
	mpz_sub_ui(*a._ZINTptr,*a._ZINTptr,b.val);
      else
	mpz_add_ui(*a._ZINTptr,*a._ZINTptr,-b.val);
      if (mpz_sizeinbase(*a._ZINTptr,2)<32){
	return a=gen(*a._ZINTptr);
      }
      return a;
    case _VECT__VECT:
      if (a.subtype==_POLY1__VECT){
	if (subgen_poly(a,b,true)._VECTptr->empty())
	  a=0;
	return a;
      }
    default:
      return a=operator_minus(a,b,contextptr);
    }    
  }

  gen operator_minus (const gen & a,const gen & b,unsigned t,GIAC_CONTEXT){
    // if (!( (++control_c_counter) & control_c_counter_mask))
      control_c();
    register mpz_t * e;
    switch ( t) {
    case _ZINT__ZINT:
      e =(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init(*e);
      mpz_sub(*e,*a._ZINTptr,*b._ZINTptr);
      return e;
    case _DOUBLE___DOUBLE_:
      return a._DOUBLE_val-b._DOUBLE_val;
    case _VECT__VECT:
      if (a.subtype==_POLY1__VECT)
	return subgen_poly(a,b);
      if (a.subtype==_PNT__VECT)
	return gen(makenewvecteur(a._VECTptr->front()-b,a._VECTptr->back()),a.subtype);
      if (a.subtype!=_POINT__VECT && equalposcomp(_GROUP__VECT_subtype,a.subtype))
	return sym_sub(a,b,contextptr);
      if (a.subtype==_POINT__VECT && b.subtype==_POINT__VECT)
	return gen(subvecteur(*a._VECTptr,*b._VECTptr),0);
      return gen(subvecteur(*a._VECTptr,*b._VECTptr),a.subtype);
    case _INT___ZINT: 
      e =(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init(*e);
      if (a.val<0)
	mpz_add_ui(*e,*b._ZINTptr,-a.val);
      else
	mpz_sub_ui(*e,*b._ZINTptr,a.val);
      mpz_neg(*e,*e);
      return(e);
    case _ZINT__INT_:
      e =(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init(*e);
      if (b.val<0)
	mpz_add_ui(*e,*a._ZINTptr,-b.val);
      else
	mpz_sub_ui(*e,*a._ZINTptr,b.val);
      return(e);
    case _INT___DOUBLE_:
      return a.val-b._DOUBLE_val;
    case _DOUBLE___INT_:
      return a._DOUBLE_val-b.val;
    case _ZINT__DOUBLE_:
      return mpz_get_d(*a._ZINTptr)-b._DOUBLE_val;
    case _DOUBLE___ZINT:
      return a._DOUBLE_val-mpz_get_d(*b._ZINTptr);
    case _DOUBLE___REAL:
      return a._DOUBLE_val-real2double(*b._REALptr);
    case _REAL__DOUBLE_:
      return real2double(*a._REALptr)-b._DOUBLE_val;
    case _CPLX__INT_: case _CPLX__ZINT: case _CPLX__DOUBLE_: case _CPLX__REAL:
      return gen(*a._CPLXptr-b,*(a._CPLXptr+1));
    case _INT___CPLX: case _ZINT__CPLX: case _DOUBLE___CPLX: case _REAL__CPLX:
      return gen(a-*b._CPLXptr,-*(b._CPLXptr+1));
    case _CPLX__CPLX:
      return gen(*a._CPLXptr - *b._CPLXptr, *(a._CPLXptr+1) - *(b._CPLXptr+1));
    case _POLY__POLY:
      return subpoly(a,b);
    case _FRAC__FRAC:
        return (*a._FRACptr)-(*b._FRACptr);
    case _SPOL1__SPOL1:
      return spsub(*a._SPOL1ptr,*b._SPOL1ptr,contextptr);
    case _EXT__EXT:
      return ext_sub(a,b,contextptr);
    case _POLY__INT_: case _POLY__ZINT: case _POLY__DOUBLE_: case _POLY__CPLX: case _POLY__MOD: case _POLY__REAL: case _POLY__USER:
      return subpoly(*a._POLYptr,b);
    case _INT___POLY: case _ZINT__POLY: case _DOUBLE___POLY: case _CPLX__POLY: case _MOD__POLY:
      return subpoly(a,*b._POLYptr);        
    case _MOD__MOD:
      return modsub(a._MODptr,b._MODptr);
    case _REAL__REAL:
      return (*a._REALptr)-(*b._REALptr);
    default:
      if (a.type==_USER)
	return (*a._USERptr)-b;
      if (b.type==_USER)
	return (-b)+a;
      if (a.type==_REAL)
	return (*a._REALptr)-b;
      if (b.type==_REAL)
	return -(*b._REALptr)+a;
      if (a.type==_STRNG)
	return a;
      return sym_sub(a,b,contextptr);
    }
  }

  gen operator_minus (const gen & a,const gen & b,GIAC_CONTEXT){
    register unsigned t=(a.type<< _DECALAGE) | b.type;
    if (!t)
      return((longlong) a.val-b.val);
    return operator_minus(a,b,t,contextptr);
  }

  gen operator - (const gen & a,const gen & b){
    register unsigned t=(a.type<< _DECALAGE) | b.type;
    if (!t)
      return((longlong) a.val-b.val);
    return operator_minus(a,b,t,context0);
  }

  gen sym_sub(const gen & a,const gen & b,GIAC_CONTEXT){
    control_c();
    if (a.is_symb_of_sommet(at_unit) || b.is_symb_of_sommet(at_unit))
      return a+(-b);
    if ( a.is_approx()){
      gen b1;
      if (has_evalf(b,b1,1,contextptr) && b!=b1)
	return a-b1;
    }
    if ( b.is_approx()){
      gen a1;
      if (has_evalf(a,a1,1,contextptr) && a!=a1)
	return a1-b;
    }
    if ( (a.type==_SYMB) && equalposcomp(plot_sommets,a._SYMBptr->sommet) ){
      if ( (b.type==_SYMB) && equalposcomp(plot_sommets,b._SYMBptr->sommet) )
	return a._SYMBptr->feuille._VECTptr->front()-b._SYMBptr->feuille._VECTptr->front();
      else
	return symbolic_plot_makevecteur(a._SYMBptr->sommet,a._SYMBptr->feuille-b,true,contextptr);
    }
    if ( (b.type==_SYMB) && equalposcomp(plot_sommets,b._SYMBptr->sommet) )
      return sym_add(-b,a,contextptr);
    if (a.type==_VECT)
      return sym_add(a,-b,contextptr);
    if (b.type==_VECT)
      return sym_add(-b,a,contextptr);
    if ((a==undef) || (b==undef))
      return undef;
    if (is_inf(a)){
      if (is_inf(b)){
	if ((a==plus_inf) && (b==minus_inf))
	  return a;
	if ((a==minus_inf) && (b==plus_inf))	
	  return a;
	return undef;
      }
      else
	return a;
    }
    if (a.type==_FRAC){
      if ( (b.type!=_SYMB) && (b.type!=_IDNT) )
        return (*a._FRACptr)-b;
      return sym_sub(_FRAC2_SYMB(a),b,contextptr);
    }
    if (b.type==_FRAC){
      if ( (a.type!=_SYMB) && (a.type!=_IDNT) )
        return a-(*b._FRACptr);
      return sym_sub(a,_FRAC2_SYMB(b),contextptr);
    }
    if (a.type==_EXT){
        if (a.is_constant() && (b.type==_POLY))
            return subpoly(a,*b._POLYptr);
        else
            return algebraic_EXTension(*a._EXTptr-b,*(a._EXTptr+1));
    }
    if (b.type==_EXT){
        if (b.is_constant() && (a.type==_POLY))
            return subpoly(*a._POLYptr,b);
        else
            return algebraic_EXTension(a-*b._EXTptr,*(b._EXTptr+1));
    }
    if (a==b)
      return chkmod(zero,a);
    if (is_inf(b))
      return -b;
    if (is_exactly_zero(b))
      return a;
    if (is_exactly_zero(a))
      return -b;
    /*
    if (a.type==_SYMB && a._SYMBptr->sommet==at_equal){
      vecteur & va=*a._SYMBptr->feuille._VECTptr;
      if (b.type==_SYMB && b._SYMBptr->sommet==at_equal){
	vecteur & vb=*b._SYMBptr->feuille._VECTptr;
	return new symbolic(at_equal,makevecteur(va.front()-vb.front(),va.back()-vb.back()));
      }
      else
	return new symbolic(at_equal,makenewvecteur(va.front()-b,va.back()-b));
    }
    if (b.type==_SYMB && b._SYMBptr->sommet==at_equal){
      vecteur & vb=*b._SYMBptr->feuille._VECTptr;
      return new symbolic(at_equal,makenewvecteur(a-vb.front(),a-vb.back()));
    }
    */
    if (is_inequality(a) || is_inequality(b))
      return a+(-b);
    if ((a.type==_SYMB) && (b.type==_SYMB)){
      if (a._SYMBptr->sommet==at_plus) {
	if (b._SYMBptr->sommet==at_plus)
	  return new symbolic(at_plus,mergevecteur(*(a._SYMBptr->feuille._VECTptr),negvecteur(*(b._SYMBptr->feuille._VECTptr))));
	else
	  return new symbolic(*a._SYMBptr,-b);
      }
      else { 
	if (b._SYMBptr->sommet==at_plus)
	  return new symbolic(*(-b)._SYMBptr,a);
	else
	  return new symbolic(at_plus,gen(makenewvecteur(a,-b)));
      }
    }
    if (b.type==_SYMB){
      if (b._SYMBptr->sommet==at_plus)
	return new symbolic(*(-b)._SYMBptr,a);
      else
	return new symbolic(at_plus,gen(makenewvecteur(a,-b)));
    }
    if (a.type==_SYMB){
      if (a._SYMBptr->sommet==at_plus)
	return new symbolic(*a._SYMBptr,-b);
      else
	return new symbolic(at_plus,gen(makenewvecteur(a,-b)));
    }
    if ((a.type==_IDNT) || (b.type==_IDNT))
      return new symbolic(at_plus,gen(makenewvecteur(a,-b)));
    if (a.type==_MOD)
      return a-makemod(b,*(a._MODptr+1));
    if (b.type==_MOD)
      return makemod(a,*(b._MODptr+1))-b;
    return new symbolic(at_plus,makenewvecteur(a,-b));
    // settypeerr("sym_sub");
  }

  vecteur negfirst(const vecteur & v){
    vecteur w(v);
    if (!w.empty())
      w.front()=-w.front();
    return w;
  }

  real_object real_object::operator -() const {
    real_object res(*this);
#ifdef HAVE_LIBMPFR
    mpfr_neg(res.inf,res.inf,GMP_RNDN);
#else
    mpf_neg(res.inf,res.inf);
#endif
    return res;
  }
    
  real_object real_object::inv() const {
    real_object res(*this);
#ifdef HAVE_LIBMPFR
    mpfr_ui_div(res.inf,1,res.inf,GMP_RNDN);
#else
    mpf_ui_div(res.inf,1,res.inf);
#endif
    return res;
  }
    
  real_object real_object::sqrt() const {
    real_object res(*this);
#ifdef HAVE_LIBMPFR
    mpfr_sqrt(res.inf,res.inf,GMP_RNDN);
#else
    mpf_sqrt(res.inf,res.inf);
#endif
    return res;
  }
    
  real_object real_object::abs() const {
#ifdef HAVE_LIBMPFR
    if (mpfr_sgn(inf)>=0)
#else
    if (mpf_sgn(inf)>=0)
#endif
      return *this;
    return -(*this);
  }

  void compile_with_mpfr(){
    setsizeerr("Compile with MPFR or USE_GMP_REPLACEMENTS if you want transcendental long float support");
  }

  real_object real_object::exp() const {
#ifdef USE_GMP_REPLACEMENTS
    *res.inf = ::exp(*res.inf);
#else
    real_object res(*this);
#ifdef HAVE_LIBMPFR
    mpfr_exp(res.inf,res.inf,GMP_RNDN);
#else
    compile_with_mpfr();
#endif
#endif
    return res;
  }

  real_object real_object::log() const {
#ifdef USE_GMP_REPLACEMENTS
    *res.inf = ::log(*res.inf);
#else
    real_object res(*this);
#ifdef HAVE_LIBMPFR
    mpfr_log(res.inf,res.inf,GMP_RNDN);
#else
    compile_with_mpfr();
#endif
#endif
    return res;
  }

  real_object real_object::sin() const {
#ifdef USE_GMP_REPLACEMENTS
    *res.inf = ::sin(*res.inf);
#else
    real_object res(*this);
#ifdef HAVE_LIBMPFR
    mpfr_sin(res.inf,res.inf,GMP_RNDN);
#else
    compile_with_mpfr();
#endif
#endif
    return res;
  }

  real_object real_object::cos() const {
#ifdef USE_GMP_REPLACEMENTS
    *res.inf = ::cos(*res.inf);
#else
    real_object res(*this);
#ifdef HAVE_LIBMPFR
    mpfr_cos(res.inf,res.inf,GMP_RNDN);
#else
    compile_with_mpfr();
#endif
#endif
    return res;
  }

  real_object real_object::tan() const {
#ifdef USE_GMP_REPLACEMENTS
    *res.inf = ::tan(*res.inf);
#else
    real_object res(*this);
#ifdef HAVE_LIBMPFR
    mpfr_tan(res.inf,res.inf,GMP_RNDN);
#else
    compile_with_mpfr();
#endif
#endif
    return res;
  }

  real_object real_object::sinh() const {
#ifdef USE_GMP_REPLACEMENTS
    *res.inf = ::sinh(*res.inf);
#else
    real_object res(*this);
#ifdef HAVE_LIBMPFR
    mpfr_sinh(res.inf,res.inf,GMP_RNDN);
#else
    compile_with_mpfr();
#endif
#endif
    return res;
  }

  real_object real_object::cosh() const {
#ifdef USE_GMP_REPLACEMENTS
    *res.inf = ::cosh(*res.inf);
#else
    real_object res(*this);
#ifdef HAVE_LIBMPFR
    mpfr_cosh(res.inf,res.inf,GMP_RNDN);
#else
    compile_with_mpfr();
#endif
#endif
    return res;
  }

  real_object real_object::tanh() const {
#ifdef USE_GMP_REPLACEMENTS
    *res.inf = ::tanh(*res.inf);
#else
    real_object res(*this);
#ifdef HAVE_LIBMPFR
    mpfr_tanh(res.inf,res.inf,GMP_RNDN);
#else
    compile_with_mpfr();
#endif
#endif
    return res;
  }

  real_object real_object::asin() const {
#ifdef USE_GMP_REPLACEMENTS
    *res.inf = ::asin(*res.inf);
#else
    real_object res(*this);
#ifdef HAVE_LIBMPFR
    mpfr_asin(res.inf,res.inf,GMP_RNDN);
#else
    compile_with_mpfr();
#endif
#endif
    return res;
  }

  real_object real_object::acos() const {
#ifdef USE_GMP_REPLACEMENTS
    *res.inf = ::acos(*res.inf);
#else
    real_object res(*this);
#ifdef HAVE_LIBMPFR
    mpfr_acos(res.inf,res.inf,GMP_RNDN);
#else
    compile_with_mpfr();
#endif
#endif
    return res;
  }

  real_object real_object::atan() const {
#ifdef USE_GMP_REPLACEMENTS
    *res.inf = ::atan(*res.inf);
#else
    real_object res(*this);
#ifdef HAVE_LIBMPFR
    mpfr_atan(res.inf,res.inf,GMP_RNDN);
#else
    compile_with_mpfr();
#endif
#endif
    return res;
  }

  real_object real_object::asinh() const {
#ifdef USE_GMP_REPLACEMENTS
    *res.inf = ::asinh(*res.inf);
#else
    real_object res(*this);
#ifdef HAVE_LIBMPFR
    mpfr_asinh(res.inf,res.inf,GMP_RNDN);
#else
    compile_with_mpfr();
#endif
#endif
    return res;
  }

  real_object real_object::acosh() const {
#ifdef USE_GMP_REPLACEMENTS
    *res.inf = ::acosh(*res.inf);
#else
    real_object res(*this);
#ifdef HAVE_LIBMPFR
    mpfr_acosh(res.inf,res.inf,GMP_RNDN);
#else
    compile_with_mpfr();
#endif
#endif
    return res;
  }

  real_object real_object::atanh() const {
#ifdef USE_GMP_REPLACEMENTS
    *res.inf = ::atanh(*res.inf);
#else
    real_object res(*this);
#ifdef HAVE_LIBMPFR
    mpfr_atanh(res.inf,res.inf,GMP_RNDN);
#else
    compile_with_mpfr();
#endif
#endif
    return res;
  }

  real_object real_interval::operator -() const {
    real_interval res(*this);
#ifdef HAVE_LIBMPFR
    mpfr_neg(res.inf,res.inf,GMP_RNDU);
#ifdef HAVE_LIBMPFI
    mpfi_neg(res.infsup,res.infsup);
#else
    mpfr_neg(res.sup,res.sup,GMP_RNDD);
    mpfr_swap(res.inf,res.sup);
#endif
#else // MPFR
    mpf_neg(res.inf,res.inf);
#ifdef HAVE_LIBMPFI
    mpfi_neg(res.infsup,res.infsup);
#else
    mpf_neg(res.sup,res.sup);
#ifdef mpf_swap
    mpf_swap(res.inf,res.sup);
#endif
#endif
#endif // MPFR
    return res;
  }

  real_object real_interval::inv() const {
    real_interval res(*this);
#ifdef HAVE_LIBMPFI
    mpfi_ui_div(res.infsup,1,res.infsup);
    mpfr_ui_div(res.inf,1,res.inf,GMP_RNDD);
#else
    // FIXME check sign
    setsizeerr();
    /* mpf_neg(res.inf,res.inf);
       mpf_neg(res.sup,res.sup);
       mpf_swap(res.inf,res.sup); */
#endif
    return res;
  }

  gen operator -(const gen & a){
    mpz_t *e ;
    switch (a.type ) {
    case _INT_: 
      return(-a.val);
    case _ZINT: 
      e=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init(*e);
      mpz_neg(*e,*a._ZINTptr);
      return(e);
    case _DOUBLE_:
      return -(a._DOUBLE_val);
    case _CPLX:
      return gen(-*a._CPLXptr,-*(a._CPLXptr+1));
    case _IDNT:
      if ((a==undef) || (a==unsigned_inf))
	return a;
      return new symbolic(at_neg,a);
    case _SYMB:
      if (a==plus_inf)
	return minus_inf;
      if (a==minus_inf)
	return plus_inf;
      if (a._SYMBptr->sommet==at_neg)
	return a._SYMBptr->feuille;
      if (a._SYMBptr->sommet==at_unit)
	return new symbolic(at_unit,makenewvecteur(-a._SYMBptr->feuille._VECTptr->front(),a._SYMBptr->feuille._VECTptr->back()));
      if (a._SYMBptr->sommet==at_plus)
	return new symbolic(at_plus,negvecteur(*a._SYMBptr->feuille._VECTptr));
      if (a._SYMBptr->sommet==at_interval && a._SYMBptr->feuille.type==_VECT && a._SYMBptr->feuille._VECTptr->size()==2){
	return new symbolic(at_interval,gen(makenewvecteur(-a._SYMBptr->feuille._VECTptr->back(),-a._SYMBptr->feuille._VECTptr->front()),_SEQ__VECT));
      }
      if (equalposcomp(plot_sommets,a._SYMBptr->sommet)){
	return symbolic_plot_makevecteur(a._SYMBptr->sommet,-a._SYMBptr->feuille,false,context0);
      }
      if (a._SYMBptr->sommet==at_equal || a._SYMBptr->sommet==at_different || a._SYMBptr->sommet==at_same)
	 return new symbolic(a._SYMBptr->sommet,makenewvecteur(-a._SYMBptr->feuille._VECTptr->front(),-a._SYMBptr->feuille._VECTptr->back())); 
      if (is_inequality(a))
	return new symbolic(a._SYMBptr->sommet,makenewvecteur(-a._SYMBptr->feuille._VECTptr->back(),-a._SYMBptr->feuille._VECTptr->front()));
      return new symbolic(at_neg,a);
    case _VECT:
      if (a.subtype==_VECTOR__VECT && a._VECTptr->size()==2)
	return gen(makenewvecteur(a._VECTptr->back(),a._VECTptr->front()),_VECTOR__VECT);
      if (a.subtype==_PNT__VECT)
	return gen(negfirst(*a._VECTptr),a.subtype);
      return gen(negvecteur(*a._VECTptr),a.subtype);
    case _POLY:
      return -(*a._POLYptr);
    case _EXT:
      return algebraic_EXTension(-(*a._EXTptr),*(a._EXTptr+1));
    case _USER:
      return -(*a._USERptr);
    case _MOD:
      return makemod(-*a._MODptr,*(a._MODptr+1));
    case _FRAC:
      return fraction(-(a._FRACptr->num),a._FRACptr->den);
    case _SPOL1:
      return spneg(*a._SPOL1ptr,context0);
    case _STRNG:
      return string2gen("-"+(*a._STRNGptr),false);
    case _REAL:
      return -*a._REALptr;
    default: 
      return new symbolic(at_neg,a);
    }
  }

  gen mulpoly(const gen & th,const gen & other){
    if ((th.type!=_POLY) || (other.type!=_POLY))
      settypeerr("mulpoly");
    vector< monomial<gen> >::const_iterator ita = th._POLYptr->coord.begin();
    vector< monomial<gen> >::const_iterator ita_end = th._POLYptr->coord.end();
    vector< monomial<gen> >::const_iterator itb = other._POLYptr->coord.begin();
    vector< monomial<gen> >::const_iterator itb_end = other._POLYptr->coord.end();
    // first some trivial cases
    if (ita==ita_end)
      return(th);
    if (itb==itb_end)
      return(other);
    if (is_one(*th._POLYptr))
      return other;
    if (is_one(*other._POLYptr))
      return th;
    // Now look if length a=1 or length b=1, happens frequently
    // think of x^3*y^2*z translated to internal form
    int c1=th._POLYptr->coord.size();
    if (c1==1)
      return other._POLYptr->shift(th._POLYptr->coord.front().index,th._POLYptr->coord.front().value);
    int c2=other._POLYptr->coord.size();
    if (c2==1)
      return th._POLYptr->shift(other._POLYptr->coord.front().index,other._POLYptr->coord.front().value);
    polynome * resptr = new polynome(th._POLYptr->dim);
    mulpoly(*th._POLYptr,*other._POLYptr,*resptr,0);
    return resptr;
  }

  vecteur multfirst(const gen & a,const vecteur & v){
    vecteur w(v);
    if (!w.empty())
      w.front()=v.front()*a;
    return w;
  }

  gen real_object::operator * (const gen & g) const{
    switch (g.type){
    case _REAL:
      return *this * *g._REALptr;
    case _INT_: case _DOUBLE_: case _ZINT: case _FRAC:
#ifdef HAVE_LIBMPFR
      return *this * real_object(g,mpfr_get_prec(inf));      
#else
      return *this * real_object(g);
#endif
    default:
      return sym_mult(*this,g,context0);
    }
  }
  
  gen real_object::operator / (const gen & g) const{
    return *this * g.inverse(context0);
  }
  
  real_object real_object::operator / (const real_object & g) const{
    return *this * g.inv();
  }
  
  real_interval mul(const real_interval & i,const real_interval & g){
    real_interval res(i);
#ifdef HAVE_LIBMPFR
    mpfr_mul(res.inf,i.inf,g.inf,GMP_RNDN);
#else
    mpf_mul(res.inf,i.inf,g.inf);
#endif
#ifdef HAVE_LIBMPFI
    mpfi_mul(res.infsup,i.infsup,g.infsup);
#else
    // FIXME: should check signs for interval arithmetic!!
    setsizeerr();
    // mpf_mul(res.sup,i.sup,g.sup); 
#endif
    return res;
  }

  real_interval real_interval::operator * (const real_interval & g) const{
    return mul(*this,g);
  }

  real_interval mul(const real_interval & i,const real_object & g){
    const real_interval * ptr=dynamic_cast<const real_interval *>(&g);
    if (ptr)
      return mul(i,*ptr);
    real_interval res(i);
#ifdef HAVE_LIBMPFR
    mpfr_mul(res.inf,i.inf,g.inf,GMP_RNDN);
#else
    mpf_mul(res.inf,i.inf,g.inf);
#endif
#ifdef HAVE_LIBMPFI
    mpfi_mul_fr(res.infsup,i.infsup,g.inf);
#else
    // FIXME: should check signs for interval arithmetic!!
    setsizeerr();
    // mpf_mul(res.sup,i.sup,g.inf);    
#endif
    return res;
  }

  real_object real_interval::operator * (const real_object & g) const{
    return mul(*this,g);
  }

  real_object real_object::operator * (const real_object & g) const{
    const real_interval * ptr=dynamic_cast<const real_interval *>(&g);
    if (ptr)
      return mul(*ptr,*this);
#ifdef HAVE_LIBMPFR
    mpfr_t sum;
    mpfr_init2(sum,min(mpfr_get_prec(this->inf),mpfr_get_prec(g.inf)));
    mpfr_mul(sum,this->inf,g.inf,GMP_RNDN);
    real_object res(sum);
    mpfr_clear(sum);
#else
    mpf_t sum;
    mpf_init(sum);
    mpf_mul(sum,this->inf,g.inf);
    real_object res(sum);
    mpf_clear(sum);
#endif
    return res;
  }

  gen multgen_poly(const gen & a,const vecteur & b,int subtype){
    gen res(vecteur(0),subtype);
    multvecteur(a,b,*res._VECTptr);
    return res;
  }

  gen multgen_poly(const vecteur & a,const vecteur & b){
    gen res(vecteur(0), _POLY1__VECT);
    operator_times(a,b,0,*res._VECTptr);
    return res;
  }

  // a*b -> tmp, modifies tmp in place
  void type_operator_times(const gen & a,const gen &b,gen & tmp){
    register unsigned t=(a.type<< _DECALAGE) | b.type;
    if (tmp.type==_DOUBLE_ && t==_DOUBLE___DOUBLE_){
      tmp._DOUBLE_val=a._DOUBLE_val*b._DOUBLE_val;
      return ;
    }
    if (!t && tmp.type==_INT_ ){
      register longlong ab=longlong(a.val)*b.val;
      tmp.val=ab;
      if (tmp.val!=ab)
	tmp=ab;
      return;
    }
    if (tmp.type==_ZINT && *tmp.ptr_val.ref_count==1){
      mpz_t * ptr=tmp._ZINTptr;
      switch (t){
      case _INT___INT_:
	tmp=longlong(a.val)*b.val;
	return;
      case _ZINT__ZINT:
	mpz_mul(*ptr,*a._ZINTptr,*b._ZINTptr);
	return ;
      case _ZINT__INT_:
	if (b.val<0){
	  mpz_mul_ui(*ptr,*a._ZINTptr,-b.val);
	  mpz_neg(*ptr,*ptr);
	}
	else
	  mpz_mul_ui(*ptr,*a._ZINTptr,b.val);
	return;
      case _INT___ZINT:
	if (a.val<0){
	  mpz_mul_ui(*ptr,*b._ZINTptr,-a.val);
	  mpz_neg(*ptr,*ptr);
	}
	else
	  mpz_mul_ui(*ptr,*b._ZINTptr,a.val);
	return;
      }
    }
    tmp=a*b;
  }

  gen operator_times (const gen & a,const gen & b,unsigned t,GIAC_CONTEXT){
    // cout << a << "*" << b << endl;
    // if (!( (++control_c_counter) & control_c_counter_mask))
      control_c();
    register  mpz_t * e;
    switch (t) {
    case _ZINT__ZINT:
      e=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init(*e);
      mpz_mul(*e,*a._ZINTptr,*b._ZINTptr);
      return e;
    case _DOUBLE___DOUBLE_:
      return a._DOUBLE_val*b._DOUBLE_val;
    case _INT___ZINT: 
      e=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init(*e);
      if (a.val<0){
	mpz_mul_ui(*e,*b._ZINTptr,-a.val);
	mpz_neg(*e,*e);
      }
      else
	mpz_mul_ui(*e,*b._ZINTptr,a.val);
      return gen(e);
    case _ZINT__INT_:
      e=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init(*e);
      if (b.val<0){
	mpz_mul_ui(*e,*a._ZINTptr,-b.val);
	mpz_neg(*e,*e);
      }
      else
	mpz_mul_ui(*e,*a._ZINTptr,b.val);
      return gen(e);
    case _INT___DOUBLE_:
      return a.val*b._DOUBLE_val;
    case _DOUBLE___INT_:
      return a._DOUBLE_val*b.val;
    case _DOUBLE___ZINT:
      return a._DOUBLE_val*mpz_get_d(*b._ZINTptr);
    case _DOUBLE___REAL:
      return a._DOUBLE_val*real2double(*b._REALptr);
    case _REAL__DOUBLE_:
      return b._DOUBLE_val*real2double(*a._REALptr);
    case _ZINT__DOUBLE_:
      return mpz_get_d(*a._ZINTptr)*b._DOUBLE_val;
    case _CPLX__INT_: case _CPLX__ZINT: case _CPLX__DOUBLE_: case _CPLX__REAL:
      return gen(*a._CPLXptr*b,*(a._CPLXptr+1)*b);
    case _INT___CPLX: case _ZINT__CPLX: case _DOUBLE___CPLX: case _REAL__CPLX:
      return gen(a*(*b._CPLXptr),a*(*(b._CPLXptr+1)));
    case _CPLX__CPLX:
      return gen(*a._CPLXptr * (*b._CPLXptr) - *(a._CPLXptr+1)* (*(b._CPLXptr+1)), (*b._CPLXptr) * (*(a._CPLXptr+1)) + *(b._CPLXptr+1) * (*a._CPLXptr));
    case _VECT__INT_: case _VECT__ZINT: case _VECT__DOUBLE_: case _VECT__CPLX: case _VECT__SYMB: case _VECT__IDNT: case _VECT__POLY: case _VECT__EXT: case _VECT__MOD: case _VECT__FRAC: case _VECT__REAL:
      if (a.subtype==_VECTOR__VECT && a._VECTptr->size()==2)
	return vector2vecteur(*a._VECTptr)*b;
      if (a.subtype==_PNT__VECT)
	return gen(multfirst(b,*a._VECTptr),_PNT__VECT);
      if (a.subtype==_POLY1__VECT && is_zero(b))
	return b;
      return multgen_poly(b,*a._VECTptr,a.subtype); // gen(multvecteur(b,*a._VECTptr),a.subtype);
    case _INT___VECT: case _ZINT__VECT: case _DOUBLE___VECT: case _CPLX__VECT: case _SYMB__VECT: case _IDNT__VECT: case _POLY__VECT: case _EXT__VECT: case _MOD__VECT: case _FRAC__VECT: case _REAL__VECT:
      if (b.subtype==_VECTOR__VECT && b._VECTptr->size()==2)
	return a*vector2vecteur(*b._VECTptr);
      if (b.subtype==_PNT__VECT)
	return gen(multfirst(a,*b._VECTptr),_PNT__VECT);
      if (b.subtype==_POLY1__VECT && is_zero(a))
	return a;
      return multgen_poly(a,*b._VECTptr,b.subtype); // gen(multvecteur(a,*b._VECTptr),b.subtype);
    case _VECT__VECT:
      if ( (a.subtype==_POLY1__VECT) || (b.subtype==_POLY1__VECT) )
	return multgen_poly(*a._VECTptr,*b._VECTptr);
      return ckmultmatvecteur(*a._VECTptr,*b._VECTptr);
    case _POLY__POLY:
      return mulpoly(a,b);
    case _FRAC__FRAC:
      if (a._FRACptr->num.type==_EXT && b._FRACptr->num.type==_EXT)
	return ((*a._FRACptr)*(*b._FRACptr)).normal();	
      return (*a._FRACptr)*(*b._FRACptr);
    case _SPOL1__SPOL1:
      return spmul(*a._SPOL1ptr,*b._SPOL1ptr,contextptr);
    case _EXT__EXT:
      return ext_mul(a,b,contextptr);
    case _POLY__INT_: case _POLY__ZINT: case _POLY__DOUBLE_: case _POLY__CPLX: case _POLY__USER: case _POLY__REAL:
      if (is_one(b))
	return a;
      return (*a._POLYptr) * b;
    case _POLY__MOD:
      return (*a._POLYptr) * b;
    case _INT___POLY: case _ZINT__POLY: case _DOUBLE___POLY: case _CPLX__POLY: case _USER__POLY: case _REAL__POLY:
      if (is_one(a))
	return b;
      return a * (*b._POLYptr);        
    case _MOD__POLY:
      return a * (*b._POLYptr);        
    case _MOD__MOD:
      return modmul(a._MODptr,b._MODptr);
    case _MOD__INT_: case _MOD__ZINT:
      return makemod(*a._MODptr*b,*(a._MODptr+1));
    case _INT___MOD: case _ZINT__MOD:
      return makemod(*b._MODptr*a,*(b._MODptr+1));
    case _REAL__REAL:
      return (*a._REALptr)*(*b._REALptr);
    default:
      if (a.type==_USER)
	return (*a._USERptr)*b;
      if (b.type==_USER)
	return (*b._USERptr)*a;      
      if (a.type==_REAL)
	return (*a._REALptr)*b;
      if (b.type==_REAL)
	return (*b._REALptr)*a;
      return sym_mult(a,b,contextptr);
    }
  }

  gen operator_times (const gen & a,const gen & b,GIAC_CONTEXT){
    register unsigned t=(a.type<< _DECALAGE) | b.type;
    if (!t)
      return gen((longlong) a.val*b.val);
    return operator_times(a,b,t,contextptr);
  }

  gen operator * (const gen & a,const gen & b){
    register unsigned t=(a.type<< _DECALAGE) | b.type;
    if (!t)
      return gen((longlong) a.val*b.val);
    return operator_times(a,b,t,context0);
  }

  bool has_i(const gen & g){
    if (g.type==_CPLX)
      return true;
    if (g.type==_VECT){
      const_iterateur it=g._VECTptr->begin(),itend=g._VECTptr->end();
      for (;it!=itend;++it){
	if (has_i(*it))
	  return true;
      }
      return false;
    }
    if (g.type!=_SYMB)
      return false;
    return has_i(g._SYMBptr->feuille);
  }

  gen giac_pow(const gen & base,const gen & exponent,GIAC_CONTEXT){
    return pow(base,exponent,contextptr);
  }

  // (-1)^n
  gen minus1pow(const gen & exponent,GIAC_CONTEXT){
    if (exponent.type==_INT_)
      return (exponent.val%2)?-1:1;
    if (exponent.type==_ZINT){
      gen q,g=irem(exponent,2,q);
      if (is_zero(g))
	return 1;
      return -1;
    }
    if (is_inf(exponent) || is_undef(exponent))
      return undef;
    if (exponent.is_symb_of_sommet(at_neg))
      return minus1pow(exponent._SYMBptr->feuille,contextptr);
    if (exponent.is_symb_of_sommet(at_plus)){
      gen res(1);
      gen & f=exponent._SYMBptr->feuille;
      if (f.type!=_VECT)
	return minus1pow(f,contextptr);
      vecteur & v = *f._VECTptr;
      int s=v.size();
      for (int i=0;i<s;++i)
	res = res * minus1pow(v[i],contextptr);
      return res;
    }
    if (exponent.is_symb_of_sommet(at_prod)){
      gen & f =exponent._SYMBptr->feuille;
      if (f.type==_VECT){
	vecteur & v = *f._VECTptr;
	int i,s=v.size();
	bool even=false;
	for (i=0;i<s;++i){
	  if (is_integer(v[i]) && is_zero(smod(v[i],2)))
	    even=true;
	  if (!is_assumed_integer(v[i],contextptr))
	    break;
	}
	if (even && i==s)
	  return 1;
      }
    }
    return new symbolic(at_pow,gen(makenewvecteur(-1,exponent),_SEQ__VECT));
  }

  gen pow(const gen & base,const gen & exponent,GIAC_CONTEXT){
    // if (!( (++control_c_counter) & control_c_counter_mask))
      control_c();
    if (is_undef(base))
      return base;
    if (base.type==_VECT && base.subtype!=_POLY1__VECT && !ckmatrix(base)){
      if (exponent.type==_VECT)
	return apply(base,exponent,contextptr,giac::giac_pow);
      return apply1st(base,exponent,contextptr,&giac::giac_pow); 
    }
    if (exponent.type==_VECT)
      return apply2nd(base,exponent,contextptr,&giac::giac_pow);
    if (exponent.type==_FRAC)
      return pow(base,new symbolic(at_prod,makenewvecteur(exponent._FRACptr->num,symb_inv(exponent._FRACptr->den))),contextptr);
    if (is_inf(base)){ 
      if (is_zero(exponent))
	return undef;
      gen d;
      bool b=has_evalf(exponent,d,1,contextptr);
      if (b && is_strictly_positive(-exponent,contextptr) )
	return 0;
      if (b && base==plus_inf &&is_strictly_positive(exponent,contextptr))
	return plus_inf;
      if ( (exponent.type==_INT_) ){
	if (exponent.val % 2)
	  return base;
	else
	  return plus_inf; // for unsigned_inf in _DOUBLE_ mode only!!
      }
      if (b && is_strictly_positive(exponent,contextptr))
	return unsigned_inf;
      return undef;
    }
    if (base.type==_SYMB){ 
      unary_function_ptr & u =base._SYMBptr->sommet;
      if (u==at_unit){
	vecteur & v=*base._SYMBptr->feuille._VECTptr;
	return new symbolic(at_unit,makenewvecteur(pow(v[0],exponent,contextptr),pow(v[1],exponent,contextptr)));
      }
      if (u==at_abs && exponent.type==_INT_ && !complex_mode(contextptr) && !has_i(base)){ 
	int n=exponent.val,m;
	if (n<0 && n%2)
	  m=n-1;
	else
	  m=(n/2)*2;
	gen basep=pow(base._SYMBptr->feuille,m);
	if (n%2)
	  return base*basep;
	else
	  return basep;
      }
      if (u==at_exp){
	// (e^a)^b=a^(a*b)
	// but we keep (e^a)^b if b is integer and e^(a*b) is not simplified
	// for rational dependance
	gen res=exp(base._SYMBptr->feuille*exponent,contextptr);
	if (exponent.type!=_INT_ || !res.is_symb_of_sommet(at_exp))
	  return res;
      }
      if (u==at_pow && !has_i(base)){
	vecteur & v=*base._SYMBptr->feuille._VECTptr;
	gen & v1=v[1];
	gen new_exp=v1*exponent;
	if (new_exp.type>_IDNT)
	  new_exp=normal(new_exp,contextptr);
	if ( v1.type==_INT_ && v1.val%2==0 && new_exp.type==_INT_ && new_exp.val%2 && !complex_mode(contextptr) ) 
	  return pow(abs(v[0],contextptr),new_exp,contextptr); 
	else 
	  return pow(v[0],new_exp,contextptr);
      }
      if (u==at_equal){
	vecteur & vb=*base._SYMBptr->feuille._VECTptr;
	return new symbolic(base._SYMBptr->sommet,makenewvecteur(pow(vb.front(),exponent,contextptr),pow(vb.back(),exponent,contextptr)));
      }
      if (exponent.type==_INT_){
	if (exponent.val==0)
	  return 1;
	if (exponent.val==1)
	  return base;
	return new symbolic(at_pow,makenewvecteur(base,exponent));
      }
    }
    switch ( (base.type<< _DECALAGE) | exponent.type ) {
    case _INT___INT_: case _ZINT__INT_: case _REAL__INT_: case _CPLX__INT_: case _IDNT__INT_: 
      return pow(base,exponent.val);
    case _DOUBLE___DOUBLE_:
      if (exponent._DOUBLE_val==int(exponent._DOUBLE_val+.5))
	return pow(base,int(exponent._DOUBLE_val+.5));
      if (base._DOUBLE_val>=0)
#ifdef _SOFTMATH_H
	return std::giac_gnuwince_pow(base._DOUBLE_val,exponent._DOUBLE_val);
#else
	return std::pow(base._DOUBLE_val,exponent._DOUBLE_val);
#endif
      else
	return exp(exponent*log(base,contextptr),contextptr);
    case _DOUBLE___INT_:
      if (base._DOUBLE_val>=0)
#ifdef _SOFTMATH_H
	return std::giac_gnuwince_pow(base._DOUBLE_val,exponent.val);
#else
      return std::pow(base._DOUBLE_val,exponent.val);
#endif
      else 
	return (exponent.val%2?-1:1)*exp(exponent*log(-base,contextptr),contextptr);
    case _INT___DOUBLE_:
#ifdef _SOFTMATH_H
      return std::giac_gnuwince_pow(base.val,exponent._DOUBLE_val);
#else
      return std::pow(base.val,exponent._DOUBLE_val);
#endif
    case _ZINT__DOUBLE_:
#ifdef _SOFTMATH_H
      return std::giac_gnuwince_pow(mpz_get_d(*base._ZINTptr),exponent._DOUBLE_val);
#else
      return std::pow(mpz_get_d(*base._ZINTptr),exponent._DOUBLE_val);
#endif
    case _DOUBLE___ZINT:
      if (base._DOUBLE_val>=0)
#ifdef _SOFTMATH_H
	return std::giac_gnuwince_pow(base._DOUBLE_val,mpz_get_d(*exponent._ZINTptr));
#else
      return std::pow(base._DOUBLE_val,mpz_get_d(*exponent._ZINTptr));
#endif
      else
	return exp(exponent*log(base,contextptr),contextptr);
    case _POLY__INT_:
      if (exponent.val<0)
	return fraction(1,pow(*base._POLYptr,-exponent.val));
      else
	return pow(*base._POLYptr,exponent.val);
    case _FRAC__INT_:
      return pow(*base._FRACptr,exponent.val);
    case _EXT__INT_: case _MOD__INT_: case _VECT__INT_: case _USER__INT_:
      return pow(base,exponent.val);
    case _MOD__ZINT:
      return makemod(powmod(*base._MODptr,exponent,*(base._MODptr+1)),*(base._MODptr+1));
    default:
      if ((base==undef) || (exponent==undef))
	return undef;
      if (is_one(base) && !is_inf(exponent))
	return base;
      if (base.type==_REAL || base.type==_DOUBLE_ || exponent.type==_REAL || exponent.type==_DOUBLE_ )
	return exp(exponent*log(base,contextptr),contextptr);
      /* 
	 if (base.is_symb_of_sommet(at_neg))
	 return minus1pow(exponent,contextptr)*pow(base._SYMBptr->feuille,exponent);
      */
      if ((base.type==_INT_) && (base.val<0)){
	if (exponent==plus_one_half)
	  return cst_i*sqrt(-base.val,contextptr);
	// if (exponent==-one_half)
	//  return rdiv(cst_i,sqrt(-base.val));
      }
      if (is_exactly_zero(base)){
	gen d;
	bool b=has_evalf(exponent,d,1,contextptr);
	if (b && is_positive(exponent,contextptr)) 
	  return base;
	if (b && is_positive(-exponent,contextptr))
	  return unsigned_inf;
	return undef;
      }
      if (is_integer(base) && is_positive(-base,contextptr))
	return minus1pow(exponent,contextptr)*pow(-base,exponent,contextptr);
      if (is_inf(exponent))
	return exp(exponent*ln(base,contextptr),contextptr);
      // extract integral powers in a product exponent
      if ((exponent.type==_SYMB) && (exponent._SYMBptr->sommet==at_prod)){
	gen subexponent_num(1),subexponent_deno(1);
	gen superexponent(1);
	const_iterateur it=exponent._SYMBptr->feuille._VECTptr->begin(),itend=exponent._SYMBptr->feuille._VECTptr->end();
	for (;it!=itend;++it){
	  if (it->type==_INT_){
	    superexponent = superexponent * (*it);
	    continue;
	  }
	  if ( (it->type==_SYMB) && (it->_SYMBptr->sommet==at_inv))
	    subexponent_deno = subexponent_deno * (it->_SYMBptr->feuille);
	  else
	    subexponent_num = subexponent_num * (*it);
	}
	if (superexponent.type!=_INT_)
	  return new symbolic(at_pow,makenewvecteur(base,exponent));
	if (subexponent_deno.type!=_INT_)
	  return new symbolic(at_pow,makenewvecteur(new symbolic(at_pow,makenewvecteur(base,_FRAC2_SYMB(subexponent_num,subexponent_deno))),superexponent));
	int q=superexponent.val / subexponent_deno.val;
	int r=superexponent.val % subexponent_deno.val;
	gen res(1);
	if (r){
	  if (fastsign(base,contextptr)==-1){ // is_strictly_positive(-base,contextptr)){
	    gen base1=-base;
	    res=exp((cst_i*cst_pi*subexponent_num)/subexponent_deno,contextptr);
	    res=new symbolic(at_pow,makenewvecteur(pow(base1,subexponent_num,contextptr),inv(subexponent_deno,contextptr)))*res;
	  }
	  else
	    res=new symbolic(at_pow,makenewvecteur(pow(base,subexponent_num,contextptr),inv(subexponent_deno,contextptr)));
	  if (r!=1)
	    res=new symbolic(at_pow,makenewvecteur(res,r));
	}
	if (!q)
	  return res;
	if (q==1)
	  return res*new symbolic(at_pow,makenewvecteur(base,subexponent_num));
	if (q==-1)
	  return res*inv(pow(base,subexponent_num,contextptr),contextptr);
	return res*new symbolic(at_pow,makenewvecteur(pow(base,subexponent_num,contextptr),q));
      }
      return new symbolic(at_pow,makenewvecteur(base,exponent));
    }  
  }

  // 0 if unknown, 1 if >0, -1 if <0
  // no test for symbolics if context_ptr=0
  int fastsign(const gen & a,GIAC_CONTEXT){
    if (is_zero(a) || is_undef(a))
      return 0;
    if (is_inf(a)){
      if (a==plus_inf)
	return 1;
      if (a==minus_inf)
	return -1;
      return 0;
    }
    switch (a.type) {
    case _INT_: 
      if (a.val>0)
	return 1;
      else
	return -1;
    case _ZINT: 
      return mpz_cmp_si(*a._ZINTptr,0);
    case _FRAC:
      return fastsign(a._FRACptr->num,contextptr)*fastsign(a._FRACptr->den,contextptr);
    case _CPLX:
      return 0;
    case _DOUBLE_:
      if (a._DOUBLE_val>0)
	return 1;
      else
	return -1;
    case _REAL:
      if (a._REALptr->is_positive())
	return 1;
      else
	return -1;
    case _SYMB:
      if (a._SYMBptr->sommet==at_abs || (a._SYMBptr->sommet==at_exp && is_real(a._SYMBptr->feuille,contextptr)))
	return 1;
    }
    if (a.type==_SYMB && a.is_symb_of_sommet(at_pow)){
      gen & f =a._SYMBptr->feuille;
      if (f.type==_VECT && f._VECTptr->size()==2){
	gen & ex = f._VECTptr->back();
	if (ex.type==_FRAC && ex._FRACptr->den.type==_INT_ && ex._FRACptr->den.val % 2 ==0 ) 
	  return 1;
      }
    }
    if (is_inf(a)){
      if (a==plus_inf)
	return 1;
      if (a==minus_inf)
	return -1;
      return 0;
    }
    gen approx;
    if (has_evalf(a,approx,1,contextptr) && a!=approx)
      return fastsign(approx,contextptr);
    // FIXME GIAC_CONTEXT??
    /*
    if (contextptr){
      gen test=superieur_strict(a,0,contextptr);
      if (test.type==_INT_){
	if (test.val)
	  return test.val;
	test=inferieur_strict(a,0,contextptr);
	if (test.type==_INT_)
	  return -test.val;
      }
    }
    */
    return 0;
  }

  bool is_greater(const gen & a,const gen &b,GIAC_CONTEXT){
    gen test=superieur_egal(a,b,contextptr);
    if ((test.type==_INT_) && (test.val==1))
      return true;
    else
      return false;    
  }

  bool is_strictly_greater(const gen & a,const gen &b,GIAC_CONTEXT){
    gen test=superieur_strict(a,b,contextptr);
    if ((test.type==_INT_) && (test.val==1))
      return true;
    else
      return false;    
  }

  bool real_object::is_positive(){
#ifdef HAVE_LIBMPFR
    if (mpfr_sgn(inf)==-1)
#else
    if (mpf_sgn(inf)==-1)
#endif
      return false;
    else
      return true;
  }

  bool is_positive(const gen & a,GIAC_CONTEXT){
    switch (a.type){
    case _INT_:
      return a.val>=0;
    case _REAL:
      return a._REALptr->is_positive();
    case _ZINT:
      if (mpz_sgn(*a._ZINTptr)==-1)
	return false;
      else
	return true;
    case _POLY:
      return is_positive(a._POLYptr->coord.front());
    case _FRAC:
      return (is_positive(a._FRACptr->num,contextptr) && is_positive(a._FRACptr->den,contextptr)) || (is_positive(-a._FRACptr->num,contextptr) && is_positive(-a._FRACptr->den,contextptr));
    case _EXT:
      return false;
    case _SYMB:
      if (a==plus_inf)
	return true;
      if (a==minus_inf)
	return false;      
      if (a._SYMBptr->sommet==at_exp)
	return true;
      if (a._SYMBptr->sommet==at_ln)
	return is_positive(a._SYMBptr->feuille-1,contextptr);
      if (a._SYMBptr->sommet==at_program)
	return true;
      return is_greater(a,0,contextptr); 
    case _FUNC:
      return true;
    default:
      return is_greater(a,0,contextptr); 
    }
  }

  bool is_strictly_positive(const gen & a,GIAC_CONTEXT){
    if (is_zero(a))
      return false;
    return is_positive(a,contextptr);
  }

  bool ck_is_greater(const gen & a,const gen &b,GIAC_CONTEXT){
    if (a==b)
      return true;
    gen test=superieur_strict(a,b,contextptr);
    if (test.type!=_INT_)
      cksignerr(test);
    if (test.val==1)
      return true;
    else
      return false;    
  }

  bool ck_is_strictly_greater(const gen & a,const gen &b,GIAC_CONTEXT){
    gen test=superieur_strict(a,b,contextptr);
    if (test.type!=_INT_)
      cksignerr(test);
    if (test.val==1)
      return true;
    else
      return false;    
  }

  bool ck_is_positive(const gen & a,GIAC_CONTEXT){
    switch (a.type){
    case _INT_:
      return a.val>=0;
    case _ZINT:
      if (mpz_sgn(*a._ZINTptr)==-1)
	return false;
      else
	return true;
    case _SYMB:
      if (a==plus_inf)
	return true;
      if (a==minus_inf)
	return false;
      if (a._SYMBptr->sommet==at_exp)
	return true;
      if (a._SYMBptr->sommet==at_ln)
	return ck_is_positive(a._SYMBptr->feuille-1,contextptr);    
      return ck_is_greater(a,0,contextptr);
    default:
      return ck_is_greater(a,0,contextptr);
    }
  }

  bool ck_is_strictly_positive(const gen & a,GIAC_CONTEXT){
    if (is_zero(a))
      return false;
    return ck_is_positive(a,contextptr);
  }

  gen min(const gen & a, const gen & b,GIAC_CONTEXT){
    if (a==b)
      return a;
    if (is_inf(a)){
      if (a==plus_inf)
	return b;
      if (a==minus_inf)
	return a;
      if (!is_inf(b))
	return undef;
    }
    if (is_inf(b)){
      if (b==plus_inf)
	return a;
      if (b==minus_inf)
	return b;
      return undef;
    }
    if (is_undef(a) || is_undef(b))
      return undef;
    gen test=superieur_strict(a,b,contextptr);
    if (test.type==_INT_){
      if (test.val==1)
	return b;
      else
	return a;
    }
    return new symbolic(at_min,makenewvecteur(a,b));
  }

  gen max(const gen & a, const gen & b,GIAC_CONTEXT){
    if (a==b)
      return a;
    if (is_inf(a)){
      if (a==plus_inf)
	return a;
      if (a==minus_inf)
	return b;
      if (!is_inf(b))
	return undef;
    }
    if (is_inf(b)){
      if (b==plus_inf)
	return b;
      if (b==minus_inf)
	return a;
      return undef;
    }
    if (is_undef(a) || is_undef(b))
      return undef;
    gen test=superieur_strict(a,b,contextptr);
    if (test.type==_INT_){
      if (test.val==1)
	return a;
      else
	return b;
    }
    return new symbolic(at_max,makenewvecteur(a,b));
  }

  bool has_evalf(const identificateur & g,int subtype,gen & res,int level,GIAC_CONTEXT){
    if (g.name==_IDNT_pi.name){
      res=m_pi(contextptr);
      return true;
    }
    gen tmp=g;
    tmp.subtype=subtype;
    tmp=tmp.evalf(level,contextptr);
    if (tmp.type==_IDNT || tmp.type==_SYMB)
      return false;
    return has_evalf(tmp,res,0,contextptr);
  }

  bool has_evalf(const gen & g,gen & res,int level,GIAC_CONTEXT){
    switch (g.type){
    case _DOUBLE_: case _REAL:
      res=g;
      return true;
    case _INT_: case _ZINT: case _CPLX:
      res=evalf(g,1,contextptr);
      return true;
    case _IDNT:
      return has_evalf(*g._IDNTptr,g.subtype,res,level,contextptr);
    case _SYMB:
      if (has_evalf(g._SYMBptr->feuille,res,level,contextptr)){
	res=g._SYMBptr->sommet(res,contextptr); 
	if (res.type==_INT_ || res.type==_ZINT)
	  res=evalf(res,1,contextptr);
	return res.type==_DOUBLE_ || res.type==_CPLX || res.type==_REAL;
      }
      else
	return false;
    }
    if (g.type==_EXT){
      gen a,b;
      if (has_evalf(*g._EXTptr,a,level,contextptr) && has_evalf(*(g._EXTptr+1),b,level,contextptr)){
	a=alg_evalf(a,b,contextptr);
	return has_evalf(a,res,level,contextptr);
      }
      return false;
    }
    if (g.type==_FRAC){
      gen num,den;
      if (has_evalf(g._FRACptr->num,num,level,contextptr) && has_evalf(g._FRACptr->den,den,level,contextptr)){
	res=num/den;
	return true;
      }
      else
	return false;
    }
    if (g.type!=_VECT)
      return false;
    if (g.subtype==_ASSUME__VECT && !g._VECTptr->empty()){
      res=g._VECTptr->back();
      return true;
    }
    vecteur v;
    const_iterateur it=g._VECTptr->begin(),itend=g._VECTptr->end();
    v.reserve(itend-it);
    for (;it!=itend;++it){
      if (!has_evalf(*it,res,level,contextptr))
	return false;
      v.push_back(res);
    }
    res=gen(v,g.subtype);
    return true;
  }

  gen sym_mult (const gen & a,const gen & b,GIAC_CONTEXT){
    control_c();
    if ((a==undef) || (b==undef))
      return undef;
    if (a.is_symb_of_sommet(at_unit)){
      vecteur & va=*a._SYMBptr->feuille._VECTptr;
      if (b.is_symb_of_sommet(at_unit)){
	vecteur & v=*b._SYMBptr->feuille._VECTptr;
	return new symbolic(at_unit,makenewvecteur(va[0]*v[0],_simplifier(va[1]*v[1],contextptr)));
      }
      else
	return new symbolic(at_unit,makenewvecteur(va[0]*b,va[1]));
    }
    if (b.is_symb_of_sommet(at_unit)){
      vecteur & v=*b._SYMBptr->feuille._VECTptr;
      return new symbolic(at_unit,makenewvecteur(a*v[0],v[1]));
    }
    if (is_inf(a)){
      if (is_exactly_zero(normal(b,contextptr)))
	return undef;
      int s=fastsign(a,contextptr)*fastsign(b,contextptr); 
      if (s==1)
	return plus_inf;
      if (s)
	return minus_inf;
      return unsigned_inf;
    }
    if (is_inf(b)){
      if (is_exactly_zero(normal(a,contextptr)))
	return undef;
      int s=fastsign(a,contextptr)*fastsign(b,contextptr); 
      if (s==1)
	return plus_inf;
      if (s)
	return minus_inf;
      return unsigned_inf;
    }
    if (a.type==_INT_ && a.val==0 )
      return a;
    if (a.type==_DOUBLE_ && a._DOUBLE_val==0 )
      return a;
    if (b.type==_INT_ && b.val==0)
      return b;
    if (b.type==_DOUBLE_ && b._DOUBLE_val==0 )
      return b;
    if (is_one(a) && ((a.type!=_MOD) || (b.type==_MOD) ))
      return b;
    if (is_one(b) && ((b.type!=_MOD) || (a.type==_MOD) ))
      return a;
    if ( a.is_approx()){
      gen b1;
      if (has_evalf(b,b1,1,contextptr)&& b!=b1)
	return a*b1;
    }
    if ( b.is_approx()){
      gen a1;
      if (has_evalf(a,a1,1,contextptr) && a!=a1)
	return a1*b;
    }
    if ((a.type==_SYMB) && equalposcomp(plot_sommets,a._SYMBptr->sommet)){
      gen tmp=remove_at_pnt(a);
      if (tmp.type==_VECT && tmp.subtype==_VECTOR__VECT){
	if (b.type==_SYMB && equalposcomp(plot_sommets,b._SYMBptr->sommet)){
	  gen tmpb=remove_at_pnt(b);
	  return dotvecteur(vector2vecteur(*tmp._VECTptr),vector2vecteur(*tmpb._VECTptr));
	}
	return _vector(vector2vecteur(*tmp._VECTptr)*b,contextptr);
      }
      return symbolic_plot_makevecteur(a._SYMBptr->sommet,a._SYMBptr->feuille*b,false,contextptr);
    }
    if ((b.type==_SYMB) && equalposcomp(plot_sommets,b._SYMBptr->sommet)){
      gen tmp=remove_at_pnt(b);
      if (tmp.type==_VECT && tmp.subtype==_VECTOR__VECT)
	return _vector(a*vector2vecteur(*tmp._VECTptr),contextptr);
      return symbolic_plot_makevecteur(b._SYMBptr->sommet,b._SYMBptr->feuille*a,false,contextptr);
    }
    if (a.type==_FRAC){
      if ( (b.type!=_SYMB) && (b.type!=_IDNT) )      
        return (*a._FRACptr)*b;
      return sym_mult(_FRAC2_SYMB(a),b,contextptr);
    }
    if (b.type==_FRAC){
      if ( (a.type!=_SYMB) && (a.type!=_IDNT) )
        return a*(*b._FRACptr);
      return sym_mult(a,_FRAC2_SYMB(b),contextptr);
    }
    if (a.type==_EXT){
        if (a.is_constant() && (b.type==_POLY))
            return a*(*b._POLYptr);
        else
            return algebraic_EXTension(*a._EXTptr*b,*(a._EXTptr+1));
    }
    if (b.type==_EXT){
        if (b.is_constant() && (a.type==_POLY))
            return (*a._POLYptr)*b;
        else
            return algebraic_EXTension(a*(*b._EXTptr),*(b._EXTptr+1));
    }
    if ( (a.type==_INT_) && (a.val<0) && (a.val!=1<<31)){
      if (b.is_symb_of_sommet(at_inv))
	return sym_mult(-a,inv(-b._SYMBptr->feuille,contextptr),contextptr);
      else
	return -sym_mult(-a,b,contextptr);
    }
    if ( (b.type==_INT_) && (b.val<0) && (b.val!=1<<31)){
      if (a.is_symb_of_sommet(at_inv))
	return sym_mult(-b,inv(-a._SYMBptr->feuille,contextptr),contextptr);
      else
	return -sym_mult(-b,a,contextptr);
    }
    if (a.type==_SYMB && a._SYMBptr->sommet==at_equal){
      vecteur & va=*a._SYMBptr->feuille._VECTptr;
      if (b.type==_SYMB && b._SYMBptr->sommet==at_equal){
	vecteur & vb=*b._SYMBptr->feuille._VECTptr;
	return new symbolic(at_equal,makenewvecteur(va.front()*vb.front(),va.back()*vb.back()));
      }
      else
	return new symbolic(at_equal,makenewvecteur(va.front()*b,va.back()*b));
    }
    if (b.type==_SYMB && b._SYMBptr->sommet==at_equal){
      vecteur & vb=*b._SYMBptr->feuille._VECTptr;
      return new symbolic(at_equal,makenewvecteur(a*vb.front(),a*vb.back()));
    }
    if ((a.type==_SYMB)&& (b.type==_SYMB)){
      if ((a._SYMBptr->sommet==at_prod) && (b._SYMBptr->sommet==at_prod))
	return new symbolic(at_prod,mergevecteur(*(a._SYMBptr->feuille._VECTptr),*(b._SYMBptr->feuille._VECTptr)));
      else {
	if (a._SYMBptr->sommet==at_prod)
	  return new symbolic(*a._SYMBptr,b);
	else {
	  if (b._SYMBptr->sommet==at_prod)
	    return new symbolic(a,b._SYMBptr->sommet,b._SYMBptr->feuille);
	  else
	    return new symbolic(at_prod,gen(makenewvecteur(a,b)));
	}
      }
    }
    if (b.type==_SYMB){
      if (b._SYMBptr->sommet==at_prod)
	return new symbolic(a,b._SYMBptr->sommet,b._SYMBptr->feuille);
      else
	return new symbolic(at_prod,gen(makenewvecteur(a,b)));
    }
    if (a.type==_SYMB){
      if (a._SYMBptr->sommet==at_prod)
	return new symbolic(*a._SYMBptr,b);
      else
	return new symbolic(at_prod,gen(makenewvecteur(a,b)));
    }
    if ((a.type==_IDNT) || (b.type==_IDNT))
      return new symbolic(at_prod,gen(makenewvecteur(a,b)));
    if (a.type==_MOD)
      return a*makemod(b,*(a._MODptr+1));
    if (b.type==_MOD)
      return b*makemod(a,*(b._MODptr+1));
    return new symbolic(at_prod,makenewvecteur(a,b));
    // settypeerr("sym_mult");
  }

  vecteur inv__VECT(const vecteur & v,GIAC_CONTEXT){
    vecteur w;
    if (is_squarematrix(v))
      w=minv(v,contextptr);
    else {
      vecteur::const_iterator it=v.begin(),itend=v.end();
      for (;it!=itend;++it)
	w.push_back(inv(*it,contextptr));
    }
    return w;
  }

  vecteur invfirst(const vecteur & v){
    vecteur w(v);
    if (!w.empty())
      w.front()=inv(w.front(),context0);
    return w;
  }

  gen inv(const gen & a,GIAC_CONTEXT){
    if ( a.type==_DOUBLE_ ? a==0 : is_exactly_zero(a))
      return unsigned_inf;
    switch (a.type ) {
    case _INT_: case _ZINT:
      if (is_one(a) || (is_minus_one(a)) )
	return a;
      else
	return fraction(1,a);
    case _REAL:
      return a._REALptr->inv();
    case _DOUBLE_:
      return 1/a._DOUBLE_val;
    case _CPLX:
      if (is_exactly_zero(*a._CPLXptr)){
	if (is_one(abs(*(a._CPLXptr+1),contextptr)))
	  return -a;
      }
      if ( a._CPLXptr->type==_DOUBLE_ || a._CPLXptr->type==_REAL || (a._CPLXptr+1)->type==_DOUBLE_ || (a._CPLXptr+1)->type==_REAL )
	return gen(rdiv(no_context_evalf(a.re(contextptr)),no_context_evalf(a.squarenorm(contextptr))),rdiv(no_context_evalf(-a.im(contextptr)),no_context_evalf(a.squarenorm(contextptr))));
      return fraction(1,a);
    case _IDNT:
      if (a==undef)
	return undef;
      if (a==unsigned_inf)
	return 0;
      return new symbolic(at_inv,a);
    case _SYMB:
      if ((a==plus_inf) || (a==minus_inf))
	return 0;
      if (a.is_symb_of_sommet(at_unit))
	return new symbolic(at_unit,makenewvecteur(inv(a._SYMBptr->feuille._VECTptr->front(),contextptr),_simplifier(inv(a._SYMBptr->feuille._VECTptr->back(),contextptr),contextptr)));
      if (equalposcomp(plot_sommets,a._SYMBptr->sommet))
	return symbolic_plot_makevecteur( a._SYMBptr->sommet,inv(a._SYMBptr->feuille,contextptr),false,contextptr);
      if (a._SYMBptr->sommet==at_inv)
	return a._SYMBptr->feuille;
      if (a._SYMBptr->sommet==at_neg)
	return -inv(a._SYMBptr->feuille,contextptr);      
      else {
	if (a._SYMBptr->sommet==at_prod)
	  return new symbolic(at_prod,inv__VECT(*(a._SYMBptr->feuille._VECTptr),contextptr));
	else
	  return new symbolic(at_inv,a);
      }
    case _VECT:
      if (a.subtype==_PNT__VECT)
	return gen(invfirst(*a._VECTptr),a.subtype);
      if (a.subtype==_POLY1__VECT)
	return fraction(gen(vecteur(1,plus_one),_POLY1__VECT),a);
      return gen(inv__VECT(*a._VECTptr,contextptr),a.subtype);
    case _EXT:
      return inv_EXT(a);
    case _USER:
      return a._USERptr->inv();
    case _MOD:
      return modinv(a);
    case _FRAC:
      if (a._FRACptr->num.type==_CPLX)
	return fraction(a._FRACptr->den,a._FRACptr->num).normal();
      return fraction(a._FRACptr->den,a._FRACptr->num);
    default: 
      return new symbolic(at_inv,a);
      // settypeerr("Inv");
    }
    
  }

  /*
  gen inv(const gen & a,GIAC_CONTEXT){
    return inv(a,context0);
  }
  */
  
  gen gen::inverse(GIAC_CONTEXT) const  { return inv(*this,contextptr); }

  void inpow(const gen & base,unsigned long int exponent,gen & res){
    if (exponent==1)
      res=base;
    else {
      inpow(base,exponent/2,res);
      res=res*res;
      if (exponent %2)
	res=res*base;
    }
  }

  gen pow(const gen & base, unsigned long int exponent){
    // if (!( (++control_c_counter) & control_c_counter_mask))
      control_c();
    mpz_t * e;
    gen res;
    switch (base.type ) {
    case _INT_: 
      if (base.val<0 && (exponent % 2))
	return(-pow(-base.val,exponent));
      else
	return(pow(absint(base.val),exponent));
    case _DOUBLE_:
#ifdef _SOFTMATH_H
      return std::giac_gnuwince_pow(base._DOUBLE_val,double(exponent));
#else
      return std::pow(base._DOUBLE_val,double(exponent));
#endif
    case _ZINT: 
      e=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init(*e);
      mpz_pow_ui(*e,*base._ZINTptr,exponent);
      return(e);
    case _CPLX: case _REAL: case _EXT: case _VECT: case _MOD: case _USER:
      // gauss integer power
      if (!exponent){
	if (ckmatrix(base))
	  return midn(base._VECTptr->size());
	return 1;
      }
      inpow(base,exponent,res);
      return(res);
    case _IDNT:
      if (is_undef(base))
	return base;
      if (!exponent)
	return 1;
      if (exponent==1)
	return base;
      return new symbolic(at_pow,makenewvecteur(base,(longlong)exponent));
    case _SYMB:
      if (!exponent)
	return 1;
      if (exponent==1)
	return base;
      if (base._SYMBptr->sommet==at_pow){
	res= (*((base._SYMBptr->feuille)._VECTptr))[1];
	return pow( (base._SYMBptr->feuille)._VECTptr->front(),gen((longlong) (exponent)) * res,context0) ;
      }
      return new symbolic(at_pow,makenewvecteur(base,(longlong) exponent));
    case _POLY:
      return pow(*base._POLYptr,(int) exponent);
    case _FRAC:
      return pow(*base._FRACptr,(int) exponent);
    default: 
      settypeerr("Pow") ;
    }
    return 0;
  }

  gen pow(const gen & base, int exponent){
    if (base==zero){
      if (exponent>0)
	return base;
      if (!exponent)
	return undef;
      if (exponent %2)
	return unsigned_inf;
      return plus_inf;
    }
    if (exponent<0)
      return inv(pow(base,-exponent),context0);
    unsigned long int expo=exponent;
    return(pow(base,expo));
  }

  gen pow(unsigned long int base, unsigned long int exponent){
    mpz_t *e=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
    mpz_init(*e);
    mpz_ui_pow_ui(*e,base,exponent);
    return(e);
  }

  void _ZINTdiv (const gen & a,const gen & b,mpz_t * & quo){
    // at least one is not an int, uncoerce remaining int
    mpz_t *aptr,*bptr;
    if (a.type!=_INT_)
      aptr=a._ZINTptr;
    else {
      aptr=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init_set_si(*aptr,a.val);
    }
    if (b.type!=_INT_)
      bptr=b._ZINTptr;
    else {
      bptr=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init_set_si(*bptr,b.val);
    }
    quo=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
    mpz_init(*quo);
    mpz_tdiv_q(*quo,*aptr,*bptr);
    if (a.type==_INT_){
      mpz_clear(*aptr);
      free(aptr);
    }
    if (b.type==_INT_){
      mpz_clear(*bptr);
      free(bptr);
    }
  }

  // a and b must be integers or Gaussian integers
  gen iquobest(const gen & a,const gen & b){
    if (is_strictly_positive(-a,0))
      return -iquobest(-a,b);
    return iquo(a+iquo(b,2),b);
  }

  // a and b must be integers or Gaussian integers
  gen iquocmplx(const gen & a,const gen & b){
    gen b2=b.squarenorm(0);
    gen ab=a*b.conj(0);
    gen res(iquobest(re(ab,context0),b2),iquobest(im(ab,context0),b2)); // ok
    return res;
  }

  // integer quotient, use rdiv for symbolic division 
  gen iquo(const gen & a,const gen & b){
    if ((b.type==_INT_)){
      switch (b.val){
      case 1:
	return a;
      case -1:
	return -a;
      case 0:
	setsizeerr("Division by 0");
      }
    }
    mpz_t * quo;
    switch ( (a.type<< _DECALAGE) | b.type ) {
    case _INT___INT_: 
      return(a.val/b.val);
    case _ZINT__ZINT: case _INT___ZINT: case _ZINT__INT_:
      _ZINTdiv(a,b,quo);
      return(quo);
    case _CPLX__INT_:  case _CPLX__ZINT:
      return gen(iquo(*a._CPLXptr,b),iquo(*(a._CPLXptr+1),b));
    case _INT___CPLX: case _ZINT__CPLX: case _CPLX__CPLX:
      return iquocmplx(a,b);
    default:
      settypeerr("iquo");
    }
    return 0;
  }

  // a and b must be integer or Gaussian integers
  gen rdivsimp(const gen & a,const gen & b){
    if (is_positive(-b,context0)) // ok
      return rdivsimp(-a,-b);
    gen c(gcd(a,b));
    if (c.type==_CPLX)
      c=gcd(c.re(context0),c.im(context0)); // ok
    return fraction(iquo(a,c),iquo(b,c));
  }

  gen divpoly(const polynome & p, const gen & e){
    if (p.coord.empty())
      return zero;
    gen d=gcd(Tcontent<gen>(p),e);
    if (is_one(d))
      return fraction(p,e);
    gen den(rdiv(e,d));
    gen iden(inv(den,context0));
    if ( (iden.type!=_SYMB) && (iden.type!=_FRAC))
      return (p/d)*iden;
    else
      return fraction(p/d,den);
  }

  gen divpoly(const gen & e,const polynome & p){
    if (is_zero(e))
      return e;
    if (Tis_constant<gen>(p)&& p.coord.front().value.type<_POLY)
      return rdiv(e,p.coord.front().value);
    gen d=gcd(Tcontent<gen>(p),e);
    gen tmp=polynome(rdiv(e,d),p.dim);
    return fraction(tmp,p/d);
  }
  
  gen divpolypoly(const gen & a,const gen &b){
    polynome ap(*a._POLYptr),bp(*b._POLYptr);
    polynome q(ap.dim),r(ap.dim);
    if (divrem1(ap,bp,q,r) && r.coord.empty())
      return q;
    return normal(fraction(a,b),context0); // ok
  }

  gen rdiv(const gen &a,const gen &b){
    // if (!( (++control_c_counter) & control_c_counter_mask))
      control_c();
    if (b.type==_MOD)
      return a*inv(b,context0);
    switch ( (a.type<< _DECALAGE) | b.type ) {
    case _INT___INT_: case _ZINT__INT_: case _ZINT__ZINT:
      if (is_exactly_zero(b)){
	if (is_exactly_zero(a))
	  return undef;
	return unsigned_inf;
      }
      if (is_exactly_zero(a%b))
	return iquo(a,b);
      else
	return rdivsimp(a,b);
    case _CPLX__INT_: case _CPLX__ZINT:
      if ( (a._CPLXptr->type==_DOUBLE_) || ((a._CPLXptr+1)->type==_DOUBLE_))
	return rdiv(no_context_evalf(a),no_context_evalf(b));
      if (a._CPLXptr->type==_REAL)
#ifdef HAVE_LIBMPFR
	return rdiv(*a._CPLXptr,b)+cst_i*rdiv(real_object(*(a._CPLXptr+1),mpfr_get_prec(a._CPLXptr->_REALptr->inf)),b);
#else
	return rdiv(*a._CPLXptr,b)+cst_i*rdiv(real_object(*(a._CPLXptr+1)),b);
#endif
      if ((a._CPLXptr+1)->type==_REAL)
#ifdef HAVE_LIBMPFR
	return rdiv(real_object(*a._CPLXptr,mpfr_get_prec((a._CPLXptr+1)->_REALptr->inf)),b)+cst_i*rdiv(*(a._CPLXptr+1),b);
#else
	return rdiv(real_object(*a._CPLXptr),b)+cst_i*rdiv(*(a._CPLXptr+1),b);
#endif
      if (is_exactly_zero(b))
	return unsigned_inf;
      if (is_exactly_zero(a%b))
	return iquo(a,b);
      else
	return rdivsimp(a,b);
    case _CPLX__CPLX: case _DOUBLE___CPLX: case _INT___CPLX: case _ZINT__CPLX: case _REAL__CPLX:
      return rdiv(a*conj(b,context0),b.squarenorm(0));
    case _DOUBLE___DOUBLE_:
      return a._DOUBLE_val/b._DOUBLE_val;
    case _DOUBLE___INT_:
      return a._DOUBLE_val/b.val;
    case _INT___DOUBLE_:
      return a.val/b._DOUBLE_val;
    case _ZINT__DOUBLE_:
      return mpz_get_d(*a._ZINTptr)/b._DOUBLE_val;
    case _CPLX__DOUBLE_: case _CPLX__REAL:
      return gen(rdiv(*a._CPLXptr,b),rdiv(*(a._CPLXptr+1),b));
    case _DOUBLE___ZINT:
      return a._DOUBLE_val/mpz_get_d(*b._ZINTptr);
      // _CPLX__DOUBLE_, _DOUBLE___CPLX, _CPLX__CPLX, _ZINT__CPLX, _INT___CPLX
    case _VECT__INT_: case _VECT__ZINT: case _VECT__DOUBLE_: case _VECT__CPLX: 
    case _VECT__SYMB: case _VECT__IDNT: case _VECT__POLY: case _VECT__EXT:
      if (a.subtype==_VECTOR__VECT)
	return a*inv(b,context0);
      return gen(divvecteur(*a._VECTptr,b),a.subtype);
    case _VECT__VECT:
      if (a.subtype==_POLY1__VECT || b.subtype==_POLY1__VECT)
	return fraction(a,b).normal();
      if (b._VECTptr->size()==1)
	return rdiv(a,b._VECTptr->front());
      return apply(a,b,rdiv);
    case _POLY__POLY:
      return divpolypoly(a,b);
    case _FRAC__FRAC:
      if (a._FRACptr->num.type==_CPLX || a._FRACptr->den.type==_CPLX ||
	  b._FRACptr->num.type==_CPLX || b._FRACptr->den.type==_CPLX)
	return (a._FRACptr->num*b._FRACptr->den)/(a._FRACptr->den*b._FRACptr->num);
      return (*a._FRACptr)/(*b._FRACptr);
    case _SPOL1__SPOL1:
      return spdiv(*a._SPOL1ptr,*b._SPOL1ptr,context0);
    case _POLY__DOUBLE_: case _POLY__REAL:
      return (*a._POLYptr)/b;
    case _POLY__INT_: case _POLY__ZINT: case _POLY__CPLX:
      return divpoly(*a._POLYptr,b);
    case _INT___POLY: case _ZINT__POLY: case _CPLX__POLY:
      return divpoly(a,*b._POLYptr);
    case _INT___VECT: case _ZINT__VECT: case _CPLX__VECT: case _SYMB__VECT:
      if (ckmatrix(b))
	return a*inv(b,context0);
      else
	return fraction(a,b);
    default:
      control_c();
      if ((a==undef) || (b==undef))
	return undef;
      if (a.is_symb_of_sommet(at_unit) || b.is_symb_of_sommet(at_unit))
	return a*inv(b,context0);
      if (is_one(b))
	return chkmod(a,b);
      if (is_minus_one(b))
	return chkmod(-a,b);
      if (is_exactly_zero(a)){
	if (!is_exactly_zero(normal(b,context0)))
	  return a;
	else
	  return undef;
      }
      if (is_inf(b)){
	if (is_inf(a))
	  return undef;
	else
	  return zero;
      }
      if (a==b)
	return chkmod(plus_one,a);
      if (a.is_approx()){
	gen b1;
	if (has_evalf(b,b1,1,context0) && b!=b1)
	  return rdiv(a,b1);
      }
      if (b.is_approx()){
	gen a1;
	if (has_evalf(a,a1,1,context0) && a!=a1)
	  return rdiv(a1,b);
      }
      if (a.type==_REAL)
	return (*a._REALptr)*inv(b,context0);
      if (b.type==_REAL)
	return a*b._REALptr->inv();
      if (a.type==_USER || b.type==_USER)
	return a*inv(b,context0);
      if (a.type==_FRAC){
	if ( (b.type!=_SYMB) && (b.type!=_IDNT) )
	  return (*a._FRACptr)/b;
	return rdiv(_FRAC2_SYMB(a),b);
      }
      if (b.type==_FRAC){
	if ( (a.type!=_SYMB) && (a.type!=_IDNT) )
	  return a/(*b._FRACptr);
	return rdiv(a,_FRAC2_SYMB(b));
      }
      if (a.type==_SYMB && a._SYMBptr->sommet==at_equal){
	vecteur & va=*a._SYMBptr->feuille._VECTptr;
	if (b.type==_SYMB && b._SYMBptr->sommet==at_equal){
	  vecteur & vb=*b._SYMBptr->feuille._VECTptr;
	  return new symbolic(at_equal,makenewvecteur(rdiv(va.front(),vb.front()),rdiv(va.back(),vb.back())));
	}
	else
	  return new symbolic(at_equal,makenewvecteur(rdiv(va.front(),b),rdiv(va.back(),b)));
      }
      if (b.type==_SYMB && b._SYMBptr->sommet==at_equal){
	vecteur & vb=*b._SYMBptr->feuille._VECTptr;
	return new symbolic(at_equal,makenewvecteur(rdiv(a,vb.front()),rdiv(a,vb.back())));
      }
      /* commented since * is not always commutative
      if (a.is_symb_of_sommet(at_prod) && a._SYMBptr->feuille.type==_VECT){
	int i=equalposcomp(*a._SYMBptr->feuille._VECTptr,b);
	if (i){
	  vecteur v(*a._SYMBptr->feuille._VECTptr);
	  v.erase(v.begin()+i-1);
	  if (v.size()==1)
	    return v.front();
	  else
	    return new symbolic(at_prod,v);
	}
      }
      */
      if ( (a.type==_SYMB) || (a.type==_IDNT) || (a.type==_FUNC) || (b.type==_SYMB) || (b.type==_IDNT) || (b.type==_FUNC) )
	return symb_prod(a,symb_inv(b));
      if (a.type==_STRNG || b.type==_STRNG)
	settypeerr();
      return fraction(a,b).normal();
    }
  }

  void _ZINTmod (const gen & a,const gen & b,mpz_t * & rem){
    // at least one is not an int, uncoerce remaining int
    mpz_t *aptr,*bptr;
    if (a.type!=_INT_)
      aptr=a._ZINTptr;
    else {
      aptr=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init_set_si(*aptr,a.val);
    }
    if (b.type!=_INT_)
      bptr=b._ZINTptr;
    else {
      bptr=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init_set_si(*bptr,b.val);
    }
    rem=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
    mpz_init(*rem);
    mpz_tdiv_r(*rem,*aptr,*bptr);
    if (a.type==_INT_){
      mpz_clear(*aptr);
      free(aptr);
    }
    if (b.type==_INT_){
      mpz_clear(*bptr);
      free(bptr);
    }
  }

  gen operator %(const gen & a,const gen & b){
    mpz_t * rem;
    switch ( (a.type<< _DECALAGE) | b.type ) {
    case _INT___INT_: 
      if (b.val)
	return(a.val % b.val);
      else
	return a.val;
    case _ZINT__ZINT: case _INT___ZINT: case _ZINT__INT_:
      _ZINTmod(a,b,rem);
      return(rem);
    case _CPLX__INT_: case _CPLX__ZINT:
      return gen(smod((*a._CPLXptr), b), smod(*(a._CPLXptr+1), b) );
    case _INT___CPLX: case _ZINT__CPLX: case _CPLX__CPLX:   
      return(a-b*iquo(a,b));
    case _VECT__VECT:
      return (*a._VECTptr)%(*b._VECTptr);
    default:
      settypeerr("%");
    }
    return 0;
  }

  void _ZINTrem(const gen & a,const gen &b,gen & q,mpz_t * & rem){
    // cout << a << " irem " << b << endl;
    mpz_t *aptr,*bptr;
    if (a.type!=_INT_)
      aptr=a._ZINTptr;
    else {
      aptr=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init_set_si(*aptr,a.val);
    }
    if (b.type!=_INT_)
      bptr=b._ZINTptr;
    else {
      bptr=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init_set_si(*bptr,b.val);
    }
    rem=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
    mpz_init(*rem);
    q.uncoerce();
    mpz_tdiv_qr(*q._ZINTptr,*rem,*aptr,*bptr);
    if (a.type==_INT_){
      mpz_clear(*aptr);
      free(aptr);
    }
    if (b.type==_INT_){
      mpz_clear(*bptr);
      free(bptr);
    }
  }

  gen irem(const gen & a,const gen & b,gen & q){
    mpz_t * rem;
    register int r;
    switch ( (a.type<< _DECALAGE) | b.type ) {
    case _INT___INT_: 
      if (!b.val)
	return a;
      r=a.val % b.val;
      /*
      if (r<0){ 
	if (b.val>0){
	  q=gen(a.val/b.val-1);
	  r += b.val;
	}
	else {
	  q=gen(a.val/b.val+1);
	  r -= b.val;
	}
      }
      else
      */
      q=gen(a.val/b.val);
      return r;
    case _ZINT__ZINT: case _INT___ZINT: case _ZINT__INT_:
      _ZINTrem(a,b,q,rem);
      return(rem);
    case _INT___CPLX: case _ZINT__CPLX: case _CPLX__CPLX: case _CPLX__INT_: case _CPLX__ZINT:
      q=iquo(a,b);
      return(a-b*q);
    default:
      settypeerr("irem");
    }
    return 0;
  }

  void _ZINTsmod(const gen & a, const gen & b, mpz_t * & rem){
    // at least one is not an int, uncoerce remaining int
    mpz_t *aptr,*bptr;
    if (a.type!=_INT_)
      aptr=a._ZINTptr;
    else {
      aptr=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init_set_si(*aptr,a.val);
    }
    if (b.type!=_INT_)
      bptr=b._ZINTptr;
    else {
      bptr=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init_set_si(*bptr,b.val);
    }
    rem=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
    mpz_t rem1,rem2,rem3;
    mpz_init(rem1); mpz_init(rem2); mpz_init(rem3);
    mpz_mod(rem1,*aptr,*bptr); // rem1 positive remainder
    if (mpz_sgn(*bptr)>0)
      mpz_sub(rem2,rem1,*bptr); // negative remainder
    else
      mpz_add(rem2,rem1,*bptr);
    // choose smallest one in abs value
    mpz_neg(rem3,rem2);
    if (mpz_cmp(rem1,rem3)>0)
      mpz_init_set(*rem,rem2);
    else
      mpz_init_set(*rem,rem1);
    if (a.type==_INT_){
      mpz_clear(*aptr);
      free(aptr);
    }
    if (b.type==_INT_){
      mpz_clear(*bptr);
      free(bptr);
    }
    mpz_clear(rem1); mpz_clear(rem2); mpz_clear(rem3); 
  }

  vecteur smod(const vecteur & v,const gen & g){
    vecteur w(v);
    iterateur it=w.begin(),itend=w.end();
    for (;it!=itend;++it)
      *it=smod(*it,g);
    return w;
  }

  gen smodSYMB(const gen & a,const gen & b){
    vecteur lv(lvar(a));
    gen n,d,f;
    f=e2r(a,lv,context0); // ok
    fxnd(f,n,d);
    n=smod(n,b);
    d=smod(d,b);
    f=n/d;
    return r2e(f,lv,context0); // ok
  }

  gen smod(const gen & a,const gen & b){
    if (b==0)
      return a;
    mpz_t * rem;
    switch ( (a.type<< _DECALAGE) | b.type ) {
    case _INT___INT_: 
      return smod(a.val,b.val);
    case _ZINT__INT_: case _INT___ZINT: case _ZINT__ZINT: 
      _ZINTsmod(a,b,rem);
      return(rem);
    case _CPLX__INT_: case _CPLX__ZINT:
      return gen(smod(*a._CPLXptr,b),smod(*(a._CPLXptr+1),b));
    case _POLY__INT_: case _POLY__ZINT:
      return smod(*a._POLYptr,b);
    case _VECT__INT_: case _VECT__ZINT:
      if (a.subtype==_POLY1__VECT)
	return gen(trim(smod(*a._VECTptr,b),0),a.subtype);
      return gen(smod(*a._VECTptr,b),a.subtype);
    default: 
      if (a.type==_SYMB)
	return smodSYMB(a,b);
      if ( (b.type==_INT_) || (b.type==_ZINT) )
	return a;
      // error, b must be _DOUBLE_
      throw(std::runtime_error("smod 2nd argument must be _DOUBLE_"));
      return a;
    }
  }

  string print_DOUBLE_(double d,GIAC_CONTEXT){
    if (my_isnan(d))
      return "undef";
    if (my_isinf(d))
      return "infinity";
    string & forme=format_double(contextptr);
    if (!forme.empty()){
      char ch=forme[0];
      if (tolower(ch)!='g' && tolower(ch)!='a' && tolower(ch)!='f' && tolower(ch)!='e' )
	ch='g';
      if (forme.size()<2)
	setsizeerr();
      if (my_isnan(d))
	return "undef";
      if (my_isinf(d))
	return "infinity";
      char s[256];
      sprintf(s,("%."+forme.substr(1,forme.size()-1)+ch).c_str(),d);
      return s;
    } 
    if (xcas_mode(contextptr)==3 && double(int(d))==d)
      return print_INT_(int(d));
    char s[256];
#ifdef SOFTMATH
    sprintf(s,"%.14g",d);
    return s;
#else
    string form("%."+print_INT_(min(decimal_digits(contextptr),14)));
    int sf=scientific_format(contextptr); 
    switch (sf){
    case 0: case 2:
      form += "g";
      break;
    case 1: 
      form += "e";
      break;
    case 3:
      form += "a";
    }
    if (sf==2){
      // engineering format
      int ndigits=int(giac_floor(std::log10(d)+0.5));
      ndigits = 3*(ndigits /3);
      sprintf(s,form.c_str(),d/std::pow(10.0,ndigits));
      return s+("e"+print_INT_(ndigits));
    }
    sprintf(s,form.c_str(),d);
    if (sf || d!=giac_floor(d) || d>=1e12 || d<=-1e12)
      return s;
    else
      return string(s)+".0";
#endif
  }

  gen maptoarray(const gen_map & m,GIAC_CONTEXT){
    vecteur res;
    gen_map::const_iterator it=m.begin(),itend=m.end();
    if (it==itend)
      setdimerr("Empty array");
    gen_map::const_reverse_iterator lastit=m.rbegin();
    // find index ranges
    gen premidx=it->first,lastidx=lastit->first;
    if (premidx.type!=_VECT || lastidx.type!=_VECT)
      settypeerr("Bad array indexes");
    vecteur & pv=*premidx._VECTptr;
    vecteur & lv=*lastidx._VECTptr;
    unsigned ps=pv.size();
    vector<int> indexes(ps);
    if (lv.size()!=ps)
      setdimerr(contextptr);
    for (unsigned i=0;i<ps;++i){
      res.push_back(symb_interval(pv[i]+(xcas_mode(contextptr)!=0),lv[i]+(xcas_mode(contextptr)!=0)));
      if (lv[i].type!=_INT_ || pv[i].type!=_INT_ || lv[i].val<pv[i].val)
	setdimerr(contextptr);
      indexes[i]=(lv[i]-pv[i]).val+1;
    }
    if (ps==1){
      vecteur tmp;
      for (;it!=itend;++it){
	tmp.push_back(it->second);
      }
      res.push_back(tmp);
    }
    else {
      vecteur tmp(indexes[ps-1]);
      for (int i=ps-2;i>=0;--i){
	vecteur newtmp;
	for (int j=0;j<indexes[i];++j)
	  newtmp.push_back(tmp);
	tmp=newtmp;
      }
      gen tab=tmp;
      for (;it!=itend;++it){
	gen * tmpptr=&tab;
	if (tmpptr->type!=_VECT)
	  setdimerr(contextptr);
	for (unsigned i=0;i<ps;++i){
	  int pos=(it->first[i]-pv[i]).val;
	  if (pos<0 || unsigned(pos)>=tmpptr->_VECTptr->size())
	    setdimerr(contextptr);
	  tmpptr=&((*tmpptr->_VECTptr)[pos]);
	}
	*tmpptr = it->second;
      }
      res.push_back(tmp);
    }
    return new symbolic(at_array,res);
  }

  string printmap(const gen_map & m){
    string s("table(\n");
    gen_map::const_iterator it=m.begin(),itend=m.end();
    for (;it!=itend;){
      s += it->first.print(context0)+ " = " + it->second.print(context0);
      ++it;
      if (it!=itend)
	s += ',';
      s += '\n';
    }
    return s+")";
  }

  std::string printmpf_t(const mpf_t & inf,GIAC_CONTEXT){
    bool negatif=mpf_sgn(inf)<0;
    mp_exp_t expo;
#ifdef VISUALC
    char * ptr=new char[decimal_digits(contextptr)+2];
#else
    char ptr[decimal_digits(contextptr)+2];
#endif
    if (negatif){
      mpf_t inf2;
      mpf_init(inf2);
      mpf_neg(inf2,inf);
      mpf_get_str(ptr,&expo,10,decimal_digits(contextptr),inf2);
      mpf_clear(inf2);
    }
    else
      mpf_get_str(ptr,&expo,10,decimal_digits(contextptr),inf);
    std::string res(ptr),reste(res.substr(1,res.size()-1));
#ifdef VISUALC
    delete [] ptr;
#endif
    res=res[0]+("."+reste);
    if (expo!=1)
      res += "e"+print_INT_(expo-1);
    if (negatif)
      return "-"+res;
    else
      return res;
  }


  string print_binary(const real_object & r){
#ifdef HAVE_LIBMPFR
    mp_exp_t expo;
    int dd=mpfr_get_prec(r.inf);
#ifdef VISUALC
    char * ptr=new char[dd+2];
#else
    char ptr[dd+2];
#endif
    if (!mpfr_get_str(ptr,&expo,2,dd,r.inf,GMP_RNDN) || !(*ptr))
      setsizeerr("MPFR print binary error "+r.print(context0));
    string res;
    if (ptr[0]=='-')
      res="-0000."+string(ptr+1);
    else
      res="0000."+string(ptr);
#ifdef VISUALC
    delete [] ptr
#endif // VISUALC
    return res+"E"+print_INT_(expo);
#else // MPFR
    setsizeerr("Error no MPFR printing "+r.print(context0));
    return "undef";
#endif
  }

  gen read_binary(const string & s,unsigned int precision){
#ifdef HAVE_LIBMPFR
    real_object r;
    mpfr_set_prec(r.inf,precision);
#ifndef HAVE_MPFR_SET_STR_RAW
    // MPFR 2.2
    mpfr_strtofr (r.inf, (char *)s.c_str(), 0, 2, GMP_RNDN); 
#else
    // FOR MPFR 2.0 use instead 
    mpfr_set_str_raw(r.inf,(char *)s.c_str());
#endif // GNUWINCE
    return r;
    setsizeerr("MPFR error reading binary "+s);
#else // HAVE_LIBMPFR
    setsizeerr("Error no MPFR reading "+s);
    return undef;
#endif // HAVE_LIBMPFR
  }

  std::string real_object::print(GIAC_CONTEXT) const{ 
#ifdef HAVE_LIBMPFR
    bool negatif=mpfr_sgn(inf)<0;
    mp_exp_t expo;
    int dd=mpfr_get_prec(inf);
    dd=bits2digits(dd);
#ifdef VISUALC
    char * ptr=new char[dd+2];
#else
    char ptr[dd+2];
#endif
    if (negatif){
      mpfr_t inf2;
      mpfr_init2(inf2,mpfr_get_prec(inf));
      mpfr_neg(inf2,inf,GMP_RNDN);
      mpfr_get_str(ptr,&expo,10,dd,inf2,GMP_RNDN);
      mpfr_clear(inf2);
    }
    else
      mpfr_get_str(ptr,&expo,10,dd,inf,GMP_RNDN);
    std::string res(ptr);
    if (expo){
      if (expo==1){
	string reste(res.substr(1,res.size()-1));
	res=res[0]+("."+reste);
      }
      else
	res = "0."+res+"e"+print_INT_(expo);
    }
    else
      res="0."+res;
#ifdef VISUALC
    delete [] ptr;
#endif
    if (negatif)
      return "-"+res;
    else
      return res;
#else
    return printmpf_t(inf,contextptr);
#endif
  }

  string print_ZINT(const mpz_t & a){
    /*
    char * s =mpz_get_str (NULL, 10,a) ;
    string res(s);
    free(s);
    return res;
    */
    size_t l=mpz_sizeinbase (a, 10) + 2;
    if (l>unsigned(MAX_PRINTABLE_ZINT))
      return "Integer_too_large_for_display";
#ifdef VISUALC
    char * s = new char[l];
#else
    char s[l];
#endif
    mpz_get_str (s, 10,a) ;
    string tmp(s);
#ifdef VISUALC
    delete [] s;
#endif
    return tmp;
  }

  string hexa_print_ZINT(const mpz_t & a){
    size_t l=mpz_sizeinbase (a, 16) + 2;
    if (l>unsigned(MAX_PRINTABLE_ZINT))
      return "Integer_too_large";
#ifdef VISUALC
    char * s = new char[l];
#else
    char s[l];
#endif
    string res("0x");
#ifdef USE_GMP_REPLACEMENTS
    if (mpz_sgn(a) == -1){
      mpz_t tmpint;
      mpz_init(tmpint);
      mpz_neg(tmpint, a);
      mpz_get_str(s,16,tmpint);
      mpz_clear(tmpint);      
      res = "-" + res + s;
    }
    else {
      mpz_get_str(s,16,a);
      res += s;
    }
#else
    mpz_get_str (s,16,a) ;
    res += s;
#endif // USE_GMP_REPLACEMENTS
#ifdef VISUALC
    delete [] s;
#endif
    return res;
  }

  string octal_print_ZINT(const mpz_t & a){
    size_t l=mpz_sizeinbase (a, 8) + 2;
    if (l>unsigned(MAX_PRINTABLE_ZINT))
      return "Integer_too_large";
#ifdef VISUALC
    char * s = new char[l];
#else
    char s[l];
#endif
    string res("0");
#ifdef USE_GMP_REPLACEMENTS
    if (mpz_sgn(a) == -1){
      mpz_t tmpint;
      mpz_init(tmpint);
      mpz_neg(tmpint, a);
      mpz_get_str(s,8,tmpint);
      mpz_clear(tmpint);
      res = "-" + res + s;
    }
    else {
      mpz_get_str(s,8,a);
      res += s;
    }
#else
    mpz_get_str (s,8,a) ;
    res += s;
#endif // USE_GMP_REPLACEMENTS
#ifdef VISUALC
    delete [] s;
#endif
    return res;
  }

  string binary_print_ZINT(const mpz_t & a){
    size_t l=mpz_sizeinbase (a, 2) + 2;
    if (l>unsigned(MAX_PRINTABLE_ZINT))
      return "Integer_too_large";
#ifdef VISUALC
    char * s = new char[l];
#else
    char s[l];
#endif
    string res("0b");
#ifdef USE_GMP_REPLACEMENTS
    if (mpz_sgn(a) == -1){
      mpz_t tmpint;
      mpz_init(tmpint);
      mpz_neg(tmpint, a);
      mpz_get_str(s,2,tmpint);
      mpz_clear(tmpint);
      res = "-" + res + s;
    }
    else {
      mpz_get_str(s,2,a);
      res += s;
    }
#else
    mpz_get_str (s,2,a) ;
    res += s;
#endif // USE_GMP_REPLACEMENTS
#ifdef VISUALC
    delete [] s;
#endif
    return res;
  }

  string printinner_VECT(const vecteur & v, int subtype,GIAC_CONTEXT){
    string s;
    vecteur::const_iterator it=v.begin(), itend=v.end();
    if (it==itend)
      return "";
    for(;;){
      if ( (subtype==_RPN_FUNC__VECT) && (it->type==_SYMB) && (it->_SYMBptr->sommet==at_quote))
	s += "'"+it->_SYMBptr->feuille.print(contextptr)+"'";
      else {
	if ( (it->type==_SYMB) && (it->_SYMBptr->sommet==at_sto) )
	  s += "("+it->print(contextptr)+")";
	else
	  s += it->print(contextptr);
      }
      ++it;
      if (it==itend){
	return s;
      }
      if ( (subtype!=_RPN_FUNC__VECT) && 
	   //	   (subtype || (!rpn_mode) ) &&
	   ( ((it-1)->type!=_SYMB) || ((it-1)->_SYMBptr->sommet!=at_comment) )
	   )
	s += ',';
      else
	s += ' ';
    }
  }

  string begin_VECT_string(int subtype,bool tex,GIAC_CONTEXT){
    string s;
    switch (subtype){
    case _SEQ__VECT:
      break;
    case _SET__VECT:
      if (xcas_mode(contextptr)>0){
	if (tex)
	  s+="\\{";
	else
	  s="{";
      }
      else
	s="set[";
      break;
    case _RPN_STACK__VECT:
      s="stack(";
      break;
    case _RPN_FUNC__VECT:
      s="<< ";
      break;
    case _GROUP__VECT:
      s="group[";
      break;
    case _LINE__VECT:
      s="line[";
      break;
    case _VECTOR__VECT:
      s="vector[";
      break;
    case _PNT__VECT:
      s="pnt[";
      break;
    case _POINT__VECT:
      s="point[";
      break;
    case _MATRIX__VECT:
      s="matrix[";
      break;
    case _POLY1__VECT:
      s="poly1[";
      break;
    case _ASSUME__VECT:
      s = "assume[";
      break;
    case _FOLDER__VECT:
      s = "folder[";
      break;
    case _POLYEDRE__VECT:
      s= "polyedre[";
      break;
    case _RGBA__VECT:
      s= "rgba[";
      break;
    default:
      s="[";
    }
    return s;
  }

  string end_VECT_string(int subtype,bool tex,GIAC_CONTEXT){
    string s;
    switch (subtype){
    case _SEQ__VECT:
      return s;
    case _SET__VECT:
      if (xcas_mode(contextptr)>0){
	if (tex)
	return "\\}";
	else
	  return "}";
      }
      else
	return "]";
    case _RPN_STACK__VECT:
      return ")";
    case _RPN_FUNC__VECT:
      return " >>";
    default:
      return "]";
    }    
  }

  string print_VECT(const vecteur & v,int subtype,GIAC_CONTEXT){
    if (v.empty()){
      switch (subtype){
      case _SEQ__VECT:
	return "NULL";
      case _SET__VECT:
	if (xcas_mode(contextptr)>0)
	  return "{ }";
	else
	  return "%{ %}";
      case _RPN_FUNC__VECT:
        return "<< >>";
      case _RPN_STACK__VECT:
        return "stack()";
      }
    }
    string s;
    if (subtype==_SPREAD__VECT){
      s = "spreadsheet[";
      int nr=v.size();
      int nc=v.front()._VECTptr->size();
      for (int i=0;;){
	int save_r,save_c;
	vecteur & w=*v[i]._VECTptr;
	s +='[';
	for (int j=0;;){
	  save_r=printcell_current_row(contextptr);
	  save_c=printcell_current_col(contextptr);
	  printcell_current_row(contextptr)=i;
	  printcell_current_col(contextptr)=j;
	  s += w[j].print(contextptr);
	  printcell_current_row(contextptr)=save_r;
	  printcell_current_col(contextptr)=save_c;
	  ++j;
	  if (j==nc)
	    break;
	  else
	    s+=',';
	}
	s+=']';
	++i;
	if (i==nr)
	  break;
	else
	  s+=',';
      }
      return s+']';
    }
    if ( ( subtype==_SEQ__VECT) && (v.size()==1) && !xcas_mode(contextptr))
      return "seq["+v.front().print(contextptr)+"]";
    else
      s=begin_VECT_string(subtype,false,contextptr);
    s += printinner_VECT(v,subtype,contextptr);
    return s+end_VECT_string(subtype,false,contextptr);
  }

  string print_SPOL1(const sparse_poly1 & p,GIAC_CONTEXT){
    if (p.empty())
      return "0";
    string s;
    vector<monome>::const_iterator it=p.begin(), itend=p.end();
    for(;;){
      s += it->print();
      ++it;
      if (it==itend)
          return s;
      s += '+';
    }
  }

  string print_the_type(int val,GIAC_CONTEXT){
    if (xcas_mode(contextptr)==1){
      switch(val){
      case _INT_:
	return "integer";
      case _DOUBLE_:
	return "float";
      case _ZINT:
	return "integer";
      case _CPLX:
	return "complex";
      case _VECT:
	return "vector";
      case _IDNT:
	return "symbol";
      case _SYMB:
	return "algebraic";
      case _FRAC:
	return "rational";
      case _MAPLE_LIST:
	return "list";
      }
    }
    switch(val){
    case _INT_:
      return "DOM_int";
    case _DOUBLE_:
      return "DOM_FLOAT";
    case _ZINT:
      return "DOM_INT";
    case _CPLX:
      return "DOM_COMPLEX";
    case _VECT:
      return "DOM_LIST";
    case _IDNT:
      return "DOM_IDENT";
    case _SYMB:
      return "DOM_SYMBOLIC";
    case _FRAC:
      return "DOM_RAT";
    case _STRNG:
      return "DOM_STRING";
    case _FUNC:
      return "DOM_FUNC";
    case _REAL:
      return "DOM_LONGFLOAT";
    }
    return print_INT_(val);
  }

  string print_STRNG(const string & s){
    string res("\"");
    int l=s.size();
    for (int i=0;i<l;++i){
      res += s[i];
      if (s[i]=='"')
	res += '"';
    }
    return res+'"';
  }

  string printi(GIAC_CONTEXT){
    if (xcas_mode(contextptr)==3)
      return "\xa1";
    if (xcas_mode(contextptr)>0)
      return "I";
    else
      return "i";
  }

  string print_EQW(const eqwdata & e){
    string s;
    s += "eqwdata(position"+print_INT_(e.x)+","+print_INT_(e.y)+",dxdy,"+print_INT_(e.dx)+","+print_INT_(e.dy);
    s += ",font,"+print_INT_(e.eqw_attributs.fontsize);
    s += ",background," +print_INT_(e.eqw_attributs.background);
    s += ",text_color," +print_INT_(e.eqw_attributs.text_color)+",";
    if (e.selected)
      s +="selected,";
    if (e.active)
      s +="active,";
    s += e.g.print(context0)+")";
    return s;
  }

  string gen::print() const{
    return print(context0);
  }

  string gen::print_universal(GIAC_CONTEXT) const{
    int lang=language(contextptr);
    language(-1,contextptr);
    string res;
    try {
      res=print(contextptr);
    }
    catch (...){
      
    }
    language(lang,contextptr);
    return res;
  }

  string gen::print(GIAC_CONTEXT) const{
    switch (type ) {
    case _INT_: 
      if (subtype==_INT_SOLVER){
	switch (val){
	case _UNFACTORED:
	  return "unfactored";
	case _BISECTION_SOLVER:
	  return "bisection_solver";
	case _FALSEPOS_SOLVER:
	  return "falsepos_solver";  
	case _BRENT_SOLVER:
	  return "brent_solver";
	case _NEWTON_SOLVER:
	  return "newton_solver";
	case _DNEWTON_SOLVER:
	  return "dnewton_solver";
	case _NEWTONJ_SOLVER:
	  return "newtonj_solver";
	case _SECANT_SOLVER:
	  return "secant_solver";
	case _STEFFENSON_SOLVER:
	  return "steffenson_solver";
	case _HYBRIDSJ_SOLVER:
	  return "hybridsj_solver";
	case _HYBRIDJ_SOLVER:
	  return "hybridj_solver";
	case _HYBRIDS_SOLVER:
	  return "hybrids_solver";
	case _HYBRID_SOLVER:
	  return "hybrid_solver";
	case _GOLUB_REINSCH_DECOMP:
	  return "golub_reinsch_decomp";
	case _GOLUB_REINSCH_MOD_DECOMP:
	  return "golub_reinsch_mode_decomp";
	case _JACOBI_DECOMP:
	  return "jacobi_decomp";
	case _MINOR_DET:
	  return "minor_det";
	case _HESSENBERG_PCAR:
	  return "hessenberg_pcar";
	case _RATIONAL_DET:
	  return "rational_det";
	case _KEEP_PIVOT:
	  return "keep_pivot";
	case _FADEEV:
	  return "fadeev";
	case _BAREISS:
	  return "bareiss";
	default:
	  return print_INT_(val);
	}
      }
      if (subtype==_INT_BOOLEAN){
	if (xcas_mode(contextptr)==2){
	  if (val)
	    return "TRUE";
	  else
	    return "FALSE";
	}
	else {
	  if (val)
	    return "true";
	  else
	    return "false";
	}
      }
      if (subtype==_INT_COLOR){
	switch (language(contextptr)){
	case 1:
	  switch (val){
	  case _BLACK:
	    return "noir";
	  case _RED:
	    return "rouge";
	  case _GREEN:
	    return "vert";
	  case _YELLOW:
	    return "jaune";
	  case _BLUE:
	    return "bleu";
	  case _MAGENTA:
	    return "magenta";
	  case _CYAN:
	    return "cyan";
	  case _WHITE:
	    return "blanc";
	  }
	  break;
	default:
	  switch (val){
	  case _BLACK:
	    return "black";
	  case _RED:
	    return "red";
	  case _GREEN:
	    return "green";
	  case _YELLOW:
	    return "yellow";
	  case _BLUE:
	    return "blue";
	  case _MAGENTA:
	    return "magenta";
	  case _CYAN:
	    return "cyan";
	  case _WHITE:
	    return "white";
	  }
	}
	switch (val){
	  case _FILL_POLYGON:
	    return "rempli";
	  case _QUADRANT2:
	    return "quadrant2";
	  case _QUADRANT3:
	    return "quadrant3";
	  case _QUADRANT4:
	    return "quadrant4";
	  case _POINT_LOSANGE:
	    return "point_losange";
	  case _POINT_CARRE:
	    return "point_carre";
	  case _POINT_PLUS:
	    return "point_plus";
	  case _POINT_TRIANGLE:
	    return "point_triangle";
	  case _POINT_ETOILE:
	    return "point_etoile";
	  case _POINT_POINT:
	    return "point_point";
	  case _POINT_INVISIBLE:
	    return "point_invisible";
	  case 49:
	    return "gomme";
	  case _LINE_WIDTH_2:
	    return "line_width_2";
	  case _LINE_WIDTH_3:
	    return "line_width_3";
	  case _LINE_WIDTH_4:
	    return "line_width_4";
	  case _LINE_WIDTH_5:
	    return "line_width_5";
	  case _LINE_WIDTH_6:
	    return "line_width_6";
	  case _LINE_WIDTH_7:
	    return "line_width_7";
	  case _LINE_WIDTH_8:
	    return "line_width_8";
	  case _POINT_WIDTH_2:
	    return "point_width_2";
	  case _POINT_WIDTH_3:
	    return "point_width_3";
	  case _POINT_WIDTH_4:
	    return "point_width_4";
	  case _POINT_WIDTH_5:
	    return "point_width_5";
	  case _POINT_WIDTH_6:
	    return "point_width_6";
	  case _POINT_WIDTH_7:
	    return "point_width_7";
	  case _POINT_WIDTH_8:
	    return "point_width_8";
	case _HIDDEN_NAME:
	  return "hidden_name";
	  case _DASH_LINE:
	    return "dash_line";
	  case _DOT_LINE:
	    return "dot_line";
	  case _DASHDOT_LINE:
	    return "dashdot_line";
	  case _DASHDOTDOT_LINE:
	    return "dashdotdot_line";
	  case _CAP_FLAT_LINE:
	    return "cap_flat_line";
	  case _CAP_ROUND_LINE:
	    return "cap_round_line";
	  case _CAP_SQUARE_LINE:
	    return "cap_square_line";
	}
      }
      if (subtype==_INT_PLOT){
	switch(val){
	case _ADAPTIVE:
	  return "adaptive";
	case _AXES:
	  return "axes";
	case _COLOR:
	  return "color";
	case _FILLED:
	  return "filled";
	case _FONT:
	  return "font";
	case _LABELS:
	  return "labels";
	case _LEGEND:
	  return "legend";
	case _LINESTYLE:
	  return "linestyle";
	case _RESOLUTION:
	  return "resolution";
	case _SAMPLE:
	  return "sample";
	case _SCALING:
	  return "scaling";
	case _STYLE:
	  return "style";
	case _SYMBOL:
	  return "symbol";
	case _SYMBOLSIZE:
	  return "symbolsize";
	case _THICKNESS:
	  return "thickness";
	case _TITLE:
	  return "title";
	case _TITLEFONT:
	  return "titlefont";
	case _VIEW:
	  return "view";
	case _AXESFONT:
	  return "axesfont";
	case _COORDS:
	  return "coords";
	case _LABELFONT:
	  return "labelfont";
	case _LABELDIRECTIONS:
	  return "labeldirections";
	case _NUMPOINTS:
	  return "numpoints";
	case _TICKMARKS:
	  return "tickmarks";
	case _XTICKMARKS:
	  return "xtickmarks";
	case _NSTEP:
	  return "nstep";
	case _XSTEP:
	  return "xstep";
	case _YSTEP:
	  return "ystep";
	case _ZSTEP:
	  return "zstep";
	case _TSTEP:
	  return "tstep";
	case _USTEP:
	  return "ustep";
	case _VSTEP:
	  return "vstep";
	case _FRAMES:
	  return "frames";
	case _GL_TEXTURE:
	  return "gl_texture";
	case _GL_LIGHT0:
	  return "gl_light0";
	case _GL_LIGHT1:
	  return "gl_light1";
	case _GL_LIGHT2:
	  return "gl_light2";
	case _GL_LIGHT3:
	  return "gl_light3";
	case _GL_LIGHT4:
	  return "gl_light4";
	case _GL_LIGHT5:
	  return "gl_light5";
	case _GL_LIGHT6:
	  return "gl_light6";
	case _GL_LIGHT7:
	  return "gl_light7";
	case _GL_AMBIENT:
	  return "gl_ambient";
	case _GL_SPECULAR:
	  return "gl_specular";
	case _GL_DIFFUSE:
	  return "gl_diffuse";
	case _GL_POSITION:
	  return "gl_position";
	case _GL_SPOT_DIRECTION:
	  return "gl_spot_direction";
	case _GL_SPOT_EXPONENT:
	  return "gl_spot_exponent";
	case _GL_SPOT_CUTOFF:
	  return "gl_spot_cutoff";
	case _GL_CONSTANT_ATTENUATION:
	  return "gl_constant_attenuation";
	case _GL_LINEAR_ATTENUATION:
	  return "gl_linear_attenuation";
	case _GL_QUADRATIC_ATTENUATION:
	  return "gl_quadratic_attenuation";
	case _GL_OPTION:
	  return "gl_option";
	case _GL_SMOOTH:
	  return "gl_smooth";
	case _GL_FLAT:
	  return "gl_flat";
	case _GL_SHININESS:
	  return "gl_shininess";
	case _GL_FRONT:
	  return "gl_front";
	case _GL_BACK:
	  return "gl_back";
	case _GL_FRONT_AND_BACK:
	  return "gl_front_and_back";
	case _GL_AMBIENT_AND_DIFFUSE:
	  return "gl_ambient_and_diffuse";
	case _GL_EMISSION:
	  return "gl_emission";
	case _GL_LIGHT_MODEL_AMBIENT:
	  return "gl_light_model_ambient";
	case _GL_LIGHT_MODEL_LOCAL_VIEWER: 
	  return "gl_light_model_local_viewer";
	case _GL_LIGHT_MODEL_TWO_SIDE:
	  return "gl_light_model_two_side";
	case _GL_LIGHT_MODEL_COLOR_CONTROL:
	  return "gl_light_model_color_control";
	case _GL_BLEND:
	  return "gl_blend";
	case _GL_SRC_ALPHA:
	  return "gl_src_alpha";
	case _GL_ONE_MINUS_SRC_ALPHA:
	  return "gl_one_minus_src_alpha";
	case _GL_SEPARATE_SPECULAR_COLOR:
	  return "gl_separate_specular_color";
	case _GL_SINGLE_COLOR:
	  return "gl_single_color";
	case _GL_MATERIAL:
	  return "gl_material";
	case _GL_COLOR_INDEXES:
	  return "gl_color_indexes";
	case _GL_LIGHT:
	  return "gl_light";
	case _GL_PERSPECTIVE:
	  return "gl_perspective";
	case _GL_ORTHO:
	  return "gl_ortho";
	case _GL_QUATERNION:
	  return "gl_quaternion";
	case _GL_ROTATION_AXIS:
	  return "gl_rotation_axis";
	case _GL_X:
	  return "gl_x";
	case _GL_Y:
	  return "gl_y";
	case _GL_Z:
	  return "gl_z";
	case _GL_XTICK:
	  return "gl_xtick";
	case _GL_YTICK:
	  return "gl_ytick";
	case _GL_ZTICK:
	  return "gl_ztick";
	case _GL_ANIMATE:
	  return "gl_animate";
	case _GL_SHOWAXES:
	  return "gl_showaxes";
	case _GL_SHOWNAMES:
	  return "gl_shownames";
	case _GL_X_AXIS_NAME:
	  return "gl_x_axis_name";
	case _GL_Y_AXIS_NAME:
	  return "gl_y_axis_name";
	case _GL_Z_AXIS_NAME:
	  return "gl_z_axis_name";
	case _GL_X_AXIS_UNIT:
	  return "gl_x_axis_unit";
	case _GL_Y_AXIS_UNIT:
	  return "gl_y_axis_unit";
	case _GL_Z_AXIS_UNIT:
	  return "gl_z_axis_unit";
	}
      }
      if (subtype==_INT_MAPLELIB){
	switch (val){
	case _LINALG:
	  return "linalg";
	case _NUMTHEORY:
	  return "numtheory";
	case _GROEBNER:
	  return "groebner";
	}
      }
      if (subtype==_INT_MAPLECONVERSION){
	switch (val){
	case _MAPLE_LIST:
	  return "list";
	case _SET__VECT:
	  return "set";
	case _MATRIX__VECT:
	  return "matrix";
	case _POLY1__VECT:
	  return "polynom";
	case _TRIG:
	  return "trig";
	case _EXPLN:
	  return "expln";
	case _PARFRAC:
	  return "parfrac";
	case _FULLPARFRAC:
	  return "fullparfrac";
	case _CONFRAC:
	  return "confrac";
	case _BASE:
	  return "base";
	case _POSINT:
	  return "posint";
	case _NEGINT:
	  return "negint";
	case _NONPOSINT:
	  return "nonposint";
	case _NONNEGINT:
	  return "nonnegint";
	}
      }
      if (subtype==_INT_MUPADOPERATOR){
	switch (val){
	case _DELETE_OPERATOR:
	  return "Delete";
	case _PREFIX_OPERATOR:
	  return "Prefix";
	case _POSTFIX_OPERATOR:
	  return "Postfix";
	case _BINARY_OPERATOR:
	  return "Binary";
	case _NARY_OPERATOR:
	  return "Nary";
	}
      }
      if (subtype==_INT_GROEBNER){
	switch (val){
	case _REVLEX_ORDER:
	  return "revlex";
	case _PLEX_ORDER:
	  return "plex";
	case _TDEG_ORDER:
	  return "tdeg";
	case _WITH_COCOA:
	  return "with_cocoa";
	case _WITH_F5:
	  return "with_f5";
	}
      }
      if (subtype!=_INT_TYPE){
	int format=integer_format(contextptr);
	switch (format){
	case 16:
	  return hexa_print_INT_(val);
	case 8:
	  return octal_print_INT_(val);
	default:
	  return print_INT_(val);
	}
      }
      return print_the_type(val,contextptr);
    case _DOUBLE_:
      return print_DOUBLE_(_DOUBLE_val,contextptr);
    case _ZINT: 
      switch (integer_format(contextptr)){
      case 16:
	return hexa_print_ZINT(*_ZINTptr);
      case 8:
	return octal_print_ZINT(*_ZINTptr);
      default:
	return print_ZINT(*_ZINTptr);
      }
    case _REAL:
      return _REALptr->print(contextptr);
    case _CPLX:
      if (is_exactly_zero(*(_CPLXptr+1)))
	return _CPLXptr->print(contextptr);
      if (is_exactly_zero(*_CPLXptr)){
	if (is_one(*(_CPLXptr+1)))
	  return printi(contextptr);
	if (is_minus_one(*(_CPLXptr+1)))
	  return "-"+printi(contextptr);
	return (_CPLXptr+1)->print(contextptr) + "*"+printi(contextptr);
      }
      if (is_one(*(_CPLXptr+1)))
	return _CPLXptr->print(contextptr) + "+"+printi(contextptr);
      if (is_minus_one(*(_CPLXptr+1)))
	return _CPLXptr->print(contextptr) + "-"+printi(contextptr);
      if (is_positive(-(*(_CPLXptr+1)),contextptr))
	return _CPLXptr->print(contextptr) + string("-") + (-(*(_CPLXptr+1))).print(contextptr) + "*"+printi(contextptr);
      return _CPLXptr->print(contextptr) + string("+") + (_CPLXptr+1)->print(contextptr) + "*"+printi(contextptr);
    case _IDNT:
      return _IDNTptr->print(contextptr);
    case _SYMB:
      if (subtype==_SPREAD__SYMB){
	if (_SYMBptr->sommet==at_sto)
	  return "=("+_SYMBptr->print(contextptr)+")";
	return "="+_SYMBptr->print(contextptr);
      }
      else
	return _SYMBptr->print(contextptr);
    case _VECT:
      return print_VECT(*_VECTptr,subtype,contextptr);
    case _POLY:
      return _POLYptr->print() ;
    case _SPOL1:
      return print_SPOL1(*_SPOL1ptr,contextptr);
    case _EXT:
      return "%%{"+_EXTptr->print(contextptr)+':'+(*(_EXTptr+1)).print(contextptr)+"%%}";
    case _USER:
      return _USERptr->print(contextptr);
    case _MOD:
      if ( (_MODptr->type==_SYMB && _MODptr->_SYMBptr->sommet!=at_pow) || (_MODptr->type==_VECT && _MODptr->subtype==_SEQ__VECT) )
	return "("+_MODptr->print(contextptr)+") % "+(*(_MODptr+1)).print(contextptr);
      return _MODptr->print(contextptr)+" % "+(*(_MODptr+1)).print(contextptr);
    case _ROOT:
      return _ROOTptr->print();
    case _FRAC:
      return _FRAC2_SYMB(*this).print(contextptr);
    case _STRNG:
      return print_STRNG(*_STRNGptr);
    case _FUNC:
      if (*this==at_return){
	if (xcas_mode(contextptr)==3)
	  return "Return";
	else
	  return "return ;";
      }
      if (rpn_mode || _FUNCptr->ptr->printsommet==&printastifunction) 
	return _FUNCptr->ptr->print(contextptr);
      else
	return "'"+_FUNCptr->ptr->print(contextptr)+"'";
    case _MAP:
      if (subtype==1)
	return maptoarray(*_MAPptr,contextptr).print(contextptr);
      else
	return printmap(*_MAPptr);
    case _EQW:
      return print_EQW(*_EQWptr);
    case _POINTER_:
      return "pointer("+hexa_print_INT_((unsigned long)_POINTER_val)+","+print_INT_(subtype)+")";
    default:
      settypeerr("print");
    }
    return "error";
  }

  void gen::dbgprint() const{    
    if (this->type==_POLY)
      _POLYptr->dbgprint();
    else
      cout << this->print(context0) << endl; 
  }

  ostream & operator << (ostream & os,const gen & a) { return os << a.print(context0); }

  string monome::print() const {
    return '<' + coeff.print(context0) + "," + exponent.print(context0) + '>' ;
  }

  void monome::dbgprint() const {
    cout << this->print();
  }

  ostream & operator << (ostream & os,const monome & m){
    return os << m.print() ;
  }

  /*
  gen string2_ZINT(string s,int l,int & pos){
    char ss[l+1];
    int neg=1;
    if (s[pos]=='-'){
      pos++;
      neg=-1;
    }
    int i=0;
    for (;(pos<l) && (s[pos]>='0') && (s[pos]<='9');pos++,i++)
      ss[i]=s[pos];
    if ((!i) && (s[pos]=='i') || (s[pos]=='I')){
      return(neg);
    }
    assert(i);
    ss[i]=char(0);
    mpz_t *mpzin = new mpz_t[1];
    mpz_init(*mpzin);
    mpz_set_str (*mpzin, ss, 10);
    if (neg>0)
      return(gen(mpzin));
    else
      return(-gen(mpzin));
  }

  istream & operator >> (istream & is,gen & a){
    string s;
    is >> s;
    int l=s.size();
    int pos=0;
    a=gen(0);
    while (pos<l){
      if ((s[pos]=='i') || (s[pos]=='I')){
	a=a+gen(0,1);
	pos++;
      }
      else {
	if (s[pos]=='+')
	  pos++;
	else {
	  gen tmp(string2_ZINT(s,l,pos));
	  if (s[pos]=='*'){
	    pos++;
	    assert( (s[pos]=='i') || (s[pos]=='I') );
	    pos++; // skip *I
	    a=a+tmp*cst_i;
	  }
	  else {
	    if ((s[pos]=='i') || (s[pos]=='I')){
	      pos++; // skip I
	      a=a+tmp*cst_i;	
	    }
	    else
	      a=a+tmp;
	  }
	}
      }
    }
    return is;
  }
  */

  istream & operator >> (istream & is,gen & a){
    string s;
    is >> s;
    a = gen(s,list_one_letter__IDNT,context0);
    return is;
  }

  gen operator !(const gen & a){
    switch (a.type){
    case _INT_: case _ZINT: case _CPLX: case _DOUBLE_:
      return (is_zero(a));
    default:
      return symb_not(a);
    }
  }

  // equality of vecteurs representing geometrical lines
  bool geo_equal(const vecteur &v,const vecteur & w,int subtype,GIAC_CONTEXT){
    int vs=v.size(),ws=w.size();
    if (vs!=ws)
      return false;
    if (v==w)
      return true;
    if ( (subtype==_LINE__VECT)  && (vs==2)){
      if (v[1]==v[0])
	return v==w;
      // v[1]!=v[0]
      if (!is_zero(im(rdiv(w[0]-v[0],v[1]-v[0]),contextptr))) 
	return false;
      if (!is_zero(im(rdiv(w[1]-v[0],v[1]-v[0]),contextptr))) 
	return false;
      return true;
    }
    if (subtype==_SET__VECT){
      vecteur w1(w),v1(v);
      sort(w1.begin(),w1.end(),islesscomplexthanf);
      sort(v1.begin(),v1.end(),islesscomplexthanf);
      return w1==v1;
    }
    return false;
  }

  bool operator_equal(const gen & a,const gen & b,GIAC_CONTEXT){
    switch ( (a.type<< _DECALAGE) | b.type ) {
    case _INT___INT_: 
      return (a.val==b.val);
    case _INT___MOD: case _ZINT__MOD:
      return a==*b._MODptr;
    case _MOD__INT_: case _MOD__ZINT:
      return b==*a._MODptr;
    case _INT___ZINT: 
      return (mpz_cmp_si(*b._ZINTptr,a.val)==0);
    case _INT___DOUBLE_:
      return double(a.val)==b._DOUBLE_val;
    case _DOUBLE___INT_:
      return a._DOUBLE_val==double(b.val);
    case _ZINT__INT_:
      return (mpz_cmp_si(*a._ZINTptr,b.val)==0);
    case _ZINT__ZINT:
      return (mpz_cmp(*a._ZINTptr,*b._ZINTptr)==0);
    case _INT___CPLX: case _ZINT__CPLX:
      return ( operator_equal(a,re(b,contextptr),contextptr) && is_zero(im(b,contextptr)));
    case _CPLX__ZINT: case _CPLX__INT_:
      return ( operator_equal(re(a,contextptr),b,contextptr) && is_zero(im(a,contextptr))) ;
    case _CPLX__CPLX:
      return( operator_equal(*a._CPLXptr,*b._CPLXptr,contextptr) && operator_equal(*(a._CPLXptr+1),*(b._CPLXptr+1),contextptr) );
    case _DOUBLE___DOUBLE_:
      if (a._DOUBLE_val==b._DOUBLE_val)
	return true;
      if (my_isnan(a._DOUBLE_val) && my_isnan(b._DOUBLE_val))
	return true; // avoid infinite loop in evalf
      return  std::abs(a._DOUBLE_val-b._DOUBLE_val)<epsilon(contextptr);
    case _IDNT__IDNT:
      // return a.subtype==b.subtype && (a._IDNTptr->name==b._IDNTptr->name || *a._IDNTptr->name==*b._IDNTptr->name);
      return a._IDNTptr->name==b._IDNTptr->name || *a._IDNTptr->name==*b._IDNTptr->name;
    case _SYMB__SYMB:
      if (a._SYMBptr==b._SYMBptr)
	return true;
      return ( (a._SYMBptr->sommet==b._SYMBptr->sommet) && (a._SYMBptr->feuille==b._SYMBptr->feuille));
    case _VECT__VECT:
      if (a._VECTptr==b._VECTptr)
	return true;
      if (a.subtype!=b.subtype)
	return false;
      if (a.subtype)
	return geo_equal(*a._VECTptr,*b._VECTptr,a.subtype,contextptr);
      return *a._VECTptr==*b._VECTptr;
    case _POLY__POLY:
      if (a._POLYptr==b._POLYptr)
	return true;
      return (a._POLYptr->dim==b._POLYptr->dim) && (a._POLYptr->coord==b._POLYptr->coord);
    case _FRAC__FRAC:
      return (a._FRACptr->num==b._FRACptr->num) && (a._FRACptr->den==b._FRACptr->den);
    case _STRNG__STRNG:
      return (a._STRNGptr==b._STRNGptr) || (*a._STRNGptr==*b._STRNGptr);
    case _FUNC__FUNC:
      return (a._FUNCptr==b._FUNCptr) || (*a._FUNCptr==*b._FUNCptr);
    case _MOD__MOD: case _EXT__EXT:
      return ( (*a._EXTptr==*b._EXTptr) && (*(a._EXTptr+1)==*(b._EXTptr+1)) );
    case _SPOL1__SPOL1:
      return *a._SPOL1ptr==*b._SPOL1ptr;
    default: // Check pointers, type subtype
      if ((a.type==b.type) && (a.subtype==b.subtype) && (a.val==b.val))
	return true;
      if (a.type<=_REAL && b.type<=_REAL)
	return is_zero(a-b);
      if (a.type==_USER)
	return *a._USERptr==b;
      if (b.type==_USER)
	return *b._USERptr==a;
      return false;
    }
  }

  bool operator ==(const gen & a,const gen & b){
    return operator_equal(a,b,context0);
  }

  bool operator !=(const gen & a,const gen & b){
    return !(a==b);
  }

  gen equal(const gen & a,const gen &b){
    return new symbolic(at_equal,makenewvecteur(a,b));
  }

  gen sign(const gen & a,GIAC_CONTEXT){
    if (is_equal(a))
      return apply_to_equal(a,sign,contextptr);
    if (is_zero(a)){
      if (a.type==_DOUBLE_)
	return 0.0;
      else
	return 0;
    }
    if (a==plus_inf)
      return 1;
    if (a==minus_inf)
      return -1;
    if (is_undef(a))
      return a;
    if (is_inf(a))
      return undef;
    switch (a.type){
    case _INT_: case _ZINT: 
      if (is_positive(a,contextptr))
	return 1;
      else
	return -1;
    case _DOUBLE_:
      if (a._DOUBLE_val>0)
	return 1.0;
      else
	return -1.0;
    }
    gen b=evalf_double(a,1,contextptr);
    if (b.type==_DOUBLE_){
      double eps=epsilon(contextptr);
      if (b._DOUBLE_val>eps)
	return plus_one;
      if (b._DOUBLE_val<-eps)
	return minus_one;
      return zero;
    }
    if (is_zero(im(a,contextptr))){
      int s=sturmsign(a,contextptr); 
      if (s)
	return s;
    }
    return new symbolic(at_sign,a);
  }
  
  gen sym_is_greater(const gen & a,const gen & b,GIAC_CONTEXT){
    if (is_undef(a) || (a==unsigned_inf) || is_undef(b) || (b==unsigned_inf) || a.type==_VECT || b.type==_VECT)
      return undef;
    if (a==b)
      return false;
    if ( (b==plus_inf) || (a==minus_inf) )
      return false;
    if ( (b==minus_inf) || (a==plus_inf) )
      return true;
    if (a.is_symb_of_sommet(at_equal) && b.is_symb_of_sommet(at_equal) ){
      gen & af=a._SYMBptr->feuille;
      gen & bf=b._SYMBptr->feuille;
      if (af.type==_VECT && bf.type==_VECT && af._VECTptr->size()==2 && bf._VECTptr->size()==2 && af._VECTptr->front()==bf._VECTptr->front())
	return sym_is_greater(af._VECTptr->back(),bf._VECTptr->back(),contextptr);
    }
    if (a.type==_USER)
      return (*a._USERptr>b);
    if (b.type==_USER)
      return (*b._USERptr<=a);
    if (a.is_symb_of_sommet(at_superieur_strict) || a.is_symb_of_sommet(at_superieur_egal) || a.is_symb_of_sommet(at_inferieur_strict) || a.is_symb_of_sommet(at_inferieur_egal) )
      return false;
    if (b.is_symb_of_sommet(at_superieur_strict) || b.is_symb_of_sommet(at_superieur_egal) || b.is_symb_of_sommet(at_inferieur_strict) || b.is_symb_of_sommet(at_inferieur_egal) )
      return false;
    gen approx;
    if (has_evalf(a-b,approx,1,contextptr) && approx.type==_DOUBLE_ )
      return (approx._DOUBLE_val>0);
    gen g=sign(a-b,contextptr); 
    if (is_one(g))
      return plus_one;
    if (is_minus_one(g))
      return false;
    return symb_superieur_strict(a,b);
  }

  gen superieur_strict(const gen & a,const gen & b,GIAC_CONTEXT){
    switch ( (a.type<< _DECALAGE) | b.type ) {
    case _INT___INT_: 
      return (a.val>b.val);
    case _INT___ZINT: 
      return (mpz_cmp_si(*b._ZINTptr,a.val)<0);
    case _ZINT__INT_:
      return (mpz_cmp_si(*a._ZINTptr,b.val)>0);
    case _ZINT__ZINT:
      return (mpz_cmp(*a._ZINTptr,*b._ZINTptr)>0);
    case _DOUBLE___DOUBLE_:
      return a._DOUBLE_val>b._DOUBLE_val;
    case _DOUBLE___INT_:
      return a._DOUBLE_val>b.val;
    case _INT___DOUBLE_:
      return a.val>b._DOUBLE_val;
    case _DOUBLE___ZINT:
      return a._DOUBLE_val>mpz_get_d(*b._ZINTptr);
    case _ZINT__DOUBLE_:
      return mpz_get_d(*a._ZINTptr)>b._DOUBLE_val;
    default:
      if (a.type<=_REAL && b.type<=_REAL)
	return is_strictly_positive(a-b,contextptr);
      return sym_is_greater(a,b,contextptr);
    }
  }

  gen inferieur_strict(const gen & a,const gen & b,GIAC_CONTEXT){
    return superieur_strict(b,a,contextptr);
  }

  gen superieur_egal(const gen & a,const gen & b,GIAC_CONTEXT){
    gen g=!superieur_strict(b,a,contextptr);
    if (g.type==_INT_)
      return g;
    return symb_superieur_egal(a,b);
  }

  gen inferieur_egal(const gen & a,const gen & b,GIAC_CONTEXT){
    return superieur_egal(b,a,contextptr);
  }

  bool is_zero__VECT(const vecteur & v){
    vecteur::const_iterator it=v.begin(),itend=v.end();
    for (;it!=itend;++it){
      if (!is_zero(*it))
	return false;
    }
    return true;
  }

  bool real_object::is_zero(){
#ifdef HAVE_LIBMPFR
    return !mpfr_sgn(inf);
#else
    return !mpf_sgn(inf);
#endif
  }

  bool real_object::is_inf(){
#ifdef HAVE_LIBMPFR
    return !mpfr_inf_p(inf);
#else
    setsizeerr();
    return false;
#endif
  }

  bool real_object::is_nan(){
#ifdef HAVE_LIBMPFR
    return !mpfr_nan_p(inf);
#else
    setsizeerr();
    return false;
#endif
  }

  bool is_zero(const gen & a){
    switch (a.type ) {
    case _INT_: 
      return !a.val; 
    case _ZINT: 
      return (!mpz_sgn(*a._ZINTptr));
    case _REAL:
      return a._REALptr->is_zero();
    case _CPLX:
      return (is_zero(*a._CPLXptr) && is_zero(*(a._CPLXptr+1)));
    case _DOUBLE_:
      return (fabs(a._DOUBLE_val)<=epsilon(context0)); 
    case _VECT:
      return is_zero__VECT(*a._VECTptr);
    case _POLY:
      return a._POLYptr->coord.empty();
    case _FRAC:
      return is_zero(a._FRACptr->num);
    case _MOD:
      return is_zero(*a._MODptr);
    case _USER:
      return a._USERptr->is_zero();
    default: 
      return false;
    }
  }

  bool is_exactly_zero(const gen & a){
    switch (a.type ) {
    case _INT_: 
      return !a.val; 
    case _ZINT: 
      return (!mpz_sgn(*a._ZINTptr));
    case _REAL:
      return a._REALptr->is_zero();
    case _CPLX:
      return (is_exactly_zero(*a._CPLXptr) && is_exactly_zero(*(a._CPLXptr+1)));
    case _DOUBLE_:
      return a._DOUBLE_val==0; 
    case _POLY:
      return a._POLYptr->coord.empty();
    case _FRAC:
      return is_exactly_zero(a._FRACptr->num);
    case _MOD:
      return is_exactly_zero(*a._MODptr);
    case _USER:
      return a._USERptr->is_zero();
    default: 
      return false;
    }
  }

  bool is_one(const gen & a){
    switch (a.type ) {
    case _INT_: 
      return a.val==1; 
    case _ZINT: 
      return (a==gen(1));
    case _CPLX:
      return (is_one(*a._CPLXptr) && is_zero(*(a._CPLXptr+1)));
    case _DOUBLE_:
      return a._DOUBLE_val==1;
    case _REAL:
      return evalf_double(a,0,0)._DOUBLE_val==1;
    case _VECT:
      return a._VECTptr->size()==1 && is_one(a._VECTptr->front());
    case _POLY:
      return Tis_constant(*a._POLYptr) && (is_one(a._POLYptr->coord.front().value));
    case _FRAC:
      return a._FRACptr->num == a._FRACptr->den;
    case _MOD:
      return is_one(*a._MODptr);
    case _USER:
      return a._USERptr->is_one();
    default: 
      return false;
    }
  }

  bool is_minus_one(const gen & a){
    switch (a.type ) {
    case _INT_: 
      return a.val==-1; 
    case _ZINT: 
      return (a==gen(-1));
    case _CPLX:
      return (is_minus_one(*a._CPLXptr) && is_zero(*(a._CPLXptr+1)));
    case _DOUBLE_:
      return a._DOUBLE_val==-1;
    case _REAL:
      return evalf_double(a,0,0)._DOUBLE_val==-1;
    case _VECT:
      return a._VECTptr->size()==1 && is_minus_one(a._VECTptr->front());
    case _POLY:
      return Tis_constant(*a._POLYptr) && (is_minus_one(a._POLYptr->coord.front().value));
    case _FRAC:
      return a._FRACptr->num == -a._FRACptr->den;
    case _MOD:
      if (*(a._MODptr+1)==plus_two)
	return is_one(*a._MODptr);
      else
	return is_minus_one(*a._MODptr);
    case _SYMB:
      return a._SYMBptr->sommet==at_neg && is_one(a._SYMBptr->feuille);
    case _USER:
      return a._USERptr->is_minus_one();
    default: 
      return false;
    }
  }

  bool is_sq_minus_one(const gen & a){
    switch (a.type ) {
    case _CPLX: case _MOD: case _USER:
      return is_minus_one(a*a);
    case _VECT:
      return a._VECTptr->size()==1 && is_sq_minus_one(a._VECTptr->front());
    case _POLY:
      return Tis_constant(*a._POLYptr) && (is_sq_minus_one(a._POLYptr->coord.front().value));
    default: 
      return false;
    }
  }

  gen _CPLXgcd(const gen & a,const gen & b){ // a & b must be gen
    if (!is_cinteger(a) || !is_cinteger(b) )
      return plus_one;
    gen acopy(a),bcopy(b),r;
    for (;;){
      if (is_zero(bcopy)){
	complex<double> c=gen2complex_d(acopy);
	double d=arg(c);
	int quadrant=int((2*d)/M_PI);
	switch (quadrant){
	case 0:
	  return acopy;
	case 1:
	  return acopy*(-cst_i);
	case -1:
	  return acopy*cst_i;
	case 2: case -2:
	  return -acopy;
	default:
	  return acopy;
	}
      }
      r=acopy%bcopy;
      acopy=bcopy;
      bcopy=r;
    }
  }

  gen rationalgcd(const gen &a,const gen & b); // in sym2poly.h

  // gcd(undef,x)=x to be used inside series
  gen symgcd(const gen & a,const gen& b){
    if (is_zero(a) || is_undef(a) || (is_one(b)))
      return b;
    if (is_one(a) || is_undef(b) || (is_zero(b)))
      return a;
    if (a==b)
      return a;
    if ( (a.type==_MOD) && (b.type==_MOD) && (a._MODptr->type<=_CPLX) && (b._MODptr->type<= _CPLX) )
      return chkmod(plus_one,a);
    if (a.type==_MOD || b.type==_MOD || a.type==_DOUBLE_ || a.type==_REAL || b.type==_DOUBLE_ || b.type==_REAL )
      return plus_one;
    if ( (a.type==_POLY) && (b.type==_POLY) )
      return gcd(*a._POLYptr,*b._POLYptr);
    if ( (a.type==_EXT) && (b.type ==_EXT) ){
      if ( (*(a._EXTptr+1)!=*(b._EXTptr+1)) || (a._EXTptr->type!=_VECT) || (b._EXTptr->type!=_VECT) )
	return plus_one;
      environment *env=new environment;
      vecteur g=gcd(*a._EXTptr->_VECTptr,*b._EXTptr->_VECTptr,env);
      delete env;
      return ext_reduce(g,*(a._EXTptr+1));
    }
    if ( (a.type==_FRAC) || (b.type==_FRAC))
        return plus_one;
    if (a.type==_EXT){
      if (a._EXTptr->type!=_VECT)
	settypeerr("symgcd");
      if ( (a._EXTptr+1)->type!=_VECT)
	return symgcd(ext_reduce(a),b);
      gen aa(lgcd(*a._EXTptr->_VECTptr));
      gen res=gcd(aa,b);
      vecteur ua((*(a._EXTptr->_VECTptr))/aa),u,uv(*((a._EXTptr+1)->_VECTptr)),v,dd;
      egcd(ua,uv,0,u,v,dd);
      gen dd0(dd.front()),b2(rdiv(b,res));
      simplify(b2,dd0);
      if (is_one(dd0))
	return res*algebraic_EXTension(ua,uv);
      return res;
    }
    if (b.type==_EXT)
      return symgcd(b,a);
    if (a.type==_POLY)
      return gcd(*a._POLYptr,polynome(b,a._POLYptr->dim));
    if (b.type==_POLY)
      return gcd(*b._POLYptr,polynome(a,b._POLYptr->dim));
    if ( (a.type!=_DOUBLE_) && (a.type!=_VECT) && (b.type!=_DOUBLE_) && (b.type!=_VECT) )
      return rationalgcd(a,b);
    return plus_one; // settypeerr("symgcd");
  }

  gen simplify(gen & n, gen & d){
    if ( (d.type==_DOUBLE_) || 
	 ( (d.type==_CPLX) && 
	   ((d._CPLXptr->type==_DOUBLE_) || ((d._CPLXptr+1)->type==_DOUBLE_)) )
	 ){
      gen dd=no_context_evalf(d);
      n=rdiv(no_context_evalf(n),dd); 
      d=plus_one;
      return dd;
    }
    if ( (n.type==_DOUBLE_) ||
	 ( (n.type==_CPLX) && 
	   ((n._CPLXptr->type==_DOUBLE_) || ((n._CPLXptr+1)->type==_DOUBLE_)) )
	 ){
      gen nn=no_context_evalf(n);
      n=plus_one;
      d=rdiv(no_context_evalf(d),nn); 
      return nn*simplify(n,d);
    }
    if (n.type==_FRAC || d.type==_FRAC)
      return plus_one;
    if (is_one(d))
      return d;
    if ( (n.type==_MOD) && (d.type!=_MOD) )
      d=makemod(d,*(n._MODptr+1));
    if (d.type==_MOD){
      if (d._MODptr->is_cinteger()){
	gen dd(d);
	n=n*inv(dd,context0);
	d=makemodquoted(plus_one,*(d._MODptr+1));
	return dd;
      }
    }
    if (is_one(n))
      return n;
    if ((n.type==_POLY) && (d.type==_POLY)){
      polynome np(*n._POLYptr),dp(*d._POLYptr);
      polynome g(np.dim);
      g=simplify(np,dp);
      n=np;
      d=dp;
      return g;
    }
    if (n.type==_VECT){
      if (d.type==_VECT){
	environment * env=new environment;
	vecteur g=gcd(*n._VECTptr,*d._VECTptr,env);
	delete env;
	n=gen(*n._VECTptr/g,_POLY1__VECT);
	d=gen(*d._VECTptr/g,_POLY1__VECT);
	return gen(g,_POLY1__VECT);
      }
      gen nd=_gcd(n);
      gen g=simplify(nd,d);
      n=divvecteur(*n._VECTptr,g);
      return g;
    }
    if (d.type==_VECT){
      gen dd=_gcd(d);
      gen g=simplify(n,dd);
      d=divvecteur(*d._VECTptr,g);
      return g;
    }
    if (d.type==_EXT){
      if ( (d._EXTptr->type==_INT_) || (d._EXTptr->type==_ZINT) ){
	n=n*inv(d,context0);
	gen d_copy=d;
	d=1;
	return d_copy;
      }
      if ( (d._EXTptr+1)->type==_EXT){
	d=ext_reduce(d);
	return simplify(n,d);
      }
      if (d._EXTptr->type==_VECT){
	vecteur u,v,dd;
	if ( (d._EXTptr+1)->type!=_VECT)
	  setsizeerr();
	egcd(*(d._EXTptr->_VECTptr),*((d._EXTptr+1)->_VECTptr),0,u,v,dd);
        gen tmp=algebraic_EXTension(u,*((d._EXTptr+1)->_VECTptr));
	n=n*tmp;
	d=d*tmp;
	return simplify(n,d)*inv_EXT(tmp);
      }
      settypeerr("simplify");
    }
    if (n.type==_EXT){
      gen n_EXT=*n._EXTptr;
      gen g=simplify(n_EXT,d);
      n=algebraic_EXTension(n_EXT,*(n._EXTptr+1));
      return g;
    }
    if (n.type==_POLY) {
      polynome np(*n._POLYptr),dp(d,n._POLYptr->dim);
      polynome g(np.dim);
      g=simplify(np,dp);
      n=np;
      d=dp;
      return g;
    }
    if (d.type==_POLY){
      polynome np(n,d._POLYptr->dim),dp(*d._POLYptr);
      polynome g(np.dim);
      g=simplify(np,dp);
      n=np;
      d=dp;
      return g;
    }
    vecteur l(lvar(n));
    lvar(d,l);
    gen num=e2r(n,l,context0),den=e2r(d,l,context0),g=gcd(num,den); // ok
    den=rdiv(den,g);
    if (is_zero(re(den,context0))){ //ok
      den=cst_i*den;
      g=-cst_i*g;
    }
    if (is_positive(-den,context0)){  // ok
      den=-den;
      g=-g;
    }
    n=r2sym(rdiv(num,g),l,context0); // ok
    d=r2sym(den,l,context0); // ok
    return r2sym(g,l,context0); // ok
  }

  gen gcd(const gen & a,const gen & b){
    mpz_t * res;
    switch ( (a.type<< _DECALAGE) | b.type ) {
    case _INT___INT_: 
      return(gcd(a.val,b.val)); 
    case _INT___ZINT: 
      if (a.val)
	return(int(mpz_gcd_ui(NULL,*b._ZINTptr,absint(a.val))));
      else
	return(b);
    case _ZINT__INT_:
      if (b.val)
	return(int(mpz_gcd_ui(NULL,*a._ZINTptr,absint(b.val))));
      else
	return(a);
    case _ZINT__ZINT:
      res=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init(*res);
      mpz_gcd(*res,*a._ZINTptr,*b._ZINTptr);
      return(res);
    case _INT___CPLX: case _ZINT__CPLX:
    case _CPLX__INT_: case _CPLX__ZINT:
    case _CPLX__CPLX:
      return _CPLXgcd(a,b);
    case _POLY__POLY:
      return gcd(*a._POLYptr,*b._POLYptr);
    case _VECT__VECT:
      return gen(gcd(*a._VECTptr,*b._VECTptr,0),_POLY1__VECT);
    case _FRAC__FRAC:
      return fraction(gcd(a._FRACptr->num,b._FRACptr->num),lcm(a._FRACptr->den,b._FRACptr->den));
    default:
      if (a.type==_FRAC)
	return fraction(gcd(a._FRACptr->num,b),a._FRACptr->den);
      if (b.type==_FRAC)
	return fraction(gcd(b._FRACptr->num,a),b._FRACptr->den);
      if (a.type==_USER)
	return a._USERptr->gcd(b);
      if (b.type==_USER)
	return b._USERptr->gcd(a);
      return symgcd(a,b);
    }
  }

  gen lcm(const gen & a,const gen & b){
    return normal(rdiv(a,gcd(a,b)),context0)*b; // ok
  }

  void ciegcd(const gen &a_orig,const gen &b_orig, gen & u,gen &v,gen &d ){
    gen a(a_orig),b(b_orig),au(plus_one),bu(zero),q,r,ru;
    while (!is_zero(b)){
      q=iquo(a,b);
      r=a-b*q;
      a=b;
      b=r;
      ru=au-bu*q;
      au=bu;
      bu=ru;
    }
    u=au;
    d=a;
    v=iquo(d-a_orig*u,b_orig);
  }

  void egcd(const gen &ac,const gen &bc, gen & u,gen &v,gen &d ){
    gen a(ac),b(bc);
    switch ( (a.type<< _DECALAGE) | b.type ) {
    case _INT___INT_: case _INT___ZINT: case _ZINT__INT_: case _ZINT__ZINT: 
      if (a.type==_INT_)
	a.uncoerce();
      if (b.type==_INT_)
	b.uncoerce();
      if (!u.type)
	u.uncoerce();
      if (!v.type)
	v.uncoerce();
      if (!d.type)
	d.uncoerce();
      mpz_gcdext(*d._ZINTptr,*u._ZINTptr,*v._ZINTptr,*a._ZINTptr,*b._ZINTptr);
      break;
    default: 
      ciegcd(a,b,u,v,d);
      break;
    }
  }

  void _ZINTinvmod(const gen & a,const gen & modulo, mpz_t * & res){
    mpz_t *aptr,*bptr;
    if (a.type!=_INT_)
      aptr=a._ZINTptr;
    else {
      aptr=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init_set_si(*aptr,a.val);
    }
    if (modulo.type)
      bptr=modulo._ZINTptr;
    else {
      bptr=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init_set_si(*bptr,modulo.val);
    }
    res=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
    mpz_init(*res);
    bool ok=mpz_invert(*res,*aptr,*bptr);
    if (a.type==_INT_){
      mpz_clear(*aptr);
      free(aptr);
    }
    if (!modulo.type){
      mpz_clear(*bptr);
      free(bptr);
    }
    if (!ok)
      setsizeerr("Not invertible "+a.print(context0)+" mod "+modulo.print(context0));
  }

  gen invmod(const gen & a,const gen & modulo){
    if (a.type==_CPLX){
      gen r=re(a,context0),i=im(a,context0); // ok
      gen n=invmod(r*r+i*i,modulo);
      return smod(r*n,modulo)-cst_i*smod(i*n,modulo);
    }
    mpz_t * res;
    switch ( (a.type<< _DECALAGE) | modulo.type) {
    case _INT___INT_: 
      return(invmod(a.val,modulo.val));
    case _INT___ZINT: case _ZINT__INT_: case _ZINT__ZINT:
      _ZINTinvmod(a,modulo,res);
      return(res);
    default: 
      settypeerr("invmod");
    }
    return 0;
  }

  gen fracmod(const gen & a_orig,const gen & modulo){
    // write a as p/q with |p| and |q|<sqrt(modulo/2)
    gen a(a_orig),m(modulo);
    if (a.type==_INT_)
      a.uncoerce();
    if (m.type==_INT_)
      m.uncoerce();
    if ( (a.type!=_ZINT) || (m.type!=_ZINT) )
      setsizeerr();
    mpz_t u,d,u1,d1,absd1,sqrtm,q,ur,r,tmp;
    mpz_init_set_si(u,0);
    mpz_init_set(d,*m._ZINTptr);
    mpz_init_set_si(u1,1);
    mpz_init_set(d1,*a._ZINTptr);
    mpz_init(absd1);
    mpz_init(sqrtm);
    mpz_init(q);
    mpz_init(ur);
    mpz_init(r);
    mpz_init(tmp);
    mpz_tdiv_q_2exp(q,*m._ZINTptr,1);
    mpz_sqrt(sqrtm,q);
    // int signe;
    for (;;){
      mpz_abs(absd1,d1);
      if (mpz_cmp(absd1,sqrtm)<=0)
	break;
      mpz_fdiv_qr(q,r,d,d1);
      // u-q*u1->ur, v-q*v1->vr
      mpz_mul(tmp,q,u1);
      mpz_sub(ur,u,tmp);
      // u1 -> u, ur -> u1 ; v1 -> v, vr -> v1, d1 -> d, r -> d1
      mpz_set(u,u1);
      mpz_set(u1,ur);
      mpz_set(d,d1);
      mpz_set(d1,r);
    }
    // u1*a+v1*m=d1 -> a=d1/u1 modulo m
    gen num(d1);
    gen den(u1);
    mpz_set(d,*m._ZINTptr);
    mpz_gcd(u,d,u1);
    if (mpz_cmp_ui(u,1)!=0)
      setsizeerr("Reconstructed denominator is not prime with modulo");
    mpz_clear(u);
    mpz_clear(d);
    mpz_clear(u1);
    mpz_clear(d1);
    mpz_clear(absd1);
    mpz_clear(sqrtm);
    mpz_clear(q);
    mpz_clear(ur);
    mpz_clear(r);
    mpz_clear(tmp);
    if (is_positive(den,context0)) // ok
      return fraction(num,den);
    else
      return fraction(-num,-den);
  }

  void _ZINTpowmod(const gen & base,const gen & expo,const gen & modulo, mpz_t * & res){
    mpz_t *aptr,*bptr;
    if (base.type)
      aptr=base._ZINTptr;
    else { 
      aptr=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init_set_si(*aptr,base.val);
    }
    if (modulo.type)
      bptr=modulo._ZINTptr;
    else {
      bptr=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init_set_si(*bptr,modulo.val);
    }
    res=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
    mpz_init(*res);
    if (!expo.type)
      mpz_powm_ui(*res,*aptr,expo.val,*modulo._ZINTptr);
    else
      mpz_powm (*res,*aptr,*expo._ZINTptr,*bptr);
    if (!base.type){
      mpz_clear(*aptr);
      free(aptr);
    }
    if (!modulo.type){
      mpz_clear(*bptr);
      free(bptr);
    }
  }

  gen powmod(const gen &base,const gen & expo,const gen & modulo){
    if (base.type==_VECT){
      const_iterateur it=base._VECTptr->begin(),itend=base._VECTptr->end();
      vecteur res;
      for (;it!=itend;++it)
	res.push_back(powmod(*it,expo,modulo));
      return gen(res,base.subtype);
    }
    if ((expo.type!=_INT_) && (expo.type!=_ZINT))
      setsizeerr(); // exponent must be a _DOUBLE_ integer
    if (!is_positive(expo,context0)) // ok 
      return(powmod(invmod(base,modulo),-expo,modulo));
    if (modulo.type==_INT_){
      // try converting base to int and expo to a long
      gen mybase(base % modulo);
      if ( (expo.type==_INT_) && (mybase.type==_INT_) ){
	unsigned long tmp=expo.val;
	return powmod(mybase.val,tmp,modulo.val);
      }
    }
    mpz_t * res;
    switch ( (base.type<< _DECALAGE) | modulo.type) {
    case _INT___INT_: case _INT___ZINT: case _ZINT__INT_: case _ZINT__ZINT:
      _ZINTpowmod(base,expo,modulo,res);
      return(res);
    default: 
      settypeerr("powmod");
    }
    return 0;
  }

  // assuming amod and bmod are prime together, find c such that
  // c = a mod amod  and c = b mod bmod
  // hence a + A*amod = b + B*bmod
  // or A*amod -B*bmod = b - a
  gen ichinrem(const gen & a,const gen &b,const gen & amod, const gen & bmod){
    gen A,B,d,q;
    egcd(amod,bmod,A,B,d);
    if (is_one(d))
      q=b-a;
    else
      if (!is_zero(irem(b-a,d,q)))
	setsizeerr("No Integer Solution");
    A=A*q;
    return smod(A*amod+a,amod*bmod);
  }

  gen isqrt(const gen & a){
    if ( (a.type!=_INT_) && (a.type!=_ZINT))
      settypeerr("isqrt");
    mpz_t *aptr;
    if (a.type!=_INT_)
      aptr=a._ZINTptr;
    else {
      aptr=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init_set_si(*aptr,a.val);
    }
    mpz_t *res=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
    mpz_init(*res);
    mpz_sqrt(*res,*aptr);
    if (a.type==_INT_){
      mpz_clear(*aptr);
      free(aptr);
    }
    return(res);
  }

  int is_perfect_square(const gen & a){
    if ( (a.type!=_INT_) && (a.type!=_ZINT))
      return false;
    mpz_t *aptr;
    if (a.type!=_INT_)
      aptr=a._ZINTptr;
    else {
      aptr=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init_set_si(*aptr,a.val);
    }
    int res= mpz_perfect_square_p(*aptr);
    if (a.type==_INT_){
      mpz_clear(*aptr);
      free(aptr);
    }
    return res;
  }

  int is_probab_prime_p(const gen & a){
    if ( (a.type!=_INT_) && (a.type!=_ZINT))
      settypeerr("is_probab_prime_p");
    mpz_t *aptr;
    if (a.type!=_INT_)
      aptr=a._ZINTptr;
    else {
      aptr=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init_set_si(*aptr,a.val);
    }
    int res= mpz_probab_prime_p(*aptr,TEST_PROBAB_PRIME);
    if (a.type==_INT_){
      mpz_clear(*aptr);
      free(aptr);
    }
    return res;
  }

  gen nextprime(const gen & a){
    if ( (a.type!=_INT_) && (a.type!=_ZINT))
      settypeerr("nextprime");
    gen res(a);
    if (is_zero(smod(res,plus_two)))
      res=res+1;
    for ( ; ; res=res+2)
      if (is_probab_prime_p(res))
	return(res);
  }

  gen prevprime(const gen & a){
    if ( (a.type!=_INT_) && (a.type!=_ZINT))
      settypeerr("prevprime");
    gen res(a);
    if (is_zero(smod(res,plus_two)))
      res=res-1;
    for ( ; ; res=res-2)
      if (is_probab_prime_p(res))
	return(res);
    return zero;
  }

  int jacobi(const gen & a, const gen &b){
    if ( (a.type!=_INT_) && (a.type!=_ZINT))
      settypeerr("jacobi");
    if ( (b.type!=_INT_) && (b.type!=_ZINT))
      settypeerr("jacobi");
    mpz_t *aptr,*bptr;
    if (a.type!=_INT_)
      aptr=a._ZINTptr;
    else {
      aptr=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init_set_si(*aptr,a.val);
    }
    if (b.type!=_INT_)
      bptr=b._ZINTptr;
    else {
      bptr=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init_set_si(*bptr,b.val);
    }
    int res= mpz_jacobi(*aptr,*bptr);
    if (a.type==_INT_){
      mpz_clear(*aptr);
      free(aptr);
    }
    if (b.type==_INT_){
      mpz_clear(*bptr);
      free(bptr);
    }
    return res;
  }

  int legendre(const gen & a, const gen & b){
    if ( (a.type!=_INT_) && (a.type!=_ZINT))
      settypeerr("legendre");
    if ( (b.type!=_INT_) && (b.type!=_ZINT))
      settypeerr("legendre");
    mpz_t *aptr,*bptr;
    if (a.type!=_INT_)
      aptr=a._ZINTptr;
    else {
      aptr=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init_set_si(*aptr,a.val);
    }
    if (b.type!=_INT_)
      bptr=b._ZINTptr;
    else {
      bptr=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
      mpz_init_set_si(*bptr,b.val);
    }
    int res=mpz_legendre(*aptr,*bptr);
    if (a.type==_INT_){
      mpz_clear(*aptr);
      free(aptr);
    }
    if (b.type==_INT_){
      mpz_clear(*bptr);
      free(bptr);
    }
    return res;
  }

  bool has_denominator(const gen & n){
    switch (n.type ) {
    case _INT_: case _ZINT: case _CPLX: case _DOUBLE_: case _IDNT: case _EXT: case _POLY: case _MOD: case _USER: case _REAL:
      return false;
    case _SYMB: case _FRAC:
      return true;
    default:
      settypeerr("has_denominator");
    }
    return 0;
  }

  // Note that this function should be optimized for large input
  string cut_string(const string & chaine,int nchar,vector<int> & ligne_end) {
    // cerr << clock() << endl;
    int pos;
    if (ligne_end.empty())
      pos=0;
    else
      pos=ligne_end.back()+1;
    int l=chaine.size();
    string res;
    for (int i=0;i<l;){
      // look for \n between i and l
      int k=chaine.find_first_of('\n',i);
      if ( (l-i<nchar) && ((k<i)||(k>=l-1)) ){
	ligne_end.push_back(pos+l);
	// cerr << clock() << endl;
	return res+chaine.substr(i,l-i);
      }
      if ((k>=i) && (k<i+nchar+4*(i==0)) ){
	ligne_end.push_back(pos+k);
	res += chaine.substr(i,k+1-i);
	i=k+1;
      }
      else {
	int j;
	int j1=chaine.find_last_of('+',i+nchar+4*(i==0));
	int j2=chaine.find_last_of('-',i+nchar+4*(i==0));
	int j3=chaine.find_last_of(',',i+nchar+4*(i==0));
	j=max(j1,max(j2,j3));
	if ((j-i)<(nchar/2))
	  j=i+nchar+4*(i==0);
	ligne_end.push_back(pos+min(j,l));
	res += chaine.substr(i,j-i);
	i=j;
	if (i<l){
	  res +="\\\n     ";
	  pos +=7;
	}
      }
    }
    // cerr << clock() << endl;
    return res;
  }

  string calc_endlines_positions(const vecteur & history_in,const vecteur & history_out,int nchar,vector<int> & endlines,vector<int> & positions){
    string res;
    endlines.clear();
    positions.clear();
    int s_in=history_in.size(),s_out=history_out.size();
    int s=max(s_in,s_out);
    for (int i=0;i<s;++i){
      string chaine;
      if (rpn_mode)
	chaine=print_INT_(s-i)+": ";
      else
	chaine=print_INT_(i)+": ";
      if (!rpn_mode){
	if (i<s_in)
	  chaine+=history_in[i].print(context0)+" = ";
      }
      else
	chaine +="   ";
      if (i<s_out)
	chaine += history_out[i].print(context0);
      if (i)
	res +='\n';
      res += cut_string(chaine,nchar,endlines);
      positions.push_back(endlines.back());
    }
    return res;
  }

  bool is_operator_char(char c){
    switch(c){
    case '+': case '-': case '*': case '/': case '^': case '%':
      return true;
    }
    return false;
  }

  bool is_operator_char(char c,char op){
    switch(c){
    case '+': case '-': 
      return true;
    case '*': case '/': case '^': case '%':
      return c==op;
    }
    return false;
  }

  bool matchpos(const string & s,int & pos){
    char c=s[pos];
    char c_orig=c;
    int l=s.size();
    int counter1=0,counter2=0,counter3=0,incr;
    if ( (c==')') || (c==']') || (c=='}') )
      incr=-1;
    else
      incr=1;
    for (;(pos>=0) && (pos<l);pos+=incr){
      switch (c=s[pos]){
      case '(':
	counter1++;
	break;
      case ')':
	counter1--;
	break;
      case '[':
	counter2++;
	break;
      case ']':
	counter2--;
	break;
      case '{':
	counter3++;
	break;
      case '}':
	counter3--;
	break;
      }
      if ( (!counter1) && (!counter2) && (!counter3) ){
	bool res=false;
	switch (c_orig){
	case '(':
	  res=c==')';
	  break;
	case '[':
	  res=c==']';
	  break;
	case '{':
	  res=c=='}';
	  break;
	case ')':
	  res=c=='(';
	  break;
	case ']':
	  res=c=='[';
	  break;
	case '}':
	  res=c=='{';
	  break;
	}
	return res;
      }
    }
    return false;
  }

  void find_left(const string & s,int & pos1,int & pos2){
    int l=s.size();
    pos1=min(max(pos1,0),l);
    int pos1orig=pos1;
    if (!pos1)
      return;
    pos2=max(min(pos2,l),0);
    int counter1=0,counter2=0;
    if (pos2==l){
      int i=pos2;
      for (;i>pos1;){
	--i;
	char ch = s[i];
	if (ch=='(')
	  ++counter1;
	if (ch==')')
	  --counter1;
	if (ch=='[')
	  ++counter2;
	if (ch==']')
	  --counter2;
      }
    }
    for (;pos1>=0;--pos1){
      char ch=s[pos1];
      if ( (!counter1) && (!counter2) && ( (ch=='(') || (ch=='[') || (ch=='+') || (ch=='-') || (ch==',')  )){
	if ( (pos1<pos1orig) && ( (ch!='(') || (s[pos2-1]!=')') ) )
	  ++pos1;
	break;
      }
      if (ch=='('){
	++counter1;
	if ( (!counter1) && (!counter2) ){
	  if (s[pos2-1]==')')
	    break;
	  if ( pos1 && isalphan(s[pos1-1])){
	    --pos1;
	    for (;pos1>=0;--pos1)
	      if (!isalphan(s[pos1]))
		break;
	    ++pos1;
	  }
	  break;
	}
      }
      if (ch==')')
	--counter1;
      if (ch=='['){
	++counter2;
	if ( (!counter1) && (!counter2) ){
	  if (s[pos2-1]==']')
	    break;
	  if ( pos1 && isalphan(s[pos1-1])){
	    --pos1;
	    for (;pos1>=0;--pos1)
	      if (!isalphan(s[pos1]))
		break;
	    ++pos1;
	  }
	  break;
	}
      }
      if (ch==']')
	--counter2;
    }
  }

  void find_right(const string & s,int & pos1,int & pos2){
    int l=s.size();
    pos1=min(max(pos1,0),l);
    pos2=max(min(pos2,l),0);
    int pos2orig=pos2;
    int counter1=0,counter2=0;
    for (int i=pos1;(i<pos2-1) && (i<l);++i){
      char ch=s[i];
      if (ch=='(')
	++counter1;
      if (ch==')')
	--counter1;
      if (ch=='[')
	++counter2;
      if (ch==']')
	--counter2;
      if ( (counter1<0) && (pos1) ){ // restart at an earlier position
	pos1=s.find_last_of('(',pos1-1);
	if (pos1<0)
	  pos1=0;
	i=pos1-1;
	counter1=0;
	counter2=0;
      }
    }      
    for (;pos2<=l;++pos2){
      char ch=s[pos2-1];
      if ( (!counter1) && (!counter2) && ( (ch==')') ||  (ch==']') || (ch=='+') || (ch=='-') || (ch==',')  ) && (pos2>pos2orig)){
	--pos2;
	break;
      }
      if (ch=='(')
	++counter1;
      if (ch==')'){
	--counter1;
	if ( (!counter1) && (!counter2) ){
	  if ( (pos1>0) && (s[pos1]=='(') && isalphan(s[pos1-1])){
	    --pos1;
	    for (;pos1>=0;--pos1)
	      if (!isalphan(s[pos1]))
		break;
	    ++pos1;
	  }
	  break;
	}
      }
      if (ch=='[')
	++counter2;
      if (ch==']'){
	--counter2;
	if ( (!counter1) && (!counter2) ){
	  if ( (pos1>0) && (s[pos1]=='[') && isalphan(s[pos1-1])){
	    --pos1;
	    for (;pos1>=0;--pos1)
	      if (!isalphan(s[pos1]))
		break;
	    ++pos1;
	  }
	  break;
	}
      }
    }
    if (pos2==l+1)
      find_left(s,pos1,pos2);
  }

  void increase_selection(const string & s,int & pos1,int& pos2){
    int l=s.size();
    int pos1_orig(pos1),pos2_orig(pos2);
    // adjust selection (does not change anything on a valid selection)
    find_left(s,pos1,pos2);
    find_right(s,pos1,pos2);
    if ( (pos1!=pos1_orig) || (pos2!=pos2_orig) )
      return;
    if (pos1 && (pos2<l) && ( (s[pos1-1]=='(') || (s[pos1-1]==',')) && (s[pos2]!=')') && (s[pos2]!=',')){
      ++pos2;
      find_right(s,pos1,pos2);
      return;
    }
    if (pos1>1){
      char op=s[pos1-1];
      --pos1;
      for (;pos1;--pos1){
	if (s[pos1]==',')
	  op=0;
	if (!is_operator_char(s[pos1],op))
	  break;
      }
      if (s[pos1]=='(' && pos1){
	--pos1;
	for (;pos1;--pos1){
	  if (!isalphan(s[pos1]))
	    break;
	}
	++pos1;
      }
      find_left(s,pos1,pos2);
      find_right(s,pos1,pos2);
      return;
    }
    pos1=0;
    ++pos2;
    find_right(s,pos1,pos2);
  }

  void decrease_selection(const string & s,int & pos1,int& pos2){
    int l=s.size();
    int pos2_orig(pos2);
    // adjust selection (does not change anything on a valid selection)
    find_left(s,pos1,pos2);
    if (pos2!=l)
      --pos2;
    if (!pos2)
      return;
    int counter1=0,counter2=0;
    char op=' ';
    if (pos2<l-1)
      op=s[pos2+1];
    for (;pos2>pos1;--pos2){
      char ch=s[pos2];
      if (ch=='('){
	++counter1;
	if ( (!counter1) && (!counter2) && pos2_orig && (s[pos2_orig-1]==')') ){
	  pos1=pos2+1;
	  pos2=pos2_orig-1;
	  return;
	}
      }
      if (ch==')')
	--counter1;
      if (ch=='[')
	++counter2;
      if (ch==']')
	--counter2;
      if (ch==',')
	op=0;
      if ( (!counter1) && (!counter2) && ( is_operator_char(ch,op) || (ch==',')) )
	return;
    }
    for (;pos1<l;++pos1){
      char ch=s[pos1];
      if ( (ch=='(') ||  (ch=='[') || (ch=='+') || (ch==','))
	break;
    }
    ++pos1;
    pos2=pos1+1;
    find_right(s,pos1,pos2);
  }

  void move_selection_right(const string & s,int & pos1, int & pos2){
    int l=s.size();
    // int pos1_orig(pos1),pos2_orig(pos2);
    // adjust selection (does not change anything on a valid selection)
    find_right(s,pos1,pos2);
    pos1=pos2;
    char op=s[pos1];
    for (;pos1<l;++pos1){
      if (s[pos1]==',')
	op=0;
      if (!is_operator_char(s[pos1],op) && (s[pos1]!=')') && (s[pos1]!=']'))
	break;
    }
    pos2=pos1+1;
    find_right(s,pos1,pos2);
  }

  void move_selection_left(const string & s,int & pos1, int & pos2){
    // int l=s.size();
    // int pos1_orig(pos1),pos2_orig(pos2);
    // adjust selection (does not change anything on a valid selection)
    find_left(s,pos1,pos2);
    pos2=pos1-1;
    char op=s[pos2];
    for (;pos2>0;--pos2){
      if (s[pos2-1]==',')
	op=0;
      if (!is_operator_char(s[pos2-1],op) && (s[pos2-1]!='(') && (s[pos2-1]!='[') )
	break;
    }
    if (pos2<=0){
      pos1=0;
      pos2=0;
      return;
    }
    pos1=pos2-1;
    find_left(s,pos1,pos2);
    find_right(s,pos1,pos2);
  }

  string remove_extension(const string & chaine){
    int s=chaine.size();
    int l=chaine.find_last_of('.',s);
    int ll=chaine.find_last_of('/',s);
    if (l>0 && l<s){
      if ( ll<=0 || ll>=s || l>ll)
	return chaine.substr(0,l);
    }
    return chaine;
  }

  //environment * env=new environment;

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC

