// -*- mode:C++ ; compile-command: "g++ -I.. -g -c gen.cc" -*-
/*
 *  Copyright (C) 2001 B. Parisse, Institut Fourier, 38402 St Martin d'Heres
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
#ifndef _GIAC_GEN_H
#define _GIAC_GEN_H
#include "first.h"

// #include <gmp.h>
#ifdef HAVE_GMPXX_H
#include <gmpxx.h>
#endif
#ifdef HAVE_LIBMPFR
#include <mpfr.h>
// #include <mpf2mpfr.h>
#endif
#ifdef HAVE_LIBMPFI
#include <mpfi.h>
#endif
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "dispatch.h"
#include "vecteur.h"
#include "fraction.h"
#include <complex>
#include <stdlib.h>

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC


  class gen ;   // a class around gmp mpz_t to simplify usage in C++
  class identificateur; // global name
  struct symbolic; // symbolics
  class unary_function_ptr; // functions
  struct real_complex_rootof; 
  template <class T> class tensor;
  typedef tensor<gen> polynome;
  typedef dbgprint_vector<polynome> vectpoly;
  typedef Tfraction<gen> fraction;
  struct eqwdata ; // type for EQW display
  struct grob ; // type for graphic object
  class gen_user ; // user defined type
  typedef std::vector< facteur< polynome > > factorization;

  // errors
  void settypeerr(GIAC_CONTEXT0);
  void setsizeerr(GIAC_CONTEXT0);
  void setdimerr(GIAC_CONTEXT0);
  void settypeerr(const std::string & s);
  void setsizeerr(const std::string & s);
  void setdimerr(const std::string & s);
  void divisionby0err(const gen &,GIAC_CONTEXT0);
  void cksignerr(const gen &,GIAC_CONTEXT0);
  void invalidserieserr(const std::string &,GIAC_CONTEXT0);
  void toofewargs(const std::string & s,GIAC_CONTEXT0);
  void toomanyargs(const std::string & s,GIAC_CONTEXT0);
  void maxordererr(GIAC_CONTEXT0);
  void setstabilityerr(GIAC_CONTEXT0);

  // short integer arithmetic
  int absint(int a);
  int min(int a,int b);
  int max(int a,int b);
  int invmod(int n,int modulo);
  unsigned invmod(unsigned a,int b);
  int invmod(longlong a,int b);
  int powmod(int a,unsigned long n,int m);
  int gcd(int a,int b);
  int smod(int a,int b); // where b is assumed to be positive
  int simplify(int & a,int & b);

  // arbitrary precision floats hierarchy (value or interval)
  std::string printmpf_t(const mpf_t & inf);
  class real_object {
  public:
#ifdef HAVE_LIBMPFR
    mpfr_t inf;
#else
    mpf_t inf;
#endif
    real_object(double d); 
#ifdef HAVE_LIBMPFR
    real_object(const mpfr_t & d); 
    real_object(const mpf_t & d); 
#else
    real_object(const mpf_t & d); 
#endif
    real_object(const gen & g);
    real_object(const gen & g,unsigned int precision);
    real_object() ;
    virtual std::string print(GIAC_CONTEXT) const;
    void dbgprint() const { std::cerr << this->print(0) << std::endl; }
    virtual ~real_object() { 
#ifdef HAVE_LIBMPFR
      mpfr_clear(inf);
#else
      mpf_clear(inf); 
#endif
    }
    virtual real_object & operator = (const real_object & g);
    real_object (const real_object & g) ;
    gen operator + (const gen & g) const;
    virtual real_object operator + (const real_object & g) const;
    gen operator * (const gen & g) const;
    virtual real_object operator * (const real_object & g) const;
    gen operator / (const gen & g) const;
    virtual real_object operator / (const real_object & g) const;
    gen operator - (const gen & g) const;
    virtual real_object operator - (const real_object & g) const;
    virtual real_object operator -() const;
    virtual real_object inv() const;
    virtual real_object sqrt() const;
    virtual real_object abs() const;
    virtual real_object exp() const;
    virtual real_object log() const;
    virtual real_object sin() const;
    virtual real_object cos() const;
    virtual real_object tan() const;
    virtual real_object sinh() const;
    virtual real_object cosh() const;
    virtual real_object tanh() const;
    virtual real_object asin() const;
    virtual real_object acos() const;
    virtual real_object atan() const;
    virtual real_object asinh() const;
    virtual real_object acosh() const;
    virtual real_object atanh() const;
    virtual bool is_zero();
    virtual bool is_inf();
    virtual bool is_nan();
    virtual bool is_positive();
    virtual double evalf_double() const;
  };
  struct ref_real_object {
    int ref_count;
    real_object r;
  };
  gen real2int(const gen & g,GIAC_CONTEXT);
  gen real2double(const gen & g);
  class real_interval : public real_object {
  public:
#ifdef HAVE_LIBMPFI
    mpfi_t infsup;
#else
#ifdef HAVE_LIBMPFR
    mpfr_t sup;
#else
    mpf_t sup;
#endif
#endif
    real_interval(const real_object & r):real_object(r) { 
#ifdef HAVE_LIBMPFI
      mpfi_init_set_fr(infsup,r.inf);
#else
#ifdef HAVE_LIBMPFR
      mpfr_init_set(sup,r.inf,GMP_RNDN); 
#else
      mpf_init_set(sup,r.inf); 
#endif
#endif
    }
    real_interval(const real_interval & r):real_object(r) { 
#ifdef HAVE_LIBMPFI
      mpfi_init_set(infsup,r.infsup);
#else
#ifdef HAVE_LIBMPFR
      mpfr_init_set(sup,r.sup,GMP_RNDN); 
#else
      mpf_init_set(sup,r.sup); 
#endif
#endif
    }
    virtual ~real_interval() { 
#ifdef HAVE_LIBMPFI
      mpfi_clear(infsup); 
#else
#ifdef HAVE_LIBMPFR
      mpfr_clear(sup); 
#else
      mpf_clear(sup); 
#endif
#endif
    }
    virtual real_object & operator = (const real_interval & g) ;
    virtual real_object & operator = (const real_object & g) ;
    virtual real_object operator + (const real_object & g) const;
    virtual real_interval operator + (const real_interval & g) const;
    virtual real_object operator * (const real_object & g) const;
    virtual real_interval operator * (const real_interval & g) const;
    virtual real_object operator - (const real_object & g) const;
    virtual real_interval operator - (const real_interval & g) const ;
    virtual real_object operator -() const;
    virtual real_object inv() const;
  };
  std::string print_binary(const real_object & r);
  gen read_binary(const std::string & s,unsigned int precision);
  // Convert g to a real or complex object of precision nbits
  gen accurate_evalf(const gen & g,int nbits);
  vecteur accurate_evalf(const vecteur & v,int nbits);

  typedef std::map<gen,gen,const std::pointer_to_binary_function < const gen &, const gen &, bool> > gen_map;
  class my_mpz;
  struct gen_ptr_t {
      // pointer types with a ref_count
      union {
	mpz_t * __ZINTptr; // long int (type _ZINT)
	real_object * __REALptr; // extended double (type _REAL)
	gen * __CPLXptr ; // complex as an gen[2] array (type _CPLX)
	identificateur * __IDNTptr; // global name identifier (type _IDNT)
	symbolic * __SYMBptr; // for symbolic objects (type _SYMB)
	gen * __MODptr;
	gen * __EXTptr; // 2 gens for alg. extension (type ext)
	// alg ext: 1st gen is a std::vector or a fraction, 2nd gen is
	// a/ a std::vector, the minimal monic polynomial (the roots are permutable)
	// b/ a real_complex_rootof given by it's min poly and 
	// c/ another type meaning that the root is expressed in terms
	//    of another rootof, in this case ext_reduce should be called
	// For 2nd order extension, X^2=d is used if d!=1 mod 4
	// X is the positive solution
	// if d=1 mod 4 the equation is X^2-X=(d-1)/4
	fraction * __FRACptr; // fraction (type _FRAC)
	polynome * __POLYptr ; // multidim. sparse polynomials (type poly)
	// _VECTosite types (std::vector<>)
	vecteur * __VECTptr ; // vecteur: std::vectors & dense_POLY1 (type _VECT)
	sparse_poly1 * __SPOL1ptr ; // std::vector<monome>: sparse 1-d poly (type _SPOL1)
	std::string * __STRNGptr;
	unary_function_ptr * __FUNCptr;
	real_complex_rootof * __ROOTptr;
	gen_user * __USERptr;
	gen_map * __MAPptr;
	eqwdata * __EQWptr;
	grob * __GROBptr;
	void * __POINTER_val;
      };
      int * ref_count;
  };

#define  _ZINTptr ptr_val.__ZINTptr
#define	 _REALptr ptr_val.__REALptr
#define  _CPLXptr ptr_val.__CPLXptr
#define  _IDNTptr ptr_val.__IDNTptr
#define  _SYMBptr ptr_val.__SYMBptr
#define  _MODptr ptr_val.__MODptr
#define  _FRACptr ptr_val.__FRACptr
#define  _EXTptr ptr_val.__EXTptr
#define  _POLYptr ptr_val.__POLYptr 
#define  _VECTptr  ptr_val.__VECTptr
#define  _SPOL1ptr ptr_val.__SPOL1ptr
#define  _STRNGptr ptr_val.__STRNGptr
#define  _FUNCptr ptr_val.__FUNCptr
#define  _ROOTptr ptr_val.__ROOTptr
#define  _USERptr ptr_val.__USERptr
#define  _MAPptr ptr_val.__MAPptr
#define  _EQWptr ptr_val.__EQWptr
#define  _GROBptr ptr_val.__GROBptr
#define  _POINTER_val ptr_val.__POINTER_val

  class gen {
  public:
    short int subtype;
    short int type;  // see dispatch.h
    union {
      // atomic types
      int val; // immediate int (type _INT_)
      double _DOUBLE_val; // immediate float (type _DOUBLE_)
      gen_ptr_t ptr_val;
    };
    gen(): subtype(0), type(_INT_),  val(0) {
#ifdef COMPILE_FOR_STABILITY
      control_c();
#endif
    };
    gen(void *ptr,short int subt): subtype(subt),type(_POINTER_) {
      _POINTER_val=ptr; ptr_val.ref_count=new int(1);
#ifdef COMPILE_FOR_STABILITY
      control_c();
#endif
    };
    gen(int i): subtype(0), type(_INT_), val(i) {
#ifdef COMPILE_FOR_STABILITY
      control_c();
#endif
    };
    gen(size_t i): subtype(0), type(_INT_), val((int)i) {
#ifdef COMPILE_FOR_STABILITY
      control_c();
#endif
    };
    gen(longlong i);
    gen(const mpz_t & m);
    // WARNING coerce *mptr to an int if possible, in this case delete mptr
    // Pls do not use this constructor unless you know exactly what you do!!
    gen(mpz_t * mptr);
    gen (double d): type(_DOUBLE_),_DOUBLE_val(d) {};
    gen(int a,int b);
    gen(double a,double b);
    gen(const gen & a,const gen & b);
    gen(const std::complex<double> & c);
    gen(const gen & e);
    gen (const identificateur & s);
    gen (const vecteur & v,short int s=0);
    gen (vecteur * vptr,short int s=0); 
    // vptr must be a pointer allocated by new, do not delete it explicitly
    gen (const symbolic & s);
    gen (symbolic * sptr);
    gen (const gen_user & g);
    gen (const real_object & g);
    // Pls do not use this constructor unless you know exactly what you do
    gen(polynome * pptr);
    gen (const polynome & p);
    gen (const fraction & p);
    gen (const real_complex_rootof & r);
    gen (const std::string & s,GIAC_CONTEXT);
    gen (const std::string & s,const vecteur & l,GIAC_CONTEXT);
    gen (const sparse_poly1 & p);
    gen (const unary_function_ptr & f,int nargs=1);
    gen (const gen_map & m);
    gen (const eqwdata & );
    gen (const grob & );
#ifdef HAVE_GMPXX_H
    gen (const mpz_class &);
#endif
    gen (const my_mpz &);
    ~gen();
    bool in_eval(int level,gen & evaled,const context * contextptr) const;
    gen eval(int level,const context * contextptr) const;
    // inline gen eval() const { return eval(DEFAULT_EVAL_LEVEL,context0); }
    bool in_evalf(int level,gen & evaled,const context * contextptr) const;
    gen evalf(int level,const context * contextptr) const;
    // inline gen evalf() const { return evalf(DEFAULT_EVAL_LEVEL,context0); }
    gen evalf_double(int level,const context * contextptr) const ;
    gen evalf2double(int level,const context * contextptr) const;
    gen & operator = (const gen & a);
    int to_int() const ;
    bool is_real(GIAC_CONTEXT) const ;
    bool is_cinteger() const ;
    bool is_integer() const ;
    bool is_constant() const;
    std::string print(GIAC_CONTEXT) const;
    std::string print_universal(GIAC_CONTEXT) const;
    std::string print() const;
    // void modify(char * & s) { *this =gen(std::string(s)); };
    void modify(int i) { *this =gen(i); };
    void dbgprint() const; 
    void uncoerce() ;
    gen conj(GIAC_CONTEXT) const;
    gen re(GIAC_CONTEXT) const ;
    gen im(GIAC_CONTEXT) const ;
    gen inverse(GIAC_CONTEXT) const;
    gen squarenorm(GIAC_CONTEXT) const;
    int bindigits() const ;
    gen operator [] (int i) const ;
    gen operator [] (const gen & i) const;
    gen operator_at(int i,GIAC_CONTEXT) const;
    gen operator_at(const gen & i,GIAC_CONTEXT) const;
    // gen & operator [] (int i) ;
    // gen & operator [] (const gen & i) ;
    gen operator () (const gen & i,GIAC_CONTEXT) const;
    bool islesscomplexthan(const gen & other) const;
    bool is_approx() const ; // true if double/real or cmplx with re/im
    int symb_size() const;
    gen change_subtype(int newsubtype);
    bool is_symb_of_sommet(const unary_function_ptr & u) const ;
    gen makegen(int i) const; // make a gen of same type as this with integer i
  };

  bool poly_is_real(const polynome & p);
  bool islesscomplexthanf(const gen & a,const gen & b);
  gen makemap(); // make a new map
  gen chartab2gen(char * & s,GIAC_CONTEXT);
  extern gen zero;
  extern gen plus_one;
  extern gen plus_two;
  extern gen plus_three;
  extern gen minus_one;
  extern gen cst_i;
  extern const vecteur null_vetor;
  extern double rad2deg_d;
  extern double deg2rad_d;
  extern gen rad2deg_g;
  extern gen deg2rad_g;

  bool is_zero(const gen & a);
  bool is_exactly_zero(const gen & a);
  bool is_one(const gen & a);
  bool is_minus_one(const gen & a);
  bool is_sq_minus_one(const gen & a);
  bool is_inf(const gen & e);
  bool is_undef(const gen & e);
  bool is_zero__VECT(const vecteur & a);
  bool has_denominator(const gen & n);
  bool has_i(const gen & g);

  // basic arithmetic
  gen operator && (const gen & a,const gen & b);
  gen operator || (const gen & a,const gen & b);
  gen operator_plus (const gen & a,const gen & b,GIAC_CONTEXT);
  gen operator + (const gen & a,const gen & b);
  gen operator_plus_eq (gen & a,const gen & b,GIAC_CONTEXT);
  inline gen operator += (gen & a,const gen & b){ 
    return operator_plus_eq(a,b,giac::context0);
  }
  Tfraction<gen> operator + (const Tfraction<gen> & a,const Tfraction<gen> & b); // specialization
  gen sym_add (const gen & a,const gen & b,GIAC_CONTEXT);
  gen operator_minus_eq (gen & a,const gen & b,GIAC_CONTEXT);
  inline gen operator -= (gen & a,const gen & b){ 
    return operator_minus_eq(a,b,giac::context0);
  }
  gen operator_minus (const gen & a,const gen & b,GIAC_CONTEXT);
  gen operator - (const gen & a,const gen & b);
  gen operator - (const gen & a);
  gen sym_sub (const gen & a,const gen & b,GIAC_CONTEXT);
  gen operator_times (const gen & a,const gen & b,GIAC_CONTEXT);
  gen operator * (const gen & a,const gen & b);
  gen sym_mult (const gen & a,const gen & b,GIAC_CONTEXT);
  gen pow(const gen & base,const gen & exponent,GIAC_CONTEXT);
  gen iquo(const gen & a,const gen & b); // same
  gen irem(const gen & a,const gen & b,gen & q); // same
  gen smod(const gen & a,const gen & b); // same
  vecteur smod(const vecteur & a,const gen & b); // same
  gen rdiv(const gen & a,const gen & b); // rational division
  inline gen operator /(const gen & a,const gen & b){ return rdiv(a,b); };
  gen operator %(const gen & a,const gen & b); // for int only
  // gen inv(const gen & a);
  gen inv(const gen & a,GIAC_CONTEXT);

  gen algebraic_EXTension(const gen & a,const gen & v);
  gen ext_reduce(const gen & a, const gen & v);
  gen maptoarray(const gen_map & m,GIAC_CONTEXT);
  gen evalf_VECT(const vecteur & v,int subtype,int level,const context * contextptr);
  gen m_gamma(int nbits); // Euler gamma constant precision nbits
  gen m_pi(int nbits); // pi precision nbits

  // a*b -> tmp, modifies tmp in place
  void type_operator_times(const gen & a,const gen &b,gen & tmp);
  inline void type_operator_plus_times(const gen & a,const gen & b,gen & c){
    gen g;
    type_operator_times(a,b,g);
    c += g;
  }

  inline void type_operator_plus_times_reduce(const gen & a,const gen & b,gen & c,int reduce){
    gen g;
    type_operator_times(a,b,g);
    c += g;
    if (reduce)
      c=smod(c,reduce);
  }

  inline void type_operator_reduce(const gen & a,const gen & b,gen & c,int reduce){
    type_operator_times(a,b,c);
    if (reduce)
      c=smod(c,reduce);
  }

  bool operator ==(const gen & a,const gen & b);
  bool operator_equal(const gen & a,const gen & b,GIAC_CONTEXT);
  bool operator !=(const gen & a,const gen & b);
  gen equal(const gen & a,const gen &b);

  gen operator !(const gen & a);

  int fastsign(const gen & a,GIAC_CONTEXT);   // 0 if unknown, 1 if >0, -1 if <0
  gen sign(const gen & a,GIAC_CONTEXT);

  // Large tests if strictly not precised, if sign is unknown return false 
  bool is_greater(const gen & a,const gen &b,GIAC_CONTEXT);
  bool is_strictly_greater(const gen & a,const gen &b,GIAC_CONTEXT);
  bool is_positive(const gen & a,GIAC_CONTEXT);
  bool is_strictly_positive(const gen & a,GIAC_CONTEXT);
  // Large tests if strictly not precised, if sign is unknown make an error
  bool ck_is_greater(const gen & a,const gen &b,GIAC_CONTEXT);
  bool ck_is_strictly_greater(const gen & a,const gen &b,GIAC_CONTEXT);
  bool ck_is_positive(const gen & a,GIAC_CONTEXT);
  bool ck_is_strictly_positive(const gen & a,GIAC_CONTEXT);

  gen superieur_strict(const gen & a,const gen & b,GIAC_CONTEXT);
  gen superieur_egal(const gen & a,const gen & b,GIAC_CONTEXT);
  gen inferieur_strict(const gen & a,const gen & b,GIAC_CONTEXT);
  gen inferieur_egal(const gen & a,const gen & b,GIAC_CONTEXT);
  bool symb_size_less(const gen & a,const gen & b);

  gen min(const gen & a, const gen & b,GIAC_CONTEXT);
  gen max(const gen & a, const gen & b,GIAC_CONTEXT=context0);
  // default context0 is required for instantiation in poly.h
  gen factorial(unsigned long int i);
  gen comb(unsigned long int i,unsigned long j);
  gen perm(unsigned long int i,unsigned long j);
  gen pow(const gen & base, unsigned long int exponent);
  gen pow(const gen & base, int exponent);
  gen pow(unsigned long int base, unsigned long int exponent);

  // more advanced arithmetic
  gen gcd(const gen & A,const gen & B);
  gen lcm(const gen & a,const gen & b);
  gen simplify(gen & n, gen & d);
  void egcd(const gen &a,const gen &b, gen & u,gen &v,gen &d );
  gen ichinrem(const gen & a,const gen &b,const gen & amod, const gen & bmod);
  gen invmod(const gen & A,const gen & modulo);
  gen fracmod(const gen & a_orig,const gen & modulo); // -> p/q=a mod modulo
  gen powmod(const gen &base,const gen & expo,const gen & modulo);
  gen isqrt(const gen & A);
  gen re(const gen & a,GIAC_CONTEXT);
  gen no_context_re(const gen & a);
  gen im(const gen & a,GIAC_CONTEXT);
  gen no_context_im(const gen & a);
  void reim(const gen & g,gen & r,gen & i,GIAC_CONTEXT);
  gen conj(const gen & a,GIAC_CONTEXT);
  gen no_context_conj(const gen & a);
  gen sq(const gen & a);
  gen abs(const gen & a,const context * contextptr=context0);
  // default context0 is required for instantiation in poly.h
  gen linfnorm(const gen & a,const context * contextptr=context0);
  // default context0 is required for instantiation in poly.h
  gen arg(const gen & a,GIAC_CONTEXT);
  gen arg_CPLX(const gen & a,GIAC_CONTEXT);
  int is_perfect_square(const gen & A);
  int is_probab_prime_p(const gen & A);
  gen nextprime(const gen & a); // more precisely next probably prime
  gen prevprime(const gen & a); // more precisely prev probably prime
  int jacobi(const gen & A, const gen &B);
  int legendre(const gen & A, const gen & B);
  vecteur pascal_next_line(const vecteur & v); 
  vecteur pascal_nth_line(int n);
  // convert a __VECTOR__VECT vecteur to a normal vecteur
  gen vector2vecteur(const vecteur & v);

  // if b is a _MOD, returns a as a b _MOD 
  gen chkmod(const gen& a,const gen & b);
  // make a _MOD a%b
  gen makemod(const gen & a,const gen & b);
  // same without evaluating %
  gen makemodquoted(const gen & a,const gen & b);

  // from a sum in x returns a list of [coeff monomial]
  // e.g. 5+2x+3*x*y -> [ [5 1] [2 x] [ 3 x*y] ]
  vecteur symbolique2liste(const gen & x,GIAC_CONTEXT);
  // v should be sorted and shrinked
  gen liste2symbolique(const vecteur & v);

  bool is_atomic(const gen & e);
  symbolic _FRAC2_SYMB(const fraction & f);
  symbolic _FRAC2_SYMB(const gen & e);
  symbolic _FRAC2_SYMB(const gen & n,const gen & d);
  gen string2gen(const std::string & ss,bool remove_ss_quotes=true);
  // by default ss is assumed to be delimited by " and "
  std::complex<double> gen2complex_d(const gen & e);
  gen eval_VECT(const vecteur & v,int subtype,int level,const context * context_ptr );
  // functional equivalent of gen methods
  inline gen eval(const gen & e,int level,const context * contextptr){ return e.eval(level,contextptr); };
  inline gen eval(const gen & e,const context * contextptr){ return e.eval(eval_level(contextptr),contextptr); };
  gen no_context_evalf(const gen & e);
  gen evalf(const gen & e,int level,const context * contextptr );
  inline gen evalf_double(const gen & e,int level,const context * contextptr){ return e.evalf_double(level,contextptr); };
  // return true if g can be converted to a double or real or complex
  bool has_evalf(const gen & g,gen & res,int level,const context * contextptr);
  inline std::string print(const gen & e,context * contextptr){ return e.print(contextptr); }
  inline bool is_real(const gen & g,GIAC_CONTEXT){ return g.is_real(contextptr); }
  inline  bool is_cinteger(const gen & g){ return g.is_cinteger();}  ;
  inline  bool is_integer(const gen & g){ return g.is_integer(); }  ;
  inline  bool is_constant(const gen & g){ return g.is_constant(); } ;
  inline bool is_approx(const gen & g){ return g.is_approx(); };
  gen aplatir_fois_plus(const gen & g);

  class gen_user{
  public:
    // You must redefine memory_alloc, copy_to and memory_free
    // memory_alloc should 
    // return dynamic_cast<gen_user *> new your_class(*this);
    // copy_to should copy your data members taking care of dynamically
    // allocated data
    // memory_free shoudl
    // delete(dynamic_cast<your_class *>(ptr));
    virtual gen_user * memory_alloc() const { return new gen_user(*this); }
    virtual void copy_to (gen_user * g) const { }
    virtual void memory_free(gen_user * ptr){
      delete(ptr);
    }
    // redefine the destructor if you allocate mem dynamically for members
    virtual ~gen_user() {}; 
    gen_user & operator = (const gen_user & g) { 
      g.copy_to(this);
      return *this;
    }
    // redefine operations if it makes sense. 
    // You can redefine gen_user + gen_user for speed
    virtual gen operator + (const gen &) const { setsizeerr("+ not redefined"); return *this; }
    virtual gen operator + (const gen_user & a) const { return (*this) + gen(a); }
    virtual gen operator - (const gen &) const { setsizeerr("Binary - not redefined"); return *this;}
    virtual gen operator - (const gen_user & a) const { return (*this) - gen(a); }
    virtual gen operator - () const { setsizeerr("Unary - not redefined"); return *this;}
    virtual gen operator * (const gen &) const { setsizeerr("Binary * not redefnied"); return *this;}
    virtual gen operator * (const gen_user & a) const { return (*this) * gen(a); }
    virtual gen operator / (const gen_user & a) const { return (*this) * a.inv(); }
    virtual bool is_zero() const { setsizeerr("==0 not redefined"); return false;}
    virtual bool is_one() const { setsizeerr("==1 not redefined"); return false;}
    virtual bool is_minus_one() const { setsizeerr("==-1 not redefined"); return false;}
    virtual gen inv () const { setsizeerr("Inv not redefined"); return *this;}
    virtual gen conj(GIAC_CONTEXT) { setsizeerr("Conj not redefnied"); return *this;}
    virtual gen re(GIAC_CONTEXT) { setsizeerr("Real part not redefnied"); return *this;}
    virtual gen im(GIAC_CONTEXT) { setsizeerr("Imaginary part not redefnied"); return *this;}
    virtual gen abs(GIAC_CONTEXT) { setsizeerr("Abs not redefnied"); return *this;}
    virtual gen arg() { setsizeerr("Arg not redefnied"); return *this;}
    virtual gen operator () (const gen &,GIAC_CONTEXT) const { setsizeerr("() not redefnied"); return *this;}
    virtual gen operator [] (const gen &) { setsizeerr("[] not redefnied"); return *this;}
    virtual bool operator == (const gen &) const { setsizeerr("== not redefnied"); return false;}
    virtual bool operator == (const gen_user & a) const { return (*this) == gen(a); }
    // must redefine > AND <= since we do not have symetrical type arguments
    virtual gen operator > (const gen &) const { setsizeerr("== not redefnied"); return *this;}
    virtual gen operator > (const gen_user & a) const { return superieur_strict(*this, gen(a),0); }
    virtual gen operator <= (const gen &) const { setsizeerr("<= not redefined"); return *this;}
    virtual gen operator <= (const gen_user & a) const { return inferieur_egal(*this, gen(a),0); }
    virtual void polygcd (const polynome &,const polynome &,polynome &) const { setsizeerr("Polynomial gcd not redefined"); }    
    virtual void polyfactor (const polynome & p,
			     factorization & f) const { 
      setsizeerr("Polynomial gcd not redefined"); 
    }    
    virtual gen gcd (const gen &) const { setsizeerr("gcd not redefined"); return *this;}    
    virtual gen gcd (const gen_user & a) const { return gcd(gen(a)); }
    virtual std::string print (GIAC_CONTEXT) const { return  "Nothing_to_print";}
    void dbgprint () const { std::cerr << this->print(0) << std::endl;}
    virtual std::string texprint (GIAC_CONTEXT) const { return "Nothing_to_print_tex"; }
    virtual gen eval(int level,const context * contextptr) const {return *this;};
    virtual gen evalf(int level,const context * contextptr) const {return *this;};
    virtual gen makegen(int i) const { return string2gen("makegen not redefined"); } ;
  };

  std::string print_the_type(int val,GIAC_CONTEXT);

  // I/O
  std::ostream & operator << (std::ostream & os,const gen & a);
  std::istream & operator >> (std::istream & is,gen & a);

  struct monome {
    gen coeff;
    gen exponent;
    monome(const gen & mycoeff) : coeff(mycoeff),exponent(zero) {};
    monome(const gen &mycoeff,const gen &myexponent) : coeff(mycoeff),exponent(myexponent) {};
    std::string print() const ;
    void dbgprint() const ;
  };
  std::ostream & operator << (std::ostream & os,const monome & m);
  inline bool operator == (const monome & a,const monome & b){ return a.coeff==b.coeff && a.exponent==b.exponent; }
  
  std::string printi(GIAC_CONTEXT);
  std::string hexa_print_ZINT(const mpz_t & a);
  std::string octal_print_ZINT(const mpz_t & a);
  std::string binary_print_ZINT(const mpz_t & a);
  std::string printinner_VECT(const vecteur & v, int subtype,GIAC_CONTEXT);
  std::string begin_VECT_string(int subtype,bool tex,GIAC_CONTEXT);
  std::string end_VECT_string(int subtype,bool tex,GIAC_CONTEXT);
  std::string print_VECT(const vecteur & v,int subtype,GIAC_CONTEXT); // subtype was 0 by default
  std::string print_SPOL1(const sparse_poly1 & p,GIAC_CONTEXT);
  // find closing or opening () [] {}
  bool matchpos(const std::string & s,int & pos);
  std::string cut_string(const std::string & chaine,int nchar,std::vector<int> & ligne_end) ;
  std::string calc_endlines_positions(const vecteur & history_in,const vecteur & history_out,int nchar,std::vector<int> & endlines,std::vector<int> & positions);
  bool is_operator_char(char c);
  void increase_selection(const std::string & s,int & pos1,int& pos2);
  void decrease_selection(const std::string & s,int & pos1,int& pos2);
  void move_selection_right(const std::string & s,int & pos1, int & pos2);
  void move_selection_left(const std::string & s,int & pos1, int & pos2);
  std::string remove_extension(const std::string & chaine);


  // This type collects global variables to enable threading
  struct environment {
    gen modulo; // characteristic
    bool moduloon; // Set to false if non modular arithmetic required
    bool complexe; // true if working on Z/pZ[i]
    gen pn; // cardinal of the field, 0 means equal to modulo
    gen coeff; // exemple of coeff, so that we can call coeff.makegen
    environment(){
      modulo=13;
      moduloon=false;
      complexe=false;
      coeff=pn=0;
    }
  };

  // extern environment * env; 

  struct attributs {
    int fontsize;
    int background;
    int text_color;
    attributs(int f,int b,int t): fontsize(f),background(b),text_color(t) {};
    attributs():fontsize(0),background(0),text_color(0) {};
  };

  // Terminal data for EQW display
  struct eqwdata {
    gen g; 
    attributs eqw_attributs;
    int x,y,dx,dy;
    bool selected;
    bool active;
    bool hasbaseline;
    bool modifiable;
    int baseline;
    eqwdata(int dxx,int dyy,int xx, int yy,const attributs & a,const gen& gg):g(gg),eqw_attributs(a),x(xx),y(yy),dx(dxx),dy(dyy),selected(false),active(false),hasbaseline(false),modifiable(true),baseline(0) {};
    eqwdata(int dxx,int dyy,int xx, int yy,const attributs & a,const gen& gg,int mybaseline):g(gg),eqw_attributs(a),x(xx),y(yy),dx(dxx),dy(dyy),selected(false),active(false),hasbaseline(true),modifiable(true),baseline(mybaseline) {};
    void dbgprint(){ std::cerr << g << ":" << dx<< ","<< dy<< "+"<<x <<","<< y<< "," << baseline << "," << eqw_attributs.fontsize << "," << eqw_attributs.background << "," << eqw_attributs.text_color << std::endl; }
  };

  // Graphic object
  struct grob {
    void (* grob_draw)(void);
    int (* grob_handle) (int);
    void * grob_data;
  };

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC

#endif // _GIAC_GEN_H
