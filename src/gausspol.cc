/* -*- mode:C++ ; compile-command: "g++-3.4 -I.. -I../include -g -c gausspol.cc -D_I386_" -*- */
#include "first.h"
/*
 *  This file implements several functions that work on univariate and
 *  multivariate polynomials and rational functions.
 *  These functions include polynomial quotient and remainder, GCD and LCM
 *  computation, factorization and rational function normalization. */

/*
 *  Copyright (C) 2000,2007 B. Parisse, Institut Fourier, 38402 St Martin d'Heres
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
#include "gausspol.h"
#include "modpoly.h"
#include "modfactor.h"
#include "solve.h" // for has_num_coeff
#include "alg_ext.h"
#include "sym2poly.h"
#include "prog.h"
#include "plot.h"
#include "modpoly.h"
#include "threaded.h"
#include <math.h>
#include <stdexcept>

#ifdef HAVE_GMPXX_H
#define myint mpz_class
#else
#define myint my_mpz
#endif

// #undef HAVE_LIBNTL

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  int primes[]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541};

  // trivial division factorization algorithm
  vector<nfactor> trivial_n_factor(gen &n){
    vector<nfactor> v;
    for (int i=0;i<100;i++){
      gen p(primes[i]);
      if (is_zero(n % p) ){
	int j=1;
	n=iquo(n,p);
	while (is_zero(n % p)){
	  n=iquo(n,p);
	  j++;
	}
	v.push_back(nfactor(p,j));
      }
      if (is_strictly_greater(p*p,n,context0))
	break;
    }
    if (n!=gen(1))
      v.push_back(nfactor(n,1));
    return v;
  }

  vecteur cyclotomic(int n){
    // Algorithm base sur les relations suivantes
    // soit p premier, si p ne divise pas n alors 
    // cyclo_{n}(X^p)=cyclo_{np}(X)*cyclo_{n}(X)
    // si p divise n alors
    // cyclo_{n}(X^p)=cyclo_{np}(X)
    //
    // poser f(x)=x-1
    // parcourir la liste des facteurs premiers de n ordonnee par ordre 
    // croissant et pour chaque facteur premier p effecter f(x)=f(x^p)/f(x).
    // Le polynome cyclotomique P_n est alors f(x^[n/pi]) ou pi est le produit
    // des facteurs premiers.
    // Par exemple
    // n=6 donne la liste { 2 3 }
    // p=2: f(x)=(x^2-1)/(x-1)=x+1
    // p=3: f(x)=(x^3+1)/(x+1)=x^2-x+1
    // n/pi=1 c'est fini
    //
    // n=100 donne { 2 5 }
    // p=2: f(x)=x+1
    // p=5: f(x)=(x^5+1)/(x+1)=x^4-x^3+x^2-x+1
    // n/pi=10 donc f(x)=x^40-x^30+x^20-x^10+1
    gen ncopy(n);
    vector<nfactor> v(trivial_n_factor(ncopy));
    vecteur res;
    res.push_back(1);
    res.push_back(-1); // res=x-1
    int pi=1;
    vector<nfactor>::const_iterator it=v.begin(),itend=v.end();
    for (;it!=itend;++it){
      if (it->fact.type!=_INT_) setsizeerr("gausspol.cc/cyclotomic");
      int p=it->fact.val;
      pi *= p;
      vecteur res_x_to_xp(x_to_xp(res,p));
      res=res_x_to_xp/res;
    }
    return x_to_xp(res,n/pi);
  }
  
  //**********************************
  // functions relative to polynomials
  //**********************************
  polynome gen2polynome(const gen & e,int dim){
    if (e.type==_POLY)
      return *e._POLYptr;
    return polynome(e,dim);
  }

  // instantiation of dbgprint for poly
  void dbg(const polynome & p){
    p.dbgprint();
  }

  bool is_one(const polynome & p){
    return Tis_one<gen>(p);
  }

  polynome firstcoeff(const polynome & p){
    return Tfirstcoeff<gen>(p);
  }

  polynome operator + (const polynome & th,const polynome & other) {
    control_c();
    // Tensor addition
    vector< monomial<gen> >::const_iterator a=th.coord.begin();
    vector< monomial<gen> >::const_iterator a_end=th.coord.end();
    if (a == a_end) {
      return other;
    }
    vector< monomial<gen> >::const_iterator b=other.coord.begin();
    vector< monomial<gen> >::const_iterator b_end=other.coord.end();
    if (b==b_end){
      return th;
    }
    polynome res(th.dim,th);
    Add<gen>(a,a_end,b,b_end,res.coord,th.is_strictly_greater);
    return res;
  }

  polynome operator - (const polynome & th,const polynome & other) {  
    control_c();
    // Tensor addition
    vector< monomial<gen> >::const_iterator a=th.coord.begin();
    vector< monomial<gen> >::const_iterator a_end=th.coord.end();
    vector< monomial<gen> >::const_iterator b=other.coord.begin();
    vector< monomial<gen> >::const_iterator b_end=other.coord.end();
    if (b == b_end) {
      return th;
    }
    polynome res(th.dim,th);
    Sub<gen>(a,a_end,b,b_end,res.coord,th.is_strictly_greater);
    return res;
  }

  polynome operator * (const gen & fact, const polynome & th){
    control_c();
    return th*fact;
  }

  void mulpoly(const polynome & th,const gen & fact0,polynome & res){
    if (&th!=&res)
      res.coord.clear();
    gen fact=fact0;
    if (fact.type!=_MOD && !th.coord.empty() && th.coord.front().value.type==_MOD){
      fact = makemod(fact,*(th.coord.front().value._MODptr+1));
    }
    if (!is_zero(fact)){
      vector< monomial<gen> >::const_iterator a = th.coord.begin();
      vector< monomial<gen> >::const_iterator a_end = th.coord.end();
      Mul<gen>(a,a_end,fact,res.coord);
    }
  }

  polynome operator * (const polynome & th, const gen & fact){
    control_c();
    // Tensor constant multiplication
    if (fact.type!=_MOD && fact==gen(1))
      return th;
    polynome res(th.dim,th);
    mulpoly(th,fact,res);
    return res;
  }
  
  ostream & operator << (ostream & os,const int_unsigned & i){
    return os << i.g << ":" << i.u ;
  }

  inline bool operator < (const int_unsigned & gu1,const int_unsigned & gu2){
    return gu1.u > gu2.u;
  }

  template <class U>
  bool convert(const polynome & p,const index_t & deg,std::vector< T_unsigned<int,U> >  & v,int reduce){
    std::vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    v.clear();
    v.reserve(itend-it);
    T_unsigned<int,U> gu;
    U u;
    index_t::const_iterator itit,ditbeg=deg.begin(),ditend=deg.end(),dit;
    gen tmp;
    for (;it!=itend;++it){
      u=0;
      itit=it->index.iptr->begin();
      for (dit=ditbeg;dit!=ditend;++itit,++dit)
	u=u*U(*dit)+U(*itit);
      gu.u=u;
      tmp=smod(it->value,reduce);
      if (tmp.type!=_INT_)
	return false;
      gu.g=tmp.val;
      v.push_back(gu);
    }
    return true;
  }

  template<class U>
  void convert(const std::vector< T_unsigned<int,U> > & v,const index_t & deg,polynome & p){
    typename std::vector< T_unsigned<int,U> >::const_iterator it=v.begin(),itend=v.end();
    index_t::const_reverse_iterator ditbeg=deg.rbegin(),ditend=deg.rend(),dit;
    p.dim=ditend-ditbeg;
    p.coord.clear();
    p.coord.reserve(itend-it);
    U u;
    index_t i(p.dim);
    int k;
    for (;it!=itend;++it){
      u=it->u;
      for (k=p.dim-1,dit=ditbeg;dit!=ditend;++dit,--k){
	i[k]=u % unsigned(*dit);
	u = u/unsigned(*dit);
      }
      p.coord.push_back(monomial<gen>(it->g,i));
    }
  }

  void smallmult(const std::vector< int_unsigned > & v1,const std::vector< int_unsigned > & v2,std::vector< int_unsigned > & v,int reduce,int possible_size=100){
#ifdef HASH_MAP_NAMESPACE
    typedef HASH_MAP_NAMESPACE::hash_map<unsigned,int> hash_prod ;
    hash_prod produit(possible_size);
    // cout << "hash " << clock() << endl;
#else
    typedef std::map<unsigned,int> hash_prod;
    hash_prod produit;
    // cout << "small map" << endl;
#endif    
    hash_prod::iterator prod_it,prod_itend;
    std::vector< int_unsigned >::const_iterator it1=v1.begin(),it1end=v1.end(),it2beg=v2.begin(),it2,it2end=v2.end();
    // FIXME if reduce is small use int for g1,g instead of longlong
    longlong g1,g;
    unsigned u1,u;
    for (;it1!=it1end;++it1){
      g1=it1->g;
      u1=it1->u;
      for (it2=it2beg;it2!=it2end;++it2){
	u=u1+it2->u;
	g=g1*it2->g ; // moved % reduce so that 1 % is done instead of 2
	prod_it=produit.find(u);
	if (prod_it==produit.end())
	  produit[u]=g % reduce;
	else {
	  int & s=prod_it->second;
	  g += s;
	  s = g % reduce;
	}
      }
    }
    int_unsigned gu;
    prod_it=produit.begin(),prod_itend=produit.end();
    v.clear();
    v.reserve(produit.size());
    for (;prod_it!=prod_itend;++prod_it){
      if (!is_zero(gu.g=prod_it->second)){
	gu.u=prod_it->first;
	v.push_back(gu);
      }
    }    
    // cout << "smallmult end " << clock() << endl;
    sort(v.begin(),v.end());
  }

  void smallmult(int x,std::vector<int_unsigned> & v,int m){
    if (!x){
      v.clear();
      return;
    }
    std::vector<int_unsigned>::iterator it=v.begin(),itend=v.end();
    for (;it!=itend;++it){
      it->g *= x;
      it->g %= m;
    }
  }

  void smalladd(const std::vector< int_unsigned > & v1,const std::vector< int_unsigned > & v2,int m,std::vector< int_unsigned > & v){
    std::vector< int_unsigned >::const_iterator it1=v1.begin(),it1end=v1.end(),it2=v2.begin(),it2end=v2.end();
    int g;
    v.clear();
    v.reserve(it1end-it1+it2end-it2); // worst case
    for (;it1!=it1end && it2!=it2end;){
      if (it1->u==it2->u){
	g=(it1->g+it2->g)%m;
	if (g)
	  v.push_back(int_unsigned(g,it1->u));
	++it1;
	++it2;
      }
      else {
	if (it1->u>it2->u){
	  v.push_back(*it1);
	  ++it1;
	}
	else {
	  v.push_back(*it2);
	  ++it2;
	}
      }
    }
    for (;it1!=it1end;++it1)
      v.push_back(*it1);
    for (;it2!=it2end;++it2)
      v.push_back(*it2);
  }

  // Poly evaluation of p at x modulo m, d is the degree in int_unsigned.u
  void peval(const vector<int_unsigned> & p,int d,int x,int m,vector<int_unsigned> & res){
    res.clear();
    vector<int_unsigned> tmp1,tmp2;
    if (p.empty())
      return;
    // cerr << p << endl;
    vector<int_unsigned>::const_iterator it=p.begin(),itend=p.end();
    int deg=d*(it->u / d),ddeg;
    for (;deg>=0;deg -=d ){ // Horner like
      // cerr << res << endl;
      smallmult(x,res,m);
      // cerr << res << endl;
      tmp2.clear();
      // Find next coeff
      for (;it!=itend;++it){
	ddeg=it->u-deg;
	if (ddeg<0)
	  break;
	tmp2.push_back(int_unsigned(it->g,ddeg));
      }
      // cerr << tmp2 << endl;
      tmp1=res;
      smalladd(tmp1,tmp2,m,res);
      // cerr << res << endl;
    }
  }

  bool degree2unsigned(index_t & deg,unsigned & u){
    u=1;
    index_t::iterator it=deg.begin(),itend=deg.end();
    for (;it!=itend;++it){
      ++(*it);
      u = u*unsigned(*it);
      if (u>RAND_MAX)
	return false;
    }
    return true;
  }

  bool peval(const polynome & p,const gen & x0,int m,polynome & g,vector<int_unsigned> * P){
    gen x1=smod(x0,m);
    if (x1.type!=_INT_)
      return false;
    int x=x1.val;
    index_t d=p.degree();
    unsigned ans;
    if (!degree2unsigned(d,ans))
      return false;
    vector<int_unsigned> Q;
    if (!P) 
      P=&Q;
    if (P->empty() && !convert(p,d,*P,m))
      return false;
    vector<int_unsigned> res;
    peval(*P,ans/d.front(),x,m,res);
    d.erase(d.begin());
    convert(res,d,g);
    return true;
  }



  void convert(const vector< T_unsigned<longlong,unsigned> > & source,vector< T_unsigned<gen,unsigned> > & target){
    target.clear();
    vector< T_unsigned<longlong,unsigned> >::const_iterator it=source.begin(),itend=source.end();
    target.reserve(itend-it);
    for (;it!=itend;++it)
      target.push_back(T_unsigned<gen,unsigned>(it->g,it->u));
  }

  gen ichrem_smod(mpz_t * Az,mpz_t * Bz,mpz_t * iz,mpz_t * tmpz,const gen & i,const gen & j){
    if (i.type==_ZINT)
      mpz_set(*iz,*i._ZINTptr);
    else
      mpz_set_si(*iz,i.val);
    // i-j
    if (j.type==_INT_){
      if (j.val>0)
	mpz_sub_ui(*tmpz,*iz,j.val);
      else
	mpz_add_ui(*tmpz,*iz,-j.val);
    }
    else
      mpz_sub(*tmpz,*iz,*j._ZINTptr);
    // times B +i
    mpz_addmul(*iz,*tmpz,*Bz);
    // mod A
    mpz_mod(*tmpz,*iz,*Az);
    // compare with *tmpz-Az
    mpz_sub(*iz,*tmpz,*Az);
    mpz_neg(*iz,*iz);
    mpz_t *res=(mpz_t *) malloc(sizeof(mpz_t));
    if (mpz_cmp(*iz,*tmpz)>=0) // use *tmpz
      mpz_init_set(*res,*tmpz);
    else {
      mpz_init_set(*res,*iz);
      mpz_neg(*res,*res);
    }
    return res;
  }

  gen ichrem_smod(mpz_t * Az,mpz_t * Bz,mpz_t * iz,mpz_t * tmpz,longlong i,longlong j){
    if (i==j)
      return i;
    longlong2mpz(i,iz);
    // longlong2mpz(i-j,tmpz); does not work since i-j might overflow
    longlong2mpz(j,tmpz);
    mpz_sub(*tmpz,*iz,*tmpz);
    // i+=B*(i-j)
    mpz_addmul(*iz,*tmpz,*Bz);
    // mod A
    mpz_mod(*tmpz,*iz,*Az);
    // compare with *tmpz-Az
    mpz_sub(*iz,*tmpz,*Az);
    mpz_neg(*iz,*iz);
    mpz_t *res=(mpz_t *) malloc(sizeof(mpz_t));
    int test=mpz_cmp(*iz,*tmpz);
    if (test>=0) // use *tmpz
      mpz_init_set(*res,*tmpz);
    else {
      mpz_init_set(*res,*iz);
      mpz_neg(*res,*res);
    }
    return res;
  }

  // smod(B*(i-j)+i,A);  
  gen ichrem_smod(const gen & A,const gen & B,const gen & i,const gen & j){
    if (i==j)
      return i;
    if (A.type!=_ZINT || B.type!=_ZINT)
      return smod(B*(i-j)+i,A);  
    mpz_t * Az=A._ZINTptr,*Bz=B._ZINTptr,iz,tmpz,*res=(mpz_t *) malloc(sizeof(mpz_t));
    mpz_init(tmpz);
    if (i.type==_ZINT)
      mpz_init_set(iz,*i._ZINTptr);
    else
      mpz_init_set_si(iz,i.val);
    // i-j
    if (j.type==_INT_){
      if (j.val>0)
	mpz_sub_ui(tmpz,iz,j.val);
      else
	mpz_add_ui(tmpz,iz,-j.val);
    }
    else
      mpz_sub(tmpz,iz,*j._ZINTptr);
    // times B +i
    mpz_addmul(iz,tmpz,*Bz);
    // mod A
    mpz_mod(tmpz,iz,*Az);
    // compare with tmpz-Az
    mpz_sub(iz,tmpz,*Az);
    mpz_neg(iz,iz);
    if (mpz_cmp(iz,tmpz)>0) // use tmpz
      mpz_init_set(*res,tmpz);
    else {
      mpz_init_set(*res,iz);
      mpz_neg(*res,*res);
    }
    mpz_clear(iz);
    mpz_clear(tmpz);
    return res;
  }

  gen ichrem_smod(const gen & A,const gen & B,longlong i,longlong j){
    // return smod((i-j)*B+i,A);  
    if (i==j)
      return i;
    return ichrem_smod(A,B,gen(i),gen(j));
  }

  template<class T,class U>
  void ichrem(const vector< T_unsigned<T,U> > & add,int addprime,const vector< T_unsigned<longlong,U> > & init,vector< T_unsigned<gen,U> > & target,gen & targetprime){
    gen A,B,d;
    egcd(addprime,targetprime,A,B,d);
    if (!is_one(d))
      setsizeerr();
    // addprime*A+targetprime*B=1
    // find c such that c=it->g mod targetprime and c=jt->g mod addprime
    // it->g + v*targetprime = jt->g + u*addprime
    // it->g - jt->g = u*addprime - v*targetprime
    // v=(jt->g-it->g)*B
    // hence c=it->g+(jt->g-it->g)*B*targetprime mod addprime*targetprime
    A=addprime*targetprime;
    B=-targetprime*B;
    mpz_t z1,z2;
    mpz_init(z1);
    mpz_init(z2);
    if (A.type!=_ZINT)
      A.uncoerce();
    if (B.type!=_ZINT)
      B.uncoerce();
    typename vector< T_unsigned<T,U> >::const_iterator jt=add.begin(),jtend=add.end();
    typename vector< T_unsigned<longlong,U> >::const_iterator kt=init.begin(),ktend=init.end();
    typename vector< T_unsigned<gen,U> >::iterator it=target.begin(),itend=target.end();
    if (it==itend){
      target.reserve(ktend-kt);
      for (;kt!=ktend && jt!=jtend;){
	if (kt->u==jt->u){
	  // it->g=smod(B*(it->g-jt->g)+it->g,A);
	  target.push_back(T_unsigned<gen,U>(ichrem_smod(A._ZINTptr,B._ZINTptr,&z1,&z2,kt->g,jt->g),jt->u));
	  ++kt; ++jt;
	}
	else {
	  if (kt->u>jt->u){
	    target.push_back(T_unsigned<gen,U>(ichrem_smod(A._ZINTptr,B._ZINTptr,&z1,&z2,kt->g,0),kt->u));
	    ++kt;
	  }
	  else {
	    target.push_back(T_unsigned<gen,U>(ichrem_smod(A._ZINTptr,B._ZINTptr,&z1,&z2,0,jt->g),jt->u));
	    ++jt;
	  }
	}
      }
      for (;jt!=jtend;++jt)
	target.push_back(T_unsigned<gen,U>(ichrem_smod(A._ZINTptr,B._ZINTptr,&z1,&z2,0,jt->g),jt->u));
      for (;kt!=ktend;++jt)
	target.push_back(T_unsigned<gen,U>(ichrem_smod(A._ZINTptr,B._ZINTptr,&z1,&z2,kt->g,0),kt->u));
    }
    else {
      for (;it!=itend && jt!=jtend;){
	if (it->u==jt->u){
	  // it->g=smod(B*(it->g-jt->g)+it->g,A);
	  it->g=ichrem_smod(A._ZINTptr,B._ZINTptr,&z1,&z2,it->g,jt->g);
	  ++it; ++jt;
	}
	else {
	  if (it->u>jt->u){
	    it->g=ichrem_smod(A._ZINTptr,B._ZINTptr,&z1,&z2,it->g,0);
	    ++it;
	  }
	  else {
	    vector< T_unsigned<gen,U> > copie(it,itend);
	    target.erase(it,itend);
	    it=copie.begin(); itend=copie.end();
	    for (;it!=itend;){
	      if (jt==jtend || it->u>jt->u){
		it->g=ichrem_smod(A._ZINTptr,B._ZINTptr,&z1,&z2,it->g,0);
		target.push_back(*it);
		++it;
	      }
	      else {
		++jt;
		if (it->u==jt->u){
		  target.push_back(T_unsigned<gen,U>(ichrem_smod(A._ZINTptr,B._ZINTptr,&z1,&z2,it->g,jt->g),it->u));
		  ++it; 
		}
		else 
		  target.push_back(T_unsigned<gen,U>(ichrem_smod(A._ZINTptr,B._ZINTptr,&z1,&z2,0,jt->g),jt->u));
	      }
	    }
	    break;
	  }
	}
      }
      for (;jt!=jtend;++jt)
	target.push_back(T_unsigned<gen,U>(ichrem_smod(A._ZINTptr,B._ZINTptr,&z1,&z2,0,jt->g),jt->u));
      for (;it!=itend;++it)
	it->g=ichrem_smod(A._ZINTptr,B._ZINTptr,&z1,&z2,it->g,0);
    } // end target empty at beginning
    targetprime = addprime*targetprime;
    mpz_clear(z1);
    mpz_clear(z2);
  }

  void smod(const vector< T_unsigned<longlong,unsigned> > & source,vector< T_unsigned<longlong,unsigned> > & target,int prime){
    if (&target==&source){
      vector< T_unsigned<longlong,unsigned> >::iterator it=target.begin(),itend=target.end();
      for (;it!=itend;++it){
	it->g %= prime;
	if (!it->g){
	  vector< T_unsigned<longlong,unsigned> > copie(target);
	  smod(copie,target,prime);
	  break;
	}
      }
      return;
    }
    target.clear();
    vector< T_unsigned<longlong,unsigned> >::const_iterator it=source.begin(),itend=source.end();
    target.reserve(itend-it);
    longlong res;
    for (;it!=itend;++it){
      res=it->g % prime;
      if (res)
	target.push_back(T_unsigned<longlong,unsigned>(res,it->u));
    }
  }

  // Fast multiplication using hash maps, might also use an int for reduction
  // but there is no garantee that res is smod-ed modulo reduce
  void mulpoly (const polynome & th, const polynome & other,polynome & res,const gen & reduce){
    control_c();
    /*
    if (th.dim==12)
      cerr << "* begin " << clock() << " " << th.coord.size() << "*" << other.coord.size() << endl;
    */
    // Multiplication
    vector< monomial<gen> >::const_iterator ita = th.coord.begin();
    vector< monomial<gen> >::const_iterator ita_end = th.coord.end();
    vector< monomial<gen> >::const_iterator itb = other.coord.begin();
    vector< monomial<gen> >::const_iterator itb_end = other.coord.end();
    //  cout << coord.size() << " " << (int) ita_end - (int) ita << " " << sizeof(monomial<gen>) << endl ;
    // first some trivial cases
    if (ita==ita_end || is_one(other)){
      res=th;
      return;
    }
    if (itb==itb_end || is_one(th)){
      res=other;
      return ;
    }
    // Now look if length a=1 or length b=1, happens frequently
    // think of x^3*y^2*z translated to internal form
    int c1=th.coord.size();
    if (c1==1){
      res=other.shift(th.coord.front().index,th.coord.front().value);
      return ;
    }
    int c2=other.coord.size();
    if (c2==1){
      res=th.shift(other.coord.front().index,other.coord.front().value);
      return;
    }
    if ( c1>50 || c2 >50 || (c1>7 && c2>7) ){
      // Degree info, try to multiply the polys using integer for the exponents
      index_t d1=th.degree(),d2=other.degree(),d(th.dim);
      ulonglong ans=1,pid1=1,pid2=1;
      for (int i=0;i<th.dim;++i){
	pid1 = pid1*unsigned(d1[i]+1);
	pid2 = pid2*unsigned(d2[i]+1);
      }
      for (int i=0;i<th.dim;++i){
	d[i]=d1[i]+d2[i]+1;
	ans = ans*unsigned(d[i]);
	if (ans/RAND_MAX>RAND_MAX)
	  break;
      }
      // ans/d[th.dim-1] == degree 1 with respect to main var and 0 for other
      // guess size of result
      // compare product of d1[i] with c1 and product of d2[i] with c2
      // for sparness factor
      double d1sparness=double(c1)/pid1;
      double d2sparness=double(c2)/pid2;
      ulonglong c1c2= ulonglong(c1)*c2;
      if (ans<c1c2)
	c1c2=ans;
      c1c2 = unsigned(std::sqrt(d1sparness*d2sparness)*c1c2);
      if (c1c2> (1<<24) )
	c1c2 = 1 << 24;
      int t1=th.coord.front().value.type,t2=other.coord.front().value.type;
      if (ans<=RAND_MAX){
	if (reduce.type==_INT_ && reduce.val<46340 && reduce.val>0){ 
	  // Modular multiplication, convert everything to integers
	  vector< int_unsigned > p1,p2,p;
	  if (convert(th,d,p1,reduce.val) && convert(other,d,p2,reduce.val)){
	    smallmult(p1,p2,p,reduce.val,c1c2);
	    convert(p,d,res);
	    return ;
	  }
	}
	if ( //false 
	     (t1==_INT_ || t1==_ZINT) && (t2==_INT_ || t2==_ZINT)
	    ){
	  longlong maxp1,maxp2;
	  // should be T_unsigned<long,unsigned>
	  // instead tmp_operator_times converts longlong args of * to long
	  vector< T_unsigned<longlong,unsigned> > p1d,p2d,pd;
	  if (convert_int(th,d,p1d,maxp1) && convert_int(other,d,p2d,maxp2) ){
	    double maxp1p2=double(maxp1)*maxp2;
	    unsigned minc1c2=min(c1,c2);
	    double un63=(ulonglong (1) << 63);
	    double res_size=double(minc1c2)*maxp1p2;
	    // Check number of required primes: 
	    res_size /= un63;
	    double nprimes=std::ceil(std::log(res_size)/std::log(2147483647.));
	    if (
#ifdef HAVE_GMPXX_H
		mpzclass_allowed?nprimes<1.5:nprimes<4.5
#else
		nprimes<4.5
#endif
		){
	      if(debug_infolevel>5) cerr << "Begin smallmult " << clock() << endl;
	      if (th.dim==1 || !threadmult<longlong,unsigned>(p1d,p2d,pd,ans/d[0],0,c1c2))
		smallmult<longlong,unsigned>(p1d,p2d,pd,0,c1c2);
	      if(debug_infolevel>5) cerr << "End smallmult " << clock() << endl;
	      unsigned pds=pd.size();
	      if ( res_size< 1 ){
		convert_from<longlong,unsigned>(pd,d,res);
		return;
	      }
	      // /*
	      else {
		if (debug_infolevel)
		  cerr << nprimes << " primes required" << endl;
		vector< T_unsigned<gen,unsigned> > target;
		// convert(pd,target);
		int prime1=2147483647;
		double primed=3037000499.; // floor(2^31.5)
		primed /= std::sqrt(double(max(minc1c2,4)));
		// the product of the two poly mod prime2 can be computed
		// without mod computation
		// FIXME Currently disabled, seems to bug
		// maybe because I should smod the product by prime2 before ichrem
		int prime2=prevprime(int(std::floor(primed))).val;
		gen targetprime = pow(plus_two,64);
		for(;res_size>=1;--nprimes){
		  bool withsmod
		    //=true;
		    =(res_size>std::pow(prime1-1000,nprimes-1)*prime2);
		  if (debug_infolevel>5) cerr << clock () << " prime " << (withsmod?prime1:prime2) << endl;
		  if (withsmod){
		    vector< int_unsigned > p1,p2,padd;
		    if (!convert(th,d,p1,prime1) || !convert(other,d,p2,prime1))
		      setsizeerr();
		    if (th.dim==1 || !threadmult<int,unsigned>(p1,p2,padd,ans/d[0],prime1,pds))
		      smallmult(p1,p2,padd,prime1,pds);
		    if (debug_infolevel>5) cerr << clock () << " ichrem " << endl;	
		    ichrem(padd,prime1,pd,target,targetprime);
		    res_size /= prime1;
		    prime1=prevprime(prime1-2).val;
		  }
		  else {
		    vector< T_unsigned<longlong,unsigned> > p1add,p2add,add;
		    smod(p1d,p1add,prime2);
		    smod(p2d,p2add,prime2);
		    if (th.dim==1 || !threadmult<longlong,unsigned>(p1add,p2add,add,ans/d[0],0,pds))
		      smallmult(p1add,p2add,add,0,pds);
		    smod(add,add,prime2);
		    if (debug_infolevel>5) cerr << clock () << " ichrem " << endl;	
		    ichrem(add,prime2,pd,target,targetprime);
		    res_size /= prime2;
		    prime2=prevprime(prime2-2).val;
		  }
		  if (debug_infolevel>5) cerr << clock () << " ichrem end" << endl;
		}
		if (debug_infolevel>5) cerr << clock() << endl;
		convert_from<gen,unsigned>(target,d,res);
		if (debug_infolevel>5) cerr << clock() << endl;
		return;
	      } // end else (some primes are required)
	      // */
	    } // end nprimes<something
	  } // end conversion possible
	} // end t1==_INT_
	if (t1==_DOUBLE_ && t2==_DOUBLE_){
	  vector< T_unsigned<double,unsigned> > p1d,p2d,pd;
	  if (convert_double(th,d,p1d) && convert_double(other,d,p2d) ){
	    if (th.dim==1 || !threadmult<double,unsigned>(p1d,p2d,pd,ans/d[0],0,c1c2))
	      smallmult<double,unsigned>(p1d,p2d,pd,0,c1c2);
	    convert_from<double,unsigned>(pd,d,res);
	    return;
	  }
	}
#ifdef HAVE_GMPXX_H
	if (t1<=_ZINT && t2<=_ZINT && mpzclass_allowed){
	  if (debug_infolevel>1)
	    cerr << "mpz mult convert begin " << clock() << endl;
	  vector< T_unsigned<myint,unsigned> > p1d,p2d,pd;
	  if (convert_myint(th,d,p1d) && convert_myint(other,d,p2d) ){
	    if (debug_infolevel>1)
	      cerr << "mpz mult begin " << clock() << endl;
	    if (th.dim==1 || !threadmult<myint,unsigned>(p1d,p2d,pd,ans/d[0],0,c1c2))
	      smallmult<myint,unsigned>(p1d,p2d,pd,0,c1c2);
	    if (debug_infolevel>1)
	      cerr << "mpz mult end " << clock() << endl;
	    convert_from<myint,unsigned>(pd,d,res);
	    return;
	  }
	}
#endif
	vector< T_unsigned<gen,unsigned> > p1,p2,p;
	convert<gen,unsigned>(th,d,p1);
	convert<gen,unsigned>(other,d,p2);
	if (th.dim==1 || !threadmult<gen,unsigned>(p1,p2,p,ans/d[0],0,c1c2))
	  smallmult<gen,unsigned>(p1,p2,p,0,c1c2);
	convert<gen,unsigned>(p,d,res);
	return ;
	// cerr << "Copy " << clock() << " " << copy_number << endl;
	// if (th.dim==12)
	//  cerr << "sort *unsigned end " << clock() << " " << res.coord.size() << endl;
	/*
	  polynome save(res);
	  sort(res.coord.begin(),res.coord.end(),th.m_is_greater); // still done
	  if (res!=save)
	  cerr << "unsorted" << endl;
	*/
	// if (res.coord.size()==1357366)
	//  cerr << "coucou" << endl;
	// if (th.dim==12){
	//  cerr << "*unsigned end " << clock() << endl;
      }
      if (ans/RAND_MAX<RAND_MAX){
	if ( (t1==_INT_ || t1==_ZINT) && (t2==_INT_ || t2==_ZINT)
	    ){
	  longlong maxp1,maxp2;
	  // should be T_unsigned<long,unsigned>
	  // instead tmp_operator_times converts longlong args of * to long
	  vector< T_unsigned<longlong,ulonglong> > p1d,p2d,pd;
	  if (convert_int(th,d,p1d,maxp1) && convert_int(other,d,p2d,maxp2) ){
	    double maxp1p2=double(maxp1)*maxp2;
	    unsigned minc1c2=min(c1,c2);
	    double un63=(ulonglong (1) << 63);
	    double res_size=double(minc1c2)*maxp1p2;
	    // Check number of required primes: 
	    res_size /= un63;
	    if ( res_size< 1 ){
	      if (th.dim==1 || !threadmult<longlong,ulonglong>(p1d,p2d,pd,ans/d[0],0,c1c2))
		smallmult<longlong,ulonglong>(p1d,p2d,pd,0,c1c2);
	      convert_from<longlong,ulonglong>(pd,d,res);
	      return;
	    }
	  }
	}
	vector< T_unsigned<gen,ulonglong> > p1,p2,p;
	convert<gen,ulonglong>(th,d,p1);
	convert<gen,ulonglong>(other,d,p2);
	smallmult<gen,ulonglong>(p1,p2,p,0,c1c2);
	convert<gen,ulonglong>(p,d,res);
	// if (th.dim==12)
	//  cerr << "sort*longlong end " << clock() << " " << res.coord.size() << endl;
	// sort(res.coord.begin(),res.coord.end(),th.m_is_greater); // still done
	// if (th.dim==12)
	//  cerr << "*longlong end " << clock() << endl;
	return ;
      }
    } // end if c1>7 && c2>7
    res.coord.clear();
    Mul<gen>(ita,ita_end,itb,itb_end,res.coord,th.is_strictly_greater,th.m_is_greater);
    // if (th.dim==12)
    //  cerr << "* end " << clock() << " " << res.coord.size() << endl;
    return ;
  }

  polynome operator * (const polynome & th, const polynome & other) {
    polynome res(th.dim,th); // reserve() is done by Mul<gen>
    mulpoly(th,other,res,0);
    return res;
  }

  /* 
    Note about Miller Pure Recurrence, see Knuth, TAOC v.2
    If   P(x) = sum_{i=0}^n p_i x^k
    Then   P(x)^m = sum_{k=0}^{m*n} a(m,k) x^k
    Where
      a(m,0) = p_0^m, 
      a(m,k) = 1/(k p_0) sum_{i=1}^max(n,k) p_i ((m+1)i-k) a(m,k-i),
    For k<=m we have a division free implementation, let
    a(m,k)=b(m,k) p_0^(m-k)
    b(m,0)=1, b(m,k)=1/k sum_{i=1}^max(n,k) p_i ((m+1)i-k) b(m,k-i) p_0^(i-1)
    But for k>m, the division by p0 must be done at each step 
    which might be too costly
    Example: P(x)=3x^2+2x+5, n=2
    m=2: P^2=9*x^4+12*x^3+34*x^2+20*x+25
    b(m,0)=1, a(m,0)=25
    b(m,1)= p_1*(3*1-1)*b(m,0)=4, a(m,1)=20
    b(m,2)=1/2*(p_1*(3*1-2)*b(m,1)+p_2*(3*2-2)*p_0*b(m,0))
          =1/2*(2*4+3*4*5)=34, a(m,2)=34
    a(m,3)=1/3/5*(p_1*(3*1-3)*a(m,2)+p_2*(3*2-3)*a(m,1))=12
    a(m,4)=1/4/5*(p_1*(3*1-4)*a(m,3)+p_2*(3*2-4)*a(m,2))=1/20*(-2*12+3*2*34)=9
    There is a case where no bad division occurs: if p_0 is a constant
    (no other variable occur) or if n==1 (binomial formula)
  */
  void powpoly (const polynome & th, int u,polynome & res){
    if (u<0)
      setsizeerr("Negative polynome power");
    if (!u){
      res= tensor<gen>(1,th.dim);
      return;
    }
    if (u==1){
      res=th;
      return;
    }
    if (u==2){
      res=th*th;
      return;
    }
    control_c();
    vector< monomial<gen> >::const_iterator ita = th.coord.begin();
    vector< monomial<gen> >::const_iterator ita_end = th.coord.end();
    int c1=ita_end-ita;
    if (c1==0){
      res=th;
      return;
    }
    if (c1==1){
      res=th;
      res.coord.front().value=pow(res.coord.front().value,u);
      res.coord.front().index = res.coord.front().index*u ;
      return ;
    }
    ulonglong ans=1,pid1=1;
    index_t d1=th.degree(),d(th.dim);
    for (int i=0;i<th.dim;i++){
      if (d1[i]==1){ // apply binomial formula
	vecteur v,w(u+1);
	polynome2poly1(th,i+1,v);
	gen a=v[0],b=v[1],bk=b;
	vecteur vbk=makevecteur(1,b);
	for (unsigned j=2;j<=u;++j){
	  bk=bk*b;
	  vbk.push_back(bk);
	}
	// (ax+b)^u=sum_j=0^u comb(u,j)*a^j*b^(u-j)*x^j
	gen aj=1,cuj=1;
	for (unsigned j=0;j<=u;++j){
	  w[u-j]=(aj*cuj)*vbk[u-j];
	  aj=aj*a;
	  cuj=(int(u-j)*cuj)/int(j+1);
	}
	poly12polynome(w,i+1,res,th.dim);
	return;
      }
    }
    for (int i=0;i<th.dim;++i){
      pid1 = pid1*unsigned(d1[i]+1);
    }
    for (int i=0;i<th.dim;++i){
      d[i]=d1[i]*u+1;
      ans = ans*unsigned(d[i]);
      if (ans/RAND_MAX>RAND_MAX)
	break;
    }
    unsigned c1c2=1; // guess number of terms for the hashmap
    if (ans<=RAND_MAX){
      vector< T_unsigned<gen,unsigned> > p1,p2,p;
      convert<gen,unsigned>(th,d,p1);
      p2=p1;
      for (unsigned i=1;i<u;++i){
	smallmult<gen,unsigned>(p1,p2,p,0,c1c2);
	p1=p;
      }
      convert<gen,unsigned>(p,d,res);
    }
    else // ans>RAND_MAX
      res=Tpow(th,u);
  }

  polynome operator - (const polynome & th) {  
    // Tensor addition
    polynome res(th.dim,th);
    vector< monomial<gen> >::const_iterator a = th.coord.begin();
    vector< monomial<gen> >::const_iterator a_end = th.coord.end();
    res.coord.reserve(a_end - a );  
    for (;a!=a_end;++a){
      res.coord.push_back(monomial<gen>(-(*a).value,(*a).index));
    }
    return res;
  }

  // exactquo==2 means we know that b divides a and we search the cofactor
  bool divrem1(const polynome & a,const polynome & b,polynome & quo,polynome & r,int exactquo) {
    quo.coord.clear();
    quo.dim=a.dim;
    r.dim=a.dim;
    r.coord.clear();
    int bs=b.coord.size();
    if ( b.dim<=1 || bs==1 || a.coord.empty() ){
      return a.TDivRem(b,quo,r,false) && (exactquo?r.coord.empty():true) ;
    }
    int bdeg=b.coord.front().index.iptr->front(),rdeg=a.lexsorted_degree(),ddeg=rdeg-bdeg;
    if (ddeg>3){
      index_t d1=a.degree(),d2=b.degree(),d3=*b.coord.front().index.iptr,d(a.dim);
      // i-th degrees of th / other in quotient and remainder
      // are <= i-th degree of th + ddeg*(i-th degree of other - i-th degree of lcoeff of other) 
      double ans=1;
      for (int i=0;i<a.dim;++i){
	if (exactquo==2)
	  d[i]=d1[i]+1;
	else
	  d[i]=d1[i]+(ddeg+1)*(d2[i]-d3[i])+1;
	int j=1;
	// round to newt power of 2
	for (;;j++){
	  if (!(d[i] >>= 1))
	    break;
	}
	d[i] = 1 << j;
	ans = ans*unsigned(d[i]);
	if (ans/RAND_MAX>RAND_MAX)
	  break;
      }
      bool doit=true;
      if (ans<RAND_MAX){
	std::vector<unsigned> vars(a.dim);
	vars[a.dim-1]=1;
	for (int i=a.dim-2;i>=0;--i){
	  vars[i]=d[i+1]*vars[i+1];
	}
	if (debug_infolevel>1)
	  std::cerr << "divrem1 convert " << clock() << std::endl;
	{
	  std::vector< T_unsigned<longlong,unsigned> > p1,p2,quot,remain;
	  longlong maxp1,maxp2;
	  doit=convert_int(a,d,p1,maxp1) && convert_int(b,d,p2,maxp2) && maxp1/RAND_MAX < RAND_MAX;
	  if (doit){
	    if (debug_infolevel>1)
	      std::cerr << "hashdivrem1 begin " << clock() << " maxp1=" << maxp1 << " maxp2=" << maxp2 << " ddeg=" << ddeg << std::endl;
	    if (hashdivrem<longlong,unsigned>(p1,p2,quot,remain,vars,/* reduce*/0,RAND_MAX/double(maxp2)/p2.size()*RAND_MAX,false)==1){
	      if (debug_infolevel>1)
		std::cerr << "hashdivrem1 end " << clock() << std::endl;
	      convert_from(quot,d,quo);
	      convert_from(remain,d,r);
	      return true;
	    }
	    doit=false;
	  }
	}
#ifdef HAVE_GMPXX_H
	if (mpzclass_allowed)
	{
	  std::vector< T_unsigned<myint,unsigned> > p1,p2,quot,remain;
	  if (debug_infolevel>1)
	    std::cerr << "divrem1mpz convert " << clock() << std::endl;
	  doit=convert_myint(a,d,p1) && convert_myint(b,d,p2);
	  if (doit){
	    if (debug_infolevel>1)
	      std::cerr << "hashdivrem1mpz begin " << clock() << " ddeg=" << ddeg << std::endl;
	    if (hashdivrem<myint,unsigned>(p1,p2,quot,remain,vars,/* reduce */ 0,/* no size check */0.0,false)==1){
	      if (debug_infolevel>1)
		std::cerr << "hashdivrem1mpz end " << clock() << std::endl;
	      convert_from(quot,d,quo);
	      convert_from(remain,d,r);
	      return true;
	    }
	  }
	}
#endif
      }
      if (doit && ans/RAND_MAX<RAND_MAX){
#ifdef __VISUALC__ // Visual C++?
	typedef unsigned __int64 ulonglong ;
#else
	typedef unsigned long long ulonglong;
#endif
	std::vector<ulonglong> vars(a.dim);
	vars[a.dim-1]=1;
	for (int i=a.dim-2;i>=0;--i){
	  vars[i]=d[i+1]*vars[i+1];
	}
	if (debug_infolevel>1)
	  std::cerr << "divrem1 convert " << clock() << std::endl;
	{
	  std::vector< T_unsigned<longlong,ulonglong> > p1,p2,quot,remain;
	  longlong maxp1,maxp2;
	  doit=convert_int(a,d,p1,maxp1) && convert_int(b,d,p2,maxp2) && maxp1/RAND_MAX < RAND_MAX;
	  if (doit){
	    if (debug_infolevel>1)
	      std::cerr << "hashdivrem1 begin " << clock() << " maxp1=" << maxp1 << " maxp2=" << maxp2 << " ddeg=" << ddeg << std::endl;
	    if (hashdivrem<longlong,ulonglong>(p1,p2,quot,remain,vars,/* reduce */0,RAND_MAX/double(maxp2)/p2.size()*RAND_MAX,false)==1){
	      if (debug_infolevel>1)
		std::cerr << "hashdivrem1 end " << clock() << std::endl;
	      convert_from(quot,d,quo);
	      convert_from(remain,d,r);
	      return true;
	    }
	  }
	}
#ifdef HAVE_GMPXX_H
	if (mpzclass_allowed)
	{
	  std::vector< T_unsigned<myint,ulonglong> > p1,p2,quot,remain;
	  longlong maxp1,maxp2;
	  if (debug_infolevel>1)
	    std::cerr << "divrem1mpz convert " << clock() << std::endl;
	  doit=convert_myint(a,d,p1) && convert_myint(b,d,p2);
	  if (doit){
	    if (debug_infolevel>1)
	      std::cerr << "hashdivrem1z begin " << clock() <<  " ddeg=" << ddeg << std::endl;
	    if (hashdivrem<myint,ulonglong>(p1,p2,quot,remain,vars,/* reduce */ 0,/* no size check */0.0,false)==1){
	      if (debug_infolevel>1)
		std::cerr << "hashdivrem1 end " << std::clock() << std::endl;
	      convert_from(quot,d,quo);
	      convert_from(remain,d,r);
	      return true;
	    }
	  }
	}
#endif
      }
    }
    return a.TDivRem1(b,quo,r,false,exactquo);
  }

  polynome operator / (const polynome & th,const polynome & other) {  
    polynome rem(th.dim,th),quo(th.dim,th);
    // if ( !(th).TDivRem1(other,quo,rem) )
    if ( !divrem1(th,other,quo,rem) )
      setsizeerr("Unable to divide, perhaps due to rounding error"+th.print()+" / "+other.print());
    return(quo);
  }

  polynome operator / (const polynome & th,const gen & fact ) {  
    if (fact==gen(1))
      return th;
    polynome res(th.dim,th);
    vector< monomial<gen> >::const_iterator a = th.coord.begin();
    vector< monomial<gen> >::const_iterator a_end = th.coord.end();
    Div(a,a_end,fact,res.coord);
    return res;
  }

  polynome operator % (const polynome & th,const polynome & other) {  
    polynome rem(th.dim,th),quo(th.dim,th);
    if ( !(th).TDivRem1(other,quo,rem) )
      setsizeerr("Unable to divide, perhaps due to rounding error"+th.print()+" % "+other.print());
    return(rem);
  }

  polynome operator % (const polynome & th, const gen & modulo) {  
    polynome res(th.dim,th);
    vector< monomial<gen> >::const_iterator a = th.coord.begin();
    vector< monomial<gen> >::const_iterator a_end = th.coord.end();
    res.coord.reserve(a_end - a );  
    for (;a!=a_end;++a){
      gen tmp((*a).value % modulo);
      if (!is_zero(tmp))
	res.coord.push_back(monomial<gen>(tmp,a->index));
    }
    return res;
  }

  polynome re(const polynome & th){
    return Tapply(th,giac::no_context_re);
  }

  polynome im(const polynome & th){
    return Tapply(th,giac::no_context_im);
  }

  polynome conj(const polynome & th){
    return Tapply(th,giac::no_context_conj);
  }

  void smod(const polynome & th, const gen & modulo,polynome & res){
    vector< monomial<gen> >::const_iterator a = th.coord.begin();
    vector< monomial<gen> >::const_iterator a_end = th.coord.end();
    res.coord.clear();
    res.coord.reserve(a_end - a );  
    for (;a!=a_end;++a){
      const gen & tmp=smod(a->value, modulo);
      if (!is_zero(tmp))
	res.coord.push_back(monomial<gen>(tmp,a->index));
    }
  }

  polynome smod(const polynome & th, const gen & modulo) {  
    polynome res(th.dim,th);
    smod(th,modulo,res);
    return res;
  }

  // var is the variable number to extract, from 1 to p.dim
  void polynome2poly1(const polynome & p,int var,vecteur & v){
    polynome pp(p);
    if (var!=1)
      pp.reorder(transposition(0,var-1,p.dim));
    v.clear();
    int current_deg=pp.lexsorted_degree();
    v.reserve(current_deg+1);
    vector< monomial<gen> >::const_iterator it=pp.coord.begin(),itend=pp.coord.end();
    for (;it!=itend;--current_deg){
      if (it->index.iptr->front()==current_deg){
	if (p.dim==1){
	  v.push_back(it->value);
	  ++it;
	}
	else
	  v.push_back(Tnextcoeff<gen>(it,itend));
      }
      else
	v.push_back(polynome(p.dim-1));
    }
    for (;current_deg>=0;--current_deg)
      v.push_back(zero);
  }

  vecteur polynome2poly1(const polynome & p,int var){
    vecteur v;
    polynome2poly1(p,var,v);
    return v;
  }

  // like polynome2poly1 for univariate p
  vecteur polynome12poly1(const polynome & p){
    int current_deg=p.lexsorted_degree();
    vecteur v;
    v.reserve(current_deg+1);
    vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    for (;it!=itend;--current_deg){
      if (it->index.iptr->front()==current_deg){
	v.push_back(it->value);
	++it;
      }
      else
	v.push_back(zero);
    }
    for (;current_deg>=0;--current_deg)
      v.push_back(zero);
    return v;
  }

  vecteur polynome2poly1(const polynome & p){
    vecteur v;
    int current_deg=p.lexsorted_degree();
    v.reserve(current_deg+1);
    vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    for (;it!=itend;--current_deg){
      if (it->index.iptr->front()==current_deg){
	v.push_back(it->value);
	++it;
      }
      else
	v.push_back(zero);
    }
    for (;current_deg>=0;--current_deg)
      v.push_back(zero);
    return v;
  }
  
  gen polynome2poly1(const gen & e,int var){
      if (e.type==_POLY)
          return polynome2poly1(*e._POLYptr,var);
      if (e.type!=_FRAC)
          return e;
      return fraction(polynome2poly1(e._FRACptr->num,var),polynome2poly1(e._FRACptr->den,var));
  }

  int inner_POLYdim(const vecteur & v){
    const_iterateur it=v.begin(),itend=v.end();
    int dim=1;
    for (;it!=itend;++it){
      if (it->type==_POLY){
	dim=it->_POLYptr->dim+1;
	break;
      }
    }
    return dim;
  }

  gen untrunc(const gen & e,int degree,int dimension){
    if (e.type==_POLY)
      return e._POLYptr->untrunc(degree,dimension);
    if (e.type==_EXT)
      return algebraic_EXTension(untrunc(*e._EXTptr,degree,dimension),untrunc(*(e._EXTptr+1),degree,dimension));
    if (e.type==_VECT){
      const_iterateur it=e._VECTptr->begin(),itend=e._VECTptr->end();
      vecteur res;
      res.reserve(itend-it);
      for (;it!=itend;++it)
	res.push_back(untrunc(*it,degree,dimension));
      return res;
    }
    if (e.type==_FRAC)
      return fraction(untrunc(e._FRACptr->num,degree,dimension),untrunc(e._FRACptr->den,degree,dimension));
    return tensor<gen>(monomial<gen>(e,degree,1,dimension));
  }

  gen vecteur2polynome(const vecteur & v,int dimension){
    const_iterateur it=v.begin(),itend=v.end();
    gen e;
    for (int d=itend-it-1;it!=itend;++it,--d){
      if (!is_zero(*it))
	e = e+untrunc(*it,d,dimension);
    }
    return e;
  }

  polynome poly12polynome(const vecteur & v){
    const_iterateur it=v.begin(),itend=v.end();
    polynome p(1);
    for (int d=itend-it-1;it!=itend;++it,--d){
      if (!is_zero(*it))
	p.coord.push_back(monomial<gen>(*it,d,1,1));
    }
    return p;
  }

  // WARNING: var begins at 1 and ends at dimension
  void poly12polynome(const vecteur & v, int var,polynome & p,int dimension){
    if (dimension)
      p.dim=dimension;
    else 
      p.dim=inner_POLYdim(v);
    p.coord.clear();
    const_iterateur it=v.begin(),itend=v.end();
    for (int d=(itend-it)-1;it!=itend;++it,--d){
      if (is_zero(*it))
	continue;
      if (it->type!=_POLY || (it->_POLYptr->dim+1)!=p.dim)
	p.coord.push_back(monomial<gen>(*it,d,1,p.dim));
      else {
	vector< monomial<gen> >::const_iterator p_it=it->_POLYptr->coord.begin(),p_itend=it->_POLYptr->coord.end();
	for (;p_it!=p_itend;++p_it)
	  p.coord.push_back(p_it->untrunc(d,p.dim));
      }
    }
    if (var!=1){
      p.reorder(transposition(0,var-1,p.dim));
    }
  }

  polynome poly1_2_polynome(const vecteur & v, int dimension){
    polynome p(dimension);
    const_iterateur it=v.begin(),itend=v.end();
    for (int d=(itend-it)-1;it!=itend;++it,--d){
      if (is_zero(*it))
	continue;
      p.coord.push_back(monomial<gen>(*it,d,1,p.dim));
    }
    return p;
  }

  polynome poly12polynome(const vecteur & v,int var,int dimension){
    polynome p(0);
    poly12polynome(v,var,p,dimension);
    return p;
  }

  // assuming pmod and qmod are prime together, find r such that
  // r = p mod pmod  and r = q mod qmod
  // hence r = p + A*pmod = q + B*qmod
  // or A*pmod -B*qmod = q - p
  // assuming u*pmod+v*pmod=d we get
  // A=u*(q-p)/d
  polynome ichinrem(const polynome &p,const polynome & q,const gen & pmod,const gen & qmod){
    gen u,v,d,tmp,pqmod(pmod*qmod);
    egcd(pmod,qmod,u,v,d);
    // cout << u << "*" << pmod << "+" << v << "*" << qmod << "=" << d << " " << u*pmod+v*qmod << endl;
    vector< monomial<gen> >::const_iterator a = p.coord.begin();
    vector< monomial<gen> >::const_iterator a_end = p.coord.end();
    vector< monomial<gen> >::const_iterator b = q.coord.begin();
    vector< monomial<gen> >::const_iterator b_end = q.coord.end();
    polynome res(p.dim);
    res.coord.reserve(a_end - a );
    for (;(a!=a_end)&&(b!=b_end);){
      if (a->index != b->index){
	if (a->index>=b->index){
	  tmp=a->value-rdiv(u*a->value,d);
	  res.coord.push_back(monomial<gen>(smod(tmp,pqmod),a->index));
	  ++a;
	}
	else {
	  tmp=rdiv(u*b->value,d);
	  res.coord.push_back(monomial<gen>(smod(tmp,pqmod),b->index));
	  ++b;
	}
      }
      else {
	tmp=a->value+rdiv(u*(b->value-a->value),d) *pmod ;
	// cout << a->value << " " << b->value << "->" << tmp << " " << pqmod << endl;
	res.coord.push_back(monomial<gen>(smod(tmp,pqmod),b->index));
	++b;
	++a;
      }
    }
    for (;a!=a_end;++a)
      res.coord.push_back(monomial<gen>(smod(a->value-rdiv(u*(a->value),d),pqmod),a->index));
    for (;b!=b_end;++b)
      res.coord.push_back(monomial<gen>(smod(rdiv(u*b->value,d),pqmod),b->index));
    return res;
  }

  bool divrem (const polynome & th, const polynome & other, polynome & quo, polynome & rem, bool allowrational ){
    return th.TDivRem(other,quo,rem,allowrational);
  }

  int inspectdim=20;
  bool exactquotient(const polynome & a,const polynome & b,polynome & quo){
    clock_t beg=clock(),delta;
    bool res= a.Texactquotient(b,quo,true);
    delta=clock()-beg;
    if (delta && debug_infolevel) // a.dim>=inspectdim
      cerr << "exactquo end " << delta << " " << res << endl;
    return res;
  }

  bool divremmod2 (const polynome & th,const polynome & other, const gen & modulo,polynome & quo, polynome & rem) {
    int asize=th.coord.size();
    if (!asize){
      quo=th;
      rem=th; 
      return true;
    }
    int bsize=other.coord.size();
    if (bsize==0)  setsizeerr("gausspol.cc/divremmod2");
    index_m a_max = th.coord.front().index;
    index_m b_max = other.coord.front().index;
    quo.coord.clear();
    quo.dim=th.dim;
    rem.dim=th.dim;
    if ( (bsize==1) && (b_max==b_max*0) ){
      rem.coord.clear();
      gen b=other.coord.front().value;
      if (b==gen(1))
	quo = th ;
      else {
	b=invmod(b,modulo);
	vector< monomial<gen> >::const_iterator itend=th.coord.end();
	for (vector< monomial<gen> >::const_iterator it=th.coord.begin();it!=itend;++it)
	  quo.coord.push_back(monomial<gen>( smod(it->value*b,modulo),it->index));
      }
    return true;
    }
    rem=th;
    if ( ! (a_max>=b_max) ){
      // test that the first power of a_max is < to that of b_max
      return (a_max.iptr->front()<b_max.iptr->front());
    }
    bool mult=is_cinteger(other.coord.front().value);
    gen b=invmod(other.coord.front().value,modulo);
    while (a_max >= b_max){
      gen q=smod(rem.coord.front().value*b, modulo);
      quo.coord.push_back(monomial<gen>(q,a_max-b_max));
      polynome temp=other.shift(a_max-b_max,q);
      rem = smod(rem-temp, modulo);
      if (rem.coord.size())
	a_max=rem.coord.front().index;
      else
	break;
    }
    return(true);
  }

  bool divremmod (const polynome & th,const polynome & other, const gen & modulo,polynome & quo, polynome & r) {
    quo.coord.clear();
    quo.dim=th.dim;
    r.dim=th.dim;
    if ( (th.dim<=1) || (th.coord.empty()) )
      return divremmod2(th,other,modulo,quo,r);
    int os=other.coord.size();
    if (!os){
      r=th;
      return true;
    }
    if (os==1){
      // Check for a division by 1
      if (is_one(other.coord.front().value) && is_zero(*other.coord.front().index.iptr)){
	quo=th;
	r.coord.clear();
	return true;
      }
      // IMPROVE Invert other.coord.front() and shift/multiply
      
    }
    std::vector< monomial<gen> >::const_iterator it=other.coord.begin();
    int bdeg=it->index.iptr->front(),rdeg=th.lexsorted_degree(),ddeg=rdeg-bdeg;
    if (//false && 
	ddeg>2 && os>10 
	){
      index_t d1=th.degree(),d2=other.degree(),d3=*other.coord.front().index.iptr,d(th.dim);
      // i-th degrees of th / other in quotient and remainder
      // are <= i-th degree of th + ddeg*(i-th degree of other - i-th degree of lcoeff of other) 
      double ans=1;
      for (int i=0;i<th.dim;++i){
	d[i]=d1[i]+(ddeg+1)*(d2[i]-d3[i])+1;
	int j=1;
	// round to 2^
	for (;;j++){
	  if (!(d[i] >>= 1))
	    break;
	}
	d[i] = 1 << j;
	ans = ans*unsigned(d[i]);
	if (ans/RAND_MAX>RAND_MAX)
	  break;
      }
      if (ans<RAND_MAX){
	if (modulo.type==_INT_ && modulo.val<46340 && modulo.val>0){ 
	  // convert everything to integers
	  vector< T_unsigned<int,unsigned> > p1,p2,quot,remain;
	  vector<unsigned> vars(th.dim);
	  vars[th.dim-1]=1;
	  for (int i=th.dim-2;i>=0;--i){
	    vars[i]=d[i+1]*vars[i+1];
	  }
	  if (debug_infolevel>1)
	    cerr << "divrem convert " << clock() << endl;
	  if (convert(th,d,p1,modulo.val) && convert(other,d,p2,modulo.val)){
	    if (debug_infolevel>1)
	      cerr << "hashdivrem begin " << clock() << endl;
	    if (hashdivrem<int,unsigned>(p1,p2,quot,remain,vars,modulo.val,0.0)==1){
	      if (debug_infolevel>1)
		cerr << "hashdivrem end " << clock() << endl;
	      convert(quot,d,quo);
	      convert(remain,d,r);
	      return true;
	    }
	    else
	      return false;
	  }
	} // end modulo.type==_INT
      } // end ans < RAND_MAX
      if (ans/RAND_MAX<RAND_MAX){
	if (modulo.type==_INT_ && modulo.val<46340 && modulo.val>0){ 
	  // convert everything to integers
	  vector< T_unsigned<int,ulonglong> > p1,p2,quot,remain;
	  vector<ulonglong> vars(th.dim);
	  vars[th.dim-1]=1;
	  for (int i=th.dim-2;i>=0;--i){
	    vars[i]=d[i+1]*vars[i+1];
	  }
	  if (convert(th,d,p1,modulo.val) && convert(other,d,p2,modulo.val)){
	    if (hashdivrem<int,ulonglong>(p1,p2,quot,remain,vars,modulo.val,0.0)==1){
	      convert(quot,d,quo);
	      convert(remain,d,r);
	      return true;
	    }
	    else
	      return false;
	  }
	} // end modulo.type==_INT
      } // end ans/RAND_MAX < RAND_MAX
    }
    tensor<gen> b0(Tnextcoeff<gen>(it,other.coord.end()));
    r=th;
    tensor<gen> q(b0.dim),q_other(th.dim);
    while ( (rdeg=r.lexsorted_degree()) >=bdeg){
      it=r.coord.begin();
      tensor<gen> a0(Tnextcoeff<gen>(it,r.coord.end())),tmp(a0.dim);
      if (!divremmod(a0,b0,modulo,q,tmp) || !tmp.coord.empty())
	return false;
      q=q.untrunc1(rdeg-bdeg);
      quo=quo+q;
      mulpoly(q,other,q_other,modulo);
      r=smod(r-q_other,modulo);
      if (r.coord.empty())
	return true;
    }
    return true;
  }


  polynome pow(const polynome & p,const gen & n){
    if (!n.is_integer())  setsizeerr("gausspol.cc/pow");
    polynome res(p.dim);
    int i=n.to_int();
    powpoly(p,i,res);
    return res;
  }

  polynome pow(const polynome & p, int  n){
    polynome res(p.dim);
    powpoly(p,n,res);
    return res;
  }

  polynome powmod1(const polynome &p,int n,const gen & modulo){
    switch (n) {
    case 0: 
      return polynome(gen(1),p.dim);
    case 1: 
      return p; 
    default: 
      polynome temp(powmod(p,n/2,modulo));
      if (n%2)
	return (temp * temp * p) % modulo;
      else
	return (temp*temp) % modulo;
    }
  }

  polynome powmod(const polynome &p,int n,const gen & modulo){
    if (p.dim<2)
      return powmod1(p,n,modulo);
    polynome res(gen(1),p.dim);
    for (int i=0;i<n;i++)
      res=(res*p) % modulo;
    return res;
  }

  polynome resultant(const polynome & p,const polynome & q){
    // FIXME polynomial subresultant does not work if p and q have approx coeff
    return Tresultant<gen>(p,q);
  }

  polynome lgcd(const polynome & p){
    return Tlgcd<gen>(p);
  }

  gen ppz(polynome & p){
    return Tppz(p);
  }

  // Find the content of p with respect to the 1st variable
  // p=p(x1,...,xn)=poly in x1 with coeff depending on x2,..,xn
  // content wrt x1 depends on x2,...,xn
  void lgcdmod(const polynome & p,const gen & modulo,polynome & pgcd){
    if (!p.dim){
      pgcd=p;
      return ;
    } 
    pgcd=pgcd.trunc1();
    vector< monomial<gen> >::const_iterator it=p.coord.begin();
    vector< monomial<gen> >::const_iterator itend=p.coord.end();
    // vector< monomial<gen> >::const_iterator itbegin=it;
    for (;it!=itend;){
      if (is_one(pgcd))
	break;
      pgcd=gcdmod(pgcd,Tnextcoeff<gen>(it,itend),modulo);
    }
    if (pgcd.coord.empty()){
      index_m i;
      for (int j=0;j<p.dim;j++)
	i.iptr->push_back(0);
      pgcd.coord.push_back(monomial<gen>(gen(1),i));
    }
    else
      pgcd=pgcd.untrunc1();
  }

  // Split a multivariate poly X_1...X_n as multivar X_dim+1...X_n
  // with coeff multivar poly of X_1..X_dim
  polynome split(const polynome & p,int inner_dim){
    int outer_dim=p.dim-inner_dim;
    polynome cur_inner(inner_dim);
    polynome res(outer_dim);
    vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    for (; it!=itend;++it){
      index_t cur_index=*(it->index.iptr);
      index_t outer_index(cur_index.begin()+inner_dim,cur_index.end());
      index_t inner_index(cur_index.begin(),cur_index.begin()+inner_dim);
      cur_inner=polynome(monomial<gen>(it->value,inner_index));
      res=res+polynome(monomial<gen>(cur_inner,outer_index));
    }
    return res;
  }

  /*
#ifdef HASH_MAP_NAMESPACE
  class hash_function_index_t {
  public:
    inline size_t operator () (const index_t & a) const { 
      size_t r=0;
      index_t::const_iterator ita=a.begin(),itaend=a.end();
      for (;ita!=itaend;++ita){
	r <<= 4;
	r += (*ita) & 0xf;
      }
      return r;
    }
    hash_function_index_t() {};
  };

  typedef HASH_MAP_NAMESPACE::hash_map<index_t,polynome,hash_function_index_t> map_index_t_polynome;
#else
  typedef std::map<index_t,polynome> map_index_t_polynome;
#endif
  */

  typedef std::map<index_t,polynome> map_index_t_polynome;

  // return true if content=1 is detected 
  bool split(const polynome & p,int inner_dim,map_index_t_polynome & res){
    // int outer_dim=p.dim-inner_dim;
    vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    for (; it!=itend;++it){
      index_t cur_index=*(it->index.iptr);
      index_t outer_index(cur_index.begin()+inner_dim,cur_index.end());
      index_t inner_index(cur_index.begin(),cur_index.begin()+inner_dim);
      map_index_t_polynome::iterator jt=res.find(outer_index),jtend=res.end();
      if (jt==jtend){
	if (is_zero(inner_index))
	  return true;
	res[outer_index]=polynome(monomial<gen>(it->value,inner_index));
      }
      else
	jt->second.coord.push_back(monomial<gen>(it->value,inner_index));
    }
    return false;
  }

  gen lcoeffn(const polynome & p){
    int dim=p.dim;
    polynome res(dim);
    vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    if (it==itend)
      return 0;
    index_t i= *it->index.iptr; 
    for (;it!=itend;++it){
      index_t & j= *it->index.iptr; 
      i[dim-1]=j[dim-1];
      if (i!=j)
	break;
      res.coord.push_back(*it);
    }
    return res;
  }

  gen lcoeff1(const polynome & p){
    if (p.coord.empty())
      return zero;
    int inner_dim=1;
    // int outer_dim=p.dim-inner_dim;
    polynome cur_inner(inner_dim);
    vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    index_t::const_iterator jt0 = it->index.iptr->begin(),jtend=it->index.iptr->end(),jt,kt0,kt;
    for (; it!=itend;++it){
      kt0 = it->index.iptr->begin();
      for (jt=jt0+inner_dim,kt=kt0+inner_dim;jt!=jtend;++kt,++jt){
	if (*kt<*jt)
	  break;
	if (*kt>*jt){
	  jt0=kt0;
	  jtend=kt0+p.dim;
	  cur_inner.coord.clear();
	  jt=jtend;
	  break;
	}
      }
      if (jt==jtend)
	cur_inner.coord.push_back(monomial<gen>(it->value,index_t(kt0,kt0+inner_dim)));
    }
    return cur_inner;
  }

  polynome content1mod(const polynome & p,const gen & modulo,bool setdim=true){
    if (p.coord.empty())
      setsizeerr();
    if (p.coord.size()==1){
      int n=p.coord.front().index.iptr->front();
      polynome c(monomial<gen>(p.coord.front().value,vector<int>(1,n)));
      if (setdim)
	change_dim(c,p.dim);
      return c;
    }
    // New code
    map_index_t_polynome m;
    polynome c(1);
    if (!split(p,1,m)){
      int c0=RAND_MAX,i0;
      map_index_t_polynome::iterator it=m.begin(),itend=m.end();
      if (m.size()==1)
	c=it->second;
      else {
	for (;c0 && it!=itend;++it){
	  if (!it->second.coord.empty() && (i0=it->second.coord.front().index.iptr->front())<c0){
	    c=it->second;
	    c0=i0;
	  }
	}
	it=m.begin();
	for (;it!=itend;++it){
	  if (!c.coord.empty() &&c.coord.front().index.iptr->front()==0 ){
	    c=polynome(plus_one,1);
	    break;
	  }
	  c=gcdmod(c,it->second,modulo);
	}
	/* Old code
	   polynome lp(split(p,1)),c(1);
	   vector< monomial<gen> >::const_iterator it=lp.coord.begin(),itend=lp.coord.end();
	   for (;it!=itend;++it){
	   if (it->value.type==_POLY)
	   c=gcdmod(c,*it->value._POLYptr,modulo);
	   else {
	   c=polynome(plus_one,1); // was c=polynome(plus_one,p.dim-1);
	   break;
	   }
	   }
	*/
      } // end if itend==it+1
    } // end if (!split()) : i.e. if the content is not trivially 1
    else
      c=polynome(plus_one,1);
    if (setdim)
      change_dim(c,p.dim);
    return c;
  }
  
  polynome pp1mod(const polynome & p,const gen & modulo){
    polynome q(p.dim),r(p.dim);
    polynome tmp(content1mod(p,modulo));
    // cerr << "pp1mod " << tmp << endl;
    divremmod(p,tmp,modulo,q,r);
    return q;
  }

  // Find non zeros coeffs of p
  int find_nonzero(const polynome & p,index_t & res){
    res.clear();
    vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    if (it==itend)
      return 0;
    int old_deg=it->index.iptr->front(),cur_deg=0;
    int nzeros=0;
    res.push_back(1);
    for (;it!=itend;++it){
      cur_deg=it->index.iptr->front();
      if (cur_deg!=old_deg){
	nzeros += old_deg - cur_deg -1 ;
	for (int i=old_deg-cur_deg;i>1;--i)
	  res.push_back(0);
	res.push_back(1);
	old_deg=cur_deg;
      }
    }
    if (cur_deg){
      nzeros += cur_deg;
      for (int i=cur_deg;i>0;--i)
	res.push_back(0);
    }
    return nzeros;
  }

  // p_orig and q_orig are primitive with respect to the main variable
  // p(x1,...,xn) q(x1,...,xn) viewed as p(x1) and q(x1) 
  // d must be the same
  void mod_gcdmod(const polynome &p_orig, const polynome & q_orig, const gen & modulo, polynome & d,int gcddeg=0){
    if (p_orig.coord.empty() || is_one(q_orig)){
      d=q_orig;
      return;
    }
    if (q_orig.coord.empty() || is_one(p_orig)){
      d=p_orig;
      return;
    }
    if (debug_infolevel)
      cerr << "gcdmod content dim " << d.dim << " " << clock() << endl;
    polynome p(p_orig.dim),q(q_orig.dim),r;
    vector<int_unsigned> pint,qint;
    index_t pintd=p_orig.degree(),qintd=q_orig.degree();
    unsigned pu,qu;
    // Make p and q primitive with respect to x2,...,xn
    // i.e. the coeff of p and q which are polynomials in x1
    // are relative prime 
    // r is the gcd of the content in this sense
    bool docontent1mod=true;
    // bug there, if false is removed lwN and lwN1 do not work
    if (false &&
	modulo.type==_INT_ && degree2unsigned(pintd,pu) && degree2unsigned(qintd,qu)){
      convert(p_orig,pintd,pint,modulo.val); 
      convert(q_orig,qintd,qint,modulo.val);
      if (is_content_trivially_1(pint,pu/pintd.front()) && is_content_trivially_1(qint,qu/qintd.front()))
	docontent1mod=false;
    }
    if (docontent1mod) {
      polynome pc(content1mod(p_orig,modulo)),qc(content1mod(q_orig,modulo));
      divremmod(p_orig,pc,modulo,p,r);
      divremmod(q_orig,qc,modulo,q,r);
      // cerr << "content end " << clock() << endl;
      change_dim(pc,1); change_dim(qc,1);
      r=gcdmod(pc,qc,modulo);
      change_dim(r,p.dim);
    }
    else {
      p=p_orig; 
      q=q_orig;
      r=polynome(plus_one,p.dim);
    }
    // Find degree of gcd with respect to x1, more precisely gcddeg>=degree/x1
    // and compute data for the sparse modular algorithm
    vector<int> vzero; // coeff of vzero correspond to zero or non zero
    int nzero=1; // Number of zero coeffs
    vecteur alphav,gcdv; // Corresponding values of alpha and gcd at alpha
    if (!gcddeg){
      vecteur b(p.dim-1);
      for (int essai=0;essai<2;++essai){
	if (essai)
	  b=vranm(p.dim-1,0,0); // find another random point
	polynome Fb(1),Gb(1);
	// Fb and Gb are p and q where x2,...,xn are evaluated at b
	if (!find_good_eval(p,q,Fb,Gb,b,(debug_infolevel>=2),modulo))
	  break;
	polynome Db(gcdmod(Fb,Gb,modulo)); // 1-d gcd wrt x1
	int Dbdeg=Db.lexsorted_degree();
	if (!Dbdeg){
	  gcddeg=0;
	  break;
	}
	if (!gcddeg){ // 1st gcd test
	  gcddeg=Dbdeg;
	  nzero=find_nonzero(Db,vzero);
	}
	else { // 2nd try
	  if (Dbdeg<gcddeg){ // 1st try unlucky, restart 1st try
	    gcddeg=Dbdeg;
	    nzero=find_nonzero(Db,vzero);
	    --essai;
	  }
	  else {
	    if (Dbdeg>gcddeg) // 2nd try unlucky, restart 2nd try
	      --essai;
	    else { // Same gcd degree for 1st and 2nd try, keep this degree
	      vector<int> tmp;
	      nzero=find_nonzero(Db,tmp);
	      if (nzero){
		vzero = vzero | tmp;
		// Recompute nzero, it is the number of 0 coeff of vzero
		vector<int>::const_iterator it=vzero.begin(),itend=vzero.end();
		for (nzero=0;it!=itend;++it){
		  if (!*it) ++nzero;
		}
	      }
	    }
	  }
	}
      }
    }
    else {
      gcddeg -= r.lexsorted_degree() ;
      nzero = 0; // No info available
    }
    if (!gcddeg){
      d=r;
      return;
    }
    d=polynome(p.dim);
    polynome interp(plus_one,p.dim);
    // gcd of leading coefficients of p and q viewed as poly in X_2...X_n
    // with coeff in Z[X_1]
    if (debug_infolevel)
      cerr << "gcdmod lcoeff1 dim " << d.dim << " " << clock() << endl;
    gen lp(lcoeff1(p)),lq(lcoeff1(q));
    polynome Delta(plus_one,p.dim);
    if ((lp.type==_POLY) && (lq.type==_POLY) )
      Delta=gcdmod(*lp._POLYptr,*lq._POLYptr,modulo);
    // we are now interpolating G=gcd(p,q)*a poly/x1
    // such that the leading coeff of G is Delta
    vector<int> pdeg(p.degree()),qdeg(q.degree()); 
    int spdeg=0,sqdeg=0;
    for (int i=1;i<p.dim;++i){
      spdeg += pdeg[i];
      sqdeg += qdeg[i];
    }
    vector<int> delta=index_min(vector<int>(pdeg.begin()+1,pdeg.end()),
			  vector<int>(qdeg.begin()+1,qdeg.end()));
    int e=0; // number of evaluations
    int alpha=0;
    if (debug_infolevel>1)
      cerr << "gcdmod find alpha dim " << d.dim << " " << clock() << endl;
    for (;;++alpha){
      vecteur valpha;
      polynome palpha(p.dim-1),qalpha(q.dim-1);
      for (;alpha<modulo.val;++alpha){
	valpha=vecteur(1,alpha);
	gen tmp(peval(p,valpha,modulo,false,&pint));
	if (is_zero(tmp))
	  continue;
	if (tmp.type!=_POLY){ 
	  if (spdeg) 
	    continue;
	  // gcd may only depend on first var
	  d=r;
	  return;
	}
	palpha=smod(*tmp._POLYptr,modulo);
	tmp=peval(q,valpha,modulo,false,&qint);
	if (is_zero(tmp)) 
	  continue;
	if (tmp.type!=_POLY){ 
	  if (sqdeg) 
	    continue;
	  d=r;
	  return;
	}
	qalpha=smod(*tmp._POLYptr,modulo);
	if ( palpha.lexsorted_degree()==pdeg[1] &&
	     qalpha.lexsorted_degree()==qdeg[1] )
	  break;
      }
      // palpha and qalpha are p and q evaluated at x1=alpha
      if (debug_infolevel>1)
	cerr << "gcdmod eval " << alpha << " dim " << d.dim << " " << clock() << endl;
      if (alpha==modulo)
	setsizeerr("Modgcd: no suitable evaluation point");
      polynome g(gcdmod(palpha,qalpha,modulo));
      vector<int> gdeg(g.degree());
      int gcd_plus_delta_deg=gcddeg+Delta.lexsorted_degree();
      if (gdeg==delta){
	// Try spmod first
	if (nzero){
	  // Add alpha,g 
	  alphav.push_back(alpha);
	  gcdv.push_back(g);
	  if (gcddeg-nzero==e){ 
	    // We have enough evaluations, let's try SPMOD
	    // Build the matrix, each line has coeffs / vzero
	    matrice m;
	    for (int j=0;j<=e;++j){
	      vector<int>::const_reverse_iterator it=vzero.rbegin(),itend=vzero.rend();
	      vecteur line;
	      for (gen p=alphav[j],pp=plus_one;it!=itend;++it,pp=smod(p*pp,modulo)){
		if (*it)
		  line.push_back( pp);
	      }
	      reverse(line.begin(),line.end());
	      line.push_back(gcdv[j]);
	      m.push_back(line);
	    }
	    // Reduce linear system modulo modulo
	    gen det; vecteur pivots; matrice mred;
	    cerr << "SPMOD " << clock() << endl;
	    modrref(m,mred,pivots,det,0,m.size(),0,m.front()._VECTptr->size()-1,true,false,modulo,false);
	    cerr << "SPMODend " << clock() << endl;
	    if (!is_zero(det)){	      
	      // Last column is the solution, it should be polynomials
	      // that must be untrunced with index = to non-0 coeff of vzero
	      polynome trygcd(p.dim);
	      vector<int>::const_iterator it=vzero.begin(),itend=vzero.end();
	      int deg=itend-it-1;
	      for (int pos=0;it!=itend;++it,--deg){
		if (!*it)
		  continue;
		gen tmp=mred[pos][e+1]; // e+1=#of points -> last col
		if (tmp.type==_POLY)
		  trygcd=trygcd+tmp._POLYptr->untrunc1(deg);
		else
		  if (!is_zero(tmp))
		    trygcd=trygcd+polynome(monomial<gen>(tmp,deg,1,p.dim));
		++pos;
	      }
	      // Check if trygcd is the gcd!
	      polynome pD(pp1mod(trygcd,modulo)),Q(p.dim),R(d.dim);
	      divremmod(p,pD,modulo,Q,R);
	      if (R.coord.empty()){
		divremmod(q,pD,modulo,Q,R);
		if (R.coord.empty()){
		  pD=pD*r;
		  d=smod(pD*invmod(pD.coord.front().value,modulo),modulo);
		  return;
		}
	      }
	    }
	    // SPMOD not successfull :-(
	    nzero=0;
	  } // end if gcddeg-nzero==e
	} // end if (nzero)
	if (debug_infolevel>1)
	  cerr << "gcdmod interp dim " << d.dim << " " << clock() << endl;
	polynome g1=(g*smod(peval(Delta,valpha,modulo),modulo))*invmod(g.coord.front().value,modulo);
	gen tmp(g1-peval(d,valpha,modulo));
	if (tmp.type==_POLY){
	  g1=smod(*tmp._POLYptr,modulo);
	  g1=g1.untrunc1();
	}
	else
	  g1=polynome(tmp,p.dim);
	d=d+g1*interp*invmod(peval(interp,valpha,modulo),modulo);
	d=smod(d,modulo);
	interp=interp*(polynome(monomial<gen>(plus_one,1,1,p.dim))-polynome(alpha,p.dim));
	++e;
	if (e>gcddeg
	    || is_zero(tmp)
	    ){
	  if (debug_infolevel)
	    cerr << "gcdmod pp1mod dim " << d.dim << " " << clock() << endl;
	  polynome pD(pp1mod(d,modulo)),Q(p.dim),R(d.dim);
	  // This removes the polynomial in x1 that we multiplied by
	  // (it was necessary to know the lcoeff of the interpolated poly)
	  if (debug_infolevel)
	    cerr << "gcdmod check dim " << d.dim << " " << clock() << endl;
	  // Now, gcd divides pD for gcddeg+1 values of x1
	  // degree(pD)<=degree(gcd)
	  divremmod(p,pD,modulo,Q,R);
	  if (debug_infolevel){
	    cerr << "test * " << clock() << endl;
	    polynome R2;
	    mulpoly(pD,Q,R2,modulo);
	    cerr << "test * end " << clock() << endl;
	  }
	  if (R.coord.empty()){
	    divremmod(q,pD,modulo,Q,R);
	    // If pD divides both P and Q, then the degree wrt variables
	    // x2,...,xn is the right one (because it is <= since pD 
	    // divides the gcd and >= since pD(x1=one of the try) was a gcd
	    // The degree in x is the right one because of the condition
	    // on the lcoeff
	    // Note that the division test might be much longer than the
	    // interpolation itself (e.g. if the degree of the gcd is small)
	    // but it seems unavoidable, for example if 
	    // P=Y-X+X(X-1)(X-2)(X-3)
	    // Q=Y-X+X(X-1)(X-2)(X-4)
	    // then gcd(P,Q)=1, but if we take Y=0, Y=1 or Y=2
	    // we get gcddeg=1 (probably degree 1 for the gcd)
	    // interpolation at X=0 and X=1 will lead to Y-X as candidate gcd
	    // and even adding X=2 will not change it
	    // We might remove division if we compute the cofactors of P and Q
	    // if P=pD*cofactor is true for degree(P) values of x1
	    // and same for Q, and the degrees wrt x1 of pD and cofactors
	    // have sum equal to degree of P or Q then pD is the gcd
	    if (R.coord.empty()){
	      pD=pD*r;
	      d=smod(pD*invmod(pD.coord.front().value,modulo),modulo);
	      if (debug_infolevel)
		cerr << "gcdmod found dim " << d.dim << " " << clock() << endl;
	      return;
	    }
	  }
	  if (debug_infolevel)
	    cerr << "Gcdmod bad guess " << endl;
	  continue;
	}
	else
	  continue;
      }
      if (gdeg[0]>delta[0]) // branch if all degree are >=
	continue;
      if (delta[0]>=gdeg[0]){ // restart with g
	gcdv=vecteur(1,g);
	alphav=vecteur(1,alpha);
	delta=gdeg;
	g=(g*smod(peval(Delta,valpha,modulo),modulo))*invmod(g.coord.front().value,modulo);
	d=g.untrunc1();
	e=1;
	interp=polynome(monomial<gen>(plus_one,1,1,p.dim))-polynome(alpha,p.dim);
	continue;
      }
    }
  }

  void psrgcdmod(polynome & a,polynome & b,const gen & modulo,polynome & prim){
    // set auxiliary polynomials g and h to 1
    polynome g(gen(1),a.dim);
    polynome h(g),quo(g),r(g);
    while (!a.coord.empty()){
      int n=b.lexsorted_degree();
      int m=a.lexsorted_degree();
      if (!n) {// if b is constant (then b!=0), gcd=original lgcdmod
	prim=polynome(gen(1),a.dim);
	return ;
      }
      int ddeg=m-n;
      if (ddeg<0)
	swap(a,b); // exchange a<->b may occur only at the beginning
      else {
	polynome b0(firstcoeff(b));
	divremmod(a*pow(b0,ddeg+1),b,modulo,quo,r); // division works always
	if (r.coord.empty())
	  break;
	// remainder is non 0, loop continue: a <- b
	a=b;
	polynome temp(powmod(h,ddeg,modulo));
	// now divides r by g*h^(m-n), result is the new b
        divremmod(r,g*temp,modulo,b,quo); // quo is the remainder here, not used
	// new g=b0 and new h=b0^(m-n)*h/temp
	if (ddeg==1) // the normal case, remainder deg. decreases by 1 each time
	  h=b0;
	else // not sure if it's better to keep temp or divide by h^(m-n+1)
	  divremmod(pow(b0,ddeg)*h,temp,modulo,h,quo);
	g=b0;
      }
    }
    // cout << "Prim" << b << endl;
    quo.coord.clear();
    lgcdmod(b,modulo,quo);
    divremmod(b,quo,modulo,prim,r);
    prim=smod(prim*invmod(prim.coord.front().value,modulo),modulo);
  }
  
  void contentgcdmod(const polynome &p, const polynome & q, const gen & modulo, polynome & cont,polynome & prim){
    if (p.coord.empty()){
      cont.coord.clear();
      lgcdmod(q,modulo,cont);
      polynome temp(cont.dim);
      divremmod(q,cont,modulo,prim,temp);
      return ;
    }
    if (q.coord.empty()){
      contentgcdmod(q,p,modulo,cont,prim);
      return;
    }
    if (p.dim!=q.dim)  setsizeerr("gausspol.cc/contentgcdmod");
    // dp and dq are the "content" of p and q w.r.t. other variables
    polynome dp(p.dim), dq(p.dim);
    // cerr << p.dim << " " << clock() << endl;
    lgcdmod(p,modulo,dp);
    lgcdmod(q,modulo,dq);
    // cerr << "End " << p.dim << " " << clock() << endl;
    cont=gcdmod(dp.trunc1(),dq.trunc1(),modulo).untrunc1();
    if (!p.dim){
      prim=polynome(gen(1),0);
      return ;
    }
    // cout << "Cont" << cont << endl; 
    polynome a(p.dim),b(p.dim),quo(p.dim),r(p.dim);
    // a and b are the primitive part of p and q
    divremmod(p,dp,modulo,a,r);
    divremmod(q,dq,modulo,b,r);
    if (modulo.val>=4*min(p.lexsorted_degree(),q.lexsorted_degree())){
      mod_gcdmod(a,b,modulo,prim);
      return ;
    }
    psrgcdmod(a,b,modulo,prim);
  }

  bool gcdmod_dim1(const polynome &p,const polynome & q,const gen & modulo,polynome & d,polynome & pcof,polynome & qcof,bool compute_cof,bool & real){
    real= poly_is_real(p) && poly_is_real(q);
    if (p.dim!=1)
      return false;
    if (q.dim!=1)
      setsizeerr();
    d.dim=pcof.dim=qcof.dim=1;
    if (real && modulo.type==_INT_ && modulo.val<46340 && gcdsmallmodpoly(p,q,modulo.val,d,pcof,qcof,compute_cof)){
      return true; 
    }
    modpoly P(polynome2poly1(p,1));
    modpoly Q(polynome2poly1(q,1));
    environment envi;
    environment * env=&envi;
    env->modulo=modulo;
    env->pn=env->modulo;
    env->moduloon=true;
    env->complexe=true;
    modpoly R,PQ,PR;
    gcdmodpoly(P,Q,env,R);
    d=poly12polynome(R);
    if (compute_cof){
      DivRem(P,R,env,PQ,PR);
      pcof=poly12polynome(PQ);
      DivRem(Q,R,env,PQ,PR);
      qcof=poly12polynome(PQ);
    }
    return true;
  }

  polynome gcdmod(const polynome &p,const polynome & q,const gen & modulo){
    if (p.dim!=q.dim)
      setsizeerr("Bug!");
    if (p==q)
      return p;
    if (p.coord.empty())
      return q;
    if (q.coord.empty())
      return p;
    if (p.dim==1){
      polynome d(1),pd(1),qd(1);
      bool estreel;
      gcdmod_dim1(p,q,modulo,d,pd,qd,false,estreel);
      return d;
    }
    // Check that there are enough points for interpolation
    // Otherwise PSR
    if (modulo.val>=4*min(p.lexsorted_degree(),q.lexsorted_degree())){
      polynome d(p.dim),pcof(p.dim),qcof(p.dim);
      if (modgcd(p,q,modulo,d,pcof,qcof,false))
	return d;
    }
    polynome a(smod(p*invmod(p.coord.front().value,modulo),modulo));
    polynome b(smod(q*invmod(q.coord.front().value,modulo),modulo));
    // Use evaluation points if enough available or modular psrh
    polynome prim(p.dim),cont(p.dim);
    contentgcdmod(a,b,modulo,prim,cont);
    if (debug_infolevel>10)
      cout << "Prim" << prim << "Cont" << cont << endl;
    return smod(prim*cont, modulo);
  }

  /*
    p and q are assumed to have integer content=1
    the leading coeff of d=gcd(p,q) divides the leading coeff of p and q
    we will therefore normalize modular gcds to have the gcd of the
    leading coeffs as leading coeff, and will try divisibility 
    of it's smodular representant after division by the content
   */
  bool gcd_modular_algo(polynome &p,polynome &q, polynome &d,bool compute_cof){
    if (p.dim==1)
      return gcd_modular_algo1(p,q,d,compute_cof);
    polynome plgcd(p.dim), qlgcd(q.dim), pp(p.dim), qq(p.dim),gcdlgcd(p.dim);
    plgcd=lgcd(p);
    qlgcd=lgcd(q);
    pp=p/plgcd; 
    qq=q/qlgcd;
    gcdlgcd=gcd(plgcd,qlgcd);
    gen gcdfirstcoeff(gcd(pp.coord.front().value, qq.coord.front().value));
    int gcddeg= min(pp.lexsorted_degree(),qq.lexsorted_degree());
    gen bound(pow(gen(2),gcddeg+1)* abs(gcdfirstcoeff,context0) * min(pp.norm(), qq.norm(),context0));
    gen modulo(nextprime(max(gcdfirstcoeff+1,gen(30000),context0))); 
    gen productmodulo(1);
    polynome currentgcd(p.dim),p_simp(p.dim),q_simp(p.dim),rem(p.dim);
    // 30000 leaves many primes below the 2^15 bound 
    for (;;modulo = nextprime(modulo+2)){
      // increment modulo to avoid modulo = 1 [4] so that it works in Z[i]
      while ( is_one(modulo % 4) || is_zero(gcdfirstcoeff % modulo))
	modulo=nextprime(modulo+2);
      polynome _gcdmod(gcdmod(smod(pp,modulo),smod(qq,modulo),modulo));
      gen adjustcoeff=gcdfirstcoeff*invmod(_gcdmod.coord.front().value,modulo);
      _gcdmod=smod((_gcdmod * adjustcoeff), modulo) ;
      int m=_gcdmod.lexsorted_degree();
      if (!m){
	p=pp*(plgcd/gcdlgcd);
	q=qq*(qlgcd/gcdlgcd);
	d=gcdlgcd;
	return true;
      }
      // combine step
      if (m<gcddeg){ // previous prime was bad
	gcddeg=m;
	currentgcd=_gcdmod;
	productmodulo=modulo;
      }
      else {
	if (m==gcddeg){ // start combine
	  if (productmodulo==gen(1)){ // no need to combine primes
	    currentgcd=_gcdmod;
	    productmodulo=modulo;
	  }
	  else {
	    //  cout << "Modulo:" << modulo << " " << _gcdmod << endl;
	    // cout << "Old gcd:" << productmodulo << " " << currentgcd << endl ;
	    currentgcd=ichinrem(_gcdmod,currentgcd,modulo,productmodulo);
	    // cout << "Combined to " << currentgcd << endl;
	    productmodulo=productmodulo*modulo;
	  }
	}
	// m>gcddeg this prime is bad, just ignore
      }
      //      if (productmodulo>bound){
      d=smod(currentgcd,productmodulo);
      ppz(d);
      //if ( pp.TDivRem1(d,p_simp,rem) && rem.coord.empty() && qq.TDivRem1(d,q_simp,rem) && rem.coord.empty() ){
      if ( divrem1(pp,d,p_simp,rem) && rem.coord.empty() && divrem1(qq,d,q_simp,rem) && rem.coord.empty() ){
	p=p_simp*(plgcd/gcdlgcd);
	q=q_simp*(qlgcd/gcdlgcd);
	d=d*gcdlgcd;
	return true;
      }
      // }
    }
    return false;
  }

  polynome pzadic(const polynome &p,const gen & n){
    monomial_v v;
    index_t i;
    for (monomial_v::const_iterator it=p.coord.begin();it!=p.coord.end();++it){
      i.clear();
      i.push_back(0);
      for (index_t::const_iterator iti=it->index.iptr->begin();iti!=it->index.iptr->end();++iti)
	i.push_back(*iti);
      gen k=it->value;
      for (int j=0;!is_zero(k);j++){
	gen r=smod(k,n.re(0));
	if (!is_zero(r)){
	  i[0]=j;
	  v.push_back(monomial<gen>(r,i));
	}
	k=iquo( (k-r),n.re(context0));
      }
    }
    // sort v
    polynome res(p.dim+1,v);
    res.tsort();
    return res;
  }

  bool listmax(const polynome &p,gen & n ){
    return Tlistmax<gen>(p,n);
  }

  polynome unmodularize(const polynome & p){
    polynome res(p.dim);
    vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    res.coord.reserve(itend-it);
    for (;it!=itend;++it){
      if (it->value.type==_MOD)
	res.coord.push_back(monomial<gen>(*it->value._MODptr,it->index));
      else
	res.coord.push_back(monomial<gen>(it->value,it->index));
    }
    return res;
  }

  void modularize(polynome & d,const gen & m){
    vector< monomial<gen> >::iterator it=d.coord.begin(),itend=d.coord.end();
    for (;it!=itend;++it)
      it->value=makemod(it->value,m);
  }

  // Find indexes of p such that p is constant, answer is in i
  void has_constant_variables(const polynome & p,index_t & i){
    i=vector<int>(p.dim,0);
    for (int j=0;j<p.dim;++j){
      i[j]=j;
    }
    vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    vector<int>::iterator iit,iitend;
    for (;it!=itend && !i.empty();++it){
      index_t & j=*it->index.iptr;
      iit=i.begin(); iitend=i.end();
      for (;iit!=iitend;++iit){
	if (j[*iit]){ // non-0 power in monomial
	  i.erase(iit);
	  iit=i.begin();
	  --iit;
	  iitend=i.end();
	}
      }
    }
  }

  // p assumed to be constant wrt variables in pi
  // vi is a vector of degree
  int extract_monomials(const polynome &p,const index_t & pi,vectpoly & vp){
    index_t pdeg=p.degree();
    // find largest degree of p with respect to these variables
    int s=pi.size(),ans=1;
    vector<int> v(s+1);
    int i=0;
    for (;i<s;++i){
      if (ans>1000) // FIXME what's the right size??
	return i;
      v[i]=pdeg[pi[i]]+1;
      ans=ans*v[i];
    }
    if (ans>10000)
      return -1;
    vp=vectpoly(ans,polynome(p.dim-s));
    if (ans==1)
      vp[0].coord.reserve(p.coord.size());
    vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    index_t::const_iterator piitbeg=pi.begin(),piit,piitend=pi.end(),vitbeg=v.begin(),vit;
    index_t::iterator iit,iti;
    int vp_pos;
    for (;it!=itend;++it){
      index_m i(p.dim-s);
      piit=piitbeg;
      vit=vitbeg;
      iit=i.iptr->begin();
      iti=it->index.iptr->begin();
      vp_pos=0;
      // construct new index without constant variables
      // and find value of index inside vp
      for (int j=0;j!=p.dim;++iti,++j){
	if (piit!=piitend && j==*piit){
	  ++piit;
	  vp_pos=vp_pos*(*vit)+(*iti);
	  ++vit;
	}
	else {
	  *iit=*iti;
	  ++iit;
	}
      }
      vp[vp_pos].coord.push_back(monomial<gen>(it->value,i));
    }
    return 0;
  }

  bool has_constant_variables_gcd(const polynome & p,const polynome & q,polynome & d){
    if (q.coord.empty()){
      d=p;
      return true;
    }
    if (p.coord.empty()){
      d=q;
      return true;
    }
    index_t pi,qi;
    has_constant_variables(p,pi);
    has_constant_variables(q,qi);
    // merge pi and qi
    index_t::iterator qit=qi.begin(),qitend=qi.end();
    for (;qit!=qitend;++qit){
      if (!equalposcomp(pi,*qit))
	pi.push_back(*qit);
    }
    if (pi.empty())
      return false;
    int s=pi.size();
    if (s==p.dim){
      gen n=gcd(Tcontent<gen>(p),Tcontent<gen>(q));
      d=polynome(monomial<gen>(n,p.dim));
      return true;
    }
    sort(pi.begin(),pi.end());
    // p or q is constant with respect to at least one variable
    // make a vector of polynomial from p and q
    vectpoly vp,vq;
    int i;
    if ( (i=extract_monomials(p,pi,vp)) ){
      pi=vector<int>(pi.begin(),pi.begin()+i);
      i=extract_monomials(p,pi,vp);
      if (i<0)
	return false;
    }
    if ( (i=extract_monomials(q,pi,vq)) ){
      pi=vector<int>(pi.begin(),pi.begin()+i);
      extract_monomials(p,pi,vp);
      i=extract_monomials(q,pi,vq);
      if (i<0)
	return false;
    }
    // find gcd of polys in vp and vq
    vectpoly::const_iterator it=vp.begin(),itend=vp.end(),jt=vq.begin(),jtend=vq.end();
    d=*jt;
    for (++jt;!is_one(d) && it!=itend;++it)
      d=gcd(d,*it);
    for (;!is_one(d) && jt!=jtend;++jt)
      d=gcd(d,*jt);
    // reconstruct gcd of p and q
    vector< monomial<gen> >::iterator dt=d.coord.begin(),dtend=d.coord.end();
    index_t::const_iterator piitbeg=pi.begin(),piit,piitend=pi.end(),dtit;
    int j;
    for (;dt!=dtend;++dt){
      index_m newi;
      newi.iptr->reserve(p.dim);
      piit=piitbeg;
      dtit=dt->index.iptr->begin();
      for (j=0;j<p.dim;++j){
	if (piit!=piitend && j==*piit){
	  newi.iptr->push_back(0);
	  ++piit;
	}
	else {
	  newi.iptr->push_back(*dtit);
	  ++dtit;
	}
      }
      dt->index=newi;
    }
    d.dim=p.dim;
    return true;
  }

  int coefftype(const polynome & p,gen & coefft){
    vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    int t;
    for (;it!=itend;++it){
      coefft=it->value;
      t=coefft.type;
      if (t==_USER)
	return t;
      if (t==_MOD)
	return t;
      if (t==_EXT)
	return t;
    }
    return 0;
  }

  bool ingcd=false;
  bool gcdheu(const polynome &p,const polynome &q, polynome & p_simp, gen & np_simp, polynome & q_simp, gen & nq_simp, polynome & d, gen & d_content,bool skip_test,bool compute_cofactors){
    // cout << "Entering gcdheu " << p.dim << endl;
    bool saveingcd=ingcd;
    ingcd=true;
    if (p.dim>=inspectdim)
      cerr << "Gcdheu begin " << p.dim << " " << clock() << " " << p.degree() << " " << p.coord.size() << " " << q.degree() << " " << q.coord.size() << endl;
    p_simp=p;
    q_simp=q;
    if (p.dim>=inspectdim)
      cerr << "Gcdheu end copy" << clock() << endl;
    // check if one coeff is a _MOD or _USER
    gen coefft;
    int pt=coefftype(p_simp,coefft);
    // If p, q have modular coeff, use modular algo
    if (!pt)
      pt=coefftype(q_simp,coefft);
    d_content=1;
    if (pt==_MOD){
      gen m=*(coefft._MODptr+1);
      if (debug_infolevel)
	cerr << "gcdmod begin " << clock() << endl;
      polynome pmod(unmodularize(p_simp));
      polynome qmod(unmodularize(q_simp));
      polynome dmod=gcdmod(pmod,qmod,m),pmodd,qmodd,tmp;
      if (debug_infolevel)
	cerr << "gcdmod end " << clock() << endl;
      d=dmod;
      modularize(d,m);
      if (compute_cofactors){
	divremmod(pmod,dmod,m,pmodd,tmp);
	divremmod(qmod,dmod,m,qmodd,tmp);
	// cerr << dmod << ":;\n" << pmodd << ":;\n" << qmodd << endl;
	p_simp=pmodd;
	modularize(p_simp,m);
	q_simp=qmodd;
	modularize(q_simp,m);
      }
      ingcd=saveingcd;
      return true;
    }
    if (pt==_USER){
      coefft._USERptr->polygcd(p_simp,q_simp,d);
      if (compute_cofactors){
	p_simp=p_simp/d;
	q_simp=q_simp/d;
      }
      ingcd=saveingcd;
      return true;
    }
    np_simp=ppz(p_simp);
    nq_simp=ppz(q_simp);
    if (pt==_EXT)
      return false;
    if (p.dim>=inspectdim)
      cerr << "Gcdheu end ppz" << clock() << " " << np_simp << " " << nq_simp << endl;
    d_content=gcd(np_simp,nq_simp);
    if (Tis_constant(p_simp) || Tis_constant(q_simp)){
      if (debug_infolevel>=2)
	cerr << "//Gcdheu p constant!" << endl;
      d=polynome(plus_one,p.dim);
      ingcd=saveingcd;
      return true;
    }
    index_t p_deg=p_simp.degree(),q_deg=q_simp.degree();
    if (!(p_deg>q_deg)){
      polynome quo(p.dim);
      if (exactquotient(q_simp,p_simp,quo)){
	d=p_simp;
	q_simp=quo;
	p_simp=polynome(monomial<gen>(plus_one,0,p.dim));
	if (is_positive(-d.coord.front())){
	  d=-d; p_simp=-p_simp; q_simp=-q_simp;
	}
	if ( (debug_infolevel>=2) || (p.dim>=inspectdim) )
	  cerr << "// End exact " << p.dim << " " << clock() << " " <<d.coord.size() << endl;
	ingcd=saveingcd;
	return true;
      }
      if ( (debug_infolevel>=2) || (p.dim>=inspectdim) )
	cerr << "//Gcdheu exact division failed! " << clock() << endl;
      if (p_simp.coord.size()==1){
	index_t i=index_gcd(*p_simp.coord.front().index.iptr,q_simp.gcddeg());
	d=polynome(monomial<gen>(plus_one,i));
	if (i!=index_t(i.size())){
	  i=-i;
	  p_simp=p_simp.shift(i);
	  q_simp=q_simp.shift(i);
	}
	ingcd=saveingcd;
	return true;
      }
    }
    if (!(q_deg>p_deg) ) {
      polynome quo(p.dim);
      if (exactquotient(p_simp,q_simp,quo)){
	d=q_simp;
	p_simp=quo;
	q_simp=polynome(monomial<gen>(plus_one,0,p.dim));
	if (is_positive(-d.coord.front())){
	  d=-d; p_simp=-p_simp; q_simp=-q_simp;
	}
	if ( (debug_infolevel>=2) || (p.dim>=inspectdim) )
	  cerr << "//End exact " << p.dim << " " << clock() << " " << d.coord.size() << endl;
	ingcd=saveingcd;
	return true;
      }
      if ( (debug_infolevel>=2) || (p.dim>=inspectdim) )
	cerr << "//Gcdheu exact division failed! " << clock() << endl;
      if (q_simp.coord.size()==1){
	index_t i=index_gcd(*q_simp.coord.front().index.iptr,p_simp.gcddeg());
	d=polynome(monomial<gen>(plus_one,i));
	if (i!=index_t(i.size())){
	  i=-i;
	  p_simp=p_simp.shift(i);
	  q_simp=q_simp.shift(i);
	}
	ingcd=saveingcd;
	return true;
      }
    } 
    if (p_simp.lexsorted_degree()==0){
      if (p.dim>=inspectdim)
	cerr << "Begin cst " << p.dim << " " << clock() << " " << d.coord.size() << endl;
      if (q_simp.lexsorted_degree()==0){
	d=gcd(p_simp.trunc1(),q_simp.trunc1()).untrunc1();
      }
      else {
	d=p_simp;
	Tlgcd<gen>(q_simp,d);
      }
      if (!is_one(d) && compute_cofactors){
	p_simp=p_simp/d;
	q_simp=q_simp/d;
      }
      if (p.dim>=inspectdim)
	cerr << "End cst " << p.dim << " " << clock() << " " << d.coord.size() << endl;
      ingcd=saveingcd;
      return true;
    }
    if (q_simp.lexsorted_degree()==0){
      if (p.dim>=inspectdim)
	cerr << "Begin cst " << p.dim << " " << clock() << " " << d.coord.size() << endl;
      d=q_simp;
      Tlgcd<gen>(p_simp,d);
      if (!is_one(d) && compute_cofactors){
	q_simp=q_simp/d;
	p_simp=p_simp/d;
      }
      if (p.dim>=inspectdim)
	cerr << "End cst " << p.dim << " " << clock() << " " << d.coord.size() << endl;
      ingcd=saveingcd;
      return true;
    }
    int Dbdeg=min(p.lexsorted_degree(),q.lexsorted_degree());
    bool est_reel=poly_is_real(p_simp) && poly_is_real(q_simp);    // FIXME: should check for extensions!
    if (debug_infolevel>=2)
      cerr << "//Gcdheu " << p.degree() << " " << p.coord.size() << " " << q.degree() << " " << q.coord.size() << endl;
   // first try evaluation for quick trivial gcd
    if (!skip_test ){
      if (p.dim>1) {
	vecteur b(p.dim-1);
	polynome Fb(1),Gb(1),Db(1);
	if ( (debug_infolevel>=2) || (p.dim>=inspectdim) )
	  cerr << "// GCD eval dimension " << p_simp.dim << " " << clock() << " " << p_simp.degree() << " " << p_simp.coord.size() << " " << q_simp.degree() << q_simp.coord.size() << " " << endl;
	gen essaimod=30011; // mod 4 = 3
	for (int essai=0;essai<2;++essai){
	  if (essai)
	    b=vranm(p.dim-1,0,0); // find another random point
	  // essaimod was est_reel?essaimod:0
	  for (;!find_good_eval(p_simp,q_simp,Fb,Gb,b,(debug_infolevel>=2) || (p.dim>=inspectdim),essaimod);){
	    for (;;){
	      essaimod=nextprime(essaimod+1);
	      if (!is_one(smod(essaimod,4)))
		break;
	    }
	  }
	  try {
	    Db=gcdmod(Fb,Gb,essaimod);
	  } catch (std::runtime_error & e){
	    Db=gcd(Fb,Gb);
	  }
	  Dbdeg=Db.lexsorted_degree();
	  if ( (debug_infolevel>=2) || (p.dim>=inspectdim) )
	    cerr << "// evaled GCD deg " << Dbdeg << endl;
	  if (!Dbdeg){
	    d.coord.clear();
	    Tcommonlgcd<gen>(p_simp,q_simp,d);
	    if ( (debug_infolevel>=2) || (p.dim>=inspectdim) )
	      cerr << "end eval " << p_simp.dim << " " << clock() << " " << d.coord.size() << endl;
	    if (compute_cofactors){
	      p_simp=p_simp/d;
	      q_simp=q_simp/d;
	    }
	    ingcd=saveingcd;
	    return true;
	  }
	  if (Dbdeg==p_simp.lexsorted_degree()){ // try p_simp/lgcd as gcd
	    if ( (debug_infolevel>=2) || (p.dim>=inspectdim) )
	      cerr << "Trying p/lgcd(p) as gcd " << p_simp.dim << " " << clock() << endl;
	    polynome p_simp_lgcd(Tlgcd<gen>(p_simp));
	    if ( (debug_infolevel>=2) || (p.dim>=inspectdim) )
	      cerr << "lgcd(p) ok " << p_simp.dim << " " << clock() << endl;
	    polynome p_simp_simp(p.dim);
	    if (!exactquotient(p_simp,p_simp_lgcd,p_simp_simp)) setsizeerr("gausspol.cc/gcdheu");
	    polynome quo(q_simp.dim);
	    if (exactquotient(q_simp,p_simp_simp,quo)){
	      if ( (debug_infolevel>=2) || (p.dim>=inspectdim) )
		cerr << "Success p/lgcd(p) as gcd " << p_simp.dim << " " << clock() << endl;
	      polynome quo_lgcd(p_simp_lgcd);
	      Tlgcd<gen>(quo,quo_lgcd);
	      d=p_simp_simp*quo_lgcd;
	      if (compute_cofactors){
		p_simp=p_simp_lgcd/quo_lgcd;
		q_simp=quo/quo_lgcd;
	      }
	      return true;
	    }
	    if ( (debug_infolevel>=2) || (p.dim>=inspectdim) )
	      cerr << "Failed p/lgcd(p) as gcd " << p_simp.dim << " " << clock() << endl;
	  }
	  if (Dbdeg==q_simp.lexsorted_degree()){ // try p_simp/lgcd as gcd
	    if ( (debug_infolevel>=2) || (p.dim>=inspectdim) )
	      cerr << "Trying q/lgcd(q) as gcd " << p_simp.dim << " " << clock() << endl;
	    polynome q_simp_lgcd(Tlgcd<gen>(q_simp));
	    if ( (debug_infolevel>=2) || (p.dim>=inspectdim) )
	      cerr << "lgcd(q) ok " << p_simp.dim << " " << clock() << endl;
	    polynome q_simp_simp(q_simp.dim);
	    if (!exactquotient(q_simp,q_simp_lgcd,q_simp_simp))  setsizeerr("gausspol.cc/gcdheu");
	    polynome quo(p_simp.dim);
	    if (exactquotient(p_simp,q_simp_simp,quo)){
	      if ( (debug_infolevel>=2) || (p.dim>=inspectdim) )
		cerr << "Success q/lgcd(q) as gcd " << p_simp.dim << " " << clock() << endl;
	      polynome quo_lgcd(q_simp_lgcd);
	      Tlgcd<gen>(quo,quo_lgcd);
	      d=q_simp_simp*quo_lgcd;
	      if (compute_cofactors){
		q_simp=q_simp_lgcd/quo_lgcd;
		p_simp=quo/quo_lgcd;
	      }
	      return true;
	    }
	    if ( (debug_infolevel>=2) || (p.dim>=inspectdim) )
	      cerr << "Failed q/lgcd(q) as gcd " << p_simp.dim << " " << clock() << endl;
	  }
	}
      }
    }
    // now work on p_simp and q_simp
    if (!p.dim){
      d=polynome(1,0);
      ingcd=saveingcd;
      return true;
    }
    gen np,nq,n;
    if (!listmax(p_simp,np)){
      ingcd=saveingcd;
      return false;
    }
    if (!listmax(q_simp,nq)){
      ingcd=saveingcd;
      return false;
    }
    if (p.dim==1){ // integer modular try, was p.dim==1 && est_reel
      environment * env= new environment;
      bool avoid_it=false;
      dense_POLY1 pp,qq;
      try {
	pp=modularize(p_simp,0,env);
	qq=modularize(q_simp,0,env);
      }
      catch (std::runtime_error & e){
	avoid_it=true;
      }
      env->moduloon = true;
      env->modulo=1001;
      env->pn=env->modulo;
      env->complexe=!est_reel;
      for (int essai=0;essai<2 && !avoid_it;++essai){
	env->modulo=nextprime(env->modulo+2);
	while ( is_one(smod(env->modulo,4)) || !is_one(gcd(gcd(env->modulo,pp.front()),qq.front())) )
	  env->modulo=nextprime(env->modulo+2);
	modpoly _gcdmod;
	gcdmodpoly(pp,qq,env,_gcdmod);
	Dbdeg=min(Dbdeg,_gcdmod.size()-1);
	if (!Dbdeg)
	  break;
      }
      delete env;
      if (!Dbdeg){
	d=polynome(1,p.dim);
	ingcd=saveingcd;
	return true;
      }	
    } 
    polynome p1(p.dim),q1(p.dim),r1(p.dim),r2(p.dim);
    gen n_2(2),n_73794(73794),n_27011(27011);
    if (is_greater(nq,np,context0))
      n=n_2*nq+n_2;
    else
      n=n_2*np+n_2;
    // PSR if gcd has a large degree, modular if low degree, else try heugcd
    // PSR complexity is proportionnal to 
    // #iteration*deg_var_n*(total_deg_other_vars*#iteration)^(2*#other_var)
    // MODGCD to product of all (part_deg_of_gcd+1+part_deg_of_gcd_lcoeff)
    int maxpqdeg0=max(p.lexsorted_degree(),q.lexsorted_degree());
    int minpqdeg0=min(p.lexsorted_degree(),q.lexsorted_degree());
    vector<int> pdeg(p.degree()),qdeg(q.degree()),maxpqdeg(p.dim);
    for (int i=0;i<p.dim;++i)
      maxpqdeg[i]=max(pdeg[i],qdeg[i]);
    double heugcddigits=1.0,maxmodop=1.0,minmodop=1.0;
    int total_deg_other_var=-maxpqdeg[0];
    for (int i=0;i<p.dim;++i){
      heugcddigits *= maxpqdeg[i] ; 
      minmodop *= Dbdeg+1; // approximation!!! should be partial degree[gcd]
      maxmodop *= min(pdeg[i],qdeg[i])+1+pdeg[i]+qdeg[i];
      total_deg_other_var += maxpqdeg[i];
    }
    double psrstep=minpqdeg0-Dbdeg;
    double psrgcdop=std::pow(psrstep,2*(p.dim-1))*(p.coord.size()*q.coord.size())/minpqdeg0, modop=std::sqrt(minmodop*maxmodop),heuop=heugcddigits*heugcddigits/10.0;
    double minop=psrgcdop; if (psrgcdop>heuop) minop=heuop; if (minop>modop) minop=modop;
    if (debug_infolevel)
      cerr << "Psr " << psrgcdop << ", Mod " << modop << ", Heu " << heuop << ", Min" << minop << endl;
    if (modop<minop) minop=modop; // was if (est_reel && modop<minop)
    if (debug_infolevel){
      if (p.dim==1)
	cerr << "GCD dim 1, n=" << n << " maxpqdeg0 " << maxpqdeg0 << "(" << maxpqdeg << ")" << endl;
      else
	cerr << "GCD dim " << p.dim << " degree " << Dbdeg << " psrgcdop " << psrgcdop << " heuop " << heuop << " modgcdop " << minmodop << "," << maxmodop << endl;
    }
    if (!skip_test){
      // int dd=p.dim*p.lexsorted_degree();
      // first try ezgcd then modgcd
      if ( // false && // uncomment to cancel EZGCD 
	  (p.dim>3) // && (Dbdeg<=maxpqdeg0/4+1) 
	   && ezgcd(p_simp,q_simp,d,true,true,0,minop)){
	if (debug_infolevel)
	  cout << "// Used EZ gcd " << endl;
	if (compute_cofactors){
	  q_simp=q_simp/d;
	  p_simp=p_simp/d;
	}
	ingcd=saveingcd;
	return true;
      }
      if ( p.dim>1 && psrgcdop< modop && psrgcdop < heuop ){
	d=gcdpsr(p_simp,q_simp,Dbdeg);
	ingcd=saveingcd;
	if (compute_cofactors){
	  q_simp=q_simp/d;
	  p_simp=p_simp/d;
	}
	return true;
      }
      if (//true  
	  modop < heuop 
	  ){ // was  if ( modop < heuop && est_reel)
	if (debug_infolevel)
	  cout << "// " << clock() << " Using modular gcd " << endl;
	ingcd=saveingcd;
	bool res=gcd_modular_algo(p_simp,q_simp,d,compute_cofactors);
	if (debug_infolevel)
	  cout << "// " << clock() << " End modular gcd " << endl;
	return res;
      }
    }
    if (debug_infolevel)
      cout << "// Using Heu gcd " << endl;
    int max_try=0;
    for (; max_try<GCDHEU_MAXTRY;max_try++){
      polynome pn(p_simp(n));
      // cout << p_simp << " pn:" << pn << " n=" << n << endl;
      polynome qn(q_simp(n));
      // NEW CODE
      if (skip_test){
	polynome pntmp,qntmp;
	gen pnsimp,qnsimp,dtmp;
	if (!gcdheu(pn,qn,pntmp,pnsimp,qntmp,qnsimp,d,dtmp,true,compute_cofactors)){
	  ingcd=saveingcd;
	  return false;
	}
	d=pzadic(d*dtmp,n);
      }
      else {
      // OLD CODE
	d=gcd(pn,qn);
	d=pzadic(d,n);
      }
      // END MODIFICATIONS
      ppz(d);
      if (!d.coord.empty()){
	gen tmp=d.coord.front().value;
	if (is_zero(re(tmp,context0))){ // checked
	  d=cst_i*d;
	  tmp=cst_i*tmp;
	}
	if (is_positive(-tmp,context0)) // checked
	  d=-d;
      }
      //      cout << "(ppz) d=" << d << endl;
      // cout << "p_simp" << p_simp << endl << "q_simp" << q_simp << endl;
      //if ( divrem(p_simp,d,p1,r1,false) && (r1.coord.size()==0) && divrem(q_simp,d,q1,r2,false) && (r2.coord.size()==0) ){
      if ( p_simp.TDivRem1(d,p1,r1) && (r1.coord.size()==0) && q_simp.TDivRem1(d,q1,r2) && (r2.coord.size()==0) ){
	//	cout << "p_simp/d" << p1 << "q_simp/d" << q1 << endl;
	p_simp=p1;
	q_simp=q1;
	// cout << "gcdheu success " << max_try << endl;
	if (p_simp.dim>=inspectdim)
	  cerr << "end gcdheu " << p_simp.dim << " " << clock() << " " << d.coord.size() << endl;
	ingcd=saveingcd;
	return true;
      }
      n=iquo(n*n_73794,n_27011);
    }
    // cout << "gcdheu failure" << endl;
    ingcd=saveingcd;
    return false; 
  }

  polynome gcdpsr(const polynome &p,const polynome &q,int gcddeg){
    if (debug_infolevel)
      cout << "// Using PSR gcd " << endl;
    if (!gcddeg && p.dim>1){ // find probable degree
      vecteur b(p.dim-1);
      polynome Fb(1),Gb(1),Db(1);
      for (int essai=0;essai<2;++essai){
	if (essai)
	  b=vranm(p.dim-1,0,0); // find another random point
	find_good_eval(p,q,Fb,Gb,b,(debug_infolevel>=2) || (p.dim>=inspectdim));
	Db=gcd(Fb,Gb);
	int Dbdeg=Db.lexsorted_degree();
	if (!Dbdeg)
	  return gcd(Tlgcd(p),Tlgcd(q));
	if (!gcddeg)
	  gcddeg=Dbdeg;
	else
	  gcddeg=min(Dbdeg,gcddeg);
      }
    }
    return Tgcdpsr<gen>(p,q,gcddeg);
  }

  // evaluate all vars but the j-th to 0
  gen peval0(const polynome & p,int j,int & total_deg){
    if (!j){
      vecteur v(p.dim-1);
      total_deg=total_degree(Tfirstcoeff<gen>(p).degree());
      return peval(p,v,0);
    }
    vecteur res;
    total_deg=0;
    int s=0,smax=0,n=p.dim,total;
    vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    int i,k=0;
    index_t::const_iterator itit;
    bool add;
    for (;it!=itend;++it){
      itit=it->index.iptr->begin();
      add=true;
      for (i=0,total=0;i<n;++i,++itit){
	total += *itit;
	if (i==j)
	  k=*itit;
	else {
	  if (*itit) 
	    add=false;
	}
      }
      if (add) {
	if (k<s)
	  res[k] = res[k]+it->value;
	else {
	  for (;s<k;++s)
	    res.push_back(0);
	  res.push_back(it->value);
	  ++s;
	}
      }
      if (k>smax){
	total_deg=total-k;
	smax=k;
      }
      if (k==smax)
	total_deg=max(total_deg,total-k);
    }
    reverse(res.begin(),res.end());
    return trim(res,0);
  }

  bool exchange_variables(polynome & p,polynome & q,index_t & permutation){
    if (p.dim<2)
      return false;
    index_t pdeg=p.degree(),qdeg=q.degree();
    int pd=pdeg.front(),qd=qdeg.front(),res=min(pd,qd),pos=0,tmp; 
    // Find first lowest degree position
    vector<int> vpos(1,0);
    for (int j=1;j<p.dim;++j){
      if ( (tmp=min(pdeg[j],qdeg[j])) <res){
	res=tmp;
	vpos=vector<int>(1,j);
      }
      if (tmp==res)
	vpos.push_back(j);
    }
    int s=vpos.size();
    // Same lowest degree, eval p at 0...0 and compare
    // (for ezgcd to find good eval: peval at zero must be non 0
    // and the lcoeff must be as small as possible)
    pos=vpos[0]; 
    if (s>1){
      int plcoeff,qlcoeff;
      gen p0=peval0(p,pos,plcoeff);
      gen q0=peval0(q,pos,qlcoeff);
      for (int j=1;j<s;++j){
	int ptmplcoeff,qtmplcoeff; 
	gen p0tmp(peval0(p,vpos[j],ptmplcoeff)),q0tmp(peval0(q,vpos[j],qtmplcoeff));
	if ( is_zero(p0tmp) || is_zero(q0tmp) )
	  continue;
	if ( is_zero(p0) || is_zero(q0) || (ptmplcoeff<plcoeff) ){
	  p0=p0tmp;
	  q0=q0tmp;
	  plcoeff=ptmplcoeff;
	  pos=vpos[j];
	}
      }
    }
    if (!pos)
      return false;
    if (p.dim>=inspectdim)
      cerr << "Exchange " << clock() << " " << p.dim << " " << p.degree() << " " << p.coord.size() << " " << q.degree() << " " << q.coord.size() << endl;
    permutation=transposition(0,pos,p.dim);
    p.reorder(permutation);
    q.reorder(permutation);
    return true;
  }

  void lcmdeno(const polynome & p, gen & res){
    vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    for (;it!=itend;++it){
      const gen & tmp = it->value;
      if (tmp.type==_FRAC)
	res=lcm(res,tmp._FRACptr->den);
    }
  }

  void simplify(polynome & p,polynome & q,polynome & p_gcd){
    if (is_one(q)){
      p_gcd=q;
      return;
    }
    if (is_one(p)){
      p_gcd=p;
      return ;
    }
    if (q.coord.empty()){
      p_gcd=polynome(gen(1),p.dim);
      swap(p_gcd.coord,p.coord);
      return ;
    }
    if (p.coord.empty()){
      p_gcd=polynome(gen(1),p.dim);
      swap(p_gcd.coord,q.coord);
      return ;
    }
    if (!p.dim){
      gen p0=p.coord.front().value,q0=q.coord.front().value;
      gen tmp=simplify(p0,q0);
      p=polynome(p0,0);
      q=polynome(q0,0);
      p_gcd=polynome(tmp,0);
      return;
    }
    if (p==q){
      p_gcd=polynome(1,p.dim);
      swap(p.coord,p_gcd.coord);
      q=p;
      return;
    }
    if (has_constant_variables_gcd(p,q,p_gcd)){
      polynome temp(p.dim);
      exactquotient(p,p_gcd,temp);
      p=temp;
      exactquotient(q,p_gcd,temp);
      q=temp;
      return ;
    }
    p_gcd.coord.clear();
    //    cout << "Simplifying " << p << endl << q << endl;
    polynome p_orig(p);
    polynome q_orig(q);
    index_t permutation;
    bool exchanged=exchange_variables(p_orig,q_orig,permutation);
    gen d_content=1,np_simp=1,nq_simp=1;
    if (gcdheu(p_orig,q_orig,p,np_simp,q,nq_simp,p_gcd,d_content)){
      p=p*rdiv(np_simp,d_content);
      q=q*rdiv(nq_simp,d_content);
      if (exchanged){
	p.reorder(permutation);
	q.reorder(permutation);
	p_gcd.reorder(permutation);
	if (!p_gcd.coord.empty() && is_strictly_positive(-p_gcd.coord.front().value,context0)){
	  p_gcd=-p_gcd;
	  p=-p;
	  q=-q;
	}
      }
      p_gcd=p_gcd*d_content;
      return ;
    }
    p_gcd=gcdpsr(p_orig,q_orig);
    polynome tmprem(p_gcd.dim);
    p_orig.TDivRem1(p_gcd,p,tmprem,true);
    q_orig.TDivRem1(p_gcd,q,tmprem,true);
    // If alg. extensions are involved, p and q may now contain fractions
    gen tmpmult(plus_one);
    lcmdeno(p,tmpmult);
    lcmdeno(q,tmpmult);
    p=p*tmpmult;
    q=q*tmpmult;
    if (exchanged){
      p.reorder(permutation);
      q.reorder(permutation);
      p_gcd.reorder(permutation);
    }
    p_gcd=p_gcd/tmpmult;
    return ;
  }

  polynome simplify(polynome &p,polynome &q){
    polynome p_gcd(p.dim);
    simplify(p,q,p_gcd);
    return p_gcd;
  }

  polynome gcd(const polynome & p,const polynome & q){
    control_c();
    if (p.coord.empty())
      return(q);
    if (q.coord.empty())
      return(p);
    /* if (p==q)
       return p; */
    if (p.dim==0){
      gen d=gcd(p.constant_term(),q.constant_term());
      index_t i;
      return polynome( monomial<gen>( d,i ) ) ;
    }
    polynome d(p.dim);
    if (has_constant_variables_gcd(p,q,d))
      return d;
    d.coord.clear();
    polynome p_orig(p),q_orig(q),p_simp(p.dim),q_simp(p.dim);
    gen d_content,np_simp,nq_simp;
    index_t permutation;
    bool exchanged=exchange_variables(p_orig,q_orig,permutation);
    if (gcdheu(p_orig,q_orig,p_simp,np_simp,q_simp,nq_simp,d,d_content,false,false)){
      if (exchanged)
	d.reorder(permutation);
      return d*d_content;
    }
    d=gcdpsr(p_simp,q_simp);
    if (exchanged)
      d.reorder(permutation);
    // if integers only, should add here gcd using modgcd
    d=d*d_content;
    if (!d.coord.empty() && d.coord.front().value.type==_MOD)
      d=d*inv(d.coord.front().value,context0);
    return d;
  }

  void egcdlgcd(const polynome &p1, const polynome & p2, polynome & u,polynome & v,polynome & d){
    TegcdTlgcd(p1,p2,u,v,d);
  }

  void egcdpsr(const polynome &p1, const polynome & p2, polynome & u,polynome & v,polynome & d){
  Tegcdpsr(p1,p2,u,v,d);
}
  
  void egcd(const polynome &p1, const polynome & p2, polynome & u,polynome & v,polynome & d){
    if (p1.dim==1)
      egcdlgcd(p1,p2,u,v,d);
    else
      egcdpsr(p1,p2,u,v,d);
  }

  bool findabcdelta(const polynome & p,polynome & a,polynome &b,polynome & c,polynome & delta){
    if (p.lexsorted_degree()!=2)
      return false;
    monomial_v::const_iterator it=p.coord.begin(),itend=p.coord.end();
    a=Tnextcoeff<gen>(it,itend);
    if (it==itend){
      b=polynome(a.dim);
      c=polynome(a.dim);
      delta=polynome(a.dim);
      return true;
    }
    if (it->index.iptr->front()==1)
      b=Tnextcoeff<gen>(it,itend);
    else
      b=polynome(a.dim);
    if (it==itend)
      c=polynome(a.dim);
    else
      c=Tnextcoeff<gen>(it,itend);
    delta=b*b-a*c*gen(4);
    return (it==itend);
  }

  bool findde(const polynome & p,polynome & d,polynome &e){
    if (p.coord.empty()){
      d=p;
      e=p;
      return true;
    }
    int n=p.lexsorted_degree();
    if (n>1)
      return false;
    monomial_v::const_iterator it=p.coord.begin(),itend=p.coord.end();
    if (!n){
      e=Tnextcoeff<gen>(it,itend);
      d=polynome(e.dim);
      return(it==itend);
    }
    d=Tnextcoeff<gen>(it,itend);
    if (it==itend)
      e=polynome(d.dim);
    else
      e=Tnextcoeff<gen>(it,itend);  
    return (it==itend);
  }

  void addtov(const polynome & tmp,vectpoly & v,bool with_sqrt,bool complexmode){
    if (!with_sqrt || tmp.lexsorted_degree()!=2 || tmp.dim>1)
      v.push_back(tmp);
    else {
      vecteur w=polynome2poly1(tmp,1);
      gen a=w.front(),b=w[1],c=w[2];
      gen delta=4*a*c-b*b,deltaf;
      if ( !complexmode && has_evalf(delta,deltaf,1,context0) && is_positive(deltaf,context0)){
	v.push_back(tmp);
	return;
      }
      gen b_over_2=rdiv(b,plus_two);
      if (b_over_2.type!=_FRAC){
	delta=a*c-b_over_2*b_over_2;
	vecteur vv(makevecteur(plus_one,rdiv(algebraic_EXTension(makevecteur(plus_one,b_over_2),makevecteur(plus_one,zero,delta)),a)));
	v.push_back(poly12polynome(vv,1));
	vv=makevecteur(a,algebraic_EXTension(makevecteur(minus_one,b_over_2),makevecteur(plus_one,zero,delta)));
	v.push_back(poly12polynome(vv,1));
      }
      else {
	vecteur vv(makevecteur(plus_one,rdiv(algebraic_EXTension(makevecteur(plus_one,b),makevecteur(plus_one,zero,delta)),plus_two*a)));
	v.push_back(poly12polynome(vv,1));
	vv=makevecteur(a,rdiv(algebraic_EXTension(makevecteur(minus_one,b),makevecteur(plus_one,zero,delta)),plus_two));
	v.push_back(poly12polynome(vv,1));
      }
    }
  }
  
  // factorize a square-free univariate polynomial
  void sqfffactor(const polynome &p, vectpoly & v,bool with_sqrt,bool test_composite,bool complexmode){
    // test if p has a numeric coeff
    if (has_num_coeff(p)){
      vecteur w=polynome2poly1(p,1);
      w=proot(w); 
      const_iterateur it=w.begin(),itend=w.end();
      polynome res(1),res2(1);
      res.coord.push_back(monomial<gen>(1,index_t(1,1)));
      res2.coord.push_back(monomial<gen>(1,index_t(1,2)));
      for (;it!=itend;++it){
	polynome copie(1);
	gen impart=im(*it,context0);
	if (!complexmode && !is_zero(impart) && (it+1)!=itend ){
	  copie = res2;
	  gen repart=re(*it,context0);
	  copie.coord.push_back(monomial<gen>(-2*repart,index_t(1,1)));
	  copie.coord.push_back(monomial<gen>(repart*repart+impart*impart,index_t(1,0)));
	  ++it;
	}
	else {
	  copie = res;
	  if (!is_zero(*it))
	    copie.coord.push_back(monomial<gen>(-*it,index_t(1,0)));
	}
	v.push_back(copie);
      }
      return ;
    }
    // find the gcd of the degrees of *it
    int d;
    if (test_composite)
      d=p.gcddeg(0);
    else
      d=1;
    if (d==1){
      // find linear factors now!
      environment * env=new environment;
      polynome temp(1);
      int ithprime=1;
      int bound=linearfind(p,env,temp,v,ithprime);
      // if degree of temp<=3, we are finished since not irred -> one fact
      // has degree 1 (hence found previously)
      int tempdeg=temp.lexsorted_degree();
      if (tempdeg<bound){
	if (tempdeg)
	  addtov(temp,v,with_sqrt,complexmode);
      }
      else {
	// find other factors
	vectpoly w;
	polynome unitaryp(1), an(0);
	int signe=1;
	if (is_positive(-temp.coord.front()))
	  signe=-1;
	unitarize(temp,unitaryp,an);
	factorunivsqff(unitaryp,env,w,ithprime,0,MODFACTOR_PRIMES);
	vectpoly::const_iterator itw=w.begin(),itwend=w.end();
	for (;itw!=itwend;++itw){
	  polynome tmp(ununitarize(*itw,an));
	  addtov(tmp,v,with_sqrt,complexmode);
	}
	if (signe==-1)
	  v.back()=-v.back();
      }
      delete env;
    }
    else { // gcddeg!=1, take the largest divisor of d
      //if (p.coord.size()==2){
      gen dd(d);
      vector<nfactor> nv(trivial_n_factor(dd));
      if (dd==gen(1))
	d=nv[nv.size()-1].fact.to_int();
      else
	d=dd.to_int();
      //}
      // use x^d as new variable, divide every degree by d
      polynome q(p.dividedegrees(d));
      vectpoly w;
      // IMPROVE: if we factor allowing 2nd order poly roots
      // we could factor bisquare poly
      // BUT that requires converting the roots to internal form
      sqfffactor(q,w,false,true,complexmode);
      vectpoly::const_iterator itw=w.begin(),itwend=w.end();
      for (;itw!=itwend;++itw)
	sqfffactor(itw->multiplydegrees(d),v,with_sqrt,false,complexmode);
    }
  }

  factorization sqff(const polynome &p ){
    return Tsqff_char0<gen>(p);
  }

  void sqff_evident(const polynome & p,factorization & f,bool with_sqrt,bool complexmode){
    // first make p primitive
    polynome pp=p/lgcd(p);
    // then square-free factorization
    factorization sqff_f(sqff(pp));
    f.clear();
    if (p.dim!=1){
      f=sqff_f;
      return ;
    }
    factorization::iterator it=sqff_f.begin();
    factorization::iterator itend=sqff_f.end();
    for (;it!=itend;++it){
      vectpoly v;
      sqfffactor(it->fact,v,with_sqrt,true,complexmode);
      vectpoly::const_iterator itv=v.begin(),itvend=v.end();
      for (;itv!=itvend;++itv)
	f.push_back(facteur<polynome>(*itv,it->mult));
    }
  }

  /* Factorization of sqff unitary polynomial with variables in reverse order
     Return number of factors, -1 if not successfull
     Might be called if polynomial is not unitary, but there is no
     proof that unlimited tries succeed in this case */
  int unitaryfactor(polynome & unitaryp, vectpoly & f,bool with_sqrt,bool complexmode){
    int dd=unitaryp.degree(unitaryp.dim-1);
    if (!dd)
      return 0; // unitaryp is cst w.r.t. x
    if (dd==1){
      f.push_back(unitaryp);
      unitaryp=unitaryp/unitaryp;
      return 1;
    }
    if (unitaryp.dim==1){
      factorization ff;
      sqff_evident(unitaryp,ff,with_sqrt,complexmode);
      if (f.empty())
	f.reserve(ff.size());
      factorization::const_iterator ff_it=ff.begin(),ff_end=ff.end();
      for (;ff_it!=ff_end;++ff_it)
	f.push_back(ff_it->fact);
      return ff.size();
    }
    ppz(unitaryp); // remove content
    gen n_2(2),np,n_73794(73794),n_27011(27011);
    polynome quo(unitaryp.dim),rem(unitaryp.dim);
    if (!listmax(unitaryp,np))
      return 0;
    gen x0(n_2*np+n_2); 
    int ntry=0;
    while (unitaryp.lexsorted_degree()){
      ntry++;
      if (ntry>GCDHEU_MAXTRY)
	return 0;
      // find evaluation point such that evaluated poly is sqff 
      // pz is unitary w.r.t. last var hence has same degree and is primitive
      polynome pz(unitaryp(x0));
      while (gcd(pz.derivative(),pz).lexsorted_degree()){
	x0=x0+gen(1);
	pz=unitaryp(x0); 
      }
      // factorization of pz
      vectpoly fz;
      int nf=unitaryfactor(pz,fz,with_sqrt,complexmode);
      if (nf==-1)
	return nf;
      if (!nf)
	return f.size();
      if (nf==1) {
	f.push_back(unitaryp);
	unitaryp=polynome(monomial<gen>(gen(1),0,unitaryp.dim));
	return f.size();
      }
      // factorization fz into factorization f
      vectpoly::iterator f_it=fz.begin(),f_itend=fz.end();
      for (;f_it!=f_itend;++f_it){
	*f_it=pzadic(*f_it,x0);
	// try division, each factor found is necessarily irreducible
	if ( (unitaryp.TDivRem1(*f_it,quo,rem)) && (rem.coord.empty())){
	  unitaryp=quo;
	  f.push_back(*f_it);
	}
      }
      x0=iquo(x0*n_73794,n_27011); // for the next try, if necessary
    }
    // factorize the cst term
    vectpoly fz;
    polynome tmp(unitaryp.trunc1());
    int nf=unitaryfactor(tmp,fz,with_sqrt,complexmode);
    if (nf==-1)
      return nf;
    if (!nf)
      return f.size();
    if (nf==1){
      f.push_back(unitaryp);
      unitaryp=polynome(monomial<gen>(gen(1),0,unitaryp.dim));
      return f.size();
    }
    vectpoly::iterator f_it=fz.begin(),f_itend=fz.end();
    for (;f_it!=f_itend;++f_it)
      f.push_back(f_it->untrunc1());
    unitaryp=polynome(monomial<gen>(gen(1),0,unitaryp.dim));
    return f.size();
  }
  
  void unitarize(const polynome &pcur, polynome &unitaryp, polynome & an){
    an=firstcoeff(pcur).trunc1();
    if (is_one(an)){
      unitaryp=pcur;
      return;
    }
    monomial_v::const_iterator it=pcur.coord.begin();
    monomial_v::const_iterator itend=pcur.coord.end();
    polynome curanpow(pow(an,0));
    int savpow=it->index.iptr->front();
    unitaryp=pow(polynome(monomial<gen>(gen(1),1,pcur.dim)),savpow);
    savpow--;
    int newpow;
    ++it;
    for (;it!=itend;){
      newpow=it->index.iptr->front();
      polynome an_1=Tnextcoeff<gen>(it,itend);
      curanpow=curanpow*pow(an,savpow-newpow);
      unitaryp=unitaryp+(an_1*curanpow).untrunc1(newpow);
      savpow=newpow;
    }
  }

  polynome ununitarize(const polynome & unitaryp, const polynome & an){
    if (is_one(an))
      return unitaryp;
    monomial_v::const_iterator it=unitaryp.coord.begin();
    monomial_v::const_iterator itend=unitaryp.coord.end();
    int curpow;
    polynome ppush(unitaryp.dim);
    for (;it!=itend;){
      curpow=it->index.iptr->front();
      polynome an_1=Tnextcoeff<gen>(it,itend);
      ppush=ppush+(an_1*pow(an,curpow)).untrunc1(curpow);
    }
    return ppush/lgcd(ppush);
  }

  // factorization over an algebraic extension
  // the main variable of G is the algebraic extension variable
  // the minimal polynomial of this variable is p_mini
  // G is assumed to be square-free
  // See algorithm 3.6.4 in Henri Cohen book starting at step 3
  bool algfactor(const polynome & G,const polynome & p_mini,int & k,factorization & f,bool complexmode){
    // search sqff norm
    polynome norme(G.dim),temp(G.dim);
    k=-1;
    for (;;) {
      ++k;
      // replace X by X-k*Y in G and _compute resultant
      if (k){
	vecteur v;
	polynome2poly1(G,2,v); // X is the second var
	polynome decal(G.dim-1);
	decal.coord.push_back(monomial<gen>(gen(-k),1,G.dim-1)); // -k*main_var
	v=taylor(v,decal);
	poly12polynome(v,2,temp,G.dim);
	norme=resultant(temp,p_mini).trunc1();
      }
      else
	norme=resultant(G,p_mini).trunc1();
      // check that norme is squarefree, first find inner dimension
      polynome dnorme=norme.derivative();
      int innerdim=0;
      vector< monomial<gen> >::const_iterator ckalg_it=norme.coord.begin(),ckalg_itend=norme.coord.end();
      for (; ckalg_it!=ckalg_itend;++ckalg_it){
	if (ckalg_it->value.type==_POLY){
	  innerdim=ckalg_it->value._POLYptr->dim;
	  break;
	}
      }
      // convert to usual multivariate polynomials
      polynome N(unsplitmultivarpoly(norme,innerdim)),Np(unsplitmultivarpoly(norme.derivative(),innerdim));
      polynome GG=gcd(N,Np);
      if (!GG.lexsorted_degree())
	break;
    }
    bool test=factor(norme,temp,f,true,false,complexmode);
    return test;
  }

  bool cfactor(const polynome & p, gen & an,factorization & f,bool with_sqrt){
    an=p.coord.front().value;
    if (has_num_coeff(p)){
      vectpoly w;
      sqfffactor(p,w,false,false,true);
      vectpoly::const_iterator itw=w.begin(),itwend=w.end();
      for (;itw!=itwend;++itw)
	f.push_back(facteur<polynome>(*itw,1));
      return true;
    }
    factorization fsqff=sqff(p);
    // factorization of each factor of fsqff
    factorization fz;
    factorization::const_iterator it=fsqff.begin(),itend=fsqff.end();
    for (;it!=itend;++it){
      polynome pcur=it->fact;
      int mult=it->mult;
      int d=pcur.lexsorted_degree();
      if (!d)
	continue;
      if (d==1){
	an=rdiv(an,pow(pcur.coord.front().value,gen(mult),context0));
	f.push_back(facteur<polynome>(pcur,mult));
	continue;
      }
      // make a polynomial with 1 more variable (i)
      polynome p_y(im(pcur).untrunc1(1)+re(pcur).untrunc1());
      polynome p_mini(p_y.dim);
      p_mini.coord.push_back(monomial<gen>(1,1,p_y.dim));
      p_mini=p_mini.multiplydegrees(2);
      p_mini.coord.push_back(monomial<gen>(1,0,p_y.dim));
      int k;
      if (!algfactor(p_y,p_mini,k,fz,false))
	return false;
      factorization::const_iterator f_it=fz.begin(),f_itend=fz.end();
      for (;f_it!=f_itend;++f_it){
	if (k){ // shift f_it->fact
	  vecteur v;
	  polynome2poly1(f_it->fact,1,v);
	  gen decal=polynome(gen(0,k),f_it->fact.dim-1);
	  v=taylor(v,decal);
	  poly12polynome(v,1,pcur,f_it->fact.dim);
	  pcur=gcd(pcur,p);
	}
	else
	  pcur=gcd(f_it->fact,p);
	an=rdiv(an,pow(pcur.coord.front().value,gen(mult),context0));
	vectpoly tmpv;
	addtov(pcur,tmpv,with_sqrt,true);
	f.push_back(facteur<polynome>(tmpv[0],mult));
	if (tmpv.size()==2)
	  f.push_back(facteur<polynome>(tmpv[1],mult));
      }
    }
    return true;   
  }

  bool ext_factor(const polynome &p,const gen & e,gen & an,polynome & p_content,factorization & f,bool complexmode){
    if (e._EXTptr->type!=_VECT)
      settypeerr("Modular factorization not yet accessible");
    gen ip=im(p,context0);
    if (!is_zero(ip)){
      // replace i by [1,0]:[1,0,1]
       gen newp=re(p,context0)+algebraic_EXTension(makevecteur(1,0),makevecteur(1,0,1))*ip;
      if (newp.type!=_POLY)
	return false;
      gen bn=1,the_ext;
      lcmdeno(*newp._POLYptr,bn);
      newp=bn*newp;
      vector< monomial<gen> >::const_iterator it=newp._POLYptr->coord.begin(),itend=newp._POLYptr->coord.end();
      for (;it!=itend;++it){
	if (it->value.type==_EXT){
	  the_ext=it->value;
	  break;
	}
      }
      if (the_ext.type!=_EXT)
	return false;
      bool res=ext_factor(*newp._POLYptr,the_ext,an,p_content,f,false);
      an=an/bn;
      return res;
    }
    an=p.coord.front().value;
    factorization fsqff=sqff(p);
    // factorization of each factor of fsqff
    factorization fz;
    factorization::const_iterator it=fsqff.begin(),itend=fsqff.end();
    for (;it!=itend;++it){
      polynome pcur=it->fact;
      int mult=it->mult;
      int d=pcur.lexsorted_degree();
      if (!d)
	continue;
      if (d==1){
	an=rdiv(an,pow(pcur.coord.front().value,gen(mult),context0));
	f.push_back(facteur<polynome>(pcur,mult));
	continue;
      }
      // make a polynomial with 1 more variable: the extension
      vecteur v_mini;
      if ((e._EXTptr+1)->type==_VECT)
	v_mini=*((e._EXTptr+1)->_VECTptr);
      else
	settypeerr("To be implemented");
      // const_iterateur v_it,v_itend=v_mini.end();
      polynome p_y(p.dim+1);
      // polynome p_mini(poly12polynome(v_mini));
      // p_mini=p_mini.untrunc(0,2);
      // p_mini.reorder(transposition(0,1,2));
      // polynome p_mini(poly12polynome(v_mini,1,p.dim+1));
      polynome p_mini(p.dim+1);
      const_iterateur tmp_it=v_mini.begin(),tmp_itend=v_mini.end();
      for (int d=(tmp_itend-tmp_it)-1;tmp_it!=tmp_itend;++tmp_it,--d){
	if (is_zero(*tmp_it))
	  continue;
	p_mini.coord.push_back(monomial<gen>(*tmp_it,d,1,p_mini.dim));
      }
      vector< monomial<gen> >::const_iterator p_it=pcur.coord.begin(),p_itend=pcur.coord.end();
      for (;p_it!=p_itend;++p_it){
	if (p_it->value.type!=_EXT){
	  p_y.coord.push_back(p_it->untrunc1());
	  continue;
	}
	if (*(p_it->value._EXTptr+1)!=*(e._EXTptr+1))
	  setsizeerr("Factor: Only one algebraic extension allowed");
	// convert the polynomial of the algebraic extension generator
	index_t ii=*(p_it->index.iptr);
	ii.insert(ii.begin(),0);
	p_y=p_y+poly1_2_polynome(*(p_it->value._EXTptr->_VECTptr),p_y.dim).shift(ii);
      }
      int k;
      if (!algfactor(p_y,p_mini,k,fz,complexmode))
	return false;
      factorization::const_iterator f_it=fz.begin(),f_itend=fz.end();
      if (f_itend-f_it==1){ // irreducible (after sqff)
	an=rdiv(an,pow(pcur.coord.front().value,gen(mult),context0));
	f.push_back(facteur<polynome>(pcur,mult));
      }
      else {
	gen bn(1);
	for (;f_it!=f_itend;++f_it){
	  if (k){ // shift f_it->fact
	    //vecteur v=polynome2poly1(f_it->fact);
	    vecteur v; polynome2poly1(f_it->fact,1,v);
	    vecteur decalv(2,zero);
	    decalv[0]=k; 
	    gen decal=algebraic_EXTension(decalv,v_mini);
	    v=taylor(v,decal);
	    // pcur=poly12polynome(v); 
	    poly12polynome(v,1,pcur,f_it->fact.dim);
	    pcur=gcd(pcur,p);
	  }
	  else
	    pcur=gcd(f_it->fact,p);
	  // unitarize pcur instead of computing bn
	  pcur=pcur/pcur.coord.front().value;
	  // bn=bn*pow(pcur.coord.front().value,gen(mult));
	  f.push_back(facteur<polynome>(pcur,mult));
	}
	an=rdiv(an,bn);
      }
    } // end for (;it!=itend;)    
    return true;   
  }

  // build a multivariate poly
  // with normal coeff from a multivariate poly with multivariate poly coeffs
  polynome unsplitmultivarpoly(const polynome & p,int inner_dim){
    polynome res(p.dim+inner_dim);
    index_t inner_index,outer_index;
    vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    for (;it!=itend;++it){
      outer_index=*(it->index.iptr);
      if (it->value.type!=_POLY){
	for (int j=0;j<inner_dim;++j)
	  outer_index.push_back(0);
	res.coord.push_back(monomial<gen>(it->value,outer_index));
      }
      else {
	vector< monomial<gen> >::const_iterator jt=it->value._POLYptr->coord.begin(),jtend=it->value._POLYptr->coord.end();
	for (;jt!=jtend;++jt){
	  inner_index=*(jt->index.iptr);
	  res.coord.push_back(monomial<gen>(jt->value,mergeindex(outer_index,inner_index)));
	}
      }
    }
    return res;
  }
  
  // build from a multivariate poly with normal coeff 
  // a multivariate poly with multivariate poly coeffs
  polynome splitmultivarpoly(const polynome & p,int inner_dim){
    int outer_dim=p.dim-inner_dim;
    index_t cur_outer;
    polynome cur_inner(inner_dim);
    polynome res(outer_dim);
    vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    for (; it!=itend;++it){
      index_t cur_index=*(it->index.iptr);
      index_t outer_index(cur_index.begin(),cur_index.begin()+outer_dim);
      index_t inner_index(cur_index.begin()+outer_dim,cur_index.end());
      if (outer_index!=cur_outer){
	if (!is_zero(cur_inner))
	  res.coord.push_back(monomial<gen>(cur_inner,cur_outer));
	cur_inner.coord.clear();
	cur_outer=outer_index;
      }
      cur_inner.coord.push_back(monomial<gen>(it->value,inner_index));
    }
    if (!is_zero(cur_inner))
      res.coord.push_back(monomial<gen>(cur_inner,cur_outer));
    return res;
  }

  // if one coeff of p is a polynomial, we must build a multivariate poly
  // with normal coeff from a multivariate poly with multivariate poly coeffs
  bool poly_factor(const polynome & p, int inner_dim,polynome & p_content,factorization & f,bool with_sqrt,bool complexmode){
    // convert p -> pp
    polynome pp(unsplitmultivarpoly(p,inner_dim)),pp_content(p.dim+inner_dim);
    // factorize pp
    if (!factor(pp,pp_content,f,false,with_sqrt,complexmode))
      return false;
    // convert back pp_content -> p_content and each term of f
    p_content=splitmultivarpoly(pp_content,inner_dim);
    factorization::iterator f_it=f.begin(),f_itend=f.end();
    for (;f_it!=f_itend;++f_it)
      f_it->fact=splitmultivarpoly(f_it->fact,inner_dim);
    return true;
  }

  // Yun algorithm in finite field of characteristic n
  // Must be called recursively since it will not detect powers multiple of n
  factorization partialsquarefree_fp(const polynome & p,unsigned n,polynome & c){
    factorization v;
    polynome y(p.derivative()),w(p);
    c=simplify(w,y);
    // If p=p_1*p_2^2*...*p_n^n, 
    // then c=gcd(p,p')=Pi_{i s.t. i%n!=0} p_i^{i-1} Pi_{i s.t. i%n==0} p_i^i
    // w=p/c=Pi_{i%n>=1} p_i, 
    // y=p'/c=Sum_{i%n>=1} ip_i'*pi_{j!=i, j%n>=1} p_j
    y=y-w.derivative(); 
    // y=Sum_{i%n>=2} (i-1)p_i'*pi_{j!=i,j%n!=0} p_j
    int k=1;
    while(!y.coord.empty()){
      // y=sum_{i%n >= k+1} (i-k) p_i' * pi_{j!=i, j>=k} p_j
      polynome g=simplify(w,y);
      if (!Tis_one(g))
	v.push_back(facteur< polynome >(g,k)); 
      // extract one time the factors of multiplicity k mod n
      c=c/w;
      // this push p_k, now w=pi_{i%n>=k+1} p_i and 
      // y=sum_{i%n>=k+1} (i-k) p_i' * pi_{j!=i, j%n>=k+1} p_j
      y=y-w.derivative();
      // y=sum_{i%n>=k+1} (i-(k+1)) p_i' * pi_{j!=i, j%n>=k+1} p_j
      k++;
    }
    if (!Tis_one(w))
      v.push_back(facteur< polynome >(w,k));
    return v; 
    // at the end c contains Pi_{i} p_i^{i-(i%n)}
  }
  
  factorization vector2factorization(const vectpoly & v){
    vectpoly::const_iterator it=v.begin(),itend=v.end();
    factorization res;
    for (int i=1;it!=itend;++it){
      if (Tis_one<gen>(*it))
	res.push_back(facteur<polynome>(*it,i));
    }
    return res;
  }

  // Yun algorithm in finite field of characteristic n
  // Requires factorization_compress after
  factorization uncompressed_squarefree_fp(const polynome & p,unsigned n,unsigned exposant){
    factorization res;
    if (Tis_one(p))
      return res;
    polynome c(p.dim);
    res=partialsquarefree_fp(p,n,c);
    if (Tis_one(c))
      return res;
    int lp=c.gcddeg(0);
    if (lp%n)
      setsizeerr("Mod factor bug!");
    polynome b(c.dividedegrees(n));
    if (exposant!=1){
      // replace all coeffs of b by coeff^(p^(n-1))
      // since in F_{p^n} we have a=(a^(p^(n-1)))^p which is not a^p 
      std::vector< monomial<gen> > ::iterator it=b.coord.begin(),itend=b.coord.end();
      int ntoexposant=pow(n,exposant-1).val;
      for (;it!=itend;++it){
	it->value=pow(it->value,ntoexposant);
      }
    }
    // Note that this is not correct, we must compute gcd of res and resn
    // that have the same residue modulo n
    // Example factor( (x+1)^3*(x-1)^4 %3 )
    // puts x-1 in res and (x^2-1)^3 in c
    // Hence we must call factorization_compress at the end
    factorization resn(uncompressed_squarefree_fp(b,n,exposant));
    factorization::const_iterator it=resn.begin(),itend=resn.end();
    for (;it!=itend;++it){
      res.push_back(facteur<polynome>(it->fact,it->mult*n));
    }
    return res;
  }

  // Compress factorization, required for sqff on finite field
  void factorization_compress(factorization & sqff_f){
    factorization sqfftmp(sqff_f);
    sqff_f.clear();
    vecteur vtmp;
    int pos;
    factorization::const_iterator it=sqfftmp.begin(),itend=sqfftmp.end();
    for (;it!=itend;++it){
      if ( (pos=equalposcomp(vtmp,it->fact)) ){
	sqff_f[pos-1].mult += it->mult;
      }
      else {
	vtmp.push_back(it->fact);
	sqff_f.push_back(*it);
      }
    }
  }

  void sqff_ffield_factor(const factorization & sqff_f,int n,environment * env,factorization & f){
    // Now factorize each factor
    factorization::const_iterator it=sqff_f.begin(),itend=sqff_f.end();
    for (;it!=itend;++it){
      const facteur<polynome> & fp=*it;
      // convert to vector 
      modpoly Qtry(modularize(env->moduloon?unmodularize(it->fact):it->fact,n,env));
      // and call sqff mod factor
      vector< facteur<modpoly> > wf;
      vector<modpoly> qmat;
      // qmatrix(Qtry,env,qmat,0);
      ddf(Qtry,qmat,env,wf);
      vector<modpoly> w;
      cantor_zassenhaus(wf,qmat,env,w);
      // put result in f
      vector<modpoly>::const_iterator jt=w.begin(),jtend=w.end();
      gen gtmp;
      for ( ;jt!=jtend;++jt){
	polynome tmp(unmodularize(*jt));
	gtmp=env->moduloon?makemod(tmp,n):tmp;
	if (gtmp.type==_POLY)
	  f.push_back(facteur<polynome>(*gtmp._POLYptr,it->mult));
      }
    }
  }

  factorization squarefree_fp(const polynome & p,unsigned n,unsigned exposant){
    factorization res(uncompressed_squarefree_fp(p,n,exposant));
    factorization_compress(res);
    return res;
  }

  // p is primitive wrt the main var
  bool mod_factor(const polynome & p_orig,polynome & p_content,int n,factorization & f){
    environment env;
    env.moduloon = true;
    env.modulo=n;
    env.pn=n;
    if (p_orig.dim!=1)
      return false; // FIXME implement multivariate fact 
    // Check that all coeff are mod
    polynome p(p_orig);
    vector< monomial<gen> >::iterator pit=p.coord.begin(),pitend=p.coord.end();
    for (;pit!=pitend;++pit){
      if (pit->value.type!=_MOD)
	pit->value=makemod(pit->value,n);
      gen & tmp = *(pit->value._MODptr+1);
      if (tmp.type!=_INT_ || tmp.val!=n)
	setsizeerr();
      gen & val = *(pit->value._MODptr);
      if (val.type==_CPLX)
	env.complexe=true;
    }
#ifdef HAVE_LIBNTL
#ifdef HAVE_LIBPTHREAD
    int locked=pthread_mutex_trylock(&ntl_mutex);
#endif // HAVE_LIBPTHREAD
    if (!locked){
      bool res=true;
      try {
	vecteur v;
	if (p.dim!=1)  setsizeerr("gausspol.cc/mod_factor");
	if (p.coord.empty())
	  setsizeerr();
	int deg=p.lexsorted_degree();
	int curpow=deg;
	v.reserve(deg+1);
	vector< monomial<gen> >::const_iterator ppit=p.coord.begin();
	vector< monomial<gen> >::const_iterator ppitend=p.coord.end();
	for (;ppit!=ppitend;++ppit){
	  int newpow=ppit->index.iptr->front();
	  for (;curpow>newpow;--curpow)
	    v.push_back(0);
	  if (ppit->value.type==_INT_)
	    v.push_back(ppit->value);
	  if (ppit->value.type==_MOD)
	    v.push_back(*ppit->value._MODptr);
	  --curpow;
	}
	for (;curpow>-1;--curpow)
	  v.push_back(0);      
	// FIXME NTL works on monic polynomials only!!
	gen v0=v.front();
	if (!is_one(v0)){
	  p_content = p_content*v0;
	  v0=invmod(v0,gen(n));
	  v = operator_times(v,v0,&env);
	}
	if (n==2){
	  NTL::GF2X ntlf(modpoly2GF2X(v));
	  NTL::vec_pair_GF2X_long fres(NTL::CanZass(ntlf,0));
	  int s=fres.length();
	  for (int i=0;i<s;i++){
	    modpoly res( GF2X2modpoly(fres[i].a));
	    f.push_back(facteur<polynome>(*makemod(unmodularize(res),2)._POLYptr,fres[i].b));
	  }
	}
	else {
	  NTL::ZZ_p::init(inttype2ZZ(n));
	  NTL::ZZ_pX ntlf(modpoly2ZZ_pX(v));
	  NTL::vec_pair_ZZ_pX_long fres(NTL::CanZass(ntlf,0));
	  int s=fres.length();
	  for (int i=0;i<s;i++){
	    modpoly res( ZZ_pX2modpoly(fres[i].a));
	    f.push_back(facteur<polynome>(*makemod(unmodularize(res),n)._POLYptr,fres[i].b));
	  }
	}
      } catch (std::runtime_error & e){
	res=false;
      }
#ifdef HAVE_LIBPTHREAD
      pthread_mutex_unlock(&ntl_mutex);
#endif
      return res;
    } // end !locked
#endif
    // sqff
    factorization sqff_f(squarefree_fp(p,n,1));
    sqff_ffield_factor(sqff_f,n,&env,f);
    return true;
  }

  bool do_factor(const polynome &p,polynome & p_content,factorization & f,bool isprimitive,bool with_sqrt,bool complexmode,const gen & divide_an_by){
    f.clear();
    if (p.coord.empty()){
      p_content=p;
      return true;
    }
    polynome p_primit(p.dim);
    if (!isprimitive){
      p_content=lgcd(p);
      if (is_strictly_positive(-p.coord.front().value,context0) && is_strictly_positive(p_content.coord.front().value,context0)) 
	p_content=-p_content;
      p_primit=p/p_content;
    }
    else
      p_primit=p;
    p_content = p_content/divide_an_by;
    if (is_one(p_primit))
      return true;
    if (!is_zero(im(divide_an_by,0)))
      complexmode=true;    
    if (!p_content.coord.empty()){
      if (!is_zero(im(p_content.coord.front().value,0)))
	complexmode=true;
      // check if one coeff is an alg. extension (only one is allowed)
      if (p_content.coord.front().value.type==_EXT){
	gen an;
	if (!ext_factor(p_primit,p_content.coord.front().value,an,p_content,f,complexmode))
	  return false;
	p_content=an*p_content;
	return true;
      }
    }
    vector< monomial<gen> >::const_iterator ckalg_it=p.coord.begin(),ckalg_itend=p.coord.end();
    for (; ckalg_it!=ckalg_itend;++ckalg_it){
      if (ckalg_it->value.type==_USER){
	ckalg_it->value._USERptr->polyfactor(p_primit,f);
	return true;
      }
      if (ckalg_it->value.type==_EXT){
	gen an;
	if (!ext_factor(p_primit,ckalg_it->value,an,p_content,f,complexmode))
	  return false;
	p_content=an*p_content;
	return true;
      }
    }
    // check if polynomial coeff are embedded inside p
    for (ckalg_it=p.coord.begin(); ckalg_it!=ckalg_itend;++ckalg_it){
      if (ckalg_it->value.type==_POLY)
	return poly_factor(p,ckalg_it->value._POLYptr->dim,p_content,f,with_sqrt,complexmode);
    }
    // check if p has modular coeff
    for (ckalg_it=p.coord.begin(); ckalg_it!=ckalg_itend;++ckalg_it){
      if (ckalg_it->value.type==_MOD){
	if ((ckalg_it->value._MODptr+1)->type!=_INT_)
	  return false;
	return mod_factor(p_primit,p_content,(ckalg_it->value._MODptr+1)->val,f);
      }
    }
    // check if one coefficient is complex
    if (complexmode || !is_zero(im(p))){
      gen an;
      bool res=cfactor(p_primit,an,f,with_sqrt);
      if (!res)
	return false;
      p_content=an*p_content;
      return true;
    }
    if (p.dim==1){
      // FIXME: if p_primit has num coeffs, we must check the leading coeff
      // and adjust p_content
      if (has_num_coeff(p_primit)){
	gen an=p_primit.coord.front().value;
	p_content=an*p_content;
	vector< monomial<gen> >::iterator it=p_primit.coord.begin(),itend=p_primit.coord.end();
	for (;it!=itend;++it)
	  it->value=evalf(it->value/an,1,context0);
      }
      sqff_evident(p_primit,f,with_sqrt,complexmode);
      return true;
    }
    // first square-free factorization
    factorization fsqff=sqff(p_primit);
    // factorization of each factor of fsqff
    factorization::const_iterator it=fsqff.begin(),itend=fsqff.end();
    for (;it!=itend;++it){
      polynome pcur=it->fact;
      int mult=it->mult;
      if (has_num_coeff(pcur)){
	f.push_back(facteur<polynome>(pcur,1));
	continue;
      }
      // try first 2 good evaluations in case pcur is irreducible
      vecteur b(pcur.dim-1);
      factorization v;
      polynome Fb(1),Gb(1);
      int essai;
      for (essai=0;essai<2;++essai){
	if (essai)
	  b=vranm(pcur.dim-1,0,0); // find another random point
	find_good_eval(pcur,pcur,Fb,Gb,b,(debug_infolevel>=2));
	factor(Fb,Gb,v,false,false,false);
	if ( (v.size()==1) && (v.front().mult==1) )
	  break;
      }
      if (essai<2){
	f.push_back(facteur<polynome>(pcur,mult));
	continue;
      }
      /* Now try heuristic factorization then call unitaryfactor
	 on each found factor */
      vectpoly fz;
      pcur.reverse();
      unitaryfactor(pcur,fz,false,false);
      pcur.reverse();
      vectpoly::iterator f_it=fz.begin(),f_itend=fz.end();
      for (;f_it!=f_itend;++f_it){
	f_it->reverse();
	// if an!=1, P(Y)=P(a_n*X) and divide by content
	f.push_back(facteur<polynome>(*f_it,mult));
      }
      if (!is_one(pcur)){
	/* now make polynomial unitary with respect to last var
	   P(x)=a_n*x^n+...+a_0, x=X/a_n,
	   P(x)=Q(X)=1/a_n^(n-1) * [ X^n+ a_{n-1}*a_n X^(n-1)+...+ a_0*a_n^{n-1}]
	*/
	fz.clear();
	polynome unitaryp(p.dim),an(p.dim-1);
	unitarize(pcur,unitaryp,an);
	// rewrite variables in inverted order
	unitaryp.reverse();
	// and call unitaryfactor
	if (unitaryfactor(unitaryp,fz,false,false)==-1)
	  return false;
	// rewrite back variables in initial order for each polynomial
	// and push back factorization
	f_it=fz.begin(),f_itend=fz.end();
	for (;f_it!=f_itend;++f_it){
	  f_it->reverse();
	  // if an!=1, P(Y)=P(a_n*X) and divide by content
	  f.push_back(facteur<polynome>(ununitarize(*f_it,an),mult));
	}
      }
    }
    return true;
  }

  bool operator < (const polynome & f,const polynome & g){
    unsigned fs=f.coord.size(),gs=g.coord.size();
    if (fs!=gs)
      return gs<gs;
    if (!fs)
      return gs;
    if (!gs)
      return false;
    vector< monomial<gen> > ::const_iterator it=f.coord.begin(),jt=g.coord.begin(),itend=f.coord.end();
    for (;it!=itend;++it,++jt){
      if (it->index!=jt->index)
	return !(jt->index <= it->index);
      if (it->value!=jt->value){
	gen a=evalf_double(it->value,1,context0),b=evalf_double(jt->value,1,context0);
	if (a.type==_DOUBLE_ && b.type==_DOUBLE_)
	  return a._DOUBLE_val<b._DOUBLE_val;
	return it->value.islesscomplexthan(jt->value);
      }
    }
    return false;
  }

  bool operator < (const facteur<polynome> & f,const facteur<polynome> & g){
    const polynome & fp=f.fact;
    const polynome & gp=g.fact;
    return fp<gp;
  }

  bool factor(const polynome &p,polynome & p_content,factorization & f,bool isprimitive,bool with_sqrt,bool complexmode,const gen & divide_an_by){
    bool res=do_factor(p,p_content,f,isprimitive,with_sqrt,complexmode,divide_an_by);
    // sort f
    sort(f.begin(),f.end());
    return res;
  }

  bool is_positive(const polynome & p){
    if (p.coord.empty())
      return true;
    return (is_positive(p.coord.front().value,context0));
  }

  void partfrac(const polynome & num, const polynome & den, const vector< facteur< polynome > > & v , vector < pf <gen> > & pfde_VECT, polynome & ipnum, polynome & ipden ){
    Tpartfrac(num,den,v,pfde_VECT,ipnum,ipden);
  }

  pf<gen> intreduce_pf(const pf<gen> & p_cst, vector< pf<gen> > & intde_VECT ){
    return Tintreduce_pf<gen>(p_cst,intde_VECT);
  }

  vecteur vector_of_polynome2vecteur(const vectpoly & v){
    vecteur res;
    vectpoly::const_iterator it=v.begin(),itend=v.end();
    res.reserve(itend-it);
    for (;it!=itend;++it)
      res.push_back(*it);
    return res;
  }

  vecteur sturm_seq(const polynome & p,polynome & cont){
    vectpoly v;
    Tsturm_seq<gen>(p,cont,v);
    return vector_of_polynome2vecteur(v);
  }

  /*
  // accumulate partial evaluation in polynomial (it,itend)
  // cur_index and nvar indicate the number of first identical eval.variables
  // vsize is the total number of eval.variables
  polynome peval(vector< monomial<gen> >::const_iterator & it,const vector< monomial<gen> >::const_iterator & itend,const vector< vecteur > & power_of_xi,index_t & cur_index,int nvar,int vsize,int var0){
    polynome res(var0);
    for (;;){
      if (it==itend)
	return res;
      index_t & it_t=*it->index.iptr;
      index_t::const_iterator it_tt=it_t.begin();
      index_t::const_iterator it_ttend=it_tt+nvar,cur_it=cur_index.begin();
      for (;it_tt!=it_ttend;++cur_it,++it_tt){
	if (*it_tt!=*cur_it)
	  return res;
      }
      // same index beginning
      if (nvar==vsize){
	res.coord.push_back(monomial<gen>(it->value,index_t(it_t.begin()+vsize,it_t.end())));
	++it;
	if (debug_infolevel)
	  cerr << "// " << itend-it << " monomials remain " << clock() << endl;
      }
      else { // go one level deeper
	cur_index.push_back(*(it->index.iptr->begin()+nvar));
	const gen & g=power_of_xi[nvar][cur_index.back()];
	if (debug_infolevel)
	  cerr << "// Enter level " << nvar+1 << " " << clock() << endl;
	if (g.type==_POLY)
	  res=res+(*g._POLYptr)*peval(it,itend,power_of_xi,cur_index,nvar+1,vsize,var0);
	else
	  res=res+g*peval(it,itend,power_of_xi,cur_index,nvar+1,vsize,var0);
	cur_index.pop_back();
	if (debug_infolevel)
	  cerr << "// Back to level " << nvar << " " << clock() << endl;
      }
    }
  }

  gen peval(const polynome & p,const vecteur & v){
    int pdim=p.dim,vsize=v.size(),var0=pdim-vsize;
    if (var0<0)
      setsizeerr("Too much substitution variables");
    polynome res(var0);
    if (p.coord.empty())
      return res;
    vecteur vnum,vden;
    gen vn,vd;
    vnum.reserve(vsize);
    vden.reserve(vsize);
    for (int i=0;i<vsize;++i){
      fxnd(v[i],vn,vd);
      vnum.push_back(vn);
      vden.push_back(vd);
    }
    vector< vecteur > power_of_xi;
    power_of_xi.reserve(pdim);
    index_t pdeg(p.degree());
    index_t deg(pdeg.begin(),pdeg.begin()+vsize);
    // compute thet table of powers
    gen global_deno(plus_one);
    for (int i=0;i<vsize;++i){
      if (debug_infolevel)
	cerr << "// Computing powers of " << i << "th var " << clock() << endl;
      // compute powers of ith component num and deno
      vecteur va(1,plus_one),vb(1,plus_one);
      gen vn(vnum[i]),vd(vden[i]),vnpow(plus_one),vdpow(plus_one);
      int degi=deg[i];
      for (int j=1;j<=degi;++j){
	vnpow=vnpow*vn;
	va.push_back(vnpow);
	vdpow=vdpow*vd;
	vb.push_back(vdpow);
      }
      // multiply in reverse order for the denominators
      for (int j=0;j<=degi;++j)
	va[j]=va[j]*vb[degi-j];
      global_deno=global_deno*vb.back();
      power_of_xi.push_back(va);
    }
    // we are now ready to evaluate the polynomial
    vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    index_t cur_index;
    return fraction(peval(it,itend,power_of_xi,cur_index,0,vsize,var0),global_deno);
  }
  */

  // a*b+c*d
  gen foisplus(const polynome & a,const polynome & b,const polynome & c,const polynome & d){
    if (a.dim>=inspectdim)
      cerr << "foisplus begin " << clock() << endl;
    index_t da=a.degree(),db=b.degree(),dc=c.degree(),dd=d.degree(),de(a.dim);
    ulonglong ans=0;
    for (int i=0;i<a.dim;++i){
      de[i]=max(da[i]+db[i]+1,dc[i]+dd[i]+1);
      ans = ans*unsigned(de[i]);
      if (ans/RAND_MAX>RAND_MAX)
	break;
    }
    if (ans<=RAND_MAX){
      polynome * res = new polynome(a.dim);
      vector< T_unsigned<gen,unsigned> > pa,pb,p,pc,pd;
      convert<gen,unsigned>(a,de,pa);
      convert<gen,unsigned>(b,de,pb);
      smallmult<gen,unsigned>(pa,pb,p,0,100);
      convert<gen,unsigned>(c,de,pc);
      convert<gen,unsigned>(d,de,pd);
      smallmult<gen,unsigned>(pc,pd,pa,0,100);
      smalladd<gen,unsigned>(p,pa,pb);
      convert<gen,unsigned>(pb,de,*res);
      if (a.dim>=inspectdim)
	cerr << "foisplus end " << clock() << endl;
      return res;
    }
    if (ans/RAND_MAX<RAND_MAX){
      polynome * res = new polynome(a.dim);
      vector< T_unsigned<gen,ulonglong> > pa,pb,p,pc,pd;
      convert<gen,ulonglong>(a,de,pa);
      convert<gen,ulonglong>(b,de,pb);
      smallmult<gen,ulonglong>(pa,pb,p,0,100);
      convert<gen,ulonglong>(c,de,pc);
      convert<gen,ulonglong>(d,de,pd);
      smallmult<gen,ulonglong>(pc,pd,pa,0,100);
      smalladd<gen,ulonglong>(p,pa,pb);
      convert<gen,ulonglong>(p,de,*res);
      return res;
    }
    return a*b+c*d;
  }

  gen pevaladd(const gen & aa,const gen & bb){
    if (!ingcd && debug_infolevel)
      cerr << "pevaladd begin " << clock() << endl;
    gen res=aa+bb;
    if (!ingcd && debug_infolevel)
      cerr << "pevaladd end " << clock() << endl;
    return res;
  }

  gen pevalmul(const gen & aa,const gen & bb,const gen & m){
    if (!ingcd && debug_infolevel)
      cerr << "pevalmul begin " << clock() << endl;
    gen res;
    if (!is_zero(m))
      res=smod(aa,m)*bb;
    else
      res=aa*bb;
    /*
    if ( (aa.type!=_FRAC) || (bb.type!=_FRAC) )
      return aa*bb;
    const Tfraction<gen> & a(*aa._FRACptr);
    const Tfraction<gen> & b(*bb._FRACptr);
    gen res(Tfraction<gen>(a.num*b.num,a.den*b.den));
    */
    if (!ingcd && debug_infolevel)
      cerr << "pevalmul end " << clock() << endl;
    return res;
  }

  // Horner like evaluation
  // m != 0 for modular evaluation
  gen peval(vector< monomial<gen> >::const_iterator & it,const vector< monomial<gen> >::const_iterator & itend,const vecteur & nums,const vecteur & dens,const index_t & deg,index_t & cur_index,int nvar,int vsize,int var0,const gen & m){
    if (it==itend)
      return zero;
    if (nvar==vsize){
      polynome res(var0);
      for (;;){
	if (it==itend)
	  return res;
	index_t & it_t=*it->index.iptr;
	index_t::const_iterator it_tt=it_t.begin();
	index_t::const_iterator it_ttend=it_tt+nvar,cur_it=cur_index.begin();
	for (;it_tt!=it_ttend;++cur_it,++it_tt){
	  if (*it_tt!=*cur_it)
	    return res;
	}
	// same main variables powers, accumulate constants
	res.coord.push_back(monomial<gen>(it->value,index_t(it_t.begin()+vsize,it_t.end())));
	++it;
      if (debug_infolevel && !ingcd)
	cerr << "// " << itend-it << " monomials remain " << clock() << endl;
      }
    } 
    // we are not at the deepest level
    gen res,tmp1,tmp2;
    int prev_power=0,cur_power=deg[nvar];
    const gen & gn=nums[nvar];
    const gen & gd=dens[nvar];
    gen cur_gd(plus_one);
    if (is_zero(gn)){ // if gn=0 we just discard monomials
      for (;;){
	if (it==itend)
	  return zero;
	index_t & it_t=*it->index.iptr;
	index_t::const_iterator it_tt=it_t.begin();
	index_t::const_iterator it_ttend=it_tt+nvar,cur_it=cur_index.begin();
	for (;it_tt!=it_ttend;++cur_it,++it_tt){
	  if (*it_tt!=*cur_it)
	    return zero;
	}
	if (!*it_tt) // break at first monomial with power = 0 at this index
	  break;
	++it;
      }
      cur_index.push_back(0);
      gen res(pow(gd,prev_power));
      if (!is_zero(m))
	res=smod(res,m);
      res=res*peval(it,itend,nums,dens,deg,cur_index,nvar+1,vsize,var0,m);
      cur_index.pop_back();
      return res;
    } // end gn==0
    for (;;){
      prev_power=cur_power;
      if (it==itend){
	if (!prev_power)
	  return res;
	else
	  return pevalmul(pow(gn,prev_power),res,m);
	// return pow(gn,prev_power)*res;
      }
      cur_power=*(it->index.iptr->begin()+nvar);
      // same powers for the beginning indices? (always true the first time)
      index_t & it_t=*it->index.iptr;
      index_t::const_iterator it_tt=it_t.begin();
      index_t::const_iterator it_ttend=it_tt+nvar,cur_it=cur_index.begin();
      for (;it_tt!=it_ttend;++cur_it,++it_tt){
	if (*it_tt!=*cur_it)
	  return pevalmul(pow(gn,prev_power),res,m);
	// return pow(gn,prev_power)*res;
      }
      // Yes: go one level deeper
      tmp1=pevalmul(pow(gn,prev_power-cur_power),res,m);
      res=zero;
      cur_index.push_back(cur_power);
      if (debug_infolevel && !ingcd)
	cerr << "// Enter level " << nvar+1 << " " << clock() << " ^ " << prev_power-cur_power << endl;
      tmp2=peval(it,itend,nums,dens,deg,cur_index,nvar+1,vsize,var0,m);
      cur_index.pop_back();
      if (debug_infolevel && !ingcd)
	cerr << "// Back to level " << nvar << " " << clock() << endl; 
      cur_gd=cur_gd*pow(gd,prev_power-cur_power);
      if (!is_zero(m))
	cur_gd=smod(cur_gd,m);
      // res=pevalmul(pow(gn,prev_power-cur_power),res)+cur_gd*peval(it,itend,nums,dens,deg,cur_index,nvar+1,vsize,var0);
      res=pevaladd(tmp1,cur_gd*tmp2);
      if (!is_zero(m))
	res=smod(res,m);
      tmp1=zero;
      tmp2=zero;
    }
  }
  
  gen peval(const polynome & p,const vecteur & v,const gen & m,bool simplify_at_end,vector<int_unsigned> * pptr){
    int pdim=p.dim,vsize=v.size(),var0=pdim-vsize;
    if (v==vecteur(vsize)){ // fast evaluation at 0
      index_t i(pdim);
      i[vsize-1]=1;
      // i=(0,0,...,0,1)
      // find the last position in p where a monomial with index i
      // could be inserted, the remaining of p truncated is the answer
      vector< monomial<gen> >::const_iterator it,itend=p.coord.end();
      it=upper_bound(p.coord.begin(),itend,monomial<gen>(plus_one,i),p.m_is_greater);
      if ( (it!=itend) && (*it->index.iptr==i))
	++it;
      polynome res(var0);
      res.coord.reserve(itend-it);
      for (;it!=itend;++it){
	index_t & itt=*it->index.iptr;
	res.coord.push_back(monomial<gen>(it->value,index_t(itt.begin()+vsize,itt.end())));
      }
      return res;
    }
    if (vsize==1 && m.type==_INT_ && m.val>0 && m.val<46340){
      polynome res;
      if (peval(p,v.front(),m.val,res,pptr))
	return res;
    }
    if (var0<0)
      setsizeerr("Too much substitution variables");
    polynome res(var0);
    if (p.coord.empty())
      return res;
    vecteur vnum,vden;
    gen vn,vd;
    vnum.reserve(vsize);
    vden.reserve(vsize);
    if (simplify_at_end){
      for (int i=0;i<vsize;++i){
	fxnd(v[i],vn,vd);
	vnum.push_back(vn);
	vden.push_back(vd);
      }
    }
    else {
      vnum=v;
      vden=vecteur(vsize,plus_one);
    }
    // we are now ready to evaluate the polynomial
    vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    index_t cur_index;
    index_t pdeg(p.degree());
    index_t deg(pdeg.begin(),pdeg.begin()+vsize);
    gen numer(peval(it,itend,vnum,vden,deg,cur_index,0,vsize,var0,m));
    if (!is_zero(m))
      numer=smod(numer,m);
    if (debug_infolevel && !ingcd){
      cerr << "// Peval end " << clock();
      if (numer.type==_POLY)
	cerr << " poly " << numer._POLYptr->coord.size();
      cerr << endl;
    }
    if ( is_zero(numer))
      return numer;
    // compute thet table of powers
    gen global_deno(plus_one);
    for (int i=0;i<vsize;++i){
      global_deno=global_deno*pow(vden[i],deg[i]);
    }
    simplify(numer,global_deno);
    return fraction(numer,global_deno);
  }

  int total_degree(const polynome & p){
    int res=0,deg;
    vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    for (;it!=itend;++it){
      deg=total_degree(*it->index.iptr);
      if (deg>res)
	res=deg;
    }
    return res;
  }


#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC
