/* -*- mode:C++ ; compile-command: "g++ -I.. -g -c threaded.cc" -*- */
/*  Multivariate GCD for large data not covered by the heuristic GCD algo
 *  Copyright (C) 2000 B. Parisse, Institut Fourier, 38402 St Martin d'Heres
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

#ifndef _GIAC_THREADED_H_
#define _GIAC_THREADED_H_
#include "first.h"
#ifndef WIN32

#if defined(__APPLE__) || defined(__FreeBSD__)
#else // was #ifndef __APPLE__
#include <sys/sysinfo.h>
#endif

#endif
#include "gausspol.h"
#include "ezgcd.h"
#ifdef HAVE_PTHREAD_H
#include <pthread.h>
#endif

#ifdef HAVE_GMPXX_H
#include <gmpxx.h>
#endif

#ifdef HAVE_GMPXX_H

  inline bool is_zero(const mpz_class & a){
    return a==0;
  }

  mpz_class invmod(const mpz_class & a,int reduce);

  mpz_class smod(const mpz_class & a,int reduce);


  inline void type_operator_times(const mpz_class & a,const mpz_class & b,mpz_class & c){
    // cerr << gen(c) << " = " << gen(a) << "*" << gen(b) << endl;
    mpz_mul(c.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t());
    // cerr << gen(c) << endl;
  }

  inline void type_operator_reduce(const mpz_class & a,const mpz_class & b,mpz_class & c,int reduce){
    mpz_mul(c.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t());
    if (reduce){
      c %= reduce;
    }
  }

  inline void type_operator_plus_times(const mpz_class & a,const mpz_class & b,mpz_class & c){
    // cerr << gen(c) << " += " << gen(a) << "*" << gen(b) << endl;
    // c+=a*b
    mpz_addmul(c.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t());
    // cerr << gen(c) << endl;
  }

  inline void type_operator_plus_times_reduce(const mpz_class & a,const mpz_class & b,mpz_class & c,int reduce){
    // c+=a*b % reduce;
    mpz_addmul(c.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t());
    if (reduce){
      c %= reduce;
    }
  }

#endif // HAVE__GMPXX_H
#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  // #define HEAP_MULT 1
  extern int heap_mult;

  inline void longlong2mpz(longlong i,mpz_t *mpzptr){
    int ii(i);
    if (i==ii)
      mpz_set_si(*mpzptr,ii);
    else {
      bool signe=(i<0);
      if (signe)
	i=-i;
      unsigned int i1=i>>32;
      unsigned int i2=i;
      mpz_init_set_ui(*mpzptr,i1);
      mpz_mul_2exp(*mpzptr,*mpzptr,32);
      mpz_add_ui(*mpzptr,*mpzptr,i2);
      if (signe)
	mpz_neg(*mpzptr,*mpzptr);
    }
  }

  // tmp is an allocated mpz_t 
  inline void mpz2longlong(mpz_t * ptr,mpz_t * tmp,longlong & ans){
    int i=mpz_sgn(*ptr);
    if (i<0)
      mpz_neg(*ptr,*ptr);
    mpz_tdiv_q_2exp(*tmp,*ptr,31);
    ans=mpz_get_ui(*tmp);
    ans <<= 31;
    mpz_tdiv_r_2exp(*tmp,*ptr,31);
    ans += mpz_get_ui(*tmp);
    if (i<0){
      mpz_neg(*ptr,*ptr);
      ans = -ans;
    }
  }

  class my_mpz {
  public:
    mpz_t ptr;
    my_mpz(){ mpz_init(ptr); }
    my_mpz(long l){ mpz_init_set_si(ptr,l); }
    my_mpz(const my_mpz & z){ mpz_init_set(ptr,z.ptr); }
    ~my_mpz() { mpz_clear(ptr); }
    my_mpz operator - () const { my_mpz tmp(*this); mpz_neg(tmp.ptr,tmp.ptr); return tmp;}
    my_mpz & operator = (const my_mpz & a) {mpz_set(ptr,a.ptr); return *this;}
  };

  my_mpz operator % (const my_mpz & a,const my_mpz & b);
  my_mpz operator + (const my_mpz & a,const my_mpz & b);
  my_mpz operator - (const my_mpz & a,const my_mpz & b);
  my_mpz operator * (const my_mpz & a,const my_mpz & b);
  my_mpz operator / (const my_mpz & a,const my_mpz & b);
  my_mpz invmod(const my_mpz & a,int reduce);
  my_mpz smod(const my_mpz & a,int reduce);
  my_mpz operator %= (my_mpz & a,const my_mpz & b);
  my_mpz operator += (my_mpz & a,const my_mpz & b);
  my_mpz operator -= (my_mpz & a,const my_mpz & b);

  template <class U>
  bool convert_myint(const polynome & p,const index_t & deg,std::vector< T_unsigned<my_mpz,U> >  & v){
    typename std::vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    v.clear();
    v.reserve(itend-it);
    U u;
    my_mpz tmp;
    index_t::const_iterator itit,ditbeg=deg.begin(),ditend=deg.end(),dit;
    T_unsigned<my_mpz,U> gu;
    for (;it!=itend;++it){
      u=0;
      itit=it->index.iptr->begin();
      for (dit=ditbeg;dit!=ditend;++itit,++dit)
	u=u*unsigned(*dit)+unsigned(*itit);
      gu.u=u;
      if (it->value.type==_ZINT)
	mpz_set(gu.g.ptr,*it->value._ZINTptr);
      else {
	if (it->value.type!=_INT_)
	  return false;
	mpz_set_si(gu.g.ptr,it->value.val);
      }
      v.push_back(gu);
    }
    return true;
  }


  inline bool is_zero(const my_mpz & a){
    return !mpz_sgn(a.ptr);
  }

  inline void type_operator_times(const my_mpz & a,const my_mpz & b,my_mpz & c){
    mpz_mul(c.ptr,a.ptr,b.ptr);
    // c=a*b;
  }

  inline void type_operator_reduce(const my_mpz & a,const my_mpz & b,my_mpz & c,int reduce){
    mpz_mul(c.ptr,a.ptr,b.ptr);
    if (reduce)
      mpz_fdiv_r_ui(c.ptr,c.ptr,reduce);
  }

  inline void type_operator_plus_times(const my_mpz & a,const my_mpz & b,my_mpz & c){
    // c+=a*b
    mpz_addmul(c.ptr,a.ptr,b.ptr);
  }

  inline void type_operator_plus_times_reduce(const my_mpz & a,const my_mpz & b,my_mpz & c,int reduce){
    // c+=a*b % reduce;
    mpz_addmul(c.ptr,a.ptr,b.ptr);
    if (reduce)
      mpz_fdiv_r_ui(c.ptr,c.ptr,reduce);
  }

#ifdef HAVE_GMPXX_H
  template <class U>
  bool convert_myint(const polynome & p,const index_t & deg,std::vector< T_unsigned<mpz_class,U> >  & v){
    typename std::vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    v.clear();
    v.reserve(itend-it);
    U u;
    index_t::const_iterator itit,ditbeg=deg.begin(),ditend=deg.end(),dit;
    for (;it!=itend;++it){
      u=0;
      itit=it->index.iptr->begin();
      for (dit=ditbeg;dit!=ditend;++itit,++dit)
	u=u*unsigned(*dit)+unsigned(*itit);
      T_unsigned<mpz_class,U> gu;
      gu.u=u;
      if (it->value.type==_ZINT){
	mpz_set(gu.g.get_mpz_t(),*it->value._ZINTptr);
      }
      else {
	if (it->value.type!=_INT_)
	  return false;
	gu.g=it->value.val;
      }
      v.push_back(gu);
    }
    return true;
  }
#endif
  
  template<class T,class U>
  struct pair_iterator {
    typename std::vector< T_unsigned<T,U> >::const_iterator it1,it2;
    U u;
    pair_iterator(typename std::vector< T_unsigned<T,U> >::const_iterator it,typename std::vector< T_unsigned<T,U> >::const_iterator jt):it1(it),it2(jt),u(it->u+jt->u) {}
    pair_iterator():u(0) {* (int *) (&it1)=0; * (int *) (&it2)=0;}
  };
  template<class T,class U>
  inline bool operator < (const pair_iterator<T,U> & p1,const pair_iterator<T,U> & p2){
    return p1.u<p2.u;
  }

  template<class T,class U>
  void smalladd(const std::vector< T_unsigned<T,U> > & v1,const std::vector< T_unsigned<T,U> > & v2,std::vector< T_unsigned<T,U> > & v){
    typename std::vector< T_unsigned<T,U> >::const_iterator it1=v1.begin(),it1end=v1.end(),it2=v2.begin(),it2end=v2.end();
    gen g;
    v.clear();
    v.reserve(it1end-it1+it2end-it2); // worst case
    for (;it1!=it1end && it2!=it2end;){
      if (it1->u==it2->u){
	g=it1->g+it2->g;
	if (!is_zero(g))
	  v.push_back(T_unsigned<T,U>(g,it1->u));
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

  template<class T,class U>
  void smallsub(const std::vector< T_unsigned<T,U> > & v1,const std::vector< T_unsigned<T,U> > & v2,std::vector< T_unsigned<T,U> > & v){
    typename std::vector< T_unsigned<T,U> >::const_iterator it1=v1.begin(),it1end=v1.end(),it2=v2.begin(),it2end=v2.end();
    gen g;
    v.clear();
    v.reserve(it1end-it1+it2end-it2); // worst case
    for (;it1!=it1end && it2!=it2end;){
      if (it1->u==it2->u){
	g=it1->g-it2->g;
	if (!is_zero(g))
	  v.push_back(T_unsigned<T,U>(g,it1->u));
	++it1;
	++it2;
      }
      else {
	if (it1->u>it2->u){
	  v.push_back(*it1);
	  ++it1;
	}
	else {
	  v.push_back(-*it2);
	  ++it2;
	}
      }
    }
    for (;it1!=it1end;++it1)
      v.push_back(*it1);
    for (;it2!=it2end;++it2)
      v.push_back(-*it2);
  }

  template<class T,class U>
  void smallmult(const T & g,const std::vector< T_unsigned<T,U> > & v1,std::vector< T_unsigned<T,U> > & v){
    typename std::vector< T_unsigned<T,U> >::const_iterator it1=v1.begin(),it1end=v1.end();
    if (&v1==&v){
      for (;it1!=it1end;++it1){
	it1->g = g*it1->g;
      }
    }
    else {
      v.clear();
      v.reserve(it1end-it1); // worst case
      for (;it1!=it1end;++it1){
	v.push_back(T_unsigned<T,U>(g*it1->g,it1->u));
      }
    }
  }

  template<class T,class U>
  void smalldiv(const std::vector< T_unsigned<T,U> > & v1,const T & g,std::vector< T_unsigned<T,U> > & v){
    typename std::vector< T_unsigned<T,U> >::const_iterator it1=v1.begin(),it1end=v1.end();
    if (&v1==&v){
      for (;it1!=it1end;++it1){
	it1->g = it1->g/g;
      }
    }
    else {
      v.clear();
      v.reserve(it1end-it1); // worst case
      for (;it1!=it1end;++it1){
	v.push_back(T_unsigned<T,U>(it1->g/g,it1->u));
      }
    }
  }


  // Possible improvement for threaded execution:
  // make each j of the k threads compute terms of the product with
  // degree wrt 1st var = j % k
  // at the end merge the results
  // For this, we might mark positions in v1 and v2 where the degree
  // wrt to the 1st var changes
  template<class T,class U>
  void smallmult(const std::vector< T_unsigned<T,U> > & v1,const std::vector< T_unsigned<T,U> > & v2,std::vector< T_unsigned<T,U> > & v,int reduce,size_t possible_size){
    typename std::vector< T_unsigned<T,U> >::const_iterator it1beg=v1.begin(),it1=v1.begin(),it1end=v1.end(),it2beg=v2.begin(),it2,it2end=v2.end();
    T g1,g2,g;
    U u1,u2,u;
    v.clear();
    if (heap_mult){
      unsigned v1s=it1end-it1beg,v2s=it2end-it2beg;
      std::vector< pair_iterator<T,U> > heap(v1s) ; // pointers to v2 monomials
      typename std::vector< pair_iterator<T,U> >::iterator heap0, heapbeg=heap.begin(),heapend=heap.begin()+v1s, heaplast=heap.begin()+v1s-1;
      // pair_iterator<T,U> * heap0, *heapbeg=heap,* heapend=heap+v1s, * heaplast=heap+v1s-1;
      for (it1=it1beg,heap0=heapbeg;heap0!=heapend;++heap0,++it1){
	*heap0=pair_iterator<T,U>(it1,it2beg);
      }
      for (;heapbeg!=heapend;){
	pop_heap(heapbeg,heapend); 
	if (!v.empty() && v.back().u==heaplast->u){
	  type_operator_plus_times_reduce(heaplast->it1->g,heaplast->it2->g,v.back().g,reduce); 
	  if (v.back().g == 0)
	    v.pop_back();
	}
	else {
	  type_operator_reduce(heaplast->it1->g,heaplast->it2->g,g,reduce); 
	  // g=0;
	  // type_operator_plus_times_reduce(heaplast->it1->g,heaplast->it2->g,g,reduce); 
	  v.push_back(T_unsigned<T,U>(g,heaplast->u));
	}
	++heaplast->it2;
	if (heaplast->it2==it2end){
	  --heaplast;
	  --heapend;
	}
	else {
	  heaplast->u = heaplast->it1->u + heaplast->it2->u;
	  push_heap(heapbeg,heapend);
	}
      }
    } // end if (heapmult)
    else {
#ifdef HASH_MAP_NAMESPACE
      typedef HASH_MAP_NAMESPACE::hash_map< U,T,hash_function_unsigned_object > hash_prod ;
      hash_prod produit(possible_size); // try to avoid reallocation
      // cout << "hash " << clock() << endl;
#else
      typedef std::map<U,T> hash_prod;
      // cout << "small map" << endl;
      hash_prod produit; 
#endif    
      typename hash_prod::iterator prod_it,prod_itend;
      for (;it1!=it1end;++it1){
	g1=it1->g;
	u1=it1->u;
	if (reduce){
	  for (it2=it2beg;it2!=it2end;++it2){
	    u=u1+it2->u;
	    prod_it=produit.find(u);
	    if (prod_it==produit.end())
	      type_operator_reduce(g1,it2->g,produit[u],reduce); // g=g1*it2->g; 
	    else 
	      type_operator_plus_times_reduce(g1,it2->g,prod_it->second,reduce); 
	  }
	}
	else {
	  for (it2=it2beg;it2!=it2end;++it2){
	    u=u1+it2->u;
	    prod_it=produit.find(u);
	    if (prod_it==produit.end()){
	      type_operator_times(g1,it2->g,produit[u]); // g=g1*it2->g; 
	    }
	    else {
	      type_operator_plus_times(g1,it2->g,prod_it->second); 
	      // g=g1*it2->g; 
	      // prod_it->second+=g;
	    }
	  }
	}
      }
      T_unsigned<T,U> gu;
      prod_it=produit.begin(),prod_itend=produit.end();
      v.reserve(produit.size());
      for (;prod_it!=prod_itend;++prod_it){
	if (!is_zero(gu.g=prod_it->second)){
	  gu.u=prod_it->first;
	  v.push_back(gu);
	}
      }    
      // cerr << "smallmult sort " << clock() << endl;
      sort(v.begin(),v.end());
      // cerr << "smallmult sort end " << clock() << endl;
    } // endif // HEAP_MULT
  }

  template<class T,class U>
  struct threadmult_t {
    const std::vector< T_unsigned<T,U> > * v1ptr ;
    std::vector< typename std::vector< T_unsigned<T,U> >::const_iterator > * v2ptr;
    std::vector< T_unsigned<T,U> > * vptr;
    U degdiv;
    unsigned current_deg;
    unsigned clock;
    int reduce;
  };

#ifdef HAVE_PTHREAD_H

  template<class T,class U> void * do_threadmult(void * ptr){
    threadmult_t<T,U> * argptr = (threadmult_t<T,U> *) ptr;
    argptr->clock=clock();
    int reduce=argptr->reduce;
    const std::vector< T_unsigned<T,U> > * v1 = argptr->v1ptr;
    std::vector< typename std::vector< T_unsigned<T,U> >::const_iterator >  * v2ptr = argptr->v2ptr;
    std::vector< T_unsigned<T,U> > & v = *argptr->vptr;
    typename std::vector< T_unsigned<T,U> >::const_iterator it1=v1->begin(),it1end=v1->end(),it2beg,it2,it2end;
    T g1,g;
    U u1,u2,u;
    int d1=0,d2,v2deg=v2ptr->size()-2,degdiv=argptr->degdiv,d=argptr->current_deg;
    // heap multiplication
#ifdef HASH_MAP_NAMESPACE
    typedef HASH_MAP_NAMESPACE::hash_map< U,T,hash_function_unsigned_object > hash_prod ;
    hash_prod produit; // try to avoid reallocation
    // cout << "hash " << clock() << endl;
#else
    typedef std::map<U,T> hash_prod;
    // cout << "small map" << endl;
    hash_prod produit; 
#endif    
    typename hash_prod::iterator prod_it,prod_itend;
    for (;it1!=it1end;++it1){
      u1=it1->u;
      d1=u1/degdiv;
      if (d1>d)
	continue;
      d2=v2deg-(d-d1);
      if (d2<0) // degree of d1 incompatible
	break;
      g1=it1->g;
      it2beg=(*v2ptr)[d2];
      it2end=(*v2ptr)[d2+1];
      if (reduce){
	for (it2=it2beg;it2!=it2end;++it2){
	  u2=it2->u;
	  u=u1+u2;
	  prod_it=produit.find(u);
	  if (prod_it==produit.end())
	    type_operator_reduce(g1,it2->g,produit[u],reduce); // g=g1*it2->g; 
	  else 
	    type_operator_plus_times_reduce(g1,it2->g,prod_it->second,reduce); 
	}
      }
      else {
	for (it2=it2beg;it2!=it2end;++it2){
	  u2=it2->u;
	  u=u1+u2;
	  prod_it=produit.find(u);
	  if (prod_it==produit.end()){
	    type_operator_times(g1,it2->g,g); // g=g1*it2->g; 
	    produit[u]=g;
	  }
	  else 
	    type_operator_plus_times(g1,it2->g,prod_it->second); 
	  // prod_it->second += g; 
	}
      }
    }
    T_unsigned<T,U> gu;
    prod_it=produit.begin(),prod_itend=produit.end();
    v.clear();
    v.reserve(produit.size());
    for (;prod_it!=prod_itend;++prod_it){
      if (!is_zero(gu.g=prod_it->second)){
	gu.u=prod_it->first;
	v.push_back(gu);
      }
    }    
    // cerr << "do_threadmult end " << clock() << endl;
    sort(v.begin(),v.end());
    // cerr << "do_threadmult sort end " << clock() << endl;
    argptr->clock = clock() - argptr->clock;
    return &v;
  }

  template<class T,class U>
  bool threadmult(const std::vector< T_unsigned<T,U> > & v1,const std::vector< T_unsigned<T,U> > & v2,std::vector< T_unsigned<T,U> > & v,U degdiv,int reduce,size_t possible_size=100){
    if (!threads_allowed)
      return false;
    if (v1.empty() || v2.empty())
      return true;
    unsigned threads_time=0;
#ifdef __APPLE__
    int nthreads=1;
#else
    int nthreads=sysconf (_SC_NPROCESSORS_ONLN);
#endif
    if (heap_mult && nthreads<2) 
      return false;
    unsigned d2=v2.front().u/degdiv,deg1v=v1.front().u/degdiv+d2;
    int cur_deg=-1,prev_deg=d2;
    // initialize iterators to the degree beginning
    std::vector< typename std::vector< T_unsigned<T,U> >::const_iterator > v2it;
    typename std::vector< T_unsigned<T,U> >::const_iterator it=v2.begin(),itend=v2.end();
    for (v2it.push_back(it);it!=itend;++it){
      cur_deg=it->u/degdiv;
      if (cur_deg==prev_deg)
	continue;
      for (int i=prev_deg-1;i>=cur_deg;--i){
	v2it.push_back(it);
	if (!i)
	  break;
      }
      prev_deg=cur_deg;
    }
    v2it.push_back(it);
    // degree of product wrt to the main variable
    // will launch deg1v+1 threads to compute each degree
    pthread_t tab[deg1v+1];
    threadmult_t<T,U> arg[deg1v+1];
    possible_size=0;
    int res=0;
    for (int I=deg1v;I>=0;I-=nthreads){
      for (int i=I;i>=0 && i>I-nthreads;--i){
	arg[i].v1ptr=&v1;
	arg[i].v2ptr=&v2it;
	arg[i].vptr = new std::vector< T_unsigned<T,U> >;
	arg[i].degdiv=degdiv;
	arg[i].current_deg=i;
	arg[i].reduce=reduce;
	if (nthreads==1)
	  do_threadmult<T,U>(&arg[i]);
	else
	  res=pthread_create(&tab[i],(pthread_attr_t *) NULL,do_threadmult<T,U>,(void *) &arg[i]);
	if (res){
	  // should cancel previous threads and delete created arg[i].vptr
	  return false;
	}
      }
      // now wait for these threads 
      for (int i=I;i>=0 && i>I-nthreads;--i){
	if (nthreads>1){
	  void * ptr;
	  pthread_join(tab[i],&ptr);
	}
	threads_time += arg[i].clock;
	possible_size += arg[i].vptr->size();
      }
    }
    // store to v
    v.reserve(possible_size);
    for (int i=deg1v;i>=0;--i){
      typename std::vector< T_unsigned<T,U> >::const_iterator it=arg[i].vptr->begin(),itend=arg[i].vptr->end();
      for (;it!=itend;++it){
	v.push_back(*it);
      }
      delete arg[i].vptr;
    }
    if (debug_infolevel)
      cerr << "Threads evaluation time " << double(threads_time)/CLOCKS_PER_SEC << endl;
    return true;
  }

#else // PTHREAD
  template<class T,class U>
  bool threadmult(const std::vector< T_unsigned<T,U> > & v1,const std::vector< T_unsigned<T,U> > & v2,std::vector< T_unsigned<T,U> > & v,U degdiv,int reduce,size_t possible_size=100){
    return false;
  }

#endif // PTHREAD


  template<class U> int coeff_type(const std::vector< T_unsigned<gen,U> > & p,unsigned & maxint){
    maxint=0;
    typename std::vector< T_unsigned<gen,U> >::const_iterator it=p.begin(),itend=p.end();
    if (it==itend)
      return -1;
    int t=it->g.type,tt;
    register int tmp;
    for (++it;it!=itend;++it){
      tt=it->g.type;
      if (tt!=t)
	return -1;
      if (!tt){
	if (it->g.val>0)
	  tmp=it->g.val;
	else
	  tmp=-it->g.val;
	if (maxint<tmp)
	  maxint=tmp;
      }
    }
    return t;
  }

  template <class U>
  bool convert_double(const polynome & p,const index_t & deg,std::vector< T_unsigned<double,U> >  & v){
    typename std::vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    v.clear();
    v.reserve(itend-it);
    T_unsigned<double,U> gu;
    U u;
    index_t::const_iterator itit,ditbeg=deg.begin(),ditend=deg.end(),dit;
    for (;it!=itend;++it){
      u=0;
      itit=it->index.iptr->begin();
      for (dit=ditbeg;dit!=ditend;++itit,++dit)
	u=u*unsigned(*dit)+unsigned(*itit);
      gu.u=u;
      if (it->value.type!=_DOUBLE_)
	return false;
      gu.g=it->value._DOUBLE_val;
      v.push_back(gu);
    }
    return true;
  }

  template <class U>
  bool convert_int(const polynome & p,const index_t & deg,std::vector< T_unsigned<longlong,U> >  & v,longlong & maxp){
    typename std::vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    v.clear();
    v.reserve(itend-it);
    T_unsigned<longlong,U> gu;
    U u;
    maxp=0;
    longlong tmp;
    mpz_t tmpz;
    mpz_init(tmpz);
    index_t::const_iterator itit,ditbeg=deg.begin(),ditend=deg.end(),dit;
    for (;it!=itend;++it){
      u=0;
      itit=it->index.iptr->begin();
      for (dit=ditbeg;dit!=ditend;++itit,++dit)
	u=u*unsigned(*dit)+unsigned(*itit);
      gu.u=u;
      if (it->value.type==_INT_)
	gu.g=it->value.val;
      else {
	if (it->value.type!=_ZINT || mpz_sizeinbase(*it->value._ZINTptr,2)>62){
	  mpz_clear(tmpz);
	  return false;
	}
	mpz2longlong(it->value._ZINTptr,&tmpz,gu.g);
      }
      tmp=gu.g>0?gu.g:-gu.g;
      if (tmp>maxp)
	maxp=tmp;
      v.push_back(gu);
    }
    mpz_clear(tmpz);
    return true;
  }

  template<class U> void convert_longlong(const std::vector< T_unsigned<gen,U> > & p,std::vector< T_unsigned<longlong,U> > & pd){
    typename std::vector< T_unsigned<gen,U> >::const_iterator it=p.begin(),itend=p.end();
    pd.reserve(itend-it);
    for (;it!=itend;++it)
      pd.push_back(T_unsigned<longlong,U>(it->g.val,it->u));
  }

  template<class T,class U> void convert_from(const std::vector< T_unsigned<T,U> > & p,std::vector< T_unsigned<gen,U> > & pd){
    typename std::vector< T_unsigned<T,U> >::const_iterator it=p.begin(),itend=p.end();
    pd.reserve(itend-it);
    for (;it!=itend;++it)
      pd.push_back(T_unsigned<gen,U>(gen(it->g),it->u));
  }

  template<class T,class U>
  void convert_from(const std::vector< T_unsigned<T,U> > & v,const index_t & deg,polynome & p){
    typename std::vector< T_unsigned<T,U> >::const_iterator it=v.begin(),itend=v.end();
    index_t::const_reverse_iterator ditbeg=deg.rbegin(),ditend=deg.rend(),dit;
    p.dim=ditend-ditbeg;
    p.coord.clear();
    p.coord.reserve(itend-it);
    U u,prevu=0;
    index_t i(p.dim);
    index_t::iterator iitbeg=i.begin(),iit,iitback=i.end()-1;
    int k;
    for (--prevu;it!=itend;++it){
      u=it->u;
      if (prevu<u+*iitback){
	*iitback -= prevu-u;
	prevu=u;
      }
      else {
	prevu=u;
	for (k=p.dim-1,dit=ditbeg;dit!=ditend;++dit,--k){
	  i[k]=u % unsigned(*dit);
	  u = u/unsigned(*dit);
	}
      }
      p.coord.push_back(monomial<gen>(gen(it->g),i));
    }
  }

  template<class T,class U>
  bool is_content_trivially_1(const typename std::vector< T_unsigned<T,U> > & v,U mainvar){
#ifdef HASH_MAP_NAMESPACE
    typedef HASH_MAP_NAMESPACE::hash_map< U,T,hash_function_unsigned_object > hash_prod ;
    hash_prod produit; // try to avoid reallocation
    // cout << "hash " << clock() << endl;
#else
    typedef std::map<U,T> hash_prod;
    // cout << "small map" << endl;
    hash_prod produit; 
#endif
    U outer_index,inner_index;
    typename std::vector< T_unsigned<T,U> >::const_iterator it=v.begin(),itend=v.end();
    for (;it!=itend;++it){
      outer_index=it->u % mainvar;
      inner_index=it->u/mainvar;
      typename hash_prod::iterator jt=produit.find(outer_index),jtend=produit.end();
      if (jt==jtend){
	if (inner_index==0)
	  return true;
	produit[outer_index]=it->g;
      }
      else
	jt->second += it->g;
    }
  }

  template<class T,class U>
  T peval_x1_xn(const typename std::vector< T_unsigned<T,U> >::const_iterator it,const typename std::vector< T_unsigned<T,U> >::const_iterator itend,const std::vector<T> & v,const T & reduce,const std::vector<U> & vars){
    if (vars.empty())
      return it->g;
    if (vars.size()!=v.size())
      setdimerr();
    U mainvar=vars.front(),deg1;
    T x=v.front();
    if (vars.size()==1){ // do Horner like eval
    }
    std::vector<T> v1(v.begin()+1,v.end());
    std::vector<U> vars1(vars.begin()+1,vars.end());
    typename std::vector< T_unsigned<T,U> >::const_iterator it2;
    T ans=0;
    for (;it!=itend;){
      deg1=(it->u/mainvar)*mainvar;
      for (it2=it;it2!=itend;++it2){
	if (it2->u<deg1)
	  break;
      }
      T tmp=peval_x1_xn(it,it2,v1,reduce,vars1);
      ans += tmp*powmod(x,deg1,reduce);
    }
    return ans;
  }

  // eval p at x2...xn
  template<class T,class U>
  void peval_x2_xn(const std::vector< T_unsigned<T,U> > & p,const std::vector<T> & v,const std::vector<U> & vars,const T & reduce,std::vector< T_unsigned<T,U> > & res){
    U mainvar=vars.front(),deg1;
    typename std::vector< T_unsigned<T,U> >::const_iterator it=p.begin(),itend=p.end(),it2;
    for (;it!=itend;){
      deg1=(it->u/mainvar)*mainvar;
      for (it2=it;it2!=itend;++it2){
	if (it2->u<deg1)
	  break;
      }
      T tmp=peval_x1_xn(it,it2,v,vars);
      if (!is_zero(tmp))
	res.push_back(T_unsigned<T,U>(tmp,deg1));
    }
  }

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // NO_NAMESPACE_GIAC

#endif // _GIAC_THREADED_H
