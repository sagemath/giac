// -*- mode:C++ -*-
/*
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
#ifndef _GIAC_MONOMIAL_H_
#define _GIAC_MONOMIAL_H_
#include "first.h"
#include <iostream>
#include <fstream>
#include <functional>
#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include "index.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  extern int debug_infolevel;

  void setsizeerr(const std::string &);

  inline bool is_zero(const longlong & a){ return a==0; }
  inline bool is_zero(const long & a){ return a==0; }
  inline bool is_zero(const int & a){ return a==0; }
  inline bool is_zero(const double & a){ return a==0; }

  inline void type_operator_times(const double & a,const double & b,double & c){
    c = a*b;
  }

  inline void type_operator_plus_times(const double & a,const double & b,double & c){
    c += a*b;
  }

  inline void type_operator_plus_times_reduce(const double & a,const double & b,double & c,int reduce){
    c += a*b;
  }

  inline void type_operator_reduce(const double & a,const double & b,double & c,int reduce){
    c = a*b;
  }

  inline void type_operator_times(const longlong & a,const longlong & b,longlong & c){
    c=a*b;
  }

  inline void type_operator_plus_times(const longlong & a,const longlong & b,longlong & c){
    c += a*b;
  }

  inline void type_operator_plus_times_reduce(const longlong & a,const longlong & b,longlong & c,int reduce){
    c += a*b;
    if (reduce)
      c %= reduce;
  }

  inline void type_operator_reduce(const longlong & a,const longlong & b,longlong & c,int reduce){
    c = a*b;
    if (reduce)
      c %= reduce;
  }

  inline void type_operator_times(const int & b,const int & c,int & a){
    a=b*c;
  }

  inline void type_operator_plus_times_reduce(const int & b,const int & c,int & a,int reduce){
    if (reduce){
#ifdef _I386_   // a<-a+b*c mod m
      if (a<0) a+=reduce;
      asm volatile("imull %%ecx; \n\t" /* b*c in edx:eax */
		   "addl %%ebx,%%eax; \n\t" /* b*c+a */
		   "adcl $0x0,%%edx; \n\t" /* b*c+a carry */
		   "idivl %%edi; \n\t"
		   :"=d"(a)
		   :"a"(b),"b"(a),"c"(c),"D"(reduce)
		   );
#else
      longlong tmp=longlong(b)*c;
      tmp+=a;
      a=tmp%reduce;
#endif
    }
    else
      a+=b*c;      
  }

  inline void type_operator_plus_times(const int & b,const int & c,int & a){
    a+=b*c;
  }

  inline void type_operator_reduce(const int & b,const int & c,int & a,int reduce){
    if (reduce){
#ifdef _I386_   // a<-b*c mod m
      asm volatile("imull %%ecx; \n\t" /* b*c in edx:eax */
		   "idivl %%edi; \n\t"
		   :"=d"(a)
		   :"a"(b),"c"(c),"D"(reduce)
		   );
#else
      longlong tmp=longlong(b)*c;
      a=tmp%reduce;
#endif
    }
    else
      a=b*c;      
  }

  inline bool has_denominator(int a){ return false; }
  inline bool has_denominator(longlong a){ return false; }

  template<class T,class U>
  struct T_unsigned {
    T g;
    U u;
    T_unsigned(const T & myg,const U & myu): g(myg),u(myu) {};
    T_unsigned(): g(0),u(0) {};
  };

  // warning, < is > so that monomial ordering is ok after back-conversion
  template <class T,class U>
  inline bool operator < (const T_unsigned<T,U> & gu1,const T_unsigned<T,U> & gu2){
    return gu1.u > gu2.u;
  }
  typedef T_unsigned<int,unsigned> int_unsigned;

  class hash_function_unsigned_object {
  public:
    // inline size_t operator () (size_t x) const { return x; }
    inline size_t operator () (int x) const { return x; }
    inline size_t operator () (unsigned x) const { return x; }
    inline size_t operator () (long long unsigned x) const { return x%12345701; }
    inline size_t operator () (long long x) const { return x%12345701; }
    hash_function_unsigned_object() {};
  };

  // a class for monomials, poly&matrices are std::vectors of couples index/monomial
  template <class T> class monomial{
  public:
    index_m index;
    T value;
    // constructors
    monomial(const monomial<T> & m) : index(m.index), value(m.value) {}
    monomial( const T & v, int dim) : value(v) {
      index.iptr->clear();
      index.iptr->reserve(dim);
      for (int i=1;i<=dim;i++)
	index.iptr->push_back(0);
    }
    monomial( const T & v, int var,int dim) : value(v) {
      index.iptr->clear();
      index.iptr->reserve(dim);
      for (int i=1;i<=dim;i++)
	index.iptr->push_back(i==var);
    }
    monomial( const T & v, int deg,int var,int dim) : value(v) {
      index.iptr->clear();
      index.iptr->reserve(dim);
      for (int i=1;i<=dim;i++)
	index.iptr->push_back(deg*(i==var));
    }
    monomial( const T & v, const index_m & i) : index(i),value(v) {}
    // unary negation
    monomial<T> operator - () const {
      return(monomial<T>(-(*this).value,(*this).index));
    }
    // members
    monomial<T> shift(index_m i,const T & fois) const {
      return monomial<T>(value*fois,i+index);
    }
    monomial<T> shift(const T & fois,index_m i) const {
      return monomial<T>(value/fois,i+index);
    }
    monomial<T> shift(index_m i) const {
      return monomial<T>(value,i+index);
    }
    void divided_by_x(){
      (*index.iptr)[0]--;
    }
    void reverse()  {
      int s=index.iptr->size();
      index_m new_i;
      new_i.iptr->reserve(s);
      index_t::const_iterator it=index.iptr->begin();
      index_t::const_iterator itend=index.iptr->end();
      --it;
      --itend;
      for (;it!=itend;--itend)
	new_i.iptr->push_back(*itend);
      index=new_i;
    }
    void reorder(const index_t & permutation)  {
      int s=index.iptr->size();
      if (unsigned(s)!=permutation.size())
	setsizeerr("Error monomial.h reorder(const index_t &)");
      index_m new_i(s);
      index_t::iterator newit=new_i.iptr->begin();
      for (int i=0;i<s;++newit,++i)
	*newit=(*index.iptr)[permutation[i]];
      index=new_i;
    }
    // truncate topmost index value (decrement by 1 the dimension)
    inline monomial<T> trunc1 () const {
      index_t & itt=*index.iptr;
#ifdef DEBUG_SUPPORT
      assert(!itt.empty());
#endif
      return monomial<T>(value,index_m(itt.begin()+1,itt.end()));
    }
    monomial<T> untrunc1 (int j=0) const {
      index_t::const_iterator it=index.iptr->begin(),itend=index.iptr->end();
      index_m new_i(itend-it+1);
      index_t::iterator newit=new_i.iptr->begin();    
      *newit=j;
      for (++newit;it!=itend;++newit,++it)
	*newit=*it;
      return monomial<T>(value,new_i);
    }
    // add a principal degree equal to j and adjust dimension to dim
    monomial<T> untrunc (int j,int dim) const {
      int s=index.iptr->size();
      assert(s<dim);
      index_m new_i(dim);
      index_t::const_iterator it=index.iptr->begin();
      index_t::iterator newit=new_i.iptr->begin();
      *newit=j;
      for (++newit,--dim;dim>s;++newit,--dim)
	*newit=0;
      for (;it!=index.iptr->end();++newit,++it)
	*newit=*it;
      return monomial<T>(value,new_i);
    }
    inline T norm () const{
      return abs(value);
    }
    std::string print() const {
      return "%%%{"+value.print()+','+print_INT_(* index.iptr)+ "%%%}";
    }
  };
  
  template <class T>
  std::ostream & operator << (std::ostream & os, const monomial<T> & m ){
    return os << m.print();
  }

  template<class T>
  bool operator == (const monomial<T>& a,const monomial<T> & b){
    return (a.value==b.value) && (a.index==b.index);
  }

  template <class T>
  monomial<T> Untrunc1(const T & t,int j=0){
    index_m new_i;
    new_i.i()->push_back(j);
    return monomial<T>(t,new_i);
  }

  // ordering monomials using index ordering
  template <class T>
  bool m_total_lex_is_greater(const monomial<T> & m1, const monomial<T> & m2){
    return(i_total_lex_is_greater(m1.index,m2.index));
  }
  template <class T>
  bool m_lex_is_greater(const monomial<T> & m1, const monomial<T> & m2){
    return(i_lex_is_greater(m1.index,m2.index));
  }

  template <class T>
  bool m_total_revlex_is_greater(const monomial<T> & m1, const monomial<T> & m2){
    return(i_total_revlex_is_greater(m1.index,m2.index));
  }

  template <class T>
  std::ostream & operator << (std::ostream & os, const  std::vector<T> v){
    typename std::vector<T>::const_iterator it=v.begin();
    typename std::vector<T>::const_iterator itend=v.end();
    os << "Vector [";
    for (;it!=itend;){
      os << *it ;
      ++it;
      if (it!=itend)
	os << ",";
    }
    os << "]" ;
    return os;
  }

  /*
    template <class T>
    std::ostream & operator << (std::ostream & os, const std::vector< monomial<T> > & v){
    for (typename std::vector< monomial<T> >::const_iterator it=v.begin();it!=v.end();++it)
    os << *it << "| " ;
    return(os);
    }
  */

  template <class T>
  void Mul ( typename std::vector< monomial<T> >::const_iterator & a,
	     typename std::vector< monomial<T> >::const_iterator & a_end,
	     const T & fact, std::vector< monomial<T> > & new_coord){
    if (new_coord.begin()==a){
      if (is_one(fact))
	return;
      typename std::vector< monomial<T> >::iterator b=new_coord.begin(),b_end=new_coord.end();
      for (;b!=b_end;++b){
	b->value = b->value * fact;
      }
    }
    else {
      new_coord.reserve(a_end - a );
      for (;a!=a_end;++a){
	new_coord.push_back(monomial<T>( ((*a).value) * fact , (*a).index) );
      }
    }
  }

  template <class T>
  std::vector< monomial<T> > operator * (const T & f,const std::vector< monomial<T> > & v){
    if (is_one(f))
      return v;
    typename std::vector< monomial<T> >::const_iterator a=v.begin(), a_end=v.end();
    std::vector< monomial<T> > res;
    if (!is_zero(f))
      Mul(a,a_end,f,res);
    return res ;
  }

  template <class T>
  std::vector< monomial<T> > operator * (const std::vector< monomial<T> > & v,const T & f){
    return f*v;
  }

  template <class T>
  void Div ( typename std::vector< monomial<T> >::const_iterator & a,
	     typename std::vector< monomial<T> >::const_iterator & a_end,
	     const T & fact, std::vector< monomial<T> > & new_coord){
    new_coord.reserve(a_end - a );
    for (;a!=a_end;++a){
      new_coord.push_back(monomial<T>( rdiv((*a).value, fact) , (*a).index) );
    }
  }

  template <class T>
  std::vector< monomial<T> > operator / (const std::vector< monomial<T> > & v,const T & f){
    if (is_one(f))
      return v;
    typename std::vector< monomial<T> >::const_iterator a=v.begin(), a_end=v.end();
    std::vector< monomial<T> > res;
    res.clear();
    Div(a,a_end,f,res);
    return res ;
  }

  template <class T>
  void Add ( typename std::vector< monomial<T> >::const_iterator & a,
	     typename std::vector< monomial<T> >::const_iterator & a_end,
	     typename std::vector< monomial<T> >::const_iterator & b,
	     typename std::vector< monomial<T> >::const_iterator & b_end,
	     std::vector< monomial<T> > & new_coord,
	     const std::pointer_to_binary_function< const index_t &, const index_t &, bool> is_strictly_greater) {
    new_coord.clear();
    new_coord.reserve( (a_end - a) + (b_end - b));
    // bool log=false;
    /* if (a!=a_end)
       log=a->index.iptr->size()>=12; 
       if (log)
       cerr << "+ begin" << clock() << endl; */
    for (;;) {
      if (a == a_end) {
	while (b != b_end) {
	  new_coord.push_back(*b);
	  ++b;
	}
	break;
      } 
      const index_m & pow_a = a->index;
    
      // If b is empty, fill up with elements from a and stop
      if (b == b_end) {
	while (a != a_end) {
	  new_coord.push_back(*a);
	  ++a;
	}
	break;
      } 

      const index_m & pow_b = b->index;
        
      // a and b are non-empty, compare powers
      if (pow_a!=pow_b){
	if (is_strictly_greater(*pow_a.iptr, *pow_b.iptr)) {
	  // a has lesser power, get coefficient from a
	  new_coord.push_back(*a);
	  ++a;
	} 
	else  {
	  // b has lesser power, get coefficient from b
	  new_coord.push_back(*b);
	  ++b;
	} 
      }
      else {
	T sum = (*a).value + (*b).value;
	if (!is_zero(sum))
	  new_coord.push_back(monomial<T>(sum,pow_a));
	++a;
	++b;
      }
    }
    //  if (log)
    //  std::cerr << "+ end " << clock() << endl;
  }

  template <class T>
  std::vector< monomial<T> > operator + (const std::vector< monomial<T> > & v,const std::vector< monomial<T> > & w){
    typename std::vector< monomial<T> >::const_iterator a=v.begin(), a_end=v.end();
    typename std::vector< monomial<T> >::const_iterator b=w.begin(), b_end=w.end();
    std::vector< monomial<T> > res;
    Add(a,a_end,b,b_end,res,lex_is_greater<int>);
    return res ;
  }

  template <class T>
  void Sub ( typename std::vector< monomial<T> >::const_iterator & a,
	     typename std::vector< monomial<T> >::const_iterator & a_end,
	     typename std::vector< monomial<T> >::const_iterator & b,
	     typename std::vector< monomial<T> >::const_iterator & b_end,
	     std::vector< monomial<T> > & new_coord,
	     const std::pointer_to_binary_function< const index_t &, const index_t &, bool> is_strictly_greater) {
    new_coord.clear();
    new_coord.reserve( (a_end - a) + (b_end - b));
    for (;;) {
      // If a is empty, fill up with elements from b and stop
      if (a == a_end) {
	while (b != b_end) {
	  new_coord.push_back(-(*b));
	  ++b;
	}
	break;
      } 
      const index_m & pow_a = a->index;
    
      // If b is empty, fill up with elements from a and stop
      if (b == b_end) {
	while (a != a_end) {
	  new_coord.push_back(*a);
	  ++a;
	}
	break;
      } 
      const index_m & pow_b = b->index;
        
      // a and b are non-empty, compare powers
      if (pow_a!=pow_b){
	if (is_strictly_greater(*pow_a.iptr, *pow_b.iptr)) {
	  // a has lesser power, get coefficient from a
	  new_coord.push_back(*a);
	  ++a;
	} 
	else  {
	  // b has lesser power, get coefficient from b
	  new_coord.push_back(-(*b));
	  ++b;
	} 
      }
      else {
	T diff = (*a).value - (*b).value;
	if (!is_zero(diff))
	  new_coord.push_back(monomial<T>(diff,pow_a));
	++a;
	++b;
      }
    }  
  }

  template <class T>
  std::vector< monomial<T> > operator - (const std::vector< monomial<T> > & v,const std::vector< monomial<T> > & w){
    typename std::vector< monomial<T> >::const_iterator a=v.begin(), a_end=v.end();
    typename std::vector< monomial<T> >::const_iterator b=w.begin(), b_end=w.end();
    std::vector< monomial<T> > res;
    Sub(a,a_end,b,b_end,res,lex_is_greater<int>);
    return res ;
  }

  template <class T>
  void addsamepower(typename std::vector< monomial<T> >::const_iterator & it,
		    typename std::vector< monomial<T> >::const_iterator & itend,
		    std::vector< monomial<T> > & new_coord){
    while (it!=itend){
      T res=(*it).value;
      index_m pow=(*it).index;
      ++it;
      while ( (it!=itend) && ((*it).index==pow)){
	res=res+(*it).value;
	++it;
      }
      if (!is_zero(res))
	new_coord.push_back(monomial<T>(res, pow));
    }
  }

  struct ltindex
  {
    std::pointer_to_binary_function < const index_t &, const index_t &, bool> is_strictly_greater;
    inline bool operator()(const index_t & s1, const index_t & s2) const
    {
      return is_strictly_greater(s1, s2) ;
    }
    ltindex(const std::pointer_to_binary_function < const index_t &, const index_t &, bool> my_is_strictly_greater) : is_strictly_greater(my_is_strictly_greater) {};

  };

  template <class T>
  void Mul ( typename std::vector< monomial<T> >::const_iterator & ita,
	     typename std::vector< monomial<T> >::const_iterator & ita_end,
	     typename std::vector< monomial<T> >::const_iterator & itb,
	     typename std::vector< monomial<T> >::const_iterator & itb_end,
	     std::vector< monomial<T> > & new_coord,
	     const std::pointer_to_binary_function < const index_t &, const index_t &, bool> is_strictly_greater,
	     const std::pointer_to_binary_function < const monomial<T> &, const monomial<T> &, bool> m_is_greater
	     ) {
    // another algorithm using a hash_map 
#ifdef HASH_MAP_NAMESPACE
    typedef HASH_MAP_NAMESPACE::hash_map< index_t,T,hash_function_object > hash_prod ;
    hash_prod produit_;
    index_t sum_(ita->index.iptr->size());
    index_t * it_aindexptr;
    index_t::const_iterator it_aindex,it_aindexbeg,it_aindexend,it_bindex;
    index_t::iterator it_sum_index,it_sum_beg=sum_.begin();
    typename hash_prod::iterator prod_it_,prod_it_end;
    typename std::vector< monomial<T> >::const_iterator it_a_cur=ita,it_b_cur;
    for ( ; it_a_cur!=ita_end; ++it_a_cur ){
      it_aindexptr= it_a_cur->index.iptr;
      it_aindexbeg=it_aindexptr->begin();
      it_aindexend=it_aindexptr->end();
      for ( it_b_cur=itb;it_b_cur!=itb_end;++it_b_cur) {
	it_bindex=it_b_cur->index.iptr->begin();
	it_sum_index=it_sum_beg;
	for (it_aindex=it_aindexbeg;it_aindex!=it_aindexend;++it_bindex,++it_sum_index,++it_aindex)
	  *it_sum_index=(*it_aindex)+(*it_bindex);
	prod_it_=produit_.find(sum_);
	if (prod_it_==produit_.end())
	  produit_[sum_]=it_a_cur->value*it_b_cur->value;
	else 	
	  prod_it_->second += it_a_cur->value*it_b_cur->value;
      }
    }
    prod_it_=produit_.begin(),prod_it_end=produit_.end();
    new_coord.reserve(produit_.size());
    for (;prod_it_!=prod_it_end;++prod_it_)
      if (!is_zero(prod_it_->second))
	new_coord.push_back(monomial<T>(prod_it_->second,prod_it_->first));
    // cerr << new_coord <<endl;
    sort(new_coord.begin(),new_coord.end(),m_is_greater);
    return ;
#endif

    /* other algorithm using a map to avoid reserving too much space */
    typedef std::map< index_t,T,const std::pointer_to_binary_function < const index_t &, const index_t &, bool> > application;
    application produit(is_strictly_greater);
    // typedef std::map<index_t,T> application;
    // application produit;
    index_t somme(ita->index.iptr->size());
    index_t * itaindexptr;
    index_t::const_iterator itaindex,itaindexbeg,itaindexend,itbindex;
    index_t::iterator itsommeindex,itsommebeg=somme.begin();
    typename application::iterator prod_it,prod_itend;
    typename std::vector< monomial<T> >::const_iterator ita_cur=ita,itb_cur;
    for ( ; ita_cur!=ita_end; ++ita_cur ){
      itaindexptr= ita_cur->index.iptr;
      itaindexbeg=itaindexptr->begin();
      itaindexend=itaindexptr->end();
      for ( itb_cur=itb;itb_cur!=itb_end;++itb_cur) {
	itbindex=itb_cur->index.iptr->begin();
	itsommeindex=itsommebeg;
	for (itaindex=itaindexbeg;itaindex!=itaindexend;++itbindex,++itsommeindex,++itaindex)
	  *itsommeindex=(*itaindex)+(*itbindex);
	prod_it=produit.find(somme);
	if (prod_it==produit.end())
	  produit[somme]=ita_cur->value*itb_cur->value;
	else 	
	  prod_it->second += ita_cur->value*itb_cur->value;
      }
    }
    prod_it=produit.begin(),prod_itend=produit.end();
    new_coord.reserve(produit.size());
    for (;prod_it!=prod_itend;++prod_it)
      if (!is_zero(prod_it->second))
	new_coord.push_back(monomial<T>(prod_it->second,prod_it->first));
    // cerr << new_coord <<endl;
    // sort(new_coord.begin(),new_coord.end(),m_is_greater);
  
    /* old algorithm
       std::vector< monomial<T> > multcoord;
       int asize=ita_end-ita,bsize=itb_end-itb;
       int d=ita->index.iptr->size();
       multcoord.reserve(asize*bsize); // correct for sparse polynomial
       typename std::vector< monomial<T> >::const_iterator ita_begin = ita,itb_begin=itb ;
       index_m old_pow=(*ita).index+(*itb).index;
       T res( 0);
       for ( ; ita!=ita_end; ++ita ){
       typename std::vector< monomial<T> >::const_iterator ita_cur=ita;
       typename std::vector< monomial<T> >::const_iterator itb_cur=itb;
       for (;itb_cur!=itb_end;--ita_cur,++itb_cur) {
       index_m cur_pow=(*ita_cur).index+(*itb_cur).index;
       if (cur_pow!=old_pow){
       if (!is_zero(res))
       multcoord.push_back( monomial<T>(res ,old_pow ));
       res=((*ita_cur).value) * ((*itb_cur).value);
       old_pow=cur_pow;
       }
       else
       res=res+((*ita_cur).value) * ((*itb_cur).value);      
       if (ita_cur==ita_begin)
       break;
       }
       }
       --ita;
       ++itb;
       for ( ; itb!=itb_end;++itb){
       typename std::vector< monomial<T> >::const_iterator ita_cur=ita;
       typename std::vector< monomial<T> >::const_iterator itb_cur=itb;
       for (;itb_cur!=itb_end;--ita_cur,++itb_cur) {
       index_m cur_pow=(*ita_cur).index+(*itb_cur).index;
       if (cur_pow!=old_pow){
       if (!is_zero(res))
       multcoord.push_back( monomial<T>(res ,old_pow ));
       res=((*ita_cur).value) * ((*itb_cur).value);
       old_pow=cur_pow;
       }
       else
       res=res+((*ita_cur).value) * ((*itb_cur).value);
    
       if (ita_cur==ita_begin)
       break;
       }
       }
       // push last monomial
       if (!is_zero(res))
       multcoord.push_back( monomial<T>(res ,old_pow ));
       // sort by asc. power
       sort( multcoord.begin(),multcoord.end(),m_is_greater);
       typename std::vector< monomial<T> >::const_iterator it=multcoord.begin();
       typename std::vector< monomial<T> >::const_iterator itend=multcoord.end();
       // adjust result size 
       // statistics about polynomial density
       // a dense poly of deg. aa and d variables has binomial(aa+d,d) monomials
       // we need to reserve at most asize*bsize
       // but less for dense polynomials since 
       // binomial(aa+d,d)*binomial(bb+d,d) > binomial(aa+bb+d,d)
       int aa=total_degree(ita_begin->index),bb=total_degree(itb_begin->index);
       double r;
       double factoriald=std::lgamma(d+1);
       // double factorialaa=std::lgamma(aa+1),factorialbb=std::lgamma(bb+1);
       // double factorialaad=std::lgamma(aa+d+1),factorialbbd=std::lgamma(bb+d+1);
       double factorialaabbd=std::lgamma(aa+bb+d+1),factorialaabb=std::lgamma(aa+bb+1);
       r=std::exp(factorialaabbd-(factorialaabb+factoriald));
       if (debug_infolevel)
       cerr << "// Mul degree " << aa << "+" << bb << " size " << asize << "*" << bsize << "=" << asize*bsize << " max " << r << endl;
       new_coord.reserve(min(int(r),itend-it));
       // add terms with same power
       addsamepower(it,itend,new_coord);
       if (debug_infolevel)
       cerr << "// Actual mul size " << new_coord.size() << endl;
    */
  }

  // polynomial multiplication
  template <class T>
  std::vector< monomial<T> > operator * (const std::vector< monomial<T> > & v,const std::vector< monomial<T> > & w){
    typename std::vector< monomial<T> >::const_iterator a=v.begin(), a_end=v.end();
    typename std::vector< monomial<T> >::const_iterator b=w.begin(), b_end=w.end();
    std::vector< monomial<T> > res;
    Mul(a,a_end,b,b_end,res,m_lex_is_greater);
    return res ;
  }

  template <class T>
  void Shift (const std::vector< monomial<T> > & v,const index_m &i, const T & fois, std::vector< monomial<T> > & new_coord){
    new_coord.clear();
    typename std::vector< monomial<T> >::const_iterator itend=v.end();
    for (typename std::vector< monomial<T> >::const_iterator it=v.begin();it!=itend;++it)
      new_coord.push_back( it->shift(i,fois) );
  }

  template <class T>
  void Shift (const std::vector< monomial<T> > & v,const T & fois, const index_m &i, std::vector< monomial<T> > & new_coord){
    new_coord.clear();
    typename std::vector< monomial<T> >::const_iterator itend=v.end();
    for (typename std::vector< monomial<T> >::const_iterator it=v.begin();it!=itend;++it)
      new_coord.push_back( it->shift(fois,i) );
  }

  template <class T>
  void Shift (const std::vector< monomial<T> > & v,const index_m &i, std::vector< monomial<T> > & new_coord){
    new_coord.clear();
    typename std::vector< monomial<T> >::const_iterator itend=v.end();
    for (typename std::vector< monomial<T> >::const_iterator it=v.begin();it!=itend;++it)
      new_coord.push_back( it->shift(i) );
  }

  template <class T>
  T Content (const std::vector< monomial<T> > & v){
    typename std::vector< monomial<T> >::const_iterator it=v.begin();
    typename std::vector< monomial<T> >::const_iterator itend=v.end();
    if (it==itend)
      return 1;
    T res=it->value;
    for (;it!=itend ;++it){
      res=gcd(res,it->value);
      if (is_one(res))
	break;
    }
    return res;
  }

  template<class T>
  T ppz(std::vector< monomial<T> > & p){
    T n=Content(p);
    p=p/n;
    return n;
  }

  template<class T>
  void Nextcoeff(typename std::vector< monomial<T> >::const_iterator & it,const typename std::vector< monomial<T> >::const_iterator & itend,std::vector< monomial<T> > & v){
    int n=it->index.iptr->front();
    int d=it->index.iptr->size();
    for (;(it!=itend) && (it->index.iptr->front()==n);++it)
      v.push_back(it->trunc1());
  }


  template <class T>
  void Trunc1(const std::vector< monomial<T> > & c,std::vector< monomial<T> > & v){
    v.reserve(c.size());
    typename std::vector< monomial<T> >::const_iterator it=c.begin();
    typename std::vector< monomial<T> >::const_iterator itend=c.end();
    for (;it!=itend;++it)
      v.push_back(it->trunc1());
  }

  template <class T>
  void Untrunc1(const std::vector< monomial<T> > & c,int j,std::vector< monomial<T> > & v){
    v.reserve(c.size());
    typename std::vector< monomial<T> >::const_iterator it=c.begin();
    typename std::vector< monomial<T> >::const_iterator itend=c.end();
    for (;it!=itend;++it)
      v.push_back(it->untrunc1(j));
  }

  template <class T>
  void Untruncn(std::vector< monomial<T> > & c,int j){
    typename std::vector< monomial<T> >::iterator it=c.begin();
    typename std::vector< monomial<T> >::iterator itend=c.end();
    index_t tmp;
    for (;it!=itend;++it){
      index_t & i=*it->index.iptr;
      if (*(it->index.ref_count)==1)
	i.push_back(j);
      else {
	tmp=i;
	tmp.push_back(j);
	it->index=tmp;
      }
    }
  }

  template <class T>
  void Untrunc(const std::vector< monomial<T> > & c,int j,int dim,std::vector< monomial<T> > & v){
    v.reserve(c.size());
    typename std::vector< monomial<T> >::const_iterator it=c.begin();
    typename std::vector< monomial<T> >::const_iterator itend=c.end();
    for (;it!=itend;++it)
      v.push_back(it->untrunc(j,dim));
  }

  template<class T>
  void Apply(typename std::vector< monomial<T> >::const_iterator & it,const typename std::vector< monomial<T> >::const_iterator & itend,T (*f)(const T &),std::vector< monomial<T> > & v ){
    v.reserve(itend-it);
    T temp;
    for (;it!=itend;++it){
      temp=f(it->value);
      if (!is_zero(temp))
	v.push_back(monomial<T>(temp,it->index));
    }
  }

  template <class T>
  bool Findpivot(std::vector< monomial<T> > & v,int rows,int cols, index_t & permut,int &pivotr,int & pivotc,std::vector< monomial<T> > & pivotline,T * pivotcol, T & pivotcoeff){
    bool found=false;
    T zero(0);
    for (int i=1;i<=rows;i++)
      pivotcol[i]=zero;
    T pivotrrefnorm;
    typename std::vector< monomial<T> >::const_iterator itend=v.end();
    for ( ; (pivotc<=cols) && (!found) ; ){
      pivotr=1;
      // search in lines not already in permut for the best pivot
      // set found to true as soon as one is non 0
      typename std::vector< monomial<T> >::const_iterator it=v.begin();
      for (;it!=itend;++it){
	if ((*it->index.iptr)[1]==pivotc){	
	  int r=it->index.iptr->front();
	  T val=it->value;
	  pivotcol[r]=val;
	  if ( !has(permut,r) ){
	    if ( (!found) || ( rrefnorm(val)>pivotrrefnorm) ) {
	      found = true;
	      pivotcoeff=val;
	      pivotr=r;
	      pivotrrefnorm=rrefnorm(pivotcoeff);
	    } 
	  }
	}
      }
      if (!found)
	pivotc++;
    }
    if (!found)
      return false;
    permut.push_back(pivotr);
    // pivot has been found, compute pivotline = line of the pivot shifted by x^
    pivotline.clear();
    typename std::vector< monomial<T> >::const_iterator it=v.begin();
    for (;it!=itend;++it){
      if (it->index.iptr->front()==pivotr)
	break;
    }
    for (;it!=itend;++it){
      if (it->index.iptr->front()!=pivotr)
	break;
      pivotline.push_back(it->trunc1());
    }
  }

  template <class T>
  void Rref (std::vector< monomial<T> > & v,int rows,int cols, index_t & permut, bool dobareiss=true){
    T bareisscoeff=1;
    permut.clear();
    permut.push_back(0); // 0 -> 0 convention for permutations
    T pivotcol[rows+1];
    int pivotr;
    for (int pivotc=1;pivotc<=cols;pivotc++){
      std::vector< monomial<T> > pivotline;
      std::vector< monomial<T> > newcoord;
      T pivotcoeff;
      if (!Findpivot(v,rows,cols,permut,pivotr,pivotc,pivotline,pivotcol,pivotcoeff))
	break;
      newcoord.reserve(v.size());
      typename std::vector< monomial<T> >::const_iterator it=v.begin(),itend=v.end(),tmpit,tmpitend;
      while (it!=itend){
	int r=it->index.iptr->front();
	std::vector< monomial<T> > temp;
	Nextcoeff(it,itend,temp);
	if (r!=pivotr){
	  // std::cout << "L" << r << "=" << pivotcoeff << "*L" << r << "-" << pivotcol[r] << "*L" << pivotr << std::endl ;
	  temp=temp*pivotcoeff-pivotline*pivotcol[r];
	  if (dobareiss)
	    temp=temp/bareisscoeff;
	}
	typename std::vector< monomial<T> >::const_iterator tmpit,tmpitend=temp.end();
	for (tmpit=temp.begin();tmpit!=tmpitend;++tmpit)
	  newcoord.push_back(tmpit->untrunc1(r));
      }
      bareisscoeff=pivotcoeff;
      v=newcoord;
      // std::cout << v << std::endl;
    }
  }

  template <class T>
  void Findbeginofrows(typename std::vector< monomial<T> >::const_iterator & it,const typename std::vector< monomial<T> >::const_iterator & itend,const int & cols,typename std::vector< monomial<T> >::const_iterator * beg){
    for (int j=0;j<=cols;j++)
      beg[j]=0;
    for (;it!=itend;){
      int row=it->index.iptr->front();
      int col=(*it->index.iptr)[1];
      beg[col]=it;
      while ((it!=itend) && (row==it->index.iptr->front()))
	++it;
    }
  }

  template <class T>
  void Findbeginofrows(typename std::vector< monomial<T> >::iterator & it,const typename std::vector< monomial<T> >::const_iterator & itend,const int & cols,typename std::vector< monomial<T> >::iterator * beg){
    for (int j=0;j<=cols;j++)
      beg[j]=0;
    for (;it!=itend;){
      int row=it->index.iptr->front();
      int col=(*it->index.iptr)[1];
      beg[col]=it;
      while ((it!=itend) && (row==it->index.iptr->front()))
	++it;
    }
  }


  template <class T>
  void Normalrref (std::vector< monomial<T> > & v,int rows,int cols, index_t & permut, bool dobareiss=true){
    Rref(v,rows,cols,permut,dobareiss);
    // divide each non-zero row by leading coeff and order
    typename std::vector< monomial<T> >::const_iterator it=v.begin(),itend=v.end();
    typename std::vector< monomial<T> >::const_iterator beg[cols+1];
    std::vector< monomial<T> > temp,newcoord;
    Findbeginofrows(it,itend,cols,beg);
    int r=1;
    for (int j=1;j<=cols;j++){
      it=beg[j];
      if (it){
	T val=it->value;
	temp.clear();
	Nextcoeff(it,itend,temp);
	Untrunc1(temp/val,r,newcoord);
	r++;
      }
      if (it==itend)
	break;
    }
    v=newcoord;
  }

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC

#endif // ndef _GIAC_MONOMIAL_H
