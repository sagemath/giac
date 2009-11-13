// -*- mode: C++ ; compile-command: "g++ -I.. -g -O2 -c index.cc" -*-
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
#ifndef _GIAC_INDEX_H_
#define _GIAC_INDEX_H_
#include "first.h"
#include <vector>
#include <iostream>
#include <string>

#if defined UNORDERED_MAP && !defined(__APPLE__)
#include <tr1/unordered_map>
#define HASH_MAP_NAMESPACE std::tr1
#define hash_map unordered_map
#else // UNORDERED_MAP

#ifdef HASH_MAP
#include <hash_map>
#ifndef HASH_MAP_NAMESPACE
#ifndef VISUALC 
#define HASH_MAP_NAMESPACE std
#endif // VISUALC
#endif // HASH_MAP_NAMESPACE
#endif

#ifdef EXT_HASH_MAP
#include <ext/hash_map>
#ifndef HASH_MAP_NAMESPACE
#define HASH_MAP_NAMESPACE __gnu_cxx
#endif
#endif

#endif // UNORDERED_MAP

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  int mygcd(int a,int b);

  // index type for tensors
  typedef std::vector<int> index_t;

  index_t operator + (const index_t & a, const index_t & b);
  index_t operator - (const index_t & a, const index_t & b);
  index_t operator | (const index_t & a, const index_t & b);
  index_t operator - (const index_t & a);
  index_t operator * (const index_t & a, int fois);
  inline index_t operator * (int fois,const index_t & a){ return a*fois; }
  index_t operator / (const index_t & a, int divisepar);
  int operator / (const index_t & a, const index_t & b);
  // >= and <= are *partial* ordering on index_t
  // they return TRUE if and only if >= or <= is true for *all* coordinates
  bool all_sup_equal (const index_t & a, const index_t & b);
  inline bool operator >= (const index_t & a, const index_t & b){ return all_sup_equal(a,b); }
  bool all_inf_equal (const index_t & a, const index_t & b);
  inline bool operator <= (const index_t & a, const index_t & b){ return all_inf_equal(a,b); }
  index_t index_gcd(const index_t & a,const index_t & b);
  index_t index_lcm(const index_t & a,const index_t & b);
  inline index_t index_min(const index_t & a,const index_t & b){ return index_gcd(a,b); }
  inline index_t index_max(const index_t & a,const index_t & b){ return index_lcm(a,b); }
  void dbgprint(const index_t & i);
  std::string print_INT_(int i);
  std::string hexa_print_INT_(int i);
  std::string octal_print_INT_(int i);
  std::string binary_print_INT_(int i);
  std::string print_INT_(const index_t & v);

  template <class T> T pow(const std::vector<T> & x, const index_t & n );

  // total degree of a std::vector
  template <class T> T total_degree(const std::vector<T> & v1);

  // two ordering function over indices: lex ordering and total order then lex
  template <class T>
  bool lex_is_greater(const std::vector<T> & v1, const std::vector<T> & v2);
  template <class T>
  bool total_revlex_is_greater(const std::vector<T> & v1, const std::vector<T> & v2);
  template <class T>
  bool total_lex_is_greater(const std::vector<T> & v1, const std::vector<T> & v2);

  index_t mergeindex(const index_t & i,const index_t & j);
  // index_t might be used for permutations
  index_t inverse(const index_t & p);
  // transposition
  index_t transposition(int i,int j,int size);
  bool has(const index_t & p,int r);
  // zero?
  bool is_zero(const index_t & p);

  template <class T> T pow(const std::vector<T> & x, const index_t & n ){
    assert(x.size()==n.size());
    typename std::vector<T>::const_iterator itx=x.begin();
    index_t::const_iterator itn=n.begin();
    T res(1);
    for (;itx!=x.end();++itx,++itn){
      res=res*pow(*itx,*itn);
    }
    return res;
  }

  template <class T> T total_degree(const std::vector<T> & v1){
    T i=0;
    for (typename std::vector<T>::const_iterator it=v1.begin();it!=v1.end();++it)
      i=i+(*it);
    return(i);
  }

  template <class T>
  bool lex_is_greater(const std::vector<T> & v1, const std::vector<T> & v2){
    assert(v1.size()==v2.size());
    typename std::vector<T>::const_iterator it1=v1.begin(),it1end=v1.end();
    typename std::vector<T>::const_iterator it2=v2.begin();
    for (;it1!=it1end;++it2,++it1){
      if ( (*it1)!=(*it2) ){
	if  ( (*it1)>(*it2))
	  return(true);
	else
	  return(false);
      }
    }
    return(true);
  }

  template <class T>
  bool lex_is_strictly_greater(const std::vector<T> & v1, const std::vector<T> & v2){
    assert(v1.size()==v2.size());
    typename std::vector<T>::const_iterator it1=v1.begin(),it1end=v1.end();
    typename std::vector<T>::const_iterator it2=v2.begin();
    for (;it1!=it1end;++it2,++it1){
      if ( (*it1)!=(*it2) ){
	if  ( (*it1)>(*it2))
	  return true;
	else
	  return false;
      }
    }
    return false;
  }


  template <class T>
  bool total_lex_is_greater(const std::vector<T> & v1, const std::vector<T> & v2){
    T d1=total_degree(v1);
    T d2=total_degree(v2);
    if (d1!=d2){
      if (d1>d2)
	return(true);
      else
	return(false);
    }
    return(lex_is_greater<T>(v1,v2));
  }

  template <class T>
  bool total_revlex_is_greater(const std::vector<T> & v1, const std::vector<T> & v2){
    T d1=total_degree(v1);
    T d2=total_degree(v2);
    if (d1!=d2){
      if (d1>d2)
	return(true);
      else
	return(false);
    }
    return(!lex_is_strictly_greater<T>(v1,v2));
  }

  //*****************************************
  // class for memory efficient indices
  //*****************************************

  class index_m {
  public:
    inline index_t * i(){ return iptr; }
    index_m(const index_m & im)  { 
      iptr=im.iptr;
      ref_count=im.ref_count;
      (*ref_count)++;
    }
    index_m(const index_t & i);
    index_m(){
      iptr=new index_t();
      ref_count=new int(1);
    }
    index_m(size_t s){
      iptr=new index_t(s);
      ref_count=new int(1);
    }
    index_m(index_t::const_iterator it,index_t::const_iterator itend){
      iptr=new index_t(it,itend);
      ref_count=new int(1);
    }
    const index_m & operator = (const index_m & other){
      (*ref_count)--;
      if (!*ref_count){
	delete iptr;
	delete ref_count;
      }
      iptr=other.iptr;
      ref_count=other.ref_count;
      (*ref_count)++;
      return *this;
    }
    ~index_m(){
      (*ref_count)--;
      if (!*ref_count){
	delete iptr;
	delete ref_count;
      }
    }
    index_t * iptr;
    int * ref_count;
    friend std::ostream & operator << (std::ostream & os, const index_m & m ){
      os << ":index_m:[ " ;
      for (index_t::const_iterator it=m.iptr->begin();it!=m.iptr->end();++it)
	os << *it << " ";
      os << "] " ;
      return(os);
    }
    void dbgprint() const {
      std::cout << *this << std::endl;
    }
    // set first index element to 0
    index_m firstzero() const {
      index_t i(*(this->iptr));
      assert(i.size());
      i[0]=0;
      return index_m(i);
    }
  };

#ifdef HASH_MAP_NAMESPACE
  inline size_t index_hash_function(const std::vector<int> & v){
    std::vector<int>::const_iterator it=v.begin(),itend=v.end();
    size_t res=0;
    if (itend-it>16)
      itend=it+16;
    if (itend-it>8){
      for (;it!=itend;++it)
	res = (res << 2) | *it;
    }
    else {
      for (;it!=itend;++it)
	res = (res << 4) | *it;
    }
    return res;
  }

  /*
  inline size_t index_hash_function(const vector<int> & v){
    vector<int>::const_iterator it=v.begin(),itend=v.end();
    if (itend-it>16)
      itend=it+16;
    size_t res=0,decal=32/(itend-it);
    for (;;){
      --itend;
      res = (res << decal) | *itend;
      if (itend==it)
	return res;
    }
  }
  */

  class hash_function_object {
  public:
    size_t operator () (const index_t & v) const { return index_hash_function(v); }
    hash_function_object() {};
  };

  typedef HASH_MAP_NAMESPACE::hash_map< index_t,index_m,hash_function_object > hash_index ;  

  extern std::vector<hash_index> global_hash_index;
  extern int copy_number;

#endif

  index_m operator + (const index_m & a, const index_m & b);
  index_m operator - (const index_m & a, const index_m & b);
  index_m operator * (const index_m & a, int fois);
  inline index_m operator * (int fois,const index_m & a){ return a*fois; }
  index_m operator / (const index_m & a, int divisepar);
  inline int operator / (const index_m & a,const index_m & b){ return *a.iptr / *b.iptr;}
  bool operator == (const index_m & i1, const index_m & i2);
  bool operator != (const index_m & i1, const index_m & i2);
  bool operator >= (const index_m & a, const index_m & b);
  bool operator <= (const index_m & a, const index_m & b);
  int total_degree(const index_m & v1);
  bool i_lex_is_greater(const index_m & v1, const index_m & v2);
  bool i_lex_is_strictly_greater(const index_m & v1, const index_m & v2);
  bool i_total_revlex_is_greater(const index_m & v1, const index_m & v2);
  bool i_total_lex_is_greater(const index_m & v1, const index_m & v2);

  template <class T> T pow(const std::vector<T> & x, const index_m & n ){
    assert(x.size()==n.iptr->size());
    typename std::vector<T>::const_iterator itx=x.begin();
    index_t::const_iterator itn=n.iptr->begin();
    T res(1);
    for (;itx!=x.end();++itx,++itn){
      res=res*pow(*itx,*itn);
    }
    return res;
  }

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC

#endif // ndef _GIAC_INDEX_H_
