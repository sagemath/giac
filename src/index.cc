// -*- mode: C++ ; compile-command: "g++ -I.. -g -O2 -c index.cc" -*-
#include "first.h"
/*
 *  Copyright (C) 2000,7 B. Parisse, Institut Fourier, 38402 St Martin d'Heres
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
#include "index.h"
#include <cmath>
#include <stdio.h>
#include <stdexcept>

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  void setsizeerr(const std::string & s);

  int mygcd(int a,int b){
    if (b)
      return mygcd(b,a%b);
    else
      return a<0?-a:a;
  }

  index_t index_gcd(const index_t & a,const index_t & b){
    index_t::const_iterator ita=a.begin(),itaend=a.end(),itb=b.begin();
    unsigned s=itaend-ita;
    index_t res(s);
    index_t::iterator itres=res.begin();  
#ifdef DEBUG_SUPPORT
    if (s!=b.size())
      setsizeerr("Error index.cc index_gcd");
#endif // DEBUG_SUPPORT
    for (;ita!=itaend;++itb,++itres,++ita)
      *itres=min(*ita,*itb);
    return res;
  }

  index_t index_lcm(const index_t & a,const index_t & b){
    index_t::const_iterator ita=a.begin(),itaend=a.end(),itb=b.begin();
    unsigned s=itaend-ita;
    index_t res(s);
    index_t::iterator itres=res.begin();  
#ifdef DEBUG_SUPPORT
    if (s!=b.size())
      setsizeerr("index.cc index_lcm");
#endif // DEBUG_SUPPORT
    for (;ita!=itaend;++itb,++itres,++ita)
      *itres=max(*ita,*itb);
    return res;
  }

  // index and monomial ordering/operations implementation

  index_t operator + (const index_t & a, const index_t & b){
    index_t::const_iterator ita=a.begin(),itaend=a.end(),itb=b.begin();
    unsigned s=itaend-ita;
    index_t res(s);
    index_t::iterator itres=res.begin();  
#ifdef DEBUG_SUPPORT
    if (s!=b.size())
      setsizeerr("index.cc operator +");
#endif // DEBUG_SUPPORT
    for (;ita!=itaend;++itb,++itres,++ita)
      *itres=(*ita)+(*itb);
    return res;
  }

  index_t operator - (const index_t & a, const index_t & b){
    index_t res;
    index_t::const_iterator ita=a.begin(),itaend=a.end(),itb=b.begin();
    unsigned s=itaend-ita;
#ifdef DEBUG_SUPPORT
    if (s!=b.size())
      setsizeerr("index.cc operator -");
#endif // DEBUG_SUPPORT
    res.reserve(s);
    for (;ita!=itaend;++ita,++itb)
      res.push_back((*ita)-(*itb));
    return res;
  }

  index_t operator | (const index_t & a, const index_t & b){
    index_t res;
    index_t::const_iterator ita=a.begin(),itaend=a.end(),itb=b.begin();
    unsigned s=itaend-ita;
#ifdef DEBUG_SUPPORT
    if (s!=b.size())
      setsizeerr("index.cc operator |");
#endif // DEBUG_SUPPORT
    res.reserve(s);
    for (;ita!=itaend;++ita,++itb)
      res.push_back((*ita) | (*itb));
    return res;
  }

  index_t operator - (const index_t & a){
    index_t res;
    index_t::const_iterator ita=a.begin(),itaend=a.end();
    int s=itaend-ita;
    res.reserve(s);
    for (;ita!=itaend;++ita)
      res.push_back(-(*ita));
    return res;
  }

  index_t operator * (const index_t & a, int fois){
    index_t res;
    index_t::const_iterator ita=a.begin(),itaend=a.end();
    res.reserve(itaend-ita);
    for (;ita!=itaend;++ita)
      res.push_back((*ita)*fois);
    return res;
  }

  index_t operator / (const index_t & a, int divisepar){
    index_t res;
    index_t::const_iterator ita=a.begin(),itaend=a.end();
    res.reserve(itaend-ita);
    for (;ita!=itaend;++ita)
      res.push_back((*ita)/divisepar);
    return res;
  }

  int operator / (const index_t & a, const index_t & b){
    index_t::const_iterator ita=a.begin(),itaend=a.end(),itb=b.begin(),itbend=b.end();
#ifdef DEBUG_SUPPORT
    if (itaend-ita!=signed(b.size()))
      setsizeerr("index.cc operator /");
#endif // DEBUG_SUPPORT
    for (;ita!=itaend;++ita,++itb){
      if (*itb)
	return *ita / *itb;
    }
    return 0;
  }

  bool all_sup_equal (const index_t & a, const index_t & b){
    index_t::const_iterator ita=a.begin(),itaend=a.end(),itb=b.begin();
#ifdef DEBUG_SUPPORT
    if (itaend-ita!=signed(b.size()))
      setsizeerr("index.cc operator >=");
#endif // DEBUG_SUPPORT
    for (;ita!=itaend;++ita,++itb){
      if ((*ita)<(*itb))
	return false;
    }
    return true;
  }

  bool all_inf_equal (const index_t & a, const index_t & b){
    index_t::const_iterator ita=a.begin(),itaend=a.end(),itb=b.begin();
#ifdef DEBUG_SUPPORT
    if (itaend-ita!=signed(b.size()))
      setsizeerr("index.cc operator <=");
#endif // DEBUG_SUPPORT
    for (;ita!=itaend;++ita,++itb){
      if ((*ita)>(*itb))
	return false;
    }
    return true;
  }

  string print_INT_(int i){
    char c[256];
    sprintf(c,"%d",i);
    return c;
  }

  string hexa_print_INT_(int i){
    char c[256];
    sprintf(c,"%X",i);
    return string("0x")+c;
  }

  string octal_print_INT_(int i){
    char c[256];
    sprintf(c,"%o",i);
    return string("0o")+c;
  }


  string binary_print_INT_(int i){
    char c[256];
    mpz_t tmp;
    mpz_init_set_ui(tmp, i);
    mpz_get_str(c, 2, tmp);
    mpz_clear(tmp);
    return string("0b")+c;
  }

  /*
  string print_INT_(int i){
    if (!i)
      return string("0");
    if (i<0)
      return string("-")+print_INT_(-i);      
    int length = (int) std::floor(std::log10((double) i));
    char s[length+2];
    s[length+1]=0;
    for (;length>-1;--length,i/=10)
      s[length]=i%10+'0';
    return s;
  }
  */

  string print_INT_(const index_t & m){
    index_t::const_iterator it=m.begin(),itend=m.end();
    if (it==itend)
      return "";
    string s("[");
    for (;;){
      s += print_INT_(*it);
      ++it;
      if (it==itend)
	return s+']';
      else
	s += ',';
    }
  }
  
  ostream & operator << (ostream & os, const index_t & m ){
    return os << ":index_t: " << print_INT_(m) << " " ;
  }

  void dbgprint(const index_t & i){
    cout << i << endl;
  }

  index_t mergeindex(const index_t & i,const index_t & j){
    index_t res(i);
    index_t::const_iterator it=j.begin(),itend=j.end();
    res.reserve(i.size()+itend-it);
    for (;it!=itend;++it)
      res.push_back(*it);
    return res;
  }

  // by convention 0 -> 0 for permutations beginning at index 1
  index_t inverse(const index_t & p){
    index_t inv(p);
    int n=p.size();
    for (int i=0;i<n;i++){
      inv[p[i]]=i; // that's the definition of inv!!
    }
    return inv;
  }

  // transposition
  index_t transposition(int i,int j,int size){
    if (i>j)
      return transposition(j,i,size);
    index_t t;
    for (int k=0;k<i;k++)
      t.push_back(k);
    t.push_back(j);
    for (int k=i+1;k<j;k++)
      t.push_back(k);
    t.push_back(i);
    for (int k=j+1;k<size;k++)
      t.push_back(k);
    return t;
  }

  bool has(const index_t & p,int r){
    index_t::const_iterator it=p.begin(),itend=p.end();
    for (;it!=itend;++it){
      if (*it==r)
	return true;
    }
    return false;
  }

  bool is_zero(const index_t & p){
    index_t::const_iterator it=p.begin(),itend=p.end();
    for (;it!=itend;++it){
      if (*it)
	return false;
    }
    return true;
  }

  // WARNING GLOBAL VARIABLE
#ifdef HASH_MAP_NAMESPACE
  std::vector<hash_index> global_hash_index;
  int copy_number=0;
#endif
  
  index_m::index_m(const index_t & i){
#ifdef HASH_MAP_NAMESPACE
#ifdef GLOBAL_HASH_MAP
    // search if i is in the hash map
    // first look if the vector of hash map is large enough
    size_t t=i.size(),th=global_hash_index.size();
    for (;th<=t;++th){
      global_hash_index.push_back(hash_index());
    }
    hash_index & hptr=global_hash_index[t];
    hash_index::iterator it=hptr.find(i),itend=hptr.end();
    if (it!=itend){
      // it's in, make a copy 
      ++copy_number;
      iptr = it->second.iptr;
      ref_count=it->second.ref_count;
      (*ref_count)++;
    }
    else { // it's not in make a new one
#endif
#endif
      iptr=new index_t(i);
      ref_count=new int(1);
#ifdef HASH_MAP_NAMESPACE
#ifdef GLOBAL_HASH_MAP
      // and copy it
      hptr[i]=*this;
    }
#endif
#endif
  }
  
  index_m operator + (const index_m & a, const index_m & b){
    index_t::const_iterator ita=a.iptr->begin();
    index_t::const_iterator itaend=a.iptr->end();
    index_t::const_iterator itb=b.iptr->begin();
    int s=itaend-ita;
#ifdef DEBUG_SUPPORT
    if (s!=signed(b.iptr->size()))
      setsizeerr("index.cc index_m operator +");
#endif // DEBUG_SUPPORT
#ifdef GLOBAL_HASH_MAP
    index_t res(s);
    index_t::iterator it=res.begin();
    for (;ita!=itaend;++ita,++itb,++it)
      (*it) = (*ita)+(*itb);
    return res;
#else
    index_m res(s);
    index_t::iterator it=res.iptr->begin();
    for (;ita!=itaend;++it,++itb,++ita)
      *it = (*ita)+(*itb);
    return res;
#endif
  }

  index_m operator - (const index_m & a, const index_m & b){
    index_t::const_iterator ita=a.iptr->begin();
    index_t::const_iterator itaend=a.iptr->end();
    index_t::const_iterator itb=b.iptr->begin();
    int s=itaend-ita;
#ifdef DEBUG_SUPPORT
    if (s!=signed(b.iptr->size()))
      setsizeerr("index.cc index_m operator -");
#endif // DEBUG_SUPPORT
#ifdef GLOBAL_HASH_MAP
    index_t res(s);
    index_t::iterator it=res.begin();
    for (;ita!=itaend;++ita,++itb,++it)
      (*it) = (*ita)-(*itb);
    return res;
#else
    index_m res(s);
    index_t::iterator it=res.iptr->begin();
    for (;ita!=itaend;++it,++itb,++ita)
      *it = (*ita)-(*itb);
    return res;
#endif
  }

  index_m operator * (const index_m & a, int fois){
    index_t::const_iterator ita=a.iptr->begin(),itaend=a.iptr->end();
#ifdef GLOBAL_HASH_MAP
    index_t res(itaend-ita);
    index_t::iterator it=res.begin();
    for (;ita!=itaend;++it,++ita)
      *it= (*ita)*fois;
    return res;
#else
    index_m res(itaend-ita);
    index_t::iterator it=res.iptr->begin();
    for (;ita!=itaend;++it,++ita)
      *it = (*ita)*fois;
    return res;
#endif
  }

  index_m operator / (const index_m & a, int divisepar){
    index_t::const_iterator ita=a.iptr->begin(),itaend=a.iptr->end();
#ifdef GLOBAL_HASH_MAP
    index_t res(itaend-ita);
    index_t::iterator it=res.begin();
    for (;ita!=itaend;++it,++ita)
      *it= (*ita)/divisepar;
    return res;
#else
    index_m res(itaend-ita);
    index_t::iterator it=res.iptr->begin();
    for (;ita!=itaend;++it,++ita)
      *it = (*ita)/divisepar;
    return res;
#endif
  }

  bool operator == (const index_m & i1, const index_m & i2){
    if (i1.iptr==i2.iptr)
      return true;
    return (*(i1.iptr)==*(i2.iptr));
  }

  bool operator != (const index_m & i1, const index_m & i2){
    return !(i1==i2);
  }

  // >= and <= are *partial* ordering on index_t
  // they return TRUE if and only if >= or <= is true for *all* coordinates
  bool operator >= (const index_m & a, const index_m & b){
    index_t::const_iterator ita=a.iptr->begin(),itaend=a.iptr->end();
    index_t::const_iterator itb=b.iptr->begin();
#ifdef DEBUG_SUPPORT
    if (itaend-ita!=signed(b.iptr->size()))
      setsizeerr("index.cc index_m operator >=");
#endif
    for (;ita!=itaend;++ita,++itb){
      if ((*ita)<(*itb))
	return false;
    }
    return true;
  }

  bool operator <= (const index_m & a, const index_m & b){
    index_t::const_iterator ita=a.iptr->begin(),itaend=a.iptr->end();
    index_t::const_iterator itb=b.iptr->begin();
#ifdef DEBUG_SUPPORT
    if (itaend-ita!=signed(b.iptr->size()))
      setsizeerr("index.cc index_m operator >=");
#endif
    for (;ita!=itaend;++ita,++itb){
      if ((*ita)>(*itb))
	return false;
    }
    return true;
  }


  int total_degree(const index_m & v1){
    int i=0;
    for (index_t::const_iterator it=v1.iptr->begin();it!=v1.iptr->end();++it)
      i=i+(*it);
    return(i);
  }


  bool i_lex_is_greater(const index_m & v1, const index_m & v2){
    index_t::const_iterator it1=v1.iptr->begin();
    index_t::const_iterator it2=v2.iptr->begin();
    index_t::const_iterator it1end=v1.iptr->end();
#ifdef DEBUG_SUPPORT
    if (it1end-it1!=signed(v2.iptr->size()))
      setsizeerr("index.cc index_m i_lex_is_greater");
#endif
    for (;it1!=it1end;++it1){
      if ( (*it1)!=(*it2) ){
	if  ( (*it1)>(*it2))
	  return(true);
	else
	  return(false);
      }
      ++it2;
    }
    return(true);
  }

  bool i_lex_is_strictly_greater(const index_m & v1, const index_m & v2){
    index_t::const_iterator it1=v1.iptr->begin();
    index_t::const_iterator it2=v2.iptr->begin();
    index_t::const_iterator it1end=v1.iptr->end();
#ifdef DEBUG_SUPPORT
    if (it1end-it1!=signed(v2.iptr->size()))
      setsizeerr("index.cc index_m i_lex_is_greater");
#endif
    for (;it1!=it1end;++it1){
      if ( (*it1)!=(*it2) ){
	if  ( (*it1)>(*it2))
	  return(true);
	else
	  return(false);
      }
      ++it2;
    }
    return(false);
  }

  /*
  bool i_revlex_is_greater(const index_m & v1, const index_m & v2){
    return revlex_is_greater(*v1.iptr,*v2.iptr);
  }
  */

  bool i_total_lex_is_greater(const index_m & v1, const index_m & v2){
    int d1=total_degree(v1);
    int d2=total_degree(v2);
    if (d1!=d2){
      if (d1>d2)
	return(true);
      else
	return(false);
    }
    return(i_lex_is_greater(v1,v2));
  }

  bool i_total_revlex_is_greater(const index_m & v1, const index_m & v2){
    int d1=total_degree(v1);
    int d2=total_degree(v2);
    if (d1!=d2){
      if (d1>d2)
	return(true);
      else
	return(false);
    }
    return !i_lex_is_strictly_greater(v1,v2);
  }

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC
