/* -*- mode:C++ ; compile-command: "g++-3.4 -I.. -I../include -g -c threaded.cc" -*- */
#include "first.h"
/*  Copyright (C) 2000,2007 B. Parisse, Institut Fourier, 38402 St Martin d'Heres
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
#include "threaded.h"
#include "sym2poly.h"
#include "gausspol.h"
#include "time.h"
#include "usual.h"

#ifdef HAVE_GMPXX_H
mpz_class invmod(const mpz_class & a,int reduce){
  mpz_class z=(a%reduce);
  int tmp=z.get_si();
  tmp=giac::invmod(tmp,reduce);
  return tmp;
}

mpz_class smod(const mpz_class & a,int reduce){
  mpz_class z=(a%reduce);
  int tmp=z.get_si();
  tmp=giac::smod(tmp,reduce);
  return tmp;
}
#endif

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  int heap_mult=0;
  gen _heap_mult(const gen & g,GIAC_CONTEXT){
    if (g.type!=_INT_)
      return heap_mult;
    return heap_mult=g.val;
  }
  const string _heap_mult_s("heap_mult");
  unary_function_eval __heap_mult(&_heap_mult,_heap_mult_s);
  unary_function_ptr at_heap_mult (&__heap_mult,0,true);

  my_mpz operator % (const my_mpz & a,const my_mpz & b){
    my_mpz tmp;
    mpz_fdiv_r(tmp.ptr,a.ptr,b.ptr);
    return tmp;
  }

  my_mpz operator %= (my_mpz & a,const my_mpz & b){
    mpz_fdiv_r(a.ptr,a.ptr,b.ptr);
    return a;
  }

  my_mpz operator += (my_mpz & a,const my_mpz & b){
    mpz_add(a.ptr,a.ptr,b.ptr);
    return a;
  }

  my_mpz operator -= (my_mpz & a,const my_mpz & b){
    mpz_sub(a.ptr,a.ptr,b.ptr);
    return a;
  }

  my_mpz operator + (const my_mpz & a,const my_mpz & b){
    my_mpz tmp;
    mpz_add(tmp.ptr,a.ptr,b.ptr);
    return tmp;
  }

  my_mpz operator - (const my_mpz & a,const my_mpz & b){
    my_mpz tmp;
    mpz_sub(tmp.ptr,a.ptr,b.ptr);
    return tmp;
  }

  my_mpz operator * (const my_mpz & a,const my_mpz & b){
    my_mpz tmp;
    mpz_mul(tmp.ptr,a.ptr,b.ptr);
    return tmp;
  }

  my_mpz operator / (const my_mpz & a,const my_mpz & b){
    my_mpz tmp;
    mpz_fdiv_q(tmp.ptr,a.ptr,b.ptr);
    return tmp;
  }

  my_mpz invmod(const my_mpz & a,int reduce){
    my_mpz z=(a%reduce);
    int tmp=mpz_get_si(z.ptr);
    tmp=invmod(tmp,reduce);
    return tmp;
  }

  my_mpz smod(const my_mpz & a,int reduce){
    my_mpz z=(a%reduce);
    int tmp=mpz_get_si(z.ptr);
    tmp=smod(tmp,reduce);
    return tmp;
  }

  void wait_1ms(context * contextptr){
    usleep(1000);
  }

  void background_callback(const gen & g,void * newcontextptr){
    if (g.type==_VECT && g._VECTptr->size()==2){
      context * cptr=(giac::context *)newcontextptr;
      if (cptr){
#ifdef HAVE_LIBPTHREAD
	pthread_mutex_lock(cptr->globalptr->_mutex_eval_status_ptr);
	sto(g._VECTptr->back(),g._VECTptr->front(),cptr);
	pthread_mutex_unlock(cptr->globalptr->_mutex_eval_status_ptr);
#endif
      }
    }
  }

#ifdef HAVE_LIBPTHREAD
  gen _background(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT || g._VECTptr->size()<2)
      setsizeerr();
    vecteur v(*g._VECTptr);
    gen target=v[0];
    gen toeval=v[1];
    int s=v.size();
    double maxsize=1e9; // in bytes
    double maxtime=1e9; // in microseconds
    int level=eval_level(contextptr);
    if (s>2){
      gen tmp=evalf_double(v[2],level,contextptr);
      if (tmp.type!=_DOUBLE_)
	settypeerr();
      maxsize=tmp._DOUBLE_val;
    }
    if (s>3){
      gen tmp=evalf_double(v[3],level,contextptr);
      if (tmp.type!=_DOUBLE_)
	settypeerr();
      maxtime=tmp._DOUBLE_val;
    }
    if (s>4 && v[4].type==_INT_)
      level=v[4].val;
    gen tmp;
    context * newcontextptr=clone_context(contextptr);
    newcontextptr->parent=contextptr;
    tmp=gen(newcontextptr,_THREAD_POINTER);
    sto(tmp,target,contextptr);
    if (!make_thread(makevecteur(symbolic(at_quote,target),toeval),level,background_callback,(void *)newcontextptr,newcontextptr)){
      sto(undef,target,contextptr);
      setsizeerr("Unable to make thread");
    }
    return tmp;
  }
  const string _background_s("background");
  unary_function_eval __background(&_background,_background_s);
  unary_function_ptr at_background (&__background,_QUOTE_ARGUMENTS,true);
#endif
  
#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC
