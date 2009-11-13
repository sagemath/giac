// -*- mode:C++ ; compile-command: "g++-3.4 -I.. -g -c series.cc" -*-
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
#ifndef _GIAC_SERIES_H
#define _GIAC_SERIES_H
#include "first.h"

/* User representation of series expansion:
   - the power series expansion as a symbolic  + O(the remainder)
   Inner representation:
   - an expansion (generalized) variable h
   - a sparse_poly1
*/

#include "gen.h"
#include "identificateur.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC
  

  bool taylor(const gen & f_x,const gen & x,const gen & lim_point,int ordre,vecteur & v,GIAC_CONTEXT);
  // The remainder term in sparse_poly1 is the term defined with coeff==undef
  // This term must be the last term in the std::vector<monome>
  gen porder(const sparse_poly1 & a); // return plus_inf if a is exact
  void padd(const sparse_poly1 & a,const sparse_poly1 &b, sparse_poly1 & res,GIAC_CONTEXT);
  sparse_poly1 spadd(const sparse_poly1 & a,const sparse_poly1 &b,GIAC_CONTEXT);
  sparse_poly1 spsub(const sparse_poly1 & a,const sparse_poly1 &b,GIAC_CONTEXT);
  void pmul(const sparse_poly1 & a,const gen & b, sparse_poly1 & res,GIAC_CONTEXT);
  void pmul(const gen & b, const sparse_poly1 & a,sparse_poly1 & res,GIAC_CONTEXT);
  void pmul(const sparse_poly1 & a,const sparse_poly1 &b, sparse_poly1 & res,bool n_truncate,const gen & n_valuation,GIAC_CONTEXT);
  sparse_poly1 spmul(const sparse_poly1 & a,const sparse_poly1 &b,GIAC_CONTEXT);
  sparse_poly1 spmul(const gen & a,const sparse_poly1 &b,GIAC_CONTEXT);
  sparse_poly1 spmul(const sparse_poly1 & a,const gen &b,GIAC_CONTEXT);
  void pneg(const sparse_poly1 & a,sparse_poly1 & res,GIAC_CONTEXT);
  sparse_poly1 spneg(const sparse_poly1 & a,GIAC_CONTEXT);
  void pshift(const sparse_poly1 & a,const gen & b, sparse_poly1 & res,GIAC_CONTEXT);
  void pdiv(const sparse_poly1 & a,const sparse_poly1 &b, sparse_poly1 & res,int ordre,GIAC_CONTEXT);
  void pdiv(const sparse_poly1 & a,const gen & b, sparse_poly1 & res,GIAC_CONTEXT);
  sparse_poly1 spdiv(const sparse_poly1 & a,const sparse_poly1 &b,GIAC_CONTEXT);  
  sparse_poly1 spdiv(const sparse_poly1 & a,const gen &b,GIAC_CONTEXT);  
  void ppow(const sparse_poly1 & base,int m,int ordre,sparse_poly1 & res,GIAC_CONTEXT); // m>=0
  void ppow(const sparse_poly1 & base,const gen & e,int ordre,int direction,sparse_poly1 & res,GIAC_CONTEXT);
  void pcompose(const vecteur & v,const sparse_poly1 & p, sparse_poly1 & res,GIAC_CONTEXT);
  void lcmdeno(vecteur &v,gen & e,GIAC_CONTEXT);
  void lcmdeno(sparse_poly1 &v,gen & e,GIAC_CONTEXT);
  
  class unary_function_ptr;
  // main series expansion C++ entry point
  gen series(const gen & e,const identificateur & x,const gen & lim_point,int ordre,int direction,GIAC_CONTEXT);
  // other series entry points, used for interactive input
  gen series(const gen & e,const gen & vars,const gen & lim_point,int ordre,int direction,GIAC_CONTEXT);
  gen series(const gen & e,const gen & vars,const gen & lim_point,const gen &ordre,GIAC_CONTEXT);
  gen _series(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_series ;

  bool contains(const gen & e,const gen & elem);
  vecteur lvarx(const gen &e,const identificateur & x,bool test=false);
  void rlvarx(const gen &e,const identificateur & x,vecteur & res);
  vecteur rlvarx(const gen &e,const identificateur & x);
  bool intersect(const vecteur & a,const vecteur &b,int & pos_a,int & pos_b);
  // void mrv(const gen & e,const identificateur & x,vecteur & faster_var,vecteur & coeff_ln, vecteur & slower_var,GIAC_CONTEXT){
  // void mrv_lead_term(const gen & e,const identificateur & x,gen & coeff, gen & mrv_var, gen & exponent,sparse_poly1 & q,int begin_ordre,GIAC_CONTEXT);

  gen limit(const gen & e,const identificateur & x,const gen & lim_point,int direction,GIAC_CONTEXT);
  gen _limit(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_limit ;

  gen _bounded_function(const gen & args,GIAC_CONTEXT);
  gen bounded_function(GIAC_CONTEXT);
  extern unary_function_ptr at_bounded_function ;

  // expansion of e at x=lim_point, order ordre, direction 0/1/-1
  // answer is in s
  void series__SPOL1(const gen & e,const identificateur & x,const gen & lim_point,int ordre,int direction,sparse_poly1 & s,GIAC_CONTEXT);
  sparse_poly1 series__SPOL1(const gen & e,const identificateur & x,const gen & lim_point,int ordre,int direction,GIAC_CONTEXT);

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC

#endif // _GIAC_SERIES_H
