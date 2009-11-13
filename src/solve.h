// -*- mode:C++ ; compile-command: "g++ -I.. -g -c solve.cc" -*-
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
#ifndef _GIAC_SOLVE_H
#define _GIAC_SOLVE_H
#include "first.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC
  extern int intvar_counter;
  extern int realvar_counter;
  gen equal2diff(const gen & g); // rewrite = as -
  vecteur protect_sort(const vecteur & res,GIAC_CONTEXT);
  vecteur find_singularities(const gen & e,const identificateur & x,int cplxmode,GIAC_CONTEXT);
  // isolate_mode & 1 is complex_mode, isolate_mode & 2 is 0 for principal sol
  vecteur solve(const gen & e,const identificateur & x,int isolate_mode,GIAC_CONTEXT);
  vecteur solve(const gen & e,const gen & x,int isolate_mode,GIAC_CONTEXT);
  vecteur solve(const vecteur & v,bool complex_mode,GIAC_CONTEXT); // v is a 1-d dense polynomial
  void solve(const gen & e,const identificateur & x,vecteur &v,int isolate_mode,GIAC_CONTEXT);
  void in_solve(const gen & e,const identificateur & x,vecteur &v,int isolate_mode,GIAC_CONTEXT);
  gen solvepostprocess(const gen & g,const gen & x,GIAC_CONTEXT);
  // convert solutions to an expression
  gen _solve(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_solve ;
  gen _fsolve(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_fsolve ;
  vecteur sxa(const vecteur & sl,const vecteur & x,GIAC_CONTEXT);
  vecteur linsolve(const vecteur & sl,const vecteur & x,GIAC_CONTEXT);
  gen symb_linsolve(const gen & syst,const gen & vars);
  gen _linsolve(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_linsolve ;
  
  /*
  gen newtona(const gen & f, const gen & x, const gen & arg,int niter1, int niter2, double eps1,double eps2,double prefact1,double prefact2, int & b);
  gen newton(const gen & f, const gen & x,const gen & guess,int niter1=5,int niter2=50,double eps1=1e-3,double eps2=1e-12,double prefact1=0.5,double prefact2=1.0);
  */
  
  gen newton(const gen & f, const gen & x,const gen & guess,int niter,double eps1,double eps2,GIAC_CONTEXT);
  gen _newton(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_newton ;

  bool has_num_coeff(const vecteur & v);
  bool has_num_coeff(const polynome & p);
  bool has_num_coeff(const gen & e);
  bool has_mod_coeff(const vecteur & v,gen & modulo);
  bool has_mod_coeff(const polynome & p,gen & modulo);
  bool has_mod_coeff(const gen & e,gen & modulo);

  polynome spoly(const polynome & p,const polynome & q,environment * env);
  polynome reduce(const polynome & p,vectpoly::const_iterator it,vectpoly::const_iterator itend,environment * env);
  polynome reduce(const polynome & p,const vectpoly & v,environment * env);
  void reduce(vectpoly & res,environment * env);
  void change_monomial_order(polynome & p,const gen & order);
  vectpoly gbasis(const vectpoly & v,const gen & order=_TDEG_ORDER,bool with_cocoa=true,bool with_f5=false,environment * env=0);
  gen remove_equal(const gen & f);
  vecteur remove_equal(const_iterateur it,const_iterateur itend);
  vecteur gsolve(const vecteur & eq_orig,const vecteur & var,bool complexmode,GIAC_CONTEXT);

  gen _greduce(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_greduce ;

  gen _gbasis(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_gbasis ;

  double nan();
  gen remove_and(const gen & g,const unary_function_ptr & u);
  vecteur solvepreprocess(const gen & args,bool complex_mode,GIAC_CONTEXT);

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC

#endif // _GIAC_SOLVE_H
