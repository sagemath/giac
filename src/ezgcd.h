/* -*- mode:C++ ; compile-command: "g++ -I.. -g -c ezgcd.cc" -*- */
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

#ifndef _GIAC_EZGCD_H_
#define _GIAC_EZGCD_H_
#include "first.h"
#include "gausspol.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  void change_dim(polynome & p,int dim);

  // Hensel quadratic lift
  // Lift the equality p(b)=qb*rb [where b is a vecteur like for peval
  // assumed to have p.dim-1 coordinates] to p=q*r mod (X-b)^deg
  // Assuming that lcoeff(q)=lcp, lcoeff(r)=lcp, lcoeff(p)=lcp^2
  // If you want to find factors of a poly P such that P(b)=Qb*Rb, 
  // if lcp is the leading coeff of P
  // then p=P*lcp, qb=Qb*lcp(b)/lcoeff(Qb), rb=Rb*lcp(b)/lcoeff(Rb)
  bool hensel_lift(const polynome & p, const polynome & lcp, const polynome & qb, const polynome & rb, const vecteur & b,polynome & q, polynome & r,bool linear_lift=true,double maxop=-1);

  bool find_good_eval(const polynome & F,const polynome & G,polynome & Fb,polynome & Gb,vecteur & b,bool debuglog=false,const gen & mod=zero);
  polynome peval_1(const polynome & p,const vecteur &v,const gen & mod);

  // max_gcddeg is used when ezgcd was not successfull to find
  // the gcd even with 2 evaluations leading to the same gcd degree
  // in this case ezgcd calls itself with a bound on the gcd degree
  // is_sqff is true if we know that F_orig or G_orig is squarefree
  // is_primitive is true if F_orig and G_orig is primitive
  bool ezgcd(const polynome & F_orig,const polynome & G_orig,polynome & GCD,bool is_sqff=false,bool is_primitive=false,int max_gcddeg=0,double maxop=-1);

  extern const std::string _ezgcd_s;
  gen _ezgcd(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_ezgcd;  

  extern const std::string _modgcd_s;
  gen _modgcd(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_modgcd;  

  extern const std::string _heugcd_s;
  gen _heugcd(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_heugcd;  

  extern const std::string _psrgcd_s;
  gen _psrgcd(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_psrgcd;  

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // NO_NAMESPACE_GIAC

#endif // _GIAC_EZGCD_H
