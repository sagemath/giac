// -*- mode:C++ ; compile-command: "g++ -I.. -g -c subst.cc" -*-
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
#ifndef _GIAC_SUBST_H
#define _GIAC_SUBST_H
#include "first.h"
#include "gen.h"
#include "identificateur.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC
  polynome gen2poly(const gen & g,int s); 
  void checkanglemode(GIAC_CONTEXT);
  gen degtorad(const gen & g,GIAC_CONTEXT);
  gen radtodeg(const gen & g,GIAC_CONTEXT);
  // find symbolic vars in g that have u has sommet
  vecteur lop(const gen & g,const unary_function_ptr & u);
  vecteur lop(const gen & g,const std::vector< unary_function_ptr > & v);
  // One substitution
  vecteur subst(const vecteur & v,const gen & i,const gen & newi,bool quotesubst,GIAC_CONTEXT);
  gen subst(const gen & e,const gen & i,const gen & newi,bool quotesubst,GIAC_CONTEXT);
  sparse_poly1 subst(const sparse_poly1 & v,const gen & i,const gen & newi,bool quotesubst,GIAC_CONTEXT);
  // Multi substitutions
  vecteur subst(const vecteur & v,const vecteur & i,const vecteur & newi,bool quotesubst,GIAC_CONTEXT);
  gen subst(const gen & e,const vecteur & i,const vecteur & ewi,bool quotesubst,GIAC_CONTEXT);
  vecteur sortsubst(const vecteur & v,const vecteur & i,const vecteur & newi,bool quotesubst,GIAC_CONTEXT); // assumes that i is sorted
  gen sortsubst(const gen & e,const vecteur & i,const vecteur & newi,bool quotesubst,GIAC_CONTEXT); // assumes that i is sorted

  gen quotesubst(const gen & e,const gen & i,const gen & newi,GIAC_CONTEXT);
  gen gen_feuille(const gen & g);

  template<class T>
    int equalposcomp(const std::vector<T> & v, const T & w){
    int n=1;
    for (typename std::vector<T>::const_iterator it=v.begin();it!=v.end();++it){
      if ((*it)==w)
	return n;
      else
	n++;
    }
    return 0;
  }
  gen subst(const gen & e,const std::vector<unary_function_ptr> & v,const std::vector< gen (*) (const gen &,GIAC_CONTEXT) > & w,bool quotesubst,GIAC_CONTEXT);

  gen subst(const gen & e,const std::vector<unary_function_ptr> & v,const std::vector< gen (*) (const gen &) > & w,bool quotesubst,GIAC_CONTEXT);

  gen halftan(const gen & e,GIAC_CONTEXT);
  extern const std::string _halftan_s;
  gen _halftan(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_halftan;

  gen hyp2exp(const gen & e,GIAC_CONTEXT);
  extern const std::string _hyp2exp_s;
  gen _hyp2exp(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_hyp2exp;

  gen sincos(const gen & e,GIAC_CONTEXT);
  extern const std::string _sincos_s;
  gen _sincos(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_sincos;

  gen trig2exp(const gen & e,GIAC_CONTEXT);
  extern const std::string _trig2exp_s;
  gen _trig2exp(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_trig2exp;

  gen halftan_hyp2exp(const gen & e,GIAC_CONTEXT);
  extern const std::string _halftan_hyp2exp_s;
  gen _halftan_hyp2exp(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_halftan_hyp2exp;

  gen rewrite_hyper(const gen & e,GIAC_CONTEXT);

  gen asin2acos(const gen & e,GIAC_CONTEXT);
  extern const std::string _asin2acos_s;
  gen _asin2acos(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_asin2acos;

  gen asin2atan(const gen & e,GIAC_CONTEXT);
  extern const std::string _asin2atan_s;
  gen _asin2atan(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_asin2atan;

  gen acos2asin(const gen & e,GIAC_CONTEXT);
  extern const std::string _acos2asin_s;
  gen _acos2asin(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_acos2asin;

  gen acos2atan(const gen & e,GIAC_CONTEXT);
  extern const std::string _acos2atan_s;
  gen _acos2atan(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_acos2atan;

  gen atan2acos(const gen & e,GIAC_CONTEXT);
  extern const std::string _atan2acos_s;
  gen _atan2acos(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_atan2acos;

  extern unary_function_ptr at_atrig2ln;

  gen atan2asin(const gen & e,GIAC_CONTEXT);
  extern const std::string _atan2asin_s;
  gen _atan2asin(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_atan2asin;

  // rewrite vars of e in terms of exp/ln if s1 resp. s2 is > 1
  // and simplify
  gen tsimplify_noexpln(const gen & e,int s1,int s2,GIAC_CONTEXT);
  gen tsimplify_common(const gen & e,GIAC_CONTEXT);

  gen tsimplify(const gen & e,GIAC_CONTEXT);
  extern const std::string _tsimplify_s;
  gen _tsimplify(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_tsimplify;

  gen simplify(const gen & e,GIAC_CONTEXT);
  extern const std::string _simplify_s;
  gen _simplify(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_simplify;

  gen trigcos(const gen & e,GIAC_CONTEXT);
  extern const std::string _trigcos_s;
  gen _trigcos(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_trigcos;

  gen trigsin(const gen & e,GIAC_CONTEXT);
  extern const std::string _trigsin_s;
  gen _trigsin(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_trigsin;

  gen trigtan(const gen & e,GIAC_CONTEXT);
  extern const std::string _trigtan_s;
  gen _trigtan(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_trigtan;

  gen tan2sincos(const gen & e,GIAC_CONTEXT);
  extern const std::string _tan2sincos_s;
  gen _tan2sincos(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_tan2sincos;

  gen tan2sincos2(const gen & e,GIAC_CONTEXT);
  extern const std::string _tan2sincos2_s;
  gen _tan2sincos2(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_tan2sincos2;

  gen tan2cossin2(const gen & e,GIAC_CONTEXT);
  extern const std::string _tan2cossin2_s;
  gen _tan2cossin2(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_tan2cossin2;

  gen tcollect(const gen & e,GIAC_CONTEXT);
  extern const std::string _tcollect_s;
  gen _tcollect(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_tcollect;

  gen lncollect(const gen & e,GIAC_CONTEXT);
  extern const std::string _lncollect_s;
  gen _lncollect(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_lncollect;

  gen _powexpand(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_powexpand;

  gen _exp2pow(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_exp2pow;
  gen _pow2exp(const gen & e,GIAC_CONTEXT);

  gen pow2expln(const gen & e,const identificateur & x,GIAC_CONTEXT);
  gen gamma2factorial(const gen & g,GIAC_CONTEXT);
  gen gammatofactorial(const gen & g,GIAC_CONTEXT);
  gen factorial2gamma(const gen & g,GIAC_CONTEXT);
  gen factorialtogamma(const gen & g,GIAC_CONTEXT);

  gen _factor_xn(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_factor_xn;

  gen sin2tan2(const gen & e,GIAC_CONTEXT);
  gen cos2tan2(const gen & e,GIAC_CONTEXT);
  gen tan2tan2(const gen & e,GIAC_CONTEXT);
  gen sinh2exp(const gen & e,GIAC_CONTEXT);
  gen cosh2exp(const gen & e,GIAC_CONTEXT);
  gen tanh2exp(const gen & e,GIAC_CONTEXT);
  gen sin2exp(const gen & e,GIAC_CONTEXT);
  gen cos2exp(const gen & e,GIAC_CONTEXT);
  gen tan2exp(const gen & e,GIAC_CONTEXT);
  gen asin2ln(const gen & e,GIAC_CONTEXT);
  gen acos2ln(const gen & e,GIAC_CONTEXT);
  gen atan2ln(const gen & e,GIAC_CONTEXT);
  gen exp2sincos(const gen & e,GIAC_CONTEXT);
  gen tantosincos(const gen & e,GIAC_CONTEXT);
  gen tantosincos2(const gen & e,GIAC_CONTEXT);
  gen tantocossin2(const gen & e,GIAC_CONTEXT);
  gen asintoacos(const gen & e,GIAC_CONTEXT);
  gen asintoatan(const gen & e,GIAC_CONTEXT);
  gen acostoasin(const gen & e,GIAC_CONTEXT);
  gen acostoatan(const gen & e,GIAC_CONTEXT);
  gen atantoasin(const gen & e,GIAC_CONTEXT);
  gen atantoacos(const gen & e,GIAC_CONTEXT);
  gen trigcospow(const gen & e,GIAC_CONTEXT);
  gen trigsinpow(const gen & e,GIAC_CONTEXT);
  gen trigtanpow(const gen & e,GIAC_CONTEXT);
  gen powtopowexpand(const gen & e,GIAC_CONTEXT);
  gen exptopower(const gen & e,GIAC_CONTEXT);
  gen Heavisidetosign(const gen & args,GIAC_CONTEXT);
  gen _Heavisidetosign(const gen & args,GIAC_CONTEXT);

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC

#endif // _GIAC_SUBST_H
