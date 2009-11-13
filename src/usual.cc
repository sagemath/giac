// -*- mode:C++ ; compile-command: "g++-3.4 -I.. -I../include -g -c usual.cc -Wall" -*-
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
#include <stdexcept>
#include <cmath>
#include <cstdlib>
#include "gen.h"
#include "identificateur.h"
#include "symbolic.h"
#include "poly.h"
#include "usual.h"
#include "series.h"
#include "modpoly.h"
#include "sym2poly.h"
#include "moyal.h"
#include "subst.h"
#include "gausspol.h"
#include "identificateur.h"
#include "ifactor.h"
#include "prog.h"
#include "rpn.h"
#include "plot.h"
#include "pari.h"
#include "tex.h"
#include "unary.h"
#include "intg.h"
#include "ti89.h"
#include "solve.h"
#include "alg_ext.h"
#include "lin.h"
#include "derive.h"
#include "series.h"
#ifdef VISUALC
#include <float.h>
#endif
#ifdef HAVE_LIBGSL
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_expint.h>
#endif

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  // must be declared before any function declaration with special handling
  vector<unary_function_ptr *> limit_tractable_functions;
  vector<gen_op_context> limit_tractable_replace;
  string messages_to_print ;

  gen frac_neg_out(const gen & g,GIAC_CONTEXT){
    if (g.type==_FRAC && g._FRACptr->num.type<=_DOUBLE_ && is_strictly_positive(-g._FRACptr->num,contextptr))
      return symbolic(at_neg,-g);
    if (g.is_symb_of_sommet(at_prod)){
      // count neg
      gen f=g._SYMBptr->feuille;
      vecteur fv(gen2vecteur(f));
      int count=0,fvs=fv.size();
      for (int i=0;i<fvs;++i){
	fv[i]=frac_neg_out(fv[i],contextptr);
	if (fv[i].is_symb_of_sommet(at_neg)){
	  ++count;
	  fv[i]=fv[i]._SYMBptr->feuille;
	}
      }
      if (fvs==1)
	f=fv[0];
      else
	f=symbolic(at_prod,fv);
      if (count%2)
	return symbolic(at_neg,f);
      else
	return f;
    }
    return g;
  }

  // utilities for trig functions
  enum { trig_deno=24 };

  bool is_multiple_of_12(const gen & k0,int & l){
    if (!k0.is_integer())
      return false;
    gen k=smod(k0,trig_deno);
    if (k.type!=_INT_)
      return false;
    l=k.val+trig_deno/2;
    return true; 
  }
  bool is_multiple_of_pi_over_12(const gen & a,int & l,bool angle_unit,GIAC_CONTEXT){
    if (is_zero(a)){
      l=0;
      return true;
    }
    gen k;
    if (angle_unit){
      if (!contains(a,_IDNT_pi))
	return false;
      k=normal(rdiv(a*gen(trig_deno/2),cst_pi),contextptr);
    }
    else 
      k=rdiv(a,15);
    return is_multiple_of_12(k,l);
  }

  // 0 = -pi, 12=0, 24=pi
  gen * table_cos[trig_deno+1]={
    &minus_one,&minus_cos_pi_12,&minus_sqrt3_2,&minus_sqrt2_2,&minus_one_half,&minus_sin_pi_12,
    &zero,&sin_pi_12,&plus_one_half,&plus_sqrt2_2,&plus_sqrt3_2,&cos_pi_12,
    &plus_one,&cos_pi_12,&plus_sqrt3_2,&plus_sqrt2_2,&plus_one_half,&sin_pi_12,
    &zero,&minus_sin_pi_12,&minus_one_half,&minus_sqrt2_2,&minus_sqrt3_2,&minus_cos_pi_12,
    &minus_one
  };
  gen * table_tan[trig_deno/2+1]={
    &zero,&tan_pi_12,&plus_sqrt3_3,&plus_one,&plus_sqrt3,&tan_5pi_12,
    &unsigned_inf,&minus_tan_5pi_12,&minus_sqrt3,&minus_one,&minus_sqrt3_3,&minus_tan_pi_12,
    &zero
  };


  bool is_rational(const gen & a,int &n,int &d){
    gen num,den;
    fxnd(a,num,den);
    if (num.type!=_INT_ || den.type!=_INT_)
      return false;
    n=num.val;
    d=den.val;
    return true;
  }
  // checking utility
  void check_2d_vecteur(const gen & args) {
    if (args.type!=_VECT)
      settypeerr("check_2d_vecteur");
    if (args._VECTptr->size()!=2)
      setsizeerr("check_2d_vecteur");
  }

  // zero arg
  unary_function_constant __1(1);
  unary_function_ptr at_one (&__1);
  unary_function_constant __0(0);
  unary_function_ptr at_zero (&__0);
  
  gen _rm_a_z(const gen & args,GIAC_CONTEXT){
    if (variables_are_files(contextptr)){
      char a_effacer[]="a.cas";
      for (;a_effacer[0]<='z';++a_effacer[0]){
	unlink(a_effacer);
      }
    }
    _purge(list_one_letter__IDNT,contextptr);
    return args;
  }
  const string _rm_a_z_s("rm_a_z");
  unary_function_eval __rm_a_z(&_rm_a_z,_rm_a_z_s);
  unary_function_ptr at_rm_a_z (&__rm_a_z,0,true);

  gen _rm_all_vars(const gen & args,const context * contextptr){
    gen g=_VARS(zero,contextptr);
    if (g.type!=_VECT)
      return g;
    vecteur & v=*g._VECTptr;
    const_iterateur it=v.begin(),itend=v.end();
    for (;it!=itend;++it){
      if (it->type==_IDNT && (*it!=cst_pi) )
	_purge(*it,contextptr);
    }
    return g;
  }
  const string _rm_all_vars_s("rm_all_vars");
  unary_function_eval __rm_all_vars(&_rm_all_vars,_rm_all_vars_s);
  unary_function_ptr at_rm_all_vars (&__rm_all_vars,0,true);

  bool is_equal(const gen & g){
    return (g.type==_SYMB) && (g._SYMBptr->sommet==at_equal);
  }

  gen apply_to_equal(const gen & g,const gen_op & f){
    if (g.type!=_SYMB || g._SYMBptr->sommet!=at_equal || g._SYMBptr->feuille.type!=_VECT)
      return f(g);
    vecteur & v=*g._SYMBptr->feuille._VECTptr;
    if (v.empty())
      setsizeerr();
    return symbolic(at_equal,makevecteur(f(v.front()),f(v.back())));
  }

  gen apply_to_equal(const gen & g,gen (* f) (const gen &, GIAC_CONTEXT),GIAC_CONTEXT){
    if (g.type!=_SYMB || g._SYMBptr->sommet!=at_equal || g._SYMBptr->feuille.type!=_VECT)
      return f(g,contextptr);
    vecteur & v=*g._SYMBptr->feuille._VECTptr;
    if (v.empty())
      setsizeerr();
    return symbolic(at_equal,makevecteur(f(v.front(),contextptr),f(v.back(),contextptr)));
  }

  // one arg
  gen _id(const gen & args){
    return args;
  }
  extern partial_derivative_onearg D_at_id;
  const string _id_s("id");
  unary_function_unary __id(&_id,&D_at_id,_id_s);
  unary_function_ptr at_id (&__id,0,true);
  partial_derivative_onearg D_at_id(at_one);


  symbolic symb_not(const gen & args){
    return symbolic(at_not,args);
  }
  gen _not(const gen & args){
    if (args.type==_VECT)
      return apply(args,_not);
    return !equaltosame(args);
  }
  const string _not_s("not");
  unary_function_unary __not(&_not,_not_s);
  unary_function_ptr at_not (&__not);

  symbolic symb_neg(const gen & args){
    return symbolic(at_neg,args);
  }
  gen _neg(const gen & args){
    return -args;
  }
  extern partial_derivative_onearg D_at_neg;
  const string _neg_s("-");
  unary_function_unary __neg(&_neg,&D_at_neg,_neg_s);
  unary_function_ptr at_neg (&__neg,0,true);
  partial_derivative_onearg D_at_neg(at_neg);

  symbolic symb_inv(const gen & a){
    return symbolic(at_inv,a);
  }
  gen _inv(const gen & args,GIAC_CONTEXT){
    if ((args.type!=_VECT) || is_squarematrix(args))
      return inv(args,contextptr);
    iterateur it=args._VECTptr->begin(), itend=args._VECTptr->end();
    gen prod(1);
    for (;it!=itend;++it)
      prod = prod * (*it);
    return inv(prod,contextptr);
  }
  const string _inv_s("inv");
  unary_function_eval __inv(&_inv,_inv_s);
  unary_function_ptr at_inv (&__inv,0,true);

  symbolic symb_ln(const gen & e){
    return symbolic(at_ln,e);
  }

  gen ln(const gen & e,GIAC_CONTEXT){
    if (e.type==_DOUBLE_){
      if (e._DOUBLE_val==0)
	return minus_inf;
      if (e._DOUBLE_val>0)
#ifdef _SOFTMATH_H
	return std::giac_gnuwince_log(e._DOUBLE_val);
#else
	return std::log(e._DOUBLE_val);
#endif
      else
#ifdef _SOFTMATH_H
	return M_PI*cst_i+std::giac_gnuwince_log(-e._DOUBLE_val);
#else
	return M_PI*cst_i+std::log(-e._DOUBLE_val);
#endif
    }
    if (e.type==_REAL){
      if (is_positive(e,contextptr))
	return e._REALptr->log();
      else
	return (-e)._REALptr->log()+cst_pi*cst_i;
    }
    if (e.type==_CPLX){ 
      if (e.subtype)
#ifdef _SOFTMATH_H
	return std::giac_gnuwince_log(gen2complex_d(e));
#else
	return std::log(gen2complex_d(e));
#endif
      if (e._CPLXptr->type==_REAL)
	return ln(abs(e,contextptr),contextptr)+cst_i*arg(e,contextptr);
      if (is_zero(*e._CPLXptr)){
	if (is_one(*(e._CPLXptr+1)))
	  return cst_i*cst_pi_over_2;
	if (is_minus_one(*(e._CPLXptr+1)))
	  return -cst_i*cst_pi_over_2;
      }
    }
    if (is_squarematrix(e))
      return analytic_apply(at_ln,*e._VECTptr,contextptr);
    if (e.type==_VECT)
      return apply(e,giac::ln,contextptr);
    if (is_zero(e))
      return minus_inf;
    if (is_one(e))
      return 0;
    if (is_minus_one(e))
      return cst_i*cst_pi;
    if ( (e==undef) || (e==unsigned_inf) || (e==plus_inf))
      return e;
    if (is_equal(e))
      return apply_to_equal(e,ln,contextptr);
    if (e.type==_SYMB){
      if (e._SYMBptr->sommet==at_inv && e._SYMBptr->feuille.type!=_VECT)
	return -ln(e._SYMBptr->feuille,contextptr);
      if (e._SYMBptr->sommet==at_exp){ 
	if (is_real(e._SYMBptr->feuille,contextptr) ) 
	  return e._SYMBptr->feuille;
      }
    }
    if (e.is_symb_of_sommet(at_pow) && e._SYMBptr->feuille.type==_VECT && e._SYMBptr->feuille._VECTptr->size()==2){
      gen a=e._SYMBptr->feuille._VECTptr->front();
      gen b=e._SYMBptr->feuille._VECTptr->back();
      // ln(a^b)
      if (is_positive(a,contextptr))
	return b*ln(a,contextptr);
    }
    return symb_ln(e);
  }
  gen log(const gen & e,GIAC_CONTEXT){
    return ln(e,contextptr);
  }
  const string _ln_s("ln"); // Using C notation, log works also for natural
  gen d_ln(const gen & args,GIAC_CONTEXT){
    return inv(args,contextptr);
  }
  partial_derivative_onearg D_at_ln(&d_ln);
  unary_function_eval __ln(&giac::ln,&D_at_ln,_ln_s);
  unary_function_ptr at_ln (&__ln,0,true);

  gen log10(const gen & e,GIAC_CONTEXT){
    if (e.type==_DOUBLE_ && e._DOUBLE_val>0 ){
#ifdef _SOFTMATH_H
      return std::giac_gnuwince_log10(e._DOUBLE_val);
#else
      return std::log10(e._DOUBLE_val);
#endif
    }
    if ( e.type==_DOUBLE_ || (e.type==_CPLX && e.subtype)){
#ifdef _SOFTMATH_H
      return std::giac_gnuwince_log(gen2complex_d(e))/std::log(10.0);
#else
      return std::log(gen2complex_d(e))/std::log(10.0);
#endif
    }
    return rdiv(ln(e,contextptr),ln(10,contextptr));
  }
  const string _log10_s("log10"); // Using C notation, log for natural
  gen d_log10(const gen & args,GIAC_CONTEXT){
    return inv(args*ln(10,contextptr),contextptr);
  }
  partial_derivative_onearg D_at_log10(&d_log10);
  unary_function_eval __log10(&giac::log10,&D_at_log10,_log10_s);
  unary_function_ptr at_log10 (&__log10,0,true);

  gen alog10(const gen & e,GIAC_CONTEXT){
    if (is_squarematrix(e))
      return analytic_apply(at_alog10,*e._VECTptr,0);
    if (e.type==_VECT)
      return apply(e,contextptr,giac::alog10);
    if (is_equal(e))
      return apply_to_equal(e,alog10,contextptr);
    return pow(gen(10),e,contextptr);
  }
  const string _alog10_s("alog10"); 
  unary_function_eval __alog10(&giac::alog10,_alog10_s);
  unary_function_ptr at_alog10 (&__alog10,0,true);

  symbolic symb_atan(const gen & e){
    return symbolic(at_atan,e);
  }
  gen atanasln(const gen & e,GIAC_CONTEXT){
    return plus_one_half*cst_i*ln(rdiv(cst_i+e,cst_i-e),contextptr);
  }
  gen atan(const gen & e0,GIAC_CONTEXT){
    gen e=frac_neg_out(e0,contextptr);
    if (e.type==_DOUBLE_){
#ifdef _SOFTMATH_H
      double res=std::giac_gnuwince_atan(e._DOUBLE_val);
#else
      double res=std::atan(e._DOUBLE_val);
#endif
      if (angle_radian(contextptr)) 
	return res;
      else
	return res*rad2deg_d;
    }
    if (e.type==_REAL)
      return e._REALptr->atan();
    if ( (e.type==_CPLX) && (e.subtype || e._CPLXptr->type==_REAL)){
      if (angle_radian(contextptr)) 
	return no_context_evalf(atanasln(e,contextptr));
      else
	return no_context_evalf(atanasln(e,contextptr))*gen(rad2deg_d);
    }
    if (is_squarematrix(e))
      return analytic_apply(at_atan,*e._VECTptr,contextptr);
    if (e.type==_VECT)
      return apply(e,giac::atan,contextptr);
    if (is_zero(e))
      return e;
    if (is_one(e)){
      if (angle_radian(contextptr)) 
	return rdiv(cst_pi,4);
      return 45;
    }
    if (is_minus_one(e)){
      if (angle_radian(contextptr)) 
	return rdiv(-cst_pi,4);
      return -45;
    }
    if (e==plus_sqrt3_3){
      if (angle_radian(contextptr)) 
	return rdiv(cst_pi,6);
      return 30;
    }
    if (e==plus_sqrt3){
      if (angle_radian(contextptr)) 
	return rdiv(cst_pi,3);
      return 60;
    }
    if (e==plus_inf){
      if (angle_radian(contextptr)) 
	return cst_pi_over_2;
      return 90;
    }
    if (e==minus_inf){
      if (angle_radian(contextptr)) 
	return -cst_pi_over_2;
      return -90;
    }
    if ((e==undef)||(e==unsigned_inf))
      return undef;
    gen tmp=evalf_double(e,0,contextptr);
    if (tmp.type==_DOUBLE_){
      gen tmp2=normal(2*e/(1-e*e),contextptr);
      if (is_one(tmp2)){
	if (angle_radian(contextptr)) 
	  return tmp._DOUBLE_val>0?rdiv(cst_pi,8):rdiv(-3*cst_pi,8);
	return tmp._DOUBLE_val>0?fraction(45,2):fraction(-145,2);
      }
      if (is_minus_one(tmp2)){
	if (angle_radian(contextptr)) 
	  return tmp._DOUBLE_val>0?rdiv(3*cst_pi,8):rdiv(-cst_pi,8);
	return tmp._DOUBLE_val>0?fraction(145,2):fraction(-45,2);
      }
      if (tmp2==plus_sqrt3_3){
	if (angle_radian(contextptr)) 
	  return tmp._DOUBLE_val>0?rdiv(cst_pi,12):rdiv(-5*cst_pi,12);
	return tmp._DOUBLE_val>0?15:-165;
      }
      if (tmp2==rdiv(minus_sqrt3,3)){
	if (angle_radian(contextptr)) 
	  return tmp._DOUBLE_val>0?rdiv(5*cst_pi,12):rdiv(-cst_pi,12);
	return tmp._DOUBLE_val>0?165:-15;
      }
    }
    if ((e.type==_SYMB) && (e._SYMBptr->sommet==at_neg))
      return -atan(e._SYMBptr->feuille,contextptr);
    if ( (e.type==_INT_) && (e.val<0) )
      return -atan(-e,contextptr);
    if (is_equal(e))
      return apply_to_equal(e,atan,contextptr);
    if (e.is_symb_of_sommet(at_tan)){
      gen tmp=cst_pi;
      if (!angle_radian(contextptr))
	tmp=180;
      return e._SYMBptr->feuille-_floor(e._SYMBptr->feuille/tmp+plus_one_half,contextptr)*tmp;
    }
    return symb_atan(e);
  }
  gen d_atan(const gen & args,GIAC_CONTEXT){
    gen g=inv(1+pow(args,2),contextptr);
    if (angle_radian(contextptr))
      return g;
    return g*deg2rad_e;
  }
  partial_derivative_onearg D_at_atan(&d_atan);
  gen taylor_atan (const gen & lim_point,const int ordre,const unary_function_ptr & f, int direction,gen & shift_coeff,GIAC_CONTEXT){
    if (ordre<0)
      return 0; // no symbolic preprocessing
    shift_coeff=0;
    if (!is_inf(lim_point))
      return taylor(lim_point,ordre,f,0,shift_coeff,contextptr);
    vecteur v;
    identificateur x(" ");
    taylor(atan(x,contextptr),x,0,ordre,v,contextptr);
    v=negvecteur(v);
    v.front()=atan(lim_point,contextptr);
    return v;
  }
  const string _atan_s("atan");
  unary_function_eval __atan(&giac::atan,&D_at_atan,&taylor_atan,_atan_s);
  unary_function_ptr at_atan (&__atan,0,true);

  symbolic symb_exp(const gen & e){
    return symbolic(at_exp,e);
  }
  gen numeric_matrix_exp(const gen & e,double eps,GIAC_CONTEXT){
    gen res=midn(e._VECTptr->size());
    gen eee(e);
    for (double j=2;j<max_numexp && linfnorm(eee,contextptr)._DOUBLE_val>eps;++j){
      res = res + eee;
      eee = gen(1/j) * eee * e ; 
    }
    return res;
  }

  gen exp(const gen & e0,GIAC_CONTEXT){
    gen e=frac_neg_out(e0,contextptr);
    if (e.type==_SPOL1)
      return symb_exp(e);
    if (e.type==_DOUBLE_){
#ifdef _SOFTMATH_H
      return std::giac_gnuwince_exp(e._DOUBLE_val);
#else
      return std::exp(e._DOUBLE_val);
#endif
    }
    if (e.type==_REAL)
      return e._REALptr->exp();
    if (e.type==_CPLX){ 
      if (e.subtype){
#ifdef _SOFTMATH_H
	return std::giac_gnuwince_exp(gen2complex_d(e));
#else
	return std::exp(gen2complex_d(e));
#endif
      }
      if (e._CPLXptr->type==_REAL){
	return exp(*e._CPLXptr,contextptr)*(cos(*(e._CPLXptr+1),contextptr)+cst_i*sin(*(e._CPLXptr+1),contextptr));
      }
    }
    if (e.type==_VECT){
      if (is_squarematrix(e)){ 
	// check for numeric entries -> numeric exp
	if (is_fully_numeric(e))
	  return numeric_matrix_exp(e,epsilon(contextptr),contextptr);
	return analytic_apply(at_exp,*e._VECTptr,contextptr);
      }
      return apply(e,contextptr,giac::exp);
    }
    if (is_zero(e))
      return 1;
    if ((e==undef) || (e==plus_inf))
      return e;
    if (e==unsigned_inf)
      return undef;
    if (e==minus_inf)
      return 0;
    if ((e.type==_SYMB) && (e._SYMBptr->sommet==at_ln))
      return e._SYMBptr->feuille;
    int k;
    if (contains(e,_IDNT_pi)){ // if (!approx_mode(contextptr)) 
      gen a,b;
      if (is_linear_wrt(e,_IDNT_pi,a,b,contextptr)){ 
	if (is_multiple_of_12(a*cst_i*gen(trig_deno/2),k) && (k%6!=1) && (k%6!=5) )
	  return (*table_cos[k]+cst_i*(*table_cos[(k+6)%24]))*exp(b,contextptr);
	else {
	  gen kk;
	  kk=normal(a*cst_i,contextptr);
	  if (is_assumed_integer(kk,contextptr)){ 
	    if (is_assumed_integer(normal(rdiv(kk,plus_two),contextptr),contextptr))
	      return exp(b,contextptr);
	    else
	      return pow(minus_one,kk,contextptr)*exp(b,contextptr);
	  }
	  int n,d,q,r;
	  if (is_rational(kk,n,d)){ 
	    q=-n/d;
	    r=-n%d;
	    if (q%2)
	      q=-1;
	    else
	      q=1;
	    if (d<0){ r=-r; d=-d; }
	    if (r<0) r += 2*d;
	    // exp(r*i*pi/d) -> use rootof([1,..,0],cyclotomic(2*d))
	    vecteur vr(r+1);
	    vr[0]=1;
	    vecteur vc(cyclotomic(2*d));
	    return q*symb_rootof(vr,vc,contextptr)*exp(b,contextptr);
	    // initially it was return q*symb_exp(r*(cst_pi*cst_i/d));
	  }
	} // end else multiple of pi/12
      } // end is_linear_wrt
    } // end if contains(e,_IDNT_pi)
    if (is_equal(e))
      return apply_to_equal(e,exp,contextptr);
    return symb_exp(e);
  }
  partial_derivative_onearg D_at_exp(&giac::exp);
  gen taylor_exp (const gen & lim_point,const int ordre,const unary_function_ptr & f,int direction,gen & shift_coeff,GIAC_CONTEXT){
    if (ordre<0)
      return 0; // no symbolic preprocessing
    shift_coeff=0;
    gen image=f(lim_point,contextptr); // should simplify if contains i*pi
    vecteur v(1,image);
    if (is_undef(image))
      return v;
    gen factorielle(1);
    for (int i=1;i<=ordre;++i,factorielle=factorielle*gen(i))
      v.push_back(rdiv(image,factorielle));
    v.push_back(undef);
    return v;
  }
  const string _exp_s("exp");
  string texprintasexp(const gen & g,const string & s,GIAC_CONTEXT){
    return "e^{"+gen2tex(g,contextptr)+"}";
  }
  unary_function_eval __exp(&giac::exp,&D_at_exp,&taylor_exp,_exp_s,0,&texprintasexp);
  unary_function_ptr at_exp (&__exp,0,true);

  symbolic symb_sqrt(const gen & e){
    return symbolic(at_sqrt,e);
  }

  void zint2simpldoublpos(const gen & e,gen & simpl,gen & doubl,bool & pos){
    gen e_copy;
    pos=ck_is_positive(e,context0); // ok
    if (!pos)
      e_copy=-e;
    else
      e_copy=e;
    simpl=1;
    doubl=1;
    vecteur u(pfacprem(e_copy));
    gen f;
    int m,k;
    const_iterateur it=u.begin(),itend=u.end();
    for (;it!=itend;++it){
      f=*it;
      ++it;
      m=it->val;
      if (m%2)
	simpl = simpl*f;
      for (k=0;k<m/2;++k)
	doubl = doubl*f;
    }
  }

  // simplified sqrt without taking care of sign
  gen sqrt_noabs(const gen & e,GIAC_CONTEXT){
    identificateur tmpx(" x");
    vecteur w=solve(tmpx*tmpx-e,tmpx,1,contextptr); 
    if (lidnt(w).empty())
      w=protect_sort(w,contextptr);
    if (w.empty())
      setsizeerr("sqrt_noabs of "+e.print(contextptr));
    return w.back();
  }

  gen sqrt(const gen & e,GIAC_CONTEXT){
    if (e.type==_DOUBLE_){
      if (e._DOUBLE_val>=0){
#ifdef _SOFTMATH_H
	return std::giac_gnuwince_sqrt(e._DOUBLE_val);
#else
	return std::sqrt(e._DOUBLE_val);
#endif
      }
      else
#ifdef _SOFTMATH_H
	return gen(0.0,std::giac_gnuwince_sqrt(-e._DOUBLE_val));
#else
	return gen(0.0,std::sqrt(-e._DOUBLE_val));
#endif
    }
    if (e.type==_REAL)
      return e._REALptr->sqrt();
    if (e.type==_CPLX){
      if (e.subtype){
#ifdef _SOFTMATH_H
	return std::giac_gnuwince_sqrt(gen2complex_d(e));
#else
	return std::sqrt(gen2complex_d(e));
#endif
      }
      // sqrt of an exact complex number
      gen x(re(e,contextptr)),y(im(e,contextptr));
      if (is_zero(y))
	return sqrt(x,contextptr);
      gen rho=sqrt(x*x+y*y,contextptr);
      gen realpart=sqrt(rdiv(x+rho,plus_two),contextptr);
      return realpart+cst_i*rdiv(y,plus_two*realpart);
    }
    if (e.type==_VECT)
      return apply(e,giac::sqrt,contextptr);
    if ((e==undef) || (e==plus_inf) || (e==unsigned_inf))
      return e;
    if (is_perfect_square(e))
      return isqrt(e);
    if (e.type==_INT_ || e.type==_ZINT){ 
      // factorization 
      gen simpl,doubl;
      bool pos;
      zint2simpldoublpos(e,simpl,doubl,pos);
      if (!pos)
	return cst_i*doubl*pow(simpl,plus_one_half,contextptr);
      else
	return doubl*pow(simpl,plus_one_half,contextptr);
    }
    if (e.type==_FRAC)
      return sqrt(e._FRACptr->num*e._FRACptr->den,contextptr)/abs(e._FRACptr->den,contextptr);
    if (e.is_symb_of_sommet(at_inv))
      return inv(sqrt(e._SYMBptr->feuille,contextptr),contextptr);
    return pow(e,plus_one_half,contextptr);
  }
  gen d_sqrt(const gen & e,GIAC_CONTEXT){
    return inv(gen(2)*sqrt(e,contextptr),contextptr);
  }
  partial_derivative_onearg D_at_sqrt(&d_sqrt);
  const string _sqrt_s("sqrt");
  string texprintassqrt(const gen & g,const string & s,GIAC_CONTEXT){
    return "\\sqrt{"+gen2tex(g,contextptr)+"}";
  }
  unary_function_eval __sqrt(&giac::sqrt,&D_at_sqrt,_sqrt_s,0,&texprintassqrt);
  unary_function_ptr at_sqrt (&__sqrt,0,true);

  gen _sq(const gen & e){
    return pow(e,2);
  }
  gen d_sq(const gen & e){
    return gen(2)*e;
  }
  partial_derivative_onearg D_at_sq(&d_sq);
  const string _sq_s("sq");
  string texprintassq(const gen & g,const string & s,GIAC_CONTEXT){
    return gen2tex(g,contextptr)+"^2";
  }
  unary_function_unary __sq(&giac::_sq,&D_at_sq,_sq_s);
  unary_function_ptr at_sq (&__sq,0,true);

  symbolic symb_cos(const gen & e){
    return symbolic(at_cos,e);
  }
  gen cos(const gen & e0,GIAC_CONTEXT){
    gen e=frac_neg_out(e0,contextptr);
    if (e.type==_SPOL1)
      return symb_cos(e);
    if (e.type==_DOUBLE_){
      double d;
      if (angle_radian(contextptr)) 
	d=e._DOUBLE_val;
      else
	d=e._DOUBLE_val*deg2rad_d;
#ifdef _SOFTMATH_H
      return std::giac_gnuwince_cos(d);
#else
      return std::cos(d);
#endif
    }	
    if (e.type==_REAL)
      return e._REALptr->cos(); // FIXME GIAC_CONTEXT
    if (e.type==_CPLX){ 
      if (e.subtype){
	complex<double> d;
	if (angle_radian(contextptr)) 
	  d=gen2complex_d(e);
	else
	  d=gen2complex_d(e)*deg2rad_d;
#ifdef _SOFTMATH_H
	return std::giac_gnuwince_cos(d);
#else
	return std::cos(d);
#endif
      }
      if (e._CPLXptr->type==_REAL){
	gen e1=e;
	if (!angle_radian(contextptr)) 
	  e1=e*deg2rad_g;
	gen e2=im(e1,contextptr);
	e1=re(e1,contextptr);
	bool tmp=angle_radian(contextptr);
	angle_radian(true,contextptr);
	e1= cos(e1,contextptr)*cosh(e2,contextptr)-cst_i*sinh(e2,contextptr)*sin(e1,contextptr);
	angle_radian(tmp,contextptr);
	return e1;
      }
    }
    if (is_squarematrix(e))
      return analytic_apply(at_cos,*e._VECTptr,contextptr);
    if (e.type==_VECT)
      return apply(e,giac::cos,contextptr);
    if (is_zero(e))
      return 1;
    if ( (e.type==_INT_) && (e.val<0) )
      return cos(-e,contextptr);
    if ((e==undef) || is_inf(e))
      return undef;
    int k;
    gen a,b;
    bool doit=false,est_multiple;
    if (angle_radian(contextptr)){
      if (contains(e,_IDNT_pi) && is_linear_wrt(e,_IDNT_pi,a,b,contextptr)){
	est_multiple=is_multiple_of_12(a*gen(trig_deno/2),k) && (k%6!=1) && (k%6!=5);
	doit=true;
      }
    }
    else {
      est_multiple=is_multiple_of_pi_over_12(e,k,angle_radian(contextptr),contextptr) && (k%6!=1) && (k%6!=5);
      doit=est_multiple;
    }
    if (doit){ 
      if (est_multiple){
	gen C=cos(b,contextptr),S=sin(b,contextptr);
	if (k%6==0 || C.type!=_SYMB || S.type!=_SYMB)	
	  return (*table_cos[k])*C+(*table_cos[(k+6)%24])*S;
      }
      else {
	if (is_assumed_integer(a,contextptr)){
	  if (is_assumed_integer(normal(rdiv(a,plus_two),contextptr),contextptr))
	    return cos(b,contextptr);
	  else
	    return pow(minus_one,a,contextptr)*cos(b,contextptr);
	}
	int n,d,q,r;
	if (is_zero(b) && is_rational(a,n,d)){
	  q=n/d;
	  r=n%d;
	  if (r>d/2){
	    r -= d;
	    ++q;
	  }
	  if (q%2)
	    q=-1;
	  else
	    q=1;
	  if (r<0)
	    r=-r;
	  if (!(d%2) && d%4){ 
	    d=d/2; // cos(r/(2*d)*pi) = sin(pi/2(1-r/d))
	    if (angle_radian(contextptr)) 
	      return -q*sin((r-d)/2*cst_pi/d,contextptr);
	    else 
	      return -q*sin(rdiv((r-d)*90,d),contextptr);
	  }
	  if (angle_radian(contextptr)) 
	    return q*symb_cos(r*cst_pi/d);
	  else
	    return q*symb_cos(rdiv(r*180,d));
	}
      }
    }
    if (e.type==_SYMB) {
      unary_function_ptr u=e._SYMBptr->sommet;
      gen f=e._SYMBptr->feuille;
      if (u==at_neg)
	return cos(f,contextptr);
      if (u==at_acos)
	return f;
      if (u==at_asin)
	return sqrt(1-pow(f,2),contextptr);
      if (u==at_atan)
	return sqrt(inv(pow(f,2)+1,contextptr),contextptr);
    }
    if (is_equal(e))
      return apply_to_equal(e,cos,contextptr);
    return symb_cos(e);
  }
  gen d_cos(const gen & e ,GIAC_CONTEXT){
    if (angle_radian(contextptr)) 
      return -(sin(e,contextptr));
    else
      return -deg2rad_e*sin(e,contextptr);
  }
  partial_derivative_onearg D_at_cos(d_cos);
  const string _cos_s("cos");
  unary_function_eval __cos(&giac::cos,&D_at_cos,_cos_s);
  unary_function_ptr at_cos (&__cos,0,true);

  symbolic symb_sin(const gen & e){
    return symbolic(at_sin,e);
  }
  gen sin(const gen & e0,GIAC_CONTEXT){
    gen e=frac_neg_out(e0,contextptr);
    if (e.type==_SPOL1)
      return symb_sin(e);
    if (e.type==_DOUBLE_){
      double d;
      if (angle_radian(contextptr)) 
	d=e._DOUBLE_val;
      else
	d=e._DOUBLE_val*deg2rad_d;
#ifdef _SOFTMATH_H
      return std::giac_gnuwince_sin(d);
#else
      return std::sin(d);
#endif
    }	
    if (e.type==_REAL)
      return e._REALptr->sin();
    if (e.type==_CPLX){ 
      if (e.subtype){
	complex<double> d;
	if (angle_radian(contextptr)) 
	  d=gen2complex_d(e);
	else
	  d=gen2complex_d(e)*deg2rad_d;
#ifdef _SOFTMATH_H
	return std::giac_gnuwince_sin(d);
#else
	return std::sin(d);
#endif
      }
      if (e._CPLXptr->type==_REAL){
	gen e1=e;
	if (!angle_radian(contextptr)) 
	  e1=e*deg2rad_g;
	gen e2=im(e1,contextptr);
	e1=re(e1,contextptr);
	bool tmp=angle_radian(contextptr);
	angle_radian(true,contextptr);
	gen res=sin(e1,contextptr)*cosh(e2,contextptr)+cst_i*sinh(e2,contextptr)*cos(e1,contextptr);
	angle_radian(tmp,contextptr);
	return res;
      }
    }
    if (is_squarematrix(e))
      return analytic_apply(at_sin,*e._VECTptr,contextptr);
    if (e.type==_VECT)
      return apply(e,giac::sin,contextptr);
    if (is_zero(e))
      return e;
    if ( (e.type==_INT_) && (e.val<0) )
      return -sin(-e,contextptr);
    if ((e==undef) || is_inf(e))
      return undef;
    int k;
    gen a,b;
    bool doit=false,est_multiple;
    if (angle_radian(contextptr)){
      if (contains(e,_IDNT_pi) && is_linear_wrt(e,_IDNT_pi,a,b,contextptr)){
	est_multiple=is_multiple_of_12(a*gen(trig_deno/2),k) && (k%6!=1) && (k%6!=5);
	doit=true;
      }
    }
    else {
      est_multiple=is_multiple_of_pi_over_12(e,k,angle_radian(contextptr),contextptr) && (k%6!=1) && (k%6!=5);
      doit=est_multiple;
    }
    if (doit){ 
      if (est_multiple){
	gen C=cos(b,contextptr),S=sin(b,contextptr);
	if (k%6==0 || C.type!=_SYMB || S.type!=_SYMB)
	  return *table_cos[(k+18)%24]*C+(*table_cos[k%24])*S;
      }
      else {
	if (is_assumed_integer(a,contextptr)){
	  if (is_assumed_integer(normal(a/2,contextptr),contextptr))
	    return sin(b,contextptr);
	  else
	    return pow(minus_one,a,contextptr)*sin(b,contextptr);
	}
	int n,d,q,r;
	if (is_zero(b) && is_rational(a,n,d)){
	  q=n/d;
	  r=n%d;
	  if (r>d/2){
	    r -= d;
	    ++q;
	  }
	  if (q%2)
	    q=-1;
	  else
	    q=1;
	  if (r<0){
	    r=-r;
	    q=-q;
	  }
	  if (!(d%2) && d%4){ 
	    d=d/2; // sin(r/(2*d)*pi) = cos(pi/2(1-r/d))
	    if (angle_radian(contextptr))
	      return q*cos((r-d)/2*cst_pi/d,contextptr);
	    else
	      return q*cos(rdiv((r-d)*90,d),contextptr);
	  }
	  if (angle_radian(contextptr)) 
	    return q*symb_sin(r*cst_pi/d);
	  else
	    return q*symb_sin(rdiv(r*180,d));
	}
      }
    }
    if (e.type==_SYMB) {
      unary_function_ptr u=e._SYMBptr->sommet;
      gen f=e._SYMBptr->feuille;
      if (u==at_neg)
	return -sin(f,contextptr);
      if (u==at_asin)
	return f;
      if (u==at_acos)
	return sqrt(1-pow(f,2),contextptr);
      if (u==at_atan)
	return rdiv(f,sqrt(pow(f,2)+1,contextptr));
    }
    if (is_equal(e))
      return apply_to_equal(e,sin,contextptr);
    return symb_sin(e);
  }
  gen d_sin(const gen & g,GIAC_CONTEXT){
    if (angle_radian(contextptr)) 
      return cos(g,contextptr);
    else
      return deg2rad_e*cos(g,contextptr);
  }
  const string _sin_s("sin");
  partial_derivative_onearg D_at_sin(&d_sin);
  unary_function_eval __sin(&giac::sin,&D_at_sin,_sin_s);
  unary_function_ptr at_sin (&__sin,0,true);

  symbolic symb_tan(const gen & e){
    return symbolic(at_tan,e);
  }
  gen tan(const gen & e0,GIAC_CONTEXT){
    gen e=frac_neg_out(e0,contextptr);
    if (e.type==_DOUBLE_){
      double d;
      if (angle_radian(contextptr)) 
	d=e._DOUBLE_val;
      else
	d=e._DOUBLE_val*deg2rad_d;
#ifdef _SOFTMATH_H
      return std::giac_gnuwince_tan(d);
#else
      return std::tan(d);
#endif
    }	
    if (e.type==_REAL)
      return e._REALptr->tan();
    if (e.type==_CPLX){ 
      if (e.subtype){
	complex<double> c(gen2complex_d(e));
	if (!angle_radian(contextptr)) 
	  c *= deg2rad_d;
#ifdef _SOFTMATH_H
	return std::giac_gnuwince_tan(c);
#else
	return std::sin(c)/std::cos(c);
#endif
      }
      if (e._CPLXptr->type==_REAL){
	gen e1=e;
	if (!angle_radian(contextptr)) 
	  e1=e*deg2rad_g;
	gen e2=im(e1,contextptr);
	e1=re(e1,contextptr);
	bool tmp=angle_radian(contextptr);
	angle_radian(true,contextptr);
	e1=tan(e1,contextptr);
	angle_radian(tmp,contextptr);
	e2=cst_i*tanh(e2,contextptr);
	return (e1+e2)/(1-e1*e2);
      }
    }
    if (is_squarematrix(e))
      return analytic_apply(at_tan,*e._VECTptr,contextptr);
    if (e.type==_VECT)
      return apply(e,contextptr,giac::tan);
    if (is_zero(e))
      return e;
    if ((e==undef) || is_inf(e))
      return undef;
    if ( (e.type==_INT_) && (e.val<0) )
      return -tan(-e,contextptr);
    int k;
    if (!approx_mode(contextptr)){ 
      if (is_multiple_of_pi_over_12(e,k,angle_radian(contextptr),contextptr)) 
	return *table_tan[(k%12)];
      else {
	gen kk;
	if (angle_radian(contextptr)) 
	  kk=normal(rdiv(e,cst_pi),contextptr);
	else
	  kk=normal(rdiv(e,180),contextptr);
	if (is_assumed_integer(kk,contextptr))
	  return zero;
	int n,d;
	if (is_rational(kk,n,d)){
	  if (angle_radian(contextptr)) 
	    return symb_tan((n%d)*inv(d,contextptr)*cst_pi);
	  else
	    return symb_tan(rdiv((n%d)*180,d));
	}
      }
    }
    if (e.type==_SYMB) {
      unary_function_ptr u=e._SYMBptr->sommet;
      gen f=e._SYMBptr->feuille;
      if (u==at_neg)
	return -tan(f,contextptr);
      if (u==at_atan)
	return f;
      if (u==at_acos)
	return rdiv(sqrt(1-pow(f,2),contextptr),f);
      if (u==at_asin)
	return rdiv(f,sqrt(1-pow(f,2),contextptr));
    }
    if (is_equal(e))
      return apply_to_equal(e,tan,contextptr);
    return symb_tan(e);
  }
  gen d_tan(const gen & e,GIAC_CONTEXT){
    if (angle_radian(contextptr)) 
      return 1+pow(tan(e,contextptr),2);
    else
      return deg2rad_e*(1+pow(tan(e,contextptr),2));
  }
  partial_derivative_onearg D_at_tan(&d_tan);
  const string _tan_s("tan");
  unary_function_eval __tan(&giac::tan,&D_at_tan,_tan_s);
  unary_function_ptr at_tan (&__tan,0,true);

  symbolic symb_asin(const gen & e){
    return symbolic(at_asin,e);
  }
  gen asinasln(const gen & x,GIAC_CONTEXT){
    return cst_i*ln(sqrt(x*x-1,contextptr)+x,contextptr)+cst_pi_over_2;
  }
  gen asin(const gen & e0,GIAC_CONTEXT){
    static gen normal_sin_pi_12=normal(sin_pi_12,contextptr);
    static gen normal_cos_pi_12=normal(cos_pi_12,contextptr);
    gen e=frac_neg_out(e0,contextptr);
    if (e.type==_DOUBLE_){
      if (e._DOUBLE_val>=-1 && e._DOUBLE_val<=1){
#ifdef _SOFTMATH_H
	double d= std::giac_gnuwince_asin(e._DOUBLE_val);
#else
	double d=std::asin(e._DOUBLE_val);
#endif
	if (angle_radian(contextptr)) 
	  return d;
	else
	  return d*rad2deg_d;
      }
    }
    if (e.type==_REAL)
      return e._REALptr->asin();
    if ( e.type==_DOUBLE_ || (e.type==_CPLX && (e.subtype || e._CPLXptr->type==_REAL)) ){
      if (angle_radian(contextptr)) 
	return no_context_evalf(asinasln(e,contextptr));
      else
	return no_context_evalf(asinasln(e,contextptr))*gen(rad2deg_d);
    }
    if (is_squarematrix(e))
      return analytic_apply(at_asin,*e._VECTptr,contextptr);
    if (e.type==_VECT)
      return apply(e,giac::asin,contextptr);
    if (is_zero(e))
      return e;
    if (is_one(e)){
      if (angle_radian(contextptr))
	return cst_pi_over_2;
      return 90;
    }
    if (e==sin_pi_12 || e==normal_sin_pi_12){
      if (angle_radian(contextptr))
	return rdiv(cst_pi,12);
      return 15;
    }
    if (e==cos_pi_12 || e==normal_cos_pi_12){
      if (angle_radian(contextptr))
	return 5*cst_pi/12;
      return 75;
    }
    if (e==plus_sqrt3_2){
      if (angle_radian(contextptr))
	return rdiv(cst_pi,3);
      return 60;
    }
    if (e==plus_sqrt2_2){
      if (angle_radian(contextptr)) 
	return rdiv(cst_pi,4);
      return 45;
    }
    if (e==plus_one_half){
      if (angle_radian(contextptr)) 
	return rdiv(cst_pi,6);
      return 30;
    }
    if (e==undef)
      return e;
    if ((e.type==_SYMB) && (e._SYMBptr->sommet==at_neg))
      return -asin(e._SYMBptr->feuille,contextptr);
    if ( (e.type==_INT_) && (e.val<0) )
      return -asin(-e,contextptr);
    if (is_equal(e))
      return apply_to_equal(e,asin,contextptr);
    if (is_positive(e*e-1,contextptr))
      return asinasln(e,contextptr);
    return symb_asin(e);
  }
  gen d_asin(const gen & args,GIAC_CONTEXT){
    gen g=inv(recursive_normal(sqrt(1-pow(args,2),contextptr),contextptr),contextptr);
    if (angle_radian(contextptr)) 
      return g;
    else
      return g*deg2rad_e;
  }
  gen taylor_asin (const gen & lim_point,const int ordre,const unary_function_ptr & f, int direction,gen & shift_coeff,GIAC_CONTEXT){
    if (ordre<0)
      return 0; // no symbolic preprocessing
    if (is_one(lim_point)){
      shift_coeff=plus_one_half;
      identificateur x(" "); vecteur v;
      taylor(pow(2+x,minus_one_half,contextptr),x,0,ordre,v,contextptr);
      // integration with shift 
      v=integrate(v,shift_coeff);
      if (!direction)
	direction=1;
      return normal((gen(-direction)*cst_i)*gen(v),contextptr);
    }
    if (is_minus_one(lim_point)){
      shift_coeff=plus_one_half;
      identificateur x(" "); vecteur v;
      taylor(pow(2-x,minus_one_half,contextptr),x,0,ordre,v,contextptr);
      // integration with shift 
      v=integrate(v,shift_coeff);
      return v;
    }
    return taylor(lim_point,ordre,f,direction,shift_coeff,contextptr);
  }
  partial_derivative_onearg D_at_asin(&d_asin);
  const string _asin_s("asin");
  unary_function_eval __asin(&giac::asin,&D_at_asin,&taylor_asin,_asin_s);
  unary_function_ptr at_asin (&__asin,0,true);

  symbolic symb_acos(const gen & e){
    return symbolic(at_acos,e);
  }
  gen acos(const gen & e0,GIAC_CONTEXT){
    gen e=frac_neg_out(e0,contextptr);
    if (e.type==_DOUBLE_){
      if (e._DOUBLE_val>=-1 && e._DOUBLE_val<=1){
#ifdef _SOFTMATH_H
	double d= std::giac_gnuwince_acos(e._DOUBLE_val);
#else
	double d=std::acos(e._DOUBLE_val);
#endif
	if (angle_radian(contextptr)) 
	  return d;
	else
	  return d*rad2deg_d;
      }
    }
    if (e.type==_REAL)
      return e._REALptr->acos();
    if ( e.type==_DOUBLE_ || (e.type==_CPLX && (e.subtype || e._CPLXptr->type==_REAL)) ){
      gen res=-cst_i*no_context_evalf(ln(sqrt(e*e-1,contextptr)+e,contextptr));
      if (angle_radian(contextptr)) 
	return res;
      else
	return res*gen(rad2deg_d);
    }
    if (is_squarematrix(e))
      return analytic_apply(at_acos,*e._VECTptr,contextptr);
    if (e.type==_VECT)
      return apply(e,giac::acos,contextptr);
    if (is_equal(e))
      return apply_to_equal(e,acos,contextptr);
    gen g=asin(e,contextptr);
    if ( (g.type==_SYMB) && (g._SYMBptr->sommet==at_asin) )
      return symb_acos(e);
    if (angle_radian(contextptr)) 
      return normal(cst_pi_over_2-asin(e,contextptr),contextptr);
    else
      return 90-asin(e,contextptr);
  }
  gen d_acos(const gen & args,GIAC_CONTEXT){
    gen g= -inv(recursive_normal(sqrt(1-pow(args,2),contextptr),contextptr),contextptr);
    if (angle_radian(contextptr)) 
      return g;
    else
      return g*deg2rad_e;
  }
  partial_derivative_onearg D_at_acos(&d_acos);
  gen taylor_acos (const gen & lim_point,const int ordre,const unary_function_ptr & f, int direction,gen & shift_coeff,GIAC_CONTEXT){
    if (ordre<0)
      return 0; // no symbolic preprocessing
    if (is_one(lim_point)){
      shift_coeff=plus_one_half;
      identificateur x(" "); vecteur v;
      taylor(pow(2+x,minus_one_half,contextptr),x,0,ordre,v,contextptr);
      // integration with shift 
      v=integrate(v,shift_coeff);
      if (!direction)
	direction=1;
      return -normal((gen(-direction)*cst_i)*gen(v),contextptr);
    }
    if (is_minus_one(lim_point)){
      shift_coeff=plus_one_half;
      identificateur x(" "); vecteur v;
      taylor(pow(2-x,minus_one_half,contextptr),x,0,ordre,v,contextptr);
      // integration with shift 
      v=integrate(v,shift_coeff);
      return -v;
    }
    return taylor(lim_point,ordre,f,direction,shift_coeff,contextptr);
  }
  const string _acos_s("acos");
  unary_function_eval __acos(&giac::acos,&D_at_acos,&taylor_acos,_acos_s);
  unary_function_ptr at_acos (&__acos,0,true);

  symbolic symb_sinh(const gen & e){
    return symbolic(at_sinh,e);
  }
  gen sinh(const gen & e0,GIAC_CONTEXT){
    gen e=frac_neg_out(e0,contextptr);
    if (e.type==_DOUBLE_){
#ifdef _SOFTMATH_H
      return std::giac_gnuwince_sinh(e._DOUBLE_val);
#else
      return std::sinh(e._DOUBLE_val);
#endif
    }
    if (e.type==_REAL)
      return e._REALptr->sinh();
    if (e.type==_CPLX){
      if (e.subtype){
#ifdef _SOFTMATH_H
      return std::giac_gnuwince_sinh(gen2complex_d(e));
#else
      return std::sinh(gen2complex_d(e));
#endif
      }
      if (e._CPLXptr->type==_REAL){
	gen g=exp(e,contextptr);
	return (g-inv(g,contextptr))/2;
      }
    }
    if (is_squarematrix(e))
      return analytic_apply(at_sinh,*e._VECTptr,contextptr);
    if (e.type==_VECT)
      return apply(e,giac::sinh,contextptr);
    if ( is_zero(e) || (e==undef) || (is_inf(e)))
      return e;
    if (is_equal(e))
      return apply_to_equal(e,sinh,contextptr);
    if (e.is_symb_of_sommet(at_neg))
      return -sinh(e._SYMBptr->feuille,contextptr);
    return symb_sinh(e);
  }
  gen d_at_sinh(const gen & e,GIAC_CONTEXT){
    return cosh(e,contextptr);
  }
  partial_derivative_onearg D_at_sinh(&d_at_sinh);
  const string _sinh_s("sinh");
  unary_function_eval __sinh(&giac::sinh,&D_at_sinh,_sinh_s);
  unary_function_ptr at_sinh (&__sinh,0,true);

  symbolic symb_cosh(const gen & e){
    return symbolic(at_cosh,e);
  }
  gen cosh(const gen & e0,GIAC_CONTEXT){
    gen e=frac_neg_out(e0,contextptr);
    if (e.type==_DOUBLE_){
#ifdef _SOFTMATH_H
      return std::giac_gnuwince_cosh(e._DOUBLE_val);
#else
      return std::cosh(e._DOUBLE_val);
#endif
    }
    if (e.type==_REAL)
      return e._REALptr->cosh();
    if (e.type==_CPLX){
      if (e.subtype){
#ifdef _SOFTMATH_H
      return std::giac_gnuwince_cosh(gen2complex_d(e));
#else
      return std::cosh(gen2complex_d(e));
#endif
      }
      if (e._CPLXptr->type==_REAL){
	gen g=exp(e,contextptr);
	return (g+inv(g,contextptr))/2;
      }
    }
    if (is_squarematrix(e))
      return analytic_apply(at_cosh,*e._VECTptr,contextptr);
    if (e.type==_VECT)
      return apply(e,giac::cosh,contextptr);
    if (is_zero(e))
      return 1;
    if (e==undef)
      return e;
    if (is_inf(e))
      return plus_inf;
    if (is_equal(e))
      return apply_to_equal(e,cosh,contextptr);
    if (e.is_symb_of_sommet(at_neg))
      return cosh(e._SYMBptr->feuille,contextptr);
    return symb_cosh(e);
  }
  partial_derivative_onearg D_at_cosh(at_sinh);
  const string _cosh_s("cosh");
  unary_function_eval __cosh(&giac::cosh,&D_at_cosh,_cosh_s);
  unary_function_ptr at_cosh (&__cosh,0,true);

  symbolic symb_tanh(const gen & e){
    return symbolic(at_tanh,e);
  }
  gen tanh(const gen & e0,GIAC_CONTEXT){
    gen e=frac_neg_out(e0,contextptr);
    if (e.type==_DOUBLE_){
#ifdef _SOFTMATH_H
      return std::giac_gnuwince_tanh(e._DOUBLE_val);
#else
      return std::tanh(e._DOUBLE_val);
#endif
    }
    if (e.type==_REAL)
      return e._REALptr->tanh();
    if (e.type==_CPLX){
      if (e.subtype){
	complex<double> c(gen2complex_d(e));
#ifdef _SOFTMATH_H
	return std::giac_gnuwince_tanh(c);
#else
	return std::sinh(c)/std::cosh(c);
#endif
      }
      if (e._CPLXptr->type==_REAL){
	gen g=exp(2*e,contextptr);
	return (g+1)/(g-1);
      }
    }
    if (is_squarematrix(e))
      return analytic_apply(at_tanh,*e._VECTptr,contextptr);
    if (e.type==_VECT)
      return apply(e,giac::tanh,contextptr);
    if (is_zero(e) || (e==undef) || (e==unsigned_inf))
      return e;
    if (e==plus_inf)
      return 1;
    if (e==minus_inf)
      return -1;
    if (is_equal(e))
      return apply_to_equal(e,tanh,contextptr);
    if (e.is_symb_of_sommet(at_neg))
      return -tanh(e._SYMBptr->feuille,contextptr);
    return symbolic(at_tanh,e);
  }
  gen d_tanh(const gen & e,GIAC_CONTEXT){
    return 1-pow(tanh(e,contextptr),2);
  }
  partial_derivative_onearg D_at_tanh(&d_tanh);
  const string _tanh_s("tanh");
  unary_function_eval __tanh(&giac::tanh,&D_at_tanh,_tanh_s);
  unary_function_ptr at_tanh (&__tanh,0,true);

  symbolic symb_asinh(const gen & e){
    return symbolic(at_asinh,e);
  }
  gen asinhasln(const gen & x,GIAC_CONTEXT){
    return ln(x+sqrt(x*x+1,contextptr),contextptr);
  }
  gen asinh(const gen & e0,GIAC_CONTEXT){
    gen e=frac_neg_out(e0,contextptr);
    if (e.type==_DOUBLE_)
      return asinhasln(e,contextptr);
    if (e.type==_REAL)
      return e._REALptr->asinh();
    if ( (e.type==_CPLX) && (e.subtype || e._CPLXptr->type==_REAL))
      return no_context_evalf(asinhasln(e,contextptr));
    if (is_squarematrix(e)){
      context tmp;
      return analytic_apply(at_asinh,*e._VECTptr,&tmp); 
    }
    if (e.type==_VECT)
      return apply(e,giac::asinh,contextptr);
    if (is_zero(e) || is_inf(e))
      return e;
    if (e==undef)
      return e;
    if (is_equal(e))
      return apply_to_equal(e,asinh,contextptr);
    return ln(e+sqrt(pow(e,2)+1,contextptr),contextptr);
    // return symbolic(at_asinh,e);
  }
  gen d_asinh(const gen & args,GIAC_CONTEXT){
    return inv(recursive_normal(sqrt(pow(args,2)+1,contextptr),contextptr),contextptr);
  }
  partial_derivative_onearg D_at_asinh(&d_asinh);
  const string _asinh_s("asinh");
  unary_function_eval __asinh(&giac::asinh,&D_at_asinh,_asinh_s);
  unary_function_ptr at_asinh (&__asinh,0,true);

  symbolic symb_acosh(const gen & e){
    return symbolic(at_cosh,e);
  }
  gen acoshasln(const gen & x,GIAC_CONTEXT){
    return ln(x+sqrt(x+1,contextptr)*sqrt(x-1,contextptr),contextptr);
  }
  gen acosh(const gen & e0,GIAC_CONTEXT){
    gen e=frac_neg_out(e0,contextptr);
    if (e.type==_DOUBLE_)
      return acoshasln(e,contextptr);
    if (e.type==_REAL)
      return e._REALptr->acosh();
    if ( (e.type==_CPLX) && (e.subtype|| e._CPLXptr->type==_REAL))
      return no_context_evalf(acoshasln(e,contextptr));
    if (is_squarematrix(e))
      return analytic_apply(at_acosh,*e._VECTptr,0);
    if (e.type==_VECT)
      return apply(e,giac::acosh,contextptr);
    if (is_one(e))
      return 0;
    if (e==plus_inf)
      return plus_inf;
    if (e==undef)
      return e;
    if (is_equal(e))
      return apply_to_equal(e,acosh,contextptr);
    return ln(e+sqrt(pow(e,2)-1,contextptr),contextptr);
    // return symbolic(at_acosh,e);
  }
  gen d_acosh(const gen & args,GIAC_CONTEXT){
    return inv(recursive_normal(sqrt(pow(args,2)-1,contextptr),contextptr),contextptr);
  }
  partial_derivative_onearg D_at_acosh(&d_acosh);
  const string _acosh_s("acosh");
  unary_function_eval __acosh(&giac::acosh,&D_at_acosh,_acosh_s);
  unary_function_ptr at_acosh (&__acosh,0,true);

  symbolic symb_atanh(const gen & e){
    return symbolic(at_atanh,e);
  }
  gen atanh(const gen & e0,GIAC_CONTEXT){
    gen e=frac_neg_out(e0,contextptr);
    if (e.type==_DOUBLE_){
#ifdef _SOFTMATH_H
      return std::giac_gnuwince_log((1+e._DOUBLE_val)/(1-e._DOUBLE_val))/2;
#else
      return std::log((1+e._DOUBLE_val)/(1-e._DOUBLE_val))/2;
#endif
    }
    if (e.type==_REAL)
      return e._REALptr->atanh();
    if ( (e.type==_CPLX) && (e.subtype || e._CPLXptr->type==_REAL))
      return no_context_evalf(rdiv(ln(rdiv(1+e,1-e),contextptr),plus_two));
    if (is_squarematrix(e))
      return analytic_apply(at_atanh,*e._VECTptr,0);
    if (e.type==_VECT)
      return apply(e,giac::atanh,contextptr);
    if (is_zero(e))
      return e;
    if (is_one(e))
      return plus_inf;
    if (is_minus_one(e))
      return minus_inf;
    if (e==undef)
      return e;
    if (is_equal(e))
      return apply_to_equal(e,atanh,contextptr);
    return rdiv(ln(rdiv(1+e,1-e),contextptr),plus_two);
    // return symbolic(at_atanh,e);
  }
  gen d_atanh(const gen & args,GIAC_CONTEXT){
    return inv(1-pow(args,2),contextptr);
  }
  partial_derivative_onearg D_at_atanh(&d_atanh);
  const string _atanh_s("atanh");
  unary_function_eval __atanh(&giac::atanh,&D_at_atanh,_atanh_s);
  unary_function_ptr at_atanh (&__atanh,0,true);

  symbolic symb_quote(const gen & arg){
    return symbolic(at_quote,arg);
  }
  gen quote(const gen & args,GIAC_CONTEXT){
    if (args.type==_VECT && args.subtype==_SEQ__VECT && !args._VECTptr->empty() && args._VECTptr->front().type==_FUNC){
      const unary_function_ptr & u =*args._VECTptr->front()._FUNCptr;
      vecteur v=vecteur(args._VECTptr->begin()+1,args._VECTptr->end());
      gen arg=eval(gen(v,_SEQ__VECT),eval_level(contextptr),contextptr);
      return symbolic(u,arg);
    }
    return args;
  }
  partial_derivative_onearg D_at_quote(&quote);
  const string _quote_s("quote");
  unary_function_eval __quote(&giac::quote,&D_at_quote,_quote_s);
  unary_function_ptr at_quote (&__quote,_QUOTE_ARGUMENTS,true);

  symbolic symb_unquote(const gen & arg){
    return symbolic(at_unquote,arg);
  }
  gen unquote(const gen & arg,GIAC_CONTEXT){
    return eval(arg,1,contextptr);
  }
  partial_derivative_onearg D_at_unquote(&unquote);
  const string _unquote_s("unquote");
  unary_function_eval __unquote(&giac::unquote,&D_at_unquote,_unquote_s);
  unary_function_ptr at_unquote (&__unquote,0,true);

  symbolic symb_order_size(const gen & e){
    return symbolic(at_order_size,e);
  }
  gen order_size(const gen & arg){
    return symb_order_size(arg);
  }
  partial_derivative_onearg D_at_order_size(&order_size);
  const string _order_size_s("order_size");
  unary_function_unary __order_size(&giac::order_size,&D_at_order_size,_order_size_s);
  unary_function_ptr at_order_size (&__order_size,0,true);

  gen re(const gen & a,GIAC_CONTEXT){
    if (is_equal(a))
      return apply_to_equal(a,re,contextptr);
    return a.re(contextptr);
  }
  const string _re_s("re");
  string texprintasre(const gen & g,const string & s,GIAC_CONTEXT){
    return "\\Re("+gen2tex(g,contextptr)+")";
  }
  unary_function_eval __re(&giac::re,_re_s,0,&texprintasre);
  unary_function_ptr at_re (&__re,0,true);

  gen im(const gen & a,GIAC_CONTEXT){
    if (is_equal(a))
      return apply_to_equal(a,im,contextptr);
    return a.im(contextptr);
  }
  const string _im_s("im");
  string texprintasim(const gen & g,const string & s,GIAC_CONTEXT){
    return "\\Im("+gen2tex(g,contextptr)+")";
  }
  unary_function_eval __im(&giac::im,_im_s,0,&texprintasim);
  unary_function_ptr at_im (&__im,0,true);

  symbolic symb_conj(const gen & e){
    return symbolic(at_conj,e);
  }
  gen conj(const gen & a,GIAC_CONTEXT){
    if (is_equal(a))
      return apply_to_equal(a,conj,contextptr);
    return a.conj(contextptr);
  }
  const string _conj_s("conj");
  string texprintasconj(const gen & g,const string & s,GIAC_CONTEXT){
    return "\\overline{"+gen2tex(g,contextptr)+"}";
  }
  unary_function_eval __conj(&giac::conj,_conj_s,0,&texprintasconj);
  unary_function_ptr at_conj (&__conj,0,true);

  gen taylor_sign (const gen & lim_point,const int ordre,const unary_function_ptr & f, int direction,gen & shift_coeff,GIAC_CONTEXT){
    if (ordre<0)
      return 0; // no symbolic preprocessing
    shift_coeff=0;
    if (is_strictly_positive(lim_point,contextptr) || (is_zero(lim_point) && direction==1))
      return makevecteur(1);
    if (is_strictly_positive(-lim_point,contextptr) || (is_zero(lim_point) && direction==-1))
      return makevecteur(-1);
    setsizeerr("Taylor sign with unsigned limit");
    return 0;
  }

  gen _sign(const gen & g,GIAC_CONTEXT){
    return apply(g,contextptr,giac::sign);
  }
  symbolic symb_sign(const gen & e){
    return symbolic(at_sign,e);
  }
  const string _sign_s("sign");
  partial_derivative_onearg D_at_sign(at_zero);
  unary_function_eval __sign(_sign,&D_at_sign,&taylor_sign,_sign_s);
  unary_function_ptr at_sign (&__sign,0,true);

  gen _abs(const gen & args,GIAC_CONTEXT){
    return apply(args,contextptr,giac::abs);
  }
  symbolic symb_abs(const gen & e){
    return symbolic(at_abs,e);
  }
  gen taylor_abs (const gen & lim_point,const int ordre,const unary_function_ptr & f, int direction,gen & shift_coeff,GIAC_CONTEXT){
    if (ordre<0)
      return 0; // no symbolic preprocessing
    shift_coeff=0;
    if (is_strictly_positive(lim_point,contextptr) || (is_zero(lim_point) && direction==1))
      return makevecteur(lim_point,1);
    if (is_strictly_positive(-lim_point,contextptr) || (is_zero(lim_point) && direction==-1))
      return makevecteur(-lim_point,-1);
    setsizeerr("Taylor abs with unsigned limit");
    return 0;
  }
  const string _abs_s("abs");
  partial_derivative_onearg D_at_abs(at_sign);
  unary_function_eval __abs(&_abs,&D_at_abs,&taylor_abs,_abs_s);
  unary_function_ptr at_abs (&__abs,0,true);

  symbolic symb_arg(const gen & e){
    return symbolic(at_arg,e);
  }
  const string _arg_s("arg");
  unary_function_eval __arg(&giac::arg,_arg_s);
  unary_function_ptr at_arg (&__arg,0,true);

  symbolic symb_cyclotomic(const gen & e){
    return symbolic(at_cyclotomic,e);
  }
  gen _cyclotomic(const gen & a){
    if (a.type!=_INT_)
      return symb_cyclotomic(a);
    return cyclotomic(a.val);
  }
  const string _cyclotomic_s("cyclotomic");
  unary_function_unary __cyclotomic(&giac::_cyclotomic,_cyclotomic_s);
  unary_function_ptr at_cyclotomic (&__cyclotomic,0,true);

  string printassto(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if ( (feuille.type!=_VECT) || (feuille._VECTptr->size()!=2) )
      return sommetstr+'('+feuille.print(contextptr)+')';
    vecteur & v=*feuille._VECTptr;
    if (xcas_mode(contextptr)==3){
      if ( (v.front().type==_SYMB) && (v.front()._SYMBptr->sommet==at_program)){
	gen & b=v.front()._SYMBptr->feuille._VECTptr->back();
	if  (b.type==_VECT || (b.type==_SYMB && (b._SYMBptr->sommet==at_local || b._SYMBptr->sommet==at_bloc))){
	  string s(v.front().print(contextptr));
	  s=s.substr(10,s.size()-10);
	  return ":"+v.back().print(contextptr)+s;
	}
	else {
	  vecteur & tmpv = *v.front()._SYMBptr->feuille._VECTptr;
	  if (tmpv[0].type==_VECT && tmpv[0].subtype==_SEQ__VECT && tmpv[0]._VECTptr->size()==1)
	    return tmpv[2].print(contextptr)+" => "+v.back().print(contextptr)+"("+tmpv[0]._VECTptr->front().print(contextptr)+")";
	  else
	    return tmpv[2].print(contextptr)+" => "+v.back().print(contextptr)+"("+tmpv[0].print(contextptr)+")";
	}
      }
      else 
	return v.front().print(contextptr)+" => "+v.back().print(contextptr);
    }
    string s(v.back().print(contextptr)+":=");
    if (v.front().type==_SEQ__VECT)
      return s+"("+v.front().print(contextptr)+")";
    else
      return s+v.front().print(contextptr);
  }
  string texprintassto(const gen & g,const string & sommetstr,GIAC_CONTEXT){
    if ( (g.type!=_VECT) || (g._VECTptr->size()!=2) )
      return sommetstr+'('+gen2tex(g,contextptr)+')';
    string s(gen2tex(g._VECTptr->back(),contextptr)+":=");
    if (g._VECTptr->front().type==_SEQ__VECT)
      return s+"("+gen2tex(g._VECTptr->front(),contextptr)+")";
    else
      return s+gen2tex(g._VECTptr->front(),contextptr);
  }

  // store a in b
  bool signal_store=true;
  gen sto(const gen & a,const gen & b,const context * contextptr){
    return sto(a,b,false,contextptr);
  }
  // in_place==true to store in vector/matrices without making a new copy
  gen sto(const gen & a,const gen & b,bool in_place,const context * contextptr_){
    // *logptr(contextptr) << "Sto " << "->" << b << endl;
    const context * contextptr=contextptr_;
    if (contextptr && contextptr->parent)
      contextptr=contextptr->parent;
    if (b.type==_SYMB){ 
      if (b._SYMBptr->sommet==at_hash && b._SYMBptr->feuille.type==_STRNG)
	return sto(a,gen(*b._SYMBptr->feuille._STRNGptr,contextptr),in_place,contextptr);
      if (b._SYMBptr->sommet==at_double_deux_points){ 
	// TI path
	gen a1,bb;
	check_binary(b._SYMBptr->feuille,a1,bb);
	gen ab=a1.eval(eval_level(contextptr),contextptr);
	if (ab.type==_VECT){
	  vecteur v=*ab._VECTptr;
	  iterateur it=v.begin(),itend=v.end();
	  for (;it!=itend;++it){
	    if (it->type!=_VECT || it->_VECTptr->size()!=2)
	      continue;
	    vecteur & w=*it->_VECTptr;
	    if (w[0]==bb)
	      w[1]=a;
	  }
	  if (it==itend)
	    v.push_back(makevecteur(bb,a));
	  return sto(gen(v,_FOLDER__VECT),a1,in_place,contextptr);
	}
	if (a1.type==_IDNT)
	  return sto(gen(vecteur(1,makevecteur(bb,a)),_FOLDER__VECT),a1,in_place,contextptr);
      } // end TI path
    }
    if (b.type==_IDNT){
      if (!contextptr){
	// Remove stale local assignements
	try {
	  b._IDNTptr->eval(1,b,contextptr); 
	} catch (std::runtime_error & e) { }
      }
      if (*b._IDNTptr==_IDNT_pi || *b._IDNTptr==_IDNT_infinity())
	setsizeerr(b.print(contextptr)+": reserved word");
      gen aa(a),ans(aa);
      if (a==b)
	return _purge(b,contextptr);
      if ( (a.type==_SYMB) && (a._SYMBptr->sommet==at_parameter)){
	gen inter=a._SYMBptr->feuille,debut,fin,saut;
	bool calc_aa=false;
	if (inter.type==_VECT){
	  vecteur & interv=*inter._VECTptr;
	  int inters=interv.size();
	  if (inters>=3){
	    debut=interv[0];
	    fin=interv[1];
	    aa=interv[2];
	    if (inters>=4)
	      saut=interv[3];
	  }
	  if (inters==2){
	    aa=interv.back();
	    inter=interv.front();
	  }
	}
	else
	  calc_aa=true;
	if ( (inter.type==_SYMB) && (inter._SYMBptr->sommet==at_interval) ){
	  debut=inter._SYMBptr->feuille._VECTptr->front();
	  fin=inter._SYMBptr->feuille._VECTptr->back();
	}
	if (calc_aa)
	  aa=rdiv(debut+fin,plus_two);
	if (is_zero(saut))
	  saut=(fin-debut)/100.;
	ans=symbolic(at_parameter,makevecteur(b,debut,fin,aa,saut));
      } // end parameter
      if (contextptr){
	const context * ptr=contextptr;
	bool done=false;
	for (;ptr->previous && ptr->tabptr;ptr=ptr->previous){
	  sym_tab::iterator it=ptr->tabptr->find(*b._IDNTptr->name),itend=ptr->tabptr->end();
	  if (it!=itend){ // found in current local context
	    // check that the current value is a thread pointer
	    if (it->second.type==_POINTER_ && it->second.subtype==_THREAD_POINTER){
	      if (it->second._POINTER_val!=(void *)contextptr_)
		settypeerr(b.print(contextptr)+" is locked by thread "+it->second.print(contextptr));
	    }
	    it->second=aa;
	    done=true;
	    break;
	  }
	}
	if (!done) {// store b globally
	  if (contains(lidnt(a),b))
	    *logptr(contextptr) << b.print(contextptr)+": recursive definition" << endl;
	  sym_tab * symtabptr=contextptr->globalcontextptr?contextptr->globalcontextptr->tabptr:contextptr->tabptr;
	  sym_tab::iterator it=symtabptr->find(*b._IDNTptr->name),itend=symtabptr->end();
	  if (it!=itend){ 
	    // check that the current value is a thread pointer
	    if (it->second.type==_POINTER_ && it->second.subtype==_THREAD_POINTER){
	      if (it->second._POINTER_val!=(void *)contextptr_)
		settypeerr(b.print(contextptr)+" is locked by thread "+it->second.print(contextptr));
	    }
	    it->second=aa;
	  }
	  else
	    (*symtabptr)[*b._IDNTptr->name]=aa;
	}
	if (!child_id && signal_store)
	  _signal(symb_quote(symbolic(at_sto,makevecteur(aa,b))),contextptr);
	return ans;
      } // end if (contextptr)
      if (contains(lidnt(a),b))
	*logptr(contextptr) << b.print(contextptr)+": recursive definition" << endl;
      if (!b._IDNTptr->localvalue->empty() && (b.subtype!=_GLOBAL__EVAL))
	b._IDNTptr->localvalue->back()=aa;
      else {
	if (current_folder_name.type==_IDNT && current_folder_name._IDNTptr->value && current_folder_name._IDNTptr->value->type==_VECT){
	  vecteur v=*current_folder_name._IDNTptr->value->_VECTptr;
	  iterateur it=v.begin(),itend=v.end();
	  for (;it!=itend;++it){
	    if (it->type!=_VECT || it->_VECTptr->size()!=2)
	      continue;
	    vecteur & w=*it->_VECTptr;
	    if (w[0]==b){
	      w[1]=aa;
	      break;
	    }
	  }
	  if (it==itend)
	    v.push_back(makevecteur(b,aa));
	  gen gv(v,_FOLDER__VECT);
	  *current_folder_name._IDNTptr->value=gv;
	  if (!child_id && signal_store)
	    _signal(symb_quote(symbolic(at_sto,makevecteur(gv,current_folder_name))),contextptr);
	  return ans;
	}
	else {
	  if (b._IDNTptr->value)
	    delete b._IDNTptr->value;
	  b._IDNTptr->value = new gen(aa);
	  if (!child_id && signal_store)
	    _signal(symb_quote(symbolic(at_sto,makevecteur(aa,b))),contextptr);
	  if (!secure_run && variables_are_files(contextptr)){
	    ofstream a_out((*b._IDNTptr->name+cas_suffixe).c_str());
	    a_out << aa << endl;
	  }
	}
      }
      return ans;
    } // end b.type==_IDNT
    if (b.type==_VECT){
      if (a.type!=_VECT)
	settypeerr();
      return apply(a,b,contextptr,giac::sto);
    }
    if ( (b.type==_SYMB) && (b._SYMBptr->sommet==at_at) ){
      // Store a in a vector or array or map
      gen destination=b._SYMBptr->feuille._VECTptr->front(); // variable name
      gen valeur;
      if (!contextptr && in_place && destination.type==_IDNT && !destination._IDNTptr->localvalue->empty() && local_eval(contextptr) )
	valeur=do_local_eval(*destination._IDNTptr,eval_level(contextptr),false);
      else
	valeur=destination.eval(in_place?1:eval_level(contextptr),contextptr);
      if ( valeur.type==_INT_ && valeur.val==0 && destination.type==_IDNT && !destination._IDNTptr->localvalue->empty() )
	valeur=destination; // non (0) initialized local var
      gen indice=b._SYMBptr->feuille._VECTptr->back().eval(eval_level(contextptr),contextptr);
      if ( (destination.type!=_IDNT) || (valeur.type!=_VECT && valeur.type!=_MAP && valeur.type!=_IDNT  && valeur.type!=_STRNG) )
	settypeerr("sto "+b.print(contextptr)+ "="+valeur.print(contextptr)+" not allowed!");
      if (valeur.type==_IDNT){ 
	// no previous vector at destination, 
	// create one in TI mode or create a map
	gen g;
	if (xcas_mode(contextptr)==3 && indice.type==_INT_ && indice.val>=0 ){
	  vecteur v(indice.val+1,zero);
	  v[indice.val]=a;
	  g=gen(v,destination.subtype);
	}
	else {
	  g=makemap();
	  (*g._MAPptr)[indice]=a;
	}
	return sto(g,destination,in_place,contextptr);
      }
      if (valeur.type==_STRNG){
	if (indice.type!=_INT_ || a.type!=_STRNG || a._STRNGptr->empty())
	  setsizeerr();
	if (indice.val<0 || indice.val>=valeur._STRNGptr->size())
	  setdimerr();
	if (in_place){
	  (*valeur._STRNGptr)[indice.val]=(*a._STRNGptr)[0];
	  return 1;
	}
	else {
	  string m(*valeur._STRNGptr);
	  m[indice.val]=(*a._STRNGptr)[0];
	  return sto(string2gen(m,false),destination,in_place,contextptr);
	}
      }
      if (valeur.type==_MAP){
	if (valeur.subtype==1){ // array
	  gen_map::iterator it=valeur._MAPptr->find(indice),itend=valeur._MAPptr->end();
	  if (it==itend)
	    setdimerr("Index outside of range");
	  if (xcas_mode(contextptr)==1)
	    in_place=true;
	}
	if (in_place){
	  (*valeur._MAPptr)[indice]=a;
	  return 1;
	}
	else {
	  gen_map m(*valeur._MAPptr);
	  m[indice]=a;
	  return sto(m,destination,in_place,contextptr);
	}
      }
      vecteur * vptr=0;
      vecteur v;
      if (in_place)
	vptr=valeur._VECTptr;
      else
	v=*valeur._VECTptr;
      if (indice.type!=_VECT){
	if (indice.type!=_INT_ || indice.val<0 )
	  settypeerr("Bad index "+indice.print(contextptr));
	// check size
	int is=in_place?vptr->size():v.size();
	for (;is<=indice.val;++is){
	  if (in_place)
	    vptr->push_back(zero);
	  else
	    v.push_back(zero);
	}
	// change indice's value
	if (in_place){
	  (*vptr)[indice.val]=a;
	  return 1;
	}
	else {
	  v[indice.val]=a;
	  return sto(gen(v,valeur.subtype),destination,in_place,contextptr);
	}
      }
      // here indice is of type _VECT, we store inside a matrix
      vecteur empile;
      iterateur it=indice._VECTptr->begin(),itend=indice._VECTptr->end();
      for (;;){
	if (it->type!=_INT_)
	  settypeerr("Bad index "+indice.print(contextptr));
	if (!in_place)
	  empile.push_back(v);
	gen tmp=in_place?(*vptr)[it->val]:v[it->val];
	++it;
	if (it==itend)
	  break;
	if (tmp.type!=_VECT)
	  settypeerr("Bad index "+indice.print(contextptr));
	if (in_place)
	  vptr= tmp._VECTptr;
	else
	  v=*tmp._VECTptr;
      }
      --itend;
      if (in_place){
	(*vptr)[itend->val]=a;
	return 1;
      }
      v[itend->val]=a;
      vecteur oldv;
      it=indice._VECTptr->begin();
      for (;;){
	if (itend==it)
	  break;
	--itend;
	empile.pop_back();
	oldv=*(empile.back()._VECTptr);
	oldv[itend->val]=v;
	v=oldv;
      }
      return sto(v,destination,in_place,contextptr);
    }
    if (b.type==_FUNC){
      string errmsg=b.print(contextptr)+ " is a reserved word, sto not allowed:";
      *logptr(contextptr) << errmsg << endl;
      return makevecteur(string2gen(errmsg,false),a);
    }
    settypeerr("sto "+b.print(contextptr)+ " not allowed!");
    return 0;
  }
  symbolic symb_sto(const gen & a,gen & b,bool in_place){
    if (in_place)
      return symbolic(at_array_sto,makevecteur(a,b));
    return symbolic(at_sto,makevecteur(a,b));
  }
  symbolic symb_sto(const gen & e){
    return symbolic(at_sto,e);
  }
  gen _sto(const gen & a,const context * contextptr){
    if (a.type!=_VECT)
      return symb_sto(a);
    if (rpn_mode){
      if (a._VECTptr->size()<2)
	toofewargs("STO");
      gen c=a._VECTptr->back();
      a._VECTptr->pop_back();
      gen b=a._VECTptr->back();
      a._VECTptr->pop_back();
      sto(b,c,contextptr);
      return gen(*a._VECTptr,_RPN_STACK__VECT);
    }
    if (a._VECTptr->size()!=2)
      settypeerr();
    return sto(a._VECTptr->front(),a._VECTptr->back(),contextptr);
  }
  const string _sto_s("sto");
  unary_function_eval __sto(&giac::_sto,_sto_s,&printassto,&texprintassto);
  unary_function_ptr at_sto (&__sto,0,true); 
  // NB argument quoting for sto is done in eval in symbolic.cc

  gen _array_sto(const gen & a,const context * contextptr){
    if (a.type!=_VECT ||a._VECTptr->size()!=2)
      settypeerr();
    gen value=a._VECTptr->front().eval(eval_level(contextptr),contextptr);
    return sto(value,a._VECTptr->back(),true,contextptr);
  }
  const string _array_sto_s("array_sto");
  unary_function_eval __array_sto(&giac::_array_sto,_array_sto_s);
  unary_function_ptr at_array_sto (&__array_sto,_QUOTE_ARGUMENTS,true);

  string printasincdec(const gen & feuille,char ch,bool tex,GIAC_CONTEXT){
    if (feuille.type!=_VECT){
      string s(tex?gen2tex(feuille,contextptr):feuille.print(contextptr));
      return xcas_mode(contextptr)?((s+":="+s+ch)+"1"):((s+ch)+ch);
    }
    vecteur & v = *feuille._VECTptr;
    if (v.size()!=2)
      setdimerr();
    gen & a=v.front();
    gen & b=v.back();
    string sa((tex?gen2tex(a,contextptr):a.print(contextptr)));
    string sb((tex?gen2tex(b,contextptr):b.print(contextptr)));
    return xcas_mode(contextptr)?sa+":="+sa+ch+sb:(sa+ch+'='+sb);
  }

  string printasincrement(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    return printasincdec(feuille,'+',false,contextptr);
  }

  string printasdecrement(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    return printasincdec(feuille,'-',false,contextptr);
  }

  string texprintasincrement(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    return printasincdec(feuille,'+',true,contextptr);
  }

  string texprintasdecrement(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    return printasincdec(feuille,'-',true,contextptr);
  }

  string printasmultcrement(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    return printasincdec(feuille,'*',false,contextptr);
  }

  string printasdivcrement(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    return printasincdec(feuille,'/',false,contextptr);
  }

  string texprintasmultcrement(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    return printasincdec(feuille,'*',true,contextptr);
  }

  string texprintasdivcrement(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    return printasincdec(feuille,'/',true,contextptr);
  }

  gen increment(const gen & var,const gen & val_orig,bool negatif,bool mult,GIAC_CONTEXT){
    if (var.type!=_IDNT)
      settypeerr("Increment");
    gen val=val_orig.eval(1,contextptr);
    if (negatif)
      val=mult?inv(val,contextptr):-val;
    if (contextptr){
      sym_tab::iterator it,itend;
      for(;contextptr;) {
	it=contextptr->tabptr->find(*var._IDNTptr->name);
	itend=contextptr->tabptr->end();
	if (it!=itend)
	  break;
	if (!contextptr->previous){
	  it=contextptr->globalcontextptr->tabptr->find(*var._IDNTptr->name);
	  if (it!=itend)
	    break;
	}
	contextptr=contextptr->previous;
      }
      if (!contextptr)
	setsizeerr("Increment");
      return it->second=mult?(it->second*val):(it->second+val);
    }
    vecteur * w=var._IDNTptr->localvalue;
    if (!w->empty() && var.subtype!=_GLOBAL__EVAL)
      return w->back()=w->back()+val;
    if (!var._IDNTptr->value)
      setsizeerr("Non assigned variable");
    return *var._IDNTptr->value=mult?(*var._IDNTptr->value*val):(*var._IDNTptr->value+val);
  }
  gen _increment(const gen & a,const context * contextptr){
    if (a.type!=_VECT)
      return increment(a,1,false,false,contextptr);
    if (a.type!=_VECT || a._VECTptr->size()!=2)
      settypeerr();
    return increment(a._VECTptr->front(),a._VECTptr->back(),false,false,contextptr);
  }
  const string _increment_s("increment");
  unary_function_eval __increment(&giac::_increment,_increment_s,&printasincrement,&texprintasincrement);
  unary_function_ptr at_increment (&__increment,_QUOTE_ARGUMENTS,true); 

  gen _decrement(const gen & a,const context * contextptr){
    if (a.type!=_VECT)
      return increment(a,1,true,false,contextptr);
    if (a._VECTptr->size()!=2)
      settypeerr();
    return increment(a._VECTptr->front(),a._VECTptr->back(),true,false,contextptr);
  }
  const string _decrement_s("decrement");
  unary_function_eval __decrement(&giac::_decrement,_decrement_s,&printasdecrement,&texprintasdecrement);
  unary_function_ptr at_decrement (&__decrement,_QUOTE_ARGUMENTS,true); 

  gen _multcrement(const gen & a,const context * contextptr){
    if (a.type!=_VECT)
      return increment(a,1,false,true,contextptr);
    if (a.type!=_VECT || a._VECTptr->size()!=2)
      settypeerr();
    return increment(a._VECTptr->front(),a._VECTptr->back(),false,true,contextptr);
  }
  const string _multcrement_s("multcrement");
  unary_function_eval __multcrement(&giac::_multcrement,_multcrement_s,&printasmultcrement,&texprintasmultcrement);
  unary_function_ptr at_multcrement (&__multcrement,_QUOTE_ARGUMENTS,true); 

  gen _divcrement(const gen & a,const context * contextptr){
    if (a.type!=_VECT)
      return increment(a,1,true,true,contextptr);
    if (a._VECTptr->size()!=2)
      settypeerr();
    return increment(a._VECTptr->front(),a._VECTptr->back(),true,true,contextptr);
  }
  const string _divcrement_s("divcrement");
  unary_function_eval __divcrement(&giac::_divcrement,_divcrement_s,&printasdivcrement,&texprintasdivcrement);
  unary_function_ptr at_divcrement (&__divcrement,_QUOTE_ARGUMENTS,true); 

  bool is_assumed_integer(const gen & g,GIAC_CONTEXT){
    if (is_integer(g))
      return true;
    if (g.type==_IDNT) {// FIXME GIAC_CONTEXT
      gen tmp=g._IDNTptr->eval(1,g,contextptr);
      if (tmp.type==_VECT && tmp.subtype==_ASSUME__VECT){
	vecteur & v = *tmp._VECTptr;
	if (!v.empty() && (v.front()==_INT_ || v.front()==_ZINT) )
	  return true;
      }
      return is_integer(tmp);
    }
    if (g.type!=_SYMB)
      return false;
    unary_function_ptr & u=g._SYMBptr->sommet;
    gen & f=g._SYMBptr->feuille;
    if ( (u==at_neg) || (u==at_abs) )
      return is_assumed_integer(f,contextptr);
    if ( (u==at_plus) || (u==at_prod) ){
      if (f.type!=_VECT)
	return is_assumed_integer(f,contextptr);
      const_iterateur it=f._VECTptr->begin(),itend=f._VECTptr->end();
      for (;it!=itend;++it)
	if (!is_assumed_integer(*it,contextptr))
	  return false;
      return true;
    }
    return false;
  }
  // v = previous assumptions, a=the real value, direction
  // is positive for [a,+inf[, negative for ]-inf,a]
  // |direction| = 1 (large) or 2 (strict) 
  // v = previous assumptions, a=the real value, direction
  // is positive for [a,+inf[, negative for ]-inf,a]
  // |direction| = 1 (large) or 2 (strict) 
  gen doubleassume_and(const vecteur & v,const gen & a,int direction,bool or_assumption,GIAC_CONTEXT){
    vecteur v_intervalle,v_excluded;
    if ( (v.size()>=3) && (v[1].type==_VECT) && (v[2].type==_VECT) ){
      v_intervalle=*v[1]._VECTptr;
      v_excluded=*v.back()._VECTptr;
    }
    gen v0=_DOUBLE_;
    v0.subtype=1;
    if (!v.empty())
      v0=v.front();
    if (!(direction %2) && !equalposcomp(v_excluded,a))
      v_excluded.push_back(a);
    if (or_assumption){ 
      // remove excluded values if they are in the interval we add
      vecteur old_v(v_excluded);
      v_excluded.clear();
      const_iterateur it=old_v.begin(),itend=old_v.end();
      for (;it!=itend;++it){
	bool a_greater_sup=ck_is_greater(a,*it,contextptr);
	if (a_greater_sup && (direction<0) )
	  continue;
	if (!a_greater_sup && (direction>0) )
	  continue;
	v_excluded.push_back(*it);
      }
    }
    if (v_intervalle.empty() || or_assumption){
      if (direction>0)
	v_intervalle.push_back(gen(makevecteur(a,plus_inf),_LINE__VECT));
      else
	v_intervalle.push_back(gen(makevecteur(minus_inf,a),_LINE__VECT));
      if (or_assumption)
	return gen(makevecteur(v0,v_intervalle,v_excluded),_ASSUME__VECT);
    }
    else { // intersection of [a.+inf[ with every interval from v_intervalle
      vecteur old_v(v_intervalle);
      v_intervalle.clear();
      const_iterateur it=old_v.begin(),itend=old_v.end();
      for (;it!=itend;++it){
	if ( (it->type!=_VECT) || (it->subtype!=_LINE__VECT) || (it->_VECTptr->size()!= 2) )
	  setsizeerr();
	gen i_inf(it->_VECTptr->front()),i_sup(it->_VECTptr->back());
	bool a_greater_sup=ck_is_greater(a,i_sup,contextptr);
	if (a_greater_sup){
	  if (direction<0)
	    v_intervalle.push_back(*it);
	  continue;
	}
	bool a_greater_inf=ck_is_greater(a,i_inf,contextptr);
	if (!a_greater_inf){
	  if (direction>0)
	    v_intervalle.push_back(*it);
	  continue;
	}
	if (direction>0)
	  v_intervalle.push_back(gen(makevecteur(a,i_sup),_LINE__VECT));
	else
	  v_intervalle.push_back(gen(makevecteur(i_inf,a),_LINE__VECT));
      }
    }
    return gen(makevecteur(v0,v_intervalle,v_excluded),_ASSUME__VECT);
  }
  // returns the assumed idnt name
  // used if assumptions are in OR conjonction
  gen assumesymbolic(const gen & a,gen idnt_must_be,GIAC_CONTEXT){
    if (a.type==_IDNT)
      return a._IDNTptr->eval(eval_level(contextptr),a,contextptr);
    if ( (a.type!=_SYMB) || (a._SYMBptr->feuille.type!=_VECT) )
      setsizeerr();
    while (idnt_must_be.type==_SYMB){
      idnt_must_be=idnt_must_be._SYMBptr->feuille;
      if ( (idnt_must_be.type==_VECT) && !(idnt_must_be._VECTptr->empty()) )
	idnt_must_be=idnt_must_be._VECTptr->front();
    }
    unary_function_ptr s(a._SYMBptr->sommet);
    vecteur v(*a._SYMBptr->feuille._VECTptr);
    int l=v.size();
    if (!l)
      setsizeerr();
    gen arg0(v.front()),arg1(v.back()),hyp(undef);
    if (s==at_sto){
      gen tmp(arg0);
      arg0=arg1;
      arg1=tmp;
    }
    if (s==at_and){
      assumesymbolic(arg0,0,contextptr);
      return assumesymbolic(arg1,0,contextptr);
    }
    if (s==at_ou){
      gen a0(assumesymbolic(arg0,0,contextptr));
      return assumesymbolic(arg1,a0,contextptr);
    }
    if (arg0.type!=_IDNT)
      arg0=arg0.eval(eval_level(contextptr),contextptr);
    if ( (arg0.type!=_IDNT) || (!is_zero(idnt_must_be) && (arg0!=idnt_must_be) ) )
      setsizeerr();
    bool or_assumption= !is_zero(idnt_must_be) && (arg0==idnt_must_be);
    vecteur last_hyp;
    arg1=arg0._IDNTptr->eval(eval_level(contextptr),arg0,contextptr);
    if ( (arg1.type!=_VECT) || (arg1.subtype!=_ASSUME__VECT) )
      last_hyp=makevecteur(vecteur(0),vecteur(0));
    else
      last_hyp=*arg1._VECTptr;
    if (l==2){
      if (s==at_sto)
	arg1=v[0].eval(eval_level(contextptr),contextptr);
      else
	arg1=v[1].eval(eval_level(contextptr),contextptr);
      gen borne_inf(gnuplot_xmin),borne_sup(gnuplot_xmax),pas;
      if ( (s==at_equal || s==at_same) || (s==at_sto) ){     
	// ex: assume(a=[1.7,1.1,2.3])
	if (arg1.type==_VECT && arg1._VECTptr->size()>=3){
	  vecteur vtmp=*arg1._VECTptr;
	  borne_inf=evalf_double(vtmp[1],eval_level(contextptr),contextptr);
	  borne_sup=evalf_double(vtmp[2],eval_level(contextptr),contextptr);
	  pas=(borne_sup-borne_inf)/100;
	  if (vtmp.size()>3)
	    pas=evalf_double(vtmp[3],eval_level(contextptr),contextptr);
	  arg1=evalf_double(vtmp[0],eval_level(contextptr),contextptr);
	}
	gen tmp=arg1.type;
	tmp.subtype=1;
	hyp=gen(makevecteur(tmp,arg1),_ASSUME__VECT);
      }
      if (s==at_inferieur_strict) // ex: assume(a<1.7)
	hyp=doubleassume_and(last_hyp,arg1,-2,or_assumption,contextptr);
      if (s==at_inferieur_egal) 
	hyp=doubleassume_and(last_hyp,arg1,-1,or_assumption,contextptr);
      if (s==at_superieur_strict)
	hyp=doubleassume_and(last_hyp,arg1,2,or_assumption,contextptr);
      if (s==at_superieur_egal) 
	hyp=doubleassume_and(last_hyp,arg1,1,or_assumption,contextptr);
      if (!is_undef(hyp)){
	sto(hyp,arg0,contextptr); 
	if ( (s==at_equal || s==at_same) || (s==at_sto) )
	  return _parameter(makevecteur(arg0,borne_inf,borne_sup,arg1,pas));
	return arg0;
      }
    }
    setsizeerr();
    return 0;
  }
  void purge_assume(const gen & a,GIAC_CONTEXT){
    if (a.type==_SYMB && (a._SYMBptr->sommet==at_and || a._SYMBptr->sommet==at_ou || a._SYMBptr->sommet==at_inferieur_strict || a._SYMBptr->sommet==at_inferieur_egal || a._SYMBptr->sommet==at_superieur_egal || a._SYMBptr->sommet==at_superieur_strict) ){
      purge_assume(a._SYMBptr->feuille,contextptr);
      return;
    }
    if (a.type==_VECT && !a._VECTptr->empty())
      purge_assume(a._VECTptr->front(),contextptr);
    else
      _purge(a,contextptr);
  }
  gen giac_assume(const gen & a,GIAC_CONTEXT){
    if ( (a.type==_VECT) && (a._VECTptr->size()==2) ){
      gen a1(a._VECTptr->front()),a2(a._VECTptr->back());
      if (a2.type==_INT_){
	// assume(a,real) for example
	a2.subtype=1;
	sto(gen(makevecteur(a2),_ASSUME__VECT),a1,contextptr);
	return a2;
      }
      if (a2==at_real){
	a2=_DOUBLE_;
	a2.subtype=1;
	sto(gen(makevecteur(a2),_ASSUME__VECT),a1,contextptr);
	return a2;
      }
      if ( (a2.type==_FUNC) && (*a2._FUNCptr==at_ou) ){
	purge_assume(a1,contextptr);
	return assumesymbolic(a1,a1,contextptr);
      }
    }
    purge_assume(a,contextptr);
    return assumesymbolic(a,0,contextptr);
  }
  const string giac_assume_s("assume");
  unary_function_eval giac__assume(&giac::giac_assume,giac_assume_s);
  unary_function_ptr at_assume (&giac__assume,_QUOTE_ARGUMENTS,true);

  gen giac_additionally(const gen & a,GIAC_CONTEXT){
    if ( (a.type==_VECT) && (a._VECTptr->size()==2) ){
      gen a1(a._VECTptr->front()),a2(a._VECTptr->back());
      if (a1.type!=_IDNT)
	setsizeerr();
      gen a1val=a1._IDNTptr->eval(1,a1,contextptr);
      if (a1val.type==_VECT && a1val.subtype==_ASSUME__VECT && !a1val._VECTptr->empty()){
	if (a2.type==_INT_){
	  // assume(a,real) for example
	  a2.subtype=1;
	  a1val._VECTptr->front()=a2;
	  return a2;
	}
	if (a2==at_real){
	  a2=_DOUBLE_;
	  a2.subtype=1;
	  a1val._VECTptr->front()=a2;
	  return a2;
	}
      }
      else
	giac_assume(a,contextptr);
    }    
    return assumesymbolic(a,0,contextptr);
  }
  const string giac_additionally_s("additionally");
  unary_function_eval giac__additionally(&giac::giac_additionally,giac_additionally_s);
  unary_function_ptr at_additionally (&giac__additionally,_QUOTE_ARGUMENTS,true);

  // multiargs
  symbolic symb_plus(const gen & a,const gen & b){
    if (a.is_symb_of_sommet(at_plus) && !is_inf(a._SYMBptr->feuille)){
      if (b.is_symb_of_sommet(at_plus) && !is_inf(b._SYMBptr->feuille))
	  return symbolic(at_plus,mergevecteur(*(a._SYMBptr->feuille._VECTptr),*(b._SYMBptr->feuille._VECTptr)));
	else
	  return symbolic(*a._SYMBptr,b);
    }
    return symbolic(at_plus,makevecteur(a,b));
  }

  inline bool plus_idnt_symb(const gen & a){
    return (a.type==_IDNT && a._IDNTptr->name!=_IDNT_undef().name && a._IDNTptr->name!=_IDNT_infinity().name) || (a.type==_SYMB && !is_inf(a) && (a._SYMBptr->sommet==at_prod || a._SYMBptr->sommet==at_pow || a._SYMBptr->sommet==at_neg));
  }
  
  inline bool idnt_symb_int(const gen & b){
    return (b.type==_INT_ && b.val!=0) || b.type==_ZINT || (b.type==_SYMB && !is_inf(b)) || (b.type==_IDNT && b._IDNTptr->name!=_IDNT_undef().name && b._IDNTptr->name!=_IDNT_infinity().name);
  }

  gen _plus(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT){
      if ((args.type==_IDNT) && (*args._IDNTptr==_IDNT_infinity()))
	return plus_inf;
      return args;
    }
    iterateur it=args._VECTptr->begin(), itend=args._VECTptr->end();
    if (it==itend)
      return zero;
    const gen & a=*it;
    ++it;
    if (itend==it)
      return a;
    const gen & b=*it;
    ++it;
    if (it==itend){
      // improve: if a is an idnt/symb and b also do not rebuild the vector
      if (idnt_symb_int(b) && plus_idnt_symb(a))
	return new symbolic(at_plus,args);
      if (idnt_symb_int(a) && plus_idnt_symb(b))
	return new symbolic(at_plus,args);
      return operator_plus(a,b,contextptr);
    }
    gen sum(operator_plus(a,b,contextptr));
    for (;it!=itend;++it){
      if (sum.type==_SYMB && sum._SYMBptr->sommet==at_plus && sum._SYMBptr->feuille.type==_VECT && sum._SYMBptr->feuille._VECTptr->size()>1){
	// Add remaining elements to the symbolic sum, check float/inf/undef
	vecteur * vptr=new vecteur(*sum._SYMBptr->feuille._VECTptr);
	vptr->reserve(vptr->size()+itend-it);
	for (;it!=itend;++it){
	  if (it->type==_USER && vptr->front().type==_USER){
	    vptr->front()=operator_plus(vptr->front(),*it,contextptr);
	    continue;
	  }
	  if (it->type<=_POLY && vptr->back().type<=_POLY)
	    vptr->back()=operator_plus(vptr->back(),*it,contextptr);
	  else {
	    if (is_inf(*it) || is_undef(*it) || (it->type==_SYMB && it->_SYMBptr->sommet==at_plus))
	      break;
	    if (!is_zero(*it))
	      vptr->push_back(*it);
	  }
	}
	if (is_zero(vptr->back()))
	  vptr->pop_back();
	if (vptr->size()==1)
	  sum=vptr->front();
	else
	  sum=symbolic(at_plus,vptr);
	if (it==itend)
	  break;
      }
      operator_plus_eq(sum ,*it,contextptr);
    }
    return sum;
  }

  unary_function_ptr _D_at_plus (int i) {
    return at_one;
  }
  partial_derivative_multiargs D_at_plus(&_D_at_plus);
  const string _plus_s("+");
  unary_function_eval __plus(&_plus,&D_at_plus,_plus_s,&printsommetasoperator);
  unary_function_ptr at_plus (&__plus);

  inline bool prod_idnt_symb(const gen & a){
    return (a.type==_IDNT && a._IDNTptr->name!=_IDNT_undef().name && a._IDNTptr->name!=_IDNT_infinity().name) || (a.type==_SYMB && !is_inf(a) && (a._SYMBptr->sommet==at_plus || a._SYMBptr->sommet==at_pow || a._SYMBptr->sommet==at_neg));
  }
  
  symbolic symb_prod(const gen & a,const gen & b){
    return symbolic(at_prod,makevecteur(a,b));
  }
  gen _prod(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      return args;
    iterateur it=args._VECTptr->begin(), itend=args._VECTptr->end();
    gen prod(1);
    /*
    if (it==itend)
      return 1;
    const gen & a=*it;
    ++it;
    if (itend==it)
      return a;
    const gen & b=*it;
    ++it;
    if (it==itend){
      // improve: if a is an idnt/symb and b also do not rebuild the vector
      if (idnt_symb_int(b) && prod_idnt_symb(a))
	return new symbolic(at_prod,args);
      if (idnt_symb_int(a) && prod_idnt_symb(b))
	return new symbolic(at_prod,args);
      return operator_plus(a,b,contextptr);
    }
    gen prod(operator_times(a,b,contextptr));
    */
    for (;it!=itend;++it){
      if ( (it->type==_SYMB) && (it->_SYMBptr->sommet==at_inv) && (it->_SYMBptr->feuille.type!=_VECT) )
	prod = rdiv(prod,it->_SYMBptr->feuille);
      else
	prod = operator_times(prod,*it,contextptr);
    }
    return prod;
  }
  unary_function_ptr _D_at_prod (int i) {
    vector<int> v;
    v.push_back(i);
    vector<unary_function_ptr> w;
    w.push_back(at_prod);
    w.push_back(unary_function_innerprod(v));
    return unary_function_compose(w);
  }
  partial_derivative_multiargs D_at_prod(&_D_at_prod);
  const string _prod_s("*");
  unary_function_eval __prod(&_prod,&D_at_prod,_prod_s,&printsommetasoperator);
  unary_function_ptr at_prod (&__prod);

  std::string cprintaspow(const gen & feuille,const string & sommetstr_orig,GIAC_CONTEXT){
    gen f(feuille);
    if (f.type==_VECT)
      f.subtype=_SEQ__VECT;
    return "pow("+f.print(contextptr)+")";
  }
  symbolic symb_pow(const gen & a,const gen & b){
    return symbolic(at_pow,makevecteur(a,b));
  }
  gen _pow(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      return args;
    vecteur & v = *args._VECTptr;
    if (v.size()!=2)
      setsizeerr("bad pow "+args.print(contextptr));
    const gen & a =v.front();
    const gen & b =v.back();
    // fast check for monomials, do not recreate the vector
    if (b.type==_INT_){
#ifdef COMPILE_FOR_STABILITY
      if (b.val > FACTORIAL_SIZE_LIMIT)
	setstabilityerr("pow");
#endif
      if (b.val==1)
	return a;
      if (a.type==_IDNT){
	if (a==undef)
	  return a;
	if (a!=unsigned_inf)
	  return b.val?symbolic(at_pow,args):gen(1);
      }
      if (a.type==_SYMB && !is_inf(a) && (a._SYMBptr->sommet==at_plus || a._SYMBptr->sommet==at_prod)){
	return b.val?symbolic(at_pow,args):gen(1);
      }
    }
    return pow(a,b,contextptr);
  }
  gen d1_pow(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      setsizeerr();
    vecteur & v=*args._VECTptr;
    if (v.size()!=2)
      setsizeerr("bad pow "+args.print(contextptr));
    if (v[1].type<=_REAL)
      return v[1]*pow(v[0],v[1]-1,contextptr);
    else
      return v[1]/v[0]*symbolic(at_pow,v);
  }
  gen d2_pow(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      setsizeerr();
    vecteur & v=*args._VECTptr;
    if (v.size()!=2)
      setsizeerr("bad pow "+args.print(contextptr));
    return ln(v[0],contextptr)*pow(v[0],v[1],contextptr);
  }
  unary_function_ptr D1_pow(unary_function_eval(&d1_pow,""));
  unary_function_ptr D2_pow(unary_function_eval(&d2_pow,""));
  unary_function_ptr d_pow(int i){
    if (i==1)
      return D1_pow;
    if (i==2)
      return D2_pow;
    setsizeerr();
    return 0;
  }
  partial_derivative_multiargs D_pow(&d_pow);
  const string _pow_s("^");
  unary_function_eval __pow(&_pow,&D_pow,_pow_s,&printsommetasoperator,0);
  unary_function_ptr at_pow (&__pow);

  // print power like a^b (args==1), pow(a,b) (args==0) or a**b (args==-1)
  gen _printpow(const gen & args,GIAC_CONTEXT){
    if (is_zero(args)){
      __pow.printsommet=&cprintaspow;
      return string2gen("pow",false);
    }
    else {
      __pow.printsommet=&printsommetasoperator;
      if (is_minus_one(args))
	__pow.s="**";
      else
	__pow.s="^";
      return string2gen(__pow.s,false);
    }
  }
  const string _printpow_s("printpow");
  unary_function_eval __printpow(&_printpow,_printpow_s);
  unary_function_ptr at_printpow (&__printpow,0,true);

  symbolic symb_powmod(const gen & a,const gen & b,const gen & n){
    return symbolic(at_powmod,makevecteur(a,b,n));
  }
  symbolic symb_powmod(const gen & a){
    return symbolic(at_powmod,a);
  }
  gen findmod(const gen & g){
    if (g.type==_MOD)
      return *(g._MODptr+1);
    if (g.type==_VECT){
      gen res;
      const_iterateur it=g._VECTptr->begin(),itend=g._VECTptr->end();
      for (;it!=itend;++it){
	res=findmod(*it);
	if (!is_zero(res))
	  return res;
      }
    }
    if (g.type==_SYMB)
      return findmod(g._SYMBptr->feuille);
    return 0;
  }
  gen _powmod(const gen & args,GIAC_CONTEXT){
    int s;
    if ( args.type!=_VECT || (s=args._VECTptr->size())<3 )
      return symb_powmod(args);
    vecteur & v = *args._VECTptr;
    gen a=v.front();
    gen n=v[1];
    if (!is_integer(n))
      return symb_powmod(args);
    gen m=v[2];
    if (s==3 && m.type!=_SYMB) // a^n mod m
      return powmod(v.front(),v[1],m);
    // powmod(a_x%m,n,p_x) or powmod(a_x,n,m,p_x,x)
    // a^n mod p,m or m,p or a^n mod p,m,x or m,p,x wrt var x
    gen var(vx_var),p;
    bool modafter=false;
    if (s==3){ // m.type==_SYMB
      p=unmod(m);
      // find m inside a
      m=findmod(a);
      a=unmod(a);
      modafter=true;
    }
    else { // s>3
      p=v[3];
      if (is_integer(p)){
	p=v[2]; m=v[3]; 
      }
    }
    if (s>=5)
      var=v[4];
    vecteur lv(1,var);
    lvar(v,lv);
    if (lv.size()!=1)
      setsizeerr("Too many variables "+gen(lv).print(contextptr));
    gen aa=e2r(a,lv,contextptr),aan,aad,bb=e2r(p,lv,contextptr),bbn,bbd;
    fxnd(aa,aan,aad);
    if ( (aad.type==_POLY) && (aad._POLYptr->lexsorted_degree() ) )
      setsizeerr();
    fxnd(bb,bbn,bbd);
    if ( (bbd.type==_POLY) && (bbd._POLYptr->lexsorted_degree() ) )
      setsizeerr(); 
    if (bbn.type!=_POLY)
      setsizeerr();
    modpoly A;
    if (aan.type==_POLY)
      A=polynome2poly1(*aan._POLYptr);
    else
      A.push_back(aan);
    modpoly B=polynome2poly1(*bbn._POLYptr);
    environment env;
    env.moduloon=true;
    env.modulo=m;
    modpoly res=powmod(A,n,B,&env);
    polynome R(poly12polynome(res));
    if (modafter)
      modularize(R,m);
    gen Res=r2e(R,lv,contextptr)/pow(r2e(aad,lv,contextptr),n,contextptr);
    return Res;
  }
  const string _powmod_s("powmod");
  unary_function_eval __powmod(&giac::_powmod,_powmod_s);
  unary_function_ptr at_powmod (&__powmod,0,true);

  symbolic symb_inferieur_strict(const gen & a,const gen & b){
    return symbolic(at_inferieur_strict,makevecteur(a,b));
  }
  symbolic symb_inferieur_strict(const gen & a){
    return symbolic(at_inferieur_strict,a);
  }
  gen _inferieur_strict(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      return symb_inferieur_strict(args);
    return inferieur_strict(args._VECTptr->front(),args._VECTptr->back(),contextptr);
  }
  const string _inferieur_strict_s("<");
  unary_function_eval __inferieur_strict(&giac::_inferieur_strict,_inferieur_strict_s,&printsommetasoperator,&texprintsommetasoperator);
  unary_function_ptr at_inferieur_strict (&__inferieur_strict);

  symbolic symb_inferieur_egal(const gen & a,const gen & b){
    return symbolic(at_inferieur_egal,makevecteur(a,b));
  }
  symbolic symb_inferieur_egal(const gen & a){
    return symbolic(at_inferieur_egal,a);
  }
  gen _inferieur_egal(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      return symb_inferieur_egal(args);
    return inferieur_egal(args._VECTptr->front(), args._VECTptr->back(),contextptr);
  }
  const string _inferieur_egal_s("<=");
  string texprintasinferieur_egal(const gen & g,const string & s,GIAC_CONTEXT){
    return texprintsommetasoperator(g,"\\leq ",contextptr);
  }
  unary_function_eval __inferieur_egal(&giac::_inferieur_egal,_inferieur_egal_s,&printsommetasoperator,&texprintasinferieur_egal);
  unary_function_ptr at_inferieur_egal (&__inferieur_egal);

  symbolic symb_superieur_strict(const gen & a,const gen & b){
    return symbolic(at_superieur_strict,makevecteur(a,b));
  }
  symbolic symb_superieur_strict(const gen & a){
    return symbolic(at_superieur_strict,a);
  }
  gen _superieur_strict(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      return symb_superieur_strict(args);
    return superieur_strict(args._VECTptr->front(),args._VECTptr->back(),contextptr);
  }
  const string _superieur_strict_s(">");
  unary_function_eval __superieur_strict(&giac::_superieur_strict,_superieur_strict_s,&printsommetasoperator,&texprintsommetasoperator);
  unary_function_ptr at_superieur_strict (&__superieur_strict);

  symbolic symb_superieur_egal(const gen & a,const gen & b){
    return symbolic(at_superieur_egal,makevecteur(a,b));
  }
  symbolic symb_superieur_egal(const gen & a){
    return symbolic(at_superieur_egal,a);
  }
  gen _superieur_egal(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      return symb_superieur_egal(args);
    return superieur_egal(args._VECTptr->front(), args._VECTptr->back(),contextptr);
  }
  const string _superieur_egal_s(">=");
  string texprintassuperieur_egal(const gen & g,const string & s,GIAC_CONTEXT){
    return texprintsommetasoperator(g,"\\geq ",contextptr);
  }
  unary_function_eval __superieur_egal(&giac::_superieur_egal,_superieur_egal_s,&printsommetasoperator,&texprintassuperieur_egal);
  unary_function_ptr at_superieur_egal (&__superieur_egal);

  string printasdifferent(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if (xcas_mode(contextptr) > 0)
      return printsommetasoperator(feuille,"<>",contextptr);
    else
      return printsommetasoperator(feuille,sommetstr,contextptr);
  }
  symbolic symb_different(const gen & a,const gen & b){
    return symbolic(at_different,makevecteur(a,b));
  }
  symbolic symb_different(const gen & a){
    return symbolic(at_different,a);
  }
  gen _different(const gen & args){
    if (args.type!=_VECT)
      return symb_different(args);
    return args._VECTptr->front() != args._VECTptr->back();
  }
  const string _different_s("!=");
  unary_function_unary __different(&giac::_different,_different_s,&printasdifferent);
  unary_function_ptr at_different (&__different);

  string printasof_(const gen & feuille,const string & sommetstr,int format,GIAC_CONTEXT){
    if ( (feuille.type!=_VECT) || (feuille._VECTptr->size()!=2) )
      return sommetstr+'('+gen2string(feuille,format,contextptr)+')';
    string s=print_with_parenthesis_if_required(feuille._VECTptr->front(),format,contextptr)+'(';
    gen & g=feuille._VECTptr->back();
    if (format==0 && g.type==_VECT && g.subtype==_SEQ__VECT)
      return s+printinner_VECT(*g._VECTptr,_SEQ__VECT,contextptr)+')';
    else
      return s+gen2string(g,format,contextptr)+')';
  }
  string texprintasof(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    return printasof_(feuille,sommetstr,1,contextptr);
  }
  string printasof(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    return printasof_(feuille,sommetstr,0,contextptr);
  }
  symbolic symb_of(const gen & a,const gen & b){
    if (b.type==_VECT && b.subtype==_SEQ__VECT && b._VECTptr->size()==1)
      return symbolic(at_of,makevecteur(a,b._VECTptr->front()));
    return symbolic(at_of,makevecteur(a,b));
  }
  symbolic symb_of(const gen & a){
    return symbolic(at_of,a);
  }
  gen _of(const gen & args,const context * contextptr){
    if (args.type!=_VECT)
      return symb_of(args);
    gen qf=args._VECTptr->front();
    gen b=args._VECTptr->back().eval(eval_level(contextptr),contextptr);
    gen f=qf.eval(eval_level(contextptr),contextptr);
    if ( (f.type==_SYMB) && (f._SYMBptr->sommet==at_program) && (qf.type==_IDNT)){
      gen tmp=f._SYMBptr->feuille;
      if (tmp.type!=_VECT)
	setsizeerr();
      (*tmp._VECTptr)[1]=b;
      return _program(tmp,qf,contextptr);
    }
    return f(b,contextptr);
  }
  const string _of_s("of");
  unary_function_eval __of(&giac::_of,_of_s,&printasof,&texprintasof);
  unary_function_ptr at_of (&__of,_QUOTE_ARGUMENTS);

  string gen2string(const gen & g,int format,GIAC_CONTEXT){
    if (format==1) 
      return gen2tex(g,contextptr); 
    else 
      return g.print(contextptr);
  }

  string print_with_parenthesis_if_required(const gen & g,int format,GIAC_CONTEXT){
    if (g.type==_SYMB || g.type==_FRAC || g.type==_CPLX || (g.type==_VECT && g.subtype==_SEQ__VECT) )
      return '('+gen2string(g,format,contextptr)+')';
    else
      return gen2string(g,format,contextptr);
  }
  
  string printasat_(const gen & feuille,const string & sommetstr,int format,GIAC_CONTEXT){
    if ( (feuille.type!=_VECT) || (feuille._VECTptr->size()!=2) )
      return sommetstr+'('+gen2string(feuille,format,contextptr)+')';
    vecteur & v=*feuille._VECTptr;
    if (xcas_mode(contextptr) > 0){
      gen indice;
      if (v.back().type==_VECT)
	indice=v.back()+vecteur(v.size(),plus_one);
      else
	indice=v.back()+plus_one;
      string s;
      return print_with_parenthesis_if_required(v.front(),format,contextptr)+'['+gen2string(indice,format,contextptr)+']';
    }
    else
      return print_with_parenthesis_if_required(feuille._VECTptr->front(),format,contextptr)+'['+gen2string(feuille._VECTptr->back(),format,contextptr)+']';
  }

  string printasat(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    return printasat_(feuille,sommetstr,0,contextptr);
  }
  string texprintasat(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    return printasat_(feuille,sommetstr,1,contextptr);
  }
  symbolic symb_at(const gen & a,const gen & b,GIAC_CONTEXT){
    if (xcas_mode(contextptr)){
      gen bb;
      if (b.type==_VECT)
	bb=b-vecteur(b._VECTptr->size(),plus_one);
      else
	bb=b-plus_one;
      return symbolic(at_at,makevecteur(a,bb));
    }
    else
      return symbolic(at_at,makevecteur(a,b));
  }
  symbolic symb_at(const gen & a){
    return symbolic(at_at,a);
  }
  gen _at(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      return symb_at(args);
    vecteur & v=*args._VECTptr;
    if (v.size()!=2)
      setsizeerr();
    gen a=v.front().eval(eval_level(contextptr),contextptr);
    gen b=v.back().eval(eval_level(contextptr),contextptr);
    if (a.type==_MAP){
      gen_map::const_iterator it=a._MAPptr->find(b),itend=a._MAPptr->end();
      if (it!=itend)
	return it->second;
      return symb_at(makevecteur(v.front(),b));
    }
    return a.operator_at(b,contextptr);
  }
  const string _at_s("at");
  unary_function_eval __at(&giac::_at,_at_s,&printasat,&texprintasat);
  unary_function_ptr at_at (&__at,_QUOTE_ARGUMENTS);

  gen _table(const gen & arg){
    vecteur v(gen2vecteur(arg));
    const_iterateur it=v.begin(),itend=v.end();
    gen_map m(ptr_fun(islesscomplexthanf));
    for (;it!=itend;++it){
      if (it->is_symb_of_sommet(at_equal)){
	gen & f =it->_SYMBptr->feuille;
	if (f.type==_VECT && f._VECTptr->size()==2){
	  vecteur & w=*f._VECTptr;
	  m[w.front()]=w.back();
	}
      }
    }
    return m;
  }
  const string _table_s("table");
  unary_function_unary __table(&giac::_table,_table_s);
  unary_function_ptr at_table (&__table,0,true);

  string printasand(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if (xcas_mode(contextptr) > 0)
      return printsommetasoperator(feuille," and ",contextptr);
    else
      return "("+printsommetasoperator(feuille,sommetstr,contextptr)+")";
  }
  string texprintasand(const gen & g,const string & s,GIAC_CONTEXT){
    return texprintsommetasoperator(g,"\\mbox{ and }",contextptr);
  }
  symbolic symb_and(const gen & a,const gen & b){
    return symbolic(at_and,makevecteur(a,b));
  }
  gen and2(const gen & a,const gen & b){
    return a && b;
  }
  gen _and(const gen & arg,GIAC_CONTEXT){
    if (arg.type==_VECT && arg.subtype==_SEQ__VECT && arg._VECTptr->size()==2)
      return apply(equaltosame(arg._VECTptr->front()).eval(eval_level(contextptr),contextptr),equaltosame(arg._VECTptr->back()).eval(eval_level(contextptr),contextptr),and2);
    gen args=apply(arg,equaltosame);
    if (args.type!=_VECT || args._VECTptr->empty())
      return args.eval(eval_level(contextptr),contextptr);
    vecteur::const_iterator it=args._VECTptr->begin(),itend=args._VECTptr->end();
    gen res=eval(*it,eval_level(contextptr),contextptr);
    ++it;
    for (;it!=itend;++it){
      if (res.type==_INT_ && res.val==0)
	return res;
      res = res && eval(*it,eval_level(contextptr),contextptr);
    }
    return res;
  }
  const string _and_s(" && ");
  unary_function_eval __and(&giac::_and,_and_s,&printasand,&texprintasand);
  unary_function_ptr at_and (&__and,_QUOTE_ARGUMENTS);

  string texprintasor(const gen & g,const string & s,GIAC_CONTEXT){
    return texprintsommetasoperator(g,"\\mbox{ or }",contextptr);
  }
  string printasor(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if (xcas_mode(contextptr) > 0)
      return printsommetasoperator(feuille," or ",contextptr);
    else
      return "("+printsommetasoperator(feuille,sommetstr,contextptr)+")";
  }
  symbolic symb_ou(const gen & a,const gen & b){
    return symbolic(at_ou,makevecteur(a,b));
  }
  gen ou2(const gen & a,const gen & b){
    return a || b;
  }
  gen _ou(const gen & arg,GIAC_CONTEXT){
    if (arg.type==_VECT && arg.subtype==_SEQ__VECT && arg._VECTptr->size()==2)
      return apply(equaltosame(arg._VECTptr->front()).eval(eval_level(contextptr),contextptr),equaltosame(arg._VECTptr->back()).eval(eval_level(contextptr),contextptr),ou2);
    gen args=apply(arg,equaltosame);
    if (args.type!=_VECT || args._VECTptr->empty())
      return eval(args,eval_level(contextptr),contextptr);
    vecteur::const_iterator it=args._VECTptr->begin(),itend=args._VECTptr->end();
    gen res=eval(*it,eval_level(contextptr),contextptr);
    ++it;
    for (;it!=itend;++it){
      if (res.type==_INT_ && res.val)
	return res;
      res = res || eval(*it,eval_level(contextptr),contextptr);
    }
    return res;
  }
  const string _ou_s(" || ");
  unary_function_eval __ou(&giac::_ou,_ou_s,&printasor,&texprintasor);
  unary_function_ptr at_ou (&__ou,_QUOTE_ARGUMENTS);

  gen xor2(const gen & a,const gen & b){
    return is_zero(a) ^ is_zero(b);
  }
  gen _xor(const gen & arg,GIAC_CONTEXT){
    if (arg.type==_VECT && arg.subtype==_SEQ__VECT && arg._VECTptr->size()==2)
      return apply(equaltosame(arg._VECTptr->front()).eval(eval_level(contextptr),contextptr),equaltosame(arg._VECTptr->back()).eval(eval_level(contextptr),contextptr),xor2);
    gen args=eval(apply(arg,equaltosame),eval_level(contextptr),contextptr);
    if (args.type!=_VECT)
      return args;
    vecteur::const_iterator it=args._VECTptr->begin(),itend=args._VECTptr->end();
    gen res=*it;
    ++it;
    for (;it!=itend;++it){
      if (is_zero(res))
	res=*it;
      else
	res = !(*it);
    }
    return res;
  }
  const string _xor_s(" xor ");
  unary_function_eval __xor(&giac::_xor,_xor_s,&printsommetasoperator,&texprintsommetasoperator);
  unary_function_ptr at_xor (&__xor,_QUOTE_ARGUMENTS);

  symbolic symb_min(const gen & a,const gen & b){
    return symbolic(at_min,makevecteur(a,b));
  }
  gen _min(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      return args;
    vecteur::const_iterator it=args._VECTptr->begin(),itend=args._VECTptr->end();
    if (ckmatrix(args))
      return apply(*it,*(it+1),contextptr,min);
    if (itend-it==2 && it->type==_VECT && (it+1)->type==_VECT )
      return matrix_apply(*it,*(it+1),contextptr,min);
    gen res=*it;
    ++it;
    for (;it!=itend;++it)
      res = min(res,*it,contextptr);
    return res;
  }
  const string _min_s("min");
  unary_function_eval giac__min(&giac::_min,_min_s);
  unary_function_ptr at_min (&giac__min,0,true);

  symbolic symb_max(const gen & a,const gen & b){
    return symbolic(at_max,makevecteur(a,b));
  }
  gen _max(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      return args;
    vecteur::const_iterator it=args._VECTptr->begin(),itend=args._VECTptr->end();
    if (itend==it)
      setdimerr();
    if (itend-it==1)
      return _max(*it,contextptr);
    if (itend-it==2 && ckmatrix(args))
      return apply(*it,*(it+1),contextptr,max);
    if (itend-it==2 && it->type==_VECT && (it+1)->type==_VECT )
      return matrix_apply(*it,*(it+1),contextptr,max);
    gen res=*it;
    ++it;
    for (;it!=itend;++it)
      res = max(res,*it,contextptr);
    return res;
  }
  const string _max_s("max");
  unary_function_eval giac__max(&giac::_max,_max_s);
  unary_function_ptr at_max (&giac__max,0,true);

  symbolic symb_gcd(const gen & a,const gen & b){
    return symbolic(at_gcd,makevecteur(a,b));
  }
  gen _gcd(const gen & args){
    if (args.type!=_VECT)
      return args;
    vecteur::const_iterator it=args._VECTptr->begin(),itend=args._VECTptr->end();
    if (ckmatrix(args) && itend-it==2)
      return apply(*it,*(it+1),gcd);
    gen res(0);
    for (;it!=itend;++it)
      res=gcd(res,*it);
    return res;
  }
  const string _gcd_s("gcd");
  unary_function_unary __gcd(&giac::_gcd,_gcd_s);
  unary_function_ptr at_gcd (&__gcd,0,true);

  symbolic symb_lcm(const gen & a,const gen & b){
    return symbolic(at_lcm,makevecteur(a,b));
  }
  gen _lcm(const gen & args){
    if (args.type!=_VECT)
      return args;
    vecteur::const_iterator it=args._VECTptr->begin(),itend=args._VECTptr->end();
    if (itend==it)
      return 1;
    if (ckmatrix(args) && itend-it==2)
      return apply(*it,*(it+1),lcm);
    gen res(*it);
    for (++it;it!=itend;++it)
      res=lcm(res,*it);
    return res;
  }
  const string _lcm_s("lcm");
  unary_function_unary __lcm(&giac::_lcm,_lcm_s);
  unary_function_ptr at_lcm (&__lcm,0,true);

  symbolic symb_egcd(const gen & a,const gen & b){
    return symbolic(at_egcd,makevecteur(a,b));
  }
  gen _egcd(const gen & args,GIAC_CONTEXT){
    if ( (args.type!=_VECT) || args._VECTptr->empty() )
      setsizeerr();
    vecteur & a = *args._VECTptr;
    if ( (a.front().type==_VECT) && (a.back().type==_VECT) ){
      vecteur u,v,d;
      egcd(*a.front()._VECTptr,*a.back()._VECTptr,0,u,v,d);
      return gen(makevecteur(u,v,d),_POLY1__VECT);
    }
    vecteur lv;
    if (a.size()==3)
      lv=vecteur(1,a[2]);
    else
      lv=vecteur(1,vx_var);
    lvar(args,lv);
    gen aa=e2r(a[0],lv,contextptr),aan,aad,bb=e2r(a[1],lv,contextptr),bbn,bbd;
    fxnd(aa,aan,aad);
    if ( (aad.type==_POLY) && (aad._POLYptr->lexsorted_degree() ) )
      setsizeerr();
    fxnd(bb,bbn,bbd);
    if ( (bbd.type==_POLY) && (bbd._POLYptr->lexsorted_degree() ) )
      setsizeerr(); 
    gen u,v,d;
    if ( (aan.type==_POLY) && (bbn.type==_POLY) ){
      polynome un(aan._POLYptr->dim),vn(aan._POLYptr->dim),dn(aan._POLYptr->dim);
      egcd(*aan._POLYptr,*bbn._POLYptr,un,vn,dn);
      u=un;
      v=vn;
      d=dn;
    }
    else {
      if (aan.type==_POLY){
	u=zero;
	v=plus_one;
	d=bbn;
      }
      else {
	u=plus_one;
	v=zero;
	d=aan;
      }
    }
    u=r2e(u*aad,lv,contextptr);
    v=r2e(v*bbd,lv,contextptr);
    d=r2e(d,lv,contextptr);
    return makevecteur(u,v,d);
  }
  const string _egcd_s("egcd");
  unary_function_eval __egcd(&giac::_egcd,_egcd_s);
  unary_function_ptr at_egcd (&__egcd,0,true);

  symbolic symb_iegcd(const gen & a,const gen & b){
    return symbolic(at_iegcd,makevecteur(a,b));
  }
  gen _iegcd(const gen & args){
    check_2d_vecteur(args);
    gen u,v,d;
    egcd(args._VECTptr->front(),args._VECTptr->back(),u,v,d);
    return makevecteur(u,v,d);
  }
  const string _iegcd_s("iegcd");
  unary_function_unary __iegcd(&giac::_iegcd,_iegcd_s);
  unary_function_ptr at_iegcd (&__iegcd,0,true);

  const string _bezout_entiers_s("bezout_entiers");
  unary_function_unary __bezout_entiers(&giac::_iegcd,_bezout_entiers_s);
  unary_function_ptr at_bezout_entiers (&__bezout_entiers,0,true);

  gen symb_equal(const gen & a,const gen & b){
    return symbolic(at_equal,gen(makevecteur(a,b),_SEQ__VECT));
  }
  gen _equal(const gen & a){
    if ((a.type!=_VECT) || (a._VECTptr->size()!=2))
      return equal(a,gen(vecteur(0),_SEQ__VECT));
    return equal( (*(a._VECTptr))[0],(*(a._VECTptr))[1] );
  }
  const string _equal_s("=");
  unary_function_unary __equal(&giac::_equal,_equal_s,&printsommetasoperator,&texprintsommetasoperator);
  unary_function_ptr at_equal (&__equal);

  string printassame(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if (xcas_mode(contextptr) > 0)
      return printsommetasoperator(feuille,"=",contextptr);
    else
      return "("+printsommetasoperator(feuille,sommetstr,contextptr)+")";
  }
  symbolic symb_same(const gen & a,const gen & b){
    return symbolic(at_same,makevecteur(a,b));
  }
  gen symb_same(const gen & a){
    return symbolic(at_same,a);
  }
  gen _same(const gen & a,GIAC_CONTEXT){
    if ((a.type!=_VECT) || (a._VECTptr->size()!=2))
      return symb_same(a);
    return operator_equal(a._VECTptr->front(),a._VECTptr->back(),contextptr);
  }
  const string _same_s("==");
  unary_function_eval __same(&giac::_same,_same_s,&printassame,&texprintsommetasoperator);
  unary_function_ptr at_same (&__same);

  // ******************
  // Arithmetic functions
  // *****************

  symbolic symb_smod(const gen & a,const gen & b){
    return symbolic(at_smod,makevecteur(a,b));
  }
  gen _smod(const gen & args){
    check_2d_vecteur(args);
    vecteur & v=*args._VECTptr;
    if (ckmatrix(v))
      return apply(v[0],v[1],smod);
    if (!is_cinteger(v.back()) )
      return symbolic(at_smod,args);
    return smod(args._VECTptr->front(),args._VECTptr->back());
  }
  const string _smod_s("smod");
  unary_function_unary __smod(&giac::_smod,_smod_s);
  unary_function_ptr at_smod (&__smod,0,true);

  symbolic symb_rdiv(const gen & a,const gen & b){
    return symbolic(at_rdiv,makevecteur(a,b));
  }
  gen _rdiv(const gen & args){
    check_2d_vecteur(args);
    return rdiv(args._VECTptr->front(),args._VECTptr->back());
  }
  const string _rdiv_s("rdiv");
  unary_function_unary __rdiv(&giac::_rdiv,_rdiv_s);
  unary_function_ptr at_rdiv (&__rdiv,0,true);

  gen unmod(const gen & g){
    if (g.type==_MOD)
      return *g._MODptr;
    if (g.type==_VECT)
      return apply(g,unmod);
    if (g.type==_SYMB)
      return symbolic(g._SYMBptr->sommet,unmod(g._SYMBptr->feuille));
    return g;
  }
  gen unmodunprod(const gen & g){
    gen h=unmod(g);
    if (h.is_symb_of_sommet(at_prod))
      h=_prod(h._SYMBptr->feuille,context0); // ok
    return h;
  }

  gen irem(const gen & a,const gen & b){
    gen q;
    return irem(a,b,q);
  }
  symbolic symb_irem(const gen & a,const gen & b){
    return symbolic(at_irem,makevecteur(a,b));
  }
  gen _normalmod(const gen & g,GIAC_CONTEXT);
  gen _irem(const gen & args,GIAC_CONTEXT){
    check_2d_vecteur(args);
    if (ckmatrix(args))
      return apply(args._VECTptr->front(),args._VECTptr->back(),irem);
    gen q;
    vecteur & v=*args._VECTptr;
    if (v.front().type==_SYMB){
      gen arg=v.front()._SYMBptr->feuille;
      if (v.front()._SYMBptr->sommet==at_pow && arg.type==_VECT && arg._VECTptr->size()==2 ){
	if (is_integer(arg._VECTptr->front()) && is_integer(arg._VECTptr->back()) )
	  return powmod(_irem(gen(makevecteur(arg._VECTptr->front(),v.back()),_SEQ__VECT),contextptr),arg._VECTptr->back(),v.back());
	return pow(_irem(gen(makevecteur(arg._VECTptr->front(),v.back()),_SEQ__VECT),contextptr),arg._VECTptr->back(),contextptr);
      }
      if (v.front()._SYMBptr->sommet==at_neg)
	return _irem(gen(makevecteur(simplifier((v.back()-1)*arg,contextptr),v.back()),_SEQ__VECT),contextptr);
      if (v.front()._SYMBptr->sommet==at_prod || v.front()._SYMBptr->sommet==at_plus){
	return v.front()._SYMBptr->sommet(_irem(gen(makevecteur(arg,v.back()),_SEQ__VECT),contextptr),contextptr);
      }
      if (v.front()._SYMBptr->sommet==at_inv){
	gen g=invmod(arg,v.back());
	if (is_positive(g,contextptr))
	  return g;
	else
	  return g+v.back();
      }
      arg=_normalmod(makevecteur(arg,v.back()),contextptr);
      return unmod(v.front()._SYMBptr->sommet(arg,contextptr));
    }
    if (v.front().type==_FRAC){
      gen g=invmod(v.front()._FRACptr->den,v.back());
      if (!is_positive(g,contextptr))
	g= g+v.back();
      return _irem(gen(makevecteur(v.front()._FRACptr->num*g,v.back()),_SEQ__VECT),contextptr);
    }
    if (v.front().type==_VECT){
      const_iterateur it=v.front()._VECTptr->begin(),itend=v.front()._VECTptr->end();
      vecteur res;
      for (;it!=itend;++it)
	res.push_back(_irem(gen(makevecteur(*it,v.back()),_SEQ__VECT),contextptr));
      return gen(res,v.front().subtype);
    }
    if (v.front().type==_IDNT)
      return v.front();
    if (!is_cinteger(v.front()) || !is_cinteger(v.back()) )
      return symbolic(at_irem,args);
    gen r=irem(v.front(),v.back(),q);
    if (is_strictly_positive(-r,contextptr)){
      if (is_strictly_positive(v.back(),contextptr)){
	r = r + v.back();
	q=q-1;
      }
      else {
	r = r - v.back();
	q=q+1;
      }
    }
    return r;
  }
  const string _irem_s("irem");
  unary_function_eval __irem(&giac::_irem,_irem_s);
  unary_function_ptr at_irem (&__irem,0,true);

  const string _mods_s("mods");
  unary_function_unary __mods(&giac::_smod,_mods_s);
  unary_function_ptr at_mods (&__mods,0,true);

  gen _quote_pow(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT || args._VECTptr->size()!=2)
      settypeerr();
    vecteur & v = *args._VECTptr;
    if (ckmatrix(v.front()))
      return pow(v.front(),v.back(),contextptr);
    return symbolic(at_pow,args);
  }
  const string _quote_pow_s("&^");
  unary_function_eval __quote_pow(&giac::_quote_pow,_quote_pow_s,&printsommetasoperator,&texprintsommetasoperator);
  unary_function_ptr at_quote_pow (&__quote_pow);

  symbolic symb_iquo(const gen & a,const gen & b){
    return symbolic(at_iquo,makevecteur(a,b));
  }
  gen Iquo(const gen & f,const gen & b){
    if (!is_cinteger(f) || !is_cinteger(b) )
      setsizeerr(); // return symbolic(at_iquo,args);
    if (is_zero(b))
      return 0;
    return (f-_irem(gen(makevecteur(f,b),_SEQ__VECT),context0))/b; // ok
  }
  gen _iquo(const gen & args){
    check_2d_vecteur(args);
    gen & f=args._VECTptr->front();
    gen & b=args._VECTptr->back();
    if (ckmatrix(args))
      return apply(f,b,iquo);
    return Iquo(f,b);
  }
  const string _iquo_s("iquo");
  unary_function_unary __iquo(&giac::_iquo,_iquo_s);
  unary_function_ptr at_iquo (&__iquo,0,true);

  vecteur iquorem(const gen & a,const gen & b){
    gen q,r;
    r=irem(a,b,q);
    return makevecteur(q,r);
  }
  symbolic symb_iquorem(const gen & a,const gen & b){
    return symbolic(at_iquorem,makevecteur(a,b));
  }
  gen _iquorem(const gen & args){
    check_2d_vecteur(args);
    vecteur & v=*args._VECTptr;
    if (!is_cinteger(v.front()) || !is_cinteger(v.back()) )
      return symbolic(at_iquorem,args);
    return iquorem(args._VECTptr->front(),args._VECTptr->back());
  }
  const string _iquorem_s("iquorem");
  unary_function_unary __iquorem(&giac::_iquorem,_iquorem_s);
  unary_function_ptr at_iquorem (&__iquorem,0,true);

  gen quorem(const gen & a,const gen & b){
    if ((a.type!=_VECT) || (b.type!=_VECT))
      return symb_quorem(a,b);
    vecteur q,r;
    environment * env=new environment;
    DivRem(*a._VECTptr,*b._VECTptr,env,q,r,true);
    delete env;
    return makevecteur(gen(q,_POLY1__VECT),gen(r,_POLY1__VECT));
  }
  symbolic symb_quorem(const gen & a,const gen & b){
    return symbolic(at_quorem,makevecteur(a,b));
  }
  gen _quorem(const gen & args,GIAC_CONTEXT){
    if ( (args.type!=_VECT) || (args._VECTptr->size()<2) )
      setsizeerr();
    vecteur & a =*args._VECTptr;
    if ( (a.front().type==_VECT) && (a.back().type==_VECT))
      return quorem(a.front(),a.back());
    if ( (a.front().type==_POLY) && (a.back().type==_POLY)){
      int dim=a.front()._POLYptr->dim;
      if (a.back()._POLYptr->dim!=dim)
	setdimerr();
      vecteur aa(polynome2poly1(*a.front()._POLYptr,1));
      vecteur bb(polynome2poly1(*a.back()._POLYptr,1));
      vecteur q,r;
      DivRem(aa,bb,0,q,r);
      return makevecteur(poly12polynome(q,1,dim),poly12polynome(r,1,dim));
    }
    vecteur lv;
    if (a.size()==3)
      lv=vecteur(1,unmodunprod(a[2]));
    else
      lv=vecteur(1,vx_var);
    lvar(args,lv);
    gen aa=e2r(a[0],lv,contextptr),aan,aad,bb=e2r(a[1],lv,contextptr),bbn,bbd;
    fxnd(aa,aan,aad);
    if ( (aad.type==_POLY) && (aad._POLYptr->lexsorted_degree() ) )
      setsizeerr();
    fxnd(bb,bbn,bbd);
    if ( (bbd.type==_POLY) && (bbd._POLYptr->lexsorted_degree() ) )
      setsizeerr();
    gen u,v;
    gen ad(r2e(aad,lv,contextptr));
    if ( (aan.type==_POLY) && (bbn.type==_POLY) ){
      vecteur aav(polynome2poly1(*aan._POLYptr,1)),bbv(polynome2poly1(*bbn._POLYptr,1)),un,vn;
      environment * env=new environment;
      DivRem(aav,bbv,env,un,vn);
      delete env;
      vecteur lvprime(lv.begin()+1,lv.end());
      u=rdiv(r2e(bbd,lv,contextptr),ad)*symb_horner(*r2e(un,lvprime,contextptr)._VECTptr,lv.front());
      v=inv(ad,contextptr)*symb_horner(*r2e(vn,lvprime,contextptr)._VECTptr,lv.front());
      return makevecteur(u,v);
    }
    else {
      if ( (bbn.type!=_POLY) || !bbn._POLYptr->lexsorted_degree() ){
	u=rdiv(aan,bbn);
	v=zero;
      }
      else {
	u=zero;
	v=aan;
      }
    }
    // aan=u*bbn+v -> aan/aad=u*bbd/aad * bbn/bbd +v/aad
    u=r2e(u*bbd,lv,contextptr);
    v=r2e(v,lv,contextptr);
    return makevecteur(rdiv(u,ad),rdiv(v,ad));
  }
  const string _quorem_s("quorem");
  unary_function_eval __quorem(&giac::_quorem,_quorem_s);
  unary_function_ptr at_quorem (&__quorem,0,true);

  symbolic symb_quo(const gen & a,const gen & b){
    return symbolic(at_quo,makevecteur(a,b));
  }
  gen _quo(const gen & args,GIAC_CONTEXT){
    return _quorem(args,contextptr)[0];
  }
  const string _quo_s("quo");
  unary_function_eval __quo(&giac::_quo,_quo_s);
  unary_function_ptr at_quo (&__quo,0,true);

  symbolic symb_rem(const gen & a,const gen & b){
    return symbolic(at_rem,makevecteur(a,b));
  }
  gen _rem(const gen & args,GIAC_CONTEXT){
    return _quorem(args,contextptr)[1];
  }
  const string _rem_s("rem");
  unary_function_eval __rem(&giac::_rem,_rem_s);
  unary_function_ptr at_rem (&__rem,0,true);

  gen double2gen(double d){
    mpz_t * m=(mpz_t *) malloc(sizeof(mpz_t)); //new mpz_t[1];
    mpz_init_set_d(*m,d);
    return m;
  }
  symbolic symb_floor(const gen & a){
    return symbolic(at_floor,a);
  }
  gen _floor(const gen & args,GIAC_CONTEXT){
    if (is_inf(args)||is_undef(args))
      return args;
    if (args.type==_VECT)
      return apply(args,contextptr,_floor);
    if (args.type==_CPLX)
      return _floor(*args._CPLXptr,contextptr)+cst_i*_floor(*(args._CPLXptr+1),contextptr);
    if ( (args.type==_INT_) || (args.type==_ZINT))
      return args;
    if (args.type==_FRAC){
      gen n=args._FRACptr->num,d=args._FRACptr->den;
      if ( ((n.type==_INT_) || (n.type==_ZINT)) && ( (d.type==_INT_) || (d.type==_ZINT)) ){
	if (is_positive(args,contextptr))
	  return iquo(n,d);
	else
	  return iquo(n,d)-1;
      }
    }
    vecteur l(lidnt(args));
    vecteur lnew=*evalf(l,1,contextptr)._VECTptr;
    gen tmp=subst(args,l,lnew,false,contextptr);
    if (tmp.type==_REAL){
      gen res=real2int(tmp,contextptr);
      if (is_strictly_positive(-tmp,contextptr) && !is_zero(res-tmp))
	return res-1;
      return res;
    }
    if (tmp.type!=_DOUBLE_)
      return symb_floor(args);
    return double2gen(giac_floor(tmp._DOUBLE_val));
  }
  gen taylor_floor (const gen & lim_point,const int ordre,const unary_function_ptr & f, int direction,gen & shift_coeff,GIAC_CONTEXT){
    if (ordre<0)
      return 0; // no symbolic preprocessing
    shift_coeff=0;
    gen l=_floor(lim_point,contextptr);
    if (l==lim_point){
      if (direction==0)
	setsizeerr("Taylor of floor with unsigned limit");
      if (direction==-1)
	l=l-1;
    }
    return is_zero(l)?vecteur(0):makevecteur(l);
  }
  const string _floor_s("floor");
  unary_function_eval __floor(&giac::_floor,&D_at_sign,&taylor_floor,_floor_s);
  unary_function_ptr at_floor (&__floor,0,true);

  symbolic symb_ceil(const gen & a){
    return symbolic(at_ceil,a);
  }
  gen _ceil(const gen & args,GIAC_CONTEXT){
    if (is_inf(args)||is_undef(args))
      return args;
    if (args.type==_VECT)
      return apply(args,contextptr,_ceil);
    if (args.type==_CPLX)
      return _ceil(*args._CPLXptr,contextptr)+cst_i*_ceil(*(args._CPLXptr+1),contextptr);
    if ( (args.type==_INT_) || (args.type==_ZINT))
      return args;
    if (args.type==_FRAC){
      gen n=args._FRACptr->num,d=args._FRACptr->den;
      if ( ((n.type==_INT_) || (n.type==_ZINT)) && ( (d.type==_INT_) || (d.type==_ZINT)) )
	return iquo(n,d)+1;
    }
    vecteur l(lidnt(args));
    vecteur lnew=*evalf(l,1,contextptr)._VECTptr;
    gen tmp=subst(args,l,lnew,false,contextptr);
    if (tmp.type==_REAL)
      return -_floor(-tmp,contextptr);
    if (tmp.type!=_DOUBLE_)
      return symb_ceil(args);
    return double2gen(giac_ceil(tmp._DOUBLE_val));
  }
  gen taylor_ceil (const gen & lim_point,const int ordre,const unary_function_ptr & f, int direction,gen & shift_coeff,GIAC_CONTEXT){
    if (ordre<0)
      return 0; // no symbolic preprocessing
    shift_coeff=0;
    gen l=_ceil(lim_point,contextptr);
    if (l==lim_point){
      if (direction==0)
	setsizeerr("Taylor of ceil with unsigned limit");
      if (direction==1)
	l=l-1;
    }
    return is_zero(l)?vecteur(0):makevecteur(l);
  }
  const string _ceil_s("ceil");
  unary_function_eval __ceil(&giac::_ceil,&D_at_sign,&taylor_ceil,_ceil_s);
  unary_function_ptr at_ceil (&__ceil,0,true);

  gen ceiltofloor(const gen & g,GIAC_CONTEXT){
    return -symbolic(at_floor,-g);
  }
  vector< unary_function_ptr > ceil_v(1,at_ceil);
  vector< gen_op_context > ceil2floor_v(1,ceiltofloor);
  gen ceil2floor(const gen & g,GIAC_CONTEXT){
    return subst(g,ceil_v,ceil2floor_v,false,contextptr);
  }


  symbolic symb_round(const gen & a){
    return symbolic(at_round,a);
  }
  gen _round(const gen & args,GIAC_CONTEXT){
    if (is_inf(args)||is_undef(args))
      return args;
    if (args.type==_VECT && args._VECTptr->size()!=2)
      return apply(args,contextptr,_round);
    if (args.type==_VECT && args.subtype==_SEQ__VECT && args._VECTptr->back().type==_INT_){
#ifdef _SOFTMATH_H
      double d=std::giac_gnuwince_pow(10.0,double(args._VECTptr->back().val));
#else
      double d=std::pow(10.0,double(args._VECTptr->back().val));
#endif
      return _round(d*args._VECTptr->front(),contextptr)/d;
    }
    if (args.type==_CPLX)
      return _round(*args._CPLXptr,contextptr)+cst_i*_round(*(args._CPLXptr+1),contextptr);
    return _floor(args+plus_one_half,contextptr);
  }
  gen taylor_round (const gen & lim_point,const int ordre,const unary_function_ptr & f, int direction,gen & shift_coeff,GIAC_CONTEXT){
    if (ordre<0)
      return 0; // no symbolic preprocessing
    shift_coeff=0;
    gen l=_round(lim_point,contextptr);
    if (is_zero(ratnormal(l-lim_point-plus_one_half))){
      if (direction==0)
	setsizeerr("Taylor of round with unsigned limit");
      if (direction==-1)
	l=l-1;
    }
    return is_zero(l)?vecteur(0):makevecteur(l);
  }
  const string _round_s("round");
  unary_function_eval __round(&giac::_round,&D_at_sign,&taylor_round,_round_s);
  unary_function_ptr at_round (&__round,0,true);

  string printasprint(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if (xcas_mode(contextptr)!=3)
      return "print("+feuille.print(contextptr)+")";
    else
      return "Disp "+feuille.print(contextptr);
  }
  symbolic symb_print(const gen & a){
    return symbolic(at_print,a);
  }
  gen _print(const gen & args,GIAC_CONTEXT){
    if ( debug_infolevel && (args.type==_IDNT) && (!args._IDNTptr->localvalue->empty()))
      *logptr(contextptr) << "Local var protected " << (*args._IDNTptr->localvalue)[args._IDNTptr->localvalue->size()-2].val << endl;
    gen tmp=args.eval(eval_level(contextptr),contextptr);
    // If giac used inside a console don't add to messages, since we print
    if (!child_id){
      if (args.type==_IDNT)
        messages_to_print += args.print(contextptr) + ":"; 
      messages_to_print += tmp.print(contextptr) +'\n';
      // *logptr(contextptr) << "Child " << messages_to_print << endl;
    }
    if (tmp.type==_VECT && !tmp._VECTptr->empty() && tmp._VECTptr->front()==gen("Unquoted",contextptr)){
      vecteur & v=*tmp._VECTptr;
      int s=v.size();
      for (int i=1;i<s;++i)
	*logptr(contextptr) << (v[i].type==_STRNG?v[i]._STRNGptr->c_str():unquote(v[i].print(contextptr)));
    }
    else {
      if (args.type==_IDNT)
	*logptr(contextptr) << args << ":";
      *logptr(contextptr) << tmp << endl;
    }
    return __interactive.op(symbolic(at_print,tmp),contextptr);
  }
  const string _print_s("print");
  unary_function_eval __print(&giac::_print,_print_s,&printasprint);
  unary_function_ptr at_print (&__print,_QUOTE_ARGUMENTS,true);

  symbolic symb_is_prime(const gen & a){
    return symbolic(at_is_prime,a);
  }
  gen _is_prime(const gen & args){
    if (!is_integer(args))
      settypeerr();
#ifdef HAVE_LIBPARI
    string s=pari_isprime(args);
    if (s.size()==1)
      return gen(s,context0);
    return string2gen(s,false);
#else
    return is_probab_prime_p(args);
#endif
  }
  const string _is_prime_s("is_prime");
  unary_function_unary __is_prime(&giac::_is_prime,_is_prime_s);
  unary_function_ptr at_is_prime (&__is_prime,0,true);

  gen _is_pseudoprime(const gen & args){
    return is_probab_prime_p(args);
  }
  const string _is_pseudoprime_s("is_pseudoprime");
  unary_function_unary __is_pseudoprime(&giac::_is_pseudoprime,_is_pseudoprime_s);
  unary_function_ptr at_is_pseudoprime (&__is_pseudoprime,0,true);

  gen nextprime1(const gen & a){
    return nextprime(a+1);
  }
  const string _nextprime_s("nextprime");
  unary_function_unary __nextprime(&giac::nextprime1,_nextprime_s);
  unary_function_ptr at_nextprime (&__nextprime,0,true);

  gen prevprime1(const gen & a){
    return prevprime(a-1);
  }
  const string _prevprime_s("prevprime");
  unary_function_unary __prevprime(&giac::prevprime1,_prevprime_s);
  unary_function_ptr at_prevprime (&__prevprime,0,true);

  symbolic symb_jacobi_symbol(const gen & a,const gen & b){
    return symbolic(at_jacobi_symbol,makevecteur(a,b));
  }
  gen _jacobi_symbol(const gen & args){
    check_2d_vecteur(args);
    return jacobi(args._VECTptr->front(),args._VECTptr->back());    
  }
  const string _jacobi_symbol_s("jacobi_symbol");
  unary_function_unary __jacobi_symbol(&giac::_jacobi_symbol,_jacobi_symbol_s);
  unary_function_ptr at_jacobi_symbol (&__jacobi_symbol,0,true);

  symbolic symb_legendre_symbol(const gen & a,const gen & b){
    return symbolic(at_legendre_symbol,makevecteur(a,b));
  }
  gen _legendre_symbol(const gen & args){
    check_2d_vecteur(args);
    return legendre(args._VECTptr->front(),args._VECTptr->back());    
  }
  const string _legendre_symbol_s("legendre_symbol");
  unary_function_unary __legendre_symbol(&giac::_legendre_symbol,_legendre_symbol_s);
  unary_function_ptr at_legendre_symbol (&__legendre_symbol,0,true);

  symbolic symb_ichinrem(const gen & a,const gen & b){
    return symbolic(at_ichinrem,makevecteur(a,b));
  }

  gen ichinrem2(const gen  & a_orig,const gen & b_orig){
    gen a=a_orig;
    gen b=b_orig;
    if (a.type==_MOD)
      a=makevecteur(*a._MODptr,*(a._MODptr+1));
    if (b.type==_MOD)
      b=makevecteur(*b._MODptr,*(b._MODptr+1));
    vecteur l(lvar(a)); lvar(b,l);
    if (l.empty()){
      check_2d_vecteur(a);
      check_2d_vecteur(b);
      vecteur & av=*a._VECTptr;
      vecteur & bv=*b._VECTptr;
      gen &ab=av.back();
      gen &bb=bv.back();
      gen &aa=av.front();
      gen &ba=bv.front();
      gen res=ichinrem(aa,ba,ab,bb);
      if (a_orig.type==_MOD)
	return makemod(res,lcm(ab,bb));
      return makevecteur(res,lcm(ab,bb));
    }
    gen x=l.front();
    if (a.type!=_VECT || b.type!=_VECT ){
      // a and b are polynomial, must have the same degrees
      // build a new polynomial calling ichinrem2 on each element
      gen ax=_e2r(makevecteur(a_orig,x),context0),bx=_e2r(makevecteur(b_orig,x),context0); // ok
      if (ax.type!=_VECT || bx.type!=_VECT )
	setsizeerr();
      int as=ax._VECTptr->size(),bs=bx._VECTptr->size();
      if (!as || !bs)
	setsizeerr("Null polynomial");
      gen a0=ax._VECTptr->front(),b0=bx._VECTptr->front(),m,n;
      if (a0.type==_MOD)
	m=*(a0._MODptr+1);
      else
	setsizeerr("Expecting modular coeff");
      if (b0.type==_MOD)
	n=*(b0._MODptr+1);
      else
	setsizeerr("Expecting modular coeff");
      gen mn=lcm(m,n);
      const_iterateur it=ax._VECTptr->begin(),itend=ax._VECTptr->end(),jt=bx._VECTptr->begin();
      vecteur res;
      for (;as>bs;--as,++it){
	res.push_back(makemod(unmod(*it),mn));
      }
      for (;bs>as;--bs,++jt){
	res.push_back(makemod(unmod(*jt),mn));
      }
      for (;it!=itend;++it,++jt)
	res.push_back(ichinrem2(makemod(unmod(*it),m),makemod(unmod(*jt),n)));
      return _r2e(makevecteur(res,x),context0); // ok
    }
    if (a.type==_VECT && a._VECTptr->size()==2 && b.type==_VECT && b._VECTptr->size()==2 ){
      // ax and bx are the polynomials, 
      gen ax=_e2r(makevecteur(a._VECTptr->front(),x),context0),bx=_e2r(makevecteur(b._VECTptr->front(),x),context0); // ok
      if (ax.type!=_VECT || bx.type!=_VECT )
	setsizeerr();
      gen m=a._VECTptr->back(),n=b._VECTptr->back(),mn=lcm(m,n);
      int as=ax._VECTptr->size(),bs=bx._VECTptr->size();
      const_iterateur it=ax._VECTptr->begin(),itend=ax._VECTptr->end(),jt=bx._VECTptr->begin();
      vecteur res;
      for (;as>bs;--as,++it){
	res.push_back(*it);
      }
      for (;bs>as;--bs,++jt){
	res.push_back(*jt);
      }
      for (;it!=itend;++it,++jt){
	gen tmp=ichinrem2(makevecteur(*it,m),makevecteur(*jt,n));
	if (tmp.type!=_VECT)
	  setsizeerr();
	res.push_back(tmp._VECTptr->front());
      }
      if (a_orig.type==_MOD)
	return makemod(_r2e(makevecteur(res,x),context0),mn); // ok
      return makevecteur(_r2e(makevecteur(res,x),context0),m*n); // ok
    }
    setsizeerr();
    return 0;
  }
  gen _ichinrem(const gen & args){
    if (args.type!=_VECT)
      settypeerr("[a % p, b % q,...]");
    vecteur & v = *args._VECTptr;
    int s=v.size();
    if (s<2)
      setdimerr();
    gen res=ichinrem2(v[0],v[1]);
    for (int i=2;i<s;++i)
      res=ichinrem2(res,v[i]);
    return res;
  }
  const string _ichinrem_s("ichinrem");
  unary_function_unary __ichinrem(&giac::_ichinrem,_ichinrem_s);
  unary_function_ptr at_ichinrem (&__ichinrem,0,true);
  
  gen _fracmod(const gen & args){
    if ( (args.type!=_VECT) || (args._VECTptr->size()!=2))
      return symbolic(at_fracmod,args);
    vecteur & v=*args._VECTptr;
    return fracmod(v[0],v[1]);
  }
  const string _fracmod_s("fracmod");
  unary_function_unary __fracmod(&giac::_fracmod,_fracmod_s);
  unary_function_ptr at_fracmod (&__fracmod,0,true);
  
  const string _iratrecon_s("iratrecon"); // maple name, fracmod takes only 2 arg
  unary_function_unary __iratrecon(&giac::_fracmod,_iratrecon_s);
  unary_function_ptr at_iratrecon (&__iratrecon,0,true);
  
  symbolic symb_chinrem(const gen & a,const gen & b){
    return symbolic(at_chinrem,makevecteur(a,b));
  }
  vecteur polyvect(const gen & a,const vecteur & v){
    if (a.type==_POLY)
      return polynome2poly1(*a._POLYptr,1);
    return vecteur(1,a);
  }
  gen _chinrem(const gen & args,GIAC_CONTEXT){
    if ( (args.type!=_VECT) || (args._VECTptr->size()<2) )
      setsizeerr();
    gen a=args._VECTptr->front();
    gen b=(*args._VECTptr)[1];
    check_2d_vecteur(a);
    check_2d_vecteur(b);
    if ((a._VECTptr->front().type!=_VECT) || (a._VECTptr->back().type!=_VECT) || (b._VECTptr->front().type!=_VECT) || (b._VECTptr->back().type!=_VECT) ){
      vecteur lv;
      if (args._VECTptr->size()==3)
	lv=vecteur(1,(*args._VECTptr)[2]);
      else
	lv=vecteur(1,vx_var);
      lvar(args,lv);
      vecteur lvprime(lv.begin()+1,lv.end());
      gen aa=e2r(a,lv,contextptr),bb=e2r(b,lv,contextptr),aan,aad,bbn,bbd;
      fxnd(aa,aan,aad);
      if (aad.type==_POLY){
	if (aad._POLYptr->lexsorted_degree() )
	  setsizeerr();
	else
	  aad=aad._POLYptr->trunc1();
      }
      fxnd(bb,bbn,bbd);
      if (bbd.type==_POLY){
	if (bbd._POLYptr->lexsorted_degree() )
	  setsizeerr();
	else
	  bbd=bbd._POLYptr->trunc1();
      }
      vecteur & aanv=*aan._VECTptr;
      vecteur & bbnv=*bbn._VECTptr;
      aanv[0]=polyvect(aanv[0],lv)/aad;
      aanv[1]=polyvect(aanv[1],lv);
      bbnv[0]=polyvect(bbnv[0],lv)/bbd;
      bbnv[1]=polyvect(bbnv[1],lv);
      vecteur res=*_chinrem(makevecteur(aanv,bbnv),contextptr)._VECTptr;
      // convert back
      res[0]=symb_horner(*r2e(res[0],lvprime,contextptr)._VECTptr,lv.front());
      res[1]=symb_horner(*r2e(res[1],lvprime,contextptr)._VECTptr,lv.front());
      return res;
    }
    modpoly produit=(*a._VECTptr->back()._VECTptr)**b._VECTptr->back()._VECTptr;
    return makevecteur(gen(chinrem(*a._VECTptr->front()._VECTptr,*b._VECTptr->front()._VECTptr,*a._VECTptr->back()._VECTptr,*b._VECTptr->back()._VECTptr,0),_POLY1__VECT),gen(produit,_POLY1__VECT));    
  }
  const string _chinrem_s("chinrem");
  unary_function_eval __chinrem(&giac::_chinrem,_chinrem_s);
  unary_function_ptr at_chinrem (&__chinrem,0,true);

  gen _factorial(const gen & args){
    if (args.type==_VECT)
      return apply(args,_factorial);
    if (args.type!=_INT_)
      return symbolic(at_factorial,args);
    if (args.val<0)
      return plus_inf;
    return factorial((unsigned long int) args.val);
  }
  const string _factorial_s("factorial");
  unary_function_unary __factorial(&giac::_factorial,_factorial_s);
  unary_function_ptr at_factorial (&__factorial,0,true);

  gen double_is_int(const gen & g,GIAC_CONTEXT){
    gen f=_floor(g,contextptr);
    gen f1=evalf(g-f,1,contextptr);
    if (f1.type==_DOUBLE_ && fabs(f1._DOUBLE_val)<epsilon(contextptr))
      return f;
    else
      return g;
  }
  gen comb(const gen & a_orig,const gen &b_orig,GIAC_CONTEXT){
    gen a=double_is_int(a_orig,contextptr);
    gen b=double_is_int(b_orig,contextptr);
    if (a.type!=_INT_ || b.type!=_INT_)
      return Gamma(a+1,contextptr)/Gamma(b+1,contextptr)/Gamma(a-b+1,contextptr);
    return comb((unsigned long int) a.val,(unsigned long int) b.val);
  }
  gen _comb(const gen & args,GIAC_CONTEXT){
    if (ckmatrix(args))
      return apply(args._VECTptr->front(),args._VECTptr->back(),contextptr,comb);
    if ( (args.type!=_VECT) || (args._VECTptr->size()!=2))
      return symbolic(at_comb,args);
    vecteur & v=*args._VECTptr;
    if (v.front().type!=_INT_ || v.back().type!=_INT_)
      return comb(v.front(),v.back(),contextptr); 
    if (v.front().val<v.back().val)
      return zero;
    if (v.front().val<0)
      return undef;
    return comb((unsigned long int) v.front().val,(unsigned long int) v.back().val);
  }
  const string _comb_s("comb");
  unary_function_eval __comb(&giac::_comb,_comb_s);
  unary_function_ptr at_comb (&__comb,0,true);

  gen perm(const gen & a,const gen &b){
    if (a.type!=_INT_ || b.type!=_INT_)
      return symbolic(at_perm,makevecteur(a,b));
    return perm((unsigned long int) a.val,(unsigned long int) b.val);
  }
  gen _perm(const gen & args){
    if (ckmatrix(args))
      return apply(args._VECTptr->front(),args._VECTptr->back(),perm);
    if ( (args.type!=_VECT) || (args._VECTptr->size()!=2) || (args._VECTptr->front().type!=_INT_) || (args._VECTptr->back().type!=_INT_) )
      return symbolic(at_perm,args);
    if (args._VECTptr->front().val<args._VECTptr->back().val)
      return zero;
    if (args._VECTptr->front().val<0)
      return undef;
    return perm((unsigned long int) args._VECTptr->front().val,(unsigned long int) args._VECTptr->back().val);
  }
  const string _perm_s("perm");
  unary_function_unary __perm(&giac::_perm,_perm_s);
  unary_function_ptr at_perm (&__perm,0,true);

  // ******************
  // Matrix functions
  // *****************

  symbolic symb_tran(const gen & a){
    return symbolic(at_tran,a);
  }
  symbolic symb_trace(const gen & a){
    return symbolic(at_trace,a);
  }
  symbolic symb_rref(const gen & a){
    return symbolic(at_rref,a);
  }
  symbolic symb_idn(const gen & e){
    return symbolic(at_idn,e);
  }
  symbolic symb_ranm(const gen & e){
    return symbolic(at_ranm,e);
  }
  symbolic symb_det(const gen & a){
    return symbolic(at_det,a);
  }
  symbolic symb_pcar(const gen & a){
    return symbolic(at_pcar,a);
  }
  symbolic symb_ker(const gen & a){
    return symbolic(at_ker,a);
  }  
  symbolic symb_image(const gen & a){
    return symbolic(at_image,a);
  }
  symbolic symb_moyal(const gen & a,const gen & b, const gen &vars,const gen & order){
    return symbolic(at_moyal,makevecteur(a,b,vars,order));
  }

  gen _evalf(const gen & a,GIAC_CONTEXT){
    if (a.is_symb_of_sommet(at_equal)&&a._SYMBptr->feuille.type==_VECT && a._SYMBptr->feuille._VECTptr->size()==2){
      vecteur & v(*a._SYMBptr->feuille._VECTptr);
      return symbolic(at_equal,makevecteur(evalf(v.front(),1,contextptr),evalf(v.back(),1,contextptr)));
    }
    if (a.type==_VECT && a.subtype==_SEQ__VECT && a._VECTptr->size()==2 && a._VECTptr->back().type==_INT_){
      int save_decimal_digits=decimal_digits(contextptr);
      int ndigits=a._VECTptr->back().val;
      set_decimal_digits(ndigits,contextptr);
      gen res=a._VECTptr->front().evalf(1,contextptr);
      if (res.type==_REAL || res.type==_CPLX)
	res=accurate_evalf(res,digits2bits(a._VECTptr->back().val));
      set_decimal_digits(save_decimal_digits,contextptr);
      if (ndigits<14)
	res=_round(gen(makevecteur(res,ndigits),_SEQ__VECT),contextptr);
      return res;
    }
    return a.evalf(1,contextptr);
  }
  const string _evalf_s("evalf");
  unary_function_eval __evalf(&giac::_evalf,_evalf_s);
  unary_function_ptr at_evalf (&__evalf,0,true);
  symbolic symb_evalf(const gen & a){
    return symbolic(at_evalf,a);
  }

  gen _eval(const gen & a,GIAC_CONTEXT){
    if (a.is_symb_of_sommet(at_equal)&&a._SYMBptr->feuille.type==_VECT && a._SYMBptr->feuille._VECTptr->size()==2){
      vecteur & v(*a._SYMBptr->feuille._VECTptr);
      return symbolic(at_equal,makevecteur(eval(v.front(),eval_level(contextptr),contextptr),eval(v.back(),eval_level(contextptr),contextptr)));
    }
    if (a.type==_VECT && a.subtype==_SEQ__VECT && a._VECTptr->size()==2){
      gen a1=a._VECTptr->front(),a2=a._VECTptr->back();
      if (a2.type==_INT_)
	return a1.eval(a2.val,contextptr);
      return _subst(gen(makevecteur(eval(a1,eval_level(contextptr),contextptr),a2),_SEQ__VECT),contextptr);
    }
    return a.eval(1,contextptr).eval(eval_level(contextptr),contextptr);
  }
  const string _eval_s("eval");
  unary_function_eval __eval(&giac::_eval,_eval_s);
  unary_function_ptr at_eval (&__eval,_QUOTE_ARGUMENTS,true);
  symbolic symb_eval(const gen & a){
    return symbolic(at_eval,a);
  }
  
  const string _evalm_s("evalm");
  unary_function_eval __evalm(&giac::_eval,_evalm_s);
  unary_function_ptr at_evalm (&__evalm,0,true);
  
  gen _ampersand_times(const gen & g){
    if (g.type!=_VECT || g._VECTptr->size()!=2)
      setsizeerr();
    return g._VECTptr->front()*g._VECTptr->back();
  }
  const string _ampersand_times_s("&*");
  unary_function_unary __ampersand_times(&giac::_ampersand_times,_ampersand_times_s,&printsommetasoperator,&texprintsommetasoperator);
  unary_function_ptr at_ampersand_times (&__ampersand_times);
  
  const string _subst_s("subst");
  gen _subst(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      settypeerr();
    vecteur & v = *args._VECTptr;
    int s=v.size();
    if (s==2){
      gen e=v.back();
      if (e.type==_VECT){
	vecteur & w = *e._VECTptr;
	vecteur vin,vout;
	const_iterateur it=w.begin(),itend=w.end();
	for (;it!=itend;++it){
	  if (it->type!=_SYMB)
	    continue;
	  if (it->_SYMBptr->sommet!=at_equal && it->_SYMBptr->sommet!=at_same)
	    continue;
	  vin.push_back(it->_SYMBptr->feuille._VECTptr->front());
	  vout.push_back(it->_SYMBptr->feuille._VECTptr->back());
	}
	gen res=subst(v.front(),vin,vout,false,contextptr);
	return res;
      }
      if (e.type!=_SYMB)
	settypeerr();
      if (e._SYMBptr->sommet!=at_equal && e._SYMBptr->sommet!=at_same)
	setsizeerr();
      return subst(v.front(),e._SYMBptr->feuille._VECTptr->front(),e._SYMBptr->feuille._VECTptr->back(),false,contextptr);
    }
    if (s<3)
      toofewargs(_subst_s);
    if (s>3)
      toomanyargs(_subst_s);
    if (v[1].is_symb_of_sommet(at_equal))
      return _subst(makevecteur(v.front(),vecteur(v.begin()+1,v.end())),contextptr);
    return subst(v.front(),v[1],v.back(),false,contextptr);
  }
  unary_function_eval __subst(&giac::_subst,_subst_s);
  unary_function_ptr at_subst (&__subst,0,true);
  symbolic symb_subst(const gen & a){
    return symbolic(at_subst,a);
  }

  string printassubs(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if (xcas_mode(contextptr)!=1 || feuille.type!=_VECT || feuille._VECTptr->size()!=2)
      return sommetstr+"("+feuille.print(contextptr)+")";
    vecteur & v=*feuille._VECTptr;
    vecteur w=mergevecteur(vecteur(1,v.back()),vecteur(v.begin(),v.end()-1));
    return sommetstr+"("+gen(w,_SEQ__VECT).print(contextptr)+")";
  }  
  gen _subs(const gen & g,GIAC_CONTEXT){
    return _subst(g,contextptr);
  }
  const string _subs_s("subs");
  unary_function_eval __subs(&_subs,_subs_s,&printassubs);
  unary_function_ptr at_subs (&__subs);

  string printasmaple_subs(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if (xcas_mode(contextptr)==1 || feuille.type!=_VECT || feuille._VECTptr->size()<2)
      return sommetstr+"("+feuille.print(contextptr)+")";
    vecteur & v=*feuille._VECTptr;
    vecteur w=mergevecteur(vecteur(1,v.back()),vecteur(v.begin(),v.end()-1));
    return sommetstr+"("+gen(w,_SEQ__VECT).print(contextptr)+")";
  }  
  gen _maple_subs(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT || g._VECTptr->size()<2)
      return _subst(g,contextptr);
    vecteur & v=*g._VECTptr;
    return _subst(makevecteur(v.back(),vecteur(v.begin(),v.end()-1)),contextptr);
  }
  const string _maple_subs_s("subs");
  unary_function_eval __maple_subs(&_maple_subs,_maple_subs_s,&printasmaple_subs);
  unary_function_ptr at_maple_subs (&__maple_subs);


  string version(){
    return string("giac ")+VERSION;
  }
  gen _version(const gen & a){
    return string2gen(version(),false);
  }
  const string _version_s("version");
  unary_function_unary __version(&giac::_version,_version_s);
  unary_function_ptr at_version (&__version,0,true);

  void prod2frac(const gen & g,vecteur & num,vecteur & den){
    num.clear();
    den.clear();
    if (g.type==_FRAC){
      vecteur num2,den2;
      prod2frac(g._FRACptr->num,num,den);
      prod2frac(g._FRACptr->den,den2,num2);
      num=mergevecteur(num,num2);
      den=mergevecteur(den,den2);
      return;      
    }
    if (g.is_symb_of_sommet(at_neg)){
      prod2frac(g._SYMBptr->feuille,num,den);
      if (!num.empty()){
	num.front()=-num.front();
	return;
      }
    }
    if ( (g.type!=_SYMB) || (g._SYMBptr->sommet!=at_prod) || (g._SYMBptr->feuille.type!=_VECT)){
      if (g.is_symb_of_sommet(at_division)){
	vecteur num2,den2;
	prod2frac(g._SYMBptr->feuille._VECTptr->front(),num,den);
	prod2frac(g._SYMBptr->feuille._VECTptr->back(),den2,num2);
	num=mergevecteur(num,num2);
	den=mergevecteur(den,den2);
	return;
      }
      if (g.is_symb_of_sommet(at_inv))
	prod2frac(g._SYMBptr->feuille,den,num);
      else
	num=vecteur(1,g);
      return;
    }
    vecteur & v=*g._SYMBptr->feuille._VECTptr;
    const_iterateur it=v.begin(),itend=v.end();
    for (;it!=itend;++it){
      if ( (it->type==_SYMB) && (it->_SYMBptr->sommet==at_inv) )
	den.push_back(it->_SYMBptr->feuille);
      else
	num.push_back(*it);
    }
  }

  gen vecteur2prod(const vecteur & num){
    if (num.empty())
      return plus_one;
    if (num.size()==1)
      return num.front();
    return symbolic(at_prod,gen(num,_SEQ__VECT));
  }

  bool need_parenthesis(const gen & g){
    if (g.type==_INT_ || g.type==_ZINT)
      return is_strictly_positive(-g,context0);  // ok
    if (g.type==_CPLX){
      gen rg=re(-g,context0),ig(im(-g,context0)); // ok
      if ( is_zero(rg))
	return is_strictly_positive(ig,context0); // ok
      if (is_zero(ig) )
	return is_strictly_positive(rg,context0); // ok
      return true;
    }
    if (g.type==_SYMB)
      return need_parenthesis(g._SYMBptr->sommet);
    if (g.type!=_FUNC)
      return false;
    unary_function_ptr & u=*g._FUNCptr;
    if (u==at_pow)
      return false;
    if (u==at_neg || u==at_minus || u==at_and || u==at_ou || u==at_xor || u==at_same || u==at_equal)
      return true;
    if (!u.ptr->printsommet)
      return false;
    return true;
  }

  gen _multistring(const gen & args,GIAC_CONTEXT){
    string res;
    if (args.type==_VECT){
      const_iterateur it=args._VECTptr->begin(),itend=args._VECTptr->end();
      for (;it!=itend;){
	if (it->type!=_STRNG)
	  break;
	res += *it->_STRNGptr;
	++it;
	if (it==itend)
	  return string2gen(res,false);;
	res += '\n';
      }
    }
    else {// newline added, otherwise Eqw_compute_size would fail
      if (args.type==_STRNG)
	res=*args._STRNGptr;
      else
	res=args.print(contextptr);
      res += '\n'; 
    }
    return string2gen(res,false);
  }
  const string _multistring_s("multistring");
  unary_function_eval __multistring(&giac::_multistring,_multistring_s);
  unary_function_ptr at_multistring (&__multistring);

  // Gamma function
  // lnGamma_minus is ln(Gamma)-(z-1/2)*ln(z)+z which is tractable at +inf
  gen taylor_lnGamma_minus(const gen & lim_point,const int ordre,const unary_function_ptr & f, int direction,gen & shift_coeff,GIAC_CONTEXT){
    if (ordre<0)
      return 0;
    if (lim_point!=plus_inf)
      setsizeerr();
    shift_coeff=0;
    vecteur v;
    // ln(Gamma(z)) = (z-1/2)*ln(z) - z +
    //                ln(2*pi)/2 + sum(B_2n /((2n)*(2n-1)*z^(2n-1)),n>=1)
    v.push_back(symbolic(at_ln,cst_two_pi)/2);
    for (int n=1;2*n<=ordre;++n){
      v.push_back(bernoulli(2*n)/(4*n*n-2*n));
      v.push_back(0);
    }
    v.push_back(undef);
    return v;
  }
  // lnGamma_minus is ln(Gamma)-(z-1/2)*ln(z)+z which is tractable at +inf
  gen d_lnGamma_minus(const gen & args,GIAC_CONTEXT){
    return Psi(args,0)+1-symbolic(at_ln,args)-(args+minus_one_half)/args;
  }
  partial_derivative_onearg D_at_lnGamma_minus(&d_lnGamma_minus);
  extern unary_function_ptr at_lnGamma_minus;
  gen _lnGamma_minus(const gen & g){
    if (is_inf(g))
      return symbolic(at_ln,cst_two_pi)/2;
    return symbolic(at_lnGamma_minus,g);
  }
  const string _lnGamma_minus_s("lnGamma_minus");
  unary_function_unary __lnGamma_minus(&_lnGamma_minus,&D_at_lnGamma_minus,&taylor_lnGamma_minus,_lnGamma_minus_s);
  unary_function_ptr at_lnGamma_minus (&__lnGamma_minus,0,true);
  // ln(Gamma) = lnGamma_minus + (z-1/2)*ln(z)-z which is tractable at +inf
  gen Gamma_replace(const gen & g,GIAC_CONTEXT){
    return symbolic(at_exp,(g+minus_one_half)*symbolic(at_ln,g)-g)*symbolic(at_exp,_lnGamma_minus(g));
  }
  gen taylor_Gamma (const gen & lim_point,const int ordre,const unary_function_ptr & f, int direction,gen & shift_coeff,GIAC_CONTEXT){
    if (ordre<0){
      limit_tractable_functions.push_back(&at_Gamma);
      limit_tractable_replace.push_back(Gamma_replace);
      return 1;
    }
    shift_coeff=0;
    if (!is_integer(lim_point) || is_strictly_positive(lim_point,contextptr))
      return taylor(lim_point,ordre,f,0,shift_coeff,contextptr);
    // Laurent series for Gamma
    if (lim_point.type!=_INT_)
      setsizeerr();
    vecteur v;
    identificateur x(" ");
    int n=-lim_point.val;
    gen decal(1);
    for (int i=1;i<=n;++i){
      decal = decal/(x-i);
    }
    taylor(decal,x,zero,ordre,v,contextptr);
    gen Psi1=taylor(1,ordre,f,0,shift_coeff,contextptr);
    shift_coeff=-1;
    if (Psi1.type!=_VECT)
      setsizeerr();
    v=operator_times(v,*Psi1._VECTptr,0);
    v=vecteur(v.begin(),v.begin()+ordre);
    v.push_back(undef);
    return v;
  }
  gen d_Gamma(const gen & args,GIAC_CONTEXT){
    return Psi(args,0)*Gamma(args,contextptr);
  }
  partial_derivative_onearg D_at_Gamma(&d_Gamma);
  gen Gamma(const gen & x,GIAC_CONTEXT){
    if (x.type==_INT_){
      if (x.val<=0)
	return plus_inf;
      return factorial(x.val-1);
    }
    if (x.type==_FRAC && x._FRACptr->den==2 && x._FRACptr->num.type==_INT_){
      int n=x._FRACptr->num.val;
      // compute Gamma(n/2)
      gen factnum=1,factden=1;
      for (;n>1;n-=2){
	factnum=(n-2)*factnum;
	factden=2*factden;
      }
      for (;n<1;n+=2){
	factnum=2*factnum;
	factden=n*factden;
      }
      return factnum/factden*sqrt(cst_pi,contextptr);
    }
#ifdef HAVE_LIBGSL
    if (x.type==_DOUBLE_)
      return gsl_sf_gamma(x._DOUBLE_val);
#endif
#ifdef HAVE_LIBMPFR
    if (x.type==_REAL && is_positive(x,contextptr)){
      mpfr_t gam;
      int prec=mpfr_get_prec(x._REALptr->inf);
      mpfr_init2(gam,prec);
      mpfr_gamma(gam,x._REALptr->inf,GMP_RNDN);
      real_object res(gam);
      mpfr_clear(gam);
      return res;
    }
#endif
#ifdef HAVE_LIBPARI
    if (x.type==_CPLX)
      return pari_gamma(x);
#endif
    if (x.type==_DOUBLE_ || x.type==_CPLX){
      if (is_positive(.5-re(x,contextptr),contextptr))
	return cst_pi / (sin(M_PI*x,contextptr)*Gamma(1-x,contextptr));
      static double p[] = {
	0.99999999999980993, 676.5203681218851, -1259.1392167224028,
	771.32342877765313, -176.61502916214059, 12.507343278686905,
	-0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
      gen z = x-1;
      gen X = p[0];
      int g=7;
      for (int i=1;i<g+2;++i)
	X += gen(p[i])/(z+i);
      gen t = z + g + 0.5;
      return sqrt(2*cst_pi,contextptr) * pow(t,z+0.5,contextptr) * exp(-t,contextptr) * X;      
    }
    return symbolic(at_Gamma,x);
  }
  gen _Gamma(const gen & args,GIAC_CONTEXT) {
    return Gamma(args,contextptr);
  }
  const string _Gamma_s("Gamma");
  unary_function_eval __Gamma(&_Gamma,&D_at_Gamma,&taylor_Gamma,_Gamma_s);
  unary_function_ptr at_Gamma (&__Gamma,0,true);

  // diGamma function
  gen taylor_Psi_minus_ln(const gen & lim_point,const int ordre,const unary_function_ptr & f, int direction,gen & shift_coeff,GIAC_CONTEXT){
    if (ordre<0)
      return 0;
    if (lim_point!=plus_inf)
      setsizeerr();
    shift_coeff=1;
    vecteur v(1,minus_one_half);
    // Psi(z)=ln(z)-1/(2*z)-sum(B_2n /(2*n*z^(2n)),n>=1)
    for (int n=2;n<=ordre;n+=2){
      v.push_back(-bernoulli(n)/(2*n));
      v.push_back(0);
    }
    v.push_back(undef);
    return v;
  }
  gen d_Psi_minus_ln(const gen & args,GIAC_CONTEXT){
    return inv(args,contextptr)-Psi(args,1,contextptr);
  }
  partial_derivative_onearg D_at_Psi_minus_ln(&d_Psi_minus_ln);
  extern unary_function_ptr at_Psi_minus_ln;
  gen _Psi_minus_ln(const gen & g){
    if (is_inf(g))
      return 0;
    return symbolic(at_Psi_minus_ln,g);
  }
  const string _Psi_minus_ln_s("Psi_minus_ln");
  unary_function_unary __Psi_minus_ln(&_Psi_minus_ln,&D_at_Psi_minus_ln,&taylor_Psi_minus_ln,_Psi_minus_ln_s);
  unary_function_ptr at_Psi_minus_ln (&__Psi_minus_ln,0,true);
  gen Psi_replace(const gen & g,GIAC_CONTEXT){
    return symbolic(at_ln,g)+_Psi_minus_ln(g);
  }
  gen taylor_Psi (const gen & lim_point,const int ordre,const unary_function_ptr & f, int direction,gen & shift_coeff,GIAC_CONTEXT){
    if (ordre<0){
      limit_tractable_functions.push_back(&at_Psi);
      limit_tractable_replace.push_back(Psi_replace);
      return 1;
    }
    shift_coeff=0;
    if (!is_integer(lim_point) || is_strictly_positive(lim_point,contextptr))
      return taylor(lim_point,ordre,f,0,shift_coeff,contextptr);
    // FIXME Laurent series for Psi
    if (lim_point.type!=_INT_)
      setsizeerr();
    vecteur v;
    identificateur x(" ");
    int n=-lim_point.val;
    gen decal;
    for (int i=0;i<n;++i){
      decal -= inv(x+i,contextptr);
    }
    taylor(decal,x,lim_point,ordre,v,contextptr);
    gen Psi1=taylor(1,ordre,f,0,shift_coeff,contextptr);
    shift_coeff=-1;
    if (Psi1.type!=_VECT)
      setsizeerr();
    v=v+*Psi1._VECTptr;
    v.insert(v.begin(),-1);
    return v;
  }
  gen d_Psi(const gen & args,GIAC_CONTEXT){
    vecteur v(gen2vecteur(args));
    if (v.size()==1)
      v.push_back(0);
    if (v.size()!=2 || v.back().type!=_INT_)
      setdimerr();
    return Psi(v.front(),v.back().val+1,contextptr);
  }
  partial_derivative_onearg D_at_Psi(&d_Psi);

  gen Psi(const gen & x,GIAC_CONTEXT){
    if (is_positive(-x,contextptr)){
      if (is_integer(x))
	return unsigned_inf;
      return Psi(ratnormal(1-x),contextptr)-cst_pi/tan(cst_pi*x,contextptr);
    }
    if (x==plus_inf)
      return x;
    if (is_inf(x) || (is_undef(x)) )
      return undef;
    if ( (x.type==_INT_) && (x.val<10000) && (x.val>=1)){
      identificateur tt(" t");
      return -cst_euler_gamma+sum_loop(inv(tt,contextptr),tt,1,x.val-1,contextptr);
    }
    if (x.type==_FRAC){
      // Psi(m/k) for 0<m<k
      // Psi(m/k) = -euler_gamma -ln(2k) - pi/2/tan(m*pi/k) +
      //    + 2 sum( cos(2 *pi*n*m/k)*ln(sin(n*pi/k)), n=1..floor (k-1)/2 )
      gen num=x._FRACptr->num,den=x._FRACptr->den;
      if (num.type==_INT_ && den.type==_INT_ && den.val<13){
	int m=num.val,k=den.val;
	gen res;
	int mk=m/k;
	for (int i=mk;i>0;--i){
	  m -= k;
	  res += inv(m,contextptr);
	}
	res = k*res - cst_euler_gamma - ln(2*k,contextptr) - cst_pi/2/tan(m*cst_pi/k,contextptr);;
	gen res1 ;
	for (int n=1;n<=(k-1)/2;n++){
	  res1 += cos(2*n*m*cst_pi/k,contextptr)*ln(sin(n*cst_pi/k,contextptr),contextptr);
	}
	return res + 2*res1;
      }
    }
#ifdef HAVE_LIBGSL
    if (x.type==_DOUBLE_)
      return gsl_sf_psi(x._DOUBLE_val);
#endif
#ifdef HAVE_LIBPARI
    if (x.type==_CPLX || x.type==_REAL)
      return pari_psi(x);
#endif
    return symbolic(at_Psi,x);
  }
  // n-th derivative of digamma function
  gen Psi(const gen & x,int n,GIAC_CONTEXT){
    if (n<-1)
      setsizeerr();
    if (n==-1)
      return Gamma(x,contextptr);
    if (is_positive(-x,contextptr))
      return unsigned_inf;
    if (is_one(x)){
      if (n%2)
	return Zeta(n+1,contextptr)*factorial(n);
      else
	return -Zeta(n+1,contextptr)*factorial(n);
    }
    if (x==plus_inf)
      return zero;
    if (is_inf(x) || (is_undef(x)) )
      return undef;
    if (!n)
      return Psi(x,contextptr);
    if ( (x.type==_INT_) && (x.val<10000) ){
      identificateur tt(" t");
      if (n%2)
	return factorial(n)*(Zeta(n+1,contextptr)-sum_loop(pow(tt,-n-1),tt,1,x.val-1,contextptr));
      else
	return -factorial(n)*(Zeta(n+1,contextptr)-sum_loop(pow(tt,-n-1),tt,1,x.val-1,contextptr));
    }
#ifdef HAVE_LIBGSL
    if (x.type==_DOUBLE_)
      return gsl_sf_psi_n(n,x._DOUBLE_val);
#endif 
    return symbolic(at_Psi,makevecteur(x,n));
  }
  gen _Psi(const gen & args,GIAC_CONTEXT) {
    if (args.type!=_VECT)
      return Psi(args,contextptr);
    if ( args._VECTptr->size()!=2 )
      return symbolic(at_Psi,args);
    gen x(args._VECTptr->front()),n(args._VECTptr->back());
    if (n.type==_REAL)
      n=n.evalf_double(1,contextptr);
    if (n.type==_DOUBLE_)
      n=int(n._DOUBLE_val);
    if (n.type!=_INT_)
      setsizeerr();
    return Psi(x,n.val,contextptr);
  }
  const string _Psi_s("Psi");
  unary_function_eval __Psi(&_Psi,&D_at_Psi,&taylor_Psi,_Psi_s);
  unary_function_ptr at_Psi (&__Psi,0,true);

  gen _normalmod(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT || g._VECTptr->size()!=2)
      setsizeerr();
    gen & f =g._VECTptr->front();
    gen res=normal(makemodquoted(f,g._VECTptr->back()),contextptr);
    if (f.type==_VECT && res.type==_VECT)
      res.subtype=f.subtype;
    return res;
  }
  const string _normalmod_s("%");
  unary_function_eval __normalmod(&_normalmod,_normalmod_s,&printsommetasoperator,&texprintsommetasoperator);
  unary_function_ptr at_normalmod (&__normalmod,0,true);

  // a=expression, x variable, n=number of terms, 
  // compute an approx value of sum((-1)^k*a(k),k,0,+infinity)
  // using Chebychev polynomials
  gen alternate_series(const gen & a,const gen & x,int n,GIAC_CONTEXT){
    gen d=normal((pow(3+2*sqrt(2,contextptr),n)+pow(3-2*sqrt(2,contextptr),n))/2,contextptr);
    gen p=1;
    gen c=d-p;
    gen S=subst(a,x,0,false,contextptr)*c;
    for (int k=1;k<n;k++) {
      p=p*gen(k+n-1)*gen(k-n-1)/gen(k-inv(2,contextptr))/gen(k);
      c=-p-c;
      S=S+subst(a,x,k,false,contextptr)*c;
    }
    return S/d;
  }

  gen Eta(const gen & s,int ndiff,GIAC_CONTEXT){
    if (s.type==_INT_ && !ndiff){
      if (s==1)
	return symbolic(at_ln,2);
      if (s%2==0)
	return (1-pow(2,1-s,contextptr))*Zeta(s,contextptr);
    }
    if (s.type==_DOUBLE_ || s.type==_REAL || (s.type==_CPLX && s.subtype==_DOUBLE_)){
      gen rx=re(s,contextptr).evalf_double(1,contextptr);
      if (rx._DOUBLE_val<0.5){
	if (ndiff){
	  identificateur id(" ");
	  gen t(id),zeta;
	  zeta=derive((1-pow(2,1-t,contextptr))*pow(2*cst_pi,t,contextptr)/cst_pi*sin(cst_pi*t/2,contextptr)*symbolic(at_Gamma,1-t)*symbolic(at_Zeta,1-t),t,ndiff,contextptr);
	  zeta=subst(zeta,t,s,false,contextptr);
	  return zeta;
	}
	gen zeta1=Eta(1-s,0,contextptr)/(1-pow(2,s,contextptr));
	gen zetas=pow(2,s,contextptr)*pow(cst_pi,s-1,contextptr)*sin(cst_pi*s/2,contextptr)*Gamma(1-s,contextptr)*zeta1;
	return (1-pow(2,1-s,contextptr))*zetas;
      }
      // find n such that 3*(1+2*|y|)*exp(|y|*pi/2)*10^ndigits < (3+sqrt(8))^n
      gen ix=im(s,contextptr).evalf_double(1,contextptr);
      if (ix.type!=_DOUBLE_)
	settypeerr();
      double y=std::abs(ix._DOUBLE_val);
      int ndigits=16; // FIXME? use decimal_digits;
      double n=(std::log10(3*(1+2*y)*std::exp(y*M_PI/2))+ndigits)/std::log10(3.+std::sqrt(8.));
      identificateur idx(" ");
      gen x(idx);
      gen res=alternate_series(inv(pow(idx+1,s,contextptr),contextptr)*pow(-ln(idx+1,contextptr),ndiff,contextptr),idx,int(std::ceil(n)),contextptr);
      return res.evalf(1,contextptr);
    }
    else {
      if (ndiff)
	return symbolic(at_Eta,gen(makevecteur(s,ndiff),_SEQ__VECT));
      else
	return symbolic(at_Eta,s);
    }
  }

  gen Eta(const gen & s0,GIAC_CONTEXT){
    gen s=s0;
    int ndiff=0;
    if (s.type==_VECT){
      if (s._VECTptr->size()!=2)
	setsizeerr();
      gen n=s._VECTptr->back();
      if (n.type==_REAL)
	n=n.evalf_double(1,contextptr);
      if (n.type==_DOUBLE_)
	n=int(n._DOUBLE_val);
      if (n.type!=_INT_)
	settypeerr();
      ndiff=n.val;
      s=s._VECTptr->front();
    }
    return Eta(s,ndiff,contextptr);
  }

  gen Zeta(const gen & x,int ndiff,GIAC_CONTEXT){
    if (!ndiff)
      return Zeta(x,contextptr);
    if (x.type==_DOUBLE_ || x.type==_REAL || (x.type==_CPLX && x.subtype==_DOUBLE_)){
      gen rex=re(x,contextptr).evalf_double(1,contextptr);
      if (rex.type!=_DOUBLE_)
	setsizeerr();
      identificateur id(" ");
      gen t(id),zeta;
      if (rex._DOUBLE_val<0.5){
	// Zeta(x)=2^x*pi^(x-1)*sin(pi*x/2)*Gamma(1-x)*zeta(1-x)
	zeta=derive(pow(2*cst_pi,t,contextptr)/cst_pi*sin(cst_pi*t/2,contextptr)*symbolic(at_Gamma,1-t)*symbolic(at_Zeta,1-t),t,ndiff,contextptr);
	zeta=subst(zeta,t,x,false,contextptr);
      }
      else {
	// Zeta=Eta(x)/(1-2^(1-x))
	zeta=derive(symbolic(at_Eta,t)/(1-pow(2,1-t,contextptr)),t,ndiff,contextptr);
	zeta=subst(zeta,t,x,false,contextptr);
      }
      return zeta;
    }
    return symbolic(at_Zeta,gen(makevecteur(x,ndiff),_SEQ__VECT));
  }
  gen Zeta(const gen & x,GIAC_CONTEXT){
    if (x.type==_VECT){
      if (x._VECTptr->size()!=2)
	setsizeerr();
      gen n=x._VECTptr->back();
      if (n.type==_REAL)
	n=n.evalf_double(1,contextptr);
      if (n.type==_DOUBLE_)
	n=int(n._DOUBLE_val);
      if (n.type!=_INT_)
	settypeerr();
      int ndiff=n.val;
      return Zeta(x._VECTptr->front(),ndiff,contextptr);
    }
    if ( (x.type==_INT_)){
      int n=x.val;
      if (!n)
	return minus_one_half;
      if (n==1)
	return plus_inf;
      if (n<0){
	if (n%2)
	  return -rdiv(bernoulli(1-n),(1-n)) ;
	else
	  return zero;
      }
      if (n%2)
	return symbolic(at_Zeta,x);
      else
	return pow(cst_pi,n)*abs(bernoulli(x),contextptr)*rdiv(pow(plus_two,n-1),factorial(n));
    }
#ifdef HAVE_LIBGSL
    if (x.type==_DOUBLE_)
      return gsl_sf_zeta(x._DOUBLE_val);
#endif // HAVE_LIBGSL
#ifdef HAVE_LIBPARI
    if (x.type==_CPLX)
      return pari_zeta(x);
#endif
#ifdef HAVE_LIBMPFR
    if (x.type==_REAL){
      mpfr_t gam;
      int prec=mpfr_get_prec(x._REALptr->inf);
      mpfr_init2(gam,prec);
      mpfr_zeta(gam,x._REALptr->inf,GMP_RNDN);
      real_object res(gam);
      mpfr_clear(gam);
      return res;
    }
#endif
    if (x.type==_CPLX || x.type==_DOUBLE_ || x.type==_REAL)
      return Eta(x,contextptr)/(1-pow(2,1-x,contextptr));
    return symbolic(at_Zeta,x);
  }
  gen _Zeta(const gen & args,GIAC_CONTEXT) {
    return Zeta(args,contextptr);
  }
  gen d_Zeta(const gen & args,GIAC_CONTEXT){
    vecteur v(gen2vecteur(args));
    if (v.size()==1)
      v.push_back(0);
    if (v.size()!=2 || v.back().type!=_INT_)
      setdimerr();
    return Zeta(v.front(),v.back().val+1,contextptr);
  }
  gen taylor_Zeta(const gen & lim_point,const int ordre,const unary_function_ptr & f, int direction,gen & shift_coeff,GIAC_CONTEXT){
    if (ordre<0)
      return 0; // no symbolic preprocessing
    if (is_one(lim_point)){
      shift_coeff=-1;
      identificateur x(" "); vecteur v,w;
      taylor(1-pow(2,1-x,contextptr),x,1,ordre+1,w,contextptr);
      w.erase(w.begin());
      reverse(w.begin(),w.end());
      if (!w.empty() && is_undef(w.front()))
	w.erase(w.begin());
      gen gw=horner(w,x);
      sparse_poly1 sp=series__SPOL1(symbolic(at_Eta,x+1)/gw,x,0,ordre,0,contextptr); 
      sparse_poly1::const_iterator it=sp.begin(),itend=sp.end();
      for (;it!=itend;++it){
	v.push_back(it->coeff); // assumes all coeffs are non zero...
      }
      return v;
    }
    return taylor(lim_point,ordre,f,direction,shift_coeff,contextptr);
  }
  partial_derivative_onearg D_at_Zeta(&d_Zeta);
  const string _Zeta_s("Zeta");
  unary_function_eval __Zeta(&_Zeta,&D_at_Zeta,&taylor_Zeta,_Zeta_s);
  unary_function_ptr at_Zeta (&__Zeta,0,true);

  gen d_Eta(const gen & args,GIAC_CONTEXT){
    vecteur v(gen2vecteur(args));
    if (v.size()==1)
      v.push_back(0);
    if (v.size()!=2 || v.back().type!=_INT_)
      setdimerr();
    return Eta(v.front(),v.back().val+1,contextptr);
  }
  partial_derivative_onearg D_at_Eta(&d_Eta);
  gen _Eta(const gen & args,GIAC_CONTEXT) {
    return Eta(args,contextptr);
  }
  const string _Eta_s("Eta");
  unary_function_eval __Eta(&_Eta,&D_at_Eta,_Eta_s);
  unary_function_ptr at_Eta (&__Eta,0,true);

  // error function
  gen taylor_erfs(const gen & lim_point,const int ordre,const unary_function_ptr & f, int direction,gen & shift_coeff,GIAC_CONTEXT){
    if (ordre<0)
      return 0;
    if (!is_inf(lim_point))
      setsizeerr();
    shift_coeff=1;
    // erfs(x)=1/sqrt(pi) * 1/x* sum( (2*k)! / (-4)^k / k! * x^(-2k) )
    gen tmp(1);
    vecteur v;
    for (int n=0;n<=ordre;){
      v.push_back(tmp);
      v.push_back(0);
      n +=2 ;
      tmp=gen(n-1)/gen(-2)*tmp;
    }
    v.push_back(undef);
    return multvecteur(inv(sqrt(cst_pi,contextptr),contextptr),v);
  }
  gen _erfs(const gen & g);
  gen d_erfs(const gen & args,GIAC_CONTEXT){
    return 2*args*_erfs(args)-gen(2)/sqrt(cst_pi,contextptr);
  }
  partial_derivative_onearg D_at_erfs(&d_erfs);
  extern unary_function_ptr at_erfs;
  gen _erfs(const gen & g){
    if (is_inf(g))
      return 0;
    return symbolic(at_erfs,g);
  }
  const string _erfs_s("erfs");
  unary_function_unary __erfs(&_erfs,&D_at_erfs,&taylor_erfs,_erfs_s);
  unary_function_ptr at_erfs (&__erfs,0,true);
  gen erf_replace(const gen & g,GIAC_CONTEXT){
    return symbolic(at_sign,g)*(1-symbolic(at_exp,-g*g)*_erfs(symbolic(at_abs,g)));
  }
  gen taylor_erf (const gen & lim_point,const int ordre,const unary_function_ptr & f, int direction,gen & shift_coeff,GIAC_CONTEXT){
    if (ordre<0){
      limit_tractable_functions.push_back(&at_erf);
      limit_tractable_replace.push_back(erf_replace);
      return 1;
    }
    shift_coeff=0;
    return taylor(lim_point,ordre,f,0,shift_coeff,contextptr);
  }
  gen d_erf(const gen & e,GIAC_CONTEXT){
    return 2*exp(-pow(e,2),contextptr)/sqrt(cst_pi,contextptr);
  }
  partial_derivative_onearg D_at_erf(d_erf);

  gen erf0(const gen & x,gen & erfc,GIAC_CONTEXT){
    if (x.type==_DOUBLE_){ 
      double absx=std::abs(x._DOUBLE_val);
      if (absx<=3){
	// numerical computation of int(exp(-t^2),t=0..x) 
	// by series expansion at x=0
	// x*sum( (-1)^n*(x^2)^n/n!/(2*n+1),n=0..inf)
	long double z=x._DOUBLE_val,z2=z*z,res=0,pi=1;
	for (int n=0;;){
	  res += pi/(2*n+1);
	  ++n;
	  pi = -pi*z2/n;
	  if (pi<1e-17 && pi>-1e-17)
	    break;
	}
	erfc=double(1-2/std::sqrt(M_PI)*z*res);
	return 2/std::sqrt(M_PI)*double(z*res);
      }
      if (absx>=6.5){
	// asymptotic expansion at infinity of int(exp(-t^2),t=x..inf)
	// z=1/x
	// z*exp(-x^2)*(1/2 - 1/4 z^2 +3/8 z^4-15/16 z^6 + ...)
	long double z=1/x._DOUBLE_val,z2=z*z/2,res=0,pi=0.5;
	for (int n=0;;++n){
	  res += pi;
	  pi = -pi*(2*n+1)*z2;
	  if (std::abs(pi)<1e-16)
	    break;
	}
	erfc=2/std::sqrt(M_PI)*double(std::exp(-1/z/z)*z*res);
	return 1-erfc;
      }
      // a:=convert(series(erfc(x)*exp(x^2),x=X,24),polynom):; b:=subst(a,x=X+h):;
      if (absx>3 && absx<=5){
	// Digits:=30; evalf(symb2poly(subst(b,X,4),h))
	long double Zl=x._DOUBLE_val-4,res=0;
	long double taberf[]={0.9323573505930262336910814663629e-18,-0.5637770672346891132663122366369e-17,0.3373969923698176600796949171416e-16,-0.1997937342757611758805760309387e-15,0.1170311628709846086671746844320e-14,-0.6779078623355796103927587022047e-14,0.3881943235655598141099274338263e-13,-0.2196789805508621713379735090290e-12,0.1228090799753488475137690971599e-11,-0.6779634525816110746734938098109e-11,0.3694326453071165814527058450923e-10,-0.1986203171147991823844885265211e-9,0.1053084120195192127202221248092e-8,-0.5503368542058483880654875851859e-8,0.2833197888944711586737808090450e-7,-0.1435964425391227330876779173688e-6,0.7160456646037012951391007806358e-6,-0.3510366649840828060143659147374e-5,0.1690564925777814684043808381146e-4,-0.7990888030555549397777128848414e-4,0.3703524689955564311420527395424e-3,-0.1681182076746114476323671722330e-2,0.7465433244975570766528102818814e-2,-0.3238350609502145478059791886069e-1,0.1369994576250613898894451230325};
	unsigned N=sizeof(taberf)/sizeof(long double);
	for (unsigned i=0;i<N;i++){
	  res *= Zl;
	  res += taberf[i];
	}
	erfc = double(std::exp(-absx*absx)*res);
	return 1-erfc;
      }
      if (absx>5 && absx<=6.5){
	// Digits:=30; evalf(symb2poly(subst(b,X,5.75),h))
	long double Zl=x._DOUBLE_val-5.75,res=0;
	long double taberf[]={-0.3899077949952308336341205103240e-12,0.2064555746398182434172952813760e-13,-0.7079917646274828801231710613009e-12,-0.2043006626755557967429543230042e-12,-0.2664588032913413248313045028978e-11,-0.3182230691773937386262907009549e-11,-0.4508687162250923186571888867300e-12,-0.2818971742901571639195611759894e-11,-0.4771270499789446447101554995178e-11,0.2345376254096117543212461524786e-11,-0.6529305258174487397807156793042e-11,0.9817004987916722489154147719630e-12,0.2085292084663647123257426988484e-10,-0.1586500138272075839895787048265e-9,0.1056533982771769784560244626854e-8,-0.6964568016562765632682760517056e-8,0.4530411628438409475101496352516e-7,-0.2918364042864784155554051827879e-6,0.1859299481340192895158490699981e-5,-0.1171241494503672776195474661763e-4,0.7292428889065898343608897828825e-4,-0.4485956983428598110336671805311e-3,0.2725273842847326036320664185043e-2,-0.1634321814380709002113440890281e-1,0.9669877816971385564543076482100e-1};
	unsigned N=sizeof(taberf)/sizeof(long double);
	for (unsigned i=0;i<N;i++){
	  res *= Zl;
	  res += taberf[i];
	}
	erfc = double(std::exp(-absx*absx)*res);
	return 1-erfc;
      }
    }
    gen z=evalf_double(abs(x,contextptr),1,contextptr);
    // loss of accuracy
    int prec=decimal_digits(contextptr);
    int newprec,nbitsz=int(z._DOUBLE_val*z._DOUBLE_val/std::log(2.)),prec2=int(prec*std::log(10.0)/std::log(2.0)+.5);
    if (nbitsz>prec2){ 
      // use asymptotic expansion at z=inf
      z = accurate_evalf(inv(x,contextptr),prec2);
      gen z2=z*z/2,res=0,pi=inv(accurate_evalf(plus_two,prec2),contextptr),eps=accurate_evalf(pow(10,-prec,contextptr),prec2)/2;
      for (int n=0;;++n){
	res += pi;
	pi = -(2*n+1)*z2*pi;
	if (is_greater(eps,abs(pi,contextptr),contextptr))
	  break;
      }
      erfc=evalf(2*inv(sqrt(cst_pi,contextptr),contextptr),1,contextptr)*exp(-inv(z*z,contextptr),contextptr)*z*res;
      return 1-erfc;
    }
    if (z._DOUBLE_val>1)
      newprec = prec2+nbitsz+int(std::log(z._DOUBLE_val))+1;
    else
      newprec = prec2+2;
    // numerical computation of int(exp(-t^2),t=0..x) 
    // by series expansion at x=0
    // x*sum( (-1)^n*(x^2)^n/n!/(2*n+1),n=0..inf)
    z=accurate_evalf(x,newprec);
    gen z2=z*z,res=0,pi=1,eps=accurate_evalf(pow(10,-prec,contextptr),prec2)/2;
    for (int n=0;;){
      res += pi/(2*n+1);
      ++n;
      pi = -pi*z2/n;
      if (is_greater(eps,abs(pi,contextptr),contextptr))
	break;
    }
    res = evalf(2*inv(sqrt(cst_pi,contextptr),contextptr),1,contextptr)*z*res;
    erfc=accurate_evalf(1-res,prec2);
    return accurate_evalf(res,prec2);
  }
  gen erf(const gen & x,GIAC_CONTEXT){
    if (x==plus_inf)
      return plus_one;
    if (x==minus_inf)
      return minus_one;
    if (is_inf(x) || x==undef)
      return undef;
    if (is_zero(x))
      return x;
    gen erfc_;
    if (x.type==_DOUBLE_ || x.type==_CPLX || x.type==_REAL)
      return erf0(x,erfc_,contextptr);
    return symbolic(at_erf,x);
    gen e=x.evalf(1,contextptr);
#ifdef HAVE_LIBGSL
    if (e.type==_DOUBLE_)
      return gsl_sf_erf(e._DOUBLE_val);
#endif
#ifdef HAVE_LIBMPFR
    if (x.type==_REAL){
      mpfr_t gam;
      int prec=mpfr_get_prec(x._REALptr->inf);
      mpfr_init2(gam,prec);
      mpfr_erf(gam,x._REALptr->inf,GMP_RNDN);
      real_object res(gam);
      mpfr_clear(gam);
      return res;
    }
#endif
    return symbolic(at_erf,x);
  }
  gen _erf(const gen & args,GIAC_CONTEXT){
    return apply(args,erf,contextptr);
  }
  const string _erf_s("erf");
  unary_function_eval __erf(&_erf,&D_at_erf,&taylor_erf,_erf_s);
  unary_function_ptr at_erf (&__erf,0,true);

  gen d_erfc(const gen & e,GIAC_CONTEXT){
    return -d_erf(e,contextptr);
  }
  partial_derivative_onearg D_at_erfc(d_erfc);
  gen erfc(const gen & x,GIAC_CONTEXT){
    gen erfc_;
    if (x.type==_DOUBLE_ || x.type==_CPLX || x.type==_REAL){
      erf0(x,erfc_,contextptr);
      return erfc_;
    }
    return 1-symbolic(at_erf,x);
    gen e=x.evalf(1,contextptr);
#ifdef HAVE_LIBGSL
    if (e.type==_DOUBLE_)
      return gsl_sf_erfc(e._DOUBLE_val);
#endif
    return 1-symbolic(at_erf,x);
  }
  gen _erfc(const gen & args,GIAC_CONTEXT){
    return apply(args,erfc,contextptr);
  }
  const string _erfc_s("erfc");
  unary_function_eval __erfc(&_erfc,&D_at_erfc,_erfc_s);
  unary_function_ptr at_erfc (&__erfc,0,true);

  // assumes z>=1
  const double exp_minus_1_over_4=std::exp(-0.25);
  void sici_fg(double z,double & fz,double & gz){
    // int([u*]exp(-u)/(u^2+z^2),0,inf)
    // #nstep=1000 in [0,1], then * e^(-1/4)
    double nstep=250,a=0;
    fz=0; gz=0;
    for (;nstep>0.25;nstep*=exp_minus_1_over_4){
      double Fz=0,Gz=0;
      int N=int(nstep+.5);
      if (N<1)
	N=1;
      // Simpson over [a,a+1]
      double t=a,tmp,z2=z*z,Ninv=1./N;
      t = a+Ninv/2.;
      double expt=std::exp(-t),expfact=std::exp(-Ninv);
      for (int i=0;i<N;++i){ // middle points
	tmp = expt/(t*t+z2); 
	Fz += tmp;
	Gz += t*tmp;
	expt *= expfact;
	t += Ninv;
      }
      Fz *= 2; Gz *= 2;
      t = a+Ninv;
      expt=std::exp(-t);
      for (int i=1;i<N;++i){ 
	tmp = expt/(t*t+z2); // endpoint
	Fz += tmp;
	Gz += t*tmp;
	expt *= expfact;
	t += Ninv;
      }
      Fz *= 2; Gz *= 2;
      tmp=std::exp(-a)/(a*a+z*z); // endpoint
      Fz += tmp;
      Gz += a*tmp;
      a++;
      tmp=std::exp(-a)/(a*a+z*z); // endpoint
      Fz += tmp;
      Gz += a*tmp;
      fz += Fz/(6*N);
      gz += Gz/(6*N);
    }
    fz *= z;
  }

  // mode=1 Si only, mode==2 Ci only
  void sici(const gen & z0,gen & siz,gen & ciz,int prec,int mode,GIAC_CONTEXT){
    gen z=evalf_double(z0,1,contextptr);
    if (z0.type==_DOUBLE_ && prec>13)
      prec=13;
    if (z.type==_DOUBLE_ && prec<=13){
      double Z=z._DOUBLE_val,fz,gz;
#if defined HAVE_LIBGSL // && 0
      if (mode==1){
	siz=gsl_sf_Si(Z);
	return;
      }
      if (mode==2){
	if (Z<0)
	  ciz=gen(gsl_sf_Ci(-Z),M_PI);
	else
	  ciz=gsl_sf_Ci(Z);
	return;
      }
#endif
      if (Z>=40 || Z<=-40){
	// use series expansion at infinity
	// Si: 1/2*PI - 1/z*cos(z) - 1/z^2*sin(z) + 2/z^3*cos(z) + 6/z^4*sin(z) - 24/z^5*cos(z) - 120/z^6*sin(z) + 720/z^7*cos(z) + 5040/z^8*sin(z) + O(1/z^9)	= 1/z^8( (pi/2*z-cz)*z ...)
	// Ci:1/z*sin(z) - 1/z^2*cos(z) - 2/z^3*sin(z) + 6/z^4*cos(z) + 24/z^5*sin(z) - 120/z^6*cos(z) - 720/z^7*sin(z) + 5040/z^8*cos(z)
	long double sz=std::sin(Z);
	long double cz=std::cos(Z);
	long double invZ=1/Z;
	long double pi=invZ;
	long double sizd=M_PI/2,cizd=0;
	for (int n=1;;++n){
	  switch (n%4){
	  case 1:
	    sizd -= pi*cz;
	    cizd += pi*sz;
	    break;
	  case 2:
	    sizd -= pi*sz;
	    cizd -= pi*cz;
	    break;
	  case 3:
	    sizd += pi*cz;
	    cizd -= pi*sz;
	    break;
	  case 0:
	    sizd += pi*sz;
	    cizd += pi*cz;
	    break;
	  }
	  if (pi<1e-16 && pi>-1e-16)
	    break;
	  pi *= n*invZ;
	}
	siz = double(sizd);
	ciz = double(cizd);
	/*
	double z8=Z*Z;
	z8*=z8;
	z8*=z8;
	siz=((((((((M_PI/2*Z-cz)*Z-sz)*Z+2*cz)*Z+6*sz)*Z-24*cz)*Z-120*sz)*Z+720*cz)*Z+5040)/z8;
	ciz=((((((((sz)*Z-cz)*Z-2*sz)*Z+6*cz)*Z+24*sz)*Z-120*cz)*Z-720*sz)*Z+5040)/z8;
	*/
	return;
      }
      bool neg=Z<0;
      if (neg) Z=-Z;
      if (Z<=8){
	long double si=1,ci=0,z2=Z*Z,pi=1;
	for (long double n=1;;n++){
	  pi = -pi*z2/2/n;
	  if (std::abs(pi)<1e-15)
	    break;
	  ci += pi/(2*n);
	  pi /= (2*n+1);
	  si += pi/(2*n+1);
	}
	siz=double(si*Z);
	ciz=double(ci)+std::log(Z)+cst_euler_gamma;
      }
      // Digits:=30; 
      // a:=convert(series(Si(x),x=X,24),polynom):; b:=subst(a,x=X+h):; 
      // c:=convert(series(Ci(x),x=X,24),polynom):; d:=subst(c,x=X+h):; 
      if (Z>8 && Z<=12){
	long double Zl=Z-10,ress=0,resc=0;
	// evalf(symb2poly(subst(b,X,10),h))
	long double tabsi[]={-0.1189416530979229549237628888274e-25,0.1535274533580010929112764295251e-23,0.5949595042899632105631168585106e-23,-0.8443813569533766615075729254715e-21,-0.2314493549701736156628360858149e-20,0.3873999840160633357506392314227e-18,0.6314325721978121436699115557946e-18,-0.1454512706077874699091092895960e-15,-0.7966560468972510895796615830931e-16,0.4362654888706402638259447137657e-13,-0.2168961216095219762135154593351e-13,-0.1013156299801104705978428503741e-10,0.1511362807589119975552764377684e-10,0.1746079592768605042057230051441e-8,-0.4215112103088864749315276215891e-8,-0.2100827763936543847478539598343e-6,0.6768578499749603731655034337648e-6,0.1604768832104729109406272571796e-4,-0.6129221770634575410555873017014e-4,-0.6629459359846050378314018176664e-3,0.2619937628043294083259734576549e-2,0.1168258324144788179314042496964e-1,-0.3923347089937577354591945908199e-1,-0.5440211108893698134047476618518e-1,1.658347594218874049330971879387};
	// evalf(symb2poly(subst(d,X,10),h))
	long double tabci[]={-0.1031363603561483377414206719978e-24,0.1612636587877966055062771505636e-24,0.3224608155602632232346779331477e-22,0.1692359618263088728225313523766e-21,-0.1902123006508838353447776226346e-19,-0.3515542663511577547829464092398e-19,0.7651991655644681341759387392123e-17,0.8949035335940300076443987350711e-17,-0.2601535596364695291068633323882e-14,0.1492297097331336833098183955735e-16,0.6873241144777191906667850002758e-12,-0.6815990267379124633604933192924e-12,-0.1385917806100726110447949877846e-9,0.2729216217490738797188995072400e-9,0.2012042413927794669971423802132e-7,-0.5698511913235728873777759634867e-7,-0.1960205632344228431569991958688e-5,0.6982250568929225912819532196480e-5,0.1127699306487082097807313715561e-3,-0.4465373163022154950275303555869e-3,-0.3158611974102019356519036678335e-2,0.1189143127195082400887895227495e-1,0.3139641318985075293153170283171e-1,-0.8390715290764524522588639478252e-1,-0.4545643300445537263453282995265e-1};
	unsigned N=sizeof(tabsi)/sizeof(long double);
	for (unsigned i=0;i<N;i++){
	  ress *= Zl;
	  ress += tabsi[i];
	  resc *= Zl;
	  resc += tabci[i];
	}
	siz = double(ress);
	ciz = double(resc);
      }
      if (Z>12 && Z<=16){
	long double Zl=Z-14,ress=0,resc=0;
	// evalf(symb2poly(subst(b,X,14),h))
	long double tabsi[]={0.4669904672048171207530336926784e-25,-0.9121786080411940331362079921979e-24,-0.2670914489796227617771058426181e-22,0.5191793400822084981037748064947e-21,0.1272661806026797328866077155393e-19,-0.2467115711789321804079825604360e-18,-0.4949985288801387245191790228644e-17,0.9598983157925822187983801948972e-16,0.1531267131408434243475192526878e-14,-0.2983793513482513341175920391396e-13,-0.3640741847215956508751181538961e-12,0.7178248462975320846576276463699e-11,0.6346880693851116027919629277690e-10,-0.1280755925139665180236732107401e-8,-0.7574841879361234814393864898316e-8,0.1596987755895423889046710494666e-6,0.5558236361973162161308791480329e-6,-0.1276894967030880532682787806556e-4,-0.2074774698099097315147678882137e-4,0.5764575129603710060263581799849e-3,0.2308201450150731015736777254153e-3,-0.1190515482960545313209358846732e-1,0.2356412497996938805126656809673e-2,0.7075766826391930770538248600113e-1,1.556211050077665053703631892805};
	// evalf(symb2poly(subst(d,X,14),h))
	long double tabci[]={0.3590233410769916717046880365555e-25,0.1141530173081541381259185089975e-23,-0.2223749846955296806257123587014e-22,-0.5971404695883450494036609310213e-21,0.1158813189176036949355044721895e-19,0.2578286744261932540063357506560e-18,-0.4996649806240905822071761629155e-17,-0.8975849912962755765138097664122e-16,0.1743615446376990117868764420851e-14,0.2446425362039524592303128232180e-13,-0.4789979360061655431004102341364e-12,-0.5015267234200784104185754752298e-11,0.9985346701639551763361641926310e-10,0.7310382166844782587890268654904e-9,-0.1502609970889405144661999732289e-7,-0.6957715037844896183079492333545e-7,0.1519752625305293353227780696896e-5,0.3762397782384872898110463611961e-5,-0.9310463095669218193046370118801e-4,-0.8685445941902185669380116006632e-4,0.2944299062831149098901196099828e-2,0.7349281020023392469738413739405e-4,-0.3572765356616331098087728605326e-1,0.9766944157702399589209205476559e-2,0.6939635592758454727438326824349e-1};
	unsigned N=sizeof(tabsi)/sizeof(long double);
	for (unsigned i=0;i<N;i++){
	  ress *= Zl;
	  ress += tabsi[i];
	  resc *= Zl;
	  resc += tabci[i];
	}
	siz = double(ress);
	ciz = double(resc);
      }
      if (Z>16 && Z<=20){
	long double Zl=Z-18,ress=0,resc=0;
	// evalf(symb2poly(subst(b,X,18),h))
	long double tabsi[]={-0.5458114686729331288343465771243e-25,-0.8535297797031477670863377809864e-25,0.3197608273507048053143023106525e-22,0.1246720240695204431925986220740e-22,-0.1566961070193290350813150654958e-19,0.1120194984050557421284180749579e-19,0.6303732529750142701943836944779e-17,-0.1093865529886308098492342051935e-16,-0.2034127294195017035840150503185e-14,0.5391577496552481117289518608625e-14,0.5113326950563000131064845650605e-12,-0.1755022270416406807283368583049e-11,-0.9642855115888351354955178064729e-10,0.3896583006063678015119758986001e-9,0.1297996094042326342264874814159e-7,-0.5741523749720224155210168829004e-7,-0.1165550911079956021485851292181e-5,0.5260587329402175094599529167340e-5,0.6336730653896082851113120333026e-4,-0.2682059741680869878136005839481e-3,-0.1788149401756335521591446074162e-2,0.6231324073036488917300258847156e-2,0.1950106172093382517043203223869e-1,-0.4172151370953756131945311580259e-1,1.536608096861185462361173893885};
	// evalf(symb2poly(subst(d,X,18),h))
	long double tabci[]={0.4605275954646862944758091918129e-26,-0.1349519299011741936395476560680e-23,-0.1307682310402805091861215069192e-23,0.7246129197780950670155731691447e-21,-0.1246271694903983843859628582537e-21,-0.3225645499768452910431606017148e-18,0.3989357101823634351098507496171e-18,0.1165949042798889880769685702457e-15,-0.2573970823221046721180373208376e-15,-0.3334412559655945734048458981847e-13,0.1020420669194920098922521974402e-12,0.7298982331099322669182907795178e-11,-0.2745289334776461767357052328935e-10,-0.1171271721735576258954174760357e-8,0.4994621772350737506214065546293e-8,0.1300541789133113153488950327832e-6,-0.5864843122180159884437825091296e-6,-0.9221666449442426870101941344121e-5,0.4080390556697611593508057267335e-4,0.3702810510394427353858531043357e-3,-0.1453024604178293371013986047647e-2,-0.6848923209258520415117450659078e-2,0.1984174958896001500414615659760e-1,0.3668426156911556360089444693281e-1,-0.4347510299950100478344114920850e-1};
	unsigned N=sizeof(tabsi)/sizeof(long double);
	for (unsigned i=0;i<N;i++){
	  ress *= Zl;
	  ress += tabsi[i];
	  resc *= Zl;
	  resc += tabci[i];
	}
	siz = double(ress);
	ciz = double(resc);
      }
      if (Z>20 && Z<=24){
	long double Zl=Z-22,ress=0,resc=0;
	// evalf(symb2poly(subst(b,X,22),h))
	long double tabsi[]={0.3374532347188926046754549614410e-25,0.9070719831057349388641148466696e-24,-0.2050470901003743543468128570765e-22,-0.4594235790411492575575248293744e-21,0.1043077621302187913107676015839e-19,0.1910932463670378689024003307389e-18,-0.4360804036884663718921536467101e-17,-0.6379227769542456535760586922387e-16,0.1464704096058402216999762277897e-14,0.1660737164480043388814198356103e-13,-0.3842072518362213364190884929092e-12,-0.3249870299707868231630032306411e-11,0.7591535353403568501001826557431e-10,0.4554066973680817003059315739181e-9,-0.1077692372064168550494073939261e-7,-0.4274452798917825988601610781581e-7,0.1030486252719749122401833140117e-5,0.2434192199039423608811016740448e-5,-0.6042868558775171616279520355687e-4,-0.7128407780774990405916452478794e-4,0.1868111001271415320776084256504e-2,0.7554565401849909491334584321080e-3,-0.2271723850350373234098156405564e-1,-0.4023322404729034509859207643533e-3,1.616083736594366543114431027190};
	// evalf(symb2poly(subst(d,X,22),h))
	long double tabci[]={-0.3765941636516127667565314425887e-25,0.8496429623966734299562910340895e-24,0.2089658856058364573723209071877e-22,-0.4733667956308362832271784627271e-21,-0.9616010057686047674500629395344e-20,0.2188568968552838678620763986263e-18,0.3594653788547495746845572205294e-17,-0.8227046014999096676822143059749e-16,-0.1063984308386789623300816541173e-14,0.2451691569424535439026086961327e-13,0.2414081477262495987003288682969e-12,-0.5610174912354327213432820922583e-11,-0.4025699100579353159380991679438e-10,0.9460092423484752010231070364358e-9,0.4662816039604651635851886085410e-8,-0.1112697436829371328490149410190e-6,-0.3461508105965530058056263010906e-6,0.8452332929171272742294005025613e-5,0.1452920166854777888625357916022e-4,-0.3688187418989882360609440618036e-3,-0.2737432060552935755550251671918e-3,0.7538061305932840665075723970069e-2,0.1234183502875539666044793790453e-2,-0.4545276483611986938428066996417e-1,0.1640691915737749726680980654224e-2};
	unsigned N=sizeof(tabsi)/sizeof(long double);
	for (unsigned i=0;i<N;i++){
	  ress *= Zl;
	  ress += tabsi[i];
	  resc *= Zl;
	  resc += tabci[i];
	}
	siz = double(ress);
	ciz = double(resc);
      }
      if (Z>24 && Z<=28){
	long double Zl=Z-26,ress=0,resc=0;
	// evalf(symb2poly(subst(b,X,26),h))
	long double tabsi[]={0.1887595283929249736721564257313e-26,-0.1139208683110194949664423805669e-23,0.1279293933254064768731875653817e-24,0.5994732455107313947881755894338e-21,-0.6938324133833276535897299455340e-21,-0.2609102593630429565352138793600e-18,0.5435603049161089391681629274496e-18,0.9202351759037677282055983009108e-16,-0.2642554055043397588710384653172e-15,-0.2565196994817985742348447471878e-13,0.8979673506711721488705329079164e-13,0.5477157384399615846158662911243e-11,-0.2160870416302321163457598611639e-10,-0.8604342982500329851062267981807e-9,0.3594354283136494037733432184009e-8,0.9424490157133423171849887794118e-7,-0.3925808731949007274318185773343e-6,-0.6671455253428628462019837043580e-5,0.2584931618742702225630929467268e-4,0.2717002055000108101387898426941e-3,-0.8869394862544894803875105833456e-3,-0.5192726828102161497754294326711e-2,0.1187673367608361403346165544076e-1,0.2932917117229241298137781896944e-1,1.544868862986338557887737260292};
	// evalf(symb2poly(subst(d,X,26),h))
	long double tabci[]={0.4649965583842175529710094230187e-25,0.2092492486523365857112019931102e-25,-0.2673016958726137563529919009844e-22,0.1734025691834661144821904431335e-22,0.1282183044327200389255292828905e-19,-0.2098089236197822391981775009963e-19,-0.5037728169303807883904931875640e-17,0.1257389975122490684981594653098e-16,0.1585106167468127577346565423777e-14,-0.5084215538447167537343863000543e-14,-0.3884144934780823340053700189354e-12,0.1455636906316532474122453840240e-11,0.7154603949717881354941340050201e-10,-0.2926061246791173118887623808558e-9,-0.9458827111752210232753218600667e-8,0.3976537604156996965393015686438e-7,0.8424410070112619898093729646951e-6,-0.3418064228754065522708009442873e-5,-0.4606856152608972144114035752184e-4,0.1664083688146614131038736660077e-3,0.1330470954446840977471854832752e-2,-0.3758634727512557357034117912936e-2,-0.1514307620917034875601537686776e-1,0.2488151239725539779697630391881e-1,0.2829515103175713190842112993963e-1};
	unsigned N=sizeof(tabsi)/sizeof(long double);
	for (unsigned i=0;i<N;i++){
	  ress *= Zl;
	  ress += tabsi[i];
	  resc *= Zl;
	  resc += tabci[i];
	}
	siz = double(ress);
	ciz = double(resc);
      }
      if (Z>28 && Z<=32){
	long double Zl=Z-30,ress=0,resc=0;
	// evalf(symb2poly(subst(b,X,30),h))
	long double tabsi[]={-0.3146268636172325352644123013773e-25,0.7254967655706034123075144269188e-24,0.1720174348300531678711926627452e-22,-0.3968569600567226603419285560888e-21,-0.7804624742706385979212734054562e-20,0.1797226284123711388607362992709e-18,0.2882267004838688287625265808775e-17,-0.6604428381371346662942676361200e-16,-0.8462601728339625717076212469620e-15,0.1921641824294710456184856743277e-13,0.1918622280824819669976604364077e-12,-0.4293140602943259867747088385429e-11,-0.3236373134197171392533933286615e-10,0.7078744765376845042933221167302e-9,0.3867645142776632643256131482060e-8,-0.8169087758905646352536906937371e-7,-0.3060269540778692505740670587310e-6,0.6120146081774563343983258719578e-5,0.1450591123346014175404634002192e-4,-0.2651270545919250186145646630982e-3,-0.3497315371034554476595429181296e-3,0.5419736490383548422011597468850e-2,0.3119763955955768506415340725399e-2,-0.3293438746976205966625829690981e-1,1.566756540030351110983731309007};
	// evalf(symb2poly(subst(d,X,30),h))
	long double tabci[]={-0.2905250884383236235247049545527e-25,-0.7522147039999421611567058749977e-24,0.1735426584530328263231524611399e-22,0.3754699661541031707840824496312e-21,-0.8657191669382300006764470296916e-20,-0.1541017153772568372597216447575e-18,0.3541395808613865455763442962592e-17,0.5090911414277503198761655016250e-16,-0.1161952195659956230574302196387e-14,-0.1318846041657448769184055840206e-13,0.2975305748369348193892361410833e-12,0.2592737707940046163669864700563e-11,-0.5742839188863978180610659896436e-10,-0.3707322050098237760097147226680e-9,0.7983406320009368395835407645496e-8,0.3641426564991240108255042664073e-7,-0.7507713655866267126314943513395e-6,-0.2264698987507883174282691040657e-5,0.4355811042213145491501868952153e-4,0.7862739829137060637224202271609e-4,-0.1341741499597397210639678491289e-2,-0.1220985799040877684843355198131e-2,0.1638149848494348313828544726235e-1,0.5141714996252801690622071553807e-2,-0.3303241728207114377922644096301e-1};
	unsigned N=sizeof(tabsi)/sizeof(long double);
	for (unsigned i=0;i<N;i++){
	  ress *= Zl;
	  ress += tabsi[i];
	  resc *= Zl;
	  resc += tabci[i];
	}
	siz = double(ress);
	ciz = double(resc);
      }
      if (Z>32 && Z<=36){
	long double Zl=Z-34,ress=0,resc=0;
	// evalf(symb2poly(subst(b,X,34),h))
	long double tabsi[]={0.3942701136493758442738611783193e-25,0.2833482488121121562882295508384e-25,-0.2240325735608330926855387141189e-22,0.7077167897405254555239797927716e-23,0.1062083455080499454129483304028e-19,-0.1296740970607253649730617065007e-19,-0.4125642950551060755694922602219e-17,0.8186731227770529065956573976157e-17,0.1284716025190091972930389282765e-14,-0.3331328985339278802607691455602e-14,-0.3121425078659463413394477246734e-12,0.9467354869456215167864364128202e-12,0.5717466510263820214761281907592e-10,-0.1880717150280359445242629986190e-9,-0.7546218335244881072513219000215e-8,0.2525355992897561612710992980882e-7,0.6743126364121930403816188890288e-6,-0.2149417013042573939863424458069e-5,-0.3721263582524267372627099974480e-4,0.1039917503622257770319033100678e-3,0.1091628590022319276482516968274e-2,-0.2344369704089296833875653244402e-2,-0.1270781662145181668046346804335e-1,0.1556125547411834767154373923623e-1,1.595256185182468624967114677624};
	// evalf(symb2poly(subst(d,X,34),h))
	long double tabci[]={-0.1989238045881366884396486216401e-26,0.9603929081698055330298655453613e-24,0.1785462486600632168852339199134e-24,-0.4994884203374138675491408864852e-21,0.3922698307156298349251565849186e-21,0.2148764364742271875693500554094e-18,-0.3483461372134677428776900288130e-18,-0.7495914823918311062321675163567e-16,0.1730727343462677961269223481865e-15,0.2069644169058104778527335300702e-13,-0.5867510845332297156886886169438e-13,-0.4387405297943553775874545604808e-11,0.1397446669253207042207953256760e-10,0.6866413079078539001893686720856e-9,-0.2296062969528945671977337418042e-8,-0.7526095978772233458886342141369e-7,0.2479954929300904010747076408018e-6,0.5360278479932581967324034110246e-5,-0.1619607535570418039830425866209e-4,-0.2210046023539758078033314699734e-3,0.5534219043710011378918782268888e-3,0.4305022897404827350706361160583e-2,-0.7413599071494898236031397261835e-2,-0.2495794925837074078235212022679e-1,0.1626491643735576698165635194377e-1};
	unsigned N=sizeof(tabsi)/sizeof(long double);
	for (unsigned i=0;i<N;i++){
	  ress *= Zl;
	  ress += tabsi[i];
	  resc *= Zl;
	  resc += tabci[i];
	}
	siz = double(ress);
	ciz = double(resc);
      }
      if (Z>36 && Z<=40){
	long double Zl=Z-38,ress=0,resc=0;
	// evalf(symb2poly(subst(b,X,38),h))
	long double tabsi[]={-0.2392585617788018962646757828598e-25,-0.6575431727115741067484987610228e-24,0.1413726322474858257725781544417e-22,0.3273737185831529184312733883924e-21,-0.6971366411141976840382089119809e-20,-0.1343722941531172977052740901750e-18,0.2818131926670801219167357117818e-17,0.4456155806753828158414772784641e-16,-0.9138685919934995126194762439509e-15,-0.1164855661769990237307373295310e-13,0.2314363358966319503074680596234e-12,0.2327043421113577898681844454632e-11,-0.4423649753121148492913255999601e-10,-0.3413423247744706973417929454576e-9,0.6100985644119355599936971303758e-8,0.3483915381970249998364550866465e-7,-0.5705724285037609250439734410714e-6,-0.2292102947589923998615699771400e-5,0.3301272634561039972744957647445e-4,0.8640908538565720173919303844142e-4,-0.1017258860929286693108579223376e-2,-0.1518531271099188526195110204881e-2,0.1246413777530741664504711921590e-1,0.7799173123931192562955174996028e-2,1.545492937235698740561891130750};
	// evalf(symb2poly(subst(d,X,38),h))
	long double tabci[]={0.2755791623901057620946655522961e-25,-0.5942947240530035802859253027110e-24,-0.1501343221496874685931095691537e-22,0.3214487908378114008198384990138e-21,0.6802473399014152745014932500302e-20,-0.1438706210063382233790032141791e-18,-0.2516685599977182931504499600334e-17,0.5224619763598943054009172687575e-16,0.7435302907498417390873416469940e-15,-0.1502856747601093160946258165352e-13,-0.1706516669026401836883958369119e-12,0.3322517447592247618873995073660e-11,0.2938003484322757578999393723335e-10,-0.5429671029677216915133045537076e-9,-0.3623244387376009573032456389105e-8,0.6223561267305907587303733846902e-7,0.3003453535122368377351197880341e-6,-0.4643099720179236845576657726334e-5,-0.1523777445337208725039230269865e-4,0.2008948838914583162973945866431e-3,0.4061768073150514090668429199326e-3,-0.4114703864552309315370227396080e-2,-0.4230290732342083420531097274815e-2,0.2513351694861302256806674303698e-1,0.7129761801971379713551376511546e-2};
	unsigned N=sizeof(tabsi)/sizeof(long double);
	for (unsigned i=0;i<N;i++){
	  ress *= Zl;
	  ress += tabsi[i];
	  resc *= Zl;
	  resc += tabci[i];
	}
	siz = double(ress);
	ciz = double(resc);
      }
      if (Z>40 && Z<40) { // not used anymore, too slow
	sici_fg(Z,fz,gz);
	siz=M_PI/2-fz*std::cos(Z)-gz*std::sin(Z);
	ciz=fz*std::sin(Z)-gz*std::cos(Z);
      }
      if (neg){
	siz = -siz;
	ciz += cst_i*cst_pi;
      }
      return;
    }
    z=evalf_double(abs(z0,contextptr),1,contextptr);
    if (z.type!=_DOUBLE_)
      settypeerr(); 
    // find number of digits that must be added to prec
    // n^n/n! equivalent to e^n*sqrt(2*pi*n)
    int newprec,nbitsz=int(z._DOUBLE_val/std::log(2.)),prec2=int(prec*std::log(10.0)/std::log(2.0)+.5);
    if (nbitsz>prec2){ 
      // use asymptotic expansion at z=inf
      z = accurate_evalf(z0,prec2);
      gen sinz=sin(z,contextptr),cosz=cos(z,contextptr);
      gen invc=1,invs=0,pi=1,eps=accurate_evalf(pow(10,-prec,contextptr),prec2)/2;
      for (int n=1;;++n){
	if (is_greater(eps,abs(pi,contextptr),contextptr))
	  break;
	pi = (n*pi)/z;
	if (n%2){
	  if (n%4==1)
	    invs += pi;
	  else
	    invs -= pi;
	}
	else {
	  if (n%4==0)
	    invc += pi;
	  else
	    invc -= pi;
	}
      }
      siz=m_pi(prec2)/2-cosz/z*invc-sinz/z*invs;
      ciz=sinz/z*invc-cosz/z*invs;
      return;
    }
    // use series expansion at z=0
    if (z._DOUBLE_val>1)
      newprec = prec2+nbitsz+int(std::log(z._DOUBLE_val)/2)+1;
    else
      newprec = prec2+2;
    z = accurate_evalf(z0,newprec);
    gen si=1,ci=0,z2=z*z,pi=1,eps=accurate_evalf(pow(10,-prec,contextptr),newprec)/2;
    for (int n=1;;n++){
      pi = pi*z2/(2*n*(2*n-1));
      if (is_greater(eps,abs(pi,contextptr),contextptr))
	break;
      if (mode!=1){
	if (n%2)
	  ci -= pi/(2*n);
	else
	  ci += pi/(2*n);
      }
      if (mode!=2){
	if (n%2)
	  si -= pi/((2*n+1)*(2*n+1));
	else
	  si += pi/((2*n+1)*(2*n+1));
      }
    }
    if (mode!=2)
      siz=si*accurate_evalf(z0,prec2);
    if (mode!=1){
      ciz=ci+ln(z,contextptr)+m_gamma(newprec);
      ciz=accurate_evalf(ciz,prec2);
    }
  }

  extern unary_function_ptr at_SiCi_f ;
  gen taylor_SiCi_f(const gen & lim_point,const int ordre,const unary_function_ptr & f, int direction,gen & shift_coeff,GIAC_CONTEXT){
    if (ordre<0)
      return 0;
    if (!is_inf(lim_point))
      return taylor(lim_point,ordre,f,0,shift_coeff,contextptr);
    shift_coeff=1;
    // f(x)=1/x* sum( +/-(2*k)!*x^(-2k) )
    gen tmp(1);
    vecteur v;
    for (int n=0;n<=ordre;){
      v.push_back(tmp);
      v.push_back(0);
      n +=2 ;
      tmp=-gen((n-1)*n)*tmp;
    }
    v.push_back(undef);
    return v;
  }
  gen _SiCi_g(const gen & args,GIAC_CONTEXT);
  gen d_SiCi_f(const gen & args,GIAC_CONTEXT){
    return -_SiCi_g(args,contextptr);
  }
  partial_derivative_onearg D_at_SiCi_f(&d_SiCi_f);
  gen _Si(const gen & args,GIAC_CONTEXT);
  gen _Ci(const gen & args,GIAC_CONTEXT);
  gen _SiCi_f(const gen & args,GIAC_CONTEXT){
    if (is_inf(args))
      return 0;
    if (is_zero(args))
      return unsigned_inf;
    if (is_undef(args))
      return args;
    if (args.type==_DOUBLE_ || args.type==_REAL)
      return _Ci(args,contextptr)*sin(args,contextptr)+(evalf(cst_pi/2,1,contextptr)-_Si(args,contextptr))*cos(args,contextptr);
    return symbolic(at_SiCi_f,args);
  }
  const string _SiCi_f_s("SiCi_f");
  unary_function_eval __SiCi_f(&_SiCi_f,&D_at_SiCi_f,&taylor_SiCi_f,_SiCi_f_s);
  unary_function_ptr at_SiCi_f (&__SiCi_f,0,true);

  extern unary_function_ptr at_SiCi_g ;
  gen taylor_SiCi_g(const gen & lim_point,const int ordre,const unary_function_ptr & f, int direction,gen & shift_coeff,GIAC_CONTEXT){
    if (ordre<0)
      return 0;
    if (!is_inf(lim_point))
      return taylor(lim_point,ordre,f,0,shift_coeff,contextptr);
    shift_coeff=2;
    // f(x)=sum( +/-(2*k+1)!*x^(-2k+2) )
    gen tmp(1);
    vecteur v;
    for (int n=1;n<=ordre;){
      v.push_back(tmp);
      v.push_back(0);
      n +=2 ;
      tmp=-gen((n-1)*n)*tmp;
    }
    v.push_back(undef);
    return v;
  }
  gen d_SiCi_g(const gen & args,GIAC_CONTEXT){
    return inv(args,contextptr)+_SiCi_f(args,contextptr);
  }
  partial_derivative_onearg D_at_SiCi_g(&d_SiCi_g);
  gen _SiCi_g(const gen & args,GIAC_CONTEXT){
    if (is_inf(args))
      return 0;
    if (is_zero(args))
      return unsigned_inf;
    if (is_undef(args))
      return args;
    if (args.type==_DOUBLE_ || args.type==_REAL)
      return -_Ci(args,contextptr)*cos(args,contextptr)+(evalf(cst_pi/2,1,contextptr)-_Si(args,contextptr))*sin(args,contextptr);
    return symbolic(at_SiCi_g,args);
  }
  const string _SiCi_g_s("SiCi_g");
  unary_function_eval __SiCi_g(&_SiCi_g,&D_at_SiCi_g,&taylor_SiCi_g,_SiCi_g_s);
  unary_function_ptr at_SiCi_g (&__SiCi_g,0,true);

  extern unary_function_ptr at_Si;
  gen Si_replace(const gen & g,GIAC_CONTEXT){
    return cst_pi_over_2-_SiCi_f(g,contextptr)*cos(g,contextptr)-_SiCi_g(g,contextptr)*sin(g,contextptr);
  }
  gen taylor_Si(const gen & lim_point,const int ordre,const unary_function_ptr & f, int direction,gen & shift_coeff,GIAC_CONTEXT){
    if (ordre<0){
      limit_tractable_functions.push_back(&at_Si);
      limit_tractable_replace.push_back(Si_replace);
      return 1;
    }
    shift_coeff=0;
    if (is_zero(lim_point)){
      shift_coeff=1;
      vecteur v;
      gen pi(1);
      for (int i=0;i<=ordre;){
	v.push_back(plus_one/pi);
	v.push_back(0);
	i += 2;
	pi = -(i*(i+1))*pi;
      }
      return v;
    }
    if (!is_inf(lim_point))
      return taylor(lim_point,ordre,f,direction,shift_coeff,contextptr);
    settypeerr();
    return undef;
  }
  gen d_Si(const gen & args,GIAC_CONTEXT){
    return sin(args,contextptr)/args;
  }
  partial_derivative_onearg D_at_Si(&d_Si);
  gen _Si(const gen & args,GIAC_CONTEXT){
    if (is_zero(args))
      return args;
    if (is_undef(args))
      return args;
    if (is_inf(args)){
      if (args==plus_inf)
	return cst_pi_over_2;
      if (args==minus_inf)
	return -cst_pi_over_2;
      return undef;
    }
    if (args.is_symb_of_sommet(at_neg))
      return -_Si(args._SYMBptr->feuille,contextptr);
    if (args.type!=_DOUBLE_ && args.type!=_REAL && args.type!=_CPLX)
      return symbolic(at_Si,args);
    gen si,ci;
    sici(args,si,ci,decimal_digits(contextptr),1,contextptr);
    return si;
  }
  const string _Si_s("Si");
  unary_function_eval __Si(&_Si,&D_at_Si,&taylor_Si,_Si_s);
  unary_function_ptr at_Si (&__Si,0,true);

  gen Ci_replace(const gen & g,GIAC_CONTEXT){
    return _SiCi_f(g,contextptr)*sin(g,contextptr)-_SiCi_g(g,contextptr)*cos(g,contextptr);
  }
  gen _Ci0(const gen &,GIAC_CONTEXT);
  gen Ci_replace0(const gen & g,GIAC_CONTEXT){
    return _Ci0(g,contextptr)+cst_euler_gamma+ln(abs(g,contextptr),contextptr);
  }
  gen taylor_Ci(const gen & lim_point,const int ordre,const unary_function_ptr & f, int direction,gen & shift_coeff,GIAC_CONTEXT){
    if (ordre<0){
      limit_tractable_functions.push_back(&at_Ci);
      limit_tractable_replace.push_back(Ci_replace);
      return 1;
    }
    shift_coeff=0;
    if (!is_inf(lim_point))
      return taylor(lim_point,ordre,f,direction,shift_coeff,contextptr);
    settypeerr();
    return undef;
  }
  gen d_Ci(const gen & args,GIAC_CONTEXT){
    return cos(args,contextptr)/args;
  }
  partial_derivative_onearg D_at_Ci(&d_Ci);
  gen _Ci(const gen & args,GIAC_CONTEXT){
    if (is_zero(args))
      return minus_inf;
    if (is_undef(args))
      return args;
    if (is_inf(args)){
      if (args==plus_inf)
	return 0;
      if (args==minus_inf)
	return cst_pi*cst_i;
      return undef;
    }
    if (args.type!=_DOUBLE_ && args.type!=_REAL && args.type!=_CPLX)
      return symbolic(at_Ci,args);
    gen si,ci;
    sici(args,si,ci,decimal_digits(contextptr),2,contextptr);
    return ci;
  }
  const string _Ci_s("Ci");
  unary_function_eval __Ci(&_Ci,&D_at_Ci,&taylor_Ci,_Ci_s);
  unary_function_ptr at_Ci (&__Ci,0,true);

  gen d_Ci0(const gen & args,GIAC_CONTEXT){
    return (cos(args,contextptr)-1)/args;
  }
  partial_derivative_onearg D_at_Ci0(&d_Ci0);
  gen taylor_Ci0(const gen & lim_point,const int ordre,const unary_function_ptr & f, int direction,gen & shift_coeff,GIAC_CONTEXT){
    if (ordre<0)
      return 0;
    if (!is_zero(lim_point))
      return taylor(lim_point,ordre,f,0,shift_coeff,contextptr);
    shift_coeff=2;
    // sum( (-1)^k/(2*k)/(2*k)! * x^(2k) )
    gen tmp(1);
    vecteur v;
    for (int n=0;n<=ordre;){
      n +=2 ;
      tmp=-gen((n-1)*n)*tmp;
      v.push_back(inv(n*tmp,contextptr));
      v.push_back(0);
    }
    v.push_back(undef);
    return v;
  }
  gen _Ci0(const gen & args,GIAC_CONTEXT){
    if (is_zero(args))
      return 0;
    if (is_undef(args))
      return args;
    if (is_inf(args))
      return minus_inf;
    if (args.type!=_DOUBLE_ && args.type!=_REAL && args.type!=_CPLX)
      return symbolic(at_Ci0,args);
    gen si,ci;
    sici(args,si,ci,decimal_digits(contextptr),2,contextptr);
    return ci-evalf(cst_euler_gamma,1,contextptr)-ln(args,contextptr);
  }
  const string _Ci0_s("Ci0");
  unary_function_eval __Ci0(&_Ci0,&D_at_Ci0,&taylor_Ci0,_Ci0_s);
  unary_function_ptr at_Ci0 (&__Ci0,0,true); /* FIXME should not registered */

  extern unary_function_ptr at_Ei ;
  extern unary_function_ptr at_Ei_f ;
  gen _Ei_f(const gen & args,GIAC_CONTEXT);
  gen taylor_Ei_f(const gen & lim_point,const int ordre,const unary_function_ptr & f, int direction,gen & shift_coeff,GIAC_CONTEXT){
    if (ordre<0)
      return 0;
    if (!is_inf(lim_point))
      setsizeerr();
    shift_coeff=1;
    // f(x)=1/x* sum( k!/x^(k) )
    gen tmp(1);
    vecteur v;
    for (int n=1;n<=ordre+1;n++){
      v.push_back(tmp);
      tmp=n*tmp;
    }
    v.push_back(undef);
    return v;
  }
  gen d_Ei_f(const gen & args,GIAC_CONTEXT){
    return -_Ei_f(args,contextptr);
  }
  partial_derivative_onearg D_at_Ei_f(&d_Ei_f);
  gen _Ei_f(const gen & args,GIAC_CONTEXT){
    if (is_inf(args))
      return 0;
    if (is_zero(args))
      return unsigned_inf;
    if (is_undef(args))
      return args;
    return symbolic(at_Ei_f,args);
  }
  const string _Ei_f_s("Ei_f");
  unary_function_eval __Ei_f(&_Ei_f,&D_at_Ei_f,&taylor_Ei_f,_Ei_f_s);
  unary_function_ptr at_Ei_f (&__Ei_f,0,true);
  gen Ei_replace(const gen & g,GIAC_CONTEXT){
    return _Ei_f(g,contextptr)*exp(g,contextptr);
  }  
  gen _Ei0(const gen & args,GIAC_CONTEXT);
  gen Ei_replace0(const gen & g,GIAC_CONTEXT){
    return _Ei0(g,contextptr)+cst_euler_gamma+ln(abs(g,contextptr),contextptr);
  }
  gen taylor_Ei(const gen & lim_point,const int ordre,const unary_function_ptr & f, int direction,gen & shift_coeff,GIAC_CONTEXT){
    if (ordre<0){
      limit_tractable_functions.push_back(&at_Ei);
      limit_tractable_replace.push_back(Ei_replace);
      return 1;
    }
    shift_coeff=0;
    if (!is_inf(lim_point))
      return taylor(lim_point,ordre,f,direction,shift_coeff,contextptr);
    settypeerr();
    return undef;
  }
  gen d_Ei(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      return exp(args,contextptr)/args;
    vecteur v=*args._VECTptr;
    if (v.size()==1)
      return exp(v.front(),contextptr)/v.front();
    setdimerr();
    return undef;
  }
  partial_derivative_onearg D_at_Ei(&d_Ei);
  gen Ei(const gen & args,GIAC_CONTEXT){
    if (is_zero(args))
      return minus_inf;
    if (args==plus_inf || is_undef(args))
      return args;
    if (args==minus_inf)
      return 0;
    if (is_inf(args))
      return undef;
    if (args.type!=_DOUBLE_ && args.type!=_REAL && args.type!=_CPLX)
      return symbolic(at_Ei,args);
    gen z=evalf_double(abs(args,contextptr),1,contextptr);
    if (z.type!=_DOUBLE_)
      settypeerr(); 
    int prec=decimal_digits(contextptr);
    if (args.type==_DOUBLE_ && prec>13)
      prec=13;
    if (args.type==_DOUBLE_ && prec<=13){
      double z=args._DOUBLE_val;
#ifdef HAVE_LIBGSL
      return gsl_sf_expint_Ei(z);
#endif
      if (z>=-4.8 && z<=40){
	// ? use __float80 or __float128
	/*
#ifdef __SSE__
#if defined __x86_64__ && defined __SSE_4_2__
	__float128 ei=0.0q,pi=1.0q;
#else
	__float80 ei=0.0w,pi=1.0w;
#endif // __SSE4_2__
	*/
	long double ei=0.0,pi=1.0,Z=z;
	for (long double n=1;;n++){
	  pi = pi*Z/n;
	  if (pi<1e-16 && pi>-1e-16)
	    break;
	  ei += pi/n;
	}
	ei=ei+std::log(std::abs(z))+0.577215664901532860610;
	return double(ei);
      }
      // a:=convert(series(Ei(x)*exp(-x)*x,x=X,24),polynom):; b:=subst(a,x=X+h):; 
      if (z>=-6.8 && z<=-4.8){
	// X:=-5.8; evalf(symb2poly(b,h),30)
	long double Z=z+5.8,res=0;
	long double tabei[]={-0.3151760388807517547897224622361e-20,-0.1956623502102099783599191666531e-19,-0.1217520387814662777242246174541e-18,-0.7595047623259136899131074978509e-18,-0.4750592523487122640844934717658e-17,-0.2979985874126626857226335496504e-16,-0.1875104560810577497563994966761e-15,-0.1183827325179695999391747354358e-14,-0.7501052898029529880772284438566e-14,-0.4771568760084635063692127230846e-13,-0.3048290706102002977129135629199e-12,-0.1956495512880190361879173037779e-11,-0.1262184835369108393063770847820e-10,-0.8188615023197271684637356284403e-10,-0.5345588144426140308176155460920e-9,-0.3513748178583494436182263562858e-8,-0.2327416524478792602782072181673e-7,-0.1554886566283128983353292920891e-6,-0.1048823536032497188023815631115e-5,-0.7151898287935501422833421755707e-5,-0.4937247597221480470432894324202e-4,-0.3456479830723309944342902721052e-3,-0.2458903778230241247990734517862e-2,-0.1781701304096371729351217642124e-1,0.8681380405349396412209368563589};
	unsigned N=sizeof(tabei)/sizeof(long double);
	for (unsigned i=0;i<N;i++){
	  res *= Z;
	  res += tabei[i];
	}
	return double(res*std::exp(z)/z);
      }
      if (z>=-10.4 && z<=-6.8){
	// X:=-8.6; evalf(symb2poly(b,h),30)
	long double Z=z+8.6,res=0;
	long double tabei[]={-0.3038274728374471199550377e-24,-0.2779136645187427028874693e-23,-0.25469559-0.3038274728374471199535898278331e-24,-0.2779136645187427028908730164655e-23,-0.2546955934813756337104395625899e-22,-0.2338906925011672805172017690719e-21,-0.2152486734154319991905416359525e-20,-0.1985479703050588294559030935503e-19,-0.1835926261877877389191827053171e-18,-0.1702095774807244266388321827485e-17,-0.1582462980029593441341590109962e-16,-0.1475687715549251315001108783383e-15,-0.1380597694124124865185429282546e-14,-0.1296174167643919255669841368436e-13,-0.1221540406726513790755510762310e-12,-0.1155953021811538418057560238246e-11,-0.1098796276688146853531628645304e-10,-0.1049579707566189469341186813219e-9,-0.1007939580405777662197671844560e-8,-0.9736450265211550752217060575644e-8,-0.9466101382146753511496418259404e-7,-0.9269139557230184816973381023220e-6,-0.9148312509513609827739907179727e-5,-0.9108785014936600974832901874747e-4,-0.9158817613130838424132243282306e-3,-0.9310767919833080253847393582574e-2,0.9041742295948504677274049567506};
	unsigned N=sizeof(tabei)/sizeof(long double);
	for (unsigned i=0;i<N;i++){
	  res *= Z;
	  res += tabei[i];
	}
	return double(res*std::exp(z)/z);
      }
      if (z>=-18 && z<=-10.4){
	// X:=-14.2; evalf(symb2poly(b,h),30)
	long double Z=z+14.2,res=0;
	long double tabei[]={-0.2146565037696152744587246594658e-29,-0.3211304301798548507083223372513e-28,-0.4810676293065620718028423299344e-27,-0.7216859970389437335874447555454e-26,-0.1084277395671112546831370934086e-24,-0.1631614493247996191078695765867e-23,-0.2459334296271721080391164633531e-22,-0.3713482299827011821020880718604e-21,-0.5617631647301656953696071542646e-20,-0.8514935739998273979648471651860e-19,-0.1293355781746933007839041344760e-17,-0.1968884006429264318605910506869e-16,-0.3004345775527340969677228988332e-15,-0.4595950524037171042443121664043e-14,-0.7049704421594200850156215866786e-13,-0.1084472101777504561037738433261e-11,-0.1673430730968429348565586095542e-10,-0.2590825472368919959973457219967e-9,-0.4025494197722650107148738975481e-8,-0.6278735181439544538790456358391e-7,-0.9834030941528033320856991985477e-6,-0.1547201663333347394190105495487e-4,-0.2446158986745476113536605827777e-3,-0.3888037771084034459793304640141e-2,0.9378427721282495585084911135452};
	unsigned N=sizeof(tabei)/sizeof(long double);
	for (unsigned i=0;i<N;i++){
	  res *= Z;
	  res += tabei[i];
	}
	return double(res*std::exp(z)/z);
      }
      if (z>=-28 && z<=-18){
	// X:=-23; evalf(symb2poly(b,h),30)
	long double Z=z+23,res=0;
	long double tabei[]={-0.2146168427075858494404136614850e-34,-0.5148499454995822704300301651629e-33,-0.1236169966295832115713010472045e-31,-0.2970790204951956634218583491140e-30,-0.7146257090942497496926778641796e-29,-0.1720741952629494954062600276217e-27,-0.4147649238148296325602078162016e-26,-0.1000823187587987074490282090878e-24,-0.2417703015904985630926544707734e-23,-0.5847381076193017965725594038649e-22,-0.1415979006710406329134293651920e-20,-0.3433325477996997527365162488189e-19,-0.8336110427404205659872823554255e-18,-0.2026897569209606939865219036707e-16,-0.4935733317936049783947792956182e-15,-0.1203807749521077525537452325013e-13,-0.2940930287406907738614747539418e-12,-0.7197370501359682859646891271814e-11,-0.1764684158114017110045008515596e-9,-0.4335208135324139738484130923065e-8,-0.1067215959531320911498195408991e-6,-0.2632980854624917892922066469020e-5,-0.6511098477481142813714679454527e-4,-0.1614117344014970753628852067852e-2,0.9598801957880143469722276499000};
	unsigned N=sizeof(tabei)/sizeof(long double);
	for (unsigned i=0;i<N;i++){
	  res *= Z;
	  res += tabei[i];
	}
	return double(res*std::exp(z)/z);
      }
      if (z>=-40 && z<=-28){
	// X:=-34; evalf(symb2poly(b,h),30)
	long double Z=z+34,res=0;
	long double tabei[]={-0.1553338170441157171980055301967e-38,-0.6629618584891807480484960239352e-37,-0.2063078177267001621688370383419e-35,-0.7866363363832491801531045634862e-34,-0.2644473347532265615154979481412e-32,-0.9583293463074797094792612009694e-31,-0.3331568830413908327708929748672e-29,-0.1185447458884996131099536033526e-27,-0.4172725550206304813111805693360e-26,-0.1478006485774725675043732714563e-24,-0.5225975670861681162526724856274e-23,-0.1851238198626063152547466792341e-21,-0.6560288456307499178019660543116e-20,-0.2327113628743444723634985600950e-18,-0.8261670399157726998435042888476e-17,-0.2935758096939480075974014556136e-15,-0.1044193898775627537692001439649e-13,-0.3717665504357675255012806230568e-12,-0.1324967081687119761655731034082e-10,-0.4727208926929774019493746342695e-9,-0.1688455962743104637225875947715e-7,-0.6037839002081685958714139659847e-6,-0.2161738491599201737451243010035e-4,-0.7749596600489697701377944551020e-3,0.9721813893840475706338481431853};
	unsigned N=sizeof(tabei)/sizeof(long double);
	for (unsigned i=0;i<N;i++){
	  res *= Z;
	  res += tabei[i];
	}
	return double(res*std::exp(z)/z);
      }
      if (z>=40 || z<=-40){
	long double ei=1,pi=1,Z=z;
	for (long double n=1;;++n){
	  if (pi<1e-16 && pi>-1e-16)
	    break;
	  pi = (n*pi)/Z;
	  ei += pi;
	}
	return double(std::exp(Z)/Z*ei);
      }
      // not used anymore, too slow
      // z<0: int(e^t/t,t,-inf,z)=e^z*int(e^(-u)/(u-z),t,0,inf)
      // z>0: Ei(9.)+int(e^t/t,t,9,z) = Ei(9.)-e^z*int(e^(-u)/(u-z),u,0,z-9)
      double nstep=400,a=0,fz=0; 
      for (;nstep>0.25;nstep*=exp_minus_1_over_4){
	double Fz=0;
	int N=int(nstep+.5);
	if (N<1)
	  N=1;
	double taille=1.0;
	if (z>0 && a+1>z-9)
	  taille=(z-9)-a;
	// Simpson over [a,a+taille]
	double t=a,tmp,Ninv=taille/N;
	t = a+Ninv/2.;
	double expt=std::exp(-t),expfact=std::exp(-Ninv);
	for (int i=0;i<N;++i){ // middle points
	  tmp = expt/(t-z); 
	  Fz += tmp;
	  expt *= expfact;
	  t += Ninv;
	}
	Fz *= 2; 
	t = a+Ninv;
	expt=std::exp(-t);
	for (int i=1;i<N;++i){ 
	  tmp = expt/(t-z); // endpoint
	  Fz += tmp;
	  expt *= expfact;
	  t += Ninv;
	}
	Fz *= 2; 
	tmp=std::exp(-a)/(a-z); // endpoint
	Fz += tmp;
	a += taille;
	tmp=std::exp(-a)/(a-z); // endpoint
	Fz += tmp;
	fz += Fz*taille/(6*N);
	if (z>0 && a>=z-9)
	  break;
      }
      fz *= std::exp(z);
      if (z<0)
	return fz;
      return 1037.878290717090-fz;
    }
    // find number of digits that must be added to prec
    // n^n/n! equivalent to e^n*sqrt(2*pi*n)
    // Note that Ei(z) might be as small as exp(-z) for relative prec
    int newprec,nbitsz=int(z._DOUBLE_val/std::log(2.)),prec2=int(prec*std::log(10.0)/std::log(2.0)+.5)+int(z._DOUBLE_val/std::log(2.0));
    if (nbitsz>prec2){ 
      // use asymptotic expansion at z=inf
      z = accurate_evalf(args,prec2);
      gen ei=1,pi=1,eps=accurate_evalf(pow(10,-prec,contextptr)*exp(-abs(z,contextptr),contextptr),prec2)/2;
      for (int n=1;;++n){
	if (is_greater(eps,abs(pi,contextptr),contextptr))
	  break;
	pi = (n*pi)/z;
	ei += pi;
      }
      return exp(z,contextptr)/z*ei;
    }
    // use series expansion at z=0
    if (z._DOUBLE_val>1)
      newprec = prec2+nbitsz+int(std::log(z._DOUBLE_val)/2)+2;
    else
      newprec = prec2+2;
    z = accurate_evalf(args,newprec);
    gen ei=0,pi=1,eps=accurate_evalf(pow(10,-prec,contextptr)*exp(-abs(z,contextptr),contextptr),newprec)/2;
    for (int n=1;;n++){
      pi = pi*z/n;
      if (is_greater(eps,abs(pi,contextptr),contextptr))
	break;
      ei += pi/n;
    }
    if (is_positive(z,contextptr))
      ei=ei+ln(z,contextptr)+m_gamma(newprec);
    else
      ei=ei+ln(-z,contextptr)+m_gamma(newprec);
    ei = accurate_evalf(ei,prec2);
    return ei;
  }
  gen Ei(const gen & args,int n,GIAC_CONTEXT){
    if (n==1)
      return -Ei(-args,contextptr);
    if (n<2)
      setdimerr();
    if (is_zero(args)){
      if (n==1)
	return plus_inf;
      return plus_one/gen(n-1);
    }
    if (args==plus_inf)
      return 0;
    if (args==minus_inf)
      return minus_inf;
    if (is_inf(args)|| is_undef(args))
      return undef;
    return (exp(-args,contextptr)-args*Ei(args,n-1,contextptr))/gen(n-1);
  }
  gen _Ei(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT){
      return Ei(args,contextptr);
    }
    if ( args._VECTptr->size()!=2 ){
      return symbolic(at_Ei,args);
    }
    gen x(args._VECTptr->front()),n(args._VECTptr->back());
    if (n.type==_REAL)
      n=n.evalf_double(1,contextptr);
    if (n.type==_DOUBLE_)
      n=int(n._DOUBLE_val);
    if (n.type!=_INT_)
      setsizeerr();
    if (n==1)
      *logptr(contextptr) << "Warning, Ei(x,1) is defined as -Ei(-x), not as Ei(x)" << endl;
    return Ei(x,n.val,contextptr);
  }
  const string _Ei_s("Ei");
  unary_function_eval __Ei(&_Ei,&D_at_Ei,&taylor_Ei,_Ei_s);
  unary_function_ptr at_Ei (&__Ei,0,true);

  gen d_Ei0(const gen & args,GIAC_CONTEXT){
    return (exp(args,contextptr)-1)/args;
  }
  partial_derivative_onearg D_at_Ei0(&d_Ei0);
  gen taylor_Ei0(const gen & lim_point,const int ordre,const unary_function_ptr & f, int direction,gen & shift_coeff,GIAC_CONTEXT){
    if (ordre<0)
      return 0;
    if (!is_zero(lim_point))
      return taylor(lim_point,ordre,f,0,shift_coeff,contextptr);
    shift_coeff=1;
    // sum( 1/(k)/(k)! * x^(k) )
    gen tmp(1);
    vecteur v;
    for (int n=0;n<=ordre;){
      n++;
      tmp=n*tmp;
      v.push_back(inv(n*tmp,contextptr));
    }
    v.push_back(undef);
    return v;
  }
  gen _Ei0(const gen & args,GIAC_CONTEXT){
    if (is_zero(args))
      return 0;
    if (is_undef(args))
      return args;
    if (is_inf(args))
      return minus_inf;
    if (args.type!=_DOUBLE_ && args.type!=_REAL && args.type!=_CPLX)
      return symbolic(at_Ei0,args);
    gen si,ci;
    sici(args,si,ci,decimal_digits(contextptr),2,contextptr);
    return ci-evalf(cst_euler_gamma,1,contextptr)-ln(args,contextptr);
  }
  const string _Ei0_s("Ei0");
  unary_function_eval __Ei0(&_Ei0,&D_at_Ei0,&taylor_Ei0,_Ei0_s);
  unary_function_ptr at_Ei0 (&__Ei0,0,true); /* FIXME should not registered */

  gen _Dirac(const gen & args,GIAC_CONTEXT){
    if (args.type==_VECT && args.subtype!=_SEQ__VECT)
      return apply(args,_Dirac,contextptr);
    gen f=args;
    if (args.type==_VECT && args.subtype==_SEQ__VECT && !args._VECTptr->empty())
      f=args._VECTptr->front();
    if (is_zero(f))
      return unsigned_inf;
    if (f.type<_IDNT)
      return 0;
    return symbolic(at_Dirac,args);
  }
  gen d_Dirac(const gen & args,GIAC_CONTEXT){
    vecteur v(gen2vecteur(args));
    if (v.size()==1)
      v.push_back(0);
    if (v.size()!=2 || v.back().type!=_INT_)
      setdimerr();
    return _Dirac(gen(makevecteur(v.front(),v.back().val+1),_SEQ__VECT),contextptr);
  }
  partial_derivative_onearg D_at_Dirac(&d_Dirac);
  const string _Dirac_s("Dirac");
  unary_function_eval __Dirac(&_Dirac,&D_at_Dirac,_Dirac_s);
  unary_function_ptr at_Dirac (&__Dirac,0,true);
  partial_derivative_onearg D_Heaviside(&_Dirac);

  gen _Heaviside(const gen & args,GIAC_CONTEXT){
    if (args.type==_VECT)
      return apply(args,_Heaviside,contextptr);
    if (is_zero(args))
      return plus_one;
    gen tmp=_sign(args,contextptr);
    if (tmp.type<=_DOUBLE_)
      return (tmp+1)/2;
    return symbolic(at_Heaviside,args);
  }
  const string _Heaviside_s("Heaviside");
  unary_function_eval __Heaviside(&_Heaviside,&D_Heaviside,_Heaviside_s);
  unary_function_ptr at_Heaviside (&__Heaviside,0,true);

  // if you add functions to solve_fcns, modify the second argument of solve_fcns_v to reflect the number of functions in the array
  unary_function_ptr solve_fcns[]={ at_exp,at_ln,at_sin,at_cos,at_tan,at_asin,at_acos,at_atan,at_sinh,at_cosh,at_tanh,at_asinh,at_acosh,at_atanh};
  vector<unary_function_ptr > solve_fcns_v(solve_fcns,solve_fcns+sizeof(solve_fcns)/sizeof(unary_function_ptr));

  // #ifndef GNUWINCE
  // #ifndef WIN32
  // #endif
  gen cst_two_pi(symbolic(at_prod,makevecteur(2,cst_pi)));
  gen cst_pi_over_2(_FRAC2_SYMB(cst_pi,2));
  gen plus_inf(symbolic(at_plus,_IDNT_infinity()));
  gen minus_inf(symbolic(at_neg,_IDNT_infinity()));
  // #ifndef GNUWINCE
  // #endif
  // #endif
  gen plus_one_half(fraction(1,2));
  gen minus_one_half(symbolic(at_neg,symb_inv(2)));
  gen plus_sqrt3(sqrt(3,context0)); // ok
  gen plus_sqrt2(sqrt(2,context0)); // ok
  gen plus_sqrt6(sqrt(6,context0)); // ok
  gen minus_sqrt6(symbolic(at_neg,plus_sqrt6));
  gen minus_sqrt3(symbolic(at_neg,plus_sqrt3));
  gen minus_sqrt2(symbolic(at_neg,plus_sqrt2));
  gen minus_sqrt3_2(_FRAC2_SYMB(minus_sqrt3,2));
  gen minus_sqrt2_2(_FRAC2_SYMB(minus_sqrt2,2));
  gen minus_sqrt3_3(_FRAC2_SYMB(minus_sqrt3,3));
  gen plus_sqrt3_2(_FRAC2_SYMB(plus_sqrt3,2));
  gen plus_sqrt2_2(_FRAC2_SYMB(plus_sqrt2,2));
  gen plus_sqrt3_3(_FRAC2_SYMB(plus_sqrt3,3));
  gen cos_pi_12(_FRAC2_SYMB(
			    symbolic(at_plus,gen(makevecteur(plus_sqrt6,plus_sqrt2),_SEQ__VECT)),
			    2));
  gen minus_cos_pi_12(_FRAC2_SYMB(
			    symbolic(at_plus,gen(makevecteur(minus_sqrt6,minus_sqrt2),_SEQ__VECT)),
			    2));
  gen sin_pi_12(_FRAC2_SYMB(
			    symbolic(at_plus,gen(makevecteur(plus_sqrt6,minus_sqrt2),_SEQ__VECT)),
			    2));
  gen minus_sin_pi_12(_FRAC2_SYMB(
			    symbolic(at_plus,gen(makevecteur(plus_sqrt2,minus_sqrt6),_SEQ__VECT)),
			    2));
  gen tan_pi_12(symbolic(at_plus,makevecteur(2,minus_sqrt3)));
  gen minus_tan_pi_12(symbolic(at_neg,tan_pi_12));
  gen tan_5pi_12(symbolic(at_plus,makevecteur(2,plus_sqrt3)));
  gen minus_tan_5pi_12(symbolic(at_neg,tan_5pi_12));
  gen rad2deg_e(_FRAC2_SYMB(180,cst_pi));
  gen deg2rad_e(_FRAC2_SYMB(cst_pi,180));
  unary_function_ptr reim_op[]={at_inv,at_exp,at_cos,at_sin,at_tan,at_cosh,at_sinh,at_tanh,at_atan,at_lnGamma_minus,at_Gamma,at_Psi_minus_ln,at_Psi,at_Zeta,at_Eta,at_sign};

  // for subst.cc
  unary_function_ptr sincostan_tab[]={at_sin,at_cos,at_tan};
  unary_function_ptr asinacosatan_tab[]={at_asin,at_acos,at_atan};
  unary_function_ptr sinhcoshtanh_tab[]={at_sinh,at_cosh,at_tanh};
  vector<unary_function_ptr> sincostan_v(sincostan_tab,sincostan_tab+3);
  vector<unary_function_ptr> asinacosatan_v(asinacosatan_tab,asinacosatan_tab+3);
  vector<unary_function_ptr> sinhcoshtanh_v(sinhcoshtanh_tab,sinhcoshtanh_tab+3);
  vector <unary_function_ptr> sincostansinhcoshtanh_v(merge(sincostan_v,sinhcoshtanh_v));
  unary_function_ptr sign_floor_ceil_round_tab[]={at_sign,at_floor,at_ceil,at_round};
  vector<unary_function_ptr> sign_floor_ceil_round_v(sign_floor_ceil_round_tab,sign_floor_ceil_round_tab+4);
  vector<unary_function_ptr> exp_v(1,at_exp);
  vector<unary_function_ptr> tan_v(1,at_tan);
  vector<unary_function_ptr> asin_v(1,at_asin);
  vector<unary_function_ptr> acos_v(1,at_acos);
  vector<unary_function_ptr> atan_v(1,at_atan);
  vector<unary_function_ptr> pow_v(1,at_pow);
  gen_op_context halftan_tab[]={sin2tan2,cos2tan2,tan2tan2};
  gen_op_context hyp2exp_tab[]={sinh2exp,cosh2exp,tanh2exp};
  gen_op_context trig2exp_tab[]={sin2exp,cos2exp,tan2exp};
  gen_op_context atrig2ln_tab[]={asin2ln,acos2ln,atan2ln};
  vector< gen_op_context > halftan_v(halftan_tab,halftan_tab+3);
  vector< gen_op_context > hyp2exp_v(hyp2exp_tab,hyp2exp_tab+3);
  vector< gen_op_context > trig2exp_v(trig2exp_tab,trig2exp_tab+3);
  vector< gen_op_context > halftan_hyp2exp_v(merge(halftan_v,hyp2exp_v));
  vector< gen_op_context > exp2sincos_v(1,exp2sincos);
  vector< gen_op_context > tan2sincos_v(1,tantosincos);
  vector< gen_op_context > tan2sincos2_v(1,tantosincos2);
  vector< gen_op_context > tan2cossin2_v(1,tantocossin2);
  vector< gen_op_context > asin2acos_v(1,asintoacos);
  vector< gen_op_context > asin2atan_v(1,asintoatan);
  vector< gen_op_context > acos2asin_v(1,acostoasin);
  vector< gen_op_context > acos2atan_v(1,acostoatan);
  vector< gen_op_context > atan2asin_v(1,atantoasin);
  vector< gen_op_context > atan2acos_v(1,atantoacos);
  vector< gen_op_context > atrig2ln_v(atrig2ln_tab,atrig2ln_tab+3);
  vector< gen_op_context > trigcos_v(1,trigcospow);
  vector< gen_op_context > trigsin_v(1,trigsinpow);
  vector< gen_op_context > trigtan_v(1,trigtanpow);
  vector< gen_op_context > powexpand_v(1,powtopowexpand);
  vector< gen_op_context > exp2power_v(1,exptopower);
  vector< unary_function_ptr > gamma_v(1,at_Gamma);
  vector< gen_op_context > gamma2factorial_v(1,gammatofactorial);
  vector< unary_function_ptr > factorial_v(1,at_factorial);
  vector< gen_op_context > factorial2gamma_v(1,factorialtogamma);

  // for integration
  unary_function_ptr primitive_tab_op[]={at_sin,at_cos,at_tan,at_exp,at_sinh,at_cosh,at_tanh,at_asin,at_acos,at_atan,at_ln,0};
  unary_function_ptr inverse_tab_op[]={at_asin,at_acos,at_atan,at_ln,at_asinh,at_acos,at_atanh,at_erf,at_erfc,at_Ei,at_Si,at_Ci,0};

  unary_function_ptr analytic_sommets[]={at_plus,at_prod,at_neg,at_inv,at_pow,at_sin,at_cos,at_tan,at_exp,at_sinh,at_cosh,at_tanh,at_asin,at_acos,at_atan,at_asinh,at_atanh,at_acosh,at_ln,at_sqrt,0};

  // test if g is < > <= >=, 
  unary_function_ptr inequality_tab[]={at_equal,at_inferieur_strict,at_inferieur_egal,at_different,at_superieur_strict,at_superieur_egal};
  vector<unary_function_ptr> inequality_sommets(inequality_tab,inequality_tab+sizeof(inequality_tab)/sizeof(unary_function_ptr));
  int is_inequality(const gen & g){
    if (g.type!=_SYMB)
      return false;
    return equalposcomp(inequality_sommets,g._SYMBptr->sommet);
  }


  string unquote(const string & s){
    int l=s.size();
    if (l>2 && s[0]=='"' && s[l-1]=='"')
      return s.substr(1,l-2);
    else
      return s;
  }


#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC
