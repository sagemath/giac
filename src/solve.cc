// -*- mode:C++ ; compile-command: "g++-3.4 -I.. -g -c solve.cc" -*-
#include "first.h"
/*
 *  Copyright (C) 2001,7 B. Parisse, R. De Graeve
 *  Institut Fourier, 38402 St Martin d'Heres
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
#include "gen.h"
#include "solve.h"
#include "modpoly.h"
#include "unary.h"
#include "symbolic.h"
#include "usual.h"
#include "sym2poly.h"
#include "subst.h"
#include "derive.h"
#include "plot.h"
#include "prog.h"
#include "series.h"
#include "alg_ext.h"
#include "intg.h"
#include "rpn.h"
#include "lin.h"
#include "misc.h"
#include "cocoa.h"
#include "ti89.h"
#ifdef HAVE_LIBGSL
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#endif

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC


  int intvar_counter=0;
  int realvar_counter=0;
  string print_intvar_counter(){
    if (intvar_counter<0)
      return print_INT_(-intvar_counter);
    string res=print_INT_(intvar_counter);
    ++intvar_counter;
    return res;
  }

  string print_realvar_counter(){
    if (realvar_counter<0)
      return print_INT_(int(-realvar_counter));
    string res=print_INT_(int(realvar_counter));
    ++realvar_counter;
    return res;
  }

  gen _reset_solve_counter(const gen & args,const context * contextptr){
    if (is_zero(args)){
      intvar_counter=0;
      return 1;
    }
    if (is_one(args)){
      realvar_counter=0;
      return 1;
    }
    if (args.type==_VECT && args._VECTptr->size()==2){
      intvar_counter=int(evalf_double(args._VECTptr->front(),1,contextptr)._DOUBLE_val);
      realvar_counter=int(evalf_double(args._VECTptr->back(),1,contextptr)._DOUBLE_val);      
    }
    else {
      intvar_counter=0;
      realvar_counter=0;
    }
    return 1;
  }
  const string _reset_solve_counter_s("reset_solve_counter");
  unary_function_eval ___reset_solve_counter(&_reset_solve_counter,_reset_solve_counter_s);
  unary_function_ptr at_reset_solve_counter (&___reset_solve_counter,0,true);

  void set_merge(vecteur & v,const vecteur & w){
    const_iterateur it=w.begin(),itend=w.end();
    for (;it!=itend;++it)
      if (!equalposcomp(v,*it))
	v.push_back(*it);
  }

  gen one_tour(GIAC_CONTEXT){
    if (angle_radian(contextptr)) 
      return cst_two_pi;
    else
      return 360;
  }
  gen one_half_tour(GIAC_CONTEXT){
    if (angle_radian(contextptr)) 
      return cst_pi;
    else
      return 180;
  }
  gen isolate_exp(const gen & e,int isolate_mode,GIAC_CONTEXT){
    if (isolate_mode &1)
      return ln(e,contextptr);
    if (e.type!=_VECT){
      if (is_strictly_positive(-e,contextptr))
	return vecteur(0);
      else
	return ln(e,contextptr);
    }
    // check in real mode for negative ln
    const_iterateur it=e._VECTptr->begin(),itend=e._VECTptr->end();
    vecteur res;
    for (;it!=itend;++it){
      if (!is_strictly_positive(-*it,contextptr))
	res.push_back(ln(*it,contextptr));
    }
    return res;
  }
  gen isolate_ln(const gen & e,int isolate_mode,GIAC_CONTEXT){
    return simplify(exp(e,contextptr),contextptr);
  }
  gen isolate_sin(const gen & e,int isolate_mode,GIAC_CONTEXT){
    gen asine=asin(e,contextptr);
    if (!(isolate_mode & 2))
      return makevecteur(asine,one_half_tour(contextptr)-asine);
    identificateur * x=new identificateur(string("n_")+print_intvar_counter());
    return makevecteur(asine+(*x)*one_tour(contextptr),one_half_tour(contextptr)-asine+(*x)*one_tour(contextptr));
  }
  gen isolate_cos(const gen & e,int isolate_mode,GIAC_CONTEXT){
    gen acose=acos(e,contextptr);
    if (!(isolate_mode & 2))
      return makevecteur(acose,-acose);
    identificateur * x=new identificateur(string("n_")+print_intvar_counter());
    return makevecteur(acose+(*x)*one_tour(contextptr),-acose+(*x)*one_tour(contextptr));
  }
  gen isolate_tan(const gen & e,int isolate_mode,GIAC_CONTEXT){
    if (!(isolate_mode & 2))
      return atan(e,contextptr);
    identificateur * x=new identificateur(string("n_")+print_intvar_counter());
    return atan(e,contextptr)+(*x)*one_half_tour(contextptr);
  }
  gen isolate_asin(const gen & e,int isolate_mode,GIAC_CONTEXT){
    return sin(e,contextptr);
  }
  gen isolate_acos(const gen & e,int isolate_mode,GIAC_CONTEXT){
    return cos(e,contextptr);
  }
  gen isolate_atan(const gen & e,int isolate_mode,GIAC_CONTEXT){
    return tan(e,contextptr);
  }

  gen isolate_asinh(const gen & e,int isolate_mode,GIAC_CONTEXT){
    return sinh(e,contextptr);
  }
  gen isolate_acosh(const gen & e,int isolate_mode,GIAC_CONTEXT){
    return cosh(e,contextptr);
  }
  gen isolate_atanh(const gen & e,int isolate_mode,GIAC_CONTEXT){
    return tanh(e,contextptr);
  }

  gen isolate_sinh(const gen & e,int isolate_mode,GIAC_CONTEXT){
    gen asine= asinh(e,contextptr);
    if (!(isolate_mode & 2))
      return asine;
    identificateur * x=new identificateur(string("n_")+print_intvar_counter());
    return makevecteur(asine+(*x)*one_tour(contextptr)*cst_i,(one_half_tour(contextptr)+(*x)*one_tour(contextptr))*cst_i-asine);
  }
  gen isolate_cosh(const gen & e,int isolate_mode,GIAC_CONTEXT){
    gen acose=acosh(e,contextptr);
    if (!(isolate_mode & 2))
      return makevecteur(acose,-acose);
    identificateur * x=new identificateur(string("n_")+print_intvar_counter());
    return makevecteur(acose+(*x)*one_tour(contextptr)*cst_i,-acose+(*x)*one_tour(contextptr)*cst_i);
  }
  gen isolate_tanh(const gen & e,int isolate_mode,GIAC_CONTEXT){
    if (!(isolate_mode & 2))
      return atanh(e,contextptr);
    identificateur * x=new identificateur(string("n_")+print_intvar_counter());
    return atanh(e,contextptr)+(*x)*one_half_tour(contextptr)*cst_i;
  }

  gen (* isolate_fcns[] ) (const gen &,int,GIAC_CONTEXT) = { isolate_exp,isolate_ln,isolate_sin,isolate_cos,isolate_tan,isolate_asin,isolate_acos,isolate_atan,isolate_sinh,isolate_cosh,isolate_tanh,isolate_asinh,isolate_acosh,isolate_atanh};

  vecteur find_excluded(const gen & g,GIAC_CONTEXT){
    if (g.type!=_IDNT)
      return vecteur(0);
    gen g2=g._IDNTptr->eval(eval_level(contextptr),g,contextptr);
    if ((g2.type==_VECT) && (g2.subtype==_ASSUME__VECT)){
      vecteur v=*g2._VECTptr;
      if ( (v.size()==3) && (v.front()==_DOUBLE_) && (v[2].type==_VECT)){
	return *v[2]._VECTptr;;
      }
    }
    return vecteur(0);
  }

  vecteur protect_sort(const vecteur & res,GIAC_CONTEXT){
    try {
      gen tmp=_sort(res,contextptr);
      if (tmp.type==_VECT){
	vecteur w=*tmp._VECTptr,res;
	iterateur it=w.begin(),itend=w.end();
	for (;it!=itend;++it){
	  if (res.empty() || *it!=res.back())
	    res.push_back(*it);
	}
	return res;
      }
    }
    catch (std::runtime_error & e){
      cerr << e.what() << endl;
    }
    return res;
  }

  vecteur find_singularities(const gen & e,const identificateur & x,int cplxmode,GIAC_CONTEXT){
    vecteur lv(lvarxpow(e,x));
    vecteur res;
    vecteur l(lvar(e));
    gen p=e2r(e,l,contextptr),n,d;
    vecteur pv=gen2vecteur(p);
    const_iterateur jt=pv.begin(),jtend=pv.end();
    for (;jt!=jtend;++jt){
      fxnd(*jt,n,d);
      if (d.type==_POLY)
	res=solve(r2e(d,l,contextptr),x,cplxmode,contextptr);
    }
    const_iterateur it=lv.begin(),itend=lv.end();
    for (;it!=itend;++it){
      if (it->type!=_SYMB)
	continue;
      unary_function_ptr & u=it->_SYMBptr->sommet;
      gen & f=it->_SYMBptr->feuille;
      res=mergevecteur(res,find_singularities(f,x,cplxmode,contextptr));
      if (u==at_ln)
	res=mergevecteur(res,solve(f,x,cplxmode,contextptr));
      if (u==at_pow && f.type==_VECT && f._VECTptr->size()==2)
	res=mergevecteur(res,solve(f._VECTptr->front(),x,cplxmode,contextptr));	
      if (u==at_tan)
	res=mergevecteur(res,solve(cos(f,contextptr),x,cplxmode,contextptr));
    }
    if (cplxmode)
      return res;
    return protect_sort(res,contextptr);
  }

  void solve_ckrange(const identificateur & x,vecteur & v,int isolate_mode,GIAC_CONTEXT){
    vecteur w,excluded(find_excluded(x,contextptr));
    find_range(x,w,contextptr);
    if (w.size()!=1 || w.front().type!=_VECT)
      return;
    w=*w.front()._VECTptr;
    if (w.size()!=2)
      return;
    gen l(w.front()),m(w.back());
    vecteur newv;
    iterateur it=v.begin(),itend=v.end();
    for (;it!=itend;++it){
      *it=simplifier(*it,contextptr);
      if (equalposcomp(excluded,*it))
	continue;
      gen sol=*it;
      if (l!=minus_inf && sign(l-sol,contextptr)==1)
	continue;
      if (m!=plus_inf && sign(sol-m,contextptr)==1)
	continue;
      sol=evalf(sol,eval_level(contextptr),contextptr);
      if (sol.type==_CPLX && !(isolate_mode &1) && !is_zero(im(sol,contextptr)))
	continue;
      if (sol.type!=_DOUBLE_){ // check for trig solutions
	newv.push_back(*it);
	vecteur lv(lvar(sol));
	if (lv.size()!=1 || l==minus_inf || m==plus_inf)
	  continue;
	gen n(lv.front()),a,b;
	// check linearity
	if (n.type!=_IDNT || n.print(contextptr).substr(0,2)!="n_" || !is_linear_wrt(*it,*n._IDNTptr,a,b,contextptr))
	  continue;
	newv.pop_back();
	a=normal(a,contextptr);
	b=normal(b,contextptr);
	if (!is_positive(a,contextptr))
	  std::swap<gen>(l,m);
	int n0(_ceil(evalf((l-b)/a,eval_level(contextptr),contextptr),contextptr).val),n1(_floor(evalf((m-b)/a,eval_level(contextptr),contextptr),contextptr).val);
	for (;n0<=n1;++n0)
	  newv.push_back(n0*a+b);
      }
      else {
	if (is_strictly_greater(l,sol,contextptr))
	  continue;
	if (is_strictly_greater(sol,m,contextptr))
	  continue;
	newv.push_back(*it);
      }
    }
    v=newv;
  }

  // Helper for the solver, make a translation using x^(n-1) coeff
  // and find gcd of deg, return true if non-trivial gcd found
  bool translate_gcddeg(const vecteur & v,vecteur & v_translated, gen & x_translation,int & gcddeg){
    int s=v.size();
    if (s<4)
      return false;
    x_translation=-v[1]/((s-1)*v[0]);
    v_translated=taylor(v,x_translation,0);
    gcddeg=0;
    for (int i=1;i<s;++i){
      if (!is_zero(v_translated[i]))
	gcddeg=gcd(gcddeg,i);
    }
    if (gcddeg<=1)
      return false;
    int newdeg=(s-1)/gcddeg+1;
    // compress v_translated, keep only terms with index multiple of gcddeg
    for (int i=1;i<newdeg;++i){
      v_translated[i]=v_translated[i*gcddeg];
    }
    v_translated=vecteur(v_translated.begin(),v_translated.begin()+newdeg);
    return true;
  }

  vecteur solve_inequation(const gen & e0,const identificateur & x,int direction,GIAC_CONTEXT);

  vecteur solve_piecewise(const gen & args,const gen & value,const identificateur & x,int isolate_mode,GIAC_CONTEXT){
    if (args.type!=_VECT)
      setsizeerr();
    vecteur & piece_args=*args._VECTptr;
    vecteur failtest; // all these tests must fail to keep solution
    gen successtest,equation; // this test must succeed
    int s=piece_args.size();
    vecteur res;
    for (int i=0;i<s;i+=2){
      if (i)
	failtest.push_back(successtest);
      if (i+1==s){
	successtest=1;
	equation=piece_args[i];
      }
      else {
	successtest=piece_args[i];
	equation=piece_args[i+1];
      }
      int fails=failtest.size();
      vecteur sol=solve(equation-value,x,isolate_mode,contextptr);
      // now test whether solutions in sol are acceptable
      const_iterateur it=sol.begin(),itend=sol.end();
      for (;it!=itend;++it){
	const gen & g=*it;
	if (g==x){
	  if (fails){
	    gen tmp=symb_not(symbolic(at_ou,gen(failtest,_SEQ__VECT)));
	    res.push_back(is_one(successtest)?tmp:symb_and(tmp,successtest));
	  }
	  else
	    res.push_back(is_one(successtest)?g:successtest);
	  continue;
	}
	if (!is_zero(derive(g,x,contextptr))){
	  if (fails){
	    gen tmp=symb_not(symbolic(at_ou,gen(failtest,_SEQ__VECT)));
	    tmp=is_one(successtest)?tmp:symb_and(tmp,successtest);
	    res.push_back(symb_and(tmp,g));
	  }
	  else
	    res.push_back(is_one(successtest)?g:symb_and(successtest,g));
	  continue;
	}
	int j;
	for (j=0;j<fails;++j){
	  if (is_one(subst(failtest[j],x,g,false,contextptr)))
	    break;
	}
	if (j==fails && is_one(subst(successtest,x,g,false,contextptr)))
	  res.push_back(g);
      }
    }
    return res;
  }

  // inner solver
  void in_solve(const gen & e,const identificateur & x,vecteur &v,int isolate_mode,GIAC_CONTEXT){
    bool complexmode=isolate_mode & 1;
    vecteur lv(lvarx(e,x));
    int s=lv.size();
    if (!s)
      return;
    if (s>1){
      for (int i=0;i<s;++i){
	gen xvar=lv[i];
	if (xvar._SYMBptr->sommet==at_sign){
	  gen new_e=subst(e,xvar,1,false,contextptr);
	  vecteur vplus;
	  in_solve(new_e,x,vplus,isolate_mode,contextptr);
	  const_iterateur it=vplus.begin(),itend=vplus.end();
	  for (;it!=itend;++it){
	    if (is_one(subst(xvar,x,*it,false,contextptr)))
	      v.push_back(*it);
	  }
	  new_e=subst(e,xvar,-1,false,contextptr);
	  vecteur vminus;
	  in_solve(new_e,x,vminus,isolate_mode,contextptr);
	  it=vminus.begin();
	  itend=vminus.end();
	  for (;it!=itend;++it){
	    if (is_one(-subst(xvar,x,*it,false,contextptr)))
	      v.push_back(*it);
	  }
	  return;
	}
      }
      throw(std::runtime_error("Unable to isolate "+x.print(contextptr)+" in "+e.print(contextptr)));
    }
    gen xvar(lv.front());
    if (xvar!=x){ // xvar must be a unary function of x
      if (xvar.type!=_SYMB)
	settypeerr();
      if (xvar._SYMBptr->sommet!=at_piecewise && xvar._SYMBptr->feuille.type==_VECT)
	throw(std::runtime_error("Unable to isolate "+x.print(contextptr)+" in "+xvar.print(contextptr)));
      if (xvar._SYMBptr->sommet==at_sign){
	gen new_e=subst(e,xvar,1,false,contextptr);
	if (is_zero(new_e)){
	  v=solve_inequation(symbolic(at_superieur_strict,makevecteur(xvar._SYMBptr->feuille,0)),x,1,contextptr);
	}
	else {
	  new_e=subst(e,xvar,-1,false,contextptr);
	  if (is_zero(new_e)){
	    v=solve_inequation(symbolic(at_inferieur_strict,makevecteur(xvar._SYMBptr->feuille,0)),x,-1,contextptr);
	  }
	}
	return;
      }
      int pos=equalposcomp(solve_fcns_v,xvar._SYMBptr->sommet);
      if (xvar._SYMBptr->sommet==at_piecewise)
	pos=-1;
      if (!pos)
	throw(std::runtime_error("Unable to isolate function "+xvar._SYMBptr->sommet.ptr->print(contextptr)));
      // solve with respect to xvar
      identificateur localt(" t");
      // ck_parameter_t();
      gen new_e=subst(e,xvar,localt,false,contextptr);
      vecteur new_v=solve(new_e,localt,isolate_mode,contextptr);
      const_iterateur it=new_v.begin(),itend=new_v.end();
      for (;it!=itend;++it){
	if (pos==-1){
	  // solve piecewise()==*it
	  set_merge(v,solve_piecewise(xvar._SYMBptr->feuille,*it,x,isolate_mode,contextptr));
	  continue;
	}
	gen res=isolate_fcns[pos-1](*it,isolate_mode,contextptr);
	if (res.type!=_VECT)
	  set_merge(v,solve(xvar._SYMBptr->feuille-res,x,isolate_mode,contextptr));
	else {
	  const_iterateur it=res._VECTptr->begin(),itend=res._VECTptr->end();
	  for (;it!=itend;++it)
	    set_merge(v,solve(xvar._SYMBptr->feuille-*it,x,isolate_mode,contextptr));
	}
      }
      solve_ckrange(x,v,isolate_mode,contextptr);
      return;
    } // end xvar!=x
    // rewrite e as a univariate polynomial, first add other vars to x
    vecteur newv;
    lv=vecteur(1,vecteur(1,x));
    alg_lvar(e,lv);
    vecteur lvrat(1,x);
    lvar(e,lvrat);
    if (lvrat==lv.front())
      lv=lvrat;
    // int lv_size=lv.size();
    gen num,den,f;
    f=e2r(e,lv,contextptr);
    fxnd(f,num,den);
    if (num.type!=_POLY)
      return;
    vecteur w=polynome2poly1(*num._POLYptr,1);
    lv.erase(lv.begin()); // remove x from lv (CDR_VECT)
    int deg;
    vecteur w_translated;
    gen delta_x;
    if (translate_gcddeg(w,w_translated,delta_x,deg)){
      // composite polynomials
      gen invdeg=inv(deg,contextptr);
      gen newe=symb_horner(*r2sym(w_translated,lv,contextptr)._VECTptr,x);
      delta_x=-r2sym(delta_x,lv,contextptr);
      vecteur vtmp;
      in_solve(newe,x,vtmp,isolate_mode,contextptr);
      vecteur unitroot(1,plus_one),munitroot;
      if (complexmode){
	for (int k=1;k<deg;++k)
	  unitroot.push_back(exp(2*k*cst_pi/deg*cst_i,contextptr));
	for (int k=0;k<deg;++k)
	  munitroot.push_back(exp((1+2*k)*cst_pi/deg*cst_i,contextptr));
      }
      const_iterateur it=vtmp.begin(),itend=vtmp.end();
      for (;it!=itend;++it){
	bool negatif=is_strictly_positive(-*it,contextptr);
	gen tmp=pow((negatif?-*it:*it),invdeg,contextptr);
	if (complexmode){
	  const_iterateur jt,jtend;
	  if (!negatif){
	    jt=unitroot.begin();
	    jtend=unitroot.end();
	  }
	  else {
	    jt=munitroot.begin();
	    jtend=munitroot.end();
	  }
	  for (;jt!=jtend;++jt)
	    newv.push_back(delta_x + (*jt) * tmp);
	}
	else {
	  if (deg%2)
	    newv.push_back(delta_x + (negatif?-tmp:tmp));
	  else {
	    if (!negatif){
	      newv.push_back(delta_x + tmp);
	      newv.push_back(delta_x - tmp);
	    }
	  }
	}
      }
      solve_ckrange(x,newv,isolate_mode,contextptr);
      v=mergevecteur(v,newv);
      return;
    }
    // if degree(w)=0, 1 or 2 solve it, otherwise error (should return ext)
    int d=w.size()-1;
    if (!d)
      return;
    if (d==1){
      gen tmp=rdiv(-r2sym(w.back(),lv,contextptr),r2sym(w.front(),lv,contextptr));
      if (!complexmode && has_i(tmp))
	return;
      newv.push_back(tmp);
      solve_ckrange(x,newv,isolate_mode,contextptr);
      v=mergevecteur(v,newv);
      return;
    }
    if (d>2){
      if (has_num_coeff(w)){
	if (complexmode)
	  newv=proot(w,epsilon(contextptr));
	else
	  newv=real_proot(w,epsilon(contextptr),contextptr);
	solve_ckrange(x,newv,isolate_mode,contextptr);
	v=mergevecteur(v,newv);
	return;
      }
      int n=is_cyclotomic(w,epsilon(contextptr));
      if (!n){
	*logptr(contextptr) << "Warning! Algebraic extension not implemented yet for poly " << r2sym(w,lv,contextptr) << endl;
	if (has_num_coeff(evalf(w,1,contextptr))){
	  if (complexmode)
	    newv=proot(w,epsilon(contextptr));
	  else
	    newv=real_proot(w,epsilon(contextptr),contextptr);
	  solve_ckrange(x,newv,isolate_mode,contextptr);
	  v=mergevecteur(v,newv);
	  return;
	}
	return;
      }
      if (complexmode){
	for (int j=1;j<=n/2;++j){
	  if (gcd(j,n)==1){
	    if (n%2){
	      newv.push_back(exp(rdiv(gen(2*j)*cst_pi*cst_i,n),contextptr));
	      newv.push_back(exp(rdiv(gen(-2*j)*cst_pi*cst_i,n),contextptr));
	    }
	    else {
	      newv.push_back(exp(rdiv(gen(j)*cst_pi*cst_i,n/2),contextptr));
	      newv.push_back(exp(rdiv(gen(-j)*cst_pi*cst_i,n/2),contextptr));
	    }
	  }
	}
      }
      solve_ckrange(x,newv,isolate_mode,contextptr);
      v=mergevecteur(v,newv);
      return ;
    }
    gen b_over_2=rdiv(w[1],plus_two);
    if (b_over_2.type!=_FRAC){
      gen a=r2sym(w.front(),lv,contextptr);
      gen minus_b_over_2=r2sym(-b_over_2,lv,contextptr);
      gen delta_prime=pow(minus_b_over_2,2)-r2sym(w.front()*w.back(),lv,contextptr);
      if (!complexmode && is_strictly_positive(-delta_prime,contextptr))
	return;
      newv.push_back(rdiv(minus_b_over_2+sqrt(delta_prime,contextptr),a));
      if (!is_zero(delta_prime))
	newv.push_back(rdiv(minus_b_over_2-sqrt(delta_prime,contextptr),a));
    }
    else {
      gen two_a=r2sym(plus_two*w.front(),lv,contextptr);
      gen minus_b=r2sym(-w[1],lv,contextptr);
      gen delta=r2sym(w[1]*w[1]-gen(4)*w.front()*w.back(),lv,contextptr);
      if (!complexmode && is_positive(-delta,contextptr))
	return;
      newv.push_back(rdiv(minus_b+sqrt(delta,contextptr),two_a));
      newv.push_back(rdiv(minus_b-sqrt(delta,contextptr),two_a));
    }
    solve_ckrange(x,newv,isolate_mode,contextptr);
    v=mergevecteur(v,newv);
  }

  // v assumed to represent an irreducible dense 1-d poly
  vecteur solve(const vecteur & v,bool complexmode,GIAC_CONTEXT){
    vecteur res;
    int d=v.size()-1;
    if (d<1)
      return res;
    if (d==1){
      res.push_back(rdiv(-v.back(),v.front()));
      return res;
    }
    if (!is_one(v.front())){
      // if v is not monic, set Y=a*X
      gen a(v.front()),puissance(plus_one);
      vecteur w;
      w.reserve(d+1);
      for (int i=0;i<=d;++i,puissance=puissance*a)
	w.push_back(v[i]*puissance);
      return divvecteur(solve(divvecteur(w,a),complex_mode(contextptr),contextptr),a);
    }
    // should call sym2rxroot for extensions of extensions
    vecteur tmp(2,zero);
    tmp.front()=plus_one;
    if (d==2){
      gen b(v[1]),c(v[2]);
      gen bprime(rdiv(b,plus_two));
      if (!has_denominator(bprime)){
	gen delta(bprime*bprime-c);
	if (!complexmode && is_positive(-delta,contextptr))
	  return res;
	vecteur w(3,zero);
	w.front()=plus_one;
	w.back()=-delta;
	tmp.back()=-bprime;
	res.push_back(algebraic_EXTension(tmp,w));
	tmp.front()=minus_one;
	tmp.back()=-bprime;
	res.push_back(algebraic_EXTension(tmp,w));	
      }
      else {
	if (!complexmode && is_positive(4*c-b*b,contextptr))
	  return res;
	tmp.back()=zero;
	res.push_back(algebraic_EXTension(tmp,v));
	tmp.front()=minus_one;
	tmp.back()=-b;
	res.push_back(algebraic_EXTension(tmp,v));
      }
      return res;
    }
    // should return a list of d algebraic extension with order number
    res.push_back(algebraic_EXTension(tmp,v));
    return res;
  }

  // works for continuous functions only
  vecteur solve_inequation(const gen & e0,const identificateur & x,int direction,GIAC_CONTEXT){
    if (has_num_coeff(e0))
      setsizeerr("Unable to solve inequations with approx coeffs "+e0.print(contextptr));
    gen e(e0._SYMBptr->feuille._VECTptr->front()-e0._SYMBptr->feuille._VECTptr->back());
    if (is_zero(ratnormal(derive(e,x,contextptr))))
      setsizeerr("Inequation is constant with respect to "+x.print(contextptr));
    vecteur veq_not_singu,veq,singu;
    singu=find_singularities(e,x,0,contextptr);
    veq_not_singu=solve(e,x,0,contextptr);
    veq=mergevecteur(veq_not_singu,singu);
    vecteur range,excluded_not_singu(find_excluded(x,contextptr));
    vecteur excluded=mergevecteur(excluded_not_singu,singu);
    find_range(x,range,contextptr); 
    // From the lower bound of range to the higher bound
    // find the sign 
    if (range.size()!=1 || range.front().type!=_VECT)
      return vecteur(0);
    vecteur veq_excluded=mergevecteur(excluded,veq);
    vecteur rangev = *range.front()._VECTptr;
    if (rangev.size()==2){
      gen &a=rangev.front();
      gen & b=rangev.back();
      veq_excluded=protect_sort(veq_excluded,contextptr);
      // keep only values inside a,b
      range=vecteur(1,a);
      const_iterateur it=veq_excluded.begin(),itend=veq_excluded.end();
      for (;it!=itend;++it){
	if (is_strictly_greater(*it,a,contextptr))
	  break;
      }
      for (;it!=itend;++it){
	if (is_greater(*it,b,contextptr))
	  break;
	range.push_back(*it);
      }
      range.push_back(b);
    }
    else {
      range=mergevecteur(rangev,veq_excluded);
      range=protect_sort(range,contextptr);
    }
    vecteur res;
    int s=range.size();
    if (s<2)
      setsizeerr();
    if (s==2 && range[0]==minus_inf && range[1]==plus_inf){
      gen test=sign(subst(e,x,0,false,contextptr),contextptr);
      if (direction<0)
	test=-test;
      if (is_one(test))
	return vecteur(1,x);
      if (is_one(-test))
	return vecteur(0);
      setsizeerr("Unable to check sign "+test.print());
    }
    for (int i=0;i<s-1;++i){
      gen l=range[i],m=range[i+1];
      if (l==m)
	continue;
      gen test;
      if (l==minus_inf)
	test=m-1;
      else {
	if (m==plus_inf)
	  test=l+1;
	else
	  test=(l+m)/2;
      }
      test=eval(subst(e0,x,test,false,contextptr),eval_level(contextptr),contextptr);
      if (test!=1)
	continue;
      gen symb_sup,symb_inf;
      if (equalposcomp(excluded_not_singu,l) || equalposcomp(singu,l) ||
	  ( !(direction %2) && equalposcomp(veq_not_singu,l)) )
	symb_inf=symb_superieur_strict(x,l);
      else
	symb_inf=symb_superieur_egal(x,l);
      if (equalposcomp(excluded_not_singu,m) || equalposcomp(singu,m) ||
	  ( !(direction %2) && equalposcomp(veq_not_singu,m)) )
	symb_sup=symb_inferieur_strict(x,m);
      else
	symb_sup=symb_inferieur_egal(x,m);
      if (l==minus_inf)
	res.push_back(symb_sup);
      else {
	if (m==plus_inf)
	  res.push_back(symb_inf);
	else
	  res.push_back(symbolic(at_and,makevecteur(symb_inf,symb_sup))); 
      }
    }
    return res;
  }

  void modsolve(const gen & e,const identificateur & x,const gen & modulo,vecteur &v,GIAC_CONTEXT){
    if (modulo.type!=_INT_)
      setdimerr("Modular equation with modulo too large");
    int m=modulo.val;
    for (int i=0;i<m;++i){
      gen tmp=subst(e,x,i,false,contextptr);
      if (is_zero(tmp.eval(eval_level(contextptr),contextptr)))
	v.push_back(i);
    }
  }

  void clean(gen & e,const identificateur & x,GIAC_CONTEXT){
    if (e.type!=_SYMB)
      return;
    if (e._SYMBptr->sommet==at_inv || (e._SYMBptr->sommet==at_pow && ck_is_positive(-e._SYMBptr->feuille._VECTptr->back(),contextptr))){
      gen ef=e._SYMBptr->feuille;
      if (e._SYMBptr->sommet==at_pow)
	ef=ef._VECTptr->front();
      // search for a tan in the variables
      vecteur lv(lvarx(e,x));
      if (lv.size()!=1)
	e=1;
      else {
	gen xvar(lv.front());
	if (!xvar.is_symb_of_sommet(at_tan))
	  e=1;
      }
      return;
    }
    if (e._SYMBptr->sommet==at_prod){
      gen ef=e._SYMBptr->feuille;
      if (ef.type!=_VECT)
	return;
      vecteur v=*ef._VECTptr;
      int vs=v.size();
      for (int i=0;i<vs;++i)
	clean(v[i],x,contextptr);
      ef=gen(v,ef.subtype);
      e=symbolic(at_prod,ef);
    }
  }
  
  // detect product and powers
  void solve(const gen & e,const identificateur & x,vecteur &v,int isolate_mode,GIAC_CONTEXT){
    if (is_zero(e)){
      v.push_back(x);
      return;
    }
    switch (e.type){
    case _IDNT:
      if (*e._IDNTptr==x && !equalposcomp(find_excluded(x,contextptr),zero))
	addtolvar(zero,v);
      return;
    case _SYMB:
      if ( e._SYMBptr->sommet==at_pow && ck_is_strictly_positive(e._SYMBptr->feuille._VECTptr->back(),contextptr) ){
	vecteur tmpv;
	solve(e._SYMBptr->feuille._VECTptr->front(),x,tmpv,isolate_mode,contextptr);
	int ncopy=1;
	// make copies of the answer (xcas_mode(contextptr)==1 compatibility)
	if (xcas_mode(contextptr)==1 && e._SYMBptr->feuille._VECTptr->back().type==_INT_)
	  ncopy=e._SYMBptr->feuille._VECTptr->back().val;
	const_iterateur it=tmpv.begin(),itend=tmpv.end();
	for (;it!=itend;++it){
	  for (int i=0;i<ncopy;++i)
	    v.push_back(*it);
	}
	return;
      }
      if (e._SYMBptr->sommet==at_prod){
	const_iterateur it=e._SYMBptr->feuille._VECTptr->begin(),itend=e._SYMBptr->feuille._VECTptr->end();
	for (;it!=itend;++it)
	  solve(*it,x,v,isolate_mode,contextptr);
	return;
      }
      if (e._SYMBptr->sommet==at_neg){
	solve(e._SYMBptr->feuille,x,v,isolate_mode,contextptr);
	return;
      }
      if (!(isolate_mode & 2) && 
	  (e._SYMBptr->sommet==at_inv || (e._SYMBptr->sommet==at_pow && ck_is_positive(-e._SYMBptr->feuille._VECTptr->back(),contextptr)))
	  ){
	gen ef=e._SYMBptr->feuille;
	if (e._SYMBptr->sommet==at_pow)
	  ef=ef._VECTptr->front();
	// search for a tan in the variables
	vecteur lv(lvarx(e,x));
	if (lv.size()!=1)
	  return;
	gen xvar(lv.front());
	if (!xvar.is_symb_of_sommet(at_tan))
	  return;
	gen arg=xvar._SYMBptr->feuille;
	// solve arg=pi/2[pi]
	in_solve(arg-isolate_tan(plus_inf,isolate_mode,contextptr),x,v,isolate_mode,contextptr);
	return;
      }
      in_solve(e,x,v,isolate_mode,contextptr);
      break;
    default:
      return;
    }
  }

  // find the arguments of sqrt inside expression e
  vecteur lvarfracpow(const gen & e){
    vecteur l0=lop(e,at_pow),l;
    const_iterateur it=l0.begin(),itend=l0.end();
    for (;it!=itend;++it){
      vecteur & arg=*it->_SYMBptr->feuille._VECTptr;
      gen g=arg[1],expnum,expden;
      if (g.type==_FRAC){
	expnum=g._FRACptr->num;
	expden=g._FRACptr->den;
      }
      else {
	if ( (g.type!=_SYMB) || (g._SYMBptr->sommet!=at_prod) )
	  continue;
	gen & arg1=g._SYMBptr->feuille;
	if (arg1.type!=_VECT)
	  continue;
	vecteur & v=*arg1._VECTptr;
	if ( (v.size()!=2) || (v[1].type!=_SYMB) || (v[1]._SYMBptr->sommet==at_inv) )
	  continue;
	expnum=v[0];
	expden=v[1]._SYMBptr->feuille;
      }
      if (expden.type!=_INT_)
	continue;
      l.push_back(arg[0]);
      l.push_back(expden.val);
      l.push_back(*it);
    }
    return l;
  }

  vecteur lvarfracpow(const gen & g,const identificateur & x,GIAC_CONTEXT){
    vecteur l0=lvarfracpow(g),l;
    const_iterateur it=l0.begin(),itend=l0.end();
    for (;it!=itend;++it){
      if (!is_zero(derive(*it,x,contextptr))){
	l.push_back(*it);
	++it;
	l.push_back(*it);
	++it;
	l.push_back(*it);
      }
      else
	it+=2;
    }
    return l;
  }

  void solve_fracpow(const gen & e,const identificateur & x,vecteur equations,const vecteur & listvars,vecteur & fullres,int isolate_mode,GIAC_CONTEXT){
    if (e.type==_IDNT){
      if (*e._IDNTptr==x && !equalposcomp(find_excluded(x,contextptr),zero)){
	addtolvar(zero,fullres);
	return;
      }
    }
    if (e.type==_SYMB){
      if ( (e._SYMBptr->sommet==at_pow) && (ck_is_positive(e._SYMBptr->feuille._VECTptr->back(),contextptr)) ){
	solve_fracpow(e._SYMBptr->feuille._VECTptr->front(),x,equations,listvars,fullres,isolate_mode,contextptr);
	return;
      }
      if (e._SYMBptr->sommet==at_prod){
	const_iterateur it=e._SYMBptr->feuille._VECTptr->begin(),itend=e._SYMBptr->feuille._VECTptr->end();
	for (;it!=itend;++it)
	  solve_fracpow(*it,x,equations,listvars,fullres,isolate_mode,contextptr);
	return;
      }
      if (e._SYMBptr->sommet==at_neg){
	solve_fracpow(e._SYMBptr->feuille,x,equations,listvars,fullres,isolate_mode,contextptr);
	return;
      }
    } // end if (e.type==_SYMB)
    /*
      // code with resultant in all var except the first one
      // disadvantage: does not check that listvars[i] are admissible
      // example assume(M<0); solve(sqrt(x)=M);
    gen expr(e);
    int s=listvars.size();
    for (int i=1;i<s;++i){
      // expr must be rationnal wrt listvars[i]
      vecteur vtmp(1,listvars[i]);
      if (listvars[i].type!=_IDNT)
	setsizeerr();
      rlvarx(expr,*listvars[i]._IDNTptr,vtmp);
      // IMPROVE: maybe a function applied to expr is rationnal
      if (vtmp.size()!=1)
	setsizeerr("Solve with fractional power:"+expr.print(contextptr)+" is not rationnal w.r.t. "+listvars[i].print(contextptr));
      if (!is_zero(derive(expr,listvars[i],contextptr)))
	expr=_resultant(makevecteur(expr,equations[i-1],listvars[i]),contextptr);
    }
    expr=ratfactor(expr,false,contextptr);
    if (is_zero(derive(expr,x,contextptr)))
      return;
    solve(expr,x,fullres,isolate_mode,contextptr);
    return;
    */
    // old code with Groebner basis
    equations.push_back(e);      
    vecteur res=gsolve(equations,listvars,complex_mode(contextptr),contextptr);
    iterateur it=res.begin(),itend=res.end();
    for (;it!=itend;++it)
      *it=(*it)[0];
    _purge(vecteur(listvars.begin()+1,listvars.end()),contextptr);
    if (listvars[0].type==_IDNT){
      fullres=mergevecteur(res,fullres);
      return;
    }
    // recursive call to solve composevar=*it with respect to x
    for (it=res.begin();it!=itend;++it){
      fullres=mergevecteur(fullres,solve(*it-listvars[0],x,isolate_mode,contextptr));
    }
  }

  vecteur solve(const gen & e,const identificateur & x,int isolate_mode,GIAC_CONTEXT){
    bool complexmode=isolate_mode & 1;
    gen expr(e);
    gen modulo;
    if (has_mod_coeff(expr,modulo)){
      vecteur v;
      modsolve(expr,x,modulo,v,contextptr);
      return v;
    }
    // Inequation?
    if (e.type==_SYMB){ 
      if (e._SYMBptr->sommet==at_inferieur_strict)
	return solve_inequation(e,x,-2,contextptr);
      if (e._SYMBptr->sommet==at_inferieur_egal)
	return solve_inequation(e,x,-1,contextptr);
      if (e._SYMBptr->sommet==at_superieur_strict)
	return solve_inequation(e,x,2,contextptr);
      if (e._SYMBptr->sommet==at_superieur_egal)
	return solve_inequation(e,x,1,contextptr);
      if (e._SYMBptr->sommet==at_equal ||e._SYMBptr->sommet==at_same)
	expr = e._SYMBptr->feuille._VECTptr->front()-e._SYMBptr->feuille._VECTptr->back();
    }
    clean(expr,x,contextptr);
    // Check for re/im/conj in complexmode
    if (complexmode){
      vecteur lc=mergevecteur(lop(expr,at_conj),mergevecteur(lop(expr,at_re),lop(expr,at_im)));
      int s=lc.size();
      for (int i=0;i<s;++i){
	gen f=lc[i]._SYMBptr->feuille;
	if (!is_zero(derive(f,x,contextptr))){
	  identificateur xrei(" x"),ximi(" y");
	  gen xre(xrei),xim(ximi);
	  bool savec=complex_mode(contextptr);
	  bool savecv=complex_variables(contextptr);
	  complex_mode(false,contextptr);
	  complex_variables(false,contextptr);
	  gen tmp=subst(e,x,xre+cst_i*xim,false,contextptr);
	  vecteur res=gsolve(makevecteur(re(tmp,contextptr),im(tmp,contextptr)),makevecteur(xre,xim),false,contextptr);
	  complex_mode(savec,contextptr);
	  complex_variables(savecv,contextptr);
	  s=res.size();
	  for (int j=0;j<s;++j){
	    if (res[j].type==_VECT && res[j]._VECTptr->size()==2){
	      gen a=res[j]._VECTptr->front();
	      gen b=res[j]._VECTptr->back();
	      if (is_zero(a))
		res[j]=cst_i*b;
	      else {
		if (is_zero(b))
		  res[j]=a;
		else
		  res[j]=symbolic(at_plus,gen(makevecteur(a,cst_i*b),_SEQ__VECT));
	      }
	    }
	  }
	  return res;
	}
      }
    }
    if ( (approx_mode(contextptr) || has_num_coeff(e)) && lidnt(e)==makevecteur(x))
      return gen2vecteur(_fsolve(gen(makevecteur(e,x),_SEQ__VECT),contextptr));
    // should rewrite e in terms of a minimal number of variables
    // first factorization of e
    // Checking for abs
    vecteur la;
    if (!complexmode)
      la=lop(expr,at_abs);
    const_iterateur itla=la.begin(),itlaend=la.end();
    for (;itla!=itlaend;++itla){
      gen g=itla->_SYMBptr->feuille;
      if (is_zero(derive(g,x,contextptr)))
	continue;
      vecteur res;
      gen ee=subst(expr,*itla,g,false,contextptr);
      vecteur v1=solve(ee,x,isolate_mode,contextptr);
      const_iterateur it=v1.begin(),itend=v1.end();
      for (;it!=itend;++it){
	gen g1=subst(g,x,*it,false,contextptr);
	if (ratnormal(abs(g1,contextptr)-g1)==0) 
	  res.push_back(*it);
      }
      ee=subst(expr,*itla,-g,false,contextptr);
      v1=solve(ee,x,isolate_mode,contextptr);
      it=v1.begin(); itend=v1.end();
      for (;it!=itend;++it){
	gen g1=subst(g,x,*it,false,contextptr);
	if (ratnormal(abs(g1,contextptr)+g1)==0) 
	  res.push_back(*it);
      }
      return res;
    }
    vecteur lv(lvarx(expr,x));
    int s=lv.size();
    if (s>1)
      expr=halftan_hyp2exp(expr,contextptr);
    // Checking for fractional power
    // Remark: algebraic extension could also be solved using resultant
    vecteur ls(lvarfracpow(expr,x,contextptr));
    if (!ls.empty()){ // Use auxiliary variables
      int s=ls.size()/3;
      vecteur substin,substout,equations,listvars(lvarx(expr,x,true));
      // remove ls from listvars, add aux var instead
      for (int i=0;i<s;++i){
	gen lsvar=ls[3*i+2];
	int j=equalposcomp(listvars,ls);
	if (j)
	  listvars.erase(listvars.begin()+j-1);
      }
      if (listvars.size()!=1)
	setsizeerr("unable to isolate "+gen(listvars).print(contextptr));
      for (int i=0;i<s;++i){
	gen lsvar=ls[3*i+2];
	substin.push_back(lsvar);
	gen tmp("c__"+print_INT_(i),contextptr);
	if (!(ls[3*i+1].val %2))
	  assumesymbolic(symb_superieur_egal(tmp,0),0,contextptr); 
	listvars.push_back(tmp);
	substout.push_back(tmp);
	equations.push_back(pow(tmp,ls[3*i+1],contextptr)-ls[3*i]);
      }
      gen expr1=subst(expr,substin,substout,false,contextptr);
      expr1=ratfactor(expr1,false,contextptr);
      vecteur fullres;
      solve_fracpow(expr1,x,equations,listvars,fullres,isolate_mode,contextptr);
      // Check that expr at x=fullres is 0
      // Only if expr1 does not depend on other variables than x
      vecteur othervar(lidnt(expr1)),res;
      if (othervar.size()==listvars.size()){
	const_iterateur it=fullres.begin(),itend=fullres.end();
	for (;it!=itend;++it){
	  vecteur algv=alg_lvar(*it);
	  if ( (!algv.empty() && algv.front().type==_VECT && !algv.front()._VECTptr->empty()) || is_zero(limit(expr,x,*it,0,contextptr)))
	    res.push_back(*it);
	}
      }
      else
	res=fullres;
      return res;
    }
    lv=lvarx(expr,x);
    if (lv.size()>1){
      gen tmp=ratfactor(simplify(expr,contextptr),false,contextptr);
      int lvs=0;
      if (tmp.is_symb_of_sommet(at_prod) && tmp._SYMBptr->feuille.type==_VECT){
	vecteur & f=*tmp._SYMBptr->feuille._VECTptr;
	int fs=f.size();
	for (int i=0;i<fs;++i){
	  lvs=lvarx(f[i],x).size();
	  if (lvs>1)
	    break;
	}
      }
      else
	lvs=lvarx(tmp,x).size();
      if (lvs<2)
	expr=tmp;
    }
    // -> exp/ln
    expr=pow2expln(expr,x,contextptr);
    bool setcplx=complexmode && complex_mode(contextptr)==false;
    if (setcplx)
      complex_mode(true,contextptr);
    expr=ratfactor(expr,false,contextptr); // factor in complex or real mode
    if (setcplx)
      complex_mode(false,contextptr);
    lv=lvarx(expr,x);
    s=lv.size();
    vecteur v;
    if (!s){
      if (is_zero(expr))
	v.push_back(x);
      return v;
    }
    solve(expr,x,v,isolate_mode,contextptr);
    if (!(isolate_mode & 2)){
      // check solutions if there is a tan inside
      for (int i=0;i<s;++i){
	if (lv[i].is_symb_of_sommet(at_tan)){
	  vecteur res;
	  const_iterateur it=v.begin(),itend=v.end();
	  for (;it!=itend;++it){
	    if (is_zero(recursive_normal(limit(_tan2sincos2(expr,contextptr),x,*it,0,contextptr),contextptr)) || is_zero(recursive_normal(limit(expr,x,*it,0,contextptr),contextptr)))
	      res.push_back(*it);
	  }
	  return res;
	}
      }
    }
    return v;
  }

  gen remove_and(const gen & g,const unary_function_ptr & u){
    if (g.type==_VECT){
      vecteur res;
      const_iterateur it=g._VECTptr->begin(),itend=g._VECTptr->end();
      for (;it!=itend;++it){
	gen tmp=remove_and(*it,u);
	if (tmp.type!=_VECT){
	  tmp=remove_and(*it,at_and);
	  res.push_back(tmp);
	}
	else
	  res=mergevecteur(res,*tmp._VECTptr);
      }
      return res;
    }
    if (!g.is_symb_of_sommet(u))
      return g;
    return remove_and(g._SYMBptr->feuille,u);
  }

  vecteur solve(const gen & e,const gen & x,int isolate_mode,GIAC_CONTEXT){
    bool complexmode=isolate_mode & 1;
    vecteur res;
    if (x.type!=_IDNT){
      if (x.type==_VECT && x._VECTptr->size()==1 && e.type==_VECT && e._VECTptr->size()==1){
	vecteur res=solve(e._VECTptr->front(),x._VECTptr->front(),isolate_mode,contextptr);
	iterateur it=res.begin(),itend=res.end();
	for (;it!=itend;++it)
	  *it=vecteur(1,*it);	
	return res;
      }
      if ( (x.type==_VECT) && (e.type==_VECT) )
	return gsolve(*e._VECTptr,*x._VECTptr,complexmode,contextptr);
      identificateur xx("x");
      res=solve(subst(e,x,xx,false,contextptr),xx,isolate_mode,contextptr);
      return res;
    }
    if (e.type==_VECT){
      const_iterateur it=e._VECTptr->begin(),itend=e._VECTptr->end();
      res.reserve(itend-it);
      for (;it!=itend;++it)
	res.push_back(solve(*it,*x._IDNTptr,isolate_mode,contextptr));
    }
    else
      res=solve(e,*x._IDNTptr,isolate_mode,contextptr);
    return res;
  }

  gen symb_solution(const gen & g,const gen & var,GIAC_CONTEXT){
    if (var.type!=_VECT){
      if (var.type==_IDNT && !lvarx(g,*var._IDNTptr).empty())
	return g;
      else
	return symbolic(at_equal,makevecteur(var,g));
    }
    vecteur v=*var._VECTptr;
    if (g.type!=_VECT || g._VECTptr->size()!=v.size())
      setsizeerr();
    iterateur it=v.begin(),itend=v.end(),jt=g._VECTptr->begin();
    for (;it!=itend;++it,++jt)
      *it=symbolic(at_equal,makevecteur(*it,*jt));
    if (xcas_mode(contextptr)==3)
      return symbolic(at_and,v);
    else
      return v;
  }

  gen quote_inferieur_strict(const gen & g,GIAC_CONTEXT){
    return symbolic(at_quote,symbolic(at_inferieur_strict,eval(g,eval_level(contextptr),contextptr)));
  }

  gen quote_superieur_strict(const gen & g,GIAC_CONTEXT){
    return symbolic(at_quote,symbolic(at_superieur_strict,eval(g,eval_level(contextptr),contextptr)));
  }

  gen quote_inferieur_egal(const gen & g,GIAC_CONTEXT){
    return symbolic(at_quote,symbolic(at_inferieur_egal,eval(g,eval_level(contextptr),contextptr)));
  }

  gen quote_superieur_egal(const gen & g,GIAC_CONTEXT){
    return symbolic(at_quote,symbolic(at_superieur_egal,eval(g,eval_level(contextptr),contextptr)));
  }

  gen quote_conj(const gen & g,GIAC_CONTEXT){
    return symbolic(at_quote,symbolic(at_conj,eval(g,eval_level(contextptr),contextptr)));
  }

  gen quote_re(const gen & g,GIAC_CONTEXT){
    return symbolic(at_quote,symbolic(at_re,eval(g,eval_level(contextptr),contextptr)));
  }

  gen quote_im(const gen & g,GIAC_CONTEXT){
    return symbolic(at_quote,symbolic(at_im,eval(g,eval_level(contextptr),contextptr)));
  }

  vecteur solvepreprocess(const gen & args,bool complexmode,GIAC_CONTEXT){
    gen g(args);
    if (g.type==_VECT && !g._VECTptr->empty() && g._VECTptr->front().is_symb_of_sommet(at_and)){
      vecteur v(*g._VECTptr);
      v.front()=remove_and(v.front(),at_and);
      g=gen(v,g.subtype);
    }
    // quote < <= > and >=
    vector<unary_function_ptr> quote_inf;
    quote_inf.push_back(at_inferieur_strict);
    quote_inf.push_back(at_inferieur_egal);
    quote_inf.push_back(at_superieur_strict);
    quote_inf.push_back(at_superieur_egal);
    if (complexmode){
      quote_inf.push_back(at_conj);
      quote_inf.push_back(at_re);
      quote_inf.push_back(at_im);
    }
    vector< gen_op_context > quote_inf_v;
    quote_inf_v.push_back(quote_inferieur_strict);
    quote_inf_v.push_back(quote_inferieur_egal);
    quote_inf_v.push_back(quote_superieur_strict);
    quote_inf_v.push_back(quote_superieur_egal);
    if (complexmode){
      quote_inf_v.push_back(quote_conj);
      quote_inf_v.push_back(quote_re);
      quote_inf_v.push_back(quote_im);
    }
    g=subst(g,quote_inf,quote_inf_v,true,contextptr);
    return plotpreprocess(g,contextptr);
  }

  gen solvepostprocess(const gen & g,const gen & x,GIAC_CONTEXT){
    if (g.type!=_VECT)
      return g;
    vecteur res=*g._VECTptr;
    // convert solution to an expression
    iterateur it=res.begin(),itend=res.end();
    if (it==itend)
      return res;
    if (x.type==_VECT || xcas_mode(contextptr)==3){
      for (;it!=itend;++it)
	*it=symb_solution(*it,x,contextptr);
    }
    if (xcas_mode(contextptr)==3)
      return symbolic(at_ou,res);
    if (xcas_mode(contextptr)==2)
      return gen(res,_SET__VECT);
    return gen(res,_SEQ__VECT);
  }

  gen _solve(const gen & args,GIAC_CONTEXT){
    int isolate_mode=complex_mode(contextptr) | (all_trig_sol(contextptr) << 1);
    vecteur v(solvepreprocess(args,complex_mode(contextptr),contextptr));
    int s=v.size();
    if (s==2 && v[1].is_symb_of_sommet(at_equal))
      return _fsolve(gen(makevecteur(v[0],v[1]._SYMBptr->feuille[0],v[1]._SYMBptr->feuille[1]),_SEQ__VECT),contextptr);
    if (s>2)
      return _fsolve(args,contextptr);
    gen arg1(v.front());
    if (arg1.type==_VECT){ // Flatten equations which are list of equations
      vecteur w;
      const_iterateur it=arg1._VECTptr->begin(),itend=arg1._VECTptr->end();
      for (;it!=itend;++it){
	gen tmp=equal2diff(*it);
	if (tmp.type==_VECT){
	  const_iterateur jt=tmp._VECTptr->begin(),jtend=tmp._VECTptr->end();
	  for (;jt!=jtend;++jt)
	    w.push_back(*jt);
	}
	else
	  w.push_back(tmp);
      }
      arg1=w;
    }
    vecteur res=solve(apply(arg1,equal2diff),v.back(),isolate_mode,contextptr);
    // if (is_fully_numeric(res))
    if (lidnt(res).empty())
      res=protect_sort(res,contextptr);
    if (!xcas_mode(contextptr))
      return res;
    return solvepostprocess(res,v[1],contextptr);
  }
  const string _solve_s("solve");
  unary_function_eval __solve(&_solve,_solve_s);
  unary_function_ptr at_solve (&__solve,_QUOTE_ARGUMENTS,true);

  double nan(){
    double x=0.0;
    return 0.0/x;
  }
#ifdef HAVE_LIBGSL
  // p should point a vector with elements the expression f(x) and x
  // OR with f(x), f'(x) and x
  double my_f (double x0, void * p) {
    gen & params = * ((gen *)p) ;
#ifdef DEBUG_SUPPORT
    if ( (params.type!=_VECT) || (params._VECTptr->size()<2))
      setsizeerr("solve.cc/my_f");
#endif	// DEBUG_SUPPORT
    gen & f=params._VECTptr->front();
    gen & x=params._VECTptr->back();
    gen res=evalf_double(subst(f,x,x0,false,context0),1,context0);
    if (res.type>_IDNT)
      setsizeerr();
    if (res.type!=_DOUBLE_)
      return nan();
    else
      return res._DOUBLE_val;
  }

  double my_df (double x0, void * p) {
    gen & params = * ((gen *)p) ;
#ifdef DEBUG_SUPPORT
    if ( (params.type!=_VECT) || (params._VECTptr->size()!=3))
      setsizeerr("solve.cc/my_df");
#endif	// DEBUG_SUPPORT
    gen & f=(*params._VECTptr)[1];
    gen & x=params._VECTptr->back();
    gen res=evalf_double(subst(f,x,x0,false,context0),1,context0);
    if (res.type>_IDNT)
      setsizeerr();
    if (res.type!=_DOUBLE_)
      return nan();
    else
      return res._DOUBLE_val;
  }

  void my_fdf (double x0, void * p,double * fx,double * dfx) {
    gen & params = * ((gen *)p) ;
#ifdef DEBUG_SUPPORT
    if ( (params.type!=_VECT) || (params._VECTptr->size()!=3))
      setsizeerr("solve.cc/my_fdf");
#endif	// DEBUG_SUPPORT
    gen & f=params._VECTptr->front();
    gen & df=(*params._VECTptr)[1];
    gen & x=params._VECTptr->back();
    gen res=evalf_double(subst(f,x,x0,false,context0),1,context0);
    if (res.type!=_DOUBLE_)
      *fx=nan();
    else
      *fx=res._DOUBLE_val;
    res=evalf_double(subst(df,x,x0,false,context0),1,context0);
    if (res.type!=_DOUBLE_)
      *dfx=nan();
    else
      *dfx=res._DOUBLE_val;
  }

  int my_F (const gsl_vector * x0, void * p,gsl_vector * F) {
    gen & params = * ((gen *)p) ;
#ifdef DEBUG_SUPPORT
    if ( (params.type!=_VECT) || (params._VECTptr->size()<2))
      setsizeerr("solve.cc/my_F");
#endif	// DEBUG_SUPPORT
    gen & f=params._VECTptr->front();
    gen & x=params._VECTptr->back();
    gen res=evalf_double(subst(f,x,gsl_vector2vecteur(x0),false,context0),1,context0);
    if (res.type!=_VECT)
      return !GSL_SUCCESS;
    return vecteur2gsl_vector(*res._VECTptr,F,context0);
  }

  int my_dF (const gsl_vector *x0, void * p,gsl_matrix * J) {
    gen & params = * ((gen *)p) ;
#ifdef DEBUG_SUPPORT
    if ( (params.type!=_VECT) || (params._VECTptr->size()!=3))
      setsizeerr("solve.cc/my_dF");
#endif	// DEBUG_SUPPORT
    gen & f=(*params._VECTptr)[1];
    gen & x=params._VECTptr->back();
    gen res=evalf_double(subst(f,x,gsl_vector2vecteur(x0),false,context0),1,context0);
    if (res.type!=_VECT)
      return !GSL_SUCCESS;
    else
      return matrice2gsl_matrix(*res._VECTptr,J,context0);
  }

  int my_FdF (const gsl_vector * x0, void * p,gsl_vector * fx,gsl_matrix * dfx) {
    gen & params = * ((gen *)p) ;
#ifdef DEBUG_SUPPORT
    if ( (params.type!=_VECT) || (params._VECTptr->size()!=3))
      setsizeerr("solve.cc/my_FdF");
#endif	// DEBUG_SUPPORT
    gen & f=params._VECTptr->front();
    gen & df=(*params._VECTptr)[1];
    gen & x=params._VECTptr->back();
    gen g0=gsl_vector2vecteur(x0);
    gen res=evalf_double(subst(f,x,g0,false,context0),1,context0);
    if (res.type!=_VECT)
      return !GSL_SUCCESS;
    int ires=vecteur2gsl_vector(*res._VECTptr,fx,context0);
    if (ires!=GSL_SUCCESS)
      return !GSL_SUCCESS;
    res=evalf_double(subst(df,x,g0,false,context0),1,context0);
    if (res.type!=_VECT)
      return !GSL_SUCCESS;
    return matrice2gsl_matrix(*res._VECTptr,dfx,context0);
  }

  gen msolve(const gen & f,const vecteur & vars,const vecteur & g,int method,double eps,GIAC_CONTEXT){
    vecteur guess(g);
    bool with_derivative=false;
    int dim=vars.size();
    switch (method){
    case _NEWTONJ_SOLVER: case _HYBRIDSJ_SOLVER: case _HYBRIDJ_SOLVER:
      with_derivative=true;
      break;
    case _DNEWTON_SOLVER: case _HYBRIDS_SOLVER: case _HYBRID_SOLVER:
      with_derivative=false;
      break;
    }
    if (with_derivative){
      gen params(makevecteur(f,mtran(*derive(f,vars,contextptr)._VECTptr),vars));
      gsl_multiroot_function_fdf FDF;
      FDF.f=&my_F;
      FDF.df=&my_dF;
      FDF.fdf=&my_FdF;
      FDF.n=dim;
      FDF.params=&params;
      const gsl_multiroot_fdfsolver_type * T=0;
      switch (method){
      case _NEWTONJ_SOLVER: 
	T=gsl_multiroot_fdfsolver_gnewton;
	break;
      case _HYBRIDSJ_SOLVER:
	T=gsl_multiroot_fdfsolver_hybridsj;
	break;
      case _HYBRIDJ_SOLVER:
	T=gsl_multiroot_fdfsolver_hybridj;
	break;
      }
      gsl_multiroot_fdfsolver * s= gsl_multiroot_fdfsolver_alloc (T, dim);
      gsl_vector * X=vecteur2gsl_vector(guess,contextptr);
      gsl_multiroot_fdfsolver_set (s, &FDF, X);
      int maxiter=SOLVER_MAX_ITERATE,res=0;
      vecteur oldguess;
      for (;maxiter;--maxiter){
	oldguess=guess;
	res=gsl_multiroot_fdfsolver_iterate(s);
	if ( (res==GSL_EBADFUNC) || (res==GSL_ENOPROG) )
	  break;
	guess=gsl_vector2vecteur(gsl_multiroot_fdfsolver_root(s));
	if (is_strictly_greater(eps,abs(guess-oldguess,contextptr),contextptr))
	  break;
      }
      gsl_multiroot_fdfsolver_free(s);
      if ( (res==GSL_EBADFUNC) || (res==GSL_ENOPROG) )
	setsizeerr("Not found");
      return guess;
    }
    else {
      gen params(makevecteur(f,vars));
      gsl_multiroot_function F;
      F.f=&my_F;
      F.n=dim;
      F.params=&params;
      const gsl_multiroot_fsolver_type * T=0;
      switch (method){
      case _DNEWTON_SOLVER: 
	T=gsl_multiroot_fsolver_dnewton;
	break;
      case _HYBRIDS_SOLVER:
	T=gsl_multiroot_fsolver_hybrids;
	break;
      case _HYBRID_SOLVER:
	T=gsl_multiroot_fsolver_hybrid;
	break;
      }
      gsl_multiroot_fsolver * s= gsl_multiroot_fsolver_alloc (T, dim);
      gsl_vector * X=vecteur2gsl_vector(guess,contextptr);
      gsl_multiroot_fsolver_set (s, &F, X);
      int maxiter=SOLVER_MAX_ITERATE,res=0;
      vecteur oldguess;
      for (;maxiter;--maxiter){
	oldguess=guess;
	res=gsl_multiroot_fsolver_iterate(s);
	if ( (res==GSL_EBADFUNC) || (res==GSL_ENOPROG) )
	  break;
	guess=gsl_vector2vecteur(gsl_multiroot_fsolver_root(s));
	if (is_strictly_greater(eps,abs(guess-oldguess,contextptr),contextptr))
	  break;
      }
      gsl_multiroot_fsolver_free(s);
      if ( (res==GSL_EBADFUNC) || (res==GSL_ENOPROG) )
	setsizeerr("Not found");
      return guess;
    }
  }
#endif // HAVE_LIBGSL

  // fsolve(expr,var[,interval/guess,method])
  gen _fsolve(const gen & args,GIAC_CONTEXT){
    vecteur v(plotpreprocess(args,contextptr));
    double gsl_eps=epsilon(contextptr);
    int s=v.size();
    if (s<2)
      toofewargs("fsolve");
    gen v0=remove_equal(v[0]);
    if (s==2 && v[1].type==_IDNT){ 
      gen v00=v0;
      // no initial guess, check for poly-like equation
      vecteur lv(lvar(v00));
      int lvs=lv.size();
      if (lv==vecteur(1,v[1])){
	gen tmp=_e2r(makevecteur(v00,v[1]),contextptr);
	if (tmp.type==_FRAC)
	  tmp=tmp._FRACptr->num;
	tmp=evalf(tmp,eval_level(contextptr),contextptr);
	if (tmp.type==_VECT)
	  return complex_mode(contextptr)?proot(*tmp._VECTptr,epsilon(contextptr)):real_proot(*tmp._VECTptr,epsilon(contextptr),contextptr);
      }
      if (lvs>1)
	v00=halftan_hyp2exp(v00,contextptr);
      lv=lvar(v00);
      lvs=lv.size();
      if (lvs==1 && lv[0].type==_SYMB && lv[0]._SYMBptr->feuille.type!=_VECT){
	int pos=equalposcomp(solve_fcns_v,lv[0]._SYMBptr->sommet);
	if (pos){
	  gen tmp=_e2r(makevecteur(v00,lv[0]),contextptr);
	  if (tmp.type==_FRAC)
	    tmp=tmp._FRACptr->num;
	  tmp=evalf(tmp,eval_level(contextptr),contextptr);
	  if (tmp.type==_VECT){
	    vecteur res0=complex_mode(contextptr)?proot(*tmp._VECTptr,epsilon(contextptr)):real_proot(*tmp._VECTptr,epsilon(contextptr),contextptr);
	    vecteur res;
	    const_iterateur it=res0.begin(),itend=res0.end();
	    for (;it!=itend;++it){
	      vecteur res0val=gen2vecteur(isolate_fcns[pos-1](*it,complex_mode(contextptr),contextptr));
	      const_iterateur jt=res0val.begin(),jtend=res0val.end();
	      for (;jt!=jtend;++jt){
		gen fs=_fsolve(gen(makevecteur(lv[0]._SYMBptr->feuille-*jt,v[1]),_SEQ__VECT),contextptr);
		if (fs.type==_VECT)
		  res=mergevecteur(res,*fs._VECTptr);
		else
		  res.push_back(fs);
	      }
	    }
	    return res;
	  }
	}
      }
    }
    gen gguess;
    if (v[1].is_symb_of_sommet(at_equal)){
      gguess=v[1]._SYMBptr->feuille[1];
      v[1]=v[1]._SYMBptr->feuille[0];
      v.insert(v.begin()+2,gguess);
      ++s;
    }
    if (s>=3)
      gguess=v[2];
    // check method
    int method=_NEWTON_SOLVER;
    //int method=0;
    if ( (s>=5) && (v[4].type==_DOUBLE_) )
      gsl_eps=v[4]._DOUBLE_val;
    if ( (s>=4) && (v[3].type==_INT_) )
      method=v[3].val;
    if (v[1].type==_VECT){
      int dim=v[1]._VECTptr->size();
      if (!dim)
	setsizeerr();
      if (s>=3){
	if (gguess.type!=_VECT)
	  setsizeerr();
	if (gguess._VECTptr->size()!=unsigned(dim))
	  setsizeerr();
      }
      else {
	gguess=vecteur(dim);
	gguess[0]=(gnuplot_xmin+gnuplot_xmax)/2;
	if (dim>1)
	  gguess[1]=(gnuplot_ymin+gnuplot_ymax)/2;
	if (dim>2)
	  gguess[2]=(gnuplot_zmin+gnuplot_zmax)/2;
	if (dim>3)
	  gguess[3]=(gnuplot_tmin+gnuplot_tmax)/2;
      }
#ifdef HAVE_LIBGSL
      if (method!=_NEWTON_SOLVER)
	return msolve(v0,*v[1]._VECTptr,*gguess._VECTptr,method,gsl_eps,contextptr);
#endif
    }
#ifdef HAVE_LIBGSL
    if (method!=_NEWTON_SOLVER){
      bool with_derivative=false;
      switch (method){
      case _BISECTION_SOLVER: case _FALSEPOS_SOLVER: case _BRENT_SOLVER:
	with_derivative=false;
	break;
      case _NEWTON_SOLVER: case _SECANT_SOLVER: case _STEFFENSON_SOLVER:
	with_derivative=true;
	break;
      }
      gen params;
      if (with_derivative){
	params= makevecteur(v0,derive(v0,v[1],contextptr),v[1]);
	double guess((gnuplot_xmin+gnuplot_xmax)/2),oldguess;
	if (s>=3){
	  gen g=evalf(gguess,eval_level(contextptr),contextptr);
	  if (gguess.type==_DOUBLE_)
	    guess=gguess._DOUBLE_val;
	}
	gsl_function_fdf FDF ;     
	FDF.f = &my_f ;
	FDF.df = &my_df ;
	FDF.fdf = &my_fdf ;
	FDF.params = &params ;
	gsl_root_fdfsolver * slv=0;
	switch (method){
	case _NEWTON_SOLVER:
	  slv=gsl_root_fdfsolver_alloc (gsl_root_fdfsolver_newton);
	  break;
	case _SECANT_SOLVER: 
	  slv=gsl_root_fdfsolver_alloc (gsl_root_fdfsolver_secant);
	  break;
	case _STEFFENSON_SOLVER:
	  slv=gsl_root_fdfsolver_alloc (gsl_root_fdfsolver_steffenson);
	  break;
	}
	if (!slv)
	  setsizeerr();
	gsl_root_fdfsolver_set(slv,&FDF,guess);
	int maxiter=SOLVER_MAX_ITERATE,res=0;
	for (;maxiter;--maxiter){
	  oldguess=guess;
	  res=gsl_root_fdfsolver_iterate(slv);
	  guess=gsl_root_fdfsolver_root(slv);
	  if ( (res==GSL_EBADFUNC) || (res==GSL_EZERODIV) )
	    break;
	  if (fabs(guess-oldguess)<gsl_eps)
	    break;
	}
	gsl_root_fdfsolver_free(slv);
	if (!maxiter)
	setsizeerr("Not found");
	if ( (res==GSL_EBADFUNC) || (res==GSL_EZERODIV) )
	  return undef;
	else
	  return guess;
      }
      else {
	params= makevecteur(v0,v[1]);
	double x_low,x_high;
	if (s>=3) {
	  vecteur w;
	  if (gguess.type==_VECT)
	    w=*gguess._VECTptr;
	  else {
	    if ( (gguess.type==_SYMB) && (gguess._SYMBptr->sommet==at_interval) )
	      w=*gguess._SYMBptr->feuille._VECTptr;
	  }
	  if (w.size()!=2)
	    settypeerr();
	  gen low=w[0].evalf(eval_level(contextptr),contextptr);
	  gen high=w[1].evalf(eval_level(contextptr),contextptr);
	  if ( (low.type!=_DOUBLE_) || (high.type!=_DOUBLE_) )
	  setsizeerr();
	  x_low=low._DOUBLE_val;
	  x_high=high._DOUBLE_val;
	}
	else {
	  x_low=gnuplot_xmin;
	  x_high=gnuplot_xmax;
	}
	if (x_low>x_high){
	  double tmp=x_low;
	  x_low=x_high;
	  x_high=tmp;
	}
	gsl_function F ;
	F.function=&my_f;
	F.params = &params ;
	gsl_root_fsolver * slv =0 ;
	switch (method){
	case  _BISECTION_SOLVER:
	  slv=gsl_root_fsolver_alloc (gsl_root_fsolver_bisection);
	  break;
	case _FALSEPOS_SOLVER: 
	  slv=gsl_root_fsolver_alloc (gsl_root_fsolver_falsepos);
	  break;
	case _BRENT_SOLVER:
	  slv=gsl_root_fsolver_alloc (gsl_root_fsolver_brent);
	  break;
	}
	if (!slv)
	  setsizeerr();
	gsl_root_fsolver_set (slv,&F,x_low,x_high);
	int res=0;
	int maxiter=SOLVER_MAX_ITERATE;
	for (;maxiter && (x_high-x_low>gsl_eps);--maxiter){
	  res=gsl_root_fsolver_iterate(slv);
	  if (res==GSL_EBADFUNC)
	    break;
	  x_low=gsl_root_fsolver_x_lower(slv);
	  x_high= gsl_root_fsolver_x_upper(slv);
	}
	gsl_root_fsolver_free (slv);
	if (res==GSL_EBADFUNC)
	  return undef;
	return makevecteur(x_low,x_high);
      } // end if derivative
    }
#else // HAVE_LIBGSL
    if (method!=_NEWTON_SOLVER)
      setsizeerr("Not linked with GSL");
#endif    // HAVE_LIBGSL
    else  // newton method, call newton
      return newton(v0,v[1],gguess,NEWTON_DEFAULT_ITERATION,gsl_eps,1e-12,contextptr);
  } // end f_solve
  const string _fsolve_s("fsolve");
  unary_function_eval __fsolve(&_fsolve,_fsolve_s);
  unary_function_ptr at_fsolve (&__fsolve,_QUOTE_ARGUMENTS,true);

  vecteur sxa(const vecteur & sl_orig,const vecteur & x,GIAC_CONTEXT){
      vecteur sl(sl_orig);
      int d;
      d=x.size();
      int de;
      de=sl.size();
      for (int i=0;i<de;i++){
              //gen e:
              //e=sl[i];    
          if ( (sl[i].type==_SYMB) && ((*sl[i]._SYMBptr).sommet==at_equal || (*sl[i]._SYMBptr).sommet==at_same)){
              sl[i]=(*sl[i]._SYMBptr).feuille[0]-(*sl[i]._SYMBptr).feuille[1];
          }
      }
      vecteur A;
      for (int i=0;i<de;i++){
          vecteur li(d+1);
          gen lo=sl[i];
          for (int j=0;j<d;j++){
              lo=subst(lo,x[j],0,false,contextptr);
              li[j]=derive(sl[i],x[j],contextptr);
          }
          li[d]=lo;
          A.push_back(li);
      }
      return(A);
  }

  vecteur linsolve(const vecteur & sl,const vecteur & x,GIAC_CONTEXT){
    vecteur A; 
    if (ckmatrix(sl)){
      unsigned int n=sl.size();
      A=mtran(sl);
      if (ckmatrix(x)){
	if (x.size()==1){
	  if (x.front()._VECTptr->size()!=n)
	    setdimerr();
	  A.push_back(-x.front());
	}
	else {
	  if (x.size()!=n)
	    setdimerr();
	  matrice xm=mtran(x);
	  if (xm.size()!=1)
	    setsizeerr();
	  A.push_back(-xm.front());
	}
      }
      else {
	if (x.size()!=n)
	  setdimerr();
	A.push_back(-x);
      }
      A=mtran(A);
      vecteur B=-mker(A,contextptr);
      if (B.empty())
	return B;
      // The last element of B must have a non-zero last component
      vecteur Bend=*B.back()._VECTptr;
      gen last=Bend.back();
      if (is_zero(last))
	return vecteur(0);
      gen R=Bend/last;
      // The solution is sum(B[k]*Ck+Blast/last)
      int s=B.size();
      for (int k=0;k<s-1;k++)
	R=R+gen("C_"+print_INT_(k),contextptr)*B[k];
      vecteur res=*R._VECTptr;
      res.pop_back();
      return res;
    }
    A=sxa(sl,x,contextptr);
    vecteur B,R(x);
    gen rep;
    B=mrref(A,contextptr);
    //cout<<B<<endl;
    int d=x.size();
    int de=sl.size();
    for (int i=0; i<de;i++){
      vecteur li(d+1);
      for(int k=0;k<d+1;k++){
	li[k]=B[i][k];
      }
      int j;
      j=i;
      while (li[j]==0 && j<d){
	j=j+1;
      }
      if (j==d && !is_zero(li[d])){
	return vecteur(0);
      } 
      else {
	if (j<d){
	  rep=-li[d];
	  for (int k=j+1;k<d;k++){
	    rep=rep-li[k]*x[k];
	  }
	  rep=rdiv(rep,li[j]);
	  R[j]=rep;
	}
      }
    }
    return R;
  }

  gen equal2diff(const gen & g){
    if ( (g.type==_SYMB) && (g._SYMBptr->sommet==at_equal || g._SYMBptr->sommet==at_same) ){
      vecteur & v=*g._SYMBptr->feuille._VECTptr;
      return v[0]-v[1];
    }
    else
      return g;
  }

  gen symb_linsolve(const gen & syst,const gen & vars){
    return symbolic(at_linsolve,makevecteur(syst,vars));
  }
 
  gen linsolve(const gen & syst,const gen & vars,GIAC_CONTEXT){
    if ((syst.type!=_VECT)||(vars.type!=_VECT))
      return symb_linsolve(syst,vars);
    return normal(linsolve(*syst._VECTptr,*vars._VECTptr,contextptr),contextptr);
  }
  
  gen _linsolve(const gen & args,GIAC_CONTEXT){
    vecteur v(plotpreprocess(args,contextptr));
    int s=v.size();
    if (s!=2)
      toomanyargs("linsolve");
    if (v[1].type==_IDNT)
      v[1]=eval(v[1],eval_level(contextptr),contextptr);
    gen syst=apply(v[0],equal2diff),vars=v[1];
    return linsolve(syst,v[1],contextptr);
  }
  const string _linsolve_s("linsolve");
  unary_function_eval __linsolve(&_linsolve,_linsolve_s);
  unary_function_ptr at_linsolve (&__linsolve,_QUOTE_ARGUMENTS,true);

  const string _resoudre_systeme_lineaire_s("resoudre_systeme_lineaire");
  unary_function_eval __resoudre_systeme_lineaire(&_linsolve,_resoudre_systeme_lineaire_s);
  unary_function_ptr at_resoudre_systeme_lineaire (&__resoudre_systeme_lineaire,_QUOTE_ARGUMENTS,true);

  /*
  gen iter(const gen & f, const gen & x,const gen & arg,int maxiter,double eps, int & b){
    gen a=arg;
    complex<double> olda;
    complex<double> ad;
    b=0;
    ad=gen2complex_d(a);
    //cout<<"a"<<a<<endl;
    //cout<<"ad"<<ad<<endl;
    for (int j=1;j<=maxiter;j++){
      olda=ad;    
      // cout << f << " " << x << " " << a << endl;
      a=subst(f,x,a).evalf();
      // cout<<"a"<<a<<endl;
      //a=a.evalf();
      //ad=a._DOUBLE_val;
      ad=gen2complex_d(a);
      // cout<<"a"<<a<<endl;
      // cout<<"ad"<<ad<<endl;
      // cout<<"j"<<j<<"abs"<<abs(ad-olda)<<endl;
      if (eps>abs(ad-olda)) {
	b=1; return(a);
      }
    } 
    return(a); 
  }
  
  gen newtona(const gen & f, const gen & x, const gen & arg,int niter1, int niter2, double eps1,double eps2,double prefact1,double prefact2, int & b){
    if (x.type!=_IDNT)
      settypeerr("2nd arg must be an identifier");
    //cout<<a<<endl;
    gen a=arg;
    gen g1;
    gen g;
    g1=x-gen(prefact1)*rdiv(f,derive(f,x));
    // sym_sub(x,sym_mult(rdiv(9,10),rdiv(f,derive(f,x))));
    try {
      a= iter(g1,x,a,niter1,eps1,b);
      g=x-gen(prefact2)*rdiv(f,derive(f,x));
      a= iter(g,x,a,niter2,eps2,b); 
    }
    catch (std::runtime_error & err){
      b=0;
    }
    return a;
  }
 
  gen newton(const gen & f, const gen & x,const gen & guess,int niter1,int niter2,double eps1,double eps2,double prefact1,double prefact2){
    bool guess_first=is_undef(guess);
    for (int j=1;j<5;j++,niter2 *=2, niter1 *=2){ 
      gen a;
      int b;
      //on prend un départ au hasard (a=x0=un _DOUBLE_)
      // a=gen(2.0);
      if (guess_first)
	a=j*4*(rand()/(RAND_MAX+1.0)-0.5);
      else {
	a=guess;
	guess_first=true;
      }
      // cout<<j<<"j"<<a<<endl; 
      gen e;
      e=newtona(f, x, a,niter1,niter2,eps1,eps2,prefact1,prefact2,b);
      if (b==1) return e;
      gen c;
      c=j*4*(rand()/(RAND_MAX+1.0)-0.5);
      // cout<<j<<"j"<<c<<endl;
      // g=x-gen(0.5)*rdiv(f,derive(f,x));
      gen ao(gen(0.0),c);
      // cout<<"ao"<<ao<<endl;
      gen e0= newtona(f, x, ao,niter1,niter2,eps1,eps2,prefact1,prefact2,b);
      if (b==1) 
	return(e0);
      gen a1(a,c);
      // cout<<j<<"j,a1"<<a1<<endl;
      e0= newtona(f, x, a1,niter1,niter2,eps1,eps2,prefact1,prefact2,b);
      if (b==1) 
	return(e0);
    }
    setsizeerr("nontrouve");
    return(0);
  }
  */

  gen newton_rand(int j){
    gen a=evalf(j*4*(gen(rand())/(gen(RAND_MAX)+1)-plus_one_half),1,0);
    if (j!=1)
      a=a+cst_i*evalf(j*4*(gen(rand())/(gen(RAND_MAX)+1)-plus_one_half),1,0);
    return a;
  }

  gen newton(const gen & f, const gen & x,const gen & guess_,int niter,double eps1, double eps2,GIAC_CONTEXT){
    gen guess(guess_);
    if (guess.is_symb_of_sommet(at_interval))
      guess=(guess._SYMBptr->feuille[0]+guess._SYMBptr->feuille[1])/2;
    gen a,b,d,fa,fb,invdf=inv(derive(f,x,contextptr),contextptr),epsg1(eps1),epsg2(eps2);
    if (ckmatrix(invdf))
      invdf=mtran(*invdf._VECTptr);
    bool guess_first=is_undef(guess);
    // Main loop with random initialization
    int j=1;
    for (;j<=5 ;j++,niter=2*niter){ 
      if (guess_first){
	if (f.type==_VECT){
	  int s=f._VECTptr->size();
	  vecteur v(s);
	  for (int i=0;i<s;++i)
	    v[i]=(newton_rand(j));
	  a=v;
	}
	else
	  a=newton_rand(j);
      }
      else {
	a=guess;
	guess_first=true;
      }
      fa=evalf(eval(subst(f,x,a,false,contextptr),eval_level(contextptr),contextptr),1,contextptr); 
      // First loop to localize the solution with prefactor
      gen lambda(1);
      int k;
      try {
	for (k=0;k<niter;++k){
	  d=-evalf(eval(subst(invdf,x,a,false,contextptr),eval_level(contextptr),contextptr)*fa,1,contextptr);
	  b=a+lambda*d;
	  gen babs=_l2norm(b,contextptr);
	  if (is_inf(babs) || is_undef(babs)){
	    guess_first=true;
	    k=niter;
	    break;
	  }
	  if (is_positive(epsg1-_l2norm(d,contextptr),contextptr)){
	    a=b;
	    break;
	  }
	  fb=evalf(eval(subst(f,x,b,false,contextptr),eval_level(contextptr),contextptr),1,contextptr);
	  if (is_positive(_l2norm(fb,contextptr)-_l2norm(fa,contextptr),contextptr)){
	    // Decrease prefactor and try again
	    lambda=evalf(plus_one_half,1,0)*lambda;
	  }
	  else {
	    // Save new value of a and increase the prefactor slightly
	    if (is_positive(lambda-0.9,contextptr))
	      lambda=1;
	    else
	      lambda=evalf(gen(12)/gen(10),1,contextptr)*lambda;
	    a=b;
	    fa=fb;
	  }
	}
	if (k==niter)
	  continue;
	// Second loop to improve precision (prefactor 1)
	for (k=0;k<niter;++k){
	  d=-evalf(subst(invdf,x,a,false,contextptr)*subst(f,x,a,false,contextptr),1,contextptr);
	  a=a+d;
	  if (is_positive(epsg2-_l2norm(d,contextptr),contextptr))
	    break;
	}
	if (k!=niter)
	  break;
      } catch (std::runtime_error & e){
	break; // start with a new initial point
      }
    } // end for
    if (j>5)
      return undef;
    return a;
  }

  gen _newton(const gen & args,GIAC_CONTEXT){
    double gsl_eps=epsilon(contextptr);
    if (args.type!=_VECT)
      return newton(args,vx_var,undef,NEWTON_DEFAULT_ITERATION,gsl_eps,1e-12,contextptr);
    vecteur v=*args._VECTptr;
    int s=v.size();
    if (s<2)
      setsizeerr();
    if (s==2){
      if (v[1].is_symb_of_sommet(at_equal))
	return newton(v[0],v[1]._SYMBptr->feuille[0],v[1]._SYMBptr->feuille[1],NEWTON_DEFAULT_ITERATION,gsl_eps,1e-12,contextptr);
      return newton(v[0],v[1],undef,NEWTON_DEFAULT_ITERATION,gsl_eps,1e-12,contextptr);
    }
    int niter=NEWTON_DEFAULT_ITERATION;
    double eps=epsilon(contextptr);
    for (int j=3;j<s;++j){
      if (v[j].type==_INT_)
	niter=v[j].val;
      else {
	gen tmp=evalf_double(v[j],1,contextptr);
	if (tmp.type==_DOUBLE_)
	  eps=tmp._DOUBLE_val;
      }
    }
    gen res=newton(v[0],v[1],v[2],niter,1e-10,eps,contextptr);
    if (debug_infolevel)
      *logptr(contextptr) << res << endl;
    return res;
    toomanyargs("newton");
  }
  const string _newton_s("newton");
  unary_function_eval __newton(&_newton,_newton_s);
  unary_function_ptr at_newton (&__newton,0,true);
  
  bool has_num_coeff(const vecteur & v){
    const_iterateur it=v.begin(),itend=v.end();
    for (;it!=itend;++it){
      if (has_num_coeff(*it))
	return true;
    }
    return false;
  }
  
  bool has_num_coeff(const polynome & p){
    vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    for (;it!=itend;++it){
      if (has_num_coeff(it->value))
	return true;
    }
    return false;
  }

  bool has_num_coeff(const gen & e){
    switch (e.type){
    case _ZINT: case _INT_: case _IDNT: case _USER:
      return false;
    case _DOUBLE_: case _REAL:
      return true;
    case _CPLX:
      return (e._CPLXptr->type==_DOUBLE_) || ((e._CPLXptr+1)->type==_DOUBLE_);
    case _SYMB:
      return has_num_coeff(e._SYMBptr->feuille);
    case _VECT:
      return has_num_coeff(*e._VECTptr);
    case _POLY:
      return has_num_coeff(*e._POLYptr);
    case _FRAC:
      return has_num_coeff(e._FRACptr->num) || has_num_coeff(e._FRACptr->den);
    default:
      return false;
    }
    return 0;
  }

  bool has_mod_coeff(const vecteur & v,gen & modulo){
    const_iterateur it=v.begin(),itend=v.end();
    for (;it!=itend;++it){
      if (has_mod_coeff(*it,modulo))
	return true;
    }
    return false;
  }

  bool has_mod_coeff(const polynome & p,gen & modulo){
    vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    for (;it!=itend;++it){
      if (has_mod_coeff(it->value,modulo))
	return true;
    }
    return false;
  }

  bool has_mod_coeff(const gen & e,gen & modulo){
    switch (e.type){
    case _MOD:
      modulo = *(e._MODptr+1);
      return true;
    case _SYMB:
      return has_mod_coeff(e._SYMBptr->feuille,modulo);
    case _VECT:
      return has_mod_coeff(*e._VECTptr,modulo);
    case _POLY:
      return has_mod_coeff(*e._POLYptr,modulo);
    default:
      return false;
    }
  }

  polynome spoly(const polynome & p,const polynome & q,environment * env){
    if (p.coord.empty() || q.coord.empty())
      setsizeerr();
    index_t & pi = *p.coord.front().index.iptr;
    index_t & qi = *q.coord.front().index.iptr;
    index_t lcm= index_lcm(pi,qi);
    polynome tmp=p.shift(lcm-pi,q.coord.front().value)-q.shift(lcm-qi,p.coord.front().value);
    return (env && env->moduloon)?smod(tmp,env->modulo):tmp;
  }


  bool GroebnerDiv (const polynome & celuici,const polynome & other, polynome & quo, polynome & rem, bool allowrational = true ) {  
    if (celuici.coord.empty()){
      quo=celuici;
      rem=celuici; 
      return true;
    }
    assert(!other.coord.empty());
    quo.coord.clear();
    rem=celuici;
    const index_m & b_max = other.coord.front().index;
    const gen & b=other.coord.front().value;
    std::vector< monomial<gen> >::const_iterator it,itend;
    for (;;){
      itend=rem.coord.end();
      it=rem.coord.begin();
      // look in rem for a monomial >= b_max
      for (; it!=itend ;++it){
	if (it->index>=b_max)
	  break;
      }
      if (it==itend) // no monomial of rem are divisible by LT(b): finished
	break;
      gen q=rdiv(it->value,b);
      if (!allowrational){ // quotient not allowed
	if ( has_denominator(q) || (!is_zero(q*b - it->value)) )
	  return false;
      }
      quo.coord.push_back(monomial<gen>(q,it->index-b_max));
      polynome temp=other.shift(it->index-b_max,q);
      rem = rem-temp;
    } // end for
    return(true);    
  }


  template <class T>
  void polynome2map(const polynome & p,map<index_t,gen,T> & pmap){
    std::vector< monomial<gen> >::const_iterator pt=p.coord.begin(),ptend=p.coord.end();
    pmap.clear();
    for (;pt!=ptend;++pt){
      pmap[*pt->index.iptr]=pt->value;
    }
  }

  // find integer content and divide
  template <class T>
  gen mapppz(map<index_t,gen,T> & pmap){
    typename std::map< index_t,gen,T>::iterator it=pmap.begin(),itend=pmap.end();
    gen d;
    for (;it!=itend;++it){
      d=gcd(d,it->second);
      if (is_one(d))
	return plus_one;
    }
    for (it=pmap.begin();it!=itend;++it)
      it->second=it->second/d;
    return d;
  }

  polynome reduce2(const polynome & p,vectpoly::const_iterator it0,vectpoly::const_iterator itend){
    if (p.coord.empty())
      return p;
    typedef std::map< index_t,gen,const std::pointer_to_binary_function < const index_t &, const index_t &, bool> > application;
    application rem(p.is_strictly_greater);
    polynome2map(p,rem);
    application::const_iterator pt,ptend;
    vectpoly::const_iterator it;
    for (;;){
      ptend=rem.end();
      // look in rem for a monomial >= to a monomial in it0, then it0+1 
      for (it=it0; it!=itend ;++it){
	for (pt=rem.begin();pt!=ptend;++pt){
	  if (pt->first>=it->coord.front().index)
	    break;
	}
	if (pt!=ptend)
	  break;
      }
      if (it==itend) // no monomial of rem are divisible by LT(b): finished
	break;
      gen a(pt->second),b(it->coord.front().value) ;
      simplify(a,b);
      if (is_minus_one(b)){
	a=-a;
	b=1;
      }
      polynome temp=it->shift(pt->first-it->coord.front().index,a);
      application::iterator pit;
      if (!is_one(b)){
	for (pit=rem.begin();pit!=ptend;++pit)
	  pit->second = b*pit->second;
      }
      // substract temp from rem
      vector< monomial<gen> >::const_iterator jt=temp.coord.begin(),jtend=temp.coord.end();
      for (;jt!=jtend;++jt){
	pit=rem.find(*jt->index.iptr);
	if (pit==ptend)
	  rem[*jt->index.iptr]=-jt->value;
	else {
	  pit->second -= jt->value;
	  if (is_zero(pit->second)){
	    rem.erase(pit->first);
	    ptend=rem.end();
	  }
	}
      }
      if (!is_one(b))
	mapppz(rem);
    }
    // convert back
    polynome res(p.dim,p);
    pt=rem.begin();
    res.coord.reserve(rem.size());
    for (;pt!=ptend;++pt){
      res.coord.push_back(monomial<gen>(pt->second,pt->first));
    }
    gen d=ppz(res);
    if (!is_one(d)) cerr << d << endl;
    return res;
  }

  polynome reducegb(const polynome & p,vectpoly::const_iterator it0,vectpoly::const_iterator itend,environment * env){
    if (p.coord.empty())
      return p;
    polynome rem(p),quo,tmp;
    for (;it0!=itend;++it0){
      if (!GroebnerDiv(rem,*it0,quo,tmp,true))
	setsizeerr();
      rem=tmp;
      ppz(rem);
    }
    return rem;
  }

  polynome reduce(const polynome & p,vectpoly::const_iterator it0,vectpoly::const_iterator itend,environment * env){
    if (p.coord.empty())
      return p;
    polynome rem(p);
    std::vector< monomial<gen> >::const_iterator pt,ptend;
    vectpoly::const_iterator it;
    for (;;){
      ptend=rem.coord.end();
      // look in rem for a monomial >= to a monomial in it0, then it0+1 
      for (it=it0; it!=itend ;++it){
	for (pt=rem.coord.begin();pt!=ptend;++pt){
	  if (pt->index>=it->coord.front().index)
	    break;
	}
	if (pt!=ptend)
	  break;
      }
      if (it==itend) // no monomial of rem are divisible by LT(b): finished
	break;
      gen a(pt->value),b(it->coord.front().value) ;
      if (env && env->moduloon){
	polynome temp=it->shift(pt->index-it->coord.front().index,a*invmod(b,env->modulo));
	rem = smod(rem - temp,env->modulo) ; // FIXME: improve!
      }
      else {
	simplify(a,b);
	polynome temp=it->shift(pt->index-it->coord.front().index,a);
	if (is_one(b))
	  rem = rem-temp;
	else {
	  rem = b*rem - temp;
	  ppz(rem);
	}
      }
    }
    if (env && env->moduloon)
      ;
    else
      ppz(rem);
    return rem;
  }

  polynome reduce(const polynome & p,const vectpoly & v,environment * env){
    vectpoly::const_iterator it=v.begin(),itend=v.end();
    return reduce(p,it,itend,env);
  }

  void reduce(vectpoly & res,environment * env){
    if (res.empty())
      return;
    sort(res.begin(),res.end(),tensor_is_greater<gen>);
    // reduce res
    int s=res.size();
    for (int i=0;i<s-1;){
      polynome & p=res[i];
      polynome pred=reduce(p,res.begin()+i+1,res.end(),env);
      if (pred.coord.empty()){
	res.erase(res.begin()+i);
	--s;
	continue;
      }
      if (pred.coord.size()==p.coord.size() && pred*p.coord.front().value==p*pred.coord.front().value){
	++i;
	continue;
      }
      res[i]=pred;
      sort(res.begin()+i,res.end(),tensor_is_greater<gen>);
      i=0;
    }
  }

  vectpoly gbasis(const vectpoly & v,const gen & order,bool with_cocoa,bool with_f5,environment * env){
    if (v.size()==1){
      return v;
    }
    vectpoly res(v);
    try {
      if (with_cocoa){
	bool ok=with_f5?f5(res,order):cocoa_gbasis(res,order);
	if (ok){
	  if (debug_infolevel)
	    res.dbgprint();
	  return res;
	}
      }
    } catch (...){
      cerr << "Unable to compute gbasis with CoCoA" << endl;
    }
    reduce(res,env);
    bool notfound=true;
    for (;notfound;){
      // cerr << res << endl;
      notfound=false;
      vectpoly::const_iterator it=res.begin(),itend=res.end(),jt;
      vectpoly newres(res);
      for (;it!=itend;++it){
	for (jt=it+1;jt!=itend;++jt){
	  polynome toadd(spoly(*it,*jt,env));
	  toadd=reduce(toadd,newres,env);
	  if (!toadd.coord.empty()){
	    newres.push_back(toadd); // should be at the right place
	    notfound=true;
	  }
	}
      }
      reduce(newres,env);
      res=newres;
    }
    return res;
  }

  gen in_ideal(const vectpoly & r,const vectpoly & v,const gen & order,bool with_cocoa,bool with_f5,environment * env){
    try {
      if (with_cocoa){
	return cocoa_in_ideal(r,v,order);
      }
    } catch (...){
     return -1;
    }
    return -1;
  }

  gen remove_equal(const gen & f){
    if ( (f.type==_SYMB) && (f._SYMBptr->sommet==at_equal || f._SYMBptr->sommet==at_same ) ){
      vecteur & v=*f._SYMBptr->feuille._VECTptr;
      return v.front()-v.back();
    }
    if (f.type==_VECT)
      return apply(f,remove_equal);
    return f;
  }

  vecteur remove_equal(const_iterateur it,const_iterateur itend){
    vecteur conditions;
    conditions.reserve(itend-it);
    for (;it!=itend;++it){
	conditions.push_back(remove_equal(*it));
    }
    return conditions;
  }

  bool vecteur2vector_polynome(const vecteur & eq_in,const vecteur & l,vectpoly & eqp){
    // remove all denominators
    const_iterateur it=eq_in.begin(),itend=eq_in.end();
    for (;it!=itend;++it){
      gen n,d;
      fxnd(*it,n,d);
      if (n.type==_POLY){
	// should reordre n with total degree+revlex order here
	eqp.push_back(*n._POLYptr);
	continue;
      }
      if (!is_zero(n))
	return false;
    }
    return true;
  }

  vecteur gsolve(const vecteur & eq_orig,const vecteur & var_orig,bool complexmode,GIAC_CONTEXT){
    // replace variables in var_orig by true identificators
    vecteur var(var_orig);
    iterateur it=var.begin(),itend=var.end();
    int s=itend-it; // # of unknowns
    bool need_subst=false;
    vector<identificateur> tab_idnt(s);
    for (int i=0;it!=itend;++it,++i){
      if (it->type!=_IDNT){
	*it=tab_idnt[i]; 
	need_subst=true;
      }
    }
    vecteur eq(remove_equal(eq_orig.begin(),eq_orig.end()));
    if (need_subst)
      eq=subst(eq,var_orig,var,false,contextptr);
    if (approx_mode(contextptr)){
#ifdef HAVE_LIBGSL
      return makevecteur(msolve(eq,var,multvecteur(zero,var),_HYBRID_SOLVER,epsilon(contextptr),contextptr));
#else
      return vecteur(1,undef);
#endif
    }
    bool convertapprox=has_num_coeff(eq);
    if (convertapprox)
      eq=*exact(evalf(eq,1,contextptr),contextptr)._VECTptr;
    // check rational
    for (it=var.begin();it!=itend;++it){
      if (it->type!=_IDNT) // should not occur!
	setsizeerr("Bad var "+it->print(contextptr));
      vecteur l(rlvarx(eq,*it->_IDNTptr));
      if (l.size()>1)
	return vecteur(1,string2gen(gen(l).print(contextptr)+" is not rational w.r.t. "+it->print(contextptr),false));
    }
    vecteur l(1,var);
    alg_lvar(eq,l);
    // convert eq to polynomial
    vecteur eq_in(*e2r(eq,l,contextptr)._VECTptr);
    vectpoly eqp;
    // remove all denominators
    it=eq_in.begin();
    itend=eq_in.end();
    for (;it!=itend;++it){
      gen n,d;
      fxnd(*it,n,d);
      if (n.type==_POLY){
	// should reordre n with total degree+revlex order here
	eqp.push_back(*n._POLYptr);
	continue;
      }
      if (!is_zero(n))
	return vecteur(0); // no solution since cst equation
    }
    vectpoly eqpr(gbasis(eqp,_PLEX_ORDER));
    // should reorder eqpr with lex order here
    // solve from right to left
    reverse(eqpr.begin(),eqpr.end());
    vecteur sols(1,vecteur(0)); // sols=[ [] ]
    vectpoly::const_iterator jt=eqpr.begin(),jtend=eqpr.end();
    for (;jt!=jtend;++jt){
      // the # of found vars is the size of sols.front()
      if (sols.empty())
	break;
      vecteur newsols;
      gen g(r2e(*jt,l,contextptr));
      const_iterateur st=sols.begin(),stend=sols.end();
      for (;st!=stend;++st){
	int foundvars=st->_VECTptr->size();
	vecteur current=*st->_VECTptr;
	gen curg=ratnormal(ratnormal(subst(g,vecteur(var.end()-foundvars,var.end()),*st,false,contextptr)));
	gen x;
	int xpos=0;
	// First search in current an identifier curg depends on
	for (;xpos<foundvars;++xpos){
	  x=current[xpos];
	  if (x==var[s-foundvars+xpos] && !is_zero(derive(curg,x,contextptr)) )
	    break;
	}
	if (xpos==foundvars){
	  xpos=0;
	  // find next var g depends on 
	  for (;foundvars<s;++foundvars){
	    x=var[s-foundvars-1];
	    current.insert(current.begin(),x);
	    if (!is_zero(derive(curg,x,contextptr)))
	      break;
	  }
	  if (s==foundvars){
	    if (is_zero(curg))
	      newsols.push_back(current);
	    continue;
	  }
	}
	// solve
	vecteur xsol(solve(curg,*x._IDNTptr,complexmode,contextptr));
	const_iterateur xt=xsol.begin(),xtend=xsol.end();
	for (;xt!=xtend;++xt){
	  // current[xpos]=*xt;
	  newsols.push_back(subst(current,*x._IDNTptr,*xt,false,contextptr));
	}
      } // end for (;st!=stend;)
      sols=newsols;
    }
    // Add var at the beginning of each solution of sols if needed
    it=sols.begin(); 
    itend=sols.end();
    for (;it!=itend;++it){
      int ss=it->_VECTptr->size();
      if (ss<s)
	*it=mergevecteur(vecteur(var.begin(),var.begin()+s-ss),*it->_VECTptr);
    }
    if (need_subst)
      sols=subst(sols,var,var_orig,false,contextptr);
    if (convertapprox)
      sols=*evalf_VECT(sols,0,1,contextptr)._VECTptr;
    return sols;
  }

  void read_gbargs(const vecteur & v,int start,int s,gen & order,bool & with_cocoa,bool & with_f5){
    for (int i=start;i<s;++i){
      if (v[i].is_symb_of_sommet(at_equal)){
	gen & tmp=v[i]._SYMBptr->feuille;
	if (tmp.type==_VECT && tmp._VECTptr->front().type==_INT_ && tmp._VECTptr->back().type==_INT_){
	  switch (tmp._VECTptr->front().val){
	  case _WITH_COCOA:
	    with_cocoa=tmp._VECTptr->back().val;
	    break;
	  case _WITH_F5:
	    with_f5=tmp._VECTptr->back().val;
	    break;
	  }
	}
      }
      if (v[i].type==_INT_ && v[i].subtype==_INT_GROEBNER){
	switch (v[i].val){
	case _WITH_COCOA:
	  with_cocoa=true;
	  break;
	case _WITH_F5:
	  with_f5=true;
	  break;
	default:
	  order=v[i].val;
	}
      }
    }
  }


  void change_monomial_order(polynome & p,const gen & order){
    switch (order.val){
      // should be strict, but does not matter since monomials are !=
    case _REVLEX_ORDER: 
      p.is_strictly_greater=std::ptr_fun(total_revlex_is_greater<int>);
      p.m_is_greater=std::ptr_fun(m_total_revlex_is_greater<gen>);
      break;
    case _TDEG_ORDER:
      p.is_strictly_greater=std::ptr_fun(total_lex_is_greater<int>);
      p.m_is_greater=std::ptr_fun(m_total_lex_is_greater<gen>);
      break;
    }
    p.tsort();
  }

  void change_monomial_order(vectpoly & eqp,const gen & order){
    // change polynomial order
    if (order.type==_INT_ && order.val){
      vectpoly::iterator it=eqp.begin(),itend=eqp.end();
      for (;it!=itend;++it){
	change_monomial_order(*it,order);
      }
    }
  }

  // gbasis([Pi],[vars]) -> [Pi']
  gen _gbasis(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      return symbolic(at_gbasis,args);
    vecteur & v = *args._VECTptr;
    int s=v.size();
    if (s<2)
      toofewargs("gbasis");
    if ( (v[0].type!=_VECT) || (v[1].type!=_VECT) )
      setsizeerr();
    vecteur l=vecteur(1,v[1]);
    alg_lvar(v[0],l);
    // v[2] will serve for ordering
    gen order=_PLEX_ORDER; // _REVLEX_ORDER;
    bool with_f5=false,with_cocoa=true;
    read_gbargs(v,2,s,order,with_cocoa,with_f5);
    // convert eq to polynomial
    vecteur eq_in(*e2r(v[0],l,contextptr)._VECTptr);
    vectpoly eqp;
    if (!vecteur2vector_polynome(eq_in,l,eqp))
      return vecteur(1,plus_one);
    gen coeff;
    environment env ;
    if (!eqp.empty() && coefftype(eqp.front(),coeff)==_MOD){
      with_cocoa = false;
      env.moduloon = true;
      env.modulo = *(coeff._MODptr+1);
      env.pn=env.modulo;
      vectpoly::iterator it=eqp.begin(),itend=eqp.end();
      for (;it!=itend;++it)
	*it=unmodularize(*it);
    }
    else
      env.moduloon = false;
    if (!with_cocoa)
      change_monomial_order(eqp,order);
    vectpoly eqpr(gbasis(eqp,order,with_cocoa,with_f5,&env));
    vecteur res;
    vectpoly::const_iterator it=eqpr.begin(),itend=eqpr.end();
    res.reserve(itend-it);
    for (;it!=itend;++it)
      res.push_back(r2e(*it,l,contextptr));
    return res;
  }
  const string _gbasis_s("gbasis");
  unary_function_eval __gbasis(&_gbasis,_gbasis_s);
  unary_function_ptr at_gbasis (&__gbasis,0,true);
  
  // greduce(P,[gbasis],[vars])
  gen _greduce(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      return symbolic(at_gbasis,args);
    vecteur & v = *args._VECTptr;
    int s=v.size();
    if (s<3)
      toofewargs("greduce");
    if ( (v[1].type!=_VECT) || (v[2].type!=_VECT) )
      setsizeerr();
    vecteur l=vecteur(1,v[2]);
    alg_lvar(v[0],l);
    alg_lvar(v[1],l);
    // v[3] will serve for ordering
    gen order=_PLEX_ORDER; // _REVLEX_ORDER;
    bool with_f5=false,with_cocoa=true;
    read_gbargs(v,3,s,order,with_cocoa,with_f5);
    gen eq(e2r(v[0],l,contextptr));
    if (eq.type!=_POLY)
      return v[0];
    gen coeff;
    environment env ;
    if (coefftype(*eq._POLYptr,coeff)==_MOD){
      with_cocoa = false;
      env.moduloon = true;
      env.modulo = *(coeff._MODptr+1);
      env.pn=env.modulo;
    }
    else
      env.moduloon = false;
    vecteur eq_in(*e2r(v[1],l,contextptr)._VECTptr);
    vectpoly eqp;
    if (!vecteur2vector_polynome(eq_in,l,eqp))
      return zero;
    change_monomial_order(eqp,order);
    polynome p(*eq._POLYptr);
    change_monomial_order(p,order);
    vectpoly rescocoa;
    if (!env.moduloon && cocoa_greduce(vectpoly(1,p),eqp,order,rescocoa))
      return r2e(rescocoa.front(),l,contextptr);
    gen C(eq._POLYptr->constant_term());
    // FIXME: get constant term, substract one to get the correct constant
    eq=eq-C+plus_one;
    // polynome res(env.moduloon?reduce(p,eqp.begin(),eqp.end(),&env):reducegb(p,eqp.begin(),eqp.end(),&env));
    polynome res(reduce(p,eqp.begin(),eqp.end(),&env));
    gen C1(res.constant_term());
    if (env.moduloon){
      res=invmod(C1,env.modulo)*res;
      modularize(res,env.modulo);
    }
    else
      res=res/C1;
    return r2e(res-plus_one,l,contextptr)+C;
  }
  const string _greduce_s("greduce");
  unary_function_eval __greduce(&_greduce,_greduce_s);
  unary_function_ptr at_greduce (&__greduce,0,true);
  
  // in_ideal([Pi],[gb],[vars]) -> true/false
  gen _in_ideal(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      setsizeerr();
    vecteur & v = *args._VECTptr;
    int s=v.size();
    if (s<3)
      toofewargs("in_ideal");
    if ( v[1].type!=_VECT || v[2].type!=_VECT )
      setsizeerr();
    vecteur atester=gen2vecteur(v[0]);
    vecteur l=vecteur(1,v[2]);
    alg_lvar(v[1],l);
    alg_lvar(v[0],l);
    gen order=_PLEX_ORDER; // _REVLEX_ORDER;
    bool with_f5=false,with_cocoa=true;
    read_gbargs(v,3,s,order,with_cocoa,with_f5);
    // convert eq to polynomial
    vecteur eq_in(*e2r(v[1],l,contextptr)._VECTptr);
    vecteur r(*e2r(atester,l,contextptr)._VECTptr);
    vectpoly eqp,eqr;
    if (!vecteur2vector_polynome(eq_in,l,eqp) || !vecteur2vector_polynome(r,l,eqr))
      setsizeerr();
    gen coeff;
    environment env ;
    if (!eqp.empty() && coefftype(eqp.front(),coeff)==_MOD){
      with_cocoa = false;
      env.moduloon = true;
      env.modulo = *(coeff._MODptr+1);
      env.pn=env.modulo;
      vectpoly::iterator it=eqp.begin(),itend=eqp.end();
      for (;it!=itend;++it)
	*it=unmodularize(*it);
    }
    else
      env.moduloon = false;
    if (!with_cocoa){
      change_monomial_order(eqp,order);
      change_monomial_order(eqr,order);
    }
    // is r in ideal eqp?
    gen res=in_ideal(eqr,eqp,order,with_cocoa,with_f5,&env);
    if (res.type==_VECT && !res._VECTptr->size()==1 && v[0].type!=_VECT)
      return res._VECTptr->front();
    return res;
  }
  const string _in_ideal_s("in_ideal");
  unary_function_eval __in_ideal(&_in_ideal,_in_ideal_s);
  unary_function_ptr at_in_ideal (&__in_ideal,0,true);

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC
