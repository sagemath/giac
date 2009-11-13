/* -*- mode:C++ ; compile-command: "g++-3.4 -I.. -g -c derive.cc" -*- */
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
#include "derive.h"
#include "usual.h"
#include "symbolic.h"
#include "unary.h"
#include "poly.h"
#include "sym2poly.h" // for equalposcomp
#include "tex.h"
#include "prog.h"
#include "intg.h"
#include "subst.h"
#include "plot.h"
#include "modpoly.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  bool depend(const gen & g,const identificateur & i){
    if (g.type==_IDNT)
      return *g._IDNTptr==i;
    if (g.type==_SYMB)
      return depend(g._SYMBptr->feuille,i);
    if (g.type!=_VECT)
      return false;
    const_iterateur it=g._VECTptr->begin(),itend=g._VECTptr->end();
    for (;it!=itend;++it){
      if (depend(*it,i))
	return true;
    }
    return false;
  }

  gen derive_SYMB(const symbolic & s,const identificateur & i,GIAC_CONTEXT){
    // if s does not depend on i return 0
    if (!depend(s,i))
      return zero;
    // rational operators are treated first for efficiency
    if (s.sommet==at_plus){
      vecteur v;
      if (s.feuille.type!=_VECT)
	return derive(s.feuille,i,contextptr);
      int taille=s.feuille._VECTptr->size();
      v.reserve(taille);
      vecteur::const_iterator iti=s.feuille._VECTptr->begin(),itend=s.feuille._VECTptr->end();
      gen e;
      for (;iti!=itend;++iti){
	e=derive(*iti,i,contextptr);
	if (!is_zero(e))
	  v.push_back(e);
      }
      if (v.size()==1)
	return v.front();
      if (v.empty())
	return zero;
      return symbolic(at_plus,v);
    }
    if (s.sommet==at_prod){
      vecteur v,w;
      if (s.feuille.type!=_VECT)
	return derive(s.feuille,i,contextptr);
      int taille=s.feuille._VECTptr->size();
      v.reserve(taille);
      w.reserve(taille);
      vecteur::const_iterator itbegin=s.feuille._VECTptr->begin(),itj,iti,itend=s.feuille._VECTptr->end();
      gen e;
      for (iti=itbegin;iti!=itend;++iti){
	w.clear();
	e=derive(*iti,i,contextptr);
	if (!is_zero(e)){
	  for (itj=itbegin;itj!=iti;++itj)
	    w.push_back(*itj);
	  w.push_back(e);
	  ++itj;
	  for (;itj!=itend;++itj)
	    w.push_back(*itj);
	  v.push_back(_prod(w,contextptr));
	}
      }
      if (v.size()==1)
	return v.front();
      if (v.empty())
	return zero;
      return symbolic(at_plus,v);
    }
    if (s.sommet==at_neg)
      return -derive(s.feuille,i,contextptr);
    if (s.sommet==at_inv){
      if (s.feuille.is_symb_of_sommet(at_pow)){
	gen & f = s.feuille._SYMBptr->feuille;
	if (f.type==_VECT && f._VECTptr->size()==2)
	  return derive(symb_pow(f._VECTptr->front(),-f._VECTptr->back()),i,contextptr);
      }
      return rdiv(-derive(s.feuille,i,contextptr),pow(s.feuille,2));
    }
    if (s.sommet==at_ln){ 
      if (s.feuille.is_symb_of_sommet(at_abs) )
	return rdiv(derive(s.feuille._SYMBptr->feuille,i,contextptr),s.feuille._SYMBptr->feuille);
      if (s.feuille.is_symb_of_sommet(at_inv))
	return -derive(symbolic(at_ln,s.feuille._SYMBptr->feuille),i,contextptr);
      if (s.feuille.is_symb_of_sommet(at_prod)){
	gen res;
	const gen &f=s.feuille._SYMBptr->feuille;
	if (f.type==_VECT){
	  const_iterateur it=f._VECTptr->begin(),itend=f._VECTptr->end();
	  for (;it!=itend;++it)
	    res=res+derive(symbolic(at_ln,*it),i,contextptr);
	  return res;
	}
      }
    }
    if (s.feuille.type==_VECT){
      vecteur v=*s.feuille._VECTptr;
      int vs=v.size();
      if (vs>=3 && s.sommet==at_ifte || s.sommet==at_when ){
	for (int j=1;j<vs;++j)
	  v[j]=derive(v[j],i,contextptr);
	return symbolic(s.sommet,v);
      }
      if (s.sommet==at_piecewise){
	for (int j=0;j<vs/2;++j)
	  v[2*j+1]=derive(v[2*j+1],i,contextptr);
	if (vs%2)
	  v[vs-1]=derive(v[vs-1],i,contextptr);
	return symbolic(s.sommet,v);
      }
    }
    // now look at other operators, first onearg operator
    if (s.sommet.ptr->D){
      if (s.feuille.type!=_VECT)
	return (*s.sommet.ptr->D)(1)(s.feuille,0)*derive(s.feuille,i,contextptr);
      // multiargs operators
      int taille=s.feuille._VECTptr->size();
      vecteur v;
      v.reserve(taille);
      vecteur::const_iterator iti=s.feuille._VECTptr->begin(),itend=s.feuille._VECTptr->end();
      gen e;
      for (int j=1;iti!=itend;++iti,++j){
	e=derive(*iti,i,contextptr);
	if (!is_zero(e))
	  v.push_back(e*(*s.sommet.ptr->D)(j)(s.feuille,0));
      }
      if (v.size()==1)
	return v.front();
      if (v.empty())
	return zero;
      return symbolic(at_plus,v);
    }
    // integrate
    if (s.sommet==at_integrate){
      if (s.feuille.type!=_VECT)
	return s.feuille;
      vecteur & v=*s.feuille._VECTptr;
      int nargs=v.size();
      if (nargs<=1)
	return s.feuille;
      gen res,newint;
      if (v[1]==i)
	res=v[0];
      else {
	res=subst(v[0],v[1],i,false,contextptr);
	newint=derive(v[0],i,contextptr);	 
	if (nargs<4)
	  newint=integrate(newint,v[1],contextptr);
      }
      if (nargs==2)
	return res+newint;
      if (nargs==3)
	return derive(v[2],i,contextptr)*subst(res,i,v[2],false,contextptr);
      if (nargs==4){
	gen a3=derive(v[3],i,contextptr);
	gen b3=is_zero(a3)?zero:limit(res,i,v[3],-1,contextptr);
	gen a2=derive(v[2],i,contextptr);
	gen b2=is_zero(a2)?zero:limit(res,i,v[2],1,contextptr);
	return a3*b3-a2*b2+_integrate(gen(makevecteur(newint,v[1],v[2],v[3]),_SEQ__VECT),contextptr);
      }
      setsizeerr();
    }
    // multi derivative and multi-indice derivatives
    if (s.sommet==at_derive){
      if (s.feuille.type!=_VECT)
	return symbolic(at_derive,gen(makevecteur(s.feuille,vx_var,2),_SEQ__VECT));
      if (s.feuille._VECTptr->size()==2){ // derive(f,x)
	gen othervar=(*s.feuille._VECTptr)[1];
	if (othervar.type!=_IDNT) setsizeerr("derive.cc/derive_SYMB");
	if (*othervar._IDNTptr==i){ // _FUNCnd derivative
	  vecteur res(*s.feuille._VECTptr);
	  symbolic sprime(s);
	  res.push_back(2);
	  return symbolic(at_derive,gen(res,_SEQ__VECT));
	}
	else {
	  vecteur var;
	  var.push_back(othervar);
	  var.push_back(i);
	  vecteur nderiv;
	  nderiv.push_back(1);
	  nderiv.push_back(1);
	  return symbolic(at_derive,gen(makevecteur((*s.feuille._VECTptr)[0],var,nderiv),_SEQ__VECT));
	}
      }
      else { // derive(f,x,n)
	if (s.feuille._VECTptr->size()!=3)  setsizeerr("derive.cc/derive_SYMB");
	gen othervar=(*s.feuille._VECTptr)[1];
	if (othervar.type==_IDNT){
	  if (*othervar._IDNTptr==i){ // n+1 derivative
	    vecteur vprime=(*s.feuille._VECTptr);
	    vprime[2] += 1;
	    return symbolic(s.sommet,gen(vprime,_SEQ__VECT));
	  }
	  else {
	    vecteur var;
	    var.push_back(othervar);
	    var.push_back(i);
	    vecteur nderiv;
	    nderiv.push_back((*s.feuille._VECTptr)[2]);
	    nderiv.push_back(1);
	    return symbolic(at_derive,gen(makevecteur((*s.feuille._VECTptr)[0],var,nderiv),_SEQ__VECT));
	  }
	} // end if othervar.type==_IDNT
	else { // othervar.type must be _VECT
	  if (othervar.type!=_VECT)  setsizeerr("derive.cc/derive_SYMB");
	  gen nder((*s.feuille._VECTptr)[2]);
	  if (nder.type!=_VECT ||
	      nder._VECTptr->size()!=othervar._VECTptr->size())  setsizeerr("derive.cc/derive_SYMB");
	  vecteur nderiv(*nder._VECTptr);
	  int pos=equalposcomp(*othervar._VECTptr,i);
	  if (pos){
	    nderiv[pos-1]=nderiv[pos-1]+1;
	  }
	  else {
	    othervar._VECTptr->push_back(i);
	    nderiv.push_back(1);
	  }
	  return symbolic(at_derive,gen(makevecteur((*s.feuille._VECTptr)[0],othervar,nderiv),_SEQ__VECT));
	}
      }
    }
    // no info about derivative
    return symbolic(at_derive,gen(makevecteur(s,i),_SEQ__VECT));
    i.dbgprint();
    s.dbgprint();
  }

  gen derive_VECT(const vecteur & v,const identificateur & i,GIAC_CONTEXT){
    vecteur w;
    w.reserve(v.size());
    vecteur::const_iterator it=v.begin(),itend=v.end();
    for (;it!=itend;++it)
      w.push_back(derive(*it,i,contextptr));
    return w;
  }

  gen derive(const gen & e,const identificateur & i,GIAC_CONTEXT){
    switch (e.type){
    case _INT_: case _DOUBLE_: case _ZINT: case _CPLX: case _MOD: case _REAL: case _USER:
      return 0;
    case _IDNT:
      if (*e._IDNTptr==i)
	return 1;
      else
	return 0;
    case _SYMB:
      return derive_SYMB(*e._SYMBptr,i,contextptr);
    case _VECT:
      return derive_VECT(*e._VECTptr,i,contextptr);
    case _FRAC:
      return fraction(derive(e._FRACptr->num,i,contextptr)*e._FRACptr->den-(e._FRACptr->num)*derive(e._FRACptr->den,i,contextptr),e._FRACptr->den);
    default:
      settypeerr();
    }
    return 0;
  }

  gen _VECTderive(const gen & e,const vecteur & v,GIAC_CONTEXT){
    vecteur w;
    w.reserve(v.size());
    vecteur::const_iterator it=v.begin(),itend=v.end();
    for (;it!=itend;++it)
      w.push_back(derive(e,*it,contextptr));
    return w;    
  }

  gen derivesymb(const gen& e,const gen & var,GIAC_CONTEXT){
    identificateur x(" x");
    gen xx(x);
    gen f=subst(e,var,xx,false,contextptr);
    f=derive(f,x,contextptr);
    f=subst(f,xx,var,false,contextptr);
    return f;
  }
  gen derive(const gen & e,const gen & vars,GIAC_CONTEXT){
          //  cout << e << " " << vars << endl;
    switch (vars.type){
    case _IDNT:
      return derive(e,*vars._IDNTptr,contextptr);
    case _VECT:
      return _VECTderive(e,*vars._VECTptr,contextptr);
    case _SYMB:
      return derivesymb(e,vars,contextptr);
    default:
      setsizeerr();
    }
    return 0;
  }

  gen derive(const gen & e,const gen & vars,const gen & nderiv,GIAC_CONTEXT){
    if (nderiv.type==_INT_){
      int n=nderiv.val;
      gen ecopie(e),eprime(e);
      int j=1;
      for (;j<=n;++j){
	eprime=derive(ecopie,vars,contextptr);
	if ( (eprime.type==_SYMB) && (eprime._SYMBptr->sommet==at_derive))
	  break;
	ecopie=eprime;
      }
      if (j==n+1)
	return eprime;
      return symbolic(at_derive,gen(makevecteur(ecopie,vars,n+1-j),_SEQ__VECT));
    }
    // multi-index derivation
    if (nderiv.type!=_VECT ||
	vars.type!=_VECT)  setsizeerr("derive.cc/derive");
    int s=nderiv._VECTptr->size();
    if (s!=signed(vars._VECTptr->size()))  setsizeerr("derive.cc/derive");
    int j=0;
    gen ecopie(e);
    for (;j<s;++j){
      ecopie=derive(ecopie,(*vars._VECTptr)[j],(*nderiv._VECTptr)[j],contextptr);
    }
    return ecopie;
  }  

  symbolic symb_derive(const gen & a,const gen & b){
    return symbolic(at_derive,gen(makevecteur(a,b),_SEQ__VECT));
  }

  gen symb_derive(const gen & a,const gen & b,const gen &c){
    if (is_zero(c))
      return a;
    if (is_one(c))
      return symb_derive(a,b);
    return symbolic(at_derive,gen(makevecteur(a,b,c),_SEQ__VECT));
  }

  // "unary" version
  gen _derive(const gen & args,GIAC_CONTEXT){
    vecteur v;
    if (args.type==_VECT && args.subtype==_POLY1__VECT)
      return gen(derivative(*args._VECTptr),_POLY1__VECT);
    if (args.type==_VECT)
      v=plotpreprocess(gen(*args._VECTptr,_SEQ__VECT),contextptr);
    else
      v=plotpreprocess(args,contextptr);
    int s=v.size();
    if (s==2){
      if (v[1].type==_VECT && v[1].subtype==_SEQ__VECT){
	vecteur & w=*v[1]._VECTptr;
	int ss=w.size();
	gen res=v[0];
	for (int i=0;i<ss;++i)
	  res=ratnormal(derive(res,w[i],contextptr));
	return res;
      }
      if (args.type!=_VECT && v[0].type==_VECT && v[0].subtype==_POLY1__VECT)
	return gen(derivative(*v[0]._VECTptr),_POLY1__VECT);
      return derive(v[0],v[1],contextptr);
    }
    if ((s==3) && (v[2].type==_INT_ || v[2].type==_VECT) )
      return derive( v[0],v[1],v[2],contextptr);    
    if (s<3)
      setsizeerr();
    const_iterateur it=v.begin()+1,itend=v.end();
    gen res(v[0]);
    for (;it!=itend;++it)
      res=ratnormal(derive(res,*it,contextptr));
    return res;
  }
  const string _derive_s("diff");
  string texprintasderive(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
      if (feuille.type!=_VECT)
          return gen2tex(feuille,contextptr)+"'";
    return "\\frac{\\partial \\left("+gen2tex(feuille._VECTptr->front(),contextptr)+"\\right)}{\\partial "+gen2tex(feuille._VECTptr->back(),contextptr)+"}";
  }
  unary_function_eval __derive(&_derive,_derive_s,0,texprintasderive);
  unary_function_ptr at_derive (&__derive,_QUOTE_ARGUMENTS,true);

  gen _grad(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT || args._VECTptr->size()!=2)
      setsizeerr();
    return _derive(args,contextptr);
  }
  const string _grad_s("grad");
  unary_function_eval __grad(&_grad,_grad_s);
  unary_function_ptr at_grad (&__grad,_QUOTE_ARGUMENTS,true);

  // FIXME: This should not use any temporary identifier
  // Should define the identity operator and write again all rules here
  // NB: It requires all D operators for functions to be functions!
  gen _function_diff(const gen & g,GIAC_CONTEXT){
    if (g.is_symb_of_sommet(at_function_diff)){
      gen & f = g._SYMBptr->feuille;
      return symbolic(at_of,makevecteur(gen(symbolic(at_composepow,makevecteur(at_function_diff,2))),f));
    }
    if (g.is_symb_of_sommet(at_of)){
      gen & f = g._SYMBptr->feuille;
      if (f.type==_VECT && f._VECTptr->size()==2){
	gen & f1=f._VECTptr->front();
	gen & f2=f._VECTptr->back();
	if (f1.is_symb_of_sommet(at_composepow)){
	  gen & f1f=f1._SYMBptr->feuille;
	  if (f1f.type==_VECT && f1f._VECTptr->size()==2 && f1f._VECTptr->front()==at_function_diff){
	    return symbolic(at_of,makevecteur(gen(symbolic(at_composepow,makevecteur(at_function_diff,f1f._VECTptr->back()+1))),f2));
	  }
	}
      }
    }
    identificateur _tmpi(" _x");
    gen _tmp(_tmpi);    
    gen dg(derive(g(_tmp,contextptr),_tmp,contextptr));
    if (lop(dg,at_derive).empty()){
      identificateur tmpi(" x");
      gen tmp(tmpi);
      gen res=symb_program(tmp,zero,quotesubst(dg,_tmp,tmp,contextptr),contextptr);
      return res;
    }
    return symbolic(at_function_diff,g);
  }
  const string _function_diff_s("function_diff");
  unary_function_eval __function_diff(&_function_diff,_function_diff_s);
  unary_function_ptr at_function_diff (&__function_diff,0,true);

  const string _fonction_derivee_s("fonction_derivee");
  unary_function_eval __fonction_derivee(&_function_diff,_fonction_derivee_s);
  unary_function_ptr at_fonction_derivee (&__fonction_derivee,0,true);

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC
