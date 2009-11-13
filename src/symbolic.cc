// -*- mode:C++ ; compile-command: "g++ -I.. -g -c symbolic.cc" -*-
#include "first.h"
/*
 *  Copyright (C) 2000, 2007 B. Parisse, Institut Fourier, 38402 St Martin d'Heres
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
#include "symbolic.h"
#include "identificateur.h"
#include "usual.h"
#include "prog.h"
#include "rpn.h"
#include "plot.h"

// NB: the operator in the symbolic (sommet) can not be replace by a pointer
// to a unary_function_ptr in the current code. Indeed, symbolic are created
// in the parser from temporary gen objects of type _FUNC which contains
// copies of unary_function_ptr, these copies are deleted once the parser
// gen object is deleted -> segfault.

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  // unary_function_ptr quoted_op[]={at_of,at_for,at_bloc,at_local,at_program,at_rpn_prog,at_ifte,at_try_catch,at_print,at_signal,at_as_function_of,at_lieu,at_legende,at_debug,at_sst,at_sst_in,at_cont,at_kill,at_halt,at_watch,at_rmwatch,at_breakpoint,at_maple2mupad,at_mupad2maple,at_maple2xcas,at_mupad2xcas,at_purge,0};
  unary_function_ptr archive_function_tab[]={at_plus,at_neg,at_binary_minus,at_prod,at_division,at_inv,at_pow,at_exp,at_ln,at_abs,at_arg,at_pnt,at_point,at_segment,at_sto,at_sin,at_cos,at_tan,at_asin,at_acos,at_atan,at_sinh,at_cosh,at_tanh,at_asinh,at_acos,at_atanh,at_interval,at_union,at_minus,at_intersect,at_not,at_and,at_ou,at_inferieur_strict,at_inferieur_egal,at_superieur_strict,at_superieur_egal,at_different,at_equal,at_rpn_prog,at_local,at_return,at_Dialog,at_double_deux_points,at_pointprod,at_pointdivision,at_pointpow,at_hash,at_pourcent,at_tilocal,at_break,at_continue,at_ampersand_times,at_maple_lib,at_unit,at_plot_style,at_xor,at_check_type,at_quote_pow,at_case,at_dollar,at_IFTE,at_RPN_CASE,at_RPN_LOCAL,at_RPN_FOR,at_RPN_WHILE,at_NOP,at_unit,0};

  int equalposcomp(unary_function_ptr tab[],const unary_function_ptr & f){
    for (int i=1;*tab!=0;++tab,++i){
      if (*tab==f)
	return i;
    }
    return 0;
  }

  symbolic::symbolic(const symbolic & mys,const gen & e): sommet(mys.sommet){
    vecteur tmp;
    if (mys.feuille.type==_VECT){
      tmp = *mys.feuille._VECTptr;
      tmp.push_back(e);
    }
    else {
      tmp.push_back(mys.feuille);
      tmp.push_back(e);
    }
    feuille = tmp;
  };
  
  symbolic::symbolic(const gen & a,const unary_function_ptr & o,const gen & b):sommet(o) {
    if (b.type==_VECT)
      feuille= mergevecteur(vecteur(1,a),*b._VECTptr);
    else
      feuille=makevecteur(a,b);
  };
  
  int symbolic::size() const {
    if (feuille.type==_SYMB)
      return 1+feuille._SYMBptr->size();
    if (feuille.type!=_VECT)
      return 2;
    int s=1;
    iterateur it=feuille._VECTptr->begin(),itend=feuille._VECTptr->end();
    for (;it!=itend;++it){
      if (it->type==_SYMB)
	s += it->_SYMBptr->size();
      else
	++s;
    }
    return s;
  }

  bool print_rewrite_prod_inv=false;
  string symbolic::print(GIAC_CONTEXT) const{
    // first look if operator has it's own printing method
    if ( (sommet.ptr->printsommet) && (sommet!=at_plus) && (sommet!=at_prod) )
      return sommet.ptr->printsommet(feuille,sommet.ptr->s,contextptr);
    // default
    if ( (feuille.type==_VECT) && (feuille._VECTptr->empty()) )
        return sommet.ptr->print(contextptr)+"(NULL)";
    if ( (feuille.type!=_VECT) || ( (sommet!=at_prod) && (feuille._VECTptr->front().type==_VECT) ) ){
      if (sommet==at_neg){
	string tmp=feuille.print(contextptr);
	if (feuille.type!=_CPLX){
	  if (feuille.type!=_SYMB || (feuille._SYMBptr->sommet==at_inv || feuille._SYMBptr->sommet==at_prod) ){
	    if (!tmp.empty() && tmp[0]=='-')
	      return tmp.substr(1,tmp.size()-1);
	    return "-"+ tmp;
	  }
	}
	return string("-(") + tmp +string(")");
      }
      if (sommet==at_inv){
	gen f = feuille;
	bool isneg=false;
	if (f.is_symb_of_sommet(at_neg) && f._SYMBptr->feuille.type!=_VECT){
	  f = f._SYMBptr->feuille;
	  isneg=true;
	}
	if (f.type<_CPLX && is_positive(-f,contextptr)){
	  f=-f;
	  isneg=!isneg;
	}
	if ((f.type!=_SYMB) && (f.type!=_CPLX) && (f.type!=_MOD))
	  return (isneg?"-1/":"1/")+f.print(contextptr);
	else
	  return (isneg?"-1/(":"1/(")+f.print(contextptr)+")";
      }
      return sommet.ptr->print(contextptr)+string("(") + feuille.print(contextptr) +string(")");
    }
    string s,tmp;
    string opstring( sommet.ptr->print(contextptr));
    int l=feuille._VECTptr->size(),tmps;
    if ( sommet==at_plus ){
      for (int i=0;i<l;++i){
	gen e((*(feuille._VECTptr))[i]);
	if ((e.type==_SYMB) && (e._SYMBptr->sommet==at_neg)){
	  tmp = e._SYMBptr->feuille.print(contextptr);
	  if (!tmp.empty() && tmp[0]=='-')
	    tmp=tmp.substr(1,tmp.size()-1);
	  else
	    tmp = "-" + tmp;
	}
	else
	  tmp = e.print(contextptr);
	if (i){
	  if (!tmp.empty() && tmp[0]=='-'){
	    if ( (tmps=tmp.size())>2 && tmp[1]=='-')
	      s+='+'+tmp.substr(2,tmps-2);
	    s+=tmp;
	  }
	  else
	    s+= "+"+tmp;
	}
	else
	  s=tmp;
      } // end_for
      return s;
    } // end at_plus
    if (sommet==at_prod) {
      gen n0,d0;
      if (print_rewrite_prod_inv &&
	  rewrite_prod_inv(feuille,n0,d0)
	  ){
	if (n0.type<_CPLX || n0.type==_IDNT)
	  s=n0.print(contextptr)+"/";
	else
	  s="("+n0.print(contextptr)+")/";
	if (d0.type<_CPLX || d0.type==_IDNT)
	  s+=d0.print(contextptr);
	else
	  s+="("+d0.print(contextptr)+")";
	return s;
      }
      for (int i=0;i<l;++i){
	gen e((*(feuille._VECTptr))[i]);
	if (e.type!=_SYMB){
	  if (i)
	    s += opstring;
	  if ( (e.type==_CPLX) || (e.type==_MOD) )
	    s += "("+e.print(contextptr)+")";
	  else
	    s +=e.print(contextptr);
	}
	else {
	  if (e._SYMBptr->sommet==at_inv){
	    gen f(e._SYMBptr->feuille);
	    if (i){
	      if ( (f.type==_CPLX) || (f.type==_MOD) ||
		   ((f.type==_SYMB) && 
		    ( (f._SYMBptr->sommet==at_plus) || (f._SYMBptr->sommet==at_prod) || (f._SYMBptr->sommet==at_neg) || f._SYMBptr->sommet==at_inv ))
		   )
		s += string("/(")+f.print(contextptr) + string(")");
	      else
		s += string("/")+f.print(contextptr);
	    }
	    else
	      s += e.print(contextptr);
	  } // end if e._SYMBptr->sommet==at_inv
	  else {
	    if (i)
	      s += opstring;
	    if ( (e._SYMBptr->sommet==at_plus) || (e._SYMBptr->sommet==at_neg) )
	      s += string("(")+e.print(contextptr)+string(")");
	    else
	      s += e.print(contextptr);
	  }
	}
      } // end_for
      return s;
    } // end if sommet_is_prod
    s = opstring + '(';
    if (feuille.subtype)
      s += begin_VECT_string(feuille.subtype,false,contextptr);
    for (int i=0;;++i){
      s += (*(feuille._VECTptr))[i].print(contextptr);
      if (i==l-1){
	break;
      }
      s += ',';
    }
    if (feuille.subtype)
      s += end_VECT_string(feuille.subtype,false,contextptr);
    return s+')';
  }

  gen symbolic::eval(int level,const context * contextptr) const {
    if (level==0)
      return *this;
    std::vector<const std::string *> & last =last_evaled_function_name(contextptr);
    last.push_back(&sommet.ptr->s);
    gen ans;
    if (sommet==at_sto && feuille.type==_VECT){ // autoname function
      gen & feuilleback=feuille._VECTptr->back();
      if ( feuilleback.type==_SYMB && (feuilleback._SYMBptr->sommet==at_unquote || feuilleback._SYMBptr->sommet==at_hash ) ){
	ans=_sto(feuille.eval(level,contextptr),contextptr);
	if (!last.empty())
	  last.pop_back();
	return ans;
      }
      bool b=show_point(contextptr);
      if (b)
	show_point(false,contextptr);
      gen e=feuille._VECTptr->front().eval(level,contextptr);
      if (b)
	show_point(b,contextptr);
      if (e.type==_SYMB && e._SYMBptr->sommet==at_pnt && e._SYMBptr->feuille.type==_VECT && e._SYMBptr->feuille._VECTptr->size()==2 && (contextptr?!contextptr->previous:!protection_level) ){
	e=new symbolic(at_pnt,gen(makevecteur(e._SYMBptr->feuille._VECTptr->front(),e._SYMBptr->feuille._VECTptr->back(),string2gen(feuilleback.print(contextptr),false)),_PNT__VECT));
	e.subtype=gnuplot_show_pnt(*e._SYMBptr,contextptr);
      }
      if ( e.type==_VECT && !e._VECTptr->empty() && e._VECTptr->back().type==_SYMB && e._VECTptr->back()._SYMBptr->sommet==at_pnt && (contextptr?!contextptr->previous:!protection_level)){
	vecteur v=*e._VECTptr;
	iterateur it=v.begin(),itend=v.end();
	gen legende;
	for (int pos=0;it!=itend;++pos,++it){
	  if ( (it->type==_SYMB) && (it->_SYMBptr->sommet==at_pnt) && (it->_SYMBptr->feuille._VECTptr->size()==2)){
	    if (feuilleback.type==_VECT && pos<feuilleback._VECTptr->size())
	      legende=(*feuilleback._VECTptr)[pos];
	    else
	      legende=feuilleback;
	    *it=new symbolic(at_pnt,gen(makevecteur(it->_SYMBptr->feuille._VECTptr->front(),it->_SYMBptr->feuille._VECTptr->back(),string2gen(legende.print(contextptr),false)),_PNT__VECT));
	    it->subtype=gnuplot_show_pnt(*it->_SYMBptr,contextptr);
	  }
	}
	e=gen(v,e.subtype);
      }
      ans=sto(e,feuilleback,contextptr);
      if (!last.empty())
	last.pop_back();
      return ans;
    } // end sommet==at_sto
    if (sommet.quoted){ 
      if (feuille.type==_SYMB){ 
	unary_function_ptr & u=feuille._SYMBptr->sommet;
	if (u==at_unquote){
	  ans=sommet(feuille.eval(level,contextptr),contextptr);
	  if (!last.empty())
	    last.pop_back();
	  return ans;
	}
	if (u==at_hash){
	  ans=sommet(gen(*feuille._SYMBptr->feuille._STRNGptr,contextptr),contextptr);
	  if (!last.empty())
	    last.pop_back();
	  return ans;
	}
      }
      int & elevel=eval_level(contextptr);
      int save_level=elevel;
      elevel=level;
      try {
	ans=sommet(feuille,contextptr);
      }
      catch (std::runtime_error & err){
	elevel=save_level;
	throw(err);
      }
      elevel=save_level;
      if (!last.empty())
	last.pop_back();
      return ans;
    } // if (sommet.quoted)
    else {
      if ((sommet==at_neg) && (feuille.type==_IDNT) && (*feuille._IDNTptr==_IDNT_infinity())){
	if (!last.empty())
	  last.pop_back();
	return minus_inf;
      }
      if (sommet==at_quote){
	if (!last.empty())
	  last.pop_back();
	return quote(feuille,contextptr);
      }
      if (feuille.in_eval(level,ans,contextptr))
	ans=(*sommet.ptr)(ans,contextptr);
      else
	ans=(*sommet.ptr)(feuille,contextptr);
      if (!last.empty())
	last.pop_back();
      return ans;
    }
  }

  bool rewrite_prod_inv(const gen & arg,gen & n,gen & d){
    n=1; d=1;
    if (arg.type==_VECT && !arg._VECTptr->empty() && arg._VECTptr->back().is_symb_of_sommet(at_inv)) {
      vecteur & uv=*arg._VECTptr;
      int tmps=uv.size(),invbegin;
      vecteur den(1,uv.back()._SYMBptr->feuille);
      // group all inv from the end to the beginning for the denominator
      for (invbegin=tmps-2;invbegin>=0;--invbegin){
	if (!uv[invbegin].is_symb_of_sommet(at_inv))
	  break;
	den.push_back(uv[invbegin]._SYMBptr->feuille);
      }
      vecteur num;
      for (int i=0;i<=invbegin;++i){
	if (uv[i].is_symb_of_sommet(at_inv) && uv[i]._SYMBptr->feuille.type<_POLY)
	  d=d*uv[i]._SYMBptr->feuille;
	else
	  num.push_back(uv[i]);
      }
      if (!is_one(d))
	den.insert(den.begin(),d);
      if (den.size()>1)
	d=new symbolic(at_prod,den);
      else
	d=den.front();
      if (!num.empty()){
	if (num.size()==1)
	  n=num.front();
	else 
	  n=new symbolic(at_prod,num);
      }
      return true;
    }
    // Group scalar denominators (warning, do not use for matrices!)
    vecteur num,den;
    prod2frac(new symbolic(at_prod,arg),num,den);
    if (!den.empty()){
      if (num.empty())
	n=plus_one;
      else {
	if (num.size()==1)
	  n=num.front();
	else
	  n=new symbolic(at_prod,num);
      }
      /* code that does not work with matrices
	 if (den.size()==1)
	 d=den.front();
	 else
	 d=symbolic(at_prod,den); 
      */
      if (den.size()==1 && den.front().type<_IDNT){
	d=den.front();
	return true;
      }
    }
    return false;
  }
    
  gen symbolic::evalf(int level,const context * contextptr) const {
    if (level==0)
      return *this;
    std::vector<const std::string *> & last =last_evaled_function_name(contextptr);
    last.push_back(&sommet.ptr->s);
    if (sommet==at_sto){ // autoname function
      gen e=feuille._VECTptr->front().evalf(level,contextptr);
      if ((e.type==_SYMB) && (e._SYMBptr->sommet==at_pnt) && (e._SYMBptr->feuille.type==_VECT) && (e._SYMBptr->feuille._VECTptr->size()==2))
	e=new symbolic(at_pnt,gen(makevecteur(e._SYMBptr->feuille._VECTptr->front(),e._SYMBptr->feuille._VECTptr->back(),string2gen(feuille._VECTptr->back().print(contextptr),false)),_PNT__VECT));
      if ( (e.type==_VECT) && (e._VECTptr->size()) && (e._VECTptr->back().type==_SYMB) && (e._VECTptr->back()._SYMBptr->sommet==at_pnt)){
	vecteur v=*e._VECTptr;
	iterateur it=v.begin(),itend=v.end();
	for (;it!=itend;++it){
	  if ( (it->type==_SYMB) && (it->_SYMBptr->sommet==at_pnt) && (it->_SYMBptr->feuille._VECTptr->size()==2))
	    *it=new symbolic(at_pnt,gen(makevecteur(it->_SYMBptr->feuille._VECTptr->front(),it->_SYMBptr->feuille._VECTptr->back(),string2gen(feuille._VECTptr->back().print(contextptr),false)),_PNT__VECT));
	}
	e=v;
      }
      if (!last.empty())
	last.pop_back();
      return sto(e,feuille._VECTptr->back(),contextptr);
    }
    gen ans;
    if (sommet.quoted && !equalposcomp(plot_sommets,sommet) ){ 
      ans=sommet(feuille,contextptr);
      if (!last.empty())
	last.pop_back();
      return ans;
    }
    else {
      if ((sommet==at_neg) && (feuille.type==_IDNT) && (*feuille._IDNTptr==_IDNT_infinity())){
	if (!last.empty())
	  last.pop_back();
	return minus_inf;
      }
      if (sommet==at_quote){
	if (!last.empty())
	  last.pop_back();
	return quote(feuille,contextptr);
      }
      if (equalposcomp(plot_sommets,sommet)){
	// bool save_is_inevalf=is_inevalf;
	// is_inevalf=true;
	ans=new symbolic(sommet,feuille.evalf(1,contextptr));
	// is_inevalf=save_is_inevalf;
	if (!last.empty())
	  last.pop_back();
	return ans;
      }
      ans=(*sommet.ptr)(feuille.evalf(level,contextptr),contextptr);
      if (!last.empty())
	last.pop_back();
      return ans;
    }
  }


  unsigned taille(const gen & g,unsigned max){
    if (g.type<=_IDNT)
      return 1;
    if (g.type==_FRAC)
      return 1+taille(g._FRACptr->num,max)+taille(g._FRACptr->den,max);
    if (g.type==_SYMB){
      if (g.is_symb_of_sommet(at_curve))
	return 10;
      return 1+taille(g._SYMBptr->feuille,max);
    }
    if (g.type==_VECT){
      unsigned res=0;
      const_iterateur it=g._VECTptr->begin(),itend=g._VECTptr->end();
      for (;it!=itend;++it){
	res += taille(*it,max);
	if (max && res>max)
	  return res;
      }
      return res;
    }
    return 2;
  }

    
  ostream & operator << (ostream & os,const symbolic & s) { return os << s.print(context0); }

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC
