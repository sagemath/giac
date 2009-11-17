// -*- mode:C++ ; compile-command: "g++ -I.. -g -c identificateur.cc" -*-
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
#include <cmath>
#include <fstream>
//#include <unistd.h> // For reading arguments from file
#include "identificateur.h"
#include "gen.h"
#include "sym2poly.h"
#include "rpn.h"
#include "prog.h"
#ifdef WIN32
#include "usual.h"
#endif

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  // bool variables_are_files=true; // FIXME -> false and change rpn.cc at_VARS
  int protection_level=0; // for local variables in null context

  string string_euler_gamma("euler_gamma");
  identificateur _IDNT_euler_gamma(string_euler_gamma,(double) .577215664901533);
  gen cst_euler_gamma(_IDNT_euler_gamma);
  string string_pi("pi");
  identificateur _IDNT_pi(string_pi,(double) M_PI);
  gen cst_pi(_IDNT_pi);
  string string_infinity("infinity");
  identificateur & _IDNT_infinity(){
    static identificateur * ans=new identificateur("infinity");
    return * ans;
  }
  gen unsigned_inf(_IDNT_infinity());
  
  /*
#ifndef GNUWINCE
#ifdef WIN32
  gen cst_two_pi(symbolic(at_prod,makevecteur(2,cst_pi)));
  gen cst_pi_over_2(_FRAC2_SYMB(cst_pi,2));
  gen plus_inf(symbolic(at_plus,_IDNT_infinity()));
  gen minus_inf(symbolic(at_neg,_IDNT_infinity()));
#endif
#endif
  */
  string string_undef("undef");
  identificateur & _IDNT_undef(){
    static identificateur * ans=new identificateur("undef");
    return * ans;
  }
  gen undef(_IDNT_undef());
  identificateur a__IDNT("a");
  gen a__IDNT_e(a__IDNT);
  identificateur b__IDNT("b");
  gen b__IDNT_e(b__IDNT);
  identificateur c__IDNT("c");
  gen c__IDNT_e(c__IDNT);
  identificateur d__IDNT("d");
  gen d__IDNT_e(d__IDNT);
  identificateur e__IDNT("e");
  gen e__IDNT_e(e__IDNT);
  identificateur f__IDNT("f");
  gen f__IDNT_e(f__IDNT);
  identificateur g__IDNT("g");
  gen g__IDNT_e(g__IDNT);
  identificateur h__IDNT("h");
  gen h__IDNT_e(h__IDNT);
  identificateur i__IDNT("i");
  gen i__IDNT_e(i__IDNT);
  identificateur I__IDNT("I");
  gen I__IDNT_e(I__IDNT);
  identificateur j__IDNT("j");
  gen j__IDNT_e(j__IDNT);
  identificateur k__IDNT("k");
  gen k__IDNT_e(k__IDNT);
  identificateur l__IDNT("l");
  gen l__IDNT_e(l__IDNT);
  identificateur m__IDNT("m");
  gen m__IDNT_e(m__IDNT);
  identificateur n__IDNT("n");
  gen n__IDNT_e(n__IDNT);
  identificateur o__IDNT("o");
  gen o__IDNT_e(o__IDNT);
  identificateur p__IDNT("p");
  gen p__IDNT_e(p__IDNT);
  identificateur q__IDNT("q");
  gen q__IDNT_e(q__IDNT);
  identificateur r__IDNT("r");
  gen r__IDNT_e(r__IDNT);
  identificateur s__IDNT("s");
  gen s__IDNT_e(s__IDNT);
  identificateur t__IDNT("t");
  gen t__IDNT_e(t__IDNT);
  identificateur u__IDNT("u");
  gen u__IDNT_e(u__IDNT);
  identificateur v__IDNT("v");
  gen v__IDNT_e(v__IDNT);
  identificateur w__IDNT("w");
  gen w__IDNT_e(w__IDNT);
  identificateur x__IDNT("x");
  gen x__IDNT_e(x__IDNT);
  gen vx_var(x__IDNT_e);
  identificateur y__IDNT("y");
  gen y__IDNT_e(y__IDNT);
  identificateur z__IDNT("z");
  gen z__IDNT_e(z__IDNT);
  identificateur CST__IDNT("CST");
  gen CST__IDNT_e(CST__IDNT);
  identificateur PICT__IDNT("PICT");
  gen PICT__IDNT_e(PICT__IDNT);
  gen tab_one_char__IDNT[]={a__IDNT_e,b__IDNT_e,c__IDNT_e,d__IDNT_e,e__IDNT_e,f__IDNT_e,g__IDNT_e,h__IDNT_e,i__IDNT_e,I__IDNT_e,j__IDNT_e,k__IDNT_e,l__IDNT_e,m__IDNT_e,n__IDNT_e,o__IDNT_e,p__IDNT_e,q__IDNT_e,r__IDNT_e,s__IDNT_e,t__IDNT_e,u__IDNT_e,v__IDNT_e,w__IDNT_e,x__IDNT_e,y__IDNT_e,z__IDNT_e,CST__IDNT_e,PICT__IDNT_e,zero};
  vecteur list_one_letter__IDNT(tab2vecteur(tab_one_char__IDNT));
  string string_break("break");
  identificateur _IDNT_break(string_break);
  string string_continue("continue");
  identificateur _IDNT_continue(string_continue);

  identificateur::identificateur(){
    ref_count = new int(1);
    value = NULL;
    quoted=new bool(false);
    localvalue = new vecteur;
    name = new string(" "+print_INT_(rand()));
  }

  identificateur::identificateur(const string & s){
    ref_count = new int(1);
    value = NULL;
    quoted=new bool(false);
    localvalue = new vecteur;
    name = new string(s);
  }

  identificateur::identificateur(const char * s){
    ref_count = new int(1);
    value = NULL;
    quoted=new bool(false);
    localvalue = new vecteur;
    name = new string(s);
  }

  identificateur::identificateur(const string & s,const gen & e){
    ref_count = new int(1);
    value = new gen(e);
    quoted=new bool(false);
    localvalue = new vecteur;
    name = new string(s);
  }

  identificateur::identificateur(const identificateur & s){
    ref_count=s.ref_count;
    if (ref_count){
      ++(*ref_count);
      value=s.value;
      quoted=s.quoted;
      localvalue=s.localvalue;
      name=s.name;
    }
  }

  identificateur::~identificateur(){
    if (ref_count){
      --(*ref_count);
      if (!(*ref_count)){
	delete ref_count;
	delete name;
        delete quoted;
	if (value)
	  delete value;
        delete localvalue;
      }
    }
  }

  identificateur & identificateur::operator =(const identificateur & s){
    if (ref_count){
      --(*ref_count);
      if (!(*ref_count)){
	delete ref_count;
	delete name;
	if (value)
	  delete value;
        delete localvalue;
      }
    }
    ref_count=s.ref_count;
    if (ref_count){
      ++(*ref_count);
      value=s.value;
      quoted=s.quoted;
      localvalue=s.localvalue;
      name=s.name;
    }
    return *this;
  }

  gen globalize(const gen & g){
    gen tmp(g);
    switch (tmp.type){
    case _IDNT:
      tmp.subtype=_GLOBAL__EVAL;
      break;
    case _VECT:
      tmp=apply(tmp,globalize);
      break;
    case _SYMB:
      if (tmp._SYMBptr->sommet!=at_program)
	tmp=symbolic(tmp._SYMBptr->sommet,globalize(tmp._SYMBptr->feuille));
      break;
    }
    return tmp;
  }

  // make g identificateurs evaluated as global
  gen global_eval(const gen & g,int level){
    if (g.type<_IDNT)
      return g;
    bool save_local_eval=local_eval(context0);
    local_eval(false,context0);
    gen tmp;
    try {
      tmp=g.eval(level,context0);
    }
    catch (std::runtime_error & e){
      cerr << e.what() << endl;
      // eval_level(level,contextptr);
    }
    local_eval(save_local_eval,context0);
    return globalize(tmp);
  }

  bool check_not_assume(const gen & not_evaled,gen & evaled, bool evalf_after,const context * contextptr);

  // make g identificateurs evaluated as global
  gen global_evalf(const gen & g,int level){
    if (g.type<_IDNT)
      return g;
    bool save_local_eval=local_eval(context0);
    local_eval(false,context0);
    gen tmp;
    try {
      tmp=g.eval(level,context0);
      if (tmp.type==_IDNT){
	gen evaled(tmp._IDNTptr->eval(level,tmp,context0));
	if (check_not_assume(tmp,evaled,true,context0))
	  tmp=evaled;
      }
    }
    catch (std::runtime_error & e){
      cerr << e.what() << endl;
      // eval_level(level,contextptr);
    }
    local_eval(save_local_eval,context0);
    return globalize(tmp);
  }

  gen identificateur::eval(int level,const gen & orig,const context * contextptr) {
    if (!ref_count)
      settypeerr();
    // cerr << "idnt::eval " << *this << " " << level << endl;
    if (level<=0){
      if (level==0) 
	return orig;
      if (contextptr){
	sym_tab::const_iterator it=contextptr->tabptr->find(*name);
	return (it==contextptr->tabptr->end())?orig:it->second;
      }
      else {
	if (localvalue->empty())
	  return orig;
	iterateur jtend=localvalue->end();
	return (protection_level>(jtend-2)->val)?localvalue->back():orig;
      }
    }
    --level;
    gen evaled;
    if (in_eval(level,orig,evaled,contextptr))
      return evaled;
    else
      return *this;
    /*
    int save_level=eval_level(contextptr);
    eval_level(level,contextptr);
    gen res=in_eval(level,contextptr);
    eval_level(save_level,contextptr);
    return res;
    */
  }

  // if globalize is true, use global value in eval
  gen do_local_eval(const identificateur & i,int level,bool globalize) {
    gen res;
    iterateur jtend=i.localvalue->end();
    if (protection_level>(jtend-2)->val)
      res=i.localvalue->back();
    else {
      for (iterateur jt=i.localvalue->begin();;){
	if (jt==jtend)
	  break;
	--jtend;
	--jtend;
	if (protection_level>jtend->val){
	  ++jtend;
	  ++jtend;
	  break;
	}
      }
      i.localvalue->erase(jtend,i.localvalue->end());
      if (!i.localvalue->empty())
	res=i.localvalue->back();
    }
    return globalize?global_eval(res,level):res; 
  }

  bool identificateur::in_eval(int level,const gen & orig,gen & evaled,const context * contextptr) {
    if (contextptr){
      const context * cur=contextptr;
      for (;cur->previous;cur=cur->previous){
	sym_tab::const_iterator it=cur->tabptr->find(*name);
	if (it!=cur->tabptr->end()){
	  if (!it->second.in_eval(level,evaled,contextptr->globalcontextptr))
	    evaled=it->second;
	  return true;
	}
      }
      // now at global level
      // check for quoted
      if (cur->quoted_global_vars && !cur->quoted_global_vars->empty() && equalposcomp(*cur->quoted_global_vars,orig)) 
	return false;
      sym_tab::const_iterator it=cur->tabptr->find(*name);
      if (it==cur->tabptr->end())
	return false;
      else {
	if (!it->second.in_eval(level,evaled,contextptr->globalcontextptr))
	  evaled=it->second;
	return true;
      }
    }
    if (local_eval(contextptr) && !localvalue->empty()){
      evaled=do_local_eval(*this,level,true);
      return true;
    }
    if (*quoted )
      return false;
    if (current_folder_name.type==_IDNT && current_folder_name._IDNTptr->value && current_folder_name._IDNTptr->value->type==_VECT){
      evaled=find_in_folder(*current_folder_name._IDNTptr->value->_VECTptr,orig);
      return (evaled!=orig);
    }
    if (value){
      evaled=value->eval(level,contextptr);
      return true;
    }
    // look in current directory for a value
    if ( secure_run || (!variables_are_files(contextptr)) || (access((*name+cas_suffixe).c_str(),R_OK))){
      evaled=orig;
      if (!local_eval(contextptr))
	evaled.subtype=_GLOBAL__EVAL;
      return true;
    }
    // set current value
    ifstream inf((*name+cas_suffixe).c_str());
    evaled=read1arg_from_stream(inf,list_one_letter__IDNT,contextptr);
    if (child_id)
      return true;
    value = new gen(evaled);
    evaled=evaled.eval(level,contextptr);
    return true;
  }

  void identificateur::push(int protection,const gen & e){
    if (!localvalue)
      settypeerr();
    localvalue->push_back(protection);
    localvalue->push_back(e);
  }

  string identificateur::print (GIAC_CONTEXT) const{
    if (!ref_count)
      return string("null__IDNT");
    if (*this==_IDNT_pi){
      switch (xcas_mode(contextptr)){
      case 1:
	return "Pi";
      case 2:
	return "PI";
      default:
	return string_pi;
      }
    }
    // index != sqrt(-1) wich has different notations
    if (xcas_mode(contextptr)==0){
      if (*name=="i")
	return "i_i";
    }
    else {
      if (*name=="I")
	return "i_i";
    }
    /*
    if (!localvalue->empty())
      return string("_") + *name ;        
    if (value)
      return string("~") + *name ;
    else
    */
    int pos=name->find(' '),taille=name->size();
    if (pos<0 || pos>=taille)
      return *name  ;
    return '`'+*name+'`';
  }

  ostream & operator << (ostream & os,const identificateur & s) { return os << s.print(context0);}

  int removecomments(const char * ss,char * ss2){
    int j=0,k=0;
    for (;ss[j];j++){
      if (ss[j]=='#'){
	ss2[k]=char(0); // end ss2 string
	break;
      }
      if (ss[j]>31){ // supress control chars
	ss2[k]=ss[j];
	k++;
      }
    }
    return k;
  }

  void identificateur::unassign(){
    if (value){
      delete(value);
      value = NULL;
    }
  }

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC
