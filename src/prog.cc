/* -*- mode:C++ ; compile-command: "g++-3.4 -I.. -I../include -g -c prog.cc -Wall" -*- */
#include "first.h"
/*
 *  Copyright (C) 2001,7 B. Parisse, Institut Fourier, 38402 St Martin d'Heres
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
#ifndef HAVE_NO_PWD_H
#include <pwd.h>
#endif
#include <stdexcept>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "prog.h"
#include "identificateur.h"
#include "symbolic.h"
#include "identificateur.h"
#include "usual.h"
#include "sym2poly.h"
#include "subst.h"
#include "plot.h"
#include "tex.h"
#include "input_parser.h"
#include "input_lexer.h"
#include "rpn.h"
#include "help.h"
#include "ti89.h" // for _unarchive_ti
#include "permu.h"
#include "modpoly.h"
#include "unary.h"
#include "input_lexer.h"
#include "maple.h"
// #include "input_parser.h"
#ifdef HAVE_LIBDL
#include <dlfcn.h>
#endif // HAVE_LIBDL

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  modules_tab giac_modules_tab;
  
  int prog_eval_level(GIAC_CONTEXT){
    if (int i=prog_eval_level_val(contextptr))
      return i;
    return std::max(1,eval_level(contextptr));
  }

  void check_secure(){
    if (secure_run)
      setsizeerr("Running in secure mode");
  }

  string indent(GIAC_CONTEXT){
    if (xcas_mode(contextptr)==3)
      return "\n:"+string(debug_ptr(contextptr)->indent_spaces,' ');
    else
      return " \n"+string(debug_ptr(contextptr)->indent_spaces,' ');
  }

  string indent2(GIAC_CONTEXT){
    return string(debug_ptr(contextptr)->indent_spaces,' ');
  }

  gen substsametoequal(const gen & g){
    return symbolic(at_subst,apply(g,sametoequal));
  }

  gen subssametoequal(const gen & g){
    return symbolic(at_subs,apply(g,sametoequal));
  }

  gen maplesubssametoequal(const gen & g){
    return symbolic(at_maple_subs,apply(g,sametoequal));
  }

  gen equaltosame(const gen & a){
    // full replacement of = by == has been commented to avoid
    // problems with tests like: if (limit(...,x=0,..))
    /*
    unary_function_ptr equaltosametab1[]={at_equal,at_subst,at_subs,at_maple_subs};
    vector<unary_function_ptr> substin(equaltosametab1,equaltosametab1+sizeof(equaltosametab1)/sizeof(unary_function_ptr));
    gen_op equaltosametab2[]={symb_same,substsametoequal,subssametoequal,maplesubssametoequal};
    vector<gen_op> substout(equaltosametab2,equaltosametab2+sizeof(equaltosametab2)/sizeof(gen_op));
    gen tmp=subst(a,substin,substout,true);
    return tmp;
    */
    if ( (a.type==_SYMB) && (a._SYMBptr->sommet==at_equal) )
      return symb_same(a._SYMBptr->feuille._VECTptr->front(),a._SYMBptr->feuille._VECTptr->back());
    else
      return a;
  }

    gen sametoequal(const gen & a){
        if ( (a.type==_SYMB) && (a._SYMBptr->sommet==at_same) )
            return equal(a._SYMBptr->feuille._VECTptr->front(),a._SYMBptr->feuille._VECTptr->back());
        else
            return a;
    }

  void increment_instruction(const const_iterateur & it0,const const_iterateur & itend,GIAC_CONTEXT){
    const_iterateur it=it0;
    for (;it!=itend;++it)
      increment_instruction(*it,contextptr);
  }

  void increment_instruction(const vecteur & v,GIAC_CONTEXT){
    const_iterateur it=v.begin(),itend=v.end();
    for (;it!=itend;++it)
      increment_instruction(*it,contextptr);
  }

  extern unary_function_ptr at_sialorssinon;
  extern unary_function_ptr at_pour;

  void increment_instruction(const gen & arg,GIAC_CONTEXT){
    // cerr << debug_ptr(contextptr)->current_instruction << " " << arg <<endl;
    ++debug_ptr(contextptr)->current_instruction;
    if (arg.type!=_SYMB)
      return;
    unary_function_ptr u=arg._SYMBptr->sommet;
    gen f=arg._SYMBptr->feuille;
    const unary_function_eval * uptr=dynamic_cast<const unary_function_eval *>(u.ptr);
    if (uptr && uptr->op==_ifte){
      --debug_ptr(contextptr)->current_instruction;
      increment_instruction(*f._VECTptr,contextptr);
      return;
    }
    if ( (u==at_local) || (uptr && uptr->op==_for) ){
      f=f._VECTptr->back();
      if (f.type!=_VECT){
	if (f.is_symb_of_sommet(at_bloc) && f._SYMBptr->feuille.type==_VECT)
	  increment_instruction(*f._SYMBptr->feuille._VECTptr,contextptr);
	else
	  increment_instruction(f,contextptr);
      }
      else 
	increment_instruction(*f._VECTptr,contextptr);
      return;
    }
    if (u==at_bloc){
      if (f.type!=_VECT)
	increment_instruction(f,contextptr);
      else
	increment_instruction(*f._VECTptr,contextptr);
      return;
    }
    if (u==at_try_catch){
      increment_instruction(f._VECTptr->front(),contextptr);
      increment_instruction(f._VECTptr->back(),contextptr);
    }
  }

  string concatenate(const vector<string> & v){
    vector<string>::const_iterator it=v.begin(),itend=v.end();
    string res;
    for (;it!=itend;++it){
      res=res + "  "+*it;
    }
    return res;
  }

  void debug_print(const vecteur & arg,vector<string> & v,GIAC_CONTEXT){
    const_iterateur it=arg.begin(),itend=arg.end();
    for (;it!=itend;++it){
      if (it->is_symb_of_sommet(at_program)){
	vector<string> tmp;
	debug_print(*it,tmp,contextptr);
	v.push_back(concatenate(tmp));
      }
      debug_print(*it,v,contextptr);
    }
  }

  void debug_print(const gen & e,vector<string>  & v,GIAC_CONTEXT){
    if (e.type!=_SYMB){
      v.push_back(indent2(contextptr)+e.print(contextptr));
      return ;
    }
    unary_function_ptr u=e._SYMBptr->sommet;
    gen f=e._SYMBptr->feuille;
    const unary_function_eval * uptr=dynamic_cast<const unary_function_eval *>(u.ptr);
    if (uptr && uptr->op==_ifte){
      string s=indent2(contextptr)+string("if(");
      vecteur w=*f._VECTptr;
      s += w.front().print(contextptr)+")";
      v.push_back(s);
      debug_ptr(contextptr)->indent_spaces += 1;
      debug_print(w[1],v,contextptr);
      debug_ptr(contextptr)->indent_spaces += 1;
      debug_print(w[2],v,contextptr);
      debug_ptr(contextptr)->indent_spaces -=2;
      return ;
    }
    if (u==at_local){
      string s(indent2(contextptr)+"local ");
      s += f._VECTptr->front().print(contextptr);
      v.push_back(s);
      debug_ptr(contextptr)->indent_spaces += 2;
      f=f._VECTptr->back();
      if (f.type!=_VECT)
	debug_print(f,v,contextptr);
      else
	debug_print(*f._VECTptr,v,contextptr);
      debug_ptr(contextptr)->indent_spaces -= 2;
      return;
    }
    if (uptr && uptr->op==_for){
      string s(indent2(contextptr)+"for(");
      vecteur w=*f._VECTptr;
      s += w[0].print(contextptr)+";"+w[1].print(contextptr)+";"+w[2].print(contextptr)+")";
      v.push_back(s);
      debug_ptr(contextptr)->indent_spaces += 2;
      f=f._VECTptr->back();
      if ((f.type==_SYMB) && (f._SYMBptr->sommet==at_bloc))
	f=f._SYMBptr->feuille;
      if (f.type!=_VECT)
	debug_print(f,v,contextptr);
      else
	debug_print(*f._VECTptr,v,contextptr);
      debug_ptr(contextptr)->indent_spaces -= 2;
      return;
    }
    if (u==at_bloc){
      v.push_back(string(indent2(contextptr)+"bloc"));
      debug_ptr(contextptr)->indent_spaces += 2;
      if (f.type!=_VECT)
	debug_print(f,v,contextptr);
      else
	debug_print(*f._VECTptr,v,contextptr);
      debug_ptr(contextptr)->indent_spaces -= 2;
      return;
    }
    if (u==at_try_catch){
      // cerr << f << endl;
      v.push_back(string(indent2(contextptr)+"try"));
      debug_ptr(contextptr)->indent_spaces += 1;
      debug_print(f._VECTptr->front(),v,contextptr);
      debug_ptr(contextptr)->indent_spaces += 1;
      debug_print(f._VECTptr->back(),v,contextptr);
      debug_ptr(contextptr)->indent_spaces -=2;
      return;
    }
    v.push_back(indent2(contextptr)+e.print(contextptr));
  }

  vecteur rm_checktype(const vecteur & v){
    vecteur addvars(v);
    iterateur it=addvars.begin(),itend=addvars.end();
    for (;it!=itend;++it){
      if (it->is_symb_of_sommet(at_check_type))
	*it=it->_SYMBptr->feuille._VECTptr->back();
      if (it->is_symb_of_sommet(at_sto))
	*it=it->_SYMBptr->feuille._VECTptr->back();	
    }
    return addvars;
  }
  // res1= list of assignation with =, res2= list of non-local vars
  void check_local_assign(const gen & g,vecteur & vars,vecteur & res1,vecteur & res2){
    if (g.is_symb_of_sommet(at_local)){
      gen &f=g._SYMBptr->feuille;
      if (f.type!=_VECT || f._VECTptr->size()!=2)
	return;
      vecteur addvars(rm_checktype(gen2vecteur(f._VECTptr->front())));
      vecteur newvars(mergevecteur(vars,addvars));
      check_local_assign(f._VECTptr->back(),newvars,res1,res2);
      return; 
    }
    if (g.is_symb_of_sommet(at_bloc) || 
	g.is_symb_of_sommet(at_for) ||
	g.is_symb_of_sommet(at_ifte)){
      check_local_assign(g._SYMBptr->feuille,vars,res1,res2);
      return;
    }
    if (g.is_symb_of_sommet(at_equal)){
      if (g._SYMBptr->feuille.type==_VECT && g._SYMBptr->feuille._VECTptr->size()==2 && g._SYMBptr->feuille._VECTptr->front().type!=_INT_ )
	res1.push_back(g);
      return;
    }
    if (g.is_symb_of_sommet(at_of)){
      gen & f=g._SYMBptr->feuille;
      if (f.type!=_VECT || f._VECTptr->size()!=2)
	return;
      check_local_assign(f._VECTptr->back(),vars,res1,res2);
      return;
    }
    if (g.type!=_VECT){
      vecteur l(*_lname(g)._VECTptr);
      const_iterateur it=l.begin(),itend=l.end();
      for (;it!=itend;++it){
	if (!equalposcomp(res2,*it) && !equalposcomp(vars,*it))
	  res2.push_back(*it);
      }
      return;
    }
    const_iterateur it=g._VECTptr->begin(),itend=g._VECTptr->end();
    for (;it!=itend;++it){
      check_local_assign(*it,vars,res1,res2);
    }
  }
  bool is_constant_idnt(const gen & g){
    return g==cst_pi || g==cst_euler_gamma;
  }
  // Return the names of variables that are not local in g
  // and the equality that are not used (warning = instead of := )
  string check_local_assign(const gen & g,GIAC_CONTEXT){
    string res;
    if (g.type==_VECT){
      const_iterateur it=g._VECTptr->begin(),itend=g._VECTptr->end();
      for (;it!=itend;++it)
	res += check_local_assign(*it,contextptr);
      return res;
    }
    if (g.is_symb_of_sommet(at_nodisp))
      return check_local_assign(g._SYMBptr->feuille,contextptr);
    if (g.is_symb_of_sommet(at_sto)){
      gen & f =g._SYMBptr->feuille;
      if (f.type!=_VECT || f._VECTptr->size()!=2)
	return res;
      res=check_local_assign(f._VECTptr->front(),contextptr);
      return res.substr(0,res.size()-1)+" compiling "+f._VECTptr->back().print(contextptr)+'\n';
    }
    if (!g.is_symb_of_sommet(at_program))
      return res;
    gen & f=g._SYMBptr->feuille;
    if (f.type!=_VECT || f._VECTptr->size()!=3)
      return "// Invalid program";
    vecteur & v =*f._VECTptr;
    vecteur vars=rm_checktype(gen2vecteur(v[0])),res1,res2(1,undef);
    gen prog=v.back();
    check_local_assign(prog,vars,res1,res2);
    int rs=res2.size();
    for (int i=0;i<rs;i++){
      if (is_constant_idnt(res2[i])){
	res2.erase(res2.begin()+i);
	--i; --rs;
      }
    }
    if (!res1.empty()){
      res="// Warning, assignation is :=, check these lines: ";
      const_iterateur it=res1.begin(),itend=res1.end();
      for (;it!=itend;++it){
	res += it->print(contextptr);
      }
      res +="\n";
    }
    if (res2.size()>1){
      res+="// Warning: ";
      const_iterateur it=res2.begin()+1,itend=res2.end();
      for (;it!=itend;++it){
	// pi already checked if (*it!=cst_pi)
	res += it->print(contextptr)+",";
      }
      res +=" declared as global variable(s)\n";
    }
    if (res.empty())
      return giac::first_error_line(contextptr)?"// Error(s)\n":"// Success\n";
    else
      return res;
  }

  string printascheck_type(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if ( (feuille.type!=_VECT) || (feuille._VECTptr->size()!=2) )
      return sommetstr+'('+feuille.print(contextptr)+')';
    string res(print_the_type(feuille._VECTptr->front().val,contextptr));
    res += ' '+feuille._VECTptr->back().print(contextptr);
    return res;
  }
  
  gen symb_check_type(const gen & args){
    return symbolic(at_check_type,args);
  }
  gen _check_type(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      return symb_check_type(args);
    if (args._VECTptr->size()!=2)
      setsizeerr("check_type must have 2 args");
    gen res=args._VECTptr->back();
    gen req=args._VECTptr->front();
    if (req.type!=_INT_) // FIXME check for matrix(...) or vector(...)
      return res;
    int type;
    switch (res.type){
    case _INT_:
      type=_ZINT;
      break;
    case _REAL:
      type=_DOUBLE_;
      break;
    default:
      type=res.type;
      break;
    }   
    if (req.val==_MAPLE_LIST){
      if (type==_VECT)
	return res;
      setsizeerr();
    }
    if (type==req.val)
      return res;
    if (type==_ZINT && type==(req.val &0xff) ){
      if (req.val==_POSINT && is_strictly_positive(res,contextptr)) 
	return res;
      if (req.val==_NEGINT && is_strictly_positive(-res,contextptr))
	return res;
      if (req.val==_NONPOSINT && is_positive(-res,contextptr))
	return res;
      if (req.val==_NONNEGINT && is_positive(res,contextptr))
	return res;
    }
    settypeerr("Argument should be of type "+print_the_type(args._VECTptr->front().val,contextptr));
    return 0;
  }
  const string _check_type_s("check_type");
  unary_function_unary __check_type(&symb_check_type,_check_type_s,&printascheck_type);
  unary_function_ptr at_check_type (&__check_type);

    gen symb_type(const gen & args){
    return symbolic(at_type,args);
  }
  gen _type(const gen & args){
    int type;
    switch (args.type){
    case _INT_:
      type=_ZINT;
      break;
    case _REAL:
      type=_DOUBLE_;
      break;
    default:
      if (args.is_symb_of_sommet(at_program))
	type=_FUNC;
      else
	type=args.type;
    }   
    gen tmp(type);
    tmp.subtype=1;
    return tmp;
  }
  const string _type_s("type");
  unary_function_unary __type(&_type,_type_s);
  unary_function_ptr at_type (&__type,0,true);

  gen _nop(const gen & a){
    if (a.type==_VECT && a.subtype==_SEQ__VECT){
      // Workaround so that sequences inside spreadsheet are saved as []
      gen tmp=a;
      tmp.subtype=0;
      return tmp;
    }
    return a;
  }
  const string _nop_s("nop");
  unary_function_unary __nop(&_nop,_nop_s);
  unary_function_ptr at_nop (&__nop,0,true);

  string printasinnerbloc(const gen & feuille,GIAC_CONTEXT){
    if ( (feuille.type==_SYMB) && feuille._SYMBptr->sommet==at_bloc)
      return printasinnerbloc(feuille._SYMBptr->feuille,contextptr);
    if (feuille.type!=_VECT)
      return indent(contextptr)+feuille.print(contextptr);
    const_iterateur it=feuille._VECTptr->begin(),itend=feuille._VECTptr->end();
    string res;
    if (it==itend)
      return res;
    for (;;){
      res += indent(contextptr)+it->print(contextptr);
      ++it;
      if (it==itend)
	return res;
      if (xcas_mode(contextptr)!=3)
	res += ";";
    }
  }

  void local_init(const gen & e,vecteur & non_init_vars,vecteur & initialisation_seq){
    vecteur v;
    if (e.type!=_VECT)
      v=vecteur(1,e);
    else
      v=*e._VECTptr;
    const_iterateur it=v.begin(),itend=v.end();
    for (;it!=itend;++it){
      if (it->type==_IDNT){
	non_init_vars.push_back(*it);
	continue;
      }
      if ( (it->type==_SYMB) && (it->_SYMBptr->sommet==at_sto)){
	non_init_vars.push_back(it->_SYMBptr->feuille._VECTptr->back());
	initialisation_seq.push_back(*it);
      }
    }
  }

  gen add_global(const gen & i,GIAC_CONTEXT){
    if (i.type==_IDNT)
      return identificateur("global_"+i.print(contextptr));
    setsizeerr("Proc Parameters");
    return 0;
  }

  string printasprogram(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if ( (feuille.type!=_VECT) || (feuille._VECTptr->size()!=3) )
      return sommetstr+'('+feuille.print(contextptr)+')';
    string res;
    if (xcas_mode(contextptr)==3)
      res="\n:lastprog";
    else
      res=" "; // was res=indent(contextptr);
    gen & feuille0=feuille._VECTptr->front();
    if (feuille0.type==_VECT && feuille0.subtype==_SEQ__VECT && feuille0._VECTptr->size()==1)
      res +="("+feuille0._VECTptr->front().print(contextptr)+")";
    else
      res +="("+feuille0.print(contextptr)+")";
    if (xcas_mode(contextptr)==3)
      res +="\n";
    else
      res += "->";
    bool test;
    string locals,inits;
    gen proc_args=feuille._VECTptr->front();
    vecteur vect,non_init_vars,initialisation_seq;
    if ((xcas_mode(contextptr)>0) && (feuille._VECTptr->back().type==_SYMB) && (feuille._VECTptr->back()._SYMBptr->sommet==at_local)){
        test=false;
        gen tmp=feuille._VECTptr->back()._SYMBptr->feuille;
	local_init(tmp._VECTptr->front(),non_init_vars,initialisation_seq);
	// For Maple add proc parameters to local vars
	if (xcas_mode(contextptr) ==1+_DECALAGE){
	  if (proc_args.type==_VECT){
	    vecteur v=*proc_args._VECTptr;
	    non_init_vars=mergevecteur(non_init_vars,v);
	    iterateur it=v.begin(),itend=v.end();
	    for (;it!=itend;++it){
	      gen tmp=add_global(*it,contextptr);
	      initialisation_seq.push_back(symb_sto(tmp,*it));
	      *it=tmp;
	    }
	    proc_args=gen(v,_SEQ__VECT);
	  }
	  else {
	    non_init_vars.push_back(proc_args);
	    gen tmp=add_global(proc_args,contextptr);
	    initialisation_seq.push_back(symb_sto(tmp,proc_args));
	    proc_args=tmp;
	  }
	}
	if (!non_init_vars.empty()){
	  if (xcas_mode(contextptr)==3)
	    locals=indent(contextptr)+"Local "+printinner_VECT(non_init_vars,_SEQ__VECT,contextptr);
	  else
	    locals=indent(contextptr)+"  local "+printinner_VECT(non_init_vars,_SEQ__VECT,contextptr)+";";
	}
	inits=printasinnerbloc(gen(initialisation_seq,_SEQ__VECT),contextptr);
	if (tmp._VECTptr->back().type==_VECT)
	  vect=*tmp._VECTptr->back()._VECTptr;
	else
	  vect=makevecteur(tmp._VECTptr->back());
    }
    else {
        test=(feuille._VECTptr->back().type!=_VECT ||feuille._VECTptr->back().subtype );
        if (!test)
            vect=*feuille._VECTptr->back()._VECTptr;
    }
    if (test){
      if (xcas_mode(contextptr)==3)
	return res+":Func "+feuille._VECTptr->back().print(contextptr)+"\n:EndFunc\n";
      return res+feuille._VECTptr->back().print(contextptr);
    }
    if (xcas_mode(contextptr)>0){
      if (xcas_mode(contextptr)==3)
	res+=":Func"+locals;
      else {
	res="proc("+proc_args.print(contextptr)+")"+locals;
	if (xcas_mode(contextptr)==2)
	  res +=indent(contextptr)+"begin ";
	if (inits.size()) 
	  res += indent(contextptr)+inits+";";
      }
    }
    else
      res +="{";
    const_iterateur it=vect.begin(),itend=vect.end();
    debug_ptr(contextptr)->indent_spaces +=2;
    for (;;){
      if (xcas_mode(contextptr)==3)
	res += indent(contextptr)+it->print(contextptr);
      else
	res += indent(contextptr)+it->print(contextptr);
      ++it;
      if (it==itend){
	debug_ptr(contextptr)->indent_spaces -=2;
	if (xcas_mode(contextptr)!=3)
	  res += "; "+indent(contextptr);
	switch (xcas_mode(contextptr)){
	case 0:
	  res += "}";
	  break;
	case 1: case 1+_DECALAGE:
	  res+=indent(contextptr)+"end;";
	  break;
	case 2:
	  return res+=indent(contextptr)+"end_proc;";
	  break;
	case 3:
	  return res+=indent(contextptr)+"EndFunc\n";
	}
	break;
      }
      else {
	if (xcas_mode(contextptr)!=3)
	  res +="; ";
      }
    }
    return res;
  }

  string texprintasprogram(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if (latex_format(contextptr)==1){
      return printasprogram(feuille,sommetstr,contextptr);
    }
    string s("\\parbox{12cm}{\\tt ");
    s += translate_underscore(printasprogram(feuille,sommetstr,contextptr));
    s+=" }";
    return s;
  }

  void replace_keywords(const gen & a,const gen & b,gen & newa,gen & newb,GIAC_CONTEXT){
    // Check that variables in a are really variables, otherwise print
    // the var and make new variables
    vecteur newv(gen2vecteur(a));
    vecteur v1,v2;
    iterateur it=newv.begin(),itend=newv.end();
    for (;it!=itend;++it){
      if (it->is_symb_of_sommet(at_sto) || it->is_symb_of_sommet(at_check_type)) // FIXME check 1st arg too
	continue;
      if (it->type!=_IDNT && it->type!=_CPLX){
	v1.push_back(*it);
	string s=gen2string(*it);
	int ss=s.size();
	if (ss>2 && s[0]=='\'' && s[ss-1]=='\'')
	  s=s.substr(1,ss-2);
	sym_tab::const_iterator i = syms().find(s);
	if (i == syms().end()) {
	  *it = *(new identificateur(s));
	  syms()[s] = *it;
	} else {
	  // std::cerr << "lexer" << s << endl;
	  *it = i->second;
	}
	v2.push_back(*it);
      }
    }
    newa=gen(newv,_SEQ__VECT);
    if (v1.empty())
      newb=b;
    else
      newb=quotesubst(b,v1,v2,contextptr);
  }

  symbolic symb_program_sto(const gen & a,const gen & b,const gen & c,const gen & d,bool embedd,GIAC_CONTEXT){
    *logptr(contextptr) << "// Parsing " << d << endl;
    gen newa,newc;
    replace_keywords(a,((embedd&&c.type==_VECT)?makevecteur(c):c),newa,newc,contextptr);
    symbolic g=symbolic(at_program,makevecteur(newa,b,newc));
    g=symbolic(at_sto,makevecteur(g,d));
    *logptr(contextptr) << check_local_assign(g,contextptr) ;
    return g;
  }
  symbolic symb_program(const gen & a,const gen & b,const gen & c,GIAC_CONTEXT){
    gen newa,newc;
    replace_keywords(a,c,newa,newc,contextptr);
    symbolic g=symbolic(at_program,makevecteur(newa,b,newc));
    *logptr(contextptr) << check_local_assign(g,contextptr) ;
    return g;
  }
  symbolic symb_program(const gen & args){
    return symbolic(at_program,args);
  }
  void lidnt_prog(const gen & g,vecteur & res);
  void lidnt_prog(const vecteur & v,vecteur & res){
    const_iterateur it=v.begin(),itend=v.end();
    for (;it!=itend;++it)
      lidnt_prog(*it,res);
  }
  void lidnt_prog(const gen & g,vecteur & res){
    switch (g.type){
    case _VECT:
      lidnt_prog(*g._VECTptr,res);
      break;
    case _IDNT:
      if (!equalposcomp(res,g))
	res.push_back(g);
      break;
    case _SYMB:
      /* if (g._SYMBptr->sommet==at_at || g._SYMBptr->sommet==at_of )
	lidnt_prog(g._SYMBptr->feuille._VECTptr->back(),res);
	else */
	lidnt_prog(g._SYMBptr->feuille,res);
      break;
    }
  }

  void local_vars(const gen & g,vecteur & v,GIAC_CONTEXT){
    if (g.type==_VECT){
      const_iterateur it=g._VECTptr->begin(),itend=g._VECTptr->end();
      for (;it!=itend;++it){
	local_vars(*it,v,contextptr);
      }
      return ;
    }
    if (g.type!=_SYMB)
      return;
    if (g._SYMBptr->sommet==at_local && g._SYMBptr->feuille.type==_VECT){
      vecteur & w = *g._SYMBptr->feuille._VECTptr;
      v=mergevecteur(v,gen2vecteur(w[0]));
      local_vars(w[1],v,contextptr);
    }
    else
      local_vars(g._SYMBptr->feuille,v,contextptr);
  }

  gen quote_program(const gen & args,GIAC_CONTEXT){
    // return symb_program(args);
    // g:=unapply(p ->translation(w,p),w);g(1)
    // Necessite d'evaluer les arguments des programmes a l'interieur d'un programme.
    // Mais il ne faut pas evaluer les variables declarees comme locales!!
    bool in_prog=debug_ptr(contextptr)->sst_at_stack.size();
    // ?? Subst all variables except arguments
    if (!in_prog || args.type!=_VECT || args._VECTptr->size()!=3)
      return symb_program(args);
    vecteur & v = *args._VECTptr;
    vecteur vars(gen2vecteur(v[0]));
    int s=vars.size(); // s vars not subst-ed
    lidnt_prog(v[2],vars);
    vars=vecteur(vars.begin()+s,vars.end());
    // Remove local variables from the list
    vecteur tmpvar,resvar;
    local_vars(v[2],tmpvar,contextptr); 
    const_iterateur it=vars.begin(),itend=vars.end();
    for (;it!=itend;++it){
      if (!equalposcomp(tmpvar,*it))
	resvar.push_back(*it);
    }
    gen tmp=gen(resvar).eval(1,contextptr);
    vecteur varsub(*tmp._VECTptr);
    return symbolic(at_program,quotesubst(args,resvar,varsub,contextptr));
  }
  bool is_return(const gen & g,gen & newres){
    if ((g.type==_SYMB) && (g._SYMBptr->sommet==at_return) ){
      gen tmp = g._SYMBptr->feuille;
      is_return(tmp,newres);
      return true;
    }
    if ( (g.type==_VECT && g.subtype ==_SEQ__VECT && g._VECTptr->size()==1) )
      return is_return(g._VECTptr->front(),newres);
    newres=g;
    return false;
  }
  void adjust_sst_at(const gen & name,GIAC_CONTEXT){
    debug_ptr(contextptr)->sst_at.clear();
    const_iterateur it=debug_ptr(contextptr)->debug_breakpoint.begin(),itend=debug_ptr(contextptr)->debug_breakpoint.end();
    for (;it!=itend;++it){
      if (it->_VECTptr->front()==name)
	debug_ptr(contextptr)->sst_at.push_back(it->_VECTptr->back().val);
    }
  }
  gen _program(const gen & args,const gen & name,const context * contextptr){
    if (args.type!=_VECT)
      return args.eval(prog_eval_level(contextptr),contextptr);
    // set breakpoints
    if ( int(debug_ptr(contextptr)->sst_at_stack.size())> MAX_RECURSION_LEVEL)
      setsizeerr("Too many recursions");
    debug_ptr(contextptr)->sst_at_stack.push_back(debug_ptr(contextptr)->sst_at);
    debug_ptr(contextptr)->sst_at.clear();
    if (name.type==_IDNT)
      adjust_sst_at(name,contextptr);
    debug_ptr(contextptr)->current_instruction_stack.push_back(debug_ptr(contextptr)->current_instruction);
    debug_ptr(contextptr)->current_instruction=0;
    bool save_sst_mode = debug_ptr(contextptr)->sst_mode ;
    // *logptr(contextptr) << "Entering prog " << args << " " << debug_ptr(contextptr)->sst_in_mode << endl;
    if (debug_ptr(contextptr)->sst_in_mode){
      debug_ptr(contextptr)->sst_in_mode=false;
      debug_ptr(contextptr)->sst_mode=true;
    }
    else
      debug_ptr(contextptr)->sst_mode=false;
    // Bind local var
    if (args._VECTptr->size()!=3)
      setsizeerr();
    gen vars=args._VECTptr->front();
    if (vars.type!=_VECT)
      vars=makevecteur(vars);
    gen values=(*(args._VECTptr))[1];
    if (values.type!=_VECT || values.subtype!=_SEQ__VECT || (vars._VECTptr->size()==1 && values._VECTptr->size()!=1))
      values=makevecteur(values);
    // *logptr(contextptr) << vars << " " << values << endl;
    debug_ptr(contextptr)->args_stack.push_back(mergevecteur(vecteur(1,name),*values._VECTptr));
    gen prog=args._VECTptr->back(),save_debug_info((*debug_ptr(contextptr)->debug_info_ptr));
    // removed sst test so that when a breakpoint is evaled
    // the correct info is displayed
    (*debug_ptr(contextptr)->debug_info_ptr)=prog;
    (*debug_ptr(contextptr)->fast_debug_info_ptr)=prog;
    int protect;
    context * newcontextptr=(context *)contextptr;
    if (!vars._VECTptr->empty())
      protect=bind(*values._VECTptr,*vars._VECTptr,newcontextptr);
    gen res,newres;
    try {
      if (prog.type!=_VECT || prog.subtype){
	++debug_ptr(newcontextptr)->current_instruction;
	if (debug_ptr(newcontextptr)->debug_mode)
	  debug_loop(res,newcontextptr);
	res=prog.eval(prog_eval_level(newcontextptr),newcontextptr);
	if (is_return(res,newres))
	  res=newres;
      }
      else {
	const_iterateur it=prog._VECTptr->begin(),itend=prog._VECTptr->end();
	bool findlabel=false;
	gen label;
	for (;it!=itend;++it){
	  ++debug_ptr(newcontextptr)->current_instruction;
	  if (debug_ptr(newcontextptr)->debug_mode)
	    debug_loop(res,newcontextptr);
	  if (!findlabel)
	    res=it->eval(prog_eval_level(newcontextptr),newcontextptr);
	  else
	    res=*it;
	  if (res.type==_SYMB){
	    unary_function_ptr & u=res._SYMBptr->sommet;
	    if (findlabel && u==at_label && label==res._SYMBptr->feuille)
	      findlabel=false;
	    if (!findlabel && u==at_goto){
	      findlabel=true;
	      label=res._SYMBptr->feuille;
	    }
	  }
	  if (findlabel && it+1==itend)
	    it=prog._VECTptr->begin()-1;
	  if (!findlabel && is_return(res,newres)){
	    res=newres;
	    break;
	  }
	}
      }
    } // end try
    catch (std::runtime_error & e){
      if (!vars._VECTptr->empty())
	leave(protect,*vars._VECTptr,newcontextptr);
      throw(std::runtime_error(e.what()));
    }
    if (!vars._VECTptr->empty())
      leave(protect,*vars._VECTptr,newcontextptr);
    debug_ptr(contextptr)->args_stack.pop_back();
    // *logptr(contextptr) << "Leaving " << args << endl;
    if (!debug_ptr(contextptr)->sst_at_stack.empty()){
      debug_ptr(contextptr)->sst_at=debug_ptr(contextptr)->sst_at_stack.back();
      debug_ptr(contextptr)->sst_at_stack.pop_back();
    }
    if (!debug_ptr(contextptr)->current_instruction_stack.empty()){
      debug_ptr(contextptr)->current_instruction=debug_ptr(contextptr)->current_instruction_stack.back();
      debug_ptr(contextptr)->current_instruction_stack.pop_back();
    }
    debug_ptr(contextptr)->sst_mode=save_sst_mode;
    if (debug_ptr(contextptr)->current_instruction_stack.empty())
      debug_ptr(contextptr)->debug_mode=false;
    (*debug_ptr(contextptr)->debug_info_ptr)=save_debug_info;
    (*debug_ptr(contextptr)->fast_debug_info_ptr)=save_debug_info;
    return res;
  }
  const string _program_s("program");
  unary_function_eval __program(&quote_program,_program_s,&printasprogram,&texprintasprogram);
  unary_function_ptr at_program (&__program,_QUOTE_ARGUMENTS);

  string printasbloc(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if ( (feuille.type!=_VECT) )
      return "{"+feuille.print(contextptr)+";}";
    const_iterateur it=feuille._VECTptr->begin(),itend=feuille._VECTptr->end();
    string res("{");
    if (xcas_mode(contextptr)>0){
      if (xcas_mode(contextptr)==3)
	res="";
      else
	res=indent(contextptr)+"begin";
    }
    debug_ptr(contextptr)->indent_spaces +=2;
    for (;;){
      if (xcas_mode(contextptr)==3)
	res += indent(contextptr)+it->print(contextptr);
      else
	res += indent(contextptr)+it->print(contextptr);
      ++it;
      if (it==itend){
	debug_ptr(contextptr)->indent_spaces -=2;
	if (xcas_mode(contextptr)==3)
	  break;
	res += "; "+indent(contextptr);
	if (xcas_mode(contextptr)>0)
	  res += indent(contextptr)+"end";
	else
	  res += "}";
	break;
      }
      else {
	if (xcas_mode(contextptr)!=3)
	  res +="; ";
      }
    }
    return res;
  }
  symbolic symb_bloc(const gen & args){
    return symbolic(at_bloc,args);
  }
  gen _bloc(const gen & prog,GIAC_CONTEXT){
    gen res,label;
    bool findlabel=false;
    if (prog.type!=_VECT){
      ++debug_ptr(contextptr)->current_instruction;
      if (debug_ptr(contextptr)->debug_mode)
	debug_loop(undef,contextptr);
      return prog.eval(eval_level(contextptr),contextptr);
    }
    else {
      const_iterateur it=prog._VECTptr->begin(),itend=prog._VECTptr->end();
      for (;it!=itend;++it){
	++debug_ptr(contextptr)->current_instruction;
	if (debug_ptr(contextptr)->debug_mode)
	  debug_loop(res,contextptr);
	if (!findlabel)
	  res=it->eval(eval_level(contextptr),contextptr);
	else 
	  res=*it;
	if (res.type==_SYMB){
	  unary_function_ptr & u=res._SYMBptr->sommet;
	  if (!findlabel && u==at_return) {
	    increment_instruction(it+1,itend,contextptr);
	    return res; // it->eval(eval_level(contextptr),contextptr);
	  }
	  if (!findlabel && u==at_goto){
	    findlabel=true;
	    label=res._SYMBptr->feuille;
	  }
	  if ( u==at_label && label==res._SYMBptr->feuille )
	    findlabel=false;
	}
	// restart the bloc if needed
	if (findlabel && it+1==itend)
	  it=prog._VECTptr->begin()-1;
      }
    }
    return res;
  }
  const string _bloc_s("bloc");
  unary_function_eval __bloc(&_bloc,_bloc_s,&printasbloc);
  unary_function_ptr at_bloc (&__bloc,_QUOTE_ARGUMENTS);

  // test
  string printasifte(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if ( (feuille.type!=_VECT) || (feuille._VECTptr->size()!=3) )
      return sommetstr+'('+feuille.print(contextptr)+')';
    const_iterateur it=feuille._VECTptr->begin();//,itend=feuille._VECTptr->end();
    string res("if ");
    if (xcas_mode(contextptr)==3)
      res="If ";
    if (xcas_mode(contextptr)>0)
      res += sametoequal(*it).print(contextptr);
    else
      res += "("+it->print(contextptr);
    ++it;
    if (xcas_mode(contextptr)>0){
      if (xcas_mode(contextptr)==3){
	if (is_undef(*(it+1)) && (it->type!=_SYMB || it->_SYMBptr->sommet!=at_bloc) )
	  return res + indent(contextptr)+"  "+it->print(contextptr);
	res += " Then ";
      }
      else
	res += " then ";
    }
    else
      res += ") ";
    debug_ptr(contextptr)->indent_spaces +=2;
    if ((xcas_mode(contextptr)>0) && (it->type==_SYMB) && (it->_SYMBptr->sommet==at_bloc))
      res += printasinnerbloc(it->_SYMBptr->feuille,contextptr);
    else
      res += it->print(contextptr) ;
    debug_ptr(contextptr)->indent_spaces -=2;
    ++it;
    while ( (xcas_mode(contextptr)>0) && (it->type==_SYMB) && (it->_SYMBptr->sommet==at_ifte) ){
      if (xcas_mode(contextptr)==3)
	res += indent(contextptr)+"Elseif ";
      else
	res+= " elif ";
      it=it->_SYMBptr->feuille._VECTptr->begin();
      res += it->print(contextptr);
      if (xcas_mode(contextptr)==3)
	res += " Then ";
      else
	res += " then ";
      ++it;
      debug_ptr(contextptr)->indent_spaces +=2;
      if ((it->type==_SYMB) && (it->_SYMBptr->sommet==at_bloc))
	res += printasinnerbloc(it->_SYMBptr->feuille,contextptr);
      else
	res += it->print(contextptr) ;
      debug_ptr(contextptr)->indent_spaces -=2;
      ++it;
    }
    if (!is_undef(*it)){
      if ((xcas_mode(contextptr)<=0) && (res[res.size()-1]=='}'))
	res += indent(contextptr)+" else ";
      else {
	if (xcas_mode(contextptr)<=0)
	  res +=";"; 
	if (xcas_mode(contextptr)==3)
	  res += indent(contextptr)+"Else ";
	else
	  res+= " else ";
      }
      debug_ptr(contextptr)->indent_spaces +=2;
      if ((xcas_mode(contextptr)>0) && (it->type==_SYMB) && (it->_SYMBptr->sommet==at_bloc))
	res += printasinnerbloc(it->_SYMBptr->feuille,contextptr);
      else {
	res += it->print(contextptr) ;
	if (!xcas_mode(contextptr))
	  res += ";";
      }
      debug_ptr(contextptr)->indent_spaces -=2;
    }
    // FIXME? NO ; AT END OF IF
    if ((xcas_mode(contextptr)<=0) && (res[res.size()-1]!='}'))
        res +=" ";
    if ( (xcas_mode(contextptr) ==1) || (xcas_mode(contextptr) == 1+_DECALAGE) )
      res += indent(contextptr)+ "fi ";
    if (xcas_mode(contextptr)==2)
      res += indent(contextptr)+ "end_if ";
    if (xcas_mode(contextptr)==3)
      res += indent(contextptr)+"EndIf ";
    return res;
  }
  symbolic symb_ifte(const gen & e){
    return symbolic(at_ifte,e);
  }
  symbolic symb_ifte(const gen & test,const gen & oui, const gen & non){
    return symbolic(at_ifte,makevecteur(test,oui,non));
  }

  gen ifte(const gen & args,bool isifte,const context * contextptr){
    if (args._VECTptr->size()!=3)
      setsizeerr("Ifte must have 3 args");
    gen test=args._VECTptr->front();
    test=equaltosame(test.eval(eval_level(contextptr),contextptr)).eval(eval_level(contextptr),contextptr);
    if (!is_integer(test)){
      test=test.evalf_double(eval_level(contextptr),contextptr);
      if ( (test.type!=_DOUBLE_) && (test.type!=_CPLX) ){
	if (isifte)
	  setsizeerr("Ifte: Unable to check test"); 
	else
	  return symbolic(at_when,args);
      }
    }
    gen clause_vraie=(*(args._VECTptr))[1];
    gen clause_fausse=args._VECTptr->back();
    // *logptr(contextptr) << "Ifte " << debug_ptr(contextptr)->current_instruction << endl ;
    gen res;
    if (is_zero(test)){ // test false, do the else part
      if (isifte){
	increment_instruction(clause_vraie,contextptr);
	// *logptr(contextptr) << "Else " << debug_ptr(contextptr)->current_instruction << endl ;
	++debug_ptr(contextptr)->current_instruction;
	if (debug_ptr(contextptr)->debug_mode)
	  debug_loop(test,contextptr);
      }
      res=clause_fausse.eval(eval_level(contextptr),contextptr);
      // *logptr(contextptr) << "Else " << debug_ptr(contextptr)->current_instruction << endl ;
    }
    else { // test true, do the then part
      if (isifte){
	++debug_ptr(contextptr)->current_instruction;
	if (debug_ptr(contextptr)->debug_mode)
	  debug_loop(test,contextptr);
      }
      res=clause_vraie.eval(eval_level(contextptr),contextptr);
      // *logptr(contextptr) << "Then " << debug_ptr(contextptr)->current_instruction << endl ;
      if (isifte)
	increment_instruction(clause_fausse,contextptr);
      // *logptr(contextptr) << "Then " << debug_ptr(contextptr)->current_instruction << endl ;
    }
    return res;
  }
  gen _ifte(const gen & args,const context * contextptr){
    if (args.type!=_VECT)
      return symb_ifte(args);
    return ifte(args,true,contextptr);
  }
  const string _ifte_s("ifte");
  unary_function_eval __ifte(&_ifte,_ifte_s,&printasifte);
  unary_function_ptr at_ifte (&__ifte,_QUOTE_ARGUMENTS);

  gen _evalb(const gen & args,GIAC_CONTEXT){
    gen test=equaltosame(args);
    test=normal(test,contextptr);
    test=test.eval(eval_level(contextptr),contextptr);
    test=test.evalf_double(1,contextptr);
    if ( (test.type!=_DOUBLE_) && (test.type!=_CPLX) )
      return symbolic(at_evalb,args);
    if (is_zero(test))
      return zero;
    return plus_one;
  }
  const string _evalb_s("evalb");
  unary_function_eval __evalb(&_evalb,_evalb_s);
  unary_function_ptr at_evalb (&__evalb,_QUOTE_ARGUMENTS,true);

  const string _maple_if_s("if");
  unary_function_eval __maple_if(&_ifte,_maple_if_s,&printasifte);
  unary_function_ptr at_maple_if (&__maple_if,_QUOTE_ARGUMENTS);

  string printaswhen(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if (xcas_mode(contextptr)||feuille.type!=_VECT || feuille._VECTptr->size()!=3)
      return sommetstr+"("+feuille.print(contextptr)+")";
    vecteur & v=*feuille._VECTptr;
    return "(("+v[0].print(contextptr)+")? "+v[1].print(contextptr)+" : "+v[2].print(contextptr)+")";
  }
  gen symb_when(const gen & t,const gen & a,const gen & b){
    return symbolic(at_when,gen(makevecteur(t,a,b),_SEQ__VECT));
  }
  gen _when(const gen & args,const context * contextptr){
    if (args.type!=_VECT)
      return symbolic(at_when,args);
    vecteur & v=*args._VECTptr;
    if (v.size()==3){
      gen res=ifte(args,false,contextptr);
      return res;
    }
    if (v.size()!=4)
      settypeerr();
    gen res=ifte(vecteur(v.begin(),v.begin()+3),false,contextptr);
    if (res.type==_SYMB && res._SYMBptr->sommet==at_when)
      return v[3];
    return res;
  }
  const string _when_s("when");
  unary_function_eval __when(&_when,_when_s,&printaswhen);
  unary_function_ptr at_when (&__when,_QUOTE_ARGUMENTS,true);

  // convert back increment and decrement to sto
  gen from_increment(const gen & g){
    int type=0;
    if (g.is_symb_of_sommet(at_increment))
      type=1;
    if (g.is_symb_of_sommet(at_decrement))
      type=-1;
    if (type){
      gen & f =g._SYMBptr->feuille;
      if (f.type!=_VECT)
	return symbolic(at_sto,makevecteur(symbolic(at_plus,makevecteur(f,type)),f));
      vecteur & v = *f._VECTptr;
      if (v.size()!=2)
	setsizeerr();
      return symbolic(at_sto,makevecteur(symbolic(at_plus,makevecteur(v[0],type*v[1])),v[0]));
    }
    return g;
  }

  // loop
  string printasfor(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if ( (feuille.type!=_VECT) || (feuille._VECTptr->size()!=4) )
      return sommetstr+'('+feuille.print(contextptr)+')';
    int maplemode=xcas_mode(contextptr) & 0x07;
    const_iterateur it=feuille._VECTptr->begin();//,itend=feuille._VECTptr->end();
    string res;
    gen inc(from_increment(*(it+2)));
    if (is_integer(*it) && is_integer(*(it+2))){
      ++it;
      if (maplemode>0){
	if (maplemode==3 && is_one(*it) ){
	  it += 2;
	  res="Loop";
	  if ((it->type==_SYMB) && (it->_SYMBptr->sommet==at_bloc))
	    res += printasinnerbloc(it->_SYMBptr->feuille,contextptr);
	  else
	    res += it->print(contextptr) ;
	  return res+indent(contextptr)+"EndLoop";
	}
	if (maplemode==3)
	  res="While "+ sametoequal(*it).print(contextptr);
	else
	  res ="while " + sametoequal(*it).print(contextptr) + " do ";
      }
      else
	res = "while("+it->print(contextptr)+")";
      ++it;
      ++it;
      debug_ptr(contextptr)->indent_spaces += 2;
      if ((maplemode>0) && (it->type==_SYMB) && (it->_SYMBptr->sommet==at_bloc))
	res += printasinnerbloc(it->_SYMBptr->feuille,contextptr)+";";
      else
	res += it->print(contextptr) ;
      debug_ptr(contextptr)->indent_spaces -= 2;
      if (maplemode==1)
	return res+indent(contextptr)+" od;";
      if (maplemode==2)
	return res+indent(contextptr)+" end_while;";
      if (maplemode==3)
	return res+indent(contextptr)+" EndWhile";
    }
    else {  
      if (maplemode>0){// pb for generic loops for Maple translation
	gen inc=from_increment(*(it+2));
	if ( (it->type!=_SYMB) || (it->_SYMBptr->sommet!=at_sto) || (inc.type!=_SYMB) || inc._SYMBptr->sommet!=at_sto || (it->_SYMBptr->feuille._VECTptr->back()!=inc._SYMBptr->feuille._VECTptr->back()) )
	  return "Maple/Mupad/TI For: unable to convert";
	gen var_name=it->_SYMBptr->feuille._VECTptr->back();
	gen step=normal(inc._SYMBptr->feuille._VECTptr->front()-var_name,contextptr);
	gen condition=*(it+1),limite=plus_inf;
	if (is_positive(-step,contextptr)) 
	  limite=minus_inf;
	bool simple_loop=false,strict=true,ascending=true;
	if (condition.type==_SYMB){
	  unary_function_ptr op=condition._SYMBptr->sommet;
	  if (condition._SYMBptr->feuille.type==_VECT){
	    if (op==at_inferieur_strict)
	      simple_loop=true;
	    if (op==at_inferieur_egal){
	      strict=false;
	      simple_loop=true;
	    }
	    if (op==at_superieur_strict){
	      simple_loop=(maplemode>=2);
	      ascending=false;
	    }
	    if (op==at_superieur_egal){
	      simple_loop=(maplemode>=2);
	      ascending=false;
	      strict=false;
	    }
	  }
	  if (simple_loop){
	    simple_loop=(condition._SYMBptr->feuille._VECTptr->front()==var_name);
	    limite=condition._SYMBptr->feuille._VECTptr->back();
	  }
	}
	if (maplemode==3)
	  res="For ";
	else
	  res ="for ";
	res += var_name.print(contextptr);
	if (maplemode==3)
	  res += ",";
	else
	  res += " from ";
	res += it->_SYMBptr->feuille._VECTptr->front().print(contextptr);
	if (maplemode==3){
	  res += ","+limite.print(contextptr);
	  if (!is_one(step))
	    res += ","+step.print(contextptr);
	  if (!simple_loop)
	    res += indent(contextptr)+"If not("+(it+1)->print(contextptr)+")"+indent(contextptr)+"Exit";
	}
	else {
	  gen absstep=step;
	  if (simple_loop){
	    absstep = abs(step,contextptr); 
	    if (ascending)
	      res += " to ";
	    else
	      res += " downto ";
	    res += limite.print(contextptr);
	    if (!strict){
	      if (ascending)
		res +="+";
	      else
		res += "-";
	      res += absstep.print(contextptr);
	      res += "/2";
	    }
	  }
	  if (!is_one(absstep)){
	    if (maplemode==2)
	      res+= " step ";
	    else
	      res += " by ";
	    res += step.print(contextptr);
	  }
	  if (!simple_loop){
	    res += " while ";
	    res += (it+1)->print(contextptr);
	  }
	  res += " do ";
	}
	it += 3;
	if ((it->type==_SYMB) && (it->_SYMBptr->sommet==at_bloc))
	  res += printasinnerbloc(it->_SYMBptr->feuille,contextptr)+";";
	else
	  res += it->print(contextptr) ;
	if (maplemode==1)
	  return res + indent(contextptr)+" od;";
	if (maplemode==2)
	  return res + indent(contextptr)+" end_for;";
	if (maplemode==3)
	  return res + indent(contextptr)+"EndFor";
      }
      res="for (";
      res += it->print(contextptr) + ';';
      ++it;
      res += it->print(contextptr) + ';';
      ++it;
      res += it->print(contextptr) + ") ";
      ++it;
      debug_ptr(contextptr)->indent_spaces += 2;
      res += it->print(contextptr) ;
      debug_ptr(contextptr)->indent_spaces -= 2;
    }
    if (res[res.size()-1]!='}')
      res += "; ";
    return res;
  }
  symbolic symb_for(const gen & e){
    return symbolic(at_for,e);
  }
  symbolic symb_for(const gen & a,const gen & b,const gen & c,const gen & d){
    return symbolic(at_for,makevecteur(a,b,c,d));
  }
  
  gen to_increment(const gen & g){
    if (!g.is_symb_of_sommet(at_sto))
      return g;
    gen & f =g._SYMBptr->feuille;
    if (f.type!=_VECT || f._VECTptr->size()!=2)
      return g;
    gen & a = f._VECTptr->front();
    gen & b = f._VECTptr->back();
    if (b.type!=_IDNT || a.type!=_SYMB)
      return g;
    gen & af=a._SYMBptr->feuille;
    if (af.type!=_VECT || af._VECTptr->empty())
      return g;
    vecteur & av= *af._VECTptr;
    int s=av.size();
    int type=0;
    if (a.is_symb_of_sommet(at_plus))
      type=1;
    // there was a wrong test with at_minus for -= (type=-1)
    if (type && av.front()==b){
      if (s==2){
	if (is_one(av.back()))
	  return symbolic(type==1?at_increment:at_decrement,b);
	if (is_minus_one(av.back()))
	  return symbolic(type==1?at_decrement:at_increment,b);
	return symbolic(type==1?at_increment:at_decrement,makevecteur(b,av.back()));
      }
      if (type)
	return symbolic(at_increment,makevecteur(b,symbolic(at_plus,vecteur(av.begin()+1,av.end()))));
    }
    return g;
  }
  bool ck_is_one(const gen & g){
    if (is_one(g))
      return true;
    if (g.type>_POLY)
      setsizeerr("Unable to eval test in loop : "+g.print());
    return false;
  }
  gen _for(const gen & args,const context * contextptr)
  {
    if (args.type!=_VECT)
      return symb_for(args);
    if (args._VECTptr->size()!=4)
      setsizeerr("For must have 4 args");
    // Initialization
    gen initialisation=args._VECTptr->front();
    // add assigned variables to be local
    bool bound=false;
    vecteur loop_var;
    int protect=0;
    context * newcontextptr=(context * ) contextptr;
    if ( (initialisation.type==_SYMB) && (initialisation._SYMBptr->sommet==at_sto)){
      gen variable=initialisation._SYMBptr->feuille._VECTptr->back();
      if (variable.type==_IDNT){
	if (contextptr==context0 && (xcas_mode(contextptr)!=1 && (variable._IDNTptr->localvalue->empty() || (*variable._IDNTptr->localvalue)[variable._IDNTptr->localvalue->size()-2].val<protection_level-1) ) ){
	  bound=true;
	  loop_var=makevecteur(variable);
	  protect=bind(makevecteur(zero),loop_var,newcontextptr);
	}
      }
      else 
	throw(std::runtime_error("Invalid loop index (hint: i=sqrt(-1)!)"));
    }
    gen test=(*(args._VECTptr))[1];
    if ((test.type==_SYMB) && (test._SYMBptr->sommet==at_equal))
      test = symb_same(test._SYMBptr->feuille._VECTptr->front(),test._SYMBptr->feuille._VECTptr->back());
    // FIXME: eval local variables in test that are not in increment and prog
    gen increment=to_increment((*(args._VECTptr))[2]);
    gen prog=(*(args._VECTptr))[3];
    if ( (prog.type==_SYMB) && (prog._SYMBptr->sommet==at_bloc))
      prog=prog._SYMBptr->feuille;
    vecteur forprog=prog.type==_VECT?*prog._VECTptr:vecteur(1,prog);
    iterateur it,itbeg=forprog.begin(),itend=forprog.end();
    for (it=itbeg;it!=itend;++it){
      *it=to_increment(*it);
    }
    gen res,oldres;
    // loop
    int save_current_instruction=debug_ptr(newcontextptr)->current_instruction;
    int eval_lev=eval_level(newcontextptr);
    debug_struct * dbgptr=debug_ptr(newcontextptr);
    try {
      bool findlabel=false;
      gen label;
      for (initialisation.eval(eval_lev,newcontextptr);ck_is_one(test.eval(eval_lev,newcontextptr).evalf(1,newcontextptr));test.val?increment.eval(eval_lev,newcontextptr):0){
	dbgptr->current_instruction=save_current_instruction;
	findlabel=false;
	// add a test for boucle of type program/composite
	// if that's the case call eval with test for break and continue
	for (it=itbeg;it!=itend;++it){
	  oldres=res;
	  ++dbgptr->current_instruction;
	  if (dbgptr->debug_mode)
	    debug_loop(res,newcontextptr);
	  if (!findlabel)
	    res=it->eval(eval_lev,newcontextptr);
	  else
	    res=*it;
	  if (res.type==_SYMB){
	    unary_function_ptr & u=res._SYMBptr->sommet;
	    if (!findlabel){ 
	      if (u==at_return) {
		increment_instruction(it+1,itend,newcontextptr);
		if (bound)
		  leave(protect,loop_var,newcontextptr);
		return res;
	      }
	      if (u==at_break){
		increment_instruction(it+1,itend,newcontextptr);
		test=zero;
		res=oldres;
		break;
	      }
	      if (u==at_continue){
		increment_instruction(it+1,itend,newcontextptr);
		res=oldres;
		break;
	      }
	    }
	    else {
	      if (u==at_label && label==res._SYMBptr->feuille)
		findlabel=false;
	    }
	    if (!findlabel && u==at_goto){
	      findlabel=true;
	      label=res._SYMBptr->feuille;
	    }
	  }
	  if (findlabel && it+1==itend)
	    it=itbeg-1;
	} // end of loop of FOR bloc instructions
      } // end of user FOR loop
      dbgptr->current_instruction=save_current_instruction;
      increment_instruction(itbeg,itend,newcontextptr);
    } // end try
    catch (std::runtime_error & e){
      if (bound)
	leave(protect,loop_var,newcontextptr);
      throw(e);
    }
    if (bound)
      leave(protect,loop_var,newcontextptr);
    return res;
  }

  const string _for_s("for");
  unary_function_eval __for(&_for,_for_s,&printasfor);
  unary_function_ptr at_for (&__for,_QUOTE_ARGUMENTS);

  int bind(const vecteur & vals,const vecteur & vars,context * & contextptr){
#ifdef DEBUG_SUPPORT
    if (vals.size()!=vars.size())
      setsizeerr(gen(vals).print(contextptr)+ " size() != " + gen(vars).print(contextptr));
#endif
    if (debug_ptr(contextptr)->debug_localvars)
      *debug_ptr(contextptr)->debug_localvars=vars;
    const_iterateur it=vals.begin(),itend=vals.end();
    const_iterateur jt=vars.begin();
    gen tmp;
    if (contextptr){
      context * newcontextptr = new context(* contextptr);
      newcontextptr->tabptr = new sym_tab;
      if (contextptr->globalcontextptr)
	newcontextptr->globalcontextptr = contextptr->globalcontextptr;
      else 
	newcontextptr->globalcontextptr = contextptr;
      newcontextptr->previous=contextptr;
      contextptr=newcontextptr;
      if (debug_ptr(contextptr))
	debug_ptr(contextptr)->debug_contextptr=contextptr;
    }
    for (;it!=itend;++it,++jt){
      if ( jt->type==_SYMB && jt->_SYMBptr->sommet==at_check_type){
	tmp=jt->_SYMBptr->feuille._VECTptr->back();
	_check_type(makevecteur(jt->_SYMBptr->feuille._VECTptr->front(),*it),contextptr);
      }
      else {
	if ( jt->type==_SYMB &&  jt->_SYMBptr->sommet==at_double_deux_points ){
	  tmp=jt->_SYMBptr->feuille._VECTptr->front();
	  _check_type(makevecteur(jt->_SYMBptr->feuille._VECTptr->back(),*it),contextptr);
	}
	else
	  tmp=*jt;
      }
      if (tmp.type==_IDNT){
	if (contextptr)
	  (*contextptr->tabptr)[*tmp._IDNTptr->name]=*it;
	else
	  tmp._IDNTptr->push(protection_level,*it);
      }
      else {
	if (tmp.type==_FUNC)
	  setsizeerr("Reserved word:"+tmp.print(contextptr));
	else
	  setsizeerr("Not bindable"+tmp.print(contextptr));
      }
    }
    if (!contextptr)
      ++protection_level;
    return protection_level-1;
  }

  bool leave(int protect,vecteur & vars,context * & contextptr){
    iterateur it=vars.begin(),itend=vars.end(),jt,jtend;
    gen tmp;
    if (contextptr){
      if (contextptr->previous){
	context * tmpptr=contextptr;
	contextptr=contextptr->previous;
	if (debug_ptr(contextptr))
	  debug_ptr(contextptr)->debug_contextptr=contextptr;
	if (tmpptr->tabptr){
	  delete tmpptr->tabptr;
	  delete tmpptr;
	  return true;
	}
      }
      return false;
    }
    for (;it!=itend;++it){
      if (it->type==_SYMB && it->_SYMBptr->sommet==at_check_type)
	tmp=it->_SYMBptr->feuille._VECTptr->back();
      else {
	if (it->type==_SYMB && it->_SYMBptr->sommet==at_double_deux_points)
	  tmp=it->_SYMBptr->feuille._VECTptr->front();
	else
	  tmp=*it;
      }
#ifdef DEBUG_SUPPORT
      if (tmp.type!=_IDNT) setsizeerr("prog.cc/leave");
#endif    
      jt=tmp._IDNTptr->localvalue->begin(),jtend=tmp._IDNTptr->localvalue->end();
      for (;;){
	if (jt==jtend)
	  break;
	--jtend;
	--jtend;
	if (protect>jtend->val){
	  ++jtend;
	  ++jtend;
	  break;
	}
      }
      tmp._IDNTptr->localvalue->erase(jtend,tmp._IDNTptr->localvalue->end());
    }
    protection_level=protect;
    return true;
  }

  string printaslocal(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if ( (feuille.type!=_VECT) || (feuille._VECTptr->size()!=2) )
      return sommetstr+'('+feuille.print(contextptr)+')';
    const_iterateur it=feuille._VECTptr->begin(),itend=feuille._VECTptr->end();
    string res;
    if (xcas_mode(contextptr)>0){
      if (xcas_mode(contextptr)==3)
	res=indent(contextptr)+"Local ";
      else
	res=indent(contextptr)+"local ";
    }
    else
      res=indent(contextptr)+"{ local ";
    if ( (it->type==_VECT) && (it->_VECTptr->size()==1) )
      res += it->_VECTptr->front().print(contextptr);
    else
      res += gen(*it->_VECTptr,_SEQ__VECT).print(contextptr);
    if (xcas_mode(contextptr)!=3)
      res +=';';
    if (xcas_mode(contextptr)>0 && xcas_mode(contextptr)!=3)
      res += indent(contextptr)+"begin ";
    debug_ptr(contextptr)->indent_spaces +=2;
    ++it;
    for ( ;;){
      if (it->type!=_VECT)
	res += indent(contextptr)+it->print(contextptr);
      else {
	const_iterateur jt=it->_VECTptr->begin(),jtend=it->_VECTptr->end();
	for (;jt!=jtend;++jt){
	  res += indent(contextptr)+jt->print(contextptr);
	  if (xcas_mode(contextptr)!=3)
	    res += "; " ;
	}
      }
      ++it;
      if (it==itend){
	debug_ptr(contextptr)->indent_spaces -= 2;
	switch (xcas_mode(contextptr)){
	case 0:
	  res += indent(contextptr)+"}";
	  break;
	case 1: case 1+_DECALAGE:
	  res+=indent(contextptr)+"end;";
	  break;
	case 2:
	  return res+=indent(contextptr)+"end_proc;";
	  break;
	}
	return res;
      }
      else
	if (xcas_mode(contextptr)!=3)
	  res +="; ";
    }
  }
  gen symb_local(const gen & args,GIAC_CONTEXT){
    if (args.type==_VECT && args._VECTptr->size()==2)
      return symb_local(args._VECTptr->front(),args._VECTptr->back(),contextptr);
    return symbolic(at_local,args);
  }

  gen symb_local(const gen & a,const gen & b,GIAC_CONTEXT){
    gen newa,newb;
    replace_keywords(a,b,newa,newb,contextptr);
    return symbolic(at_local,makevecteur(newa,newb));
  }
  gen _local(const gen & args,const context * contextptr) {
    if (args.type!=_VECT)
      return symb_local(args,contextptr);
    int s=args._VECTptr->size();
    if (s!=2)
      setsizeerr("Local must have 2 args");
    // Initialization
    gen vars=args._VECTptr->front();
    if (vars.type!=_VECT)
      vars=makevecteur(vars);
    vecteur names,values;
    iterateur it=vars._VECTptr->begin(),itend=vars._VECTptr->end();
    names.reserve(itend-it);
    values.reserve(itend-it);
    for (;it!=itend;++it){
      if (it->type==_IDNT){
	names.push_back(*it);
	values.push_back(0);
	continue;
      }
      if ( (it->type!=_SYMB) || (it->_SYMBptr->sommet!=at_sto))
	settypeerr();
      gen nom=it->_SYMBptr->feuille._VECTptr->back();
      gen val=it->_SYMBptr->feuille._VECTptr->front().eval(eval_level(contextptr),contextptr);
      if (nom.type!=_IDNT)
	settypeerr();
      names.push_back(nom);
      values.push_back(val);
    }
    context * newcontextptr = (context *) contextptr;
    int protect=bind(values,names,newcontextptr);
    gen prog=args._VECTptr->back(),res,newres;
    if (prog.type!=_VECT){
      ++debug_ptr(newcontextptr)->current_instruction;
      if (debug_ptr(newcontextptr)->debug_mode)
	debug_loop(res,newcontextptr);
      res=prog.eval(eval_level(newcontextptr),newcontextptr);
    }
    else {
      it=prog._VECTptr->begin(),itend=prog._VECTptr->end();
      bool findlabel=false;
      gen label;
      for (;it!=itend;++it){
	++debug_ptr(newcontextptr)->current_instruction;
	// cout << *it << endl;
	if (debug_ptr(newcontextptr)->debug_mode)
	  debug_loop(res,newcontextptr);
	if (!findlabel)
	  res=it->eval(eval_level(newcontextptr),newcontextptr);
	else
	  res=*it;
	if (res.type==_SYMB){
	  unary_function_ptr & u=res._SYMBptr->sommet;
	  if (findlabel && u==at_label && label==res._SYMBptr->feuille)
	    findlabel=false;
	  if (!findlabel && u==at_goto){
	    findlabel=true;
	    label=res._SYMBptr->feuille;
	  }
	}
	if (findlabel && it+1==itend)
	  it=prog._VECTptr->begin()-1;
	if (!findlabel && is_return(res,newres) ){
	  // res=newres;
	  break;
	}
      }
    }
    leave(protect,names,newcontextptr);
    return res;
  }

  const string _local_s("local");
  unary_function_eval __local(&_local,_local_s,&printaslocal);
  unary_function_ptr at_local (&__local,_QUOTE_ARGUMENTS);

  string printasreturn(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if ( (xcas_mode(contextptr)==1) || (xcas_mode(contextptr)==1+_DECALAGE) )
      return "RETURN("+feuille.print(contextptr)+")";
    if (xcas_mode(contextptr)==3)
      return "Return "+feuille.print(contextptr);
    return sommetstr+"("+feuille.print(contextptr)+")";
  }
  gen symb_return(const gen & args){
    return symbolic(at_return,args);
  }
  const string _return_s("return");
  unary_function_unary __return(&symb_return,_return_s,&printasreturn);
  unary_function_ptr at_return (&__return);

  string printastry_catch(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if ( (feuille.type!=_VECT) || (feuille._VECTptr->size()!=3) )
      return sommetstr+'('+feuille.print(contextptr)+')';
    const_iterateur it=feuille._VECTptr->begin();//,itend=feuille._VECTptr->end();
    string res;
    if (xcas_mode(contextptr)==3)
      res += "Try";
    else
      res += "try ";
    res += it->print(contextptr);
    ++it;
    if (xcas_mode(contextptr)==3){
      res += indent(contextptr)+"Else";
      ++it;
      if (!is_undef(*it))
	res += printasinnerbloc(*it,contextptr);
      res += indent(contextptr)+"EndTry";
    }
    else {
      if (res[res.size()-1]!='}')
	res += "; ";
      res += "catch(" + it->print(contextptr) + ")";
      ++it;
      res += it->print(contextptr);
      if (res[res.size()-1]!='}')
	res += "; ";
    }
    return res;
  }
  
  gen symb_try_catch(const gen & args){
    return symbolic(at_try_catch,args);
  }
  gen _try_catch(const gen & args,const context * contextptr){
    if (args.type!=_VECT)
      return symb_try_catch(args);
    if (args._VECTptr->size()!=3)
      setsizeerr("Try_catch must have 3 args");
    gen res;
    int saveprotect=protection_level;
    vector< vector<int> > save_sst_at_stack(debug_ptr(contextptr)->sst_at_stack);
    vecteur save_args_stack(debug_ptr(contextptr)->args_stack);
    vector<int> save_current_instruction_stack=debug_ptr(contextptr)->current_instruction_stack;
    int save_current_instruction=debug_ptr(contextptr)->current_instruction;
    try {
      ++debug_ptr(contextptr)->current_instruction;
      if (debug_ptr(contextptr)->debug_mode)
	debug_loop(res,contextptr);
      res=args._VECTptr->front().eval(eval_level(contextptr),contextptr);
    }
    catch (std::runtime_error & error ){
      if (!contextptr)
	protection_level=saveprotect;
      debug_ptr(contextptr)->sst_at_stack=save_sst_at_stack;
      debug_ptr(contextptr)->args_stack=save_args_stack;
      debug_ptr(contextptr)->current_instruction_stack=save_current_instruction_stack;
      gen id=(*(args._VECTptr))[1];
      string er(error.what());
      er = '"'+er+'"';
      if (id.type==_IDNT)
	sto(gen(er,contextptr),id,contextptr);
      debug_ptr(contextptr)->current_instruction=save_current_instruction;
      increment_instruction(args._VECTptr->front(),contextptr);
      ++debug_ptr(contextptr)->current_instruction;
      if (debug_ptr(contextptr)->debug_mode)
	debug_loop(res,contextptr);
      res=args._VECTptr->back().eval(eval_level(contextptr),contextptr);
    }
    debug_ptr(contextptr)->current_instruction=save_current_instruction;
    increment_instruction(args._VECTptr->front(),contextptr);
    increment_instruction(args._VECTptr->back(),contextptr);
    return res;
  }
  const string _try_catch_s("try_catch");
  unary_function_eval __try_catch(&_try_catch,_try_catch_s,&printastry_catch);
  unary_function_ptr at_try_catch (&__try_catch,_QUOTE_ARGUMENTS);

  gen feuille_(const gen & g,const gen & interval,GIAC_CONTEXT){
    vecteur v;
    if (g.type==_SYMB){
      gen & f=g._SYMBptr->feuille;
      if (f.type==_VECT)
	v=*f._VECTptr;
      else
	v=vecteur(1,f);
    }
    else {
      if (g.type==_VECT)
	v=*g._VECTptr;
      else
	v=vecteur(1,g);
    }
    int s=v.size();
    if (interval.type==_INT_){
      int i=interval.val-(xcas_mode(contextptr)!=0);
      if (i==-1 && g.type==_SYMB)
	return g._SYMBptr->sommet;
      if (i<0 || i>=s)
	setdimerr();
      return v[i];
    }
    if (interval.is_symb_of_sommet(at_interval)&& interval._SYMBptr->feuille.type==_VECT){
      vecteur & w=*interval._SYMBptr->feuille._VECTptr;
      if (w.size()!=2 || w.front().type!=_INT_ || w.back().type!=_INT_)
	settypeerr();
      int i=w.front().val,j=w.back().val;
      if (i>j)
	return gen(vecteur(0),_SEQ__VECT);
      if (xcas_mode(contextptr)){
	--i;
	--j;
      }
      if (i<0 || i>=s || j<0 || j>=s)
	setdimerr();
      return gen(vecteur(v.begin()+i,v.begin()+j+1),_SEQ__VECT);
    }
    setsizeerr();
    return 0;
  }
  gen symb_feuille(const gen & args){
    return symbolic(at_feuille,args);
  }
  gen _feuille(const gen & args,GIAC_CONTEXT){
    if (args.type==_VECT){
      if (args.subtype==_SEQ__VECT && args._VECTptr->size()==2)
	return feuille_(args._VECTptr->front(),args._VECTptr->back(),contextptr);
      return gen(*args._VECTptr,_SEQ__VECT);
    }
    if (args.type!=_SYMB)
      return symb_feuille(args);
    gen tmp=args._SYMBptr->feuille;
    if (tmp.type==_VECT)
      tmp.subtype=_SEQ__VECT;
    return tmp;
  }
  const string _feuille_s("op");
  unary_function_eval __feuille(&_feuille,_feuille_s,&printassubs);
  unary_function_ptr at_feuille (&__feuille,0,true);
  
  gen _maple_op(const gen & args,GIAC_CONTEXT){
    if (args.type==_VECT){
      vecteur & v=*args._VECTptr;
      if (args.subtype==_SEQ__VECT && v.size()>1)
	return feuille_(v.back(),v.front(),contextptr);
      return gen(v,_SEQ__VECT);
    }
    if (args.type!=_SYMB)
      return args; // was symbolic(at_maple_op,args);
    return args._SYMBptr->feuille;
  }
  const string _maple_op_s("op");
  unary_function_eval __maple_op(&_maple_op,_maple_op_s,&printasmaple_subs);
  unary_function_ptr at_maple_op (&__maple_op);
  
  gen symb_sommet(const gen & args){
    return symbolic(at_sommet,args);
  }
  gen _sommet(const gen & args){
    if (args.type!=_SYMB)
      return symb_sommet(args);
    int nargs;
    if (args._SYMBptr->feuille.type==_VECT)
        nargs=args._SYMBptr->feuille._VECTptr->size();
    else
        nargs=1;
    return gen(args._SYMBptr->sommet,nargs);
  }
  const string _sommet_s("sommet");
  unary_function_unary __sommet(&_sommet,_sommet_s);
  unary_function_ptr at_sommet (&__sommet,0,true);

  // replace in g using equalities in v
  gen subsop(const vecteur & g,const vecteur & v,const gen & sommet,GIAC_CONTEXT){
    gen newsommet=sommet;
    vecteur res(g);
    const_iterateur it=v.begin(),itend=v.end();
    for (;it!=itend;++it){
      if ( (!it->is_symb_of_sommet(at_equal) && !it->is_symb_of_sommet(at_same)) || it->_SYMBptr->feuille.type!=_VECT || it->_SYMBptr->feuille._VECTptr->size()!=2){
	*logptr(contextptr) << "Unknown subsop rule " << *it << endl;
	continue;
      }
      vecteur w=*it->_SYMBptr->feuille._VECTptr;
      if (w.front().type==_VECT){
	vecteur rec=*w.front()._VECTptr;
	if (rec.size()<1)
	  setdimerr();
	if (rec.size()==1)
	  w.front()=rec.front();
	else {
	  int i=rec.front().val;
	  if (xcas_mode(contextptr))
	    --i;
	  if (rec.front().type!=_INT_ || i<0 || i>=signed(res.size()))
	    setdimerr();
	  res[i]=subsop(res[i],vecteur(1,symbolic(at_equal,makevecteur(vecteur(rec.begin()+1,rec.end()),w.back()))),contextptr);
	  continue;
	}
      }
      if (w.front().type!=_INT_)
	continue;
      int i=w.front().val;
      if (xcas_mode(contextptr))
	--i;
      if (i==-1){
	newsommet=w.back();
	continue;
      }
      if (i<0 || i>=signed(res.size()))
	setdimerr();
      res[i]=w.back();
    }
    it=res.begin();
    itend=res.end();
    vecteur res1;
    res1.reserve(itend-it);
    for (;it!=itend;++it){
      if (it->type!=_VECT || it->subtype!=_SEQ__VECT || !it->_VECTptr->empty() )
	res1.push_back(*it);
    }
    if (newsommet.type!=_FUNC)
      return res1;
    else
      return symbolic(*newsommet._FUNCptr,res1);
  }
  gen subsop(const gen & g,const vecteur & v,GIAC_CONTEXT){
    if (g.type==_VECT)
      return subsop(*g._VECTptr,v,0,contextptr);
    if (g.type!=_SYMB)
      return g;
    vecteur w(gen2vecteur(g._SYMBptr->feuille));
    return subsop(w,v,g._SYMBptr->sommet,contextptr);
  }
  gen _maple_subsop(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      setsizeerr();
    vecteur & v=*args._VECTptr;
    int s=v.size();
    if (s<2)
      setdimerr();
    return subsop(v.back(),vecteur(v.begin(),v.end()-1),contextptr);
  }
  const string _maple_subsop_s("subsop");
  unary_function_eval __maple_subsop(&_maple_subsop,_maple_subsop_s,&printasmaple_subs);
  unary_function_ptr at_maple_subsop (&__maple_subsop);
  
  gen _subsop(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      setsizeerr();
    vecteur & v=*args._VECTptr;
    int s=v.size();
    if (s<2)
      setdimerr();
    return subsop(v.front(),vecteur(v.begin()+1,v.end()),contextptr);
  }
  const string _subsop_s("subsop");
  unary_function_eval __subsop(&_subsop,_subsop_s,&printassubs);
  unary_function_ptr at_subsop (&__subsop);
  
  gen symb_append(const gen & args){
    return symbolic(at_append,args);
  }
  gen _append(const gen & args){
    if ( args.type!=_VECT || !args._VECTptr->size() )
      setsizeerr();
    const_iterateur it=args._VECTptr->begin(),itend=args._VECTptr->end();
    if (itend-it==2 && it->type==_STRNG && (it+1)->type==_STRNG)
      return string2gen(*it->_STRNGptr+*(it+1)->_STRNGptr,false);
    if (it->type!=_VECT)
      setsizeerr();
    vecteur v(*it->_VECTptr);
    int subtype=it->subtype;
    ++it;
    for (;it!=itend;++it)
      v.push_back(*it);
    return gen(v,subtype);
  }
  const string _append_s("append");
  unary_function_unary __append(&_append,_append_s);
  unary_function_ptr at_append (&__append,0,true);

  gen symb_prepend(const gen & args){
    return symbolic(at_prepend,args);
  }
  gen _prepend(const gen & args){
    if ( (args.type!=_VECT) || (!args._VECTptr->size()) || (args._VECTptr->front().type!=_VECT) )
      setsizeerr();
    gen debut=args._VECTptr->front();
    return gen(mergevecteur(cdr_VECT(*args._VECTptr),*debut._VECTptr),debut.subtype);
  }
  const string _prepend_s("prepend");
  unary_function_unary __prepend(&_prepend,_prepend_s);
  unary_function_ptr at_prepend (&__prepend,0,true);

  gen symb_contains(const gen & args){
    return symbolic(at_contains,args);
  }
  gen _contains(const gen & args){
    if ( (args.type!=_VECT) || (args._VECTptr->size()!=2) || (args._VECTptr->front().type!=_VECT) )
      return symb_contains(args);
    return equalposcomp(*args._VECTptr->front()._VECTptr,args._VECTptr->back());
  }
  const string _contains_s("contains");
  unary_function_unary __contains(&_contains,_contains_s);
  unary_function_ptr at_contains (&__contains,0,true);

  gen select_remove(const gen & args,bool selecting,const context * contextptr){
    if ( (args.type!=_VECT) || (args._VECTptr->size()<2)){
      if (selecting)
	return symb_select(args);
      else
	return symb_remove(args);
    }
    gen v((*(args._VECTptr))[1]);
    int subtype;
    unary_function_ptr * fn=0;
    if (v.type==_SYMB){
      if (v._SYMBptr->feuille.type==_VECT)
	v=v._SYMBptr->feuille;
      else
	v=makevecteur(v._SYMBptr->feuille);
      subtype=-1;
      fn=&v._SYMBptr->sommet;
    }
    else
      subtype=v.subtype;
    if ( (v.type!=_VECT) && (v.type!=_SYMB)){
      if (selecting)
	return symb_select(args);
      else
	return symb_remove(args);
    }
    gen f(args._VECTptr->front());
    vecteur otherargs(args._VECTptr->begin()+1,args._VECTptr->end());
    const_iterateur it=v._VECTptr->begin(),itend=v._VECTptr->end();
    vecteur res;
    res.reserve(itend-it);
    if (otherargs.size()==1){
      for (;it!=itend;++it){
	if (is_zero(f(*it,contextptr))!=selecting)
	  res.push_back(*it);
      }
    }
    else {
      for (;it!=itend;++it){
	otherargs.front()=*it;
	if (is_zero(f(otherargs,contextptr))!=selecting)
	  res.push_back(*it);
      }
    }
    if (subtype<0)
      return symbolic(*fn,res);
    else
      return gen(res,subtype);
  }
  gen symb_select(const gen & args){
    return symbolic(at_select,args);
  }
  gen _select(const gen & args,const context * contextptr){
    return select_remove(args,1,contextptr);
  }
  const string _select_s("select");
  unary_function_eval __select(&_select,_select_s);
  unary_function_ptr at_select (&__select,0,true);

  gen symb_remove(const gen & args){
    return symbolic(at_remove,args);
  }
  gen _remove(const gen & args,const context * contextptr){
    return select_remove(args,0,contextptr);
  }
  const string _remove_s("remove");
  unary_function_eval __remove(&_remove,_remove_s);
  unary_function_ptr at_remove (&__remove,0,true);

  string printnostring(const gen & g,GIAC_CONTEXT){
    if (g.type==_STRNG)
      return *g._STRNGptr;
    else
      return g.print(contextptr);
  }
  gen concat(const gen & g,bool glue_lines,GIAC_CONTEXT){
    if (g.type!=_VECT)
      return symb_concat(g);
    vecteur & v=*g._VECTptr;
    if (v.size()!=2){
      if (g.subtype==_SEQ__VECT)
	return g;
      return symb_concat(g);
    }
    gen v0=v[0],v1=v[1];
    if (v0.type==_VECT && v1.type==_VECT){
      if (!glue_lines && v1.subtype!=_SEQ__VECT && ckmatrix(v0) && ckmatrix(v1) && v0._VECTptr->size()==v1._VECTptr->size() )
	return gen(mtran(mergevecteur(mtran(*v0._VECTptr),mtran(*v1._VECTptr))));
      else
	return gen(mergevecteur(*v0._VECTptr,*v1._VECTptr),v0.subtype);
    }
    if (v0.type==_VECT)
      return gen(mergevecteur(*v0._VECTptr,vecteur(1,v1)),v0.subtype);
    if (v1.type==_VECT)
      return gen(mergevecteur(vecteur(1,v0),*v1._VECTptr),v1.subtype);
    if ( (v0.type==_STRNG) || (v1.type==_STRNG) )
      return string2gen(printnostring(v0,contextptr) + printnostring(v1,contextptr),false);
    return 0;
  }
  gen symb_concat(const gen & args){
    return symbolic(at_concat,args);
  }
  gen _concat(const gen & args,GIAC_CONTEXT){
    return concat(args,false,contextptr);
  }
  const string _concat_s("concat");
  unary_function_eval __concat(&_concat,_concat_s);
  unary_function_ptr at_concat (&__concat,0,true);

  gen symb_option(const gen & args){
    return symbolic(at_option,args);
  }
  gen _option(const gen & args){
    return symb_option(args);
  }
  const string _option_s("option");
  unary_function_unary __option(&_option,_option_s);
  unary_function_ptr at_option (&__option);

  string printascase(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if ( (feuille.type!=_VECT) || (feuille._VECTptr->size()!=2) || (feuille._VECTptr->back().type!=_VECT))
      return sommetstr+'('+feuille.print(contextptr)+')';
    string res("switch (");
    res += feuille._VECTptr->front().print(contextptr);
    res += "){";
    debug_ptr(contextptr)->indent_spaces +=2;
    const_iterateur it=feuille._VECTptr->back()._VECTptr->begin(),itend=feuille._VECTptr->back()._VECTptr->end();
    for (;it!=itend;++it){
      ++it;
      if (it==itend){
	res += indent(contextptr)+"default:";
	--it;
	debug_ptr(contextptr)->indent_spaces += 2;
	res += indent(contextptr)+it->print(contextptr);
	debug_ptr(contextptr)->indent_spaces -= 2;
	break;
      }
      res += indent(contextptr)+"case "+(it-1)->print(contextptr)+":";
      debug_ptr(contextptr)->indent_spaces += 2;
      res += indent(contextptr)+it->print(contextptr);
      debug_ptr(contextptr)->indent_spaces -=2;
    }
    debug_ptr(contextptr)->indent_spaces -=2;
    res+=indent(contextptr)+"}";
    return res;
  }
  gen symb_case(const gen & args){
    return symbolic(at_case,args);
  }
  gen symb_case(const gen & a,const gen & b){
    return symbolic(at_case,makevecteur(a,b));
  }
  gen _case(const gen & args,GIAC_CONTEXT){ // FIXME DEBUGGER
    if ( (args.type!=_VECT) || (args._VECTptr->size()!=2) || (args._VECTptr->back().type!=_VECT) )
      return symb_case(args);
    gen expr=args._VECTptr->front().eval(eval_level(contextptr),contextptr),res=undef,oldres,newres;
    const_iterateur it=args._VECTptr->back()._VECTptr->begin(),itend=args._VECTptr->back()._VECTptr->end();
    for (;it!=itend;){
      if (it+1==itend){
	res=it->eval(eval_level(contextptr),contextptr);
	break;
      }
      if (expr==it->eval(eval_level(contextptr),contextptr)){
	++it;
	oldres=res;
	res=it->eval(eval_level(contextptr),contextptr);
	if (res==symbolic(at_break,zero)){
	  res=oldres;
	  break;
	}
	if (res.is_symb_of_sommet(at_return))
	  break;
      }
      else
	++it;
      if (it!=itend)
	++it;
    }
    return res;
  }
  const string _case_s("case");
  unary_function_eval __case(&_case,_case_s,&printascase);
  unary_function_ptr at_case (&__case,_QUOTE_ARGUMENTS);

  gen symb_rand(const gen & args){
    return symbolic(at_rand,args);
  }
  gen rand_interval(const vecteur & v,GIAC_CONTEXT){
    static gen randmax_plus_one=gen(RAND_MAX)+1;
    gen x1=v.front(),x2=v.back();
    if (x1==x2)
      return x1;
    if (xcas_mode(contextptr)==1 && is_integer(x1) && is_integer(x2) ){
      if (is_positive(x1-x2,contextptr))
	std::swap<gen>(x1,x2);
      int n=(x2-x1).bindigits()/gen(RAND_MAX).bindigits()+1;
      // Make n random numbers
      gen res=zero;
      for (int i=0;i<n;++i)
	res=(gen(RAND_MAX)+1)*res+giac_rand(contextptr);
      // Now res is in [0,(RAND_MAX+1)^n-1]
      // Rescale in x1..x2
      return x1+_iquo(makevecteur(res*(x2-x1),pow(gen(RAND_MAX)+1,n)));
    }
#ifdef HAVE_LIBMPFR
    if (x1.type==_REAL && x2.type==_REAL){
      int n=mpfr_get_prec(x1._REALptr->inf);
      int nr=int(n*std::log(2.0)/std::log(RAND_MAX+1.0));
      gen xr=0;
      for (int i=0;i<=nr;++i){
	xr=xr*randmax_plus_one+giac_rand(contextptr);
      }
      return x1+((x2-x1)*xr)/pow(randmax_plus_one,nr+1);
    }
#endif
    gen x=evalf_double(x1,1,contextptr),y=evalf_double(x2,1,contextptr);
    if ( (x.type==_DOUBLE_) && (y.type==_DOUBLE_) ){
      double xd=x._DOUBLE_val,yd=y._DOUBLE_val;
      double xr= (giac_rand(contextptr)/(RAND_MAX+1.0))*(yd-xd)+xd;
      return xr;
    }
    return symb_rand(gen(v,_SEQ__VECT));
  }

  gen rand_n_in_list(int n,const vecteur & v,GIAC_CONTEXT){
    n=absint(n);
    if (signed(v.size())<n)
      setdimerr();
    vecteur w(v);
    vecteur res;
    for (int i=0;i<n;++i){
      int tmp=int((double(giac_rand(contextptr))*w.size())/RAND_MAX);
      res.push_back(w[tmp]);
      w.erase(w.begin()+tmp);
    }
    return res;
  }
  gen _rand(const gen & args,GIAC_CONTEXT){
    if (args.type==_INT_){
      if (args.val<0)
	return -(xcas_mode(contextptr)==3)+int(args.val*(giac_rand(contextptr)/(RAND_MAX+1.0)));
      else
	return (xcas_mode(contextptr)==3)+int(args.val*(giac_rand(contextptr)/(RAND_MAX+1.0)));
    }
    if (args.type==_ZINT)
      return rand_interval(makevecteur(zero,args),contextptr);
    if (args.type==_VECT){ 
      if (args._VECTptr->empty())
	return giac_rand(contextptr);
      vecteur & v=*args._VECTptr;
      int s=v.size();
      if (s==2){ 
	if (v.front().type==_INT_ && v.back().type==_VECT){ // rand(n,list) choose n in list
	  return rand_n_in_list(v.front().val,*v.back()._VECTptr,contextptr);
	}
	if ( (v.back().type==_SYMB) && (v.back()._SYMBptr->sommet==at_interval) ){
	  // arg1=loi, arg2=intervalle
	}
	return rand_interval(v,contextptr);
      }
      if (s==3 && v[0].type==_INT_ && v[1].type==_INT_ && v[2].type==_INT_){ 
	// 3 integers expected, rand(n,min,max) choose n in min..max
	int n=v[0].val;
	int m=v[1].val;
	int M=v[2].val;
	if (m>M){ int tmp=m; m=M; M=tmp; }
	vecteur v;
	for (int i=m;i<=M;++i) v.push_back(i);
	return rand_n_in_list(n,v,contextptr);
      }
    }
    if ( (args.type==_SYMB) && (args._SYMBptr->sommet==at_interval) ){
      vecteur & v=*args._SYMBptr->feuille._VECTptr;
      return symb_program(vecteur(0),vecteur(0),symb_rand(gen(v,_SEQ__VECT)),contextptr);
      // return rand_interval(v);
    }
    return symb_rand(args);
  }
  const string _rand_s("rand");
  unary_function_eval __rand(&_rand,_rand_s);
  unary_function_ptr at_rand (&__rand,0,true);

  gen _srand(const gen & args,GIAC_CONTEXT){
    if (args.type==_INT_){
      srand(args.val);
      rand_seed(args.val,contextptr);
      return args;
    }
    else {
      int t=time(NULL);
      rand_seed(t,contextptr);
      srand(t);
      return t;
    }
  }
  const string _srand_s("srand");
  unary_function_eval __srand(&_srand,_srand_s);
  unary_function_ptr at_srand (&__srand);

  gen symb_char(const gen & args){
    return symbolic(at_char,args);
  }
  gen _char(const gen & args){
    string s;
    if (args.type==_INT_){
      s += args.val ;
    }
    else {
      if (args.type==_VECT){
	const_iterateur it=args._VECTptr->begin(),itend=args._VECTptr->end();
	for (;it!=itend;++it){
	  s += it->val;
	}
      }
      else return symb_char(args);
    }
    gen tmp;
    tmp.type=_STRNG;
    tmp.ptr_val.ref_count = new int(1);
    tmp._STRNGptr = new string(s);
    return tmp;
  }
  const string _char_s("char");
  unary_function_unary __char(&_char,_char_s);
  unary_function_ptr at_char (&__char,0,true);

  gen symb_asc(const gen & args){
    return symbolic(at_asc,args);
  }
  gen _asc(const gen & args){
    if (args.type==_STRNG){
      int l=args._STRNGptr->size();
      vecteur v(l);
      for (int i=0;i<l;++i)
	v[i]=int( (unsigned char) ((*args._STRNGptr)[i]));
      return v;
    }
    if (args.type==_VECT){
      if ( (args._VECTptr->size()!=2) ||(args._VECTptr->front().type!=_STRNG) || (args._VECTptr->back().type!=_INT_) )
	setsizeerr("asc");
      return int( (unsigned char) (*args._VECTptr->front()._STRNGptr)[args._VECTptr->back().val]);
    }
    else return symb_asc(args);
  }

  const string _asc_s("asc");
  unary_function_unary __asc(&_asc,_asc_s);
  unary_function_ptr at_asc (&__asc,0,true);

  gen symb_map(const gen & args){
    return symbolic(at_map,args);
  }
  gen _map(const gen & args,const context * contextptr){
    if (args.type!=_VECT)
      return symb_map(args);
    vecteur v=*args._VECTptr;
    int s=v.size();
    if (s<2)
      toofewargs("");
    gen objet=v.front();
    gen to_map=v[1];
    // FIXME: should have maple_map and mupad_map functions
    if (xcas_mode(contextptr)==1){
      objet=v[1];
      to_map=v.front();
    }
    bool matrix = ckmatrix(objet) && s>2;
    if (matrix){
      matrix=false;
      for (int i=2;i<s;++i){
	if (v[i]==at_matrix){
	  v.erase(v.begin()+i);
	  --s;
	  matrix=true;
	  break;
	}
      }
    }
    if (to_map.type==_VECT)
      setsizeerr();
    if (v.size()==2){
      if (objet.type==_SYMB){
	gen & f=objet._SYMBptr->feuille;
	gen tmp=_map(makevecteur(f,to_map),contextptr);
	if (f.type==_VECT && tmp.type==_VECT)
	  tmp.subtype=f.subtype;
	if (objet._SYMBptr->sommet==at_equal || objet._SYMBptr->sommet==at_same)
	  return symbolic(at_equal,tmp);
	return objet._SYMBptr->sommet(tmp,contextptr);
      }
      // if (to_map.type==_FUNC) return apply(objet,*to_map._FUNCptr);
      if (objet.type==_POLY){
	int dim=objet._POLYptr->dim;
	polynome * resptr=new polynome(dim);
	vector< monomial<gen> >::const_iterator it=objet._POLYptr->coord.begin(),itend=objet._POLYptr->coord.end();
	resptr->coord.reserve(itend-it);
	vecteur argv(dim+1);
	for (;it!=itend;++it){
	  argv[0]=it->value;
	  index_t * i=it->index.iptr;
	  for (int j=0;j<dim;j++)
	    argv[j+1]=(*i)[j];
	  gen g=to_map(gen(argv,_SEQ__VECT),contextptr);
	  if (!is_zero(g))
	    resptr->coord.push_back(monomial<gen>(g,it->index));
	}
	return resptr;
      }
      if (objet.type!=_VECT)
	return to_map(objet,contextptr);
      const_iterateur it=objet._VECTptr->begin(),itend=objet._VECTptr->end();
      vecteur res;
      res.reserve(itend-it);
      for (;it!=itend;++it){
	if (matrix && it->type==_VECT){
	  const vecteur & tmp = *it->_VECTptr;
	  const_iterateur jt=tmp.begin(),jtend=tmp.end();
	  vecteur tmpres;
	  tmpres.reserve(jtend-jt);
	  for (;jt!=jtend;++jt){
	    tmpres.push_back(to_map(*jt,contextptr));
	  }
	  res.push_back(tmpres);
	}
	else
	  res.push_back(to_map(*it,contextptr));
      }
      return res;
    }
    if (objet.type==_POLY){
      int dim=objet._POLYptr->dim;
      vecteur opt(v.begin()+2,v.end());
      opt=mergevecteur(vecteur(dim+1),opt);
      polynome * resptr=new polynome(dim);
      vector< monomial<gen> >::const_iterator it=objet._POLYptr->coord.begin(),itend=objet._POLYptr->coord.end();
      resptr->coord.reserve(itend-it);
      for (;it!=itend;++it){
	opt[0]=it->value;
	index_t * i=it->index.iptr;
	for (int j=0;j<dim;j++)
	  opt[j+1]=(*i)[j];
	gen g=to_map(gen(opt,_SEQ__VECT),contextptr);
	if (!is_zero(g))
	  resptr->coord.push_back(monomial<gen>(g,it->index));
      }
      return resptr;
    }
    vecteur opt(v.begin()+1,v.end());
    opt[0]=objet;
    if (objet.type!=_VECT)
      return to_map(opt,contextptr);
    const_iterateur it=objet._VECTptr->begin(),itend=objet._VECTptr->end();
    vecteur res;
    res.reserve(itend-it);
    for (;it!=itend;++it){
      if (matrix && it->type==_VECT){
	const vecteur & tmp = *it->_VECTptr;
	const_iterateur jt=tmp.begin(),jtend=tmp.end();
	vecteur tmpres;
	tmpres.reserve(jtend-jt);
	for (;jt!=jtend;++jt){
	  opt[0]=*jt;
	  tmpres.push_back(to_map(gen(opt,_SEQ__VECT),contextptr));
	}
	res.push_back(tmpres);
      }
      else {
	opt[0]=*it;
	res.push_back(to_map(gen(opt,_SEQ__VECT),contextptr));
      }
    }
    return res;
  }
  const string _map_s("map");
  unary_function_eval __map(&_map,_map_s);
  unary_function_ptr at_map (&__map,0,true);
  
  gen symb_apply(const gen & args){
    return symbolic(at_apply,args);
  }
  gen _apply(const gen & args,const context * contextptr){
    if (args.type!=_VECT)
      return symb_apply(args);
    if (args._VECTptr->empty())
      setsizeerr("apply");
    vecteur & v=*args._VECTptr;
    gen to_apply=v.front();
    int n=to_apply.subtype;
    int n2=v.size();
    if (to_apply.type!=_FUNC)
      n=n2-1;
    if (n && (n2==n+1) ){
      vecteur res;
      for (int i=0;;++i){
	vecteur tmp;
	bool finished=true;
	for (int j=1;j<=n;++j){
	  gen & g=v[j];
	  if (g.type!=_VECT)
	    tmp.push_back(g);
	  else {
	    if (signed(g._VECTptr->size())>i){
	      finished=false;
	      tmp.push_back((*g._VECTptr)[i]);
	    }
	    else
	      tmp.push_back(zero);
	  }
	}
	if (finished)
	  break;
	if (n==1)
	  res.push_back(to_apply(tmp.front(),contextptr));
	else
	  res.push_back(to_apply(tmp,contextptr));
      }
      return res;
    }
    else
      setsizeerr();
    return 0;
  }
  const string _apply_s("apply");
  unary_function_eval __apply(&_apply,_apply_s);
  unary_function_ptr at_apply (&__apply,0,true);
  
  gen symb_makelist(const gen & args){
    return symbolic(at_makelist,args);
  }
  gen _makelist(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      setsizeerr();
    vecteur v(*args._VECTptr);
    int s=v.size();
    if (s<2)
      setsizeerr();
    gen f(v[0]),debut,fin,step(1);
    if (v[1].is_symb_of_sommet(at_interval)){
      debut=v[1]._SYMBptr->feuille._VECTptr->front();
      fin=v[1]._SYMBptr->feuille._VECTptr->back();
      if (s>2)
	step=v[2];
    }
    else {
      if (s<3)
	setsizeerr();
      debut=v[1];
      fin=v[2];
      if (s>3)
	step=v[3];
    }
    if (is_zero(step))
      setsizeerr("Invalid null step");
    vecteur w;
    if (is_greater(fin,debut,contextptr)){
      step=abs(step,contextptr);
      for (gen i=debut;is_greater(fin,i,contextptr);i=i+step)
	w.push_back(f(i,contextptr));
    }
    else {
      step=-abs(step,contextptr);
      for (gen i=debut;is_greater(i,fin,contextptr);i=i+step)
	w.push_back(f(i,contextptr));
    }
    return w;
  }
  const string _makelist_s("makelist");
  unary_function_eval __makelist(&_makelist,_makelist_s);
  unary_function_ptr at_makelist (&__makelist,0,true);

  gen symb_interval(const gen & args){
    return symbolic(at_interval,args);
  }
  gen symb_interval(const gen & a,const gen & b){
    return symbolic(at_interval,makevecteur(a,b));
  }
  gen _interval(const gen & args){
    return symb_interval(args);
  }
  const string _interval_s(" .. ");
  unary_function_unary __interval(&_interval,_interval_s,&printsommetasoperator,&texprintsommetasoperator);
  unary_function_ptr at_interval (&__interval);
  
  string printascomment(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if (feuille.type!=_STRNG)
        return sommetstr+'('+feuille.print(contextptr)+')';
    string chaine=*feuille._STRNGptr;
    int s=chaine.size();
    if ( (xcas_mode(contextptr)==1) || (xcas_mode(contextptr)==1+_DECALAGE)){
        string res("# ");
        for (int i=0;i<s;++i){
            if ((i==s-1)||(chaine[i]!='\n'))
                res +=chaine[i];
            else
                res += indent(contextptr)+"# ";
        }
        return res;
    }
    int l=chaine.find_first_of('\n');
    if ((l<0)|| (l>=s))
        return "//"+chaine + indent(contextptr);
    return "/*"+chaine+"*/";
  }
  gen symb_comment(const gen & args){
    return symbolic(at_comment,args);
  }
  gen _comment(const gen & args,GIAC_CONTEXT){
    return symb_comment(args);
  }
  const string _comment_s("comment");
  unary_function_eval __comment(&_comment,_comment_s,&printascomment);
  unary_function_ptr at_comment (&__comment);

  gen symb_throw(const gen & args){
    return symbolic(at_throw,args);
  }
  gen _throw(const gen & args,GIAC_CONTEXT){
    setsizeerr(args.print(contextptr));
    return undef;
  }
  const string _throw_s("throw");
  unary_function_eval __throw(&_throw,_throw_s);
  unary_function_ptr at_throw (&__throw);

  gen symb_union(const gen & args){
    return symbolic(at_union,args);
  }
  gen _union(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      setsizeerr();
    vecteur & v =*args._VECTptr;
    if (v.empty())
      return args;
    if (v.size()==1 && v.front().type==_VECT)
      return gen(v,_SET__VECT).eval(1,contextptr);
    if (v.size()!=2)
      setsizeerr();
    gen a=v.front(),b=v.back();
    if ( (a.type!=_VECT) || (b.type!=_VECT))
      setsizeerr("Union");
    return gen(mergevecteur(*a._VECTptr,*b._VECTptr),_SET__VECT).eval(1,contextptr);
  }
  const string _union_s(" union ");
  unary_function_eval __union(&_union,_union_s,&printsommetasoperator,&texprintsommetasoperator);
  unary_function_ptr at_union (&__union);

  gen symb_intersect(const gen & args){
    return symbolic(at_intersect,args);
  }
  gen _intersect(const gen & args,GIAC_CONTEXT){
    if ((args.type!=_VECT) || (args._VECTptr->size()!=2))
      setsizeerr();
    gen a=args._VECTptr->front(),b=args._VECTptr->back();
    if ( a.type==_VECT && b.type==_VECT){
      vecteur v;
      const_iterateur it=a._VECTptr->begin(),itend=a._VECTptr->end();
      for (;it!=itend;++it){
	if (equalposcomp(*b._VECTptr,*it))
	  v.push_back(*it);
      }
      return gen(v,_SET__VECT);
    }
    setsizeerr();
    return 0;
  }
  const string _intersect_s(" intersect ");
  unary_function_eval __intersect(&_intersect,_intersect_s,&printsommetasoperator,&texprintsommetasoperator);
  unary_function_ptr at_intersect (&__intersect);

  gen _inter(const gen & args,GIAC_CONTEXT){
    vecteur attributs(1,default_color(contextptr));
    vecteur v(seq2vecteur(args));
    int s=read_attributs(v,attributs,contextptr);
    if (s!=2)
      setdimerr();
    gen a=v[0],b=v[1];
    a=gen(inter(a,b,contextptr),_GROUP__VECT);
    return put_attributs(a,attributs,contextptr);
  }
  const string _inter_s("inter");
  unary_function_eval __inter(&_inter,_inter_s);
  unary_function_ptr at_inter (&__inter);

  gen symb_minus(const gen & args){
    return symbolic(at_minus,args);
  }
  gen _minus(const gen & args){
    if ((args.type!=_VECT) || (args._VECTptr->size()!=2))
      return symb_minus(args);
    gen a=args._VECTptr->front(),b=args._VECTptr->back();
    if ( (a.type!=_VECT) || (b.type!=_VECT))
      setsizeerr("Minus");
    vecteur v;
    const_iterateur it=a._VECTptr->begin(),itend=a._VECTptr->end();
    for (;it!=itend;++it){
      if (!equalposcomp(*b._VECTptr,*it))
	v.push_back(*it);
    }
    return gen(v,_SET__VECT);
  }
  const string _minus_s(" minus ");
  unary_function_unary __minus(&_minus,_minus_s,&printsommetasoperator,&texprintsommetasoperator);
  unary_function_ptr at_minus (&__minus);

  string printasdollar(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if (feuille.type!=_VECT)
      return sommetstr+feuille.print(contextptr);
    vecteur & v=*feuille._VECTptr;
    int s=v.size();
    if (s==2)
      return printsommetasoperator(feuille,sommetstr,contextptr);
    if (s==3)
      return v[0].print(contextptr)+sommetstr+v[1].print(contextptr)+" in "+v[2].print(contextptr);
    return "error";
  }
  gen symb_dollar(const gen & args){
    return symbolic(at_dollar,args);
  }
  gen _dollar(const gen & args,GIAC_CONTEXT){
    vecteur vargs;
    if (args.type!=_VECT){
      identificateur tmp(" _t");
      vargs=makevecteur(tmp,symbolic(at_equal,makevecteur(tmp,args)));
    }
    else
      vargs=*args._VECTptr;
    int s=vargs.size();
    if (s<2)
      return symb_dollar(args);
    gen a=vargs.front(),b=vargs[1],b1=eval(b,eval_level(contextptr),contextptr);
    if (b1.type==_INT_ && b1.val>=0)
      return gen(vecteur(b1.val,eval(a,eval_level(contextptr),contextptr)),_SEQ__VECT);
    gen var,intervalle,step=1;
    if ( (b.type==_SYMB) && (b._SYMBptr->sommet==at_equal || b._SYMBptr->sommet==at_same ) ){
      var=b._SYMBptr->feuille._VECTptr->front();
      if (var.type!=_IDNT)
	setsizeerr();
      /* Commented example seq(irem(g&^((p-1)/(Div[i])),p),i=(1 .. 2))
	 bool status=*var._IDNTptr->quoted;
	 *var._IDNTptr->quoted=true;
	 a=eval(a,contextptr);
	 *var._IDNTptr->quoted=status;      
      */
      intervalle=eval(b._SYMBptr->feuille._VECTptr->back(),eval_level(contextptr),contextptr);
      if (s>=3)
	step=vargs[2];
    }
    else {
      if (s>=3){
	var=vargs[1];
	intervalle=eval(vargs[2],eval_level(contextptr),contextptr);
      }
      if (s>=4)
	step=vargs[3];
    }
    if (intervalle.type==_VECT){
      const_iterateur it=intervalle._VECTptr->begin(),itend=intervalle._VECTptr->end();
      vecteur res;
      for (;it!=itend;++it)
	res.push_back(eval(quotesubst(a,var,*it,contextptr),eval_level(contextptr),contextptr));
      return gen(res,_SEQ__VECT);
      // return gen(res,intervalle.subtype);
    }
    if ( (intervalle.type==_SYMB) && (intervalle._SYMBptr->sommet==at_interval)){
      gen c=intervalle._SYMBptr->feuille._VECTptr->front(),d=intervalle._SYMBptr->feuille._VECTptr->back();
      gen debut=c,fin=d;
      bool reverse=ck_is_greater(debut,fin,contextptr);
      step=abs(step,contextptr);
      step=eval(reverse?-step:step,eval_level(contextptr),contextptr);
      vecteur res;
      for (;;debut+=step){
	if (ck_is_strictly_greater(reverse?fin:debut,reverse?debut:fin,contextptr))
	  break;
	res.push_back(eval(quotesubst(a,var,debut,contextptr),eval_level(contextptr),contextptr));
	if (debut==fin)
	  break;
      }
      return gen(res,_SEQ__VECT);
    }
    return symb_dollar(args);    
  }
  const string _dollar_s("$");
  string texprintasdollar(const gen & g,const string & s,GIAC_CONTEXT){
    if ( (g.type==_VECT) && (g._VECTptr->size()==2))
      return gen2tex(g._VECTptr->front(),contextptr)+"\\$"+gen2tex(g._VECTptr->back(),contextptr);
    return "\\$ "+g.print(contextptr);
  }
  unary_function_eval __dollar(&_dollar,_dollar_s,&printasdollar,&texprintasdollar);
  unary_function_ptr at_dollar (&__dollar,_QUOTE_ARGUMENTS,0);

  gen symb_makemat(const gen & args){
    return symbolic(at_makemat,args);
  }
  gen _makemat(const gen & args,const context * contextptr){
    if (args.type!=_VECT)
      return symb_makemat(args);
    int s=args._VECTptr->size();
    if ( (s!=3) && (s!=2) )
        return symb_makemat(args);
    gen fonction,intervalle1,intervalle2;
    if (s==3){
        fonction=args._VECTptr->front();
        intervalle1=(*(args._VECTptr))[1];
        intervalle2=args._VECTptr->back();
    }
    else {
        intervalle1=args._VECTptr->front();
        intervalle2=args._VECTptr->back();
    }
    if (intervalle1.type==_INT_)
        intervalle1=symb_interval(makevecteur(zero,intervalle1-1));
    if (intervalle2.type==_INT_)
        intervalle2=symb_interval(makevecteur(zero,intervalle2-1));
    if ( (intervalle1.type!=_SYMB) || (intervalle1._SYMBptr->sommet!=at_interval) ||(intervalle2.type!=_SYMB) || (intervalle2._SYMBptr->sommet!=at_interval))
        setsizeerr("makemat");
    intervalle1=intervalle1._SYMBptr->feuille;
    intervalle2=intervalle2._SYMBptr->feuille;
    if ((intervalle1.type!=_VECT) || (intervalle1._VECTptr->size()!=2) || (intervalle2.type!=_VECT) || (intervalle2._VECTptr->size()!=2))
      setsizeerr("interval");
    gen debut_i=intervalle1._VECTptr->front(),fin_i=intervalle1._VECTptr->back();
    gen debut_j=intervalle2._VECTptr->front(),fin_j=intervalle2._VECTptr->back();
    if ( (debut_i.type!=_INT_) || (fin_i.type!=_INT_) || (debut_j.type!=_INT_) || (fin_j.type!=_INT_) )
      setsizeerr("Boundaries not integer");
    int di=debut_i.val,fi=fin_i.val,dj=debut_j.val,fj=fin_j.val;
    int stepi=1,stepj=1;
    if (di>fi)
        stepi=-1;
    if (dj>fj)
        stepj=-1;
    if ((fonction.type!=_SYMB) || (fonction._SYMBptr->sommet!=at_program)){
      int s=(fj-dj+1)*stepj;
      vecteur w(s,fonction);
      int t=(fi-di+1)*stepi;
      vecteur res(t);
      for (int i=0;i<t;++i)
	res[i]=w; // each element of res will be a free line, so that =< works
      return res;
    }
    vecteur v,w,a(2);
    v.reserve((fi-di)*stepi);
    w.reserve((fj-dj)*stepj);
    for (;;di+=stepi){
        a[0]=di;
        w.clear();
        for (int djj=dj;;djj+=stepj){
            a[1]=djj;
            w.push_back(fonction(gen(a,_SEQ__VECT),contextptr));
            if (djj==fj)
                break;
        }
        v.push_back(w);
        if (di==fi)
            break;
    }
    return v;
  }
  const string _makemat_s("makemat");
  unary_function_eval __makemat(&_makemat,_makemat_s);
  unary_function_ptr at_makemat (&__makemat,0,true);

  gen symb_compose(const gen & args){
    return symbolic(at_compose,args);
  }
  gen _compose(const gen & args){
    return symb_compose(args);
  }
  const string _compose_s("@");
  unary_function_unary __compose(&_compose,_compose_s,&printsommetasoperator,&texprintsommetasoperator);
  unary_function_ptr at_compose (&__compose);

  gen _composepow(const gen & args){
    return symbolic(at_composepow,args);
  }
  const string _composepow_s("@@");
  unary_function_unary __composepow(&_composepow,_composepow_s,&printsommetasoperator,&texprintsommetasoperator);
  unary_function_ptr at_composepow (&__composepow);

  gen symb_args(const gen & args){
    return symbolic(at_args,args);
  }
  gen _args(const gen & args,const context * contextptr){
    gen e;
    if (debug_ptr(contextptr)->args_stack.empty())
      e=vecteur(0);
    else
      e=debug_ptr(contextptr)->args_stack.back();
    if ( (args.type==_VECT) && (args._VECTptr->empty()))
      return e;
    else
      return e(args,contextptr);
  }
  const string _args_s("args");
  unary_function_eval __args(&_args,_args_s);
  unary_function_ptr at_args (&__args);
  
  gen symb_lname(const gen & args){
    return symbolic(at_lname,args);
  }
  void lidnt(const vecteur & v,vecteur & res){
    const_iterateur it=v.begin(),itend=v.end();
    for (;it!=itend;++it)
      lidnt(*it,res);
  }
  void lidnt(const gen & args,vecteur & res){
    switch (args.type){
    case _IDNT:
      if (!equalposcomp(res,args))
	res.push_back(args);
      break;
    case _SYMB:
      if (args._SYMBptr->sommet==at_program && args._SYMBptr->feuille.type==_VECT && args._SYMBptr->feuille._VECTptr->size()==3){
	lidnt(args._SYMBptr->feuille._VECTptr->front(),res);
	lidnt(args._SYMBptr->feuille._VECTptr->back(),res);
      }
      else
	lidnt(args._SYMBptr->feuille,res);
      break;
    case _VECT:
      lidnt(*args._VECTptr,res);
      break;
    }       
  }
  vecteur lidnt(const gen & args){
    vecteur res;
    lidnt(args,res);
    return res;
  }
  gen _lname(const gen & args){
    return lidnt(args);
  }
  const string _lname_s("lname");
  unary_function_unary __lname(&_lname,_lname_s);
  unary_function_ptr at_lname (&__lname,0,true);

  gen symb_has(const gen & args){
    return symbolic(at_has,args);
  }
  gen _has(const gen & args){
      if ( (args.type!=_VECT) || (args._VECTptr->size()!=2))
          return symb_has(args);
      return equalposcomp(*_lname(args._VECTptr->front())._VECTptr,args._VECTptr->back());
  }
  const string _has_s("has");
  unary_function_unary __has(&_has,_has_s);
  unary_function_ptr at_has (&__has,0,true);

  gen _kill(const gen & args,GIAC_CONTEXT){
    if (args.type==_VECT && args._VECTptr->empty()){
      if (!contextptr)
	protection_level=0;
      debug_ptr(contextptr)->debug_mode=false;
      debug_ptr(contextptr)->current_instruction_stack.clear();
      debug_ptr(contextptr)->sst_at_stack.clear();
      debug_ptr(contextptr)->args_stack.clear();
      setsizeerr("Program killed");
      return undef;
    }
#ifdef HAVE_LIBPTHREAD
    if (args.type==_VECT)
      return apply(args,_kill,contextptr);
    if (args.type==_POINTER_ && args.subtype==_THREAD_POINTER){
      context * cptr=(context *) args._POINTER_val;
      thread_param * tptr =thread_param_ptr(cptr);
      if (cptr && tptr->eval_thread){
	gen g=tptr->v[0];
	if (g.type==_VECT && g._VECTptr->size()==2 && g._VECTptr->front().is_symb_of_sommet(at_quote)){
	  pthread_mutex_lock(cptr->globalptr->_mutex_eval_status_ptr);
	  sto(undef,g._VECTptr->front()._SYMBptr->feuille,cptr);
	  pthread_mutex_unlock(cptr->globalptr->_mutex_eval_status_ptr);
	}
      }
      kill_thread(true,cptr);
      return 1;
    }
#endif
    settypeerr();
    return undef;
  }
  const string _kill_s("kill");
  unary_function_eval __kill(&_kill,_kill_s);
  unary_function_ptr at_kill (&__kill,0,true);

  gen _halt(const gen & args,GIAC_CONTEXT){
    if (debug_ptr(contextptr)->debug_allowed){
      debug_ptr(contextptr)->debug_mode=true;
      debug_ptr(contextptr)->sst_mode=true;
      return plus_one;
    }
    return zero;
  }
  const string _halt_s("halt");
  unary_function_eval __halt(&_halt,_halt_s);
  unary_function_ptr at_halt (&__halt,_QUOTE_ARGUMENTS,true);

  gen _debug(const gen & args,GIAC_CONTEXT){
    if (child_id && thread_eval_status(contextptr)!=1)
      return args;
    if (debug_ptr(contextptr)->debug_allowed){
      debug_ptr(contextptr)->debug_mode=true;
      debug_ptr(contextptr)->sst_in_mode=true;
      debug_ptr(contextptr)->debug_prog_name=0;
    }
    return args.eval(eval_level(contextptr),contextptr);
  }
  const string _debug_s("debug");
  unary_function_eval __debug(&_debug,_debug_s);
  unary_function_ptr at_debug (&__debug,_QUOTE_ARGUMENTS,true);

  gen _sst_in(const gen & args,GIAC_CONTEXT){
    if (child_id)
      return zero;
    if (debug_ptr(contextptr)->debug_allowed){
      debug_ptr(contextptr)->debug_mode=true;
      debug_ptr(contextptr)->sst_in_mode=true;
      return plus_one;
    }
    return zero;
  }
  const string _sst_in_s("sst_in");
  unary_function_eval __sst_in(&_sst_in,_sst_in_s);
  unary_function_ptr at_sst_in (&__sst_in,_QUOTE_ARGUMENTS,true);

  gen _sst(const gen & args,GIAC_CONTEXT){
    if (child_id)
      return args;
    if (debug_ptr(contextptr)->debug_allowed){
      debug_ptr(contextptr)->debug_mode=true;
      debug_ptr(contextptr)->sst_mode=true;
      return plus_one;
    }
    return zero;
  }
  const string _sst_s("sst");
  unary_function_eval __sst(&_sst,_sst_s);
  unary_function_ptr at_sst (&__sst,_QUOTE_ARGUMENTS,true);

  gen _cont(const gen & args,GIAC_CONTEXT){
    if (child_id)
      return args;
    if (debug_ptr(contextptr)->debug_allowed){
      debug_ptr(contextptr)->sst_mode=false;
      return plus_one;
    }
    return zero;
  }
  const string _cont_s("cont");
  unary_function_eval __cont(&_cont,_cont_s);
  unary_function_ptr at_cont (&__cont,_QUOTE_ARGUMENTS,true);

  gen watch(const gen & args,GIAC_CONTEXT){
    if (!equalposcomp(debug_ptr(contextptr)->debug_watch,args))
      debug_ptr(contextptr)->debug_watch.push_back(args);
    return args;
  }
  gen _watch(const gen & args,GIAC_CONTEXT){
    if (child_id && thread_eval_status(contextptr)!=1 )
      return args;
    if (args.type==_VECT && args._VECTptr->empty() && debug_ptr(contextptr)->debug_localvars)
      apply( *debug_ptr(contextptr)->debug_localvars,contextptr,watch);
    else
      apply(args,contextptr,watch);
    return debug_ptr(contextptr)->debug_watch;
  }
  const string _watch_s("watch");
  unary_function_eval __watch(&_watch,_watch_s);
  unary_function_ptr at_watch (&__watch,_QUOTE_ARGUMENTS,true);

  gen rmwatch(const gen & args,GIAC_CONTEXT){
    int pos;
    if (args.type==_INT_){
      pos=args.val+1;
      if (pos>signed(debug_ptr(contextptr)->debug_watch.size()))
	return debug_ptr(contextptr)->debug_watch;
    }
    else 
      pos=equalposcomp(debug_ptr(contextptr)->debug_watch,args);
    if (pos){
      debug_ptr(contextptr)->debug_watch.erase(debug_ptr(contextptr)->debug_watch.begin()+pos-1,debug_ptr(contextptr)->debug_watch.begin()+pos);
      return debug_ptr(contextptr)->debug_watch;
    }
    return zero;
  }

  gen _rmwatch(const gen & args,GIAC_CONTEXT){
    if (args.type==_VECT && args._VECTptr->empty() && debug_ptr(contextptr)->debug_localvars)
      return apply( *debug_ptr(contextptr)->debug_localvars,contextptr,rmwatch);
    else
      return apply(args,contextptr,rmwatch);
  }
  const string _rmwatch_s("rmwatch");
  unary_function_eval __rmwatch(&_rmwatch,_rmwatch_s);
  unary_function_ptr at_rmwatch (&__rmwatch,_QUOTE_ARGUMENTS,true);

  gen _breakpoint(const gen & args,GIAC_CONTEXT){
    if (child_id && thread_eval_status(contextptr)!=1)
      return args;
    if ( (args.type!=_VECT) || (args._VECTptr->size()!=2) || (args._VECTptr->front().type!=_IDNT) || (args._VECTptr->back().type!=_INT_) )
      return zero;
    if (!equalposcomp(debug_ptr(contextptr)->debug_breakpoint,args)){
      debug_ptr(contextptr)->debug_breakpoint.push_back(args);
      // FIXME should also modify debug_ptr(contextptr)->sst_at_stack if the breakpoint applies
      // to a program != current program
      if (!debug_ptr(contextptr)->args_stack.empty() && debug_ptr(contextptr)->args_stack.back().type==_VECT && debug_ptr(contextptr)->args_stack.back()._VECTptr->front()==args._VECTptr->front())
	debug_ptr(contextptr)->sst_at.push_back(args._VECTptr->back().val);
    }
    return debug_ptr(contextptr)->debug_breakpoint;
  }
  const string _breakpoint_s("breakpoint");
  unary_function_eval __breakpoint(&_breakpoint,_breakpoint_s);
  unary_function_ptr at_breakpoint (&__breakpoint,_QUOTE_ARGUMENTS,true);

  gen _rmbreakpoint(const gen & args,GIAC_CONTEXT){
    if (child_id&& thread_eval_status(contextptr)!=1)
      return args;
    int pos;
    if (args.type==_INT_){
      pos=args.val;
      if (pos<1 || pos>signed(debug_ptr(contextptr)->debug_breakpoint.size())){
	adjust_sst_at(*debug_ptr(contextptr)->debug_prog_name,contextptr);
	return debug_ptr(contextptr)->debug_breakpoint;
      }
    }
    else 
      pos=equalposcomp(debug_ptr(contextptr)->debug_breakpoint,args);
    if (pos){
      debug_ptr(contextptr)->debug_breakpoint.erase(debug_ptr(contextptr)->debug_breakpoint.begin()+pos-1,debug_ptr(contextptr)->debug_breakpoint.begin()+pos);
      adjust_sst_at(*debug_ptr(contextptr)->debug_prog_name,contextptr);
      return debug_ptr(contextptr)->debug_breakpoint;
    }
    return zero;
  }
  const string _rmbreakpoint_s("rmbreakpoint");
  unary_function_eval __rmbreakpoint(&_rmbreakpoint,_rmbreakpoint_s);
  unary_function_ptr at_rmbreakpoint (&__rmbreakpoint,_QUOTE_ARGUMENTS,true);

  void debug_loop(const gen &res,GIAC_CONTEXT){
    if (!debug_ptr(contextptr)->debug_allowed || (!debug_ptr(contextptr)->sst_mode && !equalposcomp(debug_ptr(contextptr)->sst_at,debug_ptr(contextptr)->current_instruction)) )
      return;
    // Detect thread debugging
    int thread_debug=thread_eval_status(contextptr);
    if (thread_debug>1)
      return;
    if (thread_debug==1){
      // Fill dbgptr->debug_info_ptr and fast_debug_info_ptr 
      // with debugging infos to be displayed
      debug_struct * dbgptr=debug_ptr(contextptr);
      vecteur w; 
      // w[0]=function, args,
      // w[1]=breakpoints
      // w[2] = instruction to eval or program if debugging a prog
      // w[3]= evaluation result
      // w[4]= current instruction number 
      // w[5] = watch vector, w[6] = watch values
      if (!debug_ptr(contextptr)->args_stack.empty()){
	w.push_back(debug_ptr(contextptr)->args_stack.back());
	w.push_back(vector_int_2_vecteur(debug_ptr(contextptr)->sst_at,contextptr));
      }
      else {
	w.push_back(undef);
	w.push_back(undef);
      }
      w.push_back((*debug_ptr(contextptr)->fast_debug_info_ptr));
      w.push_back(res);
      w.push_back(debug_ptr(contextptr)->current_instruction);
      vecteur dw=debug_ptr(contextptr)->debug_watch;
      if (contextptr && dw.empty()){
	// put the last 2 environments
	const context * cur=contextptr;
	sym_tab::const_iterator it=cur->tabptr->begin(),itend=cur->tabptr->end();
	for (;it!=itend;++it){
	  dw.push_back(identificateur(it->first));
	}
	if (cur->previous && cur->previous!=cur->globalcontextptr){
	  cur=cur->previous;
	  sym_tab::const_iterator it=cur->tabptr->begin(),itend=cur->tabptr->end();
	  for (;it!=itend;++it){
	    dw.push_back(identificateur(it->first));
	  }
	}
      }
      w.push_back(dw);
      // evaluate watch with debug_ptr(contextptr)->debug_allowed=false
      debug_ptr(contextptr)->debug_allowed=false;
      iterateur it=dw.begin(),itend=dw.end();
      for (;it!=itend;++it)
	*it=protecteval(*it,1,contextptr);
      w.push_back(dw);
      debug_ptr(contextptr)->debug_allowed=true;
      *dbgptr->debug_info_ptr=w;
      dbgptr->debug_refresh=false;
      // Switch to level 2, waiting for main
      thread_eval_status(2,contextptr);
      for (;;){
	// Wait until status is put back by main to level 1
	usleep(10000);
	if (thread_eval_status(contextptr)==1){
	  // the wait function of the main thread should put in debug_info_ptr
	  // the next instruction, here we check for sst/sst_in/cont/kill
	  if (dbgptr->fast_debug_info_ptr){
	    gen test=*dbgptr->fast_debug_info_ptr;
	    if (test.type==_SYMB)
	      test=test._SYMBptr->sommet;
	    if (test.type==_FUNC){
	      if (test==at_sst){
		dbgptr->sst_in_mode=false;
		dbgptr->sst_mode=true;
		return;
	      }
	      if (test==at_sst_in){
		dbgptr->sst_in_mode=true;
		dbgptr->sst_mode=true;
		return;
	      }
	      if (test==at_cont){
		dbgptr->sst_in_mode=false;
		dbgptr->sst_mode=false;
		return;
	      }
	      if (test==at_kill){
		_kill(0,contextptr);
		return;
	      }
	    } // end type _FUNC
	    // eval
	    w[2] = *dbgptr->fast_debug_info_ptr;
	    w[3] = *dbgptr->fast_debug_info_ptr = protecteval(w[2],1,contextptr);
	    *dbgptr->debug_info_ptr=w;
	    dbgptr->debug_refresh=true;
	  } // end if (*dbgptr->debug_info_ptr)
	  thread_eval_status(2,contextptr); // Back to level 2
	} // end if (thread_eval_status()==1)
      } // end endless for loop
    } // end thread debugging
#ifdef WIN32
    *logptr(contextptr) << "Sorry! Debugging requires a true operating system" << endl;
    *logptr(contextptr) << "Please try xcas on Linux or an Unix" << endl;
    return;
#else // WIN32
    if (child_id)
      return;
    vecteur w; 
    // w[0]=[function + args, breakpoints]
    // w[2]= res of last evaluation, 
    // w[3] = next instruction, w[4]=debug_ptr(contextptr)->current_instruction
    // w[5] = watch vector, w[6] = watch values
    // evaluate watch with debug_ptr(contextptr)->debug_allowed=false
    debug_ptr(contextptr)->debug_allowed=false;
    debug_ptr(contextptr)->debug_allowed=true;
    if (!debug_ptr(contextptr)->args_stack.empty()){
      w.push_back(makevecteur(debug_ptr(contextptr)->args_stack.back(),vector_int_2_vecteur(debug_ptr(contextptr)->sst_at,contextptr)));
    }
    else
      w.push_back(undef);
    w.push_back(undef);
    w.push_back(res);
    w.push_back((*debug_ptr(contextptr)->fast_debug_info_ptr));
    (*debug_ptr(contextptr)->fast_debug_info_ptr)=undef;
    w.push_back(debug_ptr(contextptr)->current_instruction);
    w.push_back(debug_ptr(contextptr)->debug_watch);
    w.push_back(undef);
    bool in_debug_loop=true;
    for (;in_debug_loop;){
      try {
	vecteur tmp=gen2vecteur(debug_ptr(contextptr)->debug_watch);
	iterateur it=tmp.begin(),itend=tmp.end();
	for (;it!=itend;++it)
	  *it=it->eval(1,contextptr);
	w[6]=tmp;
      }
      catch (std::runtime_error & error){
	w[6]=string2gen(error.what(),false);
      }
      ofstream child_out(cas_sortie_name().c_str());
      gen e(symbolic(at_debug,w));
      *logptr(contextptr) << "Archiving " << e << endl;
      archive(child_out,e,contextptr);
      archive(child_out,zero,contextptr);
      child_out << "Debugging\n" << '¤' ;
      child_out.close();
      kill_and_wait_sigusr2();
      ifstream child_in(cas_entree_name().c_str());
      w[1]= unarchive(child_in,list_one_letter__IDNT,contextptr);
      child_in.close();
      *logptr(contextptr) << "Click reads " << w[1] << endl;
      if (w[1].type==_SYMB){
	if (w[1]._SYMBptr->sommet==at_sst){
	  debug_ptr(contextptr)->sst_in_mode=false;
	  debug_ptr(contextptr)->sst_mode=true;
	  return;
	}
	if (w[1]._SYMBptr->sommet==at_sst_in){
	  debug_ptr(contextptr)->sst_in_mode=true;
	  debug_ptr(contextptr)->sst_mode=true;
	  return;
	}
	if (w[1]._SYMBptr->sommet==at_cont){
	  debug_ptr(contextptr)->sst_in_mode=false;
	  debug_ptr(contextptr)->sst_mode=false;
	  return;
	}
	if (w[1]._SYMBptr->sommet==at_kill){
	  _kill(0,contextptr);
	}
      }
      try {
	w[2] =w[1].eval(1,contextptr);
      }
      catch (std::runtime_error & error ){
	w[2]=string2gen(error.what(),false);
      }
      
    }
#endif // WIN32
  }

  string printasbackquote(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    return "`"+feuille.print(contextptr)+"`";
  }
  gen _backquote(const gen & args){
    return args;
  }
  const string _backquote_s("backquote");
  unary_function_unary __backquote(&_backquote,_backquote_s,&printasbackquote);
  unary_function_ptr at_backquote (&__backquote);

  string printasdouble_deux_points(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    gen a,b;
    check_binary(feuille,a,b);
    return a.print(contextptr)+"::"+b.print(contextptr)+ " ";
  }
  gen symb_double_deux_points(const gen & args){
    return symbolic(at_double_deux_points,args);
  }
  gen _double_deux_points(const gen & args,GIAC_CONTEXT){
    gen a,b,c;
    check_binary(args,a,b);
    c=b;
    if (b.is_symb_of_sommet(at_of))
      c=b._SYMBptr->feuille[0];
    string cs=c.print(contextptr);
    /* // following code not used since qualified names after export 
       // make b a symbolic not just the function name
    int l=cs.size(),j=0;
    for (;j<l-1;++j){
      if (cs[j]==':' && cs[j+1]==':')
	break;
    }
    if (j==l-1)
    */      
    cs=a.print(contextptr)+"::"+cs;
    sym_tab::const_iterator it=lexer_functions().find(cs);
    if (it!=lexer_functions().end()){
      c=it->second;
      if (b.is_symb_of_sommet(at_of))
	return c(b._SYMBptr->feuille[1],contextptr);
      else
	return c;
    }
    if (b.type==_FUNC) // ? should be != _IDNT 
      return b;
    if (b.type==_SYMB)
      return b.eval(eval_level(contextptr),contextptr);
    gen aa=a.eval(1,contextptr);
    if (aa.type==_VECT)
      return find_in_folder(*aa._VECTptr,b);
    return symb_double_deux_points(args);
  }
  const string _double_deux_points_s("double_deux_points");
  unary_function_eval __double_deux_points(&_double_deux_points,_double_deux_points_s,&printasdouble_deux_points);
  unary_function_ptr at_double_deux_points (&__double_deux_points,_QUOTE_ARGUMENTS);

  bool is_binary(const gen & args){
    return (args.type==_VECT) && (args._VECTptr->size()==2) ;
  }

  void check_binary(const gen & args,gen & a,gen & b){
    if ( (args.type!=_VECT) || (args._VECTptr->size()!=2) )
      setsizeerr();
    a=args._VECTptr->front();
    b=args._VECTptr->back();
  }

  bool maple2mupad(const gen & args,int in_maple_mode,int out_maple_mode,GIAC_CONTEXT){
    check_secure();
    gen a,b;
    check_binary(args,a,b);
    string as,bs;
    if (a.type==_IDNT)
      as=*a._IDNTptr->name;
    if (a.type==_STRNG)
      as=*a._STRNGptr;
    if (b.type==_IDNT)
      bs=*b._IDNTptr->name;
    if (b.type==_STRNG)
      bs=*b._STRNGptr;
    int save_maple_mode=xcas_mode(contextptr);
    xcas_mode(contextptr)=in_maple_mode;
    ifstream infile(as.c_str());
    vecteur v;
    try {
      readargs_from_stream(infile,v,list_one_letter__IDNT,contextptr);
    }
    catch (std::runtime_error & error ){
      xcas_mode(contextptr)=save_maple_mode;
      return false;
    }
    xcas_mode(contextptr)=out_maple_mode;
    ofstream outfile(bs.c_str());
    const_iterateur it=v.begin(),itend=v.end();
    for (;it!=itend;++it)
      outfile << *it << endl;
    xcas_mode(contextptr)=save_maple_mode;
    return true;
  }

  gen _maple2mupad(const gen & args,GIAC_CONTEXT){
    return maple2mupad(args,1,2,contextptr);
  }
  const string _maple2mupad_s("maple2mupad");
  unary_function_eval __maple2mupad(&_maple2mupad,_maple2mupad_s);
  unary_function_ptr at_maple2mupad (&__maple2mupad,_QUOTE_ARGUMENTS,true);

  gen _maple2xcas(const gen & args,GIAC_CONTEXT){
    return maple2mupad(args,1,0,contextptr);
  }
  const string _maple2xcas_s("maple2xcas");
  unary_function_eval __maple2xcas(&_maple2xcas,_maple2xcas_s);
  unary_function_ptr at_maple2xcas (&__maple2xcas,_QUOTE_ARGUMENTS,true);

  gen _mupad2maple(const gen & args,GIAC_CONTEXT){
    return maple2mupad(args,2,1,contextptr);
  }
  const string _mupad2maple_s("mupad2maple");
  unary_function_eval __mupad2maple(&_mupad2maple,_mupad2maple_s);
  unary_function_ptr at_mupad2maple (&__mupad2maple,_QUOTE_ARGUMENTS,true);

  gen _mupad2xcas(const gen & args,GIAC_CONTEXT){
    return maple2mupad(args,2,0,contextptr);
  }
  const string _mupad2xcas_s("mupad2xcas");
  unary_function_eval __mupad2xcas(&_mupad2xcas,_mupad2xcas_s);
  unary_function_ptr at_mupad2xcas (&__mupad2xcas,_QUOTE_ARGUMENTS,true);

  string printasvirgule(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if ( (feuille.type!=_VECT) || (feuille._VECTptr->size()!=2) )
      return sommetstr+'('+feuille.print(contextptr)+')';
    return feuille._VECTptr->front().print(contextptr)+','+feuille._VECTptr->back().print(contextptr);
  }
  gen _virgule(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      return args;
    const_iterateur it=args._VECTptr->begin(),itend=args._VECTptr->end();
    if (itend-it<2)
      return args;
    gen res=makesuite(*it,*(it+1));
    ++it;
    ++it;
    for (;it!=itend;++it)
      res=makesuite(res,*it);
    return res;
  }
  const string _virgule_s("virgule");
  unary_function_eval __virgule(&_virgule,_virgule_s,&printasvirgule);
  unary_function_ptr at_virgule (&__virgule);

  gen _pwd(const gen & args,GIAC_CONTEXT){
#ifndef HAVE_NO_CWD
    char * buffer=getcwd(0,0);
    if (buffer){
      string s(buffer);
#ifndef HAVE_LIBGC
      free(buffer);
#endif
      return string2gen(s,false);
    }
#endif
    setsizeerr();
    return 0;
  }
  const string _pwd_s("pwd");
  unary_function_eval __pwd(&_pwd,_pwd_s);
  unary_function_ptr at_pwd (&__pwd,0,true);

  gen _cd(const gen & args,GIAC_CONTEXT){
    check_secure();
    if (args.type!=_STRNG)
      settypeerr();
    int res;
    string s(*args._STRNGptr);
    string ss(*_pwd(zero,contextptr)._STRNGptr+'/'),current;
    int l=s.size();
    for (int i=0;i<=l;i++){
      if ( (i==l) || (s[i]=='/') ){
	if (i){
	  if (current==".."){
	    int t=ss.size()-2;
	    for (;t>0;--t){
	      if (ss[t]=='/')
		break;
	    }
	    if (t)
	      ss=ss.substr(0,t+1);
	    else
	      ss="/";
	  } 
	  else { // not ..
	    if (current[0]=='~'){
	      if (current.size()==1){ // uid user directory
		ss = home_directory();
	      }
	      else { // other user directory
		current=current.substr(1,current.size()-1);
#ifndef HAVE_NO_PWD_H
		passwd * p=getpwnam(current.c_str());
		if (!p)
		  setsizeerr("No such user "+current);
		ss = p->pw_dir ;
		ss +='/';
#else
		ss = "/";
#endif
	      }
	    }
	    else
	      ss+=current+"/";
	  } // end .. detection
	}
	else // i==0 / means absolute path
	  ss="/";
	current="";
      } // end / detection
      else {
	if (s[i]>' ')
	  current += s[i];
      }
    } // end for
#ifndef HAVE_NO_CWD
    res=chdir(ss.c_str());
#else
    res=-1;
#endif
    if (res)
      setsizeerr();
    gen g=symbolic(at_cd,_pwd(zero,contextptr));
    if (!child_id)
      _signal(symb_quote(g),contextptr);
    // *logptr(contextptr) << g << endl;
    return g;
  }
  const string _cd_s("cd");
  unary_function_eval __cd(&_cd,_cd_s);
  unary_function_ptr at_cd (&__cd,0,true);

  gen _scientific_format(const gen & g,GIAC_CONTEXT){
    check_secure();
    gen args(g);
    if (g.type==_DOUBLE_)
      args=int(g._DOUBLE_val);    
    if (args.type!=_INT_)
      return scientific_format(contextptr);
    scientific_format(args.val,contextptr);
    return args;
  }
  const string _scientific_format_s("scientific_format");
  unary_function_eval __scientific_format(&_scientific_format,_scientific_format_s,&printasDigits);
  unary_function_ptr at_scientific_format (&__scientific_format);

  gen _integer_format(const gen & g,GIAC_CONTEXT){
    check_secure();
    gen args(g);
    if (g.type==_DOUBLE_)
      args=int(g._DOUBLE_val);    
    if (args.type!=_INT_)
      return integer_format(contextptr);
    integer_format(args.val,contextptr);
    return args;
  }
  const string _integer_format_s("integer_format");
  unary_function_eval __integer_format(&_integer_format,_integer_format_s,&printasDigits);
  unary_function_ptr at_integer_format (&__integer_format,0,true);

  // 0: xcas, 1: maple, 2: mupad, 3: ti
  gen _xcas_mode(const gen & g,GIAC_CONTEXT){
    gen args(g);
    if (g.type==_DOUBLE_)
      args=int(g._DOUBLE_val);    
    if (args.type!=_INT_)
      return xcas_mode(contextptr);
    xcas_mode(contextptr)=args.val;
    return string2gen("Warning: some commands like subs might change arguments order",false);
  }
  const string _xcas_mode_s("xcas_mode");
  unary_function_eval __xcas_mode(&_xcas_mode,_xcas_mode_s);
  unary_function_ptr at_xcas_mode (&__xcas_mode,0,true);
  const string _maple_mode_s("maple_mode");
  unary_function_eval __maple_mode(&_xcas_mode,_maple_mode_s);
  unary_function_ptr at_maple_mode (&__maple_mode,0,true);

  gen _eval_level(const gen & g,GIAC_CONTEXT){
    check_secure();
    gen args(g);
    if (g.type==_DOUBLE_)
      args=int(g._DOUBLE_val);    
    if (args.type!=_INT_)
      return eval_level(contextptr);
    eval_level(contextptr)=args.val;
    DEFAULT_EVAL_LEVEL=args.val;
    return args;
  }
  const string _eval_level_s("eval_level");
  unary_function_eval __eval_level(&_eval_level,_eval_level_s,&printasDigits);
  unary_function_ptr at_eval_level (&__eval_level,0,true);

  gen _prog_eval_level(const gen & g,GIAC_CONTEXT){
    check_secure();
    gen args(g);
    if (g.type==_DOUBLE_)
      args=int(g._DOUBLE_val);    
    if (args.type!=_INT_)
      return prog_eval_level(contextptr);
    prog_eval_level_val(contextptr)=args.val;
    return args;
  }
  const string _prog_eval_level_s("prog_eval_level");
  unary_function_eval __prog_eval_level(&_prog_eval_level,_prog_eval_level_s,&printasDigits);
  unary_function_ptr at_prog_eval_level (&__prog_eval_level,0,true);

  gen _with_sqrt(const gen & g,GIAC_CONTEXT){
    check_secure();
    gen args(g);
    if (g.type==_DOUBLE_)
      args=int(g._DOUBLE_val);    
    if (args.type!=_INT_)
      return withsqrt(contextptr);
    withsqrt(contextptr)=args.val;
    return args;
  }
  const string _with_sqrt_s("with_sqrt");
  unary_function_eval __with_sqrt(&_with_sqrt,_with_sqrt_s,&printasDigits);
  unary_function_ptr at_with_sqrt (&__with_sqrt,0,true);

  gen _all_trig_solutions(const gen & g,GIAC_CONTEXT){
    gen args(g);
    if (g.type==_DOUBLE_)
      args=int(g._DOUBLE_val);    
    if (args.type!=_INT_)
      return all_trig_sol(contextptr);
    all_trig_sol(args.val,contextptr);
    parent_cas_setup(contextptr);
    return args;
  }
  const string _all_trig_solutions_s("all_trig_solutions");
  unary_function_eval __all_trig_solutions(&_all_trig_solutions,_all_trig_solutions_s,&printasDigits);
  unary_function_ptr at_all_trig_solutions (&__all_trig_solutions);

  gen _ntl_on(const gen & g,GIAC_CONTEXT){
    gen args(g);
    if (g.type==_DOUBLE_)
      args=int(g._DOUBLE_val);    
    if (args.type!=_INT_)
      return ntl_on(contextptr);
    ntl_on(args.val,contextptr);
    return args;
  }
  const string _ntl_on_s("ntl_on");
  unary_function_eval __ntl_on(&_ntl_on,_ntl_on_s,&printasDigits);
  unary_function_ptr at_ntl_on (&__ntl_on);

  gen _complex_mode(const gen & g,GIAC_CONTEXT){
    gen args(g);
    if (g.type==_DOUBLE_)
      args=int(g._DOUBLE_val);    
    if (args.type!=_INT_)
      return complex_mode(contextptr);
    complex_mode(args.val,contextptr);
    parent_cas_setup(contextptr);
    return args;
  }
  const string _complex_mode_s("complex_mode");
  unary_function_eval __complex_mode(&_complex_mode,_complex_mode_s,&printasDigits);
  unary_function_ptr at_complex_mode (&__complex_mode);

  gen _angle_radian(const gen & g,GIAC_CONTEXT){
    gen args(g);
    if (g.type==_DOUBLE_)
      args=int(g._DOUBLE_val);    
    if (args.type!=_INT_)
      return angle_radian(contextptr);
    angle_radian(args.val,contextptr);
    parent_cas_setup(contextptr);
    return args;
  }
  const string _angle_radian_s("angle_radian");
  unary_function_eval __angle_radian(&_angle_radian,_angle_radian_s,&printasDigits);
  unary_function_ptr at_angle_radian (&__angle_radian);

  gen _epsilon(const gen & arg,GIAC_CONTEXT){
    gen args=evalf_double(arg,0,contextptr);
    if (args.type!=_DOUBLE_)
      return epsilon(contextptr);
    epsilon(fabs(args._DOUBLE_val),contextptr);
    parent_cas_setup(contextptr);
    return args;
  }
  const string _epsilon_s("epsilon");
  unary_function_eval __epsilon(&_epsilon,_epsilon_s,&printasDigits);
  unary_function_ptr at_epsilon (&__epsilon);

  gen _proba_epsilon(const gen & arg,GIAC_CONTEXT){
    gen args=evalf_double(arg,0,contextptr);
    if (args.type!=_DOUBLE_)
      return proba_epsilon(contextptr);
    proba_epsilon(contextptr)=fabs(args._DOUBLE_val);
    parent_cas_setup(contextptr);
    return args;
  }
  const string _proba_epsilon_s("proba_epsilon");
  unary_function_eval __proba_epsilon(&_proba_epsilon,_proba_epsilon_s,&printasDigits);
  unary_function_ptr at_proba_epsilon (&__proba_epsilon);

  gen _complex_variables(const gen & g,GIAC_CONTEXT){
    gen args(g);
    if (g.type==_DOUBLE_)
      args=int(g._DOUBLE_val);    
    if (args.type!=_INT_)
      return complex_variables(contextptr);
    complex_variables(args.val,contextptr);
    parent_cas_setup(contextptr);
    return args;
  }
  const string _complex_variables_s("complex_variables");
  unary_function_eval __complex_variables(&_complex_variables,_complex_variables_s,&printasDigits);
  unary_function_ptr at_complex_variables (&__complex_variables);

  gen _approx_mode(const gen & g,GIAC_CONTEXT){
    gen args(g);
    if (g.type==_DOUBLE_)
      args=int(g._DOUBLE_val);
    if (args.type!=_INT_)
      return approx_mode(contextptr);
    approx_mode(args.val,contextptr);
    parent_cas_setup(contextptr);
    return args;
  }
  const string _approx_mode_s("approx_mode");
  unary_function_eval __approx_mode(&_approx_mode,_approx_mode_s,&printasDigits);
  unary_function_ptr at_approx_mode (&__approx_mode);

  gen _threads(const gen & g,GIAC_CONTEXT){
    gen args(g);
    if (g.type==_DOUBLE_)
      args=int(g._DOUBLE_val);    
    if (args.type!=_INT_)
      return threads;
    threads=max(absint(args.val),1);
    parent_cas_setup(contextptr);
    return args;
  }
  const string _threads_s("threads");
  unary_function_eval __threads(&_threads,_threads_s,&printasDigits);
  unary_function_ptr at_threads (&__threads);

  int digits2bits(int n){
#ifdef OLDGNUWINCE
    return (n*33)/10;
#else
    return int(floor(std::log(10.0)/std::log(2.0)*n))+1;
#endif
  }

  int bits2digits(int n){
#ifdef OLDGNUWINCE
    return (n*3)/10;
#else
    return int(floor(std::log(2.0)/std::log(10.0)*n))+1;
#endif
  }

  void set_decimal_digits(int n,GIAC_CONTEXT){
#ifdef GNUWINCE
    return undef;
#else
    decimal_digits(contextptr)=max(absint(n),1);
    // deg2rad_g=evalf(cst_pi,1,0)/180;
    // rad2deg_g=inv(deg2rad_g);
#endif
  }

  void cas_setup(const vecteur & v_orig,GIAC_CONTEXT){
    vecteur v(v_orig);
    if (v.size()<7)
      setsizeerr();
    if (*logptr(contextptr) && debug_infolevel) 
      *logptr(contextptr) << "Cas_setup " << v << char(10) << char(13) ;
    if (v[0].type==_INT_)
      approx_mode(v[0].val,contextptr);
    else {
      v[0]=evalf_double(v[0],1,contextptr);
      if (v[0].type==_DOUBLE_)
	approx_mode(int(v[0]._DOUBLE_val),contextptr);
    }
    if (v[1].type==_INT_)
      complex_variables(v[1].val,contextptr);
    else {
      v[1]=evalf_double(v[1],1,contextptr);
      if (v[1].type==_DOUBLE_)
	complex_variables(v[1]._DOUBLE_val,contextptr);
    }
    if (v[2].type==_INT_)
      complex_mode(v[2].val,contextptr);
    else {
      v[2]=evalf_double(v[2],1,contextptr);
      if (v[2].type==_DOUBLE_)
	complex_mode(v[2]._DOUBLE_val,contextptr);
    }
    if (v[3].type==_INT_)
      angle_radian(v[3].val,contextptr);
    else {
      v[3]=evalf_double(v[3],1,contextptr);
      if (v[3].type==_DOUBLE_)
	angle_radian(v[3]._DOUBLE_val,contextptr);
    }
    v[4]=evalf_double(v[4],1,contextptr);
    if (v[4].type==_DOUBLE_){
      int format=int(v[4]._DOUBLE_val);
      scientific_format(format % 16,contextptr);
      integer_format(format/16,contextptr);
    }
    v[5]=evalf_double(v[5],1,contextptr);
    if (v[5].type==_DOUBLE_)
      epsilon(fabs(v[5]._DOUBLE_val),contextptr);
    if (v[5].type==_VECT && v[5]._VECTptr->size()==2 && v[5]._VECTptr->front().type==_DOUBLE_ && v[5]._VECTptr->back().type==_DOUBLE_){
      epsilon(fabs(v[5]._VECTptr->front()._DOUBLE_val),contextptr);
      proba_epsilon(contextptr)=fabs(v[5]._VECTptr->back()._DOUBLE_val); 
    }
    if (v[6].type==_INT_)
      set_decimal_digits(v[6].val,contextptr);
    else {
      v[6]=evalf_double(v[6],1,contextptr);
      if (v[6].type==_DOUBLE_)
	set_decimal_digits(int(v[6]._DOUBLE_val),contextptr);
    }
    if (v.size()>=8){
      if (v[7].type==_VECT){
	vecteur & vv =*v[7]._VECTptr;
	if (vv.size()>=4){
	  threads=std::max(1,int(evalf_double(vv[0],1,contextptr)._DOUBLE_val));
	  MAX_RECURSION_LEVEL=std::max(int(evalf_double(vv[1],1,contextptr)._DOUBLE_val),1);
	  debug_infolevel=std::max(0,int(evalf_double(vv[2],1,contextptr)._DOUBLE_val));
	  DEFAULT_EVAL_LEVEL=std::max(1,int(evalf_double(vv[3],1,contextptr)._DOUBLE_val));
	}
      }
    }
    if (v.size()>=9){ 
      if (v[8].type==_INT_)
	increasing_power(v[8].val,contextptr);
      else {
	v[8]=evalf_double(v[8],1,contextptr);
	if (v[8].type==_DOUBLE_)
	  increasing_power(v[8]._DOUBLE_val,contextptr);
      }
    }
    if (v.size()>=10){ 
      if (v[9].type==_INT_)
	withsqrt(v[9].val,contextptr);
      else {
	v[9]=evalf_double(v[9],1,contextptr);
	if (v[9].type==_DOUBLE_)
	  withsqrt(v[9]._DOUBLE_val,contextptr);
      }
    }
    if (v.size()>=11){ 
      if (v[10].type==_INT_)
	all_trig_sol(v[10].val,contextptr);
      else {
	v[10]=evalf_double(v[10],1,contextptr);
	if (v[10].type==_DOUBLE_)
	  all_trig_sol(v[10]._DOUBLE_val,contextptr);
      }
    }
  }
  vecteur cas_setup(GIAC_CONTEXT){
    vecteur v;
    v.push_back(approx_mode(contextptr));
    v.push_back(complex_variables(contextptr));
    v.push_back(complex_mode(contextptr));
    v.push_back(angle_radian(contextptr));
    v.push_back(scientific_format(contextptr)+16*integer_format(contextptr));
    v.push_back(makevecteur(epsilon(contextptr),proba_epsilon(contextptr)));
    v.push_back(decimal_digits(contextptr));
    v.push_back(makevecteur(threads,MAX_RECURSION_LEVEL,debug_infolevel,DEFAULT_EVAL_LEVEL));
    v.push_back(increasing_power(contextptr));
    v.push_back(withsqrt(contextptr));
    v.push_back(all_trig_sol(contextptr));    
    return v;
  }
  gen _cas_setup(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      setsizeerr();
    vecteur & w=*args._VECTptr;
    cas_setup(w,contextptr);
    if (!child_id){
      _signal(symbolic(at_quote,symbolic(at_cas_setup,w)),contextptr);
    }
    return args;
  }
  const string _cas_setup_s("cas_setup");
  unary_function_eval __cas_setup(&giac::_cas_setup,_cas_setup_s);
  unary_function_ptr at_cas_setup (&__cas_setup,0,true);

  void parent_cas_setup(GIAC_CONTEXT){
    if (!child_id){
      _signal(symbolic(at_quote,symbolic(at_cas_setup,cas_setup(contextptr))),contextptr);
    }
  }

  string printasDigits(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if (feuille.type==_VECT && feuille._VECTptr->empty())
      return sommetstr;
    return sommetstr+" := "+feuille.print(contextptr);
  }
  gen _Digits(const gen & g,GIAC_CONTEXT){
    gen args(g);
    if (g.type==_DOUBLE_)
      args=int(g._DOUBLE_val);    
    if (args.type!=_INT_)
      return decimal_digits(contextptr);
    set_decimal_digits(args.val,contextptr);
    return _cas_setup(cas_setup(contextptr),contextptr);
  }
  const string _Digits_s("Digits");
  unary_function_eval __Digits(&giac::_Digits,_Digits_s,&printasDigits);
  unary_function_ptr at_Digits (&__Digits);

#ifdef GNUWINCE
  gen _xport(const gen & args){
#else
  gen _export(const gen & args){
#endif
    string libname(gen2string(args));
    std::map<std::string,std::vector<std::string> >::iterator it=library_functions().find(libname);
    if (it==library_functions().end())
      return zero;
    sort(it->second.begin(),it->second.end());
    // Add library function names to the translator
    std::vector<std::string>::iterator jt=it->second.begin(),jtend=it->second.end(),kt,ktend;
    for (;jt!=jtend;++jt){
      string tname=libname+"::"+*jt;
      // Find if the name exists in the translator base
      it=lexer_translator().find(*jt);
      if (it==lexer_translator().end())
	lexer_translator()[*jt]=vector<string>(1,tname);
      else { // Name exists, check if tname is in the vector, else push it
	kt=it->second.begin(); ktend=it->second.end();
	for (;kt!=ktend;++kt){
	  if (*kt==tname)
	    break;
	}
	if (kt!=ktend)
	  it->second.erase(kt);
	it->second.push_back(tname);
      }
    }
    return plus_one;
  }
  gen _insmod(const gen & args,GIAC_CONTEXT){
    if (args.type!=_STRNG)
#ifdef GNUWINCE
      return _xport(args);
#else
      return _export(args);
#endif
#ifdef HAVE_LIBDL
    string libname=*args._STRNGptr;
    if (libname.empty())
      return 0;
    // a way to add the current path to the search
    if (libname[0]!='/'){
      gen pwd=_pwd(0,contextptr);
      if (pwd.type==_STRNG){
	string libname1 = *pwd._STRNGptr+'/'+libname;
	if (is_file_available(libname1.c_str()))
	  libname=libname1;
      }
    }
    modules_tab::const_iterator i = giac_modules_tab.find(libname);
    if (i!=giac_modules_tab.end())
      return plus_two; // still registered
    registered_lexer_functions().clear();
    void * handle = dlopen (libname.c_str(), RTLD_LAZY);
    if (!handle) {
      setsizeerr (string(dlerror()));
    }
    // if (debug_infolevel)
    //  *logptr(contextptr) << registered_lexer_functions << endl;
    giac_modules_tab[libname]=module_info(registered_lexer_functions(),handle);
    if (!child_id)
      _signal(symb_quote(symbolic(at_insmod,args)),contextptr);
    else
      *logptr(contextptr) << "Parent insmod" <<endl;
#ifdef GNUWINCE
    _xport(args);
#else
    _export(args);
#endif
    return plus_one;
#else // HAVE_LIBDL
    return zero;
#endif // HAVE_LIBDL
  }
  const string _insmod_s("insmod");
  unary_function_eval __insmod(&_insmod,_insmod_s);
  unary_function_ptr at_insmod (&__insmod,0,true);

  gen _rmmod(const gen & args,GIAC_CONTEXT){
    if (args.type!=_STRNG)
      settypeerr();
#ifdef HAVE_LIBDL
    string libname=*args._STRNGptr;
    modules_tab::const_iterator i = giac_modules_tab.find(libname);
    if (i==giac_modules_tab.end())
      return plus_two; // not registered
    dlclose(i->second.handle);
    bool res= lexer_function_remove(i->second.registered_names);
    giac_modules_tab.erase(libname);
    if (!child_id)
      _signal(symb_quote(symbolic(at_rmmod,args)),contextptr);
    return(res);
#else // HAVE_LIBDL
    return zero;
#endif // HAVE_LIBDL
  }

  /*
  gen _rmmod(const gen & args){
    if (args.type==_VECT)
      apply(args,giac::rmmod);
    rmmod(args);    
  }
  */
  const string _rmmod_s("rmmod");
  unary_function_eval __rmmod(&_rmmod,_rmmod_s);
  unary_function_ptr at_rmmod (&__rmmod,0,true);

  gen _lsmod(const gen & args,GIAC_CONTEXT){
    vecteur v;
    modules_tab::const_iterator i = giac_modules_tab.begin(),iend=giac_modules_tab.end();
    for (;i!=iend;++i)
      v.push_back(string2gen(i->first,false));
    return v;
  }
  const string _lsmod_s("lsmod");
  unary_function_eval __lsmod(&_lsmod,_lsmod_s);
  unary_function_ptr at_lsmod (&__lsmod,0,true);

  class gen_sort {
    gen sorting_function;
    const context * contextptr;
  public:
    bool operator () (const gen & a,const gen & b){
      gen c=sorting_function(gen(makevecteur(a,b),_SEQ__VECT),contextptr);
      if (c.type!=_INT_)
	setsizeerr("Unable to sort "+c.print(contextptr));
      return !is_zero(c);
    }
    gen_sort(const gen & f,const context * ptr): sorting_function(f),contextptr(ptr) {};
  };

  /*
  gen sorting_function;
  bool sort_sort(const gen & a,const gen & b){
    gen c=sorting_function(gen(makevecteur(a,b),_SEQ__VECT),0);
    if (c.type!=_INT_)
      setsizeerr("Unable to sort "+c.print(contextptr));
    return !is_zero(c);
  }
  */

  gen simplifier(const gen & g,GIAC_CONTEXT){
    return liste2symbolique(symbolique2liste(g,contextptr));
  }
  gen _simplifier(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT)
      return simplifier(g,contextptr);
    return apply(g,_simplifier,contextptr);
  }

  gen _sort(const gen & args,GIAC_CONTEXT){
    if (args.type==_SYMB)
      return simplifier(args,contextptr);
    if (args.type!=_VECT)
      return args; // FIXME sort in additions, symbolic(at_sort,args);
    vecteur v=*args._VECTptr;
    int subtype;
    gen f;
    if ( v.size()==2 && v[0].type==_VECT ){
      f=v[1];
      subtype=v[0].subtype;
      v=*v[0]._VECTptr;
    }
    else {
      f=at_inferieur_strict;
      subtype=args.subtype;
    }
    sort(v.begin(),v.end(),gen_sort(f,contextptr));
    return gen(v,subtype);
  }
  const string _sort_s("sort");
  unary_function_eval __sort(&_sort,_sort_s);
  unary_function_ptr at_sort (&__sort,0,true);

  gen remove_nodisp(const gen & g){
    if (g.is_symb_of_sommet(at_nodisp))
      return g._SYMBptr->feuille;
    return g;
  }
  gen _ans(const gen & args,GIAC_CONTEXT){
    int s=history_out(contextptr).size();
    if (!s)
      return undef;
    int i;
    if (args.type!=_INT_)
      i=-1;
    else {
      i=args.val;
      if (xcas_mode(contextptr)==3)
	i=-i;
    }
    if (i>=0){
      if (i>=s)
	toofewargs(print_INT_(i));
      return remove_nodisp(history_out(contextptr)[i]);
    }
    if (s+i<0)
      toofewargs(print_INT_(-i));
    return remove_nodisp(history_out(contextptr)[s+i]);
  }
  const string _ans_s("ans");
  unary_function_eval __ans(&_ans,_ans_s);
  unary_function_ptr at_ans (&__ans,0,true);

  gen _quest(const gen & args,GIAC_CONTEXT){
    if (rpn_mode)
      setsizeerr();
    int s=history_in(contextptr).size();
    if (!s)
      return undef;
    int i;
    if (args.type!=_INT_)
      i=-2;
    else
      i=args.val;
    if (i>=0){
      if (i>=s)
	toofewargs(print_INT_(i));
      return remove_nodisp(history_in(contextptr)[i]);
    }
    if (s+i<0)
      toofewargs(print_INT_(-i));
    return remove_nodisp(history_in(contextptr)[s+i]);
  }
  const string _quest_s("quest");
  unary_function_eval __quest(&_quest,_quest_s);
  unary_function_ptr at_quest (&__quest,0,true);

  vector<int> float2continued_frac(double d_orig,double eps){
    double d=fabs(d_orig);
    if (d>RAND_MAX)
      setsizeerr(); 
    vector<int> v;
    double i;
    for (;;){
      i=floor(d);
      v.push_back(int(i));
      d=d-i;
      if (d<eps)
	return v;
      d=1/d;
      eps=eps*d*d;
    }
  }

  gen continued_frac2gen(vector<int> v,double d_orig,double eps,GIAC_CONTEXT){
    gen res(v.back());
    for (;;){
      v.pop_back();
      if (v.empty()){
	if (
	    !my_isnan(d_orig) &&
	    fabs(evalf_double(res-d_orig,1,contextptr)._DOUBLE_val)>eps)
	  return d_orig;
	return res;
      }
      res=inv(res,contextptr);
      res=res+v.back();
    }
    return res;
  }

  gen chk_not_unit(const gen & g){
    if (g.is_symb_of_sommet(at_unit))
      setsizeerr("Incompatible units");
    return g;
  }

  gen _convert(const gen & args,const context * contextptr){
    if (args.type!=_VECT)
      setsizeerr();
    vecteur & v=*args._VECTptr;
    int s=v.size();
    if (s<2)
      setsizeerr();
    gen & f=v[1];
    gen g=v.front();
    if (f.is_symb_of_sommet(at_unit)){
      return chk_not_unit(mksa_reduce(g/f,contextptr))*f;
    }
    if (s==3 && f.type==_INT_ ){
      if (f.val==_BASE && v.back().type==_INT_ ){
	if (is_integer(g)){
	  // convert(integer,base,integer)
	  bool positif=is_positive(g,contextptr);
	  g=abs(g,contextptr);
	  vecteur res;
	  gen q;
	  for (;!is_zero(g);){
	    res.push_back(irem(g,v.back(),q));
	    g=q;
	  }
	  // reverse(res.begin(),res.end());
	  if (positif)
	    return res;
	  return -res;
	}
	if (g.type==_VECT){
	  vecteur w(*g._VECTptr);
	  reverse(w.begin(),w.end());
	  return horner(w,v.back());
	}
      }
      if (f.val==_CONFRAC && v.back().type==_IDNT){
	g=evalf_double(g,1,contextptr);
	if (g.type==_DOUBLE_)
	  return sto(vector_int_2_vecteur(float2continued_frac(g._DOUBLE_val,epsilon(contextptr)),contextptr),v.back(),contextptr);
      }
    }
    if (s>2)
      g=gen(mergevecteur(vecteur(1,g),vecteur(v.begin()+2,v.begin()+s)),args.subtype);
    if (v[1].type==_FUNC){
      if (f==at_sincos)
	return sincos(g,contextptr);
      if (f==at_sin)
	return trigsin(g,contextptr);
      if (f==at_cos)
	return trigcos(g,contextptr);
      if (f==at_tan)
	return halftan(g,contextptr);
      if (f==at_exp || f==at_ln)
	return trig2exp(g,contextptr);
      if (f==at_string)
	return string2gen(g.print(contextptr),false);
      if (f==at_matrix || f==at_vector || f==at_array){
	g.subtype=_MATRIX__VECT;
	return g;
      }
      return f(g,contextptr);
      // setsizeerr();
    }
    if (f.type==_INT_ && f.val>=0) {
      int i=f.val;
      if (f.val==_POLY1__VECT){ // remove order_size
	vecteur l(lop(g,at_order_size));
	vecteur lp(l.size(),zero);
	g=subst(g,l,lp,false,contextptr);
	return g;
      }
      if (f.subtype==_INT_MAPLECONVERSION){
	switch (i){
	case _TRIG:
	  return sincos(g,contextptr);
	case _EXPLN:
	  return trig2exp(g,contextptr);
	case _PARFRAC: case _FULLPARFRAC:
	  return _partfrac(g,contextptr);
	case _MAPLE_LIST:
	  g.subtype=0;
	  return g;
	default:
	  setsizeerr();
	}
      }
      g.subtype=v.back().val;
      return g;
    }
    setsizeerr();
    return 0;
  }
  const string _convert_s("convert");
  unary_function_eval __convert(&_convert,_convert_s);
  unary_function_ptr at_convert (&__convert,0,true);

  const string _convertir_s("convertir");
  unary_function_eval __convertir(&_convert,_convertir_s);
  unary_function_ptr at_convertir (&__convertir,0,true);

  gen _deuxpoints(const gen & args){
    return symbolic(at_deuxpoints,args);
  }
  const string _deuxpoints_s(":");
  unary_function_unary __deuxpoints(&_deuxpoints,_deuxpoints_s,&printsommetasoperator,&texprintsommetasoperator);
  unary_function_ptr at_deuxpoints (&__deuxpoints);

  // FIXME SECURITY
  gen quote_read(const gen & args,GIAC_CONTEXT){
    if (args.type!=_STRNG)
      return symbolic(at_read,args);
    string fichier=*args._STRNGptr;
    ifstream inf(fichier.c_str());
    if (!inf)
      return undef;
#ifdef VISUALC
    char * thebuf = new char[BUFFER_SIZE];
#else
    char thebuf[BUFFER_SIZE];
#endif
    inf.getline(thebuf,BUFFER_SIZE,'\n');
    string lu(thebuf),thet;
#ifdef VISUALC
    delete [] thebuf;
#endif
    if (lu.size()>9 && lu.substr(0,9)=="{VERSION "){ // Maple Worksheet
      ofstream of("__.map");
      mws_translate(inf,of);
      of.close();
      xcas_mode(contextptr)=1;
      *logptr(contextptr) << "Running maple text translation __.map" << endl;
      fichier="__.map";
    }
    if (lu.size()>6 && lu.substr(0,6)=="**TI92"){ // TI archive
      inf.close();
      xcas_mode(contextptr)=3;
      eval(_unarchive_ti(args,contextptr),1,contextptr);
      return symbolic(at_xcas_mode,3);
    }
    if (lu=="\\START92\\\r"){ // TI text
      ofstream of("__.ti");
      ti_translate(inf,of);
      of.close();
      xcas_mode(contextptr)=3;
      *logptr(contextptr) << "Running TI89 text translation __.ti" << endl;
      fichier="__.ti";
    } // end file of type TI
    inf.close();
    ifstream inf2(fichier.c_str());
    vecteur v,l(list_one_letter__IDNT);
    readargs_from_stream(inf2,v,l,contextptr);
    return v;
  }
  gen _read(const gen & args,GIAC_CONTEXT){
    if (args.type!=_STRNG)
      return symbolic(at_read,args);
    return eval(quote_read(args,contextptr),eval_level(contextptr),contextptr);
  }
  const string _read_s("read");
  unary_function_eval __read(&_read,_read_s);
  unary_function_ptr at_read (&__read);

  gen _write(const gen & args,GIAC_CONTEXT){
    check_secure();
    if (args.type==_VECT){
      vecteur & v=*args._VECTptr;
      if (v.size()<2 || v.front().type!=_STRNG)
	setsizeerr();
      ofstream inf(v[0]._STRNGptr->c_str());
      const_iterateur it=v.begin()+1,itend=v.end();
      for (;it!=itend;++it){
	if (it->type==_IDNT){
	  gen tmp=eval(*it,1,contextptr);
	  gen tmp2=*it;
	  inf << symb_sto(tmp,tmp2) << ";" << endl;
	}
      }
      return plus_one;
    }
    if (args.type!=_STRNG)
      return symbolic(at_write,args);
    ofstream inf(args._STRNGptr->c_str());
    const_iterateur it=history_in(contextptr).begin(),itend=history_in(contextptr).end();
    if (it==itend)
      return zero;
    for (;it!=itend;++it){
      inf << *it << ";" << endl;
    }
    return plus_one;
  }
  const string _write_s("write");
  unary_function_eval __write(&_write,_write_s);
  unary_function_ptr at_write (&__write,_QUOTE_ARGUMENTS,true);

  gen _save_history(const gen & args,GIAC_CONTEXT){
    check_secure();
    if (args.type!=_STRNG)
      return symbolic(at_save_history,args);
    ofstream of(args._STRNGptr->c_str());
    vecteur v(history_in(contextptr));
    if (!v.empty() && v.back().is_symb_of_sommet(at_save_history))
      v.pop_back();
    of << gen(history_in(contextptr),_SEQ__VECT) << endl;
    return plus_one;
  }
  const string _save_history_s("save_history");
  unary_function_eval __save_history(&_save_history,_save_history_s);
  unary_function_ptr at_save_history (&__save_history,0,true);
  /*
  gen _matrix(const gen & args){
    if (!ckmatrix(args))
      return symbolic(at_matrix,args);
    gen res=args;
    res.subtype=_MATRIX__VECT;
    return res;
  }
  const string _matrix_s("matrix");
  unary_function_unary __matrix(&_matrix,_matrix_s);
  unary_function_ptr at_matrix (&__matrix);
  */

#ifdef HAVE_LIBPARI
  string pari_help(const gen & g);
#else
  string pari_help(const gen & g){
    return "Sorry! PARI support not compiled in";
  }
#endif

  gen symb_findhelp(const gen & args){
    return symbolic(at_findhelp,args);
  }
  gen _findhelp(const gen & g,GIAC_CONTEXT){
    gen args(g);
    int lang=language(contextptr);
    int helpitems = 0;
    if (g.type==_VECT && g.subtype==_SEQ__VECT && g._VECTptr->size()==2 && g._VECTptr->back().type==_INT_){
      args=g._VECTptr->front();
      lang=absint(g._VECTptr->back().val);
    }
    if (args.type==_FUNC && args._FUNCptr->ptr->s=="pari")
      return string2gen(pari_help(0),false);
    if (args.type==_SYMB && args._SYMBptr->sommet.ptr->s=="pari")
      return string2gen(pari_help(args._SYMBptr->feuille),false);
    string argss=args.print(contextptr);
    // remove space at the end if required
    while (!argss.empty() && argss[argss.size()-1]==' ')
      argss=argss.substr(0,argss.size()-1);
    if (argss.size()>5 && argss.substr(0,5)=="pari_")
      return string2gen(pari_help(string2gen(argss.substr(5,argss.size()-5),false)),false);      
    if (!vector_aide_ptr || vector_aide_ptr->empty()){
      if (!vector_aide_ptr)
	vector_aide_ptr = new vector<aide>;
      * vector_aide_ptr=readhelp("aide_cas",helpitems,false);
      if (!helpitems){
	* vector_aide_ptr=readhelp(default_helpfile,helpitems);
      }
      if (!helpitems){
	* vector_aide_ptr=readhelp((giac_aide_dir()+"aide_cas").c_str(),helpitems);
      }
    }
    if (vector_aide_ptr){
      string s=argss; // args.print(contextptr);
      int l=s.size();
      if ( (l>2) && (s[0]=='\'') && (s[l-1]=='\'') )
	s=s.substr(1,l-2);
      l=s.size();
      if (l && s[l-1]==')'){
	int i=l-1;
	for (;i;--i){
	  if (s[i]=='(')
	    break;
	}
	if (i)
	  s=s.substr(0,i);
      }
      s=writehelp(helpon(s,*vector_aide_ptr,lang,vector_aide_ptr->size()),lang);
      return string2gen(s,false);
    }
    else
      setsizeerr("No help file found");
    return 0;
  }
  const string _findhelp_s("findhelp");
  unary_function_eval __findhelp(&_findhelp,_findhelp_s);
  unary_function_ptr at_findhelp (&__findhelp,_QUOTE_ARGUMENTS,true);

  gen _member(const gen & args,GIAC_CONTEXT){
    gen g=args;
    vecteur v;
    if (args.type!=_VECT){
      g=args.eval(eval_level(contextptr),contextptr);
      if (g.type!=_VECT)
	return symbolic(at_member,args);
      v=*g._VECTptr;
    }
    else {
      v=*args._VECTptr;
      if (v.size()>1){
	v[0]=eval(v[0],eval_level(contextptr),contextptr);
	v[1]=eval(v[1],eval_level(contextptr),contextptr);
      }
    }
    int s=v.size();
    if (s<2)
      toofewargs("");
    if (v[1].type!=_VECT)
      setsizeerr();
    int i=equalposcomp(*v[1]._VECTptr,v[0]);
    if (s==3){
      if (xcas_mode(contextptr))
	sto(i,v[2],contextptr);
      else
	sto(i-1,v[2],contextptr);
    }
    return i;
  }
  const string _member_s("member");
  unary_function_eval __member(&_member,_member_s);
  unary_function_ptr at_member (&__member,_QUOTE_ARGUMENTS,true);

  // tablefunc(expression,[var,min,step])
  gen _tablefunc(const gen & args,GIAC_CONTEXT){
    gen f,x=vx_var,xstart=gnuplot_xmin,step=(gnuplot_xmax-gnuplot_xmin)/10;
    gen xmax=gnuplot_xmax;
    if (args.type==_VECT){
      vecteur & v=*args._VECTptr;
      int s=v.size();
      if (!s)
	toofewargs("");
      f=v[0];
      if (s>1)
	x=v[1];
      if (s>2)
	xstart=v[2];
      if (s>3)
	step=v[3];
      if (s>4)
	xmax=v[4];
    }
    else
      f=args;
    vecteur l0(makevecteur(x,f));
    gen graphe=symbolic(at_plotfunc,
			gen(makevecteur(_cell(makevecteur(vecteur(1,minus_one),vecteur(1,zero))),
					symb_equal(_cell(makevecteur(vecteur(1,minus_one),vecteur(1,minus_one))),symb_interval(xstart,xmax))
				    ),_SEQ__VECT));
    graphe.subtype=_SPREAD__SYMB;
    vecteur l1(makevecteur(step,graphe));
    gen l31(_cell(makevecteur(vecteur(1,minus_one),vecteur(1,zero)))+_cell(makevecteur(plus_one,vecteur(1,zero))));
    l31.subtype=_SPREAD__SYMB;
    gen l32(symb_evalf(symbolic(at_subst,makevecteur(_cell(makevecteur(zero,vecteur(1,zero))),_cell(makevecteur(zero,vecteur(1,minus_one))),_cell(makevecteur(vecteur(1,zero),vecteur(1,minus_one)))))));
    l32.subtype=_SPREAD__SYMB;
    vecteur l2(makevecteur(xstart,l32));
    vecteur l3(makevecteur(l31,l32));
    return makevecteur(l0,l1,l2,l3);
  }
  const string _tablefunc_s("tablefunc");
  unary_function_eval __tablefunc(&_tablefunc,_tablefunc_s);
  unary_function_ptr at_tablefunc (&__tablefunc,0,true);

  // tableseq(expression,[var,value])
  // var is a vector of dim the number of terms in the recurrence
  gen _tableseq(const gen & args,GIAC_CONTEXT){
    gen f,x=vx_var,uzero=zero;
    int dim=1;
    double xmin=gnuplot_xmin,xmax=gnuplot_xmax;
    if (args.type==_VECT){
      vecteur & v=*args._VECTptr;
      int s=v.size();
      if (!s)
	toofewargs("");
      f=v[0];
      if (s>1)
	x=v[1];
      if (s>2)
	uzero=evalf_double(v[2],1,contextptr);
      if (x.type==_VECT){
	dim=x._VECTptr->size();
	if (uzero.type!=_VECT)
	  settypeerr();
	if (uzero._VECTptr->front().type==_VECT)
	  uzero=uzero._VECTptr->front();
	if ( (uzero.type!=_VECT) || (signed(uzero._VECTptr->size())!=dim) )
	  setdimerr("");
      }
      else {
	if (uzero.type==_VECT && uzero._VECTptr->size()==3){
	  vecteur & uv=*uzero._VECTptr;
	  if (uv[1].type!=_DOUBLE_ || uv[2].type!=_DOUBLE_)
	    setsizeerr();
	  xmin=uv[1]._DOUBLE_val;
	  xmax=uv[2]._DOUBLE_val;
	  uzero=uv[0];
	}
      }
    }
    else
      f=args;
    vecteur res;
    res.push_back(f);
    if (x.type!=_VECT){
      res.push_back(x);
      if (dim!=1)
	res.push_back(dim);
      else {
	gen l31(symbolic(at_plotseq,
		       gen(makevecteur(
				       _cell(makevecteur(zero,vecteur(1,zero))),
				       symb_equal(_cell(makevecteur(plus_one,vecteur(1,zero))),makevecteur(_cell(makevecteur(vecteur(1,plus_one),vecteur(1,zero))),xmin,xmax)),
				       10),_SEQ__VECT
			   )
		       )
	      );
	l31.subtype=_SPREAD__SYMB;
	res.push_back(l31);
      }
      res.push_back(uzero);
      gen l51(symb_evalf(symbolic(at_subst,makevecteur(_cell(makevecteur(zero,vecteur(1,zero))),_cell(makevecteur(plus_one,vecteur(1,zero))),_cell(makevecteur(vecteur(1,minus_one),vecteur(1,zero)))))));
      l51.subtype=_SPREAD__SYMB;
      res.push_back(l51);
    }
    else {
      for (int i=0;i<dim;++i)
	res.push_back(x[i]);
      vecteur tmp1,tmp2;
      for (int i=0;i<dim;++i){
	res.push_back(uzero[i]);
	tmp1.push_back(_cell(makevecteur(i+1,vecteur(1,zero))));
	tmp2.push_back(_cell(makevecteur(vecteur(1,i-dim),vecteur(1,zero))));
      }
      gen l41(symb_eval(symbolic(at_subst,makevecteur(_cell(makevecteur(zero,vecteur(1,zero))),tmp1,tmp2))));
      l41.subtype=_SPREAD__SYMB;
      res.push_back(l41);
    }
    return mtran(vecteur(1,res));
  }
  const string _tableseq_s("tableseq");
  unary_function_eval __tableseq(&_tableseq,_tableseq_s);
  unary_function_ptr at_tableseq (&__tableseq,_QUOTE_ARGUMENTS,true);

  gen protecteval(const gen & g,int level, GIAC_CONTEXT){
    gen res;
    ctrl_c = false;
    // save cas_setup in case of an exception
    vecteur cas_setup_save = cas_setup(contextptr);
    try {
      res=approx_mode(contextptr)?g.evalf(level,contextptr):g.eval(level,contextptr);
    }
    catch (std::runtime_error & e){
      res=string2gen(e.what(),false);
      ctrl_c=false;
      // something went wrong, so restore the old cas_setup
      cas_setup(cas_setup_save, contextptr);
    }
    return res;
  }

  string printasnodisp(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    int maplemode=xcas_mode(contextptr) & 0x07;
    if (maplemode==1 || maplemode==2){
      string res=feuille.print(contextptr);
      int l=res.size(),j;
      for (j=l-1;j>=0 && res[j]==' ';--j)
	;
      if (res[j]==';')
	res[j]=':';
      else
	res += ':';
      return res;
    }
    return sommetstr+"("+feuille.print(contextptr)+")";
  }
  gen _nodisp(const gen & args){
    return string2gen("Done",false);
  }
  const string _nodisp_s("nodisp");
  unary_function_unary __nodisp(&_nodisp,_nodisp_s,&printasnodisp);
  unary_function_ptr at_nodisp (&__nodisp,0,true);

  gen _unapply(const gen & args,GIAC_CONTEXT){
    if ( (args.type!=_VECT) || args._VECTptr->empty() )
      settypeerr();
    vecteur v=*args._VECTptr,w;
    int s=v.size();
    if (s<2)
      w=vecteur(1,vx_var);
    else {
      if (s==2 && v[1].type==_VECT)
	w=*v[1]._VECTptr;
      else
	w=vecteur(v.begin()+1,v.end());
    }
    gen g=subst(v[0].eval(eval_level(contextptr),contextptr),w,w,false,contextptr);
    if (v[0].type!=_VECT && g.type==_VECT && !g.subtype)
      g=makevecteur(g);
    return symbolic(at_program,makevecteur(gen(w,_SEQ__VECT),w*zero,g));
  }
  const string _unapply_s("unapply");
  unary_function_eval __unapply(&_unapply,_unapply_s);
  unary_function_ptr at_unapply (&__unapply,_QUOTE_ARGUMENTS,true);

  gen _makevector(const gen & args){
    if (args.type!=_VECT)
      return vecteur(1,args);
    vecteur & v=*args._VECTptr;
    if (ckmatrix(args))
      return gen(v,_MATRIX__VECT);
    return v;
  }
  const string _makevector_s("makevector");
  unary_function_unary __makevector(&_makevector,_makevector_s);
  unary_function_ptr at_makevector (&__makevector,0,true);


  gen _makesuite(const gen & args){
    if (args.type!=_VECT)
      return vecteur(1,args);
    vecteur & v=*args._VECTptr;
    return gen(v,_SEQ__VECT);
  }
  const string _makesuite_s("makesuite");
  unary_function_unary __makesuite(&_makesuite,_makesuite_s);
  unary_function_ptr at_makesuite (&__makesuite,0,true);

  gen _matrix(const gen & g,const context * contextptr){
    if (g.type!=_VECT)
      settypeerr();
    vecteur v=*g._VECTptr;
    if (ckmatrix(v))
      return gen(v,_MATRIX__VECT);
    int vs=v.size();
    if (vs<2)
      settypeerr();
    if (vs==2){
      v.push_back(zero);
      ++vs;
    }
    if ( (v[0].type!=_INT_) || (v[1].type!=_INT_) )
      setsizeerr();
    int l(v[0].val),c(v[1].val);
    bool transpose=(vs>3);
    if (transpose){ // try to merge arguments there
      // v[2]..v[vs-1] represents flattened submatrices 
      vecteur v2;
      for (int i=2;i<vs;++i){
	if (v[i].type!=_VECT)
	  settypeerr();
	vecteur & w = *v[i]._VECTptr;
	int vis=w.size();
	if (vis % l)
	  setdimerr();
	int nc=vis/l;
	for (int J=0;J<nc;++J){
	  for (int I=J;I<vis;I+=nc)
	    v2.push_back(w[I]);
	}
      }
      v[2]=v2;
      std::swap<int>(l,c);
    }
    if (v[2].type==_VECT){
      vecteur w=*v[2]._VECTptr;
      int s=w.size();
      if (ckmatrix(w)){
	int ss=0;
	if (s)
	  ss=w[0]._VECTptr->size();
	int ll=min(l,s);
	for (int i=0;i<ll;++i){
	  if (ss<c)
	    w[i]=mergevecteur(*w[i]._VECTptr,vecteur(c-ss));
	  else
	    w[i]=vecteur(w[i]._VECTptr->begin(),w[i]._VECTptr->begin()+c);
	}
	if (s<l)
	  w=mergevecteur(w,vecteur(l-s,vecteur(c)));
	else
	  w=vecteur(w.begin(),w.begin()+l);
	return gen(makefreematrice(w),_MATRIX__VECT);
      }
      else {
	vecteur res;
	if (s<l*c)
	  w=mergevecteur(w,vecteur(l*c-s));
	for (int i=0;i<l;++i)
	  res.push_back(vecteur(w.begin()+i*c,w.begin()+(i+1)*c));
	if (transpose)
	  res=mtran(res);
	return gen(makefreematrice(res),_MATRIX__VECT);
      }
    }
    // v[2] as a function, should take 2 args
    gen f=v[2];
    if (!f.is_symb_of_sommet(at_program))
      return gen(vecteur(l,vecteur(c,f)),_MATRIX__VECT);
    vecteur res(l);
    int decal=(xcas_mode(contextptr)!=0);
    for (int i=decal;i<l+decal;++i){
      vecteur tmp(c);
      for (int j=decal;j<c+decal;++j)
	tmp[j-decal]=f(gen(makevecteur(i,j),_SEQ__VECT),contextptr);
      res[i-decal]=tmp;
    }
    return gen(res,_MATRIX__VECT);
  }
  const string _matrix_s("matrix");
  unary_function_eval __matrix(&_matrix,_matrix_s);
  unary_function_ptr at_matrix (&__matrix,0,true);

  string printasbreak(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if (xcas_mode(contextptr)==3)
      return "Exit ";
    else
      return sommetstr;
  }
  gen _break(const gen & args){
    return symbolic(at_break,0);
  }
  const string _break_s("break");
  unary_function_unary __break(&_break,_break_s,&printasbreak);
  unary_function_ptr at_break (&__break);

  string printascontinue(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if (xcas_mode(contextptr)==3)
      return "Cycle ";
    else
      return sommetstr;
  }
  gen _continue(const gen & args){
    return symbolic(at_continue,0);
  }
  const string _continue_s("continue");
  unary_function_unary __continue(&_continue,_continue_s,&printascontinue);
  unary_function_ptr at_continue (&__continue);

  string printaslabel(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if (xcas_mode(contextptr)==3)
      return "Lbl "+feuille.print(contextptr);
    else
      return "label "+feuille.print(contextptr);
  }
  gen _label(const gen & args){
    return symbolic(at_label,args);
  }
  const string _label_s("label");
  unary_function_unary __label(&_label,_label_s,&printaslabel);
  unary_function_ptr at_label (&__label);

  string printasgoto(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if (xcas_mode(contextptr)==3)
      return "Goto "+feuille.print(contextptr);
    else
      return "goto "+feuille.print(contextptr);
  }
  gen _goto(const gen & args){
    return symbolic(at_goto,args);
  }
  const string _goto_s("goto");
  unary_function_unary __goto(&_goto,_goto_s,&printasgoto);
  unary_function_ptr at_goto (&__goto);

  vecteur local_vars(const vecteur & v,GIAC_CONTEXT){
    const_iterateur it=v.begin(),itend=v.end();
    vecteur res;
    for (;it!=itend;++it){
      if (it->type==_IDNT && 
	  (contextptr?contextptr->tabptr->find(*it->_IDNTptr->name)==contextptr->tabptr->end():it->_IDNTptr->localvalue->empty())
	  )
	res.push_back(*it);
    }
    return res;
  }
  gen _tilocal(const gen & args,const context * contextptr){
    if (args.type!=_VECT || args._VECTptr->size()!=2)
      return symbolic(at_tilocal,args);
    vecteur & v=*args._VECTptr;
    // find local variables
    vecteur cond(gen2vecteur(v[1]));
    vecteur docond,vars;
    const_iterateur it=cond.begin(),itend=cond.end();
    for (;it!=itend;++it){
      if (it->type!=_SYMB)
	continue;
      unary_function_ptr & u=it->_SYMBptr->sommet;
      gen & g=it->_SYMBptr->feuille;
      if ( (g.type!=_VECT) || (g._VECTptr->empty()) )
	setsizeerr();
      if (u==at_equal){
	gen tmp=g._VECTptr->front();
	if (tmp.type==_IDNT){
	  gen tmp1(eval(tmp,eval_level(contextptr),contextptr));
	  if (tmp1.type==_IDNT)
	    tmp=tmp1;
	  tmp.subtype=0; // otherwise if inside a folder sto will affect tmp!
	  vars.push_back(tmp);
	}
	docond.push_back(symbolic(at_sto,makevecteur(g._VECTptr->back(),tmp)));
	continue;
      }
      if (u==at_sto){
	if (g._VECTptr->back().type==_IDNT)
	  vars.push_back(g._VECTptr->back());
	docond.push_back(*it);
	continue;
      }
      if (g._VECTptr->front().type==_IDNT)
	vars.push_back(g._VECTptr->front());
      docond.push_back(symbolic(at_assume,*it));
    }
    vecteur v0(vars.size(),zero);
    gen gv(v[0]);
    // Replace v[0] by its value if it is a global identifier
    if (gv.type==_IDNT){
      if (contextptr){
	sym_tab::const_iterator it=contextptr->tabptr->find(*gv._IDNTptr->name),itend=contextptr->tabptr->end();
	if (it!=itend)
	  gv=it->second;
      }
      else {
	if (gv._IDNTptr->value)
	  gv=*gv._IDNTptr->value;
      }
    }
    // Replace local variables by their value in gv
    vecteur vname(local_vars(*_lname(gv)._VECTptr,contextptr));
    vecteur vval(vname);
    iterateur jt=vval.begin(),jtend=vval.end();
    for (;jt!=jtend;++jt){
      *jt=contextptr?jt->_IDNTptr->eval(1,*jt,contextptr):jt->_IDNTptr->localvalue->back();
    }
    gv=quotesubst(gv,vname,vval,contextptr);
    // Replace vars global IDNT by local IDNT
    vname=vars;
    jt=vname.begin(),jtend=vname.end();
    for (;jt!=jtend;++jt)
      jt->subtype=_GLOBAL__EVAL;
    vval=vars;
    jt=vval.begin(),jtend=vval.end();
    for (;jt!=jtend;++jt)
      jt->subtype=0;
    gv=quotesubst(gv,vname,vval,contextptr);
    docond=*quotesubst(docond,vname,vval,contextptr)._VECTptr;
    gen prg=symb_program(gen(vname,_SEQ__VECT),gen(v0,_SEQ__VECT),symb_bloc(makevecteur(docond,gv)),contextptr);
    return prg(v0,contextptr);
  }
  const string _tilocal_s("|");
  unary_function_eval __tilocal(&_tilocal,_tilocal_s,&printsommetasoperator,&texprintsommetasoperator);
  unary_function_ptr at_tilocal (&__tilocal,_QUOTE_ARGUMENTS);

  string printasdialog(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    return "Dialog "+symbolic(at_bloc,feuille).print(contextptr)+indent(contextptr)+"EndDialog";
  }  
  string printasinputform(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if (xcas_mode(contextptr)==3)
      return printasdialog(feuille,sommetstr,contextptr);
    return sommetstr+"("+feuille.print(contextptr)+")";
  }  

  // Eval everything except IDNT and symbolics with
  vecteur inputform_pre_analysis(const gen & g,GIAC_CONTEXT){
    vecteur v(gen2vecteur(g));
    int s=v.size();
    for (int i=0;i<s;++i){
      if (v[i].type==_IDNT || v[i].type!=_SYMB)
	continue;
      unary_function_ptr & u =v[i]._SYMBptr->sommet;
      if ( (u==at_output) || (u==at_Text) || (u==at_Title) || (u==at_click) || (u==at_Request) || (u==at_choosebox) || (u==at_DropDown) || (u==at_Popup) )
	continue;
      v[i]=protecteval(v[i],eval_level(contextptr),contextptr);
    }
    return v;
  }
  gen inputform_post_analysis(const vecteur & v,const gen & res,GIAC_CONTEXT){
    return res.eval(eval_level(contextptr),contextptr);
  }
  // user input sent back to the parent process
  gen _inputform(const gen & args,GIAC_CONTEXT){
    string cs("inputform may be used in a window environment only");
#ifdef WIN32
    *logptr(contextptr) << cs << endl;
    return string2gen(cs,false);
#endif
    if (child_id){ 
      *logptr(contextptr) << cs << endl;
      return string2gen(cs,false);
    }
    // pre-analysis
    vecteur v(gen2vecteur(args));
    // int vs=signed(v.size());
    gen res;
    // form
    ofstream child_out(cas_sortie_name().c_str());
    gen e(symbolic(at_inputform,args));
    *logptr(contextptr) << "Archiving " << e << endl;
    archive(child_out,e,contextptr);
    archive(child_out,e,contextptr);
    if ( (args.type==_VECT) && (args._VECTptr->empty()) )
      child_out << "User input requested\n" << '¤' ;
    else
      child_out << args << '¤' ;
    child_out.close();
    kill_and_wait_sigusr2();
    ifstream child_in(cas_entree_name().c_str());
    res= unarchive(child_in,list_one_letter__IDNT,contextptr);
    child_in.close();
    *logptr(contextptr) << "Inputform reads " << res << endl;
    // post analysis
    return inputform_post_analysis(v,res,contextptr);
  }
  const string _inputform_s("inputform");
  unary_function_eval __inputform(&giac::_inputform,_inputform_s,&printasinputform);
  unary_function_ptr at_inputform (&__inputform,_QUOTE_ARGUMENTS,true);

  gen _choosebox(const gen & args,GIAC_CONTEXT){
    return __inputform.op(symbolic(at_choosebox,args),contextptr);
  }
  const string _choosebox_s("choosebox");
  unary_function_eval __choosebox(&giac::_choosebox,_choosebox_s);
  unary_function_ptr at_choosebox (&__choosebox,_QUOTE_ARGUMENTS,true);

  gen _output(const gen & args,GIAC_CONTEXT){
    return __inputform.op(symbolic(at_output,args),contextptr);
  }
  const string _output_s("output");
  unary_function_eval __output(&giac::_output,_output_s);
  unary_function_ptr at_output (&__output,_QUOTE_ARGUMENTS,true);

  gen _input(const gen & args,bool textinput,GIAC_CONTEXT){
    vecteur v(gen2vecteur(args));
    const_iterateur it=v.begin(),itend=v.end();
    if (it==itend)
      return __click.op(args,contextptr);
    gen res;
    for (;it!=itend;++it){
      if (it->type==_IDNT){
	if (textinput)
	  res=__click.op(makevecteur(string2gen(it->print(contextptr)),0,*it,1),contextptr);
	else
	  res=__click.op(makevecteur(string2gen(it->print(contextptr),false),0,*it),contextptr);
      }
      if (it+1==itend)
	break;
      if (it->type==_STRNG && (it+1)->type==_IDNT){
	if (textinput)
	  res=__click.op(makevecteur(*it,0,*(it+1),1),contextptr);
	else
	  res=__click.op(makevecteur(*it,0,*(it+1)),contextptr);
	++it;
      }
    }
    return res;
  }

  string printastifunction(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if (feuille.type==_VECT && feuille.subtype==_SEQ__VECT && feuille._VECTptr->empty())
      return sommetstr+" ";
    return sommetstr+" "+feuille.print(contextptr);
  }
  gen _Text(const gen & args,GIAC_CONTEXT){
    return __inputform.op(symbolic(at_Text,args),contextptr);
  }
  const string _Text_s("Text");
  unary_function_eval __Text(&giac::_Text,_Text_s,&printastifunction);
  unary_function_ptr at_Text (&__Text,_QUOTE_ARGUMENTS);

  gen _Title(const gen & args,GIAC_CONTEXT){
    return __inputform.op(symbolic(at_Title,args),contextptr);
  }
  const string _Title_s("Title");
  unary_function_eval __Title(&giac::_Title,_Title_s,&printastifunction);
  unary_function_ptr at_Title (&__Title,_QUOTE_ARGUMENTS);

  gen _Request(const gen & args,GIAC_CONTEXT){
    return __inputform.op(symbolic(at_Request,args),contextptr);
  }
  const string _Request_s("Request");
  unary_function_eval __Request(&giac::_Request,_Request_s,&printastifunction);
  unary_function_ptr at_Request (&__Request,_QUOTE_ARGUMENTS);

  gen _DropDown(const gen & args,GIAC_CONTEXT){
    return __inputform.op(symbolic(at_DropDown,args),contextptr);
  }
  const string _DropDown_s("DropDown");
  unary_function_eval __DropDown(&giac::_DropDown,_DropDown_s,&printastifunction);
  unary_function_ptr at_DropDown (&__DropDown,_QUOTE_ARGUMENTS);

  gen _Popup(const gen & args,GIAC_CONTEXT){
    return __inputform.op(symbolic(at_Popup,args),contextptr);
  }
  const string _Popup_s("Popup");
  unary_function_eval __Popup(&giac::_Popup,_Popup_s,&printastifunction);
  unary_function_ptr at_Popup (&__Popup,_QUOTE_ARGUMENTS);

  gen _Dialog(const gen & args,GIAC_CONTEXT){
    return __inputform.op(args,contextptr);
  }
  const string _Dialog_s("Dialog");
  unary_function_eval __Dialog(&giac::_Dialog,_Dialog_s,&printasdialog);
  unary_function_ptr at_Dialog (&__Dialog,_QUOTE_ARGUMENTS);

  gen _expr(const gen & args,GIAC_CONTEXT){
    if (args.type==_VECT && args._VECTptr->size()==2 && args._VECTptr->front().type==_STRNG && args._VECTptr->back().type==_INT_){
      int mode=args._VECTptr->back().val;
      bool rpnmode=mode<0;
      mode=absint(mode) % 256;
      if (mode>3)
	setsizeerr();
      int save_mode=xcas_mode(contextptr);
      bool save_rpnmode=rpn_mode;
      xcas_mode(contextptr)=mode;
      rpn_mode=rpnmode;
      gen res=eval(gen(*args._VECTptr->front()._STRNGptr,contextptr),eval_level(contextptr),contextptr);
      xcas_mode(contextptr)=save_mode;
      rpn_mode=save_rpnmode;
      return res;
    }
    if (args.type!=_STRNG)
      return symbolic(at_expr,args);
    return eval(gen(*args._STRNGptr,contextptr),eval_level(contextptr),contextptr);
  }
  const string _expr_s("expr");
  unary_function_eval __expr(&giac::_expr,_expr_s);
  unary_function_ptr at_expr (&__expr,0,true);

  const string _execute_s("execute");
  unary_function_eval __execute(&giac::_expr,_execute_s);
  unary_function_ptr at_execute (&__execute,0,true);

  gen _string(const gen & args,GIAC_CONTEXT){
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
    else
      res=args.print(contextptr);
    return string2gen(res,false);
  }
  const string _string_s("string");
  unary_function_eval __string(&giac::_string,_string_s);
  unary_function_ptr at_string (&__string,0,true);

  gen _part(const gen & args,GIAC_CONTEXT){
    if ( (args.type==_VECT) && args._VECTptr->size()==2 ){
      gen & i=args._VECTptr->back();
      gen & g=args._VECTptr->front();
      if (i.type!=_INT_ || i.val<=0){
	if (g.type!=_SYMB)
	  return string2gen(g.print(contextptr),false);
	else
	  return string2gen(g._SYMBptr->sommet.ptr->s,false);
      }
      else {
	if (g.type!=_SYMB){
	  if (i.val!=1)
	    setsizeerr();
	  return g;
	}
	else {
	  vecteur v(gen2vecteur(g._SYMBptr->feuille));
	  if (signed(v.size())<i.val)
	    setsizeerr();
	  return v[i.val-1];
	}
      }
    }
    if (args.type==_SYMB)
      return gen2vecteur(args._SYMBptr->feuille).size();
    return 0;
  }
  const string _part_s("part");
  unary_function_eval __part(&giac::_part,_part_s);
  unary_function_ptr at_part (&__part,0,true);

  string tiasc_translate(const string & s){
    int l=s.size();
    string t("");
    for (int i=0;i<l;++i){
      char c=s[i];
      if (c=='\r')
	continue;
      if (c=='@'){
	t += "//";
	continue;
      }
      if (c=='\\'){
	++i;
	string ti_escape("");
	for (;i<l;++i){
	  char c=s[i];
	  if (c=='\\' || c==' '){
	    break;
	  }
	  ti_escape += c;
	}
	if (i==l || c==' ')
	  return t+"::"+ti_escape;
	if (ti_escape=="->"){
	  t += "=>";
	  continue;
	}
	if (ti_escape=="(C)"){ // comment
	  t += "//";
	  continue;
	}
	if (ti_escape=="(-)"){
	  t += '-';
	  continue;
	}
	if (ti_escape=="e"){
	  t += "exp(1)";
	  continue;
	}
	if (ti_escape=="i"){
	  t += '\xa1';
	  continue;
	}
	t += ti_escape;
      }
      else
	t += c;
    }
    if (t==string(t.size(),' '))
      return "";
    return t;
  }

  gen _Pause(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT || !g._VECTptr->empty())
      __interactive.op(symbolic(at_print,g),contextptr);
    __interactive.op(symbolic(at_Pause,0),contextptr);
    return 0;
  }
  const string _Pause_s("Pause");
  unary_function_eval __Pause(&_Pause,_Pause_s,&printastifunction);
  unary_function_ptr at_Pause (&__Pause);

  const string _DelVar_s("DelVar");
  unary_function_eval __DelVar(&_purge,_DelVar_s,&printastifunction);
  unary_function_ptr at_DelVar (&__DelVar,_QUOTE_ARGUMENTS);

  gen _Row(const gen & g,GIAC_CONTEXT){
    return spread_Row(contextptr);
  }
  const string _Row_s("Row");
  unary_function_eval __Row(&_Row,_Row_s,&printastifunction);
  unary_function_ptr at_Row (&__Row);

  gen _Col(const gen & g,GIAC_CONTEXT){
    return spread_Col(contextptr);
  }
  const string _Col_s("Col");
  unary_function_eval __Col(&_Col,_Col_s,&printastifunction);
  unary_function_ptr at_Col (&__Col);

  gen matrix_apply(const gen & a,const gen & b,gen (* f) (const gen &, const gen &) ){
    if (a.type!=_VECT || b.type!=_VECT || a._VECTptr->size()!=b._VECTptr->size())
      return apply(a,b,f);
    const_iterateur it=a._VECTptr->begin(),itend=a._VECTptr->end(),jt=b._VECTptr->begin();
    vecteur res;
    res.reserve(itend-it);
    for (;it!=itend;++it,++jt){
      res.push_back(apply(*it,*jt,f));
    }
    return gen(res,a.subtype);
  }
  gen matrix_apply(const gen & a,const gen & b,GIAC_CONTEXT,gen (* f) (const gen &, const gen &,GIAC_CONTEXT) ){
    if (a.type!=_VECT || b.type!=_VECT || a._VECTptr->size()!=b._VECTptr->size())
      return apply(a,b,contextptr,f);
    const_iterateur it=a._VECTptr->begin(),itend=a._VECTptr->end(),jt=b._VECTptr->begin();
    vecteur res;
    res.reserve(itend-it);
    for (;it!=itend;++it,++jt){
      res.push_back(apply(*it,*jt,contextptr,f));
    }
    return gen(res,a.subtype);
  }
  gen prod(const gen & a,const gen &b){
    return a*b;
  }
  gen somme(const gen & a,const gen &b){
    return a+b;
  }
  gen _pointprod(const gen & g){
    gen a,b;
    check_binary(g,a,b);
    return matrix_apply(a,b,prod);
  }
  const string _pointprod_s(".*");
  unary_function_unary __pointprod(&_pointprod,_pointprod_s,&printsommetasoperator,&texprintsommetasoperator);
  unary_function_ptr at_pointprod (&__pointprod);

  gen _pointdivision(const gen & g){
    gen a,b;
    check_binary(g,a,b);
    return matrix_apply(a,b,rdiv);
  }
  const string _pointdivision_s("./");
  unary_function_unary __pointdivision(&_pointdivision,_pointdivision_s,&printsommetasoperator,&texprintsommetasoperator);
  unary_function_ptr at_pointdivision (&__pointdivision);

  gen giac_pow(const gen &,const gen &,GIAC_CONTEXT);
  gen _pointpow(const gen & g,GIAC_CONTEXT){
    gen a,b;
    check_binary(g,a,b);
    return matrix_apply(a,b,contextptr,giac_pow);
  }
  const string _pointpow_s(".^");
  unary_function_eval __pointpow(&_pointpow,_pointpow_s,&printsommetasoperator,&texprintsommetasoperator);
  unary_function_ptr at_pointpow (&__pointpow);

  string printassuffix(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    return feuille.print(contextptr)+sommetstr;
  }  
  gen _pourcent(const gen & g){
    return rdiv(g,100);
  }
  const string _pourcent_s("%");
  unary_function_unary __pourcent(&_pourcent,_pourcent_s,&printassuffix);
  unary_function_ptr at_pourcent (&__pourcent);

  gen _hash(const gen & g,GIAC_CONTEXT){
    if (g.type!=_STRNG)
      return g;
    return gen(*g._STRNGptr,contextptr);
  }
  const string _hash_s("#");
  unary_function_eval __hash(&_hash,_hash_s);
  unary_function_ptr at_hash (&__hash);

  bool user_screen=false;
  int user_screen_io_x=0,user_screen_io_y=0;
  int user_screen_fontsize=14;
  gen _interactive(const gen & args,GIAC_CONTEXT){
#ifndef WIN32
    if (child_id)
#endif
      { /* *logptr(contextptr) << "Must be under a Unix GUI" << endl; */ return 0; }
#ifndef WIN32
    ofstream child_out(cas_sortie_name().c_str());
    gen e(symbolic(at_interactive,args));
    // *logptr(contextptr) << e << endl;
    archive(child_out,e,contextptr);
    archive(child_out,e,contextptr);
    child_out.close();
    kill_and_wait_sigusr2();
    ifstream child_in(cas_entree_name().c_str());
    gen res= unarchive(child_in,list_one_letter__IDNT,contextptr);
    child_in.close();
    return res;
#endif
  }
  const string _interactive_s("interactive");
  unary_function_eval __interactive(&_interactive,_interactive_s);
  unary_function_ptr at_interactive (&__interactive,_QUOTE_ARGUMENTS,true);

  // v=[ [idnt,value] ... ]
  // search g in v if found return value
  // else return g unevaluated
  gen find_in_folder(vecteur & v,const gen & g){
    const_iterateur it=v.begin(),itend=v.end();
    for (;it!=itend;++it){
      if (it->type!=_VECT || it->_VECTptr->size()!=2)
	continue;
      vecteur & w=*it->_VECTptr;
      if (w[0]==g)
	return w[1];
    }
    return g;
  }

  gen _ti_semi(const gen & args){
    if (args.type!=_VECT || args._VECTptr->size()!=2)
      return symbolic(at_ti_semi,args);
    vecteur & v=*args._VECTptr;
    matrice m1,m2;
    if (!ckmatrix(v[0])){
      if (v[0].type==_VECT)
	m1=vecteur(1,*v[0]._VECTptr);
      else
	m1=vecteur(1,vecteur(1,v[0]));
    }
    else
      m1=*v[0]._VECTptr;
    if (!ckmatrix(v[1])){
      if (v[1].type==_VECT)
	m2=vecteur(1,*v[1]._VECTptr);
      else
	m2=vecteur(1,vecteur(1,v[1]));
    }
    else
      m2=*v[1]._VECTptr;
    // *logptr(contextptr) << m1 << " " << m2 << endl;
    return mergevecteur(m1,m2); 
  }
  const string _ti_semi_s(";");
  unary_function_unary __ti_semi(&_ti_semi,_ti_semi_s,&printsommetasoperator,&texprintsommetasoperator);
  unary_function_ptr at_ti_semi (&__ti_semi);

  gen widget_size(const gen & g){
    return zero;
  }
  const string _widget_size_s("widget_size");
  unary_function_unary __widget_size(&widget_size,_widget_size_s);
  unary_function_ptr at_widget_size (&__widget_size,0,true);

  gen keyboard(const gen & g){
    return zero;
  }
  const string _keyboard_s("keyboard");
  unary_function_unary __keyboard(&keyboard,_keyboard_s);
  unary_function_ptr at_keyboard (&__keyboard,0,true);

  gen current_sheet(const gen & g,GIAC_CONTEXT){
    return zero;
  }
  const string _current_sheet_s("current_sheet");
  unary_function_eval __current_sheet(&current_sheet,_current_sheet_s);
  unary_function_ptr at_current_sheet (&__current_sheet,_QUOTE_ARGUMENTS,true);
  
  string printasmaple_lib(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if (feuille.type!=_VECT || feuille._VECTptr->size()!=2)
      setsizeerr();
    vecteur & v=*feuille._VECTptr;
    return v[0].print(contextptr)+"["+v[1].print(contextptr)+"]";
  }
  gen maple_lib(const gen & g){
    if (g.type!=_VECT || g._VECTptr->size()!=2)
      setsizeerr();
    vecteur & v=*g._VECTptr;
    return v[1];
  }
  const string _maple_lib_s("maple_lib");
  unary_function_unary __maple_lib(&maple_lib,_maple_lib_s,&printasmaple_lib);
  unary_function_ptr at_maple_lib (&__maple_lib);

  gen window_switch(const gen & g){
    return zero; // defined by GUI handler
  }
  const string _window_switch_s("window_switch");
  unary_function_unary __window_switch(&window_switch,_window_switch_s);
  unary_function_ptr at_window_switch (&__window_switch);


  // Funcion has been changed -> simplify
  const string _simplifier_s("simplifier");
  unary_function_eval __simplifier(&_simplify,_simplifier_s);
  unary_function_ptr at_simplifier (&__simplifier,0,true);

  const string _regrouper_s("regrouper");
  unary_function_eval __regrouper(&_simplifier,_regrouper_s);
  unary_function_ptr at_regrouper (&__regrouper,0,true);

  gen find_or_make_symbol(const string & s,GIAC_CONTEXT){
    gen tmp;
    find_or_make_symbol(s,tmp,contextptr);
    return tmp;
  }

  // To each unit we associate a number and a vector of powers of kg, m, s
  std::map<std::string,gen> unit_conversion_map;
  gen mksa_register(const string & s,const gen & equiv){
    unit_conversion_map[s]=equiv;
    return find_or_make_symbol("_"+s,context0);
  }
  // fundemental metric units
  gen _m_unit(mksa_register("m",makevecteur(1,1,0,0,0)));
  gen _kg_unit(mksa_register("kg",makevecteur(1,0,1,0,0)));
  gen _s_unit(mksa_register("s",makevecteur(1,0,0,1,0)));
  gen _A_unit(mksa_register("A",makevecteur(1,0,0,0,1)));
  gen _K_unit(mksa_register("K",makevecteur(1,0,0,0,0,1))); // Kelvin
  gen _mol_unit(mksa_register("mol",makevecteur(1,0,0,0,0,0,1))); // mol
  gen _molK_unit(mksa_register("molK",makevecteur(1,0,0,0,0,1,1))); 
  gen _cd_unit(mksa_register("cd",makevecteur(1,0,0,0,0,0,0,1))); // candela
  gen _E_unit(mksa_register("E",makevecteur(1,0,0,0,0,0,0,0,1))); // euro
  // other metric units in m,kg,s,A
  gen _Bq_unit(mksa_register("Bq",makevecteur(1,0,0,-1,0)));
  gen _C_unit(mksa_register("C",makevecteur(1,0,0,1,1)));
  gen _F_unit(mksa_register("F",makevecteur(1,-2,-1,4,2)));
  gen _Gy_unit(mksa_register("Gy",makevecteur(1,2,0,-2,0)));
  gen _H_unit(mksa_register("H",makevecteur(1,2,1,-2,-2)));
  gen _Hz_unit(mksa_register("Hz",makevecteur(1,0,0,-1,0)));
  gen _J_unit(mksa_register("J",makevecteur(1,2,1,-2,0)));
  gen _mho_unit(mksa_register("mho",makevecteur(1,-2,-1,3,2)));
  gen _N_unit(mksa_register("N",makevecteur(1,1,1,-2,0)));
  gen _Ohm_unit(mksa_register("Ohm",makevecteur(1,2,1,-3,-2)));
  gen _Pa_unit(mksa_register("Pa",makevecteur(1,-1,1,-2,0)));
  gen _rad_unit(mksa_register("rad",makevecteur(1,0,0,0,0))); // radian
  gen _S_unit(mksa_register("S",makevecteur(1,-2,-1,3,2)));
  gen _st_unit(mksa_register("st",makevecteur(1,3,0,0,0)));
  gen _Sv_unit(mksa_register("Sv",makevecteur(1,2,0,-2,0)));
  gen _T_unit(mksa_register("T",makevecteur(1,0,1,-2,-1)));
  gen _V_unit(mksa_register("V",makevecteur(1,2,1,-3,-1)));
  gen _W_unit(mksa_register("W",makevecteur(1,2,1,-3,0)));
  gen _Wb_unit(mksa_register("Wb",makevecteur(1,2,1,-2,-1)));
  vecteur usual_units(mergevecteur(
				   mergevecteur(makevecteur(_Bq_unit,_C_unit,_F_unit,_Gy_unit,_H_unit,_Hz_unit,_J_unit,_mho_unit),
						makevecteur(_N_unit,_Ohm_unit,_Pa_unit,_rad_unit,_S_unit,_Sv_unit,_T_unit)),
				   makevecteur(_V_unit,_W_unit,_Wb_unit)
				  )
		      );
  // useful non metric units
  gen _a_unit(mksa_register("a",makevecteur(100,2,0,0,0)));
  gen _acre_unit(mksa_register("acre",makevecteur(4046.87260987,2,0,0,0)));
  gen _arcmin_unit(mksa_register("arcmin",makevecteur(2.90888208666)));
  gen _arcs_unit(mksa_register("arcs",makevecteur(4.8481368111)));
  gen _atm_unit(mksa_register("atm",makevecteur(101325.0,-1,1,0,2)));
  gen _au_unit(mksa_register("au",makevecteur(1.495979e11,1,0,0,0)));
  gen _Angstrom_unit(mksa_register("Angstrom",makevecteur(1e-10,1,0,0,0)));
  gen _micron_unit(mksa_register("µ",makevecteur(1e-6,1,0,0,0)));
  gen _b_unit(mksa_register("b",makevecteur(1e28,1,0,0,0)));
  gen _bar_unit(mksa_register("bar",makevecteur(1e5,-1,1,-2,0)));
  gen _bbl_unit(mksa_register("bbl",makevecteur(.158987294928,3,0,0,0)));
  gen _buUS(mksa_register("buUS",makevecteur(0.03523907,3,0,0,0)));
  gen _bu(mksa_register("bu",makevecteur(0.036368736,3,0,0,0)));
  gen _Btu_unit(mksa_register("Btu",makevecteur(1055.05585262,2,1,-2,0)));
  gen _cal_unit(mksa_register("cal",makevecteur(4.1868,2,1,-2,0)));
  gen _chain_unit(mksa_register("chain",makevecteur(20.1168402337,1,0,0,0)));
  gen _Ci_unit(mksa_register("Ci",makevecteur(3.7e10,0,0,-1,0)));
  gen _ct_unit(mksa_register("ct",makevecteur(0.0002,0,1,0,0)));
  gen _deg_unit(mksa_register("deg",makevecteur(1.74532925199e-2,0,0,0,0)));
  gen _d_unit(mksa_register("d",makevecteur(86400,0,0,1,0)));
  gen _dB_unit(mksa_register("dB",makevecteur(1,0,0,0,0)));
  gen _dyn_unit(mksa_register("dyn",makevecteur(1e-5,1,1,-2,0)));
  gen _erg_unit(mksa_register("erg",makevecteur(1e-7,2,1,-2,0)));
  gen _eV_unit(mksa_register("eV",makevecteur(1.60217733e-19,2,1,-2,0)));
  gen _degreeF_unit(mksa_register("degreeF",makevecteur(gen(5)/9,0,0,0,0,1)));
  gen _Rankine_unit(mksa_register("Rankine",makevecteur(gen(5)/9,0,0,0,0,1)));
  gen _fath_unit(mksa_register("fath",makevecteur(1.82880365761,1,0,0,0)));
  gen _fm_unit(mksa_register("fm",makevecteur(1.82880365761,1,0,0,0)));
  gen _fbm_unit(mksa_register("fbm",makevecteur(0.002359737216,3,0,0,0)));
  // gen _fc_unit(mksa_register("fc",makevecteur(10.7639104167,1,0,0,0)));
  gen _Fdy_unit(mksa_register("Fdy",makevecteur(96487,0,0,1,1)));
  gen _fermi_unit(mksa_register("fermi",makevecteur(1e-15,1,0,0,0)));
  gen _flam_unit(mksa_register("flam",makevecteur(3.42625909964,-2,0,0,0,0,0,1)));
  gen _ft_unit(mksa_register("ft",makevecteur(0.3048,1,0,0,0)));
  gen _ftUS_unit(mksa_register("ftUS",makevecteur(0.304800609601,1,0,0,0)));
  gen _Gal(mksa_register("Gal",makevecteur(0.01,1,0,-2,0)));
  gen _g_unit(mksa_register("g",makevecteur(1e-3,0,1,0,0)));
  gen _galUS_unit(mksa_register("galUS",makevecteur(0.003785411784,3,0,0,0)));
  gen _galC_unit(mksa_register("galC",makevecteur(0.00454609,3,0,0,0)));
  gen _galUK_unit(mksa_register("galUK",makevecteur(0.004546092,3,0,0,0)));
  gen _gf_unit(mksa_register("gf",makevecteur(0.00980665,1,1,-2,0)));
  gen _gmol_unit(mksa_register("gmol",makevecteur(1,0,0,0,0,0,1)));
  gen _grad_unit(mksa_register("grad",makevecteur(1.57079632679e-2)));
  gen _gon_unit(mksa_register("gon",makevecteur(1.57079632679e-2)));
  gen _grain_unit(mksa_register("grain",makevecteur(0.00006479891,0,1,0,0)));
  gen _ha_unit(mksa_register("ha",makevecteur(10000,2,0,0,0)));
  gen _h_unit(mksa_register("h",makevecteur(3600,0,0,1,0)));
  gen _hp_unit(mksa_register("hp",makevecteur(745.699871582,2,1,-3,0)));
  gen _in_unit(mksa_register("in",makevecteur(0.0254,1,0,0,0)));
  gen _inHg_unit(mksa_register("inHg",makevecteur(3386.38815789,-1,1,-2,0)));
  gen _inH2O_unit(mksa_register("inH2O",makevecteur(248.84,-1,1,-2,0)));
  gen _j_unit(mksa_register("j",makevecteur(86400,0,0,1,0)));
  gen _FF_unit(mksa_register("FF",makevecteur(.152449017237,0,0,0,0,0,0,0,1)));
  gen _kip_unit(mksa_register("kip",makevecteur(4448.22161526,1,1,-2,0)));
  gen _knot_unit(mksa_register("knot",makevecteur(0.51444444444,1,0,-1,0)));
  gen _kph_unit(mksa_register("kph",makevecteur(0.2777777777777,1,0,-1,0)));
  gen _l_unit(mksa_register("l",makevecteur(0.001,3,0,0,0)));
  gen _L_unit(mksa_register("L",makevecteur(0.001,3,0,0,0)));
  gen _lam_unit(mksa_register("lam",makevecteur(3183.09886184,-2,0,0,0,0,0,1)));
  gen _lb_unit(mksa_register("lb",makevecteur(0.45359237,0,1,0,0)));
  gen _lbf_unit(mksa_register("lbf",makevecteur(4.44922161526,1,1,-2,0)));
  gen _lbmol_unit(mksa_register("lbmol",makevecteur(453.59237,0,0,0,0,0,1)));
  gen _lbt_unit(mksa_register("lbt",makevecteur(0.3732417216,0,1,0,0)));
  gen _lyr_unit(mksa_register("lyr",makevecteur(9.46052840488e15,1,0,0,0)));
  gen _mi_unit(mksa_register("mi",makevecteur(1609.344,1,0,0,0)));
  gen _mil_unit(mksa_register("mil",makevecteur(0.0000254,1,0,0,0)));
  gen _mile_unit(mksa_register("mile",makevecteur(1609.344,1,0,0,0)));
  gen _mille_unit(mksa_register("mille",makevecteur(1852,1,0,0,0)));
  gen _mn_unit(mksa_register("mn",makevecteur(60,0,0,1,0)));
  gen _miUS_unit(mksa_register("miUS",makevecteur(1609.34721869,1,0,0,0)));
  gen _mmHg_unit(mksa_register("mmHg",makevecteur(133.322368421,-1,1,-2,0)));
  gen _mph_unit(mksa_register("mph",makevecteur(0.44704,1,0,-1,0)));
  gen _nmi_unit(mksa_register("nmi",makevecteur(1852,1,0,0,0)));
  gen _oz_unit(mksa_register("oz",makevecteur(0.028349523125,0,1,0,0)));
  gen _ozfl_unit(mksa_register("ozfl",makevecteur(2.95735295625e-5,3,0,0,0)));
  gen _ozt_unit(mksa_register("ozt",makevecteur(0.0311034768,0,1,0,0)));
  gen _ozUK_unit(mksa_register("ozUK",makevecteur(2.8413075e-5,3,0,0,0)));
  gen _P_unit(mksa_register("P",makevecteur(.1,-1,1,-1,0)));
  gen _pc_unit(mksa_register("pc",makevecteur(3.08567818585e16,1,0,0,0)));
  gen _pdl_unit(mksa_register("pdl",makevecteur(0.138254954376,1,1,-2,0)));
  gen _pk_unit(mksa_register("pk",makevecteur(0.0088097675,3,0,0,0)));
  gen _psi_unit(mksa_register("psi",makevecteur(6894.75729317,-1,1,-2,0)));
  gen _pt_unit(mksa_register("pt",makevecteur(0.000473176473,3,0,0,0)));
  gen _ptUK_unit(mksa_register("ptUK",makevecteur(0.0005682615,3,0,0,0)));
  gen _liqpt_unit(mksa_register("liqpt",makevecteur(0.000473176473,3,0,0,0)));
  gen _qt_unit(mksa_register("qt",makevecteur(0.000946359246,3,0,0,0)));
  gen _R_unit(mksa_register("R",makevecteur(0.000258,0,-1,1,1)));
  gen _rd_unit(mksa_register("rd",makevecteur(0.01,2,0,-2,0)));
  gen _rod_unit(mksa_register("rod",makevecteur(5.02921005842,1,0,0,0)));
  gen _rem_unit(mksa_register("rem",makevecteur(0.01,2,0,-2,0)));
  gen _rpm_unit(mksa_register("rpm",makevecteur(0.0166666666667,0,0,-1,0)));
  gen _sb_unit(mksa_register("sb",makevecteur(10000,-2,0,0,0,0,0,1)));
  gen _slug_unit(mksa_register("slug",makevecteur(14.5939029372,0,1,0,0)));
  gen _St_unit(mksa_register("St",makevecteur(0.0001,2,0,-1,0)));
  gen _t_unit(mksa_register("t",makevecteur(1000,0,1,0,0)));
  gen _tbsp_unit(mksa_register("tbsp",makevecteur(1.47867647813e-5,3,0,0,0)));
  gen _tex(mksa_register("tex",makevecteur(1e-6,-1,1,0,0)));
  gen _therm_unit(mksa_register("therm",makevecteur(105506000,2,1,-2,0)));
  gen _ton_unit(mksa_register("ton",makevecteur(907.18474,0,1,0,0)));
  gen _tonUK_unit(mksa_register("tonUK",makevecteur(1016.0469088,0,1,0,0)));
  gen _torr_unit(mksa_register("torr",makevecteur(133.322368421,-1,1,-2,0)));
  gen _tr_unit(mksa_register("tr",makevecteur(cst_two_pi))); // radian
  gen _u_unit(mksa_register("u",makevecteur(1.6605402e-27,0,1,0,0)));
  gen _yd_unit(mksa_register("yd",makevecteur(0.9144,1,0,0,0)));
  gen _yr_unit(mksa_register("yr",makevecteur(31556925.9747,0,0,1,0)));

  // Some hydrocarbur energy equivalent
  // tep=tonne equivalent petrole, lep litre equivalent petrole
  // toe=(metric) ton of oil equivalent
  // bblep = baril equivalent petrole, boe=baril of oil equivalent
  gen _tep_unit(mksa_register("tep",makevecteur(41.76e9,2,1,-2,0)));
  gen _toe_unit(mksa_register("toe",makevecteur(41.76e9,2,1,-2,0)));
  gen _cf_unit(mksa_register("cf",makevecteur(1.08e6,2,1,-2,0)));
  gen _tec_unit(mksa_register("tec",makevecteur(41.76e9/1.5,2,1,-2,0)));
  gen _lep_unit(mksa_register("lep",makevecteur(0.857*41.76e6,2,1,-2,0)));
  gen _bblep_unit(mksa_register("bblep",makevecteur(.158987294928*0.857*41.76e9,2,1,-2,0)));
  gen _boe_unit(mksa_register("boe",makevecteur(.158987294928*0.857*41.76e9,2,1,-2,0)));
  gen _Wh_unit(mksa_register("Wh",makevecteur(3600,2,1,-2,0)));
  // Equivalent Carbon for 1 tep, oil, gas, coal
  gen _tepC_unit(mksa_register("tepC",makevecteur(830,1,0,0,0)));
  gen _tepgC_unit(mksa_register("tepgC",makevecteur(650,1,0,0,0)));
  gen _tepcC_unit(mksa_register("tepcC",makevecteur(1000,1,0,0,0)));
  // mean PRG for HFC in kg C unit
  gen _HFCC_unit(mksa_register("HFCC",makevecteur(1400,1,0,0,0)));

  // return a vector of powers in MKSA system
  vecteur mksa_convert(const identificateur & g,GIAC_CONTEXT){
    string s=g.print(contextptr);
    // Find prefix in unit
    int exposant=0;
    int l=s.size();
    if (l>1 && s[0]=='_'){
      --l;
      s=s.substr(1,l);
    }
    else
      return makevecteur(g);
    gen res=plus_one;
    std::map<string,gen>::const_iterator it=unit_conversion_map.find(s),itend=unit_conversion_map.end();
    if (it==itend && l>1){
      switch (s[0]){
      case 'Y':
	exposant=24;
	break;
      case 'Z':
	exposant=21;
	break;
      case 'E':
	exposant=18;
	break;
      case 'P':
	exposant=15;
	break;
      case 'T':
	exposant=12;
	break;
      case 'G':
	exposant=9;
	break;
      case 'M':
	exposant=6;
	break;
      case 'K': case 'k':
	exposant=3;
	break;
      case 'H': case 'h':
	exposant=2;
	break;
      case 'D':
	exposant=1;
	break;
      case 'd':
	exposant=-1;
	break;
      case 'c':
	exposant=-2;
	break;
      case 'm':
	exposant=-3;
	break;
      case 'µ':
	exposant=-6;
	break;
      case 'n':
	exposant=-9;
	break;
      case 'p':
	exposant=-12;
	break;
      case 'f':
	exposant=-15;
	break;
      case 'a':
	exposant=-18;
	break;
      case 'z':
	exposant=-21;
	break;
      case 'y':
	exposant=-24;
	break;
      }
    }
    if (exposant!=0){
      s=s.substr(1,l-1);
      res=pow(gen(10),gen(exposant),contextptr);
      it=unit_conversion_map.find(s);
    }
    if (it==itend)
      return makevecteur(res*find_or_make_symbol("_"+s,contextptr));
    gen tmp=it->second;
    if (tmp.type!=_VECT || tmp._VECTptr->size()<4)
      return makevecteur(g);
    vecteur v=*tmp._VECTptr;
    v[0]=res*v[0];
    return v;
  }

  vecteur mksa_convert(const gen & g,GIAC_CONTEXT){
    if (g.type==_IDNT)
      return mksa_convert(*g._IDNTptr,contextptr);
    if (g.type!=_SYMB)
      return makevecteur(g);
    if (g.is_symb_of_sommet(at_unit)){
      vecteur & v=*g._SYMBptr->feuille._VECTptr;
      vecteur res=mksa_convert(v[1],contextptr);
      res[0]=v[0]*res[0];
      return res;
    }
    if (g._SYMBptr->sommet==at_inv){
      vecteur res(mksa_convert(g._SYMBptr->feuille,contextptr));
      res[0]=inv(res[0],contextptr);
      int s=res.size();
      for (int i=1;i<s;++i)
	res[i]=-res[i];
      return res;
    }
    if (g._SYMBptr->sommet==at_pow){
      gen & f=g._SYMBptr->feuille;
      if (f.type!=_VECT||f._VECTptr->size()!=2)
	setsizeerr();
      vecteur res(mksa_convert(f._VECTptr->front(),contextptr));
      gen e=f._VECTptr->back();
      res[0]=pow(res[0],e,contextptr);
      int s=res.size();
      for (int i=1;i<s;++i)
	res[i]=e*res[i];
      return res;
    }
    if (g._SYMBptr->sommet==at_prod){
      gen & f=g._SYMBptr->feuille;
      if (f.type!=_VECT)
	return mksa_convert(f,contextptr);
      vecteur & v=*f._VECTptr;
      vecteur res(makevecteur(plus_one));
      const_iterateur it=v.begin(),itend=v.end();
      for (;it!=itend;++it){
	vecteur tmp(mksa_convert(*it,contextptr));
	res[0]=res[0]*tmp[0];
	iterateur it=res.begin()+1,itend=res.end(),jt=tmp.begin()+1,jtend=tmp.end();
	for (;it!=itend && jt!=jtend;++it,++jt)
	  *it=*it+*jt;
	for (;jt!=jtend;++jt)
	  res.push_back(*jt);
      }
      return res;
    }
    return makevecteur(g);
  }

  gen unitpow(const gen & g,const gen & exponent){
    if (is_zero(exponent))
      return plus_one;
    if (is_one(exponent))
      return g;
    return symbolic(at_pow,makevecteur(g,exponent));
  }
  gen mksa_reduce(const gen & g,GIAC_CONTEXT){
    vecteur v(mksa_convert(g,contextptr));
    gen res1=v[0];
    gen res=plus_one;
    int s=v.size();
    if (s>2)
      res = res *unitpow(_kg_unit,v[2]);
    if (s>1)
      res = res *unitpow(_m_unit,v[1]);
    if (s>3)
      res = res *unitpow(_s_unit,v[3]);
    if (s>4)
      res = res * unitpow(_A_unit,v[4]);
    if (s>5)
      res = res * unitpow(_K_unit,v[5]);
    if (s>6)
      res = res * unitpow(_mol_unit,v[6]);
    if (s>7)
      res = res * unitpow(_cd_unit,v[7]);
    if (s>8)
      res = res * unitpow(_E_unit,v[8]);
    if (is_one(res))
      return res1;
    else
      return symbolic(at_unit,makevecteur(res1,res));
  }
  const string _mksa_s("mksa");
  unary_function_eval __mksa(&mksa_reduce,_mksa_s);
  unary_function_ptr at_mksa (&__mksa,0,true);
  
  gen _ufactor(const gen & g,GIAC_CONTEXT){
    if (g.type==_VECT && g.subtype==_SEQ__VECT && g._VECTptr->size()==2){
      vecteur & v=*g._VECTptr;
      return v.back()*mksa_reduce(v.front()/v.back(),contextptr);
    }
    setsizeerr();
    return 0;
  }
  const string _ufactor_s("ufactor");
  unary_function_eval __ufactor(&_ufactor,_ufactor_s);
  unary_function_ptr at_ufactor (&__ufactor,0,true);
  
  gen _usimplify(const gen & g,GIAC_CONTEXT){
    if (g.type==_VECT)
      return apply(g,_usimplify,contextptr);
    if (!g.is_symb_of_sommet(at_unit))
      return g;
    vecteur v=mksa_convert(g,contextptr);
    gen res1=v[0];
    int s=v.size();
    if (s>5)
      return g;
    for (int i=s;i<5;++i)
      v.push_back(zero);
    // look first if it's a mksa
    int pos=0;
    for (int i=1;i<5;++i){
      if (v[i]==zero)
	continue;
      if (pos){
	pos=0;
	break;
      }
      pos=i;
    }
    if (pos)
      return mksa_reduce(g,contextptr);
    v[0]=plus_one;
    const_iterateur it=usual_units.begin(),itend=usual_units.end();
    for (;it!=itend;++it){
      string s=it->print(contextptr);
      gen tmp=unit_conversion_map[s.substr(1,s.size()-1)];
      if (tmp==v)
	return _ufactor(gen(makevecteur(g,symbolic(at_unit,makevecteur(1,*it))),_SEQ__VECT),contextptr);
    }
    it=usual_units.begin();
    for (;it!=itend;++it){
      string s=it->print(contextptr);
      gen tmp=unit_conversion_map[s.substr(1,s.size()-1)];
      vecteur w(*tmp._VECTptr);
      for (int j=0;j<2;j++){
	vecteur vw;
	if (j)
	  vw=addvecteur(v,w);
	else
	  vw=subvecteur(v,w);
	for (int i=1;i<5;++i){
	  if (vw[i]==zero)
	    continue;
	  if (pos){
	    pos=0;
	    break;
	  }
	  pos=i;
	}
	if (pos){
	  if (j)
	    return _ufactor(gen(makevecteur(g,symbolic(at_unit,makevecteur(1,unitpow(*it,-1)))),_SEQ__VECT),contextptr);
	  else
	    return _ufactor(gen(makevecteur(g,symbolic(at_unit,makevecteur(1,*it))),_SEQ__VECT),contextptr);
	}
      }
    }
    return g;
  }
  const string _usimplify_s("usimplify");
  unary_function_eval __usimplify(&_usimplify,_usimplify_s);
  unary_function_ptr at_usimplify (&__usimplify,0,true);
  
  gen symb_unit(const gen & a,const gen & b,GIAC_CONTEXT){
    // Add a _ to all identifiers in b
    vecteur v(lidnt(b));
    vecteur w(v);
    iterateur it=w.begin(),itend=w.end();
    for (;it!=itend;++it){
      find_or_make_symbol("_"+it->print(contextptr),*it,contextptr);
    }
    return symbolic(at_unit,makevecteur(a,subst(b,v,w,false,contextptr)));
  }
  string printasunit(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if (feuille.type!=_VECT || feuille._VECTptr->size()!=2)
      setsizeerr();
    vecteur & v=*feuille._VECTptr;
    vecteur v1(lidnt(v[1]));
    vecteur w(v1);
    iterateur it=w.begin(),itend=w.end();
    for (;it!=itend;++it){
      string s(it->print(contextptr));
      if (!s.empty() && s[0]=='_')
	s=s.substr(1,s.size()-1);
      find_or_make_symbol(s,*it,contextptr);
    }
    string tmp(subst(v[1],v1,w,false,contextptr).print(contextptr));
    if (tmp[0]=='c' || (v[1].type==_SYMB && !v[1].is_symb_of_sommet(at_pow)) )
      tmp="_("+tmp+")";
    else
      tmp="_"+tmp;
    if (v[0].type<_POLY)
      return v[0].print(contextptr)+tmp;
    else
      return "("+v[0].print(contextptr)+")"+tmp;
  }
  gen unit(const gen & g){
    if (g.type!=_VECT || g._VECTptr->size()!=2)
      setsizeerr();
    return symbolic(at_unit,g);
  }
  const string _unit_s("_");
  unary_function_unary __unit(&unit,_unit_s,&printasunit);
  unary_function_ptr at_unit (&__unit);

  unary_function_ptr binary_op_tab[]={at_plus,at_prod,at_pow,at_and,at_ou,at_xor,at_different,at_same,at_equal,at_unit,at_compose,at_composepow,at_deuxpoints,at_tilocal,at_pointprod,at_pointdivision,at_pointpow,at_division,at_normalmod,at_minus,at_intersect,at_union,at_interval,at_inferieur_egal,at_inferieur_strict,at_superieur_egal,at_superieur_strict,0};
  // unary_function_ptr binary_op_tab[]={at_and,at_ou,at_different,at_same,0};

  // Physical constants
  identificateur _cst_hbar("_hbar_",symbolic(at_unit,makevecteur(1.05457266e-34,_J_unit*_s_unit)));
  gen cst_hbar(_cst_hbar);
  identificateur _cst_clightspeed("_c_",symbolic(at_unit,makevecteur(299792458,_m_unit/_s_unit)));
  gen cst_clightspeed(_cst_clightspeed);
  identificateur _cst_ga("_g_",symbolic(at_unit,makevecteur(9.80665,_m_unit*unitpow(_s_unit,2))));
  gen cst_ga(_cst_ga);
  identificateur _cst_IO("_IO_",symbolic(at_unit,makevecteur(1e-12,_W_unit*unitpow(_m_unit,-2))));
  // gen cst_IO("_io",context0); //  IO 1e-12W/m^2
  gen cst_IO(_cst_IO);
  identificateur _cst_epsilonox("_epsilonox_",3.9);
  gen cst_epsilonox(_cst_epsilonox); // 3.9
  identificateur _cst_epsilonsi("_epsilonsi_",11.9);
  gen cst_epsilonsi(_cst_epsilonsi); // 11.9
  identificateur _cst_qepsilon0("_qepsilon0_",symbolic(at_unit,makevecteur(1.4185979e-30,_F_unit*_C_unit/_m_unit)));
  gen cst_qepsilon0(_cst_qepsilon0); // qeps0 1.4185979e-30 F*C/m
  identificateur _cst_epsilon0q("_epsilon0q_",symbolic(at_unit,makevecteur(55263469.6,_F_unit/(_m_unit*_C_unit))));
  gen cst_epsilon0q(_cst_epsilon0q); // eps0q 55263469.6 F/(m*C)
  identificateur _cst_kq("_kq_",symbolic(at_unit,makevecteur(8.617386e-5,_J_unit/(_K_unit*_C_unit))));
  gen cst_kq(_cst_kq); // kq 8.617386e-5 J/(K*C)
  identificateur _cst_c3("_c3_",symbolic(at_unit,makevecteur(.002897756,_m_unit*_K_unit)));
  gen cst_c3(_cst_c3); // c3 .002897756m*K
  identificateur _cst_lambdac("_lambdac_",symbolic(at_unit,makevecteur( 0.00242631058e-9,_m_unit)));
  gen cst_lambdac(_cst_lambdac); // lambdac 0.00242631058 nm
  identificateur _cst_f0("_f0_",symbolic(at_unit,makevecteur(2.4179883e14,_Hz_unit)));
  gen cst_f0(_cst_f0); //  f0 2.4179883e14Hz
  identificateur _cst_lambda0("_lambda0_",symbolic(at_unit,makevecteur(1239.8425e-9,_m_unit)));
  gen cst_lambda0(_cst_lambda0); // lambda0 1239.8425_nm
  identificateur _cst_muN("_muN_",symbolic(at_unit,makevecteur(5.0507866e-27,_J_unit/_T_unit)));
  gen cst_muN(_cst_muN); // muN 5.0507866e-27_J/T
  identificateur _cst_muB("_muB_",symbolic(at_unit,makevecteur( 9.2740154e-24,_J_unit/_T_unit)));
  gen cst_muB(_cst_muB); // muB 9.2740154e-24 J/T
  identificateur _cst_a0("_a0_",symbolic(at_unit,makevecteur(.0529177249e-9,_m_unit)));
  gen cst_a0(_cst_a0); // a0 .0529177249_nm
  identificateur _cst_Rinfinity("_Rinfinity_",symbolic(at_unit,makevecteur(10973731.534,unitpow(_m_unit,-1))));
  gen cst_Rinfinity(_cst_Rinfinity); // Rinf 10973731.534 m^-1
  identificateur _cst_Faraday("_Faraday_",symbolic(at_unit,makevecteur(96485.309,_C_unit/_mol_unit)));
  gen cst_Faraday(_cst_Faraday); // F 96485.309 C/gmol
  identificateur _cst_phi("_phi_",symbolic(at_unit,makevecteur(2.06783461e-15,_Wb_unit)));
  gen cst_phi(_cst_phi); // phi 2.06783461e-15 Wb
  identificateur _cst_alpha("_alpha_",7.29735308e-3);
  gen cst_alpha(_cst_alpha); // alpha 7.29735308e-3
  identificateur _cst_mpme("_mpme_",1836.152701);
  gen cst_mpme(_cst_mpme); // mpme 1836.152701
  identificateur _cst_mp("_mp_",symbolic(at_unit,makevecteur(1.6726231e-27,_kg_unit)));
  gen cst_mp(_cst_mp); // mp 1.6726231e-27 kg
  identificateur _cst_qme("_qme_",symbolic(at_unit,makevecteur(1.75881962e11,_C_unit/_kg_unit)));
  gen cst_qme(_cst_qme); // qme 175881962000 C/kg
  identificateur _cst_me("_me_",symbolic(at_unit,makevecteur(9.1093897e-31,_kg_unit)));
  gen cst_me(_cst_me); // me 9.1093897e-31 kg
  identificateur _cst_qe("_qe_",symbolic(at_unit,makevecteur(1.60217733e-19,_C_unit)));
  gen cst_qe(_cst_qe); // q 1.60217733e-19 C
  identificateur _cst_hPlanck("_h_",symbolic(at_unit,makevecteur(6.6260755e-34,_J_unit*_s_unit)));
  gen cst_hPlanck(_cst_hPlanck); //  h 6.6260755e-34 Js
  identificateur _cst_G("_G_",symbolic(at_unit,makevecteur(6.67259e-11,unitpow(_m_unit,3)*unitpow(_s_unit,-2)*unitpow(_kg_unit,-1))));
  gen cst_G(_cst_G); // G 6.67259e-11m^3/s^2kg
  identificateur _cst_mu0("_mu0_",symbolic(at_unit,makevecteur(1.25663706144e-6,_H_unit/_m_unit)));
  gen cst_mu0(_cst_mu0); // mu0 1.25663706144e-6 H/m
  identificateur _cst_epsilon0("_epsilon0_",symbolic(at_unit,makevecteur(8.85418781761e-12,_F_unit/_m_unit)));
  gen cst_epsilon0(_cst_epsilon0); // eps0 8.85418781761e-12 F/m
  identificateur _cst_sigma("_sigma_",symbolic(at_unit,makevecteur( 5.67051e-8,_W_unit*unitpow(_m_unit,-2)*unitpow(_K_unit,-4))));
  gen cst_sigma(_cst_sigma); // sigma 5.67051e-8 W/m^2*K^4
  identificateur _cst_StdP("_StdP_",symbolic(at_unit,makevecteur(101325.0,_Pa_unit)));
  gen cst_StdP(_cst_StdP); // StdP 101.325_kPa
  identificateur _cst_StdT("_StdT_",symbolic(at_unit,makevecteur(273.15,_K_unit)));
  gen cst_StdT(_cst_StdT); // StdT 273.15_K
  identificateur _cst_Rydberg("_R_",symbolic(at_unit,makevecteur(8.31451,_J_unit/_molK_unit)));
  gen cst_Rydberg(_cst_Rydberg); // Rydberg 8.31451_J/(gmol*K)
  identificateur _cst_Vm("_Vm_",symbolic(at_unit,makevecteur(22.4141,_l_unit/_mol_unit)));
  gen cst_Vm(_cst_Vm); // Vm 22.4141_l/gmol
  identificateur _cst_kBoltzmann("_k_",symbolic(at_unit,makevecteur(1.380658e-23,_J_unit/_K_unit)));
  gen cst_kBoltzmann(_cst_kBoltzmann); // k 1.380658e-23 J/K
  identificateur _cst_NA("_NA_",symbolic(at_unit,makevecteur(6.0221367e23,unitpow(_mol_unit,-1))));
  gen cst_NA(_cst_NA); // NA 6.0221367e23 1/gmol

  gen maple_root(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT || g._VECTptr->size()!=2)
      return symbolic(at_maple_root,g);
    vecteur & v=*g._VECTptr;
    return pow(v[1],inv(v[0],contextptr),contextptr);
  }
  const string _maple_root_s("root");
  unary_function_eval __maple_root(&maple_root,_maple_root_s);
  unary_function_ptr at_maple_root (&__maple_root);

  gen symb_interrogation(const gen & e1,const gen & e3){
    if (e3.is_symb_of_sommet(at_deuxpoints)){
      gen & f =e3._SYMBptr->feuille;
      if (f.type==_VECT && f._VECTptr->size()==2)
	return symb_when(e1,f._VECTptr->front(),f._VECTptr->back());
    }
    return symb_when(e1,e3,undef);
  }

  bool first_ascend_sort(const gen & a,const gen & b){
    gen g=inferieur_strict(a[0],b[0],context0); 
    if (g.type!=_INT_)
      return a[0].islesscomplexthan(b[0]);
    return g.val==1;
  }
  bool first_descend_sort(const gen & a,const gen & b){
    gen g=superieur_strict(a[0],b[0],context0); 
    if (g.type!=_INT_)
      return !a[0].islesscomplexthan(b[0]);
    return g.val==1;
  }

  // Create an operator with a given syntax
  vector<unary_function_ptr> user_operator_list;   // GLOBAL VAR
  gen user_operator(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT || g._VECTptr->size()<3)
      setsizeerr();
    vecteur & v=*g._VECTptr;
    // int s=signed(v.size());
    if (v[0].type!=_STRNG)
      return string2gen("Operator name must be of type string",false);
    string & ss=*v[0]._STRNGptr;
    vector<unary_function_ptr>::iterator it=user_operator_list.begin(),itend=user_operator_list.end();
    for (;it!=itend;++it){
      if (it->ptr->s==ss){
	break;
      }
    }
    if (it!=itend){
      const unary_function_abstract * ptr0=it->ptr;
      const unary_function_user * ptr=dynamic_cast<const unary_function_user *>(ptr0);
      if (!ptr)
	return zero;
      if (ptr->f==v[1])
	return plus_one;
      return zero;
    }
    if (v[2].type==_INT_){ 
      int token_value=v[2].val;
      unary_function_user * uf;
      if (v[2].subtype==_INT_MUPADOPERATOR){
	switch (v[2].val){
	case _POSTFIX_OPERATOR:
	  uf= new unary_function_user (v[1],ss,0,0,0);
	  token_value=T_FACTORIAL; // like factorial
	  break;
	case _PREFIX_OPERATOR:
	  uf=new unary_function_user(v[1],ss,0,0,0);
	  token_value=T_NOT; // like not
	  break;
	case _BINARY_OPERATOR:
	  uf = new unary_function_user (v[1],ss);
	  token_value=T_FOIS; // like *
	  break;
	default:
	  return zero;
	}
      }
      else 
	// non mupad syntax, v[2] is input_parser.yy token value
	uf = new unary_function_user(v[1],ss);
      unary_function_ptr u(uf);
      // cout << symbolic(u,makevecteur(1,2)) << endl;
      user_operator_list.push_back(u);
      bool res=lexer_functions_register(u,ss,token_value);
      if (res){
	if (!child_id)
	  _signal(symb_quote(symbolic(at_user_operator,g)),contextptr);
	return plus_one;
      }
      user_operator_list.pop_back();
      delete uf;
    }
    return zero;
  }
  const string _user_operator_s("user_operator");
  unary_function_eval __user_operator(&user_operator,_user_operator_s);
  unary_function_ptr at_user_operator (&__user_operator);

  gen current_folder_name;

  gen getfold(const gen & g){
    if (is_zero(g))
      return string2gen("main",false);
    return g;
  }

  gen _SetFold(const gen & g,GIAC_CONTEXT){
    if (!is_zero(g) && g.type!=_IDNT)
      setsizeerr();
    bool ok=is_zero(g);
    if (g.type==_IDNT && g._IDNTptr->value && g._IDNTptr->value->type==_VECT && g._IDNTptr->value->subtype==_FOLDER__VECT)
      ok=true;
    if ( ok || (g.type==_IDNT && g._IDNTptr->name && (*g._IDNTptr->name=="main"|| *g._IDNTptr->name=="home") ) ){
      gen res=current_folder_name;
      current_folder_name=g;
      if (!child_id)
	_signal(symb_quote(symbolic(at_SetFold,g)),contextptr);
      return getfold(res);
    }
    setsizeerr("Non existent Folder");
    return 0;
  }
  const string _SetFold_s("SetFold");
  unary_function_eval __SetFold(&_SetFold,_SetFold_s,&printastifunction);
  unary_function_ptr at_SetFold (&__SetFold,_QUOTE_ARGUMENTS,T_RETURN); 

  gen _piecewise(const gen & g,GIAC_CONTEXT){
    // evaluate couples of condition/expression, like in a case
    if (g.type!=_VECT)
      return g;
    vecteur & v =*g._VECTptr;
    int s=v.size();
    gen test;
    for (int i=0;i<s/2;++i){
      test=v[2*i];
      test=equaltosame(test.eval(eval_level(contextptr),contextptr)).eval(eval_level(contextptr),contextptr);
      test=test.evalf_double(eval_level(contextptr),contextptr);
      if ( (test.type!=_DOUBLE_) && (test.type!=_CPLX) )
	return symbolic(at_piecewise,g.eval(eval_level(contextptr),contextptr));
      if (is_zero(test))
	continue;
      return v[2*i+1].eval(eval_level(contextptr),contextptr);
    }
    if (s%2)
      return v[s-1].eval(eval_level(contextptr),contextptr);
    setsizeerr();
    return undef;
  }
  const string _piecewise_s("piecewise");
  unary_function_eval __piecewise(&_piecewise,_piecewise_s);
  unary_function_ptr at_piecewise (&__piecewise,_QUOTE_ARGUMENTS,true);

  gen _geo2d(const gen & g,GIAC_CONTEXT){
    return g;
  }
  const string _geo2d_s("geo2d");
  unary_function_eval __geo2d(&_geo2d,_geo2d_s);
  unary_function_ptr at_geo2d (&__geo2d,0,true);

  const string _geo3d_s("geo3d");
  unary_function_eval __geo3d(&_geo2d,_geo3d_s);
  unary_function_ptr at_geo3d (&__geo3d,0,true);

  const string _spreadsheet_s("spreadsheet");
  unary_function_eval __spreadsheet(&_geo2d,_spreadsheet_s);
  unary_function_ptr at_spreadsheet (&__spreadsheet,0,true);

  std::string print_program_syntax(int maple_mode){
    string logs;
    switch (maple_mode){
    case 0:
      logs="xcas";
      break;
    case 1:
      logs="maple";
      break;
    case 2:
      logs="mupad";
      break;
    case 3:
      logs="ti";
      break;
    default:
      logs=print_INT_(maple_mode);
    }
    return logs;
  }

  gen _threads_allowed(const gen & g,GIAC_CONTEXT){
    if (is_zero(g))
      threads_allowed=false;
    else
      threads_allowed=true;
    return threads_allowed;
  }
  const string _threads_allowed_s("threads_allowed");
  unary_function_eval __threads_allowed(&_threads_allowed,_threads_allowed_s);
  unary_function_ptr at_threads_allowed (&__threads_allowed,0,true);

  gen _mpzclass_allowed(const gen & g,GIAC_CONTEXT){
    if (is_zero(g))
      mpzclass_allowed=false;
    else
      mpzclass_allowed=true;
    return mpzclass_allowed;
  }
  const string _mpzclass_allowed_s("mpzclass_allowed");
  unary_function_eval __mpzclass_allowed(&_mpzclass_allowed,_mpzclass_allowed_s);
  unary_function_ptr at_mpzclass_allowed (&__mpzclass_allowed,0,true);

  gen whentopiecewise(const gen & g,GIAC_CONTEXT){
    return symbolic(at_piecewise,g);
  }
  vector< unary_function_ptr > when_v(1,at_when);
  vector< gen_op_context > when2piecewise_v(1,whentopiecewise);
  gen when2piecewise(const gen & g,GIAC_CONTEXT){
    return subst(g,when_v,when2piecewise_v,false,contextptr);
  }

  gen piecewisetowhen(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT)
      return g;
    vecteur v = *g._VECTptr;
    int s=v.size();
    if (s==1)
      setsizeerr();
    if (s==2){
      v.push_back(undef);
      s++;
    }
    if (s==3)
      return symbolic(at_when,g);
    gen tmp=piecewisetowhen(vecteur(v.begin()+2,v.end()),contextptr);
    return symbolic(at_when,gen(makevecteur(v[0],v[1],tmp),_SEQ__VECT));
  }
  vector< unary_function_ptr > piecewise_v(1,at_piecewise);
  vector< gen_op_context > piecewise2when_v(1,piecewisetowhen);
  gen piecewise2when(const gen & g,GIAC_CONTEXT){
    return subst(g,piecewise_v,piecewise2when_v,false,contextptr);
  }

  gen whentosign(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT || g._VECTptr->size()!=3)
      setsizeerr();
    vecteur v = *g._VECTptr;
    if (v[0].is_symb_of_sommet(at_equal) || v[0].is_symb_of_sommet(at_same)){
      *logptr(contextptr) << "Assuming false condition " << v[0].print(contextptr) << endl;
      return v[2];
    }
    if (v[0].is_symb_of_sommet(at_different)){
      *logptr(contextptr) << "Assuming true condition " << v[0].print(contextptr) << endl;
      return v[1];
    }
    bool ok=false;
    if (v[0].is_symb_of_sommet(at_superieur_strict) || v[0].is_symb_of_sommet(at_superieur_egal)){
      v[0]=v[0]._SYMBptr->feuille[0]-v[0]._SYMBptr->feuille[1];
      ok=true;
    }
    if (!ok && (v[0].is_symb_of_sommet(at_inferieur_strict) || v[0].is_symb_of_sommet(at_inferieur_egal)) ){
      v[0]=v[0]._SYMBptr->feuille[1]-v[0]._SYMBptr->feuille[0];
      ok=true;
    }
    if (!ok)
      setsizeerr("Unable to handle when condition "+v[0].print(contextptr));
    return symbolic(at_sign,v[0])*(v[1]-v[2])/2+(v[1]+v[2])/2;
  }
  vector< gen_op_context > when2sign_v(1,whentosign);
  gen when2sign(const gen & g,GIAC_CONTEXT){
    return subst(g,when_v,when2sign_v,false,contextptr);
  }

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC
