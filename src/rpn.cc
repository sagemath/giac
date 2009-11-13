// -*- mode:C++ ; compile-command: "g++-3.4 -I.. -g -c rpn.cc" -*-
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
#include "rpn.h"
#include "symbolic.h"
#include "unary.h"
#include <algorithm>
#include "prog.h"
#include "usual.h"
#include "identificateur.h"
#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#include <dirent.h>
#endif
#include "input_lexer.h"
#include "plot.h"
#include "tex.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC
  bool rpn_mode=false;

  string enmajuscule(const string & s){
    string res;
    string::const_iterator it=s.begin(),itend=s.end();
    for (;it!=itend;++it){
      if ((*it>='a') && (*it<='z'))
        res += *it-32;
      else
        res += *it;
    }
    return res;
  }

  string printasconstant(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
      return sommetstr;
  }
  gen symb_rpn(const gen & args){
    return symbolic(at_rpn,args);
  }
  gen _rpn(const gen & args){
      rpn_mode=true;
      return symb_rpn(args);
  }
  const string _rpn_s("rpn");
  unary_function_unary __rpn(&_rpn,_rpn_s,&printasconstant);
  unary_function_ptr at_rpn (&__rpn);

  gen symb_alg(const gen & args){
    return symbolic(at_alg,args);
  }
  gen _alg(const gen & args){
      rpn_mode=false;
      return symb_alg(args);
  }
  const string _alg_s("alg");
  unary_function_unary __alg(&_alg,_alg_s,&printasconstant);
  unary_function_ptr at_alg (&__alg);

  void roll(int i,vecteur & v){
    if (i<2)
      return;
    iterateur it=v.begin(),itend=v.end();
    if (itend-it<i)
      return;
    it=itend-i;
    gen save=*it;
    for (;;){
      ++it;
      if (it==itend)
	break;
      *(it-1)=*it;
    }
    *(it-1)=save;
  }

  void ROLL(int i,GIAC_CONTEXT){
    roll(i,history_in(contextptr));
    roll(i,history_out(contextptr));
  }

  gen symb_ROLL(const gen & args){
    return symbolic(at_ROLL,args);
  }
  gen _ROLL(const gen & args){
    if (args._VECTptr->empty())
      return args;
    gen e=args._VECTptr->back();
    args._VECTptr->pop_back();
    if (e.type==_INT_)
      roll(e.val,*args._VECTptr);
    if (e.type==_DOUBLE_)
      roll(int(e._DOUBLE_val),*args._VECTptr);
    return args;
  }
  const string _ROLL_s("ROLL");
  unary_function_unary __ROLL(&_ROLL,_ROLL_s);
  unary_function_ptr at_ROLL (&__ROLL);

  void rolld(int i,vecteur & v){
    if (i<2)
      return;
    iterateur it=v.begin(),itend=v.end();
    if (itend-it<i)
      return;
    it=itend-i;
    --itend;
    gen save=*itend;
    for (;it!=itend;){
      --itend;
      *(itend+1)=*itend;
    }
    *it=save;
  }

  void ROLLD(int i,GIAC_CONTEXT){
    rolld(i,history_in(contextptr));
    rolld(i,history_out(contextptr));
  }
  gen symb_ROLLD(const gen & args){
    return symbolic(at_ROLLD,args);
  }
  gen _ROLLD(const gen & args){
    if (args._VECTptr->empty())
      return args;
    gen e=args._VECTptr->back();
    args._VECTptr->pop_back();
    if (e.type==_INT_)
      rolld(e.val,*args._VECTptr);
    if (e.type==_DOUBLE_)
      rolld(int(e._DOUBLE_val),*args._VECTptr);
    return args;
  }
  const string _ROLLD_s("ROLLD");
  unary_function_unary __ROLLD(&_ROLLD,_ROLLD_s);
  unary_function_ptr at_ROLLD (&__ROLLD);

  void stack_swap(vecteur & v){
    iterateur it=v.begin(),itend=v.end();
    int s=itend-it;
    if (s<2)
      return;
    --itend;
    gen tmp=*itend;
    *itend=*(itend-1);
    *(itend-1)=tmp;
  }

  void SWAP(GIAC_CONTEXT){
    stack_swap(history_in(contextptr));
    stack_swap(history_out(contextptr));
  }

  gen symb_SWAP(const gen & args){
    return symbolic(at_SWAP,args);
  }
  gen _SWAP(const gen & args){
    stack_swap(*args._VECTptr);
    return args;
  }
  const string _SWAP_s("SWAP");
  unary_function_unary __SWAP(&_SWAP,_SWAP_s);
  unary_function_ptr at_SWAP (&__SWAP);

  void dup(vecteur & v){
    if (!v.empty())
      v.push_back(v.back());
  }
  gen symb_DUP(const gen & args){
    return symbolic(at_DUP,args);
  }
  gen _DUP(const gen & args){
    dup(*args._VECTptr);
    return args;
  }
  const string _DUP_s("DUP");
  unary_function_unary __DUP(&_DUP,_DUP_s);
  unary_function_ptr at_DUP (&__DUP);

  void over(vecteur & v){
    int s=v.size();
    if (s>=2)
      v.push_back(v[s-2]);
  }
  gen symb_OVER(const gen & args){
    return symbolic(at_OVER,args);
  }
  gen _OVER(const gen & args){
    over(*args._VECTptr);
    return args;
  }
  const string _OVER_s("OVER");
  unary_function_unary __OVER(&_OVER,_OVER_s);
  unary_function_ptr at_OVER (&__OVER);

  void pick(int i,vecteur & v){
    int s=v.size();
    if ((i>=1) && (s>=i))
      v.push_back(v[s-i]);
  }
  gen symb_PICK(const gen & args){
    return symbolic(at_PICK,args);
  }
  gen _PICK(const gen & args){
    if (args._VECTptr->empty())
      return args;
    gen e=args._VECTptr->back();
    args._VECTptr->pop_back();
    if (e.type==_INT_)
      pick(e.val,*args._VECTptr);
    if (e.type==_DOUBLE_)
      pick(int(e._DOUBLE_val),*args._VECTptr);
    return args;
  }
  const string _PICK_s("PICK");
  unary_function_unary __PICK(&_PICK,_PICK_s);
  unary_function_ptr at_PICK (&__PICK);

  void drop(vecteur & v){
    if (v.empty())
      return;
    v.pop_back();
  }
  gen symb_DROP(const gen & args){
    return symbolic(at_DROP,args);
  }
  gen _DROP(const gen & args){
    drop(*args._VECTptr);
    return args;
  }
  const string _DROP_s("DROP");
  unary_function_unary __DROP(&_DROP,_DROP_s);
  unary_function_ptr at_DROP (&__DROP);

  string printasNOP(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    return sommetstr;
  }
  gen symb_NOP(const gen & args){
    return vecteur(1,symbolic(at_NOP,args));
  }
  gen _NOP(const gen & args){
    return args;
  }
  const string _NOP_s("NOP");
  unary_function_unary __NOP(&_NOP,_NOP_s,&printasNOP);
  unary_function_ptr at_NOP (&__NOP);

  string printasIFTE(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    const_iterateur it=feuille._VECTptr->begin();
    string res("IF " + printinner_VECT(*it->_VECTptr,_RPN_FUNC__VECT,contextptr));
    res += " THEN ";
    ++it;
    res += printinner_VECT(*it->_VECTptr,_RPN_FUNC__VECT,contextptr) + " ELSE ";
    ++it;
    return res + printinner_VECT(*it->_VECTptr,_RPN_FUNC__VECT,contextptr)+ " END";
  }
  gen symb_IFTE(const gen & args){
    return symbolic(at_IFTE,args);
  }
  gen _IFTE(const gen & args,const context * contextptr){
    if (args._VECTptr->size()<3)
      return args;
    gen no=args._VECTptr->back();
    args._VECTptr->pop_back();
    gen yes=args._VECTptr->back();
    args._VECTptr->pop_back();
    gen e=args._VECTptr->back();
    args._VECTptr->pop_back();
    if (e.type==_VECT){
      rpn_eval(e,*args._VECTptr,contextptr);
      if (args._VECTptr->empty())
	return args;
      e=args._VECTptr->back();
      args._VECTptr->pop_back();
    }
    if (is_zero(e))
      return rpn_eval(no,*args._VECTptr,contextptr);
    else
      return rpn_eval(yes,*args._VECTptr,contextptr);
  }
  const string _IFTE_s("IFTE");
  unary_function_eval __IFTE(&_IFTE,_IFTE_s,&printasIFTE);
  unary_function_ptr at_IFTE (&__IFTE);

  string printasRPN_LOCAL(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    string s("-> ");
    s += printinner_VECT(*feuille._VECTptr->front()._VECTptr,_RPN_FUNC__VECT,contextptr);
    gen e= feuille._VECTptr->back();
    if ( (e.type==_VECT) && (e.subtype==_RPN_FUNC__VECT))
      return s + " " +e.print(contextptr);
    else { 
      if ( (e.type==_SYMB) && (e._SYMBptr->sommet==at_quote))
	return s + " '"+e._SYMBptr->feuille.print(contextptr)+"'";
      else
	return s + " '"+e.print(contextptr)+"'";
    }
  }
  gen symb_RPN_LOCAL(const gen & a,const gen & b){
    return symbolic(at_RPN_LOCAL,makevecteur(a,b));
  }
  gen _RPN_LOCAL(const gen & args,const context * contextptr) {
    // stack level 2=symbolic names, level 1=program
    if (args.type!=_VECT)
      return symbolic(at_RPN_LOCAL,args);
    int s=args._VECTptr->size();
    if (s<3)
      toofewargs("RPN_LOCAL must have at least 3 args");
    gen prog=args._VECTptr->back();
    args._VECTptr->pop_back();
    vecteur names=*(args._VECTptr->back()._VECTptr); // must be a vector
    args._VECTptr->pop_back();
    // get values from stack
    int nvars=names.size();
    if (s-2<nvars)
      toofewargs("RPN_LOCAL");
    vecteur values(names);
    for (int j=nvars-1;j>=0;--j){
      values[j]=args._VECTptr->back();
      args._VECTptr->pop_back();
    }
    // Initialization
    context * newcontextptr = (context *) contextptr;
    int protect=bind(values,names,newcontextptr);
    vecteur res;
    if ((prog.type==_SYMB) && (prog._SYMBptr->sommet==at_quote)){
      args._VECTptr->push_back(prog._SYMBptr->feuille.eval(eval_level(contextptr),newcontextptr));
      res=*args._VECTptr;
    }
    else
      res=rpn_eval(prog,*args._VECTptr,newcontextptr);
    leave(protect,names,newcontextptr);
    return gen(res,_RPN_STACK__VECT);
  }
  const string _RPN_LOCAL_s("RPN_LOCAL");
  unary_function_eval __RPN_LOCAL(&_RPN_LOCAL,_RPN_LOCAL_s,&printasRPN_LOCAL);
  unary_function_ptr at_RPN_LOCAL (&__RPN_LOCAL);

  string printasRPN_FOR(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if ( (feuille.type!=_VECT) && (feuille._VECTptr->size()!=2))
      return "Invalid_RPN_FOR";
    string s;
    gen controle=feuille._VECTptr->front();
    gen prog=feuille._VECTptr->back();
    bool is_start= (controle[0].print(contextptr)==" j");
    if (is_start)
      s="START ";
    else
      s="FOR "+controle[0].print(contextptr)+" ";
    s += printinner_VECT(*prog._VECTptr,_RPN_FUNC__VECT,contextptr);
    if (is_one(controle[1]))
      s+=" NEXT";
    else
      s+= " "+controle[1].print(contextptr)+" STEP";
    return s;
  }
  gen symb_RPN_FOR(const gen & a,const gen & b){
    return symbolic(at_RPN_FOR,makevecteur(a,b));
  }
  gen _RPN_FOR(const gen & args,const context * contextptr) {
    // stack level 4=init
    // stack level 3=end
    // stack level 2=[init_var,step]
    // level 1=program to execute, 
    if (args.type!=_VECT)
      return symbolic(at_RPN_FOR,args);
    int s=args._VECTptr->size();
    if (s<4)
      toofewargs("RPN_FOR must have at least 4 args");
    gen prog=args._VECTptr->back();
    args._VECTptr->pop_back();
    vecteur controle=*(args._VECTptr->back()._VECTptr); // it must be a vector
    args._VECTptr->pop_back();
    gen fin=args._VECTptr->back();
    args._VECTptr->pop_back();
    gen debut=args._VECTptr->back();
    args._VECTptr->pop_back();
    // Initialization
    vecteur names(1,controle[0]);
    gen test=inferieur_egal(controle[0],fin,contextptr);
    context * newcontextptr = (context *) contextptr;
    int protect=bind(vecteur(1,debut),names,newcontextptr);
    vecteur res;
    for (;!is_zero(test.eval(eval_level(newcontextptr),newcontextptr).evalf(eval_level(contextptr),newcontextptr));sto(eval(controle[0]+controle[1],eval_level(newcontextptr),newcontextptr),controle[0],newcontextptr)){
      res=rpn_eval(prog,*args._VECTptr,newcontextptr);
    }
    leave(protect,names,newcontextptr);
    return gen(res,_RPN_STACK__VECT);
  }
  const string _RPN_FOR_s("RPN_FOR");
  unary_function_eval __RPN_FOR(&_RPN_FOR,_RPN_FOR_s,&printasRPN_FOR);
  unary_function_ptr at_RPN_FOR (&__RPN_FOR);

  string printasRPN_WHILE(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if ( (feuille.type!=_VECT) && (feuille._VECTptr->size()!=2))
      return "Invalid_RPN_WHILE";
    return "WHILE "+ printinner_VECT(*feuille._VECTptr->front()._VECTptr,_RPN_FUNC__VECT,contextptr) + " REPEAT "+printinner_VECT(*feuille._VECTptr->back()._VECTptr,_RPN_FUNC__VECT,contextptr)+ " END ";
  }
  gen symb_RPN_WHILE(const gen & a,const gen & b){
    return symbolic(at_RPN_WHILE,makevecteur(a,b));
  }
  gen _RPN_WHILE(const gen & args,const context * contextptr) {
    // stack level 2=condition
    // level 1=program to execute 
    if (args.type!=_VECT)
      return symbolic(at_RPN_WHILE,args);
    int s=args._VECTptr->size();
    if (s<2)
      toofewargs("RPN_WHILE must have at least 2 args");
    gen prog=args._VECTptr->back();
    args._VECTptr->pop_back();
    gen controle=args._VECTptr->back();
    args._VECTptr->pop_back();
    vecteur res;
    for (;;){
      res=rpn_eval(controle,*args._VECTptr,contextptr);
      if (args._VECTptr->empty())
	toofewargs("WHILE");
      gen tmp=args._VECTptr->back();
      args._VECTptr->pop_back();
      if (is_zero(tmp.eval(1,contextptr).evalf(1,contextptr)))
	break;
      res=rpn_eval(prog,*args._VECTptr,contextptr);
    }
    return gen(res,_RPN_STACK__VECT);
  }
  const string _RPN_WHILE_s("RPN_WHILE");
  unary_function_eval __RPN_WHILE(&_RPN_WHILE,_RPN_WHILE_s,&printasRPN_WHILE);
  unary_function_ptr at_RPN_WHILE (&__RPN_WHILE);

  string printasRPN_CASE(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if (feuille.type!=_VECT)
      return "Invalid_RPN_CASE";
    vecteur v=*feuille._VECTptr;
    if ((v.size()!=1) ||(v.front().type!=_VECT))
      return "Invalid_RPN_CASE";
    string res("CASE ");
    const_iterateur it=v.front()._VECTptr->begin(),itend=v.front()._VECTptr->end();
    for (;it!=itend;){
      res += printinner_VECT(*it->_VECTptr,_RPN_FUNC__VECT,contextptr);
      ++it;
      if (it==itend)
	break;
      res += " THEN " + printinner_VECT(*it->_VECTptr,_RPN_FUNC__VECT,contextptr) + " END ";
      ++it;
    }
    return res+" END ";
  }
  gen symb_RPN_CASE(const gen & a){
    return symbolic(at_RPN_CASE,vecteur(1,a));
  }
  gen _RPN_CASE(const gen & args,const context * contextptr) {
    // level 1=[case1, prg1, case2,prg2,..., [default]]
    if (args.type!=_VECT)
      return symbolic(at_RPN_CASE,args);
    int s=args._VECTptr->size();
    if (s<1)
      toofewargs("RPN_CASE must have at least 1 arg");
    vecteur controle=*args._VECTptr->back()._VECTptr;
    args._VECTptr->pop_back();
    const_iterateur it=controle.begin(),itend=controle.end();
    vecteur res;
    for (;it!=itend;){
      res=rpn_eval(*it,*args._VECTptr,contextptr);
      if (args._VECTptr->empty())
	toofewargs("CASE");
      ++it;
      if (it==itend) // default of the case struct
	break;
      gen test=args._VECTptr->back();
      args._VECTptr->pop_back();
      if (!(is_zero(test.eval(1,contextptr).evalf(1,contextptr)))){
	res=rpn_eval(*it,*args._VECTptr,contextptr);
	break;
      }
      ++it;
    }
    return gen(*args._VECTptr,_RPN_STACK__VECT);
  }
  const string _RPN_CASE_s("RPN_CASE");
  unary_function_eval __RPN_CASE(&_RPN_CASE,_RPN_CASE_s,&printasRPN_CASE);
  unary_function_ptr at_RPN_CASE (&__RPN_CASE);

  gen symb_RCL(const gen & a){
    return symbolic(at_RCL,a);
  }
  gen _RCL(const gen & args,const context * contextptr) {
    // stack level 2=condition
    // level 1=program to execute 
    if (args.type!=_IDNT)
      return symbolic(at_RCL,args);
    return args._IDNTptr->eval(1,args,contextptr);
  }
  const string _RCL_s("RCL");
  unary_function_eval __RCL(&_RCL,_RCL_s);
  unary_function_ptr at_RCL (&__RCL);

  gen symb_VARS(const gen & a){
    return symbolic(at_VARS,a);
  }
#if defined(__APPLE__) || defined(__FreeBSD__)
  static int int_one (struct dirent *unused){
#else
  static int int_one (const struct dirent *unused){
#endif
    return 1;
  }
  gen _VARS(const gen & args,const context * contextptr) {
    vecteur res;
    if (contextptr){
      if (contextptr->globalcontextptr && contextptr->globalcontextptr->tabptr){
	sym_tab::const_iterator it=contextptr->globalcontextptr->tabptr->begin(),itend=contextptr->globalcontextptr->tabptr->end();
	for (;it!=itend;++it)
	  res.push_back(identificateur(it->first));
      }
      return res;
    }
    if (!variables_are_files(contextptr)){
      sym_tab::const_iterator it=syms().begin(),itend=syms().end();
      for (;it!=itend;++it){
	gen id=it->second;
	if (id.type==_IDNT && id._IDNTptr->value){
	  res.push_back(id);
	}
      }
      if (is_one(args) && current_folder_name.type==_IDNT && current_folder_name._IDNTptr->value){ // add variables of the current folder
	gen & tmp=*current_folder_name._IDNTptr->value;
	if (tmp.type==_VECT){
	  vecteur v=*current_folder_name._IDNTptr->value->_VECTptr;
	  iterateur it=v.begin(),itend=v.end();
	  for (;it!=itend;++it){
	    if (it->type!=_VECT || it->_VECTptr->size()!=2)
	      continue;
	    vecteur & w=*it->_VECTptr;
	    res.push_back(w[0]);
	  }
	}
      }
      return res;
    }
#ifdef VISUALC
    return undef;
#else
    struct dirent **eps;
    int n;
    n = scandir ("./", &eps, int_one, alphasort);
    if (n >= 0){
      string s;
      int cnt;
      for (cnt = 0; cnt < n; ++cnt){
	s=string(eps[cnt]->d_name);
	unsigned k=0;
	for (;k<s.size();++k){
	  if (!isalphan(s[k]))
	    break;
	}
	if ( (k==s.size()-4) && (s[k]=='.') && (s.substr(k+1,3)=="cas") )
	  res.push_back(identificateur(s.substr(0,k)));
      }
      if ( rpn_mode && (args.type==_VECT) && (args.subtype==_RPN_STACK__VECT) ){
	args._VECTptr->push_back(res);
	return args;
      }
      return res;
    }
    else
      settypeerr ("Couldn't open the directory");
    return 0;
#endif
  }

  const string _VARS_s("VARS");
  unary_function_eval __VARS(&_VARS,_VARS_s);
  unary_function_ptr at_VARS (&__VARS);

  gen symb_purge(const gen & a){
    return symbolic(at_purge,a);
  }
  gen _purge(const gen & args,const context * contextptr) {
    if (rpn_mode && (args.type==_VECT)){
      if (!args._VECTptr->size())
	toofewargs("purge");
      gen apurger=args._VECTptr->back();
      _purge(apurger,contextptr);
      args._VECTptr->pop_back();
      return gen(*args._VECTptr,_RPN_STACK__VECT);
    }
    if (args.type==_VECT)
      return apply(args,contextptr,giac::_purge);
    if (args.is_symb_of_sommet(at_at)){
      gen & f = args._SYMBptr->feuille;
      if (f.type==_VECT && f._VECTptr->size()==2){
	gen m = eval(f._VECTptr->front(),eval_level(contextptr),contextptr);
	gen indice=eval(f._VECTptr->back(),eval_level(contextptr),contextptr);
	if (m.type==_MAP){
	  gen_map::iterator it=m._MAPptr->find(indice),itend=m._MAPptr->end();
	  if (it==itend)
	    setsizeerr("Bad index"+indice.print(contextptr));
	  m._MAPptr->erase(it);
	  return 1;
	}
      }
    }
    if (args.type!=_IDNT)
      return symbolic(at_purge,args);
    // REMOVED! args.eval(eval_level(contextptr),contextptr); 
    if (contextptr){
      if (contextptr->globalcontextptr!=contextptr){ 
	// purge a local variable = set it to assume(DOM_SYMBOLIC)
	gen a2(_SYMB);
	a2.subtype=1;
	return sto(gen(makevecteur(a2),_ASSUME__VECT),args,contextptr);
      }
      // purge a global variable
      sym_tab::iterator it=contextptr->tabptr->find(*args._IDNTptr->name),itend=contextptr->tabptr->end();
      if (it==itend)
	return string2gen("No such variable "+args.print(contextptr),false);
      gen res=it->second;
      if (it->second.type==_POINTER_ && it->second.subtype==_THREAD_POINTER)
	settypeerr(args.print(contextptr)+" is locked by thread "+it->second.print(contextptr));
      contextptr->tabptr->erase(it);
      return res;
    }
    if (current_folder_name.type==_IDNT && current_folder_name._IDNTptr->value && current_folder_name._IDNTptr->value->type==_VECT){
      vecteur v=*current_folder_name._IDNTptr->value->_VECTptr;
      iterateur it=v.begin(),itend=v.end();
      gen val;
      for (;it!=itend;++it){
	if (it->type!=_VECT || it->_VECTptr->size()!=2)
	  continue;
	vecteur & w=*it->_VECTptr;
	if (w[0]==args){
	  val=w[1];
	  break;
	}
      }
      if (it!=itend){
	v.erase(it);
	gen res=gen(v,_FOLDER__VECT);
	*current_folder_name._IDNTptr->value=res;
	if (!child_id && signal_store)
	  _signal(symb_quote(symbolic(at_sto,makevecteur(res,current_folder_name))),contextptr);
	return val;
      }
    }
    if (args._IDNTptr->value){
      if (variables_are_files(contextptr))
	unlink((*args._IDNTptr->name+cas_suffixe).c_str());
      gen res=*args._IDNTptr->value;
      if (res.type==_VECT && res.subtype==_FOLDER__VECT){
	if (res._VECTptr->size()!=1)
	  setsizeerr("Non-empty folder");
      }
      delete args._IDNTptr->value;
      args._IDNTptr->value=0;
      if (!child_id && signal_store)
	_signal(symb_quote(symb_purge(args)),contextptr);
      return res;
    }
    else
      return string2gen(args.print(contextptr)+" not assigned",false);
  }
  const string _purge_s("purge");
  unary_function_eval __purge(&_purge,_purge_s);
  unary_function_ptr at_purge (&__purge,_QUOTE_ARGUMENTS);

  string printasRPN_UNTIL(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if ( (feuille.type!=_VECT) && (feuille._VECTptr->size()!=2))
      return "Invalid_RPN_UNTIL";
    return "DO "+ printinner_VECT(*feuille._VECTptr->front()._VECTptr,_RPN_FUNC__VECT,contextptr) + " UNTIL "+printinner_VECT(*feuille._VECTptr->back()._VECTptr,_RPN_FUNC__VECT,contextptr)+ " END ";
  }
  gen symb_RPN_UNTIL(const gen & a,const gen & b){
    return symbolic(at_RPN_UNTIL,makevecteur(a,b));
  }
  gen _RPN_UNTIL(const gen & args,const context * contextptr) {
    // stack level 2=program
    // level 1=condition 
    if (args.type!=_VECT)
      return symbolic(at_RPN_UNTIL,args);
    int s=args._VECTptr->size();
    if (s<2)
      toofewargs("RPN_UNTIL must have at least 2 args");
    gen controle=args._VECTptr->back();
    args._VECTptr->pop_back();
    gen prog=args._VECTptr->back();
    args._VECTptr->pop_back();
    vecteur res;
    for (;;){
      res=rpn_eval(prog,*args._VECTptr,contextptr);
      res=rpn_eval(controle,*args._VECTptr,contextptr);
      if (args._VECTptr->empty())
	toofewargs("UNTIL");
      gen tmp=args._VECTptr->back();
      args._VECTptr->pop_back();
      if (!is_zero(tmp.eval(eval_level(contextptr),contextptr).evalf(eval_level(contextptr),contextptr)))
	break;
    }
    return gen(res,_RPN_STACK__VECT);
  }
  const string _RPN_UNTIL_s("RPN_UNTIL");
  unary_function_eval __RPN_UNTIL(&_RPN_UNTIL,_RPN_UNTIL_s,&printasRPN_UNTIL);
  unary_function_ptr at_RPN_UNTIL (&__RPN_UNTIL);

  // RPN evaluation loop (_VECTEVAL), no return stack currently
  vecteur rpn_eval(const vecteur & prog,vecteur & pile,const context * contextptr){
    const_iterateur it=prog.begin(),itend=prog.end();
    for (;it!=itend;++it){
      if (it->type==_FUNC){
	// test nargs with subtype
	int nargs=it->subtype;
	if (nargs>signed(pile.size()))
	  toofewargs(it->print(contextptr)+": stack "+gen(pile).print(contextptr));
	if (nargs==1){
	  pile.back()=(*it->_FUNCptr)(pile.back(),contextptr);
	}
	else {
	  if (nargs){
	    vecteur v(nargs);
	    for (int k=nargs-1;k>=0;--k){
	      v[k]=pile.back();
	      pile.pop_back();
	    }
	    pile.push_back((*it->_FUNCptr)(v,contextptr));
	  }
	  else {
	    gen res;
	    if (*it->_FUNCptr==at_eval){ // eval stack level 1
	      if (pile.empty())
		toofewargs("EVAL");
	      res=pile.back();
	      pile.pop_back();
	      if ( (res.type==_SYMB) && (res._SYMBptr->sommet==at_rpn_prog))
		res=res._SYMBptr->feuille;
	      res=rpn_eval(res,pile,contextptr);
	    }
	    else
	      res=(*it->_FUNCptr)(pile,contextptr);
	    if ( (res.type==_VECT) && (res.subtype=_RPN_STACK__VECT) )
	      pile= *res._VECTptr;
	    else
	      pile= vecteur(1,res);
	  }
	}
      }
      else {
	// test for special symbolic (control struct)
	unary_function_ptr control_op[]={at_RPN_LOCAL,at_RPN_FOR,at_IFTE,at_RPN_CASE,at_RPN_WHILE,at_RPN_UNTIL,0};
	if ( (it->type==_SYMB) && equalposcomp(control_op,it->_SYMBptr->sommet)){
	  // push args of it to the stack and call sommet on the stack
	  if (it->_SYMBptr->feuille.type!=_VECT) // should not happen!
	    pile.push_back(it->_SYMBptr->feuille);
	  else
	    pile=mergevecteur(pile,*it->_SYMBptr->feuille._VECTptr);
	  gen res=it->_SYMBptr->sommet(pile,contextptr);
	  if ( (res.type==_VECT) && (res.subtype=_RPN_STACK__VECT) )
	    pile= *res._VECTptr;
	  else
	    pile= vecteur(1,res);	  
	}
	else {
	  if ( (it->type!=_VECT) 
	       // || (it->subtype==_RPN_FUNC__VECT) 
	       ){
	    gen res=it->eval(1,contextptr);
	    if ( (res.type==_VECT) && (res.subtype==_RPN_STACK__VECT))
	      pile=*res._VECTptr;
	    else
	      pile.push_back(res);
	  }
	  else
	    pile.push_back(*it);
	}
      }
    }
    return pile;
  }

  vecteur rpn_eval(const gen & prog,vecteur & pile,const context * contextptr){
    if (prog.type!=_VECT)
      return rpn_eval(vecteur(1,prog),pile,contextptr);
    else
      return rpn_eval(*prog._VECTptr,pile,contextptr);
  }

  string printasrpn_prog(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if (feuille.type!=_VECT)
      return "<< "+feuille.print(contextptr)+" >>";
    return "<< "+printinner_VECT(*feuille._VECTptr,_RPN_FUNC__VECT,contextptr)+" >>";
  }
  gen symb_rpn_prog(const gen & args){
    return symbolic(at_rpn_prog,args);
  }
  gen _rpn_prog(const gen & args,const context * contextptr){
    if (!rpn_mode || (args.type!=_VECT))
      return symbolic(at_rpn_prog,args);
    vecteur pile(history_out(contextptr));
    *logptr(contextptr) << pile << " " << args << endl;
    return gen(rpn_eval(*args._VECTptr,pile,contextptr),_RPN_STACK__VECT);
  }
  const string _rpn_prog_s("rpn_prog");
  unary_function_eval __rpn_prog(&_rpn_prog,_rpn_prog_s,&printasrpn_prog);
  unary_function_ptr at_rpn_prog (&__rpn_prog,_QUOTE_ARGUMENTS);

  string texprintasdivision(const gen & feuille,const string & s,GIAC_CONTEXT){
    if (feuille.type!=_VECT || feuille._VECTptr->size()!=2)
      return "invalid /";
    return "\\frac{"+gen2tex(feuille._VECTptr->front(),contextptr)+"}{"+gen2tex(feuille._VECTptr->back(),contextptr)+"}";
  }
  gen symb_division(const gen & a,const gen & b){
    return symbolic(at_division,makevecteur(a,b));
  }
  gen symb_division(const gen & args){
    return symbolic(at_division,args);
  }
  gen _division(const gen & args){
    if ( (args.type!=_VECT) || (args._VECTptr->size()!=2) )
      return symb_division(args);
    return rdiv(args._VECTptr->front(),args._VECTptr->back());
  }
  const string _division_s("/");
  unary_function_unary __division(&_division,_division_s,&printsommetasoperator,&texprintasdivision);
  unary_function_ptr at_division (&__division);

  gen symb_binary_minus(const gen & a,const gen & b){
    return symbolic(at_binary_minus,makevecteur(a,b));
  }
  gen symb_binary_minus(const gen & args){
    return symbolic(at_binary_minus,args);
  }
  gen _binary_minus(const gen & args){
    if ( (args.type!=_VECT) || (args._VECTptr->size()!=2) )
      return symb_binary_minus(args);
    return args._VECTptr->front()-args._VECTptr->back();
  }
  const string _binary_minus_s("-");
  unary_function_unary __binary_minus(&_binary_minus,_binary_minus_s,&printsommetasoperator,&texprintsommetasoperator);
  unary_function_ptr at_binary_minus (&__binary_minus);

  vecteur tab2vecteur(gen tab[]){
    vecteur res;
    for (;!is_zero(*tab);++tab)
      res.push_back(*tab);
    return res;
  }

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC

