// -*- mode:C++ ; compile-command: "g++ -I.. -g -c unary.cc" -*-
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
#include "unary.h"
#include "gen.h"
#include "usual.h"
#include "rpn.h"
#include "tex.h"
#include "input_lexer.h"
#include "symbolic.h"
#include "input_lexer.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC
  // unary_function_ptr
  gen unary_function_ptr::operator () (const gen & arg,const context * context_ptr) const{
    return (*ptr)(arg,context_ptr);
    this->dbgprint();
  }

  void unary_function_ptr::dbgprint() const {
    cerr << ptr->s << endl; 
  }

  string unary_function_abstract::print(GIAC_CONTEXT) const { 
    int lang=language(contextptr);
    multimap<string,localized_string>::iterator it=back_lexer_localization_map.find(s),backend=back_lexer_localization_map.end(),itend=back_lexer_localization_map.upper_bound(s);
    if (it!=backend){
      for (;it!=itend;++it){
	if (it->second.language==lang)
	  return it->second.chaine;
      }
    }
    return s;
    /*
    if (rpn_mode) 
      return enmajuscule(s); 
    else 
      return s; 
    */
  }

  unary_function_ptr::unary_function_ptr(const unary_function_ptr & myptr):ptr(myptr.ptr),ref_count(myptr.ref_count),quoted(myptr.quoted){
    if (ref_count)
      ++(*ref_count);
  }

  // dynamic unary_function_abstract
  unary_function_ptr::unary_function_ptr(const unary_function_abstract & myptr): quoted(0) {
    ref_count = new int(1);
    ptr = myptr.recopie();
  }

  unary_function_ptr::unary_function_ptr(const unary_function_abstract & myptr,int myquoted,int parser_token): quoted(myquoted) {
    ref_count = new int(1);
    ptr = myptr.recopie();
    if (parser_token)
      if (!lexer_functions_register(*this,ptr->s,parser_token))
	setsizeerr("Unable to register "+ptr->s);
  }
  // global unary_function_abstract pointer, no reference count here
  unary_function_ptr::unary_function_ptr(const unary_function_abstract * myptr): ptr(myptr),ref_count(0),quoted(0){}

  unary_function_ptr::unary_function_ptr(const unary_function_abstract * myptr,int myquoted,int parser_token): ptr(myptr),ref_count(0),quoted(myquoted){
    if (parser_token)
      if (!lexer_functions_register(*this,ptr->s,parser_token))
	setsizeerr("Unable to register "+ptr->s);
  }

  // copy constructor
  unary_function_ptr & unary_function_ptr::operator = (const unary_function_ptr & acopier){
    if (ref_count){
      if (!((*ref_count)--)){
	delete ptr;
	delete ref_count;
      }
    }
    ptr = acopier.ptr;
    ref_count=acopier.ref_count;
    quoted = acopier.quoted;
    if (ref_count)
      ++(*ref_count);
    return *this;
  }

  unary_function_ptr::~unary_function_ptr(){
    if (ref_count){
      if (!((*ref_count)--)){
	delete ptr;
	delete ref_count;
      }
    }
  }

  gen apply(const gen & e,const unary_function_ptr & f){
    if (e.type!=_VECT)
      return f(e,0);
    const_iterateur it=e._VECTptr->begin(),itend=e._VECTptr->end();
    vecteur v;
    v.reserve(itend-it);
    for (;it!=itend;++it)
      v.push_back(f(*it,0));
    return gen(v,e.subtype);
  }

  gen apply(const gen & e, gen (* f) (const gen &) ){
    if (e.type!=_VECT)
      return f(e);
    const_iterateur it=e._VECTptr->begin(),itend=e._VECTptr->end();
    vecteur v;
    v.reserve(itend-it);
    for (;it!=itend;++it)
      v.push_back(f(*it));
    return gen(v,e.subtype);
  }

  gen apply(const gen & e, gen (* f) (const gen &,const context *),GIAC_CONTEXT ){
    if (e.type!=_VECT)
      return f(e,contextptr);
    const_iterateur it=e._VECTptr->begin(),itend=e._VECTptr->end();
    vecteur v;
    v.reserve(itend-it);
    for (;it!=itend;++it)
      v.push_back(f(*it,contextptr));
    return gen(v,e.subtype);
  }

  gen apply(const gen & e1, const gen & e2,gen (* f) (const gen &, const gen &) ){
    if ((e1.type!=_VECT) && (e2.type!=_VECT))
      return f(e1,e2);
    if (e1.type!=_VECT){
      const_iterateur it=e2._VECTptr->begin(),itend=e2._VECTptr->end();
      vecteur v;
      v.reserve(itend-it);
      for (;it!=itend;++it)
	v.push_back(f(e1,*it));
      return gen(v,e2.subtype);
    }
    if (e2.type!=_VECT){
      const_iterateur it=e1._VECTptr->begin(),itend=e1._VECTptr->end();
      vecteur v;
      v.reserve(itend-it);
      for (;it!=itend;++it)
	v.push_back(f(*it,e2));
      return gen(v,e1.subtype);
    }
    const_iterateur it1=e1._VECTptr->begin(),it1end=e1._VECTptr->end();
    const_iterateur it2=e2._VECTptr->begin(),it2end=e2._VECTptr->end();
    if (it2end-it2!=it1end-it1)
      setdimerr();
    vecteur v;
    v.reserve(it1end-it1);
    for (;it1!=it1end;++it1,++it2)
      v.push_back(f(*it1,*it2));
    return gen(v,e1.subtype);
  }

  gen apply(const gen & e, const context * contextptr,const gen_op_context & f ){
    if (e.type!=_VECT)
      return f(e,contextptr);
    const_iterateur it=e._VECTptr->begin(),itend=e._VECTptr->end();
    vecteur v;
    v.reserve(itend-it);
    for (;it!=itend;++it)
      v.push_back(f(*it,contextptr));
    return gen(v,e.subtype);
  }

  gen apply(const gen & e1, const gen & e2,const context * contextptr,gen (* f) (const gen &, const gen &,const context *) ){
    if ((e1.type!=_VECT) && (e2.type!=_VECT))
      return f(e1,e2,contextptr);
    if (e1.type!=_VECT){
      const_iterateur it=e2._VECTptr->begin(),itend=e2._VECTptr->end();
      vecteur v;
      v.reserve(itend-it);
      for (;it!=itend;++it)
	v.push_back(f(e1,*it,contextptr));
      return gen(v,e2.subtype);
    }
    if (e2.type!=_VECT){
      const_iterateur it=e1._VECTptr->begin(),itend=e1._VECTptr->end();
      vecteur v;
      v.reserve(itend-it);
      for (;it!=itend;++it)
	v.push_back(f(*it,e2,contextptr));
      return gen(v,e1.subtype);
    }
    const_iterateur it1=e1._VECTptr->begin(),it1end=e1._VECTptr->end();
    const_iterateur it2=e2._VECTptr->begin(),it2end=e2._VECTptr->end();
    if (it2end-it2!=it1end-it1)
      setdimerr();
    vecteur v;
    v.reserve(it1end-it1);
    for (;it1!=it1end;++it1,++it2)
      v.push_back(f(*it1,*it2,contextptr));
    return gen(v,e1.subtype);
  }

  gen apply1st(const gen & e1, const gen & e2,gen (* f) (const gen &, const gen &) ){
    if (e1.type!=_VECT)
      return f(e1,e2);
    const_iterateur it=e1._VECTptr->begin(),itend=e1._VECTptr->end();
    vecteur v;
    v.reserve(itend-it);
    for (;it!=itend;++it)
      v.push_back(f(*it,e2));
    return gen(v,e1.subtype);
  }

  gen apply1st(const gen & e1, const gen & e2,const context * contextptr, gen (* f) (const gen &, const gen &,const context *) ){
    if (e1.type!=_VECT)
      return f(e1,e2,contextptr);
    const_iterateur it=e1._VECTptr->begin(),itend=e1._VECTptr->end();
    vecteur v;
    v.reserve(itend-it);
    for (;it!=itend;++it)
      v.push_back(f(*it,e2,contextptr));
    return gen(v,e1.subtype);
  }

  gen apply2nd(const gen & e1, const gen & e2,gen (* f) (const gen &, const gen &) ){
    if (e2.type!=_VECT)
      return f(e1,e2);
    const_iterateur it=e2._VECTptr->begin(),itend=e2._VECTptr->end();
    vecteur v;
    v.reserve(itend-it);
    for (;it!=itend;++it)
      v.push_back(f(e1,*it));
    return gen(v,e2.subtype);
  }

  gen apply2nd(const gen & e1, const gen & e2,const context * contextptr, gen (* f) (const gen &, const gen &,const context *) ){
    if (e2.type!=_VECT)
      return f(e1,e2,contextptr);
    const_iterateur it=e2._VECTptr->begin(),itend=e2._VECTptr->end();
    vecteur v;
    v.reserve(itend-it);
    for (;it!=itend;++it)
      v.push_back(f(e1,*it,contextptr));
    return gen(v,e2.subtype);
  }

  
  // unary_function_abstract
  unary_function_abstract * unary_function_abstract::recopie() const{
    unary_function_abstract * ptr=new unary_function_abstract(s);
    ptr->D = D;
    return ptr;
  }

  unary_function_unary * unary_function_unary::recopie() const{
    unary_function_unary * ptr=new unary_function_unary(op,s);
    ptr->D = D;
    return ptr;
  }

  unary_function_eval * unary_function_eval::recopie() const{
    unary_function_eval * ptr=new unary_function_eval(op,s);
    ptr->D = D;
    return ptr;
  }

  unary_function_compose * unary_function_compose::recopie() const{
    unary_function_compose * ptr=new unary_function_compose(op_v);
    ptr->D = D;
    return ptr;
  }

  unary_function_list * unary_function_list::recopie() const{
    unary_function_list * ptr=new unary_function_list(op_l);
    ptr->D = D;
    return ptr;
  }

  unary_function_constant * unary_function_constant::recopie() const{
    unary_function_constant * ptr=new unary_function_constant(constant);
    ptr->D = D;
    return ptr;
  }

  unary_function_innerprod * unary_function_innerprod::recopie() const{
    unary_function_innerprod * ptr=new unary_function_innerprod(i);
    ptr->D = D;
    return ptr;
  }
  
  unary_function_user * unary_function_user::recopie() const{
    unary_function_user * ptr=new unary_function_user(f,s,printsommet,texprint,cprint);
    ptr->D = D;
    return ptr;
  }
  

  // unary_function_compose related
  gen unary_function_compose::operator () (const gen & arg,const context * context_ptr) const{
    vector<unary_function_ptr>::const_iterator it=op_v.begin(),itend=op_v.end();
    gen res(arg);
    for (;it!=itend;++it)
      res=(*it)(res,context_ptr);
    return res;
  }

  unary_function_compose::unary_function_compose(const vector<unary_function_ptr> & myop_v) : unary_function_abstract(":: ") {
    vector<unary_function_ptr>::const_iterator it=myop_v.begin(),itend=myop_v.end();
    for (;it!=itend;++it){
      op_v.push_back( (*it) );
      s += it->ptr->s;
      s +=" ";
    }
    s += string(";");
  }

  gen unary_function_list::operator () (const gen & arg,const context * context_ptr) const{
    vector<unary_function_ptr>::const_iterator it=op_l.begin(),itend=op_l.end();
    vecteur res;
    for (;it!=itend;++it)
      res.push_back( (*it)(arg,context_ptr));
    return res;
  }

  unary_function_list::unary_function_list(const vector<unary_function_ptr> & myop_v) : unary_function_abstract("{ ") {
    vector<unary_function_ptr>::const_iterator it=myop_v.begin(),itend=myop_v.end();
    for (;it!=itend;++it){
      op_l.push_back( (*it) );
      s += it->ptr->s;
      s +=" ";
    }
    s += string("}");
  }


  gen unary_function_innerprod::operator () (const gen & arg,const context * contextptr) const{
    if (arg.type!=_VECT)
      setsizeerr(arg.print(contextptr)+ " should be of type _VECT (unary.cc)");
    vecteur res;
    // remove i indices from arg
    vecteur::const_iterator jt=arg._VECTptr->begin(),jtend=arg._VECTptr->end();
    vector<int>::const_iterator it=i.begin(), itend=i.end();
    for (int j=0;jt!=jtend;++jt,++j){
      if (it==itend)
	break;
      else {
	if (j!=*it)
	  res.push_back(*jt);
	else
	  ++it;
      }
    }
    for (;jt!=jtend;++jt)
      res.push_back(*jt);
    return res;
  }

  // I/O
  ostream & operator << (ostream & os,const unary_function_abstract & o){ return os << o.s; }
  ostream & operator << (ostream & os,const unary_function_unary & o) { return os << o.s ; }
  ostream & operator << (ostream & os,const unary_function_eval & o) { return os << o.s ; }
  ostream & operator << (ostream & os,const unary_function_compose & p){ return os << p.s;} 
  ostream & operator << (ostream & os,const unary_function_list & p){ return os<< p.s; }
  ostream & operator << (ostream & os,const unary_function_constant & c){ return os<< c.s; }
  ostream & operator << (ostream & os,const unary_function_innerprod & i){ return os<< i.s; }


  string printsommetasoperator(const gen & feuille,const string & sommetstr_orig,GIAC_CONTEXT){
    if (feuille.type!=_VECT)
      return feuille.print(contextptr);
    string sommetstr(sommetstr_orig);
    if (isalpha(sommetstr[0]) || sommetstr[0]=='%')
      sommetstr=' '+sommetstr+' ';
    if (sommetstr==_pow_s) {
      gen pui=feuille._VECTptr->back();
      gen arg=feuille._VECTptr->front();
      if ( pui==plus_one_half )
	return "sqrt("+arg.print(contextptr)+')';
      if ( pui==minus_one_half  || pui==fraction(minus_one,plus_two) )
	return "1/sqrt("+arg.print(contextptr)+')';
      if (arg.type==_SYMB && arg._SYMBptr->sommet!=at_neg && !arg._SYMBptr->sommet.ptr->printsommet){
	bool puisymb=pui.type==_SYMB;
	sommetstr=arg.print(contextptr)+"^";
	if (puisymb)
	  sommetstr += "(";
	sommetstr += pui.print(contextptr);
	if (puisymb)
	  sommetstr += ")";
	return sommetstr;
      }
    }
    vecteur::const_iterator itb=feuille._VECTptr->begin(),itend=feuille._VECTptr->end();
    if (itb==itend)
      return "";
    string s;
    if (itb->type==_FRAC && sommetstr!="=")
      s='('+itb->print(contextptr)+")";
    else {
      if ( sommetstr=="=" || itb->type==_IDNT || (itb->type<=_CPLX && is_positive(*itb,contextptr)) )
	s=itb->print(contextptr);
      else
	s='('+itb->print(contextptr)+")";
    }
    ++itb;
    for (;;){
      if (itb==itend)
	return s;
      if ( itb->type==_SYMB || itb->type==_FRAC || itb->type==_CPLX || (itb->type==_VECT && itb->subtype==_SEQ__VECT) )
	s += sommetstr + '('+itb->print(contextptr)+")";
      else
	s += sommetstr + itb->print(contextptr);
      ++itb;
    }
  }
    
  string texprintsommetasoperator(const gen & feuille,const string & sommetstr,GIAC_CONTEXT){
    if (feuille.type!=_VECT)
      return feuille.print(contextptr);
    vecteur::const_iterator itb=feuille._VECTptr->begin(),itend=feuille._VECTptr->end();
    if (itb==itend)
      return "";
    string s;
    if (itb->type==_FRAC)
      s="("+gen2tex(*itb,contextptr)+")";
    else {
      if ( sommetstr=="=" || itb->type==_IDNT || (itb->type<=_CPLX && is_positive(*itb,contextptr)) )
	s=gen2tex(*itb,contextptr);
      else
	s="("+gen2tex(*itb,contextptr)+")";
    }
    ++itb;
    for (;;){
      if (itb==itend)
	return s;
      if ( itb->type==_SYMB || itb->type==_FRAC || itb->type==_CPLX || (itb->type==_VECT && itb->subtype==_SEQ__VECT) )
	s += sommetstr + '('+gen2tex(*itb,contextptr)+")";
      else
	s += sommetstr + gen2tex(*itb,contextptr);
      ++itb;
    }
  }

  partial_derivative_onearg::partial_derivative_onearg(gen (* mydf) (const gen & args) ) :     df(unary_function_ptr(unary_function_unary(mydf,""))) {}

  partial_derivative_onearg::partial_derivative_onearg(gen (* mydf) (const gen & args,const context * contextptr) ) :     df(unary_function_ptr(unary_function_eval(mydf,""))) {}

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC
