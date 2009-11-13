// -*- mode:C++ ; compile-command: "g++ -I.. -g -c identificateur.cc" -*-
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
#ifndef _GIAC_IDENTIFICATEUR_H
#define _GIAC_IDENTIFICATEUR_H
#include "first.h"
#include <string>
#include <iostream>
#include "global.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  extern int protection_level; // for local vars

  const int MAXLENSIZE = 1000000; // max size of a line in files

  class gen;

  class identificateur {
  public:
    int * ref_count;
    gen * value;
    std::string * name;
    vecteur * localvalue;
    // value / localvalue might be an assumption if it's a vecteur 
    // of subtype _ASSUME__VECT
    // The first gen of an assumption vecteur is the type (_FRAC for rational)
    // If the type is _REAL, the vecteur has 2 other elements
    // * an interval or a _SET_VECT of intervals 
    //   where interval=vecteur of length 2 of subtype _LINE__VECT
    // * a list of excluded particular values
    // If the type is _DOUBLE_ the variable will be evalf-ed but not eval-ed
    // This is useful in geometry to make figures and get exact results
    // If the type is _INT_ it 
    bool * quoted;
    identificateur();
    explicit identificateur(const std::string & s);
    explicit identificateur(const char * s);
    identificateur(const std::string & s,const gen & e);
    identificateur(const identificateur & s);
    ~identificateur();
    identificateur & operator =(const identificateur & s);
    gen eval(int level,const gen & orig,const context * context_ptr) ;
    bool in_eval(int level,const gen & orig,gen & evaled,const context * context_ptr) ;
    std::string print(const context * context_ptr) const ;
    void dbgprint() const { std::cout << this->print(context0); }
    void unassign() ;
    void push(int protection,const gen & e);
  };

  // make g identificateurs evaluated as global in null context
  gen global_eval(const gen & g,int level);
  gen global_evalf(const gen & g,int level);
  // return the local value of i, if globalize is true, replace idnt with
  // global idnt in returned value
  gen do_local_eval(const identificateur & i,int level,bool globalize);

  std::ostream & operator << (std::ostream & os,const identificateur & s);
  
  extern std::string string_euler_gamma;
  extern identificateur _IDNT_euler_gamma;
  extern gen cst_euler_gamma;
  extern std::string string_pi;
  extern identificateur _IDNT_pi;
  extern gen cst_pi;
  extern std::string string_infinity;
  identificateur & _IDNT_infinity();
  extern gen unsigned_inf;
  extern std::string string_undef;
  identificateur & _IDNT_undef();
  extern gen undef;
  extern identificateur a__IDNT;
  extern gen a__IDNT_e;
  extern identificateur b__IDNT;
  extern gen b__IDNT_e;
  extern identificateur c__IDNT;
  extern gen c__IDNT_e;
  extern identificateur d__IDNT;
  extern gen d__IDNT_e;
  extern identificateur e__IDNT;
  extern gen e__IDNT_e;
  extern identificateur f__IDNT;
  extern gen f__IDNT_e;
  extern identificateur g__IDNT;
  extern gen g__IDNT_e;
  extern identificateur h__IDNT;
  extern gen h__IDNT_e;
  extern identificateur i__IDNT;
  extern gen i__IDNT_e;
  extern identificateur I__IDNT;
  extern gen I__IDNT_e;
  extern identificateur j__IDNT;
  extern gen j__IDNT_e;
  extern identificateur k__IDNT;
  extern gen k__IDNT_e;
  extern identificateur l__IDNT;
  extern gen l__IDNT_e;
  extern identificateur m__IDNT;
  extern gen m__IDNT_e;
  extern identificateur n__IDNT;
  extern gen n__IDNT_e;
  extern identificateur o__IDNT;
  extern gen o__IDNT_e;
  extern identificateur p__IDNT;
  extern gen p__IDNT_e;
  extern identificateur q__IDNT;
  extern gen q__IDNT_e;
  extern identificateur r__IDNT;
  extern gen r__IDNT_e;
  extern identificateur s__IDNT;
  extern gen s__IDNT_e;
  extern identificateur t__IDNT;
  extern gen t__IDNT_e;
  extern identificateur u__IDNT;
  extern gen u__IDNT_e;
  extern identificateur v__IDNT;
  extern gen v__IDNT_e;
  extern identificateur w__IDNT;
  extern gen w__IDNT_e;
  extern identificateur x__IDNT;
  extern gen x__IDNT_e;
  extern gen vx_var;
  extern identificateur y__IDNT;
  extern gen y__IDNT_e;
  extern identificateur z__IDNT;
  extern gen z__IDNT_e;
  extern vecteur list_one_letter__IDNT;
  extern identificateur CST__IDNT;
  extern gen CST__IDNT_e;
  extern identificateur _IDNT_break;
  extern identificateur _IDNT_continue;

  // small utility to remove #...
  int removecomments(const char * ss,char * ss2);

#ifndef NO_NAMESPACE_GIAC
}
#endif // ndef NO_NAMESPACE_GIAC

#endif // _GIAC_IDENTIFICATEUR_H
