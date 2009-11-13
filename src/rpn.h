// -*- mode:C++ ; compile-command: "g++ -I.. -g -c rpn.cc" -*-
/*
 *  Copyright (C) 2001 B. Parisse, Institut Fourier, 38402 St Martin d'Heres
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
#ifndef _GIAC_RPN_H
#define _GIAC_RPN_H
#include "first.h"

#include "gen.h"
#include "vecteur.h"
#include <string>
#include <ctype.h>

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  std::string printasconstant(const gen & feuille,const std::string & sommetstr,GIAC_CONTEXT);
  extern bool rpn_mode;
  std::string enmajuscule(const std::string & s);
  void roll(int i,vecteur & v);
  void ROLL(int i,GIAC_CONTEXT);
  gen symb_ROLL(const gen & args);
  extern const std::string _ROLL_s;
  gen _ROLL(const gen & args);
  extern unary_function_ptr at_ROLL;

  void rolld(int i,vecteur & v);
  void ROLLD(int i,GIAC_CONTEXT);
  gen symb_ROLLD(const gen & args);
  extern const std::string _ROLLD_s;
  gen _ROLLD(const gen & args);
  extern unary_function_ptr at_ROLLD;

  void stack_swap(vecteur & v);
  void SWAP(GIAC_CONTEXT);
  gen symb_SWAP(const gen & args);
  extern const std::string _SWAP_s;
  gen _SWAP(const gen & args);
  extern unary_function_ptr at_SWAP;

  void dup(vecteur & v);
  gen symb_DUP(const gen & args);
  extern const std::string _DUP_s;
  gen _DUP(const gen & args);
  extern unary_function_ptr at_DUP;

  void over(vecteur & v);
  gen symb_OVER(const gen & args);
  extern const std::string _OVER_s;
  gen _OVER(const gen & args);
  extern unary_function_ptr at_OVER;

  void pick(int i,vecteur & v);
  gen symb_PICK(const gen & args);
  extern const std::string _PICK_s;
  gen _PICK(const gen & args);
  extern unary_function_ptr at_PICK;

  void drop(vecteur & v);
  gen symb_DROP(const gen & args);
  extern const std::string _DROP_s;
  gen _DROP(const gen & args);
  extern unary_function_ptr at_DROP;

  gen symb_NOP(const gen & args);
  extern const std::string _NOP_s;
  gen _NOP(const gen & args);
  extern unary_function_ptr at_NOP;

  gen symb_IFTE(const gen & args);
  extern const std::string _IFTE_s;
  gen _IFTE(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_IFTE;

  gen symb_RPN_LOCAL(const gen & a,const gen & b);
  extern const std::string _RPN_LOCAL_s;
  gen _RPN_LOCAL(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_RPN_LOCAL;

  gen symb_RPN_FOR(const gen & a,const gen & b);
  extern const std::string _RPN_FOR_s;
  gen _RPN_FOR(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_RPN_FOR;

  gen symb_RPN_WHILE(const gen & a,const gen & b);
  extern const std::string _RPN_WHILE_s;
  gen _RPN_WHILE(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_RPN_WHILE;

  gen symb_RPN_UNTIL(const gen & a,const gen & b);
  extern const std::string _RPN_UNTIL_s;
  gen _RPN_UNTIL(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_RPN_UNTIL;

  gen symb_RPN_CASE(const gen & a);
  extern const std::string _RPN_CASE_s;
  gen _RPN_CASE(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_RPN_CASE;

  gen symb_RCL(const gen & a,const gen & b);
  extern const std::string _RCL_s;
  gen _RCL(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_RCL;

  gen symb_VARS(const gen & a,const gen & b);
  extern const std::string _VARS_s;
  gen _VARS(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_VARS;

  gen symb_purge(const gen & a,const gen & b);
  extern const std::string _purge_s;
  gen _purge(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_purge;

  gen symb_rpn(const gen & args);
  extern const std::string _rpn_s;
  gen _rpn(const gen & args);
  extern unary_function_ptr at_rpn;

  gen symb_alg(const gen & args);
  extern const std::string _alg_s;
  gen _alg(const gen & args);
  extern unary_function_ptr at_alg;

  gen symb_rpn_prog(const gen & args);
  extern const std::string _rpn_prog_s;
  gen _rpn_prog(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_rpn_prog;
  vecteur rpn_eval(const vecteur & prog,vecteur & pile,GIAC_CONTEXT);
  vecteur rpn_eval(const gen & prog,vecteur & pile,GIAC_CONTEXT);

  gen symb_division(const gen & a,const gen & b);
  gen symb_division(const gen & args);
  extern const std::string _division_s;
  gen _division(const gen & args);
  extern unary_function_ptr at_division;

  gen symb_binary_minus(const gen & a,const gen & b);
  gen symb_binary_minus(const gen & args);
  extern const std::string _binary_minus_s;
  gen _binary_minus(const gen & args);
  extern unary_function_ptr at_binary_minus;
  vecteur tab2vecteur(gen tab[]);

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC

#endif // _GIAC_RPN_H
