// -*- mode:C++ ; compile-command: "g++ -I.. -g -c permu.cc" -*-
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

#ifndef _GIAC_PERMU_H_
#define _GIAC_PERMU_H_
#include "first.h"
#include "global.h"
#include "gen.h"
#include "unary.h"
#include "symbolic.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC
  vecteur vector_int_2_vecteur(const std::vector<int> & v,GIAC_CONTEXT);
  std::vector<int> vecteur_2_vector_int(const vecteur & v,GIAC_CONTEXT);
  std::vector< std::vector<int> > vecteur_2_vectvector_int(const vecteur & v,GIAC_CONTEXT);
  vecteur vectvector_int_2_vecteur(const std::vector< std::vector<int> > & v,GIAC_CONTEXT);
  gen square_hadamard_bound(const matrice & m);

  std::vector<int> randperm(const int & n);
  gen _randperm(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_randperm;

  bool is_permu(const vecteur &p,std::vector<int> & p1,GIAC_CONTEXT);
  gen _is_permu(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_is_permu;

  bool is_cycle(const vecteur & c,std::vector<int> & c1,GIAC_CONTEXT);
  gen _is_cycle(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_is_cycle;

  std::vector< std::vector<int> > permu2cycles(const std::vector<int> & p) ;
  gen _permu2cycles(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_permu2cycles;

  std::vector<int> cycle2permu(const std::vector< std::vector<int> > & c);
  gen _cycle2permu(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_cycle2permu;

  int signature(const std::vector<int> & p);
  gen _signature(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_signature;

  std::vector<int> p1op2(const std::vector<int> & p1,const std::vector<int> & p2);
  gen _p1op2(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_p1op2;
  
  std::vector<int> c1oc2(const std::vector<int> & c1,const std::vector<int> & c2);
  gen _c1oc2(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_c1oc2;

  std::vector<int> c1op2(const std::vector<int> & c1, const std::vector<int> & p2);
  gen _c1op2(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_c1op2;

  std::vector<int> p1oc2(const std::vector<int> & p1, const std::vector<int> & c2);
  gen _p1oc2(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_p1oc2;

  // arithmetic mean column by column
  vecteur mean(const matrice & m,bool column=true);
  vecteur stddev(const matrice & m,bool column=true,int variance=1);
  matrice ascsort(const matrice & m,bool column=true);

  gen _divergence(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_divergence;

  gen _curl(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_curl;

  gen _permu2mat(const gen & args,GIAC_CONTEXT); // permutation vector -> matrix

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // NO_NAMESPACE_GIAC

#endif // _GIAC_PERMU_H
