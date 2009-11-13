// -*- mode:C++ ; compile-command: "g++ -I.. -g -c moyal.cc" -*-
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
#ifndef _GIAC_MOYAL_H
#define _GIAC_MOYAL_H
#include "first.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC
  class gen;
  gen moyal(const gen & a,const gen & b,const gen & vars,const gen & order);
  gen _moyal(const gen & args);
  extern unary_function_ptr at_moyal ;

  gen Airy_Ai(const gen & a,const gen & b,const gen & vars,const gen & order);
  gen _Airy_Ai(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_Airy_Ai ;

  gen Airy_Bi(const gen & a,const gen & b,const gen & vars,const gen & order);
  gen _Airy_Bi(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_Airy_Bi ;

  gen randNorm();
  gen _randNorm(const gen & args);
  extern unary_function_ptr at_randNorm ;

  gen _UTPN(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_UTPN ;

  gen _UTPC(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_UTPC ;

  gen _UTPT(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_UTPT ;

  gen _UTPF(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_UTPF ;

  gen _binomial(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_binomial ;

  gen _binomial_cdf(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_binomial_cdf ;

  gen _binomial_icdf(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_binomial_icdf ;

  gen _poisson(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_poisson ;

  gen _poisson_cdf(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_poisson_cdf ;

  gen _poisson_icdf(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_poisson_icdf ;

  gen _normal_cdf(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_normal_cdf ;

  gen _normal_icdf(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_normal_icdf ;

  gen _student(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_student ;

  gen _student_cdf(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_student_cdf ;

  gen _student_icdf(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_student_icdf ;

  gen _chisquare(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_chisquare ;

  gen _chisquare_cdf(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_chisquare_cdf ;

  gen _chisquare_icdf(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_chisquare_icdf ;

  gen _snedecor(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_snedecor ;

  gen _snedecor_cdf(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_snedecor_cdf ;

  gen _snedecor_icdf(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_snedecor_icdf ;

  gen _Beta(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_Beta ;

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC


#endif // _GIAC_MOYAL_H
