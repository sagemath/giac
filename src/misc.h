/* -*- mode:C++ ; compile-command: "g++ -I.. -g -c misc.cc" -*-
 *
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

#ifndef _GIAC_MISC_H_
#define _GIAC_MISC_H_
#include "first.h"
#include "global.h"
#include "gen.h"
#include "unary.h"
#include "symbolic.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC
  gen _normalize(const gen & a,GIAC_CONTEXT);
  extern unary_function_ptr at_normalize;
  gen _degree(const gen & args,GIAC_CONTEXT);
  gen _tcoeff(const gen & args,GIAC_CONTEXT);
  gen _lcoeff(const gen & args,GIAC_CONTEXT);
  gen _primpart(const gen & g,GIAC_CONTEXT);
  void aplatir(const matrice & m,vecteur & v);

  gen _dfc2f(const gen & g,GIAC_CONTEXT);
  vecteur gen2continued_fraction(const gen & g,int n,GIAC_CONTEXT);
  gen _float2rational(const gen & g,GIAC_CONTEXT);
  gen float2rational(double d_orig,double eps,GIAC_CONTEXT);

  gen _preval(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_preval;

  gen _dotprod(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_dotprod;

  gen _mean(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_mean;

  gen _median(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_median;

  gen _stddev(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_stddev;

  gen _variance(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_variance;

  vecteur divided_differences(const vecteur & x,const vecteur & y);
  gen _lagrange(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_lagrange;

  gen _reorder(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_reorder;

  gen _adjoint_matrix(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_adjoint_matrix;

  gen _equal2diff(const gen & args);
  extern unary_function_ptr at_equal2diff;

  gen _rank(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_rank;

  gen _diag(const gen & args);
  extern unary_function_ptr at_diag;

  gen _sec(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_sec;

  gen _csc(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_csc;

  gen _cot(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_cot;

  gen _asec(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_asec;

  gen _acsc(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_acsc;

  gen _acot(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_acot;

  gen _ibpu(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_ibpu;

  gen _changebase(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_changebase;

  gen _epsilon2zero(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_epsilon2zero;

  gen _suppress(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_suppress;

  gen _froot(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_froot;

  gen _fcoeff(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_fcoeff;

  gen _truncate(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_truncate;

  gen _divpc(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_divpc;

  gen _ptayl(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_ptayl;

  gen _float2rational(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_float2rational;

  gen _gramschmidt(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_gramschmidt;

  gen _pmin(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_pmin;

  gen _potential(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_potential;

  gen _vpotential(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_vpotential;

  gen _symb2poly1(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_symb2poly1;

  gen _poly12symb(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_poly12symb;

  gen _exp2trig(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_exp2trig;

  gen _nrows(const gen & args);
  extern unary_function_ptr at_nrows;

  gen _ncols(const gen & args);
  extern unary_function_ptr at_ncols;

  gen _l2norm(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_l2norm;

  gen _input(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_input;
  class unary_function_unary;
  extern unary_function_eval __input;

  gen _histogram(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_histogram;
  matrice effectifs(const vecteur & v,double class_minimum,double class_size,GIAC_CONTEXT);

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // NO_NAMESPACE_GIAC

#endif // _GIAC_MISC_H
