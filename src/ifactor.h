// -*- mode:C++ -*-
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

#ifndef _GIAC_IFACTOR_H_
#define _GIAC_IFACTOR_H_
#include "first.h"
#include "global.h"
#include "gen.h"
#include "unary.h"
#include "symbolic.h"
#include "identificateur.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC
  symbolic symb_ifactor(const gen & args);
  vecteur ifactors(const gen & n0);
  gen ifactors(const gen & args,int maplemode);
  extern unary_function_ptr at_ifactors;
  extern unary_function_ptr at_maple_ifactors;

  vecteur factors(const gen & g,const gen & x,GIAC_CONTEXT);
  extern unary_function_ptr at_factors;

  gen _divis(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_divis ;

  gen _idivis(const gen & args);
  extern unary_function_ptr at_idivis ;

  vecteur pfacprem(gen & n,bool addlast=true);

  gen _ifactor(const gen & args);
  extern unary_function_ptr at_ifactor ;

  gen _euler(const gen & args);
  symbolic symb_euler(const gen & args);
  extern unary_function_ptr at_euler;

  gen _pa2b2(const gen & args);
  extern unary_function_ptr at_pa2b2 ;
 
  gen _propfrac(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_propfrac ;

  gen _iabcuv(const gen & args);
  extern unary_function_ptr at_iabcuv ;

  gen _abcuv(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_abcuv ;

  gen _simp2(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_simp2 ;

  gen _fxnd(const gen & args);
  extern unary_function_ptr at_fxnd ;

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // NO_NAMESPACE_GIAC

#endif // _GIAC_IFACTOR_H
