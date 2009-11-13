// -*- mode:C++ ; compile-command: "g++ -I.. -g -c maple.cc" -*-
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
#ifndef _GIAC_MAPLE_H
#define _GIAC_MAPLE_H
#include "first.h"
#ifdef HAVE_LIBPNG
#include <png.h>
#endif

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC
  class gen;

  gen _about(const gen & g,GIAC_CONTEXT);
  gen _zip(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_zip ;

  // product(P,n=a..b) where the first variable in v is n
  gen product(const polynome & P,const vecteur & v,const gen & n,const gen & a,const gen & b,GIAC_CONTEXT);

  gen _accumulate_head_tail(const gen & args);
  extern unary_function_ptr at_accumulate_head_tail ;
  gen fft(const gen & g_orig,int direct,GIAC_CONTEXT);
  gen _evalc(const gen & g,GIAC_CONTEXT);
  extern unary_function_ptr at_gcdex ;
  extern unary_function_ptr at_seqsolve ;
  extern unary_function_ptr at_array ;

#ifdef HAVE_LIBPNG
  int write_png(const char *file_name, png_bytep *rows, int w, int h, int colortype, int bitdepth);
#endif

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC


#endif // _GIAC_MAPLE_H
