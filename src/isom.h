// -*- mode:C++ ; compile-command: "g++ -I.. -g -c isom.cc " -*- 
/*
 *  Copyright (C) 2001 R. De Graeve, Institut Fourier, 38402 St Martin d'Heres
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

#ifndef _GIAC_ISOM_H
#define _GIAC_ISOM_H
#include "first.h"
#include "gen.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  vecteur isom(const vecteur & M,GIAC_CONTEXT);
  extern const std::string _isom_s;
  gen symb_isom(const gen & args);
  extern unary_function_ptr at_isom ;

  vecteur mkisom(const gen & n,int b,GIAC_CONTEXT);
  gen symb_mkisom(const gen & args);
  gen symb_mkisom(const gen & q,const gen & x);
  extern const std::string _mkisom_s;
  extern unary_function_ptr at_mkisom;

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC

#endif // _GIAC_ISOM_H
