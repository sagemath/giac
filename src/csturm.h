// -*- mode:C++ ; compile-command: "g++ -I.. -g -c csturm.cc" -*-
/*
 *  Copyright (C) 2007 B. Parisse, Institut Fourier, 38402 St Martin d'Heres
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
#ifndef _GIAC_CSTURM_H
#define _GIAC_CSTURM_H
#include "first.h"
#include <complex>
#include <iostream>
#include "modpoly.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  // P(x)->P(a*x)
  void change_scale(modpoly & p,const gen & l);

  // compute Sturm sequence of r0 and r1,
  // returns gcd (without content)
  // and compute list of quotients, coeffP, coeffR
  // such that coeffR*r_(k+2) = Q_k*r_(k+1) - coeffP_k*r_k
  gen csturm_seq(modpoly & r0,modpoly & r1,vecteur & listquo,vecteur & coeffP, vecteur & coeffR,GIAC_CONTEXT);

  // number of complex roots of p in a..b, returns -1 if a factor
  // has been found
  int csturm_square(const gen & p,const gen & a,const gen & b,gen& pgcd,GIAC_CONTEXT);

  // Find complex roots of P in a0,b0 -> a1,b1
  // Returns false if a factor of P has been found (->in pgcd)
  bool complex_roots(const modpoly & P,const gen & a0,const gen & b0,const gen & a1,const gen & b1,vecteur & roots,gen & pgcd,double eps);



#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC

#endif // _GIAC_CSTURM_H
