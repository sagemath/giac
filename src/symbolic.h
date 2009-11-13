// -*- mode:C++ ; compile-command: "g++ -I.. -g -c symbolic.cc" -*-
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
#ifndef _GIAC_SYMBOLIC_H
#define _GIAC_SYMBOLIC_H
#include "first.h"
#include <iostream>
#include <string>
#include <vector>
#include "unary.h"
#include "gen.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  extern unary_function_ptr archive_function_tab[];
  struct symbolic {
    unary_function_ptr sommet; 
    gen feuille;
    symbolic(const unary_function_ptr & o,const gen & e): sommet(o),feuille(e){};
    symbolic(const unary_function_ptr & o,const gen & e1,const gen &e2): sommet(o), feuille(makevecteur(e1,e2)) {};
    symbolic(const unary_function_ptr & o,const gen & e1,const gen &e2,const gen & e3): sommet(o), feuille(makevecteur(e1,e2,e3)) {};
    symbolic(const unary_function_ptr & o,const gen & e1,const gen &e2,const gen & e3,const gen & e4): sommet(o), feuille(makevecteur(e1,e2,e3,e4)) {};
    symbolic(const symbolic & mys) : sommet(mys.sommet),feuille(mys.feuille) {};
    symbolic(const symbolic & mys,const gen & e);
    symbolic(const gen & a,const unary_function_ptr & o,const gen & b);
    std::string print(GIAC_CONTEXT) const;
    void dbgprint() const{ std::cout << this->print(context0) << std::endl; }
    gen eval(int level,const context * context_ptr) const;
    gen evalf(int level,const context * context_ptr) const;
    int size() const;
  };

  struct ref_symbolic {
    int ref_count;
    symbolic s;
  };
  
  std::ostream & operator << (std::ostream & os,const symbolic & s);
  
  int equalposcomp(unary_function_ptr tab[],const unary_function_ptr & f);

  // find the "size" of g but limited by max
  unsigned taille(const gen & g,unsigned max);
  extern bool print_rewrite_prod_inv;
  // try to rewrite arg the argument of a product as a fraction n/d
  bool rewrite_prod_inv(const gen & arg,gen & n,gen & d);

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC

#endif // _GIAC_SYMBOLIC_H
