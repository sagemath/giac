// -*- mode:C++ ; compile-command: "g++ -I.. -fPIC -DPIC -g -c renee.cc -o renee.lo && ln -sf renee.lo renee.o && gcc -shared renee.lo -lc  -Wl,-soname -Wl,librenee.so.0 -o librenee.so.0.0.0 && ln -sf librenee.so.0.0.0 librenee.so.0 && ln -sf librenee.so.0.0.0 librenee.so" -*-
/*
 *  Copyright (C) 2002 Renee De Graeve, Inst. Fourier, 38402 St Martin d'Heres
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
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
using namespace std;
#include <stdexcept>
#include <cmath>
#include <cstdlib>
#include <giac/giac.h>
//#include "renee.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

#include <stdexcept>
#include <cmath>
#include <cstdlib>
#include <giac/giac.h>
  gen _split(const gen & args){
    //renvoie [ax,ay] si les arg sont g1=ax*an (sans denominateur) et g2=[x,y]
    if ( (args.type!=_VECT)  || (args._VECTptr->size()!=2) )
      settypeerr();
    vecteur v(plotpreprocess(args));
    gen g1=v.front(),g2=v.back();
    if (g2.type!=_VECT)
      settypeerr();
    vecteur v2(*g2._VECTptr);
    int n=v2.size();
    if (n!=2) setsizeerr();
    vecteur fa;
    fa=factors(g1);
    int l=fa.size();
    gen ax=1;
    gen ay=1;
    for (int k=0;k<l;k=k+2){
      gen f=fa[k];
      if (derive(f,v2[0])==0) {
        ay=ay*pow(f,fa[k+1]);}
      else {if (derive(f,v2[1])==0){
        ax=ax*pow(f,fa[k+1]);}
      else {vecteur res(1);return (res);}
      }
    }
    vecteur res(2);
    res[0]=ax;
    res[1]=ay;
    return res;
  }
 
  const string _split_s("split");
  unary_function_unary __split(&_split,_split_s);
  unary_function_ptr at_split (&__split,_QUOTE_ARGUMENTS,true);
 
 gen _sum_riemann(const gen & args){
   //sum de k=1 a n de g1 avec g2=[n,k])                                            if ( (args.type!=_VECT)  || (args._VECTptr->size()!=2) )
      settypeerr();
    gen g1=args[0],g2=args[1];
    if (g2.type!=_VECT)
      settypeerr();
    vecteur v2(*g2._VECTptr);
    int sv2=v2.size();
    if (sv2!=2) setsizeerr();
    identificateur x("_x");
    gen a=normal(v2[0]*subst(g1,v2[1],x*v2[0]));
    gen nda=_fxnd(a);
    vecteur nada(*nda._VECTptr);
    gen na=nada[0];
    gen da=nada[1];
    //gen na;
    //gen da;
    //fxnd(a,na,da);
    vecteur var(3);
    var[0]=na;
    var[1]=x;
    var[2]=v2[0];
    //var=[x,n]
    gen b0=_split(var);
    vecteur vb0(*b0._VECTptr);
    if ((vb0.size()==1)&& (vb0[0]==0)) {settypeerr("ce n'est pas une
somme de riemann");}
    var[0]=da;
    gen b1=_split(var);
    vecteur vb1(*b1._VECTptr);
    if ((vb1.size()==1)&& (vb1[0]==0)) {settypeerr ("ce n'est pas une
somme de riemann");}
    gen ax=b0[0]/b1[0];
    gen an=b0[1]/b1[1];
    return
(_integrate(makevecteur(ax,x,0,1))*_limit(makevecteur(an,v2[0],plus_inf)));
  }
  const string _sum_riemann_s("sum_riemann");
  unary_function_unary __sum_riemann(&_sum_riemann,_sum_riemann_s);
  unary_function_ptr at_sum_riemann(&__sum_riemann,_QUOTE_ARGUMENTS,true);
#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC
