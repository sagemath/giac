// -*- mode:C++ ; compile-command: "g++ -I.. -fPIC -DPIC -g -c pgcd.cpp -o pgcd.lo && ln -sf pgcd.lo pgcd.o && gcc -shared pgcd.lo -lc  -Wl,-soname -Wl,libpgcd.so.0 -o libpgcd.so.0.0.0 && ln -sf libpgcd.so.0.0.0 libpgcd.so.0 && ln -sf libpgcd.so.0.0.0 libpgcd.so" -*-
using namespace std;
#include <stdexcept>
#include <cmath>
#include <cstdlib>
#include <giac/giac.h>
//#include "pgcd.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  gen monpgcd(gen a,gen b){
    gen q,r;
    for (;b!=0;){
      r=irem(a,b,q);
      a=b;
      b=r;
    }
    return a;
  }
  gen _monpgcd(const gen & args){
    if ( (args.type!=_VECT) || (args._VECTptr->size()!=2))
      setsizeerr();
    vecteur &v=*args._VECTptr;
    return monpgcd(v[0],v[1]);
  }
  const string _monpgcd_s("monpgcd");
  unary_function_unary __monpgcd(&_monpgcd,_monpgcd_s);
  unary_function_ptr at_monpgcd (&__monpgcd,0,true);
  

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC
