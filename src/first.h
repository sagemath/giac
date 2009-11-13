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
#ifndef _GIAC_FIRST_H_
#define _GIAC_FIRST_H_
#include "config.h"

#ifdef USE_GMP_REPLACEMENTS
#include "gmp_replacements.h"
#else
#include <gmp.h>
#endif

#include <assert.h>
#ifdef HAVE_LIBGC
#define GC_DEBUG
#include <gc_cpp.h>
void* operator new( size_t size ); 
  
void operator delete( void* obj ); 
  
void* operator new[]( size_t size );
  
void operator delete[]( void* obj );

void * RS_gmpalloc(size_t a);

void * RS_gmprealloc(void * old_p,size_t old_size,size_t new_size);

void RS_gmpfree(void * old_p,size_t old_size);
class _init_gmp_memory {
  int done;
 public:
  _init_gmp_memory(): done(1){
    // bind GMP allocation to GC
    mp_set_memory_functions(RS_gmpalloc,RS_gmprealloc,RS_gmpfree);
  }
};
#endif // HAVE_LIBGC

#ifdef __VISUALC__ // Visual C++?
typedef long pid_t;
typedef __int64 longlong ;
typedef unsigned __int64 ulonglong ;
#else
typedef long long longlong;
typedef unsigned long long ulonglong;
#endif



#endif // _GIAC_FIRST_H_
