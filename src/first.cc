/* -*- compile-command: "g++ -g -c -I.. first.cc" -*-
 *  Copyright (C) 2000,7 B. Parisse, Institut Fourier, 38402 St Martin d'Heres
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
#include "first.h"

#ifdef HAVE_LIBGC
void* operator new( size_t size ) {
    return GC_MALLOC_UNCOLLECTABLE( size );}
  
void operator delete( void* obj ) {
    GC_FREE( obj );}
  
void* operator new[]( size_t size ) {
    return GC_MALLOC_UNCOLLECTABLE( size );}
  
void operator delete[]( void* obj ) {
    GC_FREE( obj );}

void * RS_gmpalloc(size_t a)
{
  return(GC_malloc_atomic(a));
  //return(GC_MALLOC_UNCOLLECTABLE(a));
}

void * RS_gmprealloc(void * old_p,size_t old_size,size_t new_size)
{
  void * tmp=GC_realloc(old_p,new_size);
  return(tmp);
}

void RS_gmpfree(void * old_p,size_t old_size)
{
  // GC_free(old_p);
}

_init_gmp_memory __giac_init_gnp_memory;

#endif

