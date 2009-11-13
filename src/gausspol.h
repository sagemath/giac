/* -*- mode:C++ ; compile-command: "g++ -I.. -g -c gausspol.cc" -*- */
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
#ifndef _GIAC_GAUSSPOL_H_
#define _GIAC_GAUSSPOL_H_
#include "first.h"
#include "poly.h"
#include "gen.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  extern int primes[100] ;

  class nfactor {
  public:
    gen fact;
    int mult;
    nfactor(const nfactor &n) : fact(n.fact),mult(n.mult) {}
    nfactor(const gen & n, int m) : fact(n),mult(m) {}
  };
  std::vector<nfactor> trivial_n_factor(gen &n);
  vecteur cyclotomic(int n);

  // gen pow(const gen & n,int k);
  typedef std::vector< monomial<gen> > monomial_v;
  typedef tensor<gen> polynome;

  // function on polynomials
  polynome gen2polynome(const gen & e,int dim);
  // check type of coefficients
  int coefftype(const polynome & p,gen & coefft);
  // remove MODulo coefficients
  polynome unmodularize(const polynome & p);
  void modularize(polynome & d,const gen & m);
  // arithmetic
  bool is_one(const polynome & p);
  polynome firstcoeff(const polynome & p);
  polynome operator + (const polynome & th,const polynome & other);
  polynome operator - (const polynome & th,const polynome & other);
  // Fast multiplication using hash maps, might also use an int for reduction
  // but there is no garantee that res is smod-ed modulo reduce
  // Use reduce=0 for non modular
  void mulpoly (const polynome & th, const polynome & other,polynome & res,const gen & reduce);
  polynome operator * (const polynome & th, const polynome & other) ;
  void mulpoly(const polynome & th,const gen & fact,polynome & res);
  polynome operator * (const polynome & th, const gen & fact) ;
  polynome operator * (const gen & fact, const polynome & th) ;
  // a*b+c*d
  gen foisplus(const polynome & a,const polynome & b,const polynome & c,const polynome & d);
  polynome operator - (const polynome & th) ;
  polynome operator / (const polynome & th,const polynome & other);
  polynome operator / (const polynome & th,const gen & fact );
  polynome operator % (const polynome & th,const polynome & other);
  polynome operator % (const polynome & th, const gen & modulo);
  polynome re(const polynome & th);
  polynome im(const polynome & th);
  polynome conj(const polynome & th);
  void polynome2poly1(const polynome & p,int var,vecteur & v);
  vecteur polynome12poly1(const polynome & p);
  int inner_POLYdim(const vecteur & v);
  vecteur polynome2poly1(const polynome & p,int var);
  vecteur polynome2poly1(const polynome & p); // for algebraic ext.
  gen polynome2poly1(const gen & e,int var);
  void poly12polynome(const vecteur & v, int var,polynome & p,int dimension=0);
  polynome poly12polynome(const vecteur & v,int var,int dimension=0);
  polynome poly12polynome(const vecteur & v);
  gen vecteur2polynome(const vecteur & v,int dimension);
  bool divrem1(const polynome & a,const polynome & b,polynome & quo,polynome & r,int exactquo=0) ;
  bool divrem (const polynome & th, const polynome & other, polynome & quo, polynome & rem, bool allowrational = false );
  bool divremmod (const polynome & th,const polynome & other, const gen & modulo,polynome & quo, polynome & rem);
  polynome pow(const polynome & th,int n);
  bool is_positive(const polynome & p);
  polynome pow(const polynome & p,const gen & n);
  polynome powmod(const polynome &p,int n,const gen & modulo);
  polynome gcd(const polynome & p,const polynome & q);
  void lcmdeno(const polynome & p, gen & res);
  polynome ichinrem(const polynome &p,const polynome & q,const gen & pmod,const gen & qmod);
  polynome resultant(const polynome & p,const polynome & q);
  polynome lgcd(const polynome & p);
  gen ppz(polynome & p);
  void lgcdmod(const polynome & p,const gen & modulo,polynome & pgcd);
  polynome gcdmod(const polynome &p,const polynome & q,const gen & modulo);
  // modular gcd via PSR, used when not enough eval points available
  // a and b must be primitive and will be scratched
  void psrgcdmod(polynome & a,polynome & b,const gen & modulo,polynome & prim);
  // Find non zeros coeffs of p, res contains the positions of non-0 coeffs
  int find_nonzero(const polynome & p,index_t & res);
  polynome pzadic(const polynome &p,const gen & n);
  bool gcd_modular_algo(polynome &p,polynome &q, polynome &d,bool compute_cof);
  bool listmax(const polynome &p,gen & n );
  bool gcdheu(const polynome &p,const polynome &q, polynome & p_simp, gen & np_simp, polynome & q_simp, gen & nq_simp, polynome & d, gen & d_content ,bool skip_test=false,bool compute_cofactors=true);
  polynome gcdpsr(const polynome &p,const polynome &q,int gcddeg=0);
  void simplify(polynome & p,polynome & q,polynome & p_gcd);
  polynome simplify(polynome &p,polynome &q);
  void egcd(const polynome &p1, const polynome & p2, polynome & u,polynome & v,polynome & d);
  bool findabcdelta(const polynome & p,polynome & a,polynome &b,polynome & c,polynome & delta);
  bool findde(const polynome & p,polynome & d,polynome &e);
  factorization sqff(const polynome &p );
  void sqff_evident(const polynome & p,factorization & f,bool withsqrt,bool complexmode);
  // factorization over Z[i]
  bool cfactor(const polynome & p, gen & an,factorization & f,bool withsqrt);
  // sqff factorization over a finite field
  factorization squarefree_fp(const polynome & p,unsigned n,unsigned exposant);
  // univariate factorization over a finite field, once sqff
  void sqff_ffield_factor(const factorization & sqff_f,int n,environment * env,factorization & f);

  // factorization over Z[e] where e is an algebraic extension
  bool ext_factor(const polynome &p,const gen & e,gen & an,polynome & p_content,factorization & f,bool complexmode);
  // factorization over Z[coeff_of_p]
  bool factor(const polynome &p,polynome & p_content,factorization & f,bool isprimitive,bool withsqrt,bool complexmode,const gen & divide_by_an=1);
  void unitarize(const polynome &pcur, polynome &unitaryp, polynome & an);
  polynome ununitarize(const polynome & unitaryp, const polynome & an);
  void partfrac(const polynome & num, const polynome & den, const std::vector< facteur< polynome > > & v , std::vector < pf <gen> > & pfde_VECT, polynome & ipnum, polynome & ipden );
  pf<gen> intreduce_pf(const pf<gen> & p_cst, std::vector< pf<gen> > & intde_VECT );
  // Sturm sequences
  vecteur vector_of_polynome2vecteur(const vectpoly & v);
  vecteur sturm_seq(const polynome & p,polynome & cont);

  // prototype of factorization of univariate sqff unitary polynomial
  // provided e.g. by smodular
  void factorunivsqff(const polynome & q,environment * env,vectpoly & v,int & ithprime,int debug,int modfactor_primes);
  // find linear factor only 
  int linearfind(const polynome & q,environment * env,polynome & qrem,vectpoly & v,int & ithprime);
  // prototype of modular 1-d gcd algorithm
  bool gcd_modular_algo1(polynome &p,polynome &q,polynome &d,bool compute_cof);
  polynome smod(const polynome & th, const gen & modulo);
  void smod(const polynome & th, const gen & modulo,polynome & res);
  bool gcdmod_dim1(const polynome &p,const polynome & q,const gen & modulo,polynome & d,polynome & pcof,polynome & qcof,bool compute_cof,bool & real);

  // evaluate p at v by replacing the last variables of p with values of v
  gen peval(const polynome & p,const vecteur & v,const gen &m,bool simplify_at_end=false,std::vector<int_unsigned> * pptr=0);
  int total_degree(const polynome & p);

  // build a multivariate poly
  // with normal coeff from a multivariate poly with multivariate poly coeffs
  polynome unsplitmultivarpoly(const polynome & p,int inner_dim);
  polynome splitmultivarpoly(const polynome & p,int inner_dim);
  polynome split(const polynome & p,int inner_dim);

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC

#endif // _GIAC_GAUSSPOL_H_
