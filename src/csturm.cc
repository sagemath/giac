// -*- mode:C++ ; compile-command: "g++ -I.. -I../include -g -c -Wall csturm.cc" -*-
#include "first.h"
/*
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
using namespace std;
#include <cmath>
#include <stdexcept>
#include <map>
#include "gen.h"
#include "csturm.h"
#include "vecteur.h"
#include "modpoly.h"
#include "unary.h"
#include "symbolic.h"
#include "usual.h"
#include "sym2poly.h"
#include "solve.h"
#include "prog.h"
#include "subst.h"
#include "permu.h"
#include "series.h"
#include "alg_ext.h"
#include "ti89.h"
#include "plot.h"
#include "modfactor.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  // compute Sturm sequence of r0 and r1,
  // returns gcd (without content)
  // and compute list of quotients, coeffP, coeffR
  // such that coeffR*r_(k+2) = Q_k*r_(k+1) - coeffP_k*r_k
  gen csturm_seq(modpoly & r0,modpoly & r1,vecteur & listquo,vecteur & coeffP, vecteur & coeffR,GIAC_CONTEXT){
    listquo.clear();
    coeffP.clear();
    coeffR.clear();
    if (r0.empty())
      return r1;
    if (r1.empty())
      return r0;
    gen tmp;
    lcmdeno(r0,tmp,contextptr);
    if (ck_is_positive(-tmp,contextptr))
      r0=-r0;
    r0=r0/abs(lgcd(r0),contextptr);
    lcmdeno(r1,tmp,contextptr);
    if (ck_is_positive(-tmp,contextptr))
      r1=-r1;
    r1=r1/abs(lgcd(r0),contextptr);
    // set auxiliary constants g and h to 1
    gen g(1),h(1);
    modpoly a(r0),b(r1),quo,r;
    gen b0(1);
    for (int loop_counter=0;;++loop_counter){
      int m=a.size()-1;
      int n=b.size()-1;
      int ddeg=m-n; // should be 1 generically
      if (!n) { // if b is constant, gcd=1
	return 1;
      }
      if (ddeg%2) // ddeg should be even if b0 is a _POLY1
	b0=b.front();
      else
	b0=abs(b.front(),contextptr); 
      coeffP.push_back(pow(b0,ddeg+1));
      DivRem(coeffP.back()*a,b,0,quo,r);
      listquo.push_back(quo);
      coeffR.push_back(g*pow(h,ddeg));
      if (r.empty()){
	return b/abs(lgcd(b),contextptr);
      }
      // remainder is non 0, loop continue: a <- b
      a=b;
      // now divides r by g*h^(m-n) and change sign, result is the new b
      b= -r/coeffR.back();
      g=b0;
      h=pow(b0,ddeg)/pow(h,ddeg-1);
    } // end while loop
  }

  int csturm_vertex_ab(const modpoly & r0,const modpoly & r1,const vecteur & listquo,const vecteur & coeffP, const vecteur & coeffR,const gen & a,int start,GIAC_CONTEXT){
    int n=listquo.size(),j,k,res=0;
    vecteur R(n+2);
    R[0]=horner(r0,a);
    R[1]=horner(r1,a);
    for (j=0;j<n;j++)
      R[j+2]=(-coeffP[j]*R[j]+R[j+1]*horner(listquo[j],a))/coeffR[j];
    // signes
    for (j=start;j<n+2;j++){
      if (R[j]!=0) break;
    }
    for (k=j+1;k<n+2;k++){
      if (is_zero(R[k])) continue;
      if (is_positive(-R[j]*R[k],contextptr)){
	res++;
	j=k;
      }
    }
    return res;
  }

  // compute vertex index at a (==0 unless s(a)==0) 
  int csturm_vertex_a(const modpoly & s,const modpoly & r,const gen & a,int direction,GIAC_CONTEXT){
    int j;
    modpoly s1,s2;
    gen sa=horner(s,a,0,s1);
    if (!is_zero(sa)) return 0;
    for (j=1;;j++){
      sa=horner(s1,a,0,s2);
      if (!is_zero(sa))
	break;
      s1=s2;
    }
    if (direction==1) j=0;
    gen tmp=sign(sa,contextptr)*sign(horner(r,a),contextptr);
    return tmp.val*((j%2)?-1:1);
  }

  void change_scale(modpoly & p,const gen & l){
    int n=p.size();
    gen lton(l);
    for (int i=n-2;i>=0;--i){
      p[i] = p[i] * lton;
      lton = lton * l;
    }
  }

  // p(x)->p(a*x+b)
  modpoly linear_changevar(const modpoly & p,const gen & a,const gen & b){
    modpoly res(taylor(p,b));
    change_scale(res,a);
    return res;
  }

  // p(a*x+b)->p(x)
  // t=a*x+b -> pgcd(t)=g((t-b)/a)
  modpoly inv_linear_changevar(const modpoly & p,const gen & a,const gen & b){
    gen A=inv(a,context0); 
    gen B=-b/a;
    modpoly res(taylor(p,B));
    change_scale(res,A);
    return res;
  }

  // Find roots of R, S=R' at precision eps, returns number of roots
  // if eps==0 does not compute intervals for roots
  int csturm_realroots(const modpoly & S,const modpoly & R,const vecteur & listquo,const vecteur & coeffP, const vecteur & coeffR,const gen & a,const gen & b,const gen & t0, const gen & t1,vecteur & realroots,double eps,GIAC_CONTEXT){
    int n1=csturm_vertex_ab(S,R,listquo,coeffP,coeffR,t0,1,contextptr);
    int n2=csturm_vertex_ab(S,R,listquo,coeffP,coeffR,t1,1,contextptr);
    int n=(n2-n1);
    if (!eps || !n)
      return n;
    if (is_strictly_greater(eps,(t1-t0)*abs(b,contextptr),contextptr)){
      realroots.push_back(makevecteur(makevecteur(a+t0*b,a+t1*b),n));
      return n;
    }
    gen t01=(t0+t1)/2;
    csturm_realroots(S,R,listquo,coeffP,coeffR,a,b,t0,t01,realroots,eps,contextptr);
    csturm_realroots(S,R,listquo,coeffP,coeffR,a,b,t01,t1,realroots,eps,contextptr);
    return n;
  }

  // Find complex sturm sequence for P(a+(b-a)*x)
  // If P is "pseudo"-real on [a,b] and eps>0 put roots in [a,b]
  // at precision eps inside realroots
  // returns a,b,R,S,g,listquo,coeffP,coeffR
  // and typeseq=0 (complex Sturm) or 1 (limit)
  // If b-a is real and horiz_sturm is not empty, it tries to replace
  // the variable by im(a)*i in horiz_sturm and if no quotient in horiz_sturm
  // has a leading 0 coefficient, 
  // it returns im(a)*i,im(a)*i+1,R,S,g,listquo,coeffP,coeffR
  // If b-a is pure imaginary and vert_sturm is not empty, it tries to replace
  // the variable by re(a) and returns re(a),re(a)+i,R,S,g,listquo,coeffP,coeffR
  vecteur csturm_segment_seq(const modpoly & P,const gen & a,const gen & b,vecteur & realroots,double eps,vecteur & horiz_sturm,vecteur & vert_sturm,GIAC_CONTEXT){
    // try with horiz_sturm and vert_sturm
    gen ab(b-a);
    // /*
    if (is_zero(re(ab,contextptr))){ // b-a is pure imaginary
      if (vert_sturm.empty()){
	gen A=gen(makevecteur(1,0),_POLY1__VECT);
	vert_sturm.push_back(undef);
	vecteur tmp;
	vert_sturm=csturm_segment_seq(P,A,A+cst_i,tmp,eps,horiz_sturm,vert_sturm,contextptr);
      }
      if (vert_sturm.size()==9){
	vecteur res(vert_sturm);
	gen A=re(a,contextptr);
	res[0]=A; // re(a)
	res[1]=A+cst_i; // re(a)+i
	res[2]=apply1st(res[2],A,horner); // R 
	res[3]=apply1st(res[3],A,horner); // S
	res[4]=horner(res[4],A); // g
	vecteur tmp(*res[5]._VECTptr);
	int tmps=tmp.size();
	for (int j=0;j<tmps;++j)
	  tmp[j]=apply1st(tmp[j],A,horner); 
	res[5]=tmp; // listquo
	res[6]=apply1st(res[6],A,horner); // coeffP
	res[7]=apply1st(res[7],A,horner); // coeffR
	if (res[6].type==_VECT && !equalposcomp(*res[6]._VECTptr,0))
	  return res;
	else
	  cerr << "list of quotients is not regular" << endl;
      }
    }
    // */
    modpoly Q(taylor(P,a));
    change_scale(Q,b-a);
    // now x is in 0..1
    gen gtmp=apply(Q,re,contextptr);
    if (gtmp.type!=_VECT)
      setsizeerr();
    modpoly R=trim(*gtmp._VECTptr,0);
    gtmp=apply(Q,im,contextptr);
    if (gtmp.type!=_VECT)
      setsizeerr();
    modpoly S=trim(*gtmp._VECTptr,0);
    modpoly listquo,coeffP,coeffR;
    gen g=csturm_seq(S,R,listquo,coeffP,coeffR,contextptr);
    int typeseq=-1;
    if (debug_infolevel)
      *logptr(contextptr)  << "segment " << a << ".." << b << ", im/re:" << S << "|" << R << ", gcd:" << g << endl;
    if (g.type==_VECT && g._VECTptr->size()==P.size()){
      // if g==P (up to a constant), use real Sturm sequences
      if (debug_infolevel)
	*logptr(contextptr)  << "Real-kind roots: " << g << endl;
      R=*g._VECTptr;
      S=derivative(R);
      g=csturm_seq(S,R,listquo,coeffP,coeffR,contextptr);
      typeseq=csturm_realroots(S,R,listquo,coeffP,coeffR,a,b-a,0,1,realroots,eps,contextptr);
    }
    if (g.type==_VECT)
      g=inv_linear_changevar(*g._VECTptr,b-a,a);
    vecteur res= makevecteur(a,b,R,S,g,listquo,coeffP,coeffR,typeseq);
    return res;
  }

  // index for segment a,b (2* number of roots when summed over a closed
  // polygon). Note that if S=ImP along the segment is 0 we remove
  // the roots on [a,b] using real Sturm sequences
  // If S=0 at a or b, this is simply ignored
  // Indeed the computed index is then the same as if S was of the
  // sign of R, and since R!=0 if S is 0 this is a property of the vertex
  // not of the segment (note that contrary to counting real roots
  // on an interval, S can vanish as many times as long as R keeps
  // the same sign, without modifying the algebraic number of Im=0
  // cuts if S has the same sign on both end)
  int csturm_segment(const vecteur & seq,const gen & a,const gen & b,GIAC_CONTEXT){
    gen t0,t1;
    if (seq.size()!=9)
      setsizeerr();
    gen aseq=seq[0];
    gen bseq=seq[1];
    gen directeur=(b-a)/(bseq-aseq);
    t0=(a-aseq)/(bseq-aseq);
    if ( !is_zero(im(directeur,contextptr)) || !is_zero(im(t0,contextptr)) )
      setsizeerr();
    t0=re(t0,contextptr); // t0=normal(t0);
    t1=re(t0+directeur,contextptr); // t1=normal(t0+directeur);
    int signe=1;
    if (is_strictly_greater(t0,t1,contextptr)){
      signe=-1;
      std::swap<gen>(t0,t1);
    } 
    const modpoly & R=*seq[2]._VECTptr;
    const modpoly & S=*seq[3]._VECTptr;
    gen g=seq[4];
    const modpoly & listquo=*seq[5]._VECTptr;
    const modpoly & coeffP=*seq[6]._VECTptr;
    const modpoly & coeffR=*seq[7]._VECTptr;
    int debut=(seq[8].val==-1)?0:1;
    int tmp = csturm_vertex_ab(S,R,listquo,coeffP,coeffR,t0,debut,contextptr);
    int res = tmp;
    tmp = csturm_vertex_ab(S,R,listquo,coeffP,coeffR,t1,debut,contextptr);
    res -= tmp;
    // tmp = (-csturm_vertex_a(S,R,t0,1,contextptr)+csturm_vertex_a(S,R,t1,-1,contextptr));
    // res += tmp;
    res=(debut?1:signe)*res;
    if (debug_infolevel)
      *logptr(contextptr)  << "segment " << a << ".." << b << " index contribution " << res << endl;
    return res;
  }

  bool csturm_square_seq(const modpoly & P,const gen & a0,const gen & b0,const gen & a1,const gen & b1,gen & pgcd,vecteur & realroots,double eps,vecteur & seq1,vecteur & seq2,vecteur & seq3,vecteur & seq4,vecteur & horiz_sturm,vecteur & vert_sturm,GIAC_CONTEXT){
    gen A=a0+cst_i*b0,B=a1+cst_i*b0;
    vecteur rroots;
    seq1=csturm_segment_seq(P,A,B,rroots,eps,horiz_sturm,vert_sturm,contextptr);
    pgcd=seq1[4];
    if (!is_one(pgcd)){
      return false;
    }
    A=a1+cst_i*b0; B=a1+cst_i*b1;
    seq2=csturm_segment_seq(P,A,B,rroots,eps,horiz_sturm,vert_sturm,contextptr);
    pgcd=seq2[4];
    if (!is_one(pgcd)){
      return false;
    }
    A=a1+cst_i*b1; B=a0+cst_i*b1;
    seq3=csturm_segment_seq(P,A,B,rroots,eps,horiz_sturm,vert_sturm,contextptr);
    pgcd=seq3[4];
    if (!is_one(pgcd)){
      return false;
    }
    A=a0+cst_i*b1; B=a0+cst_i*b0;
    seq4=csturm_segment_seq(P,A,B,rroots,eps,horiz_sturm,vert_sturm,contextptr);
    pgcd=seq4[4];
    if (!is_one(pgcd)){
      return false;
    }
    realroots=mergevecteur(realroots,rroots);
    return true;
  }

  // find 2* number of roots of P inside the square of vertex of affixes a,b
  // roots on the square are not counted. P must not vanish at the vertices.
  // The complex Sturm sequences must be known
  int csturm_square(const modpoly & P,const gen & a0,const gen & b0,const gen & a1,const gen & b1,const vecteur & seq1,const vecteur & seq2,const vecteur & seq3,const vecteur & seq4,GIAC_CONTEXT){
    int ind,tmp;
    ind = 0;
    gen A=a0+cst_i*b0,B=a1+cst_i*b0;
    tmp = csturm_segment(seq1,A,B,contextptr);
    ind += tmp;
    A=a1+cst_i*b0; B=a1+cst_i*b1;
    tmp = csturm_segment(seq2,A,B,contextptr);
    ind += tmp;
    A=a1+cst_i*b1; B=a0+cst_i*b1;
    tmp = csturm_segment(seq3,A,B,contextptr);
    ind += tmp;
    A=a0+cst_i*b1; B=a0+cst_i*b0;
    tmp = csturm_segment(seq4,A,B,contextptr);
    ind += tmp;
    return ind;
  }

  void csturm_normalize(modpoly & p,const gen & a0,const gen & b0,const gen & a1,const gen & b1,vecteur & roots){
    int n=p.size()-1;
    // Make sure that x->a+i*x does not return a multiple 
    // of a real polynomial with the multiple non real
    // If degree of p is even the multiple will be a real (because of lcoeff)
    if (n%2){
      // If degree is odd then look at q=p(x-a_{n-1}/n*an)
      // it has the same property
      // if its cst coeff is zero remove
      gen an=p.front();
      gen b=p[1];
      gen shift=-b/n/an;
      modpoly q(taylor(p,shift));
      gen q0;
      // remove valuation
      int qs=q.size();
      int n1=0;
      for (;qs>0;--qs,++n1){
	if (!is_zero(q0=q[qs-1]))
	  break;
      }
      if (is_zero(re(q0,context0))){
	q=cst_i*q;
	p=cst_i*p;
      }
      if (n1){
	q=modpoly(q.begin(),q.begin()+qs);
	gen a=re(shift,context0),b=im(shift,context0);
	if (is_greater(a,a0,context0) && is_greater(b,b0,context0) && is_greater(a1,a,context0) && is_greater(b1,b,context0))
	  roots.push_back(makevecteur(shift,n1));
	p=taylor(q,-shift);
      }
    }
  }

  void ab2a0b0a1b1(const gen & a,const gen & b,gen & a0,gen & b0,gen & a1,gen & b1,GIAC_CONTEXT){
    a0=re(a,contextptr); b0=im(a,contextptr);
    a1=re(b,contextptr); b1=im(b,contextptr);
    if (ck_is_greater(a0,a1,contextptr)) std::swap<gen>(a0,a1);
    if (ck_is_greater(b0,b1,contextptr)) std::swap<gen>(b0,b1);
  }

  // find 2* number of roots of P inside the square of vertex of affixes a,b
  // excluding those on the square
  // returns -1 on error
  int csturm_square(const gen & p,const gen & a,const gen & b,gen& pgcd,GIAC_CONTEXT){
    if (p.type==_POLY){
      int res=0;
      factorization f(sqff(*p._POLYptr));
      factorization::const_iterator it=f.begin(),itend=f.end();
      for (;it!=itend;++it){
	int tmp=csturm_square(polynome2poly1(it->fact),a,b,pgcd,contextptr);
	if (tmp==-1)
	  return -1;
	res += it->mult*tmp;
      }
      return res;
    }
    if (p.type!=_VECT)
      return 0;
    modpoly P=*p._VECTptr;
    vecteur realroots;
    gen a0,b0,a1,b1;
    ab2a0b0a1b1(a,b,a0,b0,a1,b1,contextptr);
    csturm_normalize(P,a0,b0,a1,b1,realroots);
    int evident=0;
    if (!realroots.empty()){
      gen r=realroots.front();
      if (r.type==_VECT && r._VECTptr->size()==2)
	r=r._VECTptr->front();
      gen rx=re(r,contextptr),ry=im(r,contextptr);
      if ( ( is_zero(ry) && (rx==a0 || rx==a1) ) ||
	   ( is_zero(rx) && (ry==b0 || ry==b1) ) )
	;
      else
	evident=1;
    }
    if (P.size()<2)
      return evident;
    vecteur seq1,seq2,seq3,seq4,horiz_seq,vert_seq;    
    if (!csturm_square_seq(P,a0,b0,a1,b1,pgcd,realroots,0.0,seq1,seq2,seq3,seq4,horiz_seq,vert_seq,contextptr)){
      if (pgcd.type!=_VECT)	      
	return -1;
      modpoly g=(*pgcd._VECTptr)/pgcd[0];
      // true factorization found, restart with each factor
      modpoly p1=P/g;
      int n1=csturm_square(p1,a,b,pgcd,contextptr);
      if (n1==-1)
	return -1;
      int n2=csturm_square(g,a,b,pgcd,contextptr);
      if (n2==-1)
	return -1;
      return evident+n1+n2;
    }
    return evident+csturm_square(P,a0,b0,a1,b1,seq1,seq2,seq3,seq4,contextptr);
  }

  void complex_roots(const modpoly & P,const gen & a0,const gen & b0,const gen & a1,const gen & b1,vecteur & realroots,vecteur & complexroots,double eps);

  void complex_roots_split(const modpoly & P,const gen & pgcd,const gen & a0,const gen & b0,const gen & a1,const gen & b1,vecteur & realroots,vecteur & complexroots,double eps){
    if (pgcd.type!=_VECT)
      setsizeerr("csturm.cc"+pgcd.print());
    modpoly g=(*pgcd._VECTptr)/pgcd[0];
    // true factorization found, restart with each factor
    modpoly p1=P/g;
    csturm_normalize(p1,a0,b0,a1,b1,realroots);
    csturm_normalize(g,a0,b0,a1,b1,realroots);
    complex_roots(p1,a0,b0,a1,b1,realroots,complexroots,eps);
    complex_roots(g,a0,b0,a1,b1,realroots,complexroots,eps);
  }

  // Find complex roots of P in a0,b0 -> a1,b1
  void complex_roots(const modpoly & P,const gen & a0,const gen & b0,const gen & a1,const gen & b1,const vecteur & seq1,const vecteur & seq2,const vecteur & seq3,const vecteur & seq4,vecteur & realroots,vecteur & complexroots,double eps,vecteur & horiz_sturm,vecteur & vert_sturm){
    int n=csturm_square(P,a0,b0,a1,b1,seq1,seq2,seq3,seq4,context0);
    if (debug_infolevel && n)
      cerr << a0 << "," << b0 << ".." << a1 << "," << b1 << ":" << n/2 << endl;
    if (!n)
      return;
    if (eps<=0)
      setsizeerr("Bad precision "+print_DOUBLE_(eps,14));
    if (is_strictly_greater(eps,a1-a0,context0) && is_strictly_greater(eps,b1-b0,context0)){
      gen r(makevecteur(a0+cst_i*b0,a1+cst_i*b1));
      complexroots.push_back(makevecteur(r,gen(n)/2));
      return;
    }
    gen a01=(a0+a1)/2,b01=(b0+b1)/2,pgcd;
    vecteur seqvert,seqhoriz;
    gen A=a0+cst_i*b01,B=a1+cst_i*b01;
    seqhoriz=csturm_segment_seq(P,A,B,realroots,eps,horiz_sturm,vert_sturm,context0);
    pgcd=seqhoriz[4];
    if (is_one(pgcd)){
      A=a01+cst_i*b0; B=a01+cst_i*b1;
      seqvert=csturm_segment_seq(P,A,B,realroots,eps,horiz_sturm,vert_sturm,context0);
      pgcd=seqvert[4];
    }
    if (!is_one(pgcd)){
      complex_roots_split(P,pgcd,a0,b0,a1,b1,realroots,complexroots,eps);
      return;
    }
    /*
      (a0,b1)  - (a01,b1)  -  (a1,b1)         seq3                seq3
         |     n4   |     n3    |        seq4   n4      seqvert    n3     seq2
      (a0,b01) - (a01,b01) -  (a1,b01)        seqhoriz           seqhoriz
         |     n1   |     n2    |        seq4   n1      seqvert   n2      seq2
      (a0,b0)  - (a01,b0)  -  (a1,b0)         seq1                seq1
     */
    complex_roots(P,a0,b0,a01,b01,seq1,seqvert,seqhoriz,seq4,realroots,complexroots,eps,horiz_sturm,vert_sturm);
    complex_roots(P,a01,b0,a1,b01,seq1,seq2,seqhoriz,seqvert,realroots,complexroots,eps,horiz_sturm,vert_sturm);
    complex_roots(P,a01,b01,a1,b1,seqhoriz,seq2,seq3,seqvert,realroots,complexroots,eps,horiz_sturm,vert_sturm);
    complex_roots(P,a0,b01,a01,b1,seqhoriz,seqvert,seq3,seq4,realroots,complexroots,eps,horiz_sturm,vert_sturm);
  }

  // Find complex roots of P in a0,b0 -> a1,b1
  void complex_roots(const modpoly & P,const gen & a0,const gen & b0,const gen & a1,const gen & b1,vecteur & realroots,vecteur & complexroots,double eps){
    if (P.size()<2)
      return;
    vecteur Seq1,Seq2,Seq3,Seq4,horiz_sturm,vert_sturm;
    gen pgcd;
    if (!csturm_square_seq(P,a0,b0,a1,b1,pgcd,realroots,eps,Seq1,Seq2,Seq3,Seq4,horiz_sturm,vert_sturm,context0))
      complex_roots_split(P,pgcd,a0,b0,a1,b1,realroots,complexroots,eps);
    else
      complex_roots(P,a0,b0,a1,b1,Seq1,Seq2,Seq3,Seq4,realroots,complexroots,eps,horiz_sturm,vert_sturm);
  }

  // Find complex roots of P in a0,b0 -> a1,b1
  bool complex_roots(const modpoly & P,const gen & a0,const gen & b0,const gen & a1,const gen & b1,gen & pgcd,vecteur & roots,double eps){
    vecteur realroots,complexroots;
    complex_roots(P,a0,b0,a1,b1,realroots,complexroots,eps);
    roots=mergevecteur(roots,mergevecteur(realroots,complexroots));
    return true;
  }

  vecteur crationalroot(polynome & p,bool complexe){
    vectpoly v;
    int i=1;
    polynome qrem;
    environment * env= new environment;
    env->complexe=complexe || !is_zero(im(p,context0));
    vecteur w;
    do_linearfind(p,env,qrem,v,w,i);
    delete env;
    p=qrem;
    return w;
  }

  vecteur keep_in_rectangle(const vecteur & croots,const gen A0,const gen & B0,const gen & A1,const gen & B1,bool embed,GIAC_CONTEXT){
    vecteur roots;
    const_iterateur it=croots.begin(),itend=croots.end();
    for (;it!=itend;++it){
      gen a=re(*it,contextptr),b=im(*it,contextptr);
      if (is_greater(a,A0,contextptr)&&is_greater(A1,a,contextptr)&&is_greater(b,B0,contextptr)&&is_greater(B1,b,contextptr))
	roots.push_back(embed?makevecteur(*it,1):*it);
    }
    return roots;
  }

  // find roots of polynomial P at precision eps using complex Sturm sequences
  // P must have numeric coefficients, in Q[i]
  vecteur complex_roots(const modpoly & P,const gen & a0,const gen & b0,const gen & a1,const gen & b1,bool complexe,double eps){
    if (P.empty())
      return P;
    eps=std::abs(eps);
    bool aplati=(a0==a1) && (b0==b1);
    if (!aplati && complexe && (a0==a1 || b0==b1) )
      setsizeerr("Square is flat!");
    gen A0(a0),B0(b0),A1(a1),B1(b1);
    // initial rectangle: |roots|< 1+ max(|a_i|)/|a_n|
    if (aplati){
      gen maxai=_max(*apply(P,abs,context0)._VECTptr,context0);
      gen tmp=1+maxai/abs(P.front(),context0);
      A0=-tmp;
      B0=-tmp;
      A1=tmp;
      B1=tmp;
    }
    gen tmp;
    modpoly p(*apply(P,exact,context0)._VECTptr);
    lcmdeno(p,tmp,context0);
    polynome pp(poly12polynome(p));
    if (!complexe){
      gen tmp=gcd(re(pp,context0),im(pp,context0));
      if (tmp.type!=_POLY)
	return vecteur(0);
      pp=*tmp._POLYptr;
    }
    vecteur croots=crationalroot(pp,complexe);
    vecteur roots=keep_in_rectangle(croots,A0,B0,A1,B1,true,context0);
    p=polynome2poly1(pp);
    gen an=p.front();
    if (!is_zero(im(an,context0)))
      p=conj(p.front(),context0)*p;
    if (!complexe){ // real root isolation
      modpoly R=p;
      modpoly S=derivative(R);
      vecteur listquo,coeffP,coeffR;
      csturm_seq(S,R,listquo,coeffP,coeffR,context0);
      csturm_realroots(S,R,listquo,coeffP,coeffR,0,1,A0,A1,roots,eps,context0);
      return roots;
    }
    csturm_normalize(p,A0,B0,A1,B1,roots);
    gen pgcd;
    if (!complex_roots(p,A0,B0,A1,B1,pgcd,roots,eps))
      setsizeerr();
    return roots;
  }


  gen complexroot(const gen & g,bool complexe,GIAC_CONTEXT){
    if (g.type!=_VECT || (g._VECTptr->size()<2) )
      setsizeerr();
    vecteur & v(*g._VECTptr);
    gen p=v.front(),prec=evalf_double(v[1],1,contextptr);
    if (prec.type!=_DOUBLE_)
      settypeerr();
    double eps=prec._DOUBLE_val;
    unsigned vs=v.size();
    gen A(0),B(0);
    if (vs>3){
      A=v[2];
      B=v[3];
    }
    gen a0=re(A,contextptr),b0=im(A,contextptr),a1=re(B,contextptr),b1=im(B,contextptr);
    if (is_greater(a0,a1,contextptr))
      std::swap<gen>(a0,a1);
    if (is_greater(b0,b1,contextptr))
      std::swap<gen>(b0,b1);
    if (p.type==_VECT)
      return complex_roots(*p._VECTptr,a0,b0,a1,b1,complexe,eps);
    vecteur l;
    lvar(p,l);
    if (l.size()!=1)
      settypeerr();
    gen px=_e2r(makevecteur(p,l),contextptr);
    if (px.type==_FRAC)
      px=px._FRACptr->num;
    if (px.type!=_POLY)
      return vecteur(0);
    factorization f(sqff(*px._POLYptr));
    factorization::const_iterator it=f.begin(),itend=f.end();
    vecteur res;
    for (;it!=itend;++it){
      vecteur tmp=complex_roots(polynome2poly1(it->fact),a0,b0,a1,b1,complexe,eps);
      iterateur jt=tmp.begin(),jtend=tmp.end();
      for (;jt!=jtend;++jt){
	if (jt->type==_VECT && jt->_VECTptr->size()==2)
	  jt->_VECTptr->back()=it->mult*jt->_VECTptr->back();
      }
      res=mergevecteur(res,tmp);
    }
    return res;
  }

  gen _complexroot(const gen & g,GIAC_CONTEXT){
    return complexroot(g,true,contextptr);
  }
  const string _complexroot_s("complexroot");
  unary_function_eval __complexroot(&giac::_complexroot,_complexroot_s);
  unary_function_ptr at_complexroot (&__complexroot,0,true);

  gen _realroot(const gen & g,GIAC_CONTEXT){
    return complexroot(g,false,contextptr);
  }
  const string _realroot_s("realroot");
  unary_function_eval __realroot(&giac::_realroot,_realroot_s);
  unary_function_ptr at_realroot (&__realroot,0,true);

  vecteur crationalroot(const gen & g0,bool complexe){
    gen g(g0),a,b;
    if (g.type==_VECT){
      if (g.subtype==_SEQ__VECT){
	vecteur & tmp=*g._VECTptr;
	if (tmp.size()!=3)
	  setdimerr();
	g=tmp[0];
	a=tmp[1];
	b=tmp[2];
      }
      else {
	g=poly12polynome(*g._VECTptr);
      }
    }
    gen a0,b0,a1,b1;
    ab2a0b0a1b1(a,b,a0,b0,a1,b1,context0);
    vecteur l;
    lvar(g,l);
    if (l.size()!=1)
      settypeerr();
    gen px=_e2r(makevecteur(g,l),context0);
    if (px.type==_FRAC)
      px=px._FRACptr->num;
    if (px.type!=_POLY)
      return vecteur(0);
    factorization f(sqff(*px._POLYptr));
    factorization::const_iterator it=f.begin(),itend=f.end();
    vecteur res;
    for (;it!=itend;++it){
      polynome p=it->fact;
      vecteur tmp=crationalroot(p,complexe);
      res=mergevecteur(res,tmp);
    }
    if (a0!=a1 || b0!=b1)
      res=keep_in_rectangle(res,a0,b0,a1,b1,false,context0);
    return res;
  }
  gen _crationalroot(const gen & g,GIAC_CONTEXT){
    return crationalroot(g,true);
  }
  const string _crationalroot_s("crationalroot");
  unary_function_eval __crationalroot(&giac::_crationalroot,_crationalroot_s);
  unary_function_ptr at_crationalroot (&__crationalroot,0,true);

  gen _rationalroot(const gen & g,GIAC_CONTEXT){
    return crationalroot(g,false);
  }
  const string _rationalroot_s("rationalroot");
  unary_function_eval __rationalroot(&giac::_rationalroot,_rationalroot_s);
  unary_function_ptr at_rationalroot (&__rationalroot,0,true);

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC
