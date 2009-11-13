/* -*- mode:C++ ; compile-command: "g++-3.4 -I.. -I../include -g -c ezgcd.cc" -*- */
#include "first.h"
/*  Multivariate GCD for large data not covered by the heuristic GCD algo
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
#include "ezgcd.h"
#include "sym2poly.h"
#include "gausspol.h"
#include "time.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  void add_dim(monomial<gen> & m,int d){
    index_t i(*m.index.iptr);
    for (int j=0;j<d;++j)
      i.push_back(0);
    m.index=i;
  }

  void change_dim(polynome & p,int dim){
    vector< monomial<gen> >::iterator it=p.coord.begin(),itend=p.coord.end();
    if (p.dim>=dim){
      p.dim=dim;
      for (;it!=itend;++it){
	index_t & i=*it->index.iptr;
	it->index=index_t(i.begin(),i.begin()+dim);
      }
      return;
    }
    int delta_dim=dim-p.dim;
    p.dim=dim;
    for (;it!=itend;++it)
      add_dim(*it,delta_dim);
  }

  // returns q such that p=q [degree] and q has only terms of degree<degree
  // p=q[N] means that p-q vanishes at v at order N
  polynome reduce(const polynome & p,const vecteur & v,int degree){
    int vsize=v.size();
    if (!vsize)
      return p;
    if (v==vecteur(vsize)){ 
      // trivial reduction, remove all terms of total deg >= degree
      polynome res(p.dim);
      vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
      for (;it!=itend;++it){
	if (total_degree(*it->index.iptr)<degree)
	  res.coord.push_back(*it);
      }
      return res;
    }
    if (degree<=1){
      gen res=peval(p,v,0);
      if (is_zero(res))
	return polynome(p.dim);
      if (res.type==_POLY){
	polynome resp(*res._POLYptr);
	change_dim(resp,p.dim);
	return resp;
      }
      else
	return polynome(monomial<gen>(res,0,p.dim));
    }
    polynome pcur(p);
    polynome y(monomial<gen>(plus_one,1,1,p.dim));
    if (!is_zero(v.front()))
      y.coord.push_back(monomial<gen>(-v.front(),0,1,p.dim));
    polynome quo(y.dim),rem(y.dim);
    pcur.TDivRem1(y,quo,rem);
    rem=reduce(rem.trunc1(),vecteur(v.begin()+1,v.end()),degree);
    quo=reduce(quo,v,degree-1);
    return quo*y+rem.untrunc1();
  }

  // Same as reduce but do it for every coefficient of p with
  // respect to the main variable
  polynome reduce_poly(const polynome & p,const vecteur & v,int degree){
    vector< monomial<gen> >::const_iterator it=p.coord.begin(),itend=p.coord.end();
    polynome res(p.dim);
    for (;it!=itend;){
      int d=it->index.iptr->front();
      polynome tmp(Tnextcoeff<gen>(it,itend));
      res=res+reduce(tmp,v,degree).untrunc1(d);
    }
    return res;
  }

  // reduce_divrem does a mixed division: euclidean w.r.t. the first var
  // and ascending power of X-b for the other vars
  // FIXME: this implementation does not work currently, except if other
  // depends only on the first var
  bool reduce_divrem2(const polynome & a,const polynome & other,const vecteur & v,int n,polynome & quo,polynome & rem,bool allowrational=false) {
    int asize=(a).coord.size();
    if (!asize){
      quo=a;
      rem=a; 
      return true;
    }
    int bsize=other.coord.size();
    if (bsize==0)  setsizeerr("ezgcd.cc/reduce_divrem2");
    index_m a_max = a.coord.front().index;
    index_m b_max = other.coord.front().index;
    quo.coord.clear();
    quo.dim=a.dim;
    rem.dim=a.dim;
    if ( (bsize==1) && (b_max==b_max*0) ){
      rem.coord.clear();
      gen b=other.coord.front().value;
      if (is_one(b))
	quo = a ;
      else {
	std::vector< monomial<gen> >::const_iterator itend=a.coord.end();
	for (std::vector< monomial<gen> >::const_iterator it=a.coord.begin();it!=itend;++it)
	  quo.coord.push_back(monomial<gen>(rdiv(it->value,b),it->index)); 
      }
      return true;
    }
    rem=a;
    if ( ! (a_max>=b_max) ){
      // test that the first power of a_max is < to that of b_max
      return (a_max.iptr->front()<b_max.iptr->front());
    }
    gen b(other.coord.front().value);
    while (a_max >= b_max){
      // errors should be trapped here and false returned if error occured
      gen q(rdiv(rem.coord.front().value,b));
      if (!allowrational){
	if ( has_denominator(q) || 
	     (!is_zero(q*b - rem.coord.front().value)) )
	  return false;
      }
      // end error trapping
      quo.coord.push_back(monomial<gen>(q,a_max-b_max));
      tensor<gen> temp=reduce_poly(other.shift(a_max-b_max,q),v,n);
      rem = rem-temp;
      if (rem.coord.size())
	a_max=rem.coord.front().index;
      else
	break;
    }
    return(true);    
  }

  bool reduce_divrem(const polynome & a,const polynome & other,const vecteur & v,int n,polynome & quo,polynome & rem) {
    quo.coord.clear();
    quo.dim=a.dim;
    rem.dim=a.dim;
    // if ( (a.dim<=1) || (a.coord.empty()) )
      return reduce_divrem2(a,other,v,n,quo,rem);
    std::vector< monomial<gen> >::const_iterator it=other.coord.begin();
    int bdeg=it->index.iptr->front(),rdeg;
    tensor<gen> b0(Tnextcoeff<gen>(it,other.coord.end()));
    tensor<gen> r(a),q(b0.dim);
    while ( (rdeg=r.lexsorted_degree()) >=bdeg){
      it=r.coord.begin();
      tensor<gen> a0(Tnextcoeff<gen>(it,r.coord.end())),tmp(a0.dim);
      // FIXME: should make ascending power division
      if (!reduce_divrem(a0,b0,v,n,q,tmp) || !tmp.coord.empty())
	return false;
      q=q.untrunc1(rdeg-bdeg);
      quo=quo+q;
      r=r-reduce_poly(q*other,v,n);
      if (r.coord.empty())
	return true;
    }
    return true;
  }

  // Hensel linear or quadratic lift
  // FIXME Quadratic lift currently works only if lcp is constant
  // Lift the equality p(b)=qb*rb [where b is a vecteur like for peval
  // assumed to have p.dim-1 coordinates] to p=q*r mod (X-b)^deg
  // Assuming that lcoeff(q)=lcp, lcoeff(r)=lcp, lcoeff(p)=lcp^2
  // If you want to find factors of a poly P such that P(b)=Qb*Rb, 
  // if lcp is the leading coeff of P
  // then p=P*lcp, qb=Qb*lcp(b)/lcoeff(Qb), rb=Rb*lcp(b)/lcoeff(Rb)
  bool hensel_lift(const polynome & p, const polynome & lcp, const polynome & qb, const polynome & rb, const vecteur & b,polynome & q, polynome & r,bool linear_lift,double maxop){
    if (maxop)
      linear_lift=true; // otherwise please adjust number of operations to do!
    double nop=0;
    int dim=p.dim;
    int deg=total_degree(p);
    if ( (qb.dim!=1) || (rb.dim!=1) || (dim==1) )
      setsizeerr("Bad dimension for qb or rb or b or degrees");
    polynome qu(1),ru(1),qbd(1);
    egcd(qb,rb,qu,ru,qbd);
    if (!Tis_constant(qbd))
      setsizeerr("qb and rb not prime together!");
    gen qrd(qbd.coord.front().value);
    // now we have qu*qb+ru*rb=qrd with 1-d polynomials
    change_dim(qu,dim);
    change_dim(ru,dim);
    // adjust dim & leading coeff of q and r by removing current leading coeff
    // and replace by lcp
    q=qb;
    r=rb;
    change_dim(q,dim);
    change_dim(r,dim);
    polynome q0(q),r0(r);
    vector<int> qshift(q.dim);
    qshift[0]=q.lexsorted_degree();
    q=q+(lcp-Tfirstcoeff<gen>(q)).shift(qshift);
    qshift[0]=r.lexsorted_degree();
    r=r+(lcp-Tfirstcoeff<gen>(r)).shift(qshift);    
    polynome p_qr(dim);
    for (int n=1;;){
      // qu*q+ru*r=qrd [n] (it's exact at the loop begin)
      // p=q*r [n] where [n] means of total valuation >= n
      // at the beginning n=1
      // enhanced at order 2*n by adding q',r' of valuation >=n
      // p-(q+q')*(r+r')=p-q*r - (r'q+q'r)-q'*r'
      // hence if we put r', q' such that p-q*r=(r'q+q'r) [2n]
      // we are done. Since p-q*r is of order [n], we get the solution
      // r'=qu*(p-qr)/qrd and q'=ru*(p-qr)/qrd
      if (debug_infolevel)
	cerr << "// Hensel " << n << " -> " << deg << endl;
      if (n>deg)
	return false;
      if (linear_lift)
	++n;
      else
	n=2*n;
      if (maxop>0){
	nop += double(q.coord.size())*r.coord.size();
	if (debug_infolevel)
	  cerr << "EZGCD " << nop << ":" << maxop << endl;
	if (nop>maxop)
	  return false;
      }
      p_qr=reduce_poly(p-q*r,b,deg);
      if (is_zero(p_qr))
	return true;
      if (n>deg)
	n=deg;
      p_qr=reduce_poly(p_qr,b,n);
      polynome qprime(reduce_poly(ru*p_qr,b,n)),qquo(qprime.dim),qrem(qprime.dim);
      polynome rprime(reduce_poly(qu*p_qr,b,n)),rquo(rprime.dim),rrem(qprime.dim);
      // reduction of qprime and rprime with respect to the main variable
      // we know that
      // (*) degree(p_qr) < degree(qr)
      // where degree is the degree wrt the main variable
      // since the leading coeffs of q and r are still adjusted
      // Then there is a unique solution to (*) with
      // degree(qprime)<degree(q), degree(rprime)<degree(r)
      if (linear_lift){
	reduce_divrem(qprime,q0,b,n,qquo,qrem);
	reduce_divrem(rprime,r0,b,n,rquo,rrem);
      }
      else {
	reduce_divrem(qprime,q,b,n,qquo,qrem);
	reduce_divrem(rprime,r,b,n,rquo,rrem);
      }
      // reduction of qprime and rprime with respect to the other variables
      // maybe we should check that q and r below have integer coeff
      q=q+inv(qrd,context0)*qrem; 
      r=r+inv(qrd,context0)*rrem;
      if (!linear_lift && (n<=deg/2)){
	// Now we modify qu and ru so that
	// (qu+qu')*q+(ru+ru')*r=qrd [2n]
	// therefore qu'*q+ru'*r=qrd-(qu*q+ru*r) [2n]
	// hence qu'= qu*[qrd-(qu*q+ru*r)]/qrd, ru'=ru*[qrd-(qu*q+ru*r)]/qrd
	p_qr=polynome(monomial<gen>(qrd,0,dim))-reduce_poly(qu*q+ru*r,b,n);
	qprime=reduce_poly(qu*p_qr,b,n);
	rprime=reduce_poly(ru*p_qr,b,n);
	reduce_divrem(qprime,r,b,n,qquo,qrem);
	reduce_divrem(rprime,q,b,n,rquo,rrem);
	qu=qu+inv(qrd,context0)*qrem; // should check that qu and ru have integer coeff
	ru=ru+inv(qrd,context0)*rrem;
      }
    }
  }

  // Replace the last coordinates of p with b instead of the first
  gen peval_back(const polynome & p,const vecteur & b){
    int pdim=p.dim,bdim=b.size();
    index_t cycle(pdim);
    int deltad=pdim-bdim;
    for (int i=0;i<bdim;++i)
      cycle[i]=i+deltad;
    for (int i=bdim;i<pdim;++i)
      cycle[i]=i-bdim;
    polynome pp(p);
    pp.reorder(cycle);
    int save=debug_infolevel;
    if (debug_infolevel)
      --debug_infolevel;
    gen res(peval(pp,b,0));
    debug_infolevel=save;
    return res;
  }

  polynome peval_1(const polynome & p,const vecteur &v,const gen & mod){
    if (p.dim!=signed(v.size()+1))
      setsizeerr();
    polynome res(1);
    index_t i(1);
    std::vector< monomial<gen> >::const_iterator it=p.coord.begin();
    std::vector< monomial<gen> >::const_iterator itend=p.coord.end();
    for (;it!=itend;){
      i[0]=it->index.iptr->front();
      polynome pactuel(Tnextcoeff<gen>(it,itend));
      gen g(peval(pactuel,v,mod));
      if ( (g.type==_POLY) && (g._POLYptr->dim==0) )
	g=g._POLYptr->coord.empty()?0:g._POLYptr->coord.front().value;
      if (!is_zero(g))
	res.coord.push_back(monomial<gen>(g,i));
    }
    return res;
  }

  // return true if a good eval point has been found
  bool find_good_eval(const polynome & F,const polynome & G,polynome & Fb,polynome & Gb,vecteur & b,bool debuglog,const gen & mod){
    int Fdeg=F.lexsorted_degree(),Gdeg=G.lexsorted_degree(),nvars=b.size();
    gen Fg,Gg;
    int essai=0;
    for (;;++essai){
      if (!is_zero(mod) && essai>mod.val)
	return false;
      if (debuglog)
	cerr << "Find_good_eval " << clock() << " " << b << endl;
      Fb=peval_1(F,b,mod);
      if (debuglog)
	cerr << "Fb= " << clock() << " " << gen(Fb) << endl;
      if (&F==&G)
	Gb=Fb;
      else {
	Gb=peval_1(G,b,mod);
      }
      if (debuglog)
	cerr << "Gb= " << clock() << " " << gen(Gb) << endl;
      if ( (Fb.lexsorted_degree()==Fdeg) && (Gb.lexsorted_degree()==Gdeg) ){
	if (debuglog)
	  cerr << "FOUND good eval" << clock() << " " << b << endl;
	return true;
      }
      b=vranm(nvars,0,0); // find another random point
    }
  }

  // It is probably required that 0 is a good evaluation point to
  // have an efficient algorithm
  // max_gcddeg is used when ezgcd was not successfull to find
  // the gcd even with 2 evaluations leading to the same gcd degree
  // in this case ezgcd calls itself with a bound on the gcd degree
  // is_sqff is true if we know that F_orig or G_orig is squarefree
  // is_primitive is true if F_orig and G_orig is primitive
  bool ezgcd(const polynome & F_orig,const polynome & G_orig,polynome & GCD,bool is_sqff,bool is_primitive,int max_gcddeg,double maxop){
    if (debug_infolevel)
      cerr << "// Starting EZGCD dimension " << F_orig.dim << endl;
    if (F_orig.dim<2)
      setsizeerr("Args must be multivariate polynomials");
    int Fdeg=F_orig.lexsorted_degree(),Gdeg=G_orig.lexsorted_degree();
    polynome F(F_orig.dim),G(F_orig.dim),cF(F_orig.dim),cG(F_orig.dim),cFG(F_orig.dim);
    if (is_primitive){
      cFG=polynome(monomial<gen>(plus_one,0,F_orig.dim));
      cF=cFG;
      cG=cFG;
      F=F_orig;
      G=G_orig;
    }
    else {
      cF=Tlgcd(F_orig);
      cG=Tlgcd(G_orig);
      cFG=gcd(cF.trunc1(),cG.trunc1()).untrunc1();
      F=F_orig/cF;
      G=G_orig/cG;
    }
    if (Tis_constant(F) || Tis_constant(G) ){
      GCD=cFG;
      return true;
    }
    polynome lcF(Tfirstcoeff(F)),lcG(Tfirstcoeff(G));
    double nop=lcF.coord.size()*F.coord.size()+lcG.coord.size()*G.coord.size();
    if (maxop>0){
      if (maxop<nop/10)
	return false;
    }
    vecteur b(F.dim-1);
    polynome Fb(1),Gb(1),Db(1);
    int old_gcddeg;
    for (;;){
      if (debug_infolevel)
	cerr << "// Back to EZGCD dimension " << F_orig.dim << endl;
      find_good_eval(F,G,Fb,Gb,b);
      Db=gcd(Fb,Gb);
      old_gcddeg=Db.lexsorted_degree();
      if (debug_infolevel)
	cerr << "// Eval at " << b << " gcd  degree " << old_gcddeg << endl;
      if (!old_gcddeg){
	GCD=cFG;
	return true;
      }
      if ( (!max_gcddeg) || (old_gcddeg<max_gcddeg) )
	break;
    }
    polynome new_Fb(1),new_Gb(1),quo(F.dim),rem(F.dim);
    for (;;){
      vecteur new_b(vranm(F.dim-1,0,0));
      find_good_eval(F,G,new_Fb,new_Gb,new_b);
      if (b==new_b)
	continue;
      polynome new_Db(gcd(new_Fb,new_Gb));
      int new_gcddeg=new_Db.lexsorted_degree();
      if (debug_infolevel)
	cerr << "// Eval at " << new_b << " gcd  degree " << new_gcddeg << endl;
      if (!new_gcddeg){
	GCD=cFG;
	return true;
      }
      if (new_gcddeg>old_gcddeg) // bad evaluation point
	continue;
      if (new_gcddeg==old_gcddeg) // might be a good guess!
	break;
      old_gcddeg=new_gcddeg;
      Db=new_Db;
      Fb=new_Fb;
      Gb=new_Gb;
      b=new_b;
    }
    // Found two times the same degree, try to lift!
    if ( (Fdeg<=Gdeg) && (old_gcddeg==Fdeg) ){
      if (G.TDivRem1(F,quo,rem) && rem.coord.empty()){
	GCD= F*cFG;
	return true;
      }
    }
    if ( (Gdeg<Fdeg) && (old_gcddeg==Gdeg)  ){
      if (G.TDivRem1(F,quo,rem) && rem.coord.empty()){
	GCD=G*cFG;
	return true;
      }
    }
    if (debug_infolevel)
      cerr << "// EZGCD degree " << old_gcddeg << endl;
    if ( (old_gcddeg==Fdeg) || (old_gcddeg==Gdeg) )
      return false;
    // this algo is fast if 0 is a good eval & the degree of the gcd is small
    if (!is_zero(b))
      return false;
    //if ( (old_gcddeg>4) && (old_gcddeg>Fdeg/4) && (old_gcddeg>Gdeg/4) )
    //  return false;
    polynome cofacteur(Fb/Db);
    if (Tis_constant(gcd(cofacteur,Db))){
      // lift Fb/Db *Db, more precisely insure that lc of each factor
      // is lcF(b)
      gen lcFb(peval_back(lcF,b));
      if (lcFb.type==_POLY)
	lcFb=lcFb._POLYptr->coord.front().value;
      Db=(lcFb*Db)/Db.coord.front().value;
      cofacteur=(lcFb*cofacteur)/cofacteur.coord.front().value;
      polynome liftF(F*lcF);
      polynome D(F_orig.dim),cofacteur_F(F_orig.dim),quo,rem;
      if (hensel_lift(liftF,lcF,cofacteur,Db,b,cofacteur_F,D,!Tis_constant(lcF),maxop) ){
	D=D/Tlgcd(D);
	if (F.TDivRem1(D,quo,rem) && is_zero(rem) && G.TDivRem1(D,quo,rem) && is_zero(rem)){
	  GCD=D*cFG;
	  return true;
	}
      }
      return false;
    }
    cofacteur=Gb/Db;
    if (Tis_constant(gcd(cofacteur,Db))){
      // lift Gb/Db *Db, more precisely insure that lc of each factor
      // is lcG(b)
      gen lcGb(peval_back(lcG,b));
      if (lcGb.type==_POLY)
	lcGb=lcGb._POLYptr->coord.front().value;
      Db=(lcGb*Db)/Db.coord.front().value;
      cofacteur=(lcGb*cofacteur)/cofacteur.coord.front().value;
      polynome liftG(G*lcG);
      polynome D(G_orig.dim),cofacteur_G(G_orig.dim),quo,rem;
      if (hensel_lift(liftG,lcG,cofacteur,Db,b,cofacteur_G,D,!Tis_constant(lcG),maxop) ){
	D=D/Tlgcd(D);
	if (F.TDivRem1(D,quo,rem) && is_zero(rem) && G.TDivRem1(D,quo,rem) && is_zero(rem)){
	  GCD=D*cFG;
	  return true;
	}
      }
      return false;
    }
    // FIXME find an integer j such that (F+jG)/D_b is coprime with D_b
    return false;
  }

  // algorithm=0 for HEUGCD, 1 for PRS, 2 for EZGCD, 3 for MODGCD
  gen heugcd_psrgcd_ezgcd_modgcd(const gen & args,int algorithm,GIAC_CONTEXT){
    vecteur & v=*args._VECTptr;
    gen p1(v[0]),p2(v[1]),n1,n2,d1,d2;
    vecteur lv;
    if ( (v.size()==3) && (v[2].type==_VECT) )
      lv=*v[2]._VECTptr;
    lvar(p1,lv);
    lvar(p2,lv);
    p1=e2r(p1,lv,contextptr);
    fxnd(p1,n1,d1);
    p2=e2r(p2,lv,contextptr);
    fxnd(p2,n2,d2);
    gen res,np_simp,nq_simp,d_content;
    polynome p,q,p_gcd;
    if ( (n1.type!=_POLY) || (n2.type!=_POLY) )
      res=gcd(n1,n2);
    else {
      polynome pres;
      bool result=false;
      switch(algorithm){
      case 0:
	p_gcd.dim=n1._POLYptr->dim;
	result=gcdheu(*n1._POLYptr,*n2._POLYptr,p,np_simp,q,nq_simp,p_gcd,d_content,true);
	pres=p_gcd*d_content;
	break;
      case 1:
	pres=gcdpsr(*n1._POLYptr,*n2._POLYptr);
	result=true;
	break;
      case 2:
	result=ezgcd(*n1._POLYptr,*n2._POLYptr,pres);
	break;
      case 3:
	result=gcd_modular_algo(*n1._POLYptr,*n2._POLYptr,pres,false);
	break;
      }
      if (result)
	res=pres;
      else
	setsizeerr("GCD not successfull");
    }
    return r2e(res,lv,contextptr);
  }

  gen _ezgcd(const gen & args,GIAC_CONTEXT){
    if ( (args.type!=_VECT) || (args._VECTptr->size()<2) )
      return symbolic(at_ezgcd,args);
    return heugcd_psrgcd_ezgcd_modgcd(args,2,contextptr);
  }
  const string _ezgcd_s("ezgcd");
  unary_function_eval __ezgcd(&giac::_ezgcd,_ezgcd_s);
  unary_function_ptr at_ezgcd (&__ezgcd,0,true);
  
  gen _modgcd(const gen & args,GIAC_CONTEXT){
    if ( (args.type!=_VECT) || (args._VECTptr->size()<2) )
      return symbolic(at_modgcd,args);
    return heugcd_psrgcd_ezgcd_modgcd(args,3,contextptr);
  }
  const string _modgcd_s("modgcd");
  unary_function_eval __modgcd(&giac::_modgcd,_modgcd_s);
  unary_function_ptr at_modgcd (&__modgcd,0,true);
  
  gen _heugcd(const gen & args,GIAC_CONTEXT){
    if ( (args.type!=_VECT) || (args._VECTptr->size()<2) )
      return symbolic(at_heugcd,args);
    return heugcd_psrgcd_ezgcd_modgcd(args,0,contextptr);
  }
  const string _heugcd_s("heugcd");
  unary_function_eval __heugcd(&giac::_heugcd,_heugcd_s);
  unary_function_ptr at_heugcd (&__heugcd,0,true);

  gen _psrgcd(const gen & args,GIAC_CONTEXT){
    if ( (args.type!=_VECT) || (args._VECTptr->size()<2) )
      return symbolic(at_psrgcd,args);
    return heugcd_psrgcd_ezgcd_modgcd(args,1,contextptr);
  }
  const string _psrgcd_s("psrgcd");
  unary_function_eval __psrgcd(&giac::_psrgcd,_psrgcd_s);
  unary_function_ptr at_psrgcd (&__psrgcd,0,true);


#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC
