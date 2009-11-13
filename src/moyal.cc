// -*- mode:C++ ; compile-command: "g++-3.4 -I.. -g -c moyal.cc" -*-
#include "first.h"
/*
 *  Copyright (C) 2000,2007 B. Parisse, Institut Fourier, 38402 St Martin d'Heres
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
#include <stdexcept>
#include <cmath>
#include <cstdlib>
#include "sym2poly.h"
#include "usual.h"
#include "moyal.h"
#include "solve.h"
#include "intg.h"
#ifdef HAVE_LIBGSL
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_sf_erf.h>
#endif

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

#ifdef HAVE_SSCL
  gen moyal(const gen & a,const gen & b, const gen &vars,const gen & order){
    return symb_moyal(a,b,vars,order);
  }

#else // HAVE_SSCL
  gen moyal(const gen & a,const gen & b, const gen &vars,const gen & order){
    return symb_moyal(a,b,vars,order);
  }

#endif // HAVE_SSCL

  // "unary" version
  gen _moyal(const gen & args){
    int s=args._VECTptr->size();
    if (s!=4) setsizeerr("moyal.cc/_moyal");
    return moyal( (*(args._VECTptr))[0],(*(args._VECTptr))[1],(*(args._VECTptr))[2],(*(args._VECTptr))[3]);
  }
  const string _moyal_s("moyal");
  string texprintasmoyal(const gen & g,const string & s,GIAC_CONTEXT){
    return texprintsommetasoperator(g,"#",contextptr);
  }
  unary_function_unary __moyal(&_moyal,_moyal_s,0,&texprintasmoyal);
  unary_function_ptr at_moyal (&__moyal,0,true);

  gen Beta(const gen & a,const gen& b,GIAC_CONTEXT){
    return Gamma(a,contextptr)*Gamma(b,contextptr)/Gamma(a+b,contextptr);
  }
  gen _Beta(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      return symbolic(at_Beta,args);
    vecteur & v=*args._VECTptr;
    int s=v.size();
    if (s!=2)
      return symbolic(at_Beta,args);
    return Beta(v[0],v[1],contextptr);
  }
  const string _Beta_s("Beta");
  unary_function_eval __Beta(&_Beta,_Beta_s);
  unary_function_ptr at_Beta (&__Beta,0,true);

  gen Airy_Ai(const gen & x,GIAC_CONTEXT){
    gen e=x.evalf(1,contextptr);
#ifdef HAVE_LIBGSL
    if (e.type==_DOUBLE_)
      return gsl_sf_airy_Ai(e._DOUBLE_val,GSL_PREC_DOUBLE);
#endif
    return symbolic(at_Airy_Ai,x);
  }
  gen _Airy_Ai(const gen & args,GIAC_CONTEXT){
    return apply(args,Airy_Ai,contextptr);
  }
  const string _Airy_Ai_s("Airy_Ai");
  unary_function_eval __Airy_Ai(&_Airy_Ai,_Airy_Ai_s);
  unary_function_ptr at_Airy_Ai (&__Airy_Ai,0,true);

  gen Airy_Bi(const gen & x,GIAC_CONTEXT){
    gen e=x.evalf(1,contextptr);
#ifdef HAVE_LIBGSL
    if (e.type==_DOUBLE_)
      return gsl_sf_airy_Bi(e._DOUBLE_val,GSL_PREC_DOUBLE);
#endif
    return symbolic(at_Airy_Bi,x);
  }
  gen _Airy_Bi(const gen & args,GIAC_CONTEXT){
    return apply(args,Airy_Bi,contextptr);
  }
  const string _Airy_Bi_s("Airy_Bi");
  unary_function_eval __Airy_Bi(&_Airy_Bi,_Airy_Bi_s);
  unary_function_ptr at_Airy_Bi (&__Airy_Bi,0,true);

  gen _UTPN(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      return erfc(args/symbolic(at_sqrt,2),contextptr)/2;
    vecteur & v=*args._VECTptr;
    int s=v.size();
    if (s!=3 || is_zero(v[1]))
      setsizeerr();
    return erfc((v[2]-v[0])/sqrt(2*v[1],contextptr),contextptr)/2;
  }
  const string _UTPN_s("UTPN");
  unary_function_eval __UTPN(&_UTPN,_UTPN_s);
  unary_function_ptr at_UTPN (&__UTPN,0,true);

  gen randNorm(){
    /*
    double d=rand()/(RAND_MAX+1.0);
    d=2*d-1;
    identificateur x(" x");
    return newton(erf(x)-d,x,d);
    */
    double u=rand()/(RAND_MAX+1.0);
    double d=rand()/(RAND_MAX+1.0);
    return std::sqrt(-2*std::log(u))*std::cos(2*M_PI*d);
  }
  gen _randNorm(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT || args._VECTptr->size()!=2)
      setsizeerr();
    vecteur & v=*args._VECTptr;
    return evalf(v[0]+v[1]*randNorm(),1,contextptr);
  }
  const string _randNorm_s("randNorm");
  unary_function_eval __randNorm(&_randNorm,_randNorm_s);
  unary_function_ptr at_randNorm (&__randNorm,0,true);

  gen normald(const gen & m,const gen & s,const gen & x,GIAC_CONTEXT){
    gen v(s*s);
    return inv(sqrt(2*cst_pi*v,contextptr),contextptr)*exp(-pow(x-m,2)/(2*v),contextptr);
  }
  gen _normald(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT)
      return normald(0,1,g,contextptr);
    vecteur & v=*g._VECTptr;
    int s=v.size();
    if (s==3)
      return normald(v[0],v[1],v[2],contextptr);
    setsizeerr();
    return 0;
  }
  const string _normald_s("normald");
  unary_function_eval __normald(&_normald,_normald_s);
  unary_function_ptr at_normald (&__normald,0,true);

  // Normal cumulative distribution function
  // proba that X<x for X following a normal distrib of mean mean and dev dev
  // arg = vector [mean,dev,x] or x alone (mean=0, dev=1)
  gen normal_cdf(const gen & g,GIAC_CONTEXT){
    return rdiv(erf(plus_sqrt2_2*g,contextptr)+plus_one,2);
  }
  gen _normal_cdf(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT)
      return normal_cdf(g,contextptr);
    vecteur & v=*g._VECTptr;
    int s=v.size();
    if (s==2)
      return normal_cdf(v[1],contextptr)-normal_cdf(v[0],contextptr); 
    if (s==3)
      return normal_cdf((v[2]-v[0])/v[1],contextptr);
    if (s==4)
      return normal_cdf((v[3]-v[0])/v[1],contextptr)-normal_cdf((v[2]-v[0])/v[1],contextptr);
    setsizeerr();
    return 0;
  }
  const string _normal_cdf_s("normal_cdf");
  unary_function_eval __normal_cdf(&_normal_cdf,_normal_cdf_s);
  unary_function_ptr at_normal_cdf (&__normal_cdf,0,true);

  const string _normald_cdf_s("normald_cdf");
  unary_function_eval __normald_cdf(&_normal_cdf,_normald_cdf_s);
  unary_function_ptr at_normald_cdf (&__normald_cdf,0,true);

  // returns x s.t. UTPN(x) ~ y
  // ref Abramowitz & Stegun equation 26.2.22
  double utpn_initial_guess(double y){
    double t=std::sqrt(-2*std::log(y));
    t=t-(2.30753+.27061*t)/(1+0.99229*t+0.04481*t*t);
    return t;
  }
  gen utpn_inverse(double y,GIAC_CONTEXT){
    identificateur x(" x");
    return newton(erf(x/std::sqrt(2.0),contextptr)+2*y-1,x,utpn_initial_guess(y),NEWTON_DEFAULT_ITERATION,1e-5,1e-12,contextptr);
  }
  gen normal_icdf(const gen & g_orig,GIAC_CONTEXT){
    gen g=evalf(g_orig,1,contextptr);
    if (g.type!=_DOUBLE_ || g._DOUBLE_val<0 || g._DOUBLE_val>1)
      setsizeerr();
    if (g._DOUBLE_val==0)
      return minus_inf;
    if (g._DOUBLE_val==1)
      return plus_inf;
    return utpn_inverse(1-g._DOUBLE_val,contextptr);
    // identificateur x(" x");
    // plus_sqrt2*newton(erf(x)+1-2*g,x,0);
  }
  gen _normal_icdf(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT)
      return normal_icdf(g,contextptr);
    vecteur & v=*g._VECTptr;
    if (v.size()!=3)
      setsizeerr();
    return v[0]+v[1]*normal_icdf(v[2],contextptr);
  }
  const string _normal_icdf_s("normal_icdf");
  unary_function_eval __normal_icdf(&_normal_icdf,_normal_icdf_s);
  unary_function_ptr at_normal_icdf (&__normal_icdf,0,true);

  const string _normald_icdf_s("normald_icdf");
  unary_function_eval __normald_icdf(&_normal_icdf,_normald_icdf_s);
  unary_function_ptr at_normald_icdf (&__normald_icdf,0,true);

  gen binomial(const gen & n,const gen & k,const gen & p,GIAC_CONTEXT){
    return comb(n,k,contextptr)*pow(p,k,contextptr)*pow(1-p,n-k,contextptr);
  }
  gen _binomial(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT)
      setsizeerr();
    vecteur & v=*g._VECTptr;
    int s=v.size();
    if (s==2)
      return comb(v[0],v[1],contextptr);
    if (s==3)
      return binomial(v[0],v[1],v[2],contextptr);
    setsizeerr();
    return 0;
  }
  const string _binomial_s("binomial");
  unary_function_eval __binomial(&_binomial,_binomial_s);
  unary_function_ptr at_binomial (&__binomial,0,true);

  gen binomial_cdf(const gen & n,const gen &p,const gen & x,GIAC_CONTEXT){
    identificateur k(" k");
    return sum(binomial(n,k,p,contextptr),k,0,_floor(x,contextptr),contextptr);
  }
  gen _binomial_cdf(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT)
      setsizeerr();
    vecteur & v=*g._VECTptr;
    int s=v.size();
    if (s==3)
      return binomial_cdf(v[0],v[1],v[2],contextptr);
    if (s==4)
      return binomial_cdf(v[0],v[1],v[3],contextptr)-binomial_cdf(v[0],v[1],v[2],contextptr);
    setsizeerr();
    return 0;
  }
  const string _binomial_cdf_s("binomial_cdf");
  unary_function_eval __binomial_cdf(&_binomial_cdf,_binomial_cdf_s);
  unary_function_ptr at_binomial_cdf (&__binomial_cdf,0,true);

  gen binomial_icdf(const gen & n,const gen &p,const gen & x_orig,GIAC_CONTEXT){
    gen x=evalf(x_orig,1,contextptr);
    if (x._DOUBLE_val==0)
      return zero;
    if (x._DOUBLE_val==1)
      return n;
    if (n.type!=_INT_ || p.type!=_DOUBLE_ || x.type!=_DOUBLE_ || x._DOUBLE_val<0 || x._DOUBLE_val>1 )
      return symbolic(at_binomial_icdf,makevecteur(n,p,x));
    int k=0;
    gen b;
    for (;k<n.val;++k){
      b=b+binomial(n,k,p,contextptr);
      if (!ck_is_strictly_greater(x,b,contextptr))
	return k;
    }
    return n;
  }
  gen _binomial_icdf(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT)
      setsizeerr();
    vecteur & v=*g._VECTptr;
    int s=v.size();
    if (s==3)
      return binomial_icdf(v[0],v[1],v[2],contextptr);
    if (s==4)
      return binomial_icdf(v[0],v[1],v[3],contextptr)-binomial_icdf(v[0],v[1],v[2],contextptr);
    setsizeerr();
    return 0;
  }
  const string _binomial_icdf_s("binomial_icdf");
  unary_function_eval __binomial_icdf(&_binomial_icdf,_binomial_icdf_s);
  unary_function_ptr at_binomial_icdf (&__binomial_icdf,0,true);

  gen poisson(const gen & m,const gen & k,GIAC_CONTEXT){
    return exp(-m,contextptr)*pow(m,k,contextptr)/_factorial(k);
  }
  gen _poisson(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT)
      setsizeerr();
    vecteur & v=*g._VECTptr;
    int s=v.size();
    if (s==2)
      return poisson(v[0],v[1],contextptr);
    setsizeerr();
    return 0;
  }
  const string _poisson_s("poisson");
  unary_function_eval __poisson(&_poisson,_poisson_s);
  unary_function_ptr at_poisson (&__poisson,0,true);

  gen poisson_cdf(const gen & n,const gen & x,GIAC_CONTEXT){
    identificateur k(" k");
    return sum(poisson(n,k,contextptr),k,0,_floor(x,contextptr),contextptr);
  }
  gen _poisson_cdf(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT)
      setsizeerr();
    vecteur & v=*g._VECTptr;
    int s=v.size();
    if (s==2)
      return poisson_cdf(v[0],v[1],contextptr);
    if (s==3)
      return poisson_cdf(v[0],v[2],contextptr)-poisson_cdf(v[0],v[1],contextptr);
    setsizeerr();
    return 0;
  }
  const string _poisson_cdf_s("poisson_cdf");
  unary_function_eval __poisson_cdf(&_poisson_cdf,_poisson_cdf_s);
  unary_function_ptr at_poisson_cdf (&__poisson_cdf,0,true);

  gen poisson_icdf(const gen & m,const gen & t_orig,GIAC_CONTEXT){
    gen t=evalf(t_orig,1,contextptr);
    if (t._DOUBLE_val==0)
      return zero;
    if (t._DOUBLE_val==1)
      return plus_inf;
    if (t.type!=_DOUBLE_ || t._DOUBLE_val<0 || t._DOUBLE_val>1)
      return symbolic(at_poisson_icdf,makevecteur(m,t));
    int k=0;
    gen b;
    for (;;++k){
      b=b+poisson(m,k,contextptr);
      if (!ck_is_strictly_greater(t,b,contextptr))
	return k;
    }
    return t;
  }
  gen _poisson_icdf(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT)
      setsizeerr();
    vecteur & v=*g._VECTptr;
    int s=v.size();
    if (s==2)
      return poisson_icdf(v[0],v[1],contextptr);
    if (s==3)
      return poisson_icdf(v[0],v[2],contextptr)-poisson_icdf(v[0],v[1],contextptr);
    setsizeerr();
    return 0;
  }
  const string _poisson_icdf_s("poisson_icdf");
  unary_function_eval __poisson_icdf(&_poisson_icdf,_poisson_icdf_s);
  unary_function_ptr at_poisson_icdf (&__poisson_icdf,0,true);

  gen student(const gen & n,const gen & x,GIAC_CONTEXT){
    return Gamma(rdiv(n+1,2),contextptr)/Gamma(rdiv(n,2),contextptr)/sqrt(n*cst_pi,contextptr)*pow((1+pow(x,2)/n),-rdiv(n+1,2),contextptr);
  }
  gen _student(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT)
      setsizeerr();
    vecteur v=*g._VECTptr;
    int s=v.size();
    if (s==2){
      if (v[1].type==_DOUBLE_)
	return evalf(student(v[0],v[1],contextptr),1,contextptr);
      return student(v[0],v[1],contextptr);
    }
    setsizeerr();
    return 0;
  }
  const string _student_s("student");
  unary_function_eval __student(&_student,_student_s);
  unary_function_ptr at_student (&__student,0,true);

  gen zero_function(const gen & e ){
    return zero;
  }
  gen d2_UTPT(const gen & e,GIAC_CONTEXT ){
    return -_student(e,contextptr);
  }
  partial_derivative_onearg D_at_UTPT(d2_UTPT);
  unary_function_ptr _zero_function(unary_function_unary(&zero_function,""));
  unary_function_ptr D2_UTPT(unary_function_eval(&d2_UTPT,""));
  unary_function_ptr d_UTPT(int i){
    if (i==1)
      return _zero_function;
    if (i==2)
      return D2_UTPT;
    setsizeerr();
    return 0;
  }
  partial_derivative_multiargs D_UTPT(&d_UTPT);

  double FTS(int ndf,double cs2,double term, int j,double sum){
    for(;;){
      term=term*cs2*(ndf+j)/(j+2);
      double oldsum=sum;
      sum=sum+term;
      if (oldsum==sum)
	return sum;
      j +=2;
    }
  }
  double FCS(double fac,double term,int j,double res){
    for (;;){
      j -= 2;
      if (j<=1)
	return res;
      res = ((j-1)*fac*res)/j + term;
    }
  }

  double FSS(double ofs,double term,double fac,int j,double sum){
    for (;;){
      j -= 2;
      if (j<=1)
	return sum;
      sum=(ofs+j-2)/j*fac*sum+term;
    }
  }
  double FSS2(double ofs,double term,double fac,int j,double sum){
    for (;;){
      j -= 2;
      if (j<=0)
	return sum;
      sum=(ofs+j)/j*fac*sum+term;
    }
  }
  double TTS(int dof,int j,double term,double res,double cs2){
    for (;;){
      j +=2;
      term= (term*cs2*(j-1))/j;
      double ores= res;
      res= res + term/(dof+j);
      if (ores==res) return res;
    }
  }
  gen UTPT(const gen & n_orig,const gen & x0,GIAC_CONTEXT){
    gen n=_floor(n_orig,0);
    if (x0==plus_inf)
      return 0;
    if (x0==minus_inf)
      return 1;
    gen x1=evalf(x0,1,contextptr);
    if (n.type!=_INT_ || x1.type!=_DOUBLE_)
      return symbolic(at_UTPT,makevecteur(n,x0));
    int dof=n.val;
    double x=x1._DOUBLE_val,x2=x*x,y2= x2/dof;
    if (dof>=100){
      double y=std::log(y2)+1, a=dof-0.5, b=48*a*a;
      y=a*y;
      double res = (((((-.4*y - 3.3)*y -24)*y - 85.5)/(.8*y*y + 100 + b)+ y + 3)/b + 1)*std::sqrt(y);
      if (x<0)
	res=-res;
      return _UTPN(res,contextptr);
    }
    double y=std::sqrt(y2),b= 1+y2,cs2=1/b;
    if (x2<25){
      double res;
      if (dof==1)
	res=0;
      else
	res=FCS(cs2,y,dof,y);
      if (dof %2)
	res=2/M_PI*(std::atan(y)+cs2*res);
      else
	res=res*std::sqrt(cs2);
      if (x>0)
	return (1-res)/2;
      else
	return (1+res)/2;
    }
    else {
      double res= TTS(dof,0,dof,1,cs2);
      res= FCS(cs2,0,dof+2,res);
      if (dof %2)
	res= 2/M_PI*std::sqrt(cs2)*res;
      res /=2;
      if (x<0)
	res=1-res;
      return res;
    }
  }
  gen _UTPT(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      setsizeerr();
    vecteur & v=*args._VECTptr;
    int s=v.size();
    if (s!=2)
      setsizeerr();
    return UTPT(v[0],v[1],contextptr);
  }
  const string _UTPT_s("UTPT");
  unary_function_eval __UTPT(&_UTPT,&D_UTPT,_UTPT_s);
  unary_function_ptr at_UTPT (&__UTPT,0,true);

  // 26.7.5 in Abramowitz & Stegun
  double utpt_initial_guess(int n,double y,GIAC_CONTEXT){
    // double xp=utpn_initial_guess(y);
    double xp=utpn_inverse(y,contextptr)._DOUBLE_val;
    double xp2=xp*xp;
    double g1xp=xp*(xp2+1)/4;
    double g2xp=((5*xp2+16)*xp2+3)*xp/96;
    xp=xp+g1xp/n+g2xp/(n*n);
    return xp;
  }

  // dof=degree of freedom
  gen student_cdf(const gen & dof,const gen & x1,const gen & x2,GIAC_CONTEXT){
    gen x=evalf(x2,1,contextptr);
    return UTPT(dof,x1,contextptr)-UTPT(dof,x2,contextptr);
  }
  gen _student_cdf(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT)
      setsizeerr();
    vecteur & v=*g._VECTptr;
    int s=v.size();
    if (s==2)
      return student_cdf(v[0],minus_inf,v[1],contextptr);
    if (s==3)
      return student_cdf(v[0],v[1],v[2],contextptr);
    setsizeerr();
    return 0;
  }
  const string _student_cdf_s("student_cdf");
  unary_function_eval __student_cdf(&_student_cdf,_student_cdf_s);
  unary_function_ptr at_student_cdf (&__student_cdf,0,true);

  gen student_icdf(const gen & m,const gen & t_orig,GIAC_CONTEXT){
    gen t=evalf(t_orig,1,contextptr);
    if (t._DOUBLE_val==0)
      return zero;
    if (t._DOUBLE_val==1)
      return plus_inf;
    if (m.type!=_INT_ || t.type!=_DOUBLE_ || t._DOUBLE_val<0 || t._DOUBLE_val>1)
      return symbolic(at_student_icdf,makevecteur(m,t));
    double y=t._DOUBLE_val;
    double x0=utpt_initial_guess(m.val,1-y,contextptr);
    // return x0;
    // FIXME: use an iterative method to improve the initial guess
    identificateur x(" x");
    return newton(_student_cdf(makevecteur(m,x),contextptr)-y,x,x0,NEWTON_DEFAULT_ITERATION,1e-5,1e-12,contextptr);
  }
  gen _student_icdf(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT)
      setsizeerr();
    vecteur & v=*g._VECTptr;
    int s=v.size();
    if (s==2)
      return student_icdf(v[0],v[1],contextptr);
    if (s==3)
      return student_icdf(v[0],v[2],contextptr)-student_icdf(v[0],v[1],contextptr);
    setsizeerr();
    return 0;
  }
  const string _student_icdf_s("student_icdf");
  unary_function_eval __student_icdf(&_student_icdf,_student_icdf_s);
  unary_function_ptr at_student_icdf (&__student_icdf,0,true);

  gen chisquare(const gen & n,const gen & x,GIAC_CONTEXT){
    gen n2=n/2;
    return rdiv(pow(x,n2-1,contextptr)*exp(-x/2,contextptr),Gamma(n2,contextptr)*pow(2,n2,contextptr));
  }
  gen _chisquare(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT)
      setsizeerr();
    vecteur & v=*g._VECTptr;
    int s=v.size();
    if (s==2)
      return chisquare(v[0],v[1],contextptr);
    setsizeerr();
    return 0;
  }
  const string _chisquare_s("chisquare");
  unary_function_eval __chisquare(&_chisquare,_chisquare_s);
  unary_function_ptr at_chisquare (&__chisquare,0,true);

  gen d2_UTPC(const gen & e,GIAC_CONTEXT){
    return -_chisquare(e,contextptr);
  }
  partial_derivative_onearg D_at_UTPC(d2_UTPC);
  unary_function_ptr D2_UTPC(unary_function_eval(&d2_UTPC,""));
  unary_function_ptr d_UTPC(int i){
    if (i==1)
      return _zero_function;
    if (i==2)
      return D2_UTPC;
    setsizeerr();
    return 0;
  }
  partial_derivative_multiargs D_UTPC(&d_UTPC);
  gen UTPC(const gen & n_orig,const gen & x0,GIAC_CONTEXT){
    gen dof=_floor(n_orig,0);
    if (x0==plus_inf)
      return 0;
    if (is_zero(x0))
      return 1;
    gen x1=evalf(x0,1,contextptr);
    if (dof.type!=_INT_ || x1.type!=_DOUBLE_)
      return symbolic(at_UTPC,makevecteur(dof,x0));
    int n=dof.val;
    double x=x1._DOUBLE_val;
    if (x<0)
      return 1;
    if (x>10000)
      return 0.0;
    if (n<1)
      setsizeerr();
    if (n==1)
      return 2*_UTPN(sqrt(x,contextptr),contextptr);
    if (n>100){
    }
    double res=1;
    if (x>2){
      int r=n%2+2;
      res=std::exp(-x/2);
      double term = res;
      for (;r<n;r += 2){
	term = term*x/r;
	res += term;
      }
    }
    else {
      int r=n-2;
      for (;r>1;r-=2){
	res = res*x/r+1;
      }
      res *= std::exp(-x/2);
    }
    if (n%2)
      return sqrt(2*x/M_PI,contextptr)*res+2*_UTPN(sqrt(x,contextptr),contextptr);
    else
      return res;
  }
  gen _UTPC(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      setsizeerr();
    vecteur & v=*args._VECTptr;
    int s=v.size();
    if (s!=2)
      setsizeerr();
    return UTPC(v[0],v[1],contextptr);
  }
  const string _UTPC_s("UTPC");
  unary_function_eval __UTPC(&_UTPC,&D_UTPC,_UTPC_s);
  unary_function_ptr at_UTPC (&__UTPC,0,true);

  gen chisquare_cdf(const gen & dof,const gen & x1,const gen & x2,GIAC_CONTEXT){
    return UTPC(dof,x1,contextptr)-UTPC(dof,x2,contextptr);
  }
  gen _chisquare_cdf(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT)
      setsizeerr();
    vecteur & v=*g._VECTptr;
    int s=v.size();
    if (s==2)
      return chisquare_cdf(v[0],0,v[1],contextptr);
    if (s==3)
      return chisquare_cdf(v[0],v[1],v[2],contextptr);
    setsizeerr();
    return 0;
  }
  const string _chisquare_cdf_s("chisquare_cdf");
  unary_function_eval __chisquare_cdf(&_chisquare_cdf,_chisquare_cdf_s);
  unary_function_ptr at_chisquare_cdf (&__chisquare_cdf,0,true);

  // Abramowitz & Stegun 26.4.17
  double utpc_initial_guess(int n,double y,GIAC_CONTEXT){
    // double xp=utpn_initial_guess(y);
    if (n==2)
      return -2*std::log(y);
    if (n==1)
      y=y/2;
    double xp=utpn_inverse(y,contextptr)._DOUBLE_val;
    if (n==1)
      return xp*xp;
    double d=2/(9.0*n);
    d=1+xp*std::sqrt(d)-d;
    return n*d*d*d;
  }

  gen chisquare_icdf(const gen & m,const gen & t_orig,GIAC_CONTEXT){
    gen t=evalf(t_orig,1,contextptr);
    if (t._DOUBLE_val==0)
      return zero;
    if (t._DOUBLE_val==1)
      return plus_inf;
    if (m.type!=_INT_ || t.type!=_DOUBLE_ || t._DOUBLE_val<0 || t._DOUBLE_val>1)
      return symbolic(at_chisquare_icdf,makevecteur(m,t));
    // return utpc_initial_guess(m.val,1-t._DOUBLE_val);
    double x0=utpc_initial_guess(m.val,1-t._DOUBLE_val,contextptr);
    // FIXME
    identificateur x(" z");
    return newton(1-UTPC(m,x,contextptr)-t,x,x0,NEWTON_DEFAULT_ITERATION,1e-5,1e-12,contextptr);   
  }
  gen _chisquare_icdf(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT)
      setsizeerr();
    vecteur & v=*g._VECTptr;
    int s=v.size();
    if (s==2)
      return chisquare_icdf(v[0],v[1],contextptr);
    if (s==3)
      return chisquare_icdf(v[0],v[2],contextptr)-chisquare_icdf(v[0],v[1],contextptr);
    setsizeerr();
    return 0;
  }
  const string _chisquare_icdf_s("chisquare_icdf");
  unary_function_eval __chisquare_icdf(&_chisquare_icdf,_chisquare_icdf_s);
  unary_function_ptr at_chisquare_icdf (&__chisquare_icdf,0,true);

  gen snedecor(const gen & a,const gen & b,const gen & x,GIAC_CONTEXT){
    if (is_positive(-x,contextptr))
      return zero;
    return pow(a/b,a/2,contextptr)/Beta(a/2,b/2,contextptr) * pow(x,a/2-1,contextptr) * pow(1+a/b*x,-(a+b)/2,contextptr);
  }
  gen _snedecor(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT)
      setsizeerr();
    vecteur & v=*g._VECTptr;
    int s=v.size();
    if (s==3)
      return snedecor(v[0],v[1],v[2],contextptr);
    setsizeerr();
    return 0;
  }
  const string _snedecor_s("snedecor");
  unary_function_eval __snedecor(&_snedecor,_snedecor_s);
  unary_function_ptr at_snedecor (&__snedecor,0,true);

  const string _fisher_s("fisher");
  unary_function_eval __fisher(&_snedecor,_fisher_s);
  unary_function_ptr at_fisher (&__fisher,0,true);

  gen d2_UTPF(const gen & e ){
    return -_snedecor(e,0);
  }
  partial_derivative_onearg D_at_UTPF(d2_UTPF);
  unary_function_ptr D2_UTPF(unary_function_unary(&d2_UTPF,""));
  unary_function_ptr d_UTPF(int i){
    if (i<3)
      return _zero_function;
    if (i==3)
      return D2_UTPF;
    setsizeerr();
    return 0;
  }
  partial_derivative_multiargs D_UTPF(&d_UTPF);
  gen UTPF(const gen & num,const gen & den,const gen & x0,GIAC_CONTEXT){
    gen gndf=_floor(num,0),gddf=_floor(den,0),gx=evalf(x0,1,contextptr);
    if (gndf.type!=_INT_ || gddf.type!=_INT_ || gx.type!=_DOUBLE_)
      return symbolic(at_UTPF,makevecteur(num,den,x0));
    if (gx._DOUBLE_val<=0)
      return plus_one;
    int ndf=gndf.val,ddf=gddf.val;
    double x=gx._DOUBLE_val;
    if (ndf<1 || ddf <1 || ndf>300 || ddf>300)
      setdimerr();
    if (ndf==1)
      return 2*UTPT(ddf,std::sqrt(x),contextptr);
    double y= (x*ndf)/ddf, sn2= y/(1+y), cs2= 1/(1+y),sum;
    y=std::sqrt(y);
    if (ndf%2){ 
      if (ddf%2){ // ndf && ddf odd
	if (y<1){
	  if (ddf==1) sum=0; else sum=1;
	  sum=y*cs2*FCS(cs2,1,ddf,sum)+std::atan(y);
	  sum=1-2/M_PI*sum;
	  double sumB=ddf*FSS(ddf,1,sn2,ndf,1);
	  sumB=FCS(cs2,0,ddf+2,sumB);
	  return sum+2/M_PI*cs2*y*sumB;
	}
	else { // y>=1
	  sum= FTS(ndf,cs2,1,ddf,1);
	  sum= FSS2(ddf,0,sn2,ndf,sum);
	  sum= FCS(cs2,0,ddf+2,sum);
	  return 2/M_PI*y*cs2*sum;
	}
      } // end ndf odd , ddf odd
      else { // ndf odd, ddf even
	if (y<1)
	  return 1-FSS(ndf,1,cs2,ddf,1)*std::pow(sn2,ndf/2.0);
	else {
	  sum= FTS(ndf,cs2,1,ddf,1);
	  sum= FSS(ndf,0,cs2,ddf+2,sum);
	  return sum*std::pow(sn2,ndf/2.0);
	}
      }
    } // end ndf odd
    else { // ndf even
      return FSS(ddf,1,sn2,ndf,1)*std::pow(cs2,ddf/2.0);
    }
  }
  gen _UTPF(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      setsizeerr();
    vecteur & v=*args._VECTptr;
    int s=v.size();
    if (s!=3)
      setsizeerr();
    return UTPF(v[0],v[1],v[2],contextptr);
  }
  const string _UTPF_s("UTPF");
  unary_function_eval __UTPF(&_UTPF,&D_UTPF,_UTPF_s);
  unary_function_ptr at_UTPF (&__UTPF,0,true);

  gen snedecor_cdf(const gen & ndof,const gen & ddof,const gen & x,GIAC_CONTEXT){
    return 1-UTPF(ndof,ddof,x,contextptr);
  }
  gen _snedecor_cdf(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT)
      setsizeerr();
    vecteur & v=*g._VECTptr;
    int s=v.size();
    if (s==3)
      return snedecor_cdf(v[0],v[1],v[2],contextptr);
    setsizeerr();
    return 0;
  }
  const string _snedecor_cdf_s("snedecor_cdf");
  unary_function_eval __snedecor_cdf(&_snedecor_cdf,_snedecor_cdf_s);
  unary_function_ptr at_snedecor_cdf (&__snedecor_cdf,0,true);
  const string _fisher_cdf_s("fisher_cdf");
  unary_function_eval __fisher_cdf(&_snedecor_cdf,_fisher_cdf_s);
  unary_function_ptr at_fisher_cdf (&__fisher_cdf,0,true);

  // Abramowitz & Stegun 26.6.16
  double utpf_initial_guess(int num,int den,double y,GIAC_CONTEXT){
    if (num==1){
      double xp=utpt_initial_guess(den,y/2,contextptr);
      return xp*xp;
    }
    if (den==1){
      return y-0.5;
    }
    double xp=utpn_inverse(y,contextptr)._DOUBLE_val;
    double lambda=(xp*xp-3)/6;
    double h=2/fabs(1.0/(num-1)+1.0/(den-1)); // harmonic
    double w=xp*std::sqrt(h+lambda)/h-(lambda+5.0/6.0-2/(3*h))*fabs(1.0/(num-1)-1.0/(den-1));
    return std::exp(2*w);
  }

  gen snedecor_icdf(const gen & num,const gen & den,const gen & t_orig,GIAC_CONTEXT){
    gen t=evalf(t_orig,1,contextptr);
    if (t._DOUBLE_val==0)
      return zero;
    if (t._DOUBLE_val==1)
      return plus_inf;
    if (num.type!=_INT_ || den.type!=_INT_ || t.type!=_DOUBLE_ || t._DOUBLE_val<0 || t._DOUBLE_val>1)
      return symbolic(at_snedecor_icdf,makevecteur(num,den,t));
    // return utpf_initial_guess(num.val,den.val,1-t._DOUBLE_val);
    double x0=utpf_initial_guess(num.val,den.val,1-t._DOUBLE_val,contextptr);
    // FIXME
    identificateur x(" z");
    return newton(1-UTPF(num,den,x,contextptr)-t,x,x0,NEWTON_DEFAULT_ITERATION,1e-5,1e-12,contextptr);   
  }
  gen _snedecor_icdf(const gen & g,GIAC_CONTEXT){
    if (g.type!=_VECT)
      setsizeerr();
    vecteur & v=*g._VECTptr;
    int s=v.size();
    if (s==3)
      return snedecor_icdf(v[0],v[1],v[2],contextptr);
    if (s==4)
      return snedecor_icdf(v[0],v[1],v[3],contextptr)-snedecor_icdf(v[0],v[1],v[2],contextptr);
    setsizeerr();
    return 0;
  }
  const string _snedecor_icdf_s("snedecor_icdf");
  unary_function_eval __snedecor_icdf(&_snedecor_icdf,_snedecor_icdf_s);
  unary_function_ptr at_snedecor_icdf (&__snedecor_icdf,0,true);

  const string _fisher_icdf_s("fisher_icdf");
  unary_function_eval __fisher_icdf(&_snedecor_icdf,_fisher_icdf_s);
  unary_function_ptr at_fisher_icdf (&__fisher_icdf,0,true);


#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC
