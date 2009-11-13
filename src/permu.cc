/* -*- mode:C++ ; compile-command: "g++-3.4 -I.. -g -c permu.cc" -*- */
#include "first.h"
/*
 *  Copyright (C) 2005, 2007 R. De Graeve & B. Parisse, 
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
#include <fstream>
#include "permu.h"
#include "usual.h"
#include "sym2poly.h"
#include "rpn.h"
#include "prog.h"
#include "derive.h"
#include "subst.h"
#include "misc.h"
#include "plot.h"
#include "intg.h"
#include "ifactor.h"
#include "lin.h"
#include "modpoly.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  vecteur vector_int_2_vecteur(const vector<int> & v,GIAC_CONTEXT){
    //transforme un vector<int> en vecteur 
    vector<int>::const_iterator it=v.begin(),itend=v.end();
    vecteur res;
    res.reserve(itend-it);
    if (xcas_mode(contextptr)){
      for (;it!=itend;++it)
	res.push_back(*it+1);
    }
    else {
      for (;it!=itend;++it)
	res.push_back(*it);
    }
    return res;
  } 

  vector<int> vecteur_2_vector_int(const vecteur & v,GIAC_CONTEXT){
    //transforme un vecteur en vector<int>  
    vecteur::const_iterator it=v.begin(),itend=v.end();
    vector<int> res;
    res.reserve(itend-it);
    if (xcas_mode(contextptr)){
      for (;it!=itend;++it)
	if ((*it).type==_INT_) 
	  res.push_back((*it).val-1);
	else 
	  settypeerr();
    }
    else {
      for (;it!=itend;++it)
	if ((*it).type==_INT_) 
	  res.push_back((*it).val); 
	else 
	  settypeerr();
    }
    return res;
  } 

  vector< vector<int> > vecteur_2_vectvector_int(const vecteur & v,GIAC_CONTEXT){
    //transforme un vecteur en vector< vector<int> >  
    vecteur::const_iterator it=v.begin(),itend=v.end();
    vector< vector<int> > res;
    res.reserve(itend-it);
    for (;it!=itend;++it){
      if (it->type!=_VECT)
	setsizeerr();
      res.push_back(vecteur_2_vector_int(*it->_VECTptr,contextptr));
    }
    return res;
  }

  vecteur vectvector_int_2_vecteur(const vector< vector<int> > & v,GIAC_CONTEXT){
    //transforme un vector< vector<int> > en vecteur  
    int s=v.size();
    vecteur res;
    res.reserve(s);
    for (int i=0;i<s;++i)
      res.push_back(vector_int_2_vecteur(v[i],contextptr));
    return res;
  }

  vector<int> sizes(const vector< vector<int> > & v){
    //donne la liste des tailles des vecteurs qui forment v
    int s=v.size();
    vector<int> res(s);
    for (int i=0;i<s;i++){
      vector<int> vi;
      vi=v[i];
      res[i]=vi.size();
      //res.push_back(vi.size());pourqoi?
    }
    return res;
  }

  gen _sizes(const gen & args){
    if (args.type!=_VECT) 
      settypeerr();
    vecteur v(*args._VECTptr); 
    vecteur res;
    vecteur::const_iterator it=v.begin(),itend=v.end();
    res.reserve(itend-it);
    for (;it!=itend;++it){
      if (it->type!=_VECT)
	setsizeerr();
      res.push_back(int(it->_VECTptr->size()));
    }
    return res;
  } 
  const string _sizes_s("sizes");
  unary_function_unary __sizes(&_sizes,_sizes_s);
  unary_function_ptr at_sizes (&__sizes,0,true);

  int cyclesorder(const vector< vector<int> > & v,GIAC_CONTEXT){
    vector<int> ss;
    ss=sizes(v);   
    return(_lcm(vector_int_2_vecteur(ss,contextptr)).val);
  } 
  int permuorder(const vector<int> & p,GIAC_CONTEXT){
    vector< vector<int> > c;  
    c=permu2cycles(p);
    return(_lcm(vector_int_2_vecteur(sizes(c),contextptr)).val);
  }

  gen _permuorder(const gen & args,GIAC_CONTEXT){
    if  (args.type!=_VECT) 
      setsizeerr();
    vecteur v(*args._VECTptr); 
    vector<int> p;
    if (!is_permu(v,p,contextptr))
      setsizeerr();
    return (permuorder(p,contextptr));
  }
  const string _permuorder_s("permuorder");
  unary_function_eval __permuorder(&_permuorder,_permuorder_s);
  unary_function_ptr at_permuorder (&__permuorder,0,true);

  vector<int> randperm(const int & n){
    //renvoie une permutation au hasard de long n
    vector<int> p(n);
    vector<int> temp(n);
    for (int k=0;k<n;k++) {temp[k]=k;}
    //on chosit au hasard h et alors p[k]=temp[h]
    int m;
    m=n;
    for (int k=0;k<n;k++) {
      int h;
      h=int(m*(rand()/(RAND_MAX*1.0)));
      p[k]=temp[h]; m=m-1;
      //mise a jour de temp :il faut supprimer temp[h]
      for (int j=h;j<m;j++) {temp[j]=temp[j+1];
      }
    }    
    return(p); 
  }
  gen _randperm(const gen & args,GIAC_CONTEXT){
    if (!is_integer(args)) 
      setsizeerr();
    gen n=args;
    return vector_int_2_vecteur(randperm(n.val),contextptr);
  }
  const string _randperm_s("randperm");
  unary_function_eval __randperm(&_randperm,_randperm_s);
  unary_function_ptr at_randperm (&__randperm,0,true);
  
  bool is_permu(const vecteur &p,vector<int> & p1,GIAC_CONTEXT) {
    //renvoie true si p est une perm et transforme p en le vector<int> p1  
    int n;
    n=p.size();
    vector<int> p2(n);
    p1=p2;
    vector<int> temp(n);
   
    for (int j=0;j<n;j++){ if (p[j].type!=_INT_){return(false);}}
     
    for (int j=0;j<n;j++){
      if (xcas_mode(contextptr)>0) 
	p1[j]=p[j].val-1; 
      else 
	p1[j]=p[j].val;
      if ((n<=p1[j])|| (p1[j])<0) {
	return(false);
      }
    }
    int k;
    k=0;
    while (k<n) {
      int p1k=p1[k];
      if (p1k<0 || p1k>=n) {return(false);}
      if (temp[p1k]) {
	return(false);} 
      else {temp[p1k]=1;}
      k=k+1;
    }
    return(true);
  } 
  gen _is_permu(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT) 
      setsizeerr();
    vecteur v(*args._VECTptr);
    vector<int> p1;
    return is_permu(v,p1,contextptr);
  }
  const string _is_permu_s("is_permu");
  unary_function_eval __is_permu(&_is_permu,_is_permu_s);
  unary_function_ptr at_is_permu (&__is_permu,0,true);


  bool is_cycle(const vecteur & c,vector<int> & c1,GIAC_CONTEXT) {
    //renvoie true si c est un cycle et transf c vecteur en le cycle c1 vector<int>
    int n1;
    n1=c.size();
    //vector<int> p;
    vector<int> c2(n1);
    c1=c2;
    for (int j=0;j<n1;j++){
      if (xcas_mode(contextptr)>0) 
	c1[j]=c[j].val-1; 
      else 
	c1[j]=c[j].val;
    }   
    for (int k=0;k<n1-1;k++) {
      int ck=c1[k];
      if (ck<0) {return(false);}
      for (int j=k+1;j<n1;j++){
	if (ck==c1[j]) {return (false);}
      }
    }
    return (true);
  } 
  gen _is_cycle(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT) 
      setsizeerr();
    vecteur v(*args._VECTptr);
    vector<int> c1;
    return is_cycle(v,c1,contextptr);
  }
  const string _is_cycle_s("is_cycle");
  unary_function_eval __is_cycle(&_is_cycle,_is_cycle_s);
  unary_function_ptr at_is_cycle (&__is_cycle,0,true);
 
  vector<int> cycle2perm(const vector<int> & c) {
    //transforme c en la permu p et renvoie p
    int n1;
    n1=c.size();
    //vector<int> c1(n1);
    int n;
    n=c[0];
    for (int k=1;k<n1;k++) {
      if (n<c[k]) {
	n=c[k];
      }
    }
    n=n+1;   
    vector<int> p(n);
    for (int k=0;k<n;k++) {
      p[k]=k;}
    for (int k=0;k<n1-1;k++) {p[c[k]]=c[k+1];}
    p[c[n1-1]]=c[0];
    return(p);
  } 
  
  gen _cycle2perm(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT) 
      setsizeerr();
    vecteur v(*args._VECTptr); 
    vector<int> c;
    if (!is_cycle(v,c,contextptr))
      setsizeerr();
    return vector_int_2_vecteur(cycle2perm(c),contextptr);
  }
  const string _cycle2perm_s("cycle2perm");
  unary_function_eval __cycle2perm(&_cycle2perm,_cycle2perm_s);
  unary_function_ptr at_cycle2perm (&__cycle2perm,0,true);
 
 
  vector<int> p1op2(const vector<int> & p1,const vector<int> & p2) {
    //composition de 2 perm p1 et p2 de long n1 et n2 : a pour long max(n1,n2)
    int n1;
    n1=p1.size();
    int n2;
    n2=p2.size();
    vector<int> p3;
    vector<int> p4;
    p3=p1;
    p4=p2;   
    if (n1>n2) {
      for (int k=n2;k<n1;k++) {p4.push_back(k);n2=n1;} 
    } else {
      for (int k=n1;k<n2;k++) {p3.push_back(k);n1=n2;}
    }
    ;
    vector<int> p(n1);     
    for (int k=0;k<n1;k++) {
      p[k]=p3[p4[k]];
    }      
    return(p);
  }
 
  gen _p1op2(const gen & args,GIAC_CONTEXT){
    if ( (args.type!=_VECT)  || (args._VECTptr->size()!=2) )
      settypeerr();
    vecteur v(*args._VECTptr);
    gen v1=v.front(),v2=v.back();
    if ( (v1.type!=_VECT) || (v2.type!=_VECT))
      settypeerr();
    vector<int> p1,p2;   
    if (!is_permu(*v1._VECTptr,p1,contextptr) || !is_permu(*v2._VECTptr,p2,contextptr))
      setsizeerr();
    return vector_int_2_vecteur(p1op2(p1,p2),contextptr);
  }
  const string _p1op2_s("p1op2");
  unary_function_eval __p1op2(&_p1op2,_p1op2_s);
  unary_function_ptr at_p1op2 (&__p1op2,0,true);

  vector<int> c1oc2(const vector<int> & c1,const vector<int> & c2) {
    //composition de 2 cycles en une perm de long min
    vector<int> p1;
    p1=cycle2perm(c1);
    vector<int> p2;
    p2=cycle2perm(c2);
    int n1;
    n1=p1.size();
    int n2;
    n2=p2.size();
    int n;
    if (n1>n2) {
      n=n1;
      for (int k=n2;k<n;k++) {p2.push_back(k);} 
    } 
    else {
      n=n2;
      for (int k=n1;k<n;k++) {p1.push_back(k);} 
    }   
    vector<int> p(n);
    for (int k=0;k<n;k++) {
      p[k]=p1[p2[k]];
    }
    return(p);
  }
  
  gen _c1oc2(const gen & args,GIAC_CONTEXT){
    if ( (args.type!=_VECT)  || (args._VECTptr->size()!=2) )
      settypeerr();
    vecteur v(*args._VECTptr);
    gen v1=v.front(),v2=v.back();
    if ( (v1.type!=_VECT) || (v2.type!=_VECT))
      settypeerr();
    vector<int> c1,c2;
    //c1=vecteur_2_vector_int(*v1._VECTptr);
    //c2=vecteur_2_vector_int(*v2._VECTptr);
    if (!is_cycle(*v1._VECTptr,c1,contextptr) || !is_cycle(*v2._VECTptr,c2,contextptr))
      setsizeerr();
    return vector_int_2_vecteur(c1oc2(c1,c2),contextptr);
  }
  const string _c1oc2_s("c1oc2");
  unary_function_eval __c1oc2(&_c1oc2,_c1oc2_s);
  unary_function_ptr at_c1oc2 (&__c1oc2,0,true);
 
  vector<int> c1op2(const vector<int> & c1, const vector<int> & p2) {
    //composition d'un cycle et d'une perm en une perm de long min
    vector<int> p1,p3;    
    p1=cycle2perm(c1);
    int n1;
    n1=p1.size();
    int n2;
    n2=p2.size();
    p3=p2;
    int n;
    if (n1>n2) {
      n=n1;
      for (int k=n2;k<n;k++) {p3.push_back(k);} 
    } 
    else {
      n=n2;
      for (int k=n1;k<n;k++) {p1.push_back(k);} 
    }
    vector<int> p(n);
    for (int k=0;k<n;k++) {
      p[k]=p1[p3[k]];
    }
    return(p);
  }
  
  gen _c1op2(const gen & args,GIAC_CONTEXT){
    if ( (args.type!=_VECT)  || (args._VECTptr->size()!=2) )
      settypeerr();
    vecteur v(*args._VECTptr);
    gen v1=v.front(),v2=v.back();
    if ( (v1.type!=_VECT) || (v2.type!=_VECT))
      settypeerr();
    vector<int> c1,p2;
    if (!is_cycle(*v1._VECTptr,c1,contextptr) || !is_permu(*v2._VECTptr,p2,contextptr))
      setsizeerr();
    return vector_int_2_vecteur(c1op2(c1,p2),contextptr);
  }
  const string _c1op2_s("c1op2");
  unary_function_eval __c1op2(&_c1op2,_c1op2_s);
  unary_function_ptr at_c1op2 (&__c1op2,0,true);

  vector<int> p1oc2(const vector<int> & p1, const vector<int> & c2) {
    //composition d'une perm et d'un cycle en une perm de long min
    vector<int> p2,p3;
    p2=cycle2perm(c2);
    int n2;
    n2=p2.size();
    int n1;
    n1=p1.size();
    p3=p1;
    int n;
    if (n1>n2) {n=n1;
    for (int k=n2;k<n;k++) {p2.push_back(k);} 
    } else {n=n2;
    for (int k=n1;k<n;k++) {p3.push_back(k);} 
    }
    vector<int> p(n);
    for (int k=0;k<n;k++) {
      p[k]=p3[p2[k]];
    }
    return(p);
  }
  
  gen _p1oc2(const gen & args,GIAC_CONTEXT){
    if ( (args.type!=_VECT)  || (args._VECTptr->size()!=2) )
      settypeerr();
    vecteur v(*args._VECTptr);
    gen v1=v.front(),v2=v.back();
    if ( (v1.type!=_VECT) || (v2.type!=_VECT))
      settypeerr();
    vector<int> p1,c2;
    if (!is_cycle(*v2._VECTptr,c2,contextptr) || !is_permu(*v1._VECTptr,p1,contextptr))
      setsizeerr();
    return vector_int_2_vecteur(p1oc2(p1,c2),contextptr);
  }
  const string _p1oc2_s("p1oc2");
  unary_function_eval __p1oc2(&_p1oc2,_p1oc2_s);
  unary_function_ptr at_p1oc2 (&__p1oc2,0,true);

  vector<int> cycles2permu(const vector< vector<int> > & c) {
    //transforme une liste de cycles en la permutation produit
    int n;
    n=c.size();
    vector<int> pk;
    vector<int> p;
    vector<int> ck;
    vector<int> c1;
    ck=c[n-1];    
    vector<int> c0(1);
    c0[0]=0;
    p=c1oc2(ck,c0);
    for (int k=n-2;k>=0;k--){
      ck=c[k];
      p=c1op2(ck,p);
    } 
    return(p);
  }

  gen _cycles2permu(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT) 
      settypeerr();
    vecteur v(*args._VECTptr); 
    vecteur::const_iterator it=v.begin(),itend=v.end();
    for (;it!=itend;++it){
    vector<int> c;
    if (!is_cycle(*it->_VECTptr,c,contextptr))
      setsizeerr();
    }
    return vector_int_2_vecteur(cycles2permu(vecteur_2_vectvector_int(v,contextptr)),contextptr);
  }
  const string _cycles2permu_s("cycles2permu");
  unary_function_eval __cycles2permu(&_cycles2permu,_cycles2permu_s);
  unary_function_ptr at_cycles2permu (&__cycles2permu,0,true);

  vector< vector<int> > permu2cycles(const vector<int> & p) {
    //transforme la permutation p en une liste de cycles (p= produit de ces cycles)
    int l=p.size();
    int n=0;
    int deb;
    vector<int> p1(l);
    p1=p;
    vector<int>  temp(l+1);
    vector< vector<int> > c;
    if (p1[l-1]==l-1) {
      vector<int> c1;  
      c1.push_back(l-1);
      c.push_back(c1); l=l-1;
    }
    temp[l]=0;
    for (int k=0;k<l;k++) temp[k]=p1[k];
    while (n<l){ 
      vector<int> v;
      v.push_back(n);deb=n;
      while (p1[n]!=deb){
	v.push_back(p1[n]);
	temp[n]=0;
	n=p1[n];
      }
      if (n!=deb) {c.push_back(v);}
      temp[n]=0;
      n=deb+1;
      while ((n<l)&&(temp[n]==0)) n++;
    }
    return(c);
  }

  gen _permu2cycles(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT) 
      setsizeerr();
    vecteur v(*args._VECTptr);
    vector<int> p;
    if (!is_permu(v,p,contextptr))
      setsizeerr();
    return vectvector_int_2_vecteur(permu2cycles(p),contextptr);
  }

  const string _permu2cycles_s("permu2cycles");
  unary_function_eval __permu2cycles(&_permu2cycles,_permu2cycles_s);
  unary_function_ptr at_permu2cycles (&__permu2cycles,0,true);

  vector<int> perminv(const vector<int> & p){
   int n;
   n=p.size();
   vector<int> p1(n);
   for (int j=0;j<n;j++){
     p1[p[j]]=j;
   }
   return p1;
 }
  gen _perminv(const gen & args,GIAC_CONTEXT){
    if  (args.type!=_VECT) 
      setsizeerr();
    vecteur v(*args._VECTptr); 
    vector<int> p;
    if (!is_permu(v,p,contextptr))
      setsizeerr();
    return vector_int_2_vecteur(perminv(p),contextptr);
  }
  const string _perminv_s("perminv");
  unary_function_eval __perminv(&_perminv,_perminv_s);
  unary_function_ptr at_perminv (&__perminv,0,true);
  
 vector<int> cycleinv(const vector<int> & c){
   int n;
   n=c.size();
   vector<int> c1(n);
   for (int j=0;j<n;j++){
     c1[j]=c[n-j-1];
   }
   return c1;
 }
  gen _cycleinv(const gen & args,GIAC_CONTEXT){
    if  (args.type!=_VECT) 
      setsizeerr();
    vecteur v(*args._VECTptr);
    vector<int> c;
    if (!is_cycle(v,c,contextptr))
      setsizeerr();
    return vector_int_2_vecteur(cycleinv(c),contextptr);
  }
  const string _cycleinv_s("cycleinv");
  unary_function_eval __cycleinv(&_cycleinv,_cycleinv_s);
  unary_function_ptr at_cycleinv (&__cycleinv,0,true);
  

  int signature(const vector<int> & p) {
    //renvoie la signature de la permutation p
    int s; 
    s=1;     
    vector< vector<int> > c;
    c=permu2cycles(p);
    int l=c.size();
    for (int k=0;k<l;k++){
      int lk=c[k].size()-1;
      if (lk%2) s=s*-1;
    }
    return(s);
  }

  gen _signature(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT) 
      setsizeerr();
    vecteur v(*args._VECTptr);
    return signature(vecteur_2_vector_int(v,contextptr));
  }

  const string _signature_s("signature");
  unary_function_eval __signature(&_signature,_signature_s);
  unary_function_ptr at_signature (&__signature,0,true);

  vecteur vector_double_2_vecteur(const vector<double> & v){
    //transforme un vector<double> en vecteur 
    vector<double>::const_iterator it=v.begin(),itend=v.end();
    vecteur res;
    res.reserve(itend-it);
    for (;it!=itend;++it)
      res.push_back(*it);
    return res;
  } 

  vecteur vectvector_double_2_vecteur(const vector< vector<double> > & v){
    //transforme un vector< vector<double> > en vecteur  
    int s=v.size();
    vecteur res;
    res.reserve(s);
    for (int i=0;i<s;++i)
      res.push_back(vector_double_2_vecteur(v[i]));
    return res;
  }  
  
  gen _hilbert(const gen & args){
    int n,p;
    if (args.type==_INT_) {
      n=args.val;
      p=args.val;
    }
    else {
      if ( (args.type!=_VECT)  || (args._VECTptr->size()!=2) )
	settypeerr();
      vecteur v(*args._VECTptr);
      gen v1=v.front(),v2=v.back();
      n=v1.val;
      p=v2.val; 
    }   
    vecteur c;
    for (int k=0;k<n;k++){
      vecteur l(p);
      for (int j=0;j<p;j++){
	l[j]=rdiv(1,k+j+1);
      }
      c.push_back(l);
    } 
    return c;
  }

  const string _hilbert_s("hilbert");
  unary_function_unary __hilbert(&_hilbert,_hilbert_s);
  unary_function_ptr at_hilbert (&__hilbert,0,true);

  gen l2norm2(const gen & g){
    if (g.type!=_VECT)
      return g*g;
    const_iterateur it=g._VECTptr->begin(),itend=g._VECTptr->end();    
    gen res(0);
    for (;it!=itend;++it)
      res = res + (*it)*(*it);
    return res;
  }

  gen square_hadamard_bound(const matrice & m){
    const_iterateur it=m.begin(),itend=m.end();
    gen prod(1);
    for (;it!=itend;++it)
      prod=prod*l2norm2(*it);
    return prod;
  }

  gen _hadamard(const gen & args,GIAC_CONTEXT){
    if (ckmatrix(args) || args[0][0].type!=_VECT){
      // Hadamard bound on det(args)
      matrice & m=*args._VECTptr;
      return sqrt(min(square_hadamard_bound(m),square_hadamard_bound(mtran(m)),contextptr),contextptr);
    }
    if ( (args.type!=_VECT)  || (args._VECTptr->size()!=2) )
      settypeerr();
    vecteur v(*args._VECTptr);
    gen g1=v.front(),g2=v.back();
     
    if ((g1.type!=_VECT) ||(g2.type!=_VECT))
      settypeerr();
    vecteur v1(*g1._VECTptr); 
    vecteur v2(*g2._VECTptr);
    if (v1.size()!=v2.size()) setsizeerr();
    int n=v1.size();
    vecteur c;
    for (int k=0;k<n;k++){
      if ((v1[k].type!=_VECT) ||(v2[k].type!=_VECT)) settypeerr(); 
      vecteur l1(*(v1[k])._VECTptr);
      vecteur l2(*(v2[k])._VECTptr);
      if (l1.size()!=l2.size()) setsizeerr();
      int p=l1.size();
      vecteur l(p);
      for (int j=0;j<p;j++){
	l[j]=l1[j]*l2[j];
      }
      c.push_back(l);
    } 
    return c;
  }
  const string _hadamard_s("hadamard");
  unary_function_eval __hadamard(&_hadamard,_hadamard_s);
  unary_function_ptr at_hadamard (&__hadamard,0,true);

  gen _trn(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      settypeerr();
    vecteur v(*args._VECTptr);
    int n=v.size();
    vecteur c;
    for (int k=0;k<n;k++){
      vecteur lc(n);
      if (v[k].type!=_VECT) settypeerr(); 
      vecteur l(*(v[k])._VECTptr);
      if (l.size()!=unsigned(n)) setsizeerr();
      for (int j=0;j<n;j++){
	lc[j]=conj(l[j],contextptr);
      }

      c.push_back(lc);
    } 
    c=mtran(c);
    return c;
  }

  const string _trn_s("trn");
  unary_function_eval __trn(&_trn,_trn_s);
  unary_function_ptr at_trn (&__trn,0,true);

  gen _syst2mat(const gen & args,GIAC_CONTEXT){
    if ( (args.type!=_VECT)  || (args._VECTptr->size()!=2) )
      settypeerr();
    vecteur v(*args._VECTptr);
    gen g1=_equal2diff(v.front()),g2=v.back();
    
    if ((g1.type!=_VECT) ||(g2.type!=_VECT))
      settypeerr();
    vecteur v1(*g1._VECTptr);
    vecteur v2(*g2._VECTptr);
    
    int n=v2.size(),m=v1.size();
    vecteur c;
    for (int k=0;k<m;k++){
      vecteur l(n+1);
      gen ln=v1[k];
      for (int j=0;j<n;j++){
	l[j]=derive(v1[k],v2[j],contextptr);
	ln=subst(ln,v2[j],0,false,contextptr);
      }
      l[n]=ln;
      c.push_back(l);
    } 
    return c;
  }
  const string _syst2mat_s("syst2mat");
  unary_function_eval __syst2mat(&_syst2mat,_syst2mat_s);
  unary_function_ptr at_syst2mat (&__syst2mat,0,true);

  gen _vandermonde(const gen & args){
    if (args.type!=_VECT)  
      settypeerr();
    vecteur v(*args._VECTptr);
    int n=v.size();
    vecteur c; 
    vecteur l(n); 
    for (int j=0;j<n;j++){
      l[j]=1;
    }
    c.push_back(l);
    for (int k=1;k<n;k++){
      for (int j=0;j<n;j++){
	l[j]=l[j]*v[j];
      }
      c.push_back(l);
    }
    return c;
  }
  const string _vandermonde_s("vandermonde");
  unary_function_unary __vandermonde(&_vandermonde,_vandermonde_s);
  unary_function_ptr at_vandermonde (&__vandermonde,0,true);


  gen _laplacian(const gen & args,GIAC_CONTEXT){
    if ( (args.type!=_VECT)  || (args._VECTptr->size()!=2) )
      settypeerr();    
    vecteur v(plotpreprocess(args,contextptr));
    gen g1=v.front(),g2=v.back();
    if (g2.type!=_VECT) 
      settypeerr();
    vecteur v2(*g2._VECTptr);
    int n=v2.size();
    gen la;
    la=0;
    for (int k=0;k<n;k++){
      la=la+derive(derive(g1,v2[k],contextptr),v2[k],contextptr);
    } 
    return normal(la,contextptr);
  }
  const string _laplacian_s("laplacian");
  unary_function_eval __laplacian(&_laplacian,_laplacian_s);
  unary_function_ptr at_laplacian (&__laplacian,_QUOTE_ARGUMENTS,true);

  gen _hessian(const gen & args,GIAC_CONTEXT){
    if ( (args.type!=_VECT)  || (args._VECTptr->size()!=2) )
      settypeerr();        
    vecteur v(plotpreprocess(args,contextptr));
    gen g1=v.front(),g2=v.back();
    if (g2.type!=_VECT) 
      settypeerr();
    vecteur v2(*g2._VECTptr);
    int n=v2.size();
    vecteur he;    
    for (int k=0;k<n;k++){
      vecteur l(n);
      for (int j=0;j<n;j++){
	l[j]=derive(derive(g1,v2[k],contextptr),v2[j],contextptr);
      }
      he.push_back(l);
    }
    return (he);
  }    
  const string _hessian_s("hessian");
  unary_function_eval __hessian(&_hessian,_hessian_s);
  unary_function_ptr at_hessian (&__hessian,_QUOTE_ARGUMENTS,true);

  gen _divergence(const gen & args,GIAC_CONTEXT){
    if ( (args.type!=_VECT)  || (args._VECTptr->size()!=2) )
      settypeerr();    
    vecteur v(plotpreprocess(args,contextptr));
    gen g1=v.front(),g2=v.back();
    if ((g1.type!=_VECT) ||(g2.type!=_VECT))
      settypeerr();
    vecteur v1(*g1._VECTptr);
    vecteur v2(*g2._VECTptr);
    int n=v2.size();
    gen di;
    di=0;
    for (int k=0;k<n;k++){
      di=di+derive(v1[k],v2[k],contextptr);
    } 
    return normal(di,contextptr);
  }    
  const string _divergence_s("divergence");
  unary_function_eval __divergence(&_divergence,_divergence_s);
  unary_function_ptr at_divergence (&__divergence,_QUOTE_ARGUMENTS,true);

  gen _curl(const gen & args,GIAC_CONTEXT){
    if ( (args.type!=_VECT)  || (args._VECTptr->size()!=2) )
      settypeerr();
    vecteur v(plotpreprocess(args,contextptr));
    gen g1=v.front(),g2=v.back();
    if ((g1.type!=_VECT) ||(g2.type!=_VECT))
      settypeerr();
    vecteur v1(*g1._VECTptr);
    vecteur v2(*g2._VECTptr);
    int n=v2.size();
    if (n!=3) setsizeerr();
    vecteur rot(3);
    rot[0]=derive(v1[2],v2[1],contextptr)-derive(v1[1],v2[2],contextptr);
    rot[1]=derive(v1[0],v2[2],contextptr)-derive(v1[2],v2[0],contextptr);
    rot[2]=derive(v1[1],v2[0],contextptr)-derive(v1[0],v2[1],contextptr);
   
    return rot;
  }
  const string _curl_s("curl");
  unary_function_eval __curl(&_curl,_curl_s);
  unary_function_ptr at_curl (&__curl,_QUOTE_ARGUMENTS,true);

  void find_n_x(const gen & args,int & n,gen & x,gen & a){
    if (args.type==_INT_){
      n=args.val;
      x=vx_var;
      a=a__IDNT_e;
      return;
    }
    if (args.type!=_VECT || args._VECTptr->size()<2)
      setsizeerr();
    vecteur & v=*args._VECTptr;
    if (v[0].type==_INT_){
      n=v[0].val;
      x=v[1];
    }
    else
      setsizeerr();
    if (v.size()>2)
      a=v[2];
    else
      a=a__IDNT_e;
  }

  vecteur hermite(int n){
    vecteur v(n+1);
    v[0]=pow(plus_two,n);
    for (int k=2;k<=n;k+=2){
      v[k]=-((n+2-k)*(n+1-k)*v[k-2])/(2*k);
    }
    return v;
  }

  gen _hermite(const gen & args,GIAC_CONTEXT){
    int n;
    gen a,x;
    find_n_x(args,n,x,a);
    return r2e(hermite(n),x,contextptr);
  }
  const string _hermite_s("hermite");
  unary_function_eval __hermite(&_hermite,_hermite_s);
  unary_function_ptr at_hermite (&__hermite,0,true);

  gen _laguerre(const gen & args,GIAC_CONTEXT){
    int n;
    gen a,x;
    find_n_x(args,n,x,a);
    gen p0,p1,p2;
    p0=1;
    p1=1+a-x;
    if (n==0) return p0;
    if (n==1) return p1;
    for (int k=2;k<=n;k++){
      //p2=rdiv(2*k+a-1-x,k)*p1-rdiv(k+a-1,k)*p0;
      p2=(2*k+a-1-x)*p1-(k-1)*(k+a-1)*p0;
      p0=p1;
      p1=p2;  
    } 
    //return normal(p2,contextptr);
    return normal(rdiv(p2,factorial(n)),contextptr);
  }
    
  const string _laguerre_s("laguerre");
  unary_function_eval __laguerre(&_laguerre,_laguerre_s);
  unary_function_ptr at_laguerre (&__laguerre,0,true);

  // Improved one
  vecteur tchebyshev1(int n){
    if (n==0) return vecteur(1,1);
    vecteur v(n+1);
    v[0]=pow(gen(2),n-1);
    if (n==1) return v;
    for (int k=2;k<=n;k+=2){
      v[k]=-((n-k+2)*(n-k+1)*v[k-2])/(2*k*(n-k/2));
    }
    return v;
  }
  gen _tchebyshev1(const gen & args,GIAC_CONTEXT){
    int n;
    gen a,x;
    find_n_x(args,n,x,a);
    return r2e(tchebyshev1(n),x,contextptr);
  }
  const string _tchebyshev1_s("tchebyshev1");
  unary_function_eval __tchebyshev1(&_tchebyshev1,_tchebyshev1_s);
  unary_function_ptr at_tchebyshev1 (&__tchebyshev1,0,true);

  // Improved one
  vecteur tchebyshev2(int n){
    vecteur v(n+1);
    v[0]=pow(gen(2),n);
    for (int k=1;k<=n/2;++k){
      v[2*k]=-(n+2-2*k)*(n+1-2*k)*v[2*k-2]/(4*k*(n+1-k));
    }
    return v;
  }
  gen _tchebyshev2(const gen & args,GIAC_CONTEXT){
    gen p0,p1,p2;
    int n;
    gen a,x;
    find_n_x(args,n,x,a);
    return r2e(tchebyshev2(n),x,contextptr);
  }
  const string _tchebyshev2_s("tchebyshev2");
  unary_function_eval __tchebyshev2(&_tchebyshev2,_tchebyshev2_s);
  unary_function_ptr at_tchebyshev2 (&__tchebyshev2,0,true);

  // Legendre is the vecteur returned divided by n!
  vecteur legendre(int n){
    vecteur v0,v1,vtmp1,vtmp2;
    v0.push_back(1);
    v1.push_back(1);
    v1.push_back(0);
    if (!n) return v0;
    if (n==1) return v1;
    for (int k=2;k<=n;k++){
      multvecteur(2*k-1,v1,vtmp1);
      vtmp1.push_back(0); // (2k-1)*x*p1
      multvecteur((k-1)*(k-1),v0,vtmp2); // (k-1)^2*p0
      vtmp1=vtmp1-vtmp2; // p2=(2*k-1)*x*p1-(k-1)*(k-1)*p0;
      v0=v1;
      v1=vtmp1;
    } 
    return v1; 
  }

  gen _legendre(const gen & args,GIAC_CONTEXT){ 
    int n;
    gen a,x;
    find_n_x(args,n,x,a);
    vecteur v=multvecteur(inv(factorial(n),contextptr),legendre(n));
    return r2e(v,x,contextptr);
  }
  const string _legendre_s("legendre");
  unary_function_eval __legendre(&_legendre,_legendre_s);
  unary_function_ptr at_legendre (&__legendre,0,true);

  // arithmetic mean column by column
  vecteur mean(const matrice & m,bool column){
    matrice mt;
    if (column)
      mt=mtran(m);
    else
      mt=m;
    vecteur res;
    const_iterateur it=mt.begin(),itend=mt.end();
    for (;it!=itend;++it){
      const gen & g =*it;
      if (g.type!=_VECT){
	res.push_back(g);
	continue;
      }
      vecteur & v=*g._VECTptr;
      if (v.empty()){
	res.push_back(undef);
	continue;
      }
      const_iterateur jt=v.begin(),jtend=v.end();
      int s=jtend-jt;
      gen somme(0);
      for (;jt!=jtend;++jt){
	//somme = somme + evalf(*jt);
	somme = somme + *jt;
      }
      res.push_back(rdiv(somme,s));
    }
    return res;
  }

  vecteur stddev(const matrice & m,bool column,int variance){
    matrice mt;
    if (column)
      mt=mtran(m);
    else
      mt=m; 
    vecteur moyenne(mean(mt,false));
    vecteur res;
    const_iterateur it=mt.begin(),itend=mt.end();
    for (int i=0;it!=itend;++it,++i){
      const gen & g =*it;
      if (g.type!=_VECT){
	res.push_back(0);
	continue;
      }
      vecteur & v=*g._VECTptr;
      if (v.empty()){
	res.push_back(undef);
	continue;
      }
      const_iterateur jt=v.begin(),jtend=v.end();
      int s=jtend-jt;
      gen somme(0);
      for (;jt!=jtend;++jt){
	// somme = somme + evalf((*jt)*(*jt));
	somme = somme + (*jt)*(*jt);
      }
      if (variance!=3)
	res.push_back(sqrt(rdiv(somme-s*moyenne[i]*moyenne[i],s-(variance==2)),context0));
      else
	res.push_back(rdiv(somme,s)-moyenne[i]*moyenne[i]);
    }
    return res;
  }

  matrice ascsort(const matrice & m,bool column){
    matrice mt;
    if (column)
      mt=mtran(m);
    else
      mt=m;
    iterateur it=mt.begin(),itend=mt.end();
    gen tmp;
    for (;it!=itend;++it){
      gen & g =*it;
      if (g.type!=_VECT)
	continue;
      vecteur v=*g._VECTptr;
      if (v.empty())
	continue;
      const_iterateur jt=v.begin(),jtend=v.end();
      int n=jtend-jt;
      vector<double> vv(n);
      for (int j=0;jt!=jtend;++jt,++j){
	if ( (jt->type==_VECT) && (jt->_VECTptr->size()==3) )
	  tmp=(*jt->_VECTptr)[1];
	else
	  tmp=*jt;
	tmp=evalf(tmp,1,0);
	if (tmp.type!=_DOUBLE_)
	  vv[j]=0;
	else
	  vv[j]=tmp._DOUBLE_val;
      }
      sort(vv.begin(),vv.end());
      for (int j=0;j<n;++j)
	v[j]=vv[j];
      *it=v;
    }
    return mt;
  }

 gen _permu2mat(const gen & args,GIAC_CONTEXT){
   //transforme une permutation en une matrice obtenue en permutant les lignes de la matrice identite  
    if (args.type!=_VECT)  
      settypeerr(); 
    vector<int> p1;
    vecteur p(*args._VECTptr);
    if (!(is_permu(p,p1,contextptr)))  
      settypeerr();
    int n=p.size();
    vecteur c; 
    vecteur l(n); 
    for (int k=0;k<n;k++){
      for (int j=0;j<n;j++){
	if (p[k]==j+(xcas_mode(contextptr)!=0)) {
	  l[j]=1;
	  } else {
	    l[j]=0;
	  }
      }
      c.push_back(l);
    }
    return c;
  }

  const string _permu2mat_s("permu2mat");
  unary_function_eval __permu2mat(&_permu2mat,_permu2mat_s);
  unary_function_ptr at_permu2mat (&__permu2mat,0,true);

  bool est_dans(const vector<int> & a , const int n, vector< vector<int> > s) {
    //teste si a est egal a l'un des n premiers elements de s 
    bool cont=true;
    int j=0;
    while (j<=n && cont) {
      if (a==s[j]) {
	cont=false;
      }
      j=j+1;
    }
    return (! cont);
  } 

 vector< vector<int> > groupermu(const vector<int> & p1,const vector<int> & p2) {
    //groupe engendre par de 2 perm p1 et p2 de long n1 et n2 : a pour long max(n1,n2)
    int n1;
    n1=p1.size();
    int n2;
    n2=p2.size();
    vector<int> a;
    vector<int> b;
    a=p1;
    b=p2;
    if (n1>n2) {
      for (int k=n2;k<n1;k++) {b.push_back(k);n2=n1;} 
    } else {
      for (int k=n1;k<n2;k++) {a.push_back(k);n1=n2;}
    }
    ;
    vector< vector<int> > s(2);
    s[0]=a;
    s[1]=b;
    int p=0;
    int q=1;
    bool pfini=true;
    while (pfini) {
      int k=0;
      for (int j=p;j<=q;j++){
	vector<int> na;
	vector<int> nb;
	na=p1op2(a,s[j]);
	if (!(est_dans(na,q+k,s))){
	  k=k+1;
	  s.push_back(na);
	}
	nb=p1op2(b,s[j]);
	if (!(est_dans(nb,q+k,s))){
	  k=k+1;
	  s.push_back(nb);
	}
      }
      if (k!=0) {
	p=q+1;
	q=q+k;
      } else {
	pfini=false;
      }
    } 
    //q est l'ordre du groupe = size(s)     
    return(s);
 }
 
  gen _groupermu(const gen & args,GIAC_CONTEXT){
    if ( (args.type!=_VECT)  || (args._VECTptr->size()!=2) )
      settypeerr();
    vecteur v(*args._VECTptr);
    gen v1=v.front(),v2=v.back();
    if ( (v1.type!=_VECT) || (v2.type!=_VECT))
      settypeerr();
    vector<int> p1,p2;   
    if (!is_permu(*v1._VECTptr,p1,contextptr) || !is_permu(*v2._VECTptr,p2,contextptr))
      setsizeerr();
    return vectvector_int_2_vecteur(groupermu(p1,p2),contextptr);
  }
  const string _groupermu_s("groupermu");
  unary_function_eval __groupermu(&_groupermu,_groupermu_s);
  unary_function_ptr at_groupermu (&__groupermu,0,true);

  gen _nextperm(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      settypeerr();
    const vecteur & v1=*args._VECTptr;
    vector<int> p1;
    if (!is_permu(v1,p1,contextptr))
      setsizeerr();
    if (next_permutation(p1.begin(),p1.end()))
      return vector_int_2_vecteur(p1,contextptr);
    else
      return undef;
  }
  const string _nextperm_s("nextperm");
  unary_function_eval __nextperm(&_nextperm,_nextperm_s);
  unary_function_ptr at_nextperm (&__nextperm,0,true);

  gen _prevperm(const gen & args,GIAC_CONTEXT){
    if (args.type!=_VECT)
      settypeerr();
    const vecteur & v1=*args._VECTptr;
    vector<int> p1;
    if (!is_permu(v1,p1,contextptr))
      setsizeerr();
    if (prev_permutation(p1.begin(),p1.end()))
      return vector_int_2_vecteur(p1,contextptr);
    else
      return undef;
  }
  const string _prevperm_s("prevperm");
  unary_function_eval __prevperm(&_prevperm,_prevperm_s);
  unary_function_ptr at_prevperm (&__prevperm,0,true);

  gen _split(const gen & args,GIAC_CONTEXT){
    //renvoie [ax,ay] si les arg sont g1=ax*an (sans denominateur) et g2=[x,y]
    //sinon renvoie [0]
    if ( (args.type!=_VECT)  || (args._VECTptr->size()!=2) )
      settypeerr();
    vecteur v(*args._VECTptr);
    gen g1=v.front(),g2=v.back();
    if (g2.type!=_VECT)
      settypeerr();
    vecteur v2(*g2._VECTptr);
    int n=v2.size();
    if (n!=2) setsizeerr();
    vecteur fa;
    fa=factors(g1,vx_var,contextptr);
    int l=fa.size();
    gen ax=1;
    gen ay=1;
    for (int k=0;k<l;k=k+2){
      gen f=fa[k];
      if (derive(f,v2[0],contextptr)==0) {
        ay=ay*pow(f,fa[k+1],contextptr);
      }
      else {
	if (derive(f,v2[1],contextptr)==0){
	  ax=ax*pow(f,fa[k+1],contextptr);
	}
	else {vecteur res(1);return (res);}
      }
    }
    vecteur res(2);
    res[0]=ax;
    res[1]=ay;
    return res;
  }
 
  const string _split_s("split");
  unary_function_eval __split(&_split,_split_s);
  unary_function_ptr at_split (&__split,0,true);
 
 gen _sum_riemann(const gen & args,GIAC_CONTEXT){
  //renvoie l'equivalent pour n=infini de sum de k=1 a n de g1 avec g2=[n,k])  
    if ( (args.type!=_VECT)  || (args._VECTptr->size()!=2) )
      settypeerr();
    gen g1=args[0],g2=args[1];
    if (g2.type!=_VECT)
      settypeerr();
    vecteur v2(*g2._VECTptr);
    int sv2=v2.size();
    if (sv2!=2) setsizeerr();
    identificateur x("_x");
    //on pose k=n*x
    //mettre que v2[0]=n et v2[1]=k sont ds N
    //_assume(makevecteur(is_plus(v2[0]))); 
    //_assume(makevecteur(is_plus(v2[1])));
    gen a=recursive_normal(v2[0]*subst(g1,v2[1],x*v2[0],false,contextptr),contextptr);
    gen nda=_fxnd(a);
    vecteur nada(*nda._VECTptr);
    //na numerateur de a et da denominateur de a ou k=n*x
    gen na=nada[0];
    gen da=nada[1];
    vecteur var(2);
    var[0]=na;
    v2[1]=x;
    var[1]=v2;
    //var=[na,[n,x]]
    gen b0=_split(var,contextptr);
    //on separe les variables du numerateur
    vecteur vb0(*b0._VECTptr);
    if ((vb0.size()==1)&& (vb0[0]==0))
      return string2gen("ce n'est probablement pas une somme de riemann",false);
    var[0]=da;
    //var=[da,[n,x]]
    gen b1=_split(var,contextptr);
    //on separe les variables du denominateur
    vecteur vb1(*b1._VECTptr);
    if ((vb1.size()==1)&& (vb1[0]==0)) 
      return string2gen("ce n'est probablement pas une somme de riemann",false);
    gen an=vb0[0]/vb1[0];
    gen ax=vb0[1]/vb1[1];
    gen tmp=_integrate(makevecteur(ax,x,0,1),contextptr);
    gen tmp2=_limit(makevecteur(an,v2[0],plus_inf),contextptr);
    //tmp n'a pas ete calcule sum_riemann(pi/(2*n)*log(sin(pi*k/(2*n))),[n,k])
    //if ((tmp.type==_SYMB)&&(tmp._SYMBptr->sommet==at_integrate))
    //return _limit(makevecteur(tmp*an,v2[0],plus_inf));
    //tmp ne doit pas etre infini qd tmp2 est nul et reciproquement
    if (is_inf(tmp)&& is_zero(tmp2))
      return string2gen("ce n'est probablement pas une somme de riemann",false);
    if (is_zero(tmp)&& is_inf(tmp2))
      return string2gen("ce n'est probablement pas une somme de riemann",false);
    return recursive_normal(tmp*_series(makevecteur(an,v2[0],plus_inf),contextptr),contextptr);
  }
  const string _sum_riemann_s("sum_riemann");
  unary_function_eval __sum_riemann(&_sum_riemann,_sum_riemann_s);
  unary_function_ptr at_sum_riemann(&__sum_riemann,0,true);

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC
