// -*- mode:C++ ; compile-command: "g++ -I.. -g -c usual.cc" -*-
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
#ifndef _GIAC_USUAL_H
#define _GIAC_USUAL_H
#include "first.h"
#include <string>

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  // Global vectors of tractable functions and equivalents
  extern std::vector<unary_function_ptr *> limit_tractable_functions;
  typedef gen ( * gen_op ) (const gen & arg);
  extern std::vector<gen_op_context> limit_tractable_replace;

  extern std::string messages_to_print ;
  extern bool signal_store; // if false then child sto is not signal to parent
  // declare here pointers to rational operators 
  // These are global variables but that's normal for functional obj

  class gen;
  class unary_function_ptr;
  class unary_function_eval;
  class partial_derivative_onearg;

  // convert a gen to a string, format=0 (normal), 1 (tex)
  std::string gen2string(const gen & g,int format,GIAC_CONTEXT);

  // Demodularize
  gen unmod(const gen & g);
  gen unmodunprod(const gen & g);

  gen double2gen(double d);
  bool is_equal(const gen & g);
  gen apply_to_equal(const gen & g,const gen_op & f);
  gen apply_to_equal(const gen & g,gen (* f) (const gen &, GIAC_CONTEXT),GIAC_CONTEXT);

  std::string print_with_parenthesis_if_required(const gen & g,int format,GIAC_CONTEXT);

  // usual zero args
  extern unary_function_ptr at_zero;
  extern unary_function_ptr at_one;
  extern unary_function_ptr at_id;
  extern const std::string _rm_a_z_s;
  gen _rm_a_z(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_rm_a_z;
  gen _rm_all_vars(const gen & args,const context * contextptr);
  extern unary_function_ptr at_rm_all_vars;
  extern unary_function_ptr reim_op[];

  // usual unary function related declarations (extended unary)
  gen _neg(const gen & args);
  extern unary_function_ptr at_neg ;
  symbolic symb_inv(const gen & a);
  gen _inv(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_inv ;
  symbolic symb_not(const gen & args);
  gen _not(const gen & args);
  extern unary_function_ptr at_not;

  symbolic symb_exp(const gen & e);
  gen exp(const gen & e,GIAC_CONTEXT);
  extern unary_function_ptr at_exp ;

  symbolic symb_ln(const gen & e);
  gen log(const gen & e,GIAC_CONTEXT);
  gen ln(const gen & e,GIAC_CONTEXT);
  extern unary_function_ptr at_ln ;

  symbolic symb_log10(const gen & e);
  gen log10(const gen & e,GIAC_CONTEXT);
  extern unary_function_ptr at_log10 ;

  symbolic symb_alog10(const gen & e,GIAC_CONTEXT);
  gen alog10(const gen & e,GIAC_CONTEXT);
  extern unary_function_ptr at_alog10 ;

  symbolic symb_atan(const gen & e);
  gen atan(const gen & e,GIAC_CONTEXT);
  extern unary_function_ptr at_atan ;

  // e = +/- simpl*doubl^2
  void zint2simpldoublpos(const gen & e,gen & simpl,gen & doubl,bool & pos);
  symbolic symb_sqrt(const gen & e);
  gen sqrt_noabs(const gen & e,GIAC_CONTEXT);
  gen sqrt(const gen & e,GIAC_CONTEXT);
  extern unary_function_ptr at_sqrt ;

  gen sq(const gen & e);
  extern unary_function_ptr at_sq ;

  symbolic symb_sin(const gen & e);
  gen sin(const gen & e,GIAC_CONTEXT);
  extern unary_function_ptr at_sin ;

  symbolic symb_cos(const gen & e);
  gen cos(const gen & e,GIAC_CONTEXT);
  extern unary_function_ptr at_cos ;

  symbolic symb_tan(const gen & e);
  gen tan(const gen & e,GIAC_CONTEXT);
  extern unary_function_ptr at_tan ;

  symbolic symb_asin(const gen & e);
  gen asin(const gen & e,GIAC_CONTEXT);
  extern unary_function_ptr at_asin ;

  symbolic symb_acos(const gen & e);
  gen acos(const gen & e,GIAC_CONTEXT);
  extern unary_function_ptr at_acos ;

  symbolic symb_sinh(const gen & e);
  gen sinh(const gen & e,GIAC_CONTEXT);
  extern unary_function_ptr at_sinh ;

  symbolic symb_cosh(const gen & e);
  gen cosh(const gen & e,GIAC_CONTEXT);
  extern unary_function_ptr at_cosh ;

  symbolic symb_tanh(const gen & e);
  gen tanh(const gen & e,GIAC_CONTEXT);
  extern unary_function_ptr at_tanh ;

  symbolic symb_asinh(const gen & e);
  gen asinh(const gen & e,GIAC_CONTEXT);
  extern unary_function_ptr at_asinh ;

  symbolic symb_acosh(const gen & e);
  gen acosh(const gen & e,GIAC_CONTEXT);
  extern unary_function_ptr at_acosh ;

  symbolic symb_atanh(const gen & e);
  gen atanh(const gen & e,GIAC_CONTEXT);
  extern unary_function_ptr at_atanh ;

  symbolic symb_quote(const gen & arg);
  gen quote(const gen & e,GIAC_CONTEXT);
  extern unary_function_ptr at_quote ;

  symbolic symb_unquote(const gen & arg);
  gen unquote(const gen & e,GIAC_CONTEXT);
  extern unary_function_ptr at_unquote ;

  symbolic symb_order_size(const gen & e);
  gen order_size(const gen & e);
  extern unary_function_ptr at_order_size ;

  symbolic symb_and(const gen & a,const gen & b);
  gen _and(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_and;

  symbolic symb_ou(const gen & a,const gen & b);
  gen _ou(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_ou;

  gen _xor(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_xor;

  symbolic symb_min(const gen & a,const gen & b);
  gen _min(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_min;

  symbolic symb_max(const gen & a,const gen & b);
  gen _max(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_max;

  symbolic symb_gcd(const gen & a,const gen & b);
  gen _gcd(const gen & args);
  extern unary_function_ptr at_gcd;

  symbolic symb_lcm(const gen & a,const gen & b);
  gen _lcm(const gen & args);
  extern unary_function_ptr at_lcm;

  symbolic symb_egcd(const gen & a,const gen & b);
  gen _egcd(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_egcd;

  symbolic symb_iegcd(const gen & a,const gen & b);
  gen _iegcd(const gen & args);
  extern unary_function_ptr at_iegcd;

  symbolic symb_iquo(const gen & a,const gen & b);
  gen Iquo(const gen & a,const gen & b);
  gen _iquo(const gen & args);
  extern unary_function_ptr at_iquo;

  symbolic symb_irem(const gen & a,const gen & b);
  gen _irem(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_irem;

  gen _mods(const gen & args);
  extern unary_function_ptr at_mods;

  gen _quote_pow(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_quote_pow;

  symbolic symb_iquorem(const gen & a,const gen & b);
  gen _iquorem(const gen & args);
  extern unary_function_ptr at_iquorem;

  symbolic symb_smod(const gen & a,const gen & b);
  gen _smod(const gen & args);
  extern unary_function_ptr at_smod;

  symbolic symb_rdiv(const gen & a,const gen & b);
  gen _rdiv(const gen & args);
  extern unary_function_ptr at_rdiv;

  symbolic symb_is_prime(const gen & a);
  gen _is_prime(const gen & args);
  extern unary_function_ptr at_is_prime;

  gen _nextprime(const gen & args);
  extern unary_function_ptr at_nextprime;

  gen _prevprime(const gen & args);
  extern unary_function_ptr at_prevprime;

  symbolic symb_floor(const gen & a);
  gen _floor(const gen & args,GIAC_CONTEXT);  
  extern const std::string _floor_s;
  extern unary_function_ptr at_floor;

  symbolic symb_ceil(const gen & a);
  gen _ceil(const gen & args,GIAC_CONTEXT);  
  extern const std::string _ceil_s;
  extern unary_function_ptr at_ceil;
  gen ceil2floor(const gen & g,GIAC_CONTEXT);

  symbolic symb_round(const gen & a);
  gen _round(const gen & args,GIAC_CONTEXT);  
  extern const std::string _round_s;
  extern unary_function_ptr at_round;

  symbolic symb_print(const gen & a);
  gen _print(const gen & args,GIAC_CONTEXT);  
  extern const std::string _print_s;
  extern unary_function_ptr at_print;
  extern unary_function_eval __print;

  symbolic symb_jacobi_symbol(const gen & a,const gen & b);
  gen _jacobi_symbol(const gen & args);
  extern unary_function_ptr at_jacobi_symbol;

  symbolic symb_legendre_symbol(const gen & a,const gen & b);
  gen _legendre_symbol(const gen & args);
  extern unary_function_ptr at_legendre_symbol;

  symbolic symb_ichinrem(const gen & a,const gen & b);
  gen _ichinrem(const gen & args);
  extern unary_function_ptr at_ichinrem;

  gen _fracmod(const gen & args);
  extern unary_function_ptr at_fracmod;

  gen _factorial(const gen & args);
  extern unary_function_ptr at_factorial;

  gen _perm(const gen & args);
  extern unary_function_ptr at_perm;

  gen comb(const gen & n,const gen &k,GIAC_CONTEXT);
  gen _comb(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_comb;

  symbolic symb_chinrem(const gen & a,const gen & b);
  gen _chinrem(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_chinrem;

  extern unary_function_ptr at_re;

  extern unary_function_ptr at_im;

  symbolic symb_conj(const gen & e);
  extern unary_function_ptr at_conj ;

  symbolic symb_abs(const gen & e);
  extern unary_function_ptr at_abs ;

  symbolic symb_arg(const gen & e);
  extern unary_function_ptr at_arg ;

  symbolic symb_sign(const gen & e);
  extern unary_function_ptr at_sign;

  symbolic symb_cyclotomic(const gen & e);
  extern unary_function_ptr at_cyclotomic ;

  symbolic symb_quo(const gen & a,const gen & b);
  extern unary_function_ptr at_quo ;

  symbolic symb_rem(const gen & a,const gen & b);
  extern unary_function_ptr at_rem ;

  symbolic symb_quorem(const gen & a,const gen & b);
  gen _quorem(const gen & args,GIAC_CONTEXT);
  gen _quo(const gen & args,GIAC_CONTEXT);
  gen _rem(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_quorem ;
  extern unary_function_ptr at_normalmod;

  symbolic symb_sto(const gen & a,gen & b,bool in_place=false);
  symbolic symb_sto(const gen & e);
  extern unary_function_ptr at_sto ;
  extern unary_function_ptr at_array_sto ;
  extern unary_function_ptr at_increment;
  extern unary_function_ptr at_decrement;
  extern unary_function_ptr at_multcrement;
  extern unary_function_ptr at_divcrement;
  gen sto(const gen & a,const gen & b,GIAC_CONTEXT);
  gen sto(const gen & a,const gen & b,bool in_place,GIAC_CONTEXT);  
  gen _sto(const gen & g,const context * contextptr);

  bool is_assumed_integer(const gen & g,GIAC_CONTEXT);

  // assume format _VECT of subtype _ASSUME__VECT
  // [ _DOUBLE_ , list of intervals, excluded ]
  // [ _INT_, ]
  extern unary_function_ptr at_assume ;
  gen giac_assume(const gen & a,GIAC_CONTEXT);
  // returns the assumed idnt name
  // used if assumptions are in OR conjonction
  gen assumesymbolic(const gen & a,gen idnt_must_be,GIAC_CONTEXT);
  // v = previous assumptions, a=the real value, direction
  // is positive for [a,+inf[, negative for ]-inf,a]
  // |direction| = 1 (large) or 2 (strict) 
  gen doubleassume_and(const vecteur & v,const gen & a,int direction,bool or_assumption,GIAC_CONTEXT);

  gen _equal(const gen & args);
  gen symb_equal(const gen & a,const gen & b);
  extern unary_function_ptr at_equal;
  gen symb_same(const gen & a);
  symbolic symb_same(const gen & a,const gen & b);
  gen _same(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_same;

  symbolic symb_inferieur_strict(const gen & a);
  symbolic symb_inferieur_strict(const gen & a,const gen & b);
  gen _inferieur_strict(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_inferieur_strict;

  symbolic symb_inferieur_egal(const gen & a);
  symbolic symb_inferieur_egal(const gen & a,const gen & b);
  gen _inferieur_egal(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_inferieur_egal;

  symbolic symb_superieur_strict(const gen & a);
  symbolic symb_superieur_strict(const gen & a,const gen & b);
  gen _superieur_strict(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_superieur_strict;

  symbolic symb_superieur_egal(const gen & a);
  symbolic symb_superieur_egal(const gen & a,const gen & b);
  gen _superieur_egal(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_superieur_egal;

  int is_inequality(const gen & g);
  extern std::vector<unary_function_ptr> inequality_sommets;

  symbolic symb_different(const gen & a);
  symbolic symb_different(const gen & a,const gen & b);
  gen _different(const gen & args);
  extern unary_function_ptr at_different;

  symbolic symb_of(const gen & a);
  symbolic symb_of(const gen & a,const gen & b);
  gen _of(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_of;

  symbolic symb_at(const gen & a);
  symbolic symb_at(const gen & a,const gen & b,GIAC_CONTEXT);
  gen _at(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_at;
  
  gen _table(const gen & args);
  extern unary_function_ptr at_table;
  
  // usual multiargs
  // for multiargs we use _name for the corresponding "unary" function
  // to avoid confusion of pointers since name is used
  symbolic symb_plus(const gen & a,const gen & b);
  gen _plus(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_plus ;

  symbolic symb_prod(const gen & a,const gen & b);
  gen _prod(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_prod ;

  symbolic symb_pow(const gen & a,const gen & b);
  gen _pow(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_pow ;
  extern const std::string _pow_s;

  symbolic symb_powmod(const gen & a,const gen & b,const gen &c);
  symbolic symb_powmod(const gen & a);
  gen _powmod(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_powmod ;
  extern const std::string _powmod_s;

  symbolic symb_tran(const gen & a);
  symbolic symb_trace(const gen & a);
  symbolic symb_rref(const gen & a);
  symbolic symb_idn(const gen & e);
  symbolic symb_ranm(const gen & e);
  symbolic symb_det(const gen & a);
  symbolic symb_pcar(const gen & a);
  symbolic symb_ker(const gen & a);
  symbolic symb_image(const gen & a);
  symbolic symb_moyal(const gen & a,const gen & b, const gen &vars,const gen & order);

  symbolic symb_evalf(const gen & e);
  extern unary_function_ptr at_evalf;
  gen _evalf(const gen &,GIAC_CONTEXT);

  symbolic symb_eval(const gen & e);
  extern unary_function_ptr at_eval;
  extern unary_function_ptr at_evalm;
  
  symbolic symb_subst(const gen & e);
  gen _subst(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_subst;
  std::string printassubs(const gen & feuille,const std::string & sommetstr,GIAC_CONTEXT);
  std::string printasmaple_subs(const gen & feuille,const std::string & sommetstr,GIAC_CONTEXT);
  gen _subs(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_subs;

  gen _maple_subs(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_maple_subs;
  extern unary_function_ptr at_ampersand_times;

  extern unary_function_ptr solve_fcns[];
  extern std::vector<unary_function_ptr> solve_fcns_v;

  symbolic symb_version(const gen & e);
  extern unary_function_ptr at_version;
  std::string version();
  gen _version(const gen & a);
  
  gen Gamma(const gen & x,GIAC_CONTEXT);
  extern const std::string _Gamma_s;
  gen _Gamma(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_Gamma;

  gen Psi(const gen & x,GIAC_CONTEXT);
  gen Psi(const gen & x,int n,GIAC_CONTEXT);
  extern const std::string _Psi_s;
  gen _Psi(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_Psi;

  gen Zeta(const gen & x,GIAC_CONTEXT);
  extern const std::string _Zeta_s;
  gen _Zeta(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_Zeta;
  extern unary_function_ptr at_Eta;

  gen _erf(const gen & args,GIAC_CONTEXT);
  gen erf(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_erf ;

  gen _erfc(const gen & args,GIAC_CONTEXT);
  gen erfc(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_erfc ;

  extern unary_function_ptr at_Si;
  extern unary_function_ptr at_Ci;
  extern unary_function_ptr at_Ci0;
  gen Ci_replace0(const gen & g,GIAC_CONTEXT);
  gen _Ci(const gen & g,GIAC_CONTEXT);
  gen _Si(const gen & g,GIAC_CONTEXT);
  gen _Ei(const gen & g,GIAC_CONTEXT);
  extern unary_function_ptr at_Ei;
  extern unary_function_ptr at_Ei0;
  gen Ei_replace0(const gen & g,GIAC_CONTEXT);
  extern unary_function_ptr at_Heaviside ;
  extern unary_function_ptr at_Dirac ;
  gen _Heaviside(const gen & args,GIAC_CONTEXT);
  gen _Dirac(const gen & args,GIAC_CONTEXT);

  extern gen plus_one_half;
  extern gen minus_one_half;
  extern gen minus_sqrt3_2;
  extern gen minus_sqrt2_2;
  extern gen plus_sqrt2_2;
  extern gen plus_sqrt3_2;
  extern gen minus_sqrt3_3;
  extern gen minus_sqrt3;
  extern gen minus_sqrt6;
  extern gen plus_sqrt3;
  extern gen plus_sqrt2;
  extern gen plus_sqrt6;
  extern gen plus_sqrt3_3;
  extern gen tan_pi_12;
  extern gen minus_tan_pi_12;
  extern gen tan_5pi_12;
  extern gen cos_pi_12;
  extern gen sin_pi_12;
  extern gen minus_cos_pi_12;
  extern gen minus_sin_pi_12;
  extern gen minus_tan_5pi_12;
  extern gen cst_two_pi;
  extern gen cst_pi_over_2;
  extern gen plus_inf;
  extern gen minus_inf;
  extern gen rad2deg_e;
  extern gen deg2rad_e;

  // for subst.cc
  extern unary_function_ptr sincostan_tab[];
  extern unary_function_ptr asinacosatan_tab[];
  extern unary_function_ptr sinhcoshtanh_tab[];
  extern std::vector<unary_function_ptr> sincostan_v;
  extern std::vector<unary_function_ptr> asinacosatan_v,sinhcoshtanh_v,sincostansinhcoshtanh_v;
  extern std::vector<unary_function_ptr> exp_v;
  extern std::vector<unary_function_ptr> tan_v;
  extern std::vector<unary_function_ptr> asin_v;
  extern std::vector<unary_function_ptr> acos_v;
  extern std::vector<unary_function_ptr> atan_v;
  extern std::vector<unary_function_ptr> pow_v;
  extern gen_op_context halftan_tab[];
  extern gen_op_context hyp2exp_tab[];
  extern gen_op_context trig2exp_tab[];
  extern gen_op_context atrig2ln_tab[];
  extern std::vector< gen_op_context > halftan_v;
  extern std::vector< gen_op_context > hyp2exp_v;
  extern std::vector< gen_op_context > trig2exp_v;
  extern std::vector< gen_op_context > halftan_hyp2exp_v;
  extern std::vector< gen_op_context > exp2sincos_v;
  extern std::vector< gen_op_context > tan2sincos_v;
  extern std::vector< gen_op_context > tan2sincos2_v;
  extern std::vector< gen_op_context > tan2cossin2_v;
  extern std::vector< gen_op_context > asin2acos_v;
  extern std::vector< gen_op_context > asin2atan_v;
  extern std::vector< gen_op_context > acos2asin_v;
  extern std::vector< gen_op_context > acos2atan_v;
  extern std::vector< gen_op_context > atan2asin_v;
  extern std::vector< gen_op_context > atan2acos_v;
  extern std::vector< gen_op_context > atrig2ln_v;
  extern std::vector<gen_op_context> trigcos_v;
  extern std::vector<gen_op_context> trigsin_v;
  extern std::vector<gen_op_context> trigtan_v;
  extern std::vector< gen_op_context > powexpand_v;
  extern std::vector< gen_op_context > exp2power_v;
  extern std::vector<unary_function_ptr> gamma_v;
  extern std::vector< gen_op_context > gamma2factorial_v;
  extern std::vector< unary_function_ptr > factorial_v;
  extern std::vector< gen_op_context > factorial2gamma_v;
  extern std::vector<unary_function_ptr> sign_floor_ceil_round_v;
  extern unary_function_ptr analytic_sommets[];

  bool need_parenthesis(const gen & g);
  void prod2frac(const gen & g,vecteur & num,vecteur & den);
  gen vecteur2prod(const vecteur & num);

  gen _multistring(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_multistring;

  std::string unquote(const std::string & s);

  
#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC

#endif // _GIAC_USUAL_H
