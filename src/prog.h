/* -*- mode:C++; compile-command: "g++ -I.. -g -c prog.cc" -*- */
/*
 *  Copyright (C) 2001 B. Parisse, Institut Fourier, 38402 St Martin d'Heres
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
#ifndef _GIAC_PROG_H
#define _GIAC_PROG_H
#include "first.h"
#include <vector>
#include <string>
#include <map>
#include "gen.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC
  void check_secure(); // in secure mode error
  struct user_function;
  struct module_info {
    std::vector<user_function> registered_names;
    void * handle;
    module_info():handle(0){};
    module_info(const std::vector<user_function> & r,void * h):registered_names(r),handle(h){};
  } ;
  typedef std::map< std::string, module_info> modules_tab;
  extern modules_tab giac_modules_tab;
  void set_decimal_digits(int n,GIAC_CONTEXT);
  int digits2bits(int n);
  int bits2digits(int n);
  // debug_info should be a vecteur containing
  // w[0]=function + args, w[2]= res of last evaluation, 
  // w[3] = source, w[4]=current_instruction
  // w[5] = watch vecteur, w[6] = watch values
  gen equaltosame(const gen & a);
  gen sametoequal(const gen & a);    
  int bind(const vecteur & vals,const vecteur & vars,context * & contextptr);
  bool leave(int protect,vecteur & vars,context * & contextptr);

  void increment_instruction(const vecteur & v,GIAC_CONTEXT);
  void increment_instruction(const gen & arg,GIAC_CONTEXT);
  void debug_print(const vecteur & arg,std::vector<std::string> & v,GIAC_CONTEXT);
  void debug_print(const gen & e,std::vector<std::string>  & v,GIAC_CONTEXT);
  std::string printasinnerbloc(const gen & feuille,GIAC_CONTEXT);
  std::string indent(GIAC_CONTEXT);
  // Return the names of variables that are not local in g
  // and the equality that are not used (warning = instead of := )
  std::string check_local_assign(const gen & g,GIAC_CONTEXT);
  extern const std::string _program_s;
  symbolic symb_program_sto(const gen & a,const gen & b,const gen & c,const gen & d,bool embedd=false,GIAC_CONTEXT=context0);
  symbolic symb_program(const gen & a,const gen & b,const gen & c,GIAC_CONTEXT);
  symbolic symb_program(const gen & args);
  gen _program(const gen & args,const gen & name,GIAC_CONTEXT);
  extern unary_function_ptr at_program ;

  extern const std::string _bloc_s;
  symbolic symb_bloc(const gen & args);
  extern unary_function_ptr at_bloc ;

  std::string printasifte(const gen & feuille,const std::string & sommetstr,GIAC_CONTEXT);
  extern const std::string _ifte_s;
  symbolic symb_ifte(const gen & test,const gen & oui, const gen & non);
  symbolic symb_ifte(const gen & e);
  gen _ifte(const gen & e,GIAC_CONTEXT);
  extern unary_function_ptr at_ifte ;
  gen _evalb(const gen & e,GIAC_CONTEXT);
  extern unary_function_ptr at_evalb ;
  gen _maple_if(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_maple_if ;

  gen symb_when(const gen & t,const gen & a,const gen & b);
  gen _when(const gen & e,GIAC_CONTEXT);
  extern unary_function_ptr at_when ;

  extern const std::string _for_s;
  symbolic symb_for(const gen & e);
  symbolic symb_for(const gen & a,const gen & b,const gen & c,const gen & d);
  gen _for(const gen & e,GIAC_CONTEXT);
  extern unary_function_ptr at_for ;

  gen symb_local(const gen & args,GIAC_CONTEXT);
  gen symb_local(const gen & a,const gen & b,GIAC_CONTEXT);
  extern const std::string _local_s;
  gen _local(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_local;

  gen symb_return(const gen & args);
  extern const std::string _return_s;
  gen _return(const gen & args);
  extern unary_function_ptr at_return;  

  gen symb_try_catch(const gen & args);
  extern const std::string _try_catch_s;
  gen _try_catch(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_try_catch;  

  gen symb_check_type(const gen & args);
  extern const std::string _check_type_s;
  gen _check_type(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_check_type;  

  gen symb_type(const gen & args);
  extern const std::string _type_s;
  gen _type(const gen & args);
  extern unary_function_ptr at_type;  

  gen symb_feuille(const gen & args);
  extern const std::string _feuille_s;
  gen _feuille(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_feuille;  
  extern const std::string _maple_op_s;
  gen _maple_op(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_maple_op;  

  gen symb_sommet(const gen & args);
  extern const std::string _sommet_s;
  gen _sommet(const gen & args);
  extern unary_function_ptr at_sommet;  

  gen subsop(const gen & g,const vecteur & v,GIAC_CONTEXT);
  gen subsop(const vecteur & g,const vecteur & v,const gen & sommet,GIAC_CONTEXT);
  extern const std::string _maple_subsop_s;
  gen _maple_subsop(const gen & args);
  extern unary_function_ptr at_maple_subsop;  

  extern const std::string _subsop_s;
  gen _subsop(const gen & args);
  extern unary_function_ptr at_subsop;  

  gen symb_append(const gen & args);
  extern const std::string _append_s;
  gen _append(const gen & args);
  extern unary_function_ptr at_append;  

  gen symb_prepend(const gen & args);
  extern const std::string _prepend_s;
  gen _prepend(const gen & args);
  extern unary_function_ptr at_prepend;  

  gen concat(const gen & g,bool glue_lines,GIAC_CONTEXT);
  gen symb_concat(const gen & args);
  extern const std::string _concat_s;
  gen _concat(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_concat;  
  
  gen symb_contains(const gen & args);
  extern const std::string _contains_s;
  gen _contains(const gen & args);
  extern unary_function_ptr at_contains;  

  gen symb_select(const gen & args);
  extern const std::string _select_s;
  gen _select(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_select;  

  gen symb_remove(const gen & args);
  extern const std::string _remove_s;
  gen _remove(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_remove;  

  gen symb_option(const gen & args);
  extern const std::string _option_s;
  gen _option(const gen & args);
  extern unary_function_ptr at_option;  

  gen symb_case(const gen & args);
  gen symb_case(const gen & a,const gen & b);
  extern const std::string _case_s;
  gen _case(const gen & args);
  extern unary_function_ptr at_case;  

  gen symb_rand(const gen & args);
  extern const std::string _rand_s;
  gen _rand(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_rand;  

  extern const std::string _srand_s;
  gen _srand(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_srand;  

  gen symb_char(const gen & args);
  extern const std::string _char_s;
  gen _char(const gen & args);
  extern unary_function_ptr at_char;  

  gen symb_asc(const gen & args);
  extern const std::string _asc_s;
  gen _asc(const gen & args);
  extern unary_function_ptr at_asc;  

  gen symb_map(const gen & args);
  extern const std::string _map_s;
  gen _map(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_map;  
  
  gen symb_apply(const gen & args);
  extern const std::string _apply_s;
  gen _apply(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_apply;  
  
  gen symb_makelist(const gen & args);
  extern const std::string _makelist_s;
  gen _makelist(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_makelist;  
  
  gen symb_interval(const gen & args);
  gen symb_interval(const gen & a,const gen & b);
  extern const std::string _interval_s;
  gen _interval(const gen & args);
  extern unary_function_ptr at_interval;
  gen symb_interval(const gen & a,const gen & b);
  
  gen symb_comment(const gen & args);
  extern const std::string _comment_s;
  gen _comment(const gen & args);
  extern unary_function_ptr at_comment;

  gen symb_throw(const gen & args);
  extern const std::string _throw_s;
  gen _throw(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_throw;

  gen symb_union(const gen & args);
  extern const std::string _union_s;
  gen _union(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_union;

  gen symb_intersect(const gen & args);
  extern const std::string _intersect_s;
  gen _intersect(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_intersect;
  extern const std::string _inter_s;
  gen _inter(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_inter;

  gen symb_minus(const gen & args);
  extern const std::string _minus_s;
  gen _minus(const gen & args);
  extern unary_function_ptr at_minus;

  gen symb_dollar(const gen & args);
  extern const std::string _dollar_s;
  gen _dollar(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_dollar;

  gen symb_makemat(const gen & args);
  extern const std::string _makemat_s;
  gen _makemat(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_makemat;

  gen symb_compose(const gen & args);
  extern const std::string _compose_s;
  gen _compose(const gen & args);
  extern unary_function_ptr at_compose;

  gen _composepow(const gen & args);
  extern unary_function_ptr at_composepow;

  gen symb_has(const gen & args);
  extern const std::string _has_s;
  gen _has(const gen & args);
  extern unary_function_ptr at_has;

  gen symb_args(const gen & args);
  extern const std::string _args_s;
  gen _args(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_args;

  gen symb_lname(const gen & args);
  extern const std::string _lname_s;

  vecteur lidnt(const gen & args);
  void lidnt(const gen & args,vecteur & res);

  gen _lname(const gen & args);
  extern unary_function_ptr at_lname;
  
  extern const std::string _halt_s;
  gen _halt(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_halt;

  extern const std::string _kill_s;
  gen _kill(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_kill;

  extern const std::string _cont_s;
  gen _cont(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_cont;

  extern const std::string _sst_s;
  gen _sst(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_sst;

  extern const std::string _sst_in_s;
  gen _sst_in(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_sst_in;

  extern const std::string _debug_s;
  gen _debug(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_debug;

  extern const std::string _watch_s;
  gen _watch(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_watch;

  extern const std::string _rmwatch_s;
  gen _rmwatch(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_rmwatch;

  extern const std::string _breakpoint_s;
  gen _breakpoint(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_breakpoint;

  extern const std::string _rmbreakpoint_s;
  gen _rmbreakpoint(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_rmbreakpoint;

  void debug_loop(const gen &res,GIAC_CONTEXT);

  extern const std::string _backquote_s;
  gen _backquote(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_backquote;

  extern const std::string _double_deux_points_s;
  gen _double_deux_points(const gen & args);
  gen symb_double_deux_points(const gen & args);
  extern unary_function_ptr at_double_deux_points;

  gen _maple2mupad(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_maple2mupad;

  gen _maple2xcas(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_maple2xcas;

  gen _mupad2maple(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_mupad2maple;

  gen _mupad2xcas(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_mupad2xcas;

  gen _cd(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_cd;

  gen _pwd(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_pwd;

  gen _scientific_format(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_scientific_format;

  gen _xcas_mode(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_xcas_mode;
  extern unary_function_ptr at_maple_mode;

  gen _all_trig_solutions(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_all_trig_solutions;

  gen _ntl_on(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_ntl_on;

  gen _complex_mode(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_complex_mode;

  gen _angle_radian(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_angle_radian;
 
  gen _epsilon(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_epsilon;

  gen _proba_epsilon(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_proba_epsilon;

  gen _complex_variables(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_complex_variables;

  gen _approx_mode(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_approx_mode;

  gen _threads(const gen & args);
  extern unary_function_ptr at_threads;

  gen _cas_setup(const gen & args,GIAC_CONTEXT);
  void parent_cas_setup(GIAC_CONTEXT); // send current cas_setup to parent
  extern unary_function_ptr at_cas_setup;
  void cas_setup(const vecteur & v,GIAC_CONTEXT);
  vecteur cas_setup(GIAC_CONTEXT);

  gen _Digits(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_Digits;

  gen _insmod(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_insmod;

  gen _rmmod(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_rmmod;

  gen _lsmod(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_lsmod;

  gen _virgule(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_virgule;

  bool is_binary(const gen & args);
  void check_binary(const gen & args,gen & a,gen & b);

  gen _sort(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_sort;

  gen _ans(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_ans;

  gen _quest(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_quest;

  std::vector<int> float2continued_frac(double d_orig,double eps);
  gen continued_frac2gen(std::vector<int> v,double d_orig,double eps,GIAC_CONTEXT);
  gen _convert(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_convert;

  gen _deuxpoints(const gen & args);
  extern unary_function_ptr at_deuxpoints;

  gen quote_read(const gen & args,GIAC_CONTEXT); // read in a file and return non evaled
  gen _read(const gen & args,GIAC_CONTEXT); // read in a file and return evaled
  extern unary_function_ptr at_read;

  gen _write(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_write;

  gen _save_history(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_save_history;

  gen symb_findhelp(const gen & args);
  gen _findhelp(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_findhelp;

  gen _member(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_member;

  gen _tablefunc(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_tablefunc;

  gen _tableseq(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_tableseq;

  gen protecteval(const gen & g,int level,GIAC_CONTEXT);

  gen _nodisp(const gen & args);
  extern unary_function_ptr at_nodisp;

  gen _unapply(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_unapply;

  gen _makevector(const gen & args);
  extern unary_function_ptr at_makevector;

  gen _matrix(const gen & args);
  extern unary_function_ptr at_matrix;

  gen _makesuite(const gen & args);
  extern unary_function_ptr at_makesuite;

  gen _break(const gen & args);
  extern unary_function_ptr at_break;

  gen _continue(const gen & args);
  extern unary_function_ptr at_continue;

  gen _label(const gen & args);
  extern unary_function_ptr at_label;

  gen _goto(const gen & args);
  extern unary_function_ptr at_goto;

  gen _tilocal(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_tilocal;

  gen inputform_post_analysis(const vecteur & v,const gen & res,GIAC_CONTEXT);
  vecteur inputform_pre_analysis(const gen & g,GIAC_CONTEXT);
  gen _inputform(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_inputform;
  class unary_function_unary;
  class unary_function_eval;
  extern unary_function_eval __inputform;

  gen _choosebox(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_choosebox;

  gen _output(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_output;
  gen _input(const gen & args,bool textinput,GIAC_CONTEXT);

  extern const std::string _nop_s;
  gen _nop(const gen & args);
  extern unary_function_ptr at_nop;

  std::string printastifunction(const gen & feuille,const std::string & sommetstr,GIAC_CONTEXT);

  gen _Dialog(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_Dialog;

  gen _Title(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_Title;

  gen _Text(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_Text;

  gen _Request(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_Request;

  gen _DropDown(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_DropDown;

  gen _Popup(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_Popup;

  gen _expr(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_expr;

  gen _string(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_string;

  gen _part(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_part;

  gen _Pause(const gen & args);
  extern unary_function_ptr at_Pause;

  gen _Row(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_Row;
  extern unary_function_eval __Row;

  gen _Col(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_Col;
  extern unary_function_eval __Col;

  gen _DelVar(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_DelVar;

  gen prod(const gen &,const gen &);
  gen somme(const gen &,const gen &);
  gen _pointprod(const gen & args);
  extern unary_function_ptr at_pointprod;

  gen _pointdivision(const gen & args);
  extern unary_function_ptr at_pointdivision;

  gen _pointpow(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_pointpow;

  gen _pourcent(const gen & args);
  extern unary_function_ptr at_pourcent;

  gen _hash(const gen & args);
  extern unary_function_ptr at_hash;

  // used to update IO screen and graph inside progs
  gen _interactive(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_interactive;
  extern unary_function_eval __interactive;
  extern bool user_screen; 
  extern int user_screen_io_x,user_screen_io_y,user_screen_fontsize;

  std::string printassuffix(const gen & feuille,const std::string & sommetstr,GIAC_CONTEXT);

  // translate TI escape sequence
  std::string tiasc_translate(const std::string & s);
  gen matrix_apply(const gen & a,const gen & b,gen (* f) (const gen &, const gen &) );
  gen matrix_apply(const gen & a,const gen & b,GIAC_CONTEXT,gen (* f) (const gen &, const gen &,GIAC_CONTEXT) );

  // v=[ [idnt,value] ... ]
  // search g in v if found return value
  // else return g evaluated
  // and add g to the list according to add_to_folder
  gen find_in_folder(vecteur & v,const gen & g);
  extern gen current_folder_name; // must be an idnt (or a path)
  gen getfold(const gen & g); // translate 0 to "main"

  gen _ti_semi(const gen & args);
  extern unary_function_ptr at_ti_semi;

  extern unary_function_unary __keyboard;
  extern unary_function_ptr at_keyboard;
  extern unary_function_unary __widget_size;
  extern unary_function_ptr at_widget_size;

  extern unary_function_eval __current_sheet;
  extern unary_function_ptr at_current_sheet;

  extern unary_function_unary __window_switch;
  extern unary_function_ptr at_window_switch;

  extern unary_function_unary __maple_lib;
  extern unary_function_ptr at_maple_lib;
  extern unary_function_ptr at_unit;
  extern unary_function_ptr at_maple_root;
  gen symb_unit(const gen & a,const gen & b,GIAC_CONTEXT);
  gen symb_interrogation(const gen & e1,const gen & e3);
  std::string printasDigits(const gen & feuille,const std::string & sommetstr,GIAC_CONTEXT);
  bool first_ascend_sort(const gen & a,const gen & b);
  bool first_descend_sort(const gen & a,const gen & b);

  extern unary_function_ptr at_user_operator;
  gen user_operator(const gen & g,GIAC_CONTEXT);
  gen _SetFold(const gen & g,GIAC_CONTEXT);
  extern unary_function_ptr at_SetFold;

  gen simplifier(const gen & g,GIAC_CONTEXT);
  gen _simplifier(const gen & g,GIAC_CONTEXT);
  // Unit management
  gen unitpow(const gen & g,const gen & exponent);
  gen mksa_reduce(const gen & g,GIAC_CONTEXT);
  gen chk_not_unit(const gen & g);
  gen find_or_make_symbol(const std::string & s,GIAC_CONTEXT);
  extern gen _m_unit;
  extern gen _kg_unit;
  extern gen _s_unit;
  extern gen _A_unit;
  extern gen _K_unit;
  extern gen _mol_unit;
  extern gen _molK_unit;
  extern gen _cd_unit;
  extern gen _E_unit;
  // other metric units in m,kg,s,A
  extern gen _Bq_unit;
  extern gen _C_unit;
  extern gen _F_unit;
  extern gen _Gy_unit;
  extern gen _H_unit;
  extern gen _Hz_unit;
  extern gen _J_unit;
  extern gen _mho_unit;
  extern gen _N_unit;
  extern gen _Ohm_unit;
  extern gen _Pa_unit;
  extern gen _r_unit;
  extern gen _S_unit;
  extern gen _st_unit;
  extern gen _Sv_unit;
  extern gen _T_unit;
  extern gen _V_unit;
  extern gen _W_unit;
  extern gen _Wb_unit;
  // useful non metric units
  extern gen _a_unit;
  extern gen _acre_unit;
  extern gen _arcmin_unit;
  extern gen _arcs_unit;
  extern gen _atm_unit;
  extern gen _au_unit;
  extern gen _angstrom_unit;
  extern gen _b_unit;
  extern gen _bar_unit;
  extern gen _bbl_unit;
  extern gen _Btu_unit;
  extern gen _cal_unit;
  extern gen _chain_unit;
  extern gen _Ci_unit;
  extern gen _ct_unit;
  // extern gen _°_unit;
  extern gen _d_unit;
  extern gen _dB_unit;
  extern gen _dyn_unit;
  extern gen _erg_unit;
  extern gen _eV_unit;
  // extern gen _°F_unit;
  extern gen _fath_unit;
  extern gen _fbm_unit;
  // extern gen _fc_unit;
  extern gen _Fdy_unit;
  extern gen _fermi_unit;
  extern gen _flam_unit;
  extern gen _ft_unit;
  extern gen _ftUS_unit;
  extern gen _g_unit;
  extern gen _gal_unit;
  extern gen _galC_unit;
  extern gen _galUK_unit;
  extern gen _gf_unit;
  extern gen _gmol_unit;
  extern gen _grad_unit;
  extern gen _grain_unit;
  extern gen _ha_unit;
  extern gen _h_unit;
  extern gen _hp_unit;
  extern gen _in_unit;
  extern gen _inHg_unit;
  extern gen _inH2O_unit;
  extern gen _FF_unit;
  extern gen _kip_unit;
  extern gen _knot_unit;
  extern gen _kph_unit;
  extern gen _l_unit;
  extern gen _lam_unit;
  extern gen _lb_unit;
  extern gen _lbf_unit;
  extern gen _lbmol_unit;
  extern gen _lbt_unit;
  extern gen _lyr_unit;
  extern gen _mi_unit;
  extern gen _mil_unit;
  extern gen _min_unit;
  extern gen _miUS_unit;
  extern gen _mmHg_unit;
  extern gen _mph_unit;
  extern gen _nmi_unit;
  extern gen _oz_unit;
  extern gen _ozfl_unit;
  extern gen _ozt_unit;
  extern gen _ozUK_unit;
  extern gen _P_unit;
  extern gen _pc_unit;
  extern gen _pdl_unit;
  extern gen _pk_unit;
  extern gen _psi_unit;
  extern gen _pt_unit;
  extern gen _qt_unit;
  extern gen _R_unit;
  extern gen _rad_unit;
  extern gen _rd_unit;
  extern gen _rem_unit;
  extern gen _rpm_unit;
  extern gen _sb_unit;
  extern gen _slug_unit;
  extern gen _St_unit;
  extern gen _t_unit;
  extern gen _tbsp_unit;
  extern gen _therm_unit;
  extern gen _ton_unit;
  extern gen _tonUK_unit;
  extern gen _torr_unit;
  extern gen _u_unit;
  extern gen _yd_unit;
  extern gen _yr_unit;
  // Physical constants
  extern gen cst_hbar;
  extern gen cst_clightspeed;
  extern gen cst_ga;
  extern gen cst_IO;
  extern gen cst_epsilonox;
  extern gen cst_epsilonsi;
  extern gen cst_qepsilon0;
  extern gen cst_epsilon0q;
  extern gen cst_kq;
  extern gen cst_c3;
  extern gen cst_lambdac;
  extern gen cst_f0;
  extern gen cst_lambda0;
  extern gen cst_muN;
  extern gen cst_muB;
  extern gen cst_a0;
  extern gen cst_Rinfinity;
  extern gen cst_Faraday;
  extern gen cst_phi;
  extern gen cst_alpha;
  extern gen cst_mpme;
  extern gen cst_mp;
  extern gen cst_qme;
  extern gen cst_me;
  extern gen cst_qe;
  extern gen cst_hPlanck;
  extern gen cst_G;
  extern gen cst_mu0;
  extern gen cst_epsilon0;
  extern gen cst_sigma;
  extern gen cst_StdP;
  extern gen cst_StdT;
  extern gen cst_Rydberg;
  extern gen cst_Vm;
  extern gen cst_kBoltzmann;
  extern gen cst_NA;
  extern unary_function_ptr binary_op_tab[];

  extern unary_function_ptr at_piecewise;

  extern unary_function_ptr at_geo2d ;
  extern unary_function_ptr at_geo3d ;
  extern unary_function_ptr at_spreadsheet ;

  std::string print_program_syntax(int maple_mode);
  gen when2piecewise(const gen & g,GIAC_CONTEXT);
  gen when2sign(const gen & g,GIAC_CONTEXT);
  gen piecewise2when(const gen & g,GIAC_CONTEXT);

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC

#endif // _GIAC_PROG_H

