 /* -*- mode: C++; compile-command: "flex input_lexer.ll && g++ -g -I.. -c input_lexer.cc -Wall" -*- */
/** @file input_lexer.ll
 *
 *  Lexical analyzer definition for reading expressions.
 *  Note Maple input should be processed replacing # with // and { } for set
 *  This file must be processed with flex. */

/*
 *  Early version modified from GiNaC by B. Parisse (C) 2001, 7
 *  GiNaC Copyright (C) 1999-2000 Johannes Gutenberg University Mainz, Germany
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


/*
 * The lexer will first check for static patterns and strings (defined below)
 * If a match is not found, it calls find_or_make_symbol
 * This function looks first if the string should be translated
 * (e.g. add a prefix from the export table)
 * then look in lexer_functions for a match, then look in sym_tab
 * if not found in sym_tab, a new identificateur is created & added in sym_tab
 * Functions in lexer_functions are added during the construction
 * of the corresponding unary_functions using lexer_functions_register
 */


/*
 *  Definitions
 */

%{
#include "first.h"
#include <iostream>
#include <stdexcept>

#include "input_lexer.h"
#include "help.h"
#include "gen.h"
#include "identificateur.h"
#include "usual.h"
#include "derive.h"
#include "series.h"
#include "intg.h"
#include "sym2poly.h"
#include "moyal.h"
#include "subst.h"
#include "vecteur.h"
#include "modpoly.h"
#include "lin.h"
#include "solve.h"
#include "ifactor.h"
#include "alg_ext.h"
#include "gauss.h"
#include "isom.h"
#include "plot.h"
#include "prog.h"
#include "rpn.h"
#include "ezgcd.h"
#include "tex.h"
#include "risch.h"
#include "input_parser.h"    


  using namespace std;
  using namespace giac;
#ifndef NO_NAMESPACE_GIAC
  namespace giac {
#endif // ndef NO_NAMESPACE_GIAC
    sym_tab & syms(){
      static sym_tab * ans=new sym_tab;
      return * ans;
    }

    sym_tab & lexer_functions(){
      static sym_tab * ans=new sym_tab;
      return * ans;
    }

    std::vector<int> lexer_localization_vector;
    std::map<std::string,std::string> & lexer_localization_map(){
      static std::map<std::string,std::string> * ans = new std::map<std::string,std::string>;
      return * ans;
    }
    std::multimap<std::string,giac::localized_string> back_lexer_localization_map;
    // lexer_localization_vector is the list of languages currently translated
    // lexer_localization_map translates keywords from the locale to giac 
    // back_lexer_localization_map lists for a giac keyword the translations

    std::map<std::string,std::vector<std::string> > & lexer_translator (){
      static std::map<std::string,std::vector<std::string> > * ans = new std::map<std::string,std::vector<std::string> >;
      return * ans;
    }
    // lexer_translator will be updated when export/with is called
    // To each string (w/o ::) in a given library, 
    // If it exists, we push_back the full string (with ::)
    // If not we create a vector with the full string
    // If a library is unexported we remove the corresponding entry in the 
    // vector and remove the entry if the vector is empty
    std::map<std::string,std::vector<std::string> > & library_functions (){
      static std::map<std::string,std::vector<std::string> > * ans=new std::map<std::string,std::vector<std::string> >;
      return *ans;
    }
    // First string is the library name, second is the vector of function names
    // User defined relations
    vector<user_function> & registered_lexer_functions(){
      static vector<user_function> * ans = new vector<user_function>;
      return * ans;
    }


#ifndef NO_NAMESPACE_GIAC
  } // namespace giac
#endif // ndef NO_NAMESPACE_GIAC

%}

%option reentrant bison-bridge
%option outfile="input_lexer.cc"
%option header-file="lexer.h"
%option noyywrap
%option prefix="giac_yy"

	/* Abbreviations */
D	[0-9]
E	[eE][eE]?[-+]?{D}+
A	[a-zA-Z~µ·‡‰„‚ÈËÍÎÓÌÏÔÛÚˆÙı˙˘¸˚˝ˇÒ¿¡√¬ƒ…» ÀÕÃŒœ”“‘’÷⁄Ÿ€‹›ÁΩ]
AN	[0-9a-zA-Z_.~µ·‡‰„‚ÈËÍÎÓÌÏÔÛÚˆÙı˙˘¸˚˝ˇÒ¿¡√¬ƒ…» ÀÕÃŒœ”“‘’÷⁄Ÿ€‹›ÁΩ] 
        /* If changed, modify isalphan in help.cc FIXME is . allowed inside alphanumeric */

%x comment
%x comment_hash
%x str
%x backquote
/*
 *  Lexical rules
 */

%%

[ \t\\]+			/* skip whitespace */
\n                increment_lexer_line_number(yyextra); //cerr << "Scanning line " << lexer_line_number(yyextra) << endl;
  /* Strings */
  /* \"[^\"]*\"        yylval = string2gen( giac_yytext); return T_STRING; */
\"                BEGIN(str); comment_s("",yyextra);
<str>\"\"         increment_comment_s('"',yyextra);
<str>\"        {  index_status(yyextra)=1; BEGIN(INITIAL); 
                  (*yylval)=string2gen(comment_s(yyextra),false); 
                  return T_STRING; }
<str>\n        increment_comment_s('\n',yyextra); increment_lexer_line_number(yyextra);
<str>\\[0-7]{1,3} {
                   /* octal escape sequence */
                   int result;
                   (void) sscanf( yytext + 1, "%o", &result );
                   increment_comment_s(char(result & 0xff),yyextra);
                   }
<str>\\[0-9]+      {
                   /* generate error - bad escape sequence; something
                    * like '\48' or '\0777777'
                    */
                   }
<str>\\n  increment_comment_s('\n',yyextra);
<str>\\t  increment_comment_s('\t',yyextra);
<str>\\r  increment_comment_s('\r',yyextra);
<str>\\b  increment_comment_s('\b',yyextra);
<str>\\f  increment_comment_s('\f',yyextra);
<str>\\(.|\n)  increment_comment_s(yytext[1],yyextra);
<str>[^\\\n\"]+ increment_comment_s(yytext,yyextra);
`                  if (rpn_mode){ index_status(yyextra)=0; return T_ACCENTGRAVE; } else { BEGIN(backquote); comment_s("",yyextra); }
<backquote>\n      increment_comment_s('\n',yyextra); increment_lexer_line_number(yyextra);
<backquote>[^\n`]+       increment_comment_s(yytext,yyextra);
<backquote>`       {  index_status(yyextra)=1; BEGIN(INITIAL); 
     return find_or_make_symbol(comment_s(yyextra),(*yylval),yyextra); }

"://"[^\n]*\n      index_status(yyextra)=0; increment_lexer_line_number(yyextra);
"//"[^\n]*\n      index_status(yyextra)=0; increment_lexer_line_number(yyextra);/* (*yylval) = string2gen('"'+string(giac_yytext).substr(2,string(giac_yytext).size()-3)+'"');   return T_COMMENT; */
"/*"              BEGIN(comment); comment_s(yyextra)="";

<comment>[^*\n]*        comment_s(yyextra)+=yytext; /* eat anything that's not a '*' */
<comment>"*"+[^*/\n]*   comment_s(yyextra)+=yytext; /* eat up '*'s not followed by '/'s */
<comment>\n             comment_s(yyextra) += '\n'; increment_lexer_line_number(yyextra); cerr << "(Comment) scanning line " << lexer_line_number(yyextra) << endl;
<comment>"*"+"/"        BEGIN(INITIAL); index_status(yyextra)=0; /* (*yylval) = string2gen(comment_s(yyextra),false); return T_COMMENT; */
"#++"[^*]*"++#"         index_status(yyextra)=0; /* (*yylval) = string2gen('"'+string(yytext).substr(3,string(yytext).size()-6)+'"'); return T_COMMENT; */
"#--"[^*]*"--#"         index_status(yyextra)=0; /* (*yylval) = string2gen('"'+string(yytext).substr(3,string(yytext).size()-6)+'"'); return T_COMMENT; */

"?"                     if (index_status(yyextra)) return T_INTERROGATION; else return T_HELP;
"_"                     return T_UNIT;
"'"                     if (opened_quote(yyextra)) { opened_quote(yyextra)=0; return T_QUOTE; } if (index_status(yyextra) && !in_rpn(yyextra) && xcas_mode(yyextra)!= 1) return T_PRIME; opened_quote(yyextra)=1; return T_QUOTE;
";"			index_status(yyextra)=0; if (xcas_mode(yyextra)==3) return TI_SEMI; (*yylval)=0; return T_SEMI;
"ß"                  index_status(yyextra)=0; if (xcas_mode(yyextra)==3) return T_SEMI; return TI_SEMI;
":"			if (spread_formula(yyextra)) return T_DEUXPOINTS; if ( xcas_mode(yyextra)==3 ) { index_status(yyextra)=0; return TI_DEUXPOINTS; }  index_status(yyextra)=0; if (xcas_mode(yyextra)>0) { (*yylval)=1; return T_SEMI; } else return T_DEUXPOINTS;
":;"                    (*yylval)=1; return T_SEMI;
"::"                    return T_DOUBLE_DEUX_POINTS;

			/* special values */
"i"			index_status(yyextra)=1; if (xcas_mode(yyextra) > 0 ) { (*yylval)=i__IDNT_e; return T_SYMBOL; } else { (*yylval) = cst_i; return T_LITERAL;};
\xa1                     index_status(yyextra)=1; (*yylval) = cst_i; return T_LITERAL;
"I"                     index_status(yyextra)=1; if (xcas_mode(yyextra)==0 || xcas_mode(yyextra)==3 || rpn_mode) { (*yylval)=I__IDNT_e; return T_SYMBOL; } else { (*yylval) = cst_i; return T_LITERAL; };
"%i"			index_status(yyextra)=1; (*yylval) = cst_i; return T_LITERAL;
"%e"			index_status(yyextra)=1; (*yylval) = symbolic(at_exp,1); return T_LITERAL;
"%pi"			index_status(yyextra)=1; (*yylval) = cst_pi; return T_LITERAL;
"pi"			index_status(yyextra)=1; (*yylval) = cst_pi; return T_LITERAL;
"Pi"			index_status(yyextra)=1; (*yylval) = cst_pi; return T_LITERAL;
"PI"			index_status(yyextra)=1; (*yylval) = cst_pi; return T_LITERAL;
"euler_gamma"		index_status(yyextra)=1; (*yylval) = cst_euler_gamma; return T_LITERAL;
"infinity"		index_status(yyextra)=1; (*yylval) = unsigned_inf; return T_LITERAL;
"inf"		index_status(yyextra)=1; (*yylval) = plus_inf; return T_LITERAL;
"unsigned_inf"		index_status(yyextra)=1; (*yylval) = unsigned_inf; return T_LITERAL;
"plus_inf"		index_status(yyextra)=1; (*yylval) = plus_inf; return T_LITERAL;
"minus_inf"		index_status(yyextra)=1; (*yylval) = minus_inf; return T_LITERAL;
"undef"		        index_status(yyextra)=1; (*yylval) = undef; return T_LITERAL;
"§"                     return T_END_INPUT;

               /* integer values */
"DOM_int"                  index_status(yyextra)=0; (*yylval) = _INT_; (*yylval).subtype=_INT_TYPE; return T_TYPE_ID;
"DOM_INT"               index_status(yyextra)=0; (*yylval) = _ZINT; (*yylval).subtype=_INT_TYPE; return T_TYPE_ID;
"integer"               index_status(yyextra)=0; (*yylval) = _ZINT; (*yylval).subtype=_INT_TYPE; return T_TYPE_ID;
"posint"               index_status(yyextra)=0; (*yylval) = _POSINT; (*yylval).subtype=_INT_MAPLECONVERSION; return T_TYPE_ID;
"negint"               index_status(yyextra)=0; (*yylval) = _NEGINT; (*yylval).subtype=_INT_MAPLECONVERSION; return T_TYPE_ID;
"nonposint"               index_status(yyextra)=0; (*yylval) = _NONPOSINT; (*yylval).subtype=_INT_MAPLECONVERSION; return T_TYPE_ID;
"nonnegint"               index_status(yyextra)=0; (*yylval) = _NONNEGINT; (*yylval).subtype=_INT_MAPLECONVERSION; return T_TYPE_ID;
"DOM_COMPLEX"               index_status(yyextra)=0; (*yylval) = _CPLX; (*yylval).subtype=_INT_TYPE; return T_TYPE_ID;
"complex"               index_status(yyextra)=0; (*yylval) = _CPLX; (*yylval).subtype=_INT_TYPE; return T_TYPE_ID;
"DOM_FLOAT"                index_status(yyextra)=0; (*yylval) = _DOUBLE_; (*yylval).subtype=_INT_TYPE; return T_TYPE_ID;
"float"                index_status(yyextra)=0; (*yylval) = _DOUBLE_; (*yylval).subtype=_INT_TYPE; return T_TYPE_ID;
"DOM_LIST"              index_status(yyextra)=0; (*yylval) = _VECT; (*yylval).subtype=_INT_TYPE; return T_TYPE_ID;
"DOM_MATRIX"              index_status(yyextra)=0; (*yylval) = _VECT; (*yylval).subtype=_INT_TYPE; return T_TYPE_ID;
"DOM_IDENT"            index_status(yyextra)=0; (*yylval) = _IDNT; (*yylval).subtype=_INT_TYPE; return T_TYPE_ID;
"DOM_SYMBOLIC"              index_status(yyextra)=0; (*yylval) = _SYMB; (*yylval).subtype=_INT_TYPE; return T_TYPE_ID;
"algebraic"              index_status(yyextra)=0; (*yylval) = _SYMB; (*yylval).subtype=_INT_TYPE; return T_TYPE_ID;
"DOM_RAT"              index_status(yyextra)=0; (*yylval) = _FRAC; (*yylval).subtype=_INT_TYPE; return T_TYPE_ID;
"rational"              index_status(yyextra)=0; (*yylval) = _FRAC; (*yylval).subtype=_INT_TYPE; return T_TYPE_ID;
"DOM_STRING"                index_status(yyextra)=0; (*yylval) = _STRNG; (*yylval).subtype=_INT_TYPE; return T_TYPE_ID;
"DOM_FUNC"                index_status(yyextra)=0; (*yylval) = _FUNC; (*yylval).subtype=_INT_TYPE; return T_TYPE_ID;
"DOM_LONGFLOAT"                index_status(yyextra)=0; (*yylval) = _REAL; (*yylval).subtype=_INT_TYPE; return T_TYPE_ID;
"operator"              if (xcas_mode(yyextra)==2){ (*yylval) = gen(at_user_operator,6); index_status(yyextra)=0; return T_UNARY_OP; }  index_status(yyextra)=0; (*yylval) = _FUNC; (*yylval).subtype=_INT_TYPE; return T_TYPE_ID;
"unfactored"         index_status(yyextra)=1; (*yylval) = _UNFACTORED; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"sans_factoriser"         index_status(yyextra)=1; (*yylval) = _UNFACTORED; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"newton_solver"         index_status(yyextra)=1; (*yylval) = _NEWTON_SOLVER; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"dnewton_solver"         index_status(yyextra)=1; (*yylval) = _DNEWTON_SOLVER; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"newtonj_solver"         index_status(yyextra)=1; (*yylval) = _NEWTONJ_SOLVER; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"hybridj_solver"         index_status(yyextra)=1; (*yylval) = _HYBRIDJ_SOLVER; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"hybridsj_solver"         index_status(yyextra)=1; (*yylval) = _HYBRIDSJ_SOLVER; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"hybrid_solver"         index_status(yyextra)=1; (*yylval) = _HYBRID_SOLVER; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"hybrids_solver"         index_status(yyextra)=1; (*yylval) = _HYBRIDS_SOLVER; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"secant_solver"         index_status(yyextra)=1; (*yylval) = _SECANT_SOLVER; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"steffenson_solver"         index_status(yyextra)=1; (*yylval) = _STEFFENSON_SOLVER; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"bisection_solver"         index_status(yyextra)=1; (*yylval) = _BISECTION_SOLVER; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"falsepos_solver"         index_status(yyextra)=1; (*yylval) = _FALSEPOS_SOLVER; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"brent_solver"         index_status(yyextra)=1; (*yylval) = _BRENT_SOLVER; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"golub_reinsch_decomp"         index_status(yyextra)=1; (*yylval) = _GOLUB_REINSCH_DECOMP; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"golub_reinsch_mod_decomp"         index_status(yyextra)=1; (*yylval) = _GOLUB_REINSCH_MOD_DECOMP; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"jacobi_decomp"         index_status(yyextra)=1; (*yylval) = _JACOBI_DECOMP; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"minor_det"         index_status(yyextra)=1; (*yylval) = _MINOR_DET; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"bareiss"         index_status(yyextra)=1; (*yylval) = _BAREISS; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"rational_det"         index_status(yyextra)=1; (*yylval) = _RATIONAL_DET; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"fadeev"         index_status(yyextra)=1; (*yylval) = _FADEEV; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"keep_pivot"         index_status(yyextra)=1; (*yylval) = _KEEP_PIVOT; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"trapezoid"         index_status(yyextra)=1; (*yylval) = _TRAPEZE; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"trapeze"         index_status(yyextra)=1; (*yylval) = _TRAPEZE; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"middle_point"         index_status(yyextra)=1; (*yylval) = _POINT_MILIEU; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"point_milieu"         index_status(yyextra)=1; (*yylval) = _POINT_MILIEU; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"right_rectangle"         index_status(yyextra)=1; (*yylval) = _RECTANGLE_DROIT; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"rectangle_droit"         index_status(yyextra)=1; (*yylval) = _RECTANGLE_DROIT; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"left_rectangle"         index_status(yyextra)=1; (*yylval) = _RECTANGLE_GAUCHE; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"rectangle_gauche"         index_status(yyextra)=1; (*yylval) = _RECTANGLE_GAUCHE; (*yylval).subtype=_INT_SOLVER;return T_NUMBER;
"filled"         index_status(yyextra)=1; (*yylval) = _FILL_POLYGON; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"rempli"         index_status(yyextra)=1; (*yylval) = _FILL_POLYGON ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"line_width_1"         index_status(yyextra)=1; (*yylval) = 0 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"line_width_2"         index_status(yyextra)=1; (*yylval) = _LINE_WIDTH_2 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"line_width_3"         index_status(yyextra)=1; (*yylval) = _LINE_WIDTH_3 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"line_width_4"         index_status(yyextra)=1; (*yylval) = _LINE_WIDTH_4 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"line_width_5"         index_status(yyextra)=1; (*yylval) = _LINE_WIDTH_5 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"line_width_6"         index_status(yyextra)=1; (*yylval) = _LINE_WIDTH_6 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"line_width_7"         index_status(yyextra)=1; (*yylval) = _LINE_WIDTH_7 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"line_width_8"         index_status(yyextra)=1; (*yylval) = _LINE_WIDTH_8 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"epaisseur_ligne_1"         index_status(yyextra)=1; (*yylval) = 0 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"epaisseur_ligne_2"         index_status(yyextra)=1; (*yylval) = _LINE_WIDTH_2 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"epaisseur_ligne_3"         index_status(yyextra)=1; (*yylval) = _LINE_WIDTH_3 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"epaisseur_ligne_4"         index_status(yyextra)=1; (*yylval) = _LINE_WIDTH_4 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"epaisseur_ligne_5"         index_status(yyextra)=1; (*yylval) = _LINE_WIDTH_5 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"epaisseur_ligne_6"         index_status(yyextra)=1; (*yylval) = _LINE_WIDTH_6 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"epaisseur_ligne_7"         index_status(yyextra)=1; (*yylval) = _LINE_WIDTH_7 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"epaisseur_ligne_8"         index_status(yyextra)=1; (*yylval) = _LINE_WIDTH_8 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"point_width_1"         index_status(yyextra)=1; (*yylval) = 0 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"point_width_2"         index_status(yyextra)=1; (*yylval) = _POINT_WIDTH_2 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"point_width_3"         index_status(yyextra)=1; (*yylval) = _POINT_WIDTH_3 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"point_width_4"         index_status(yyextra)=1; (*yylval) = _POINT_WIDTH_4 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"point_width_5"         index_status(yyextra)=1; (*yylval) = _POINT_WIDTH_5 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"point_width_6"         index_status(yyextra)=1; (*yylval) = _POINT_WIDTH_6 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"point_width_7"         index_status(yyextra)=1; (*yylval) = _POINT_WIDTH_7 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"point_width_8"         index_status(yyextra)=1; (*yylval) = _POINT_WIDTH_8 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"epaisseur_point_1"         index_status(yyextra)=1; (*yylval) = 0 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"epaisseur_point_2"         index_status(yyextra)=1; (*yylval) = _POINT_WIDTH_2 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"epaisseur_point_3"         index_status(yyextra)=1; (*yylval) = _POINT_WIDTH_3 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"epaisseur_point_4"         index_status(yyextra)=1; (*yylval) = _POINT_WIDTH_4 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"epaisseur_point_5"         index_status(yyextra)=1; (*yylval) = _POINT_WIDTH_5 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"epaisseur_point_6"         index_status(yyextra)=1; (*yylval) = _POINT_WIDTH_6 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"epaisseur_point_7"         index_status(yyextra)=1; (*yylval) = _POINT_WIDTH_7 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"epaisseur_point_8"         index_status(yyextra)=1; (*yylval) = _POINT_WIDTH_8 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"hidden_name"         index_status(yyextra)=1; (*yylval) = _HIDDEN_NAME ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"nom_cache"         index_status(yyextra)=1; (*yylval) = _HIDDEN_NAME ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"quadrant1"         index_status(yyextra)=1; (*yylval) = _QUADRANT1 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"quadrant2"         index_status(yyextra)=1; (*yylval) = _QUADRANT2 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"quadrant3"         index_status(yyextra)=1; (*yylval) = _QUADRANT3 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"quadrant4"         index_status(yyextra)=1; (*yylval) = _QUADRANT4 ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"rhombus_point"         index_status(yyextra)=1; (*yylval) = _POINT_LOSANGE ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"point_losange"         index_status(yyextra)=1; (*yylval) = _POINT_LOSANGE ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"plus_point"         index_status(yyextra)=1; (*yylval) = _POINT_PLUS ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"point_plus"         index_status(yyextra)=1; (*yylval) = _POINT_PLUS ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"square_point"         index_status(yyextra)=1; (*yylval) = _POINT_CARRE  ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"point_carre"         index_status(yyextra)=1; (*yylval) = _POINT_CARRE  ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"cross_point"         index_status(yyextra)=1; (*yylval) = 0  ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"point_croix"         index_status(yyextra)=1; (*yylval) = 0  ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"triangle_point"         index_status(yyextra)=1; (*yylval) = _POINT_TRIANGLE ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"point_triangle"         index_status(yyextra)=1; (*yylval) = _POINT_TRIANGLE ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"star_point"         index_status(yyextra)=1; (*yylval) = _POINT_ETOILE ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"point_etoile"         index_status(yyextra)=1; (*yylval) = _POINT_ETOILE ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"point_point"         index_status(yyextra)=1; (*yylval) = _POINT_POINT ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"invisible_point"         index_status(yyextra)=1; (*yylval) = _POINT_INVISIBLE ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"point_invisible"         index_status(yyextra)=1; (*yylval) = _POINT_INVISIBLE ; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"black"         index_status(yyextra)=1; (*yylval) = _BLACK; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"gomme"         index_status(yyextra)=1; (*yylval) = _WHITE; /* was 49 */ (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"noir"         index_status(yyextra)=1; (*yylval) = _BLACK; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"red"         index_status(yyextra)=1; (*yylval) = _RED; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"rouge"         index_status(yyextra)=1; (*yylval) = _RED; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"green"         index_status(yyextra)=1; (*yylval) = _GREEN; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"vert"         index_status(yyextra)=1; (*yylval) = _GREEN; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"yellow"         index_status(yyextra)=1; (*yylval) = _YELLOW; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"jaune"         index_status(yyextra)=1; (*yylval) = _YELLOW; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"blue"         index_status(yyextra)=1; (*yylval) = _BLUE; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"bleu"         index_status(yyextra)=1; (*yylval) = _BLUE; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"magenta"         index_status(yyextra)=1; (*yylval) = _MAGENTA; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"cyan"         index_status(yyextra)=1; (*yylval) = _CYAN; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"white"         index_status(yyextra)=1; (*yylval) = _WHITE; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"blanc"         index_status(yyextra)=1; (*yylval) = _WHITE; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"solid_line"         index_status(yyextra)=1; (*yylval) = 0; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"ligne_trait_plein"         index_status(yyextra)=1; (*yylval) = 0; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"dash_line"         index_status(yyextra)=1; (*yylval) = _DASH_LINE; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"ligne_tiret"         index_status(yyextra)=1; (*yylval) = _DASH_LINE; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"dot_line"         index_status(yyextra)=1; (*yylval) = _DOT_LINE; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"ligne_point"         index_status(yyextra)=1; (*yylval) = _DOT_LINE; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"dashdot_line"         index_status(yyextra)=1; (*yylval) = _DASHDOT_LINE; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"ligne_tiret_point"         index_status(yyextra)=1; (*yylval) = _DASHDOT_LINE; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"dashdotdot_line"         index_status(yyextra)=1; (*yylval) = _DASHDOTDOT_LINE; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"ligne_tiret_pointpoint"         index_status(yyextra)=1; (*yylval) = _DASHDOTDOT_LINE; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"cap_flat_line"         index_status(yyextra)=1; (*yylval) = _CAP_FLAT_LINE; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"ligne_chapeau_plat"         index_status(yyextra)=1; (*yylval) = _CAP_FLAT_LINE; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"cap_round_line"         index_status(yyextra)=1; (*yylval) = _CAP_ROUND_LINE; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"ligne_chapeau_rond"         index_status(yyextra)=1; (*yylval) = _CAP_ROUND_LINE; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"cap_square_line"         index_status(yyextra)=1; (*yylval) = _CAP_SQUARE_LINE; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"ligne_chapeau_carre"         index_status(yyextra)=1; (*yylval) = _CAP_SQUARE_LINE; (*yylval).subtype=_INT_COLOR ;return T_NUMBER;
"adaptive"         index_status(yyextra)=1; (*yylval) =_ADAPTIVE ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"axes"         index_status(yyextra)=1; (*yylval) =_AXES ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"frames"         index_status(yyextra)=1; (*yylval) =_FRAMES ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"trames"         index_status(yyextra)=1; (*yylval) =_FRAMES ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"nstep"         index_status(yyextra)=1; (*yylval) = _NSTEP ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"tstep"         index_status(yyextra)=1; (*yylval) = _TSTEP ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"ustep"         index_status(yyextra)=1; (*yylval) = _USTEP ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"vstep"         index_status(yyextra)=1; (*yylval) = _VSTEP ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"xstep"         index_status(yyextra)=1; (*yylval) = _XSTEP ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"ystep"         index_status(yyextra)=1; (*yylval) = _YSTEP ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"zstep"         index_status(yyextra)=1; (*yylval) = _ZSTEP ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"font"         index_status(yyextra)=1; (*yylval) =_FONT ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"labels"         index_status(yyextra)=1; (*yylval) = _LABELS; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"linestyle"         index_status(yyextra)=1; (*yylval) = _LINESTYLE; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"resolution"         index_status(yyextra)=1; (*yylval) =_RESOLUTION ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"sample"         index_status(yyextra)=1; (*yylval) = _SAMPLE; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"scaling"         index_status(yyextra)=1; (*yylval) = _SCALING; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"style"         index_status(yyextra)=1; (*yylval) = _STYLE; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"symbol"         index_status(yyextra)=0; (*yylval) = _SYMBOL; (*yylval).subtype=_INT_PLOT ;return T_TYPE_ID;
"symbolsize"         index_status(yyextra)=1; (*yylval) =_SYMBOLSIZE ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"thickness"         index_status(yyextra)=1; (*yylval) = _THICKNESS; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"epaisseur"         index_status(yyextra)=1; (*yylval) = _THICKNESS; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"title"         index_status(yyextra)=1; (*yylval) = _TITLE; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"titre"         index_status(yyextra)=1; (*yylval) = _TITLE; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"titlefont"         index_status(yyextra)=1; (*yylval) = _TITLEFONT; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"view"         index_status(yyextra)=1; (*yylval) = _VIEW; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"axesfont"         index_status(yyextra)=1; (*yylval) =_AXESFONT ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"coords"         index_status(yyextra)=1; (*yylval) = _COORDS; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"labelfont"         index_status(yyextra)=1; (*yylval) = _LABELFONT; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"labeldirections"         index_status(yyextra)=1; (*yylval) = _LABELDIRECTIONS; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"numpoints"         index_status(yyextra)=1; (*yylval) = _NUMPOINTS; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"tickmarks"         index_status(yyextra)=1; (*yylval) =_TICKMARKS ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"xtickmarks"         index_status(yyextra)=1; (*yylval) =_XTICKMARKS ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_texture"         index_status(yyextra)=1; (*yylval) =_GL_TEXTURE ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_light0"         index_status(yyextra)=1; (*yylval) =_GL_LIGHT0 ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_light1"         index_status(yyextra)=1; (*yylval) =_GL_LIGHT1 ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_light2"         index_status(yyextra)=1; (*yylval) =_GL_LIGHT2 ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_light3"         index_status(yyextra)=1; (*yylval) =_GL_LIGHT3 ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_light4"         index_status(yyextra)=1; (*yylval) =_GL_LIGHT4 ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_light5"         index_status(yyextra)=1; (*yylval) =_GL_LIGHT5 ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_light6"         index_status(yyextra)=1; (*yylval) =_GL_LIGHT6 ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_light7"         index_status(yyextra)=1; (*yylval) =_GL_LIGHT7 ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_ambient"       index_status(yyextra)=1; (*yylval) =_GL_AMBIENT ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_specular"       index_status(yyextra)=1; (*yylval) =_GL_SPECULAR ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_diffuse"       index_status(yyextra)=1; (*yylval) =_GL_DIFFUSE ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_position"       index_status(yyextra)=1; (*yylval) =_GL_POSITION ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_spot_direction"       index_status(yyextra)=1; (*yylval) =_GL_SPOT_DIRECTION ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_spot_exponent"       index_status(yyextra)=1; (*yylval) =_GL_SPOT_EXPONENT ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_spot_cutoff"       index_status(yyextra)=1; (*yylval) =_GL_SPOT_CUTOFF ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_constant_attenuation"       index_status(yyextra)=1; (*yylval) =_GL_CONSTANT_ATTENUATION ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_linear_attenuation"       index_status(yyextra)=1; (*yylval) =_GL_LINEAR_ATTENUATION ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_quadratic_attenuation"       index_status(yyextra)=1; (*yylval) =_GL_QUADRATIC_ATTENUATION ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_option"       index_status(yyextra)=1; (*yylval) =_GL_OPTION ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_smooth"       index_status(yyextra)=1; (*yylval) =_GL_SMOOTH ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_flat"       index_status(yyextra)=1; (*yylval) =_GL_FLAT ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_shininess"       index_status(yyextra)=1; (*yylval) =_GL_SHININESS ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_front"       index_status(yyextra)=1; (*yylval) =_GL_FRONT ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_back"       index_status(yyextra)=1; (*yylval) =_GL_BACK ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_front_and_back"       index_status(yyextra)=1; (*yylval) =_GL_FRONT_AND_BACK ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_ambient_and_diffuse"       index_status(yyextra)=1; (*yylval) =_GL_AMBIENT_AND_DIFFUSE ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_emission"       index_status(yyextra)=1; (*yylval) =_GL_EMISSION ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_light_model_ambient"       index_status(yyextra)=1; (*yylval) =_GL_LIGHT_MODEL_AMBIENT ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_light_model_local_viewer"       index_status(yyextra)=1; (*yylval) =_GL_LIGHT_MODEL_LOCAL_VIEWER ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_light_model_two_side"       index_status(yyextra)=1; (*yylval) =_GL_LIGHT_MODEL_TWO_SIDE ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_light_model_color_control"       index_status(yyextra)=1; (*yylval) =_GL_LIGHT_MODEL_COLOR_CONTROL ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_blend"       index_status(yyextra)=1; (*yylval) =_GL_BLEND ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_separate_specular_color"       index_status(yyextra)=1; (*yylval) =_GL_SEPARATE_SPECULAR_COLOR ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_single_color"       index_status(yyextra)=1; (*yylval) =_GL_SINGLE_COLOR ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_src_alpha"       index_status(yyextra)=1; (*yylval) =_GL_SRC_ALPHA ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_one_minus_src_alpha"       index_status(yyextra)=1; (*yylval) =_GL_ONE_MINUS_SRC_ALPHA ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_material"       index_status(yyextra)=1; (*yylval) =_GL_MATERIAL ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_color_indexes"       index_status(yyextra)=1; (*yylval) =_GL_COLOR_INDEXES ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_light"       index_status(yyextra)=1; (*yylval) =_GL_LIGHT ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_perspective"       index_status(yyextra)=1; (*yylval) =_GL_PERSPECTIVE ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_quaternion"       index_status(yyextra)=1; (*yylval) =_GL_QUATERNION ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_rotation_axis"       index_status(yyextra)=1; (*yylval) =_GL_ROTATION_AXIS ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_ortho"       index_status(yyextra)=1; (*yylval) =_GL_ORTHO ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_x"       index_status(yyextra)=1; (*yylval) =_GL_X ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_y"       index_status(yyextra)=1; (*yylval) =_GL_Y ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_z"       index_status(yyextra)=1; (*yylval) =_GL_Z ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_xtick"       index_status(yyextra)=1; (*yylval) =_GL_XTICK ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_ytick"       index_status(yyextra)=1; (*yylval) =_GL_YTICK ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_ztick"       index_status(yyextra)=1; (*yylval) =_GL_ZTICK ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_animate"       index_status(yyextra)=1; (*yylval) =_GL_ANIMATE ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_showaxes"       index_status(yyextra)=1; (*yylval) =_GL_SHOWAXES ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_shownames"       index_status(yyextra)=1; (*yylval) =_GL_SHOWNAMES ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;
"gl_x_axis_name"       index_status(yyextra)=1; (*yylval) =_GL_X_AXIS_NAME ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER; 
"gl_y_axis_name"       index_status(yyextra)=1; (*yylval) =_GL_Y_AXIS_NAME ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER; 
"gl_z_axis_name"       index_status(yyextra)=1; (*yylval) =_GL_Z_AXIS_NAME ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER; 
"gl_x_axis_unit"       index_status(yyextra)=1; (*yylval) =_GL_X_AXIS_UNIT ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER; 
"gl_y_axis_unit"       index_status(yyextra)=1; (*yylval) =_GL_Y_AXIS_UNIT ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER; 
"gl_z_axis_unit"       index_status(yyextra)=1; (*yylval) =_GL_Z_AXIS_UNIT ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER; 
"gl_x_axis_color"       index_status(yyextra)=1; (*yylval) =_GL_X_AXIS_COLOR ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER; 
"gl_y_axis_color"       index_status(yyextra)=1; (*yylval) =_GL_Y_AXIS_COLOR ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER; 
"gl_z_axis_color"       index_status(yyextra)=1; (*yylval) =_GL_Z_AXIS_COLOR ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER; 
  /* "gl_"       index_status(yyextra)=1; (*yylval) =_GL_ ; (*yylval).subtype=_INT_PLOT ;return T_NUMBER;  */
"set"         index_status(yyextra)=0; (*yylval) =_SET__VECT ; (*yylval).subtype=_INT_MAPLECONVERSION ;return T_TYPE_ID;
"list"         if (xcas_mode(yyextra)==3) { index_status(yyextra)=1; return find_or_make_symbol(yytext,(*yylval),yyextra); } index_status(yyextra)=0; (*yylval) = _MAPLE_LIST ; (*yylval).subtype=_INT_MAPLECONVERSION ;return T_TYPE_ID;
"polynom"         index_status(yyextra)=0; (*yylval) =_POLY1__VECT ; (*yylval).subtype=_INT_MAPLECONVERSION ;return T_TYPE_ID;
"trig"         index_status(yyextra)=1; (*yylval) =_TRIG ; (*yylval).subtype=_INT_MAPLECONVERSION ;return T_NUMBER;
"expln"         index_status(yyextra)=1; (*yylval) =_EXPLN ; (*yylval).subtype=_INT_MAPLECONVERSION ;return T_NUMBER;
"parfrac"         index_status(yyextra)=1; (*yylval) =_PARFRAC ; (*yylval).subtype=_INT_MAPLECONVERSION ;return T_NUMBER;
"base"         index_status(yyextra)=1; (*yylval) =_BASE ; (*yylval).subtype=_INT_MAPLECONVERSION ;return T_NUMBER;
"confrac"         index_status(yyextra)=1; (*yylval) =_CONFRAC ; (*yylval).subtype=_INT_MAPLECONVERSION ;return T_NUMBER;
"fullparfrac"         index_status(yyextra)=1; (*yylval) =_FULLPARFRAC ; (*yylval).subtype=_INT_MAPLECONVERSION ;return T_NUMBER;
"Delete"         index_status(yyextra)=1; (*yylval) =_DELETE_OPERATOR ; (*yylval).subtype=_INT_MUPADOPERATOR ;return T_NUMBER;
"Prefix"         index_status(yyextra)=1; (*yylval) =_PREFIX_OPERATOR ; (*yylval).subtype=_INT_MUPADOPERATOR ;return T_NUMBER;
"Postfix"         index_status(yyextra)=1; (*yylval) =_POSTFIX_OPERATOR ; (*yylval).subtype=_INT_MUPADOPERATOR ;return T_NUMBER;
"Binary"         index_status(yyextra)=1; (*yylval) =_BINARY_OPERATOR ; (*yylval).subtype=_INT_MUPADOPERATOR ;return T_NUMBER;
"Nary"         index_status(yyextra)=1; (*yylval) =_NARY_OPERATOR ; (*yylval).subtype=_INT_MUPADOPERATOR ;return T_NUMBER;
"revlex"         index_status(yyextra)=1; (*yylval) =_REVLEX_ORDER ; (*yylval).subtype=_INT_GROEBNER ;return T_NUMBER;
"plex"         index_status(yyextra)=1; (*yylval) =_PLEX_ORDER ; (*yylval).subtype=_INT_GROEBNER ;return T_NUMBER;
"tdeg"         index_status(yyextra)=1; (*yylval) =_TDEG_ORDER ; (*yylval).subtype=_INT_GROEBNER ;return T_NUMBER;
"with_cocoa"         index_status(yyextra)=1; (*yylval) =_WITH_COCOA ; (*yylval).subtype=_INT_GROEBNER ;return T_NUMBER;
"with_f5"         index_status(yyextra)=1; (*yylval) =_WITH_F5 ; (*yylval).subtype=_INT_GROEBNER ;return T_NUMBER;


    /* vector/polynom/matrice delimiters */
"seq["              (*yylval) = _SEQ__VECT; return T_VECT_DISPATCH;
"set["              (*yylval) = _SET__VECT; return T_VECT_DISPATCH;
"rpn_func["         (*yylval) = _RPN_FUNC__VECT; return T_VECT_DISPATCH;
"group["            (*yylval) = _GROUP__VECT; return T_VECT_DISPATCH;
"line["             (*yylval) = _LINE__VECT; return T_VECT_DISPATCH;
"vector["           (*yylval) = _VECTOR__VECT; return T_VECT_DISPATCH;
"matrix["           (*yylval) = _MATRIX__VECT; return T_VECT_DISPATCH;
"pnt["              (*yylval) = _PNT__VECT; return T_VECT_DISPATCH;
"point["            (*yylval) = _POINT__VECT; return T_VECT_DISPATCH;
"curve["            (*yylval) = _CURVE__VECT; return T_VECT_DISPATCH;
"halfline["         (*yylval) = _HALFLINE__VECT; return T_VECT_DISPATCH;
"poly1["            (*yylval) = _POLY1__VECT; return T_VECT_DISPATCH;
"assume["           (*yylval) = _ASSUME__VECT; return T_VECT_DISPATCH;
"spreadsheet["      (*yylval) = _SPREAD__VECT; return T_VECT_DISPATCH;
"folder["      (*yylval) = _FOLDER__VECT; return T_VECT_DISPATCH;
"polyedre["      (*yylval) = _POLYEDRE__VECT; return T_VECT_DISPATCH;
"rgba["      (*yylval) = _RGBA__VECT; return T_VECT_DISPATCH;
"<"                     index_status(yyextra)=0; (*yylval)=gen(at_inferieur_strict,2);  return T_TEST_EQUAL;
">"                     index_status(yyextra)=0; (*yylval)=gen(at_superieur_strict,2); return T_TEST_EQUAL;
","                     index_status(yyextra)=0; return T_VIRGULE;
",,"                     index_status(yyextra)=0; return T_VIRGULE;
"("                     index_status(yyextra)=0; return T_BEGIN_PAR;
")"                     index_status(yyextra)=1; return T_END_PAR;
\[			if (index_status(yyextra)) { index_status(yyextra)=0; return T_INDEX_BEGIN; } else { (*yylval) = 0; return T_VECT_DISPATCH; } ;
\]			index_status(yyextra)=1; return T_VECT_END;
"%["                    index_status(yyextra)=0; return T_POLY1_BEGIN;
"%]"                    index_status(yyextra)=1; return T_POLY1_END;
"%%["                   index_status(yyextra)=0; return T_MATRICE_BEGIN;
"%%]"                   index_status(yyextra)=1; return T_MATRICE_END;
"%%%["                  index_status(yyextra)=0; return T_ASSUME_BEGIN;
"%%%]"                  index_status(yyextra)=1; return T_ASSUME_END;
    /* geometric delimiters */
"%("                    index_status(yyextra)=0; return T_GROUPE_BEGIN;
"%)"                    index_status(yyextra)=1; return T_GROUPE_END;
"%%("                   index_status(yyextra)=0; return T_LINE_BEGIN;
"%%)"                   index_status(yyextra)=1; return T_LINE_END;
"%%%("                  index_status(yyextra)=0; return T_VECTOR_BEGIN;
"%%%)"                  index_status(yyextra)=1; return T_VECTOR_END;
"%%%%("                 index_status(yyextra)=0; return T_CURVE_BEGIN;
"%%%%)"                 index_status(yyextra)=1; return T_CURVE_END;
    /* gen delimiters */
"{"                     index_status(yyextra)=0; if (rpn_mode) { (*yylval)=0; return T_VECT_DISPATCH; } if (xcas_mode(yyextra)==3) return T_VECT_BEGIN; if (xcas_mode(yyextra) > 0 ) return T_SET_BEGIN; else return T_BLOC_BEGIN;
"}"                     index_status(yyextra)=1; if (rpn_mode) return T_VECT_END; if (xcas_mode(yyextra)==3) return T_VECT_END; if (xcas_mode(yyextra) > 0) return T_SET_END; else return T_BLOC_END;
"%{"                    index_status(yyextra)=0; return T_SET_BEGIN;
"%}"                    index_status(yyextra)=1; return T_SET_END;
"%%{"                   index_status(yyextra)=0; return T_ROOTOF_BEGIN;
"%%}"                   index_status(yyextra)=1; return T_ROOTOF_END;
"%%%{"                  index_status(yyextra)=0; return T_SPOLY1_BEGIN;
"%%%}"                  index_status(yyextra)=1; return T_SPOLY1_END;
"<<"                    index_status(yyextra)=0; ++in_rpn(yyextra); return T_RPN_BEGIN;
">>"                    index_status(yyextra)=0; --in_rpn(yyextra); return T_RPN_END;

    /* Maple libraries names */
"linalg"                index_status(yyextra)=1; (*yylval)=_LINALG; (*yylval).subtype=_INT_MAPLELIB; return T_MAPLELIB;
"numtheory"                index_status(yyextra)=1; (*yylval)=_NUMTHEORY; (*yylval).subtype=_INT_MAPLELIB; return T_MAPLELIB;
"groebner"                index_status(yyextra)=1; (*yylval)=_GROEBNER; (*yylval).subtype=_INT_MAPLELIB; return T_MAPLELIB;

    /* binary operators */
"->"                    index_status(yyextra)=0; return T_MAPSTO;
"-<"                    (*yylval) = gen(at_couleur,2); index_status(yyextra)=0; return T_INTERVAL;
"=="			index_status(yyextra)=0; (*yylval)=gen(at_same,2); return T_TEST_EQUAL;
"'=='"                  index_status(yyextra)=0; (*yylval)=gen(at_same,2); return T_QUOTED_BINARY;
"_equal"                  index_status(yyextra)=0; (*yylval)=gen(at_same,2); return T_QUOTED_BINARY;
"!="			index_status(yyextra)=0; (*yylval)=gen(at_different,2); return T_TEST_EQUAL;
"'!='"                  index_status(yyextra)=0; (*yylval)=gen(at_different,2); return T_QUOTED_BINARY;
"<>"			index_status(yyextra)=0; (*yylval)=gen(at_different,2); return T_TEST_EQUAL;
"'<>'"                  index_status(yyextra)=0; (*yylval)=gen(at_different,2); return T_QUOTED_BINARY;
"_unequal"                  index_status(yyextra)=0; (*yylval)=gen(at_different,2); return T_QUOTED_BINARY;
"'<='"                  index_status(yyextra)=0; (*yylval)=gen(at_inferieur_egal,2); return T_QUOTED_BINARY;
"_leequal"                  index_status(yyextra)=0; (*yylval)=gen(at_inferieur_egal,2); return T_QUOTED_BINARY;
"<="			index_status(yyextra)=0; (*yylval)=gen(at_inferieur_egal,2); return T_TEST_EQUAL;
"'<'"                  index_status(yyextra)=0; (*yylval)=gen(at_inferieur_strict,2); return T_QUOTED_BINARY;
"_less"                  index_status(yyextra)=0; (*yylval)=gen(at_inferieur_strict,2); return T_QUOTED_BINARY;
"'>'"                  index_status(yyextra)=0; (*yylval)=gen(at_superieur_strict,2); return T_QUOTED_BINARY;
">="			index_status(yyextra)=0; (*yylval)=gen(at_superieur_egal,2); return T_TEST_EQUAL;
"'>='"                  index_status(yyextra)=0; (*yylval)=gen(at_superieur_egal,2); return T_QUOTED_BINARY;
"="                     spread_formula(yyextra)=!index_status(yyextra); index_status(yyextra)=0; (*yylval)=gen(at_equal,2); return T_EQUAL;
"'='"                  index_status(yyextra)=0; (*yylval)=gen(at_equal,2); return T_QUOTED_BINARY;
"$"                     index_status(yyextra)=0; (*yylval)=gen(at_dollar,2); if (xcas_mode(yyextra)>0) return T_DOLLAR_MAPLE; else return T_DOLLAR;
"%$"                   index_status(yyextra)=0; (*yylval)=gen(at_dollar,2); return T_DOLLAR_MAPLE;
"'$'"                  index_status(yyextra)=0; (*yylval)=gen(at_dollar,2); return T_QUOTED_BINARY;
"_seqgen"                  index_status(yyextra)=0; (*yylval)=gen(at_dollar,2); return T_QUOTED_BINARY;
":="			index_status(yyextra)=0; (*yylval)=gen(at_sto,2); return T_AFFECT;
"':='"                  index_status(yyextra)=0; (*yylval)=gen(at_sto,2); return T_QUOTED_BINARY;
"_assign"                  index_status(yyextra)=0; (*yylval)=gen(at_sto,2); return T_QUOTED_BINARY;
"=>"                    index_status(yyextra)=0; (*yylval)=gen(at_sto,2); return TI_STO;
"=<"                    index_status(yyextra)=0; (*yylval)=gen(at_array_sto,2); return T_AFFECT;
"@"{D}+                   index_status(yyextra)=1; yytext[0]='0'; (*yylval) = symb_double_deux_points(makevecteur(_IDNT_break,chartab2gen(yytext,yyextra))); return T_SYMBOL;
"@"                     if (xcas_mode(yyextra)!=3) {index_status(yyextra)=0; (*yylval)=gen(at_compose,2); return T_COMPOSE; } BEGIN(comment_hash);
"@@"                     index_status(yyextra)=0; (*yylval)=gen(at_composepow,2); return T_POW;
"'@@'"                     index_status(yyextra)=0; (*yylval)=gen(at_composepow,2); return T_QUOTED_BINARY;
"_fnest"                     index_status(yyextra)=0; (*yylval)=gen(at_composepow,2); return T_QUOTED_BINARY;
"'@'"                  index_status(yyextra)=0; (*yylval)=gen(at_compose,2); return T_QUOTED_BINARY;
"_fconcat"                  index_status(yyextra)=0; (*yylval)=gen(at_compose,2); return T_QUOTED_BINARY;
"&&"                    index_status(yyextra)=0; (*yylval)=gen(at_and,2); return T_AND_OP;
"and"                   index_status(yyextra)=0; (*yylval)=gen(at_and,2); return T_AND_OP;
"AND"                   index_status(yyextra)=0; (*yylval)=gen(at_and,2); return T_AND_OP;
"'&&'"                  index_status(yyextra)=0; (*yylval)=gen(at_and,2); return T_QUOTED_BINARY;
"'and'"                 index_status(yyextra)=0; (*yylval)=gen(at_and,2); return T_QUOTED_BINARY;
"_and"                 index_status(yyextra)=0; (*yylval)=gen(at_and,2); return T_QUOTED_BINARY;
"|"                     index_status(yyextra)=0; (*yylval)=gen(at_tilocal,2); return T_PIPE;
"||"                    index_status(yyextra)=0; (*yylval)=gen(at_ou,2); return T_AND_OP;
"'||'"                  index_status(yyextra)=0; (*yylval)=gen(at_ou,2); return T_QUOTED_BINARY;
"or"                    index_status(yyextra)=0; (*yylval)=gen(at_ou,2); return T_AND_OP;
"'or'"                  index_status(yyextra)=0; (*yylval)=gen(at_ou,2); return T_QUOTED_BINARY;
"_or"                  index_status(yyextra)=0; (*yylval)=gen(at_ou,2); return T_QUOTED_BINARY;
"OR"                    index_status(yyextra)=0; (*yylval)=gen(at_ou,2); return T_AND_OP;
"xor"                    index_status(yyextra)=0; (*yylval)=gen(at_xor,2); return T_AND_OP;
"_xor"                  index_status(yyextra)=0; (*yylval)=gen(at_xor,2); return T_QUOTED_BINARY;
"'xor'"                  index_status(yyextra)=0; (*yylval)=gen(at_xor,2); return T_QUOTED_BINARY;
"XOR"                    index_status(yyextra)=0; (*yylval)=gen(at_xor,2); return T_AND_OP;
".."                    index_status(yyextra)=0; (*yylval)=gen(at_interval,2); return T_INTERVAL;
"'..'"                  index_status(yyextra)=0; (*yylval)=gen(at_interval,2); return T_QUOTED_BINARY;
"_range"                  index_status(yyextra)=0; (*yylval)=gen(at_interval,2); return T_QUOTED_BINARY;
"!"                     if (xcas_mode(yyextra) || index_status(yyextra)) { (*yylval)=gen(at_factorial); return T_FACTORIAL; } else { index_status(yyextra)=0; (*yylval)=gen(at_not,1); return T_NOT; }

    /* standard functions */
"+"                     index_status(yyextra)=0; (*yylval)=gen(at_plus,2); return T_PLUS;
"++"                    index_status(yyextra)=0; (*yylval)=gen(at_increment,1); return T_FACTORIAL;
"+="                    index_status(yyextra)=0; (*yylval)=gen(at_increment,1); return T_PLUS;
"--"                    index_status(yyextra)=0; (*yylval)=gen(at_decrement,1); return T_FACTORIAL;
"-="                    index_status(yyextra)=0; (*yylval)=gen(at_decrement,1); return T_PLUS;
".+"                    index_status(yyextra)=0; (*yylval)=gen(at_plus,2); return T_PLUS;
"&"                     index_status(yyextra)=0; (*yylval)=gen(at_plus,2); return T_PLUS;
"≤"                     index_status(yyextra)=0; (*yylval)=2; return T_SQ;
"≥"                     index_status(yyextra)=0; (*yylval)=3; return T_SQ;
  /* "','"                   index_status(yyextra)=0; (*yylval)=gen(at_makevector,2); return T_QUOTED_BINARY; commented because of f('a','b') */
"'+'"                   index_status(yyextra)=0; (*yylval)=gen(at_plus,2); return T_QUOTED_BINARY;
"_plus"                   index_status(yyextra)=0; (*yylval)=gen(at_plus,2); return T_QUOTED_BINARY;
"-"                     index_status(yyextra)=0; (*yylval)=gen(at_binary_minus,2); return T_MOINS;
".-"                     index_status(yyextra)=0; (*yylval)=gen(at_binary_minus,2); return T_MOINS;
"'-'"                   index_status(yyextra)=0; (*yylval)=gen(at_binary_minus,2); return T_QUOTED_BINARY;
"_subtract"                   index_status(yyextra)=0; (*yylval)=gen(at_binary_minus,2); return T_QUOTED_BINARY;
"*"                     index_status(yyextra)=0; (*yylval)=gen(at_prod,2); return T_FOIS;
"*="                    index_status(yyextra)=0; (*yylval)=gen(at_multcrement,1); return T_FOIS;
"."                     index_status(yyextra)=0; (*yylval)=gen(at_prod,2); return T_FOIS;
"&*"                     index_status(yyextra)=0; (*yylval)=gen(at_ampersand_times,2); return T_FOIS;
"&^"                     index_status(yyextra)=0; (*yylval)=gen(at_quote_pow,2); return T_POW;
".*"                     index_status(yyextra)=0; (*yylval)=gen(at_pointprod,2); return T_FOIS;
"'*'"                   index_status(yyextra)=0; (*yylval)=gen(at_prod,2); return T_QUOTED_BINARY;
"_mult"                   index_status(yyextra)=0; (*yylval)=gen(at_prod,2); return T_QUOTED_BINARY;
"/"                     index_status(yyextra)=0; (*yylval)=gen(at_division,2); return T_DIV;
"/="                    index_status(yyextra)=0; (*yylval)=gen(at_divcrement,1); return T_DIV;
"./"                     index_status(yyextra)=0; (*yylval)=gen(at_pointdivision,2); return T_DIV;
"'/'"                   index_status(yyextra)=0; (*yylval)=gen(at_division,2); return T_QUOTED_BINARY;
"_divide"                   index_status(yyextra)=0; (*yylval)=gen(at_division,2); return T_QUOTED_BINARY;
"%"                     index_status(yyextra)=0; if (xcas_mode(yyextra)==3) { (*yylval)=gen(at_pourcent); return T_FACTORIAL; } if (xcas_mode(yyextra)==1) { (*yylval)=symbolic(at_ans,vecteur(0)); return T_NUMBER; }  if (xcas_mode(yyextra)) (*yylval)=gen(at_irem,2); else (*yylval)=0; return T_MOD;
"'%'"                  index_status(yyextra)=0; (*yylval)=gen(at_irem,2); return T_QUOTED_BINARY;
"mod"                   index_status(yyextra)=0; if (xcas_mode(yyextra)==3) { (*yylval)=gen(at_irem,2); return T_UNARY_OP; } else { if (xcas_mode(yyextra)) (*yylval)=gen(at_irem,2); else (*yylval)=0; return T_MOD; }
"'mod'"                  index_status(yyextra)=0; (*yylval)=gen(at_irem,2); return T_QUOTED_BINARY;
"_mod"                  index_status(yyextra)=0; (*yylval)=gen(at_irem,2); return T_QUOTED_BINARY;
"MOD"                   index_status(yyextra)=0; return T_MOD;
"^"                     index_status(yyextra)=0; (*yylval)=gen(at_pow,2); return T_POW;
"pow"		         (*yylval) = gen(at_pow,2); index_status(yyextra)=0; return T_UNARY_OP;
"**"                     index_status(yyextra)=0; (*yylval)=gen(at_pow,2); return T_POW;
".^"                     index_status(yyextra)=0; (*yylval)=gen(at_pointpow,2); return T_POW;
"'^'"                   index_status(yyextra)=0; (*yylval)=gen(at_pow,2); return T_QUOTED_BINARY;
"_power"                   index_status(yyextra)=0; (*yylval)=gen(at_pow,2); return T_QUOTED_BINARY;
"SWAP"                  (*yylval) = gen(at_SWAP,0); index_status(yyextra)=0; return T_RPN_OP;
"DROP"                  (*yylval) = gen(at_DROP,0); index_status(yyextra)=0; return T_RPN_OP;
"DUP"                   (*yylval) = gen(at_DUP,0); index_status(yyextra)=0; return T_RPN_OP;
"ROLL"                  (*yylval) = gen(at_ROLL,0); index_status(yyextra)=0; return T_RPN_OP;
"ROLLD"                 (*yylval) = gen(at_ROLLD,0); index_status(yyextra)=0; return T_RPN_OP;
"PICK"                  (*yylval) = gen(at_PICK,0); index_status(yyextra)=0; return T_RPN_OP;
"Digits"                  (*yylval) = gen(at_Digits,0); index_status(yyextra)=0; return T_DIGITS;
"DIGITS"                  (*yylval) = gen(at_Digits,0); index_status(yyextra)=0; return T_DIGITS;
"arccos"		(*yylval) = gen(at_acos,1); index_status(yyextra)=0; return T_UNARY_OP;
"randnorm"		(*yylval) = gen(at_randNorm,1); index_status(yyextra)=0; return T_UNARY_OP;
"arccosh"		(*yylval) = gen(at_acosh,1); index_status(yyextra)=0; return T_UNARY_OP;
"alg"                   (*yylval) = gen(at_alg,0); index_status(yyextra)=0; return T_UNARY_OP;
"angle_radian"		(*yylval) = gen(at_angle_radian,0); index_status(yyextra)=0; return T_DIGITS;
"approx_mode"		(*yylval) = gen(at_approx_mode,0); index_status(yyextra)=0; return T_DIGITS;
"args"			index_status(yyextra)=1; return T_ARGS;
"'args'"                index_status(yyextra)=0; (*yylval)=gen(at_args,0); return T_QUOTED_BINARY;
"posons"		(*yylval) = gen(at_assume,1); index_status(yyextra)=0; return T_UNARY_OP;
"arcsin"		(*yylval) = gen(at_asin,1); index_status(yyextra)=0; return T_UNARY_OP;
"arcsinh"		(*yylval) = gen(at_asinh,1); index_status(yyextra)=0; return T_UNARY_OP;
"at"			(*yylval) = gen(at_at,2); index_status(yyextra)=0; return T_UNARY_OP;
"arctan"		(*yylval) = gen(at_atan,1); index_status(yyextra)=0; return T_UNARY_OP;
"arctanh"		(*yylval) = gen(at_atanh,1); index_status(yyextra)=0; return T_UNARY_OP;
"backquote"		(*yylval) = gen(at_backquote,1); index_status(yyextra)=0; return T_UNARY_OP;
"bloc"  		(*yylval) = gen(at_bloc,1); index_status(yyextra)=0; return T_UNARY_OP;
"begin"                 index_status(yyextra)=0; return T_BLOC_BEGIN;
"BEGIN"                 index_status(yyextra)=0; return T_BLOC_BEGIN;
"break"  		index_status(yyextra)=0; (*yylval)=gen(at_break,0); return T_BREAK;
"BREAK"  		index_status(yyextra)=0; (*yylval)=gen(at_break,0); return T_BREAK;
"by"                    index_status(yyextra)=0; return T_BY;
"case"                  index_status(yyextra)=0; return T_CASE;
"CASE"                  index_status(yyextra)=0; return T_CASE;
"catch"                 index_status(yyextra)=0; return T_CATCH;
"continue"  		index_status(yyextra)=0; return T_CONTINUE;
"click"		        (*yylval) = gen(at_click,1); index_status(yyextra)=0; return T_UNARY_OP;
"comment"		(*yylval) = gen(at_comment,1); index_status(yyextra)=0; return T_UNARY_OP;
"all_trig_solutions"		(*yylval) = gen(at_all_trig_solutions,1); index_status(yyextra)=0; return T_DIGITS;
"ntl_on"		(*yylval) = gen(at_ntl_on,1); index_status(yyextra)=0; return T_DIGITS;
"complex_mode"		(*yylval) = gen(at_complex_mode,1); index_status(yyextra)=0; return T_DIGITS;
"complex_variables"	(*yylval) = gen(at_complex_variables,0); index_status(yyextra)=0; return T_DIGITS;
"CONT"		        (*yylval) = gen(at_cont,1); index_status(yyextra)=0; return T_UNARY_OP;
"DEBUG"		        (*yylval) = gen(at_debug,1); index_status(yyextra)=0; return T_UNARY_OP;
"default"               index_status(yyextra)=1; return T_DEFAULT;
"derive"		(*yylval) = gen(at_derive,2); index_status(yyextra)=0; return T_UNARY_OP;
"do"			index_status(yyextra)=0; (*yylval)=gen(at_for,4); return T_DO;
"DO"			index_status(yyextra)=0; (*yylval)=gen(at_for,4); return T_DO;
"D"  		if (xcas_mode(yyextra)==1 || xcas_mode(yyextra)==2) { (*yylval) = gen(at_function_diff,1); index_status(yyextra)=1; return T_UNARY_OP;} else { index_status(yyextra)=1; return find_or_make_symbol(yytext,(*yylval),yyextra); }
"e"                     if (xcas_mode(yyextra)==1 || xcas_mode(yyextra)==2) { (*yylval)=e__IDNT_e; }else (*yylval)=symbolic(at_exp,1); return T_NUMBER;
"else"			index_status(yyextra)=0; (*yylval)=gen(at_ifte,3); return T_ELSE;
"ELSE"			index_status(yyextra)=0; (*yylval)=gen(at_ifte,3); return T_ELSE;
"elif"                  index_status(yyextra)=0; return T_ELIF;
"end"			index_status(yyextra)=0; return T_BLOC_END;
"END"			index_status(yyextra)=0; return T_BLOC_END;
"end_for"               index_status(yyextra)=0; return T_BLOC_END;
"end_if"                index_status(yyextra)=0; (*yylval)=gen(at_ifte,3); return T_BLOC_END;
"end if"                index_status(yyextra)=0; (*yylval)=gen(at_ifte,3); return T_BLOC_END;
"end do"                index_status(yyextra)=0; (*yylval)=gen(at_ifte,3); return T_BLOC_END;
"end proc"                index_status(yyextra)=0; (*yylval)=gen(at_ifte,3); return T_BLOC_END;
"end_while"             index_status(yyextra)=0; return T_BLOC_END;
"end_case"              index_status(yyextra)=0; return T_ENDCASE;
"end_proc"              index_status(yyextra)=0; return T_BLOC_END;
"epsilon"               (*yylval) = gen(at_epsilon,0); index_status(yyextra)=0; return T_DIGITS;
"proba_epsilon"               (*yylval) = gen(at_proba_epsilon,0); index_status(yyextra)=0; return T_DIGITS;
"equal"			(*yylval) = gen(at_equal,2); index_status(yyextra)=0; return T_UNARY_OP;
"error"		        index_status(yyextra)=0; (*yylval)=gen(at_throw,1); return T_RETURN;
"erase"                 (*yylval) = gen(at_erase,0); index_status(yyextra)=0; return T_UNARY_OP;
"ERROR"		        index_status(yyextra)=0; (*yylval)=gen(at_throw,1); return T_RETURN;
"esac"                  index_status(yyextra)=0; return T_ENDCASE;
"expand"                if (xcas_mode(yyextra)==3) (*yylval)=gen(at_partfrac); else (*yylval) = gen(at_expand,1); index_status(yyextra)=0; return T_UNARY_OP;
"export"		(*yylval) = gen(at_insmod,1); index_status(yyextra)=0; return T_UNARY_OP;
"false"                  index_status(yyextra)=0; (*yylval) = zero; (*yylval).subtype=_INT_BOOLEAN; return T_NUMBER;
"faux"                  index_status(yyextra)=0; (*yylval) = zero; (*yylval).subtype=_INT_BOOLEAN; return T_NUMBER;
"FALSE"                  index_status(yyextra)=0; (*yylval) = zero; (*yylval).subtype=_INT_BOOLEAN; return T_NUMBER;
"fdistrib"		(*yylval) = gen(at_expand,1); index_status(yyextra)=0; return T_UNARY_OP;
"for"			index_status(yyextra)=0; (*yylval)=gen(at_for,4); return T_FOR;
"FOR"			index_status(yyextra)=0; (*yylval)=gen(at_for,4); return T_FOR;
"from"			index_status(yyextra)=0; return T_FROM;
"fi"                    index_status(yyextra)=0; return T_BLOC_END;
"global"                (*yylval)=1; index_status(yyextra)=0; return T_LOCAL;
"HALT"		        (*yylval) = gen(at_halt,1); index_status(yyextra)=0; return T_UNARY_OP;
"if"      		index_status(yyextra)=0; (*yylval)=gen(at_ifte,3); return T_IF;
"IF"      		index_status(yyextra)=0; (*yylval)=gen(at_ifte,3); return T_RPN_IF;
"ifte"      		index_status(yyextra)=0; (*yylval)=gen(at_ifte,3); return T_IFTE;
"IFTE"      		index_status(yyextra)=0; (*yylval)=gen(at_ifte,3); return T_IFTE;
"'ifte'"                index_status(yyextra)=0; (*yylval)=gen(at_ifte,3); return T_QUOTED_BINARY;
"'if'"                index_status(yyextra)=0; (*yylval)=gen(at_ifte,3); return T_QUOTED_BINARY;
"_if"                index_status(yyextra)=0; (*yylval)=gen(at_ifte,3); return T_QUOTED_BINARY;
"ifactors"      	if (xcas_mode(yyextra)==1) (*yylval) = gen(at_maple_ifactors); else (*yylval) = gen(at_ifactors,1); index_status(yyextra)=0; return T_UNARY_OP;
"in"                    index_status(yyextra)=0; return T_IN;
"INTVX"			(*yylval) = gen(at_integrate,1); index_status(yyextra)=0; return T_UNARY_OP;
"inter"			(*yylval) = gen(at_inter,2); index_status(yyextra)=0; return T_UNARY_OP;
"intersect"             index_status(yyextra)=0; return T_INTERSECT;
"'intersect'"                  index_status(yyextra)=0; (*yylval)=gen(at_intersect,2); return T_QUOTED_BINARY;
"_intersect"                  index_status(yyextra)=0; (*yylval)=gen(at_intersect,2); return T_QUOTED_BINARY;
"isom"		        (*yylval) = gen(at_isom,1); index_status(yyextra)=0; return T_UNARY_OP;
"KILL"		        (*yylval) = gen(at_kill,1); index_status(yyextra)=0; return T_UNARY_OP;
"local"			(*yylval)=0; index_status(yyextra)=0; return T_LOCAL;
"localbloc"		index_status(yyextra)=0; return T_LOCALBLOC;
"LOCAL"			(*yylval)=0; index_status(yyextra)=0; return T_LOCAL;
"log"			(*yylval) = gen(at_ln,1); index_status(yyextra)=1; return T_UNARY_OP; /* index_status(yyextra)=1 to accept log[] for a basis log */
"mkisom"		(*yylval) = gen(at_mkisom,2); index_status(yyextra)=0; return T_UNARY_OP;
"minus"                 index_status(yyextra)=0; return T_MINUS;
"'minus'"                  index_status(yyextra)=0; (*yylval)=gen(at_minus,2); return T_QUOTED_BINARY;
"_minus"                  index_status(yyextra)=0; (*yylval)=gen(at_minus,2); return T_QUOTED_BINARY;
"next"		        index_status(yyextra)=0; return T_CONTINUE;
"NEXT"		        index_status(yyextra)=0; return T_CONTINUE;
"not"                 (*yylval) = gen(at_not,1); if (xcas_mode(yyextra)) return T_NOT;  index_status(yyextra)=0; return T_UNARY_OP;
"neg"		(*yylval) = gen(at_neg,1); index_status(yyextra)=0; return T_UNARY_OP;
"'not'"                  index_status(yyextra)=0; (*yylval)=gen(at_not,1); return T_QUOTED_BINARY;
"_not"                  index_status(yyextra)=0; (*yylval)=gen(at_not,1); return T_QUOTED_BINARY;
"NULL"                  index_status(yyextra)=1; return T_NULL;
"normalf"		(*yylval) = gen(at_greduce,1); index_status(yyextra)=0; return T_UNARY_OP;
"od"                    index_status(yyextra)=0; return T_BLOC_END;
"of"			index_status(yyextra)=0; return T_OF;
"'of'"                  index_status(yyextra)=0; (*yylval)=gen(at_of,2); return T_QUOTED_BINARY;
"op"                    if (xcas_mode(yyextra)==1) (*yylval) = gen(at_maple_op,1); else (*yylval) = gen(at_feuille,1); index_status(yyextra)=0; return T_UNARY_OP;
"feuille"               (*yylval) = gen(at_feuille,1); index_status(yyextra)=0; return T_UNARY_OP;
"option"                (*yylval)=2; index_status(yyextra)=0; return T_LOCAL;
"otherwise"             index_status(yyextra)=0; return T_DEFAULT;
"pcoef"		        (*yylval) = gen(at_pcoeff,1); index_status(yyextra)=0; return T_UNARY_OP;
"plotfunc2d" 		(*yylval) = gen(at_funcplot,2); index_status(yyextra)=0; return T_UNARY_OP;
"user_operator" 		(*yylval) = gen(at_user_operator,6); index_status(yyextra)=0; return T_UNARY_OP;
"Ox" 		        (*yylval) = _droite(makevecteur(zero,plus_one),0); return T_EXPRESSION;
"Oy" 		        (*yylval) = _droite(makevecteur(zero,cst_i),0); return T_EXPRESSION;
"program"		index_status(yyextra)=0; return T_PROGRAM;
"proc"                  index_status(yyextra)=0; return T_PROC;
"purge"                 if (rpn_mode) {(*yylval)=gen(at_purge,0); index_status(yyextra)=0; return T_RPN_OP;} else {(*yylval) = gen(at_purge,1); index_status(yyextra)=0; return T_UNARY_OP;};
"unassign"                 if (rpn_mode) {(*yylval)=gen(at_purge,0); index_status(yyextra)=0; return T_RPN_OP;} else {(*yylval) = gen(at_purge,1); index_status(yyextra)=0; return T_UNARY_OP;};
"PURGE"                 if (rpn_mode) {(*yylval)=gen(at_purge,0); index_status(yyextra)=0; return T_RPN_OP;} else {(*yylval) = gen(at_purge,1); index_status(yyextra)=0; return T_UNARY_OP;};
"randseed"			(*yylval) = gen(at_srand,1); index_status(yyextra)=0; return T_RETURN;
"rcl"			(*yylval) = gen(at_RCL,1); index_status(yyextra)=0; return T_UNARY_OP;
"RCL"			(*yylval) = gen(at_RCL,1); index_status(yyextra)=0; return T_UNARY_OP;
"read"			(*yylval) = gen(at_read,1); index_status(yyextra)=0; return T_RETURN;
"repeat"                (*yylval) = gen(at_for,1) ; index_status(yyextra)=0; return T_REPEAT;
"repeter"                (*yylval) = gen(at_for,1) ; index_status(yyextra)=0; return T_REPEAT;
"REPEAT"                (*yylval) = gen(at_for,1) ;index_status(yyextra)=0; return T_REPEAT;
"return"		(*yylval) = gen(at_return,1) ; index_status(yyextra)=0; return T_RETURN;
"retourne"		(*yylval) = gen(at_return,1) ; index_status(yyextra)=0; return T_RETURN;
"RETURN"		(*yylval) = gen(at_return,1) ; index_status(yyextra)=0; return T_RETURN;
"'return'"              (*yylval) = gen(at_return,1) ; index_status(yyextra)=0; return T_QUOTED_BINARY;
"risch"		(*yylval) = gen(at_risch,4); index_status(yyextra)=0; return T_UNARY_OP;
"root"                (*yylval) = gen(at_maple_root,1); index_status(yyextra)=1; return T_UNARY_OP;
"rpn"                   (*yylval) = gen(at_rpn,0); index_status(yyextra)=0; return T_UNARY_OP;
"same"                  (*yylval) = gen(at_same,1); index_status(yyextra)=0; return T_UNARY_OP;
"signal"                (*yylval) = gen(at_signal,1); index_status(yyextra)=0; return T_UNARY_OP;
"SST"		        (*yylval) = gen(at_sst,1); index_status(yyextra)=0; return T_UNARY_OP;
"SST_IN"		(*yylval) = gen(at_sst_in,1); index_status(yyextra)=0; return T_UNARY_OP;
"srand"			(*yylval) = gen(at_srand,1); index_status(yyextra)=0; return T_RETURN;
"START"                 index_status(yyextra)=0; return T_START;
"stack"                 index_status(yyextra)=0; return T_STACK;
"step"                  index_status(yyextra)=0; return T_BY;
"STEP"                  index_status(yyextra)=0; return T_BY;
"STO"			(*yylval) = gen(at_sto,0); index_status(yyextra)=0; return T_UNARY_OP;
"subs"			if (xcas_mode(yyextra)==1) (*yylval) = gen(at_maple_subs,2); else (*yylval) = gen(at_subs,2); index_status(yyextra)=0; return T_UNARY_OP;
"subsop"		if (xcas_mode(yyextra)==1) (*yylval) = gen(at_maple_subsop,2); else (*yylval) = gen(at_subsop,2); index_status(yyextra)=0; return T_UNARY_OP;
"switch"                index_status(yyextra)=0; return T_SWITCH;
"SWITCH"                index_status(yyextra)=0; return T_SWITCH;
"tanh"			(*yylval) = gen(at_tanh,1); index_status(yyextra)=0; return T_UNARY_OP;
"then"		        index_status(yyextra)=0; (*yylval)=gen(at_ifte,3); return T_THEN;
"THEN"		        index_status(yyextra)=0; (*yylval)=gen(at_ifte,3); return T_THEN;
"throw"		        index_status(yyextra)=0; (*yylval)=gen(at_throw,1); return T_RETURN;
"to"                    (*yylval)=1; return T_TO;
"TO"                    (*yylval)=1; return T_TO;
"downto"                (*yylval)=-1; return T_TO;
"DOWNTO"                (*yylval)=-1; return T_TO;
"true"                  index_status(yyextra)=0; (*yylval) = plus_one; (*yylval).subtype=_INT_BOOLEAN; return T_NUMBER;
"vrai"                  index_status(yyextra)=0; (*yylval) = plus_one; (*yylval).subtype=_INT_BOOLEAN; return T_NUMBER;
"TRUE"                  index_status(yyextra)=0; (*yylval) = plus_one; (*yylval).subtype=_INT_BOOLEAN; return T_NUMBER;
"try"                   index_status(yyextra)=0; return T_TRY;
"try_catch"             index_status(yyextra)=0; return T_TRY_CATCH;
"union"                 index_status(yyextra)=0; return T_UNION;
"'union'"                  index_status(yyextra)=0; (*yylval)=gen(at_union,2); return T_QUOTED_BINARY;
"_union"                  index_status(yyextra)=0; (*yylval)=gen(at_union,2); return T_QUOTED_BINARY;
"until"                 index_status(yyextra)=0; return T_UNTIL;
"jusqu_a"                 index_status(yyextra)=0; return T_UNTIL;
"jusqua"                 index_status(yyextra)=0; return T_UNTIL;
"UNTIL"                 index_status(yyextra)=0; return T_UNTIL;
"virgule"               (*yylval) = gen(at_virgule,2); index_status(yyextra)=0; return T_UNARY_OP;
"VARS"                  (*yylval) = gen(at_VARS,0); index_status(yyextra)=0; return T_UNARY_OP;
"scientific_format"		        (*yylval) = gen(at_scientific_format,0); index_status(yyextra)=0; return T_DIGITS;
"while"                 index_status(yyextra)=0; (*yylval)=gen(at_for,4); if (xcas_mode(yyextra)==3) return TI_WHILE; if (xcas_mode(yyextra)!=0) return T_MUPMAP_WHILE; return T_WHILE;
"WHILE"                 index_status(yyextra)=0; (*yylval)=gen(at_for,4); return T_RPN_WHILE;
"Text"                  (*yylval) = gen(at_Text,1); index_status(yyextra)=0; return T_RETURN;
"DropDown"                  (*yylval) = gen(at_DropDown,1); index_status(yyextra)=0; return T_RETURN;
"Popup"                  (*yylval) = gen(at_Popup,1); index_status(yyextra)=0; return T_RETURN;
"Request"                  (*yylval) = gen(at_Request,1); index_status(yyextra)=0; return T_RETURN;
"Title"                  (*yylval) = gen(at_Title,1); index_status(yyextra)=0; return T_RETURN;
"Define"               (*yylval)=0; index_status(yyextra)=0; return TI_DEFINE;
":Prgm"                 (*yylval)=0; index_status(yyextra)=0; return TI_PRGM;
"Prgm"                 (*yylval)=0; index_status(yyextra)=0; return TI_PRGM;
"EndPrgm"              (*yylval)=0; index_status(yyextra)=0; return T_BLOC_END;
"Endprgm"              (*yylval)=0; index_status(yyextra)=0; return T_BLOC_END;
":Func"                 (*yylval)=0; index_status(yyextra)=0; return TI_PRGM;
"Func"                 (*yylval)=0; index_status(yyextra)=0; return TI_PRGM;
"func"                 (*yylval)=0; index_status(yyextra)=0; return TI_PRGM;
":func"                 (*yylval)=0; index_status(yyextra)=0; return TI_PRGM;
"EndFunc"              (*yylval)=0; index_status(yyextra)=0; return T_BLOC_END;
"endfunc"              (*yylval)=0; index_status(yyextra)=0; return T_BLOC_END;
"Endfunc"              (*yylval)=0; index_status(yyextra)=0; return T_BLOC_END;
"Local"                (*yylval)=0; index_status(yyextra)=0; return TI_LOCAL;
"If"      	       index_status(yyextra)=0; (*yylval)=gen(at_ifte,3); return T_IF;
"Then"                 index_status(yyextra)=0; (*yylval)=gen(at_ifte,3); return T_THEN;
"Else"                 index_status(yyextra)=0; (*yylval)=gen(at_ifte,3); return T_ELSE;
"EndIf"                index_status(yyextra)=0; (*yylval)=gen(at_ifte,3); return T_BLOC_END;
"Endif"                index_status(yyextra)=0; (*yylval)=gen(at_ifte,3); return T_BLOC_END;
"Return"	       (*yylval) = gen(at_return,1) ; index_status(yyextra)=0; return T_RETURN;
"ElseIf"               (*yylval) = gen(at_ifte,3) ;  index_status(yyextra)=0; return T_ELIF;
"Elseif"               (*yylval) = gen(at_ifte,3) ;  index_status(yyextra)=0; return T_ELIF;
"Exit"  	       index_status(yyextra)=0; (*yylval)=gen(at_breakpoint,0); return T_BREAK;
"Loop"                  index_status(yyextra)=0; (*yylval)=gen(at_for,0); return TI_LOOP;
"EndLoop"               index_status(yyextra)=0; (*yylval)=gen(at_for,0); return T_BLOC_END;
"Endloop"               index_status(yyextra)=0; (*yylval)=gen(at_for,0); return T_BLOC_END;
"For"                   index_status(yyextra)=0; (*yylval)=gen(at_for,0); return TI_FOR;
"EndFor"               index_status(yyextra)=0; (*yylval)=gen(at_for,0); return T_BLOC_END;
"Endfor"               index_status(yyextra)=0; (*yylval)=gen(at_for,0); return T_BLOC_END;
"While"                   index_status(yyextra)=0; (*yylval)=gen(at_for,0); return TI_WHILE;
"EndWhile"               index_status(yyextra)=0; (*yylval)=gen(at_for,0); return T_BLOC_END;
"Endwhile"               index_status(yyextra)=0; (*yylval)=gen(at_for,0); return T_BLOC_END;
"Cycle"               index_status(yyextra)=0; (*yylval)=gen(at_for,0); return T_CONTINUE;
"Try"           index_status(yyextra)=0; (*yylval)=gen(at_for,0); return TI_TRY;
"EndTry"           index_status(yyextra)=0; (*yylval)=gen(at_for,0); return T_BLOC_END;
"Endtry"           index_status(yyextra)=0; (*yylval)=gen(at_for,0); return T_BLOC_END;
"Disp"	       (*yylval) = gen(at_print,1) ; index_status(yyextra)=0; return T_RETURN;
"Pause"	       (*yylval) = gen(at_Pause,1) ; index_status(yyextra)=0; return T_RETURN;
"DelVar"       (*yylval) = gen(at_DelVar,1) ; index_status(yyextra)=0; return T_RETURN;
"label"	       (*yylval) = gen(at_label,1) ; index_status(yyextra)=0; return T_RETURN;
"Lbl"	       (*yylval) = gen(at_label,1) ; index_status(yyextra)=0; return T_RETURN;
"goto"	       (*yylval) = gen(at_goto,1) ; index_status(yyextra)=0; return T_RETURN;
"Goto"	       (*yylval) = gen(at_goto,1) ; index_status(yyextra)=0; return T_RETURN;
"Dialog"       (*yylval) = gen(at_Dialog,1) ; index_status(yyextra)=0; return TI_DIALOG; 
"EndDlog"    (*yylval) = gen(at_Dialog,1) ; index_status(yyextra)=0; return T_BLOC_END;
"Row"	       (*yylval) = gen(at_Row,0) ; index_status(yyextra)=0; return T_DIGITS;
"Col"	       (*yylval) = gen(at_Col,0) ; index_status(yyextra)=0; return T_DIGITS;
"threads"      (*yylval) = gen(at_threads,0) ; index_status(yyextra)=0; return T_DIGITS;
"_hbar_"        (*yylval) = symbolic(at_unit,makevecteur(1.05457266e-34,_J_unit*_s_unit)); index_status(yyextra)=0; return T_SYMBOL;
"_c_"        (*yylval) = symbolic(at_unit,makevecteur(299792458,_m_unit/_s_unit)); index_status(yyextra)=0; return T_SYMBOL;
"_g_"        (*yylval) = symbolic(at_unit,makevecteur(9.80665,_m_unit*unitpow(_s_unit,2))); index_status(yyextra)=0; return T_SYMBOL;
"_IO_" (*yylval) = symbolic(at_unit,makevecteur(1e-12,_W_unit*unitpow(_m_unit,-2))); index_status(yyextra)=0; return T_SYMBOL; 
"_epsilonox_" (*yylval) = 3.9; index_status(yyextra)=0; return T_SYMBOL; 
"_epsilonsi_" (*yylval) = 11.9; index_status(yyextra)=0; return T_SYMBOL; 
"_qepsilon0_" (*yylval) = symbolic(at_unit,makevecteur(1.4185979e-30,_F_unit*_C_unit/_m_unit)); index_status(yyextra)=0; return T_SYMBOL; 
"_epsilon0q_" (*yylval) = symbolic(at_unit,makevecteur(55263469.6,_F_unit/(_m_unit*_C_unit))); index_status(yyextra)=0; return T_SYMBOL; 
"_kq_" (*yylval) = symbolic(at_unit,makevecteur(8.617386e-5,_J_unit/(_K_unit*_C_unit))); index_status(yyextra)=0; return T_SYMBOL; 
"_c3_" (*yylval) = symbolic(at_unit,makevecteur(.002897756,_m_unit*_K_unit)); index_status(yyextra)=0; return T_SYMBOL; 
"_lambdac_" (*yylval) = symbolic(at_unit,makevecteur( 0.00242631058e-9,_m_unit)); index_status(yyextra)=0; return T_SYMBOL; 
"_f0_" (*yylval) = symbolic(at_unit,makevecteur(2.4179883e14,_Hz_unit)); index_status(yyextra)=0; return T_SYMBOL; 
"_lambda0_" (*yylval) = symbolic(at_unit,makevecteur(1239.8425e-9,_m_unit)); index_status(yyextra)=0; return T_SYMBOL; 
"_muN_" (*yylval) = symbolic(at_unit,makevecteur(5.0507866e-27,_J_unit/_T_unit)); index_status(yyextra)=0; return T_SYMBOL; 
"_muB_" (*yylval) = symbolic(at_unit,makevecteur( 9.2740154e-24,_J_unit/_T_unit)); index_status(yyextra)=0; return T_SYMBOL; 
"_a0_" (*yylval) = symbolic(at_unit,makevecteur(.0529177249e-9,_m_unit)); index_status(yyextra)=0; return T_SYMBOL; 
"_Rinfinity_" (*yylval) = symbolic(at_unit,makevecteur(10973731.534,unitpow(_m_unit,-1))); index_status(yyextra)=0; return T_SYMBOL; 
"_Faraday_" (*yylval) = symbolic(at_unit,makevecteur(96485.309,_C_unit/_mol_unit)); index_status(yyextra)=0; return T_SYMBOL; 
"_phi_" (*yylval) = symbolic(at_unit,makevecteur(2.06783461e-15,_Wb_unit)); index_status(yyextra)=0; return T_SYMBOL; 
"_alpha_" (*yylval) = 7.29735308e-3; index_status(yyextra)=0; return T_SYMBOL; 
"_mpme_" (*yylval) = 1836.152701; index_status(yyextra)=0; return T_SYMBOL; 
"_mp_" (*yylval) = symbolic(at_unit,makevecteur(1.6726231e-27,_kg_unit)); index_status(yyextra)=0; return T_SYMBOL; 
"_qme_" (*yylval) = symbolic(at_unit,makevecteur(1.75881962e11,_C_unit/_kg_unit)); index_status(yyextra)=0; return T_SYMBOL; 
"_me_" (*yylval) = symbolic(at_unit,makevecteur(9.1093897e-31,_kg_unit)); index_status(yyextra)=0; return T_SYMBOL; 
"_qe_" (*yylval) = symbolic(at_unit,makevecteur(1.60217733e-19,_C_unit)); index_status(yyextra)=0; return T_SYMBOL; 
"_h_" (*yylval) = symbolic(at_unit,makevecteur(6.6260755e-34,_J_unit*_s_unit)); index_status(yyextra)=0; return T_SYMBOL; 
"_G_" (*yylval) = symbolic(at_unit,makevecteur(6.67259e-11,unitpow(_m_unit,3)*unitpow(_s_unit,-2)*unitpow(_kg_unit,-1))); index_status(yyextra)=0; return T_SYMBOL; 
"_mu0_" (*yylval) = symbolic(at_unit,makevecteur(1.25663706144e-6,_H_unit/_m_unit)); index_status(yyextra)=0; return T_SYMBOL; 
"_epsilon0_" (*yylval) = symbolic(at_unit,makevecteur(8.85418781761e-12,_F_unit/_m_unit)); index_status(yyextra)=0; return T_SYMBOL; 
"_sigma_" (*yylval) = symbolic(at_unit,makevecteur( 5.67051e-8,_W_unit*unitpow(_m_unit,-2)*unitpow(_K_unit,-4))); index_status(yyextra)=0; return T_SYMBOL; 
"_StdP_" (*yylval) = symbolic(at_unit,makevecteur(101325.0,_Pa_unit)); index_status(yyextra)=0; return T_SYMBOL; 
"_StdT_" (*yylval) = symbolic(at_unit,makevecteur(273.15,_K_unit)); index_status(yyextra)=0; return T_SYMBOL; 
"_R_" (*yylval) = symbolic(at_unit,makevecteur(8.31451,_J_unit/_molK_unit)); index_status(yyextra)=0; return T_SYMBOL; 
"_Vm_" (*yylval) = symbolic(at_unit,makevecteur(22.4141,_l_unit/_mol_unit)); index_status(yyextra)=0; return T_SYMBOL; 
"_k_" (*yylval) = symbolic(at_unit,makevecteur(1.380658e-23,_J_unit/_K_unit)); index_status(yyextra)=0; return T_SYMBOL; 
"_NA_" (*yylval) = symbolic(at_unit,makevecteur(6.0221367e23,unitpow(_mol_unit,-1))); index_status(yyextra)=0; return T_SYMBOL; 
"_in"        (*yylval) = symbolic(at_unit,makevecteur(1,_in_unit)); index_status(yyextra)=0; return T_SYMBOL;
 
			/* numbers */
{D}+			|
"#"[0-7]+"o"		|
"#"[0-1]+"b"		|
"#"[0-9a-fA-F]+"h"	|
"#o"[0-7]+		|
"#b"[0-1]+		|
"#x"[0-9a-fA-F]+	|
"0o"[0-7]+		|
"0b"[0-1]+		|
"0x"[0-9a-fA-F]+	|
{D}+"."{D}*({E})?	|
{D}*"."{D}+({E})?	|
{D}+{E}			{ 
  index_status(yyextra)=1;
  int l=strlen(yytext);
  char ch,ch2;
  if (l>2 && yytext[1]!='x' && (yytext[l-1]=='o' || yytext[l-1]=='b' || yytext[l-1]=='h') ){
    char base=yytext[l-1];
    for (int i=l-1;i>1;--i){
      yytext[i]=yytext[i-1];
    }
    if (base=='h')
      base='x';
    yytext[1]=base;
  }
  else {
    for (l=0;(ch=*(yytext+l));++l){
      if (ch=='x')
	break;
      if (ch=='e' || ch=='E'){
	if ( (ch2=*(yytext+l+1)) && (ch2=='e' || ch2=='E')){
	  ++l;
	  for (;(ch=*(yytext+l));++l)
	    *(yytext+l-1)=ch;
	  *(yytext+l-1)=0;
	  --l;
	  break;
	}
      }
    }
  }
  (*yylval) = chartab2gen(yytext,yyextra); return T_NUMBER; 
}

			/* symbols */
{A}{AN}*		{
  index_status(yyextra)=1; 
  return find_or_make_symbol(yytext,(*yylval),yyextra);
}
"#"                     if (!xcas_mode(yyextra) || xcas_mode(yyextra)==3) { 
  // cerr << "hash" << endl;
  (*yylval)=gen(at_hash,1); return TI_HASH; 
} else BEGIN(comment_hash);
<comment_hash>[^*\n]*\n BEGIN(INITIAL); index_status(yyextra)=0; increment_lexer_line_number(yyextra);  /* comment_s(yyextra)=string(yytext); (*yylval)=string2gen(comment_s(yyextra).substr(0,comment_s(yyextra).size()-1),false); return T_COMMENT; */
			/* everything else */
.			(*yylval)=string2gen(string(yytext),false); return T_STRING;

%%

/*
 *  Routines
 */
#ifndef NO_NAMESPACE_GIAC
  namespace giac {
#endif // ndef NO_NAMESPACE_GIAC
    
    void update_lexer_localization(const std::vector<int> & v,std::map<std::string,std::string> &lexer_map,std::multimap<std::string,giac::localized_string> &back_lexer_map){
      lexer_map.clear();
      back_lexer_map.clear();
      int s=v.size();
      for (int i=0;i<s;++i){
	int lang=v[i];
	if (lang>=1 && lang<=3){
	  std::string doc=find_doc_prefix(lang);
	  std::string file=giac::giac_aide_dir()+doc+"keywords";
	  std::string giac_kw,local_kw;
	  size_t l;
	  char * line = (char *)malloc(1024);
	  ifstream f(file.c_str());
	  if (f){
	    cerr << "// Using keyword file " << file << endl;
	    for (;;){
	      f.getline(line,1023,'\n');
	      l=strlen(line);
	      if (f.eof()){
		f.close();
		break;
	      }
	      if (l>3 && line[0]!='#'){
		if (line[l-1]=='\n')
		  --l;
		// read giac keyword
		size_t j;
		giac_kw="";
		for (j=0;j<l;++j){
		  if (line[j]==' ')
		    break;
		  giac_kw += line[j];
		}
		// read corresponding local keywords
		local_kw="";
		for (++j;j<l;++j){
		  if (line[j]==' '){
		    lexer_map[local_kw]=giac_kw;
		    back_lexer_map.insert(pair<string,localized_string>(giac_kw,localized_string(lang,local_kw)));
		    local_kw="";
		  }
		  else
		    local_kw += line[j];
		}
		if (!local_kw.empty()){
		  lexer_map[local_kw]=giac_kw;
		  back_lexer_map.insert(pair<string,localized_string>(giac_kw,localized_string(lang,local_kw)));
		}
	      }
	    }
	    free(line);
	  } // if (f)
	  else
	    cerr << "// Unable to find keyword file " << file << endl;
	}
      }
    }

    bool has_special_syntax(const string & s){
      sym_tab::const_iterator i = lexer_functions().find(s);
      if (i==lexer_functions().end())
	return false;
      return (i->second.subtype!=T_UNARY_OP);
    }
    
    bool lexer_functions_register(const unary_function_ptr & u,const string & s,int parser_token){
      sym_tab::const_iterator i = lexer_functions().find(s);
      if (i!=lexer_functions().end())
	return false;
      registered_lexer_functions().push_back(user_function(s,parser_token));
      lexer_functions()[s] = gen(u);
      if (parser_token==1)
	lexer_functions()[s].subtype=T_UNARY_OP;
      else
	lexer_functions()[s].subtype=parser_token;
      // If s is a library function name (with ::), update the library
      int ss=s.size(),j=0;
      for (;j<ss-1;++j){
	if (s[j]==':' && s[j+1]==':')
	  break;
      }
      if (j<ss-1){
	string libname=s.substr(0,j);
	string funcname=s.substr(j+2,ss-j-2);
	std::map<std::string,std::vector<string> >::iterator it=library_functions().find(libname);
	if (it!=library_functions().end())
	  it->second.push_back(funcname);
	else
	  library_functions()[libname]=vector<string>(1,funcname);
      }
      return true;
    }

    bool lexer_function_remove(const vector<user_function> & v){
      vector<user_function>::const_iterator it=v.begin(),itend=v.end();
      sym_tab::const_iterator i,iend;
      bool ok=true;
      for (;it!=itend;++it){
	i = lexer_functions().find(it->s);
	iend=lexer_functions().end();
	if (i==iend)
	  ok=false;
	else
	  lexer_functions().erase(it->s);
      }
      return ok;
    }

    int find_or_make_symbol(const string & s,gen & res,GIAC_CONTEXT){
      if (s.size()==1){
	switch (s[0]){
	case '+':
	  res=at_plus;
	  return T_UNARY_OP;
	case '-':
	  res=at_neg;
	  return T_UNARY_OP;
	case '*':
	  res=at_prod;
	  return T_UNARY_OP;
	case '/':
	  res=at_division;
	  return T_UNARY_OP;
	case '^':
	  res=at_pow;
	  return T_UNARY_OP;
	}
      }
      string ts(s);
      std::map<std::string,std::string>::const_iterator trans=lexer_localization_map().find(ts);
      if (trans!=lexer_localization_map().end())
	ts=trans->second;
      std::map<std::string,std::vector<string> >::const_iterator j=lexer_translator().find(ts);
      if (j!=lexer_translator().end() && !j->second.empty())
	ts=j->second.back();
      sym_tab::const_iterator i = lexer_functions().find(ts);
      if (i!=lexer_functions().end()){
	if (i->second.subtype==T_TO)
	  res=plus_one;
	else
	  res = i->second;
	res.subtype=1;
	index_status(contextptr)=(i->second.subtype==T_UNARY_OP);
	return i->second.subtype ;
      }
      i = syms().find(s);
      if (i == syms().end()) {
	// std::cerr << "lexer new" << s << endl;
	res = *(new identificateur(s));
	syms()[s] = res;
      } else {
	// std::cerr << "lexer" << s << endl;
	res = i->second;
      }
      return T_SYMBOL;
    }

  // Add to the list of predefined symbols
  void set_lexer_symbols(const vecteur & l,GIAC_CONTEXT){
    if (initialisation_done(contextptr) && (l==list_one_letter__IDNT) )
      return;
    initialisation_done(contextptr)=true;
    const_iterateur it=l.begin(),itend=l.end();
    for (; it!=itend; ++it) {
      if (it->type!=_IDNT)
	continue;
      sym_tab::const_iterator i = syms().find(* (it->_IDNTptr->name));
      if (i==syms().end())
	syms()[* (it->_IDNTptr->name)] = *it;
    }
  }


    // Set the input string
    YY_BUFFER_STATE set_lexer_string(const std::string &s_orig,yyscan_t & scanner,GIAC_CONTEXT){
      string s(s_orig),lexer_string;
      bool instring=false;
      // stupid match of bracket then parenthesis
      int l=s.size(),nb=0,np=0;
      int i=0;
      for (;i<l;++i){
	if (s[i]==92){
	  i += 2;
	  if (i>=l)
	    break;
	}
	if (instring){
	  if (s[i]=='"')
	    instring=false;
	}
	else {
	  switch (s[i]){
	  case '"':
	    instring=true;
	    break;
	  case '(':
	    ++np;
	    break;
	  case ')':
	    --np;
	    break;
	  case '[':
	    ++nb;
	    break;
	  case ']':
	    --nb;
	    break;
	  }
	}
      }
      while (np<0 && i>=0 && s[i-1]==')'){
	--i;
	++np;
      }
      while (nb<0 && i>=0 && s[i-1]==']'){
	--i;
	++nb;
      }
      s=s.substr(0,i);
      if (nb<0)
	cerr << "Too many ]" << endl;
      if (np<0)
	cerr << "Too many )" << endl;
      if (nb>0)
	s=s+string(nb,']');
      if (np>0)
	s=s+string(np,')');
      index_status(contextptr)=0;
      opened_quote(contextptr)=0;
      in_rpn(contextptr)=0;
      lexer_line_number(contextptr)=1;
      first_error_line(contextptr)=0;
      spread_formula(contextptr)=0;
      l=s.size();
      for (;l;l--){
	if (s[l-1]!=' ')
	  break;
      }
      while (l>=4 && s[l-1]==';' && s[l-2]==':' && s[l-3]==';'){
	if (s[l-4]==':')
	  l -= 2;
	else {
	  s[l-3]=':';
	  s[l-2]=';';
	  l--;
	}
      }
      s=s.substr(0,l);
      /* if (l && ( (s[l-1]==';') || (s[l-1]==':')))
	 l--; */
      string ss;
      for (int i=0;i<l;++i){
	if (s[i]=='.'){
	  if ( i && (i<l-1) && (s[i-1]!=' ') && (s[i+1]=='.') ){
	    ss+= " ..";
	    ++i;
	  }
	  else
	    ss+='.';
	}
	else {
	  if (xcas_mode(contextptr) > 0 && xcas_mode(contextptr) !=3){
	    if (s[i]=='#')
	      ss += "//";
	    else
	      ss += s[i];
	  }
	  else
	    ss+=s[i];
	}
      }
      lexer_string = ss+"\n§";
      yylex_init(&scanner);
      yyset_extra(contextptr, scanner);
      YY_BUFFER_STATE state=yy_scan_string(lexer_string.c_str(),scanner);
      return state;
    }

    int delete_lexer_string(YY_BUFFER_STATE & state,yyscan_t & scanner){
      yy_delete_buffer(state,scanner);
      yylex_destroy(scanner);
      return 1;
    }

#ifndef NO_NAMESPACE_GIAC
  } // namespace giac
#endif // ndef NO_NAMESPACE_GIAC
  
