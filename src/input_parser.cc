/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.3"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 1

/* Using locations.  */
#define YYLSP_NEEDED 0

/* Substitute the variable and function names.  */
#define yyparse giac_yyparse
#define yylex   giac_yylex
#define yyerror giac_yyerror
#define yylval  giac_yylval
#define yychar  giac_yychar
#define yydebug giac_yydebug
#define yynerrs giac_yynerrs


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     T_NUMBER = 258,
     T_SYMBOL = 259,
     T_LITERAL = 260,
     T_DIGITS = 261,
     T_STRING = 262,
     T_END_INPUT = 263,
     T_EXPRESSION = 264,
     T_UNARY_OP = 265,
     T_OF = 266,
     T_NOT = 267,
     T_TYPE_ID = 268,
     T_VIRGULE = 269,
     T_AFFECT = 270,
     T_MAPSTO = 271,
     T_BEGIN_PAR = 272,
     T_END_PAR = 273,
     T_PLUS = 274,
     T_MOINS = 275,
     T_FOIS = 276,
     T_DIV = 277,
     T_MOD = 278,
     T_POW = 279,
     T_QUOTED_BINARY = 280,
     T_QUOTE = 281,
     T_PRIME = 282,
     T_TEST_EQUAL = 283,
     T_EQUAL = 284,
     T_INTERVAL = 285,
     T_UNION = 286,
     T_INTERSECT = 287,
     T_MINUS = 288,
     T_AND_OP = 289,
     T_COMPOSE = 290,
     T_DOLLAR = 291,
     T_DOLLAR_MAPLE = 292,
     T_INDEX_BEGIN = 293,
     T_VECT_BEGIN = 294,
     T_VECT_DISPATCH = 295,
     T_VECT_END = 296,
     T_SET_BEGIN = 297,
     T_SET_END = 298,
     T_SEMI = 299,
     T_DEUXPOINTS = 300,
     T_DOUBLE_DEUX_POINTS = 301,
     T_IF = 302,
     T_RPN_IF = 303,
     T_ELIF = 304,
     T_THEN = 305,
     T_ELSE = 306,
     T_IFTE = 307,
     T_SWITCH = 308,
     T_CASE = 309,
     T_DEFAULT = 310,
     T_ENDCASE = 311,
     T_FOR = 312,
     T_FROM = 313,
     T_TO = 314,
     T_DO = 315,
     T_BY = 316,
     T_WHILE = 317,
     T_MUPMAP_WHILE = 318,
     T_RPN_WHILE = 319,
     T_REPEAT = 320,
     T_UNTIL = 321,
     T_IN = 322,
     T_START = 323,
     T_BREAK = 324,
     T_CONTINUE = 325,
     T_TRY = 326,
     T_CATCH = 327,
     T_TRY_CATCH = 328,
     T_PROC = 329,
     T_BLOC = 330,
     T_BLOC_BEGIN = 331,
     T_BLOC_END = 332,
     T_RETURN = 333,
     T_LOCAL = 334,
     T_LOCALBLOC = 335,
     T_NAME = 336,
     T_PROGRAM = 337,
     T_NULL = 338,
     T_ARGS = 339,
     T_FACTORIAL = 340,
     T_RPN_OP = 341,
     T_RPN_BEGIN = 342,
     T_RPN_END = 343,
     T_STACK = 344,
     T_GROUPE_BEGIN = 345,
     T_GROUPE_END = 346,
     T_LINE_BEGIN = 347,
     T_LINE_END = 348,
     T_VECTOR_BEGIN = 349,
     T_VECTOR_END = 350,
     T_CURVE_BEGIN = 351,
     T_CURVE_END = 352,
     T_ROOTOF_BEGIN = 353,
     T_ROOTOF_END = 354,
     T_SPOLY1_BEGIN = 355,
     T_SPOLY1_END = 356,
     T_POLY1_BEGIN = 357,
     T_POLY1_END = 358,
     T_MATRICE_BEGIN = 359,
     T_MATRICE_END = 360,
     T_ASSUME_BEGIN = 361,
     T_ASSUME_END = 362,
     T_HELP = 363,
     TI_DEUXPOINTS = 364,
     TI_LOCAL = 365,
     TI_LOOP = 366,
     TI_FOR = 367,
     TI_WHILE = 368,
     TI_STO = 369,
     TI_TRY = 370,
     TI_DIALOG = 371,
     T_PIPE = 372,
     TI_DEFINE = 373,
     TI_PRGM = 374,
     TI_SEMI = 375,
     TI_HASH = 376,
     T_ACCENTGRAVE = 377,
     T_MAPLELIB = 378,
     T_INTERROGATION = 379,
     T_UNIT = 380,
     T_BIDON = 381,
     T_LOGO = 382,
     T_SQ = 383,
     NEG = 384
   };
#endif
/* Tokens.  */
#define T_NUMBER 258
#define T_SYMBOL 259
#define T_LITERAL 260
#define T_DIGITS 261
#define T_STRING 262
#define T_END_INPUT 263
#define T_EXPRESSION 264
#define T_UNARY_OP 265
#define T_OF 266
#define T_NOT 267
#define T_TYPE_ID 268
#define T_VIRGULE 269
#define T_AFFECT 270
#define T_MAPSTO 271
#define T_BEGIN_PAR 272
#define T_END_PAR 273
#define T_PLUS 274
#define T_MOINS 275
#define T_FOIS 276
#define T_DIV 277
#define T_MOD 278
#define T_POW 279
#define T_QUOTED_BINARY 280
#define T_QUOTE 281
#define T_PRIME 282
#define T_TEST_EQUAL 283
#define T_EQUAL 284
#define T_INTERVAL 285
#define T_UNION 286
#define T_INTERSECT 287
#define T_MINUS 288
#define T_AND_OP 289
#define T_COMPOSE 290
#define T_DOLLAR 291
#define T_DOLLAR_MAPLE 292
#define T_INDEX_BEGIN 293
#define T_VECT_BEGIN 294
#define T_VECT_DISPATCH 295
#define T_VECT_END 296
#define T_SET_BEGIN 297
#define T_SET_END 298
#define T_SEMI 299
#define T_DEUXPOINTS 300
#define T_DOUBLE_DEUX_POINTS 301
#define T_IF 302
#define T_RPN_IF 303
#define T_ELIF 304
#define T_THEN 305
#define T_ELSE 306
#define T_IFTE 307
#define T_SWITCH 308
#define T_CASE 309
#define T_DEFAULT 310
#define T_ENDCASE 311
#define T_FOR 312
#define T_FROM 313
#define T_TO 314
#define T_DO 315
#define T_BY 316
#define T_WHILE 317
#define T_MUPMAP_WHILE 318
#define T_RPN_WHILE 319
#define T_REPEAT 320
#define T_UNTIL 321
#define T_IN 322
#define T_START 323
#define T_BREAK 324
#define T_CONTINUE 325
#define T_TRY 326
#define T_CATCH 327
#define T_TRY_CATCH 328
#define T_PROC 329
#define T_BLOC 330
#define T_BLOC_BEGIN 331
#define T_BLOC_END 332
#define T_RETURN 333
#define T_LOCAL 334
#define T_LOCALBLOC 335
#define T_NAME 336
#define T_PROGRAM 337
#define T_NULL 338
#define T_ARGS 339
#define T_FACTORIAL 340
#define T_RPN_OP 341
#define T_RPN_BEGIN 342
#define T_RPN_END 343
#define T_STACK 344
#define T_GROUPE_BEGIN 345
#define T_GROUPE_END 346
#define T_LINE_BEGIN 347
#define T_LINE_END 348
#define T_VECTOR_BEGIN 349
#define T_VECTOR_END 350
#define T_CURVE_BEGIN 351
#define T_CURVE_END 352
#define T_ROOTOF_BEGIN 353
#define T_ROOTOF_END 354
#define T_SPOLY1_BEGIN 355
#define T_SPOLY1_END 356
#define T_POLY1_BEGIN 357
#define T_POLY1_END 358
#define T_MATRICE_BEGIN 359
#define T_MATRICE_END 360
#define T_ASSUME_BEGIN 361
#define T_ASSUME_END 362
#define T_HELP 363
#define TI_DEUXPOINTS 364
#define TI_LOCAL 365
#define TI_LOOP 366
#define TI_FOR 367
#define TI_WHILE 368
#define TI_STO 369
#define TI_TRY 370
#define TI_DIALOG 371
#define T_PIPE 372
#define TI_DEFINE 373
#define TI_PRGM 374
#define TI_SEMI 375
#define TI_HASH 376
#define T_ACCENTGRAVE 377
#define T_MAPLELIB 378
#define T_INTERROGATION 379
#define T_UNIT 380
#define T_BIDON 381
#define T_LOGO 382
#define T_SQ 383
#define NEG 384




/* Copy the first part of user declarations.  */
#line 25 "input_parser.yy"

         #define YYPARSE_PARAM scanner
         #define YYLEX_PARAM   scanner
	 #line 34 "input_parser.yy"

#include "first.h"
#include <stdexcept>
#include <cstdlib>
#include "index.h"
#include "gen.h"
#define YYSTYPE giac::gen
#define YY_EXTRA_TYPE  const giac::context *
#include "lexer.h"
#include "input_lexer.h"
#include "usual.h"
#include "derive.h"
#include "sym2poly.h"
#include "vecteur.h"
#include "modpoly.h"
#include "alg_ext.h"
#include "prog.h"
#include "rpn.h"
#include "intg.h"
#include "plot.h"

using namespace std;
#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

// It seems there is a bison bug when it reallocates space for the stack
// therefore I redefine YYINITDEPTH to 1000 (max size is YYMAXDEPTH)
// instead of 200
// Feel free to change if you need but then readjust YYMAXDEPTH
#ifdef GNUWINCE
#define YYINITDEPTH 1000
#else
#define YYINITDEPTH 20000
#define YYMAXDEPTH 100000
#define YYERROR_VERBOSE 1
#endif
// Note that the compilation by bison with -v option generates a file y.output
// to debug the grammar, compile input_parser.yy with bison
// then add yydebug=1 in input_parser.cc at the beginning of yyparse (
// #define YYDEBUG 1


gen polynome_or_sparse_poly1(const gen & coeff, const gen & index){
  if (index.type==_VECT){
    index_t i;
    const_iterateur it=index._VECTptr->begin(),itend=index._VECTptr->end();
    i.reserve(itend-it);
    for (;it!=itend;++it){
      if (it->type!=_INT_)
         settypeerr();
      i.push_back(it->val);
    }
    monomial<gen> m(coeff,i);
    return polynome(m);
  }
  else {
    sparse_poly1 res;
    res.push_back(monome(coeff,index));
    return res;
  }
}


/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif

#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef int YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 440 "y.tab.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int i)
#else
static int
YYID (i)
    int i;
#endif
{
  return i;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  183
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   11267

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  130
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  29
/* YYNRULES -- Number of rules.  */
#define YYNRULES  260
/* YYNRULES -- Number of states.  */
#define YYNSTATES  592

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   384

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   109,   110,   111,   112,   113,   114,
     115,   116,   117,   118,   119,   120,   121,   122,   123,   124,
     125,   126,   127,   128,   129
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     5,     8,    12,    16,    18,    21,    26,
      30,    33,    35,    37,    44,    51,    58,    65,    69,    83,
      96,   109,   118,   122,   125,   128,   136,   150,   160,   165,
     170,   172,   177,   179,   181,   183,   187,   191,   195,   198,
     202,   206,   210,   214,   218,   222,   226,   230,   234,   237,
     240,   246,   250,   256,   258,   262,   265,   270,   275,   277,
     282,   287,   291,   293,   296,   299,   303,   311,   320,   327,
     335,   342,   347,   352,   358,   363,   365,   370,   372,   376,
     380,   385,   387,   392,   394,   397,   400,   402,   406,   409,
     411,   415,   417,   419,   429,   440,   445,   453,   463,   473,
     481,   491,   493,   499,   505,   511,   516,   520,   527,   533,
     537,   545,   550,   552,   558,   563,   570,   576,   584,   589,
     594,   596,   600,   603,   609,   613,   617,   620,   624,   628,
     632,   636,   640,   642,   646,   650,   657,   662,   666,   670,
     674,   678,   682,   686,   690,   694,   698,   702,   706,   710,
     714,   716,   718,   720,   723,   727,   730,   734,   737,   739,
     741,   745,   748,   751,   752,   755,   758,   763,   766,   770,
     774,   776,   780,   782,   786,   790,   792,   796,   798,   800,
     801,   803,   804,   806,   808,   811,   814,   815,   818,   822,
     824,   828,   830,   832,   834,   837,   841,   845,   847,   849,
     851,   853,   855,   857,   859,   861,   863,   865,   867,   869,
     871,   873,   875,   877,   879,   883,   885,   889,   893,   895,
     901,   909,   913,   917,   922,   927,   933,   939,   943,   947,
     952,   954,   958,   962,   964,   967,   968,   971,   972,   975,
     976,   980,   983,   987,   992,   994,   998,  1004,  1011,  1013,
    1016,  1018,  1021,  1022,  1028,  1029,  1033,  1039,  1040,  1043,
    1049
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int16 yyrhs[] =
{
     131,     0,    -1,   132,    -1,   133,     8,    -1,   133,    44,
       8,    -1,   133,    44,   132,    -1,     3,    -1,     3,     4,
      -1,     3,     4,    24,   133,    -1,     3,     4,   128,    -1,
       3,     5,    -1,     7,    -1,     9,    -1,   134,    17,   142,
      18,    15,   151,    -1,   134,    17,   142,    18,    15,   133,
      -1,   133,   114,   134,    17,   142,    18,    -1,   133,   114,
     134,    38,   133,    41,    -1,   133,   114,   134,    -1,   109,
     134,    17,   142,    18,   119,   143,   109,   110,   142,   109,
     143,   153,    -1,   109,   134,    17,   142,    18,   119,   143,
     110,   142,   109,   143,   153,    -1,   109,   134,    17,   142,
      18,   119,   109,   110,   142,   109,   143,   153,    -1,   109,
     134,    17,   142,    18,   119,   143,   153,    -1,   116,   143,
     153,    -1,   116,   151,    -1,   109,   133,    -1,   118,   134,
      17,   142,    18,    29,   133,    -1,   118,   134,    17,   142,
      18,    29,   119,   109,   110,   142,   109,   143,   153,    -1,
     118,   134,    17,   142,    18,    29,   119,   143,   153,    -1,
     134,    17,   142,    18,    -1,   133,    17,   142,    18,    -1,
     123,    -1,   123,    38,   133,    41,    -1,   134,    -1,     5,
      -1,     6,    -1,     6,    15,   133,    -1,   133,    28,   133,
      -1,   133,    29,   133,    -1,    29,   133,    -1,   133,    19,
     133,    -1,   133,    20,   133,    -1,   133,    21,   133,    -1,
     133,    22,   133,    -1,   133,    24,   133,    -1,   133,    23,
     133,    -1,   133,    30,   133,    -1,   133,    34,   133,    -1,
     133,    45,   133,    -1,    20,   133,    -1,    19,   133,    -1,
     100,   133,    14,   133,   101,    -1,    98,   133,    99,    -1,
      98,   133,    14,   133,    99,    -1,    11,    -1,   133,    15,
     133,    -1,    12,   133,    -1,    84,    17,   133,    18,    -1,
      84,    38,   133,    41,    -1,    84,    -1,    10,    17,   133,
      18,    -1,    10,    38,   133,    41,    -1,    10,    17,    18,
      -1,    10,    -1,   133,    27,    -1,   133,    85,    -1,    87,
     144,    88,    -1,    74,    17,   142,    18,   135,   143,    77,
      -1,    74,    17,   142,    18,   135,    76,   143,    77,    -1,
      47,    17,   133,    18,   151,   150,    -1,    47,    17,   133,
      18,   133,    44,   150,    -1,    47,   133,    50,   151,    51,
     151,    -1,    47,   133,    50,   151,    -1,    47,   133,   109,
     133,    -1,    47,   133,    50,   143,   152,    -1,    52,    17,
     133,    18,    -1,    52,    -1,    82,    17,   133,    18,    -1,
      82,    -1,   133,    16,   151,    -1,   133,    16,   133,    -1,
      75,    17,   133,    18,    -1,    75,    -1,    80,    17,   133,
      18,    -1,    80,    -1,    78,   109,    -1,    78,   133,    -1,
      78,    -1,    26,    78,    26,    -1,   127,   133,    -1,   127,
      -1,   127,    17,    18,    -1,    69,    -1,    70,    -1,    57,
      17,   141,    44,   141,    44,   141,    18,   151,    -1,    57,
      17,   141,    44,   141,    44,   141,    18,   133,    44,    -1,
      57,    17,   133,    18,    -1,    57,   134,    67,   133,    60,
     143,    77,    -1,    57,   134,   149,    59,   133,   148,    60,
     143,    77,    -1,    57,   134,   149,   148,    59,   133,    60,
     143,    77,    -1,    57,   134,   149,   148,    60,   143,    77,
      -1,    57,   134,   149,   148,    63,   133,    60,   143,    77,
      -1,    57,    -1,   112,   142,   109,   143,   153,    -1,   113,
     133,   109,   143,   153,    -1,    62,    17,   133,    18,   151,
      -1,    65,   143,    66,   133,    -1,    60,   143,    77,    -1,
      62,    17,   133,    18,   133,    44,    -1,    63,   133,    60,
     143,    77,    -1,   111,   143,   153,    -1,    71,   151,    72,
      17,   133,    18,   151,    -1,    73,    17,   133,    18,    -1,
      73,    -1,   115,   143,    51,   143,   153,    -1,   115,   143,
      51,   153,    -1,   115,   143,   109,    51,   143,   153,    -1,
     115,   143,   109,    51,   153,    -1,    53,    17,   133,    18,
      76,   156,    77,    -1,    54,    17,     4,    18,    -1,    54,
     133,   157,    56,    -1,    13,    -1,    26,    13,    26,    -1,
      37,   133,    -1,   133,    37,   134,    67,   133,    -1,   133,
      37,   133,    -1,   133,    36,   133,    -1,    36,     4,    -1,
     133,    35,   133,    -1,   133,    31,   133,    -1,   133,    32,
     133,    -1,   133,    33,   133,    -1,   133,   117,   133,    -1,
      25,    -1,    26,   133,    26,    -1,   122,   145,   122,    -1,
      17,   133,    18,    17,   142,    18,    -1,   133,    38,   133,
      41,    -1,    17,   133,    18,    -1,    40,   142,    41,    -1,
      39,   142,    41,    -1,    42,   142,    43,    -1,    90,   142,
      91,    -1,    92,   142,    93,    -1,    94,   142,    95,    -1,
      96,   142,    97,    -1,   102,   142,   103,    -1,   104,   142,
     105,    -1,   106,   142,   107,    -1,   133,    14,   133,    -1,
     133,   120,   133,    -1,    83,    -1,   136,    -1,    86,    -1,
     108,   133,    -1,   133,   124,   133,    -1,   125,   133,    -1,
     133,   125,   133,    -1,   133,   128,    -1,     1,    -1,     4,
      -1,     4,    46,   133,    -1,    13,     4,    -1,   121,   133,
      -1,    -1,   135,   137,    -1,   138,   135,    -1,    89,    17,
     133,    18,    -1,    89,    83,    -1,    79,   139,    44,    -1,
      81,   133,    44,    -1,   140,    -1,   140,    14,   139,    -1,
     134,    -1,     4,    15,   133,    -1,    17,   140,    18,    -1,
      10,    -1,    10,    46,   133,    -1,    13,    -1,     3,    -1,
      -1,   133,    -1,    -1,   133,    -1,   133,    -1,   143,   133,
      -1,   143,   158,    -1,    -1,   145,   144,    -1,   145,    14,
     144,    -1,    10,    -1,    26,    10,    26,    -1,     3,    -1,
     134,    -1,     7,    -1,   125,   145,    -1,     3,   125,   145,
      -1,    40,   144,    41,    -1,    19,    -1,    20,    -1,    22,
      -1,    21,    -1,    24,    -1,    29,    -1,    23,    -1,    30,
      -1,    34,    -1,    28,    -1,    11,    -1,    36,    -1,    35,
      -1,    31,    -1,    32,    -1,    33,    -1,    86,    -1,    26,
      86,    26,    -1,    25,    -1,    87,   144,    88,    -1,    26,
     133,    26,    -1,    52,    -1,    48,   144,    50,   144,    77,
      -1,    48,   144,    50,   144,    51,   144,    77,    -1,    68,
     144,    61,    -1,    68,   144,    70,    -1,    57,   134,   144,
      61,    -1,    57,   134,   144,    70,    -1,    64,   144,    65,
     144,    77,    -1,    60,   144,    66,   144,    77,    -1,    16,
     147,   146,    -1,    54,   155,    77,    -1,    54,   155,   144,
      77,    -1,   136,    -1,    87,   144,    88,    -1,    26,   133,
      26,    -1,   134,    -1,   147,   134,    -1,    -1,    61,   133,
      -1,    -1,    58,   133,    -1,    -1,   154,   133,    44,    -1,
     154,   151,    -1,    76,   143,    77,    -1,    76,   135,   143,
      77,    -1,   153,    -1,   154,   143,   153,    -1,    49,   133,
      50,   143,   152,    -1,   109,    49,   133,    50,   143,   152,
      -1,    77,    -1,   109,    77,    -1,    51,    -1,   109,    51,
      -1,    -1,   155,   144,    50,   144,    77,    -1,    -1,    55,
      45,   151,    -1,    54,     3,    45,   151,   156,    -1,    -1,
      55,   143,    -1,    11,     3,    60,   143,   157,    -1,    44,
      -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   173,   173,   181,   182,   183,   186,   187,   188,   189,
     190,   191,   192,   195,   196,   197,   198,   202,   203,   207,
     211,   215,   218,   219,   220,   221,   222,   226,   227,   228,
     229,   230,   231,   232,   233,   234,   235,   236,   237,   243,
     245,   246,   247,   248,   249,   250,   253,   254,   255,   260,
     266,   267,   272,   274,   275,   276,   277,   278,   279,   280,
     290,   294,   301,   304,   305,   306,   307,   311,   315,   318,
     330,   331,   332,   333,   336,   337,   338,   339,   340,   344,
     351,   352,   353,   357,   359,   360,   361,   362,   363,   364,
     365,   367,   368,   369,   370,   371,   372,   379,   387,   395,
     396,   397,   398,   417,   421,   425,   430,   431,   434,   435,
     436,   437,   438,   439,   440,   441,   442,   443,   444,   445,
     446,   450,   451,   452,   453,   454,   455,   456,   457,   458,
     459,   460,   463,   464,   469,   470,   471,   475,   476,   481,
     482,   483,   484,   485,   486,   487,   488,   489,   490,   495,
     496,   497,   498,   499,   500,   501,   505,   509,   510,   519,
     526,   533,   538,   542,   543,   544,   548,   549,   552,   555,
     558,   559,   566,   567,   568,   569,   570,   571,   579,   589,
     590,   593,   594,   597,   599,   604,   607,   608,   609,   612,
     613,   614,   615,   616,   617,   621,   625,   626,   627,   628,
     629,   630,   631,   632,   633,   634,   635,   636,   637,   638,
     639,   640,   641,   642,   643,   644,   645,   646,   647,   648,
     649,   650,   651,   652,   653,   654,   655,   656,   657,   658,
     659,   662,   663,   666,   667,   670,   671,   674,   675,   678,
     679,   680,   684,   687,   694,   695,   698,   701,   706,   707,
     710,   711,   714,   715,   722,   723,   724,   727,   728,   729,
     732
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "T_NUMBER", "T_SYMBOL", "T_LITERAL",
  "T_DIGITS", "T_STRING", "T_END_INPUT", "T_EXPRESSION", "T_UNARY_OP",
  "T_OF", "T_NOT", "T_TYPE_ID", "T_VIRGULE", "T_AFFECT", "T_MAPSTO",
  "T_BEGIN_PAR", "T_END_PAR", "T_PLUS", "T_MOINS", "T_FOIS", "T_DIV",
  "T_MOD", "T_POW", "T_QUOTED_BINARY", "T_QUOTE", "T_PRIME",
  "T_TEST_EQUAL", "T_EQUAL", "T_INTERVAL", "T_UNION", "T_INTERSECT",
  "T_MINUS", "T_AND_OP", "T_COMPOSE", "T_DOLLAR", "T_DOLLAR_MAPLE",
  "T_INDEX_BEGIN", "T_VECT_BEGIN", "T_VECT_DISPATCH", "T_VECT_END",
  "T_SET_BEGIN", "T_SET_END", "T_SEMI", "T_DEUXPOINTS",
  "T_DOUBLE_DEUX_POINTS", "T_IF", "T_RPN_IF", "T_ELIF", "T_THEN", "T_ELSE",
  "T_IFTE", "T_SWITCH", "T_CASE", "T_DEFAULT", "T_ENDCASE", "T_FOR",
  "T_FROM", "T_TO", "T_DO", "T_BY", "T_WHILE", "T_MUPMAP_WHILE",
  "T_RPN_WHILE", "T_REPEAT", "T_UNTIL", "T_IN", "T_START", "T_BREAK",
  "T_CONTINUE", "T_TRY", "T_CATCH", "T_TRY_CATCH", "T_PROC", "T_BLOC",
  "T_BLOC_BEGIN", "T_BLOC_END", "T_RETURN", "T_LOCAL", "T_LOCALBLOC",
  "T_NAME", "T_PROGRAM", "T_NULL", "T_ARGS", "T_FACTORIAL", "T_RPN_OP",
  "T_RPN_BEGIN", "T_RPN_END", "T_STACK", "T_GROUPE_BEGIN", "T_GROUPE_END",
  "T_LINE_BEGIN", "T_LINE_END", "T_VECTOR_BEGIN", "T_VECTOR_END",
  "T_CURVE_BEGIN", "T_CURVE_END", "T_ROOTOF_BEGIN", "T_ROOTOF_END",
  "T_SPOLY1_BEGIN", "T_SPOLY1_END", "T_POLY1_BEGIN", "T_POLY1_END",
  "T_MATRICE_BEGIN", "T_MATRICE_END", "T_ASSUME_BEGIN", "T_ASSUME_END",
  "T_HELP", "TI_DEUXPOINTS", "TI_LOCAL", "TI_LOOP", "TI_FOR", "TI_WHILE",
  "TI_STO", "TI_TRY", "TI_DIALOG", "T_PIPE", "TI_DEFINE", "TI_PRGM",
  "TI_SEMI", "TI_HASH", "T_ACCENTGRAVE", "T_MAPLELIB", "T_INTERROGATION",
  "T_UNIT", "T_BIDON", "T_LOGO", "T_SQ", "NEG", "$accept", "input",
  "correct_input", "exp", "symbol", "entete", "stack", "local", "nom",
  "suite_symbol", "affectable_symbol", "exp_or_empty", "suite",
  "prg_suite", "rpn_suite", "rpn_token", "rpn_sub_prog", "symbol_suite",
  "step", "from", "else", "bloc", "elif", "ti_bloc_end", "ti_else",
  "rpn_case", "switch", "case", "semi", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   343,   344,
     345,   346,   347,   348,   349,   350,   351,   352,   353,   354,
     355,   356,   357,   358,   359,   360,   361,   362,   363,   364,
     365,   366,   367,   368,   369,   370,   371,   372,   373,   374,
     375,   376,   377,   378,   379,   380,   381,   382,   383,   384
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,   130,   131,   132,   132,   132,   133,   133,   133,   133,
     133,   133,   133,   133,   133,   133,   133,   133,   133,   133,
     133,   133,   133,   133,   133,   133,   133,   133,   133,   133,
     133,   133,   133,   133,   133,   133,   133,   133,   133,   133,
     133,   133,   133,   133,   133,   133,   133,   133,   133,   133,
     133,   133,   133,   133,   133,   133,   133,   133,   133,   133,
     133,   133,   133,   133,   133,   133,   133,   133,   133,   133,
     133,   133,   133,   133,   133,   133,   133,   133,   133,   133,
     133,   133,   133,   133,   133,   133,   133,   133,   133,   133,
     133,   133,   133,   133,   133,   133,   133,   133,   133,   133,
     133,   133,   133,   133,   133,   133,   133,   133,   133,   133,
     133,   133,   133,   133,   133,   133,   133,   133,   133,   133,
     133,   133,   133,   133,   133,   133,   133,   133,   133,   133,
     133,   133,   133,   133,   133,   133,   133,   133,   133,   133,
     133,   133,   133,   133,   133,   133,   133,   133,   133,   133,
     133,   133,   133,   133,   133,   133,   133,   133,   133,   134,
     134,   134,   134,   135,   135,   135,   136,   136,   137,   138,
     139,   139,   140,   140,   140,   140,   140,   140,   140,   141,
     141,   142,   142,   143,   143,   143,   144,   144,   144,   145,
     145,   145,   145,   145,   145,   145,   145,   145,   145,   145,
     145,   145,   145,   145,   145,   145,   145,   145,   145,   145,
     145,   145,   145,   145,   145,   145,   145,   145,   145,   145,
     145,   145,   145,   145,   145,   145,   145,   145,   145,   145,
     145,   146,   146,   147,   147,   148,   148,   149,   149,   150,
     150,   150,   151,   151,   152,   152,   152,   152,   153,   153,
     154,   154,   155,   155,   156,   156,   156,   157,   157,   157,
     158
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     2,     3,     3,     1,     2,     4,     3,
       2,     1,     1,     6,     6,     6,     6,     3,    13,    12,
      12,     8,     3,     2,     2,     7,    13,     9,     4,     4,
       1,     4,     1,     1,     1,     3,     3,     3,     2,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     2,     2,
       5,     3,     5,     1,     3,     2,     4,     4,     1,     4,
       4,     3,     1,     2,     2,     3,     7,     8,     6,     7,
       6,     4,     4,     5,     4,     1,     4,     1,     3,     3,
       4,     1,     4,     1,     2,     2,     1,     3,     2,     1,
       3,     1,     1,     9,    10,     4,     7,     9,     9,     7,
       9,     1,     5,     5,     5,     4,     3,     6,     5,     3,
       7,     4,     1,     5,     4,     6,     5,     7,     4,     4,
       1,     3,     2,     5,     3,     3,     2,     3,     3,     3,
       3,     3,     1,     3,     3,     6,     4,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       1,     1,     1,     2,     3,     2,     3,     2,     1,     1,
       3,     2,     2,     0,     2,     2,     4,     2,     3,     3,
       1,     3,     1,     3,     3,     1,     3,     1,     1,     0,
       1,     0,     1,     1,     2,     2,     0,     2,     3,     1,
       3,     1,     1,     1,     2,     3,     3,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     3,     1,     3,     3,     1,     5,
       7,     3,     3,     4,     4,     5,     5,     3,     3,     4,
       1,     3,     3,     1,     2,     0,     2,     0,     2,     0,
       3,     2,     3,     4,     1,     3,     5,     6,     1,     2,
       1,     2,     0,     5,     0,     3,     5,     0,     2,     5,
       1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint16 yydefact[] =
{
       0,   158,     6,   159,    33,    34,    11,    12,    62,    53,
       0,   120,     0,     0,     0,   132,     0,     0,     0,     0,
       0,     0,     0,     0,    75,     0,     0,   101,     0,     0,
       0,     0,    91,    92,     0,   112,     0,    81,     0,    83,
      77,   150,    58,   152,   186,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    30,     0,     0,     0,     2,
       0,    32,   151,     7,    10,     0,     0,     0,     0,    55,
     161,     0,    49,    48,   120,     0,     0,    38,   126,   122,
     182,     0,     0,     0,     0,     0,     0,     0,     0,   257,
       0,     0,   237,   183,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    85,     0,     0,     0,     0,   191,
     193,   189,   207,     0,   197,   198,   200,   199,   203,   201,
     215,     0,   206,   202,   204,   210,   211,   212,   205,   209,
     208,   186,   186,   218,   252,     0,   186,   186,   186,   213,
     186,     0,   192,   230,     0,   186,     0,   167,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   153,    24,    32,
       0,     0,     0,     0,     0,    23,     0,   162,     0,     0,
     155,     0,    88,     1,     3,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    63,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    64,
       0,     0,     0,     0,     0,   157,     0,     0,     9,   160,
      35,    61,     0,     0,   137,   121,     0,   133,   139,   138,
     140,     0,     0,     0,     0,     0,   159,     0,     0,     0,
     180,     0,     0,     0,   235,   260,   106,   184,   185,     0,
       0,     0,     0,     0,   163,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   233,     0,    62,   152,     0,
       0,     0,   186,   186,     0,     0,     0,     0,   194,    65,
     186,   187,     0,   141,   142,   143,   144,     0,    51,     0,
     145,   146,   147,     0,   248,     0,   109,     0,     0,     0,
       0,    22,     0,   134,     0,    90,   148,    54,    79,    78,
       0,    39,    40,    41,    42,    44,    43,    36,    37,    45,
     128,   129,   130,    46,   127,   125,   124,    32,     0,     4,
       5,    47,    17,   131,   149,   154,   156,     0,     8,    59,
      60,     0,     0,     0,    71,    72,    74,     0,   118,     0,
       0,   119,    95,     0,   238,     0,     0,     0,     0,     0,
       0,   105,     0,     0,   164,     0,   165,   242,     0,   111,
     163,    80,    82,    76,    56,    57,   195,     0,   186,   234,
     227,   190,   214,   217,   196,   186,   228,     0,     0,   186,
     186,   221,   222,   216,   188,   166,   148,   148,     0,   249,
       0,     0,     0,   114,     0,     0,    31,    29,     0,   136,
       0,     0,    28,     0,     0,     0,   239,     0,   250,     0,
      73,   244,     0,     0,   254,     0,   180,     0,     0,   235,
     236,     0,     0,     0,     0,   104,   108,   169,   178,   159,
     175,   177,     0,   172,     0,   170,   243,     0,     0,     0,
       0,     0,   186,   229,   223,   224,     0,     0,    52,    50,
      28,   102,   103,   113,     0,   116,     0,   123,     0,     0,
       0,   135,     0,   239,     0,    68,     0,     0,     0,   251,
       0,    70,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   107,     0,     0,     0,   168,     0,     0,     0,
       0,   232,   231,   186,   219,     0,   226,   225,     0,   115,
       0,    15,    16,    14,    13,    69,     0,   241,     0,     0,
     245,     0,     0,   117,    53,   259,     0,    96,     0,     0,
      99,     0,   173,   176,   174,   171,   110,     0,    66,     0,
     253,     0,     0,     0,    25,   240,     0,     0,     0,   255,
       0,     0,     0,     0,    67,   220,     0,     0,     0,    21,
       0,     0,   246,     0,   254,     0,    93,    97,    98,   100,
       0,     0,     0,     0,    27,   247,   256,    94,     0,     0,
       0,     0,     0,     0,     0,     0,    20,     0,    19,     0,
      18,    26
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    68,    69,   103,    71,   253,    72,   364,   254,   444,
     445,   241,   413,   104,   154,   155,   380,   266,   358,   244,
     475,   109,   420,   421,   422,   272,   484,   239,   248
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -522
static const yytype_int16 yypact[] =
{
    7209,  -522,    40,    11,  -522,    91,  -522,  -522,    70,  -522,
    7209,    89,  7209,  7209,  7209,  -522,  7334,  7209,   118,  7209,
    4334,  4334,  4459,  7459,   106,   111,  7584,    20,  7209,   128,
    7209,  7209,  -522,  -522,    54,   133,   135,   136,  1071,   138,
     140,  -522,    83,  -522,  8353,   -11,  4584,  4709,  4834,  4959,
    7209,  7209,  5084,  5209,  5334,  7209,  7209,  7209,  7209,  7209,
    7209,  5459,    28,  7209,  8353,   120,  7209,  1199,   164,  -522,
    8435,   151,  -522,   -19,  -522,  7209,  7209,  5584,  7209,   545,
    -522,  8623,   545,   545,    69,  1581,  8650, 11080,  -522, 11139,
   10581,   141,   142,   152,  7209,  8508,  7209,  7209,  7709,  8471,
      89,  5709,    16, 10581,  2084,  7209,  8696,  2209,  2334,   112,
    7209,  5834,  7209,  1325, 10581,  7209,  7209,  7209,  7209,    53,
    -522,  -522,  -522,    28,  -522,  -522,  -522,  -522,  -522,  -522,
    -522,  7834,  -522,  -522,  -522,  -522,  -522,  -522,  -522,  -522,
    -522,  8353,  8353,  -522,  -522,    28,  8353,  8353,  8353,  -522,
    8353,  8353,  -522,  -522,   121,  8209,  7209,  -522,   119,   125,
     124,   123,  8769, 10640,   122,   116,   115, 10581, 10581,   206,
    2459,   131,  8811,  2584,  2459,  -522,   213,    64,   109,  7209,
      44,  5959, 11080,  -522,  -522,  7209,  7209,  5459,  5834,  7209,
    7209,  7209,  7209,  7209,  7209,  -522,  7209,  7209,  7209,  7209,
    7209,  7209,  7209,  7209,  7209,  7209,  7209,  6084,  7209,  -522,
      28,  7209,  7209,  7209,  7209,  -522,  5834,  7209,  -522,    64,
   10758,  -522,  8926,  8953,   215,  -522,   943,  -522,  -522,  -522,
    -522,  8986,  5459,  7209,  9101,  9134,    18,   234,  7209,   197,
    9175,   210,  7209,  7209,    66,  -522,  -522, 10581,  -522,  9249,
    7209,  7209,  7209,  6209,   175,  2709,   240,  9290,   241,  9323,
    9364,  9438,  9479,  9512,  8353,  -522,     5,   137,   232,  9553,
     220,   217,  8281,  8353,   196,   198,    29,   177,  -522,  -522,
    8353,  -522,  9627,  -522,  -522,  -522,  -522,  7209,  -522,  7209,
    -522,  -522,  -522,  5834,  -522,  6334,  -522,  7209,  7209,  6459,
    6584,  -522,  5834,  -522,  9668,  -522, 11080, 10758, 11112,  -522,
     252, 10849, 10849,   212,   419, 10908,   211,   687, 11139, 10790,
     439,  9784,   717,   295,    64,   247, 11139,    17,  9701,  -522,
    -522, 11080,    86, 10699, 10817, 10876,    44,   254,   211,  -522,
    -522,  5834,  1453,   781,   222, 10581,  -522,   199,  -522,   218,
    2834,  -522,  -522,  5709, 10581,  9742,  7209,  7209,   117,  5459,
    2959, 10581,  9816,    -2,  -522,  3084,   200,  -522,  7209,  -522,
     175,  -522,  -522,  -522,  -522,  -522,  -522,  7209,  8353,  -522,
    -522,  -522,  -522,  -522,  -522,  8353,  -522,   -15,    78,  8353,
    8353,  -522,  -522,  -522,  -522,  -522, 10994, 11053,   259,  -522,
    2459,  2459,  2459,  -522,  6459,   262,  -522,  -522,  7209,  -522,
    5834,  7209,   266,   268,  5834,  9858,   -44,  7209,  -522,  1709,
    -522,  -522,  7209,    54,    34,  7209, 10581,   246,  7209,  9931,
   10581,  7209,  7209,  7209,  9990,  -522,  -522,  -522,  -522,    24,
     245,    89,    -2,  -522,   250,   282,  -522, 10049,  3209, 10109,
     216,    -8,  8353,  -522,  -522,  -522,   221,   226,  -522,  -522,
      13,  -522,  -522,  -522,  2459,  -522,   270, 10581,   287, 10168,
    5459,  -522,  8623,   -44,   255,  -522,  5459, 10227,  7209,  -522,
    2459,  -522,   304,   263,   233,  1834,  6709,  3334,   249, 10286,
    3459, 10345,  -522,  7209,  7209,   302,  -522,    -2,    54,  7209,
    3584,  -522,  -522,  8353,  -522,   236,  -522,  -522,  7959,  -522,
    6834,  -522,  -522, 10758,  -522,  -522, 10404,  -522,  7209, 10463,
    -522,   276,    54,  -522,   234,  -522,   316,  -522,  7209,  7209,
    -522,  7209, 10581, 10935,  -522,  -522,  -522,  3709,  -522,   258,
    -522,  6959,  1959,  8084, 11139,  -522,   781,  7209,    54,  -522,
    5459,  3834,  3959,  4084,  -522,  -522,  7209,  4209,  7209,  -522,
    7084,  2459,  -522,   781,    34, 10522,  -522,  -522,  -522,  -522,
     227,  7209,   229,  7209,  -522,  -522,  -522,  -522,  7209,   235,
    7209,   237,  2459,  7209,  2459,  7209,  -522,  2459,  -522,  2459,
    -522,  -522
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -522,  -522,   144,     0,   385,  -250,    -4,  -522,  -522,  -152,
     -95,  -350,   113,   184,   248,   -37,  -522,  -522,   -81,  -522,
    -124,    79,  -521,   114,  -387,  -522,  -212,  -131,  -522
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -259
static const yytype_int16 yytable[] =
{
      70,   438,   439,   427,   366,   217,   156,   418,   440,     3,
      79,   441,    81,    82,    83,   442,    86,    87,   100,    89,
      90,    90,    90,    95,     3,   562,    99,   178,   470,   476,
     106,   377,     3,   100,   216,   452,   348,   101,   114,   493,
     153,   100,   575,   503,    73,    74,    90,    90,    90,    90,
     162,   163,    90,    90,    90,   167,   168,    75,    90,   172,
     153,   188,   453,   177,    75,   474,   180,   182,   194,   504,
      75,   195,   157,    80,   242,   219,   220,   222,   223,   203,
     204,   188,   206,   243,   408,   114,   476,    77,   482,   483,
     391,   195,   378,    80,   231,   225,   234,   235,    81,   392,
     117,   240,   206,   410,   247,   249,    76,   247,    78,   218,
     257,    90,   259,   168,   278,   260,   261,   262,   263,    63,
     448,   118,    88,    96,   411,   356,    63,   357,    97,   209,
     108,   269,   508,    91,    92,    93,   526,   153,   153,   454,
     175,    63,   153,   153,   153,   105,   153,   153,   455,    63,
     110,   153,   111,   112,    77,   115,   282,   116,   179,   158,
     159,   160,   161,   381,   183,   164,   165,   166,   216,  -259,
     247,   171,   215,   247,   247,    78,   431,   432,   264,   304,
     433,    81,   228,   229,   256,   306,   307,   308,    90,   311,
     312,   313,   314,   315,   316,   230,   317,   318,   319,   320,
     321,   322,   323,   324,   325,   326,   328,    70,   331,   279,
     283,   333,   334,   335,   336,   107,    90,   338,   284,   285,
     286,   291,   292,   293,   258,   290,    86,   376,   188,   188,
     302,   303,   341,   345,   192,   194,   194,   349,   195,   195,
     297,   170,   354,   355,   173,   174,   203,   203,   204,   206,
     206,   361,   362,   351,   353,   247,   252,   368,   382,   370,
     153,   384,   389,   390,   188,   393,   309,   385,   153,   153,
     407,   194,   412,   423,   195,   424,   153,   460,   425,   363,
     466,   470,   203,  -259,   296,   206,   471,   396,   301,   397,
     486,   494,   255,    90,   496,   168,   497,   209,   506,   510,
     168,   310,    90,   507,   502,   511,   479,   521,   522,   528,
     523,   344,   188,   540,   189,   190,   191,   192,   193,   194,
     534,   548,   195,   196,   197,   198,   199,   200,   201,   337,
     203,   204,   205,   206,   550,   555,   578,   214,   580,   215,
     215,    90,   415,   247,   583,   535,   585,   495,   488,   515,
     247,   330,   576,   426,   525,     0,   429,   430,     0,   434,
     247,     0,     0,     0,     0,   247,     0,     0,   447,     0,
       0,     0,     0,     0,   153,   215,     0,   449,     0,     0,
     209,   153,     0,     0,     0,   153,   153,     0,     0,   270,
     271,     0,     0,     0,   274,   275,   276,     0,   277,     0,
     247,   247,   247,   281,     0,     0,   398,     0,   467,     0,
      90,   469,   102,   403,   472,   405,   343,   477,     0,   168,
     214,   416,   350,   215,     0,     0,     0,     0,     0,   152,
       0,   489,     0,   491,   360,     0,   188,   365,   435,     0,
       0,   169,     0,   194,     0,     0,   195,   176,   153,   152,
       0,     0,     0,     0,   203,   204,   188,   206,   189,   190,
     191,   192,   193,   194,   247,     0,   195,     0,     0,   198,
     513,   200,   201,     0,   203,   204,   516,   206,   519,     0,
     247,   400,   401,   402,     0,   247,   426,   247,     0,     0,
     247,     0,     0,   532,   533,     0,     0,     0,   169,   153,
     247,     0,   481,     0,   209,     0,     0,     0,   265,     0,
     544,     0,     0,     0,   461,   462,   463,     0,   465,     0,
     387,   388,     0,   468,   209,     0,   152,   152,   394,     0,
     273,   152,   152,   152,     0,   152,   152,   247,     0,     0,
     152,   168,   247,     0,   214,     0,   247,   215,     0,   514,
     565,   247,   247,   247,     0,   517,    90,   168,    90,     0,
     168,   247,   188,   247,   214,     0,     0,   215,     0,   194,
       0,    90,   195,    90,     0,     0,     0,   536,   509,     0,
     203,   204,   247,   206,   247,     0,     0,   247,   464,   247,
     327,     0,     0,     0,   520,   332,     0,     0,     0,     0,
       0,   549,     0,     0,     0,     0,   480,     0,     0,   485,
       0,     0,   487,     0,     0,     0,   490,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   450,   564,     0,   566,
     209,     0,   500,   451,     0,     0,     0,   456,   457,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   152,
       0,   379,     0,     0,     0,     0,   559,   152,   152,     0,
       0,     0,     0,     0,     0,   152,     0,     0,     0,   570,
       0,   572,     0,   215,     0,   574,     0,     0,     0,     0,
     169,     0,     0,   537,   579,   169,   581,     0,     0,     0,
       0,     0,   542,     0,     0,     0,   586,     0,   588,     0,
     505,   590,   546,   591,   188,     0,   189,   190,   191,   192,
     193,   194,   551,   552,   195,   553,     0,   198,   199,   200,
     201,     0,   203,   204,     0,   206,     0,   561,     0,     0,
       0,   563,     0,     0,   188,     0,   189,   190,   191,   192,
     193,   194,     0,     0,   195,     0,     0,   198,   443,   200,
    -259,   539,   203,   204,     0,   206,     0,     0,     0,     0,
       0,     0,   582,   152,   584,     0,     0,   587,     0,   589,
     152,     0,   209,     0,   152,   152,     0,     0,     0,     0,
       0,     0,     1,     0,     2,     3,     4,     5,     6,     0,
       7,     8,     9,    10,    11,     0,     0,     0,    12,     0,
      13,    14,   209,     0,   169,     0,    15,    16,     0,     0,
      17,     0,   214,     0,     0,   215,     0,    18,    19,     0,
      20,    21,     0,    22,     0,   245,     0,   443,    23,     0,
     417,     0,   418,    24,    25,    26,     0,   152,    27,     0,
       0,    28,   214,    29,    30,   215,    31,     0,     0,     0,
      32,    33,    34,     0,    35,    36,    37,     0,   294,    38,
       0,    39,     0,    40,    41,    42,     0,    43,    44,     0,
      45,    46,     0,    47,     0,    48,     0,    49,     0,    50,
       0,    51,   443,    52,     0,    53,     0,    54,   152,    55,
     419,     0,    57,    58,    59,     0,    60,    61,     0,    62,
       0,     0,    63,    64,    65,     0,    66,     0,    67,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   169,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   169,     0,     1,   169,     2,     3,     4,     5,
       6,   -87,     7,     8,     9,    10,    84,   -87,   -87,   -87,
      12,   -87,    13,    14,   -87,   -87,   -87,   -87,    15,    16,
     -87,   -87,    17,   -87,   -87,   -87,   -87,   -87,   -87,    18,
      19,   -87,    20,    21,   -87,    22,   -87,   -87,   -87,     0,
      23,   -87,   -87,   -87,   -87,    24,    25,    26,   -87,   -87,
      27,   -87,   -87,    28,   -87,    29,    30,   -87,    31,   -87,
     -87,   -87,    32,    33,    34,     0,    35,    36,    37,     0,
     -87,    85,     0,    39,     0,    40,    41,    42,   -87,    43,
      44,   -87,    45,    46,   -87,    47,   -87,    48,   -87,    49,
     -87,    50,   -87,    51,   -87,    52,   -87,    53,   -87,    54,
     -87,    55,    56,   -87,    57,    58,    59,   -87,    60,    61,
     -87,    62,     0,   -87,    63,    64,    65,   -87,    66,     0,
      67,   -87,     1,     0,     2,     3,     4,     5,     6,   -86,
       7,     8,     9,    10,    11,   -86,   -86,   -86,    12,   -86,
      13,    14,   -86,   -86,   -86,   -86,    15,    16,   -86,   -86,
      17,   -86,   -86,   -86,   -86,   -86,   -86,    18,    19,   -86,
      20,    21,   -86,    22,   -86,   -86,   -86,     0,    23,   -86,
     -86,   -86,   -86,    24,    25,    26,   -86,   -86,    27,   -86,
     -86,    28,   -86,    29,    30,   -86,    31,   -86,   -86,   -86,
      32,    33,    34,     0,    35,    36,    37,     0,   -86,     0,
       0,    39,     0,    40,    41,    42,   -86,    43,    44,   -86,
      45,    46,   -86,    47,   -86,    48,   -86,    49,   -86,    50,
     -86,    51,   -86,    52,   -86,    53,   -86,    54,   -86,    55,
     -86,   -86,    57,    58,    59,   -86,    60,    61,   -86,    62,
       0,   -86,    63,    64,    65,   -86,    66,     0,    67,   -86,
       1,     0,     2,     3,     4,     5,     6,   -89,     7,     8,
       9,    10,    11,   -89,   -89,   -89,   181,   -89,    13,    14,
     -89,   -89,   -89,   -89,    15,    16,   -89,   -89,    17,   -89,
     -89,   -89,   -89,   -89,   -89,    18,    19,   -89,    20,    21,
     -89,    22,   -89,   -89,   -89,     0,    23,   -89,   -89,   -89,
     -89,    24,    25,    26,   -89,   -89,    27,   -89,   -89,    28,
     -89,    29,    30,   -89,    31,   -89,   -89,   -89,    32,    33,
      34,     0,    35,    36,    37,     0,   -89,   -89,     0,    39,
       0,    40,    41,    42,   -89,    43,    44,   -89,    45,    46,
     -89,    47,   -89,    48,   -89,    49,   -89,    50,   -89,    51,
     -89,    52,   -89,    53,   -89,    54,   -89,    55,   -89,   -89,
      57,    58,    59,   -89,    60,    61,   -89,    62,     0,   -89,
      63,    64,    65,   -89,    66,     0,     1,   -89,     2,     3,
       4,     5,     6,   -84,     7,     8,     9,    10,    11,   -84,
     -84,   -84,    12,   -84,    13,    14,   -84,   -84,   -84,   -84,
      15,    16,   -84,   -84,    17,   -84,   -84,   -84,   -84,   -84,
     -84,    18,    19,   -84,    20,    21,   -84,    22,   -84,   -84,
     -84,     0,    23,   -84,   -84,   -84,   -84,    24,    25,    26,
     -84,   -84,    27,   -84,   -84,    28,   -84,    29,    30,   -84,
      31,   -84,   -84,   -84,    32,    33,    34,     0,    35,    36,
      37,     0,   -84,    38,     0,    39,     0,    40,    41,    42,
     -84,    43,    44,   -84,    45,    46,   -84,    47,   -84,    48,
     -84,    49,   -84,    50,   -84,    51,   -84,    52,   -84,    53,
     -84,    54,   -84,    55,     0,   -84,    57,    58,    59,   -84,
      60,    61,   -84,    62,     0,   -84,    63,    64,    65,   -84,
      66,     0,    67,   -84,     1,     0,     2,     3,     4,     5,
       6,     0,     7,     8,     9,    10,    11,  -137,  -137,  -137,
     414,     0,    13,    14,  -137,  -137,  -137,  -137,    15,    16,
    -137,  -137,    17,  -137,  -137,  -137,  -137,  -137,  -137,    18,
      19,  -137,    20,    21,     0,    22,     0,     0,  -137,     0,
      23,     0,     0,  -137,     0,    24,    25,    26,     0,     0,
      27,     0,     0,    28,     0,    29,    30,     0,    31,     0,
       0,     0,    32,    33,    34,     0,    35,    36,    37,   108,
       0,    38,     0,    39,     0,    40,    41,    42,  -137,    43,
      44,     0,    45,    46,     0,    47,     0,    48,     0,    49,
       0,    50,     0,    51,     0,    52,     0,    53,     0,    54,
       0,    55,    56,     0,    57,    58,    59,  -137,    60,    61,
    -137,    62,     0,  -137,    63,    64,    65,  -137,    66,     0,
      67,  -137,     1,     0,     2,     3,     4,     5,     6,     0,
       7,     8,     9,    10,    11,   -86,   -86,   -86,    12,     0,
      13,    14,   -86,   -86,   -86,   -86,    15,   226,   -86,   -86,
      17,   -86,   -86,   -86,   -86,   -86,   -86,    18,    19,   -86,
      20,    21,     0,    22,     0,     0,   -86,     0,    23,     0,
       0,     0,     0,    24,    25,    26,     0,     0,    27,     0,
       0,    28,     0,    29,    30,     0,    31,     0,     0,     0,
      32,    33,    34,     0,    35,    36,    37,     0,     0,    38,
       0,    39,     0,    40,    41,    42,   -86,    43,    44,     0,
      45,    46,     0,    47,     0,    48,     0,    49,     0,    50,
       0,    51,     0,    52,     0,    53,     0,    54,     0,    55,
     113,     0,    57,    58,    59,   -86,    60,    61,   -86,    62,
       0,   -86,    63,    64,    65,   -86,    66,     0,    67,   -86,
       1,     0,     2,     3,     4,     5,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,     0,    13,    14,
       0,     0,     0,     0,    15,    16,     0,     0,    17,     0,
       0,     0,     0,     0,     0,    18,    19,     0,    20,    21,
       0,    22,     0,     0,     0,     0,    23,     0,   478,     0,
     479,    24,    25,    26,     0,     0,    27,     0,     0,    28,
       0,    29,    30,     0,    31,     0,     0,     0,    32,    33,
      34,     0,    35,    36,    37,     0,   399,    38,     0,    39,
       0,    40,    41,    42,     0,    43,    44,     0,    45,    46,
       0,    47,     0,    48,     0,    49,     0,    50,     0,    51,
       0,    52,     0,    53,     0,    54,     0,    55,    56,     0,
      57,    58,    59,     0,    60,    61,     0,    62,     0,     0,
      63,    64,    65,     0,    66,     1,    67,     2,     3,     4,
       5,     6,     0,     7,     8,   524,    10,    11,     0,     0,
       0,    12,     0,    13,    14,     0,     0,     0,     0,    15,
      16,     0,     0,    17,     0,     0,     0,     0,     0,     0,
      18,    19,     0,    20,    21,     0,    22,     0,   245,     0,
       0,    23,     0,     0,     0,     0,    24,    25,    26,   238,
    -257,    27,     0,     0,    28,     0,    29,    30,     0,    31,
       0,     0,     0,    32,    33,    34,     0,    35,    36,    37,
       0,     0,    38,     0,    39,     0,    40,    41,    42,     0,
      43,    44,     0,    45,    46,     0,    47,     0,    48,     0,
      49,     0,    50,     0,    51,     0,    52,     0,    53,     0,
      54,     0,    55,    56,     0,    57,    58,    59,     0,    60,
      61,     0,    62,     0,     0,    63,    64,    65,     0,    66,
       1,    67,     2,     3,     4,     5,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,     0,    13,    14,
       0,     0,     0,     0,    15,    16,     0,     0,    17,     0,
       0,     0,     0,     0,     0,    18,    19,     0,    20,    21,
       0,    22,     0,   245,     0,     0,    23,     0,     0,     0,
       0,    24,    25,    26,     0,     0,    27,     0,     0,    28,
       0,    29,    30,     0,    31,     0,     0,     0,    32,    33,
      34,     0,    35,    36,    37,     0,   294,    38,     0,    39,
       0,    40,    41,    42,     0,    43,    44,     0,    45,    46,
       0,    47,     0,    48,     0,    49,     0,    50,     0,    51,
       0,    52,     0,    53,     0,    54,     0,    55,   557,   558,
      57,    58,    59,     0,    60,    61,     0,    62,     0,     0,
      63,    64,    65,     0,    66,     1,    67,     2,     3,     4,
       5,     6,     0,     7,     8,     9,    10,    11,     0,     0,
       0,    12,     0,    13,    14,     0,     0,     0,     0,    15,
      16,     0,     0,    17,     0,     0,     0,     0,     0,     0,
      18,    19,     0,    20,    21,     0,    22,     0,   245,     0,
       0,    23,     0,     0,     0,     0,    24,    25,    26,     0,
       0,    27,     0,     0,    28,     0,    29,    30,     0,    31,
       0,     0,     0,    32,    33,    34,     0,    35,    36,    37,
       0,   246,    38,     0,    39,     0,    40,    41,    42,     0,
      43,    44,     0,    45,    46,     0,    47,     0,    48,     0,
      49,     0,    50,     0,    51,     0,    52,     0,    53,     0,
      54,     0,    55,    56,     0,    57,    58,    59,     0,    60,
      61,     0,    62,     0,     0,    63,    64,    65,     0,    66,
       1,    67,     2,     3,     4,     5,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,     0,    13,    14,
       0,     0,     0,     0,    15,    16,     0,     0,    17,     0,
       0,     0,     0,     0,     0,    18,    19,     0,    20,    21,
       0,    22,     0,   245,     0,     0,    23,     0,     0,     0,
       0,    24,    25,    26,     0,     0,    27,     0,     0,    28,
       0,    29,    30,     0,    31,   251,     0,     0,    32,    33,
      34,     0,    35,    36,    37,     0,     0,    38,     0,    39,
       0,    40,    41,    42,     0,    43,    44,     0,    45,    46,
       0,    47,     0,    48,     0,    49,     0,    50,     0,    51,
       0,    52,     0,    53,     0,    54,     0,    55,    56,     0,
      57,    58,    59,     0,    60,    61,     0,    62,     0,     0,
      63,    64,    65,     0,    66,     1,    67,     2,     3,     4,
       5,     6,     0,     7,     8,     9,    10,    11,     0,     0,
       0,    12,     0,    13,    14,     0,     0,     0,     0,    15,
      16,     0,     0,    17,     0,     0,     0,     0,     0,     0,
      18,    19,     0,    20,    21,     0,    22,     0,     0,     0,
       0,    23,     0,     0,     0,     0,    24,    25,    26,     0,
       0,    27,     0,     0,    28,     0,    29,    30,     0,    31,
       0,     0,     0,    32,    33,    34,     0,    35,    36,    37,
       0,     0,    38,  -163,    39,   252,    40,    41,    42,     0,
      43,    44,     0,    45,    46,     0,    47,     0,    48,     0,
      49,     0,    50,     0,    51,     0,    52,     0,    53,     0,
      54,     0,    55,    56,     0,    57,    58,    59,     0,    60,
      61,     0,    62,     0,     0,    63,    64,    65,     0,    66,
       1,    67,     2,     3,     4,     5,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,     0,    13,    14,
       0,     0,     0,     0,    15,    16,     0,     0,    17,     0,
       0,     0,     0,     0,     0,    18,    19,     0,    20,    21,
       0,    22,     0,   245,     0,     0,    23,     0,     0,     0,
       0,    24,    25,    26,     0,     0,    27,     0,     0,    28,
       0,    29,    30,     0,    31,     0,     0,     0,    32,    33,
      34,     0,    35,    36,    37,     0,   294,    38,     0,    39,
       0,    40,    41,    42,     0,    43,    44,     0,    45,    46,
       0,    47,     0,    48,     0,    49,     0,    50,     0,    51,
       0,    52,     0,    53,     0,    54,     0,    55,   295,     0,
      57,    58,    59,     0,    60,    61,     0,    62,     0,     0,
      63,    64,    65,     0,    66,     1,    67,     2,     3,     4,
       5,     6,     0,     7,     8,     9,    10,    11,     0,     0,
       0,    12,     0,    13,    14,     0,     0,     0,     0,    15,
      16,     0,     0,    17,     0,     0,     0,     0,     0,     0,
      18,    19,     0,    20,    21,     0,    22,     0,   245,     0,
       0,    23,     0,     0,     0,   299,    24,    25,    26,     0,
       0,    27,     0,     0,    28,     0,    29,    30,     0,    31,
       0,     0,     0,    32,    33,    34,     0,    35,    36,    37,
       0,     0,    38,     0,    39,     0,    40,    41,    42,     0,
      43,    44,     0,    45,    46,     0,    47,     0,    48,     0,
      49,     0,    50,     0,    51,     0,    52,     0,    53,     0,
      54,     0,    55,   300,     0,    57,    58,    59,     0,    60,
      61,     0,    62,     0,     0,    63,    64,    65,     0,    66,
       1,    67,     2,     3,     4,     5,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,     0,    13,    14,
       0,     0,     0,     0,    15,    16,     0,     0,    17,     0,
       0,     0,     0,     0,     0,    18,    19,     0,    20,    21,
       0,    22,     0,   245,     0,     0,    23,     0,     0,     0,
       0,    24,    25,    26,     0,     0,    27,     0,     0,    28,
       0,    29,    30,     0,    31,     0,     0,     0,    32,    33,
      34,     0,    35,    36,    37,     0,   367,    38,     0,    39,
       0,    40,    41,    42,     0,    43,    44,     0,    45,    46,
       0,    47,     0,    48,     0,    49,     0,    50,     0,    51,
       0,    52,     0,    53,     0,    54,     0,    55,    56,     0,
      57,    58,    59,     0,    60,    61,     0,    62,     0,     0,
      63,    64,    65,     0,    66,     1,    67,     2,     3,     4,
       5,     6,     0,     7,     8,     9,    10,    11,     0,     0,
       0,    12,     0,    13,    14,     0,     0,     0,     0,    15,
      16,     0,     0,    17,     0,     0,     0,     0,     0,     0,
      18,    19,     0,    20,    21,     0,    22,     0,   245,     0,
       0,    23,     0,     0,     0,     0,    24,    25,    26,     0,
    -258,    27,     0,     0,    28,     0,    29,    30,     0,    31,
       0,     0,     0,    32,    33,    34,     0,    35,    36,    37,
       0,     0,    38,     0,    39,     0,    40,    41,    42,     0,
      43,    44,     0,    45,    46,     0,    47,     0,    48,     0,
      49,     0,    50,     0,    51,     0,    52,     0,    53,     0,
      54,     0,    55,    56,     0,    57,    58,    59,     0,    60,
      61,     0,    62,     0,     0,    63,    64,    65,     0,    66,
       1,    67,     2,     3,     4,     5,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,     0,    13,    14,
       0,     0,     0,     0,    15,    16,     0,     0,    17,     0,
       0,     0,     0,     0,     0,    18,    19,     0,    20,    21,
       0,    22,     0,   245,     0,     0,    23,     0,     0,     0,
       0,    24,    25,    26,     0,     0,    27,     0,     0,    28,
       0,    29,    30,     0,    31,     0,     0,     0,    32,    33,
      34,     0,    35,    36,    37,     0,   436,    38,     0,    39,
       0,    40,    41,    42,     0,    43,    44,     0,    45,    46,
       0,    47,     0,    48,     0,    49,     0,    50,     0,    51,
       0,    52,     0,    53,     0,    54,     0,    55,    56,     0,
      57,    58,    59,     0,    60,    61,     0,    62,     0,     0,
      63,    64,    65,     0,    66,     1,    67,     2,     3,     4,
       5,     6,     0,     7,     8,     9,    10,    11,     0,     0,
       0,    12,     0,    13,    14,     0,     0,     0,     0,    15,
      16,     0,     0,    17,     0,     0,     0,     0,     0,     0,
      18,    19,     0,    20,    21,     0,    22,     0,   245,     0,
       0,    23,     0,     0,     0,     0,    24,    25,    26,     0,
       0,    27,     0,     0,    28,     0,    29,    30,     0,    31,
       0,     0,     0,    32,    33,    34,     0,    35,    36,    37,
       0,   446,    38,     0,    39,     0,    40,    41,    42,     0,
      43,    44,     0,    45,    46,     0,    47,     0,    48,     0,
      49,     0,    50,     0,    51,     0,    52,     0,    53,     0,
      54,     0,    55,    56,     0,    57,    58,    59,     0,    60,
      61,     0,    62,     0,     0,    63,    64,    65,     0,    66,
       1,    67,     2,     3,     4,     5,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,     0,    13,    14,
       0,     0,     0,     0,    15,    16,     0,     0,    17,     0,
       0,     0,     0,     0,     0,    18,    19,     0,    20,    21,
       0,    22,     0,     0,     0,     0,    23,     0,     0,     0,
       0,    24,    25,    26,     0,     0,    27,     0,     0,    28,
       0,    29,    30,     0,    31,     0,     0,     0,    32,    33,
      34,     0,    35,    36,    37,   499,     0,    38,   363,    39,
       0,    40,    41,    42,     0,    43,    44,     0,    45,    46,
       0,    47,     0,    48,     0,    49,     0,    50,     0,    51,
       0,    52,     0,    53,     0,    54,     0,    55,    56,     0,
      57,    58,    59,     0,    60,    61,     0,    62,     0,     0,
      63,    64,    65,     0,    66,     1,    67,     2,     3,     4,
       5,     6,     0,     7,     8,     9,    10,    11,     0,     0,
       0,    12,     0,    13,    14,     0,     0,     0,     0,    15,
      16,     0,     0,    17,     0,     0,     0,     0,     0,     0,
      18,    19,     0,    20,    21,     0,    22,     0,   245,     0,
       0,    23,     0,     0,     0,     0,    24,    25,    26,     0,
       0,    27,     0,     0,    28,     0,    29,    30,     0,    31,
       0,     0,     0,    32,    33,    34,     0,    35,    36,    37,
       0,   527,    38,     0,    39,     0,    40,    41,    42,     0,
      43,    44,     0,    45,    46,     0,    47,     0,    48,     0,
      49,     0,    50,     0,    51,     0,    52,     0,    53,     0,
      54,     0,    55,    56,     0,    57,    58,    59,     0,    60,
      61,     0,    62,     0,     0,    63,    64,    65,     0,    66,
       1,    67,     2,     3,     4,     5,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,     0,    13,    14,
       0,     0,     0,     0,    15,    16,     0,     0,    17,     0,
       0,     0,     0,     0,     0,    18,    19,     0,    20,    21,
       0,    22,     0,   245,     0,     0,    23,     0,     0,     0,
       0,    24,    25,    26,     0,     0,    27,     0,     0,    28,
       0,    29,    30,     0,    31,     0,     0,     0,    32,    33,
      34,     0,    35,    36,    37,     0,   530,    38,     0,    39,
       0,    40,    41,    42,     0,    43,    44,     0,    45,    46,
       0,    47,     0,    48,     0,    49,     0,    50,     0,    51,
       0,    52,     0,    53,     0,    54,     0,    55,    56,     0,
      57,    58,    59,     0,    60,    61,     0,    62,     0,     0,
      63,    64,    65,     0,    66,     1,    67,     2,     3,     4,
       5,     6,     0,     7,     8,     9,    10,    11,     0,     0,
       0,    12,     0,    13,    14,     0,     0,     0,     0,    15,
      16,     0,     0,    17,     0,     0,     0,     0,     0,     0,
      18,    19,     0,    20,    21,     0,    22,     0,   245,     0,
       0,    23,     0,     0,     0,     0,    24,    25,    26,     0,
       0,    27,     0,     0,    28,     0,    29,    30,     0,    31,
       0,     0,     0,    32,    33,    34,     0,    35,    36,    37,
       0,   538,    38,     0,    39,     0,    40,    41,    42,     0,
      43,    44,     0,    45,    46,     0,    47,     0,    48,     0,
      49,     0,    50,     0,    51,     0,    52,     0,    53,     0,
      54,     0,    55,    56,     0,    57,    58,    59,     0,    60,
      61,     0,    62,     0,     0,    63,    64,    65,     0,    66,
       1,    67,     2,     3,     4,     5,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,     0,    13,    14,
       0,     0,     0,     0,    15,    16,     0,     0,    17,     0,
       0,     0,     0,     0,     0,    18,    19,     0,    20,    21,
       0,    22,     0,   245,     0,     0,    23,     0,     0,     0,
       0,    24,    25,    26,     0,     0,    27,     0,     0,    28,
       0,    29,    30,     0,    31,     0,     0,     0,    32,    33,
      34,     0,    35,    36,    37,     0,   554,    38,     0,    39,
       0,    40,    41,    42,     0,    43,    44,     0,    45,    46,
       0,    47,     0,    48,     0,    49,     0,    50,     0,    51,
       0,    52,     0,    53,     0,    54,     0,    55,    56,     0,
      57,    58,    59,     0,    60,    61,     0,    62,     0,     0,
      63,    64,    65,     0,    66,     1,    67,     2,     3,     4,
       5,     6,     0,     7,     8,     9,    10,    11,     0,     0,
       0,    12,     0,    13,    14,     0,     0,     0,     0,    15,
      16,     0,     0,    17,     0,     0,     0,     0,     0,     0,
      18,    19,     0,    20,    21,     0,    22,     0,   245,     0,
       0,    23,     0,     0,     0,     0,    24,    25,    26,     0,
       0,    27,     0,     0,    28,     0,    29,    30,     0,    31,
       0,     0,     0,    32,    33,    34,     0,    35,    36,    37,
       0,   567,    38,     0,    39,     0,    40,    41,    42,     0,
      43,    44,     0,    45,    46,     0,    47,     0,    48,     0,
      49,     0,    50,     0,    51,     0,    52,     0,    53,     0,
      54,     0,    55,    56,     0,    57,    58,    59,     0,    60,
      61,     0,    62,     0,     0,    63,    64,    65,     0,    66,
       1,    67,     2,     3,     4,     5,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,     0,    13,    14,
       0,     0,     0,     0,    15,    16,     0,     0,    17,     0,
       0,     0,     0,     0,     0,    18,    19,     0,    20,    21,
       0,    22,     0,   245,     0,     0,    23,     0,     0,     0,
       0,    24,    25,    26,     0,     0,    27,     0,     0,    28,
       0,    29,    30,     0,    31,     0,     0,     0,    32,    33,
      34,     0,    35,    36,    37,     0,   568,    38,     0,    39,
       0,    40,    41,    42,     0,    43,    44,     0,    45,    46,
       0,    47,     0,    48,     0,    49,     0,    50,     0,    51,
       0,    52,     0,    53,     0,    54,     0,    55,    56,     0,
      57,    58,    59,     0,    60,    61,     0,    62,     0,     0,
      63,    64,    65,     0,    66,     1,    67,     2,     3,     4,
       5,     6,     0,     7,     8,     9,    10,    11,     0,     0,
       0,    12,     0,    13,    14,     0,     0,     0,     0,    15,
      16,     0,     0,    17,     0,     0,     0,     0,     0,     0,
      18,    19,     0,    20,    21,     0,    22,     0,   245,     0,
       0,    23,     0,     0,     0,     0,    24,    25,    26,     0,
       0,    27,     0,     0,    28,     0,    29,    30,     0,    31,
       0,     0,     0,    32,    33,    34,     0,    35,    36,    37,
       0,   569,    38,     0,    39,     0,    40,    41,    42,     0,
      43,    44,     0,    45,    46,     0,    47,     0,    48,     0,
      49,     0,    50,     0,    51,     0,    52,     0,    53,     0,
      54,     0,    55,    56,     0,    57,    58,    59,     0,    60,
      61,     0,    62,     0,     0,    63,    64,    65,     0,    66,
       1,    67,     2,     3,     4,     5,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,     0,    13,    14,
       0,     0,     0,     0,    15,    16,     0,     0,    17,     0,
       0,     0,     0,     0,     0,    18,    19,     0,    20,    21,
       0,    22,     0,     0,     0,     0,    23,     0,     0,     0,
       0,    24,    25,    26,     0,     0,    27,     0,     0,    28,
       0,    29,    30,     0,    31,     0,     0,     0,    32,    33,
      34,     0,    35,    36,    37,     0,   399,    38,     0,    39,
       0,    40,    41,    42,     0,    43,    44,     0,    45,    46,
       0,    47,     0,    48,     0,    49,     0,    50,     0,    51,
       0,    52,     0,    53,     0,    54,     0,    55,    56,   571,
      57,    58,    59,     0,    60,    61,     0,    62,     0,     0,
      63,    64,    65,     0,    66,     1,    67,     2,     3,     4,
       5,     6,     0,     7,     8,     9,    10,    11,     0,     0,
       0,    12,     0,    13,    14,     0,     0,     0,     0,    15,
      16,     0,     0,    17,     0,     0,     0,     0,     0,     0,
      18,    19,     0,    20,    21,  -181,    22,     0,     0,     0,
       0,    23,     0,     0,     0,     0,    24,    25,    26,     0,
       0,    27,     0,     0,    28,     0,    29,    30,     0,    31,
       0,     0,     0,    32,    33,    34,     0,    35,    36,    37,
       0,     0,    38,     0,    39,     0,    40,    41,    42,     0,
      43,    44,     0,    45,    46,     0,    47,     0,    48,     0,
      49,     0,    50,     0,    51,     0,    52,     0,    53,     0,
      54,     0,    55,    56,     0,    57,    58,    59,     0,    60,
      61,     0,    62,     0,     0,    63,    64,    65,     0,    66,
       1,    67,     2,     3,     4,     5,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,     0,    13,    14,
       0,     0,     0,     0,    15,    16,     0,     0,    17,     0,
       0,     0,     0,     0,     0,    18,    19,     0,    20,    21,
       0,    22,  -181,     0,     0,     0,    23,     0,     0,     0,
       0,    24,    25,    26,     0,     0,    27,     0,     0,    28,
       0,    29,    30,     0,    31,     0,     0,     0,    32,    33,
      34,     0,    35,    36,    37,     0,     0,    38,     0,    39,
       0,    40,    41,    42,     0,    43,    44,     0,    45,    46,
       0,    47,     0,    48,     0,    49,     0,    50,     0,    51,
       0,    52,     0,    53,     0,    54,     0,    55,    56,     0,
      57,    58,    59,     0,    60,    61,     0,    62,     0,     0,
      63,    64,    65,     0,    66,     1,    67,     2,     3,     4,
       5,     6,     0,     7,     8,     9,    10,    11,     0,     0,
       0,    12,     0,    13,    14,     0,     0,     0,     0,    15,
      16,     0,     0,    17,     0,     0,     0,     0,     0,     0,
      18,    19,     0,    20,    21,     0,    22,     0,     0,     0,
       0,    23,     0,     0,     0,     0,    24,    25,    26,     0,
       0,    27,     0,     0,    28,     0,    29,    30,     0,    31,
       0,     0,     0,    32,    33,    34,     0,    35,    36,    37,
       0,     0,    38,     0,    39,     0,    40,    41,    42,     0,
      43,    44,     0,    45,    46,  -181,    47,     0,    48,     0,
      49,     0,    50,     0,    51,     0,    52,     0,    53,     0,
      54,     0,    55,    56,     0,    57,    58,    59,     0,    60,
      61,     0,    62,     0,     0,    63,    64,    65,     0,    66,
       1,    67,     2,     3,     4,     5,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,     0,    13,    14,
       0,     0,     0,     0,    15,    16,     0,     0,    17,     0,
       0,     0,     0,     0,     0,    18,    19,     0,    20,    21,
       0,    22,     0,     0,     0,     0,    23,     0,     0,     0,
       0,    24,    25,    26,     0,     0,    27,     0,     0,    28,
       0,    29,    30,     0,    31,     0,     0,     0,    32,    33,
      34,     0,    35,    36,    37,     0,     0,    38,     0,    39,
       0,    40,    41,    42,     0,    43,    44,     0,    45,    46,
       0,    47,  -181,    48,     0,    49,     0,    50,     0,    51,
       0,    52,     0,    53,     0,    54,     0,    55,    56,     0,
      57,    58,    59,     0,    60,    61,     0,    62,     0,     0,
      63,    64,    65,     0,    66,     1,    67,     2,     3,     4,
       5,     6,     0,     7,     8,     9,    10,    11,     0,     0,
       0,    12,     0,    13,    14,     0,     0,     0,     0,    15,
      16,     0,     0,    17,     0,     0,     0,     0,     0,     0,
      18,    19,     0,    20,    21,     0,    22,     0,     0,     0,
       0,    23,     0,     0,     0,     0,    24,    25,    26,     0,
       0,    27,     0,     0,    28,     0,    29,    30,     0,    31,
       0,     0,     0,    32,    33,    34,     0,    35,    36,    37,
       0,     0,    38,     0,    39,     0,    40,    41,    42,     0,
      43,    44,     0,    45,    46,     0,    47,     0,    48,  -181,
      49,     0,    50,     0,    51,     0,    52,     0,    53,     0,
      54,     0,    55,    56,     0,    57,    58,    59,     0,    60,
      61,     0,    62,     0,     0,    63,    64,    65,     0,    66,
       1,    67,     2,     3,     4,     5,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,     0,    13,    14,
       0,     0,     0,     0,    15,    16,     0,     0,    17,     0,
       0,     0,     0,     0,     0,    18,    19,     0,    20,    21,
       0,    22,     0,     0,     0,     0,    23,     0,     0,     0,
       0,    24,    25,    26,     0,     0,    27,     0,     0,    28,
       0,    29,    30,     0,    31,     0,     0,     0,    32,    33,
      34,     0,    35,    36,    37,     0,     0,    38,     0,    39,
       0,    40,    41,    42,     0,    43,    44,     0,    45,    46,
       0,    47,     0,    48,     0,    49,  -181,    50,     0,    51,
       0,    52,     0,    53,     0,    54,     0,    55,    56,     0,
      57,    58,    59,     0,    60,    61,     0,    62,     0,     0,
      63,    64,    65,     0,    66,     1,    67,     2,     3,     4,
       5,     6,     0,     7,     8,     9,    10,    11,     0,     0,
       0,    12,     0,    13,    14,     0,     0,     0,     0,    15,
      16,     0,     0,    17,     0,     0,     0,     0,     0,     0,
      18,    19,     0,    20,    21,     0,    22,     0,     0,     0,
       0,    23,     0,     0,     0,     0,    24,    25,    26,     0,
       0,    27,     0,     0,    28,     0,    29,    30,     0,    31,
       0,     0,     0,    32,    33,    34,     0,    35,    36,    37,
       0,     0,    38,     0,    39,     0,    40,    41,    42,     0,
      43,    44,     0,    45,    46,     0,    47,     0,    48,     0,
      49,     0,    50,     0,    51,     0,    52,  -181,    53,     0,
      54,     0,    55,    56,     0,    57,    58,    59,     0,    60,
      61,     0,    62,     0,     0,    63,    64,    65,     0,    66,
       1,    67,     2,     3,     4,     5,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,     0,    13,    14,
       0,     0,     0,     0,    15,    16,     0,     0,    17,     0,
       0,     0,     0,     0,     0,    18,    19,     0,    20,    21,
       0,    22,     0,     0,     0,     0,    23,     0,     0,     0,
       0,    24,    25,    26,     0,     0,    27,     0,     0,    28,
       0,    29,    30,     0,    31,     0,     0,     0,    32,    33,
      34,     0,    35,    36,    37,     0,     0,    38,     0,    39,
       0,    40,    41,    42,     0,    43,    44,     0,    45,    46,
       0,    47,     0,    48,     0,    49,     0,    50,     0,    51,
       0,    52,     0,    53,  -181,    54,     0,    55,    56,     0,
      57,    58,    59,     0,    60,    61,     0,    62,     0,     0,
      63,    64,    65,     0,    66,     1,    67,     2,     3,     4,
       5,     6,     0,     7,     8,     9,    10,    11,     0,     0,
       0,    12,     0,    13,    14,     0,     0,     0,     0,    15,
      16,     0,     0,    17,     0,     0,     0,     0,     0,     0,
      18,    19,     0,    20,    21,     0,    22,     0,     0,     0,
       0,    23,     0,     0,     0,     0,    24,    25,    26,     0,
       0,    27,     0,     0,    28,     0,    29,    30,     0,    31,
       0,     0,     0,    32,    33,    34,     0,    35,    36,    37,
       0,     0,    38,     0,    39,     0,    40,    41,    42,     0,
      43,    44,     0,    45,    46,     0,    47,     0,    48,     0,
      49,     0,    50,     0,    51,     0,    52,     0,    53,     0,
      54,  -181,    55,    56,     0,    57,    58,    59,     0,    60,
      61,     0,    62,     0,     0,    63,    64,    65,     0,    66,
       1,    67,     2,     3,     4,     5,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,     0,    13,    14,
       0,     0,     0,     0,    15,    16,     0,     0,    17,     0,
       0,     0,     0,     0,     0,    18,    19,     0,    20,    21,
       0,    22,     0,     0,     0,     0,    23,     0,     0,     0,
       0,    24,    25,    26,     0,     0,    27,     0,     0,    28,
       0,    29,    30,     0,    31,     0,     0,     0,    32,    33,
      34,     0,    35,    36,    37,   108,     0,    38,     0,    39,
       0,    40,    41,    42,     0,    43,    44,     0,    45,    46,
       0,    47,     0,    48,     0,    49,     0,    50,     0,    51,
       0,    52,     0,    53,     0,    54,     0,    55,    56,     0,
      57,    58,    59,     0,    60,    61,     0,    62,     0,     0,
      63,    64,    65,     0,    66,     1,    67,     2,     3,     4,
       5,     6,     0,     7,     8,     9,    10,    11,     0,     0,
       0,    12,   221,    13,    14,     0,     0,     0,     0,    15,
      16,     0,     0,    17,     0,     0,     0,     0,     0,     0,
      18,    19,     0,    20,    21,     0,    22,     0,     0,     0,
       0,    23,     0,     0,     0,     0,    24,    25,    26,     0,
       0,    27,     0,     0,    28,     0,    29,    30,     0,    31,
       0,     0,     0,    32,    33,    34,     0,    35,    36,    37,
       0,     0,    38,     0,    39,     0,    40,    41,    42,     0,
      43,    44,     0,    45,    46,     0,    47,     0,    48,     0,
      49,     0,    50,     0,    51,     0,    52,     0,    53,     0,
      54,     0,    55,    56,     0,    57,    58,    59,     0,    60,
      61,     0,    62,     0,     0,    63,    64,    65,     0,    66,
       1,    67,     2,     3,     4,     5,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,     0,    13,    14,
       0,     0,     0,     0,    15,    16,     0,     0,    17,     0,
       0,     0,     0,     0,     0,    18,    19,     0,    20,    21,
       0,    22,     0,  -179,     0,     0,    23,     0,     0,     0,
       0,    24,    25,    26,     0,     0,    27,     0,     0,    28,
       0,    29,    30,     0,    31,     0,     0,     0,    32,    33,
      34,     0,    35,    36,    37,     0,     0,    38,     0,    39,
       0,    40,    41,    42,     0,    43,    44,     0,    45,    46,
       0,    47,     0,    48,     0,    49,     0,    50,     0,    51,
       0,    52,     0,    53,     0,    54,     0,    55,    56,     0,
      57,    58,    59,     0,    60,    61,     0,    62,     0,     0,
      63,    64,    65,     0,    66,     1,    67,     2,     3,     4,
       5,     6,     0,     7,     8,     9,    10,    11,     0,     0,
       0,    12,  -181,    13,    14,     0,     0,     0,     0,    15,
      16,     0,     0,    17,     0,     0,     0,     0,     0,     0,
      18,    19,     0,    20,    21,     0,    22,     0,     0,     0,
       0,    23,     0,     0,     0,     0,    24,    25,    26,     0,
       0,    27,     0,     0,    28,     0,    29,    30,     0,    31,
       0,     0,     0,    32,    33,    34,     0,    35,    36,    37,
       0,     0,    38,     0,    39,     0,    40,    41,    42,     0,
      43,    44,     0,    45,    46,     0,    47,     0,    48,     0,
      49,     0,    50,     0,    51,     0,    52,     0,    53,     0,
      54,     0,    55,    56,     0,    57,    58,    59,     0,    60,
      61,     0,    62,     0,     0,    63,    64,    65,     0,    66,
       1,    67,     2,     3,     4,     5,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,   305,    13,    14,
       0,     0,     0,     0,    15,    16,     0,     0,    17,     0,
       0,     0,     0,     0,     0,    18,    19,     0,    20,    21,
       0,    22,     0,     0,     0,     0,    23,     0,     0,     0,
       0,    24,    25,    26,     0,     0,    27,     0,     0,    28,
       0,    29,    30,     0,    31,     0,     0,     0,    32,    33,
      34,     0,    35,    36,    37,     0,     0,    38,     0,    39,
       0,    40,    41,    42,     0,    43,    44,     0,    45,    46,
       0,    47,     0,    48,     0,    49,     0,    50,     0,    51,
       0,    52,     0,    53,     0,    54,     0,    55,    56,     0,
      57,    58,    59,     0,    60,    61,     0,    62,     0,     0,
      63,    64,    65,     0,    66,     1,    67,     2,     3,     4,
       5,     6,   329,     7,     8,     9,    10,    11,     0,     0,
       0,    12,     0,    13,    14,     0,     0,     0,     0,    15,
      16,     0,     0,    17,     0,     0,     0,     0,     0,     0,
      18,    19,     0,    20,    21,     0,    22,     0,     0,     0,
       0,    23,     0,     0,     0,     0,    24,    25,    26,     0,
       0,    27,     0,     0,    28,     0,    29,    30,     0,    31,
       0,     0,     0,    32,    33,    34,     0,    35,    36,    37,
       0,     0,    38,     0,    39,     0,    40,    41,    42,     0,
      43,    44,     0,    45,    46,     0,    47,     0,    48,     0,
      49,     0,    50,     0,    51,     0,    52,     0,    53,     0,
      54,     0,    55,    56,     0,    57,    58,    59,     0,    60,
      61,     0,    62,     0,     0,    63,    64,    65,     0,    66,
       1,    67,     2,     3,     4,     5,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,     0,    13,    14,
       0,     0,     0,     0,    15,    16,     0,     0,    17,     0,
       0,     0,     0,     0,     0,    18,    19,     0,    20,    21,
       0,    22,     0,     0,     0,     0,    23,     0,     0,     0,
       0,    24,    25,    26,     0,     0,    27,     0,     0,    28,
       0,    29,    30,     0,    31,     0,     0,     0,    32,    33,
      34,     0,    35,    36,    37,     0,     0,    38,   363,    39,
       0,    40,    41,    42,     0,    43,    44,     0,    45,    46,
       0,    47,     0,    48,     0,    49,     0,    50,     0,    51,
       0,    52,     0,    53,     0,    54,     0,    55,    56,     0,
      57,    58,    59,     0,    60,    61,     0,    62,     0,     0,
      63,    64,    65,     0,    66,     1,    67,     2,     3,     4,
       5,     6,     0,     7,     8,     9,    10,    11,     0,     0,
       0,    12,     0,    13,    14,     0,     0,     0,     0,    15,
      16,     0,     0,    17,     0,     0,     0,     0,     0,     0,
      18,    19,     0,    20,    21,     0,    22,     0,     0,     0,
       0,    23,     0,     0,     0,     0,    24,    25,    26,     0,
       0,    27,     0,     0,    28,     0,    29,    30,     0,    31,
       0,     0,     0,    32,    33,    34,     0,    35,    36,    37,
       0,   399,    38,     0,    39,     0,    40,    41,    42,     0,
      43,    44,     0,    45,    46,     0,    47,     0,    48,     0,
      49,     0,    50,     0,    51,     0,    52,     0,    53,     0,
      54,     0,    55,    56,     0,    57,    58,    59,     0,    60,
      61,     0,    62,     0,     0,    63,    64,    65,     0,    66,
       1,    67,     2,     3,     4,     5,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,     0,    13,    14,
       0,     0,     0,     0,    15,    16,     0,     0,    17,     0,
       0,     0,     0,     0,     0,    18,    19,     0,    20,    21,
       0,    22,     0,     0,     0,     0,    23,     0,     0,     0,
       0,    24,    25,    26,     0,     0,    27,     0,     0,    28,
       0,    29,    30,     0,    31,     0,     0,     0,    32,    33,
      34,     0,    35,    36,    37,     0,   294,    38,     0,    39,
       0,    40,    41,    42,     0,    43,    44,     0,    45,    46,
       0,    47,     0,    48,     0,    49,     0,    50,     0,    51,
       0,    52,     0,    53,     0,    54,     0,    55,   295,     0,
      57,    58,    59,     0,    60,    61,     0,    62,     0,     0,
      63,    64,    65,     0,    66,     1,    67,     2,     3,     4,
       5,     6,     0,     7,     8,     9,    10,    11,     0,     0,
       0,    12,     0,    13,    14,     0,     0,     0,     0,    15,
      16,     0,     0,    17,     0,     0,     0,     0,     0,     0,
      18,    19,     0,    20,    21,     0,    22,     0,     0,     0,
       0,    23,     0,     0,     0,   404,    24,    25,    26,     0,
       0,    27,     0,     0,    28,     0,    29,    30,     0,    31,
       0,     0,     0,    32,    33,    34,     0,    35,    36,    37,
       0,     0,    38,     0,    39,     0,    40,    41,    42,     0,
      43,    44,     0,    45,    46,     0,    47,     0,    48,     0,
      49,     0,    50,     0,    51,     0,    52,     0,    53,     0,
      54,     0,    55,    56,     0,    57,    58,    59,     0,    60,
      61,     0,    62,     0,     0,    63,    64,    65,     0,    66,
       1,    67,     2,     3,     4,     5,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,  -179,    13,    14,
       0,     0,     0,     0,    15,    16,     0,     0,    17,     0,
       0,     0,     0,     0,     0,    18,    19,     0,    20,    21,
       0,    22,     0,     0,     0,     0,    23,     0,     0,     0,
       0,    24,    25,    26,     0,     0,    27,     0,     0,    28,
       0,    29,    30,     0,    31,     0,     0,     0,    32,    33,
      34,     0,    35,    36,    37,     0,     0,    38,     0,    39,
       0,    40,    41,    42,     0,    43,    44,     0,    45,    46,
       0,    47,     0,    48,     0,    49,     0,    50,     0,    51,
       0,    52,     0,    53,     0,    54,     0,    55,    56,     0,
      57,    58,    59,     0,    60,    61,     0,    62,     0,     0,
      63,    64,    65,     0,    66,     1,    67,     2,     3,     4,
       5,     6,     0,     7,     8,     9,    10,    11,     0,     0,
       0,    12,     0,    13,    14,     0,     0,     0,     0,    15,
      16,     0,     0,    17,     0,     0,     0,     0,     0,     0,
      18,    19,     0,    20,    21,     0,    22,     0,     0,     0,
       0,    23,     0,     0,     0,     0,    24,    25,    26,     0,
       0,    27,     0,     0,    28,     0,    29,    30,     0,    31,
       0,     0,     0,    32,    33,    34,     0,    35,    36,    37,
       0,     0,    38,     0,    39,     0,    40,    41,    42,     0,
      43,    44,     0,    45,    46,     0,    47,     0,    48,     0,
      49,     0,    50,     0,    51,     0,    52,     0,    53,     0,
      54,     0,    55,    56,     0,    57,    58,    59,     0,    60,
      61,     0,    62,   543,     0,    63,    64,    65,     0,    66,
       1,    67,     2,     3,     4,     5,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,     0,    13,    14,
       0,     0,     0,     0,    15,    16,     0,     0,    17,     0,
       0,     0,     0,     0,     0,    18,    19,     0,    20,    21,
       0,    22,     0,     0,     0,     0,    23,     0,     0,     0,
       0,    24,    25,    26,     0,     0,    27,     0,     0,    28,
       0,    29,    30,     0,    31,     0,     0,     0,    32,    33,
      34,     0,    35,    36,    37,     0,     0,    38,     0,    39,
       0,    40,    41,    42,     0,    43,    44,     0,    45,    46,
       0,    47,     0,    48,     0,    49,     0,    50,     0,    51,
       0,    52,     0,    53,     0,    54,     0,    55,    56,   556,
      57,    58,    59,     0,    60,    61,     0,    62,     0,     0,
      63,    64,    65,     0,    66,     1,    67,     2,     3,     4,
       5,     6,     0,     7,     8,     9,    10,    11,     0,     0,
       0,    12,     0,    13,    14,     0,     0,     0,     0,    15,
      16,     0,     0,    17,     0,     0,     0,     0,     0,     0,
      18,    19,     0,    20,    21,     0,    22,     0,     0,     0,
       0,    23,     0,     0,     0,     0,    24,    25,    26,     0,
       0,    27,     0,     0,    28,     0,    29,    30,     0,    31,
       0,     0,     0,    32,    33,    34,     0,    35,    36,    37,
       0,     0,    38,     0,    39,     0,    40,    41,    42,     0,
      43,    44,     0,    45,    46,     0,    47,     0,    48,     0,
      49,     0,    50,     0,    51,     0,    52,     0,    53,     0,
      54,     0,    55,    56,   573,    57,    58,    59,     0,    60,
      61,     0,    62,     0,     0,    63,    64,    65,     0,    66,
       1,    67,     2,     3,     4,     5,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,     0,    13,    14,
       0,     0,     0,     0,    15,    16,     0,     0,    17,     0,
       0,     0,     0,     0,     0,    18,    19,     0,    20,    21,
       0,    22,     0,     0,     0,     0,    23,     0,     0,     0,
       0,    24,    25,    26,     0,     0,    27,     0,     0,    28,
       0,    29,    30,     0,    31,     0,     0,     0,    32,    33,
      34,     0,    35,    36,    37,     0,     0,    38,     0,    39,
       0,    40,    41,    42,     0,    43,    44,     0,    45,    46,
       0,    47,     0,    48,     0,    49,     0,    50,     0,    51,
       0,    52,     0,    53,     0,    54,     0,    55,    56,     0,
      57,    58,    59,     0,    60,    61,     0,    62,     0,     0,
      63,    64,    65,     0,    66,     1,    67,     2,     3,     4,
       5,     6,     0,     7,     8,     9,    10,    84,     0,     0,
       0,    12,     0,    13,    14,     0,     0,     0,     0,    15,
      16,     0,     0,    17,     0,     0,     0,     0,     0,     0,
      18,    19,     0,    20,    21,     0,    22,     0,     0,     0,
       0,    23,     0,     0,     0,     0,    24,    25,    26,     0,
       0,    27,     0,     0,    28,     0,    29,    30,     0,    31,
       0,     0,     0,    32,    33,    34,     0,    35,    36,    37,
       0,     0,    85,     0,    39,     0,    40,    41,    42,     0,
      43,    44,     0,    45,    46,     0,    47,     0,    48,     0,
      49,     0,    50,     0,    51,     0,    52,     0,    53,     0,
      54,     0,    55,    56,     0,    57,    58,    59,     0,    60,
      61,     0,    62,     0,     0,    63,    64,    65,     0,    66,
       1,    67,     2,     3,     4,     5,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    94,     0,    13,    14,
       0,     0,     0,     0,    15,    16,     0,     0,    17,     0,
       0,     0,     0,     0,     0,    18,    19,     0,    20,    21,
       0,    22,     0,     0,     0,     0,    23,     0,     0,     0,
       0,    24,    25,    26,     0,     0,    27,     0,     0,    28,
       0,    29,    30,     0,    31,     0,     0,     0,    32,    33,
      34,     0,    35,    36,    37,     0,     0,    38,     0,    39,
       0,    40,    41,    42,     0,    43,    44,     0,    45,    46,
       0,    47,     0,    48,     0,    49,     0,    50,     0,    51,
       0,    52,     0,    53,     0,    54,     0,    55,    56,     0,
      57,    58,    59,     0,    60,    61,     0,    62,     0,     0,
      63,    64,    65,     0,    66,     1,    67,     2,     3,     4,
       5,     6,     0,     7,     8,     9,    10,    11,     0,     0,
       0,    98,     0,    13,    14,     0,     0,     0,     0,    15,
      16,     0,     0,    17,     0,     0,     0,     0,     0,     0,
      18,    19,     0,    20,    21,     0,    22,     0,     0,     0,
       0,    23,     0,     0,     0,     0,    24,    25,    26,     0,
       0,    27,     0,     0,    28,     0,    29,    30,     0,    31,
       0,     0,     0,    32,    33,    34,     0,    35,    36,    37,
       0,     0,    38,     0,    39,     0,    40,    41,    42,     0,
      43,    44,     0,    45,    46,     0,    47,     0,    48,     0,
      49,     0,    50,     0,    51,     0,    52,     0,    53,     0,
      54,     0,    55,    56,     0,    57,    58,    59,     0,    60,
      61,     0,    62,     0,     0,    63,    64,    65,     0,    66,
       1,    67,     2,   236,     4,     5,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,     0,    13,    14,
       0,     0,     0,     0,    15,    16,     0,     0,    17,     0,
       0,     0,     0,     0,     0,    18,    19,     0,    20,    21,
       0,    22,     0,     0,     0,     0,    23,     0,     0,     0,
       0,    24,    25,    26,     0,     0,    27,     0,     0,    28,
       0,    29,    30,     0,    31,     0,     0,     0,    32,    33,
      34,     0,    35,    36,    37,     0,     0,    38,     0,    39,
       0,    40,    41,    42,     0,    43,    44,     0,    45,    46,
       0,    47,     0,    48,     0,    49,     0,    50,     0,    51,
       0,    52,     0,    53,     0,    54,     0,    55,    56,     0,
      57,    58,    59,     0,    60,    61,     0,    62,     0,     0,
      63,    64,    65,     0,    66,     1,    67,     2,     3,     4,
       5,     6,     0,     7,   267,     9,    10,    11,     0,     0,
       0,    12,     0,    13,    14,     0,     0,     0,     0,    15,
      16,     0,     0,    17,     0,     0,     0,     0,     0,     0,
      18,    19,     0,    20,    21,     0,    22,     0,     0,     0,
       0,    23,     0,     0,     0,     0,    24,    25,    26,     0,
       0,    27,     0,     0,    28,     0,    29,    30,     0,    31,
       0,     0,     0,    32,    33,    34,     0,    35,    36,    37,
       0,     0,    38,     0,    39,     0,    40,    41,    42,     0,
     268,    44,     0,    45,    46,     0,    47,     0,    48,     0,
      49,     0,    50,     0,    51,     0,    52,     0,    53,     0,
      54,     0,    55,    56,     0,    57,    58,    59,     0,    60,
      61,     0,    62,     0,     0,    63,    64,    65,     0,    66,
       1,    67,     2,     3,     4,     5,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,     0,    13,    14,
       0,     0,     0,     0,    15,    16,     0,     0,    17,     0,
       0,     0,     0,     0,     0,    18,    19,     0,    20,    21,
       0,    22,     0,     0,     0,     0,    23,     0,     0,     0,
       0,    24,    25,    26,     0,     0,    27,     0,     0,    28,
       0,    29,    30,     0,    31,     0,     0,     0,    32,    33,
      34,     0,    35,    36,    37,     0,     0,    38,     0,    39,
       0,    40,    41,    42,     0,    43,    44,     0,    45,    46,
       0,    47,     0,    48,     0,    49,     0,    50,     0,    51,
       0,    52,     0,    53,     0,    54,     0,    55,   541,     0,
      57,    58,    59,     0,    60,    61,     0,    62,     0,     0,
      63,    64,    65,     0,    66,     1,    67,     2,     3,     4,
       5,     6,     0,     7,     8,     9,    10,    11,     0,     0,
       0,    12,     0,    13,    14,     0,     0,     0,     0,    15,
      16,     0,     0,    17,     0,     0,     0,     0,     0,     0,
      18,    19,     0,    20,    21,     0,    22,     0,     0,     0,
       0,    23,     0,     0,     0,     0,    24,    25,    26,     0,
       0,    27,     0,     0,    28,     0,    29,    30,     0,    31,
       0,     0,     0,    32,    33,    34,     0,    35,    36,    37,
       0,     0,    38,     0,    39,     0,    40,    41,    42,     0,
      43,    44,     0,    45,    46,     0,    47,     0,    48,     0,
      49,     0,    50,     0,    51,     0,    52,     0,    53,     0,
      54,     0,    55,   560,     0,    57,    58,    59,     0,    60,
      61,     0,    62,     0,     0,    63,    64,    65,     0,    66,
       0,    67,   119,     3,     0,     0,   120,     0,     0,   121,
     122,     0,   100,   280,     0,   123,     0,     0,   124,   125,
     126,   127,   128,   129,   130,   131,     0,   132,   133,   134,
     135,   136,   137,   138,   139,   140,     0,     0,     0,   141,
       0,     0,     0,     0,     0,     0,     0,   142,     0,     0,
       0,   143,     0,   144,     0,     0,   145,     0,     0,   146,
       0,     0,     0,   147,     0,     0,     0,   148,     0,     0,
       0,     0,     0,     0,   119,     3,     0,     0,   120,     0,
       0,   121,   122,     0,   100,   149,   150,   123,    45,     0,
     124,   125,   126,   127,   128,   129,   130,   131,     0,   132,
     133,   134,   135,   136,   137,   138,   139,   140,     0,     0,
       0,   141,     0,     0,     0,     0,     0,     0,     0,   142,
      63,     0,     0,   143,   151,   144,     0,     0,   145,     0,
       0,   146,     0,     0,     0,   147,     0,     0,     0,   148,
       0,     0,     0,     0,     0,     0,   119,     3,   386,     0,
     120,     0,     0,   121,   122,     0,   100,   149,   150,   123,
      45,     0,   124,   125,   126,   127,   128,   129,   130,   131,
       0,   132,   133,   134,   135,   136,   137,   138,   139,   140,
       0,     0,     0,   141,     0,     0,     0,     0,     0,     0,
       0,   142,    63,     0,     0,   143,   151,   144,     0,     0,
     145,     0,     0,   146,     0,     0,     0,   147,     0,     0,
       0,   148,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   149,
     150,     0,    45,   184,     0,     0,     0,     0,     0,   185,
     186,   187,   188,     0,   189,   190,   191,   192,   193,   194,
       0,     0,   195,   196,   197,   198,   199,   200,   201,   202,
     203,   204,   205,   206,    63,     0,     0,     0,   151,   207,
     208,     0,   237,     0,     0,   185,   186,   187,   188,     0,
     189,   190,   191,   192,   193,   194,     0,     0,   195,   196,
     197,   198,   199,   200,   201,   202,   203,   204,   205,   206,
       0,     0,     0,     0,     0,     0,   208,     0,     0,     0,
     209,     0,   185,   186,   187,   188,   238,   189,   190,   191,
     192,   193,   194,     0,     0,   195,   196,   197,   198,   199,
     200,   201,   202,   203,   204,   205,   206,     0,     0,   210,
       0,     0,   211,   208,     0,   212,   209,     0,   232,   213,
     214,     0,     0,   215,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   210,     0,     0,   211,     0,
       0,   212,     0,   209,     0,   213,   214,     0,     0,   215,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   233,     0,     0,
       0,     0,   210,     0,     0,   211,     0,     0,   212,     0,
       0,     0,   213,   214,     0,     0,   215,   185,   186,   187,
     188,   224,   189,   190,   191,   192,   193,   194,     0,     0,
     195,   196,   197,   198,   199,   200,   201,   202,   203,   204,
     205,   206,     0,     0,   185,   186,   187,   188,   208,   189,
     190,   191,   192,   193,   194,     0,   227,   195,   196,   197,
     198,   199,   200,   201,   202,   203,   204,   205,   206,     0,
       0,     0,     0,     0,     0,   208,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   209,     0,
     185,   186,   187,   188,     0,   189,   190,   191,   192,   193,
     194,     0,     0,   195,   196,   197,   198,   199,   200,   201,
     202,   203,   204,   205,   206,   209,     0,   210,     0,     0,
     211,   208,     0,   212,     0,     0,     0,   213,   214,     0,
       0,   215,     0,     0,     0,     0,   250,     0,     0,     0,
       0,     0,     0,     0,   210,     0,     0,   211,     0,     0,
     212,     0,     0,     0,   213,   214,     0,     0,   215,     0,
       0,   209,     0,   287,   186,   187,   188,     0,   189,   190,
     191,   192,   193,   194,     0,     0,   195,   196,   197,   198,
     199,   200,   201,   202,   203,   204,   205,   206,     0,     0,
     210,     0,     0,   211,   208,     0,   212,     0,     0,     0,
     213,   214,     0,     0,   215,   185,   186,   187,   188,     0,
     189,   190,   191,   192,   193,   194,     0,     0,   195,   196,
     197,   198,   199,   200,   201,   202,   203,   204,   205,   206,
       0,     0,     0,     0,   209,     0,   208,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   288,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   210,     0,     0,   211,     0,     0,   212,
       0,     0,     0,   213,   214,     0,   209,   215,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     298,     0,     0,     0,     0,   210,     0,     0,   211,     0,
       0,   212,     0,     0,     0,   213,   214,     0,     0,   215,
     185,   186,   187,   188,   339,   189,   190,   191,   192,   193,
     194,     0,     0,   195,   196,   197,   198,   199,   200,   201,
     202,   203,   204,   205,   206,     0,     0,   185,   186,   187,
     188,   208,   189,   190,   191,   192,   193,   194,     0,     0,
     195,   196,   197,   198,   199,   200,   201,   202,   203,   204,
     205,   206,     0,     0,   340,     0,     0,     0,   208,     0,
     185,   186,   187,   188,   342,   189,   190,   191,   192,   193,
     194,   209,     0,   195,   196,   197,   198,   199,   200,   201,
     202,   203,   204,   205,   206,     0,     0,     0,     0,     0,
       0,   208,     0,     0,     0,     0,     0,     0,   209,     0,
     210,     0,     0,   211,     0,     0,   212,     0,     0,     0,
     213,   214,     0,     0,   215,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   210,     0,     0,
     211,   209,     0,   212,     0,     0,     0,   213,   214,     0,
       0,   215,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     210,     0,     0,   211,     0,     0,   212,     0,     0,     0,
     213,   214,     0,     0,   215,   185,   186,   187,   188,   346,
     189,   190,   191,   192,   193,   194,     0,     0,   195,   196,
     197,   198,   199,   200,   201,   202,   203,   204,   205,   206,
       0,     0,     0,     0,     0,     0,   208,     0,   185,   186,
     187,   188,   347,   189,   190,   191,   192,   193,   194,     0,
       0,   195,   196,   197,   198,   199,   200,   201,   202,   203,
     204,   205,   206,     0,     0,     0,     0,     0,     0,   208,
       0,     0,     0,     0,     0,     0,   209,     0,     0,   185,
     186,   187,   188,   352,   189,   190,   191,   192,   193,   194,
       0,     0,   195,   196,   197,   198,   199,   200,   201,   202,
     203,   204,   205,   206,     0,   210,     0,     0,   211,   209,
     208,   212,     0,     0,     0,   213,   214,     0,     0,   215,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   210,     0,
       0,   211,     0,     0,   212,     0,     0,     0,   213,   214,
     209,     0,   215,   185,   186,   187,   188,   359,   189,   190,
     191,   192,   193,   194,     0,     0,   195,   196,   197,   198,
     199,   200,   201,   202,   203,   204,   205,   206,     0,   210,
       0,     0,   211,     0,   208,   212,     0,     0,     0,   213,
     214,     0,     0,   215,   185,   186,   187,   188,   369,   189,
     190,   191,   192,   193,   194,     0,     0,   195,   196,   197,
     198,   199,   200,   201,   202,   203,   204,   205,   206,     0,
       0,     0,     0,     0,   209,   208,     0,   185,   186,   187,
     188,   371,   189,   190,   191,   192,   193,   194,     0,     0,
     195,   196,   197,   198,   199,   200,   201,   202,   203,   204,
     205,   206,     0,   210,     0,     0,   211,     0,   208,   212,
       0,     0,     0,   213,   214,   209,     0,   215,   185,   186,
     187,   188,   372,   189,   190,   191,   192,   193,   194,     0,
       0,   195,   196,   197,   198,   199,   200,   201,   202,   203,
     204,   205,   206,     0,   210,     0,     0,   211,   209,   208,
     212,     0,     0,     0,   213,   214,     0,     0,   215,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   210,     0,     0,
     211,     0,     0,   212,     0,     0,     0,   213,   214,   209,
       0,   215,   185,   186,   187,   188,   373,   189,   190,   191,
     192,   193,   194,     0,     0,   195,   196,   197,   198,   199,
     200,   201,   202,   203,   204,   205,   206,     0,   210,     0,
       0,   211,     0,   208,   212,     0,     0,     0,   213,   214,
       0,     0,   215,   185,   186,   187,   188,   374,   189,   190,
     191,   192,   193,   194,     0,     0,   195,   196,   197,   198,
     199,   200,   201,   202,   203,   204,   205,   206,     0,     0,
       0,     0,     0,   209,   208,     0,   185,   186,   187,   188,
       0,   189,   190,   191,   192,   193,   194,     0,     0,   195,
     196,   197,   198,   199,   200,   201,   202,   203,   204,   205,
     206,     0,   210,   375,     0,   211,     0,   208,   212,     0,
       0,     0,   213,   214,   209,     0,   215,   185,   186,   187,
     188,     0,   189,   190,   191,   192,   193,   194,     0,   383,
     195,   196,   197,   198,   199,   200,   201,   202,   203,   204,
     205,   206,     0,   210,     0,     0,   211,   209,   208,   212,
       0,     0,     0,   213,   214,     0,     0,   215,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   210,     0,     0,   211,
       0,     0,   212,     0,     0,     0,   213,   214,   209,     0,
     215,   185,   186,   187,   188,   395,   189,   190,   191,   192,
     193,   194,     0,     0,   195,   196,   197,   198,   199,   200,
     201,   202,   203,   204,   205,   206,     0,   210,     0,     0,
     211,     0,   208,   212,     0,     0,     0,   213,   214,     0,
       0,   215,   185,   186,   187,   188,     0,   189,   190,   191,
     192,   193,   194,     0,     0,   195,   196,   197,   198,   199,
     200,   201,   202,   203,   204,   205,   206,     0,     0,   406,
       0,     0,   209,   208,     0,   185,   186,   187,   188,     0,
     189,   190,   191,   192,   193,   194,     0,     0,   195,   196,
     197,   198,   199,   200,   201,   202,   203,   204,   205,   206,
       0,   210,   409,     0,   211,     0,   208,   212,     0,     0,
       0,   213,   214,   209,     0,   215,   185,   186,   187,   188,
       0,   189,   190,   191,   192,   193,   194,     0,     0,   195,
     196,   197,   198,   199,   200,   201,   202,   203,   204,   205,
     206,     0,   210,     0,     0,   211,   209,   208,   212,     0,
       0,     0,   213,   214,     0,     0,   215,     0,     0,     0,
       0,   188,   428,   189,   190,   191,   192,   193,   194,     0,
       0,   195,     0,     0,   198,   210,     0,     0,   211,   203,
     204,   212,   206,     0,     0,   213,   214,   209,     0,   215,
     185,   186,   187,   188,     0,   189,   190,   191,   192,   193,
     194,     0,     0,   195,   196,   197,   198,   199,   200,   201,
     202,   203,   204,   205,   206,     0,   210,     0,     0,   211,
     437,   208,   212,     0,     0,     0,   213,   214,     0,   209,
     215,     0,   185,   186,   187,   188,     0,   189,   190,   191,
     192,   193,   194,     0,     0,   195,   196,   197,   198,   199,
     200,   201,   202,   203,   204,   205,   206,     0,     0,     0,
       0,   209,   473,   208,     0,     0,     0,     0,     0,   214,
       0,     0,   215,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     210,     0,     0,   211,     0,     0,   212,     0,     0,     0,
     213,   214,     0,   209,   215,   185,   186,   187,   188,     0,
     189,   190,   191,   192,   193,   194,     0,     0,   195,   196,
     197,   198,   199,   200,   201,   202,   203,   204,   205,   206,
       0,     0,   210,     0,     0,   211,   208,     0,   212,     0,
       0,     0,   213,   214,     0,     0,   215,     0,     0,     0,
       0,     0,   357,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   185,   186,   187,   188,     0,   189,
     190,   191,   192,   193,   194,     0,   209,   195,   196,   197,
     198,   199,   200,   201,   202,   203,   204,   205,   206,     0,
       0,     0,     0,     0,   492,   208,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   210,     0,     0,   211,     0,
       0,   212,     0,     0,     0,   213,   214,     0,     0,   215,
       0,     0,     0,   185,   186,   187,   188,   498,   189,   190,
     191,   192,   193,   194,     0,   209,   195,   196,   197,   198,
     199,   200,   201,   202,   203,   204,   205,   206,     0,     0,
       0,     0,     0,     0,   208,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   210,     0,     0,   211,     0,     0,
     212,     0,     0,     0,   213,   214,     0,     0,   215,     0,
       0,     0,     0,   185,   186,   187,   188,     0,   189,   190,
     191,   192,   193,   194,   209,   501,   195,   196,   197,   198,
     199,   200,   201,   202,   203,   204,   205,   206,     0,     0,
       0,     0,     0,     0,   208,     0,     0,     0,     0,     0,
       0,     0,     0,   210,     0,     0,   211,     0,     0,   212,
       0,     0,     0,   213,   214,     0,     0,   215,     0,     0,
       0,     0,   185,   186,   187,   188,     0,   189,   190,   191,
     192,   193,   194,     0,   209,   195,   196,   197,   198,   199,
     200,   201,   202,   203,   204,   205,   206,     0,     0,   512,
       0,     0,     0,   208,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   210,     0,     0,   211,     0,     0,   212,
       0,     0,     0,   213,   214,     0,     0,   215,     0,     0,
       0,   185,   186,   187,   188,     0,   189,   190,   191,   192,
     193,   194,     0,   209,   195,   196,   197,   198,   199,   200,
     201,   202,   203,   204,   205,   206,     0,     0,     0,     0,
       0,     0,   208,     0,     0,     0,     0,   518,     0,     0,
       0,     0,   210,     0,     0,   211,     0,     0,   212,     0,
       0,     0,   213,   214,     0,     0,   215,     0,     0,     0,
     185,   186,   187,   188,     0,   189,   190,   191,   192,   193,
     194,     0,   209,   195,   196,   197,   198,   199,   200,   201,
     202,   203,   204,   205,   206,     0,     0,     0,     0,     0,
       0,   208,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   210,     0,     0,   211,     0,   529,   212,     0,     0,
       0,   213,   214,     0,     0,   215,     0,     0,     0,   185,
     186,   187,   188,     0,   189,   190,   191,   192,   193,   194,
       0,   209,   195,   196,   197,   198,   199,   200,   201,   202,
     203,   204,   205,   206,     0,     0,     0,     0,     0,     0,
     208,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     210,     0,     0,   211,     0,   531,   212,     0,     0,     0,
     213,   214,     0,     0,   215,     0,     0,     0,   185,   186,
     187,   188,     0,   189,   190,   191,   192,   193,   194,     0,
     209,   195,   196,   197,   198,   199,   200,   201,   202,   203,
     204,   205,   206,     0,     0,     0,     0,     0,   545,   208,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   210,
       0,     0,   211,     0,     0,   212,     0,     0,     0,   213,
     214,     0,     0,   215,     0,     0,     0,   185,   186,   187,
     188,     0,   189,   190,   191,   192,   193,   194,     0,   209,
     195,   196,   197,   198,   199,   200,   201,   202,   203,   204,
     205,   206,     0,     0,     0,     0,     0,     0,   208,     0,
       0,     0,     0,   547,     0,     0,     0,     0,   210,     0,
       0,   211,     0,     0,   212,     0,     0,     0,   213,   214,
       0,     0,   215,     0,     0,     0,   185,   186,   187,   188,
       0,   189,   190,   191,   192,   193,   194,     0,   209,   195,
     196,   197,   198,   199,   200,   201,   202,   203,   204,   205,
     206,     0,     0,     0,     0,     0,   577,   208,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   210,     0,     0,
     211,     0,     0,   212,     0,     0,     0,   213,   214,     0,
       0,   215,     0,     0,     0,   185,   186,   187,   188,     0,
     189,   190,   191,   192,   193,   194,     0,   209,   195,   196,
     197,   198,   199,   200,   201,   202,   203,   204,   205,   206,
       0,     0,     0,     0,     0,     0,   208,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   210,     0,     0,   211,
       0,     0,   212,     0,     0,     0,   213,   214,     0,     0,
     215,     0,     0,     0,   289,   186,   187,   188,     0,   189,
     190,   191,   192,   193,   194,     0,   209,   195,   196,   197,
     198,   199,   200,   201,   202,   203,   204,   205,   206,     0,
       0,     0,     0,     0,     0,   208,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   210,     0,     0,   211,     0,
       0,   212,     0,     0,     0,   213,   214,     0,     0,   215,
       0,     0,     0,   185,   186,   187,   188,     0,   189,   190,
     191,   192,   193,   194,     0,   209,   195,   196,   197,   198,
     199,   200,   201,   202,   203,   204,   205,   206,     0,     0,
       0,     0,     0,     0,   208,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   210,     0,     0,   211,     0,     0,
     212,     0,     0,     0,   213,   214,     0,     0,   215,     0,
       0,     0,   185,  -259,   187,   188,     0,   189,   190,   191,
     192,   193,   194,     0,   209,   195,   196,   197,   198,   199,
     200,   201,   202,   203,   204,   205,   206,     0,     0,     0,
       0,     0,     0,   208,     0,     0,     0,   188,     0,   189,
     190,   191,   192,   193,   194,     0,  -259,   195,     0,   212,
       0,     0,     0,   213,   214,   203,   204,   215,   206,     0,
       0,   185,     0,   187,   188,     0,   189,   190,   191,   192,
     193,   194,     0,   209,   195,   196,   197,   198,   199,   200,
     201,   202,   203,   204,   205,   206,     0,     0,     0,     0,
       0,     0,   208,     0,     0,     0,   188,     0,     0,     0,
     191,   192,   193,   194,     0,   209,   195,     0,   212,     0,
       0,     0,   213,   214,   203,   204,   215,   206,     0,     0,
     185,     0,   187,   188,     0,   189,   190,   191,   192,   193,
     194,     0,   209,   195,   196,   197,   198,   199,   200,   201,
     202,   203,   204,   205,   206,   214,     0,     0,   215,     0,
       0,   208,     0,     0,     0,   188,     0,     0,     0,   191,
     192,  -259,   194,     0,   209,   195,     0,     0,     0,     0,
       0,   213,   214,   203,   204,   215,   206,     0,     0,     0,
     186,   187,   188,     0,   189,   190,   191,   192,   193,   194,
       0,   209,   195,   196,   197,   198,   199,   200,   201,   202,
     203,   204,   205,   206,   214,     0,     0,   215,     0,     0,
     208,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   209,     0,     0,     0,     0,     0,     0,
    -259,   214,     0,     0,   215,     0,     0,     0,     0,     0,
     187,   188,     0,   189,   190,   191,   192,   193,   194,     0,
     209,   195,   196,   197,   198,   199,   200,   201,   202,   203,
     204,   205,   206,   214,     0,     0,   215,     0,     0,   208,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   210,
       0,     0,   211,     0,     0,   212,     0,     0,     0,   213,
     214,     0,     0,   215,     0,     0,     0,     0,     0,   187,
     188,     0,   189,   190,   191,   192,   193,   194,     0,   209,
     195,   196,   197,   198,   199,   200,   201,   202,   203,   204,
     205,   206,     0,   458,     0,     0,   187,   188,   208,   189,
     190,   191,   192,   193,   194,     0,     0,   195,   196,   197,
     198,   199,   200,   201,   202,   203,   204,   205,   206,   214,
       0,     0,   215,     0,     0,   208,     0,     0,  -259,   188,
       0,   189,   190,   191,   192,   193,   194,     0,   209,   195,
     196,   197,   198,   199,   200,   201,   202,   203,   204,   205,
     206,     0,     0,     0,   459,     0,   188,     0,   189,   190,
     191,   192,   193,   194,     0,   209,   195,   196,   197,   198,
     199,   200,   201,     0,   203,   204,     0,   206,   214,     0,
       0,   215,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   209,     0,     0,
       0,     0,     0,     0,     0,   214,     0,     0,   215,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   209,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   214,     0,     0,
     215,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   214,     0,     0,   215
};

static const yytype_int16 yycheck[] =
{
       0,     3,     4,   353,   254,    24,    17,    51,    10,     4,
      10,    13,    12,    13,    14,    17,    16,    17,    13,    19,
      20,    21,    22,    23,     4,   546,    26,    64,    15,   416,
      30,    26,     4,    13,    17,    50,    18,    17,    38,    15,
      44,    13,   563,    51,     4,     5,    46,    47,    48,    49,
      50,    51,    52,    53,    54,    55,    56,    46,    58,    59,
      64,    17,    77,    63,    46,   109,    66,    67,    24,    77,
      46,    27,    83,     4,    58,    75,    76,    77,    78,    35,
      36,    17,    38,    67,    67,    85,   473,    17,    54,    55,
      61,    27,    87,     4,    94,    26,    96,    97,    98,    70,
      17,   101,    38,    17,   104,   105,    15,   107,    38,   128,
     110,   111,   112,   113,   151,   115,   116,   117,   118,   121,
     370,    38,     4,    17,    38,    59,   121,    61,    17,    85,
      76,   131,   119,    20,    21,    22,   486,   141,   142,    61,
      61,   121,   146,   147,   148,    17,   150,   151,    70,   121,
      17,   155,    17,    17,    17,    17,   156,    17,    38,    46,
      47,    48,    49,    26,     0,    52,    53,    54,    17,   125,
     170,    58,   128,   173,   174,    38,    59,    60,   125,   179,
      63,   181,    41,    41,    72,   185,   186,   187,   188,   189,
     190,   191,   192,   193,   194,    43,   196,   197,   198,   199,
     200,   201,   202,   203,   204,   205,   206,   207,   208,    88,
      91,   211,   212,   213,   214,    31,   216,   217,    93,    95,
      97,   105,   107,    17,   111,   103,   226,   264,    17,    17,
      17,   122,    17,   233,    22,    24,    24,     3,    27,    27,
     109,    57,   242,   243,    60,    61,    35,    35,    36,    38,
      38,   251,   252,    56,    44,   255,    81,    17,    26,    18,
     264,    41,    66,    65,    17,    88,   187,    50,   272,   273,
      18,    24,    18,    51,    27,    76,   280,    18,    60,    79,
      18,    15,    35,    36,   170,    38,    18,   287,   174,   289,
      44,    46,   108,   293,    44,   295,    14,    85,    77,    29,
     300,   188,   302,    77,    88,    18,    51,     3,    45,    60,
      77,   232,    17,    77,    19,    20,    21,    22,    23,    24,
      18,    45,    27,    28,    29,    30,    31,    32,    33,   216,
      35,    36,    37,    38,    18,    77,   109,   125,   109,   128,
     128,   341,   342,   343,   109,   497,   109,   442,   429,   473,
     350,   207,   564,   353,   485,    -1,   356,   357,    -1,   359,
     360,    -1,    -1,    -1,    -1,   365,    -1,    -1,   368,    -1,
      -1,    -1,    -1,    -1,   378,   128,    -1,   377,    -1,    -1,
      85,   385,    -1,    -1,    -1,   389,   390,    -1,    -1,   141,
     142,    -1,    -1,    -1,   146,   147,   148,    -1,   150,    -1,
     400,   401,   402,   155,    -1,    -1,   293,    -1,   408,    -1,
     410,   411,    27,   299,   414,   302,   232,   417,    -1,   419,
     125,   342,   238,   128,    -1,    -1,    -1,    -1,    -1,    44,
      -1,   431,    -1,   433,   250,    -1,    17,   253,   359,    -1,
      -1,    56,    -1,    24,    -1,    -1,    27,    62,   452,    64,
      -1,    -1,    -1,    -1,    35,    36,    17,    38,    19,    20,
      21,    22,    23,    24,   464,    -1,    27,    -1,    -1,    30,
     470,    32,    33,    -1,    35,    36,   476,    38,   478,    -1,
     480,   297,   298,   299,    -1,   485,   486,   487,    -1,    -1,
     490,    -1,    -1,   493,   494,    -1,    -1,    -1,   113,   503,
     500,    -1,   423,    -1,    85,    -1,    -1,    -1,   123,    -1,
     510,    -1,    -1,    -1,   400,   401,   402,    -1,   404,    -1,
     272,   273,    -1,   410,    85,    -1,   141,   142,   280,    -1,
     145,   146,   147,   148,    -1,   150,   151,   537,    -1,    -1,
     155,   541,   542,    -1,   125,    -1,   546,   128,    -1,   470,
     550,   551,   552,   553,    -1,   476,   556,   557,   558,    -1,
     560,   561,    17,   563,   125,    -1,    -1,   128,    -1,    24,
      -1,   571,    27,   573,    -1,    -1,    -1,   498,   464,    -1,
      35,    36,   582,    38,   584,    -1,    -1,   587,   404,   589,
     205,    -1,    -1,    -1,   480,   210,    -1,    -1,    -1,    -1,
      -1,   522,    -1,    -1,    -1,    -1,   422,    -1,    -1,   425,
      -1,    -1,   428,    -1,    -1,    -1,   432,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   378,   548,    -1,   550,
      85,    -1,   448,   385,    -1,    -1,    -1,   389,   390,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   264,
      -1,   266,    -1,    -1,    -1,    -1,   542,   272,   273,    -1,
      -1,    -1,    -1,    -1,    -1,   280,    -1,    -1,    -1,   556,
      -1,   558,    -1,   128,    -1,   561,    -1,    -1,    -1,    -1,
     295,    -1,    -1,   499,   571,   300,   573,    -1,    -1,    -1,
      -1,    -1,   508,    -1,    -1,    -1,   582,    -1,   584,    -1,
     452,   587,   518,   589,    17,    -1,    19,    20,    21,    22,
      23,    24,   528,   529,    27,   531,    -1,    30,    31,    32,
      33,    -1,    35,    36,    -1,    38,    -1,   543,    -1,    -1,
      -1,   547,    -1,    -1,    17,    -1,    19,    20,    21,    22,
      23,    24,    -1,    -1,    27,    -1,    -1,    30,   363,    32,
      33,   503,    35,    36,    -1,    38,    -1,    -1,    -1,    -1,
      -1,    -1,   578,   378,   580,    -1,    -1,   583,    -1,   585,
     385,    -1,    85,    -1,   389,   390,    -1,    -1,    -1,    -1,
      -1,    -1,     1,    -1,     3,     4,     5,     6,     7,    -1,
       9,    10,    11,    12,    13,    -1,    -1,    -1,    17,    -1,
      19,    20,    85,    -1,   419,    -1,    25,    26,    -1,    -1,
      29,    -1,   125,    -1,    -1,   128,    -1,    36,    37,    -1,
      39,    40,    -1,    42,    -1,    44,    -1,   442,    47,    -1,
      49,    -1,    51,    52,    53,    54,    -1,   452,    57,    -1,
      -1,    60,   125,    62,    63,   128,    65,    -1,    -1,    -1,
      69,    70,    71,    -1,    73,    74,    75,    -1,    77,    78,
      -1,    80,    -1,    82,    83,    84,    -1,    86,    87,    -1,
      89,    90,    -1,    92,    -1,    94,    -1,    96,    -1,    98,
      -1,   100,   497,   102,    -1,   104,    -1,   106,   503,   108,
     109,    -1,   111,   112,   113,    -1,   115,   116,    -1,   118,
      -1,    -1,   121,   122,   123,    -1,   125,    -1,   127,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   541,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   557,    -1,     1,   560,     3,     4,     5,     6,
       7,     8,     9,    10,    11,    12,    13,    14,    15,    16,
      17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    43,    44,    45,    -1,
      47,    48,    49,    50,    51,    52,    53,    54,    55,    56,
      57,    58,    59,    60,    61,    62,    63,    64,    65,    66,
      67,    68,    69,    70,    71,    -1,    73,    74,    75,    -1,
      77,    78,    -1,    80,    -1,    82,    83,    84,    85,    86,
      87,    88,    89,    90,    91,    92,    93,    94,    95,    96,
      97,    98,    99,   100,   101,   102,   103,   104,   105,   106,
     107,   108,   109,   110,   111,   112,   113,   114,   115,   116,
     117,   118,    -1,   120,   121,   122,   123,   124,   125,    -1,
     127,   128,     1,    -1,     3,     4,     5,     6,     7,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,    45,    -1,    47,    48,
      49,    50,    51,    52,    53,    54,    55,    56,    57,    58,
      59,    60,    61,    62,    63,    64,    65,    66,    67,    68,
      69,    70,    71,    -1,    73,    74,    75,    -1,    77,    -1,
      -1,    80,    -1,    82,    83,    84,    85,    86,    87,    88,
      89,    90,    91,    92,    93,    94,    95,    96,    97,    98,
      99,   100,   101,   102,   103,   104,   105,   106,   107,   108,
     109,   110,   111,   112,   113,   114,   115,   116,   117,   118,
      -1,   120,   121,   122,   123,   124,   125,    -1,   127,   128,
       1,    -1,     3,     4,     5,     6,     7,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    39,    40,
      41,    42,    43,    44,    45,    -1,    47,    48,    49,    50,
      51,    52,    53,    54,    55,    56,    57,    58,    59,    60,
      61,    62,    63,    64,    65,    66,    67,    68,    69,    70,
      71,    -1,    73,    74,    75,    -1,    77,    78,    -1,    80,
      -1,    82,    83,    84,    85,    86,    87,    88,    89,    90,
      91,    92,    93,    94,    95,    96,    97,    98,    99,   100,
     101,   102,   103,   104,   105,   106,   107,   108,   109,   110,
     111,   112,   113,   114,   115,   116,   117,   118,    -1,   120,
     121,   122,   123,   124,   125,    -1,     1,   128,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    -1,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    -1,    73,    74,
      75,    -1,    77,    78,    -1,    80,    -1,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,    -1,   110,   111,   112,   113,   114,
     115,   116,   117,   118,    -1,   120,   121,   122,   123,   124,
     125,    -1,   127,   128,     1,    -1,     3,     4,     5,     6,
       7,    -1,     9,    10,    11,    12,    13,    14,    15,    16,
      17,    -1,    19,    20,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    -1,    42,    -1,    -1,    45,    -1,
      47,    -1,    -1,    50,    -1,    52,    53,    54,    -1,    -1,
      57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,    -1,
      -1,    -1,    69,    70,    71,    -1,    73,    74,    75,    76,
      -1,    78,    -1,    80,    -1,    82,    83,    84,    85,    86,
      87,    -1,    89,    90,    -1,    92,    -1,    94,    -1,    96,
      -1,    98,    -1,   100,    -1,   102,    -1,   104,    -1,   106,
      -1,   108,   109,    -1,   111,   112,   113,   114,   115,   116,
     117,   118,    -1,   120,   121,   122,   123,   124,   125,    -1,
     127,   128,     1,    -1,     3,     4,     5,     6,     7,    -1,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    -1,
      19,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    -1,    42,    -1,    -1,    45,    -1,    47,    -1,
      -1,    -1,    -1,    52,    53,    54,    -1,    -1,    57,    -1,
      -1,    60,    -1,    62,    63,    -1,    65,    -1,    -1,    -1,
      69,    70,    71,    -1,    73,    74,    75,    -1,    -1,    78,
      -1,    80,    -1,    82,    83,    84,    85,    86,    87,    -1,
      89,    90,    -1,    92,    -1,    94,    -1,    96,    -1,    98,
      -1,   100,    -1,   102,    -1,   104,    -1,   106,    -1,   108,
     109,    -1,   111,   112,   113,   114,   115,   116,   117,   118,
      -1,   120,   121,   122,   123,   124,   125,    -1,   127,   128,
       1,    -1,     3,     4,     5,     6,     7,    -1,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    17,    -1,    19,    20,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    36,    37,    -1,    39,    40,
      -1,    42,    -1,    -1,    -1,    -1,    47,    -1,    49,    -1,
      51,    52,    53,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    62,    63,    -1,    65,    -1,    -1,    -1,    69,    70,
      71,    -1,    73,    74,    75,    -1,    77,    78,    -1,    80,
      -1,    82,    83,    84,    -1,    86,    87,    -1,    89,    90,
      -1,    92,    -1,    94,    -1,    96,    -1,    98,    -1,   100,
      -1,   102,    -1,   104,    -1,   106,    -1,   108,   109,    -1,
     111,   112,   113,    -1,   115,   116,    -1,   118,    -1,    -1,
     121,   122,   123,    -1,   125,     1,   127,     3,     4,     5,
       6,     7,    -1,     9,    10,    11,    12,    13,    -1,    -1,
      -1,    17,    -1,    19,    20,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    -1,    39,    40,    -1,    42,    -1,    44,    -1,
      -1,    47,    -1,    -1,    -1,    -1,    52,    53,    54,    55,
      56,    57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,
      -1,    -1,    -1,    69,    70,    71,    -1,    73,    74,    75,
      -1,    -1,    78,    -1,    80,    -1,    82,    83,    84,    -1,
      86,    87,    -1,    89,    90,    -1,    92,    -1,    94,    -1,
      96,    -1,    98,    -1,   100,    -1,   102,    -1,   104,    -1,
     106,    -1,   108,   109,    -1,   111,   112,   113,    -1,   115,
     116,    -1,   118,    -1,    -1,   121,   122,   123,    -1,   125,
       1,   127,     3,     4,     5,     6,     7,    -1,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    17,    -1,    19,    20,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    36,    37,    -1,    39,    40,
      -1,    42,    -1,    44,    -1,    -1,    47,    -1,    -1,    -1,
      -1,    52,    53,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    62,    63,    -1,    65,    -1,    -1,    -1,    69,    70,
      71,    -1,    73,    74,    75,    -1,    77,    78,    -1,    80,
      -1,    82,    83,    84,    -1,    86,    87,    -1,    89,    90,
      -1,    92,    -1,    94,    -1,    96,    -1,    98,    -1,   100,
      -1,   102,    -1,   104,    -1,   106,    -1,   108,   109,   110,
     111,   112,   113,    -1,   115,   116,    -1,   118,    -1,    -1,
     121,   122,   123,    -1,   125,     1,   127,     3,     4,     5,
       6,     7,    -1,     9,    10,    11,    12,    13,    -1,    -1,
      -1,    17,    -1,    19,    20,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    -1,    39,    40,    -1,    42,    -1,    44,    -1,
      -1,    47,    -1,    -1,    -1,    -1,    52,    53,    54,    -1,
      -1,    57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,
      -1,    -1,    -1,    69,    70,    71,    -1,    73,    74,    75,
      -1,    77,    78,    -1,    80,    -1,    82,    83,    84,    -1,
      86,    87,    -1,    89,    90,    -1,    92,    -1,    94,    -1,
      96,    -1,    98,    -1,   100,    -1,   102,    -1,   104,    -1,
     106,    -1,   108,   109,    -1,   111,   112,   113,    -1,   115,
     116,    -1,   118,    -1,    -1,   121,   122,   123,    -1,   125,
       1,   127,     3,     4,     5,     6,     7,    -1,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    17,    -1,    19,    20,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    36,    37,    -1,    39,    40,
      -1,    42,    -1,    44,    -1,    -1,    47,    -1,    -1,    -1,
      -1,    52,    53,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    62,    63,    -1,    65,    66,    -1,    -1,    69,    70,
      71,    -1,    73,    74,    75,    -1,    -1,    78,    -1,    80,
      -1,    82,    83,    84,    -1,    86,    87,    -1,    89,    90,
      -1,    92,    -1,    94,    -1,    96,    -1,    98,    -1,   100,
      -1,   102,    -1,   104,    -1,   106,    -1,   108,   109,    -1,
     111,   112,   113,    -1,   115,   116,    -1,   118,    -1,    -1,
     121,   122,   123,    -1,   125,     1,   127,     3,     4,     5,
       6,     7,    -1,     9,    10,    11,    12,    13,    -1,    -1,
      -1,    17,    -1,    19,    20,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    -1,    39,    40,    -1,    42,    -1,    -1,    -1,
      -1,    47,    -1,    -1,    -1,    -1,    52,    53,    54,    -1,
      -1,    57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,
      -1,    -1,    -1,    69,    70,    71,    -1,    73,    74,    75,
      -1,    -1,    78,    79,    80,    81,    82,    83,    84,    -1,
      86,    87,    -1,    89,    90,    -1,    92,    -1,    94,    -1,
      96,    -1,    98,    -1,   100,    -1,   102,    -1,   104,    -1,
     106,    -1,   108,   109,    -1,   111,   112,   113,    -1,   115,
     116,    -1,   118,    -1,    -1,   121,   122,   123,    -1,   125,
       1,   127,     3,     4,     5,     6,     7,    -1,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    17,    -1,    19,    20,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    36,    37,    -1,    39,    40,
      -1,    42,    -1,    44,    -1,    -1,    47,    -1,    -1,    -1,
      -1,    52,    53,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    62,    63,    -1,    65,    -1,    -1,    -1,    69,    70,
      71,    -1,    73,    74,    75,    -1,    77,    78,    -1,    80,
      -1,    82,    83,    84,    -1,    86,    87,    -1,    89,    90,
      -1,    92,    -1,    94,    -1,    96,    -1,    98,    -1,   100,
      -1,   102,    -1,   104,    -1,   106,    -1,   108,   109,    -1,
     111,   112,   113,    -1,   115,   116,    -1,   118,    -1,    -1,
     121,   122,   123,    -1,   125,     1,   127,     3,     4,     5,
       6,     7,    -1,     9,    10,    11,    12,    13,    -1,    -1,
      -1,    17,    -1,    19,    20,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    -1,    39,    40,    -1,    42,    -1,    44,    -1,
      -1,    47,    -1,    -1,    -1,    51,    52,    53,    54,    -1,
      -1,    57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,
      -1,    -1,    -1,    69,    70,    71,    -1,    73,    74,    75,
      -1,    -1,    78,    -1,    80,    -1,    82,    83,    84,    -1,
      86,    87,    -1,    89,    90,    -1,    92,    -1,    94,    -1,
      96,    -1,    98,    -1,   100,    -1,   102,    -1,   104,    -1,
     106,    -1,   108,   109,    -1,   111,   112,   113,    -1,   115,
     116,    -1,   118,    -1,    -1,   121,   122,   123,    -1,   125,
       1,   127,     3,     4,     5,     6,     7,    -1,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    17,    -1,    19,    20,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    36,    37,    -1,    39,    40,
      -1,    42,    -1,    44,    -1,    -1,    47,    -1,    -1,    -1,
      -1,    52,    53,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    62,    63,    -1,    65,    -1,    -1,    -1,    69,    70,
      71,    -1,    73,    74,    75,    -1,    77,    78,    -1,    80,
      -1,    82,    83,    84,    -1,    86,    87,    -1,    89,    90,
      -1,    92,    -1,    94,    -1,    96,    -1,    98,    -1,   100,
      -1,   102,    -1,   104,    -1,   106,    -1,   108,   109,    -1,
     111,   112,   113,    -1,   115,   116,    -1,   118,    -1,    -1,
     121,   122,   123,    -1,   125,     1,   127,     3,     4,     5,
       6,     7,    -1,     9,    10,    11,    12,    13,    -1,    -1,
      -1,    17,    -1,    19,    20,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    -1,    39,    40,    -1,    42,    -1,    44,    -1,
      -1,    47,    -1,    -1,    -1,    -1,    52,    53,    54,    -1,
      56,    57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,
      -1,    -1,    -1,    69,    70,    71,    -1,    73,    74,    75,
      -1,    -1,    78,    -1,    80,    -1,    82,    83,    84,    -1,
      86,    87,    -1,    89,    90,    -1,    92,    -1,    94,    -1,
      96,    -1,    98,    -1,   100,    -1,   102,    -1,   104,    -1,
     106,    -1,   108,   109,    -1,   111,   112,   113,    -1,   115,
     116,    -1,   118,    -1,    -1,   121,   122,   123,    -1,   125,
       1,   127,     3,     4,     5,     6,     7,    -1,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    17,    -1,    19,    20,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    36,    37,    -1,    39,    40,
      -1,    42,    -1,    44,    -1,    -1,    47,    -1,    -1,    -1,
      -1,    52,    53,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    62,    63,    -1,    65,    -1,    -1,    -1,    69,    70,
      71,    -1,    73,    74,    75,    -1,    77,    78,    -1,    80,
      -1,    82,    83,    84,    -1,    86,    87,    -1,    89,    90,
      -1,    92,    -1,    94,    -1,    96,    -1,    98,    -1,   100,
      -1,   102,    -1,   104,    -1,   106,    -1,   108,   109,    -1,
     111,   112,   113,    -1,   115,   116,    -1,   118,    -1,    -1,
     121,   122,   123,    -1,   125,     1,   127,     3,     4,     5,
       6,     7,    -1,     9,    10,    11,    12,    13,    -1,    -1,
      -1,    17,    -1,    19,    20,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    -1,    39,    40,    -1,    42,    -1,    44,    -1,
      -1,    47,    -1,    -1,    -1,    -1,    52,    53,    54,    -1,
      -1,    57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,
      -1,    -1,    -1,    69,    70,    71,    -1,    73,    74,    75,
      -1,    77,    78,    -1,    80,    -1,    82,    83,    84,    -1,
      86,    87,    -1,    89,    90,    -1,    92,    -1,    94,    -1,
      96,    -1,    98,    -1,   100,    -1,   102,    -1,   104,    -1,
     106,    -1,   108,   109,    -1,   111,   112,   113,    -1,   115,
     116,    -1,   118,    -1,    -1,   121,   122,   123,    -1,   125,
       1,   127,     3,     4,     5,     6,     7,    -1,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    17,    -1,    19,    20,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    36,    37,    -1,    39,    40,
      -1,    42,    -1,    -1,    -1,    -1,    47,    -1,    -1,    -1,
      -1,    52,    53,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    62,    63,    -1,    65,    -1,    -1,    -1,    69,    70,
      71,    -1,    73,    74,    75,    76,    -1,    78,    79,    80,
      -1,    82,    83,    84,    -1,    86,    87,    -1,    89,    90,
      -1,    92,    -1,    94,    -1,    96,    -1,    98,    -1,   100,
      -1,   102,    -1,   104,    -1,   106,    -1,   108,   109,    -1,
     111,   112,   113,    -1,   115,   116,    -1,   118,    -1,    -1,
     121,   122,   123,    -1,   125,     1,   127,     3,     4,     5,
       6,     7,    -1,     9,    10,    11,    12,    13,    -1,    -1,
      -1,    17,    -1,    19,    20,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    -1,    39,    40,    -1,    42,    -1,    44,    -1,
      -1,    47,    -1,    -1,    -1,    -1,    52,    53,    54,    -1,
      -1,    57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,
      -1,    -1,    -1,    69,    70,    71,    -1,    73,    74,    75,
      -1,    77,    78,    -1,    80,    -1,    82,    83,    84,    -1,
      86,    87,    -1,    89,    90,    -1,    92,    -1,    94,    -1,
      96,    -1,    98,    -1,   100,    -1,   102,    -1,   104,    -1,
     106,    -1,   108,   109,    -1,   111,   112,   113,    -1,   115,
     116,    -1,   118,    -1,    -1,   121,   122,   123,    -1,   125,
       1,   127,     3,     4,     5,     6,     7,    -1,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    17,    -1,    19,    20,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    36,    37,    -1,    39,    40,
      -1,    42,    -1,    44,    -1,    -1,    47,    -1,    -1,    -1,
      -1,    52,    53,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    62,    63,    -1,    65,    -1,    -1,    -1,    69,    70,
      71,    -1,    73,    74,    75,    -1,    77,    78,    -1,    80,
      -1,    82,    83,    84,    -1,    86,    87,    -1,    89,    90,
      -1,    92,    -1,    94,    -1,    96,    -1,    98,    -1,   100,
      -1,   102,    -1,   104,    -1,   106,    -1,   108,   109,    -1,
     111,   112,   113,    -1,   115,   116,    -1,   118,    -1,    -1,
     121,   122,   123,    -1,   125,     1,   127,     3,     4,     5,
       6,     7,    -1,     9,    10,    11,    12,    13,    -1,    -1,
      -1,    17,    -1,    19,    20,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    -1,    39,    40,    -1,    42,    -1,    44,    -1,
      -1,    47,    -1,    -1,    -1,    -1,    52,    53,    54,    -1,
      -1,    57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,
      -1,    -1,    -1,    69,    70,    71,    -1,    73,    74,    75,
      -1,    77,    78,    -1,    80,    -1,    82,    83,    84,    -1,
      86,    87,    -1,    89,    90,    -1,    92,    -1,    94,    -1,
      96,    -1,    98,    -1,   100,    -1,   102,    -1,   104,    -1,
     106,    -1,   108,   109,    -1,   111,   112,   113,    -1,   115,
     116,    -1,   118,    -1,    -1,   121,   122,   123,    -1,   125,
       1,   127,     3,     4,     5,     6,     7,    -1,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    17,    -1,    19,    20,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    36,    37,    -1,    39,    40,
      -1,    42,    -1,    44,    -1,    -1,    47,    -1,    -1,    -1,
      -1,    52,    53,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    62,    63,    -1,    65,    -1,    -1,    -1,    69,    70,
      71,    -1,    73,    74,    75,    -1,    77,    78,    -1,    80,
      -1,    82,    83,    84,    -1,    86,    87,    -1,    89,    90,
      -1,    92,    -1,    94,    -1,    96,    -1,    98,    -1,   100,
      -1,   102,    -1,   104,    -1,   106,    -1,   108,   109,    -1,
     111,   112,   113,    -1,   115,   116,    -1,   118,    -1,    -1,
     121,   122,   123,    -1,   125,     1,   127,     3,     4,     5,
       6,     7,    -1,     9,    10,    11,    12,    13,    -1,    -1,
      -1,    17,    -1,    19,    20,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    -1,    39,    40,    -1,    42,    -1,    44,    -1,
      -1,    47,    -1,    -1,    -1,    -1,    52,    53,    54,    -1,
      -1,    57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,
      -1,    -1,    -1,    69,    70,    71,    -1,    73,    74,    75,
      -1,    77,    78,    -1,    80,    -1,    82,    83,    84,    -1,
      86,    87,    -1,    89,    90,    -1,    92,    -1,    94,    -1,
      96,    -1,    98,    -1,   100,    -1,   102,    -1,   104,    -1,
     106,    -1,   108,   109,    -1,   111,   112,   113,    -1,   115,
     116,    -1,   118,    -1,    -1,   121,   122,   123,    -1,   125,
       1,   127,     3,     4,     5,     6,     7,    -1,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    17,    -1,    19,    20,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    36,    37,    -1,    39,    40,
      -1,    42,    -1,    44,    -1,    -1,    47,    -1,    -1,    -1,
      -1,    52,    53,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    62,    63,    -1,    65,    -1,    -1,    -1,    69,    70,
      71,    -1,    73,    74,    75,    -1,    77,    78,    -1,    80,
      -1,    82,    83,    84,    -1,    86,    87,    -1,    89,    90,
      -1,    92,    -1,    94,    -1,    96,    -1,    98,    -1,   100,
      -1,   102,    -1,   104,    -1,   106,    -1,   108,   109,    -1,
     111,   112,   113,    -1,   115,   116,    -1,   118,    -1,    -1,
     121,   122,   123,    -1,   125,     1,   127,     3,     4,     5,
       6,     7,    -1,     9,    10,    11,    12,    13,    -1,    -1,
      -1,    17,    -1,    19,    20,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    -1,    39,    40,    -1,    42,    -1,    44,    -1,
      -1,    47,    -1,    -1,    -1,    -1,    52,    53,    54,    -1,
      -1,    57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,
      -1,    -1,    -1,    69,    70,    71,    -1,    73,    74,    75,
      -1,    77,    78,    -1,    80,    -1,    82,    83,    84,    -1,
      86,    87,    -1,    89,    90,    -1,    92,    -1,    94,    -1,
      96,    -1,    98,    -1,   100,    -1,   102,    -1,   104,    -1,
     106,    -1,   108,   109,    -1,   111,   112,   113,    -1,   115,
     116,    -1,   118,    -1,    -1,   121,   122,   123,    -1,   125,
       1,   127,     3,     4,     5,     6,     7,    -1,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    17,    -1,    19,    20,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    36,    37,    -1,    39,    40,
      -1,    42,    -1,    -1,    -1,    -1,    47,    -1,    -1,    -1,
      -1,    52,    53,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    62,    63,    -1,    65,    -1,    -1,    -1,    69,    70,
      71,    -1,    73,    74,    75,    -1,    77,    78,    -1,    80,
      -1,    82,    83,    84,    -1,    86,    87,    -1,    89,    90,
      -1,    92,    -1,    94,    -1,    96,    -1,    98,    -1,   100,
      -1,   102,    -1,   104,    -1,   106,    -1,   108,   109,   110,
     111,   112,   113,    -1,   115,   116,    -1,   118,    -1,    -1,
     121,   122,   123,    -1,   125,     1,   127,     3,     4,     5,
       6,     7,    -1,     9,    10,    11,    12,    13,    -1,    -1,
      -1,    17,    -1,    19,    20,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    -1,    39,    40,    41,    42,    -1,    -1,    -1,
      -1,    47,    -1,    -1,    -1,    -1,    52,    53,    54,    -1,
      -1,    57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,
      -1,    -1,    -1,    69,    70,    71,    -1,    73,    74,    75,
      -1,    -1,    78,    -1,    80,    -1,    82,    83,    84,    -1,
      86,    87,    -1,    89,    90,    -1,    92,    -1,    94,    -1,
      96,    -1,    98,    -1,   100,    -1,   102,    -1,   104,    -1,
     106,    -1,   108,   109,    -1,   111,   112,   113,    -1,   115,
     116,    -1,   118,    -1,    -1,   121,   122,   123,    -1,   125,
       1,   127,     3,     4,     5,     6,     7,    -1,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    17,    -1,    19,    20,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    36,    37,    -1,    39,    40,
      -1,    42,    43,    -1,    -1,    -1,    47,    -1,    -1,    -1,
      -1,    52,    53,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    62,    63,    -1,    65,    -1,    -1,    -1,    69,    70,
      71,    -1,    73,    74,    75,    -1,    -1,    78,    -1,    80,
      -1,    82,    83,    84,    -1,    86,    87,    -1,    89,    90,
      -1,    92,    -1,    94,    -1,    96,    -1,    98,    -1,   100,
      -1,   102,    -1,   104,    -1,   106,    -1,   108,   109,    -1,
     111,   112,   113,    -1,   115,   116,    -1,   118,    -1,    -1,
     121,   122,   123,    -1,   125,     1,   127,     3,     4,     5,
       6,     7,    -1,     9,    10,    11,    12,    13,    -1,    -1,
      -1,    17,    -1,    19,    20,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    -1,    39,    40,    -1,    42,    -1,    -1,    -1,
      -1,    47,    -1,    -1,    -1,    -1,    52,    53,    54,    -1,
      -1,    57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,
      -1,    -1,    -1,    69,    70,    71,    -1,    73,    74,    75,
      -1,    -1,    78,    -1,    80,    -1,    82,    83,    84,    -1,
      86,    87,    -1,    89,    90,    91,    92,    -1,    94,    -1,
      96,    -1,    98,    -1,   100,    -1,   102,    -1,   104,    -1,
     106,    -1,   108,   109,    -1,   111,   112,   113,    -1,   115,
     116,    -1,   118,    -1,    -1,   121,   122,   123,    -1,   125,
       1,   127,     3,     4,     5,     6,     7,    -1,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    17,    -1,    19,    20,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    36,    37,    -1,    39,    40,
      -1,    42,    -1,    -1,    -1,    -1,    47,    -1,    -1,    -1,
      -1,    52,    53,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    62,    63,    -1,    65,    -1,    -1,    -1,    69,    70,
      71,    -1,    73,    74,    75,    -1,    -1,    78,    -1,    80,
      -1,    82,    83,    84,    -1,    86,    87,    -1,    89,    90,
      -1,    92,    93,    94,    -1,    96,    -1,    98,    -1,   100,
      -1,   102,    -1,   104,    -1,   106,    -1,   108,   109,    -1,
     111,   112,   113,    -1,   115,   116,    -1,   118,    -1,    -1,
     121,   122,   123,    -1,   125,     1,   127,     3,     4,     5,
       6,     7,    -1,     9,    10,    11,    12,    13,    -1,    -1,
      -1,    17,    -1,    19,    20,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    -1,    39,    40,    -1,    42,    -1,    -1,    -1,
      -1,    47,    -1,    -1,    -1,    -1,    52,    53,    54,    -1,
      -1,    57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,
      -1,    -1,    -1,    69,    70,    71,    -1,    73,    74,    75,
      -1,    -1,    78,    -1,    80,    -1,    82,    83,    84,    -1,
      86,    87,    -1,    89,    90,    -1,    92,    -1,    94,    95,
      96,    -1,    98,    -1,   100,    -1,   102,    -1,   104,    -1,
     106,    -1,   108,   109,    -1,   111,   112,   113,    -1,   115,
     116,    -1,   118,    -1,    -1,   121,   122,   123,    -1,   125,
       1,   127,     3,     4,     5,     6,     7,    -1,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    17,    -1,    19,    20,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    36,    37,    -1,    39,    40,
      -1,    42,    -1,    -1,    -1,    -1,    47,    -1,    -1,    -1,
      -1,    52,    53,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    62,    63,    -1,    65,    -1,    -1,    -1,    69,    70,
      71,    -1,    73,    74,    75,    -1,    -1,    78,    -1,    80,
      -1,    82,    83,    84,    -1,    86,    87,    -1,    89,    90,
      -1,    92,    -1,    94,    -1,    96,    97,    98,    -1,   100,
      -1,   102,    -1,   104,    -1,   106,    -1,   108,   109,    -1,
     111,   112,   113,    -1,   115,   116,    -1,   118,    -1,    -1,
     121,   122,   123,    -1,   125,     1,   127,     3,     4,     5,
       6,     7,    -1,     9,    10,    11,    12,    13,    -1,    -1,
      -1,    17,    -1,    19,    20,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    -1,    39,    40,    -1,    42,    -1,    -1,    -1,
      -1,    47,    -1,    -1,    -1,    -1,    52,    53,    54,    -1,
      -1,    57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,
      -1,    -1,    -1,    69,    70,    71,    -1,    73,    74,    75,
      -1,    -1,    78,    -1,    80,    -1,    82,    83,    84,    -1,
      86,    87,    -1,    89,    90,    -1,    92,    -1,    94,    -1,
      96,    -1,    98,    -1,   100,    -1,   102,   103,   104,    -1,
     106,    -1,   108,   109,    -1,   111,   112,   113,    -1,   115,
     116,    -1,   118,    -1,    -1,   121,   122,   123,    -1,   125,
       1,   127,     3,     4,     5,     6,     7,    -1,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    17,    -1,    19,    20,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    36,    37,    -1,    39,    40,
      -1,    42,    -1,    -1,    -1,    -1,    47,    -1,    -1,    -1,
      -1,    52,    53,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    62,    63,    -1,    65,    -1,    -1,    -1,    69,    70,
      71,    -1,    73,    74,    75,    -1,    -1,    78,    -1,    80,
      -1,    82,    83,    84,    -1,    86,    87,    -1,    89,    90,
      -1,    92,    -1,    94,    -1,    96,    -1,    98,    -1,   100,
      -1,   102,    -1,   104,   105,   106,    -1,   108,   109,    -1,
     111,   112,   113,    -1,   115,   116,    -1,   118,    -1,    -1,
     121,   122,   123,    -1,   125,     1,   127,     3,     4,     5,
       6,     7,    -1,     9,    10,    11,    12,    13,    -1,    -1,
      -1,    17,    -1,    19,    20,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    -1,    39,    40,    -1,    42,    -1,    -1,    -1,
      -1,    47,    -1,    -1,    -1,    -1,    52,    53,    54,    -1,
      -1,    57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,
      -1,    -1,    -1,    69,    70,    71,    -1,    73,    74,    75,
      -1,    -1,    78,    -1,    80,    -1,    82,    83,    84,    -1,
      86,    87,    -1,    89,    90,    -1,    92,    -1,    94,    -1,
      96,    -1,    98,    -1,   100,    -1,   102,    -1,   104,    -1,
     106,   107,   108,   109,    -1,   111,   112,   113,    -1,   115,
     116,    -1,   118,    -1,    -1,   121,   122,   123,    -1,   125,
       1,   127,     3,     4,     5,     6,     7,    -1,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    17,    -1,    19,    20,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    36,    37,    -1,    39,    40,
      -1,    42,    -1,    -1,    -1,    -1,    47,    -1,    -1,    -1,
      -1,    52,    53,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    62,    63,    -1,    65,    -1,    -1,    -1,    69,    70,
      71,    -1,    73,    74,    75,    76,    -1,    78,    -1,    80,
      -1,    82,    83,    84,    -1,    86,    87,    -1,    89,    90,
      -1,    92,    -1,    94,    -1,    96,    -1,    98,    -1,   100,
      -1,   102,    -1,   104,    -1,   106,    -1,   108,   109,    -1,
     111,   112,   113,    -1,   115,   116,    -1,   118,    -1,    -1,
     121,   122,   123,    -1,   125,     1,   127,     3,     4,     5,
       6,     7,    -1,     9,    10,    11,    12,    13,    -1,    -1,
      -1,    17,    18,    19,    20,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    -1,    39,    40,    -1,    42,    -1,    -1,    -1,
      -1,    47,    -1,    -1,    -1,    -1,    52,    53,    54,    -1,
      -1,    57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,
      -1,    -1,    -1,    69,    70,    71,    -1,    73,    74,    75,
      -1,    -1,    78,    -1,    80,    -1,    82,    83,    84,    -1,
      86,    87,    -1,    89,    90,    -1,    92,    -1,    94,    -1,
      96,    -1,    98,    -1,   100,    -1,   102,    -1,   104,    -1,
     106,    -1,   108,   109,    -1,   111,   112,   113,    -1,   115,
     116,    -1,   118,    -1,    -1,   121,   122,   123,    -1,   125,
       1,   127,     3,     4,     5,     6,     7,    -1,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    17,    -1,    19,    20,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    36,    37,    -1,    39,    40,
      -1,    42,    -1,    44,    -1,    -1,    47,    -1,    -1,    -1,
      -1,    52,    53,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    62,    63,    -1,    65,    -1,    -1,    -1,    69,    70,
      71,    -1,    73,    74,    75,    -1,    -1,    78,    -1,    80,
      -1,    82,    83,    84,    -1,    86,    87,    -1,    89,    90,
      -1,    92,    -1,    94,    -1,    96,    -1,    98,    -1,   100,
      -1,   102,    -1,   104,    -1,   106,    -1,   108,   109,    -1,
     111,   112,   113,    -1,   115,   116,    -1,   118,    -1,    -1,
     121,   122,   123,    -1,   125,     1,   127,     3,     4,     5,
       6,     7,    -1,     9,    10,    11,    12,    13,    -1,    -1,
      -1,    17,    18,    19,    20,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    -1,    39,    40,    -1,    42,    -1,    -1,    -1,
      -1,    47,    -1,    -1,    -1,    -1,    52,    53,    54,    -1,
      -1,    57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,
      -1,    -1,    -1,    69,    70,    71,    -1,    73,    74,    75,
      -1,    -1,    78,    -1,    80,    -1,    82,    83,    84,    -1,
      86,    87,    -1,    89,    90,    -1,    92,    -1,    94,    -1,
      96,    -1,    98,    -1,   100,    -1,   102,    -1,   104,    -1,
     106,    -1,   108,   109,    -1,   111,   112,   113,    -1,   115,
     116,    -1,   118,    -1,    -1,   121,   122,   123,    -1,   125,
       1,   127,     3,     4,     5,     6,     7,    -1,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    17,    18,    19,    20,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    36,    37,    -1,    39,    40,
      -1,    42,    -1,    -1,    -1,    -1,    47,    -1,    -1,    -1,
      -1,    52,    53,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    62,    63,    -1,    65,    -1,    -1,    -1,    69,    70,
      71,    -1,    73,    74,    75,    -1,    -1,    78,    -1,    80,
      -1,    82,    83,    84,    -1,    86,    87,    -1,    89,    90,
      -1,    92,    -1,    94,    -1,    96,    -1,    98,    -1,   100,
      -1,   102,    -1,   104,    -1,   106,    -1,   108,   109,    -1,
     111,   112,   113,    -1,   115,   116,    -1,   118,    -1,    -1,
     121,   122,   123,    -1,   125,     1,   127,     3,     4,     5,
       6,     7,     8,     9,    10,    11,    12,    13,    -1,    -1,
      -1,    17,    -1,    19,    20,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    -1,    39,    40,    -1,    42,    -1,    -1,    -1,
      -1,    47,    -1,    -1,    -1,    -1,    52,    53,    54,    -1,
      -1,    57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,
      -1,    -1,    -1,    69,    70,    71,    -1,    73,    74,    75,
      -1,    -1,    78,    -1,    80,    -1,    82,    83,    84,    -1,
      86,    87,    -1,    89,    90,    -1,    92,    -1,    94,    -1,
      96,    -1,    98,    -1,   100,    -1,   102,    -1,   104,    -1,
     106,    -1,   108,   109,    -1,   111,   112,   113,    -1,   115,
     116,    -1,   118,    -1,    -1,   121,   122,   123,    -1,   125,
       1,   127,     3,     4,     5,     6,     7,    -1,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    17,    -1,    19,    20,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    36,    37,    -1,    39,    40,
      -1,    42,    -1,    -1,    -1,    -1,    47,    -1,    -1,    -1,
      -1,    52,    53,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    62,    63,    -1,    65,    -1,    -1,    -1,    69,    70,
      71,    -1,    73,    74,    75,    -1,    -1,    78,    79,    80,
      -1,    82,    83,    84,    -1,    86,    87,    -1,    89,    90,
      -1,    92,    -1,    94,    -1,    96,    -1,    98,    -1,   100,
      -1,   102,    -1,   104,    -1,   106,    -1,   108,   109,    -1,
     111,   112,   113,    -1,   115,   116,    -1,   118,    -1,    -1,
     121,   122,   123,    -1,   125,     1,   127,     3,     4,     5,
       6,     7,    -1,     9,    10,    11,    12,    13,    -1,    -1,
      -1,    17,    -1,    19,    20,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    -1,    39,    40,    -1,    42,    -1,    -1,    -1,
      -1,    47,    -1,    -1,    -1,    -1,    52,    53,    54,    -1,
      -1,    57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,
      -1,    -1,    -1,    69,    70,    71,    -1,    73,    74,    75,
      -1,    77,    78,    -1,    80,    -1,    82,    83,    84,    -1,
      86,    87,    -1,    89,    90,    -1,    92,    -1,    94,    -1,
      96,    -1,    98,    -1,   100,    -1,   102,    -1,   104,    -1,
     106,    -1,   108,   109,    -1,   111,   112,   113,    -1,   115,
     116,    -1,   118,    -1,    -1,   121,   122,   123,    -1,   125,
       1,   127,     3,     4,     5,     6,     7,    -1,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    17,    -1,    19,    20,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    36,    37,    -1,    39,    40,
      -1,    42,    -1,    -1,    -1,    -1,    47,    -1,    -1,    -1,
      -1,    52,    53,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    62,    63,    -1,    65,    -1,    -1,    -1,    69,    70,
      71,    -1,    73,    74,    75,    -1,    77,    78,    -1,    80,
      -1,    82,    83,    84,    -1,    86,    87,    -1,    89,    90,
      -1,    92,    -1,    94,    -1,    96,    -1,    98,    -1,   100,
      -1,   102,    -1,   104,    -1,   106,    -1,   108,   109,    -1,
     111,   112,   113,    -1,   115,   116,    -1,   118,    -1,    -1,
     121,   122,   123,    -1,   125,     1,   127,     3,     4,     5,
       6,     7,    -1,     9,    10,    11,    12,    13,    -1,    -1,
      -1,    17,    -1,    19,    20,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    -1,    39,    40,    -1,    42,    -1,    -1,    -1,
      -1,    47,    -1,    -1,    -1,    51,    52,    53,    54,    -1,
      -1,    57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,
      -1,    -1,    -1,    69,    70,    71,    -1,    73,    74,    75,
      -1,    -1,    78,    -1,    80,    -1,    82,    83,    84,    -1,
      86,    87,    -1,    89,    90,    -1,    92,    -1,    94,    -1,
      96,    -1,    98,    -1,   100,    -1,   102,    -1,   104,    -1,
     106,    -1,   108,   109,    -1,   111,   112,   113,    -1,   115,
     116,    -1,   118,    -1,    -1,   121,   122,   123,    -1,   125,
       1,   127,     3,     4,     5,     6,     7,    -1,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    17,    18,    19,    20,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    36,    37,    -1,    39,    40,
      -1,    42,    -1,    -1,    -1,    -1,    47,    -1,    -1,    -1,
      -1,    52,    53,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    62,    63,    -1,    65,    -1,    -1,    -1,    69,    70,
      71,    -1,    73,    74,    75,    -1,    -1,    78,    -1,    80,
      -1,    82,    83,    84,    -1,    86,    87,    -1,    89,    90,
      -1,    92,    -1,    94,    -1,    96,    -1,    98,    -1,   100,
      -1,   102,    -1,   104,    -1,   106,    -1,   108,   109,    -1,
     111,   112,   113,    -1,   115,   116,    -1,   118,    -1,    -1,
     121,   122,   123,    -1,   125,     1,   127,     3,     4,     5,
       6,     7,    -1,     9,    10,    11,    12,    13,    -1,    -1,
      -1,    17,    -1,    19,    20,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    -1,    39,    40,    -1,    42,    -1,    -1,    -1,
      -1,    47,    -1,    -1,    -1,    -1,    52,    53,    54,    -1,
      -1,    57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,
      -1,    -1,    -1,    69,    70,    71,    -1,    73,    74,    75,
      -1,    -1,    78,    -1,    80,    -1,    82,    83,    84,    -1,
      86,    87,    -1,    89,    90,    -1,    92,    -1,    94,    -1,
      96,    -1,    98,    -1,   100,    -1,   102,    -1,   104,    -1,
     106,    -1,   108,   109,    -1,   111,   112,   113,    -1,   115,
     116,    -1,   118,   119,    -1,   121,   122,   123,    -1,   125,
       1,   127,     3,     4,     5,     6,     7,    -1,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    17,    -1,    19,    20,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    36,    37,    -1,    39,    40,
      -1,    42,    -1,    -1,    -1,    -1,    47,    -1,    -1,    -1,
      -1,    52,    53,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    62,    63,    -1,    65,    -1,    -1,    -1,    69,    70,
      71,    -1,    73,    74,    75,    -1,    -1,    78,    -1,    80,
      -1,    82,    83,    84,    -1,    86,    87,    -1,    89,    90,
      -1,    92,    -1,    94,    -1,    96,    -1,    98,    -1,   100,
      -1,   102,    -1,   104,    -1,   106,    -1,   108,   109,   110,
     111,   112,   113,    -1,   115,   116,    -1,   118,    -1,    -1,
     121,   122,   123,    -1,   125,     1,   127,     3,     4,     5,
       6,     7,    -1,     9,    10,    11,    12,    13,    -1,    -1,
      -1,    17,    -1,    19,    20,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    -1,    39,    40,    -1,    42,    -1,    -1,    -1,
      -1,    47,    -1,    -1,    -1,    -1,    52,    53,    54,    -1,
      -1,    57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,
      -1,    -1,    -1,    69,    70,    71,    -1,    73,    74,    75,
      -1,    -1,    78,    -1,    80,    -1,    82,    83,    84,    -1,
      86,    87,    -1,    89,    90,    -1,    92,    -1,    94,    -1,
      96,    -1,    98,    -1,   100,    -1,   102,    -1,   104,    -1,
     106,    -1,   108,   109,   110,   111,   112,   113,    -1,   115,
     116,    -1,   118,    -1,    -1,   121,   122,   123,    -1,   125,
       1,   127,     3,     4,     5,     6,     7,    -1,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    17,    -1,    19,    20,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    36,    37,    -1,    39,    40,
      -1,    42,    -1,    -1,    -1,    -1,    47,    -1,    -1,    -1,
      -1,    52,    53,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    62,    63,    -1,    65,    -1,    -1,    -1,    69,    70,
      71,    -1,    73,    74,    75,    -1,    -1,    78,    -1,    80,
      -1,    82,    83,    84,    -1,    86,    87,    -1,    89,    90,
      -1,    92,    -1,    94,    -1,    96,    -1,    98,    -1,   100,
      -1,   102,    -1,   104,    -1,   106,    -1,   108,   109,    -1,
     111,   112,   113,    -1,   115,   116,    -1,   118,    -1,    -1,
     121,   122,   123,    -1,   125,     1,   127,     3,     4,     5,
       6,     7,    -1,     9,    10,    11,    12,    13,    -1,    -1,
      -1,    17,    -1,    19,    20,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    -1,    39,    40,    -1,    42,    -1,    -1,    -1,
      -1,    47,    -1,    -1,    -1,    -1,    52,    53,    54,    -1,
      -1,    57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,
      -1,    -1,    -1,    69,    70,    71,    -1,    73,    74,    75,
      -1,    -1,    78,    -1,    80,    -1,    82,    83,    84,    -1,
      86,    87,    -1,    89,    90,    -1,    92,    -1,    94,    -1,
      96,    -1,    98,    -1,   100,    -1,   102,    -1,   104,    -1,
     106,    -1,   108,   109,    -1,   111,   112,   113,    -1,   115,
     116,    -1,   118,    -1,    -1,   121,   122,   123,    -1,   125,
       1,   127,     3,     4,     5,     6,     7,    -1,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    17,    -1,    19,    20,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    36,    37,    -1,    39,    40,
      -1,    42,    -1,    -1,    -1,    -1,    47,    -1,    -1,    -1,
      -1,    52,    53,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    62,    63,    -1,    65,    -1,    -1,    -1,    69,    70,
      71,    -1,    73,    74,    75,    -1,    -1,    78,    -1,    80,
      -1,    82,    83,    84,    -1,    86,    87,    -1,    89,    90,
      -1,    92,    -1,    94,    -1,    96,    -1,    98,    -1,   100,
      -1,   102,    -1,   104,    -1,   106,    -1,   108,   109,    -1,
     111,   112,   113,    -1,   115,   116,    -1,   118,    -1,    -1,
     121,   122,   123,    -1,   125,     1,   127,     3,     4,     5,
       6,     7,    -1,     9,    10,    11,    12,    13,    -1,    -1,
      -1,    17,    -1,    19,    20,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    -1,    39,    40,    -1,    42,    -1,    -1,    -1,
      -1,    47,    -1,    -1,    -1,    -1,    52,    53,    54,    -1,
      -1,    57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,
      -1,    -1,    -1,    69,    70,    71,    -1,    73,    74,    75,
      -1,    -1,    78,    -1,    80,    -1,    82,    83,    84,    -1,
      86,    87,    -1,    89,    90,    -1,    92,    -1,    94,    -1,
      96,    -1,    98,    -1,   100,    -1,   102,    -1,   104,    -1,
     106,    -1,   108,   109,    -1,   111,   112,   113,    -1,   115,
     116,    -1,   118,    -1,    -1,   121,   122,   123,    -1,   125,
       1,   127,     3,     4,     5,     6,     7,    -1,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    17,    -1,    19,    20,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    36,    37,    -1,    39,    40,
      -1,    42,    -1,    -1,    -1,    -1,    47,    -1,    -1,    -1,
      -1,    52,    53,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    62,    63,    -1,    65,    -1,    -1,    -1,    69,    70,
      71,    -1,    73,    74,    75,    -1,    -1,    78,    -1,    80,
      -1,    82,    83,    84,    -1,    86,    87,    -1,    89,    90,
      -1,    92,    -1,    94,    -1,    96,    -1,    98,    -1,   100,
      -1,   102,    -1,   104,    -1,   106,    -1,   108,   109,    -1,
     111,   112,   113,    -1,   115,   116,    -1,   118,    -1,    -1,
     121,   122,   123,    -1,   125,     1,   127,     3,     4,     5,
       6,     7,    -1,     9,    10,    11,    12,    13,    -1,    -1,
      -1,    17,    -1,    19,    20,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    -1,    39,    40,    -1,    42,    -1,    -1,    -1,
      -1,    47,    -1,    -1,    -1,    -1,    52,    53,    54,    -1,
      -1,    57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,
      -1,    -1,    -1,    69,    70,    71,    -1,    73,    74,    75,
      -1,    -1,    78,    -1,    80,    -1,    82,    83,    84,    -1,
      86,    87,    -1,    89,    90,    -1,    92,    -1,    94,    -1,
      96,    -1,    98,    -1,   100,    -1,   102,    -1,   104,    -1,
     106,    -1,   108,   109,    -1,   111,   112,   113,    -1,   115,
     116,    -1,   118,    -1,    -1,   121,   122,   123,    -1,   125,
       1,   127,     3,     4,     5,     6,     7,    -1,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    17,    -1,    19,    20,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    36,    37,    -1,    39,    40,
      -1,    42,    -1,    -1,    -1,    -1,    47,    -1,    -1,    -1,
      -1,    52,    53,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    62,    63,    -1,    65,    -1,    -1,    -1,    69,    70,
      71,    -1,    73,    74,    75,    -1,    -1,    78,    -1,    80,
      -1,    82,    83,    84,    -1,    86,    87,    -1,    89,    90,
      -1,    92,    -1,    94,    -1,    96,    -1,    98,    -1,   100,
      -1,   102,    -1,   104,    -1,   106,    -1,   108,   109,    -1,
     111,   112,   113,    -1,   115,   116,    -1,   118,    -1,    -1,
     121,   122,   123,    -1,   125,     1,   127,     3,     4,     5,
       6,     7,    -1,     9,    10,    11,    12,    13,    -1,    -1,
      -1,    17,    -1,    19,    20,    -1,    -1,    -1,    -1,    25,
      26,    -1,    -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    -1,    39,    40,    -1,    42,    -1,    -1,    -1,
      -1,    47,    -1,    -1,    -1,    -1,    52,    53,    54,    -1,
      -1,    57,    -1,    -1,    60,    -1,    62,    63,    -1,    65,
      -1,    -1,    -1,    69,    70,    71,    -1,    73,    74,    75,
      -1,    -1,    78,    -1,    80,    -1,    82,    83,    84,    -1,
      86,    87,    -1,    89,    90,    -1,    92,    -1,    94,    -1,
      96,    -1,    98,    -1,   100,    -1,   102,    -1,   104,    -1,
     106,    -1,   108,   109,    -1,   111,   112,   113,    -1,   115,
     116,    -1,   118,    -1,    -1,   121,   122,   123,    -1,   125,
      -1,   127,     3,     4,    -1,    -1,     7,    -1,    -1,    10,
      11,    -1,    13,    14,    -1,    16,    -1,    -1,    19,    20,
      21,    22,    23,    24,    25,    26,    -1,    28,    29,    30,
      31,    32,    33,    34,    35,    36,    -1,    -1,    -1,    40,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    48,    -1,    -1,
      -1,    52,    -1,    54,    -1,    -1,    57,    -1,    -1,    60,
      -1,    -1,    -1,    64,    -1,    -1,    -1,    68,    -1,    -1,
      -1,    -1,    -1,    -1,     3,     4,    -1,    -1,     7,    -1,
      -1,    10,    11,    -1,    13,    86,    87,    16,    89,    -1,
      19,    20,    21,    22,    23,    24,    25,    26,    -1,    28,
      29,    30,    31,    32,    33,    34,    35,    36,    -1,    -1,
      -1,    40,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    48,
     121,    -1,    -1,    52,   125,    54,    -1,    -1,    57,    -1,
      -1,    60,    -1,    -1,    -1,    64,    -1,    -1,    -1,    68,
      -1,    -1,    -1,    -1,    -1,    -1,     3,     4,    77,    -1,
       7,    -1,    -1,    10,    11,    -1,    13,    86,    87,    16,
      89,    -1,    19,    20,    21,    22,    23,    24,    25,    26,
      -1,    28,    29,    30,    31,    32,    33,    34,    35,    36,
      -1,    -1,    -1,    40,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    48,   121,    -1,    -1,    52,   125,    54,    -1,    -1,
      57,    -1,    -1,    60,    -1,    -1,    -1,    64,    -1,    -1,
      -1,    68,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    86,
      87,    -1,    89,     8,    -1,    -1,    -1,    -1,    -1,    14,
      15,    16,    17,    -1,    19,    20,    21,    22,    23,    24,
      -1,    -1,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,   121,    -1,    -1,    -1,   125,    44,
      45,    -1,    11,    -1,    -1,    14,    15,    16,    17,    -1,
      19,    20,    21,    22,    23,    24,    -1,    -1,    27,    28,
      29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      -1,    -1,    -1,    -1,    -1,    -1,    45,    -1,    -1,    -1,
      85,    -1,    14,    15,    16,    17,    55,    19,    20,    21,
      22,    23,    24,    -1,    -1,    27,    28,    29,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    -1,    -1,   114,
      -1,    -1,   117,    45,    -1,   120,    85,    -1,    50,   124,
     125,    -1,    -1,   128,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   114,    -1,    -1,   117,    -1,
      -1,   120,    -1,    85,    -1,   124,   125,    -1,    -1,   128,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   109,    -1,    -1,
      -1,    -1,   114,    -1,    -1,   117,    -1,    -1,   120,    -1,
      -1,    -1,   124,   125,    -1,    -1,   128,    14,    15,    16,
      17,    18,    19,    20,    21,    22,    23,    24,    -1,    -1,
      27,    28,    29,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    -1,    -1,    14,    15,    16,    17,    45,    19,
      20,    21,    22,    23,    24,    -1,    26,    27,    28,    29,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    -1,
      -1,    -1,    -1,    -1,    -1,    45,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    85,    -1,
      14,    15,    16,    17,    -1,    19,    20,    21,    22,    23,
      24,    -1,    -1,    27,    28,    29,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    85,    -1,   114,    -1,    -1,
     117,    45,    -1,   120,    -1,    -1,    -1,   124,   125,    -1,
      -1,   128,    -1,    -1,    -1,    -1,    60,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   114,    -1,    -1,   117,    -1,    -1,
     120,    -1,    -1,    -1,   124,   125,    -1,    -1,   128,    -1,
      -1,    85,    -1,    14,    15,    16,    17,    -1,    19,    20,
      21,    22,    23,    24,    -1,    -1,    27,    28,    29,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    -1,    -1,
     114,    -1,    -1,   117,    45,    -1,   120,    -1,    -1,    -1,
     124,   125,    -1,    -1,   128,    14,    15,    16,    17,    -1,
      19,    20,    21,    22,    23,    24,    -1,    -1,    27,    28,
      29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      -1,    -1,    -1,    -1,    85,    -1,    45,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    99,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   114,    -1,    -1,   117,    -1,    -1,   120,
      -1,    -1,    -1,   124,   125,    -1,    85,   128,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     109,    -1,    -1,    -1,    -1,   114,    -1,    -1,   117,    -1,
      -1,   120,    -1,    -1,    -1,   124,   125,    -1,    -1,   128,
      14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    -1,    -1,    27,    28,    29,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    -1,    -1,    14,    15,    16,
      17,    45,    19,    20,    21,    22,    23,    24,    -1,    -1,
      27,    28,    29,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    -1,    -1,    41,    -1,    -1,    -1,    45,    -1,
      14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    85,    -1,    27,    28,    29,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    -1,    -1,    -1,    -1,    -1,
      -1,    45,    -1,    -1,    -1,    -1,    -1,    -1,    85,    -1,
     114,    -1,    -1,   117,    -1,    -1,   120,    -1,    -1,    -1,
     124,   125,    -1,    -1,   128,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   114,    -1,    -1,
     117,    85,    -1,   120,    -1,    -1,    -1,   124,   125,    -1,
      -1,   128,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     114,    -1,    -1,   117,    -1,    -1,   120,    -1,    -1,    -1,
     124,   125,    -1,    -1,   128,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,    -1,    -1,    27,    28,
      29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      -1,    -1,    -1,    -1,    -1,    -1,    45,    -1,    14,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    -1,
      -1,    27,    28,    29,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    -1,    -1,    -1,    -1,    -1,    -1,    45,
      -1,    -1,    -1,    -1,    -1,    -1,    85,    -1,    -1,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      -1,    -1,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    -1,   114,    -1,    -1,   117,    85,
      45,   120,    -1,    -1,    -1,   124,   125,    -1,    -1,   128,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   114,    -1,
      -1,   117,    -1,    -1,   120,    -1,    -1,    -1,   124,   125,
      85,    -1,   128,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,    -1,    -1,    27,    28,    29,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    -1,   114,
      -1,    -1,   117,    -1,    45,   120,    -1,    -1,    -1,   124,
     125,    -1,    -1,   128,    14,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    -1,    -1,    27,    28,    29,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    -1,
      -1,    -1,    -1,    -1,    85,    45,    -1,    14,    15,    16,
      17,    18,    19,    20,    21,    22,    23,    24,    -1,    -1,
      27,    28,    29,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    -1,   114,    -1,    -1,   117,    -1,    45,   120,
      -1,    -1,    -1,   124,   125,    85,    -1,   128,    14,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    -1,
      -1,    27,    28,    29,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    -1,   114,    -1,    -1,   117,    85,    45,
     120,    -1,    -1,    -1,   124,   125,    -1,    -1,   128,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   114,    -1,    -1,
     117,    -1,    -1,   120,    -1,    -1,    -1,   124,   125,    85,
      -1,   128,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    -1,    -1,    27,    28,    29,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    -1,   114,    -1,
      -1,   117,    -1,    45,   120,    -1,    -1,    -1,   124,   125,
      -1,    -1,   128,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,    -1,    -1,    27,    28,    29,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    -1,    -1,
      -1,    -1,    -1,    85,    45,    -1,    14,    15,    16,    17,
      -1,    19,    20,    21,    22,    23,    24,    -1,    -1,    27,
      28,    29,    30,    31,    32,    33,    34,    35,    36,    37,
      38,    -1,   114,    41,    -1,   117,    -1,    45,   120,    -1,
      -1,    -1,   124,   125,    85,    -1,   128,    14,    15,    16,
      17,    -1,    19,    20,    21,    22,    23,    24,    -1,    26,
      27,    28,    29,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    -1,   114,    -1,    -1,   117,    85,    45,   120,
      -1,    -1,    -1,   124,   125,    -1,    -1,   128,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   114,    -1,    -1,   117,
      -1,    -1,   120,    -1,    -1,    -1,   124,   125,    85,    -1,
     128,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    -1,    -1,    27,    28,    29,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    -1,   114,    -1,    -1,
     117,    -1,    45,   120,    -1,    -1,    -1,   124,   125,    -1,
      -1,   128,    14,    15,    16,    17,    -1,    19,    20,    21,
      22,    23,    24,    -1,    -1,    27,    28,    29,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    -1,    -1,    41,
      -1,    -1,    85,    45,    -1,    14,    15,    16,    17,    -1,
      19,    20,    21,    22,    23,    24,    -1,    -1,    27,    28,
      29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      -1,   114,    41,    -1,   117,    -1,    45,   120,    -1,    -1,
      -1,   124,   125,    85,    -1,   128,    14,    15,    16,    17,
      -1,    19,    20,    21,    22,    23,    24,    -1,    -1,    27,
      28,    29,    30,    31,    32,    33,    34,    35,    36,    37,
      38,    -1,   114,    -1,    -1,   117,    85,    45,   120,    -1,
      -1,    -1,   124,   125,    -1,    -1,   128,    -1,    -1,    -1,
      -1,    17,    60,    19,    20,    21,    22,    23,    24,    -1,
      -1,    27,    -1,    -1,    30,   114,    -1,    -1,   117,    35,
      36,   120,    38,    -1,    -1,   124,   125,    85,    -1,   128,
      14,    15,    16,    17,    -1,    19,    20,    21,    22,    23,
      24,    -1,    -1,    27,    28,    29,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    -1,   114,    -1,    -1,   117,
      44,    45,   120,    -1,    -1,    -1,   124,   125,    -1,    85,
     128,    -1,    14,    15,    16,    17,    -1,    19,    20,    21,
      22,    23,    24,    -1,    -1,    27,    28,    29,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    -1,    -1,    -1,
      -1,    85,    44,    45,    -1,    -1,    -1,    -1,    -1,   125,
      -1,    -1,   128,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     114,    -1,    -1,   117,    -1,    -1,   120,    -1,    -1,    -1,
     124,   125,    -1,    85,   128,    14,    15,    16,    17,    -1,
      19,    20,    21,    22,    23,    24,    -1,    -1,    27,    28,
      29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      -1,    -1,   114,    -1,    -1,   117,    45,    -1,   120,    -1,
      -1,    -1,   124,   125,    -1,    -1,   128,    -1,    -1,    -1,
      -1,    -1,    61,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    14,    15,    16,    17,    -1,    19,
      20,    21,    22,    23,    24,    -1,    85,    27,    28,    29,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    -1,
      -1,    -1,    -1,    -1,    44,    45,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   114,    -1,    -1,   117,    -1,
      -1,   120,    -1,    -1,    -1,   124,   125,    -1,    -1,   128,
      -1,    -1,    -1,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,    -1,    85,    27,    28,    29,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    -1,    -1,
      -1,    -1,    -1,    -1,    45,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   114,    -1,    -1,   117,    -1,    -1,
     120,    -1,    -1,    -1,   124,   125,    -1,    -1,   128,    -1,
      -1,    -1,    -1,    14,    15,    16,    17,    -1,    19,    20,
      21,    22,    23,    24,    85,    26,    27,    28,    29,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    -1,    -1,
      -1,    -1,    -1,    -1,    45,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   114,    -1,    -1,   117,    -1,    -1,   120,
      -1,    -1,    -1,   124,   125,    -1,    -1,   128,    -1,    -1,
      -1,    -1,    14,    15,    16,    17,    -1,    19,    20,    21,
      22,    23,    24,    -1,    85,    27,    28,    29,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    -1,    -1,    41,
      -1,    -1,    -1,    45,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   114,    -1,    -1,   117,    -1,    -1,   120,
      -1,    -1,    -1,   124,   125,    -1,    -1,   128,    -1,    -1,
      -1,    14,    15,    16,    17,    -1,    19,    20,    21,    22,
      23,    24,    -1,    85,    27,    28,    29,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    -1,    -1,    -1,    50,    -1,    -1,
      -1,    -1,   114,    -1,    -1,   117,    -1,    -1,   120,    -1,
      -1,    -1,   124,   125,    -1,    -1,   128,    -1,    -1,    -1,
      14,    15,    16,    17,    -1,    19,    20,    21,    22,    23,
      24,    -1,    85,    27,    28,    29,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    -1,    -1,    -1,    -1,    -1,
      -1,    45,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   114,    -1,    -1,   117,    -1,    60,   120,    -1,    -1,
      -1,   124,   125,    -1,    -1,   128,    -1,    -1,    -1,    14,
      15,    16,    17,    -1,    19,    20,    21,    22,    23,    24,
      -1,    85,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    -1,    -1,    -1,    -1,    -1,    -1,
      45,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     114,    -1,    -1,   117,    -1,    60,   120,    -1,    -1,    -1,
     124,   125,    -1,    -1,   128,    -1,    -1,    -1,    14,    15,
      16,    17,    -1,    19,    20,    21,    22,    23,    24,    -1,
      85,    27,    28,    29,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    -1,    -1,    -1,    -1,    -1,    44,    45,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   114,
      -1,    -1,   117,    -1,    -1,   120,    -1,    -1,    -1,   124,
     125,    -1,    -1,   128,    -1,    -1,    -1,    14,    15,    16,
      17,    -1,    19,    20,    21,    22,    23,    24,    -1,    85,
      27,    28,    29,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    -1,    -1,    -1,    -1,    -1,    -1,    45,    -1,
      -1,    -1,    -1,    50,    -1,    -1,    -1,    -1,   114,    -1,
      -1,   117,    -1,    -1,   120,    -1,    -1,    -1,   124,   125,
      -1,    -1,   128,    -1,    -1,    -1,    14,    15,    16,    17,
      -1,    19,    20,    21,    22,    23,    24,    -1,    85,    27,
      28,    29,    30,    31,    32,    33,    34,    35,    36,    37,
      38,    -1,    -1,    -1,    -1,    -1,    44,    45,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   114,    -1,    -1,
     117,    -1,    -1,   120,    -1,    -1,    -1,   124,   125,    -1,
      -1,   128,    -1,    -1,    -1,    14,    15,    16,    17,    -1,
      19,    20,    21,    22,    23,    24,    -1,    85,    27,    28,
      29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      -1,    -1,    -1,    -1,    -1,    -1,    45,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   114,    -1,    -1,   117,
      -1,    -1,   120,    -1,    -1,    -1,   124,   125,    -1,    -1,
     128,    -1,    -1,    -1,    14,    15,    16,    17,    -1,    19,
      20,    21,    22,    23,    24,    -1,    85,    27,    28,    29,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    -1,
      -1,    -1,    -1,    -1,    -1,    45,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   114,    -1,    -1,   117,    -1,
      -1,   120,    -1,    -1,    -1,   124,   125,    -1,    -1,   128,
      -1,    -1,    -1,    14,    15,    16,    17,    -1,    19,    20,
      21,    22,    23,    24,    -1,    85,    27,    28,    29,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    -1,    -1,
      -1,    -1,    -1,    -1,    45,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   114,    -1,    -1,   117,    -1,    -1,
     120,    -1,    -1,    -1,   124,   125,    -1,    -1,   128,    -1,
      -1,    -1,    14,    15,    16,    17,    -1,    19,    20,    21,
      22,    23,    24,    -1,    85,    27,    28,    29,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    -1,    -1,    -1,
      -1,    -1,    -1,    45,    -1,    -1,    -1,    17,    -1,    19,
      20,    21,    22,    23,    24,    -1,   117,    27,    -1,   120,
      -1,    -1,    -1,   124,   125,    35,    36,   128,    38,    -1,
      -1,    14,    -1,    16,    17,    -1,    19,    20,    21,    22,
      23,    24,    -1,    85,    27,    28,    29,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    -1,    -1,    -1,    17,    -1,    -1,    -1,
      21,    22,    23,    24,    -1,    85,    27,    -1,   120,    -1,
      -1,    -1,   124,   125,    35,    36,   128,    38,    -1,    -1,
      14,    -1,    16,    17,    -1,    19,    20,    21,    22,    23,
      24,    -1,    85,    27,    28,    29,    30,    31,    32,    33,
      34,    35,    36,    37,    38,   125,    -1,    -1,   128,    -1,
      -1,    45,    -1,    -1,    -1,    17,    -1,    -1,    -1,    21,
      22,    23,    24,    -1,    85,    27,    -1,    -1,    -1,    -1,
      -1,   124,   125,    35,    36,   128,    38,    -1,    -1,    -1,
      15,    16,    17,    -1,    19,    20,    21,    22,    23,    24,
      -1,    85,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,   125,    -1,    -1,   128,    -1,    -1,
      45,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    85,    -1,    -1,    -1,    -1,    -1,    -1,
     124,   125,    -1,    -1,   128,    -1,    -1,    -1,    -1,    -1,
      16,    17,    -1,    19,    20,    21,    22,    23,    24,    -1,
      85,    27,    28,    29,    30,    31,    32,    33,    34,    35,
      36,    37,    38,   125,    -1,    -1,   128,    -1,    -1,    45,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   114,
      -1,    -1,   117,    -1,    -1,   120,    -1,    -1,    -1,   124,
     125,    -1,    -1,   128,    -1,    -1,    -1,    -1,    -1,    16,
      17,    -1,    19,    20,    21,    22,    23,    24,    -1,    85,
      27,    28,    29,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    -1,    99,    -1,    -1,    16,    17,    45,    19,
      20,    21,    22,    23,    24,    -1,    -1,    27,    28,    29,
      30,    31,    32,    33,    34,    35,    36,    37,    38,   125,
      -1,    -1,   128,    -1,    -1,    45,    -1,    -1,    16,    17,
      -1,    19,    20,    21,    22,    23,    24,    -1,    85,    27,
      28,    29,    30,    31,    32,    33,    34,    35,    36,    37,
      38,    -1,    -1,    -1,   101,    -1,    17,    -1,    19,    20,
      21,    22,    23,    24,    -1,    85,    27,    28,    29,    30,
      31,    32,    33,    -1,    35,    36,    -1,    38,   125,    -1,
      -1,   128,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    85,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   125,    -1,    -1,   128,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    85,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   125,    -1,    -1,
     128,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   125,    -1,    -1,   128
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     1,     3,     4,     5,     6,     7,     9,    10,    11,
      12,    13,    17,    19,    20,    25,    26,    29,    36,    37,
      39,    40,    42,    47,    52,    53,    54,    57,    60,    62,
      63,    65,    69,    70,    71,    73,    74,    75,    78,    80,
      82,    83,    84,    86,    87,    89,    90,    92,    94,    96,
      98,   100,   102,   104,   106,   108,   109,   111,   112,   113,
     115,   116,   118,   121,   122,   123,   125,   127,   131,   132,
     133,   134,   136,     4,     5,    46,    15,    17,    38,   133,
       4,   133,   133,   133,    13,    78,   133,   133,     4,   133,
     133,   142,   142,   142,    17,   133,    17,    17,    17,   133,
      13,    17,   134,   133,   143,    17,   133,   143,    76,   151,
      17,    17,    17,   109,   133,    17,    17,    17,    38,     3,
       7,    10,    11,    16,    19,    20,    21,    22,    23,    24,
      25,    26,    28,    29,    30,    31,    32,    33,    34,    35,
      36,    40,    48,    52,    54,    57,    60,    64,    68,    86,
      87,   125,   134,   136,   144,   145,    17,    83,   142,   142,
     142,   142,   133,   133,   142,   142,   142,   133,   133,   134,
     143,   142,   133,   143,   143,   151,   134,   133,   145,    38,
     133,    17,   133,     0,     8,    14,    15,    16,    17,    19,
      20,    21,    22,    23,    24,    27,    28,    29,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    44,    45,    85,
     114,   117,   120,   124,   125,   128,    17,    24,   128,   133,
     133,    18,   133,   133,    18,    26,    26,    26,    41,    41,
      43,   133,    50,   109,   133,   133,     4,    11,    55,   157,
     133,   141,    58,    67,   149,    44,    77,   133,   158,   133,
      60,    66,    81,   135,   138,   143,    72,   133,   142,   133,
     133,   133,   133,   133,   125,   134,   147,    10,    86,   133,
     144,   144,   155,   134,   144,   144,   144,   144,   145,    88,
      14,   144,   133,    91,    93,    95,    97,    14,    99,    14,
     103,   105,   107,    17,    77,   109,   153,   109,   109,    51,
     109,   153,    17,   122,   133,    18,   133,   133,   133,   151,
     142,   133,   133,   133,   133,   133,   133,   133,   133,   133,
     133,   133,   133,   133,   133,   133,   133,   134,   133,     8,
     132,   133,   134,   133,   133,   133,   133,   142,   133,    18,
      41,    17,    18,   143,   151,   133,    18,    18,    18,     3,
     143,    56,    18,    44,   133,   133,    59,    61,   148,    18,
     143,   133,   133,    79,   137,   143,   135,    77,    17,    18,
      18,    18,    18,    18,    18,    41,   145,    26,    87,   134,
     146,    26,    26,    26,    41,    50,    77,   144,   144,    66,
      65,    61,    70,    88,   144,    18,   133,   133,   142,    77,
     143,   143,   143,   153,    51,   142,    41,    18,    67,    41,
      17,    38,    18,   142,    17,   133,   151,    49,    51,   109,
     152,   153,   154,    51,    76,    60,   133,   141,    60,   133,
     133,    59,    60,    63,   133,   151,    77,    44,     3,     4,
      10,    13,    17,   134,   139,   140,    77,   133,   135,   133,
     144,   144,    50,    77,    61,    70,   144,   144,    99,   101,
      18,   153,   153,   153,   143,   153,    18,   133,   142,   133,
      15,    18,   133,    44,   109,   150,   154,   133,    49,    51,
     143,   151,    54,    55,   156,   143,    44,   143,   148,   133,
     143,   133,    44,    15,    46,   140,    44,    14,    18,    76,
     143,    26,    88,    51,    77,   144,    77,    77,   119,   153,
      29,    18,    41,   133,   151,   150,   133,   151,    50,   133,
     153,     3,    45,    77,    11,   157,   141,    77,    60,    60,
      77,    60,   133,   133,    18,   139,   151,   143,    77,   144,
      77,   109,   143,   119,   133,    44,   143,    50,    45,   151,
      18,   143,   143,   143,    77,    77,   110,   109,   110,   153,
     109,   143,   152,   143,   151,   133,   151,    77,    77,    77,
     142,   110,   142,   110,   153,   152,   156,    44,   109,   142,
     109,   142,   143,   109,   143,   109,   153,   143,   153,   143,
     153,   153
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (scanner, YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (&yylval, YYLEX_PARAM)
#else
# define YYLEX yylex (&yylval)
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value, scanner); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, void * scanner)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep, scanner)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
    void * scanner;
#endif
{
  if (!yyvaluep)
    return;
  YYUSE (scanner);
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, void * scanner)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep, scanner)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
    void * scanner;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep, scanner);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *bottom, yytype_int16 *top)
#else
static void
yy_stack_print (bottom, top)
    yytype_int16 *bottom;
    yytype_int16 *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule, void * scanner)
#else
static void
yy_reduce_print (yyvsp, yyrule, scanner)
    YYSTYPE *yyvsp;
    int yyrule;
    void * scanner;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      fprintf (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       , scanner);
      fprintf (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule, scanner); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, void * scanner)
#else
static void
yydestruct (yymsg, yytype, yyvaluep, scanner)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
    void * scanner;
#endif
{
  YYUSE (yyvaluep);
  YYUSE (scanner);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void * scanner);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */






/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void * scanner)
#else
int
yyparse (scanner)
    void * scanner;
#endif
#endif
{
  /* The look-ahead symbol.  */
int yychar;

/* The semantic value of the look-ahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;

  int yystate;
  int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int yytoken = 0;
#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  yytype_int16 yyssa[YYINITDEPTH];
  yytype_int16 *yyss = yyssa;
  yytype_int16 *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  YYSTYPE *yyvsp;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     look-ahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to look-ahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a look-ahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
#line 173 "input_parser.yy"
    {   const giac::context * contextptr = giac_yyget_extra(scanner);
			    if ((yyvsp[(1) - (1)])._VECTptr->size()==1)
			     parsed_gen((yyvsp[(1) - (1)])._VECTptr->front(),contextptr);
                          else
			     parsed_gen(gen(*(yyvsp[(1) - (1)])._VECTptr,_SEQ__VECT),contextptr);
			 }
    break;

  case 3:
#line 181 "input_parser.yy"
    { (yyval)=vecteur(1,(yyvsp[(1) - (2)])); }
    break;

  case 4:
#line 182 "input_parser.yy"
    { if ((yyvsp[(2) - (3)]).val==1) (yyval)=vecteur(1,symbolic(at_nodisp,(yyvsp[(1) - (3)]))); else (yyval)=vecteur(1,(yyvsp[(1) - (3)])); }
    break;

  case 5:
#line 183 "input_parser.yy"
    { if ((yyvsp[(2) - (3)]).val==1) (yyval)=mergevecteur(makevecteur(symbolic(at_nodisp,(yyvsp[(1) - (3)]))),*(yyvsp[(3) - (3)])._VECTptr); else (yyval)=mergevecteur(makevecteur((yyvsp[(1) - (3)])),*(yyvsp[(3) - (3)])._VECTptr); }
    break;

  case 6:
#line 186 "input_parser.yy"
    {(yyval) = (yyvsp[(1) - (1)]);}
    break;

  case 7:
#line 187 "input_parser.yy"
    {(yyval)=symbolic(at_prod,makevecteur((yyvsp[(1) - (2)]),(yyvsp[(2) - (2)])));}
    break;

  case 8:
#line 188 "input_parser.yy"
    {(yyval)=symbolic(at_prod,makevecteur((yyvsp[(1) - (4)]),symb_pow((yyvsp[(2) - (4)]),(yyvsp[(4) - (4)]))));}
    break;

  case 9:
#line 189 "input_parser.yy"
    {(yyval)=symbolic(at_prod,makevecteur((yyvsp[(1) - (3)]),symb_pow((yyvsp[(2) - (3)]),(yyvsp[(3) - (3)]))));}
    break;

  case 10:
#line 190 "input_parser.yy"
    {(yyval)=symbolic(at_prod,makevecteur((yyvsp[(1) - (2)]),(yyvsp[(2) - (2)])));}
    break;

  case 11:
#line 191 "input_parser.yy"
    { (yyval)=(yyvsp[(1) - (1)]); }
    break;

  case 12:
#line 192 "input_parser.yy"
    { if ((yyvsp[(1) - (1)]).type==_FUNC) (yyval)=symbolic(*(yyvsp[(1) - (1)])._FUNCptr,gen(vecteur(0),_SEQ__VECT)); else (yyval)=(yyvsp[(1) - (1)]); }
    break;

  case 13:
#line 195 "input_parser.yy"
    {(yyval) = symb_program_sto((yyvsp[(3) - (6)]),(yyvsp[(3) - (6)])*zero,(yyvsp[(6) - (6)]),(yyvsp[(1) - (6)]),false,giac_yyget_extra(scanner));}
    break;

  case 14:
#line 196 "input_parser.yy"
    {(yyval) = symb_program_sto((yyvsp[(3) - (6)]),(yyvsp[(3) - (6)])*zero,(yyvsp[(6) - (6)]),(yyvsp[(1) - (6)]),true,giac_yyget_extra(scanner));}
    break;

  case 15:
#line 197 "input_parser.yy"
    {(yyval) = symb_program_sto((yyvsp[(5) - (6)]),(yyvsp[(5) - (6)])*zero,(yyvsp[(1) - (6)]),(yyvsp[(3) - (6)]),false,giac_yyget_extra(scanner));}
    break;

  case 16:
#line 198 "input_parser.yy"
    { 
         const giac::context * contextptr = giac_yyget_extra(scanner);
         gen g=symb_at((yyvsp[(3) - (6)]),(yyvsp[(5) - (6)]),contextptr); (yyval)=symb_sto((yyvsp[(1) - (6)]),g); 
        }
    break;

  case 17:
#line 202 "input_parser.yy"
    { (yyval)=symb_sto((yyvsp[(1) - (3)]),(yyvsp[(3) - (3)])); }
    break;

  case 18:
#line 203 "input_parser.yy"
    { 
          const giac::context * contextptr = giac_yyget_extra(scanner);
          (yyval)=symb_program_sto((yyvsp[(4) - (13)]),(yyvsp[(4) - (13)])*zero,symb_local((yyvsp[(10) - (13)]),mergevecteur(*(yyvsp[(7) - (13)])._VECTptr,*(yyvsp[(12) - (13)])._VECTptr),contextptr),(yyvsp[(2) - (13)]),false,contextptr); 
	}
    break;

  case 19:
#line 207 "input_parser.yy"
    { 
          const giac::context * contextptr = giac_yyget_extra(scanner);
	(yyval)=symb_program_sto((yyvsp[(4) - (12)]),(yyvsp[(4) - (12)])*zero,symb_local((yyvsp[(9) - (12)]),mergevecteur(*(yyvsp[(7) - (12)])._VECTptr,*(yyvsp[(11) - (12)])._VECTptr),contextptr),(yyvsp[(2) - (12)]),false,contextptr); 
	}
    break;

  case 20:
#line 211 "input_parser.yy"
    { 
          const giac::context * contextptr = giac_yyget_extra(scanner);
	(yyval)=symb_program_sto((yyvsp[(4) - (12)]),(yyvsp[(4) - (12)])*zero,symb_local((yyvsp[(9) - (12)]),(yyvsp[(11) - (12)]),contextptr),(yyvsp[(2) - (12)]),false,contextptr); 
	}
    break;

  case 21:
#line 215 "input_parser.yy"
    { 
	(yyval)=symb_program_sto((yyvsp[(4) - (8)]),(yyvsp[(4) - (8)])*zero,symb_bloc((yyvsp[(7) - (8)])),(yyvsp[(2) - (8)]),false,giac_yyget_extra(scanner)); 
	}
    break;

  case 22:
#line 218 "input_parser.yy"
    { (yyval)=symbolic(*(yyvsp[(1) - (3)])._FUNCptr,(yyvsp[(2) - (3)])); }
    break;

  case 23:
#line 219 "input_parser.yy"
    { (yyval)=symbolic(*(yyvsp[(1) - (2)])._FUNCptr,(yyvsp[(2) - (2)])); }
    break;

  case 24:
#line 220 "input_parser.yy"
    { (yyval)=(yyvsp[(2) - (2)]); }
    break;

  case 25:
#line 221 "input_parser.yy"
    { (yyval)=symb_program_sto((yyvsp[(4) - (7)]),(yyvsp[(4) - (7)])*zero,(yyvsp[(7) - (7)]),(yyvsp[(2) - (7)]),false,giac_yyget_extra(scanner));}
    break;

  case 26:
#line 222 "input_parser.yy"
    { 
          const giac::context * contextptr = giac_yyget_extra(scanner);
          (yyval)=symb_program_sto((yyvsp[(4) - (13)]),(yyvsp[(4) - (13)])*zero,symb_local((yyvsp[(10) - (13)]),(yyvsp[(12) - (13)]),contextptr),(yyvsp[(2) - (13)]),false,contextptr);
        }
    break;

  case 27:
#line 226 "input_parser.yy"
    { (yyval)=symb_program_sto((yyvsp[(4) - (9)]),(yyvsp[(4) - (9)])*zero,symb_bloc((yyvsp[(8) - (9)])),(yyvsp[(2) - (9)]),false,giac_yyget_extra(scanner)); }
    break;

  case 28:
#line 227 "input_parser.yy"
    {(yyval) = symb_of((yyvsp[(1) - (4)]),(yyvsp[(3) - (4)]));}
    break;

  case 29:
#line 228 "input_parser.yy"
    {(yyval) = symb_of((yyvsp[(1) - (4)]),(yyvsp[(3) - (4)]));}
    break;

  case 30:
#line 229 "input_parser.yy"
    { (yyval)=(yyvsp[(1) - (1)]); }
    break;

  case 31:
#line 230 "input_parser.yy"
    { (yyval)=symbolic(at_maple_lib,makevecteur((yyvsp[(1) - (4)]),(yyvsp[(3) - (4)]))); }
    break;

  case 32:
#line 231 "input_parser.yy"
    {(yyval) = (yyvsp[(1) - (1)]);}
    break;

  case 33:
#line 232 "input_parser.yy"
    {(yyval) = (yyvsp[(1) - (1)]);}
    break;

  case 34:
#line 233 "input_parser.yy"
    {(yyval) = (yyvsp[(1) - (1)]);}
    break;

  case 35:
#line 234 "input_parser.yy"
    {(yyval) = symbolic(*(yyvsp[(1) - (3)])._FUNCptr,(yyvsp[(3) - (3)]));}
    break;

  case 36:
#line 235 "input_parser.yy"
    {(yyval) = symbolic(*(yyvsp[(2) - (3)])._FUNCptr,gen(makevecteur((yyvsp[(1) - (3)]),(yyvsp[(3) - (3)])),_SEQ__VECT));}
    break;

  case 37:
#line 236 "input_parser.yy"
    {(yyval) = equal((yyvsp[(1) - (3)]),(yyvsp[(3) - (3)])); }
    break;

  case 38:
#line 237 "input_parser.yy"
    { 
	if ((yyvsp[(2) - (2)]).type==_SYMB) (yyval)=(yyvsp[(2) - (2)]); else (yyval)=symbolic(at_nop,(yyvsp[(2) - (2)])); 
	(yyval).change_subtype(_SPREAD__SYMB); 
        const giac::context * contextptr = giac_yyget_extra(scanner);
        spread_formula(false,contextptr); 
	}
    break;

  case 39:
#line 243 "input_parser.yy"
    { /* if ($2==at_plus) $$=symb_plus($1,$3); else */
  (yyval) =symbolic(*(yyvsp[(2) - (3)])._FUNCptr,makevecteur((yyvsp[(1) - (3)]),(yyvsp[(3) - (3)])));}
    break;

  case 40:
#line 245 "input_parser.yy"
    {(yyval) = symb_plus((yyvsp[(1) - (3)]),symbolic(at_neg,(yyvsp[(3) - (3)])));}
    break;

  case 41:
#line 246 "input_parser.yy"
    {(yyval) =symbolic(*(yyvsp[(2) - (3)])._FUNCptr,gen(makevecteur((yyvsp[(1) - (3)]),(yyvsp[(3) - (3)])),_SEQ__VECT));}
    break;

  case 42:
#line 247 "input_parser.yy"
    {(yyval) =symbolic(*(yyvsp[(2) - (3)])._FUNCptr,makevecteur((yyvsp[(1) - (3)]),(yyvsp[(3) - (3)])));}
    break;

  case 43:
#line 248 "input_parser.yy"
    {(yyval) =symbolic(*(yyvsp[(2) - (3)])._FUNCptr,makevecteur((yyvsp[(1) - (3)]),(yyvsp[(3) - (3)])));}
    break;

  case 44:
#line 249 "input_parser.yy"
    {if ((yyvsp[(2) - (3)]).type==_FUNC) (yyval)=symbolic(*(yyvsp[(2) - (3)])._FUNCptr,gen(makevecteur((yyvsp[(1) - (3)]),(yyvsp[(3) - (3)])),_SEQ__VECT)); else (yyval) = symbolic(at_normalmod,gen(makevecteur((yyvsp[(1) - (3)]),(yyvsp[(3) - (3)])),_SEQ__VECT));}
    break;

  case 45:
#line 250 "input_parser.yy"
    {(yyval) = symbolic(*(yyvsp[(2) - (3)])._FUNCptr,makevecteur((yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]))); }
    break;

  case 46:
#line 253 "input_parser.yy"
    {(yyval) = symbolic(*(yyvsp[(2) - (3)])._FUNCptr,gen(makevecteur((yyvsp[(1) - (3)]),(yyvsp[(3) - (3)])),_SEQ__VECT));}
    break;

  case 47:
#line 254 "input_parser.yy"
    {(yyval)= symbolic(at_deuxpoints,makevecteur((yyvsp[(1) - (3)]),(yyvsp[(3) - (3)])));}
    break;

  case 48:
#line 255 "input_parser.yy"
    { 
					if ((yyvsp[(2) - (2)])==unsigned_inf)
						(yyval) = minus_inf;
					else { if ((yyvsp[(2) - (2)]).type==_INT_) (yyval)=(-(yyvsp[(2) - (2)]).val); else (yyval)=symbolic(at_neg,(yyvsp[(2) - (2)])); }
				}
    break;

  case 49:
#line 260 "input_parser.yy"
    {
					if ((yyvsp[(2) - (2)])==unsigned_inf)
						(yyval) = plus_inf;
					else
						(yyval) = symbolic(at_plus,(yyvsp[(2) - (2)]));
				}
    break;

  case 50:
#line 266 "input_parser.yy"
    {(yyval) = polynome_or_sparse_poly1((yyvsp[(2) - (5)]),(yyvsp[(4) - (5)]));}
    break;

  case 51:
#line 267 "input_parser.yy"
    { 
           if ( ((yyvsp[(2) - (3)]).type==_SYMB) && ((yyvsp[(2) - (3)])._SYMBptr->sommet==at_deuxpoints) )
             (yyval) = algebraic_EXTension((yyvsp[(2) - (3)])._SYMBptr->feuille._VECTptr->front(),(yyvsp[(2) - (3)])._SYMBptr->feuille._VECTptr->back());
           else (yyval)=(yyvsp[(2) - (3)]);
        }
    break;

  case 52:
#line 272 "input_parser.yy"
    {if ((yyvsp[(2) - (5)]).type==_VECT)
     (yyval) = real_complex_rootof(*(yyvsp[(2) - (5)])._VECTptr,(yyvsp[(4) - (5)])); else (yyval)=zero;}
    break;

  case 53:
#line 274 "input_parser.yy"
    { (yyval)=gen(at_of,2); }
    break;

  case 54:
#line 275 "input_parser.yy"
    {(yyval) = symb_sto((yyvsp[(3) - (3)]),(yyvsp[(1) - (3)]),(yyvsp[(2) - (3)])==at_array_sto); if ((yyvsp[(3) - (3)]).is_symb_of_sommet(at_program)) *logptr(giac_yyget_extra(scanner))<<"// End defining "<<(yyvsp[(1) - (3)])<<endl;}
    break;

  case 55:
#line 276 "input_parser.yy"
    { (yyval) = symbolic(*(yyvsp[(1) - (2)])._FUNCptr,(yyvsp[(2) - (2)]));}
    break;

  case 56:
#line 277 "input_parser.yy"
    {(yyval) = symb_args((yyvsp[(3) - (4)]));}
    break;

  case 57:
#line 278 "input_parser.yy"
    {(yyval) = symb_args((yyvsp[(3) - (4)]));}
    break;

  case 58:
#line 279 "input_parser.yy"
    { (yyval)=symb_args(vecteur(0)); }
    break;

  case 59:
#line 280 "input_parser.yy"
    {
	(yyval) = symbolic(*(yyvsp[(1) - (4)])._FUNCptr,(yyvsp[(3) - (4)]));
        const giac::context * contextptr = giac_yyget_extra(scanner);
	if (*(yyvsp[(1) - (4)])._FUNCptr==at_maple_mode ||*(yyvsp[(1) - (4)])._FUNCptr==at_xcas_mode ){
          xcas_mode(contextptr)=(yyvsp[(3) - (4)]).val;
        }
	if (*(yyvsp[(1) - (4)])._FUNCptr==at_user_operator){
          user_operator((yyvsp[(3) - (4)]),contextptr);
        }
	}
    break;

  case 60:
#line 290 "input_parser.yy"
    { 
          const giac::context * contextptr = giac_yyget_extra(scanner);
          (yyval)=symb_at((yyvsp[(1) - (4)]),(yyvsp[(3) - (4)]),contextptr);
        }
    break;

  case 61:
#line 294 "input_parser.yy"
    {
	(yyval) = symbolic(*(yyvsp[(1) - (3)])._FUNCptr,gen(vecteur(0),_SEQ__VECT));
	if (*(yyvsp[(1) - (3)])._FUNCptr==at_rpn)
          rpn_mode=1;
	if (*(yyvsp[(1) - (3)])._FUNCptr==at_alg)
          rpn_mode=0;
	}
    break;

  case 62:
#line 301 "input_parser.yy"
    {
	(yyval) = (yyvsp[(1) - (1)]);
	}
    break;

  case 63:
#line 304 "input_parser.yy"
    {(yyval) = symbolic(at_derive,(yyvsp[(1) - (2)]));}
    break;

  case 64:
#line 305 "input_parser.yy"
    { (yyval)=symbolic(*(yyvsp[(2) - (2)])._FUNCptr,(yyvsp[(1) - (2)])); }
    break;

  case 65:
#line 306 "input_parser.yy"
    { (yyval)=symb_rpn_prog((yyvsp[(2) - (3)])); }
    break;

  case 66:
#line 307 "input_parser.yy"
    { 
          const giac::context * contextptr = giac_yyget_extra(scanner);
           (yyval)=symb_program((yyvsp[(3) - (7)]),zero*(yyvsp[(3) - (7)]),symb_local((yyvsp[(5) - (7)]),(yyvsp[(6) - (7)]),contextptr),contextptr); 
        }
    break;

  case 67:
#line 311 "input_parser.yy"
    { 
          const giac::context * contextptr = giac_yyget_extra(scanner);
         (yyval)=symb_program((yyvsp[(3) - (8)]),zero*(yyvsp[(3) - (8)]),symb_local((yyvsp[(5) - (8)]),(yyvsp[(7) - (8)]),contextptr),contextptr); 
        }
    break;

  case 68:
#line 315 "input_parser.yy"
    {
	(yyval) = symbolic(*(yyvsp[(1) - (6)])._FUNCptr,makevecteur(equaltosame((yyvsp[(3) - (6)])),symb_bloc((yyvsp[(5) - (6)])),(yyvsp[(6) - (6)])));
	}
    break;

  case 69:
#line 318 "input_parser.yy"
    {
        vecteur v=makevecteur(equaltosame((yyvsp[(3) - (7)])),(yyvsp[(5) - (7)]),(yyvsp[(7) - (7)]));
	// *logptr(giac_yyget_extra(scanner)) << v << endl;
	(yyval) = symbolic(*(yyvsp[(1) - (7)])._FUNCptr,v);
	}
    break;

  case 70:
#line 330 "input_parser.yy"
    {(yyval) = symbolic(*(yyvsp[(1) - (6)])._FUNCptr,makevecteur(equaltosame((yyvsp[(2) - (6)])),symb_bloc((yyvsp[(4) - (6)])),symb_bloc((yyvsp[(6) - (6)]))));}
    break;

  case 71:
#line 331 "input_parser.yy"
    {(yyval) = symbolic(*(yyvsp[(1) - (4)])._FUNCptr,makevecteur(equaltosame((yyvsp[(2) - (4)])),(yyvsp[(4) - (4)]),undef));}
    break;

  case 72:
#line 332 "input_parser.yy"
    {(yyval) = symbolic(*(yyvsp[(1) - (4)])._FUNCptr,makevecteur(equaltosame((yyvsp[(2) - (4)])),(yyvsp[(4) - (4)]),undef));}
    break;

  case 73:
#line 333 "input_parser.yy"
    {
	(yyval) = symbolic(*(yyvsp[(1) - (5)])._FUNCptr,makevecteur(equaltosame((yyvsp[(2) - (5)])),symb_bloc((yyvsp[(4) - (5)])),(yyvsp[(5) - (5)])));
	}
    break;

  case 74:
#line 336 "input_parser.yy"
    {(yyval) = symbolic(*(yyvsp[(1) - (4)])._FUNCptr,(yyvsp[(3) - (4)]));}
    break;

  case 75:
#line 337 "input_parser.yy"
    {(yyval) = (yyvsp[(1) - (1)]);}
    break;

  case 76:
#line 338 "input_parser.yy"
    {(yyval) = symb_program((yyvsp[(3) - (4)]));}
    break;

  case 77:
#line 339 "input_parser.yy"
    {(yyval) = gen(at_program,3);}
    break;

  case 78:
#line 340 "input_parser.yy"
    {
          const giac::context * contextptr = giac_yyget_extra(scanner);
         (yyval) = symb_program((yyvsp[(1) - (3)]),zero*(yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]),contextptr);
        }
    break;

  case 79:
#line 344 "input_parser.yy"
    {
          const giac::context * contextptr = giac_yyget_extra(scanner);
             if ((yyvsp[(3) - (3)]).type==_VECT) 
                (yyval) = symb_program((yyvsp[(1) - (3)]),zero*(yyvsp[(1) - (3)]),symb_bloc(makevecteur(at_nop,(yyvsp[(3) - (3)]))),contextptr); 
             else 
                (yyval) = symb_program((yyvsp[(1) - (3)]),zero*(yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]),contextptr);
		}
    break;

  case 80:
#line 351 "input_parser.yy"
    {(yyval) = symb_bloc((yyvsp[(3) - (4)]));}
    break;

  case 81:
#line 352 "input_parser.yy"
    {(yyval) = at_bloc;}
    break;

  case 82:
#line 353 "input_parser.yy"
    {
          const giac::context * contextptr = giac_yyget_extra(scanner);
          (yyval) = symb_local((yyvsp[(3) - (4)]),contextptr);
        }
    break;

  case 83:
#line 357 "input_parser.yy"
    {(yyval) = gen(at_local,2);}
    break;

  case 84:
#line 359 "input_parser.yy"
    {(yyval) = gen(*(yyvsp[(1) - (2)])._FUNCptr,0);}
    break;

  case 85:
#line 360 "input_parser.yy"
    { (yyval)=symbolic(*(yyvsp[(1) - (2)])._FUNCptr,(yyvsp[(2) - (2)])); }
    break;

  case 86:
#line 361 "input_parser.yy"
    {(yyval) = gen(*(yyvsp[(1) - (1)])._FUNCptr,0);}
    break;

  case 87:
#line 362 "input_parser.yy"
    { (yyval)=(yyvsp[(2) - (3)]);}
    break;

  case 88:
#line 363 "input_parser.yy"
    { (yyval)=symbolic(*(yyvsp[(1) - (2)])._FUNCptr,(yyvsp[(2) - (2)])); }
    break;

  case 89:
#line 364 "input_parser.yy"
    {(yyval) = symbolic(*(yyvsp[(1) - (1)])._FUNCptr,gen(vecteur(0),_SEQ__VECT));}
    break;

  case 90:
#line 365 "input_parser.yy"
    {(yyval) = symbolic(*(yyvsp[(1) - (3)])._FUNCptr,gen(vecteur(0),_SEQ__VECT));}
    break;

  case 91:
#line 367 "input_parser.yy"
    {(yyval) = symbolic(at_break,zero);}
    break;

  case 92:
#line 368 "input_parser.yy"
    {(yyval) = symbolic(at_continue,zero);}
    break;

  case 93:
#line 369 "input_parser.yy"
    {(yyval) = symbolic(*(yyvsp[(1) - (9)])._FUNCptr,makevecteur((yyvsp[(3) - (9)]),equaltosame((yyvsp[(5) - (9)])),(yyvsp[(7) - (9)]),symb_bloc((yyvsp[(9) - (9)]))));}
    break;

  case 94:
#line 370 "input_parser.yy"
    {(yyval) = symbolic(*(yyvsp[(1) - (10)])._FUNCptr,makevecteur((yyvsp[(3) - (10)]),equaltosame((yyvsp[(5) - (10)])),(yyvsp[(7) - (10)]),(yyvsp[(9) - (10)])));}
    break;

  case 95:
#line 371 "input_parser.yy"
    {(yyval) = symbolic(*(yyvsp[(1) - (4)])._FUNCptr,makevecteur((yyvsp[(3) - (4)])));}
    break;

  case 96:
#line 372 "input_parser.yy"
    { 
	  gen kk(identificateur("index"));
	  vecteur v(*(yyvsp[(6) - (7)])._VECTptr);
          const giac::context * contextptr = giac_yyget_extra(scanner);
	  v.insert(v.begin(),symb_sto(symb_at((yyvsp[(4) - (7)]),kk,contextptr),(yyvsp[(2) - (7)])));
	  (yyval)=symbolic(*(yyvsp[(1) - (7)])._FUNCptr,makevecteur(symb_sto(xcas_mode(contextptr)!=0,kk),symb_inferieur_strict(kk,symb_size((yyvsp[(4) - (7)]))+(xcas_mode(contextptr)!=0)),symb_sto(symb_plus(kk,plus_one),kk),symb_bloc(v))); 
	  }
    break;

  case 97:
#line 379 "input_parser.yy"
    { 
          gen tmp; 
          const giac::context * contextptr = giac_yyget_extra(scanner);
          if (!has_evalf((yyvsp[(6) - (9)]),tmp,1,contextptr) || is_positive(tmp,contextptr)) 
             (yyval)=symbolic(*(yyvsp[(1) - (9)])._FUNCptr,makevecteur(symb_sto((yyvsp[(3) - (9)]),(yyvsp[(2) - (9)])),symb_inferieur_egal((yyvsp[(2) - (9)]),(yyvsp[(5) - (9)])),symb_sto(symb_plus((yyvsp[(2) - (9)]),(yyvsp[(6) - (9)])),(yyvsp[(2) - (9)])),symb_bloc((yyvsp[(8) - (9)])))); 
          else 
            (yyval)=symbolic(*(yyvsp[(1) - (9)])._FUNCptr,makevecteur(symb_sto((yyvsp[(3) - (9)]),(yyvsp[(2) - (9)])),symb_superieur_egal((yyvsp[(2) - (9)]),(yyvsp[(5) - (9)])),symb_sto(symb_plus((yyvsp[(2) - (9)]),(yyvsp[(6) - (9)])),(yyvsp[(2) - (9)])),symb_bloc((yyvsp[(8) - (9)])))); 
        }
    break;

  case 98:
#line 387 "input_parser.yy"
    { 
         gen tmp; 
          const giac::context * contextptr = giac_yyget_extra(scanner);
         if (!has_evalf((yyvsp[(6) - (9)]),tmp,1,contextptr) || is_positive(tmp,contextptr)) 
           (yyval)=symbolic(*(yyvsp[(1) - (9)])._FUNCptr,makevecteur(symb_sto((yyvsp[(3) - (9)]),(yyvsp[(2) - (9)])),symb_inferieur_egal((yyvsp[(2) - (9)]),(yyvsp[(6) - (9)])),symb_sto(symb_plus((yyvsp[(2) - (9)]),(yyvsp[(4) - (9)])),(yyvsp[(2) - (9)])),symb_bloc((yyvsp[(8) - (9)])))); 
         else 
           (yyval)=symbolic(*(yyvsp[(1) - (9)])._FUNCptr,makevecteur(symb_sto((yyvsp[(3) - (9)]),(yyvsp[(2) - (9)])),symb_superieur_egal((yyvsp[(2) - (9)]),(yyvsp[(6) - (9)])),symb_sto(symb_plus((yyvsp[(2) - (9)]),(yyvsp[(4) - (9)])),(yyvsp[(2) - (9)])),symb_bloc((yyvsp[(8) - (9)])))); 
        }
    break;

  case 99:
#line 395 "input_parser.yy"
    { (yyval)=symbolic(*(yyvsp[(1) - (7)])._FUNCptr,makevecteur(symb_sto((yyvsp[(3) - (7)]),(yyvsp[(2) - (7)])),plus_one,symb_sto(symb_plus((yyvsp[(2) - (7)]),(yyvsp[(4) - (7)])),(yyvsp[(2) - (7)])),symb_bloc((yyvsp[(6) - (7)])))); }
    break;

  case 100:
#line 396 "input_parser.yy"
    { (yyval)=symbolic(*(yyvsp[(1) - (9)])._FUNCptr,makevecteur(symb_sto((yyvsp[(3) - (9)]),(yyvsp[(2) - (9)])),(yyvsp[(6) - (9)]),symb_sto(symb_plus((yyvsp[(2) - (9)]),(yyvsp[(4) - (9)])),(yyvsp[(2) - (9)])),symb_bloc((yyvsp[(8) - (9)])))); }
    break;

  case 101:
#line 397 "input_parser.yy"
    {(yyval) = gen(*(yyvsp[(1) - (1)])._FUNCptr,4);}
    break;

  case 102:
#line 398 "input_parser.yy"
    {
           vecteur & v=*(yyvsp[(2) - (5)])._VECTptr;
           if ( (v.size()<3) || v[0].type!=_IDNT){
             *logptr(giac_yyget_extra(scanner)) << "Syntax For name,begin,end[,step]" << endl;
             (yyval)=undef;
           }
           else {
             gen pas(plus_one);
             if (v.size()==4)
               pas=v[3];
             gen condition;
             if (is_positive(pas,0))
               condition=symb_inferieur_egal(v[0],v[2]);
             else
               symb_superieur_egal(v[0],v[2]);
             vecteur w=makevecteur(symb_sto(v[1],v[0]),condition,symb_sto(symb_plus(v[0],pas),v[0]),symb_bloc((yyvsp[(4) - (5)])));
             (yyval)=symbolic(*(yyvsp[(1) - (5)])._FUNCptr,w);
           }
	}
    break;

  case 103:
#line 417 "input_parser.yy"
    { 
	vecteur v=makevecteur(zero,equaltosame((yyvsp[(2) - (5)])),zero,symb_bloc((yyvsp[(4) - (5)])));
	(yyval)=symbolic(*(yyvsp[(1) - (5)])._FUNCptr,v); 
	}
    break;

  case 104:
#line 421 "input_parser.yy"
    { 
	vecteur v=makevecteur(zero,equaltosame((yyvsp[(3) - (5)])),zero,symb_bloc((yyvsp[(5) - (5)])));
	(yyval)=symbolic(*(yyvsp[(1) - (5)])._FUNCptr,v); 
	}
    break;

  case 105:
#line 425 "input_parser.yy"
    { 
        vecteur v=gen2vecteur((yyvsp[(2) - (4)]));
        v.push_back(symb_ifte(equaltosame((yyvsp[(4) - (4)])),symbolic(at_break,zero),undef));
	(yyval)=symbolic(*(yyvsp[(1) - (4)])._FUNCptr,makevecteur(zero,1,zero,symb_bloc(v))); 
	}
    break;

  case 106:
#line 430 "input_parser.yy"
    { vecteur v=makevecteur(zero,plus_one,zero,symb_bloc((yyvsp[(2) - (3)]))); (yyval)=symbolic(*(yyvsp[(1) - (3)])._FUNCptr,v); }
    break;

  case 107:
#line 431 "input_parser.yy"
    { 
	(yyval)=symbolic(*(yyvsp[(1) - (6)])._FUNCptr,makevecteur(zero,equaltosame((yyvsp[(3) - (6)])),zero,(yyvsp[(5) - (6)]))); 
	}
    break;

  case 108:
#line 434 "input_parser.yy"
    { (yyval)=symbolic(*(yyvsp[(1) - (5)])._FUNCptr,makevecteur(zero,equaltosame((yyvsp[(2) - (5)])),zero,symb_bloc((yyvsp[(4) - (5)])))); }
    break;

  case 109:
#line 435 "input_parser.yy"
    { (yyval)=symbolic(*(yyvsp[(1) - (3)])._FUNCptr,makevecteur(zero,plus_one,zero,symb_bloc((yyvsp[(2) - (3)])))); }
    break;

  case 110:
#line 436 "input_parser.yy"
    { (yyval)=symb_try_catch(makevecteur(symb_bloc((yyvsp[(2) - (7)])),(yyvsp[(5) - (7)]),symb_bloc((yyvsp[(7) - (7)]))));}
    break;

  case 111:
#line 437 "input_parser.yy"
    {(yyval)=symb_try_catch(makevecteur((yyvsp[(3) - (4)])));}
    break;

  case 112:
#line 438 "input_parser.yy"
    {(yyval)=gen(at_try_catch,3);}
    break;

  case 113:
#line 439 "input_parser.yy"
    { (yyval)=symb_try_catch(makevecteur(symb_bloc((yyvsp[(2) - (5)])),_IDNT_break,symb_bloc((yyvsp[(4) - (5)])))); }
    break;

  case 114:
#line 440 "input_parser.yy"
    { (yyval)=symb_try_catch(makevecteur(symb_bloc((yyvsp[(2) - (4)])),_IDNT_break,undef)); }
    break;

  case 115:
#line 441 "input_parser.yy"
    { (yyval)=symb_try_catch(makevecteur(symb_bloc((yyvsp[(2) - (6)])),_IDNT_break,symb_bloc((yyvsp[(5) - (6)])))); }
    break;

  case 116:
#line 442 "input_parser.yy"
    { (yyval)=symb_try_catch(makevecteur(symb_bloc((yyvsp[(2) - (5)])),_IDNT_break,undef)); }
    break;

  case 117:
#line 443 "input_parser.yy"
    { (yyval)=symb_case((yyvsp[(3) - (7)]),(yyvsp[(6) - (7)])); }
    break;

  case 118:
#line 444 "input_parser.yy"
    { (yyval) = symb_case((yyvsp[(3) - (4)])); }
    break;

  case 119:
#line 445 "input_parser.yy"
    { (yyval)=symb_case((yyvsp[(2) - (4)]),(yyvsp[(3) - (4)])); }
    break;

  case 120:
#line 446 "input_parser.yy"
    { 
	(yyval)=(yyvsp[(1) - (1)]); 
	// $$.subtype=1; 
	}
    break;

  case 121:
#line 450 "input_parser.yy"
    { (yyval)=(yyvsp[(2) - (3)]); /* $$.subtype=1; */ }
    break;

  case 122:
#line 451 "input_parser.yy"
    { (yyval) = symb_dollar((yyvsp[(2) - (2)])); }
    break;

  case 123:
#line 452 "input_parser.yy"
    {(yyval)=symb_dollar(makevecteur((yyvsp[(1) - (5)]),(yyvsp[(3) - (5)]),(yyvsp[(5) - (5)])));}
    break;

  case 124:
#line 453 "input_parser.yy"
    { (yyval) = symb_dollar(makevecteur((yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]))); }
    break;

  case 125:
#line 454 "input_parser.yy"
    { (yyval) = symb_dollar(makevecteur((yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]))); }
    break;

  case 126:
#line 455 "input_parser.yy"
    { (yyval)=symb_dollar((yyvsp[(2) - (2)])); }
    break;

  case 127:
#line 456 "input_parser.yy"
    { (yyval) = symb_compose(makevecteur((yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]))); }
    break;

  case 128:
#line 457 "input_parser.yy"
    { (yyval) = symb_union(makevecteur((yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]))); }
    break;

  case 129:
#line 458 "input_parser.yy"
    { (yyval) = symb_intersect(makevecteur((yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]))); }
    break;

  case 130:
#line 459 "input_parser.yy"
    { (yyval) = symb_minus(makevecteur((yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]))); }
    break;

  case 131:
#line 460 "input_parser.yy"
    { 
	(yyval)=symbolic(*(yyvsp[(2) - (3)])._FUNCptr,makevecteur((yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]))); 
	}
    break;

  case 132:
#line 463 "input_parser.yy"
    { (yyval) = (yyvsp[(1) - (1)]); }
    break;

  case 133:
#line 464 "input_parser.yy"
    {if ((yyvsp[(2) - (3)]).type==_FUNC) (yyval)=(yyvsp[(2) - (3)]); else { 
          // const giac::context * contextptr = giac_yyget_extra(scanner);
          (yyval)=symb_quote((yyvsp[(2) - (3)]));
          } 
        }
    break;

  case 134:
#line 469 "input_parser.yy"
    { (yyval)=(yyvsp[(2) - (3)]); }
    break;

  case 135:
#line 470 "input_parser.yy"
    {(yyval) = symb_of((yyvsp[(2) - (6)]),(yyvsp[(5) - (6)]));}
    break;

  case 136:
#line 471 "input_parser.yy"
    {
          const giac::context * contextptr = giac_yyget_extra(scanner);
	  (yyval) = symb_at((yyvsp[(1) - (4)]),(yyvsp[(3) - (4)]),contextptr);
        }
    break;

  case 137:
#line 475 "input_parser.yy"
    {(yyval) = (yyvsp[(2) - (3)]);}
    break;

  case 138:
#line 476 "input_parser.yy"
    { (yyval) = gen(*((yyvsp[(2) - (3)])._VECTptr),(yyvsp[(1) - (3)]).val);
	if ((yyvsp[(2) - (3)])._VECTptr->size()==1 && (yyvsp[(2) - (3)])._VECTptr->front().is_symb_of_sommet(at_ti_semi) ) {
	(yyval)=(yyvsp[(2) - (3)])._VECTptr->front();
  }
}
    break;

  case 139:
#line 481 "input_parser.yy"
    {(yyval) = gen(*((yyvsp[(2) - (3)])._VECTptr));}
    break;

  case 140:
#line 482 "input_parser.yy"
    {(yyval) = gen(*((yyvsp[(2) - (3)])._VECTptr),_SET__VECT);}
    break;

  case 141:
#line 483 "input_parser.yy"
    {(yyval) = gen(*((yyvsp[(2) - (3)])._VECTptr),_GROUP__VECT);}
    break;

  case 142:
#line 484 "input_parser.yy"
    {(yyval) = gen(*((yyvsp[(2) - (3)])._VECTptr),_LINE__VECT);}
    break;

  case 143:
#line 485 "input_parser.yy"
    {(yyval) = gen(*((yyvsp[(2) - (3)])._VECTptr),_VECTOR__VECT);}
    break;

  case 144:
#line 486 "input_parser.yy"
    {(yyval) = gen(*((yyvsp[(2) - (3)])._VECTptr),_CURVE__VECT);}
    break;

  case 145:
#line 487 "input_parser.yy"
    {(yyval) = gen(*((yyvsp[(2) - (3)])._VECTptr),_POLY1__VECT);}
    break;

  case 146:
#line 488 "input_parser.yy"
    {(yyval) = gen(*((yyvsp[(2) - (3)])._VECTptr),_MATRIX__VECT);}
    break;

  case 147:
#line 489 "input_parser.yy"
    {(yyval) = gen(*((yyvsp[(2) - (3)])._VECTptr),_ASSUME__VECT);}
    break;

  case 148:
#line 490 "input_parser.yy"
    { 
         if ((yyvsp[(1) - (3)]).type==_VECT && (yyvsp[(1) - (3)]).subtype==_SEQ__VECT && !((yyvsp[(3) - (3)]).type==_VECT && (yyvsp[(2) - (3)]).subtype==_SEQ__VECT)){ (yyval)=(yyvsp[(1) - (3)]); (yyval)._VECTptr->push_back((yyvsp[(3) - (3)])); }
	 else
           (yyval) = makesuite((yyvsp[(1) - (3)]),(yyvsp[(3) - (3)])); 
        }
    break;

  case 149:
#line 495 "input_parser.yy"
    { vecteur v1(gen2vecteur((yyvsp[(1) - (3)]))),v3(gen2vecteur((yyvsp[(3) - (3)]))); (yyval)=symbolic(at_ti_semi,makevecteur(v1,v3)); }
    break;

  case 150:
#line 496 "input_parser.yy"
    { (yyval)=gen(vecteur(0),_SEQ__VECT); }
    break;

  case 151:
#line 497 "input_parser.yy"
    { (yyval)=(yyvsp[(1) - (1)]); }
    break;

  case 152:
#line 498 "input_parser.yy"
    { (yyval)=(yyvsp[(1) - (1)]); }
    break;

  case 153:
#line 499 "input_parser.yy"
    {(yyval)=symb_findhelp((yyvsp[(2) - (2)]));}
    break;

  case 154:
#line 500 "input_parser.yy"
    { (yyval)=symb_interrogation((yyvsp[(1) - (3)]),(yyvsp[(3) - (3)])); }
    break;

  case 155:
#line 501 "input_parser.yy"
    {
          const giac::context * contextptr = giac_yyget_extra(scanner);
          (yyval)=symb_unit(plus_one,(yyvsp[(2) - (2)]),contextptr); 
        }
    break;

  case 156:
#line 505 "input_parser.yy"
    {
          const giac::context * contextptr = giac_yyget_extra(scanner);
          (yyval)=symb_unit((yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]),contextptr); 
        }
    break;

  case 157:
#line 509 "input_parser.yy"
    { (yyval)=symb_pow((yyvsp[(1) - (2)]),(yyvsp[(2) - (2)])); }
    break;

  case 158:
#line 510 "input_parser.yy"
    { 
        const giac::context * contextptr = giac_yyget_extra(scanner);
	messages_to_print += parser_filename(contextptr) + parser_error(contextptr); 
	/* *logptr(giac_yyget_extra(scanner)) << messages_to_print; */
	(yyval)=undef;
        spread_formula(false,contextptr); 
	}
    break;

  case 159:
#line 519 "input_parser.yy"
    { (yyval)=(yyvsp[(1) - (1)]); }
    break;

  case 160:
#line 526 "input_parser.yy"
    { 
        if ((yyvsp[(3) - (3)]).type==_INT_ && (yyvsp[(3) - (3)]).subtype==_INT_TYPE){
	   (yyval)=symb_check_type(makevecteur((yyvsp[(3) - (3)]),(yyvsp[(1) - (3)]))); 
        }
        else
	  (yyval)=symb_double_deux_points(makevecteur((yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]))); 
	}
    break;

  case 161:
#line 533 "input_parser.yy"
    { 
	  gen tmp((yyvsp[(1) - (2)])); 
	  // tmp.subtype=1; 
	  (yyval)=symb_check_type(makevecteur(tmp,(yyvsp[(2) - (2)]))); 
	  }
    break;

  case 162:
#line 538 "input_parser.yy"
    {(yyval)=symbolic(*(yyvsp[(1) - (2)])._FUNCptr,(yyvsp[(2) - (2)])); }
    break;

  case 163:
#line 542 "input_parser.yy"
    { (yyval)=vecteur(0); }
    break;

  case 164:
#line 543 "input_parser.yy"
    { (yyval)=mergevecteur(*(yyvsp[(1) - (2)])._VECTptr,*(yyvsp[(2) - (2)])._VECTptr); }
    break;

  case 165:
#line 544 "input_parser.yy"
    { (yyval)=(yyvsp[(2) - (2)]); }
    break;

  case 166:
#line 548 "input_parser.yy"
    { if ((yyvsp[(3) - (4)]).type==_VECT) (yyval)=gen(*(yyvsp[(3) - (4)])._VECTptr,_RPN_STACK__VECT); else (yyval)=gen(vecteur(1,(yyvsp[(3) - (4)])),_RPN_STACK__VECT); }
    break;

  case 167:
#line 549 "input_parser.yy"
    { (yyval)=gen(vecteur(0),_RPN_STACK__VECT); }
    break;

  case 168:
#line 552 "input_parser.yy"
    { if (!(yyvsp[(1) - (3)]).val) (yyval)=(yyvsp[(2) - (3)]); else (yyval)=vecteur(0);}
    break;

  case 169:
#line 555 "input_parser.yy"
    { (yyval)=(yyvsp[(2) - (3)]); }
    break;

  case 170:
#line 558 "input_parser.yy"
    { (yyval)=gen(vecteur(1,(yyvsp[(1) - (1)])),_SEQ__VECT); }
    break;

  case 171:
#line 559 "input_parser.yy"
    { 
	       vecteur v=*(yyvsp[(3) - (3)])._VECTptr;
	       v.insert(v.begin(),(yyvsp[(1) - (3)]));
	       (yyval)=gen(v,_SEQ__VECT);
	     }
    break;

  case 172:
#line 566 "input_parser.yy"
    { (yyval)=(yyvsp[(1) - (1)]); }
    break;

  case 173:
#line 567 "input_parser.yy"
    { (yyval)=symb_sto((yyvsp[(3) - (3)]),(yyvsp[(1) - (3)]),(yyvsp[(2) - (3)])==at_array_sto); }
    break;

  case 174:
#line 568 "input_parser.yy"
    { (yyval)=(yyvsp[(2) - (3)]); }
    break;

  case 175:
#line 569 "input_parser.yy"
    { (yyval)=(yyvsp[(1) - (1)]); *logptr(giac_yyget_extra(scanner)) << "Error: reserved word "<< (yyvsp[(1) - (1)]) <<endl;}
    break;

  case 176:
#line 570 "input_parser.yy"
    { (yyval)=symb_double_deux_points(makevecteur((yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]))); *logptr(giac_yyget_extra(scanner)) << "Error: reserved word "<< (yyvsp[(1) - (3)]) <<endl; }
    break;

  case 177:
#line 571 "input_parser.yy"
    { 
  const giac::context * contextptr = giac_yyget_extra(scanner);
  (yyval)=string2gen("_"+(yyvsp[(1) - (1)]).print(contextptr),false); 
  if (!giac::first_error_line(contextptr)){
    giac::first_error_line(giac::lexer_line_number(contextptr),contextptr);
    giac:: error_token_name((yyvsp[(1) - (1)]).print(contextptr)+ " (reserved word)",contextptr);
  }
}
    break;

  case 178:
#line 579 "input_parser.yy"
    { 
  const giac::context * contextptr = giac_yyget_extra(scanner);
  (yyval)=string2gen("_"+(yyvsp[(1) - (1)]).print(contextptr),false);
  if (!giac::first_error_line(contextptr)){
    giac::first_error_line(giac::lexer_line_number(contextptr),contextptr);
    giac:: error_token_name((yyvsp[(1) - (1)]).print(contextptr)+ " reserved word",contextptr);
  }
}
    break;

  case 179:
#line 589 "input_parser.yy"
    { (yyval)=plus_one;}
    break;

  case 180:
#line 590 "input_parser.yy"
    { (yyval)=(yyvsp[(1) - (1)]); }
    break;

  case 181:
#line 593 "input_parser.yy"
    { (yyval)=vecteur(0,_SEQ__VECT); }
    break;

  case 182:
#line 594 "input_parser.yy"
    { (yyval)=makesuite((yyvsp[(1) - (1)])); }
    break;

  case 183:
#line 597 "input_parser.yy"
    { (yyval) = makevecteur((yyvsp[(1) - (1)])); }
    break;

  case 184:
#line 599 "input_parser.yy"
    { vecteur v(1,(yyvsp[(1) - (2)])); 
			  if ((yyvsp[(1) - (2)]).type==_VECT) v=*((yyvsp[(1) - (2)])._VECTptr); 
			  v.push_back((yyvsp[(2) - (2)])); 
			  (yyval) = v;
			}
    break;

  case 185:
#line 604 "input_parser.yy"
    { (yyval) = (yyvsp[(1) - (2)]);}
    break;

  case 186:
#line 607 "input_parser.yy"
    { (yyval)=vecteur(0); }
    break;

  case 187:
#line 608 "input_parser.yy"
    { (yyval)=mergevecteur(vecteur(1,(yyvsp[(1) - (2)])),*((yyvsp[(2) - (2)])._VECTptr));}
    break;

  case 188:
#line 609 "input_parser.yy"
    { (yyval)=mergevecteur(vecteur(1,(yyvsp[(1) - (3)])),*((yyvsp[(3) - (3)])._VECTptr));}
    break;

  case 189:
#line 612 "input_parser.yy"
    { (yyval)=(yyvsp[(1) - (1)]); }
    break;

  case 190:
#line 613 "input_parser.yy"
    { (yyval)=(yyvsp[(2) - (3)]); }
    break;

  case 191:
#line 614 "input_parser.yy"
    {(yyval) = (yyvsp[(1) - (1)]);}
    break;

  case 192:
#line 615 "input_parser.yy"
    {(yyval) = (yyvsp[(1) - (1)]);}
    break;

  case 193:
#line 616 "input_parser.yy"
    { (yyval)=(yyvsp[(1) - (1)]); }
    break;

  case 194:
#line 617 "input_parser.yy"
    { 
            const giac::context * contextptr = giac_yyget_extra(scanner);
            (yyval)=symb_unit(plus_one,(yyvsp[(2) - (2)]),contextptr);
           }
    break;

  case 195:
#line 621 "input_parser.yy"
    { 
             const giac::context * contextptr = giac_yyget_extra(scanner);
             (yyval)=symb_unit((yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]),contextptr);
           }
    break;

  case 196:
#line 625 "input_parser.yy"
    { (yyval)=(yyvsp[(2) - (3)]); }
    break;

  case 197:
#line 626 "input_parser.yy"
    { (yyval)=gen(at_plus,2); }
    break;

  case 198:
#line 627 "input_parser.yy"
    { (yyval)=gen(at_binary_minus,2); }
    break;

  case 199:
#line 628 "input_parser.yy"
    { (yyval)=gen(at_division,2); }
    break;

  case 200:
#line 629 "input_parser.yy"
    { (yyval)=gen(at_prod,2); }
    break;

  case 201:
#line 630 "input_parser.yy"
    { (yyval)=gen(at_pow,2); }
    break;

  case 202:
#line 631 "input_parser.yy"
    { (yyval)=gen(at_equal); }
    break;

  case 203:
#line 632 "input_parser.yy"
    { (yyval)=gen(*(yyvsp[(1) - (1)])._FUNCptr,2); }
    break;

  case 204:
#line 633 "input_parser.yy"
    { (yyval)=gen(at_interval,2); }
    break;

  case 205:
#line 634 "input_parser.yy"
    {(yyval) = gen(at_and,2);}
    break;

  case 206:
#line 635 "input_parser.yy"
    { (yyval)=(yyvsp[(1) - (1)]); }
    break;

  case 207:
#line 636 "input_parser.yy"
    { (yyval)=gen(at_of,2); }
    break;

  case 208:
#line 637 "input_parser.yy"
    { (yyval) = gen(at_dollar,2); }
    break;

  case 209:
#line 638 "input_parser.yy"
    { (yyval) = gen(at_compose,2); }
    break;

  case 210:
#line 639 "input_parser.yy"
    { (yyval) = gen(at_union,2); }
    break;

  case 211:
#line 640 "input_parser.yy"
    { (yyval) = gen(at_intersect,2); }
    break;

  case 212:
#line 641 "input_parser.yy"
    { (yyval) = gen(at_minus,2); }
    break;

  case 213:
#line 642 "input_parser.yy"
    { (yyval)=(yyvsp[(1) - (1)]); }
    break;

  case 214:
#line 643 "input_parser.yy"
    { (yyval)=(yyvsp[(2) - (3)]); }
    break;

  case 215:
#line 644 "input_parser.yy"
    { (yyval)=(yyvsp[(1) - (1)]); }
    break;

  case 216:
#line 645 "input_parser.yy"
    { (yyval)=gen(*(yyvsp[(2) - (3)])._VECTptr,_RPN_FUNC__VECT); }
    break;

  case 217:
#line 646 "input_parser.yy"
    {(yyval) = symb_quote((yyvsp[(2) - (3)]));}
    break;

  case 218:
#line 647 "input_parser.yy"
    {(yyval) = gen(at_IFTE,3);}
    break;

  case 219:
#line 648 "input_parser.yy"
    { (yyval)=symb_IFTE(makevecteur((yyvsp[(2) - (5)]),(yyvsp[(4) - (5)]),symb_NOP(vecteur(0)))); }
    break;

  case 220:
#line 649 "input_parser.yy"
    { (yyval)=symb_IFTE(makevecteur((yyvsp[(2) - (7)]),(yyvsp[(4) - (7)]),(yyvsp[(6) - (7)]))); }
    break;

  case 221:
#line 650 "input_parser.yy"
    { vecteur v=*(yyvsp[(2) - (3)])._VECTptr; gen step(plus_one); if (!v.empty()) { step=v.back(); v.pop_back();} (yyval)=symb_RPN_FOR(makevecteur(identificateur(" j"),step),gen(v,_RPN_FUNC__VECT)); }
    break;

  case 222:
#line 651 "input_parser.yy"
    { (yyval)=symb_RPN_FOR(makevecteur(identificateur(" j"),plus_one),(yyvsp[(2) - (3)])); }
    break;

  case 223:
#line 652 "input_parser.yy"
    { vecteur v=*(yyvsp[(3) - (4)])._VECTptr; gen step(plus_one); if (!v.empty()) { step=v.back(); v.pop_back();} (yyval)=symb_RPN_FOR(makevecteur((yyvsp[(2) - (4)]),step),gen(v,_RPN_FUNC__VECT)); }
    break;

  case 224:
#line 653 "input_parser.yy"
    { (yyval)=symb_RPN_FOR(makevecteur((yyvsp[(2) - (4)]),plus_one),(yyvsp[(3) - (4)])); }
    break;

  case 225:
#line 654 "input_parser.yy"
    { (yyval)=symb_RPN_WHILE((yyvsp[(2) - (5)]),(yyvsp[(4) - (5)]));}
    break;

  case 226:
#line 655 "input_parser.yy"
    { (yyval)=symb_RPN_UNTIL((yyvsp[(2) - (5)]),(yyvsp[(4) - (5)])); }
    break;

  case 227:
#line 656 "input_parser.yy"
    { (yyval)=symb_RPN_LOCAL((yyvsp[(2) - (3)]),(yyvsp[(3) - (3)])); }
    break;

  case 228:
#line 657 "input_parser.yy"
    { (yyval)=symb_RPN_CASE((yyvsp[(2) - (3)])); }
    break;

  case 229:
#line 658 "input_parser.yy"
    { vecteur v(*(yyvsp[(2) - (4)])._VECTptr); v.push_back((yyvsp[(3) - (4)])); (yyval)=symb_RPN_CASE(v); }
    break;

  case 230:
#line 659 "input_parser.yy"
    { (yyval)=(yyvsp[(1) - (1)]); }
    break;

  case 231:
#line 662 "input_parser.yy"
    { (yyval)=gen(*(yyvsp[(2) - (3)])._VECTptr,_RPN_FUNC__VECT); }
    break;

  case 232:
#line 663 "input_parser.yy"
    {(yyval) = symb_quote((yyvsp[(2) - (3)]));}
    break;

  case 233:
#line 666 "input_parser.yy"
    { (yyval)=vecteur(1,(yyvsp[(1) - (1)])); }
    break;

  case 234:
#line 667 "input_parser.yy"
    { vecteur v=*(yyvsp[(1) - (2)])._VECTptr; v.push_back((yyvsp[(2) - (2)])); (yyval)=v; }
    break;

  case 235:
#line 670 "input_parser.yy"
    { (yyval)=plus_one; }
    break;

  case 236:
#line 671 "input_parser.yy"
    { (yyval)=(yyvsp[(2) - (2)]); }
    break;

  case 237:
#line 674 "input_parser.yy"
    { (yyval)=plus_one; }
    break;

  case 238:
#line 675 "input_parser.yy"
    { (yyval)=(yyvsp[(2) - (2)]); }
    break;

  case 239:
#line 678 "input_parser.yy"
    { (yyval)=undef; }
    break;

  case 240:
#line 679 "input_parser.yy"
    { (yyval)=(yyvsp[(2) - (3)]); }
    break;

  case 241:
#line 680 "input_parser.yy"
    { (yyval)=symb_bloc((yyvsp[(2) - (2)])); }
    break;

  case 242:
#line 684 "input_parser.yy"
    { 
	(yyval) = (yyvsp[(2) - (3)]);
	}
    break;

  case 243:
#line 687 "input_parser.yy"
    {
          const giac::context * contextptr = giac_yyget_extra(scanner);
          (yyval) = symb_local((yyvsp[(2) - (4)]),(yyvsp[(3) - (4)]),contextptr);
         }
    break;

  case 244:
#line 694 "input_parser.yy"
    { (yyval)=undef; }
    break;

  case 245:
#line 695 "input_parser.yy"
    {
	(yyval)=symb_bloc((yyvsp[(2) - (3)])); 
	}
    break;

  case 246:
#line 698 "input_parser.yy"
    { 
	  (yyval)=symb_ifte(equaltosame((yyvsp[(2) - (5)])),symb_bloc((yyvsp[(4) - (5)])),(yyvsp[(5) - (5)]));
	  }
    break;

  case 247:
#line 701 "input_parser.yy"
    { 
	  (yyval)=symb_ifte(equaltosame((yyvsp[(3) - (6)])),symb_bloc((yyvsp[(5) - (6)])),(yyvsp[(6) - (6)]));
	  }
    break;

  case 248:
#line 706 "input_parser.yy"
    { (yyval)=undef; }
    break;

  case 249:
#line 707 "input_parser.yy"
    { (yyval)=undef; }
    break;

  case 250:
#line 710 "input_parser.yy"
    { (yyval)=undef; }
    break;

  case 251:
#line 711 "input_parser.yy"
    { (yyval)=undef; }
    break;

  case 252:
#line 714 "input_parser.yy"
    { (yyval)=vecteur(0); }
    break;

  case 253:
#line 715 "input_parser.yy"
    { 
	  vecteur v(*(yyvsp[(1) - (5)])._VECTptr); 
	  v.push_back((yyvsp[(2) - (5)])); 
	  v.push_back((yyvsp[(4) - (5)])); (yyval)=v; 
	  }
    break;

  case 254:
#line 722 "input_parser.yy"
    { (yyval)=vecteur(0); }
    break;

  case 255:
#line 723 "input_parser.yy"
    { (yyval)=makevecteur(symb_bloc((yyvsp[(3) - (3)])));}
    break;

  case 256:
#line 724 "input_parser.yy"
    { (yyval)=mergevecteur(makevecteur((yyvsp[(2) - (5)]),symb_bloc((yyvsp[(4) - (5)]))),*((yyvsp[(5) - (5)])._VECTptr));}
    break;

  case 257:
#line 727 "input_parser.yy"
    { (yyval)=vecteur(0); }
    break;

  case 258:
#line 728 "input_parser.yy"
    { (yyval)=vecteur(1,symb_bloc((yyvsp[(2) - (2)]))); }
    break;

  case 259:
#line 729 "input_parser.yy"
    { (yyval)=mergevecteur(makevecteur((yyvsp[(2) - (5)]),symb_bloc((yyvsp[(4) - (5)]))),*((yyvsp[(5) - (5)])._VECTptr));}
    break;

  case 260:
#line 732 "input_parser.yy"
    { (yyval)=(yyvsp[(1) - (1)]); }
    break;


/* Line 1267 of yacc.c.  */
#line 5838 "y.tab.c"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (scanner, YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (scanner, yymsg);
	  }
	else
	  {
	    yyerror (scanner, YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval, scanner);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse look-ahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp, scanner);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  *++yyvsp = yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (scanner, YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEOF && yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval, scanner);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp, scanner);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}


#line 739 "input_parser.yy"


#ifndef NO_NAMESPACE_GIAC
} // namespace giac


#endif // ndef NO_NAMESPACE_GIAC

// Error print routine (store error string in parser_error)
int giac_yyerror(yyscan_t scanner,const char *s)
{
  const giac::context * contextptr = giac_yyget_extra(scanner);
  if ( (*giac_yyget_text( scanner )) && (*giac_yyget_text( scanner )!='¤')){
    std::string txt=giac_yyget_text( scanner );
    parser_error( ":" + giac::print_INT_(giac::lexer_line_number(contextptr)) + ": " + string(s) + " line " + giac::print_INT_(giac::lexer_line_number(contextptr)) + " at " + txt +"\n",contextptr);
     giac::parsed_gen(giac::string2gen(txt,false),contextptr);
  }
  else {
    parser_error(":" + giac::print_INT_(giac::lexer_line_number(contextptr)) + ": " +string(s) + " at end of input\n",contextptr);
    giac::parsed_gen(giac::undef,contextptr);
  }
  if (!giac::first_error_line(contextptr)){
    giac::first_error_line(giac::lexer_line_number(contextptr),contextptr);
    std::string s=string(giac_yyget_text( scanner ));
    if (s=="¤")
      s="end of input";
    giac:: error_token_name(s,contextptr);
  }
  return giac::lexer_line_number(contextptr);
}

