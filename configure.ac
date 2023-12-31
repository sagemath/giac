dnl Run autoreconf -iv to produce the configure script and associated files.

dnl Giac version number
m4_define([giac_major_version], [1])
m4_define([giac_minor_version], [9])
m4_define([giac_micro_version], [0])
m4_define([giac_version], [giac_major_version.giac_minor_version.giac_micro_version])

dnl Attempt to change archive name made by make dist 
dnl REVNUM=`date +%Y%m%d`
dnl AC_SUBST(REVNUM)


dnl Giac library ABI versioning information 
dnl (NOT to be confused with version number) 
m4_define([giac_lt_current], [0])
m4_define([giac_lt_age], [0])
m4_define([giac_lt_revision], [0])

AC_INIT([giac], giac_version)
AC_PREREQ(2.59)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([src/index.h])
dnl XXX: gettext hard-codes location of config.h. f@#$!
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([gnu 1.8 dist-bzip2])

dnl AM_MAINTAINER_MODE

dnl AC_GNU_SOURCE
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_CXX
AC_PROG_CXXCPP
AC_LANG([C++])
AC_PROG_LIBTOOL
AC_PROG_YACC
AM_PROG_LEX
AC_C_BIGENDIAN

dnl Check for standard C+headers
GINAC_STD_CXX_HEADERS
AC_CHECK_HEADERS(regex)
AC_CHECK_FUNCS(tgammaf)

dnl CXXFLAGS="$CXXFLAGS -std=c++0x"
dnl avoid "security" checks for vectors
CXXFLAGS="$CXXFLAGS -U_GLIBCXX_ASSERTIONS"

dnl Define DOUBLEVAL if bigendian
if test "x$ac_cv_c_bigendian" = "xyes"; then
	CXXFLAGS="$CXXFLAGS -DDOUBLEVAL"
fi

dnl Check for GMP
CONFIG_GMPXX="yes"
CONFIG_TOMMATH="no"
AC_ARG_ENABLE([tommath],
	[AS_HELP_STRING([--enable-tommath], [Use libtommath instead of GMP [[default=no]]])],
	[ if test "x$enableval" = "xyes"; then CONFIG_TOMMATH="yes"; fi ], [])

CONFIG_LIBBF="no"
AC_ARG_ENABLE([libbf],
	[AS_HELP_STRING([--enable-libbf], [Use libbf instead of GMP [[default=no]]])],
	[ if test "x$enableval" = "xyes"; then CONFIG_LIBBF="yes"; fi ], [])

if test "$CONFIG_LIBBF" = "yes" ; then 
	AC_DEFINE([BF2GMP], [], [Replace GMP by libbf])
else 
if test "$CONFIG_TOMMATH" = "yes" ; then 
	AC_DEFINE([USE_GMP_REPLACEMENTS], [], [Replace GMP by tommath])
	AC_CHECK_HEADERS(tommath.h,AC_CHECK_LIB(tommath,mp_init))
else
  AC_CHECK_HEADERS(gmp.h,AC_CHECK_LIB(gmp,main))
  dnl Check for MPFR
  AC_CHECK_HEADERS(mpfr.h,AC_CHECK_LIB(mpfr,main))
  AC_CHECK_FUNCS(mpfr_set_str_raw)
  dnl Check for MPFI
  AC_CHECK_HEADERS(mpfi.h,AC_CHECK_LIB(mpfi,main))
  AC_CHECK_HEADERS(mps/mps.h,AC_CHECK_LIB(mps,main))
  AC_ARG_ENABLE([gmpxx],
	[AS_HELP_STRING([--enable-gmpxx], [Use GMP C+wrappers [[default=yes]]])],
	[if test "$enableval" = "no"; then CONFIG_GMPXX="no"; fi], [])
  if test "$CONFIG_GMPXX" = "yes"; then 
	AC_CHECK_HEADERS(gmpxx.h, [], [CONFIG_GMPXX="no"])
    fi
  if test "$CONFIG_GMPXX" = "yes"; then 
	save_LIBS="$LIBS"
	AC_CHECK_LIB(gmpxx, main, [], [CONFIG_GMPXX="no"])
	LIBS="$save_LIBS"
    fi
  GMPXX_LIBS=""
  if test "$CONFIG_GMPXX" = "yes"; then
	GMPXX_LIBS="-lgmpxx"
    fi
  AC_SUBST(CONFIG_GMPXX)
  AC_SUBST(GMPXX_LIBS)

  dnl AC_CHECK_LIB(ecm, main)
  CONFIG_ECM="yes"
  AC_ARG_ENABLE([ecm],
	[AS_HELP_STRING([--enable-ecm], [Use ECM [[default=yes]]])],
	[if test "$enableval" = "no"; then CONFIG_ECM="no"; fi], [])

  if test "$CONFIG_ECM" = "yes"; then
	AC_CHECK_HEADER(ecm.h, [], [CONFIG_ECM="no"])
    fi
  if test "$CONFIG_ECM" = "yes"; then
	AC_CHECK_LIB(ecm, main, [], [CONFIG_ECM="no"])
    fi
  AC_SUBST(CONFIG_ECM)
  AC_SUBST(ECM_LIBS)

  dnl AC_CHECK_LIB(bernmm, main)
  CONFIG_BERNMM="yes"
  AC_ARG_ENABLE([bernmm],
	[AS_HELP_STRING([--enable-bernmm], [Use BERNMM [[default=yes]]])],
	[if test "$enableval" = "no"; then CONFIG_BERNMM="no"; fi], [])

  if test "$CONFIG_BERNMM" = "yes"; then
	AC_CHECK_HEADER(bern_rat.h, [], [CONFIG_BERNMM="no"])
    fi
  if test "$CONFIG_BERNMM" = "yes"; then
	AC_CHECK_LIB(bernmm, main, [], [CONFIG_BERNMM="no"])
    fi
  AC_SUBST(CONFIG_BERNMM)
  AC_SUBST(BERNMM_LIBS)  

  dnl Checking for NTL
  CONFIG_NTL="yes"
  AC_ARG_ENABLE([ntl],
	[AS_HELP_STRING([--enable-ntl], [Use NTL [[default=yes]]])],
	[if test "$enableval" = "no"; then CONFIG_NTL="no"; fi], [])

  if test "$CONFIG_NTL" = "yes"; then
	AC_CHECK_HEADER(NTL/ZZ.h, [], [CONFIG_NTL="no"])
    fi
  if test "$CONFIG_NTL" = "yes"; then
	save_LIBS="$LIBS"
	AC_CHECK_LIB(ntl, main, [], [CONFIG_NTL="no"])
	LIBS="$save_LIBS"
    fi
  NTL_LIBS=""
  if test "$CONFIG_NTL" = "yes"; then
	NTL_LIBS="-lntl"
    fi
  AC_SUBST(CONFIG_NTL)
  AC_SUBST(NTL_LIBS)

  CONFIG_PARI="yes"
  AC_ARG_ENABLE([pari],
	[AS_HELP_STRING([--enable-pari], [Use PARI [[default=yes]]])],
	[ if test "$enableval" = "no"; then CONFIG_PARI="no"; fi], [])
  
  if test "$CONFIG_PARI" = "yes"; then
	AC_CHECK_HEADERS(pari/pari.h, [], [CONFIG_PARI="no"])
    fi
  if test "$CONFIG_PARI" = "yes"; then
	save_LIBS="$LIBS"
	AC_CHECK_LIB(pari, main, [], [CONFIG_PARI="no"])
	LIBS="$save_LIBS"
    fi

fi
fi
AM_CONDITIONAL(CONFIG_NTL, [test "$CONFIG_NTL" = "yes"])
AM_CONDITIONAL(CONFIG_ECM, [test "$CONFIG_ECM" = "yes"])
AM_CONDITIONAL(CONFIG_BERNMM, [test "$CONFIG_BERNMM" = "yes"])
AM_CONDITIONAL(CONFIG_GMPXX, [test "$CONFIG_GMPXX" = "yes"])
  
dnl Optimization or debug
CONFIG_DEBUG_P="yes"
AC_ARG_ENABLE([debug],
	[AS_HELP_STRING([--enable-debug], [enables debug support [[default=no]] ])],
	[ if test "x$enableval" = "xno"; then CONFIG_DEBUG_P="no"; fi ], [])

if test "$CONFIG_DEBUG_P" = "yes"; then
	AC_DEFINE(DEBUG_SUPPORT, ,[Set if debugging is enabled])
  fi
  
  dnl Garbage collector
CONFIG_GC="no"
AC_ARG_ENABLE([gc],
	[AS_HELP_STRING([--enable-gc], [Use Boehm garbage collector [[default=no]]])],
	[ if test "$enableval" = "yes"; then CONFIG_GC="yes"; fi ], [])

if test "$CONFIG_GC" = "xyes" ; then 
	AC_CHECK_LIB(gc,GC_malloc_atomic)
  fi
  
CONFIG_SSCL="no"
  dnl Include semi-classical
AC_ARG_ENABLE([sscl], [
	AS_HELP_STRING([--enable-sscl],
		[enable semi-classical algorithms (not working yet) [[default=no]]])],
	[ if "x$enableval" = "xyes"; then CONFIG_SSCL="yes"; fi ], [])
  
if test "x$CONFIG_SSCL" = "xyes"; then 
	AC_DEFINE(HAVE_SSCL, ,[Set to enable semi-classical routines])
fi
  
AC_CHECK_LIB(m, main)

CONFIG_PNG="yes"
AC_ARG_ENABLE(png,
	[AS_HELP_STRING([--enable-png], [Enable PNG library])],
	[ if test "x$enableval" = "xno"; then CONFIG_PNG="no"; fi], [])
  
if test "x$CONFIG_PNG" = "xyes"; then 
   AC_CHECK_HEADERS(png.h, AC_SEARCH_LIBS(png_sig_cmp,[png16 png14 png12 png]))
   dnl AC_CHECK_HEADERS(png.h, AC_CHECK_LIB(png,main)) dnl AC_CHECK_LIBS(main,[png16,png14,png12,png]))
fi

  
dnl Checking for dl
  
CONFIG_DL="yes"
AC_ARG_ENABLE(dl,
	[AS_HELP_STRING([--enable-dl], [Enable dynamically loaded modules])],
	[ if test "x$enableval" = "xno"; then CONFIG_DL="no"; fi], [])
  
if test "x$CONFIG_DL" = "xyes"; then 
	AC_CHECK_LIB(dl, dlopen)
  fi
  
AC_CHECK_SIZEOF(void *)
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(long)
AC_CHECK_SIZEOF(long long)
  
dnl Checking for Gnu Sci Lib
CONFIG_GSL="yes"
AC_ARG_ENABLE([gsl],[
	AS_HELP_STRING([--enable-gsl], [Use GNU scientific library [[default=yes]]])],
	[ if test "x$enableval" = "xno"; then CONFIG_GSL="no"; fi], [])

if test "$CONFIG_GSL" = "yes"; then
	AC_CHECK_HEADERS(gsl/gsl_blas.h)
	AC_CHECK_HEADERS(gsl/gsl_eigen.h)
	if test "$ac_cv_header_gsl_gsl_blas_h" != "yes" -o "$ac_cv_header_gsl_gsl_eigen_h" != "yes"; then
		CONFIG_GSL="no"
	fi
fi
if test "$CONFIG_GSL" = "yes"; then
	save_LIBS="$LIBS"
	AC_CHECK_LIB(gslcblas, main, [], [CONFIG_GSL="no"])
	AC_CHECK_LIB(gsl, gsl_sf_gamma, [], [CONFIG_GSL="no"])
	LIBS="$save_LIBS"
  fi
GSL_LIBS=""
if test "$CONFIG_GSL" = "yes"; then
	GSL_LIBS="-lgsl -lgslcblas"
  fi
AC_SUBST(CONFIG_GSL)
AC_SUBST(GSL_LIBS)
AM_CONDITIONAL(CONFIG_GSL, [test "$CONFIG_GSL" = "yes"])
  
CONFIG_LAPACK="yes"
AC_ARG_ENABLE([lapack],
	[AS_HELP_STRING([--enable-lapack], [Use LAPACK [[default=yes]]])],
	[ if test "$enableval" = "no"; then CONFIG_LAPACK="no"; fi], [])
  
if test "$CONFIG_LAPACK" = "yes"; then
#	AC_CHECK_LIB(f2c, main, [], [CONFIG_LAPACK="no"])
#	AC_CHECK_LIB(blas, main, [], [CONFIG_LAPACK="no"])
#	AC_CHECK_LIB(tmglib, main, [], [CONFIG_LAPACK="no"])
	AC_CHECK_LIB(gfortran, main)
	AC_CHECK_LIB(blas, main)
	AC_CHECK_LIB(lapack, main, [], [CONFIG_LAPACK="no"])
#        AX_BLAS([have_blas=yes],[have_blas=no])
#        AX_LAPACK([have_lapack=yes],[have_lapack=no])
fi
  
dnl if test "$CONFIG_PARI" = "yes"; then
dnl                      AC_MSG_CHECKING([whether PARI is configured with pthreads])
dnl                      AC_RUN_IFELSE([AC_LANG_PROGRAM([
dnl                          [#include <pari/pari.h>
dnl                           #include <string.h>]],
dnl                         [[return strcmp(PARI_MT_ENGINE, "pthread") != 0]])],
dnl                        [AC_MSG_RESULT([yes. Good])],
dnl                        [AC_MSG_RESULT([no. Not good])
dnl                         CONFIG_PARI="no"])
dnl  fi

PARI_LIBS=""
if test "$CONFIG_PARI" = "yes"; then
	PARI_LIBS="-lpari"
fi
AC_SUBST(PARI_LIBS)
AC_SUBST(CONFIG_PARI)
AM_CONDITIONAL(CONFIG_PARI, [test "$CONFIG_PARI" = "yes"])
  


CONFIG_COCOA="yes"
AC_ARG_ENABLE([cocoa],
	[AS_HELP_STRING([--enable-cocoa], [Use CoCoA library [[default=autodetect]]])],
	[if test "$enableval" = "no"; then CONFIG_COCOA="no"; fi], [])

if test "$CONFIG_COCOA" = "yes"; then
	AC_CHECK_HEADER(CoCoA/io.H, [], [CONFIG_COCOA="no"])
  fi
if test "$CONFIG_COCOA" = "yes"; then
	save_LIBS="$LIBS"
	AC_CHECK_LIB(cocoa, main, [], [CONFIG_COCOA="no"])
	LIBS="$save_LIBS"
fi
COCOA_LIBS=""
if test "$CONFIG_COCOA" = "yes"; then
	COCOA_LIBS="-lcocoa"
fi
AC_SUBST(CONFIG_COCOA)
AC_SUBST(COCOA_LIBS)
AM_CONDITIONAL(CONFIG_COCOA, [test "$CONFIG_COCOA" = "yes"])

dnl AC_CHECK_LIB(ao, main)
CONFIG_AO="yes"
AC_ARG_ENABLE([ao],
	[AS_HELP_STRING([--enable-ao], [Use AO [[default=yes]]])],
	[if test "$enableval" = "no"; then CONFIG_AO="no"; fi], [])

if test "$CONFIG_AO" = "yes"; then
	AC_CHECK_HEADER(ao/ao.h, [], [CONFIG_AO="no"])
  fi
if test "$CONFIG_AO" = "yes"; then
	AC_CHECK_LIB(ao, main, [], [CONFIG_AO="no"])
  fi
AC_SUBST(CONFIG_AO)
AC_SUBST(AO_LIBS)
AM_CONDITIONAL(CONFIG_AO, [test "$CONFIG_AO" = "yes"])

dnl AC_CHECK_LIB(glpk, main)
CONFIG_GLPK="yes"
AC_ARG_ENABLE([glpk],
	[AS_HELP_STRING([--enable-glpk], [Use GLPK [[default=yes]]])],
	[if test "$enableval" = "no"; then CONFIG_GLPK="no"; fi], [])

if test "$CONFIG_GLPK" = "yes"; then
	AC_CHECK_HEADER(glpk.h, [], [CONFIG_GLPK="no"])
  fi
if test "$CONFIG_GLPK" = "yes"; then
	AC_CHECK_LIB(glpk, main, [], [CONFIG_GLPK="no"])
  fi
AC_SUBST(CONFIG_GLPK)
AC_SUBST(GLPK_LIBS)
AM_CONDITIONAL(CONFIG_GLPK, [test "$CONFIG_GLPK" = "yes"])

dnl AC_CHECK_LIB(samplerate, main)
CONFIG_SAMPLERATE="yes"
AC_ARG_ENABLE([samplerate],
	[AS_HELP_STRING([--enable-samplerate], [Use SAMPLERATE [[default=yes]]])],
	[if test "$enableval" = "no"; then CONFIG_SAMPLERATE="no"; fi], [])

if test "$CONFIG_SAMPLERATE" = "yes"; then
	AC_CHECK_HEADER(samplerate.h, [], [CONFIG_SAMPLERATE="no"])
  fi
if test "$CONFIG_SAMPLERATE" = "yes"; then
	AC_CHECK_LIB(samplerate, main, [], [CONFIG_SAMPLERATE="no"])
  fi
AC_SUBST(CONFIG_SAMPLERATE)
AC_SUBST(SAMPLERATE_LIBS)
AM_CONDITIONAL(CONFIG_SAMPLERATE, [test "$CONFIG_SAMPLERATE" = "yes"])

dnl AC_CHECK_LIB(curl, main)
CONFIG_CURL="yes"
AC_ARG_ENABLE([curl],
	[AS_HELP_STRING([--enable-curl], [Use CURL [[default=yes]]])],
	[if test "$enableval" = "no"; then CONFIG_CURL="no"; fi], [])

if test "$CONFIG_CURL" = "yes"; then
	AC_CHECK_HEADER(curl/curl.h, [], [CONFIG_CURL="no"])
  fi
if test "$CONFIG_CURL" = "yes"; then
	AC_CHECK_LIB(curl, main, [], [CONFIG_CURL="no"])
  fi
AC_SUBST(CONFIG_CURL)
AC_SUBST(CURL_LIBS)
AM_CONDITIONAL(CONFIG_CURL, [test "$CONFIG_CURL" = "yes"])

  
AC_CANONICAL_HOST
dnl Check for stuff needed for building the Giac interactive shell (cas).
GINAC_READLINE

dnl Ugh...
CONFIG_GUI="yes"
AC_ARG_ENABLE([gui],
	[AS_HELP_STRING([--enable-gui], [The disable-gui option has no more effect because FLTK 1.3.8 is now included in Giac/Xcas. If FLTK compilation fails or if your local FLTK install is incompatible, you should configure with --disable-fltk])],
	[ if test "x$enableval" = "xno"; then CONFIG_GUI="no"; fi], [])

CONFIG_FLTK="yes"
AC_ARG_ENABLE([fltk],
	[AS_HELP_STRING([--enable-fltk], [enable FLTK [[default=yes]]])],
	[ if test "x$enableval" = "xno"; then CONFIG_FLTK="no"; fi], [])
	
save_LIBS="$LIBS"
if test "$CONFIG_FLTK" = "yes"; then 

dnl Check for X11 anf fltk or cygwin and fltk...
	case $host_os in
		*cygwin*)
		AC_CHECK_LIB(GL, main)
dnl		AC_CHECK_LIB(GLU, main)
		AC_CHECK_LIB(fltk_gl, main)
		AC_CHECK_LIB(fltk, main)
		AC_CHECK_LIB(z, main)
		AC_CHECK_LIB(jpeg, main)
		AC_CHECK_LIB(fltk_images, main)
		if test "$HAVE_LIBFLTK" = "1" ; then
			AC_DEFINE([HAVE_LIBFLVW],1,[Now defined if fltk is available])
		fi
		LIBS="$LIBS -mwindows -lole32 -luuid -lcomctl32 -lwsock32"
		if test "$HAVE_LIBFLTK_GL" = "1" ; then
		    AC_DEFINE([USE_OPENGL32],1,[Defined if fltk has OpenGL enabled])
		    LIBS="$LIBS -lglu32 -lopengl32 -lfltk_png"
		fi

		CXXFLAGS="$CXXFLAGS -DWIN32"
		AC_CHECK_LIB(intl.dll, main)
		AC_CHECK_LIB(intl, main)
		;;

		*mingw*)
		CXXFLAGS="$CXXFLAGS -D__MINGW_H"
		CFLAGS="$CFLAGS -D__MINGW_H"
		LIBS="$LIBS -lopengl32"
		AC_CHECK_LIB(fltk_gl, main)
		AC_CHECK_LIB(fltk, main)
		if test "$HAVE_LIBFLTK_GL" = "1" ; then
		    AC_DEFINE([USE_OPENGL32],1,[Defined if fltk has OpenGL enabled])
		fi		
		AC_CHECK_LIB(z, main)
		AC_CHECK_LIB(jpeg, main)
		AC_CHECK_LIB(fltk_images, main)
		;;

		*darwin*)
		dnl FLTK does not need X libraries on this platform
		;;
		*)
		AC_PATH_XTRA
		AC_CHECK_LIB(GL, main)
dnl		AC_CHECK_LIB(GLU, main)
		;;
	esac
	dnl force FLTK static compilation if FLTK not present
	AC_CHECK_HEADERS(FL/Fl.H,[],[tar xfj fltk-1.3.8-source.tar.bz2 ; cd fltk-1.3.8 && ./configure --disable-shared && make && cd .. && CXXFLAGS="$CXXFLAGS -I../fltk-1.3.8 -I./fltk-1.3.8" && LDFLAGS="$LDFLAGS -Lfltk-1.3.8/lib "] ) 
	dnl X11, check for opengl, fltk
	AC_CHECK_LIB(X11,main)
	AC_CHECK_LIB(Xcursor,main)
	AC_CHECK_LIB(Xfixes,main)
	AC_CHECK_LIB(Xinerama, main)
	AC_CHECK_LIB(fontconfig, main)
	AC_CHECK_LIB(Xft, main,[HAVE_LIBXFT=1]) dnl  ,[], [AC_MSG_ERROR(exiting, as libXft cannot be found. Hint: you can run ./configure --disable-fltk to compile without GUI)])
        if test "$HAVE_LIBXFT" = "1" ; then
           AC_CHECK_LIB(Xft, main)
	   AC_CHECK_LIB(Xext, main)
	   AC_CHECK_LIB(Xrender, main)
	   AC_CHECK_LIB(fltk_gl, main)
	   AC_CHECK_LIB(fltk, main)
	   if test "$HAVE_LIBFLTK" = "1"; then
	      AC_DEFINE([HAVE_LIBFLVW],1, [Now defined if fltk is available])
	      fi
	      AC_CHECK_LIB(z,main)
	      AC_CHECK_HEADERS(png.h, AC_CHECK_LIB(png,main))
	      AC_CHECK_LIB(jpeg, main)
	      AC_CHECK_LIB(fltk_images, main)
	      FLTK_LIBS="$LIBS $X_LIBS"
	      AC_SUBST(FLTK_LIBS)
              cd src && ln -sf ../fltk-1.3.8 && cd ..
           fi
        fi
LIBS="$save_LIBS"

AC_CHECK_HEADERS(dfu_lib.h) 
AC_CHECK_LIB(usb-1.0,main) 
AC_CHECK_LIB(dfu,dfu_main) 
AC_CHECK_LIB(cliquer,main) 
AC_CHECK_LIB(nauty,main)
AC_CHECK_HEADERS(nauty/naututil.h)
dnl https://github.com/vectorclass by Agner Fog
AC_CHECK_HEADERS(vcl1/vectorclass.h)
  
AC_CHECK_LIB(pthread,main)
AC_CHECK_HEADERS(pthread.h)
AC_CHECK_HEADERS(malloc.h,AC_CHECK_FUNCS(mallinfo))
AC_CHECK_FUNCS(sysconf)
AC_CHECK_LIB(rt,clock_gettime)
  
  dnl Checking for hash_map
AC_CHECK_HEADER(unordered_map, AC_DEFINE(C11_UNORDERED_MAP,1, [Define if <unordered_map> header is available]))
AC_CHECK_HEADER(ext/hash_map, AC_DEFINE(EXT_HASH_MAP,1, [Define if <ext/hash_map> header is aviailable]))
AC_CHECK_HEADER(tr1/unordered_map, AC_DEFINE(UNORDERED_MAP,1, [Define if <tr1/unordered_map> header is available]))
AC_CHECK_HEADER(hash_map,AC_DEFINE(HASH_MAP,1, [Define if <hash_map> header is aviailable]))
AC_CHECK_HEADERS(pwd.h, ,AC_DEFINE(HAVE_NO_PWD_H,1, [Set if pwd.h is NOT available]))
AC_CHECK_HEADERS(locale.h)
AC_CHECK_HEADERS(unistd.h)
AC_CHECK_HEADERS(sys/time.h)
AC_CHECK_HEADERS(sys/times.h, , AC_DEFINE(HAVE_NO_SYS_TIMES_H,1, [Set if <sys/times.h> is NOT available]))
AC_CHECK_HEADERS(sys/resource.h, , AC_DEFINE(HAVE_NO_SYS_RESOURCE_WAIT_H,1, [Set if <sys/resource.h> is NOT available]))
AC_CHECK_HEADERS(sys/types.h)
AC_CHECK_HEADERS(signal.h)
AC_CHECK_FUNCS(getpid, , AC_DEFINE(HAVE_NO_SIGNAL_H,1, [Set if getpid() function is NOT available]))
AC_CHECK_FUNCS(getcwd, , AC_DEFINE(HAVE_NO_CWD,1, [Set if getcwd() function is NOT available]))
AC_CHECK_FUNCS(getpwuid, , AC_DEFINE(HAVE_NO_HOME_DIRECTORY,1, [Set if getpwuid() function is NOT available]))
AC_CHECK_FUNCS(system, , AC_DEFINE(HAVE_NO_SYSTEM, 1, [Set if system() function is NOT available]))
  
ALL_LINGUAS="es fr el pt it en zh de"
AM_GNU_GETTEXT
AM_GNU_GETTEXT_VERSION([0.14.5])
dnl for Slowlaris comment next line
dnl AC_CHECK_LIB(intl, main) 
dnl auto-check will work if the function checked is alone in a file
dnl and independant from the whole micropython library
dnl otherwise it will fail because it depends on giac
dnl AC_CHECK_LIB(micropython,mp_version)

AC_PROG_EGREP
CONFIG_DOC_PROGS="yes"
AC_PATH_PROG(LATEX, latex, "")
if test -z "$LATEX"; then CONFIG_DOC_PROGS = "no"; fi

if test "$CONFIG_DOC_PROGS" = "yes"; then
	AC_PATH_PROG(MAKEINDEX, makeindex, "")
	if test -z "$MAKEINDEX"; then CONFIG_DOC_PROGS = "no"; fi
fi

if test "$CONFIG_DOC_PROGS" = "yes"; then
	AC_PATH_PROG(PDFLATEX, pdflatex, "")
	if test -z "$PDFLATEX"; then CONFIG_DOC_PROGS = "no"; fi
fi
if test "$CONFIG_DOC_PROGS" = "yes"; then
	AC_PATH_PROG(DVIPS, dvips, "")
	if test -z "$DVIPS"; then CONFIG_DOC_PROGS = "no"; fi
fi
dnl if test "$CONFIG_DOC_PROGS" = "yes"; then
dnl 	AC_PATH_PROG(LATEX2HTML, latex2html, "")
dnl	if test -z "$LATEX2HTML"; then CONFIG_DOC_PROGS = "no"; fi
dnl fi
AC_SUBST(CONFIG_DOC_PROGS)
AM_CONDITIONAL(DOC_PROGS, [test "$CONFIG_DOC_PROGS" = "yes"])

dnl Checking for sstream for mathml
AC_CHECK_HEADERS(sstream, AC_DEFINE(HAVE_SSTREAM,1), )

dnl CXXFLAGS="$CXXFLAGS -fno-strict-aliasing -DSMARTPTR64"
CXXFLAGS="$CXXFLAGS -fno-strict-aliasing -DGIAC_GENERIC_CONSTANTS -DTIMEOUT"

dnl enables MicroPython, -fpermissive to avoid failure on Debian 10
CONFIG_MICROPY="yes"
AC_ARG_ENABLE([micropy],
	[AS_HELP_STRING([--enable-micropy], [enable MICROPY [[default=yes]]])],
	[ if test "x$enableval" = "xno"; then CONFIG_MICROPY="no"; fi], [])
	
if test "$CONFIG_MICROPY" = "yes"; then 
  LDFLAGS="$LDFLAGS -L. -L.. -lmicropython"
  CXXFLAGS="$CXXFLAGS -DHAVE_LIBMICROPYTHON -fpermissive"
fi
AM_CONDITIONAL([BUILD_MICROPY], [test "$CONFIG_MICROPY" = "yes"])

dnl enables QuickJS, 
CONFIG_QUICKJS="yes"
AC_ARG_ENABLE([quickjs],
	[AS_HELP_STRING([--enable-quickjs], [enable QUICKJS [[default=yes]]])],
	[ if test "x$enableval" = "xno"; then CONFIG_QUICKJS="no"; fi], [])
	
if test "$CONFIG_QUICKJS" = "yes"; then 
  CXXFLAGS="$CXXFLAGS -DQUICKJS"
  LDFLAGS="$LDFLAGS"
  CFLAGS="$CFLAGS -D_GNU_SOURCE -DQUICKJS -DCONFIG_BIGNUM -DCONFIG_VERSION=\\\"2020-11-08\\\""
fi
AM_CONDITIONAL([BUILD_QUICKJS], [test "$CONFIG_QUICKJS" = "yes"])

dnl Output makefiles etc.
AC_OUTPUT([
Makefile
src/Makefile
check/Makefile
doc/Makefile
doc/en/Makefile
doc/es/Makefile
doc/fr/Makefile
doc/el/Makefile
doc/el/tutoriel/Makefile
doc/el/casinter/Makefile
doc/local/Makefile
doc/pt/Makefile
doc/zh/Makefile
doc/de/Makefile
intl/Makefile
po/Makefile.in
examples/Makefile
debianold/Makefile
debian/Makefile
micropython-1.12/Makefile
micropython-1.12/docs/Makefile
micropython-1.12/drivers/Makefile
micropython-1.12/examples/Makefile
micropython-1.12/extmod/Makefile
micropython-1.12/lib/Makefile
micropython-1.12/logo/Makefile
micropython-1.12/ports/Makefile
micropython-1.12/py/Makefile
micropython-1.12/tests/Makefile
micropython-1.12/tools/Makefile
micropython-1.12/unix/Makefile
micropython-1.12/ulab/Makefile
micropython-1.12/ulab/numpy/Makefile
micropython-1.12/ulab/numpy/fft/Makefile
micropython-1.12/ulab/numpy/ndarray/Makefile
micropython-1.12/ulab/numpy/linalg/Makefile
micropython-1.12/ulab/scipy/Makefile
micropython-1.12/ulab/scipy/linalg/Makefile
micropython-1.12/ulab/scipy/optimize/Makefile
micropython-1.12/ulab/scipy/signal/Makefile
micropython-1.12/ulab/scipy/special/Makefile
micropython-1.12/ulab/user/Makefile
micropython-1.12/ulab/utils/Makefile
])

echo "Adding link . to giac in src"
rm -f src/giac && ln -s . src/giac
rm -f src/config.h && cp config.h src
rm -f path.h
echo "#define giac_html_location \"$prefix/share/giac/doc/\"" >> path.h
echo "#define giac_aide_location \"$prefix/share/giac/aide_cas\"" >> path.h
echo "#ifdef __APPLE__" >> path.h
echo "#define giac_locale_location \"/Applications/usr/share/locale/\"" >> path.h
echo "#else" >> path.h
echo "#define giac_locale_location \"$prefix/share/locale/\"" >> path.h
echo "#endif" >> path.h
echo "#define giac_gnuplot_location \"$GNUPLOT\"" >> path.h
rm -f src/path.h && cp path.h src
