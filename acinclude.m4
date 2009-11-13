dnl old names
AM_LC_MESSAGES
AM_WITH_NLS
AM_GNU_GETTEXT
AM_PATH_PROG_WITH_TEST
dnl AC_DEFUN([AM_PROG_LIBTOOL],   [AC_PROG_LIBTOOL])
dnl AC_DEFUN([AM_ENABLE_SHARED],  [AC_ENABLE_SHARED($@)])
dnl AC_DEFUN([AM_ENABLE_STATIC],  [AC_ENABLE_STATIC($@)])
dnl AC_DEFUN([AM_DISABLE_SHARED], [AC_DISABLE_SHARED($@)])
dnl AC_DEFUN([AM_DISABLE_STATIC], [AC_DISABLE_STATIC($@)])
dnl AC_DEFUN([AM_PROG_LD],        [AC_PROG_LD])
dnl AC_DEFUN([AM_PROG_NM],        [AC_PROG_NM])

dnl Additional macros used by configure.
dnl Usage: GIAC_RLVERSION
dnl The maintainers of libreadline are complete morons: they don't care a shit
dnl about compatiblilty (which is not so bad by itself) and at the same time 
dnl they don't export the version to the preprocessor so we could kluge around 
dnl incomatiblities.  The only reliable way to figure out the version is by 
dnl checking the extern variable rl_library_version at runtime.  &#@$%*!
AC_DEFUN(GIAC_LIB_READLINE_VERSION,
[AC_CACHE_CHECK([for version of libreadline], giac_cv_rlversion, [
AC_TRY_RUN([
#include <stdio.h>
#include <sys/types.h>
#include <readline/readline.h>

main()
{
    FILE *fd;
    fd = fopen("conftest.out", "w");
    fprintf(fd, "%s\n", rl_library_version);
    fclose(fd);
    return 0;
}], giac_cv_rlversion=`cat 'conftest.out'`, giac_cv_rlversion='unknown', giac_cv_rlversion='4.2')])
if test "x${giac_cv_rlversion}" != "xunknown"; then
  RL_VERSION_MAJOR=`echo ${giac_cv_rlversion} | sed -e 's/\([[0-9]]*\)\.\([[0-9]]*\).*/\1/'`
  AC_DEFINE_UNQUOTED(GIAC_RL_VERSION_MAJOR, $RL_VERSION_MAJOR, [Major version of installed readline library.])
  RL_VERSION_MINOR=`echo ${giac_cv_rlversion} | sed -e 's/\([[0-9]]*\)\.\([[0-9]]*\).*/\2/'`
  AC_DEFINE_UNQUOTED(GIAC_RL_VERSION_MINOR, $RL_VERSION_MINOR, [Minor version of installed readline library.])
else
  GIAC_WARNING([I could not run a test of libreadline (needed for building ginsh).])
fi
])

dnl Usage: GIAC_TERMCAP
dnl libreadline is based on the termcap functions.
dnl Some systems have tgetent(), tgetnum(), tgetstr(), tgetflag(), tputs(),
dnl tgoto() in libc, some have it in libtermcap, some have it in libncurses.
dnl When both libtermcap and libncurses exist, we prefer the latter, because
dnl libtermcap is being phased out.
AC_DEFUN(GIAC_TERMCAP,
[LIBTERMCAP=
AC_CHECK_FUNCS(tgetent)
if test "x$ac_cv_func_tgetent" = "xyes"; then
    :
else
    AC_CHECK_LIB(ncurses, tgetent, LIBTERMCAP="-lncurses")
    if test -z "$LIBTERMCAP"; then
        AC_CHECK_LIB(termcap, tgetent, LIBTERMCAP="-ltermcap")
    fi
fi
AC_SUBST(LIBTERMCAP)
])


# Is the gmp header file new enough? (should be implemented with an argument)
AC_DEFUN(GIAC_GMP_H_VERSION,
[AC_CACHE_CHECK([for recent enough gmp.h], cl_cv_new_gmp_h, [
  AC_TRY_CPP([#include <gmp.h>
#if !defined(__GNU_MP_VERSION) || (__GNU_MP_VERSION < 2)
 #error "ancient gmp.h"
#endif],
cl_cv_new_gmp_h="yes",
cl_cv_new_gmp_h="no"; GIAC_ERROR([The header file for GMP version 2 or above could not be found.]))
])])dnl

# Does libgmp provide some functionality introduced in version 3.0?
AC_DEFUN(GIAC_GMP_CHECK,
[AC_CACHE_CHECK([for working libgmp], cl_cv_new_libgmp, [
    SAVELIBS=$LIBS
    LIBS="$LIBS -lgmp"
    AC_TRY_LINK([#include <gmp.h>],[mpz_sqrt(0,0)],
cl_cv_new_libgmp="yes", 
cl_cv_new_libgmp="no"; LIBS=$SAVELIBS; GIAC_ERROR([A test program could not be linked against GMP version 2 or above.]))
])])

dnl Usage: GIAC_ERROR(message)
dnl This macro displays the warning "message" and sets the flag giac_error
dnl to yes.
AC_DEFUN(GIAC_ERROR,[
giac_error_txt="$giac_error_txt
** $1
"
giac_error=yes])

dnl Usage: GIAC_WARNING(message)
dnl This macro displays the warning "message" and sets the flag giac_warning
dnl to yes.
AC_DEFUN(GIAC_WARNING,[
giac_warning_txt="$giac_warning_txt
== $1
"
giac_warning=yes])

dnl Usage: GIAC_CHECK_ERRORS
dnl (preferably to be put at end of configure.in)
dnl This macro displays a warning message if GIAC_ERROR or GIAC_WARNING
dnl has occured previously.
AC_DEFUN(GIAC_CHECK_ERRORS,[
if test "x${giac_error}" = "xyes"; then
    echo "**** The following problems have been detected by configure."
    echo "**** Please check the messages below before running \"make\"."
    echo "**** (see the section 'Common Problems' in the INSTALL file)"
    echo "$giac_error_txt"
    if test "x${giac_warning_txt}" != "x"; then
        echo "${giac_warning_txt}"
    fi
    echo "deleting cache ${cache_file}"
    rm -f $cache_file
    else
        if test x$giac_warning = xyes; then
            echo "==========================================================="
            echo "=== configure has detected the following install options."
            echo "=== Please check the messages below before running \"make\"."
            echo "=== If configure was wrong, you can modify config.h and"
            echo "=== add manually libraries -Lpath -llibname in src/Makefile"
            echo "=== at the line LIBS ="
            echo "=== If built seems too long, edit src/Makefile and remove "
            echo "=== -O2 at the line CXXFLAGS ="
            echo "=== (for more details see the INSTALL file)"
            echo "$giac_warning_txt"
        fi
    echo "Configuration done. Now type \"make\"."
fi])
