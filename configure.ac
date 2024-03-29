
dnl /**-------------------------------------------------------------------**
dnl  **                              CLooG                                **
dnl  **-------------------------------------------------------------------**
dnl  **                           configure.in                            **
dnl  **-------------------------------------------------------------------**
dnl  **                   First version: august 7th 2002                  **
dnl  **-------------------------------------------------------------------**/
dnl
dnl Input file for autoconf to build a configuration shellscript.
dnl To build the configure script from the CLooG's top-level directory, use
dnl autoconf -l autoconf autoconf/configure.in > configure
dnl if it doesn't work (invalid option -l) try -I instead
dnl autoconf -I autoconf autoconf/configure.in > configure

dnl /**************************************************************************
dnl *               CLooG : the Chunky Loop Generator (experimental)          *
dnl ***************************************************************************
dnl *                                                                         *
dnl * Copyright (C) 2001 Cedric Bastoul                                       *
dnl *                                                                         *
dnl * This is free software; you can redistribute it and/or modify it under   *
dnl * the terms of the GNU General Public License as published by the Free    *
dnl * Software Foundation; either version 2 of the License, or (at your       *
dnl * option) any later version.                                              *
dnl *                                                                         *
dnl * This software is distributed in the hope that it will be useful, but    *
dnl * WITHOUT ANY WARRANTY; without even the implied warranty of              *
dnl * MERCHANTABILITYor FITNESS FOR A PARTICULAR PURPOSE.  See the GNU        *
dnl * General Public License for more details.                                *
dnl *                                                                         *
dnl * You should have received a copy of the GNU General Public License along *
dnl * with software; if not, write to the Free Software Foundation, Inc.,     *
dnl * 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA                  *
dnl *                                                                         *
dnl * CLooG, the Chunky Loop Generator                                        *
dnl * Written by Cedric Bastoul, Cedric.Bastoul@inria.fr                      *
dnl *                                                                         *
dnl ***************************************************************************/

AC_PREREQ(2.13)
AC_INIT(source/cloog.c)
AC_CONFIG_AUX_DIR(autoconf)

dnl Every other copy of the package version number gets its value from here.
AM_INIT_AUTOMAKE(cloog, 0.15)
AC_CONFIG_HEADER(include/cloog/cloog-config.h)

dnl default version
BITS="64"

AC_SUBST(VERSION)

# Version number machinery.
changequote(<<, >>)dnl
if test -n "`expr $VERSION : '\([0-9]*\)\.[0-9]*\.[0-9]*'`"
then
  CLOOG_VERSION_MAJOR=`expr $VERSION : '\([0-9]*\)\.[0-9]*\.[0-9]*'`
  CLOOG_VERSION_MINOR=`expr $VERSION : '[0-9]*\.\([0-9]*\)\.[0-9]*'`
  CLOOG_VERSION_REVISION=`expr $VERSION : '[0-9]*\.[0-9]*\.\([0-9]*\)'`
  CLOOG_VERSION_BETA=`expr $VERSION : '[0-9]*\.[0-9]*\.[0-9]*pre\([0-9]*\)'`
else
  CLOOG_VERSION_MAJOR=`expr $VERSION : '\([0-9]*\)\.[0-9]*'`
  CLOOG_VERSION_MINOR=`expr $VERSION : '[0-9]*\.\([0-9]*\)'`
  CLOOG_VERSION_REVISION=0
  CLOOG_VERSION_BETA=`expr $VERSION : '[0-9]*\.[0-9]*pre\([0-9]*\)'`
fi
if test -z "$CLOOG_VERSION_BETA"
then
  CLOOG_VERSION_BETA=0
fi
changequote([, ])dnl
AC_SUBST(CLOOG_VERSION_MAJOR)
AC_SUBST(CLOOG_VERSION_MINOR)
AC_SUBST(CLOOG_VERSION_REVISION)
AC_SUBST(CLOOG_VERSION_BETA)


dnl /**************************************************************************
dnl  *                              Checking                                  *
dnl  **************************************************************************/

dnl Checks for programs.
AC_PROG_CC
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_CHECK_PROG(CD, cd)
dnl Configure needs an empty install.sh file with this, i HATE that...
AC_PROG_INSTALL
AC_PROG_LIBTOOL


dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T



dnl Checks for header files.
AC_HEADER_STDC


dnl Checks for library functions.
AC_CHECK_FUNCS(strtol)

AC_CHECK_HEADER(sys/resource.h,[AC_DEFINE([HAS_SYS_RESOURCE_H],1,[presence of sys/resource.h])],
	        [AC_MSG_RESULT([no])])

dnl /**************************************************************************
dnl  *                             Option setting                             *
dnl  **************************************************************************/

dnl Some default values cause I'm not sure whether autoconf set them, while
dnl documentation says it does...
gmp_package="yes"
gmp_include_package="yes"
gmp_library_package="yes"

MP_BITS="MP"
AC_DEFINE([LINEAR_VALUE_IS_MP], 1, [Empty comment])
NEED_MP="no"

dnl --with-ppl=ppl-path
AC_ARG_WITH(ppl,
	    [AS_HELP_STRING([--with-ppl=DIR],
			    [Prefix of Parma Polyhedral Library installation])])

dnl --with-polylib=polylib-path
AC_ARG_WITH(polylib_prefix,
	    [AS_HELP_STRING([--with-polylib-prefix=DIR],
			    [Prefix of PolyLib installation])])
AC_ARG_WITH(polylib_exec_prefix,
	    [AS_HELP_STRING([--with-polylib-exec-prefix=DIR],
			    [Exec prefix of PolyLib installation])])
AC_ARG_WITH(polylib_builddir,
	    [AS_HELP_STRING([--with-polylib-builddir=DIR],
			    [Location of PolyLib builddir])])
AC_ARG_WITH(bits,
	    [  --with-bits=32|64|gmp   PolyLib library type to use],
	    [ echo "polylib$withval"; BITS=$withval ])

dnl --with-gmp=gmp-path
AC_ARG_WITH(gmp,
        [  --with-gmp=DIR          DIR where the gmp package is installed],
	[ echo "Package gmp : $withval" &&
	  gmp_package=$withval &&
          GMP_INC=$gmp_package/include &&
          GMP_LIB=$gmp_package/lib &&
          NEED_MP="yes"])

AC_ARG_WITH(gmp-include,
        [  --with-gmp-include=DIR  DIR where gmp.h is installed],
	[ echo "Package gmp-include : $withval" &&
	  gmp_include_package=$withval &&
          GMP_INC=$gmp_include_package &&
	  NEED_MP="yes"])

AC_ARG_WITH(gmp-library,
        [  --with-gmp-library=DIR  DIR where the gmp library is installed],
	[ echo "Package gmp-library : $withval" &&
	  gmp_library_package=$withval &&
          GMP_LIB=$gmp_library_package &&
	  NEED_MP="yes"])


dnl /**************************************************************************
dnl  *                            Where is GMP?                               *
dnl  **************************************************************************/


dnl Checking for gmp
AC_MSG_CHECKING(whether gmp works)
if test "$gmp_package" = "no"; then
  echo "GMP package not defined"
  AC_MSG_RESULT(no)
  TO_BUILD_MP=""
else
  if test "$NEED_MP" = "no"; then
    echo "Mode normal GMP"
    TO_BUILD="$TO_BUILD MP"
    AC_CHECK_HEADER(gmp.h,
                    [AC_CHECK_LIB(gmp,
		                  __gmpz_init,
		                  [LIBS="$LIBS -lgmp"],
				  [echo "Can't find gmp library." &&
				   echo "MP version will not be builded." &&
				   TO_BUILD_MP=""])],
		    [echo "Can't find gmp headers." &&
		     echo "MP version will not be builded." &&
		     TO_BUILD_MP=""])
  else
    dnl Default given by --with-X is "yes", --without-X is "no". We also
    dnl initialized manually all gmp_package* variables to "yes" (thus they are
    dnl supposed to be "yes" except if the user set them himself).
    
    if test "$gmp_package" != "yes" ; then
      echo "(GMP path has been set by user)"
      GMP_DIR=$gmp_package
      dnl Useful for AC_CHECK_X to find what we want.
      CPPFLAGS="-I$GMP_DIR/include $CPPFLAGS"
      LDFLAGS="-L$GMP_DIR/lib $LDFLAGS"
    fi
    
    if test "$gmp_include_package" != "yes" ; then
      CPPFLAGS="-I$GMP_INC $CPPFLAGS"
    fi
    
    if test "$gmp_library_package" != "yes" ; then
      LDFLAGS="-L$GMP_LIB $LDFLAGS"
    fi
     
    AC_CHECK_HEADER(gmp.h,
                    [],
                    [AC_MSG_ERROR(Can't find gmp headers.)])
    AC_CHECK_LIB(gmp,
                 __gmpz_init,
                 [LIBS="$LIBS -lgmp"],
                 [AC_MSG_ERROR(Can't find gmp library.)])
    
    AC_MSG_RESULT(yes)
  fi
fi

dnl /**************************************************************************
dnl  *                          Where is PolyLib?                             *
dnl  **************************************************************************/

polylibs="64 32"
case "$BITS" in
    32|64|gmp)
	polylibs=$BITS
    ;;
esac

dnl Checking for PolyLib
AC_SUBST(cl_cv_polylib)
AC_SUBST(POLYLIB_CPPFLAGS)
AC_SUBST(POLYLIB_LDFLAGS)
AC_SUBST(POLYLIB_LIBS)

if test "x$with_polylib_prefix" != "x" -a "x$with_polylib_exec_prefix" = "x"
then
	with_polylib_exec_prefix=$with_polylib_prefix
fi
if test "x$with_polylib_prefix" != "x" -a "x$with_polylib_builddir" != "x"; then
	AC_MSG_ERROR(
	[--with-polylib-prefix and --with-polylib-builddir are mutually exclusive])
fi

AC_MSG_CHECKING(for location of PolyLib)
if test "x$with_polylib_builddir" != "x"; then
	AC_MSG_RESULT(built in $with_polylib_builddir)
	polylib_srcdir=`echo @srcdir@ | $with_polylib_builddir/config.status --file=-`
	AC_MSG_NOTICE(PolyLib sources in $polylib_srcdir)
	polylib_polylibs=`echo @polylibs@ | $with_polylib_builddir/config.status --file=-`
	AC_MSG_NOTICE(Configured polylibs: $polylib_polylibs)
	CPPFLAGS="-I$with_polylib_builddir/include -I$polylib_srcdir/include $CPPFLAGS"
	for BITS in $polylibs; do
		cl_cv_polylib=missing
		for PL in $polylib_polylibs; do
			if test "libpolylib$BITS.la" = $PL; then
				cl_cv_polylib="polylib$BITS"
				break;
			fi
		done
		if test "$cl_cv_polylib" != "missing"; then
			LIBS="$with_polylib_builddir/lib$cl_cv_polylib.la $LIBS"
			break;
		fi
	done
else
	if test "x$with_polylib_prefix" != "x"; then
		AC_MSG_RESULT(installed in $with_polylib_exec_prefix)
		CPPFLAGS="-I$with_polylib_prefix/include $CPPFLAGS"
		LDFLAGS="-L$with_polylib_exec_prefix/lib $LDFLAGS"
		LD_LIBRARY_PATH="$with_polylib_exec_prefix/lib:$LD_LIBRARY_PATH"
		export LD_LIBRARY_PATH
	else
		AC_MSG_RESULT(installed in standard location)
	fi
	for BITS in $polylibs; do
		cl_cv_polylib="polylib$BITS"
		AC_CHECK_LIB($cl_cv_polylib,PolyhedronTSort,[
			LIBS="-l$cl_cv_polylib $LIBS"
			break
		])
		cl_cv_polylib=missing
	done
fi
if test "$cl_cv_polylib" = "missing"; then
	AC_MSG_RESULT(Can't find PolyLib.)
else
	AC_CHECK_HEADER(polylib/$cl_cv_polylib.h,[],
		        [AC_MSG_ERROR(Can't find PolyLib headers.)])
fi

AC_SUBST(GIT_INDEX)
if test -f $srcdir/.git/HEAD; then
	GIT_INDEX="\$(top_srcdir)/.git/index"
fi


dnl /**************************************************************************
dnl  *                     Where is Parma Polyhedral Library?                 *
dnl  **************************************************************************/

# The default back-end for polyhedral operations is PolyLib, unless
# -with-ppl is used.
AC_SUBST(POLYHEDRAL_BACKEND)
POLYHEDRAL_BACKEND=polylib

dnl Checking for PPL
AC_SUBST(PPL_CPPFLAGS)
AC_SUBST(PPL_LDFLAGS)
AC_SUBST(PPL_LIBS)

AC_MSG_CHECKING(for Parma Polyhedral Library (PPL))
if test "x$with_ppl" != "x" -a "x$with_ppl" != "xno"; then

   	if test "x$with_polylib_prefix" != "x" -o "x$with_polylib_exec_prefix" != "x" -o "x$with_polylib_builddir" != "x"; then
	   	AC_MSG_ERROR([--with-polylib and --with-ppl are mutually exclusive])
	fi

	if test "x$with_ppl" != "xyes" ; then
	    ppl_prefix=$with_ppl
	    AC_MSG_RESULT(installed in $ppl_prefix)
	    CPPFLAGS="-I$ppl_prefix/include $CPPFLAGS"
	    LDFLAGS="-L$ppl_prefix/lib $LDFLAGS"
	    LD_LIBRARY_PATH="$ppl_prefix/lib:$LD_LIBRARY_PATH"
	    export LD_LIBRARY_PATH
	else
	    AC_MSG_RESULT(installed at system prefix)
	fi
	POLYHEDRAL_BACKEND=ppl

	AC_CHECK_HEADER(ppl_c.h,[],
		        [AC_MSG_ERROR(Can't find PPL headers.)])

	AC_MSG_CHECKING([for version 0.10 (or later revision) of PPL])
	AC_TRY_COMPILE([#include "ppl_c.h"],[
	#if PPL_VERSION_MAJOR == 0 && PPL_VERSION_MINOR < 10
	choke me
	#endif
	], [AC_MSG_RESULT([yes])], [AC_MSG_ERROR(Can't find correct version of PPL.) ])

	LIBS="$LIBS -lppl_c -lppl -lgmpxx"
	AC_DEFINE([CLOOG_PPL_BACKEND], 1, [Use the PPL backend])

else
	AC_MSG_RESULT(not using PPL)
fi

# Allow host libstdc++ to be specified for static linking with PPL.
AC_ARG_WITH(host-libstdcxx, [  --with-host-libstdcxx=L Use linker arguments L to link with libstdc++
                          when linking with PPL])

case $with_host_libstdcxx in
  no|yes) AC_MSG_ERROR([-with-host-libstdcxx needs an argument]) ;;
  *) LIBS="$LIBS $with_host_libstdcxx" ;;
esac

dnl /**************************************************************************
dnl  *                            Substitutions                               *
dnl  **************************************************************************/


dnl Substitutions to do in Makefile.in.
AC_SUBST(target)
AC_SUBST(CC)
AC_SUBST(LN_S)
AC_SUBST(SHEXT)

AC_SUBST(prefix)
AC_SUBST(exec_prefix)
AC_SUBST(INSTALL)

AC_SUBST(BITS)

AC_CONFIG_FILES(Makefile test/Makefile include/cloog/cloog.h autoconf/Doxyfile)
AC_CONFIG_FILES(doc/Makefile)
dnl AC_CONFIG_FILES(include/cloog/version.h)
AC_CONFIG_FILES([genversion.sh], [chmod +x genversion.sh])
AC_CONFIG_COMMANDS([version.h],
	[echo '#define CLOOG_HEAD "'`./genversion.sh`'"' > version.h])
AC_CONFIG_COMMANDS([doc/gitversion.texi],
	[echo '@set VERSION '`./genversion.sh`'' > doc/gitversion.texi])
AC_OUTPUT


echo "             /*-----------------------------------------------*"
echo "              *           CLooG configuration is OK           *"
echo "              *-----------------------------------------------*/"
echo "It appears that your system is OK to start CLooG compilation. You need"
echo "now to type \"make\". After compilation, you should check CLooG by typing"
echo "\"make check\". If no problem occur, you can type \"make uninstall\" if"
echo "you are upgrading an old version. Lastly type \"make install\" to install"
echo "CLooG on your system (log as root if necessary)."
