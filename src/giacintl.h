#ifndef _GIACINTL_H
#define _GIACINTL_H

#include "first.h"
#include "config.h"

#if defined(__APPLE__) || defined(__FreeBSD__)
#include <libintl.h>
#endif

#ifdef HAVE_GETTEXT
#include <libintl.h>
#else

#ifndef _LIBINTL_H
#define _LIBINTL_H      1
inline const char * gettext(const char * s) { return s; };
#endif // _LIBINTL_H

#endif // HAVE_GETTEXT
#endif // _GIACINTL_H
