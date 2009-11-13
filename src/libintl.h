#include "config.h"
#ifdef __APPLE__
#include "/usr/local/include/libintl.h"
#else // APPLE

#ifdef __FreeBSD__
# define gettext(Msgid) (Msgid)
# define dgettext(Domainname, Msgid) (Msgid)
# define dcgettext(Domainname, Msgid, Category) (Msgid)
# define textdomain(Domainname) ((char *) Domainname)
# define bindtextdomain(Domainname, Dirname) ((char *) Dirname)
#else // __FreeBSD__

#ifdef HAVE_GETTEXT
#include "/usr/include/libintl.h"
#else

#ifndef _LIBINTL_H
#define _LIBINTL_H      1
inline const char * gettext(const char * s) { return s; };
#endif // libintl_h

#endif // have_gettext

#endif // __FreeBSD__

#endif //APPLE
