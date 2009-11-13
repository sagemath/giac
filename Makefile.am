PACKAGE = giac
VERSION = 0.8.4
SUBDIRS = 	src check intl po debianold debian doc examples 
DIST_SUBDIRS = 	src check intl po debianold debian doc examples 
EXTRA_DIST = autogen.sh setup inst acconfig.h TROUBLES config.h config.h.in Init.scm config.h.win config.h.ipaq config.h.gc config.h.osx compte INSTALACION LEAME.WIN LEAME.ipaq libtool giac.specs giac.specs32 fltk.1.1.changes README.MACOSX pariinl.h changeversion mkosx config.h.wince depcomp configure.in.debian mkdebian mkdebian2 win2unix.cc addu.cc xcas.nsi xcas.desktop xcas.xpm xcas.applications README.visualcpp

# Rule to build tar-bzipped distribution package
$(PACKAGE)-$(VERSION).tar.gz: dist

install-exec-hook:
	cp config.h $(includedir)/giac
	./mkinstalldirs $(prefix)/share/
	./mkinstalldirs $(prefix)/share/pixmaps
	cp xcas.xpm $(prefix)/share/pixmaps/xcas.xpm
	./mkinstalldirs $(prefix)/share/applications
	cp xcas.desktop $(prefix)/share/applications/xcas.desktop
	./mkinstalldirs $(prefix)/share/application-registry
	cp xcas.applications $(prefix)/share/application-registry/xcas.applications
