SUBDIRS = 	en es fr local
DIST_SUBDIRS = 	en es fr local
EXTRA_DIST = aide_cas  README index.html dxcas/* pari/*
install-exec-hook: 
	../mkinstalldirs $(prefix)/share/giac
	../mkinstalldirs $(prefix)/share/giac/doc
	../mkinstalldirs $(prefix)/share/giac/doc/dxcas
	../mkinstalldirs $(prefix)/share/giac/doc/pari
	cp -f aide_cas $(prefix)/share/giac
	cp -f README $(prefix)/share/giac/doc
	cp -f index.html $(prefix)/share/giac/doc
	cp -f dxcas/* $(prefix)/share/giac/doc/dxcas
	cp -f pari/* $(prefix)/share/giac/doc/pari
