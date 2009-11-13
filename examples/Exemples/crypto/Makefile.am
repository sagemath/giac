EXTRA_DIST = cesar.xws codage.xws rsa.xws reed_s.xws

install-exec-hook:
	../../../mkinstalldirs $(prefix)/share/giac
	../../../mkinstalldirs $(prefix)/share/giac/examples
	../../../mkinstalldirs $(prefix)/share/giac/examples/Exemples
	../../../mkinstalldirs $(prefix)/share/giac/examples/Exemples/crypto
	cp -f *.xws $(prefix)/share/giac/examples/Exemples/crypto
