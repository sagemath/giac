EXTRA_DIST = glace.xws soleil.xws terre.xws terre.jpg insol.xws vostok.xws  co2.xws energie.xws co2t150.xws edc.xws global_temp.xws

install-exec-hook:
	../../../mkinstalldirs $(prefix)/share/giac
	../../../mkinstalldirs $(prefix)/share/giac/examples
	../../../mkinstalldirs $(prefix)/share/giac/examples/Exemples
	../../../mkinstalldirs $(prefix)/share/giac/examples/Exemples/climat
	cp -f *.xws *.jpg $(prefix)/share/giac/examples/Exemples/climat
