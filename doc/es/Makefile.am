info_TEXINFOS = giac_es.texinfo
EXTRA_DIST = xcasmenu xcasex keywords cas_es.tex html_* 
html:	giac_es.texinfo
	touch giac_es.html && makeinfo --html --no-split giac_es.texinfo
dvi: cas_es.tex 
	latex cas_es
install-exec-hook: html
	../../mkinstalldirs $(prefix)/share/giac/doc/es
	cp -f *.html $(prefix)/share/giac/doc/es
	cp -f html_* $(prefix)/share/giac/doc/es/
	ln -sf ../en/casinter $(prefix)/share/giac/doc/es/casinter
	ln -sf ../en/cascmd_en $(prefix)/share/giac/doc/es/cascmd_es
	ln -sf ../fr/castor $(prefix)/share/giac/doc/es/castor
	ln -sf ../fr/casexo $(prefix)/share/giac/doc/es/casexo
	ln -sf ../fr/casrouge $(prefix)/share/giac/doc/es/casrouge
	ln -sf ../fr/casgeo $(prefix)/share/giac/doc/es/casgeo
	ln -sf ../fr/casexo $(prefix)/share/giac/doc/es/casexo
	ln -sf ../fr/cascas $(prefix)/share/giac/doc/es/cascas
	ln -sf ../fr/cassim $(prefix)/share/giac/doc/es/cassim
	ln -sf ../fr/algo.html $(prefix)/share/giac/doc/es/
	cp -f xcasmenu xcasex keywords $(prefix)/share/giac/doc/es/
