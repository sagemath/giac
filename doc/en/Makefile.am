info_TEXINFOS = giac_us.texinfo
EXTRA_DIST = xcasmenu xcasex keywords ifpdf.sty cas.tex slipres.tex pgcd.cc pgcd.cpp .gdbinit giacxcas.tex casinter.tex cascmd_en.tex html_* 
html:	giac_us.texinfo casinter.tex cascmd_en.tex
	touch giac_us.html && makeinfo --html --no-split giac_us.texinfo
	latex casinter
	makeindex casinter.idx
	latex casinter.tex
	latex2html -no_math -html_version 3.2,math casinter.tex --address '<A HREF="http://www-fourier.ujf-grenoble.fr/~parisse/giac.html">giac</A> documentation written by Bernard Parisse'
	latex cascmd_en.tex
	makeindex cascmd_en.idx
	latex cascmd_en.tex
	latex2html -local_icons -no_math -html_version 3.2,math cascmd_en.tex --address 'Documentation of <A HREF="http://www-fourier.ujf-grenoble.fr/~parisse/giac.html">giac</A> written by Ren&eacute;e De Graeve and Bernard Parisse'
dvi: cas.tex # giacxcas.tex
	latex cas
#	latex giacxcas
install-exec-hook: html
	../../mkinstalldirs $(prefix)/share/giac/doc/en
	cp -f *.html $(prefix)/share/giac/doc/en
	../../mkinstalldirs $(prefix)/share/giac/doc/en
	../../mkinstalldirs $(prefix)/share/giac/doc/en/casinter
	../../mkinstalldirs $(prefix)/share/giac/doc/en/cascmd_en
	rm -f casinter/*.old
	cp -f casinter/* $(prefix)/share/giac/doc/en/casinter
	rm -f cascmd_en/*.old
	cp -f cascmd_en/* $(prefix)/share/giac/doc/en/cascmd_en
	cp -f html_* $(prefix)/share/giac/doc/en/
	ln -sf ../fr/castor $(prefix)/share/giac/doc/en/castor
	ln -sf ../fr/casexo $(prefix)/share/giac/doc/en/casexo
	ln -sf ../fr/casrouge $(prefix)/share/giac/doc/en/casrouge
	ln -sf ../fr/casgeo $(prefix)/share/giac/doc/en/casgeo
	ln -sf ../fr/cassim $(prefix)/share/giac/doc/en/cassim
	ln -sf ../fr/casexo $(prefix)/share/giac/doc/en/casexo
	ln -sf ../fr/cascas $(prefix)/share/giac/doc/en/cascas
	ln -sf ../fr/algo.html $(prefix)/share/giac/doc/en/
	cp -f xcasmenu xcasex keywords $(prefix)/share/giac/doc/en/
