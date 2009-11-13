EXTRA_DIST = xcasmenu xcasex keywords ifpdf.sty castor.tex casrouge.tex cascmd_fr.tex casgeo.tex casexo.tex cascas.tex sy.eps rot.eps syrot.eps codage.eps cassim.tex permutation.tex tutoriel.tex windowinter.* demarr1.eps demarr1.png demarr1.fig demarr.eps figtortue/* figgeo/* casinter.tex algo.tex algo.html fltkcol.tex figure? html_* mkindex th1968* windowinter* undemix* cesaro* morleyp* baton* camembert* plotlist* nuage* nuee* quadrill* etoile* inscrit* sphere* nico* suite*  trialea* factalea*
html:	tutoriel.tex casrouge.tex cascmd_fr.tex casgeo.tex cassim.tex castor.tex casinter.tex algo.html
	latex casrouge.tex
	makeindex casrouge.idx
	latex casrouge.tex
	latex2html -local_icons -no_math -html_version 3.2,math casrouge.tex --address 'Documentation de <A HREF="http://www-fourier.ujf-grenoble.fr/~parisse/giac_fr.html">giac</A> &eacute;crite par Ren&eacute;e De Graeve'
	latex castor.tex
	makeindex castor.idx
	latex castor.tex
	latex2html -local_icons -no_math -html_version 3.2,math castor.tex --address 'Documentation de <A HREF="http://www-fourier.ujf-grenoble.fr/~parisse/giac_fr.html">giac</A> &eacute;crite par Ren&eacute;e De Graeve'
	latex cascmd_fr.tex
	makeindex cascmd_fr.idx
	latex cascmd_fr.tex
	latex2html -local_icons -no_math -html_version 3.2,math cascmd_fr.tex --address 'Documentation de <A HREF="http://www-fourier.ujf-grenoble.fr/~parisse/giac_fr.html">giac</A> &eacute;crite par Ren&eacute;e De Graeve'
	latex casgeo.tex
	makeindex casgeo.idx
	latex casgeo.tex
	latex2html -local_icons -no_math -html_version 3.2,math casgeo.tex --address 'Documentation de <A HREF="http://www-fourier.ujf-grenoble.fr/~parisse/giac_fr.html">giac</A> &eacute;crite par Ren&eacute;e De Graeve'
	latex cassim.tex
	makeindex cassim.idx
	latex cassim.tex
	latex2html -local_icons -no_math -html_version 3.2,math cassim.tex --address 'Documentation de <A HREF="http://www-fourier.ujf-grenoble.fr/~parisse/giac_fr.html">giac</A> &eacute;crite par Ren&eacute;e De Graeve'
	latex tutoriel.tex
	makeindex tutoriel.idx
	latex tutoriel.tex
	latex2html -local_icons -no_math -html_version 3.2,math tutoriel.tex --address 'Documentation de <A HREF="http://www-fourier.ujf-grenoble.fr/~parisse/giac_fr.html">giac</A> &eacute;crite par Ren&eacute;e De Graeve, Bernard Parisse et Bernard Ycart'
	latex casinter.tex
	latex casinter.tex
	latex2html -local_icons -no_math -html_version 3.2,math casinter.tex --address 'Documentation de <A HREF="http://www-fourier.ujf-grenoble.fr/~parisse/giac_fr.html">giac</A> &eacute;crite par Ren&eacute;e De Graeve'
	latex casexo.tex
	latex casexo.tex
	latex2html -local_icons -no_math -html_version 3.2,math casexo.tex --address 'Documentation de <A HREF="http://www-fourier.ujf-grenoble.fr/~parisse/giac_fr.html">giac</A> &eacute;crite par Ren&eacute;e De Graeve'
	latex cascas.tex
	latex cascas.tex
	latex2html -local_icons -no_math -html_version 3.2,math cascas.tex --address 'Documentation de <A HREF="http://www-fourier.ujf-grenoble.fr/~parisse/giac_fr.html">giac</A> &eacute;crite par Ren&eacute;e De Graeve'
dvi: casrouge.tex cascmd_fr.tex casgeo.tex permutation.tex castor.tex
	latex casrouge.tex
	latex cascmd_fr.tex
	latex casgeo.tex
	latex castor.tex
	latex permutation
install-exec-hook: html 
	../../mkinstalldirs $(prefix)/share/giac/doc/fr
	../../mkinstalldirs $(prefix)/share/giac/doc/fr/tutoriel
	../../mkinstalldirs $(prefix)/share/giac/doc/fr/casinter
	../../mkinstalldirs $(prefix)/share/giac/doc/fr/casrouge
	../../mkinstalldirs $(prefix)/share/giac/doc/fr/casexo
	../../mkinstalldirs $(prefix)/share/giac/doc/fr/cascas
	../../mkinstalldirs $(prefix)/share/giac/doc/fr/casgeo
	../../mkinstalldirs $(prefix)/share/giac/doc/fr/cascmd_fr
	../../mkinstalldirs $(prefix)/share/giac/doc/fr/cassim
	../../mkinstalldirs $(prefix)/share/giac/doc/fr/castor
	rm -f tutoriel/*.old
	rm -f casinter/*.old
	rm -f casrouge/*.old
	rm -f casgeo/*.old
	rm -f cascmd_fr/*.old
	rm -f cassim/*.old
	rm -f castor/*.old
	cp -f tutoriel/* $(prefix)/share/giac/doc/fr/tutoriel
	cp -f casinter/* $(prefix)/share/giac/doc/fr/casinter
	cp -f casrouge/* $(prefix)/share/giac/doc/fr/casrouge
	cp -f casgeo/* $(prefix)/share/giac/doc/fr/casgeo
	cp -f cascmd_fr/* $(prefix)/share/giac/doc/fr/cascmd_fr
	cp -f cassim/* $(prefix)/share/giac/doc/fr/cassim
	cp -f castor/* $(prefix)/share/giac/doc/fr/castor
	cp -f casexo/* $(prefix)/share/giac/doc/fr/casexo
	cp -f cascas/* $(prefix)/share/giac/doc/fr/cascas
	cp -f algo.html $(prefix)/share/giac/doc/fr/
	cp -f html_* $(prefix)/share/giac/doc/fr/
	cp -f xcasmenu xcasex keywords $(prefix)/share/giac/doc/fr/
#	../../mkinstalldirs /usr/
#	../../mkinstalldirs /usr/local
#	../../mkinstalldirs /usr/local/share
#	rm -rf /usr/local/share/giac
#	ln -s $(prefix)/share/giac /usr/local/share/giac