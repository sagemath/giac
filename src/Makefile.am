## Process this file with automake to produce Makefile.in

lib_LTLIBRARIES = libgiac.la libxcas.la
libgiac_la_SOURCES = sym2poly.cc gausspol.cc threaded.cc \
  moyal.cc maple.cc ti89.cc mathml.cc\
  misc.cc permu.cc quater.cc desolve.cc input_parser.yy \
  symbolic.cc index.cc modpoly.cc modfactor.cc \
  ezgcd.cc derive.cc solve.cc intg.cc intgab.cc risch.cc lin.cc series.cc \
  subst.cc vecteur.cc csturm.cc tex.cc global.cc \
  ifactor.cc alg_ext.cc gauss.cc isom.cc plot.cc plot3d.cc help.cc rpn.cc prog.cc \
  pari.cc cocoa.cc unary.cc usual.cc identificateur.cc gen.cc input_lexer.ll first.cc \
  TmpLESystemSolver.C TmpFGLM.C
## does not work libgiac_la_LIBADD = @LIBOBJS@
giacincludedir = $(includedir)/giac
giacinclude_HEADERS = dispatch.h fraction.h gen.h desolve.h misc.h ti89.h \
  mathml.h \
  identificateur.h unary.h symbolic.h index.h modpoly.h modfactor.h permu.h \
  monomial.h poly.h gausspol.h threaded.h sym2poly.h ezgcd.h derive.h \
  solve.h quater.h \
  subst.h usual.h intg.h intgab.h risch.h lin.h series.h \
  vecteur.h csturm.h moyal.h tex.h global.h \
  ifactor.h alg_ext.h gauss.h isom.h plot.h plot3d.h help.h rpn.h prog.h \
  pari.h cocoa.h giac.h lexer.h\
  input_lexer.h input_parser.h first.h maple.h giacintl.h path.h libintl.h \
  TmpFGLM.H TmpLESystemSolver.H

libxcas_la_SOURCES = History.cc Input.cc Xcas1.cc Equation.cc Print.cc \
	Tableur.cc Editeur.cc Graph.cc Graph3d.cc Help1.cc Cfg.cc \
	Flv_CStyle.cc \
	Flve_Check_Button.cc Flve_Input.cc Flv_Style.cc Flv_Data_Source.cc \
	Flve_Combo.cc Flv_List.cc Flv_Table.cc gl2ps.c 
xcasincludedir = $(includedir)/xcas
xcasinclude_HEADERS = History.h Input.h Xcas1.h Equation.h Print.h Tableur.h \
	Editeur.h Graph.h Graph3d.h Help1.h Cfg.h \
	Flv_CStyle.H  Flve_Check_Button.H  Flve_Input.H  Flv_Style.H \
	Flv_Data_Source.H  Flve_Combo.H Flv_List.H Flv_Table.H gl2ps.h

bin_PROGRAMS = icas aide xcas 
icas_SOURCES = icas.cc 
icas_LDADD = -lxcas -lgiac -lgmp #-lcurses
icas_LDFLAGS = -L../src/.libs 
#mycas_SOURCES = mycas.cc 
#mycas_LDFLAGS = -s -L../src/.libs 
#mycas_LDADD = -lgiac -lgmp #-lcurses
aide_SOURCES = help.cc aide.cc
aide_LDADD = -lgiac -lgmp#-lcurses
aide_LDFLAGS = -L../src/.libs 
# If you want a static binary remove the hash below and comment -lcurses
#xcas_SOURCES = xcas.cc plotfltk.cc sheet.cc Eqw.cc Flv_CStyle.cc Flve_Check_Button.cc Flve_Input.cc Flv_Style.cc Flv_Data_Source.cc Flve_Combo.cc Flv_List.cc Flv_Table.cc #$(libgiac_la_SOURCES) 
#xcas_LDADD = -lgiac -lgmp # /usr/local/lib/libfltk.a  -lcurses
#xcas_LDFLAGS = -L../src/.libs  #-s --static 
xcas_SOURCES = xcas.cxx 
xcas_LDADD = -lxcas -lgiac -lgmp # /usr/local/lib/libfltk.a  -lcurses
xcas_LDFLAGS = -L../src/.libs  #-s --static 

INCLUDES = -I$(srcdir) -I$(top_srcdir)
#LFLAGS = -Pgiac_yy -olex.yy.c
YFLAGS = -p giac_yy -d
EXTRA_DIST = .gdbinit xgiac pgiac prime.table gen_user \
	input_parser.h softmath.h softmath.cc\
	mkxcas mkwcas mkgcas mkaout Makefile.wince\
	Makefile.rpm Makefile.ipaq Makefilenew.ipaq ipaq_usr.tgz \
	Makefile.fast Makefilenew.win Makefile.O3 Makefile.64 \
	normalize.cc integrate.cc partfrac.cc factor.cc xcas.cc \
	hist.fl hist.h hist.cxx  \
	progfr.cc renee.cc renee.h cas2.cc README.ipaq README.WIN   \
	readme.doc latex.bat xdvi.bat dvips.bat xcasfr.bat xcases.bat \
	xcasfr.sh xcases.sh xcasen.sh mkstat mkhist updatecas update.bat \
	giac.pl cas2html.cc all_global_var find_global_var.cc static_init \
	Fl_GDI_Printer.cxx  Fl_PS_Printer.cxx \
	casce.cc xcasce.cc xcasctrl.cc casctrl.cc config.h.vcc Makefile.vcc giac_oo.cpp giaclib.vcproj

install-exec-hook:
	../mkinstalldirs $(prefix)/share/giac
	cp -f updatecas $(bindir)
	cp -f xgiac $(bindir)
	cp -f pgiac $(bindir)
	ln -sf aide $(bindir)/cas_help
	ln -sf aide $(bindir)/en_cas_help
	ln -sf aide $(bindir)/es_cas_help
	ln -sf aide $(bindir)/fr_cas_help
	ln -sf icas $(bindir)/giac
	ln -sf xcas $(bindir)/xcasnew
	../mkinstalldirs $(prefix)/share/giac/src
	cp -f $(libgiac_la_SOURCES) $(prefix)/share/giac/src
	../mkinstalldirs $(includedir)
	../mkinstalldirs $(includedir)/giac
	cp -f ../config.h $(includedir)/giac
