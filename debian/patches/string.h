Index: giac/src/hist.cxx
===================================================================
--- giac.orig/src/hist.cxx	2009-11-24 21:58:08.000000000 +0100
+++ giac/src/hist.cxx	2009-11-24 22:00:00.000000000 +0100
@@ -9,6 +9,9 @@
 #ifdef HAVE_MALLOC_H
 #include <malloc.h>
 #endif
+#ifdef HAVE_STRING_H
+#include <string.h>
+#endif
 static char ** xcas_argv; 
 static int xcas_argc,xcas_user_level; 
 static giac::vecteur rpn_menu; 
