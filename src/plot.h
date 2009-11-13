/* -*- mode:C++ ; compile-command: "g++ -I.. -g -c plot.cc" -*- */
#ifndef _GIAC_PLOT_H
#define _GIAC_PLOT_H
/*
 *  Copyright (C) 2000 B. Parisse, Institut Fourier, 38402 St Martin d'Heres
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
#include "first.h"
#include <stdexcept>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include "gen.h"
#include "plot3d.h"
#ifdef HAVE_SIGNAL_H
#include <signal.h>
#endif


#ifdef HAVE_LIBFLTK
#include <FL/Enumerations.H>
#else
enum Fl_Color {	// standard colors
  FL_BLACK		= 0,
  FL_RED		= 1,
  FL_GREEN		= 2,
  FL_YELLOW		= 3,
  FL_BLUE		= 4,
  FL_MAGENTA		= 5,
  FL_CYAN		= 6,
  FL_WHITE		= 7,
  FL_INACTIVE_COLOR	= 8,
  FL_SELECTION_COLOR	= 15,

  FL_FREE_COLOR		= 16,
  FL_NUM_FREE_COLOR	= 16,

  FL_GRAY_RAMP		= 32,

  // boxtypes limit themselves to these colors so whole ramp is not allocated:
  FL_GRAY0		= 32,	// 'A'
  FL_DARK3		= 39,	// 'H'
  FL_DARK2		= 45,   // 'N'
  FL_DARK1		= 47,	// 'P'
  FL_GRAY		= 49,	// 'R' default color
  FL_LIGHT1		= 50,	// 'S'
  FL_LIGHT2		= 52,	// 'U'
  FL_LIGHT3		= 54,	// 'W'

  FL_COLOR_CUBE		= 56
};
#endif

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  std::string print_DOUBLE_(double d,unsigned ndigits);

  void local_sto_double(double value,const identificateur & i,GIAC_CONTEXT);
  void local_sto_double_increment(double value,const identificateur & i,GIAC_CONTEXT);
  vecteur quote_eval(const vecteur & v,const vecteur & quoted,GIAC_CONTEXT);

  // if v is the vector argument of pnt, get_style returns attributs and
  // set legende to the legende to be displayed
  vecteur get_style(const vecteur & v,std::string & legende);
  gen readvar(const gen & g);
  void read_tmintmaxtstep(vecteur & vargs,gen & t,int vstart,double &tmin,double & tmax,double &tstep,bool & tminmax_defined,bool & tstep_defined,GIAC_CONTEXT);
  int read_attributs(const vecteur & v,vecteur & attributs,GIAC_CONTEXT);
  void read_option(const vecteur & v,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,vecteur & attributs, int & nstep,int & jstep,int & kstep,GIAC_CONTEXT);
  gen curve_surface_apply(const gen & elem,const gen & b,gen (* func) (const gen &, const gen &,const context *),GIAC_CONTEXT);
  gen apply3d(const gen & e1, const gen & e2,const context * contextptr, gen (* f) (const gen &, const gen &,const context *) );

  gen mkrand2d3d(int dim,int nargs,gen (* f)(const gen &,const context *),GIAC_CONTEXT);
  gen droite_by_equation(const vecteur & v,bool est_plan,GIAC_CONTEXT);
  gen abs_norm(const gen & g,GIAC_CONTEXT);
  gen abs_norm2(const gen & g,GIAC_CONTEXT);
  gen dotvecteur(const gen & a,const gen & b,GIAC_CONTEXT);
  bool check3dpoint(const gen & g);
  vecteur remove_not_in_segment(const gen & a,const gen & b,int subtype,const vecteur & v,GIAC_CONTEXT);
  vecteur interpolygone(const vecteur & p,const gen & bb,GIAC_CONTEXT);

  bool set_turtle_state(const vecteur & v,GIAC_CONTEXT);
  gen turtle2gen(const logo_turtle & turtle);
  vecteur turtlevect2vecteur(const std::vector<logo_turtle> & v);
  std::vector<logo_turtle> vecteur2turtlevect(const vecteur & v);
  logo_turtle vecteur2turtle(const vecteur & v);

  // File "Out.txt" for gnuwince
#ifdef GNUWINCE
  extern std::ofstream * outptr;
#endif

  vecteur plotpreprocess(const gen & args,GIAC_CONTEXT);
  vecteur gen2vecteur(const gen & arg);
  bool chk_double_interval(const gen & g,double & inf,double & sup,GIAC_CONTEXT);
  bool readrange(const gen & g,double defaultxmin,double defaultxmax,gen & x, double & xmin, double & xmax,GIAC_CONTEXT);
  void ck_parameter_x(GIAC_CONTEXT);
  void ck_parameter_y(GIAC_CONTEXT);
  void ck_parameter_z(GIAC_CONTEXT);
  void ck_parameter_t(GIAC_CONTEXT);
  void ck_parameter_u(GIAC_CONTEXT);
  void ck_parameter_v(GIAC_CONTEXT);
  void ck_parameter(const gen & ,GIAC_CONTEXT);
  extern int LEGENDE_SIZE;
  extern int COORD_SIZE;
  extern int PARAM_STEP;

  void autoname_plus_plus(std::string & autoname);
  int erase3d();
  int erase_pos(GIAC_CONTEXT);
  int erase_pos(int current,GIAC_CONTEXT);
  bool is_segment(const gen & e);
  gen remove_at_pnt(const gen & e);
  gen remove_sto(const gen & e);
  vecteur selection2vecteur(const std::vector<int> & selected,GIAC_CONTEXT);
  vecteur selection2vecteureval(const std::vector<int> & selected,GIAC_CONTEXT);
  // find best int in selected (and modify selected)
  // p is the pointed mouse point, eps the precision
  // try_perp=-1 if no try of perp line, =an history_position otherwise
  bool find_best(std::vector<int> & selected,const gen & p,double eps,int try_perp_history_pos,int & pnt_pos,int & history_position,gen & res,GIAC_CONTEXT);
  int findfirstcercle(const vecteur & v);
  int findfirstpoint(const vecteur & v);

  extern bool run_modif;  // obsolete, used by make_child and co
  extern int run_modif_pos;
  extern bool fastcurveprint;
  void rewrite_with_t_real(gen & eq,const gen & t,GIAC_CONTEXT);
  extern int _GROUP__VECT_subtype[];
  extern unary_function_ptr plot_sommets[];
  extern unary_function_ptr not_point_sommets[];
  extern unary_function_ptr notexprint_plot_sommets[];
  extern pid_t gnuplot_pid;
  int set_nonblock_flag (int desc, int value);
  extern std::string gnuplot_name; // name of the program gnuplot
  extern std::string gnuplot_filename; // name of files where we save plots
  extern int gnuplot_fileno; // current index in save plot files
  extern int gnuplot_pixels_per_eval;
  extern double gnuplot_xmin,gnuplot_xmax,gnuplot_ymin,gnuplot_ymax,gnuplot_zmin,gnuplot_zmax,gnuplot_tmin,gnuplot_tmax,gnuplot_tstep,global_window_xmin,global_window_xmax,global_window_ymin,global_window_ymax,x_tick,y_tick; // ranges
  extern double class_minimum,class_size; // histogram
  extern bool autoscale;
  extern bool has_gnuplot;
  extern bool win9x; // True for windows (win95/98/Me can not pipe to gnuplot)
  // runs gnuplot if necessary, returns the FD of the pipe to write to gnuplot
  int run_gnuplot(int & r);
  void gnuplot_wait(int handle,FILE * gnuplot_out_readstream,int ngwait=0);
  FILE * open_gnuplot(bool & clrplot,FILE * & gnuplot_out_readstream,int & r);
  // Plot a function or a set of functions
  void reset_gnuplot_hidden3d(FILE * stream);
  extern bool gnuplot_do_splot;
  std::string gnuplot_traduit(const gen & g);
  void win9x_gnuplot(FILE * stream);
  extern int gnuplot_wait_times;

  // return parametrization for a parametric curve and translate
  // ellipsis/hyperbola to a rational parametrization
  // m will contain the complex depending on gen_t 
  bool find_curve_parametrization(const gen & geo_obj,gen & m,const gen & gen_t,double T,gen & tmin,gen & tmax,gen & tstep,GIAC_CONTEXT);

  gen plotfunc(const gen & f,const gen & vars,const vecteur & attributs,bool clrplot,double function_xmin,double function_xmax,double function_ymin,double function_ymax,double function_zmin, double function_zmax,int nstep,int jstep,bool showeq,GIAC_CONTEXT);
  // return a vector of values with simple decimal representation
  // between xmin/xmax or including xmin/xmax (if bounds is true)
  vecteur ticks(double xmin,double xmax,bool bounds);
  gen plotcontour(const gen & f0,bool contour,GIAC_CONTEXT);
  gen plot_array(const std::vector< std::vector< double> > & fij,int imax,int jmax,double xmin,double xmax,double dx,double ymin,double ymax,double dy,const vecteur & lz,const vecteur & attributs,bool contour,GIAC_CONTEXT);
  symbolic symb_plotfunc(const gen & a,const gen & b);
  symbolic symb_plotfunc(const gen & a);
  bool latex_replot(FILE * stream,const std::string & s);
  bool png_replot(int i);
  bool png_replot(const std::string & s);
  bool terminal_replot(const char * terminal,int i,const char * file_extension);
  bool terminal_replot(const char * terminal,const std::string & s);
  extern const std::string _plotfunc_s;
  gen _aire(const gen & args,GIAC_CONTEXT);
  gen funcplotfunc(const gen & args,bool densityplot,const context * contextptr);
  gen _plotfunc(const gen &,GIAC_CONTEXT);
  extern unary_function_ptr at_plotfunc;
  extern unary_function_ptr at_funcplot;
  extern const std::string _plot_s;
  // gen _plot(const gen &);
  extern unary_function_ptr at_plot;
  void kill_gnuplot();
  gen remove_at_pnt(const gen & e);

  symbolic symb_erase(const gen & args);
  extern const std::string _erase_s;
  gen _erase(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_erase;

  extern const std::string _pixon_s;
  gen _pixon(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_pixon;

  extern const std::string _pixoff_s;
  gen _pixoff(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_pixoff;

  extern const std::string _droite_s;
  gen _droite(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_droite;

  extern const std::string _demi_droite_s;
  gen _demi_droite(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_demi_droite;

  extern const std::string _segment_s;
  gen _segment(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_segment;

  extern unary_function_ptr at_inter_unique;
  extern unary_function_ptr at_polygone_ouvert;

  // segment x->y with attributs c
  gen symb_segment(const gen & x,const gen & y,const vecteur & ,int ,GIAC_CONTEXT);
  // point x color c name nom
  gen symb_pnt_name(const gen & x,const gen & c,const gen & nom,GIAC_CONTEXT);
  // point x and color c
  gen symb_pnt(const gen & x,const gen & c,GIAC_CONTEXT);
  gen pnt_attrib(const gen & point,const vecteur & attributs,GIAC_CONTEXT);
  // point x with default color FL_BLACK
  gen symb_pnt(const gen & x,GIAC_CONTEXT);
  extern const std::string _pnt_s;
  gen _pnt(const gen & args);
  extern unary_function_ptr at_pnt;
  extern unary_function_ptr at_animation;
  int animations(const gen & g); // number of animations inside g
  gen get_animation_pnt(const gen & g,int pos);

  symbolic symb_point(const gen & x,const gen & y);
  extern const std::string _point_s;
  gen _point(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_point;

  extern const std::string _affixe_s;
  gen _affixe(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_affixe;

  gen _abscisse(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_abscisse;

  gen _ordonnee(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_ordonnee;

  extern const std::string _cercle_s;
  gen _cercle(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_cercle;
  bool centre_rayon(const gen & cercle,gen & centre,gen & rayon,bool absrayon, GIAC_CONTEXT);
  extern const std::string _centre_s;
  gen _centre(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_centre;
  extern const std::string _rayon_s;
  gen _rayon(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_rayon;

  symbolic symb_milieu(const gen & x,const gen & y);
  symbolic symb_milieu(const gen & x);
  extern const std::string _milieu_s;
  gen _milieu(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_milieu;

  extern const std::string _mediatrice_s;
  gen _mediatrice(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_mediatrice;

  gen bissectrice(const gen & args,bool interieur,GIAC_CONTEXT);
  extern const std::string _bissectrice_s;
  gen _bissectrice(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_bissectrice;

  extern const std::string _exbissectrice_s;
  gen _exbissectrice(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_exbissectrice;

  extern const std::string _mediane_s;
  gen _mediane(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_mediane;

  extern const std::string _circonscrit_s;
  gen _circonscrit(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_circonscrit;

  extern const std::string _inscrit_s;
  gen _inscrit(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_inscrit;

  extern const std::string _exinscrit_s;
  gen _exinscrit(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_exinscrit;

  extern const std::string _isobarycentre_s;
  gen _isobarycentre(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_isobarycentre;

  symbolic symb_perpendiculaire(const gen & a);
  extern const std::string _perpendiculaire_s;
  gen _perpendiculaire(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_perpendiculaire;

  extern const std::string _parallele_s;
  gen _parallele(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_parallele;

  gen distance2pp(const gen & ee,const gen & ff,GIAC_CONTEXT);
  gen distance2(const gen & f1,const gen & f2,GIAC_CONTEXT);
  extern const std::string _longueur2_s;
  gen _longueur2(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_longueur2;

  gen longueur(const gen & f1,const gen & f2,GIAC_CONTEXT);
  extern const std::string _longueur_s;
  gen _longueur(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_longueur;

  gen angle(const gen & f1,const gen & f2,GIAC_CONTEXT);
  extern const std::string _angle_s;
  gen _angle(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_angle;

  gen scalar_product(const gen & a,const gen & b,GIAC_CONTEXT);
  // return t such that ta+(1-t)b is the projection of c on [a,b] 
  gen projection(const gen & a,const gen & b,const gen & c,GIAC_CONTEXT);
  // projection of p on a parametric curve
  // e=symb_cercle or line, returns t
  gen projection(const gen & e,const gen & p,GIAC_CONTEXT);
  gen parameter2point(const vecteur & v,GIAC_CONTEXT);
  
  std::vector<int> nearest_point(const vecteur & v,const gen & p,double eps,GIAC_CONTEXT);

  vecteur inter(const gen & a,const gen & b,GIAC_CONTEXT);

  extern const std::string _signal_s;
  gen _signal(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_signal;

  extern const std::string _click_s;
  gen _click(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_click;
  class unary_function_eval;
  extern unary_function_eval __click;

  extern const std::string _element_s;
  gen _element(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_element;

  extern const std::string _as_function_of_s;
  gen _as_function_of(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_as_function_of;

  extern const std::string _lieu_s;
  gen _lieu(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_lieu;

  extern const std::string _head_s;
  gen _head(const gen & args);
  extern unary_function_ptr at_head;

  extern const std::string _tail_s;
  gen _tail(const gen & args);
  extern unary_function_ptr at_tail;

  extern const std::string _sommets_s;
  gen _sommets(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_sommets;

  extern const std::string _symetrie_s;
  gen _symetrie(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_symetrie;

  extern const std::string _rotation_s;
  gen _rotation(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_rotation;

  extern const std::string _projection_s;
  gen _projection(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_projection;

  extern const std::string _homothetie_s;
  gen _homothetie(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_homothetie;

  extern const std::string _est_aligne_s;
  gen _est_aligne(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_est_aligne;

  extern const std::string _est_cocyclique_s;
  gen _est_cocyclique(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_est_cocyclique;

  extern const std::string _est_parallele_s;
  gen _est_parallele(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_est_parallele;

  extern const std::string _est_perpendiculaire_s;
  gen _est_perpendiculaire(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_est_perpendiculaire;

  gen _est_element(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_est_element;

  extern const std::string _inversion_s;
  gen _inversion(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_inversion;

  extern const std::string _similitude_s;
  gen _similitude(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_similitude;

  extern const std::string _translation_s;
  gen translation(const gen & a,const gen & bb,GIAC_CONTEXT);
  gen _translation(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_translation;

  symbolic symb_curve(const gen & source,const gen & plot);
  extern const std::string _curve_s;
  gen _curve(const gen & args);
  extern unary_function_ptr at_curve;

  extern const std::string _plotparam_s;
  gen _plotparam(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_plotparam;
  extern unary_function_ptr at_paramplot;
  gen paramplotparam(const gen & args,bool clrplot,const context * contextptr);

  extern const std::string _plotpolar_s;
  gen _plotpolar(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_plotpolar;
  extern unary_function_ptr at_polarplot;

  extern const std::string _parameq_s;
  gen _parameq(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_parameq;

  extern const std::string _equation_s;
  gen _equation(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_equation;

  extern const std::string _tangent_s;
  gen _tangent(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_tangent;

  extern const std::string _ellipse_s;
  gen _ellipse(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_ellipse;

  extern const std::string _hyperbole_s;
  gen _hyperbole(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_hyperbole;

  extern const std::string _parabole_s;
  gen _parabole(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_parabole;

  extern const std::string _legende_s;
  gen _legende(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_legende;

  extern const std::string _couleur_s;
  gen _couleur(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_couleur;
  extern unary_function_ptr at_display ;

  extern const std::string _parameter_s;
  gen _parameter(const gen & args);
  extern unary_function_ptr at_parameter;

  extern const std::string _hauteur_s;
  gen _hauteur(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_hauteur;

  extern const std::string _triangle_s;
  gen _triangle(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_triangle;

  extern const std::string _triangle_rectangle_s;
  gen _triangle_rectangle(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_triangle_rectangle;

  extern const std::string _triangle_isocele_s;
  gen _triangle_isocele(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_triangle_isocele;

  extern const std::string _triangle_equilateral_s;
  gen _triangle_equilateral(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_triangle_equilateral;

  extern const std::string _parallelogramme_s;
  gen _parallelogramme(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_parallelogramme;

  extern const std::string _carre_s;
  gen _carre(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_carre;

  extern const std::string _quadrilatere_s;
  gen _quadrilatere(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_quadrilatere;

  extern const std::string _rectangle_s;
  gen _rectangle(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_rectangle;

  extern const std::string _losange_s;
  gen _losange(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_losange;

  extern const std::string _polygone_s;
  gen _polygone(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_polygone;

  extern const std::string _plotfield_s;
  gen _plotfield(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_plotfield;
  extern unary_function_ptr at_fieldplot;

  extern const std::string _interactive_plotode_s;
  gen _interactive_plotode(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_interactive_plotode;
  extern unary_function_ptr at_interactive_odeplot;

  extern const std::string _plotode_s;
  gen _plotode(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_plotode;
  extern unary_function_ptr at_odeplot;

  std::ostream & archive(std::ostream & os,const gen & e,GIAC_CONTEXT);
  gen unarchive(std::istream & is,const vecteur & l,GIAC_CONTEXT);
  gen archive_session(bool save_history,std::ostream & os,GIAC_CONTEXT);
  gen archive_session(bool save_history,const std::string & s,GIAC_CONTEXT);
  std::string archive_session(bool save_history,GIAC_CONTEXT);
  gen unarchive_session(std::istream & is,int level, const gen & replace,GIAC_CONTEXT);
  gen unarchive_session(const std::string & s,int level, const gen & replace,GIAC_CONTEXT);
  gen unarchive_session_string(const std::string & s,int level, const gen & replace,GIAC_CONTEXT);

  extern const std::string _archive_s;
  gen _archive(bool save_history,const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_archive;

  extern const std::string _unarchive_s;
  gen _unarchive(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_unarchive;

  extern const std::string _xyztrange_s;
  void geo_setup(const vecteur & w,GIAC_CONTEXT);
  gen xyztrange(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,double tmin,double tmax,double wxmin,double wxmax,double wymin, double wymax, int axes,double class_minimum,double class_size,bool gnuplot_hidden3d);
  gen _xyztrange(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_xyztrange;
  class unary_function_unary;
  extern unary_function_eval __xyztrange;

  gen _switch_axes(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_switch_axes;

  gen plotseq(const gen& f,const gen&x,double x0,double xmin,double xmax,int niter,GIAC_CONTEXT);
  gen _plotseq(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_plotseq;
  extern unary_function_ptr at_seqplot;

  gen _plotimplicit(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_plotimplicit;
  extern unary_function_ptr at_implicitplot;

  gen _plotcontour(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_plotcontour;
  extern unary_function_ptr at_contourplot;

  gen _bitmap(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_bitmap;

  gen _Pictsize(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_Pictsize;

  extern unary_function_ptr implicittex_plot_sommets[];

  gen _plot_style(const gen & args);
  extern unary_function_ptr at_plot_style ;

  gen _DrawInv(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_DrawInv ;

  extern unary_function_ptr at_Graph;
  extern unary_function_ptr at_DrawFunc;
  extern unary_function_ptr at_DrawParm;
  extern unary_function_ptr at_DrawPol;
  extern unary_function_ptr at_DrwCtour;

  extern const std::string _est_isocele_s;
  gen _est_isocele(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_est_isocele;

  extern const std::string _est_equilateral_s;
  gen _est_equilateral(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_est_equilateral;

  extern const std::string _est_carre_s;
  gen _est_carre(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_est_carre;

  extern const std::string _est_losange_s;
  gen _est_losange(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_est_losange;

  extern const std::string _est_parallelogramme_s;
  gen _est_parallelogramme(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_est_parallelogramme;

  extern const std::string _est_rectangle_s;
  gen _est_rectangle(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_est_rectangle;
 
  extern const std::string _est_harmonique_s;
  gen _est_harmonique(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_est_harmonique;
  
  extern const std::string _div_harmonique_s;
  gen _div_harmonique(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_div_harmonique;
  
  extern const std::string _point_div_s;
  gen _point_div(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_point_div;
  
  extern const std::string _birapport_s;
  gen _birapport(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_birapport;

  extern const std::string _est_harmonique_s;
  gen _est_harmonique(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_est_harmonique;

  extern const std::string _div_harmonique_s;
  gen _div_harmonique(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_div_harmonique;

  extern const std::string _conj_harmonique_s;
  gen _conj_harmonique(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_conj_harmonique;

  extern const std::string _conj_harmoniques_s;
  gen _conj_harmoniques(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_conj_harmoniques;

  extern const std::string _point_div_s;
  gen _point_div(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_point_div;

  extern const std::string _birapport_s;
  gen _birapport(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_birapport;

  extern const std::string _puissance_s;
  gen _puissance(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_puissance;

  extern const std::string _axe_radical_s;
  gen _axe_radical(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_axe_radical;

  extern const std::string _polaire_s;
  gen _polaire(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_polaire;

  extern const std::string _pole_s;
  gen _pole(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_pole;

  extern const std::string _polaire_reciproque_s;
  gen _polaire_reciproque(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_polaire_reciproque;

  extern const std::string _est_orthogonal_s;
  gen _est_orthogonal(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_est_orthogonal;

  extern const std::string _est_conjugue_s;
  gen _est_conjugue(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_est_conjugue;

  extern const std::string _est_faisceau_cercle_s;
  gen _est_faisceau_cercle(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_est_faisceau_cercle;

  extern const std::string _est_faisceau_droite_s;
  gen _est_faisceau_droite(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_est_faisceau_droite;

  extern const std::string _enveloppe_s;
  gen _enveloppe(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_enveloppe;

  extern bool gnuplot_hidden3d,gnuplot_pm3d;
  void gnuplot_set_hidden3d(bool hidden);
  void gnuplot_set_pm3d(bool b);

  int gnuplot_show_pnt(const symbolic & e,GIAC_CONTEXT);
  int graph_output_type(const giac::gen & g);
  gen put_attributs(const gen & lieu_geo,const vecteur & attributs,GIAC_CONTEXT);
  vecteur seq2vecteur(const gen & g);

  int est_aligne(const gen & a,const gen & b,const gen & c,GIAC_CONTEXT);
  bool est_coplanaire(const gen & a,const gen & b,const gen & c,const gen & d,GIAC_CONTEXT);
  bool est_cocyclique(const gen & a,const gen & b,const gen & c,const gen & d,GIAC_CONTEXT);
  // True if a=coeff*b
  bool est_parallele_vecteur(const vecteur & a,const vecteur &b,gen & coeff,GIAC_CONTEXT);
  bool est_parallele_vecteur(const vecteur & a,const vecteur &b,GIAC_CONTEXT);
  bool est_parallele(const gen & a,const gen & b,GIAC_CONTEXT);
  bool est_perpendiculaire(const gen & a,const gen & b,GIAC_CONTEXT);
  // check if a belongs to b, a must be a complex, b a line or circle or curve
  int est_element(const gen & a_orig,const gen & b_orig,GIAC_CONTEXT);
  bool est_carre(const gen & a,const gen & b,const gen & c,const gen & d,GIAC_CONTEXT);
  int est_isocele(const gen & a,const gen & b,const gen & c,GIAC_CONTEXT);
  bool est_equilateral(const gen & a,const gen & b,const gen & c,GIAC_CONTEXT);
  int est_trianglerect(const gen & a,const gen & b,const gen & c,GIAC_CONTEXT);
  //teste si deux cercles C1 centre c1 rayon R1 et C2  centre c2 rayon R2
  //sont orthogonaux
  bool est_orthogonal(const gen & c1,const gen & R1,const gen & c2,const gen & R2,GIAC_CONTEXT);
  //teste si 4 points forment une division harmonique
  bool est_harmonique(const gen & a,const gen & b,const gen & c,const gen & d,GIAC_CONTEXT);

  gen _vector(const gen & args,GIAC_CONTEXT);
  extern unary_function_ptr at_vector;

  extern unary_function_ptr point_sommet_tab_op[];

  /* following obsolete declaration that will be remove in the future
   when all fork/child etc. will be removed */
  extern vecteur plot_instructions;
  // extern int plot_instructionsh,plot_instructionsw;
  extern std::string PICTautoname;
  void PICTautoname_plus_plus();

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC

#endif // _GIAC_PLOT_H
