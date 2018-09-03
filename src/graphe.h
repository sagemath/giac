/*
 * graphe.h
 *
 * (c) 2018 Luka Marohnić
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __GRAPHE_H
#define __GRAPHE_H
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "first.h"
#include "gen.h"
#include "unary.h"
#include "moyal.h"
#include <string>
#include <iostream>
#include <fstream>
#include <queue>
#include <stack>
#include <set>
#ifdef HAVE_LIBGLPK
#include <glpk.h>
#endif

#ifndef DBL_MAX
#define DBL_MAX 1.79769313486e+308
#endif
#define PLASTIC_NUMBER 1.32471795724
#define PLASTIC_NUMBER_SQUARED 1.75487766625
#define PLASTIC_NUMBER_CUBED 2.32471795724
#define MARGIN_FACTOR 0.139680581996 // pow(PLASTIC_NUMBER,-7)
#define NEGLIGIBILITY_FACTOR 0.0195106649868 // MARGIN_FACTOR^2

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

typedef unsigned long ulong;

enum gt_dot_token_type {
    _GT_DOT_TOKEN_TYPE_IDENTIFIER = 1,
    _GT_DOT_TOKEN_TYPE_NUMBER = 2,
    _GT_DOT_TOKEN_TYPE_OPERATOR = 3,
    _GT_DOT_TOKEN_TYPE_STRING = 4,
    _GT_DOT_TOKEN_TYPE_DELIMITER = 5
};

enum gt_attribute {
    _GT_ATTRIB_LABEL,
    _GT_ATTRIB_WEIGHT,
    _GT_ATTRIB_COLOR,
    _GT_ATTRIB_SHAPE,
    _GT_ATTRIB_STYLE,
    _GT_ATTRIB_DIRECTED,
    _GT_ATTRIB_WEIGHTED,
    _GT_ATTRIB_POSITION,
    _GT_ATTRIB_NAME,
    //add more here
    _GT_ATTRIB_USER  // this one must be the last
};

enum gt_layout_style {
    _GT_STYLE_DEFAULT,
    _GT_STYLE_SPRING,
    _GT_STYLE_PLANAR,
    _GT_STYLE_3D,
    _GT_STYLE_CIRCLE,
    _GT_STYLE_TREE,
    _GT_STYLE_BIPARTITE
};

class graphe {
public:
    typedef std::vector<int> ivector;
    typedef std::vector<ivector> ivectors;
    typedef std::map<int,gen> attrib;
    typedef std::pair<int,int> ipair;
    typedef std::vector<ipair> ipairs;
    typedef std::pair<double,double> dpair;
    typedef std::vector<dpair> dpairs;
    typedef std::vector<double> point;
    typedef std::vector<point> layout;
    typedef std::map<int,std::map<int,double> > sparsemat;
    typedef std::set<ipair> edgeset;
    typedef std::vector<std::map<int,int> > edgemap;

    class vertex {
        int m_subgraph;
        bool m_sorted;
        // used for traversing
        bool m_visited;
        bool m_on_stack;
        int m_low;
        int m_disc;
        int m_ancestor;
        int m_color;
        // used for drawing trees
        int m_position;
        int m_gaps;
        double m_prelim;
        double m_modifier;
        int m_children;
        // used for planar embedding
        bool m_embedded;
        int m_number;
        std::map<int,int> m_edge_faces;
        // used for grid drawing
        int m_left;
        int m_right;
        int m_x_offset;
        int m_y;
        // *
        attrib m_attributes;
        std::map<int,attrib> m_neighbor_attributes;
        ivector m_neighbors;
    public:
        vertex();
        vertex(const vertex &other);
        vertex& operator =(const vertex &other);
        void assign(const vertex &other);
        inline bool sorted() const { return m_sorted; }
        gen label() const;
        inline void set_label(const gen &s) { m_attributes[_GT_ATTRIB_LABEL]=s; }
        inline void set_subgraph(int s) { m_subgraph=s; }
        inline int subgraph() const { return m_subgraph; }
        inline void set_embedded(bool yes) { m_embedded=yes; }
        inline bool is_embedded() const { return m_embedded; }
        inline void set_number(int n) { m_number=n; }
        inline int number() const { return m_number; }
        inline void set_visited(bool yes) { m_visited=yes; }
        inline bool is_visited() const { return m_visited; }
        inline bool is_on_stack() const { return m_on_stack; }
        inline void set_on_stack(bool yes) { m_on_stack=yes; }
        inline void set_low(int l) { m_low=l; }
        inline int low() const { return m_low; }
        inline void set_disc(int t) { m_disc=t; }
        inline int disc() const { return m_disc; }
        inline void set_ancestor(int i) { m_ancestor=i; }
        inline void unset_ancestor() { m_ancestor=-1; }
        inline int ancestor() const { return m_ancestor; }
        inline void set_color(int c) { m_color=c; }
        inline int color() const { return m_color; }
        inline void inc_children() { ++m_children; }
        inline void clear_children() { m_children=0; }
        inline int children() const { return m_children; }
        inline bool is_leaf() const { return m_children==0; }
        inline void set_position(int p) { m_position=p; }
        inline int position() const { return m_position; }
        inline void set_gaps(int n) { m_gaps=n; }
        inline int gaps() const { return m_gaps; }
        inline void set_prelim(double val) { m_prelim=val; }
        inline double prelim() const { return m_prelim; }
        inline void set_modifier(double val) { m_modifier=val; }
        inline double modifier() const { return m_modifier; }
        inline void set_left(int l) { m_left=l; }
        inline int left() const { return m_left; }
        inline void set_right(int r) { m_right=r; }
        inline int right() const { return m_right; }
        inline void set_x_offset(int dx) { m_x_offset=dx; }
        inline int x_offset() const { return m_x_offset; }
        inline void set_y(int y) { m_y=y; }
        inline int y() const { return m_y; }
        inline const attrib &attributes() const { return m_attributes; }
        inline attrib &attributes() { return m_attributes; }
        inline void set_attribute(int key,const gen &val) { m_attributes[key]=val; }
        inline void set_attributes(const attrib &attr) { copy_attributes(attr,m_attributes); }
        inline const ivector &neighbors() const { return m_neighbors; }
        void add_neighbor(int i,const attrib &attr=attrib());
        attrib &neighbor_attributes(int i);
        const attrib &neighbor_attributes(int i) const;
        bool has_neighbor(int i,bool include_temp_edges=true) const;
        void move_neighbor(int i,int j,bool after=true);
        void remove_neighbor(int i);
        template<class Compare>
        inline void sort_neighbors(Compare comp) { sort(m_neighbors.begin(),m_neighbors.end(),comp); m_sorted=true; }
        inline void sort_neighbors() { sort(m_neighbors.begin(),m_neighbors.end()); m_sorted=true; }
        inline void clear_neighbors() { m_neighbors.clear(); m_neighbor_attributes.clear(); m_sorted=false; }
        void incident_faces(ivector &F) const;
        inline void add_edge_face(int nb,int f) { assert(m_edge_faces.find(nb)==m_edge_faces.end()); m_edge_faces[nb]=f+1; }
        inline void clear_edge_faces() { m_edge_faces.clear(); }
        inline int edge_face(int nb) { return m_edge_faces[nb]-1; }
    };

    class dotgraph {
        int m_index;
        attrib vertex_attr;
        attrib edge_attr;
        attrib chain_attr;
        ivector m_chain;
        int pos;
    public:
        dotgraph();
        dotgraph(const dotgraph &other);
        dotgraph(int i);
        dotgraph& operator =(const dotgraph &other);
        void assign(const dotgraph &other);
        inline int index() const { return m_index; }
        inline void set_index(int i) { m_chain[pos]=i; }
        inline const attrib &vertex_attributes() const { return vertex_attr; }
        inline const attrib &edge_attributes() const { return edge_attr; }
        inline const attrib &chain_attributes() const { return chain_attr; }
        inline attrib &vertex_attributes() { return vertex_attr; }
        inline attrib &edge_attributes() { return edge_attr; }
        inline attrib &chain_attributes() { return chain_attr; }
        inline const ivector &chain() const { return m_chain; }
        inline ivector &chain() { return m_chain; }
        inline int position() const { return pos; }
        inline void incr() { ++pos; if (int(m_chain.size())<=pos) m_chain.resize(pos+1,0); }
        inline void clear_chain() { pos=0; m_chain.resize(1); m_chain.front()=0; chain_attr.clear(); }
        inline bool chain_completed() { return m_chain.back()!=0; }
        inline bool chain_empty() { return pos==0 && m_chain.front()==0; }
    };

    class matching_maximizer {
        graphe *G;
        std::map<int,ivector> blossoms;
        std::map<int,int> forest;
        int mate(const ipairs &matching,int v);
        int find_root(int k);
        int root_distance(std::map<int,int>::const_iterator it);
        int root_distance(int v);
        int find_base(int v,int w);
        bool tree_path(int v,int w,ivector &path);
        std::map<int,ivector>::iterator in_blossom(int v);
        std::map<int,ivector>::iterator is_blossom_base(int v);
        void append_non_blossom_adjacents(int v,std::map<int,ivector>::const_iterator bit,ivector &lst);
        ivector adjacent(int v);
        inline ipair make_edge(int i,int j) { return std::make_pair(i<j?i:j,i<j?j:i); }
    public:
        matching_maximizer(graphe *gr);
        bool find_augmenting_path(const ipairs &matching,ivector &path);
        void find_maximum_matching(ipairs &matching);
    };

    class tree_node_positioner {
        graphe *G;
        layout *x;
        double hsep;
        double vsep;
        ivectors levels;
        ivector node_counters;
        ivector gap_counters;
        std::queue<int> placed;
        int depth;
        void walk(int i,int pass,int level=0,double modsum=0);
        void process_level(int i);
    public:
        tree_node_positioner(graphe *gr,layout *ly,double hs,double vs);
        void positioning(int apex);
    };

#ifdef HAVE_LIBGLPK
    class painter {
        graphe *G;
        ivectors values;
        ivector cover_number;
        ivector clique;
        ipairs col2ij;
        ivector branch_candidates;
        int lb;
        int ub;
        int nxcols;
        glp_prob *mip;
        double timeout;
        void compute_bounds(int max_colors);
        void make_values();
        void formulate_mip();
        static void mip_callback(glp_tree *tree,void *info);
    public:
        painter(graphe *gr,double tm=5.0) { G=gr; timeout=tm; }
        int color_vertices(ivector &colors,int max_colors=0);
        int select_branching_variable(glp_tree *tree);
    };
#endif

    class rectangle {
        double m_x;
        double m_y;
        double m_width;
        double m_height;
        bool m_locked_above;
        bool m_locked_right;
        layout *L;
    public:
        struct comparator {
            bool operator()(const rectangle &r1,const rectangle &r2) const {
                return r1.height()>r2.height();
            }
        };
        rectangle();
        rectangle(double X,double Y,double W,double H,layout *ly=NULL);
        rectangle(const rectangle &rect);
        rectangle& operator =(const rectangle &other);
        void assign(const rectangle &other);
        inline void set_anchor(double x,double y) { m_x=x; m_y=y; }
        inline double x() const { return m_x; }
        inline double y() const { return m_y; }
        inline double width() const { return m_width; }
        inline double height() const { return m_height; }
        inline void set_locked_above(bool yes) { m_locked_above=yes; }
        inline void set_locked_right(bool yes) { m_locked_right=yes; }
        inline bool is_locked_above() const { return m_locked_above; }
        inline bool is_locked_right() const { return m_locked_right; }
        bool intersects(const rectangle &other) const;
        bool intersects(const std::vector<rectangle> &rectangles) const;
        bool intersects(std::vector<rectangle>::const_iterator first,std::vector<rectangle>::const_iterator last) const;
        layout *get_layout() const { return L; }
    };

    class union_find {
        struct element {
            int id;
            int parent;
            int rank;
            element() { id=-1; }
        };
        std::vector<element> elements;
    public:
        union_find(int n) { elements.resize(n); }
        void make_set(int id);
        bool is_stored(int id);
        int find(int id);
        void unite(int id1,int id2);
    };

    class ostergard {
        graphe *G;
        int maxsize;
        bool found;
        double timeout; // in seconds
        ivector c;
        ivector incumbent;
        ivector clique_nodes;
        clock_t start;
        bool timed_out;
        void recurse(ivector &U,int size);
    public:
        ostergard(graphe *gr,double max_time=0) { G=gr; timeout=max_time; }
        int maxclique(ivector &clique);
    };

    struct ipairs_comparator {
        bool operator()(const ipair &a,const ipair &b) {
            return a.second<b.second;
        }
    };

    struct edges_comparator {
        graphe *G;
        bool operator()(const ipair &a,const ipair &b) {
            return G->weight(a)<G->weight(b);
        }
        edges_comparator(graphe *gr) { G=gr; }
    };

    struct degree_comparator {
        graphe *G;
        bool asc;
        bool operator()(int v,int w) {
            return (asc && G->degree(v)<G->degree(w)) ||
                    (!asc && G->degree(v)>G->degree(w));
        }
        degree_comparator(graphe *gr,bool ascending=true) { G=gr; asc=ascending; }
    };

    struct ivectors_degree_comparator {
        graphe *G;
        bool operator()(const ivector &a,const ivector &b) {
            int deg_a,deg_b;
            for (ivector_iter it=a.begin();it!=a.end();++it) {
                deg_a+=G->degree(*it);
            }
            for (ivector_iter it=b.begin();it!=b.end();++it) {
                deg_b+=G->degree(*it);
            }
            return deg_a*b.size()>=deg_b*a.size();
        }
        ivectors_degree_comparator(graphe *gr) { G=gr; }
    };

    typedef std::vector<vertex>::const_iterator node_iter;
    typedef std::map<int,attrib>::const_iterator neighbor_iter;
    typedef attrib::const_iterator attrib_iter;
    typedef std::vector<point>::const_iterator layout_iter;
    typedef ivector::const_iterator ivector_iter;
    typedef ivectors::const_iterator ivectors_iter;
    typedef ipairs::const_iterator ipairs_iter;
    typedef edgeset::const_iterator edgeset_iter;
    typedef edgemap::const_iterator edgemap_iter;
    typedef vertex* vptr;
    typedef std::vector<vptr> vpointers;
    static const gen FAUX;
    static const gen VRAI;
    static bool verbose;
    static int default_vertex_color;
    static int default_edge_color;
    static int default_vertex_label_color;
    static int default_highlighted_edge_color;
    static int default_highlighted_vertex_color;
    static int default_edge_width;
    static int bold_edge_width;
    // special graphs
    static const int clebsch_graph[];
    static const char* coxeter_graph[];
    static const int dodecahedron_graph[];
    static const int dyck_graph[];
    static const int grinberg_graph[];
    static const int grotzsch_graph[];
    static const int harries_graph_lcf[];
    static const int harries_wong_graph_lcf[];
    static const int heawood_graph[];
    static const int herschel_graph[];
    static const int mcgee_graph[];
    static const int pappus_graph[];
    static const int robertson_graph[];
    static const int soccer_ball_graph[];
    static const int tetrahedron_graph[];
    static const int octahedron_graph[];
    static const int icosahedron_graph[];
    static const int levi_graph[];
    static const int ljubljana_graph_lcf[];

private:
    const context *ctx;
    std::vector<vertex> nodes;
    attrib attributes;
    std::vector<std::string> user_tags;
    ivector marked_nodes;
    int disc_time;
    ivector discovered_nodes;
    std::stack<ipair> edge_stack;
    std::stack<int> node_stack;
    std::queue<int> node_queue;
    ivectors visited_edges;
    void clear_node_stack();
    void clear_node_queue();
    void message(const char *str) const;
    void message(const char *format,int a) const;
    void message(const char *format,int a,int b) const;
    void message(const char *format,int a,int b,int c) const;
    std::string giac_version() const;
    inline vertex &node(int i) { return nodes[i]; }
    bool dot_parse_attributes(std::ifstream &dotfile,attrib &attr);
    static bool insert_attribute(attrib &attr,int key,const gen &val,bool overwrite=true);
    static bool remove_attribute(attrib &attr,int key);
    static bool genmap2attrib(const gen_map &m,attrib &attr);
    static void attrib2genmap(const attrib &attr,gen_map &m);
    static void copy_attributes(const attrib &src,attrib &dest);
    void write_attrib(std::ofstream &dotfile,const attrib &attr) const;
    static int sets_union(const ivector &A,const ivector &B,ivector &U);
    static int sets_intersection(const ivector &A,const ivector &B,ivector &I);
    static int sets_difference(const ivector &A,const ivector &B,ivector &D);
    static point make_point(double x,double y) { point p(2,x); p.back()=y; return p; }
    static point make_point(double x,double y,double z) { point p(3,x); p[1]=y; p.back()=z; return p; }
    static void add_point(point &a,const point &b);
    static void subtract_point(point &a,const point &b);
    static void scale_point(point &p,double s);
    static double point_vecprod2d(const point &v,const point &w);
    static double point_dotprod(const point &p,const point &q);
    static void clear_point_coords(point &p);
    static double point_displacement(const point &p,bool sqroot=true);
    static double point_distance(const point &p,const point &q,point &pq);
    static void point_mirror(double a,double b,double c,const point &src,point &dest);
    static void point_lincomb(const point &p,const point &q,double d1,double d2,point &res);
    static void copy_point(const point &src,point &dest);
    void rand_point(point &p,double radius=1.0);
    static vecteur point2vecteur(const point &p);
    static bool points_coincide(const point &p,const point &q,double tol);
    static void copy_layout(const layout &src,layout &dest);
    static void rotate_layout(layout &x,double phi);
    static double layout_min(const layout &x,int d);
    static double layout_diameter(const layout &x);
    static void point2polar(point &p,double &r,double &phi);
    static bool sparse_matrix_element(const sparsemat &A,int i,int j,double &val);
    static void multiply_sparse_matrices(const sparsemat &A,const sparsemat &B,sparsemat &P,int ncols,bool symmetric=false);
    static void transpose_sparsemat(const sparsemat &A,sparsemat &T);
    void multilevel_recursion(layout &x,int d,double R,double K,double tol,int depth=0);
    int mdeg(const ivector &V,int i) const;
    void coarsening(graphe &G,const sparsemat &P,const ivector &V) const;
    void tomita(ivector &R,ivector &P,ivector &X,std::map<int,int> &m,bool store_matching) const;
    int cp_maxclique(ivector &clique);
    void cp_recurse(ivector &C,ivector &P,ivector &incumbent);
    int ost_maxclique(ivector &clique);
    void ost_recursive(ivector &U,int size,int &maxsize,ivector &incumbent,bool &found);
    void remove_isolated_node(int i);
    void find_cut_vertices_dfs(int i,std::set<int> &ap,int sg);
    void find_blocks_dfs(int i,std::vector<ipairs> &blocks,int sg);
    void find_bridges_dfs(int i,ipairs &B,int sg);
    int find_cycle_dfs(int i,int sg);
    bool find_path_dfs(int dest,int i,int sg,bool skip_embedded);
    static void sort_rectangles(std::vector<rectangle> &rectangles);
    static bool embed_rectangles(std::vector<rectangle> &rectangles,double maxheight);
    static bool segments_crossing(const point &p,const point &r,const point &q,const point &s,point &crossing);
    static bool point2segment_projection(const point &p,const point &q,const point &r,point &proj);
    void force_directed_placement(layout &x,double K,double R=DBL_MAX,double tol=0.01,bool ac=true);
    static bool get_node_position(const attrib &attr,point &p);
    void coarsening_mis(const ivector &V,graphe &G,sparsemat &P) const;
    void coarsening_ec(const ipairs &M,graphe &G,sparsemat &P) const;
    int best_quadrant(const point &p,const layout &x) const;
    void append_segment(vecteur &drawing, const point &p,const point &q,int color,int width,int style,bool arrow=false) const;
    void append_node(vecteur &drawing,const point &p,int color,int width,int shape) const;
    void append_label(vecteur &drawing,const point &p,const gen &label,int quadrant,int color=_BLACK) const;
    static int face_has_edge(const ivector &face,int i,int j);
    int first_neighbor_from_subgraph(const vertex &v,int sg) const;
    void set_nodes_embedded(const ivector &v,bool yes=true);
    void unembed_all_nodes();
    int planar_embedding(ivectors &faces);
    void set_embedding(const ivectors &faces);
    void clear_embedding();
    int choose_embedding_face(const ivectors &faces,int v);
    int choose_outer_face(const ivectors &faces);
    static void make_regular_polygon_layout(layout &x,const ivector &v,double R=1.0);
    bool edges_crossing(const ipair &e1,const ipair &e2,const layout &x,point &crossing) const;
    static void build_block_tree(int i,ivectors &blocks);
    static int common_element(const ivector &v1,const ivector &v2,int offset=0);
    void embed_children_blocks(int i,ivectors &block_tree,std::vector<ivectors> &blocks_faces);
    void periphericity(const ivector &outer_face,ivector &p);
    static void tree_layout2polar(layout &ly);
    void tree_height_dfs(int i,int level,int &depth);
    void make_product_nodes(const graphe &G,graphe &P) const;
    static void extract_path_from_cycle(const ivector &cycle,int i,int j,ivector &path);
    static void generate_nk_sets(int n,int k,std::vector<ulong> &v);
    void strongconnect_dfs(ivectors &components,int i,int sg);
    bool degrees_equal(const ivector &v,int deg=0) const;
    void lca_recursion(int u,const ipairs &p,ivector &lca_recursion,union_find &ds);
    void pathfinder(int i,ivector &path);
    void rdfs(int i,ivector &d,bool rec,int sg,bool skip_embedded);
    bool is_descendant(int v,int anc) const;
    static int pred(int i,int n);
    static int succ(int i,int n);
    static void arc_path(int i,int j,const ivector &cycle,ivector &path);
    void fold_face(const ivector &face, bool subdivide,int &label);
    void find_chords(const ivector &face,ipairs &chords);
    void augment(const ivectors &faces,int outer_face,bool subdivide=false);
    int saturation_degree(const vertex &v,std::set<int> &colors) const;
    int uncolored_degree(const vertex &v) const;
    bool is_partially_colored() const;
    void remove_maximal_clique(ivector &V) const;
    bool bipartite_matching_bfs(ivector &dist);
    bool bipartite_matching_dfs(int u,ivector &dist);

public:
    graphe(const context *contextptr=context0);
    graphe(const graphe &G);
    graphe(const std::string &name,const context *contextptr=context0);
    inline int rand_integer(int n) const { return giac::giac_rand(ctx)%n; }
    inline double rand_uniform() const { return giac::giac_rand(ctx)/(RAND_MAX+1.0); }
    inline double rand_normal() const { return giac::randNorm(ctx); }
    static bool is_real_number(const gen &g);
    static gen to_binary(int number,int chars);
    inline const context *giac_context() const { return ctx; }
    static gen make_idnt(const char* name,int index=-1,bool intern=true);
    void make_default_labels(vecteur &labels,int n,int n0=0,int offset=-1) const;
    static gen boole(bool b) { return b?VRAI:FAUX; }
    static gen word2gen(const std::string &word);
    static gen str2gen(const std::string &str,bool isstring=false);
    static std::string genstring2str(const gen &g);
    static std::string gen2str(const gen &g);
    static gen plusinf();
    void ivectors2vecteur(const ivectors &v,vecteur &res,bool sort_all=false) const;
    bool read_gen(const gen &g);
    void read_special(const int *special_graph);
    void read_special(const char **special_graph);
    void copy(graphe &G) const;
    inline void copy_nodes(const std::vector<vertex> &V) { nodes=std::vector<vertex>(V.begin(),V.end()); }
    void join_edges(const graphe &G);
    void clear();
    int tag2index(const std::string &tag);
    std::string index2tag(int index) const;
    int register_user_tag(const std::string &tag);
    void register_user_tags(const std::vector<std::string> &tags);
    inline const ivector &get_marked_nodes() const { return marked_nodes; }
    void get_marked_nodes(vecteur &V) const;
    void get_marked_nodes_in_subgraph(int s,ivector &m) const;
    inline void copy_marked_nodes(const ivector &mv) { marked_nodes=ivector(mv.begin(),mv.end()); }
    void mark_node(int v);
    void mark_node(const gen &v) { mark_node(node_index(v)); }
    bool unmark_node(int v);
    inline bool unmark_node(const gen &v) { return unmark_node(node_index(v)); }
    inline void unmark_all_nodes() { marked_nodes.clear(); }
    inline void sort_marked_nodes() { std::sort(marked_nodes.begin(),marked_nodes.end()); }
    void set_edge_visited(int i,int j);
    inline void set_edge_visited(const ipair &e) { set_edge_visited(e.first,e.second); }
    bool is_edge_visited(int i,int j) const;
    inline bool is_edge_visited(const ipair &e) const { return is_edge_visited(e.first,e.second); }
    inline void unvisit_all_edges() { visited_edges.clear(); }
    void move_neighbor(int v,int w,int ref=-1,bool after=true);
    template<class Compare>
    inline void sort_neighbors(int v,Compare comp) { node(v).sort_neighbors(comp); }
    gen to_gen() const;
    bool write_latex(const std::string &filename,const gen &drawing) const;
    bool write_dot(const std::string &filename) const;
    bool read_dot(const std::string &filename);
    inline bool is_empty() const { return nodes.empty(); }
    matrice weight_matrix() const;
    gen weight(int i,int j) const;
    inline gen weight(const ipair &edge) const { return weight(edge.first,edge.second); }
    int edge_count() const;
    inline int node_count() const { return nodes.size(); }
    vecteur vertices(int sg=-1) const;
    void unvisit_all_nodes(int sg=-1);
    void unset_all_ancestors(int sg=-1);
    void uncolor_all_nodes(int base_color=0,int sg=-1);
    void set_node_color(int i,int c) { node(i).set_color(c); }
    void dfs(int root,bool rec=true,bool clr=true,ivector *D=NULL,int sg=-1,bool skip_embedded=false);
    void bfs(int root,bool rec=true,bool clr=true,ivector *D=NULL,int sg=-1,bool skip_embedded=false);
    inline const ivector &get_discovered_nodes() const { return discovered_nodes; }
    bool is_connected(int sg=-1);
    bool is_biconnected(int sg=-1);
    bool is_triconnected(int sg=-1);
    void adjacent_nodes(int i,ivector &adj,bool include_temp_edges=true) const;
    void translate_indices_to(const graphe &G,const ivector &indices,ivector &dest) const;
    void translate_indices_from(const graphe &G,const ivector &indices,ivector &dest) const;
    void get_edges_as_pairs(ipairs &E, bool include_temp_edges=true,int sg=-1) const;
    vecteur edges(bool include_weights,int sg=-1) const;
    int add_node(const gen &v);
    inline int add_node(const gen &v,const attrib &attr) { int i=add_node(v); nodes[i].set_attributes(attr); return i; }
    void add_nodes(const vecteur &v);
    bool remove_node(int i);
    bool remove_node(const gen &v);
    void remove_nodes(const ivector &V);
    void remove_nodes(const vecteur &V);
    inline const vertex &node(int i) const { return nodes[i]; }
    inline const gen node_label(int i) const { assert(i>=0 && i<node_count()); return nodes[i].label(); }
    vecteur get_node_labels(const ivector &v) const;
    int node_index(const gen &v) const;
    int largest_integer_label() const;
    void set_subgraph(const ivector &v,int s);
    void merge_subgraphs(int s,int t);
    int max_subgraph_index() const;
    inline const attrib &graph_attributes() const { return attributes; }
    inline const attrib &node_attributes(int i) const { assert(i>=0 && i<node_count()); return node(i).attributes(); }
    const attrib &edge_attributes(int i,int j) const;
    attrib &edge_attributes(int i,int j);
    void attrib2vecteurs(const attrib &attr,vecteur &tags,vecteur &values) const;
    void add_edge(int i,int j,const gen &w=gen(1));
    void add_edge(int i,int j,const attrib &attr);
    inline void add_edge(const ipair &edge) { add_edge(edge.first,edge.second); }
    inline void add_edge(const ipair &edge,const attrib &attr) { add_edge(edge.first,edge.second,attr); }
    ipair add_edge(const gen &v,const gen &w,const gen &weight=gen(1));
    ipair add_edge(const gen &v,const gen &w,const attrib &attr);
    void add_temporary_edge(int i,int j);
    void remove_temporary_edges();
    bool remove_edge(int i,int j);
    inline bool remove_edge(const ipair &p) { return remove_edge(p.first,p.second); }
    bool has_edge(int i,int j) const;
    inline bool has_edge(ipair p) const { return has_edge(p.first,p.second); }
    ipair make_edge(const vecteur &v) const;
    bool edges2ipairs(const vecteur &E,ipairs &ev,bool &notfound) const;
    static void ipairs2edgeset(const ipairs &E,edgeset &Eset);
    bool nodes_are_adjacent(int i,int j) const;
    int in_degree(int index,bool count_temp_edges=true) const;
    int out_degree(int index,bool count_temp_edges=true) const;
    int degree(int index,bool count_temp_edges=true) const;
    vecteur degree_sequence(int sg=-1) const;
    void adjacency_matrix(matrice &m) const;
    void adjacency_sparse_matrix(sparsemat &sm) const;
    matrice incidence_matrix() const;
    inline void set_graph_attribute(int key,const gen &val) { attributes[key]=val; }
    inline void set_graph_attributes(const attrib &attr) { copy_attributes(attr,attributes); }
    void set_node_attribute(int index,int key,const gen &val);
    void set_edge_attribute(int i,int j,int key,const gen &val);
    bool get_graph_attribute(int key,gen &val) const;
    bool get_node_attribute(int i,int key,gen &val) const;
    bool get_edge_attribute(int i,int j,int key,gen &val) const;
    void discard_graph_attribute(int key);
    void discard_node_attribute(int i,int key);
    void discard_edge_attribute(int i,int j,int key);
    inline void set_name(const std::string &str) { set_graph_attribute(_GT_ATTRIB_NAME,str2gen(str,true)); }
    inline std::string name() const { gen s; if (get_graph_attribute(_GT_ATTRIB_NAME,s)) return genstring2str(s); else return ""; }
    bool is_directed() const;
    bool is_weighted() const;
    inline void set_directed(bool yes) { set_graph_attribute(_GT_ATTRIB_DIRECTED,boole(yes)); }
    inline void set_weighted(bool yes) { set_graph_attribute(_GT_ATTRIB_WEIGHTED,boole(yes)); }
    void make_weighted(const matrice &m);
    void make_directed() { set_directed(true); }
    void make_unweighted();
    void randomize_edge_weights(double a,double b,bool integral_weights=false);
    bool is_regular(int d) const;
    bool is_equal(const graphe &G) const;
    void underlying(graphe &G) const;
    void complement(graphe &G) const;
    bool isomorphic_copy(graphe &G,const ivector &sigma);
    bool relabel_nodes(const vecteur &labels);
    void induce_subgraph(const ivector &vi,graphe &G,bool copy_attrib=true) const;
    void subgraph(const ipairs &E,graphe &G,bool copy_attrib=true) const;
    bool is_subgraph(const graphe &G) const;
    void maximal_independent_set(ivector &ind) const;
    void maximize_matching(ipairs &matching);
    void find_maximal_matching(ipairs &matching) const;
    bool trail(const vecteur &v);
    bool demoucron(ivectors &faces);
    void create_random_layout(layout &x,int dim);
    void make_spring_layout(layout &x,int d,double tol=0.001);
    void make_circular_layout(layout &x,const ivector &hull,double A=0,double tol=0.005);
    void make_tutte_layout(layout &x,const ivector &outer_face);
    bool make_planar_layout(layout &x);
    void make_tree_layout(layout &x,double sep,int apex=0);
    void make_bipartite_layout(layout &x,const ivector &p1,const ivector &p2);
    void layout_best_rotation(layout &x);
    static rectangle layout_bounding_rect(layout &ly,double padding=0);
    static void pack_rectangles(std::vector<rectangle> &rectangles);
    static gen point2gen(const point &p,bool vect=false);
    static bool gen2point(const gen &g,point &p);
    static point layout_center(const layout &x);
    static void scale_layout(layout &x,double diam);
    inline bool is_tree() { return !is_directed() && edge_count()+1==node_count() && is_connected(); }
    bool is_forest();
    bool is_tournament() const;
    bool is_planar();
    bool is_clique() const;
    bool is_triangle_free() const;
    int tree_height(int root);
    void clique_stats(std::map<int,int> &m,bool store_matching=false);
    int maximum_clique(ivector &clique);
    void greedy_neighborhood_clique_cover_numbers(ivector &cover_numbers);
    bool clique_cover(ivectors &cover,int k=0);
    int maximum_independent_set(ivector &v) const;
    int girth(bool odd=false,int sg=-1);
    bool hakimi(const ivector &L);
    void make_plane_dual(const ivectors &faces);
    void make_lcf_graph(const ivector &jumps,int e);
    void make_lcf_graph(const int *j,int e);
    void make_sierpinski_graph(int n,int k,bool triangle);
    void make_shrikhande_graph();
    void make_complete_graph();
    void make_complete_multipartite_graph(const ivector &partition_sizes,layout *x=NULL);
    void make_petersen_graph(int n,int k,layout *x=NULL);
    bool make_kneser_graph(int n,int k);
    void make_path_graph();
    void make_cycle_graph();
    void make_grid_graph(int m,int n,bool torus=false);
    void make_web_graph(int n,int m,layout *x=NULL);
    void make_wheel_graph(int n,layout *x=NULL);
    void make_antiprism_graph(int n,layout *x=NULL);
    void make_complete_kary_tree(int k,int d);
    void make_random_tree(const vecteur &V,int maxd,bool addnodes=true);
    void make_random_planar(double p,int connectivity);
    void make_random(bool dir,const vecteur &V,double p);
    void make_random_bipartite(const vecteur &V,const vecteur &W,double p);
    void make_random_regular(const vecteur &V,int d,bool connected);
    void make_random_flow_network(const vecteur &V,int capacity);
    static void translate_layout(layout &x,const point &dx);
    void cartesian_product(const graphe &G,graphe &P) const;
    void tensor_product(const graphe &G,graphe &P) const;
    void connected_components(ivectors &components,int sg=-1,bool skip_embedded=false,int *count=NULL);
    int connected_components_count(int sg=-1);
    void strongly_connected_components(ivectors &components,int sg=-1);
    bool has_cut_vertex(int sg=-1,int i=0);
    void find_cut_vertices(ivector &articulation_points,int sg=-1);
    void find_blocks(std::vector<ipairs> &blocks,int sg=-1);
    void find_bridges(ipairs &B,int sg=-1);
    bool find_cycle(ivector &cycle,int sg=-1);
    bool find_path(int i,int j,ivector &path,int sg=-1,bool skip_embedded=false);
    bool find_eulerian_path(ivector &path);
    int eulerian_path_start(bool &iscycle) const;
    bool fleury(int start,ivector &path);
    void hierholzer(ivector &path);
    void contract_edge(int i,int j);
    void incident_edges(const ivector &V,edgeset &E);
    static bool edges_incident(const ipair &e1,const ipair &e2);
    void convex_hull(const layout &x,layout &hull);
    void edge_labels_placement(const layout &x);
    void draw_edges(vecteur &drawing,const layout &x);
    void draw_nodes(vecteur &drawing,const layout &x) const;
    void draw_labels(vecteur &drawing,const layout &x) const;
    void highlight_edges(const ipairs &E,int color);
    void highlight_nodes(const ivector &V,int color);
    void distance(int i,const ivector &J,ivector &dist,ivectors *shortest_paths=NULL);
    void allpairs_distance(matrice &m) const;
    void dijkstra(int src,const ivector &dest,vecteur &path_weights,ivectors *cheapest_paths=NULL);
    bool topologic_sort(ivector &ordering);
    bool is_arborescence() const;
    void reverse(graphe &G) const;
    void spanning_tree(int i,graphe &T,int sg=-1);
    void minimal_spanning_tree(graphe &T,int sg=-1);
    void lowest_common_ancestors(int root,const ipairs &p,ivector &lca_recursion);
    int lowest_common_ancestor(int i,int j,int root);
    void compute_st_numbering(int s,int t);
    vecteur get_st_numbering() const;
    void greedy_vertex_coloring_biggs(ivector &ordering);
    int greedy_vertex_coloring(const ivector &p);
    int exact_vertex_coloring(int max_colors=0);
    void get_node_colors(ivector &colors);
    bool is_bipartite(ivector &V1,ivector &V2,int sg=-1);
    bool is_vertex_colorable(int k);
    void dsatur();
    int color_count() const;
    ipair adjacent_color_count(int i) const;
    ipair chromatic_number_bounds();
    void store_layout(const layout &x);
    bool has_stored_layout(layout &x) const;
    int bipartite_matching(const ivector &p1,const ivector &p2,ipairs &matching);
    void line_graph(graphe &G) const;
    void transitive_closure(graphe &G,bool weighted=false);
    graphe &operator =(const graphe &other) { nodes.clear(); other.copy(*this); return *this; }
};

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC
#endif // __GRAPHE_H