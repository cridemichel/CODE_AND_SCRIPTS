#define MAXN 25000 /* maximum number of molecules to plot */
#define NUMBW 256
#define MAXAT 16 /* maximum number of aroms per molecule allowed
		       (equal to the number of greys) */

#define PI 2.0*acos(0.0)
#define STACKS 10
#define SLIDES 10
#define MGL_NO_FADE  1
#define MGL_FADE_LIN 2
#define MGL_FADE_QUAD 3
#define MGL_MAX_FRAMES 4

#define Sqr(x) (x)*(x)
struct colStruct 
{
  float rgba[4];
  char name[32];
} *mgl_col;
enum atom_types {MGL_ATOM_SPHERE, MGL_ATOM_DISK, MGL_ATOM_CYLINDER};
typedef enum atom_types atom_types_e;

struct atom_common {
  atom_types_e type;
  double rx;
  double ry;
  double rz;
  int atcol;
  int greyLvl;
};
struct atom_sphere
{
  struct atom_common common;
  double radius;
};
struct atom_disk
{
  struct atom_common common;
  double nx; 
  double ny; 
  double nz;
  double toprad;
  double radius;
  double height;
};
struct atom_cylinder
{
  struct atom_common common;
  double nx; 
  double ny; 
  double nz;
  double toprad;
  double botrad;
  double height;
};
union atom 
{
  struct atom_common common; 
  struct atom_disk disk;
  struct atom_sphere sphere;
  struct atom_cylinder cylinder;
};
enum bond_types {NONE, MGL_BOND_WIRE, MGL_BOND_CYLINDER};
typedef enum bond_types bond_types_e;
typedef union atom atom_u;
struct bond 
{
  /* indici relativi ai due atomi tra cui ci deve essere un bond */
  int from;
  int to;
  /* tipo del bond: ossia come deve essere rappresentato graficamente */
  bond_types_e type;
  double thickness; /* spessore del bond */
  int color; /* colore del bond */ 
};
typedef union atom atom_s; 
typedef struct bond bond_s;
struct molecule 
{
  atom_s *atom;
  bond_s *bond; 
  int numbonds;
  int numat;
};
struct global_settings
{
  int saved_counter;
  /*int saveimage=0, savedimg=0;*/
  int saveandquit;
  char *savefile;
  int drawcube;
  double *sig;
  /*double *height;*/
  int numAt;
  int setdiameter;
  int setheight;
  int *NumMols;
  int Width;
  int Height;
  int infos;
  int frameNo; 
  int fadeMode; /* default = linear fading */
  int NA; 
  int axon;/* perspective by default */
  int bw;
  double L; /* lenght of an edge of the cubic box */
  float degx;
  float degy;
  float degz;
  float deginc;
  int *colIdxCol;
  int *colIdxBw;
  double viewangle;
  double near;
  double far;
  int setvp;
  double diameter;
  double height;
  double  dist;
  double ivpx;
  double ivpy;
  double ivpz;
  int default_bw;
  int default_col;
  double defbondthick;
  int defbondcol;
};

