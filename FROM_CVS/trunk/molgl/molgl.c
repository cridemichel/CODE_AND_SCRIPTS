#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "molgl.h"

const int NUMCOLS = 746;
float mgl_bw[NUMBW][4];
char *mglrgb[]={
#include "mglrgb.h"
};
struct molecule *mols = NULL;
struct global_settings globset;

/* array con i valori di default */
char inputFile[512];
GLuint *atomsList = NULL;

/* mgl_*[colIdx*[j]] is the color of the j-th atom */ 
void setLight(void)
{
  GLfloat light_ambient0[] = { 1.0, 1.0, 1.0, 0.15 };
  GLfloat light_diffuse0[] = { 1.0, 1.0, 1.0, 0.15 };
  GLfloat light_specular0[] = { 1.0, 1.0, 1.0, 0.15 };
  GLfloat light_ambient1[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat light_diffuse1[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat light_specular1[] = { 1.0, 1.0, 1.0, 1.0 }; 
 
  /*	light_position is NOT default value	*/
  GLfloat light_position0[] = { 10.0, 10.0, 10.0, 0.0 };
  GLfloat light_position1[] = { -10.0, -10.0, -10.0, 0.0};
  GLfloat lmodel_ambient[] = { 0.2, 0.2, 0.2, 0.15 };
  GLfloat local_view[] = { 0.0 };
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
  glLightModelfv(GL_LIGHT_MODEL_LOCAL_VIEWER, local_view);
  
  //glLightModeli(GL_LIGHT_MODEL, GL_LIGHT_TWO_SIDE);

  glLightfv (GL_LIGHT0, GL_AMBIENT, light_ambient0);
  glLightfv (GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
  glLightfv (GL_LIGHT0, GL_SPECULAR, light_specular0);
  glLightfv (GL_LIGHT0, GL_POSITION, light_position0);
  glLightfv (GL_LIGHT1, GL_AMBIENT, light_ambient1);
  glLightfv (GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
  glLightfv (GL_LIGHT1, GL_SPECULAR, light_specular1);
  glLightfv (GL_LIGHT1, GL_POSITION, light_position1);

}

/* =============================== >>> setBW <<< ===========================*/
void setBW(void)
{
  int nc;
  for (nc = 0; nc < NUMBW; ++nc)
    {
      mgl_bw[nc][0] = ((float) nc) / 255.0;
      mgl_bw[nc][1] = ((float) nc) / 255.0;
      mgl_bw[nc][2] = ((float) nc) / 255.0;
      mgl_bw[nc][3] = 1.0; 
    }  
  /* 255 levels of gray */
}

/*  Initialize material property and light source. */
void myinit (void)
{
    setLight();
    //glFrontFace (GL_CW);
    glEnable (GL_LIGHTING);
    glEnable (GL_LIGHT0);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    
#if 0
    glEnable(GL_AUTO_NORMAL);
    glEnable(GL_NORMALIZE);
#endif
    glDepthFunc(GL_LESS);
    glEnable(GL_DEPTH_TEST);
    glClearColor(1.0, 1.0, 1.0, 0.0);
}

/* ========================== >>> setColor <<< =============================*/
void setColor(float col[4], double ff)
{
  /* col is the specular color, diffusie and the ambient are calculated
     scaling this by df and af respectively */
  float mat[4]; /* atoms color */
  float df = 0.95, af = 0.2;
#if 0
  int i;
#endif
#if 0
    mat[0] = ff*af*col[0]; mat[1] = ff*af*col[1]; 
    mat[2] = ff*af*col[2]; mat[3] = ff*col[3];
    glMaterialfv (GL_FRONT, GL_AMBIENT, mat);
    mat[0] = ff*df*col[0]; mat[1] = ff*df*col[1]; mat[2] = ff*df*col[2];	
    glMaterialfv (GL_FRONT, GL_DIFFUSE, mat);
    mat[0] = ff*col[0]; mat[1] = ff*col[1]; mat[2] = ff*col[2];
    glMaterialfv (GL_FRONT, GL_SPECULAR, mat);
    glMaterialf (GL_FRONT, GL_SHININESS, 0.9*128.0);
#endif
  mat[0] = af*col[0]; mat[1] = af*col[1]; 
  mat[2] = af*col[2]; mat[3] = ff*col[3];
  glMaterialfv (GL_FRONT, GL_AMBIENT, mat);
  mat[0] = df*col[0]; mat[1] = df*col[1]; mat[2] = df*col[2];	
  glMaterialfv (GL_FRONT, GL_DIFFUSE, mat);
  mat[0] = col[0]; mat[1] = col[1]; mat[2] = col[2];
  glMaterialfv (GL_FRONT, GL_SPECULAR, mat);
  glMaterialf (GL_FRONT, GL_SHININESS, 0.9*128.0);
}
/* ======================== >>> calcFadeFact <<< ===========================*/
double  calcFadeFact(int mode, int nf)
{
  double ff;
  if (mode == MGL_FADE_LIN)
    {
      ff =  1.0 - ((double)nf) / ((double)globset.frameNo) ;
    }
  else if(mode == MGL_FADE_QUAD)
    {
      ff = 1.0 / Sqr(((double)nf) / ((double)globset.frameNo));
    }
  else ff = 1;

  return ff;
}
/* ========================== >>> displayMol <<< ===========================*/
void displayAtom(int nf, int nm, int na)
{
  float fadeFact;
  GLUquadricObj* ss;
  atom_s *atom;
  glPushMatrix();
  atom = &mols[nm].atom[na];
  glTranslatef(atom->common.rx,atom->common.ry,atom->common.rz);/* 1st atom */ 
  
  fadeFact = calcFadeFact(globset.fadeMode, nf);
  
  if (atom->common.greyLvl)
    {
      setColor(mgl_bw[atom->common.greyLvl], fadeFact);
    }
  else 
    {
      if (globset.bw)
	setColor(mgl_bw[globset.colIdxBw[na]], fadeFact);
      else
	{
	  if (atom->common.atcol>=0 && atom->common.atcol<NUMCOLS)
	    {
	      setColor(mgl_col[atom->common.atcol].rgba, fadeFact);
	    }
	  else
	    setColor(mgl_col[globset.colIdxCol[na]].rgba, fadeFact);
	}
    }
  if (atom->common.type==MGL_ATOM_SPHERE)
    {
      glutSolidSphere (atom->sphere.radius, STACKS, SLIDES);
    }
  else if (atom->common.type==MGL_ATOM_DISK)
    {
      ss = gluNewQuadric();
      gluCylinder(ss, atom->disk.radius, 
		      atom->disk.radius, 
		      atom->disk.height, 
		      STACKS, SLIDES);
    }
  else if (atom->common.type==MGL_ATOM_CYLINDER)
    {
      ss = gluNewQuadric();
      gluCylinder(ss, atom->cylinder.toprad, 
		      atom->cylinder.botrad, 
		      atom->cylinder.height, 
		      STACKS, SLIDES);
    }
 /*ss = gluNewQuadric();
    gluSphere(ss, sig[na], 12, 12);
  */
  glPopMatrix();

}
/* ========================= >>> buildAtomsList <<< ======================= */
void buildAtomsList()
{
  int i, j, nf;
  /*float col[4];*/

  /* NOTE: frameno is the number of frames read from postions file */
  for (nf = 0; nf < globset.frameNo; ++nf )
    {
      atomsList[nf] = glGenLists(1);
      glNewList(atomsList[nf], GL_COMPILE);
      for(i = 0; i < globset.NumMols[nf]; ++i)
	{
	  displayAtom(nf, i, j);
	}
      glEndList();
    }
}

void myReshape(int w, int h);

/* ======================== >>> printStr <<< =============================== */
void printStr(int x, int y, const char* text)
{
  /* print a string at the position (x,y) */
  int i;
  
  glRasterPos2i(x, y);
  for (i = 0; i < strlen(text); ++i)
    {
      glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, text[i]);
    }
}
/* ======================== >>> onScreenInfo <<< ========================== */
void onScreenInfo()
{
  char text[255];
  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0, globset.Width, 0, globset.Height);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glColor3f(0.8, 0.0, 0.0);
  sprintf(text, "degx: %.1f degy: %.1f degz: %.1f", globset.degx, globset.degy, globset.degz);
  printStr(10, globset.Height - 20, text);
  sprintf(text, "L: %.1f deginc: %.1f", globset.L, globset.deginc);
  printStr(10, globset.Height - 45, text);
  myReshape(globset.Width, globset.Height);
  glPopAttrib();
  glPopMatrix();
}

#include<png.h>
/*extern PNG_EXPORT(void,png_set_zbuf_size);*/
void save_image(void)
{
   FILE *fp;
   char *fn;
   char numstr[300];
   png_structp png_ptr;
   png_infop info_ptr;
   png_bytep *row_pointer;
   unsigned char *image;
   int width,height;
   unsigned int NbBytes;
   int i, j ;
   int pixel_size=3;
   /* Allocate our buffer for the image */
   width= glutGet(GLUT_WINDOW_WIDTH);
   height=glutGet(GLUT_WINDOW_HEIGHT) ;
   NbBytes= 3*width*height*sizeof(char);
   if ((image = malloc(NbBytes)) == NULL) {
      fprintf(stderr,"Failed to allocate memory for image\n");
      return;
   }

   glPixelStorei(GL_PACK_ALIGNMENT,1);

   /* Open the file */
   fn = malloc(sizeof(char)*(strlen(globset.savefile)+257));
   strcpy(fn,globset.savefile);
   if (globset.saved_counter > 0)
     {
       snprintf(numstr, 256, "%d",globset.saved_counter);
       strcat(fn, numstr);
     }
   if ((fp = fopen(fn,"wb")) == NULL) {
      fprintf(stderr,"Failed to open file for window dump\n");
      free(fn);
      free(image);
      return;
   }

   png_ptr = png_create_write_struct
       (PNG_LIBPNG_VER_STRING, NULL,
	NULL, NULL);
   if (!png_ptr)
     {
       free(image);
       free(fn);
       return;
     }
   info_ptr = png_create_info_struct(png_ptr);
   if (!info_ptr)
     {
       png_destroy_write_struct(&png_ptr,
			 	(png_infopp)NULL);
       free(image);
       free(fn);
       return;
     }
   if (setjmp(png_jmpbuf(png_ptr)))
    {
       png_destroy_write_struct(&png_ptr, &info_ptr);
       free(image);
       free(fn);
       fclose(fp);
       return;
    }
   png_init_io(png_ptr, fp);
    /* turn on or off filtering, and/or choose
       specific filters.  You can use either a single
       PNG_FILTER_VALUE_NAME or the logical OR of one
       or more PNG_FILTER_NAME masks. */
   /*
   png_set_filter(png_ptr, 0,
       PNG_FILTER_NONE  | PNG_FILTER_VALUE_NONE |
       PNG_FILTER_SUB   | PNG_FILTER_VALUE_SUB  |
       PNG_FILTER_UP    | PNG_FILTER_VALUE_UP   |
       PNG_FILTER_AVE   | PNG_FILTER_VALUE_AVE  |
       PNG_FILTER_PAETH | PNG_FILTER_VALUE_PAETH|
       PNG_ALL_FILTERS);*/
   /* set the zlib compression level */
   png_set_compression_level(png_ptr,
		      	     Z_BEST_COMPRESSION);
   
   /* set other zlib parameters */
   png_set_compression_mem_level(png_ptr, 8);
   png_set_compression_strategy(png_ptr,
				Z_DEFAULT_STRATEGY);
   png_set_compression_window_bits(png_ptr, 15);
   png_set_compression_method(png_ptr, 8);
   png_set_compression_buffer_size(png_ptr, 8192);
   png_set_IHDR(png_ptr, info_ptr, width, height,
       8,  PNG_COLOR_TYPE_RGB,  PNG_INTERLACE_NONE,
       PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT); 

   /* Copy the image into our buffer */
   /*glReadBuffer(GL_BACK_LEFT);*/
   glReadPixels(0,0,width,height,GL_RGB,GL_UNSIGNED_BYTE,image);
   row_pointer=(png_bytep*)malloc(sizeof(png_bytep)*height);
   row_pointer = png_malloc(png_ptr,
      height*sizeof(png_bytep));
   for (i=0; i<height; i++)
      row_pointer[i]=png_malloc(png_ptr,
         width*pixel_size);
   /* Write the raw file */
   /* fprintf(fptr,"P6\n%d %d\n255\n",width,height); for ppm */
   for (j=height-1;j>=0;j--) {
      for (i=0;i<width;i++) {
         row_pointer[height-j-1][3*i+0]=image[3*j*width+3*i+0];
         row_pointer[height-j-1][3*i+1]=image[3*j*width+3*i+1];
         row_pointer[height-j-1][3*i+2]=image[3*j*width+3*i+2];
      }
   }
   png_set_rows(png_ptr, info_ptr, row_pointer);
   png_write_png(png_ptr, info_ptr,  PNG_TRANSFORM_IDENTITY, NULL);
   
   fclose(fp);
   /* Clean up */
   free(fn);
   free(image);
}

/* ======================== >>> display <<< ===============================*/
void display (void)
{
  int j, nf;
  /*
   * double fadeFact;
   */
  if (globset.saveandquit)
    glutHideWindow();


  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  glLoadIdentity();
  
  if (!globset.axon) 
    glTranslatef(0.0, 0.0, -(globset.near+globset.L/2.0));
  else
    glTranslatef(0.0, 0.0,-globset.L/2.0);
  glTranslatef(0.0, 0.0, -globset.dist);
  glPushMatrix ();
  /* NOTE:
     First arg is the degrees of the rotation, the others are the component
     of the vector around which we perfomr the rotation */
  glRotatef(-globset.degz, 0.0, 0.0, 1.0);
  glRotatef(-globset.degy, 0.0, 1.0, 0.0);
  glRotatef(-globset.degx, 1.0, 0.0, 0.0);

  setColor(mgl_bw[0], 1.0);
  if (globset.drawcube)
    glutWireCube(globset.L);
 
  for (nf = 0; nf < globset.frameNo; ++nf)
    {
      if (nf == 0)
	{
	  glDepthMask(GL_TRUE);
	}
      else 
	{
	  glDepthMask(GL_FALSE);
	}
      
      glCallList(atomsList[nf]);
	
    }
  glPopMatrix ();
  
  if (globset.infos) onScreenInfo();
  /*glFlush ();*/
  glutSwapBuffers();
  if (globset.saveandquit==1)
    {
      save_image();
      exit(0);
    }
}


/* =========================== >>> setproj <<< ============================*/
void setproj(void)
{
  double L;
  L = globset.L;
  /*int w = Width, h = Height;*/
  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  if (globset.axon)
    {
      /*printf("orto\n");*/
      glOrtho (-L, L, -L, 
      	       L, -2.0*L, 2.0*L);
      /*if (w <= h) 
	glOrtho (-L/2.0, L/2.0, -L*(GLfloat)h/(GLfloat)w/2.0, 
		 L*(GLfloat)h/(GLfloat)w/2.0, -L/2.0, L/2.0);
      else 
	glOrtho (-L*(GLfloat)w/(GLfloat)h/2.0, 
	L*(GLfloat)w/(GLfloat)h/2.0, -L/2.0, L/2.0, -L/2.0, L/2.0);*/
    }
  else 
    {
      globset.near = L / tan(PI*globset.viewangle/360.0) / 2.0;
      globset.far = 100*globset.near;
      gluPerspective(globset.viewangle, (GLdouble) globset.Width / (GLdouble) globset.Height, 
		     globset.near, globset.far);
    }

  glMatrixMode (GL_MODELVIEW);
}

/* ========================== >>> myReshape <<< ========================*/
void myReshape(int w, int h)
{
  globset.Width = w;
  globset.Height = h;
  glViewport (0, 0, w, h);
  setproj();
  glMatrixMode (GL_MODELVIEW);
  glLoadIdentity();
}
void print_usage(void)
{
  printf ("USAGE:\n");
  printf("molgl [-h/--help | --saveandquit/-sq | --pngfile/-f <filename> \n");
  printf("| --viewpoint/-vp (x,y,z) | --diameter/-d <atoms_diameter> | --noinfos/-ni\n");
  printf("| --nobox|-nb ] <input_file> \n");
}
/* ============================= >>> args <<< ============================= */
void args(int argc, char* argv[])
{
  int i=1;
  int ii;

  if (argc == 1)
    {
      print_usage();
      exit(-1);
    }
  else if (argc > 2)
    {
      while (i < argc-1)
	{
	  if (!strcmp(argv[i], "-h") || !strcmp(argv[i],"--help"))
	    {
	      print_usage();
	      exit(-1);
	    }
	  else if (!strcmp(argv[i],"--saveandquit") || !strcmp(argv[i],"-sq"))
	    {
	      globset.saveandquit = 1;
	    }
	  else if (!strcmp(argv[i],"--pngfile") || !strcmp(argv[i],"-f"))
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: You must supply a file name!\n");
		  exit(-1);
		}
	      globset.savefile = malloc(sizeof(char)*(strlen(argv[i])+1));
	      strcpy(globset.savefile,argv[i]);
	    }
	  else if (!strcmp(argv[i],"--viewpoint")||!strcmp(argv[i],"-vp"))
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: You must supply the viewpoint (x,y,z)!\n");
		  exit(-1);
		}
	      sscanf(argv[i], "(%lf,%lf,%lf)", &globset.ivpx, &globset.ivpy, &globset.ivpz);
	      globset.setvp = 1;
	    }
	  else if (!strcmp(argv[i],"--diameter")|| !strcmp(argv[i],"-d"))
	    {
	      i++;
	      globset.setdiameter = 1;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: You must supply the viewpoint (x,y,z)!\n");
		  exit(-1);
		}
	      globset.diameter = atof(argv[i]);
	    }
	  else if (!strcmp(argv[i],"--nobox")|| !strcmp(argv[i],"-nb"))
	    {
	      globset.drawcube = 0;
	    }
	  else if (!strcmp(argv[i],"--noinfos")|| !strcmp(argv[i],"-ni"))
	    {
	      globset.infos = 0;
	    }
	  else
	    {
	      fprintf(stderr, "ERROR: Invalid argumet!\n");
	      exit(-1);
	    }
      	  i++;
	}
    }

  if (i == argc)
    {
      fprintf(stderr, "ERROR: You must supply an input file!\n");
      exit(-1);
    }
  
      
  if (globset.savefile == NULL)
    globset.savefile = "molglimg.png";

  strcpy(inputFile, argv[i]);
}
void dropSpaces(char *S);
int getColByName(const char* name);

int parsecol(char *str)
{
  int colNum;
  char* eptr;
  
  colNum = (int) strtod(str, &eptr);
  if (eptr == str) /* not a number */
    {
      dropSpaces(str);
      return  getColByName(str);
      /* Find the number of the color named 's1'*/
      /*printf("col:%s:, %d\n", S, colIdxCol[j]);
      */
    }
  else
    {
      /*printf("Color Number: %d\n", colNum);*/
      return colNum;
    }
}
/* ========================== >>> assignAtom <<< ===========================*/
void assignAtom(int nf, int i, int j, const char* L)
{
  char s1[128], s2[128], s3[128], s4[128], s5[128], s6[128], s7[128], s8[128];

  if (sscanf(L,"%s %s %s %s %s %s", s1, s2, s3, s4, s5, s6) == 6)
    {
      /*printf("Uso il raggio specificato per l'atomo [%d][%d]\n", i, j);
      printf("Uso il livello di grigio: %d per l'atomo [%d][%d]",
	     atoi(s5),i, j);*/
      rx[j][i][nf] = atof(s1);
      ry[j][i][nf] = atof(s2);
      rz[j][i][nf] = atof(s3);
      /* nx, ny, nz sono le componenti del vettore normale al dischetto */
      nx[j][i][nf] = atof(s4);
      ny[j][i][nf] = atof(s5);
      nz[j][i][nf] = atof(s6);
      radius[j][i] = globset.sig[j];
      greyLvl[j][i] = 0; /*colIdxBW[j];// default value of grey level */
      atcol[j][i]  = -1;
    }
  else if (sscanf(L,"%s %s %s %s %s %s @ %s C[%[^]]", s1, s2, s3, s4, s5, s6, s7, s8) == 8)
    {
      /*printf("Uso il raggio specificato per l'atomo [%d][%d]\n", i, j);
      printf("Uso il livello di grigio: %d per l'atomo [%d][%d]",
	     atoi(s5),i, j);*/
      rx[j][i][nf] = atof(s1);
      ry[j][i][nf] = atof(s2);
      rz[j][i][nf] = atof(s3);
      /* nx, ny, nz sono le componenti del vettore normale al dischetto */
      nx[j][i][nf] = atof(s4);
      ny[j][i][nf] = atof(s5);
      nz[j][i][nf] = atof(s6);
      radius[j][i] = atof(s7);
      greyLvl[j][i] = 0; /*colIdxBW[j];// default value of grey level */
      atcol[j][i]  = parsecol(s8);
    }
  else if (sscanf(L,"%s %s %s @ %s $ %s ", s1, s2, s3, s4, s5) == 5 )
    {
      /*printf("Uso il raggio specificato per l'atomo [%d][%d]\n", i, j);
      printf("Uso il livello di grigio: %d per l'atomo [%d][%d]",
	     atoi(s5),i, j);*/
      rx[j][i][nf]= atof(s1);
      ry[j][i][nf]= atof(s2);
      rz[j][i][nf]= atof(s3);
      /*greylLvl[j][i] = colIdxBW[j];*/
      /* default value of grey level*/
      radius[j][i] = atof(s4);
      greyLvl[j][i] = atoi(s5);
      atcol[j][i]  = -1;
    }
  if (sscanf(L,"%s %s %s C[%[^]]", s1, s2, s3, s4) == 4 )
    {
      rx[j][i][nf] = atof(s1);
      ry[j][i][nf] = atof(s2);
      rz[j][i][nf] = atof(s3);
      radius[j][i] = sig[j];
      greyLvl[j][i] = 0;
      atcol[j][i] = parsecol(s4);
    }
  else if (sscanf(L,"%s %s %s @ %s C[%[^]]", s1, s2, s3, s4, s5) == 5)
    {
      /* printf("Uso il raggio specificato per l'atomo [%d][%d]\n", i, j);
      */
      rx[j][i][nf]= atof(s1);
      ry[j][i][nf]= atof(s2);
      rz[j][i][nf]= atof(s3);
      //greylLvl[j][i] = colIdxBW[j];// default value of grey level
      radius[j][i] = atof(s4);
      greyLvl[j][i] = 0;
      atcol[j][i] = parsecol(s5);
    }
  else if (sscanf(L,"%s %s %s @ %s ", s1, s2, s3, s4) == 4 )
    {
      /* printf("Uso il raggio specificato per l'atomo [%d][%d]\n", i, j);
      */
      rx[j][i][nf]= atof(s1);
      ry[j][i][nf]= atof(s2);
      rz[j][i][nf]= atof(s3);
      //greylLvl[j][i] = colIdxBW[j];// default value of grey level
      radius[j][i] = atof(s4);
      greyLvl[j][i] = 0;
      atcol[j][i]  = -1;
    }
  else if (sscanf(L,"%s %s %s ", s1, s2, s3) == 3 )
    {
      rx[j][i][nf]= atof(s1);
      ry[j][i][nf]= atof(s2);
      rz[j][i] [nf]= atof(s3);
      //greylLvl[j][i] = colIdxBW[j];// default value of grey level
      radius[j][i] = sig[j];
      greyLvl[j][i] = 0;
      atcol[j][i]  = -1;
    }
  else
    {
      printf("Line: %s\n", L);
      printf("ERROR: Bad format in input file! Exiting...\n");
      exit(-1);
    }
} 

/* ========================== >>> getColByName <<< ========================= */
int getColByName(const char* name)
{
  /* NOTE: Not optimized search !!!! */
  int nc;
  for (nc = 0; nc < NUMCOLS; ++nc)
    {
      if (!strcmp(mgl_col[nc].name, name))
	return nc;
    }
  printf("ERROR: Unrecognized color %s!\n", name);
  exit(-1);
}

void add_atom(int cf, int curmol, int curat)
{
  mols[cf][curmol].atom = realloc(mols[cf][curmol].atom, sizeof(atom_u)*curat);
  mols[cf][curmol].atom[curat-1]
}

void add_mol(int cf, int curmol)
{ 
  mols[cf] = realloc(mols[cf], sizeof(struct molecules)*curmol);
  mols[cf][curmol-1].atom = NULL;
  mols[cf][curmol-1].bond = NULL;
}

void add_frame(int curframes)
{
  int nf;
  nf = curframes + 1;
  globset.NumMols = realloc(globset.NumMols, sizeof(int)*nf);
  globset.NumMols[nf] = 0;
  atomsList = realloc(atomsList, sizeof(GLuint)*nf);
  mols = realloc(mols, sizeof(struct molecule*)*nf);
  mols[nf-1] = NULL;
}

/* =========================== >>> readLine <<< =============================*/
void readLine(FILE* stream, char* L)
{
  if (fscanf(stream, "%[^\n]\n", L) < 1)
    {
      fprintf(stderr,"ERROR: Void line...aborting\n");
      exit(-1);
    }
}

/* ========================== >>> dropSpaces <<< ==========================*/
void dropSpaces(char *S)
{
  char s1[512];
  int i;

  sscanf(S, " %[^\n]", s1);/* drop initial spaces */
//  printf("----->>>>>%s<<<\n", s1);
  /* and now final spaces */
  for (i = strlen(s1) - 1; (i >= 0) && (s1[i] == ' '); --i)
    {
      s1[i] = '\0';  
    }
  strcpy(S, s1);
}

/* ========================== >>> assignCol <<< ===========================*/
void assignCol(char* S, int j)
{
  int colNum;
  char* eptr;

  colNum = (int) strtod(S, &eptr);
  if (eptr == S) /* not a number */
    {
      dropSpaces(S);
      globset.colIdxCol[j] = getColByName(S);
      /* Find the number of the color named 's1'*/
      /*printf("col:%s:, %d\n", S, colIdxCol[j]);
      */
    }
  else
    {
      //printf("Color Number: %d\n", colNum);
      globset.colIdxCol[j] = colNum;
    }
}

/* ========================== >>> pareseLine <<< =========================== */
int parseLine(const char* Line, int* nf, int* i, int alloc)
{
  /*
    return true if this is a parameter, 0 otherwise
  */
  char parName[128], parVal[512], s1[512];
  int lett, j, jj, ii;
  /* Syntax:
     <parname> : <value>
     where <parname> is of the form ".<name>"
   */
  
  lett =  sscanf(Line, " %[^:#\n] : %[^#\n]", parName, parVal);
  //printf("lett: %d pn:%s|pv: %s\n", lett, parName, parVal);
  if (lett == 0) return 1; /* 2 = comment only */
  if (parName[0] != '.') return 0; // no a parameter line!
  
  /* new frames */
  if (!strcmp(parName, ".newmol"))
    {
      ++(*i);
      if (alloc)
	add_mol(*i);
    }
  else if (!strcmp(parName, ".newframe"))
    {
      if (!((*nf == 0) && (*i == 0)))
	{
	  NumMols[nf] = i;
	  ++(*nf);
	  *i = 0;
	  if (alloc)
      	    /* allocate a new frame here */
	    add_frame(*nf);
	}
      return 1;
    }
      
  /* Atoms radi */
  if (!strcmp(parName, ".atomRad"))
    {
      /* Build a string of this type: "%f , %f , ..." with NA '%f' */
      //strcpy(s1, strtok(parVal, ","));
      //printf("radius[0]: %s\n", s1);
      ns = strtok(parVal, ",");
      a = 0;
      if (globset.sig)
	free(globset.sig);
      while(ns)
	{
	  //strcpy(s1, strtok(NULL, ","));
	  globset.sig = realloc(globset.sig,a+1);
	  globset.sig[a] = atof(ns);
	  a++;
	  ns = strtok(NULL, ",");
	}
    
      //printf("---->%f %f\n", sig[0], sig[1]);
      return 1;
    }
  if (!strcmp(parName, ".Vol"))
    {
      globset.L = atof(parVal);
      globset.L = cbrt(globset.L);
      return 1;
    }
  
  if (!strcmp(parName, ".fadeMode"))
    {
      globset.fadeMode = atoi(parVal);
      return 1;
    }
  
  if (!strcmp(parName, ".atomCol"))
    {
      ns = strtok(parVal, ",");
      j = 0; 
      if (globset.colIdxCol)
	free(globset.colIdxCol);
      while(ns)
	{
	  strcpy(s1, ns);
	  globset.colIdxCol = realloc(globset.colIdxCol, j+1);
	  assignCol(s1, j);
	  ns = strtok(NULL, ",");
	  j++;
	}
    
      return 1;
    }
  if (!strcmp(parName, ".atomBw"))
    {
      /* Build a string of this type: "%f , %f , ..." with NA '%f' */
      ns = strtok(parVal, ",");
      if (globset.colIdxBw)
	free(globset.colIdxBw);
      a = 0;
      while (ns)
	{
	  globset.colIdxBw = realloc(globset.colIdxBw, a+1);
	  globset.colIdxBw[a] = atoi(ns);
	  ns = strtok(NULL, ",");
	  a++;
	}
      return 1;
    }

  if (!strcmp(parName, ".numAt"))
    {
      globset.NA = atoi(parVal);
      return 1;
    }
  printf("ERROR: Invalid parameter %s!\n", parName);
  exit(-1);
}
void setdefaults_after_fakeread(void)
{
  int a;
  if (globset.NA)
    {
      for (a = 0; a < globset.NA; ++a)
	{
	  /*  sig[i] = 0.5;*/
	  globset.colIdxCol[a] = a;
	  /* first (and only those ones) 25 grey are quite different in this way */
	  if (i < 25)
	    globset.colIdxBw[a] = a*10;
	  else 
	    globset.colIdxBw[a] = a;
	}
      globset.colIdxCol[0] = 466; /* red1 */
      if (globset.NA > 1)
	globset.colIdxCol[1] = 126; /* green */
      globset.colIdxBw[0] = 120;
      if (globset.NA > 1 )
	globset.colIdxBw[1] = 220;
    }
  if (globset.setdiameter && globset.NA)
    {
      for (a = 0; a < globset.NA; ++a)
	{
	  globset.sig[a] = globset.diameter;
	} 
    }
  
}
/* ========================== >>> loadAtomPos <<< ===========================*/
void loadAtomPos(void)
{
  /* DESCRIPTION:
     The file should be an ascii file of the form:
     x0 y0 z0 
     x1 y1 z1 
     .
     .
     .
     where 0 refers to atom 0 and 1 refers to atom 1*/
  FILE* ifs;
  int i, a, nf, first=1;
  char line[1024];
  long fpostmp, fpos=0;
  /*printf("Loading file %s\n", inputFile);*/
  i = 0;
  nf = 0;
  a = 0;
  add_one_frame(0);
  ifs = fopen(inputFile, "r"); 
  if (!ifs)
    {
      fprintf(stderr, "ERROR: invalid input file or wrong arguments!\n");
      exit(-1);
    }
  /* fake reading to get number of particles etc. */
  while(!feof(ifs))
    {
#if 0
      fpostmp=ftell(ifs);
#endif
      readLine(ifs, line);

      if (parseLine(line, &nf, &i, 1)) 
	continue;
      if (globset.NA && a >= globset.NA) 
	{
	  a = 0;
	  add_mol(i);
	}
#if 0  
      if (first)
	{
	  first=0;
	  fpos = fpostmp;  
	}
#endif
      /*assignAtom(nf, i, 0, line);*/
      
      /*readLine(ifs, line);*/
      /*assignAtom(nf, i, j, line);*/
      
      ++a;
      add_atom(i, a);
    }

  /*printf("NumMols:%d", NumMols);*/
  
  if (nf == 0 && globset.NumMols[0] == 0 && mols == NULL)
    {
      printf("ERROR: Presumbly the input file you supplied is void!\n");
      exit(-1);
    }
  if (globset.NA)
    globset.sig = malloc(sizeof(double)*globset.NA);
  setdefaults_after_fakeread();
  globset.frameNo = ++nf;
  globset.NumMols[nf] = i;
  /* printf("Read %i molecule\n", NumMols);
     printf("Number of frames: %d\n", frameNo);
  */
  nf = i = 0;
#if 0
  if (fseek(ifs, fpos, SEEK_SET))
    {
      perror("Error seeking the positions file:");
      exit(-1);
    }
#else
  rewind(ifs);
#endif
  while(!feof(ifs))
    {
      readLine(ifs, line);

      assignAtom(nf, i, a, line);
      
      if (parseLine(line, &nf, &i, 0)) 
	continue;
      ++a;
      if (NA && i >= NA)
	{
	  a = 0;
	  i++;
	}
    }
 fclose(ifs);
  /*
     printf("File with atoms positions successfully read.\n");*/
}


/* ======================== >>> readRGB <<< ==============================*/
void readRGB(void)
{
  int i, rc, gc, bc;
  /* determine the number of color to redad */
  mgl_col = (struct colStruct*) malloc(NUMCOLS * sizeof(struct colStruct));
  for (i=0; i < NUMCOLS; i++)
    {
      sscanf(mglrgb[i], "%d %d %d %[^\n]\n", &rc, &gc, &bc, 
	     mgl_col[i].name);
      //printf("rgb(%d, %d, %d)\n", rc, gc, bc);
      mgl_col[i].rgba[3] = 1.0;
      /* float normalized to [0.0, 1.0] */
      mgl_col[i].rgba[0] = ((float) rc) / 255.0; 
      mgl_col[i].rgba[1] = ((float) gc) / 255.0; 
      mgl_col[i].rgba[2] = ((float) bc) / 255.0; 
    }
}

/* ======================== >>> default_pars <<< ==========================*/
void default_pars(void)
{
  int i;

  globset.L = 14.0;
  globset.saved_counter = 0;
  globset saveandquit = 0;
  globset.savefile = NULL;
  globset.drawcube = 1;
  globset.sig = NULL;
  globset.height = NULL;
  globset.setdiameter = 0;
  globset.numAt = 0; /* atomi per molecola 0=illimitati a meno che non si usi .newmol*/
  globset.setdiameter = 0;
  globset.NumMols = NULL;
  globset.Width = 500;
  globset.Height = 500;
  globset.infos = 1;
  globset.frameNo = 1;
  globset.fadeMode = 2;
  globset.NA = 1;
  globset.axon = 0;
  globset.bw = 0;
  globset.degx = 0.0;
  globset.degy = 0.0;
  globset.degz = 0.0;
  globset.deginc = 5.0;
  globset.viewangle=45.0;
  globset.setvp = 0;
  globset.diameter = 1.0;
  globset.dist = 1.0;
 
  readRGB();

  setBW();

}


/* =========================== >>> rotatex <<< =========================== */
void special(int k, int x, int y)
{
  switch (k) {
  case GLUT_KEY_UP:
    globset.degx += globset.deginc;
    //printf("degx: %f\n", degx);
    break;
  case GLUT_KEY_DOWN:
    globset.degx -= globset.deginc;
    //printf("degx: %f\n", degx);
    break;
  case GLUT_KEY_LEFT:
    globset.degy -= globset.deginc;
    //printf("degy: %f\n", degy);
    break;
  case GLUT_KEY_RIGHT:
    globset.degy += globset.deginc;
    //printf("degy: %f\n", degy);
    break;
  case GLUT_KEY_PAGE_UP:
  case GLUT_KEY_F1:
    globset.dist += globset.L / 10.0;
    break;
  case GLUT_KEY_PAGE_DOWN:
  case GLUT_KEY_F2:
    globset.dist -= globset.L / 10.0;
    break;

  default:
    return;
  }
  glutPostRedisplay();
}


/* ======================== >>> degadd <<< ================================ */
void key(unsigned char k, int x, int y)
 {
   switch(k) 
     {
     case '+':
       if (globset.deginc > 179.0) break;/* Maximu increment is 180.0 degrees */
       globset.deginc += 1.0; 
       //printf("degree increment: %.0f\n", deginc);
       break;
     case '-':
       if (globset.deginc < 2.0) break; /* Minimum increment is 1.0 degree */
       globset.deginc -= 1.0;
       //printf("deg increment: %.0f\n", deginc);
       break;
     case 'a':
       globset.L += 1.0;
       //printf("L=%.0f\n", L);
       setproj();
       break;
     case 'z':
       if (globset.L > 1.0)
	 globset.L -= 1.0;
       //printf("L=%.0f\n", L);
       setproj();
       break;
     case 'i':
       globset.deginc = -globset.deginc;
       //printf("deginc inverted %.0f to %.0f\n", -deginc, deginc);
       break;
     case 'r':
       globset.degx=0.0;
       globset.degy=0.0;
       globset.degz=0.0;
       //printf("Reset rotations to (0,0,0)\n");
       break;
     case 27: // escape  
       exit(0);
       break;
     case 'p':
       globset.infos = !globset.infos;
       break;
     case 'c':
       globset.bw = !globset.bw;
       //printf("bw: %d\n", bw);     
       break;
     case 'f':
       ++globset.fadeMode;
       if (globset.fadeMode > 3) globset.fadeMode = 1;
       break;
     case 'v':
       globset.axon = !globset.axon;
       myReshape(globset.Width, globset.Height);
       break;
     case 's':
       /* save the current image here */
       save_image();
       if (globset.saved_counter < 256000) 
	 globset.saved_counter++;
     break;
     default: return;

     }
   glutPostRedisplay();
}
/*  Main Loop

 *  Open window with initial window size, title bar, 
 *  RGBA display mode, and handle input events.
 */
int main(int argc, char** argv)
{
  default_pars();
  args(argc, argv);
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowPosition(100, 100);
  glutInitWindowSize(globset.Width,globset.Height);
  glutCreateWindow("MOLGL by Cristiano De Michele (C) 1998-2004");
  myinit();
  loadAtomPos();
  buildAtomsList();
  glutDisplayFunc(display);
  glutReshapeFunc(myReshape);
  glutKeyboardFunc(key);
  glutSpecialFunc(special);
  /*  glutVisibilityFunc(visible);*/
  glutPostRedisplay();
  glutMainLoop();
  return 0;             /* ANSI C requires main to return int. */
}
