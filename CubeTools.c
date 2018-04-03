#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glib.h>
#include <glib/gprintf.h>
/* PLplot */
#include <plcdemos.h>
/* Imagemagick */
#include <wand/MagickWand.h>
/* GLIB */
/* https://www.ibm.com/developerworks/linux/tutorials/l-glib/ */
#include <glib.h>
/* GSL */
#include <gsl/gsl_histogram.h>
GSList* history = NULL;
/* SDL */
#include <SDL/SDL.h>
/* openGL */
#include <GL/gl.h>
#include <GL/glu.h>


#define a0bohr 0.529177249
#define X 0
#define Y 1
#define Z 2
#define MINI 0
#define MAXI 1
#define DATA 3
#define NCHAR 4096
#define DELIM   " \n" 
#define MAXWORD 50
#define MAXLEN  512
#define True 1
#define False 0
typedef struct {
  double* data;
  int n;
} LDOS;
typedef struct {
  char filename[NCHAR];
  int natom;
  double O[3];
  double L[3][3];
  int N[3];
  double **q;
  int *elt;
  double ****data;
  double DIM[3];
  double MC[3];
  double atom_lim[3][2];
  int       idx_atom_lim[3][2];
  double LDOS;			/* current LDOS corresponding to LDOSmap */
  double **LDOSmap;		/* 2D grid which contains the z coordinate at LDOS=Ctes */
  double **zmap;		/* 2D grid which contains the LDOS at z */
  double **ymap;		/* LDOS 2D grid x-z map at y=cte */
  double **xmap;		/* LDOS 2D grid x-z map at y=cte */
  char *pdffile;
} CubeFile;
typedef struct {
  char line[NCHAR];
  int ntokens;
  char words[MAXWORD][MAXLEN];
} Line;
/* ********************************************************************* */
void distrib(double **map,int nx,int ny,double histmin,double histmax,int histlevel,double histymin,double histymax){
  fprintf(stdout,"### histmin= %g histmax= %g nhist= %d\n",histmin,histmax,histlevel);

  gsl_histogram * h = gsl_histogram_alloc (nx*ny);
  gsl_histogram_set_ranges_uniform (h, histmin, histmax);
  /* double *data; */
  /* data=malloc(sizeof(double)); */
  int i,j;
  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
      gsl_histogram_increment (h, map[i][j]);
      /* data[k]=map[i][j]; */
      /* k++; */
      /* data=realloc(data,(k+1)*sizeof(double)); */
    }
  }
  FILE *out;
  out=fopen("hist.dat","w+");
  gsl_histogram_fprintf (out, h, "%g", "%g");
  fclose(out);
  system("xmgrace -nxy hist.dat&");
  gsl_histogram_free (h);
  /* plcol0( 1 ); */
  

  /* plenv( histmin,histmax, histymin,histymax,1, 0 ); */
  /* plhist( k, data, histmin,histmax, histlevel, PL_HIST_IGNORE_OUTLIERS ); */
  /* //plhist( k,data, histmin,histmax, histlevel, PL_HIST_NOSCALING ); */
  /* plcol0( 2 ); */
  /* pllab( "#frValue", "#frFrequency", */
  /* 	 "Distribution" ); */
  /* free(data); */
}
/* ********************************************************************* */
void plot_distribution(int n, double* data, int histmin,int histmax, int nhistlevel){
  plsdev( "pdfcairo" );
  plsfnam( "essai.pdf" ); 
  //plspal0( "cmap0_black_on_white.pal" );
  plspal1( "cmap1_gray.pal", 1 );
  plscmap0n( 3 );
  plinit();
  plcol0( 1 );
  plhist(n,data, histmin,histmax, nhistlevel, PL_HIST_IGNORE_OUTLIERS );
  plcol0( 2 );
  pllab( "#frValue", "#frFrequency","#Titre" );
  plend();
}
/* ********************************************************************* */
double **matrix(int n,int m){
  double **mat;
  int i;
  mat=malloc(n*sizeof(double*));
  for(i=0;i<n;i++) mat[i]=malloc(m*sizeof(double));
  return mat;
}
/* ********************************************************************* */
double ****tensor4(int n,int m,int l,int p){
  double ****mat;
  int i,j,k;
  mat=malloc(n*sizeof(double***));
  for(i=0;i<n;i++) {
    mat[i]=malloc(m*sizeof(double**));
    for(j=0;j<m;j++) {
      mat[i][j]=malloc(l*sizeof(double*));
      for(k=0;k<l;k++){
	mat[i][j][k]=malloc(p*sizeof(double));
      }
    }
  }
  return mat;
}
/* ********************************************************************* */
int *ivector(int n){
  int *ivec;
  ivec=malloc(n*sizeof(int));
  return ivec;
}
/* -------------------------------------------------------------------------------------------------------------- */
Line *read_line(char *line){
  Line *sline;sline=malloc(sizeof(Line));   sprintf(sline->line,"%s",line);
  char *sep=NULL;  sep=strtok(line,DELIM);

  sline->ntokens=0;
  while (sep != NULL){
    strcpy(sline->words[sline->ntokens++], sep);
    sep = strtok(NULL, DELIM);
  }

  
  free(sep);
  return sline;
}

/* ********************************************************************* */
CubeFile *NewCubeFile(CubeFile *fileini){
  
  int ix,iy,iz,i;
  CubeFile *file;file=malloc(sizeof(CubeFile));
  sprintf(file->filename,"%s","new");
  
  


  file->natom=fileini->natom;
  file->O[X]=fileini->O[X];
  file->O[Y]=fileini->O[Y];
  file->O[Z]=fileini->O[Z];
  file->N[X]=fileini->N[X];
  file->L[X][X]=fileini->L[X][X];  file->L[X][Y]=fileini->L[X][Y]; file->L[X][Z]=fileini->L[X][Z];
  file->N[Y]=fileini->N[Y];
  file->L[Y][X]=fileini->L[Y][X];  file->L[Y][Y]=fileini->L[Y][Y];  file->L[Y][Z]=fileini->L[Y][Z];
  file->N[Z]=fileini->N[Z];
  file->L[Z][X]=fileini->L[Z][X];  file->L[Z][Y]=fileini->L[Z][Y];  file->L[Z][Z]=fileini->L[Z][Z];

  file->DIM[X]=fileini->DIM[X];
  file->DIM[Y]=fileini->DIM[Y];
  file->DIM[Z]=fileini->DIM[Z];
  file->MC[X]=fileini->MC[X] ;
  file->MC[Y]=fileini->MC[Y] ;
  file->MC[Z]=fileini->MC[Z];
  
  file->q=matrix(file->natom,3);
  file->elt=ivector(file->natom);

  /* atoms */
  for(i=0;i<file->natom;i++){
    file->elt[i]=fileini->elt[i];
    file->q[i][X]=fileini->q[i][X];
    file->q[i][Y]=fileini->q[i][Y];
    file->q[i][Z]=fileini->q[i][Z];
  }
  file->data=tensor4(4,file->N[X],file->N[Y],file->N[Z]);
  for (ix=0;ix<file->N[X];ix++) {
    for (iy=0;iy<file->N[Y];iy++) {
      for (iz=0;iz<file->N[Z];iz++) {
	file->data[X][ix][iy][iz]=fileini->data[X][ix][iy][iz];
	file->data[Y][ix][iy][iz]=fileini->data[Y][ix][iy][iz];
	file->data[Z][ix][iy][iz]=fileini->data[Z][ix][iy][iz];
	file->data[DATA][ix][iy][iz]=0.0;
      }
    }
  }
  
  file->atom_lim[X][MINI]=fileini->atom_lim[X][MINI];
  file->atom_lim[Y][MINI]=  fileini->atom_lim[Y][MINI];
  file->atom_lim[Z][MINI]=  fileini->atom_lim[Z][MINI];
  file->atom_lim[X][MAXI]=  fileini->atom_lim[X][MAXI];
  file->atom_lim[Y][MAXI]=  fileini->atom_lim[Y][MAXI];
  file->atom_lim[Z][MAXI]=  fileini->atom_lim[Z][MAXI];

  file->idx_atom_lim[X][MINI]=  fileini->idx_atom_lim[X][MINI];
  file->idx_atom_lim[Y][MINI]=  fileini->idx_atom_lim[Y][MINI];
  file->idx_atom_lim[Z][MINI]=  fileini->idx_atom_lim[Z][MINI];
  file->idx_atom_lim[X][MAXI]=  fileini->idx_atom_lim[X][MAXI];
  file->idx_atom_lim[Y][MAXI]=  fileini->idx_atom_lim[Y][MAXI];
  file->idx_atom_lim[Z][MAXI]=  fileini->idx_atom_lim[Z][MAXI];

  file->LDOS=0.0;
  file->LDOSmap=matrix(file->N[X],file->N[Y]);
  file->zmap=matrix(file->N[X],file->N[Y]);
  file->ymap=matrix(file->N[X],file->N[Z]);
  file->xmap=matrix(file->N[Y],file->N[Z]);
  file->pdffile=malloc(NCHAR*sizeof(char));sprintf(file->pdffile,"%s","essai.pdf");
  
  return file;
}
/* ********************************************************************* */
CubeFile *ReadCubeFile(char *filename){
  /* WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING */
  /* WARNING                               data in cube file are in a. u. but they are in angstroem in ReadCubeFile         WARNING */
  /* WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING */



  
  int ix,iy,iz,i;
  CubeFile *file;file=malloc(sizeof(CubeFile));
  sprintf(file->filename,"%s",filename);
  
  
  FILE *in;
  in=fopen(file->filename,"r");   if(in==NULL){    fprintf(stdout,"### ERR ### impossible to open %s\n",filename);    exit(0);  }
  char line[NCHAR];
  //cdump=
    fgets(line,NCHAR,in);    fprintf(stdout,"### MSG ### %s\n",line);
    //cdump=
    fgets(line,NCHAR,in);    fprintf(stdout,"### MSG ### %s\n",line);
  fscanf(in,"%d %lg %lg %lg\n",&file->natom,&file->O[X],&file->O[Y],&file->O[Z]);
  file->O[X]*=a0bohr;    file->O[Y]*=a0bohr;    file->O[Z]*=a0bohr;
  fprintf(stdout,"# Origin @ (%g,%g,%g)\n",file->O[X],file->O[Y],file->O[Z]);
  fscanf(in,"%d %lg %lg %lg\n",&file->N[X],&file->L[X][X],&file->L[X][Y],&file->L[X][Z]);
  fscanf(in,"%d %lg %lg %lg\n",&file->N[Y],&file->L[Y][X],&file->L[Y][Y],&file->L[Y][Z]);
  fscanf(in,"%d %lg %lg %lg\n",&file->N[Z],&file->L[Z][X],&file->L[Z][Y],&file->L[Z][Z]);
  file->L[X][X]*=a0bohr;    file->L[X][Y]*=a0bohr;    file->L[X][Z]*=a0bohr;
  file->L[Y][X]*=a0bohr;    file->L[Y][Y]*=a0bohr;    file->L[Y][Z]*=a0bohr;
  file->L[Z][X]*=a0bohr;    file->L[Z][Y]*=a0bohr;    file->L[Z][Z]*=a0bohr;

  file->DIM[X]=file->N[X]*file->L[X][X];
  file->DIM[Y]=file->N[Y]*file->L[Y][Y];
  file->DIM[Z]=file->N[Z]*file->L[Z][Z];
  
  file->q=matrix(file->natom,3);
  file->elt=ivector(file->natom);
  file->MC[X]=0.0 ; file->MC[Y]=0.0 ; file->MC[Z]=0.0;
  for(i=0;i<file->natom;i++){
    fscanf(in,"%d %*g %lg %lg %lg\n",&file->elt[i],&file->q[i][X],&file->q[i][Y],&file->q[i][Z]);
    file->q[i][X]*=a0bohr;file->q[i][Y]*=a0bohr;file->q[i][Z]*=a0bohr;
    file->MC[X]+=file->q[i][X];    file->MC[Y]+=file->q[i][Y];    file->MC[Z]+=file->q[i][Z];
  }
  file->atom_lim[X][MINI]=file->q[0][X];
  file->atom_lim[Y][MINI]=file->q[0][Y];
  file->atom_lim[Z][MINI]=file->q[0][Z];
  file->atom_lim[X][MAXI]=file->q[0][X];
  file->atom_lim[Y][MAXI]=file->q[0][Y];
  file->atom_lim[Z][MAXI]=file->q[0][Z];
  for(i=1;i<file->natom;i++){
    if(file->atom_lim[X][MINI]>file->q[i][X]) file->atom_lim[X][MINI]=file->q[i][X] ;
    if(file->atom_lim[Y][MINI]>file->q[i][Y]) file->atom_lim[Y][MINI]=file->q[i][Y] ;
    if(file->atom_lim[Z][MINI]>file->q[i][Z]) file->atom_lim[Z][MINI]=file->q[i][Z] ;
    if(file->atom_lim[X][MAXI]<file->q[i][X]) file->atom_lim[X][MAXI]=file->q[i][X] ;
    if(file->atom_lim[Y][MAXI]<file->q[i][Y]) file->atom_lim[Y][MAXI]=file->q[i][Y] ;
    if(file->atom_lim[Z][MAXI]<file->q[i][Z]) file->atom_lim[Z][MAXI]=file->q[i][Z] ;
  }

  file->MC[X]/=file->natom;    file->MC[Y]/=file->natom;      file->MC[Z]/=file->natom;



  file->data=tensor4(4,file->N[X],file->N[Y],file->N[Z]);
  
  for (ix=0;ix<file->N[X];ix++) {
    for (iy=0;iy<file->N[Y];iy++) {
      for (iz=0;iz<file->N[Z];iz++) {
	fscanf(in,"%lg ",&file->data[DATA][ix][iy][iz]);
	file->data[X][ix][iy][iz]=file->O[X]+ix*file->L[X][X];
	file->data[Y][ix][iy][iz]=file->O[Y]+iy*file->L[Y][Y];
	file->data[Z][ix][iy][iz]=file->O[Z]+iz*file->L[Z][Z];
	/* file->data[X][ix][iy][iz]=file->O[X]-.5*file->N[X]*file->L[X][X]+ix*file->L[X][X]; */
	/* file->data[Y][ix][iy][iz]=file->O[Y]-.5*file->N[Y]*file->L[Y][Y]+iy*file->L[Y][Y]; */
	/* file->data[Z][ix][iy][iz]=file->O[Z]-.5*file->N[Z]*file->L[Z][Z]+iz*file->L[Z][Z]; */
	//	fprintf(stdout,"%g %g %g %g\n ",file->data[X][ix][iy][iz],file->data[Y][ix][iy][iz],file->data[Z][ix][iy][iz],file->data[DATA][ix][iy][iz]);

      }
    }
  }

  file->idx_atom_lim[X][MINI]=0;  file->idx_atom_lim[Y][MINI]=0;  file->idx_atom_lim[Z][MINI]=0;
  file->idx_atom_lim[X][MAXI]=file->N[X]-1;  file->idx_atom_lim[Y][MAXI]=file->N[Y]-1;  file->idx_atom_lim[Z][MAXI]=file->N[Z]-1;
  ix=1;while(file->atom_lim[X][MINI]>file->data[X][ix][0][0]) ix++; file->idx_atom_lim[X][MINI]=ix-1;
  ix=1;while(file->atom_lim[Y][MINI]>file->data[Y][0][ix][0]) ix++; file->idx_atom_lim[Y][MINI]=ix-1;
  ix=1;while(file->atom_lim[Z][MINI]>file->data[Z][0][0][ix]) ix++; file->idx_atom_lim[Z][MINI]=ix-1;
  ix=file->N[X]-2;while(file->atom_lim[X][MAXI]<file->data[X][ix][0][0]) ix--; file->idx_atom_lim[X][MAXI]=ix+1;
  ix=file->N[Y]-2;while(file->atom_lim[Y][MAXI]<file->data[Y][0][ix][0]) ix--; file->idx_atom_lim[Y][MAXI]=ix+1;
  ix=file->N[Z]-2;while(file->atom_lim[Z][MAXI]<file->data[Z][0][0][ix]) ix--; file->idx_atom_lim[Z][MAXI]=ix+1;


  fclose(in);
  file->LDOS=0.0;
  file->LDOSmap=matrix(file->N[X],file->N[Y]);
  file->zmap=matrix(file->N[X],file->N[Y]);
  file->ymap=matrix(file->N[X],file->N[Z]);
  file->xmap=matrix(file->N[Y],file->N[Z]);
  file->pdffile=malloc(NCHAR*sizeof(char));sprintf(file->pdffile,"%s","essai.pdf");

  
  return file;
}
PLFLT tr[6] ={ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 };
// pltr_data is unused so mark it with the PL_UNUSED macro
static void mypltr( PLFLT x, PLFLT y, PLFLT *tx, PLFLT *ty, void * PL_UNUSED( pltr_data ) ) {
    *tx = tr[0] * x + tr[1] * y + tr[2];
    *ty = tr[3] * x + tr[4] * y + tr[5];
}
/* -------------------------------------------------------------------------------------------------------------- */
void WriteCubeFile(CubeFile *file,char *filename){
  FILE *out;out=fopen(filename,"w+");
  int i,ix,iy,iz;
  fprintf(out,"comment\ncomment\n");
  fprintf(out,"%5d%12.6f%12.6f%12.6f\n",file->natom,file->O[X]/a0bohr,file->O[Y]/a0bohr,file->O[Z]/a0bohr);
  fprintf(out,"%5d%12.6f%12.6f%12.6f\n",file->N[X],file->L[X][X]/a0bohr,file->L[X][Y]/a0bohr,file->L[X][Z]/a0bohr);
  fprintf(out,"%5d%12.6f%12.6f%12.6f\n",file->N[Y],file->L[Y][X]/a0bohr,file->L[Y][Y]/a0bohr,file->L[Y][Z]/a0bohr);
  fprintf(out,"%5d%12.6f%12.6f%12.6f\n",file->N[Z],file->L[Z][X]/a0bohr,file->L[Z][Y]/a0bohr,file->L[Z][Z]/a0bohr);
  for(i=0;i<file->natom;i++){
    fprintf(out,"%5d%12.6f%12.6f%12.6f%12.6f\n",file->elt[i],file->elt[i]*1.0,file->q[i][X]/a0bohr,file->q[i][Y]/a0bohr,file->q[i][Z]/a0bohr);
  }
  for (ix=0;ix<file->N[X];ix++) {
    for (iy=0;iy<file->N[Y];iy++) {
      for (iz=0;iz<file->N[Z];iz++) {
	fprintf(out,"%13.5E",file->data[DATA][ix][iy][iz]);
	if (iz % 6 == 5)
	  fprintf(out,"\n");
      }
      fprintf(out,"\n");
    }
  }
  
  
  fclose(out);
}
/* -------------------------------------------------------------------------------------------------------------- */
void execute(char *filename){
  CubeFile **file;
  file=malloc(sizeof(CubeFile));
  


  FILE *param;
  char * line = NULL;
  int i;i=0;
  param=fopen(filename,"r");
  while(feof(param)==0){
    line=malloc(NCHAR*sizeof(char));
    fgets(line,NCHAR,param);

    Line *sline=NULL;
    sline=read_line(line);
    fprintf(stdout,"ntokens=%d\n",sline->ntokens);
    if(strcmp(sline->words[0],"read")==0)  {
      file=realloc(file,(i+1)*sizeof(CubeFile));
      fprintf(stdout,"\treading %s...\n",sline->words[1]);
      file[i]=ReadCubeFile(sline->words[1]);
      i++;
    }
    
    free(line);
    free(sline);
  }
  fclose(param);

}
/* -------------------------------------------------------------------------------------------------------------- */
char * read_line_from_console(void) {
  char * line = malloc(100), * linep = line;
  size_t lenmax = 100, len = lenmax;
  int c;
  
  if(line == NULL)
    return NULL;
  
  for(;;) {
    c = fgetc(stdin);
    if(c == EOF)
      break;
    
    if(--len == 0) {
      len = lenmax;
      char * linen = realloc(linep, lenmax *= 2);
      
      if(linen == NULL) {
	free(linep);
	return NULL;
      }
      line = linen + (line - linep);
      linep = linen;
    }
    
    if((*line++ = c) == '\n')
      break;
  }
  *line = '\0';
  return linep;
}
/* ******************************************************************************* */
static void setup_light( CubeFile *file ) {
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  /* glEnable(GL_LIGHT1); */
  //glMatrixMode(GL_MODELVIEW);
  //glLoadIdentity( );
  float lightpos[] = {1.0,1.0,1.0, 1.};
  glLightfv(GL_LIGHT0, GL_POSITION, lightpos);
  static float blanc[] = { 1.0F,1.0F,1.0F,1.0F };
  static float bleu[] = { 0.0F,0.0F,1.0F,1.0F };
  static float gris[] = { 0.25F,0.25F,0.25F,1.0F };
  glLightfv(GL_LIGHT0,GL_DIFFUSE,gris);
  glEnable(GL_AUTO_NORMAL);
  glEnable(GL_NORMALIZE);
  glDepthFunc(GL_LESS);
  glEnable(GL_DEPTH_TEST);
  glMaterialf(GL_FRONT,GL_SHININESS,32.0F);
  /* float lightpos2[] = {file->MC[X],file->MC[Y],file->MC[Z], 1.}; */
  /* glLightfv(GL_LIGHT1, GL_POSITION, lightpos2); */

}

/* ******************************************************************************* */
static void setup_opengl( int width, int height ) {
  
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_COLOR_MATERIAL);

  float ratio = (float) width / (float) height;
  /* Our shading model--Gouraud (smooth). */
  glShadeModel( GL_SMOOTH );
  
  /* Culling. */
  glCullFace( GL_BACK );
  glFrontFace( GL_CCW );
  glEnable( GL_CULL_FACE );
  
  /* Set the clear color. */
  glClearColor( 1.0 , 0, 0, 0 );
  
  /* Setup our viewport. */
  glViewport( 0, 0, width, height );
  
  /*
   * Change to the projection matrix and set
   * our viewing volume.
   */
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity( );
  /*
   * EXERCISE:
   * Replace this with a call to glFrustum.
   */
  gluPerspective( 70.0, ratio, 1.0, 1024.0 );
}
/* ******************************************************************************* */
static void unit_cell(CubeFile *file,GLfloat tx,GLfloat ty,GLfloat tz) {
  glPushMatrix() ;
  /* Move down the z-axis. */
  //  glTranslatef( tx,ty,tz);
  /* Send our triangle data to the pipeline. */
  glBegin( GL_LINES );
  glVertex3f(0.0f, 0.0f, 0.0f);
  glVertex3f(file->DIM[X], 0.0f, 0.0f);

  glVertex3f(0.0f, 0.0f, 0.0f);
  glVertex3f(0.0f, file->DIM[Y], 0.0f);

  glVertex3f(0.0f, 0.0f, 0.0f);
  glVertex3f(0.0f, 0.0f, file->DIM[Z]);

  glVertex3f(file->DIM[X], 0.0f, 0.0f);
  glVertex3f(file->DIM[X], 0.0f, file->DIM[Z]);

  glVertex3f(file->DIM[X], 0.0f, file->DIM[Z]);
  glVertex3f(0.0f, 0.0f, file->DIM[Z]);

  glVertex3f(0.0f, file->DIM[Y], 0.0f);
  glVertex3f(0.0f, file->DIM[Y], file->DIM[Z]);

  glVertex3f(0.0f, file->DIM[Y], file->DIM[Z]);
  glVertex3f(0.0f, 0.0f, file->DIM[Z]);

  glVertex3f(0.0f, file->DIM[Y], 0.0f);
  glVertex3f(file->DIM[X], file->DIM[Y], 0.0f);
  
  glVertex3f(file->DIM[X], file->DIM[Y], 0.0f);
  glVertex3f(file->DIM[X], file->DIM[Y], file->DIM[Z]);

  glVertex3f(file->DIM[X], file->DIM[Y], file->DIM[Z]);
  glVertex3f(0.0f, file->DIM[Y], file->DIM[Z]);

  glVertex3f(file->DIM[X], file->DIM[Y], file->DIM[Z]);
  glVertex3f(file->DIM[X], 0.0f, file->DIM[Z]);

  glVertex3f(file->DIM[X], file->DIM[Y], 0.0f);
  glVertex3f(file->DIM[X], 0.0f, 0.0f);
  
  glEnd();
  glPopMatrix();
}
/* ******************************************************************************* */
static void atom(GLfloat tx,GLfloat ty,GLfloat tz) {
  /* We don't want to modify the projection matrix. */
  //glMatrixMode( GL_MODELVIEW );
  //glLoadIdentity( );

  glPushMatrix() ;
  /* Move down the z-axis. */
  glTranslatef( tx,ty,tz);
  /* Send our triangle data to the pipeline. */
  GLfloat radius=.7;
  GLUquadricObj *quad= gluNewQuadric();
  gluSphere( quad , radius , 36 , 18 );
  glPopMatrix();
}
/* ******************************************************************************* */
static void draw_screen( CubeFile *file) {
  /* Clear the color and depth buffers. */
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  int i;
  fprintf(stdout,"draw\n");




  /* glPushMatrix() ; */
  /* /\* gluLookAt(file->MC[X]+1,file->MC[Y]+1,file->MC[Z], *\/ */
  /* /\* 	    file->MC[X],file->MC[Y],file->MC[Z], *\/ */
  /* /\* 	    0,1,0); *\/ */
  /* glMatrixMode( GL_MODELVIEW );  */
  /* glLoadIdentity( ); */
  /* gluPerspective( 90.0, 1.0, 1.0, 1024.0 ); */
  /* gluLookAt(3,4,4, */
  /* 	    0,0,0, */
  /* 	    0,0,1); */
  gluLookAt(10,50,40,
  	    0,0,30,
  	    0,0,1.5);
  	    
  	    
  /* glPopMatrix(); */


  unit_cell(file,0.0,0.0,0.0);


  for(i=0;i<file->natom;i++)  atom(file->q[i][X],file->q[i][Y],file->q[i][Z]);
  
  //atom(0.0,0.0,0.0);

  glPushMatrix() ;

  /* glClear(GL_COLOR_BUFFER_BIT); */
  /* glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); */
  /* Move down the z-axis. */
  //  glTranslatef( tx,ty,tz);
  /* Send our triangle data to the pipeline. */
  glEnable(GL_BLEND); glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 
  glBegin( GL_POLYGON );
  glColor4d(0.5,0.5,0.5,0.5);
  glVertex3f(0.0f, 0.0f, file->data[Z][0][0][270]-file->data[Z][0][0][0]);
  glColor4d(0.5,0.5,0.0,0.5);
  glVertex3f(file->DIM[X], 0.0f, file->data[Z][0][0][270]);
  glVertex3f(file->DIM[X], file->DIM[Y], file->data[Z][0][0][270]);
  glVertex3f(0.0f, file->DIM[Y], file->data[Z][0][0][270]);
  glEnd();
  glDisable(GL_BLEND);
  glPopMatrix();
  fprintf(stdout,"%g %g %g\n",file->data[Z][0][0][0],file->data[Z][0][0][270],file->data[Z][0][0][449]);


  
  glFlush(); 



  
  /*
  * Swap the buffers. This this tells the driver to
  * render the next frame from the contents of the
  * back-buffer, and to set all rendering operations
  * to occur on what was the front-buffer.
  *
  * Double buffering prevents nasty visual tearing
  * from the application drawing on areas of the
  * screen that are being updated at the same time.
  */
  SDL_GL_SwapBuffers( );
  
}

/* -------------------------------------------------------------------------------------------------------------- */
int main(int nargument,char **argument) {
  

  if (g_file_test (argument[1], G_FILE_TEST_EXISTS)) {
    g_printf ("Reading file %s ...\n", argument[1]);
    execute(argument[1]);
  } else {
    g_printf ("Cannot open file %s. Command line mode\n", argument[1]);
    char *line;
    int nfile=0;    CubeFile **file;    file=malloc(sizeof(CubeFile));
    int nLDOSfile=0;    LDOS *LDOSfile; LDOSfile=malloc(sizeof(LDOS));
    for(;;){
      printf("CubeTools > ");
      line=read_line_from_console();
      history = g_slist_append(history, line);
      printf("%s",line);
      if(strcmp(line,"quit\n")==0) break;
      else {

	Line *sline=NULL;
	sline=read_line(line);
	fprintf(stdout,"ntokens=%d\n",sline->ntokens);
	/* ################################################################################# */
	if(strcmp(sline->words[0],"povray")==0)  {
	  /* ################################################################################# */
	  int idxfile=0;
	  FILE *out;
	  out=fopen("CubeTools.pov","w+");

	  fprintf(out,"global_settings {\n");
	  fprintf(out,"  ambient_light rgb <0.200000002980232, 0.200000002980232, 0.200000002980232>\n");
	  fprintf(out,"  max_trace_level 15\n");
	  fprintf(out,"}\n");
	  
	  fprintf(out,"background { color rgb <1,1,1> }\n");
	  fprintf(out,"#declare DY1 = 8;\n");
	  fprintf(out,"#declare DY2 = 2;\n");
	  fprintf(out,"#declare DZ = 3;\n");
	  fprintf(out,"#declare ALPHA = .7;\n");
	  fprintf(out,"camera {\n");
	    
	  double camera[3]={-2,-25,10};
	  int i;
	  for(i=0;i<sline->ntokens;i++){
	    if(strcmp(sline->words[i],"camera")==0) 	{
	      camera[X]=atof(sline->words[i+1]);
	      camera[Y]=atof(sline->words[i+2]);
	      camera[Z]=atof(sline->words[i+3]);
	    }
	  }


	  double angle=40;
	  
	  double direction[3];
	  direction[X]=file[idxfile]->MC[X]-camera[X];
	  direction[Y]=file[idxfile]->MC[Y]-camera[Y];
	  direction[Z]=file[idxfile]->MC[Z]-camera[Z];
	  double look_at[3];
	  look_at[X]=file[idxfile]->MC[X];
	  look_at[Y]=file[idxfile]->MC[Y];
	  look_at[Z]=file[idxfile]->MC[Z];
	
	  fprintf(out," location <%g,%g,%g>\n",camera[X],camera[Y],camera[Z]);
	  fprintf(out,"  angle %g\n",angle);
	  fprintf(out,"  look_at <%g,%g,%g> \n",look_at[X],look_at[Y],look_at[Z]);
	  fprintf(out,"}\n");

	  double light1[3]={100.0,0.0,0.0};
	  double light2[3]={100.0,0.0,0.0};
	  double light3[3]={100.0,0.0,0.0};
	  for(i=0;i<sline->ntokens;i++){
	    if(strcmp(sline->words[i],"light1")==0) 	{
	      light1[X]=atof(sline->words[i+1]);
	      light1[Y]=atof(sline->words[i+2]);
	      light1[Z]=atof(sline->words[i+3]);
	    }
	    if(strcmp(sline->words[i],"light2")==0) 	{
	      light2[X]=atof(sline->words[i+1]);
	      light2[Y]=atof(sline->words[i+2]);
	      light1[Z]=atof(sline->words[i+3]);
	    }
	    if(strcmp(sline->words[i],"light3")==0) 	{
	      light3[X]=atof(sline->words[i+1]);
	      light3[Y]=atof(sline->words[i+2]);
	      light3[Z]=atof(sline->words[i+3]);
	    }
	  }
	  fprintf(out,"light_source {\n");
	  fprintf(out,"  <%g,%g,%g>\n",light1[X],light1[Y],light1[Z]);
	  fprintf(out,"  color rgb <1, 1, 1>\n");
	  fprintf(out,"  fade_distance 171.046048057225\n");
	  fprintf(out,"  fade_power 0\n");
	  fprintf(out,"  parallel\n");
	  fprintf(out,"  point_at <%g,%g,%g>\n",file[idxfile]->MC[X],file[idxfile]->MC[Y],file[idxfile]->MC[Z]);
	  fprintf(out,"}\n");
	  fprintf(out,"light_source {\n");
	  fprintf(out,"  <%g,%g,%g>\n",light2[X],light2[Y],light2[Z]);
	  fprintf(out,"  color rgb <1, 1, 1>\n");
	  fprintf(out,"  fade_distance 171.046048057225\n");
	  fprintf(out,"  fade_power 0\n");
	  fprintf(out,"  parallel\n");
	  fprintf(out,"  point_at <%g,%g,%g>\n",file[idxfile]->MC[X],file[idxfile]->MC[Y],file[idxfile]->MC[Z]);
	  fprintf(out,"}\n");
	  fprintf(out,"light_source {\n");
	  fprintf(out,"  <%g,%g,%g>\n",light3[X],light3[Y],light3[Z]);
	  fprintf(out,"  color rgb <1, 1, 1>\n");
	  fprintf(out,"  fade_distance 171.046048057225\n");
	  fprintf(out,"  fade_power 0\n");
	  fprintf(out,"  parallel\n");
	  fprintf(out,"  point_at <%g,%g,%g>\n",file[idxfile]->MC[X],file[idxfile]->MC[Y],file[idxfile]->MC[Z]);
	  fprintf(out,"}\n");
	  
	  /* fprintf(out,"light_source {\n"); */
	  /* fprintf(out,"  <-84.098973324669, -13.6489285743175, 53.2414161510446>\n"); */
	  /* fprintf(out,"  color rgb <0.300000011920929, 0.300000011920929, 0.300000011920929>\n"); */
	  /* fprintf(out,"  fade_distance 171.046048057225\n"); */
	  /* fprintf(out,"  fade_power 0\n"); */
	  /* fprintf(out,"  parallel\n"); */
	  /* fprintf(out,"  point_at <84.098973324669, 13.6489285743175, -53.2414161510446>\n"); */
	  /* fprintf(out,"}\n"); */


	  /* fprintf(out,"light_source {\n"); */
	  /* fprintf(out,"  <-2,-25,10>\n"); */
	  /* fprintf(out,"  color rgb <0.300000011920929, 0.300000011920929, 0.300000011920929>\n"); */
	  /* fprintf(out,"  fade_distance 171.046048057225\n"); */
	  /* fprintf(out,"  fade_power 0\n"); */
	  /* fprintf(out,"  parallel\n"); */
	  /* fprintf(out,"  point_at <%g,%g,%g>\n",look_at[0],look_at[1],look_at[2]); */
	  /* fprintf(out,"}\n"); */

        
	  fprintf(out,"#default {\n");
	  fprintf(out,"  finish {ambient .8 diffuse 1 specular 1 roughness .005 metallic 0.5}\n");
	  fprintf(out,"}\n");

	  fprintf(out,"//---------------------------- macro Vector(Start,End,Radius)\n");
	  fprintf(out,"#macro Vector(P_s,P_e, Rv)\n");
	  fprintf(out,"  union{\n");
	  fprintf(out,"    cylinder{ P_s, P_e - (vnormalize(P_e - P_s)*9.5*Rv), Rv }\n");
	  fprintf(out,"    cone    { P_e - (vnormalize(P_e - P_s)*10*Rv), 3*Rv, P_e,0}\n");
	  fprintf(out,"  }// end of union\n");
	  fprintf(out,"#end //----------------------------------------- end of macro\n");
	  /* fprintf(out,"#declare Vector_plus_Texture =\n"); */
	  /* fprintf(out,"  texture{ pigment{ color rgb<0.2,0.5,0.0>}\n"); */
	  /* fprintf(out,"    finish { phong 1} \n"); */
	  /* fprintf(out,"  }\n"); */
	  /* fprintf(out,"#declare Vector_moins_Texture =\n"); */
	  /* fprintf(out,"  texture{ pigment{ color rgb<1.0,0.,0.0>}\n"); */
	  /* fprintf(out,"    finish { phong 1} \n"); */
	  /* fprintf(out,"  }\n"); */


	  
	  /* fprintf(out,"//---------------------------- macro Vector(Start,End,Radius)\n"); */
	  /* fprintf(out,"#macro Vector(P_s,P_e, Rv)\n"); */
	  /* fprintf(out,"  union{\n"); */
	  /* fprintf(out,"    cylinder{ P_s, P_e - (vnormalize(P_e - P_s)*9.5*Rv), Rv }\n"); */
	  /* fprintf(out,"    cone    { P_e - (vnormalize(P_e - P_s)*10*Rv), 3*Rv, P_e,0}\n"); */
	  /* fprintf(out,"  }// end of union\n"); */
	  /* fprintf(out,"#end //----------------------------------------- end of macro\n"); */
	  fprintf(out,"#declare Vector_Texture_X =\n");
	  fprintf(out,"  texture{ pigment{ color rgb<1.0,0.0,0.0>}\n");
	  fprintf(out,"    finish { phong 1} \n");
	  fprintf(out,"  }\n");
	  fprintf(out,"#declare Vector_Texture_Y =\n");
	  fprintf(out,"  texture{ pigment{ color rgb<0.0,1.0,0.0>}\n");
	  fprintf(out,"    finish { phong 1} \n");
	  fprintf(out,"  }\n");
	  fprintf(out,"#declare Vector_Texture_Z =\n");
	  fprintf(out,"  texture{ pigment{ color rgb<0.0,0.0,1.0>}\n");
	  fprintf(out,"    finish { phong 1} \n");
	  fprintf(out,"  }\n");


	  fprintf(out,"object{\n");
	  fprintf(out,"Vector(\n");
	  fprintf(out,"<%g,%g,%g>,\n",file[idxfile]->O[X],file[idxfile]->O[Y],file[idxfile]->O[Z]);
	  fprintf(out,"<%g,%g,%g>,\n",file[idxfile]->O[X]+file[idxfile]->DIM[X],file[idxfile]->O[Y],file[idxfile]->O[Z]);
	  fprintf(out,"0.05)\n");
	  fprintf(out,"texture{ Vector_Texture_X }\n");
	  fprintf(out,"}\n");
	  fprintf(out,"object{\n");
	  fprintf(out,"Vector(\n");
	  fprintf(out,"<%g,%g,%g>,\n",file[idxfile]->O[X],file[idxfile]->O[Y],file[idxfile]->O[Z]);
	  fprintf(out,"<%g,%g,%g>,\n",file[idxfile]->O[X],file[idxfile]->O[Y]+file[idxfile]->DIM[Y],file[idxfile]->O[Z]);
	  fprintf(out,"0.05)\n");
	  fprintf(out,"texture{ Vector_Texture_Y }\n");
	  fprintf(out,"}\n");
	  fprintf(out,"object{\n");
	  fprintf(out,"Vector(\n");
	  fprintf(out,"<%g,%g,%g>,\n",file[idxfile]->O[X],file[idxfile]->O[Y],file[idxfile]->O[Z]);
	  fprintf(out,"<%g,%g,%g>,\n",file[idxfile]->O[X],file[idxfile]->O[Y],file[idxfile]->O[Z]+file[idxfile]->DIM[Z]);
	  fprintf(out,"0.05)\n");
	  fprintf(out,"texture{ Vector_Texture_Z }\n");
	  fprintf(out,"}\n");



	  double re=1.6,d;
	  double resqr=re*re;
	  int j;
	  for(i=0;i<file[idxfile]->natom;i++){
	    for(j=0;j<file[idxfile]->natom;j++){
	      if (j>i){
		d=pow(file[idxfile]->q[i][X]-file[idxfile]->q[j][X],2)+pow(file[idxfile]->q[i][Y]-file[idxfile]->q[j][Y],2)+pow(file[idxfile]->q[i][Z]-file[idxfile]->q[j][Z],2);
		if (d<= resqr){
		  fprintf(out,"cylinder {\n");
		  fprintf(out,"<%f,%f,%f>,\n",file[idxfile]->q[i][X],file[idxfile]->q[i][Y],file[idxfile]->q[i][Z]);
		  fprintf(out,"<%f,%f,%f>,\n",file[idxfile]->q[j][X],file[idxfile]->q[j][Y],file[idxfile]->q[j][Z]);
		  fprintf(out,"0.1\n");
		  fprintf(out,"pigment { rgbt <0.75, 0.75, 0.75, 0> }\n");
		  fprintf(out,"}\n");
		}
	      }
	    }
	  }
	    
	  
	  double color[3]={0.4,0.4,0.4},radius=0.51;  // C
	  for(i=0;i<file[idxfile]->natom;i++){
            fprintf(out,"sphere { \n");
	    fprintf(out,"    <%g,%g,%g>, %g\n",file[idxfile]->q[i][X],file[idxfile]->q[i][Y],file[idxfile]->q[i][Z],radius);
	    fprintf(out,"    pigment { rgbt <%g,%g,%g,0> }\n",color[0],color[1],color[2]);
	    fprintf(out,"}\n");

	  }




	  
	  double P1[3]={0,0,25};
	  double P2[3]={20,0,25};
	  double P3[3]={0,20,25};
	  double P4[3]={20,20,25};
	  double transpa=0.5;
	  double colo[3]={0.5,0.5,0.5};
	  for(i=0;i<sline->ntokens;i++){
	    if(strcmp(sline->words[i],"transpa")==0) 		      transpa=atof(sline->words[i+1]);
	    if(strcmp(sline->words[i],"color")==0){
	      colo[0]=atof(sline->words[i+1]);
	      colo[1]=atof(sline->words[i+2]);
	      colo[2]=atof(sline->words[i+3]);
	    }
	    if(strcmp(sline->words[i],"plane")==0) 	{
	      P1[X]=atof(sline->words[i+1]);
	      P1[Y]=atof(sline->words[i+2]);
	      P1[Z]=atof(sline->words[i+3]);
	      P2[X]=atof(sline->words[i+4]);
	      P2[Y]=atof(sline->words[i+5]);
	      P2[Z]=atof(sline->words[i+6]);
	      P3[X]=atof(sline->words[i+7]);
	      P3[Y]=atof(sline->words[i+8]);
	      P3[Z]=atof(sline->words[i+9]);
	      P4[X]=atof(sline->words[i+10]);
	      P4[Y]=atof(sline->words[i+11]);
	      P4[Z]=atof(sline->words[i+12]);

	    }
	  }

	  fprintf(out,"triangle{\n");
	  fprintf(out,"  <%g,%g,%g>,\n",P1[X],P1[Y],P1[Z]);
	  fprintf(out,"  <%g,%g,%g>,\n",P2[X],P2[Y],P2[Z]);
	  fprintf(out,"  <%g,%g,%g>\n",P3[X],P3[Y],P3[Z]);
	  fprintf(out,"  texture{pigment{color rgbt <%g,%g,%g,%g>}finish {phong 1}}\n",colo[0],colo[1],colo[2],transpa);
	  fprintf(out,"}\n");
	  fprintf(out,"triangle{\n");
	  fprintf(out,"  <%g,%g,%g>,\n",P1[X],P1[Y],P1[Z]);
	  fprintf(out,"  <%g,%g,%g>,\n",P3[X],P3[Y],P3[Z]);
	  fprintf(out,"  <%g,%g,%g>\n",P4[X],P4[Y],P4[Z]);
	  fprintf(out,"  texture{pigment{color rgbt <%g,%g,%g,%g>}finish {phong 1}}\n",colo[0],colo[1],colo[2],transpa);
	  fprintf(out,"}\n");
	  

	  fclose(out);
	  system("povray -D CubeTools.pov ; display CubeTools.png &");
	}
	  /* ################################################################################# */
	if(strcmp(sline->words[0],"read")==0)  {
	  /* ################################################################################# */
	  /* 
	     Read the data files. Two kinds of files
	      - cube file
	      - LDOS file (output from QE)
	   */
	  if(strcmp(sline->words[1],"cube")==0)  {
	    file=realloc(file,(nfile+1)*sizeof(CubeFile));
	    file[nfile]=ReadCubeFile(sline->words[2]);
	    printf("### Cubefile with index %d created\n",nfile);
	    nfile++;
	  }
	  if(strcmp(sline->words[1],"LDOS")==0)  {
	    LDOSfile=realloc(LDOSfile,(nLDOSfile+1)*sizeof(LDOS));
	    LDOSfile[nLDOSfile].data=malloc(sizeof(double));
	    LDOSfile[nLDOSfile].n=0;
	    FILE *in;
	    in=fopen(sline->words[2],"r");
	    while(feof(in)==0){
	      LDOSfile[nLDOSfile].data=realloc(LDOSfile[nLDOSfile].data,(LDOSfile[nLDOSfile].n+1)*sizeof(double));
	      fscanf(in,"%lg ",&LDOSfile[nLDOSfile].data[LDOSfile[nLDOSfile].n]);
	      //	      fprintf(stdout,"%d %g \n",LDOSfile[nLDOSfile].n,LDOSfile[nLDOSfile].data[LDOSfile[nLDOSfile].n]);
	      LDOSfile[nLDOSfile].n++;
	    }
	    fclose(in);
	    fprintf(stdout,"### LDOS file %s (idx %d)--> %d records\n",sline->words[2],nLDOSfile,LDOSfile[nLDOSfile].n);

	    nLDOSfile++;
	  }

	}
	/* ################################################################################# */
	if(strcmp(sline->words[0],"pdfname")==0)  {
	  /* ################################################################################# */
	  sprintf(file[atoi(sline->words[1])]->pdffile,"%s",sline->words[2]);
	}
	/* ############## ########################################## */
	if(strcmp(sline->words[0],"plot")==0)  {
	  /* ################################################################################# */
	  if(strcmp(sline->words[1],"LDOS")==0)  {
	    int idxfile=0,i;
	    for(i=0;i<sline->ntokens;i++){
	      if(strcmp(sline->words[i],"file")==0) 	  idxfile=atoi(sline->words[i+1]);
	    }
	    int nhistlevel=100;
	    double histmin,histmax;
	    histmin=LDOSfile[idxfile].data[0];
	    histmax=LDOSfile[idxfile].data[LDOSfile[idxfile].n-1];
	    fprintf(stdout,"### min= %g max= %g\n",histmin,histmax);
	    for(i=0;i<sline->ntokens;i++){
	      if(strcmp(sline->words[i],"histmax")==0) 	  histmax=atof(sline->words[i+1]);
	      if(strcmp(sline->words[i],"histmin")==0) 	  histmin=atof(sline->words[i+1]);
	      if(strcmp(sline->words[i],"histlevel")==0) 	  nhistlevel=atoi(sline->words[i+1]);
	    }
	    
	    plot_distribution(LDOSfile[idxfile].n, LDOSfile[idxfile].data, histmin,histmax, nhistlevel);
	  }
	}
	/* ############## ########################################## */
	if(strcmp(sline->words[0],"info_file")==0)  {
	/* ################################################################################# */
	  int idxfile; idxfile=atoi(sline->words[1]);
	  fprintf(stdout,"### NX=%d\nNY=%d\nNZ=%d\n",file[idxfile]->N[X],file[idxfile]->N[Y],file[idxfile]->N[Z]);
	  fprintf(stdout,"### DIMX=%g\n### DIMY=%g\n### DIMZ=%g\n",file[idxfile]->DIM[X],file[idxfile]->DIM[Y],file[idxfile]->DIM[Z]);
	  fprintf(stdout,"### Origin=(%g, %g, %g)\n",file[idxfile]->O[X],file[idxfile]->O[Y],file[idxfile]->O[Z]);
	  fprintf(stdout,"###        %g %g %g\n",file[idxfile]->L[X][X],file[idxfile]->L[X][Y],file[idxfile]->L[X][Z]);
	  fprintf(stdout,"### L=     %g %g %g\n",file[idxfile]->L[Y][X],file[idxfile]->L[Y][Y],file[idxfile]->L[Y][Z]);
	  fprintf(stdout,"###        %g %g %g\n",file[idxfile]->L[Z][X],file[idxfile]->L[Z][Y],file[idxfile]->L[Z][Z]);
	  fprintf(stdout,"### Limits atoms [%g,%g,%g] --> [%g,%g,%g]\n",
		  file[idxfile]->atom_lim[X][MINI],
		  file[idxfile]->atom_lim[Y][MINI],
		  file[idxfile]->atom_lim[Z][MINI],
		  file[idxfile]->atom_lim[X][MAXI],
		  file[idxfile]->atom_lim[Y][MAXI],
		  file[idxfile]->atom_lim[Z][MAXI]);
	  fprintf(stdout,"###              [%d,%d,%d] --> [%d,%d,%d]\n",
		  file[idxfile]->idx_atom_lim[X][MINI],
		  file[idxfile]->idx_atom_lim[Y][MINI],
		  file[idxfile]->idx_atom_lim[Z][MINI],
		  file[idxfile]->idx_atom_lim[X][MAXI],
		  file[idxfile]->idx_atom_lim[Y][MAXI],
		  file[idxfile]->idx_atom_lim[Z][MAXI]);
	  int i;
	  for(i=0;i<sline->ntokens;i++){
	    if(strcmp(sline->words[i],"point")==0) 	  {
	      int ix,iy,iz;
	      ix=atoi(sline->words[i+1]);
	      iy=atoi(sline->words[i+2]);
	      iz=atoi(sline->words[i+3]);
	      fprintf(stdout,"### DATA point (%d,%d,%d): x=%g y=%g z=%g data=%g\n",ix,iy,iz,
		      file[idxfile]->data[X][ix][iy][iz],file[idxfile]->data[Y][ix][iy][iz],file[idxfile]->data[Z][ix][iy][iz],file[idxfile]->data[DATA][ix][iy][iz]);
	    }
	  }
	}
	/* ############## ########################################## */
	if(strcmp(sline->words[0],"new_cubefile")==0)  {
	/* ################################################################################# */
	  int idxfileini=0;
	  if(strcmp(sline->words[1],"init_file")==0)  {
	    idxfileini=atoi(sline->words[2]);
	  }
	  file=realloc(file,(nfile+1)*sizeof(CubeFile));
	  file[nfile]=NewCubeFile(file[idxfileini]);
	  fprintf(stdout,"### Cubefile with index %d created\n",nfile);
	  nfile++;
	  fprintf(stdout,"### %d file in CubeTools\n",nfile);
	}
	/* ############## ########################################## */
	if(strcmp(sline->words[0],"xmap")==0)  {
	/* ################################################################################# */
	  int idxfile=0;
	  int ix0=0;
	  int i;
	  for(i=0;i<sline->ntokens;i++)	    {
	    if(strcmp(sline->words[i],"file")==0) 	  idxfile=atoi(sline->words[i+1]);
	    if(strcmp(sline->words[i],"ix0")==0) 	  ix0=atoi(sline->words[i+1]);
	  }
	  fprintf(stdout,"### x= %g\n",file[idxfile]->data[X][ix0][0][0]);
	  /* compute ymap */
	  int ix,iy,iz;
	  ix=ix0;
	  for (iz=0;iz<file[idxfile]->N[Z];iz++) {
	    for (iy=0;iy<file[idxfile]->N[Y];iy++) {
	      file[idxfile]->xmap[iy][iz]=file[idxfile]->data[DATA][ix][iy][iz];
	    }
	  }

	  double min,max;
	  min=file[idxfile]->xmap[0][0]; max=file[idxfile]->xmap[0][0];
	  for (iz=0;iz<file[idxfile]->N[Z];iz++) {
	    for (iy=0;iy<file[idxfile]->N[X];iy++) {
	      if(min>file[idxfile]->xmap[iy][iz]) min=file[idxfile]->xmap[iy][iz];
	      if(max<file[idxfile]->xmap[iy][iz]) max=file[idxfile]->xmap[iy][iz];
	    }
	  }
	  fprintf(stdout,"### min= %g max= %g\n",min,max);
	  
	  int nlevel=50;
	  for(i=0;i<sline->ntokens;i++){
	    if(strcmp(sline->words[i],"nlevel")==0) 	  nlevel=atof(sline->words[i+1]);
	    if(strcmp(sline->words[i],"min")==0) 	  min=atof(sline->words[i+1]);
	    if(strcmp(sline->words[i],"max")==0) 	  max=atof(sline->words[i+1]);
	  }
	  double d;
	  d=(max-min)/(nlevel+1);
	  PLFLT *clevel; clevel=malloc(nlevel*sizeof(PLFLT));
	  for(i=0;i<nlevel;i++) 	clevel[i]=min+(i+1)*d;

	
	  /* ************************ */
	  /* initialisation of PLPLOT */
	  plsdev( "pdfcairo" );
	  plsfnam( file[idxfile]->pdffile ); 
	  plspal1( "cmap1_gray.pal", 1 );	  //plspal0( "cmap0_black_on_white.pal" );
	  plscmap0n( 3 );
	  plinit();
	  pl_setcontlabelformat( 4, 3 );
	  pl_setcontlabelparam( 0.006, 0.3, 0.1, 1 );
	  PLFLT         fill_width = 2., cont_width = 0.;
	  PLINT         cont_color = 0;
	  /* set up the limits of the map */
	  PLFLT plenvxmin,plenvxmax,plenvymin,plenvymax;
	  plenvxmin=0.0;
	  plenvxmax=14.0;
	  plenvymin=4.0;
	  plenvymax=18.0;
	  for(i=0;i<sline->ntokens;i++){
	    if(strcmp(sline->words[i],"plenvxmin")==0) 	  plenvxmin=atof(sline->words[i+1]);
	    if(strcmp(sline->words[i],"plenvxmax")==0) 	  plenvxmax=atof(sline->words[i+1]);
	    if(strcmp(sline->words[i],"plenvymin")==0) 	  plenvymin=atof(sline->words[i+1]);
	    if(strcmp(sline->words[i],"plenvymax")==0) 	  plenvymax=atof(sline->words[i+1]);
	  }
	  plenv( plenvxmin,plenvxmax, plenvymin,plenvymax,1, 0 );
	  /* transformation matrix for the coordinate */
	  tr[0]=file[idxfile]->L[Y][Y];     tr[1]= file[idxfile]->L[Z][Y];    tr[2]=file[idxfile]->O[Y];
	  tr[3]=file[idxfile]->L[Y][Z];     tr[4]= file[idxfile]->L[Z][Z];    tr[5]=file[idxfile]->O[Z];

	  plcol0( 2 );
	  plshades( (const PLFLT * const *) file[idxfile]->xmap,
		    file[idxfile]->N[Y],file[idxfile]->N[Z], NULL,
		    1, file[idxfile]->N[Y], 1, file[idxfile]->N[Z],
		    clevel, nlevel,
		    fill_width,cont_color,cont_width,plfill,1,mypltr, NULL );
	  pllab( "x (ang.)", "y (ang.)", "Constant height LUMO+1 E" );

	  
	  int distribution=False;
	  for(i=0;i<sline->ntokens;i++){
	    if(strcmp(sline->words[i],"distrib")==0) 	  distribution=True;
	  }
	  if(distribution==True){
	    plcol0( 2 );
	    double histmin,histmax,histymin=0.0,histymax=1.0;
	    
	    int histlevel;
	    histmin=min;
	    histmax=max;
	    histlevel=100;
	    for(i=0;i<sline->ntokens;i++){
	      if(strcmp(sline->words[i],"histymax")==0) 	  histymax=atof(sline->words[i+1]);
	      if(strcmp(sline->words[i],"histymin")==0) 	  histymin=atof(sline->words[i+1]);
	      if(strcmp(sline->words[i],"histmax")==0) 	  histmax=atof(sline->words[i+1]);
	      if(strcmp(sline->words[i],"histmin")==0) 	  histmin=atof(sline->words[i+1]);
	      if(strcmp(sline->words[i],"histlevel")==0) 	  histlevel=atoi(sline->words[i+1]);
	    }
	    distrib(file[idxfile]->xmap,file[idxfile]->N[Y],file[idxfile]->N[Z],histmin,histmax,histlevel,histymin,histymax);
	  }
	  plend();
	  free(clevel);
	}
	/* ############## #################################################################### */
	if(strcmp(sline->words[0],"ymap")==0)  {
	  /*                                                                                   */
	  /* ymap is a 2D x-z map of LDOS                                                      */
	  /*                                                                                   */
	  /* ################################################################################# */
	  int idxfile=0;
	  int idxref=0;
	  int iy0=0;
	  int i;
	  for(i=0;i<sline->ntokens;i++)	    {
	    if(strcmp(sline->words[i],"file")==0) 	  idxfile=atoi(sline->words[i+1]);
	    if(strcmp(sline->words[i],"ref")==0) 	  idxref=atoi(sline->words[i+1]);
	    if(strcmp(sline->words[i],"iy0")==0) 	  iy0=atoi(sline->words[i+1]);
	  }
	  fprintf(stdout,"### y[%d]=%g\n",iy0,file[idxfile]->data[Y][0][iy0][0]);
	  /* compute ymap */
	  int ix,iy,iz;
	  iy=iy0;
	  for (iz=0;iz<file[idxfile]->N[Z];iz++) {
	    for (ix=0;ix<file[idxfile]->N[X];ix++) {
	      file[idxfile]->ymap[ix][iz]=file[idxfile]->data[DATA][ix][iy][iz];
	    }
	  }

	  double min,max;
	  min=file[idxfile]->ymap[0][0]; max=file[idxfile]->ymap[0][0];
	  for (iz=0;iz<file[idxfile]->N[Z];iz++) {
	    for (ix=0;ix<file[idxfile]->N[X];ix++) {
	      if(min>file[idxfile]->ymap[ix][iz]) min=file[idxfile]->ymap[ix][iz];
	      if(max<file[idxfile]->ymap[ix][iz]) max=file[idxfile]->ymap[ix][iz];
	    }
	  }
	  fprintf(stdout,"### min= %g max= %g\n",min,max);
	  
	  int nlevel=50;
	  for(i=0;i<sline->ntokens;i++){
	    if(strcmp(sline->words[i],"nlevel")==0) 	  nlevel=atoi(sline->words[i+1]);
	    if(strcmp(sline->words[i],"min")==0) 	  min=atof(sline->words[i+1]);
	    if(strcmp(sline->words[i],"max")==0) 	  max=atof(sline->words[i+1]);
	  }
	  double d;
	  d=(max-min)/(nlevel+1);
	  PLFLT *clevel; clevel=malloc(nlevel*sizeof(PLFLT));
	  for(i=0;i<nlevel;i++) 	clevel[i]=min+(i+1)*d;

	
	  /* ************************ */
	  /* initialisation of PLPLOT */
	  plsdev( "pdfcairo" );
	  plsfnam( file[idxfile]->pdffile ); 
	  plspal1( "cmap1_gray.pal", 1 );	  //plspal0( "cmap0_black_on_white.pal" );
	  //plscmap0n( 3 );
	  plinit();
	  pl_setcontlabelformat( 4, 3 );
	  pl_setcontlabelparam( 0.006, 0.3, 0.1, 1 );
	  PLFLT         fill_width = 2., cont_width = 0.;
	  PLINT         cont_color = 0;
	  /* set up the limits of the map */
	  PLFLT plenvxmin,plenvxmax,plenvymin,plenvymax;
	  plenvxmin=file[idxfile]->O[X];
	  plenvxmax=file[idxfile]->O[X]+file[idxfile]->DIM[X];
	  plenvymin=file[idxfile]->O[Z];
	  plenvymax=file[idxfile]->O[Z]+file[idxfile]->DIM[Z];
	  for(i=0;i<sline->ntokens;i++){
	    if(strcmp(sline->words[i],"plenvxmin")==0) 	  plenvxmin=atof(sline->words[i+1]);
	    if(strcmp(sline->words[i],"plenvxmax")==0) 	  plenvxmax=atof(sline->words[i+1]);
	    if(strcmp(sline->words[i],"plenvymin")==0) 	  plenvymin=atof(sline->words[i+1]);
	    if(strcmp(sline->words[i],"plenvymax")==0) 	  plenvymax=atof(sline->words[i+1]);
	  }
	  plenv( plenvxmin,plenvxmax, plenvymin,plenvymax,1, 0 );
	  /* transformation matrix for the coordinate */
	  tr[0]=file[idxfile]->L[X][X];     tr[1]= file[idxfile]->L[Z][X];    tr[2]=file[idxfile]->O[X];
	  tr[3]=file[idxfile]->L[X][Z];     tr[4]= file[idxfile]->L[Z][Z];    tr[5]=file[idxfile]->O[Z];

	  plcol0( 2 );
	  plshades( (const PLFLT * const *) file[idxfile]->ymap,
		    file[idxfile]->N[X],file[idxfile]->N[Z], NULL,
		    1, file[idxfile]->N[X], 1, file[idxfile]->N[Z],
		    clevel, nlevel,
		    fill_width,cont_color,cont_width,plfill,1,mypltr, NULL );

	  /* -------------------- */
	  /* for plotting contour */
	  /* -------------------- */
	  int contour=False;
	  for(i=0;i<sline->ntokens;i++){
	    if(strcmp(sline->words[i],"contour")==0) contour=True;
	  }
	  if(contour==True){
	    int nlevelcont=5;
	    for(i=0;i<sline->ntokens;i++){
	      if(strcmp(sline->words[i],"nlevelcont")==0) nlevelcont=atoi(sline->words[i+1]);
	    }
	    PLFLT *clevelcont; clevelcont=malloc(nlevelcont*sizeof(PLFLT));

	    for(i=0;i<sline->ntokens;i++){
	      if(strcmp(sline->words[i],"automatic")==0) {
		double contmin,contmax; contmin=min;contmax=max;
		for(i=0;i<sline->ntokens;i++){
		  if(strcmp(sline->words[i],"contmin")==0) 	  contmin=atof(sline->words[i+1]);
		  if(strcmp(sline->words[i],"contmax")==0) 	  contmax=atof(sline->words[i+1]);
		}
		d=(contmax-contmin)/(nlevelcont+1);
		for(i=0;i<nlevelcont;i++) 	clevelcont[i]=contmin+(i+1)*d;
	      }
	      if(strcmp(sline->words[i],"man")==0) {
		int ii;
		for(ii=0;ii<nlevelcont;ii++){
		  clevelcont[ii]=atof(sline->words[i+ii+1]);
		}
	      }
	    }
	    pl_setcontlabelparam( 0.006, 0.3, 0.1, 0);
	    //plschr(0.0,0.0);
	    plcol0( 2 );
	    plcont( (const PLFLT * const *) file[idxfile]->ymap,
		    file[idxfile]->N[X],file[idxfile]->N[Z],
		    1, file[idxfile]->N[X], 1, file[idxfile]->N[Z],
		    clevelcont, nlevelcont,
		    mypltr, NULL );
	    pllab( "x (ang.)", "z (ang.)", "ymap" );
	    free(clevelcont);
	  }

	  /* -------- */
	  /* Atom  */
	  /* -------- */
	  double *posx,*posz;
	  int ii;
	  double datom=.5;
	  for(ii=0;ii<sline->ntokens;ii++){
	    if(strcmp(sline->words[ii],"atom")==0) {
	      datom=atof(sline->words[ii+1]);
	      fprintf(stdout,"### datom= %g\n",datom);
	    }
	  }
	  posx=malloc(file[idxfile]->natom*sizeof(double));
	  posz=malloc(file[idxfile]->natom*sizeof(double));
	  plcol0(4);
	  
	  for(i=0;i<file[idxfile]->natom;i++){
	    if(fabs(file[idxfile]->q[i][Y]-file[idxfile]->data[Y][0][iy0][0])<datom){
	      //posx[i]=file[idxfile]->q[i][X]-file[idxfile]->MC[X]+0.5*file[idxfile]->DIM[X]/2;
	      //posz[i]=file[idxfile]->q[i][Z]-file[idxfile]->MC[Z]+0.5*file[idxfile]->DIM[Z]/2;
	      posx[i]=file[idxfile]->q[i][X]+file[idxfile]->data[Y][0][iy0][0]/sqrt(3.0);
	      posz[i]=file[idxfile]->q[i][Z];
	      fprintf(stdout,"%g %g %g\n",posx[i],file[idxfile]->q[i][Y],posz[i]); 
	    }
	    /* posx[i]=tr[0]*file[idxfile]->q[i][X]+tr[1]*file[idxfile]->q[i][Z]+tr[3]+plenvxmin; */
	    /* posz[i]=tr[3]*file[idxfile]->q[i][X]+tr[4]*file[idxfile]->q[i][Z]+tr[5]+plenvymin; */

	    /* plptex(posx[i],posz[i],1.0,0,0.5,text); */
	  }
	  /* plschr(30.0,1.0); */
	  /* plwidth(10.0); */
	  plssym(0.0,5.0);
	  plcol0(5);
	  plpoin(file[idxfile]->natom,posx,posz,17);
	  pllab( "x (ang.)", "z (ang.)", "ymap" );


	  /* --------------------------------------------------------- 
	     For plotting a line on the ymap and then plot a LDOS=f(x) 
	     Note that the trajectory plotted on the map can be different  
	        from the trajectory computed from idxfile
	     --------------------------------------------------------- */
	  for(i=0;i<sline->ntokens;i++) {
	    if(strcmp(sline->words[i],"line")==0) {
	      int iz0; iz0=atoi(sline->words[i+2]);
	      int nline;nline=file[idxfile]->N[X];
	      double *linex,*liney;
	      linex=malloc(nline*sizeof(double)); liney=malloc(nline*sizeof(double));
	      iz=0;ix=0;iy=0;
	      for(ix=0;ix<file[idxfile]->N[X];ix++) {
		linex[ix]=file[idxfile]->data[X][ix][0][0];
		if(strcmp(sline->words[i+1],"LDOS")==0)
		  liney[ix]=file[idxref]->LDOSmap[ix][iy0];
		if(strcmp(sline->words[i+1],"z")==0)
		  liney[ix]=file[idxfile]->data[Z][0][0][iz0];
	      }
	      plcol0(4);
	      plwidth(3.0);
	      plline(nline,linex,liney);
	      plwidth(1.0);
	      free(linex);free(liney);

	      /* linescan plot */
	      {
		double *datax;datax=malloc(file[idxfile]->N[X]*sizeof(double));
		double *datay;datay=malloc(file[idxfile]->N[X]*sizeof(double));
		if(strcmp(sline->words[i+1],"z")==0){
		  /* for a constant height */
		  for(ix=0;ix<file[idxfile]->N[X];ix++){
		    datax[ix]=file[idxfile]->data[X][ix][iy0][iz0];
		    datay[ix]=file[idxfile]->data[DATA][ix][iy0][iz0];
		  }
		}
		if(strcmp(sline->words[i+1],"LDOS")==0){
		  /* for LDOS */
		  /* file[idxref]->LDOSmap[ix][iy0] give the z coordinate  */
		  double ztarget,z0,z1,I1,I0;
		  for(ix=0;ix<file[idxfile]->N[X];ix++){
		    ztarget=file[idxref]->LDOSmap[ix][iy0];
		    iz=file[idxfile]->N[Z]-1;
		    while(file[idxfile]->data[Z][ix][iy0][iz]>ztarget && iz>=0) iz--;
		    z0=file[idxfile]->data[Z][ix][iy0][iz];
		    z1=file[idxfile]->data[Z][ix][iy0][iz+1];
		    I0=file[idxfile]->data[DATA][ix][iy0][iz];
		    I1=file[idxfile]->data[DATA][ix][iy0][iz+1];
		    datax[ix]=file[idxfile]->data[X][ix][iy0][iz];
		    datay[ix]=((I1-I0)*ztarget+z1*I0-z0*I1)/(z1-z0);
		  }
		}
		
		double miny,maxy,minx,maxx;
		miny=datay[0];maxy=datay[0];
		for(ix=0;ix<file[idxfile]->N[X];ix++){
		  if(miny>datay[ix]) miny=datay[ix];
		  if(maxy<datay[ix]) maxy=datay[ix]; 
		}
		//pladv( 0 );

		minx=datax[0];
		maxx=datax[file[idxfile]->N[X]-1];
		int ii;
		for(ii=0;ii<sline->ntokens;ii++) {
		  if(strcmp(sline->words[ii],"minx")==0) minx=atof(sline->words[ii+1]);
		  if(strcmp(sline->words[ii],"maxx")==0) maxx=atof(sline->words[ii+1]);
		  if(strcmp(sline->words[ii],"miny")==0) miny=atof(sline->words[ii+1]);
		  if(strcmp(sline->words[ii],"maxy")==0) maxy=atof(sline->words[ii+1]);
		}
		plenv(minx,maxx,miny,maxy,0, 0 );
		
		plline(file[idxfile]->N[X],datax,datay);
		pllab( "x (ang.)", "LDOS", "" );		

		{
		  char *name;name=malloc(NCHAR*sizeof(char));
		  sprintf(name,"%s.LDOS.dat",file[idxfile]->pdffile);
		  FILE *out; out=fopen(name,"w+");
		  fprintf(stdout,"# %s %s\n","x(ang.)","LDOS");
		  for(ix=0;ix<file[idxfile]->N[X];ix++){
		    fprintf(out,"%g %g\n",datax[ix],datay[ix]);
		  }
		  fclose(out);
		  free(name);
		}
		/* ---------------------------------------- */
		//pladv( 0 );
		for(ix=0;ix<file[idxfile]->N[X];ix++){
		  datay[ix]=file[idxfile]->LDOSmap[ix][iy0];
		}
		miny=datay[0];maxy=datay[0];
		for(ix=0;ix<file[idxfile]->N[X];ix++){
		  if(miny>datay[ix]) miny=datay[ix];
		  if(maxy<datay[ix]) maxy=datay[ix]; 
		}
		minx=datax[0];
		maxx=datax[file[idxfile]->N[X]-1];
		for(ii=0;ii<sline->ntokens;ii++) {
		  if(strcmp(sline->words[ii],"zminx")==0) minx=atof(sline->words[ii+1]);
		  if(strcmp(sline->words[ii],"zmaxx")==0) maxx=atof(sline->words[ii+1]);
		  if(strcmp(sline->words[ii],"zminy")==0) miny=atof(sline->words[ii+1]);
		  if(strcmp(sline->words[ii],"zmaxy")==0) maxy=atof(sline->words[ii+1]);
		}
		plenv(minx,maxx,miny,maxy,0, 0 );




		plline(file[idxfile]->N[X],datax,datay);
		pllab( "x (ang.)", "z (ang.)", "" );
		{
		  char *name;name=malloc(NCHAR*sizeof(char));
		  sprintf(name,"%s.z.dat",file[idxfile]->pdffile);
		  FILE *out; out=fopen(name,"w+");
		  fprintf(stdout,"# %s %s\n","x(ang.)","z(ang.)");
		  for(ix=0;ix<file[idxfile]->N[X];ix++){
		    fprintf(out,"%g %g\n",datax[ix],datay[ix]);
		  }
		  fclose(out);
		  free(name);
		}


		
		free(datax);free(datay);
		
	      }
	    }
	  }
	  free(posx);free(posz);
	  
	  /* ------------------------------------------- */
	  /* for computing and plotting the distribution */
	  /* ------------------------------------------- */
	  int distribution=False;
	  for(i=0;i<sline->ntokens;i++){
	    if(strcmp(sline->words[i],"distrib")==0) 	  distribution=True;
	  }
	  if(distribution==True){
	    double histmin,histmax,histymin=0.0,histymax=1.0;
	    
	    int histlevel;
	    histmin=min;
	    histmax=max;
	    histlevel=100;
	    for(i=0;i<sline->ntokens;i++){
	      if(strcmp(sline->words[i],"histymax")==0) 	  histymax=atof(sline->words[i+1]);
	      if(strcmp(sline->words[i],"histymin")==0) 	  histymin=atof(sline->words[i+1]);
	      if(strcmp(sline->words[i],"histmax")==0) 	  histmax=atof(sline->words[i+1]);
	      if(strcmp(sline->words[i],"histmin")==0) 	  histmin=atof(sline->words[i+1]);
	      if(strcmp(sline->words[i],"histlevel")==0) 	  histlevel=atoi(sline->words[i+1]);
	    }
	    distrib(file[idxfile]->ymap,file[idxfile]->N[X],file[idxfile]->N[Z],histmin,histmax,histlevel,histymin,histymax);

	  }
	  plend();
	  free(clevel);
	}
	/* ############## ########################################## */
	if(strcmp(sline->words[0],"constant_height")==0)  {
	/* ################################################################################# */
	  int idxfile=0;
	  int iz0=0;
	  int i;
	  for(i=0;i<sline->ntokens;i++)	    {
	    if(strcmp(sline->words[i],"file")==0) 	  idxfile=atoi(sline->words[i+1]);
	    if(strcmp(sline->words[i],"iz0")==0) 	  iz0=atoi(sline->words[i+1]);
	  }

	  fprintf(stdout,"NX=%d\nNY=%d\nNZ=%d\n",file[idxfile]->N[X],file[idxfile]->N[Y],file[idxfile]->N[Z]);

	  static int nx      ;        // Default number of data points in x
	  static int ny      ;        // Default number of data points in y
	  nx=file[idxfile]->N[X];
	  ny=file[idxfile]->N[Y];
	  PLFLT        **z;
	  plsdev( "pdfcairo" );
	  plsfnam( file[idxfile]->pdffile ); 
	  //plspal0( "cmap0_black_on_white.pal" );
	  plspal1( "cmap1_gray.pal", 1 );
	  plscmap0n( 3 );
	  plinit();
	  plAlloc2dGrid( &z, nx, ny );
	  
	  int ix,iy,iz;
	  double zmin,zmax,dz;
	  PLFLT plenvxmin,plenvxmax,plenvymin,plenvymax;
	  ix=0;
	  iy=0;
	  iz=iz0;
	  fprintf(stdout,"z=%g zatom=%g dz=%g angstroems\n",
		  file[idxfile]->data[Z][ix][iy][iz0],
		  file[idxfile]->atom_lim[MAXI][Z],
		  file[idxfile]->data[Z][ix][iy][iz0]- file[idxfile]->atom_lim[MAXI][Z]);
	  zmin=file[idxfile]->data[DATA][ix][iy][iz];
	  zmax=file[idxfile]->data[DATA][ix][iy][iz];
	  
	  tr[0]=file[idxfile]->L[X][X];       tr[1]=file[idxfile]->L[Y][X];     tr[2]=file[idxfile]->O[X];
	  tr[3]=file[idxfile]->L[X][Y];       tr[4]=-file[idxfile]->L[Y][Y];    tr[5]=file[idxfile]->O[Y]+file[idxfile]->DIM[Y];
	  for (ix=0;ix<file[idxfile]->N[X];ix++) {
	    for (iy=0;iy<file[idxfile]->N[Y];iy++) {
	      z[ix][iy]=file[idxfile]->data[DATA][ix][iy][iz];
	      if(z[ix][iy]>zmax) zmax=z[ix][iy];
	      if(z[ix][iy]<zmin) zmin=z[ix][iy];
	    }
	  }
	  //define the range and scale of the graph, and draw labels, axes, etc.
	  plenvxmin=file[idxfile]->O[X];
	  plenvxmax=file[idxfile]->DIM[X]*.5;
	  plenvymin=file[idxfile]->O[Y]+file[idxfile]->DIM[Y]*.3;
	  plenvymax=file[idxfile]->DIM[Y]*.7;
	  plenvxmin=0.0;
	  plenvxmax=14.0;
	  plenvymin=4.0;
	  plenvymax=18.0;
	  
	  int nlevel=50;
	  dz=(zmax-zmin)/(nlevel+1);
	  fprintf(stdout,"# zmin= %g zmax= %g dz= %g nlevel= %d\n",zmin,zmax,dz,nlevel);
	  PLFLT *clevel; clevel=malloc(nlevel*sizeof(PLFLT));
	  /* clevel[0]=1.0e-5; */
	  /* clevel[1]=1.0e-4; */
	  /* clevel[2]=1.0e-3; */
	  /* clevel[3]=1.0e-2; */

	  for(i=0;i<nlevel;i++) {
	    clevel[i]=zmin+(i+1)*dz;
	    //fprintf(stdout,"# clevel[%d]= %g\n",i,clevel[i]);
	  }
	  
	  
	  pl_setcontlabelformat( 4, 3 );
	  pl_setcontlabelparam( 0.006, 0.3, 0.1, 1 );
	  //double alpha;alpha=0.3;
	  //      plenv( 0.0, 1.0*nx, 0.0, 1.0*ny, 2, 0 );
	  //plenv( alpha*nx, (1.0-alpha)*nx, alpha*ny, (1.0-alpha)*ny, 2, 0 );
	  plenv( plenvxmin,plenvxmax, plenvymin,plenvymax,1, 0 );
	  plcol0( 2 );
	  
	  
	  PLFLT         fill_width = 2., cont_width = 0.;
	  //      PLFLT         colorbar_width, colorbar_height;
	  PLINT         cont_color = 0;
	  plshades( (const PLFLT * const *) z, nx, ny, NULL,1, nx, 1, ny, clevel, nlevel, fill_width,cont_color,cont_width,plfill,1,mypltr, NULL );
	  pllab( "x (ang.)", "y (ang.)", "Constant height LUMO+1 E" );


	  

	  
	  plFree2dGrid( z, nx,ny);


	  int line=False;
	  for(i=0;i<sline->ntokens;i++){
	    if(strcmp(sline->words[i],"line")==0) 	  {
	      line=True;
	      double linex[2],liney[2];
	      iz=0;ix=0;iy=0;
	      linex[0]=file[idxfile]->data[X][atoi(sline->words[i+1])][0][0]+file[idxfile]->O[X];
	      liney[0]=file[idxfile]->data[Y][0][atoi(sline->words[i+2])][0]+file[idxfile]->O[Y];
	      linex[1]=file[idxfile]->data[X][atoi(sline->words[i+3])][0][0]+file[idxfile]->O[X];
	      liney[1]=file[idxfile]->data[Y][0][atoi(sline->words[i+4])][0]+file[idxfile]->O[Y];
	      plcol0( 2 );
	      plline(2,linex,liney);

	    }

	  }

	  
	  plend();
  
	  
	}
	/* ################################################################################# */
	if(strcmp(sline->words[0],"save")==0)  {
	  /* ################################################################################# */
	  /* 
	   */
	  int idxfiletosave=-1,i;
	  char name[NCHAR]; strcpy(name,"default.cube");
	  for(i=0;i<sline->ntokens;i++)	    {
	    if(strcmp(sline->words[i],"idx")==0) 	  {
	      idxfiletosave=atoi(sline->words[i+1]);
	    }
	    if(strcmp(sline->words[i],"name")==0) 	  {
	      strcpy(name,sline->words[i+1]);
	    }
	  }
	  if(idxfiletosave>=0){
	    WriteCubeFile(file[idxfiletosave],name);
	  }
	}
	/* ############## ########################################## */
	if(strcmp(sline->words[0],"diff")==0)  {
	/* ################################################################################# */
	  int i,idx1=-1,idx2=-1,idxfiletostock=-1;
	  double normalize=1.0;
	  for(i=0;i<sline->ntokens;i++)	    {
	    if(strcmp(sline->words[i],"to")==0) 	  {
	      idxfiletostock=atoi(sline->words[i+1]);
	    }
	    if(strcmp(sline->words[i],"idx1")==0) 	  {
	      idx1=atoi(sline->words[i+1]);
	    }
	    if(strcmp(sline->words[i],"idx2")==0) 	  {
	      idx2=atoi(sline->words[i+1]);
	    }
	    if(strcmp(sline->words[i],"normalize")==0) 	  {
	      normalize=0.7071067;
	    }
	  }
	  if(idxfiletostock>=0 && idx1>=0 && idx2>=0){
	    int ix,iy,iz;
	    for (ix=0;ix<file[idxfiletostock]->N[X];ix++) {
	      for (iy=0;iy<file[idxfiletostock]->N[Y];iy++) {
		for (iz=0;iz<file[idxfiletostock]->N[Z];iz++) {
		  file[idxfiletostock]->data[DATA][ix][iy][iz]=normalize*(file[idx1]->data[DATA][ix][iy][iz]-file[idx2]->data[DATA][ix][iy][iz]);
		}
	      }
	    }
	  } else {
	    printf("### !!! ERROR !!! in diff %d %d %d\n",idxfiletostock,idx1,idx2);
	  }
	  
	}
	/* ############## ########################################## */
	if(strcmp(sline->words[0],"transform")==0)  {
	/* ################################################################################# */
	  /*
	    transform idx 0 by sqrt to 2
	   */
	  char op[8];
	  int i,idxfiletostock=-1,idxfiletotransform=-1;
	  for(i=0;i<sline->ntokens;i++)	    {
	    if(strcmp(sline->words[i],"to")==0) 	  {
	      idxfiletostock=atoi(sline->words[i+1]);
	    }
	    if(strcmp(sline->words[i],"idx")==0) 	  {
	      idxfiletotransform=atoi(sline->words[i+1]);
	    }
	    if(strcmp(sline->words[i],"by")==0) 	  {
	      strcpy(op,sline->words[i+1]);
	    }
	  }
	  printf("%s\n",op);
	  if(idxfiletostock>=0 && idxfiletotransform>=0){
	    int ix,iy,iz;
	    double sign;
	    if(strcmp(op,"sqrt")==0){
	      for (ix=0;ix<file[idxfiletostock]->N[X];ix++) {
		for (iy=0;iy<file[idxfiletostock]->N[Y];iy++) {
		  for (iz=0;iz<file[idxfiletostock]->N[Z];iz++) {
		    if(file[idxfiletotransform]->data[DATA][ix][iy][iz]<0.0) {
		      sign=-1.0;
		    }else{
		      sign=1.0;
		    }
		    file[idxfiletostock]->data[DATA][ix][iy][iz]=sign*pow(fabs(file[idxfiletotransform]->data[DATA][ix][iy][iz]),.5);
		  }
		}
	      }
	    } else if(strcmp(op,"sqr")==0){
	      for (ix=0;ix<file[idxfiletostock]->N[X];ix++) {
		for (iy=0;iy<file[idxfiletostock]->N[Y];iy++) {
		  for (iz=0;iz<file[idxfiletostock]->N[Z];iz++) {
		    file[idxfiletostock]->data[DATA][ix][iy][iz]=pow(fabs(file[idxfiletotransform]->data[DATA][ix][iy][iz]),2);
		  }
		}
	      }
	    }	    else {
	      printf("### !!! ERROR !!! %s unknown operator in transform\n",op);
	    }
	  } else {
	    printf("### !!! ERROR !!! in transform %d %d\n",idxfiletotransform,idxfiletostock);
	  }
	}
	/* ############## ########################################## */
	if(strcmp(sline->words[0],"add")==0)  {
	/* ################################################################################# */
	  int i,idxfiletostock;
	  int nfiletoadd,*listfiletoadd; nfiletoadd=sline->ntokens-1;
	  int normalize=FALSE;
	  double norm=1.0;
	  for(i=0;i<sline->ntokens;i++)	    {
	    if(strcmp(sline->words[i],"to")==0) 	  {
	      idxfiletostock=atoi(sline->words[i+1]);
	      nfiletoadd-=2;
	    }
	    if(strcmp(sline->words[i],"normalize")==0) 	  {
	      nfiletoadd--;
	      normalize=TRUE;
	    }
	  }
	  int j;
	  for(i=0;i<sline->ntokens;i++)	    {
	    if(strcmp(sline->words[i],"list")==0) 	  {
	      nfiletoadd--;
	      listfiletoadd=malloc(nfiletoadd*sizeof(int));
	      for(j=0;j<nfiletoadd;j++) listfiletoadd[j]=atoi(sline->words[i+1+j]);
	    }
	  }
	  if(normalize==TRUE) norm=pow(nfiletoadd,-.5);
	  fprintf(stdout,"### %d file to add\n### list:",nfiletoadd);
	  for(j=0;j<nfiletoadd;j++) fprintf(stdout,"### %d ",listfiletoadd[j]);
	  fprintf(stdout,"\n");
	  int ix,iy,iz;
	  for (ix=0;ix<file[idxfiletostock]->N[X];ix++) {
	    for (iy=0;iy<file[idxfiletostock]->N[Y];iy++) {
	      for (iz=0;iz<file[idxfiletostock]->N[Z];iz++) {
		for(i=0;i<nfiletoadd;i++){
		  j=listfiletoadd[i];
		  file[idxfiletostock]->data[DATA][ix][iy][iz]+=norm*file[j]->data[DATA][ix][iy][iz]; 
		}
	      }
	    }
	  }



	  if(listfiletoadd) free(listfiletoadd);
	}
	/* ############## ########################################## */
	if(strcmp(sline->words[0],"average")==0)  {
	/* ################################################################################# */
	  int i,idxfiletostock;
	  int nfiletoadd,*listfiletoadd; nfiletoadd=sline->ntokens-1;
	  for(i=0;i<sline->ntokens;i++)	    {
	    if(strcmp(sline->words[i],"to")==0) 	  {
	      idxfiletostock=atoi(sline->words[i+1]);
	      nfiletoadd-=2;
	    }
	  }
	  int j;
	  for(i=0;i<sline->ntokens;i++)	    {
	    if(strcmp(sline->words[i],"list")==0) 	  {
	      nfiletoadd--;
	      listfiletoadd=malloc(nfiletoadd*sizeof(int));
	      for(j=0;j<nfiletoadd;j++) listfiletoadd[j]=atoi(sline->words[i+1+j]);
	    }
	  }
	  fprintf(stdout,"### %d file to add\n### list:",nfiletoadd);
	  for(j=0;j<nfiletoadd;j++) fprintf(stdout,"### %d ",listfiletoadd[j]);
	  fprintf(stdout,"\n");
	  int ix,iy,iz;
	  for (ix=0;ix<file[idxfiletostock]->N[X];ix++) {
	    for (iy=0;iy<file[idxfiletostock]->N[Y];iy++) {
	      for (iz=0;iz<file[idxfiletostock]->N[Z];iz++) {
		for(i=0;i<nfiletoadd;i++){
		  j=listfiletoadd[i];
		  file[idxfiletostock]->data[DATA][ix][iy][iz]+=file[j]->data[DATA][ix][iy][iz]; 
		}
		file[idxfiletostock]->data[DATA][ix][iy][iz]/=nfiletoadd;
	      }
	    }
	  }



	  if(listfiletoadd) free(listfiletoadd);
	}
	/* ############## ########################################## */
	if(strcmp(sline->words[0],"3D")==0)  {
	  SDL_Init(SDL_INIT_VIDEO);// Démarrage de la SDL (ici : chargement du système vidéo)
	  SDL_WM_SetCaption("CubeTools 3D view !",NULL);
	  SDL_SetVideoMode(640, 640, 32, SDL_OPENGL);
	  
	  
	  SDL_GL_SetAttribute( SDL_GL_RED_SIZE, 5 );
	  SDL_GL_SetAttribute( SDL_GL_GREEN_SIZE, 5 );
	  SDL_GL_SetAttribute( SDL_GL_BLUE_SIZE, 5 );
	  SDL_GL_SetAttribute( SDL_GL_DEPTH_SIZE, 16 );
	  SDL_GL_SetAttribute( SDL_GL_DOUBLEBUFFER, 1 );
	  
	  double width = 640;
	  double height = 640;
	  
	  setup_opengl( width, height );
	  setup_light(file[0]);
	  draw_screen(file[0]);
	  SDL_Event event;
	  //for(;;){
	  int cont=True;
	  while(cont==True){
	    SDL_WaitEvent(&event);
	    switch(event.type) {
	    case SDL_QUIT:
	      cont=False;
	      break;
	      
	    }

	    //	draw_screen(file[0]);
	    
	  }

	  
	  SDL_Quit();

	/* ################################################################################# */
	}
	/* ############## ########################################## */
	if(strcmp(sline->words[0],"projection")==0)  {
	/* ################################################################################# */
	  int i,ix,iy,iz;
	  double x1,y1,x0,y0;
	  int idxfile=0,idxref=0;
	  for(i=0;i<sline->ntokens;i++)	    {
	    if(strcmp(sline->words[i],"file")==0) 	  idxfile=atoi(sline->words[i+1]);
	    if(strcmp(sline->words[i],"ref")==0) 	  idxref=atoi(sline->words[i+1]);
	  }

	  PLFLT histmin,histmax;
	  int histlevel;
	  histmin=0.0;
	  histmax=1.0;
	  histlevel=100;
	  for(i=0;i<sline->ntokens;i++){
	    if(strcmp(sline->words[i],"histmax")==0) 	  histmax=atof(sline->words[i+1]);
	    if(strcmp(sline->words[i],"histmin")==0) 	  histmin=atof(sline->words[i+1]);
	    //if(strcmp(sline->words[i],"mapmax")==0) 	  max=atof(sline->words[i+1]);
	    //if(strcmp(sline->words[i],"mapmin")==0) 	  min=atof(sline->words[i+1]);
	    if(strcmp(sline->words[i],"histlevel")==0) 	  histlevel=atoi(sline->words[i+1]);
	  }




	  double *data;
	  /* PLINT NPTS;NPTS=(PLINT) file[idxfile]->N[X]*file[idxfile]->N[Y]; */
	  /* data=malloc(NPTS*sizeof(PLFLT)); */
	  data=malloc(sizeof(double));
	  i=0;
	  /* LDOSmap contains the z coordinate  */
	  for (iy=0;iy<file[idxfile]->N[Y];iy++) {
	    for (ix=0;ix<file[idxfile]->N[X];ix++) {
	      file[idxfile]->zmap[ix][iy]=0.0;
	      if(file[idxref]->LDOSmap[ix][iy]>0.0){
		iz=file[idxfile]->N[Z]-1;
		while(file[idxfile]->data[Z][ix][iy][iz]>file[idxref]->LDOSmap[ix][iy] && iz>0) iz--;
		if(iz>0){
		  x0=file[idxfile]->data[Z][ix][iy][iz];
		  y0=file[idxfile]->data[DATA][ix][iy][iz];
		  x1=file[idxfile]->data[Z][ix][iy][iz+1];
		  y1=file[idxfile]->data[DATA][ix][iy][iz+1];
		  file[idxfile]->zmap[ix][iy]=(file[idxref]->LDOSmap[ix][iy]*(y1-y0)-y1*x0+y0*x1)/(x1-x0);
		  
		  data[i]=file[idxfile]->zmap[ix][iy];
		  //		  fprintf(stdout,"%g %g %g %g %g %g\n",file[idxref]->LDOSmap[ix][iy],x1,y1,x0,y0,data[i]);
		  i++;
		  data=realloc(data,(i+1)*sizeof(double));

		}
	      }
	    }
	  }
	  PLINT NPTS;NPTS=(PLINT)i;

	  plsdev( "pdfcairo" );
	  plsfnam( file[idxfile]->pdffile );
	  //plspal0( "cmap0_black_on_white.pal" );
	  plspal1( "cmap1_gray.pal", 1 );
	  plscmap0n( 3 );
	  plinit();
	  tr[0]=file[idxfile]->L[X][X];       tr[1]=file[idxfile]->L[Y][X];       tr[2]=file[idxfile]->O[X];
	  tr[3]=file[idxfile]->L[X][Y];       tr[4]=-file[idxfile]->L[Y][Y];       tr[5]=file[idxfile]->O[Y]+file[idxfile]->DIM[Y];
	  PLFLT plenvxmin,plenvxmax,plenvymin,plenvymax;
	  plenvxmin=0.0;
	  plenvxmax=14.0;
	  plenvymin=4.0;
	  plenvymax=18.0;
	  
	  double min,max;
	  min=data[0];max=data[0];
	  for(i=1;i<NPTS;i++) {
	    if(min>data[i]) min=data[i];
	    if(max<data[i]) max=data[i];
	  }
	  for(i=0;i<sline->ntokens;i++){
	    if(strcmp(sline->words[i],"mapmax")==0) 	  max=atof(sline->words[i+1]);
	    if(strcmp(sline->words[i],"mapmin")==0) 	  min=atof(sline->words[i+1]);
	  }




	  int nlevel=50;
	  double d;
	  d=(max-min)/(nlevel+1);
	  PLFLT *clevel; clevel=malloc(nlevel*sizeof(PLFLT));
	  for(i=0;i<nlevel;i++) 	clevel[i]=min+(i+1)*d;
	  pl_setcontlabelformat( 4, 3 );
	  pl_setcontlabelparam( 0.006, 0.3, 0.1, 1 );
	  plenv( plenvxmin,plenvxmax, plenvymin,plenvymax,1, 0 );
	  plcol0( 2 );
	  PLFLT         fill_width = 2., cont_width = 0.;
	  //      PLFLT         colorbar_width, colorbar_height;
	  PLINT         cont_color = 0;
	  static int nx      ;        // Default number of data points in x
	  static int ny      ;        // Default number of data points in y
	  nx=file[idxfile]->N[X];
	  ny=file[idxfile]->N[Y];

	  plshades( (const PLFLT * const *) file[idxfile]->zmap, nx, ny, NULL,1, nx, 1, ny, clevel, nlevel, fill_width,cont_color,cont_width,plfill,1,mypltr, NULL );
	  char *label;label=malloc(NCHAR*sizeof(char));
	  sprintf(label,"Constant LDOS (%g)",file[idxfile]->LDOS);
	  pllab( "x (ang.)", "y (ang.)", label );
	  free(label);




	  plcol0( 1 );
	  plhist( NPTS, data, histmin,histmax, histlevel, PL_HIST_IGNORE_OUTLIERS );
	  plcol0( 2 );
	  pllab( "#frValue", "#frFrequency",
		 "#frPLplot Example 5 - Probability function of Oscillator" );
	  free(data);
	  plend();
	  

	}
	/* ############## ########################################## */
	if(strcmp(sline->words[0],"shell")==0)  {
	  /* ################################################################################# */
	  char **string;string=malloc(2*sizeof(char*));
	  string[0]=malloc(NCHAR*sizeof(char));
	  string[1]=malloc(NCHAR*sizeof(char));
	  
	  int i,icur,iold,itmp;
	  sprintf(string[0],"%s ",sline->words[1]);
	  fprintf(stdout,"%s\n",string[0]);
	  icur=1;iold=0;
	  for(i=2;i<sline->ntokens;i++){
	    sprintf(string[icur],"%s %s ",string[iold],sline->words[i]);
	    itmp=iold;iold=icur;icur=itmp;
	    
	  }
	  system(string[iold]);
	  free(string[0]);free(string[1]);free(string);
	}
	/* ################################################################################### */
	if(strcmp(sline->words[0],"constant_LDOS")==0)  {
	  /* Compute file[idxfile]->LDOSmap which is a 2D x-y map containing the z coordinate  
	             at file[idxfile]->LDOS=Cst                                                */
	  /* ################################################################################# */
	  int i;
	  int idxfile=0;
	  int nlevel=50;
	  for(i=0;i<sline->ntokens;i++)	    {
	    if(strcmp(sline->words[i],"file")==0) 	  idxfile=atoi(sline->words[i+1]);
	    if(strcmp(sline->words[i],"nlevel")==0) 	  nlevel=atoi(sline->words[i+1]);
	  }
	  file[idxfile]->LDOS=1.0e-4;
	  double zmin,zmax,dz;
	  zmin=file[idxfile]->atom_lim[MINI][Z];
	  zmax=file[idxfile]->atom_lim[MAXI][Z];
	  PLFLT histmin,histmax;
	  histmin=file[idxfile]->O[Z];
	  histmax=file[idxfile]->O[Z]+file[idxfile]->DIM[Z];
	  for(i=0;i<sline->ntokens;i++){
	    if(strcmp(sline->words[i],"zmin")==0) 	  zmin=atof(sline->words[i+1]);
	    if(strcmp(sline->words[i],"zmax")==0) 	  zmax=atof(sline->words[i+1]);
	    if(strcmp(sline->words[i],"histmax")==0) 	  histmax=atof(sline->words[i+1]);
	    if(strcmp(sline->words[i],"histmin")==0) 	  histmin=atof(sline->words[i+1]);
	    if(strcmp(sline->words[i],"LDOS")==0) 	  file[idxfile]->LDOS=atof(sline->words[i+1]);
	  }

	  fprintf(stdout,"### FILE %d\n",idxfile);
	  fprintf(stdout,"### zmin %g\n",zmin);
	  fprintf(stdout,"### zmax %g\n",zmax);
	  fprintf(stdout,"### histmin %g\n",histmin);
	  fprintf(stdout,"### histmax %g\n",histmax);
	  fprintf(stdout,"### nlevel %d\n",nlevel);
	  fprintf(stdout,"NX=%d\nNY=%d\nNZ=%d\n",file[idxfile]->N[X],file[idxfile]->N[Y],file[idxfile]->N[Z]);
	  static int nx      ;        // Default number of data points in x
	  static int ny      ;        // Default number of data points in y
	  nx=file[idxfile]->N[X];
	  ny=file[idxfile]->N[Y];


	  // Transformation function
	  tr[0]=file[idxfile]->L[X][X];tr[1]=file[idxfile]->L[Y][X];tr[2]=file[idxfile]->O[X];
	  tr[3]=file[idxfile]->L[X][Y];tr[4]=-file[idxfile]->L[Y][Y];tr[5]=file[idxfile]->O[Y]+file[idxfile]->DIM[Y];
      
	  plsdev( "pdfcairo" );
	  plsfnam(file[idxfile]->pdffile ); 
	  //plspal0( "cmap0_black_on_white.pal" );
	  plspal1( "cmap1_gray.pal", 1 );
	  plscmap0n( 3 );
	  plinit();
	  int ix,iy,iz;

	  PLFLT *data;
	  PLINT NPTS;NPTS=(PLINT) file[idxfile]->N[X]*file[idxfile]->N[Y];
	  data=malloc(NPTS*sizeof(PLFLT));

	  FILE *out;
	  out=fopen("tmp","w+");
	  double x0,x1,y0,y1;
	  i=0;
	  /* the principle of constant_LDOS is to build a 2D map of z of constant LDOS
	     The search part consists in starting from the upper z, to descent until the
	     target LDOS (file[idxfile]->LDOS) value is reached.
	   */
	  for (iy=0;iy<file[idxfile]->N[Y];iy++) {
	    for (ix=0;ix<file[idxfile]->N[X];ix++) {
	      iz=file[idxfile]->N[Z]-1;
	      /*  */
	      file[idxfile]->LDOSmap[ix][iy]=0.0;
	      while(file[idxfile]->data[DATA][ix][iy][iz]<file[idxfile]->LDOS && iz>0) iz--;
	      if(iz>0) {
		x0=file[idxfile]->data[Z][ix][iy][iz-1];
		y0=file[idxfile]->data[DATA][ix][iy][iz-1];
		x1=file[idxfile]->data[Z][ix][iy][iz];
		y1=file[idxfile]->data[DATA][ix][iy][iz];
		file[idxfile]->LDOSmap[ix][iy]=(file[idxfile]->LDOS*(x1-x0)-x1*y0+x0*y1)/(y1-y0);
		data[i]=file[idxfile]->LDOSmap[ix][iy];i++;
	      }
	      fprintf(out,"%g %g %g\n",
		      file[idxfile]->data[X][ix][iy][iz],
		      file[idxfile]->data[Y][ix][iy][iz],
		      file[idxfile]->LDOSmap[ix][iy]);
	    }
	  }
	  fclose(out);
	  dz=(zmax-zmin)/(nlevel+1);
	  fprintf(stdout,"### dz %g\n",dz);
	  PLFLT *clevel; clevel=malloc(nlevel*sizeof(PLFLT));
	  for(i=0;i<nlevel;i++) 	clevel[i]=zmin+(i+1)*dz;
	  
	  pl_setcontlabelformat( 4, 3 );
	  pl_setcontlabelparam( 0.006, 0.3, 0.1, 1 );
	  /* set up the limits of the map */
	  PLFLT plenvxmin,plenvxmax,plenvymin,plenvymax;
	  plenvxmin=0.0;
	  plenvxmax=14.0;
	  plenvymin=4.0;
	  plenvymax=18.0;
	  for(i=0;i<sline->ntokens;i++){
	    if(strcmp(sline->words[i],"plenvxmin")==0) 	  plenvxmin=atof(sline->words[i+1]);
	    if(strcmp(sline->words[i],"plenvxmax")==0) 	  plenvxmax=atof(sline->words[i+1]);
	    if(strcmp(sline->words[i],"plenvymin")==0) 	  plenvymin=atof(sline->words[i+1]);
	    if(strcmp(sline->words[i],"plenvymax")==0) 	  plenvymax=atof(sline->words[i+1]);
	  }
	  plenv( plenvxmin,plenvxmax, plenvymin,plenvymax,1, 0 );

	  plcol0( 2 );
	  PLFLT         fill_width = 2., cont_width = 0.;
	  //      PLFLT         colorbar_width, colorbar_height;
	  PLINT         cont_color = 0;
	  plshades( (const PLFLT * const *) file[idxfile]->LDOSmap, nx, ny, NULL,1, nx, 1, ny, clevel, nlevel, fill_width,cont_color,cont_width,plfill,1,mypltr, NULL );

	  int scan=False;
	  double scanx0,scany0,scanx1,scany1;
	    int npts=5;
	    for(i=0;i<sline->ntokens;i++){
	    if(strcmp(sline->words[i],"scan")==0) 	  {
	      scan=True;
	      scanx0=atof(sline->words[i+1]);
	      scany0=atof(sline->words[i+2]);
	      scanx1=atof(sline->words[i+3]);
	      scany1=atof(sline->words[i+4]);
	      npts=atoi(sline->words[i+5]);
	    }
	  }

	  if(scan==True){
	    double *scanx,*scany,d,dcur0,dcur1;
	    int nscan=2;
	    scanx=malloc(nscan*sizeof(double));
	    scany=malloc(nscan*sizeof(double));
	    iz=0;ix=0;iy=0;
	    

	    

	    scanx[0]= file[idxfile]->data[X][ix][iy][iz];
	    scany[0]= file[idxfile]->data[Y][ix][iy][iz];
	    dcur0=sqrt(pow(scanx[0]-scanx0,2)+pow(scany[0]-scany0,2));
	    scanx[1]= file[idxfile]->data[X][ix][iy][iz];
	    scany[1]= file[idxfile]->data[Y][ix][iy][iz];
	    dcur1=sqrt(pow(scanx[1]-scanx1,2)+pow(scany[1]-scany1,2));
	    
	    for(ix=0;ix<file[idxfile]->N[X];ix++){
	      for(iy=0;iy<file[idxfile]->N[Y];iy++){
		d=sqrt(
		       pow(file[idxfile]->data[X][ix][iy][iz]-scanx0,2)+
		       pow(file[idxfile]->data[Y][ix][iy][iz]-scany0,2));
		if(d<dcur0){
		  scanx[0]= file[idxfile]->data[X][ix][iy][iz];
		  scany[0]= file[idxfile]->data[Y][ix][iy][iz];
		  dcur0=d;
		}
		d=sqrt(
		       pow(file[idxfile]->data[X][ix][iy][iz]-scanx1,2)+
		       pow(file[idxfile]->data[Y][ix][iy][iz]-scany1,2));
		if(d<dcur1){
		  scanx[1]= file[idxfile]->data[X][ix][iy][iz];
		  scany[1]= file[idxfile]->data[Y][ix][iy][iz];
		  dcur1=d;
		}
	      }
	    }
	    
	    fprintf(stdout,"line (%g %g) to ( %g %g)\n",scanx[0],scany[0],scanx[1],scany[1]);
	    plcol0( 2 );
	    //	    plline(nscan,scanx,scany);

	    double x,y;
	    double *px,*py;
	    px=malloc(npts*sizeof(double));
	    py=malloc(npts*sizeof(double));
	    dcur0=sqrt(pow(scanx[1]-scanx[0],2)+pow(scany[1]-scany[0],2));
	    d=dcur0/(npts-1);
	    px[0]=scanx[0];py[0]=scany[0];
	    for(i=1;i<npts-1;i++){
	      px[i]=px[0]+i*d*(scanx[1]-scanx[0])/dcur0;
	      py[i]=py[0]+i*d*(scany[1]-scany[0])/dcur0;
	      fprintf(stdout,"%d %g %g\n",i,px[i],py[i]);
	    }
	    px[npts-1]=scanx[1];py[npts-1]=scany[1];
	    plcol0( 2 );
	    plpoin(npts,px,py,9);

	    
	    double *ppx,*ppy;
	    int *ppix,*ppiy;
	    ppx=malloc(npts*sizeof(double));
	    ppy=malloc(npts*sizeof(double));
	    ppix=malloc(npts*sizeof(int));
	    ppiy=malloc(npts*sizeof(int));
	    for(i=0;i<npts;i++){
	      ix=0;iy=0;iz=0;
	      ppix[i]=ix;	      ppiy[i]=iy;
	      ppx[i]= file[idxfile]->data[X][ix][iy][iz];
	      ppy[i]= file[idxfile]->data[Y][ix][iy][iz];
	      dcur0=sqrt(pow(ppx[i]-px[i],2)+pow(ppy[i]-py[i],2));
	      for(ix=0;ix<file[idxfile]->N[X];ix++){
		for(iy=0;iy<file[idxfile]->N[Y];iy++){
		  d=sqrt(
			 pow(file[idxfile]->data[X][ix][iy][iz]-px[i],2)+
			 pow(file[idxfile]->data[Y][ix][iy][iz]-py[i],2));
		  if(d<dcur0){
		    ppx[i]= file[idxfile]->data[X][ix][iy][iz];
		    ppy[i]= file[idxfile]->data[Y][ix][iy][iz];
		    ppix[i]=ix;	      ppiy[i]=iy;
		    dcur0=d;
		  }
		}
	      }
	    }
	    for(i=0;i<npts;i++){
	      fprintf(stdout,"%g %g %g %d %d\n",
		      ppx[i],ppy[i],
		      file[idxfile]->LDOSmap[ppix[i]][ppiy[i]],
		      ppix[i],ppiy[i]);
	    }

	    plcol0( 1 );
	    plpoin(npts, ppx,ppy,9);

	    
	    
	  }
	  char *label;label=malloc(NCHAR*sizeof(char));
	  sprintf(label,"Constant LDOS (%g)",file[idxfile]->LDOS);
	  pllab( "x (ang.)", "y (ang.)", label );
	  free(label);


	  plcol0( 1 );
	  plhist( NPTS, data, histmin,histmax, 44, PL_HIST_IGNORE_OUTLIERS );
	  plcol0( 2 );
	  pllab( "#frValue", "#frFrequency",
		 "#frPLplot Example 5 - Probability function of Oscillator" );

	  
	  //	  plFree2dGrid( file[idxfile]->LDOSmap, nx,ny);



	  free(data);
	  plend();
	}
      }
      free(line);
    }
  }

  /* FILE *out; */
  /* out=fopen("CubeTools.hist","w+"); */
  /* int i;GSList *node; */
  /* for(i=0;node=g_slist_nth(history,i);i++){ */
  /*   g_print("%s\n",(char *)node->data); */
  /* } */
  /* fclose(out); */
  return EXIT_SUCCESS;
}
