PLPLOTINC=-I/home/bulou/ownCloud/src/PLplot/plplot-5.10.0/examples/c -I/usr/include/plplot
PLPLOTLIB=-lplplotd 
GLIB=`pkg-config --cflags --libs glib-2.0`
IMAGEMAGICK=`pkg-config --cflags --libs MagickWand`
GSLLIBPATH=/usr/local/lib
GSLLIBPATH=/usr/lib/x86_64-linux-gnu
GSLLIB=$(GSLLIBPATH)/libgslcblas.a $(GSLLIBPATH)/libgsl.a
SDLLIB=-lGL -lGLU -lSDL -lglut


INCLUDE=-I/usr/local/include $(PLPLOTINC) 
LIBRARY= $(GSLLIB)  $(GLIB)  $(PLPLOTLIB) $(IMAGEMAGICK) $(SDLLIB) -lm 
all: CubeTools.c
	gcc -Wall CubeTools.c $(INCLUDE) -o CubeTools $(LIBRARY) 
