/*
 * Compile with:
 * 
 *      3D with MPI: CC99='mpicc -std=c99' qcc -Wall -O2 -grid=octree -D_MPI=1 droplet.c -o droplet -lm -L$BASILISK/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11
 *      
 *      2D with MPI: CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 droplet.c -o droplet -lm -L$BASILISK/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11
 * 
 *      2D without MPI: qcc -Wall -O2 droplet.c -o droplet -lm
 * 
 * Warning: mkdir ./results
 * Run with:
 * 
 *      with MPI: mpirun -np 2 ./droplet
 *      without MPI: ./droplet
 * 
 * Optional:
 * 
 * qcc -source droplet.c
 * mpicc -O2 -Wall -std=c99 -D_MPI=1 _droplet.c -o droplet -lm
 * mpirun -np 2 ./droplet
 * 
 */

// Basilisk provides and extensive toolbox to do interface reconstruction. On this page we check the convergence properties for determining the surface area of an iso-surface countour of an analytical function.

// The sphere area is A = 4 PI R²
// is we define a sphere A = PI
// then the radius should be:
// PI = 4 PI R²
// 1/4 / R²
// R = 1/2

#include "grid/octree.h"
#include "fractions.h"
#include "utils.h"

# define func (sq(x-xo)+sq(y-yo)+sq(z-zo)-sq(R))

double xo=0.1,yo=M_PI/10.0;
double zo=-1./2.44,R=0.5;
scalar f[];
int main(){
  L0=5.;
  X0=Y0=Z0=-L0/2;
  init_grid(1<<4);
  int i = 0;
  static FILE * fp1 = fopen("interface.dat","w");
  for (i=1;i<8;i++){
    double s = 1./pow(2.,(double)i); 
    refine(fabs(func)<(s) && level < i+4);
    fraction(f,func);
    int cells=0;
    foreach()
      if (f[]>0. && f[]<1.)
	cells++;
    fprintf(fp1,"%d\t%d\t%g\n",i,cells,100*fabs(interface_area(f)-M_PI)/M_PI);
//     static FILE * fp = popen("gfsview-batch3D interface.interface.gfv | ppm2gif --delay 200 > surface.gif","w");
//     output_gfs(fp);
//     fprintf(fp, "Save stdout { format = PPM width = 600 height = 600}\n");
  }
  fclose(fp1);
}