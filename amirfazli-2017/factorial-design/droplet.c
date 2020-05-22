/*
 * Factorial Design
 * 
 * Droplet-Particle impact - Immersed Boundary: Imposed q - conservative approach
 *
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


#include "momentum.h"
#include "tension.h"
#include "view.h"

#if dimension == 3
#include "lambda2.h"
#endif

#define cf


// TODO Check mesh size
// Refinement level (level, n): (5, 32) ... (10, 1024)
int LEVEL = 5;
int LEVEL_MAX = 10;

// 28 × 8 mm domain
#define LARGER 2
#define DOMAIN_LENGTH LARGER*28.e-3

// Smallest grid size for domain length 2 * 28.e-3 mm 
// (level, micrometers): (9, 109) (10, 54) (11, 27)
// Particle diameter deviation of 10 micrometers would require (13, 6.8)
// Finest experimental film thickness of about 1 micrometer would require (15, 1.7)

// The front view images was set on a resolutionof 640 × 480 pixels
// Precision for the length measurement was determined to be 0.03 mm, it would require (11, 27.3)
// In the benchmark simulation, a 28 × 8 mm domain meshed with initial cells of 0.02 mm size in the impact area was used

// For this problem SMALL is considered a microsecond or a micrometer
#define SMALL 1.e-6

// Gravity acceleration
#define GRAVITY 9.81
#define PRESSURE 101000

// Fluid flow properties - Liquid: water, Gas: air, Solid: glass
#define RHO_L 1000.
#define RHO_G 1.
#define RHO_S 500.
#define MU_L 1.e-3
#define MU_G 1.e-6

// Glass beads with a diameter of 2 ± 0.01 mm were used as target particles. 
// In the numerical simulation, a 2 mm particle was placed at the centre of a 28 × 8 mm domain
scalar particle[];
// Nominal diameter
#define Dp 2.e-3
// Little changes on particle position could be used as experimental replicates,
// however it is not known this experimental uncertainties
#define Xp 0.
#define HIGHER 1
#define Yp HIGHER*0.6*DOMAIN_LENGTH
#define Up 0
#define Vp 0
#if dimension == 3
#define Zp 0
#define Wp 0
#endif

// Diameter of the droplet just before the impact
#define D0 3.3e-3
// Little changes on droplet shape is used as experimental replicates
// Ellipsoid factor
#define EF 0.9
// Reference length: droplet diameter = 3.3 mm
#define Dd D0
#define Xd 0.
const double Yd = Yp + 0.5*Dp + 0.55*Dd; // additional gap of 5% x Dd
#define Ud 0
#if dimension == 3
#define Zd 0
#define Wd 0
#endif

// TODO Check number of replicates and variations
#define replicates 3
#define variations 2

// Design parameters
double Reynolds = 0.;
double Weber = 0.;
int Replicate;

// Ellipsoid radius
double ra = 0.5*Dd;
double rb = 0.5*Dd;
#if dimension == 3
double rc = 0.5*Dd;
#endif

// Initial droplet volume and interfacial area
double initial_droplet_volume = 0.;
double initial_droplet_area = 0.;

// Initial droplet velocity
double Vd = 0.;

// Surface tension coefficient
double sigma = 0.;

// Both cameras were synchronized to capture 7000 fps
// The front view images was set an exposure time of 97 µs
// Plot and simulation time 
// TODO Check times
const double checkat = 1000.e-6;
const double duration = 100000.e-6;


static FILE * factorialFile;

int main() {
  
  L0 = DOMAIN_LENGTH;

  #if dimension == 2    
  origin(-0.5*DOMAIN_LENGTH, 0);
  #else
  origin(-0.5*DOMAIN_LENGTH, 0, -0.5*DOMAIN_LENGTH);
  #endif    
    
  rho1 = RHO_L;
  rho2 = RHO_G;  
  mu1 = MU_L;
  mu2 = MU_G;
  
  // Decreasing the tolerance on the Poisson solve improves the results 
  TOLERANCE = 1e-4;
        
  // Initial uniform mesh size 
  N = 1 << LEVEL;
  init_grid(N);
  
  // TODO Check ellipsoid factors, Reynolds and Weber numbers, and dimensionless timestep for plotting
  const double ef[replicates] = {1., 0.9, 0.8};
  const double reynolds[variations] = {1.e3, 1.e4};
  const double weber[variations] = {1.e2, 1.e3};
    
  for (int i=0; i<replicates; i++) {
    
    Replicate = i;
    
    rb = ef[i]*0.5*Dd; 
    
    // Relation amining at sphere-ellipsoid conserving volume 
    ra = sqrt(ef[i])*0.5*Dd;
    
    // Sphere volume and area
    #if dimension == 2
    initial_droplet_volume = pi*ra*rb;
    initial_droplet_area = pi*Dd;
    #else
    rc = ra;
    initial_droplet_volume = (4/3)*pi*ra*rb*rc;
    initial_droplet_area = pi*Dd*Dd;
    #endif
          
    for (int j=0; j<variations; j++) {
      
      Reynolds = reynolds[j];
      Vd = Reynolds*MU_L/(RHO_L*Dd);
      
      for (int k=0; k<variations; k++) {
        
        Weber = weber[k];
        sigma = RHO_L*Dd*Vd*Vd/Weber;
        f.sigma = sigma;
       
        (pid() == 0 ? printf("\nReplicate: %d   Re: %d   We: %d\n", Replicate, (int)(Reynolds), (int)(Weber)) : 0);
        
        run();
        
      }
      
    }
    
  }
  
}

// Boundary conditions
q.n[bottom] = neumann(0);
p[bottom] = dirichlet(PRESSURE);

q.n[top] = neumann(0);
p[top] = dirichlet(PRESSURE);

q.n[right] = neumann(0);
p[right] = dirichlet(PRESSURE);

q.n[left] = neumann(0);
p[left] = dirichlet(PRESSURE);

#if dimension == 3
q.n[front] = neumann(0);
p[front] = dirichlet(PRESSURE);

q.n[back] = neumann(0);
p[back] = dirichlet(PRESSURE);
#endif


event init(i = 0) {
  
  // TODO Check restart file dump
  if (!restore (file = "dump")) {  
    
    #if dimension == 2    
    coord vd = {Ud, -Vd};    
    refine ( ( sq(x-Xp) + sq(y-Yp) < sq(1.1*0.5*Dp) ||  sq(x-Xd)/sq(1.1*ra) + sq(y-Yd)/sq(1.1*rb) < 1. ) && level < LEVEL_MAX );
    fraction ( particle, sq(0.5*Dp) - sq(x-Xp) - sq(y-Yp) );
    fraction ( f, 1. - sq(x-Xd)/sq(ra) - sq(y-Yd)/sq(rb) );    
    #else    
    coord vd = {Ud, -Vd, Wd};    
    refine ( ( sq(x-Xp) + sq(y-Yp) + sq(z-Zp) < sq(1.1*0.5*Dp) ||  sq(x-Xd)/sq(1.1*ra) + sq(y-Yd)/sq(1.1*rb) + sq(z-Zd)/sq(1.1*rc) < 1. ) && level < LEVEL_MAX );    
    fraction ( particle, sq(0.5*Dp) - sq(x-Xp) - sq(y-Yp) - sq(z-Zp) );
    fraction ( f, 1. - sq(x-Xd)/sq(ra) - sq(y-Yd)/sq(rb) - sq(z-Zd)/sq(rc) );    
    #endif
    
    foreach() 
      foreach_dimension()      
        q.x[] = f[]*rho1*vd.x;
    boundary ((scalar *){q}); 
    
#ifndef cf
    foreach()
      cf[] = f[];
#endif
    
  }
  
}


event stationary_particle (i++) {
  
  #if dimension == 2  
  coord vp = {Up, Vp};  
  fraction ( particle, sq(0.5*Dp) - sq(x-Xp) - sq(y-Yp) );  
  #else  
  coord vp = {Up, Vp, Wp};  
  fraction ( particle, sq(0.5*Dp) - sq(x-Xp) - sq(y-Yp) - sq(z-Zp) );  
  #endif

  foreach()
    foreach_dimension()
      q.x[] = particle[]*RHO_S*vp.x + (1. - particle[])*q.x[];
  boundary ((scalar *){q});
  
}


event acceleration (i++) {  
  face vector av = a;
  foreach_face(y)
    av.y[] -= GRAVITY;   
}


event adapt (i++) {

  vector u[];
  foreach()
    foreach_dimension()
      u.x[] = q.x[]/rho[];
  boundary((scalar *){u});
  
#if dimension == 2
  adapt_wavelet({f, particle, u}, (double[]){1.e-3, 1.e-3, 1.e-2, 1.e-2, 1.e-2}, LEVEL_MAX);
#else
  adapt_wavelet({f, particle, u}, (double[]){1.e-3, 1.e-3, 1.e-2, 1.e-2, 1.e-2}, LEVEL_MAX);  
#endif  
  
  //unrefine (y<(1/2)*DOMAIN_LENGTH);
  
}


// Factorial Design File
event factorial (t=0; t<duration; t+=checkat) {
  
  double area = interface_area(f);
  
  // Header: Replicate, Re, We, t*, A*
  
  if (pid() == 0) {
    
    if (i == 0)
      factorialFile = fopen ("factorialDesign.csv", "a+");
    
    if ( t == 0 || fabs(t-0.5*Dd/Vd) < checkat || fabs(t-1.0*Dd/Vd) < checkat || fabs(t-1.5*Dd/Vd) < checkat ) {
      fprintf (factorialFile, "%d, %d, %d, %1.3f, %g\n", Replicate, (int)(Reynolds), (int)(Weber), t*Vd/Dd, area/initial_droplet_area);
      fflush (factorialFile);
    }
    
  }
  
  if (t > 1.51*Dd/Vd)
    return 1;
  
}


event snapshot (t=0; t<duration; t+=checkat) {
  
  if ( t == 0 || fabs(t-0.5*Dd/Vd) < checkat || fabs(t-1.0*Dd/Vd) < checkat || fabs(t-1.5*Dd/Vd) < checkat ) {
  
    // We dump a snapshot which can be used to restart the simulation
    char dump_name[80];
    sprintf (dump_name, "./dumps/dump-R%d-Re%d-We%d-%1.1f", Replicate, (int)(Reynolds), (int)(Weber), t*Vd/Dd);
    dump (file = dump_name);
    
    char ppm_name_vof[80];
    sprintf(ppm_name_vof, "./figs/vof-R%d-Re%d-We%d-%1.1f.ppm", Replicate, (int)(Reynolds), (int)(Weber), t*Vd/Dd); 

    // VoF cut-off
    scalar ff[];
    foreach() {
      ff[] = f[] < 1e-4 ? 0 : f[] > 1. - 1e-4 ? 1. : f[];
    }
    boundary({ff});
      
    // Cells for which m is negative will be black in the picture  
    scalar m[];
    foreach()
      m[] = 0.5 - particle[];
    boundary ({m});
        
    output_ppm(ff, file=ppm_name_vof, mask=m, n=1<<LEVEL_MAX);
      
// Doesn't work remotely    
//     #if dimension == 3
//     char name_iso[80]; 
//     sprintf(name_iso, "./figs/iso-R%d-Re%d-We%d-%1.1f.png", Replicate, (int)(Reynolds), (int)(Weber), t*Vd/Dd);
//       
//     char name_front[80]; 
//     sprintf(name_front, "./figs/front-R%d-Re%d-We%d-%1.1f.png", Replicate, (int)(Reynolds), (int)(Weber), t*Vd/Dd);
//       
//     char name_bottom[80]; 
//     sprintf(name_bottom, "./figs/bottom-R%d-Re%d-We%d-%1.1f.png", Replicate, (int)(Reynolds), (int)(Weber), t*Vd/Dd);
//       
//     
//     clear (); 
//     view (ty = -.35, fov = 10, width = 2048, height = 2048, camera = "iso", bg = {1, 1, 1});
//     draw_vof ("f", color = "f", min = 0, max = 10);
//     draw_vof ("particle", color = "f", min = -10, max = 0);
//     save (name_iso);
//           
//     clear ();    
//     view (ty = -.5, fov = 25, width = 2048, height = 2048, camera = "front", bg = {1, 1, 1});
//     box (notics=true);
//     cells ();
//     draw_vof ("f", color = "f", min = 0, max = 10);
//     draw_vof ("particle", color = "f", min = -10, max = 0);
//     save (name_front);
//       
//     clear ();    
//     view (ty = -0.01, fov = 20, width = 2048, height = 2048, camera = "bottom", bg = {1, 1, 1});
//     box (notics=true);
//     cells ();
//     draw_vof ("f", color = "f", min = 0, max = 10);
//     draw_vof ("particle", color = "f", min = -10, max = 0);
//     save (name_bottom);
//           
//     #endif
    
  }
  
}
   