/*
 * Plot and Save from DUMP files using BVIEW
 * 
 * You may want to build an animation further
 * $ convert -delay 10 -loop 0 *.png animation.gif
 */


/* Compile with:
 * 
 *      3D with MPI: CC99='mpicc -std=c99' qcc -Wall -O2 -grid=octree -D_MPI=1 plot3D.c -o plot3D -lm -L$BASILISK/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11
 *      
 *      2D with MPI: CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 plot.c -o plot -lm -L$BASILISK/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11
 * 
 *      2D without MPI: qcc -Wall -O2 plot.c -o plot -lm
 * 
 * Run with:
 * 
 *      with MPI: mpirun -np 2 ./plot3D
 *      without MPI: ./plot3D
 * 
 * Optional:
 * 
 * qcc -source plot.c
 * mpicc -O2 -Wall -std=c99 -D_MPI=1 _plot.c -o plot -lm
 * mpirun -np 2 ./plot
 * 
 */


#include "run.h"
#include "view.h"

#define mesh 512

int main() {
  
  run();
  
}

// TODO Automaticaly loop over dump files
#define numberOfFiles 10


event init(i = 0) {

  scalar f[];
  //scalar particle[];
  
  scalar concentration[];
  
  const double dump_number[numberOfFiles] = {0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0};
  char dump_name[20], image_name[20];
  
  for (int i=0; i<numberOfFiles; i++) {    
     
    sprintf(dump_name, "dump-%2.1f", dump_number[i]);
    printf("%s\n", dump_name);
    restore (file = dump_name, list = all); //list = {f, particle}); //list = all);
    
//     vector u[];
//     foreach() { 
//       u.x[] = q.x[]/rho[];
//       u.y[] = q.y[]/rho[];
//       #if dimension == 3
//       u.z[] = q.z[]/rho[];
//       #endif  
//     }
//     boundary((scalar *){u});

    #if dimension == 3  
    
    clear ();   
           
    /* The block following this command will be drawn in a translated coordinate system.
      struct _translate { x, y, z }; */
    //translate (y = 28.e-3);
    
    /* tx, ty: shifts the camera center point.
      fov: changes the field-of-view (angle of view [°]). The larger the field of view the more of the scene you can see, but each bit will be seen at a lower resolution.
      quat[]: the quaternion defining the camera angles.
      sx, sy, sz: stretch factors for each axis.
      width, height, samples: the width and height (in pixels) of the image to render (default is 800 x 800). The image can optionally be generated by first rendering an image with samples times more pixels in each direction followed by subsampling. This provides a form of antialiasing. Default is four samples.
      bg[]: an array of red, green, blue values between 0 and 1 which defines the background color.
      theta, phi, psi: Euler-like angles (in radians), used (instead of quat[]) to define the camera angle.
      relative: whether the theta and phi angles are absolute or relative to the current position (i.e. increments of the current angles).
      camera: predefined camera angles: “left”, “right”, “top”, “bottom”, “front”, “back” and “iso”.
      map: an optional coordinate mapping function. */     
     view (ty = -.01, fov = 5, width = mesh, height = mesh, camera = "iso", bg = {1, 1, 1});
    //view (ty = -0.0, fov = 45, width = mesh, height = mesh, camera = "iso", bg = {1, 1, 1});
    
    /* The block following this command will be drawn in a coordinate system symmetric relative to the given plane. The plane is given by n and α.
    sstruct _mirror { coord n; double alpha; }; */ 
    //mirror ();
        
    /* Displays box boundaries and axis coordinates
    notics: do not draw tick marks (default is false).
    lc[]: an array of red, green, blue values between 0 and 1 which defines the line color.
    lw: the line width. */  
    box (notics=true);
    //box();
       
    /* c: the name (as a string) of the Volume-Of-Fluid field.
      s: the (optional) name of the face fraction field.
      edges: whether to display the edges of the facets.
      larger: makes each cell larger by this factor. This helps close the gaps in the VOF interface representation. Default is 1.1 in 3D and when edges are not displayed, otherwise it is 1.
      filled: in 2D, whether to fill the inside (1) or outside (-1).
      color: use this field to color each interface fragment.
      min, max: the minimum and maximum values to use for color mapping.
      spread: the “spread factor” to use if min and max are not defined. The maximum and minimum values will be taken as the average plus or minus spread times the standard deviation. Default is 5. If negative, the minimum and maximum values of the field are used.
      linear: if true the color will be linearly interpolated for each vertex of the facet.
      map: the colormap to use. Default is jet. Options: https://www.npmjs.com/package/colormap
      fc[]: an array of red, green, blue values between 0 and 1 which defines the facet color.
      lc[]: an array of red, green, blue values between 0 and 1 which defines the line color.
      lw: the line width. */
//     draw_vof ("concentration", color = "concentration");
    draw_vof ("f", color = "temperature", min = 0, max = t);
//     draw_vof ("particle", color = "f", min = -10, max = 0);
      
    /* In 3D the intersections of the cells with a plane are displayed. The default plane is z=0. 
      This can be changed by setting n and alpha which define the plane: nx*x + ny*y + nz*z = alpha
      struct _cells { coord n; double alpha; float lc[3], lw; // the line color and width }; */    
    //cells ();
        
    /* The field name is given by color. The min, max, spread, map etc. arguments work as described in draw_vof().
      In 3D the intersections of the field with a plane are displayed. The default plane is z=0. 
      struct _squares {
        char * color;
        double min, max, spread;
        bool linear;
        colormap map;
        float fc[3], lc[3];  
        coord n;
        double α;
      }; */  
    //squares ("f", min = 0, max = 1);
    //squares ("particle", min = 0, max = 1);
   
    
     /* displays an isosurface of a field
      struct _isosurface {
        char * f; // the name (as a string) of the field.
        double v; // the value of the isosurface.

        char * color; // use this field to color each interface fragment.
        double min, max, spread;
        bool linear;
        colormap map; 
        float fc[3], lc[3];
      }; */ 
    //#if dimension == 3
    //  scalar l2[];
    //  lambda2 (u, l2);
    //  isosurface ("l2", -0.0001);
    //#endif

    /* Moves the camera to a different viewpoint
      start: starting time of the camera motion.
      end: time at which the viewpoint should be reached.
      tx, ty, quat, fov: definition of the target viewpoint. */
    //travelling ();

    /* Prints facets as log files */
    // output_facets (f, stderr);    

    /* Draws strings on a separate layer (for annotations)
      str: string to display.
      pos: position: “0” bottom left, “1” top left, “2” top right and “3” bottom right (default 0).
      size: the size of the text, given as the number of characters which can fit within the width of the screen. Default is 40.
      lc[]: an array of red, green, blue values between 0 and 1 which defines the text color.
      lw: the line width. */
    //draw_string ()
    
    sprintf (image_name, "iso-%2.1f.png", dump_number[i]);
    save (image_name);
        
        
    
    clear ();    
    //view (ty = -.01, fov = 25, width = mesh, height = mesh, camera = "front", bg = {1, 1, 1});
    view (ty = -.01, fov = 10, width = mesh, height = mesh, camera = "front", bg = {1, 1, 1});
    //view (ty = -0.0, fov = 45, width = mesh, height = mesh, camera = "front", bg = {1, 1, 1});
   box (notics=true);
    //cells ();    
     squares ("concentration", min=0.2, max=1);
//     draw_vof ("concentration", color = "concentration");
    //draw_vof ("f", color = "f", min = -10, max = 0);
     draw_vof ("f", edges=true);
    //draw_vof ("particle", color = "f", min = -10, max = 0);
    
    sprintf (image_name, "front-%2.1f.png", dump_number[i]);
    save (image_name);
    
    
    clear ();    
    view (ty = -0.01, fov = 25, width = mesh, height = mesh, camera = "bottom", bg = {1, 1, 1});
    box (notics=true);
    draw_vof ("f", color="temperature", min=0, max=1);
    cells(n = {0,1,0});    
    squares ("temperature", min=0, max=1, n = {0,1,0});
    
//     draw_vof ("concentration", color = "concentration");
    //draw_vof ("f", color = "f", min = -10, max = 0);
//    draw_vof ("particle", color = "f", min = -10, max = 0);
    
    sprintf (image_name, "bottom-%2.1f.png", dump_number[i]);
    save (image_name);
        
    #endif
           
  } 
  
}
  
  