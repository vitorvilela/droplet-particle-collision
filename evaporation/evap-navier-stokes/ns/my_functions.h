#ifndef BASILISK_HEADER_15
#define BASILISK_HEADER_15
#line 1 "./../../my_functions.h"
coord normal (Point point, scalar c) {
  coord n = mycs (point, c);
  double nn = 0.;
  foreach_dimension()
    nn += sq(n.x);
  nn = sqrt(nn);
  foreach_dimension()
    n.x /= nn;
  return n;
}


void compute_normal (scalar f, vector normal_vector) {
  foreach() {
    coord n = normal (point, f);
    foreach_dimension() 
      normal_vector.x[] = n.x;
  }
  boundary((scalar*){normal_vector});
}


void magnet (scalar f, double error) {
  foreach() {
    f[] = clamp(f[], 0., 1.);
    f[] = (f[] < error ? 0. : (f[] > 1. - error ? 1. : f[]));
  }
  boundary ({f});
}


void my_laplacian (scalar f, scalar l, face vector D) {
  boundary({f, D});
  foreach() {
    l[] = 0.;
    foreach_dimension()
      l[] += (f[1] - f[0])*D.x[1] - (f[] - f[-1])*D.x[];
    l[] /= sq(Delta);
  }
  boundary({l});
}


double interface_length (Point point, scalar c)
{
  coord n = mycs (point, c);
  double alpha = line_alpha (c[], n);
  coord coord_centroid = {0, 0};
  return line_length_center(n, alpha, &coord_centroid);
}



#endif
