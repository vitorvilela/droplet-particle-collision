#!/bin/bash
rm droplet
rm *.png
rm *.ppm
rm *.csv
CC99='mpicc -std=c99' qcc -Wall -O2 -grid=octree -D_MPI=1 droplet.c -o droplet -lm -L$BASILISK/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11
mpirun -np 6 ./droplet &

