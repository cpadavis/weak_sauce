#############################################################################
#
#	Makefile
#
#	libr3d.a - a static library containing just three simple functions
#	to voxelize tetrahedra with up to 2nd-order (quadratic) density fields.
#
#	See Readme.md and r3d.h for usage
#
#############################################################################

######## User options #########

# Use an octree to efficiently search enclosed regions
#OPT += -DUSE_TREE

# Turn off cell shifting (for precision checks)
#OPT += -DNO_SHIFTING

# Turn off clipping (for timing checks)
#OPT += -DNO_CLIPPING

# Turn off computation of vertex terms (for timing checks)
#OPT += -DNO_REDUCTION

# For making POV-Ray figures
#OPT += -DPRINT_POVRAY_TRAVERSAL
#OPT += -DPRINT_POVRAY_TREE
#OPT += -DPRINT_POVRAY_DECOMP

###############################

CC = gcc
CFLAGS = -Wall -W -g -fPIC

# super lazy let's just do this my way

all: r2d_lib

r2d_lib: r2d.o
	$(CC) $(CFLAGS) -shared -Wl,-install_name,r2d_lib.so -o r2d_lib.so r2d.o

r2d.o: r2d.c
	$(CC) $(CFLAGS) -c -o r2d.o r2d.c

clean:
	rm -f *.so *.o
