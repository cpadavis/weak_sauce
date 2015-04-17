/**************************************************************
 *
 *
 *  	r2d.h
 *
 *		Devon Powell
 *		8 February 2015
 *
 * 		Routines for computing exact intersection volumes
 * 		and integrals between tetrahedra and regular meshes
 *
 *  	Copyright (C) 2014 Stanford University.
 *  	See License.txt for more information.
 *
 *      Some important usage notes:
 *
 *      - patch is a row-major array which must be allocated and zeroed prior 
 *        to calling this function, i.e. with calloc(). It will be filled in with
 *        the appropriate values. Just be sure to free() it when done.
 *
 *
 *************************************************************/

#ifndef _R2D_H_
#define _R2D_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

/*
 * The precision to be used in calculations
 */
#ifdef SINGLE_PRECISION
typedef float r2d_real;
#else
typedef double r2d_real;
#endif

/*
 * Integer types used for grid indexing and bit flags
 */
typedef int32_t r2d_int;
typedef int64_t r2d_long;

/*
 * Vector struct for vertex data
 */
typedef struct {
	r2d_real x, y;
} r2d_rvec2;

/*
 * Integer vector struct for passing grid indices
 */
typedef struct {
	r2d_int x, y;
} r2d_dvec2;

/*
 * Struct representing a plane
 */
typedef struct {
	r2d_rvec2 n;
	r2d_real d;
} r2d_plane;


/*
 * Struct representing a doubly-linked vertex
 */
typedef struct {
	r2d_rvec2 pos;
	unsigned char pnbrs[2];
	r2d_real fdist[4];
	unsigned char fflags;
} r2d_vertex;

/*
 * return type for voxelization functions
 * can be easily modified for debugging purposes
 */
typedef struct {
	r2d_int good;
	char* errmsg;
	// user added fields:
	r2d_real vtot;
	r2d_real vox_min;
	r2d_real vox_max;
} r2d_info;

/*
 * r2d
 *
 * Functions for voxelizing 2D primitives and convex polygons
 *
 */

r2d_info r2d_init(r2d_int bufsz); 
r2d_info r2d_set_dest_grid(r2d_real *dest, r2d_dvec2 dest_dims, r2d_rvec2 dest_window);
void r2d_finalize();
r2d_info r2d_rasterize_quad(r2d_real C, r2d_plane faces[4], r2d_dvec2 ibounds[2]);

/*
 * r2du
 *
 * Utility functions for r2d, including basic clipping operations
 *
 */

void r2du_clip_quad(r2d_vertex* vertbuffer, r2d_int* nverts, unsigned char andcmp);
void r2du_reduce(r2d_vertex* vertbuffer, r2d_int* nverts, r2d_real C, r2d_long flatind);

void r2du_init_box(r2d_vertex* vertbuffer, r2d_int* nverts, r2d_rvec2 rbounds[2]);
void r2du_faces_from_verts(r2d_rvec2* verts, r2d_int nverts, r2d_plane* faces);
void r2du_get_box_index_bounds(r2d_dvec2* ibounds,
				r2d_real xmin, r2d_real ymin,
				r2d_real xmax, r2d_real ymax);


// TO BE IMPLEMENTED
// typedef struct {
// 	r2d_rvec2 pos;
// 	r2d_int pnbrs[2];
// 	r2d_real fdist[64];
// 	r2d_ulong fflags;
// } r2d_vertex64;
//
// typedef uint64_t r2d_ulong;
//
// r2d_info r2d_rasterize_tri(r2d_real C, r2d_plane faces[4], r2d_dvec2 ibounds[2]);
// void r2du_clip_tri(r2d_vertex* vertbuffer, r2d_int* nverts, unsigned char andcmp);
// void r2du_clip_n(r2d_vertex* vertbuffer, r2d_int* nverts, r2d_plane* faces, r2d_int nfaces);

double
TestWoot(double x);

#endif // _R2D_H_
