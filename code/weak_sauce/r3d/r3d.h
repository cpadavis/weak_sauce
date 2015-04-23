/**************************************************************
 *
 *    r3d.h
 *    
 *    Devon Powell
 *    8 February 2015
 *    
 *    Routines for computing exact intersection volumes
 *    and integrals between tetrahedra and regular meshes
 *    
 *    Copyright (C) 2014 Stanford University.
 *    See License.txt for more information.
 *
 *************************************************************/

/**
 * \file r3d.h
 * \author Devon Powell
 * \date 8 Feb 2015
 * \brief Interface for r3d
 */

#ifndef _R3D_H_
#define _R3D_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

/**
 * \brief Real type specifying the precision to be used in calculations
 *
 * Default is `double` (recommended). `float` precision is enabled by 
 * compiling with `-DSINGLE_PRECISION`.
 */
#ifdef SINGLE_PRECISION
typedef float r3d_real;
#else 
typedef double r3d_real;
#endif

/** \struct r3d_rvec3
 *  \brief Vector struct.
 */
typedef struct {
	r3d_real x, /*!< \f$x\f$-component. */
			 y, /*!< \f$y\f$-component. */
			 z; /*!< \f$z\f$-component. */
} r3d_rvec3;

/** \struct r3d_plane
 *  \brief A plane.
 */
typedef struct {
	r3d_rvec3 n; /*!< Unit-length normal vector. */
	r3d_real d; /*!< Signed perpendicular distance to the origin. */
} r3d_plane;

/** \struct r3d_orientation
 *  \brief Perpendicular distances and bit flags for up to 6 faces.
 */
typedef struct {
	//r3d_real fdist[6];
	r3d_real fdist[4];
	unsigned char fflags;
} r3d_orientation;

/* \struct r3d_vertex
 * \brief A doubly-linked vertex.
 */
typedef struct {
	r3d_rvec3 pos;
	unsigned char pnbrs[3];
	r3d_orientation orient;
} r3d_vertex;

// to be implemented?
#if 0
/* \struct r3d_poly
 * \brief A polyhedron.
 */
typedef struct {
	r3d_vertex* verts;
	r3d_orientation* orient;
	r3d_int nverts;
	r3d_int funclip;
} r3d_poly;
#endif

/**
 * \brief Integer type used for grid indexing and bit flags
 */
typedef int32_t r3d_int;

/**
 * \brief Long integer type used for grid indexing
 */
typedef int64_t r3d_long;

/** \struct r3d_dvec3
 *  \brief Integer vector struct for grid indexing.
 */
typedef struct {
	r3d_int i, /*!< \f$x\f$-index. */
			j, /*!< \f$y\f$-index. */
			k; /*!< \f$z\f$-index. */
} r3d_dvec3;

/** \struct r3d_dest_grid
 *  \brief Destination grid information.
 */
typedef struct {
	r3d_real* moments[10]; /*!< Gridded moments in row-major (C) order.
							   Must be at least n.x*n.y*n.z in size. */ 
	r3d_int polyorder; /*!< Polynomial order (0 for constant, 1 for linear, 2 for quadratic) to be voxelized. */
	r3d_orientation* orient; /*!< Buffer for gridpoint orientation checks.
							   Must be at least (n.x+1)*(n.y+1)*(n.z+1) in size.
							   Unused if compiled with -DUSE_TREE.*/ 
	r3d_long bufsz; /*!< Allocated size of moments and orient buffers. */
	r3d_dvec3 n; /*!< Grid dimensions (cells per coordinate). */
	r3d_rvec3 d; /*!< Grid cell size. */
} r3d_dest_grid;

/**
 * \brief Tabulated number of moments needed for a given polynomial order.
 */
const static r3d_int r3d_num_moments[3] = {1, 4, 10};

/**
 * \brief Voxelize a tetrahedron to the destination grid
 *
 * Conservatively voxelizes a tetrahedron to the destination grid (set with r3d_set_dest_grid()).
 *
 * \param [in] faces
 * The four faces of the tetrahedron to be voxelized.
 *
 */
void r3d_voxelize_tet(r3d_plane faces[4], r3d_dest_grid grid);


void r3d_clip_tet(r3d_vertex* vertbuffer, r3d_int* nverts, unsigned char andcmp);


/**
 * \brief Integrate a polynomial density over a convex polyhedron using simplicial decomposition
 *
 * \param [in] vertbuffer
 * Coefficients for the input density, which is a polynomial of up to order 2.
 *
 * \param [in] polyorder
 * Order of the polynomial density specified in coeffs. 0 for constant (1 coefficient), 1 for linear
 * (4 coefficients), 2 for quadratic (10 coefficients).
 *
 * \param [in] faces
 * The four faces of the tetrahedron to be voxelized.
 *
 * \param [in] ibounds
 * Index range of the destination grid to be considered in the voxelization process, specified as
 * the lower and upper corners of a box. Indices must lie inside the destination grid, as described
 * in r3d_set_dest_grid(). A convenience function is provided in r3du_get_box_index_bounds().
 *
 * \return r3d_info containing the current state of r3d
 *
 */

void r3d_reduce(r3d_vertex* vertbuffer, r3d_int* nverts, r3d_real* moments, r3d_int polyorder);

/*
 * r3du
 *
 * Utility functions for r3d.
 *
 */

void r3du_init_box(r3d_vertex* vertbuffer, r3d_int* nverts, r3d_rvec3 rbounds[2]);
void r3du_tet_faces_from_verts(r3d_rvec3* verts, r3d_plane* faces);
r3d_real r3du_orient(r3d_rvec3 pa, r3d_rvec3 pb, r3d_rvec3 pc, r3d_rvec3 pd);

#endif // _R3D_H_
