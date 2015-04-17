/*
 *  
 *
 *  	r2d.c
 *
 *  	See Readme.md and r2d.h for usage.
 * 
 *  	Copyright (C) 2014 Stanford University.
 *  	See License.txt for more information.
 *
 *
 */

#include "r2d.h"

// tells us which bit signals a clipped vertex
// the last one - leaving seven to flag faces
#define CLIP_MASK 0x80 
#define CLIP_MASK_64 0x8000000000000000 

// macros for vector manipulation
#define dot(va, vb) (va.x*vb.x + va.y*vb.y)
#define wav(va, wa, vb, wb, vr) {			\
	vr.x = (wa*va.x + wb*vb.x)/(wa + wb);	\
	vr.y = (wa*va.y + wb*vb.y)/(wa + wb);	\
}
#define norm(v) {					\
	r2d_real tmplen = sqrt(dot(v, v));	\
	v.x /= (tmplen + 1.0e-299);		\
	v.y /= (tmplen + 1.0e-299);		\
}

double
TestWoot(double x)
{
    double x2;
    x2 = x * 2;
    return x2 + 1;
}

// dest grid (global for now)
static r2d_real* r2d_dest_grid;
static r2d_dvec2 r2d_dest_dims;
static r2d_rvec2 r2d_d;

// return information (global)
static r2d_info current_info;

//////////////// r2du: utility functions for r2d ///////////////////
//////////////////////////////////////////////////////////////////

inline void r2du_clip_quad(r2d_vertex* vertbuffer, r2d_int* nverts, unsigned char andcmp) {

	unsigned char v, f, ff;
	unsigned char vstart, vnext, vcur, searching;
	unsigned char fmask, ffmask;

	// loop over faces
	for(f = 0; f < 4; ++f) {

		fmask = (1 << f);
		if(andcmp & fmask) continue;

		// find the first vertex lying outside of the face
		// only need to find one (taking advantage of convexity)
		// TODO: add a smarter graph traversal here!
		vcur = 127;
		for(v = 0; vcur >= 127 && v < *nverts; ++v) 
			if(!(vertbuffer[v].fflags & (CLIP_MASK | fmask))) vcur = v;
		if(vcur >= 127) continue;
				
/*
		// traverse to the vertex lying oustide of the clip plane
		while(r2d_vertbuffer[vcur].fflags & fmask) {
			if(r2d_vertbuffer[vcur].fdist)
			if(r2d_vertbuffer[r2d_vertbuffer[vcur].pnbrs[0]].fdist[f] < r2d_vertbuffer[r2d_vertbuffer[vcur].pnbrs[1]].fdist[f]) {
				vcur = r2d_vertbuffer[vcur].pnbrs[0];
			}
			else {
				vcur = r2d_vertbuffer[vcur].pnbrs[1];
			}

		}
*/

		// traverse to the right
		vstart = vcur;
		searching = 1;
		while(searching) { 
			vnext = vertbuffer[vcur].pnbrs[1];
			if(fmask & vertbuffer[vnext].fflags) {
				// vnext is inside the face

				// compute the intersection point using a weighted
				// average of perpendicular distances to the plane
				wav(vertbuffer[vcur].pos, vertbuffer[vnext].fdist[f],
					vertbuffer[vnext].pos, -vertbuffer[vcur].fdist[f],
					vertbuffer[*nverts].pos);
				
				vertbuffer[*nverts].pnbrs[1] = vnext;

				// reciprocal connetivity
				vertbuffer[vnext].pnbrs[0] = *nverts;

				// do face intersections and flags
				vertbuffer[*nverts].fflags = 0x00;
				for(ff = f + 1; ff < 4; ++ff) {

					// skip if all initial verts are inside ff
					ffmask = (1 << ff); 
					if(andcmp & ffmask) continue;

					// weighted average keeps us in a relative coordinate system
					
					vertbuffer[*nverts].fdist[ff] = (vertbuffer[vcur].fdist[ff]*vertbuffer[vnext].fdist[f] 
							- vertbuffer[vnext].fdist[ff]*vertbuffer[vcur].fdist[f])
							/(vertbuffer[vnext].fdist[f] - vertbuffer[vcur].fdist[f]);
					if(vertbuffer[*nverts].fdist[ff] > 0.0) vertbuffer[*nverts].fflags |= ffmask;
				}
				++(*nverts);
				searching = 0;
			}
			vertbuffer[vcur].fflags |= CLIP_MASK;
			vcur = vnext;
		}

		// traverse to the left
		vcur = vstart;
		searching = 1;
		while(searching) { 
			vnext = vertbuffer[vcur].pnbrs[0];
			if(fmask & vertbuffer[vnext].fflags) {
				// vnext is inside the face

				// compute the intersection point using a weighted
				// average of perpendicular distances to the plane
				wav(vertbuffer[vcur].pos, vertbuffer[vnext].fdist[f],
					vertbuffer[vnext].pos, -vertbuffer[vcur].fdist[f],
					vertbuffer[*nverts].pos);
				
				vertbuffer[*nverts].pnbrs[0] = vnext;

				// reciprocal connetivity
				vertbuffer[vnext].pnbrs[1] = *nverts;

				// do face intersections and flags
				vertbuffer[*nverts].fflags = 0x00;
				for(ff = f + 1; ff < 4; ++ff) {

					// skip if all verts are inside ff
					ffmask = (1 << ff); 
					if(andcmp & ffmask) continue;

					// weighted average keeps us in a relative coordinate system
					vertbuffer[*nverts].fdist[ff] = (vertbuffer[vcur].fdist[ff]*vertbuffer[vnext].fdist[f] 
							- vertbuffer[vnext].fdist[ff]*vertbuffer[vcur].fdist[f])
							/(vertbuffer[vnext].fdist[f] - vertbuffer[vcur].fdist[f]);
					if(vertbuffer[*nverts].fdist[ff] > 0.0) vertbuffer[*nverts].fflags |= ffmask;

				}
				++(*nverts);
				searching = 0;
			}
			vertbuffer[vcur].fflags |= CLIP_MASK;
			vcur = vnext;
		}
		vertbuffer[*nverts-2].pnbrs[0] = *nverts-1; 
		vertbuffer[*nverts-1].pnbrs[1] = *nverts-2;
	}
}

inline void r2du_reduce(r2d_vertex* vertbuffer, r2d_int* nverts, r2d_real C, r2d_long flatind) {

	r2d_real tetvol, locvol;
	unsigned char v;
	r2d_rvec2 v0, v1; 

	// iterate over vertices and compute sum for this voxel
	locvol = 0.0;

	// finally, get the moments of each decomposed simplex
	for(v = 0; v < *nverts; ++v) {
		if(vertbuffer[v].fflags & CLIP_MASK) continue;
		v0 = vertbuffer[v].pos;
		v1 = vertbuffer[vertbuffer[v].pnbrs[1]].pos;

		tetvol = 0.5*(v0.x*v1.y - v0.y*v1.x); 
		locvol += tetvol; 

		/*
		 * HIGHER_ORDER MOMENTS, TO BE IMPLEMENTED.
		*/
	}

	r2d_dest_grid[flatind] += C*locvol;

	current_info.vtot += locvol;
	if(locvol < current_info.vox_min) current_info.vox_min = locvol;
	if(locvol > current_info.vox_max) current_info.vox_max = locvol;

}

inline void r2du_init_box(r2d_vertex* vertbuffer, r2d_int* nverts, r2d_rvec2 rbounds[2]) {
	*nverts = 4;
	vertbuffer[0].pnbrs[0] = 3;	
	vertbuffer[0].pnbrs[1] = 1;	
	vertbuffer[1].pnbrs[0] = 0;	
	vertbuffer[1].pnbrs[1] = 2;	
	vertbuffer[2].pnbrs[0] = 1;	
	vertbuffer[2].pnbrs[1] = 3;	
	vertbuffer[3].pnbrs[0] = 2;	
	vertbuffer[3].pnbrs[1] = 0;	
	vertbuffer[0].pos.x = rbounds[1].x; 
	vertbuffer[0].pos.y = rbounds[0].y; 
	vertbuffer[1].pos.x = rbounds[1].x; 
	vertbuffer[1].pos.y = rbounds[1].y; 
	vertbuffer[2].pos.x = rbounds[0].x; 
	vertbuffer[2].pos.y = rbounds[1].y; 
	vertbuffer[3].pos.x = rbounds[0].x; 
	vertbuffer[3].pos.y = rbounds[0].y; 
}



inline void r2du_faces_from_verts(r2d_rvec2* verts, r2d_int nverts, r2d_plane* faces) {

	// TODO: warn for convexity?

	// compute unit face normals and distances to origin
	r2d_rvec2 v0, v1, tmpcent;
	r2d_int f;

	// Assumes vertices are CCW
	for(f = 0; f < nverts; ++f) {
		v0 = verts[f];
		v1 = verts[(f+1)%nverts];
		faces[f].n.x = -(v1.y - v0.y);
		faces[f].n.y = (v1.x - v0.x);
		norm(faces[f].n);
		tmpcent.x = 0.5*(v0.x + v1.x);
		tmpcent.y = 0.5*(v0.y + v1.y);
		faces[f].d = -dot(faces[f].n, tmpcent);
	}
}

inline void r2du_get_box_index_bounds(r2d_dvec2* ibounds, r2d_real xmin, r2d_real ymin,
										r2d_real xmax, r2d_real ymax) {
	// get integer bounds
	ibounds[0].x = floor(xmin/r2d_d.x);
	ibounds[0].y = floor(ymin/r2d_d.y);
	ibounds[1].x = ceil(xmax/r2d_d.x);
	ibounds[1].y = ceil(ymax/r2d_d.y);

	// and clamp to the projection window
	if(ibounds[0].x < 0) ibounds[0].x = 0;
	if(ibounds[0].y < 0) ibounds[0].y = 0;
	if(ibounds[1].x > r2d_dest_dims.x) ibounds[1].x = r2d_dest_dims.x;
	if(ibounds[1].y > r2d_dest_dims.y) ibounds[1].y = r2d_dest_dims.y;

}

//////////////////////////////////////////
/// r2d //////////////////////////////////
//////////////////////////////////////////

typedef struct {
	r2d_real fdist[4];
	unsigned char fflags;
} r2d_gridpoint;

// typedef struct {
// 	r2d_real fdist[64];
// 	r2d_ulong fflags;
// } r2d_gridpoint64;

#ifdef USE_TREE
// tree node for recursive splitting
typedef struct {
	r2d_int imin, jmin;
	r2d_int ioff, joff;
	r2d_gridpoint gridpt[4];
} r2d_treenode;

// long version for many clip faces
// typedef struct {
// 	r2d_int imin, jmin;
// 	r2d_int ioff, joff;
// 	r2d_gridpoint gridpt64[4];
// } r2d_treenode64;

#else // !USE_TREE
static r2d_long r2d_gblen;
static r2d_gridpoint* r2d_gridbuffer;
// static r2d_gridpoint64* r2d_gridbuffer64;
#endif // USE_TREE 

r2d_info r2d_init(r2d_int bufsz) {

	// start keeping information
	current_info.good = 1;

#ifndef USE_TREE
	// initial allocation. This buffer is resized on the fly if needed
	r2d_gblen = bufsz; 
	r2d_gridbuffer = (r2d_gridpoint*) malloc(r2d_gblen*sizeof(r2d_gridpoint));
//	r2d_gridbuffer64 = (r2d_gridpoint64*) malloc(r2d_gblen*sizeof(r2d_gridpoint64));
	if(!r2d_gridbuffer /* || !r2d_gridbuffer64 */) {
		current_info.good = 0;
		current_info.errmsg = "Bad allocation of grid buffer";
	}
#endif
	
	return current_info;

}

void r2d_finalize() {
#ifndef USE_TREE
	if(r2d_gridbuffer) free(r2d_gridbuffer);
	//if(r2d_gridbuffer64) free(r2d_gridbuffer64);
#endif
}


r2d_info r2d_set_dest_grid(r2d_real* dest, r2d_dvec2 dims, r2d_rvec2 window) {
	r2d_dest_grid = dest;
	r2d_dest_dims = dims;
	r2d_d.x = window.x/r2d_dest_dims.x; 
	r2d_d.y = window.y/r2d_dest_dims.y; 
	if(dims.x <= 0 || dims.y <= 0 || window.x <= 0.0 || window.y <= 0.0) {
		current_info.good = 0;
		current_info.errmsg = "Invalid grid dimensions";
	}
	return current_info;
}


///// WARNING! MODIFIES THE INPUT FACES! ///////
r2d_info r2d_rasterize_quad(r2d_real C, r2d_plane faces[4], r2d_dvec2 ibounds[2]) {

	// variables used in this function
	r2d_real gor, locvol;
	r2d_rvec2 gpt;
	r2d_int i, j, nx, ny;
	unsigned char v, f;
	unsigned char orcmp, andcmp;

// macros for grid access
#define vind(i, j) (r2d_dest_dims.y*(i + ibounds[0].x) + (j + ibounds[0].y))
#ifdef USE_TREE
	// stack for the tree
	r2d_treenode treestack[256];
	r2d_int ntreestack;
	r2d_treenode curnode;
#else
#define gind(i, j) ((ny+1)*(i) + (j))
	r2d_long vv[4];
	r2d_long vv0;
#endif

	// voxel bounds in shifted coordinates
	r2d_rvec2 rbounds[2] = {
		{-0.5*r2d_d.x, -0.5*r2d_d.y}, {0.5*r2d_d.x, 0.5*r2d_d.y}
	};
	// vertex buffer
	r2d_int nverts;
	r2d_vertex vertbuffer[128];

	// start keeping track for this quad
	current_info.good = 1;
	current_info.vtot = 0.0;
	current_info.vox_min = 1.0e99;
	current_info.vox_max = -1.0e99;

	
	// shift the faces to align with
	// the working range of indices
	// THIS CHANGES THE INPUT FACES
	for(f = 0; f < 4; ++f) {
		faces[f].d += ibounds[0].x*r2d_d.x*faces[f].n.x
							+ ibounds[0].y*r2d_d.y*faces[f].n.y;
	}
	nx = ibounds[1].x - ibounds[0].x; 
	ny = ibounds[1].y - ibounds[0].y; 

#ifdef USE_TREE
	
	// get the initial face orientations for each corner of the node
	gpt.x = nx*r2d_d.x;
	gpt.y = 0.0;
	curnode.gridpt[0].fflags = 0x00;
	for(f = 0; f < 4; ++f) {
		gor = faces[f].d + dot(gpt, faces[f].n);
		if(gor > 0.0) curnode.gridpt[0].fflags |= (1 << f);
		curnode.gridpt[0].fdist[f] = gor;
	}
	gpt.y = ny*r2d_d.y;
	curnode.gridpt[1].fflags = 0x00;
	for(f = 0; f < 4; ++f) {
		gor = faces[f].d + dot(gpt, faces[f].n);
		if(gor > 0.0) curnode.gridpt[1].fflags |= (1 << f);
		curnode.gridpt[1].fdist[f] = gor;
	}
	gpt.x = 0.0;
	curnode.gridpt[2].fflags = 0x00;
	for(f = 0; f < 4; ++f) {
		gor = faces[f].d + dot(gpt, faces[f].n);
		if(gor > 0.0) curnode.gridpt[2].fflags |= (1 << f);
		curnode.gridpt[2].fdist[f] = gor;
	}
	gpt.y = 0.0;
	curnode.gridpt[3].fflags = 0x00;
	for(f = 0; f < 4; ++f) {
		gor = faces[f].d + dot(gpt, faces[f].n);
		if(gor > 0.0) curnode.gridpt[3].fflags |= (1 << f);
		curnode.gridpt[3].fdist[f] = gor;
	}

	curnode.imin = 0;
	curnode.jmin = 0;
	curnode.ioff = nx;
	curnode.joff = ny;

	ntreestack = 0;
	treestack[ntreestack++] = curnode;
	while(ntreestack > 0) {

		// pop the top node
		curnode = treestack[--ntreestack];

		orcmp = 0x00;
        andcmp = 0x0f;
		for(v = 0; v < 4; ++v) {
			orcmp |= curnode.gridpt[v].fflags; 
			andcmp &= curnode.gridpt[v].fflags; 
		}

		if(andcmp == 0x0f) {
			
			// all cells in this leaf are fully contained
#ifndef NO_REDUCTION
			locvol = r2d_d.x*r2d_d.y;
			
			for(i = curnode.imin; i < curnode.imin + curnode.ioff; ++i)
			for(j = curnode.jmin; j < curnode.jmin + curnode.joff; ++j) {

				r2d_dest_grid[vind(i, j)] += C*locvol;

				current_info.vtot += locvol;
				if(locvol < current_info.vox_min) current_info.vox_min = locvol;
				if(locvol > current_info.vox_max) current_info.vox_max = locvol;
			
			}
#endif // NO_REDUCTION
			continue;
		}
		if(orcmp != 0x0f) {
			// the leaf lies entirely outside the tet
			// skip it
			continue;
		}
		if(curnode.ioff == 1 && curnode.joff == 1) {
			// we've reached a single cell that straddles the tet boundary 
			// Clip and voxelize it.

			// initialize the unit cube connectivity 
			// copy all flags from tree traversal
			
			// if cell shifting is off, use the absolute vertex coordinates

#ifdef NO_SHIFTING
			// voxel bounds in shifted coordinates
			rbounds[0].x = curnode.imin*r2d_d.x;
			rbounds[0].y = curnode.jmin*r2d_d.y;
			rbounds[1].x = (curnode.imin+1)*r2d_d.x;
			rbounds[1].y = (curnode.jmin+1)*r2d_d.y;
#endif // NO_SHIFTING

			r2du_init_box(vertbuffer, &nverts, rbounds);
			for(v = 0; v < nverts; ++v) {
				vertbuffer[v].fflags = curnode.gridpt[v].fflags;
				for(f = 0; f < 4; ++f)
					vertbuffer[v].fdist[f] = curnode.gridpt[v].fdist[f];
			}

			// Clipping
#ifndef NO_CLIPPING
			r2du_clip_quad(vertbuffer, &nverts, andcmp);
#endif // NO_CLIPPING

			// Reducetion
#ifndef NO_REDUCTION
			r2du_reduce(vertbuffer, &nverts, C, vind(curnode.imin, curnode.jmin));
#endif // NO_REDUCTION
			continue;	
		}

		// else, split the node along its longest dimension
		// and push to the the stack
		i = curnode.ioff/2;
		j = curnode.joff/2;
		if(i >= j) {

			// LEFT NODE
			treestack[ntreestack].imin = curnode.imin;
			treestack[ntreestack].ioff = i;
			treestack[ntreestack].jmin = curnode.jmin;
			treestack[ntreestack].joff = curnode.joff;
			treestack[ntreestack].gridpt[2] = curnode.gridpt[2];
			treestack[ntreestack].gridpt[3] = curnode.gridpt[3];

			// RIGHT NODE
			treestack[ntreestack+1].imin = curnode.imin + i;
			treestack[ntreestack+1].ioff = curnode.ioff - i;
			treestack[ntreestack+1].jmin = curnode.jmin;
			treestack[ntreestack+1].joff = curnode.joff;
			treestack[ntreestack+1].gridpt[0] = curnode.gridpt[0];
			treestack[ntreestack+1].gridpt[1] = curnode.gridpt[1];

			// FILL IN COMMON POINTS
			gpt.x = r2d_d.x*(curnode.imin + i);
			gpt.y = r2d_d.y*curnode.jmin;
			treestack[ntreestack].gridpt[0].fflags = 0x00;
			treestack[ntreestack+1].gridpt[3].fflags = 0x00;
			for(f = 0; f < 4; ++f) {
				gor = faces[f].d + dot(gpt, faces[f].n);
				treestack[ntreestack].gridpt[0].fdist[f] = gor;
				treestack[ntreestack+1].gridpt[3].fdist[f] = gor;
				if(gor > 0.0) {
					treestack[ntreestack].gridpt[0].fflags |= (1 << f);
					treestack[ntreestack+1].gridpt[3].fflags |= (1 << f);
				}
			}
			gpt.y = r2d_d.y*(curnode.jmin + curnode.joff);
			treestack[ntreestack].gridpt[1].fflags = 0x00;
			treestack[ntreestack+1].gridpt[2].fflags = 0x00;
			for(f = 0; f < 4; ++f) {
				gor = faces[f].d + dot(gpt, faces[f].n);
				treestack[ntreestack].gridpt[1].fdist[f] = gor;
				treestack[ntreestack+1].gridpt[2].fdist[f] = gor;
				if(gor > 0.0) {
					treestack[ntreestack].gridpt[1].fflags |= (1 << f);
					treestack[ntreestack+1].gridpt[2].fflags |= (1 << f);
				}
			}
			ntreestack += 2;
			continue;
		}
		else { // j > i

			// LEFT NODE
			treestack[ntreestack].imin = curnode.imin;
			treestack[ntreestack].ioff = curnode.ioff;
			treestack[ntreestack].jmin = curnode.jmin;
			treestack[ntreestack].joff = j;
			treestack[ntreestack].gridpt[0] = curnode.gridpt[0];
			treestack[ntreestack].gridpt[3] = curnode.gridpt[3];

			// RIGHT NODE
			treestack[ntreestack+1].imin = curnode.imin;
			treestack[ntreestack+1].ioff = curnode.ioff;
			treestack[ntreestack+1].jmin = curnode.jmin + j;
			treestack[ntreestack+1].joff = curnode.joff - j;
			treestack[ntreestack+1].gridpt[1] = curnode.gridpt[1];
			treestack[ntreestack+1].gridpt[2] = curnode.gridpt[2];

			// FILL IN COMMON POINTS
			gpt.x = r2d_d.x*curnode.imin;
			gpt.y = r2d_d.y*(curnode.jmin + j);
			treestack[ntreestack].gridpt[2].fflags = 0x00;
			treestack[ntreestack+1].gridpt[3].fflags = 0x00;
			for(f = 0; f < 4; ++f) {
				gor = faces[f].d + dot(gpt, faces[f].n);
				treestack[ntreestack].gridpt[2].fdist[f] = gor;
				treestack[ntreestack+1].gridpt[3].fdist[f] = gor;
				if(gor > 0.0) {
					treestack[ntreestack].gridpt[2].fflags |= (1 << f);
					treestack[ntreestack+1].gridpt[3].fflags |= (1 << f);
				}
			}
			gpt.x = r2d_d.x*(curnode.imin + curnode.ioff);
			treestack[ntreestack].gridpt[1].fflags = 0x00;
			treestack[ntreestack+1].gridpt[0].fflags = 0x00;
			for(f = 0; f < 4; ++f) {
				gor = faces[f].d + dot(gpt, faces[f].n);
				treestack[ntreestack].gridpt[1].fdist[f] = gor;
				treestack[ntreestack+1].gridpt[0].fdist[f] = gor;
				if(gor > 0.0) {
					treestack[ntreestack].gridpt[1].fflags |= (1 << f);
					treestack[ntreestack+1].gridpt[0].fflags |= (1 << f);
				}
			}
			ntreestack += 2;
			continue;
		}
	}

#else // !USE_TREE

	// make sure the grid buffer is large enough
	if(r2d_gblen < (nx+1)*(ny+1)) {
		r2d_gblen = (nx+1)*(ny+1)*2; // extra factor of two to minimize reallocations 
		r2d_gridpoint* ngb = (r2d_gridpoint*) realloc((void*) r2d_gridbuffer, r2d_gblen*sizeof(r2d_gridpoint));
		if(ngb == r2d_gridbuffer) {
			current_info.errmsg = "Bad reallocation of grid buffer";
			current_info.good = 0;
			return current_info;
		}
		r2d_gridbuffer = ngb;//(r2d_gridpoint*) realloc((void*) r2d_gridbuffer, r2d_gblen*sizeof(r2d_gridpoint));
	}

	// check all grid vertices in the patch against each tet face
	for(i = 0; i <= nx; ++i)
	for(j = 0; j <= ny; ++j) {

		gpt.x = i*r2d_d.x; gpt.y = j*r2d_d.y;
		vv0 = gind(i, j);
		r2d_gridbuffer[vv0].fflags = 0x00;

		// flag the vertex for each face it lies inside
		// also save its distance from each face
		for(f = 0; f < 4; ++f) {
			gor = faces[f].d + dot(gpt, faces[f].n);
			if(gor > 0.0) r2d_gridbuffer[vv0].fflags |= (1 << f);
			r2d_gridbuffer[vv0].fdist[f] = gor;
		}
	}

	// iterate over all voxels in the patch
	for(i = 0; i < nx; ++i)
	for(j = 0; j < ny; ++j) {

		// precompute flattened grid indices
		vv[0] = gind(i+1, j);
		vv[1] = gind(i+1, j+1);
		vv[2] = gind(i, j+1);
		vv[3] = gind(i, j);
	
		// check inclusion of each voxel within the tet
		orcmp = 0x00;
        andcmp = 0x0f;
		for(v = 0; v < 4; ++v) {
			orcmp |= r2d_gridbuffer[vv[v]].fflags; 
			andcmp &= r2d_gridbuffer[vv[v]].fflags; 
		}

		if(andcmp == 0x0f) {
			// the voxel is entirely inside the quad

#ifndef NO_REDUCTION
			
			locvol = r2d_d.x*r2d_d.y;

			r2d_dest_grid[vind(i, j)] += C*locvol;

			current_info.vtot += locvol;
			if(locvol < current_info.vox_min) current_info.vox_min = locvol;
			if(locvol > current_info.vox_max) current_info.vox_max = locvol;

#endif // NO_REDUCTION
		}	
		else if(orcmp == 0x0f) {
			// the voxel crosses the boundary of the tet
			// need to process further
			
			// initialize the unit cube connectivity 
			// if shifting is turned off, use absolute coordinates

#ifdef NO_SHIFTING
			// voxel bounds in shifted coordinates
			rbounds[0].x = i*r2d_d.x;
			rbounds[0].y = j*r2d_d.y;
			rbounds[1].x = (i+1)*r2d_d.x;
			rbounds[1].y = (j+1)*r2d_d.y;
#endif // NO_SHIFTING

			r2du_init_box(vertbuffer, &nverts, rbounds);
			for(v = 0; v < nverts; ++v) {
				vertbuffer[v].fflags = r2d_gridbuffer[vv[v]].fflags;
				for(f = 0; f < 4; ++f)
					vertbuffer[v].fdist[f] = r2d_gridbuffer[vv[v]].fdist[f];
			}

			// Clipping
#ifndef NO_CLIPPING
		r2du_clip_quad(vertbuffer, &nverts, andcmp);
		/*r2du_clip_n(vertbuffer, &nverts, faces, 4);*/
#endif // NO_CLIPPING

#ifndef NO_REDUCTION
		r2du_reduce(vertbuffer, &nverts, C, vind(i, j));
#endif // NO_REDUCTION
		}
	}
#endif // USE_TREE

	return current_info;
}



#if 0
///// WARNING! MODIFIES THE INPUT FACES! ///////
r2d_info r2d_rasterize_n(r2d_real C, r2d_plane* faces, r2d_int nfaces, r2d_dvec2 ibounds[2]) {

	// variables used in this function
	r2d_real gor, locvol;
	r2d_rvec2 gpt;
	r2d_int i, j, nx, ny;
	unsigned char v, f;
	unsigned char orcmp, andcmp;

// macros for grid access
#define vind(i, j) (r2d_dest_dims.y*(i + ibounds[0].x) + (j + ibounds[0].y))
#ifdef USE_TREE
	// stack for the tree
	r2d_treenode treestack[256]; // longer
	r2d_int ntreestack;
	r2d_treenode curnode;
#else
#define gind(i, j) ((ny+1)*(i) + (j))
	r2d_long vv[4];
	r2d_long vv0;
#endif

	// voxel bounds in shifted coordinates
	r2d_rvec2 rbounds[2] = {
		{-0.5*r2d_d.x, -0.5*r2d_d.y}, {0.5*r2d_d.x, 0.5*r2d_d.y}
	};
	// vertex buffer
	r2d_int nverts;
	r2d_vertex vertbuffer[128]; // longer

	// start keeping track for this quad
	current_info.good = 1;
	current_info.vtot = 0.0;
	current_info.vox_min = 1.0e99;
	current_info.vox_max = -1.0e99;

	
	// shift the faces to align with
	// the working range of indices
	// THIS CHANGES THE INPUT FACES
	for(f = 0; f < nfaces; ++f) {
		faces[f].d += ibounds[0].x*r2d_d.x*faces[f].n.x
							+ ibounds[0].y*r2d_d.y*faces[f].n.y;
	}
	nx = ibounds[1].x - ibounds[0].x; 
	ny = ibounds[1].y - ibounds[0].y; 

#ifdef USE_TREE
	
	// get the initial face orientations for each corner of the node
	gpt.x = nx*r2d_d.x;
	gpt.y = 0.0;
	curnode.gridpt[0].fflags = 0x00;
	for(f = 0; f < nfaces; ++f) {
		gor = faces[f].d + dot(gpt, faces[f].n);
		if(gor > 0.0) curnode.gridpt[0].fflags |= (1 << f);
		curnode.gridpt[0].fdist[f] = gor;
	}
	gpt.y = ny*r2d_d.y;
	curnode.gridpt[1].fflags = 0x00;
	for(f = 0; f < nfaces; ++f) {
		gor = faces[f].d + dot(gpt, faces[f].n);
		if(gor > 0.0) curnode.gridpt[1].fflags |= (1 << f);
		curnode.gridpt[1].fdist[f] = gor;
	}
	gpt.x = 0.0;
	curnode.gridpt[2].fflags = 0x00;
	for(f = 0; f < nfaces; ++f) {
		gor = faces[f].d + dot(gpt, faces[f].n);
		if(gor > 0.0) curnode.gridpt[2].fflags |= (1 << f);
		curnode.gridpt[2].fdist[f] = gor;
	}
	gpt.y = 0.0;
	curnode.gridpt[3].fflags = 0x00;
	for(f = 0; f < nfaces; ++f) {
		gor = faces[f].d + dot(gpt, faces[f].n);
		if(gor > 0.0) curnode.gridpt[3].fflags |= (1 << f);
		curnode.gridpt[3].fdist[f] = gor;
	}

	curnode.imin = 0;
	curnode.jmin = 0;
	curnode.ioff = nx;
	curnode.joff = ny;

	ntreestack = 0;
	treestack[ntreestack++] = curnode;
	while(ntreestack > 0) {

		// pop the top node
		curnode = treestack[--ntreestack];

		orcmp = 0x00;
        andcmp = 0x0f;
		for(f = 0; f < nfaces; ++f) {
			orcmp |= curnode.gridpt[f].fflags; 
			andcmp &= curnode.gridpt[f].fflags; 
		}

		if(andcmp == 0x0f) { // longer
			
			// all cells in this leaf are fully contained
#ifndef NO_REDUCTION
#ifndef NO_SHIFTING
			locvol = r2d_d.x*r2d_d.y;
#endif // NO_SHIFTING
			for(i = curnode.imin; i < curnode.imin + curnode.ioff; ++i)
			for(j = curnode.jmin; j < curnode.jmin + curnode.joff; ++j) {
#ifdef NO_SHIFTING
				/*
				locvol = 0.0;
				locvol -= (i*r2d_d.x)*(j*r2d_d.y)*(k*r2d_d.z);
				locvol += ((i+1)*r2d_d.x)*(j*r2d_d.y)*(k*r2d_d.z);
				locvol -= ((i+1)*r2d_d.x)*((j+1)*r2d_d.y)*(k*r2d_d.z);
				locvol += (i*r2d_d.x)*((j+1)*r2d_d.y)*(k*r2d_d.z);
				locvol += (i*r2d_d.x)*(j*r2d_d.y)*((k+1)*r2d_d.z);
				locvol -= ((i+1)*r2d_d.x)*(j*r2d_d.y)*((k+1)*r2d_d.z);
				locvol += ((i+1)*r2d_d.x)*((j+1)*r2d_d.y)*((k+1)*r2d_d.z);
				locvol -= (i*r2d_d.x)*((j+1)*r2d_d.y)*((k+1)*r2d_d.z);
				*/
#endif
				r2d_dest_grid[vind(i, j)] += C*locvol;

				current_info.vtot += locvol;
				if(locvol < current_info.vox_min) current_info.vox_min = locvol;
				if(locvol > current_info.vox_max) current_info.vox_max = locvol;
			
			}
#endif // NO_REDUCTION
			continue;
		}
		if(orcmp != 0x0f) { // longer
			// the leaf lies entirely outside the tet
			// skip it
			continue;
		}
		if(curnode.ioff == 1 && curnode.joff == 1) {
			// we've reached a single cell that straddles the tet boundary 
			// Clip and voxelize it.

			// initialize the unit cube connectivity 
			// copy all flags from tree traversal
			
			// if cell shifting is off, use the absolute vertex coordinates
#ifdef NO_SHIFTING
			rbounds[2] = {
				{i*r2d_d.x, j*r2d_d.y}, {(i+1)*r2d_d.x, (j+1)*r2d_d.y}
			};
#endif // NO_SHIFTING
			r2du_init_box(vertbuffer, &nverts, rbounds);
			for(v = 0; v < nverts; ++v) {
				vertbuffer[v].fflags = curnode.gridpt[v].fflags;
				for(f = 0; f < nfaces; ++f)
					vertbuffer[v].fdist[f] = curnode.gridpt[v].fdist[f];
			}

			// Clipping
#ifndef NO_CLIPPING
			r2du_clip_quad(vertbuffer, &nverts, andcmp);
#endif // NO_CLIPPING

			// Reducetion
#ifndef NO_REDUCTION
			r2du_reduce(vertbuffer, &nverts, C, vind(curnode.imin, curnode.jmin));
#endif // NO_REDUCTION
			continue;	
		}

		// else, split the node along its longest dimension
		// and push to the the stack
		i = curnode.ioff/2;
		j = curnode.joff/2;
		if(i >= j) {

			// LEFT NODE
			treestack[ntreestack].imin = curnode.imin;
			treestack[ntreestack].ioff = i;
			treestack[ntreestack].jmin = curnode.jmin;
			treestack[ntreestack].joff = curnode.joff;
			treestack[ntreestack].gridpt[2] = curnode.gridpt[2];
			treestack[ntreestack].gridpt[3] = curnode.gridpt[3];

			// RIGHT NODE
			treestack[ntreestack+1].imin = curnode.imin + i;
			treestack[ntreestack+1].ioff = curnode.ioff - i;
			treestack[ntreestack+1].jmin = curnode.jmin;
			treestack[ntreestack+1].joff = curnode.joff;
			treestack[ntreestack+1].gridpt[0] = curnode.gridpt[0];
			treestack[ntreestack+1].gridpt[1] = curnode.gridpt[1];

			// FILL IN COMMON POINTS
			gpt.x = r2d_d.x*(curnode.imin + i);
			gpt.y = r2d_d.y*curnode.jmin;
			treestack[ntreestack].gridpt[0].fflags = 0x00;
			treestack[ntreestack+1].gridpt[3].fflags = 0x00;
			for(f = 0; f < nfaces; ++f) {
				gor = faces[f].d + dot(gpt, faces[f].n);
				treestack[ntreestack].gridpt[0].fdist[f] = gor;
				treestack[ntreestack+1].gridpt[3].fdist[f] = gor;
				if(gor > 0.0) {
					treestack[ntreestack].gridpt[0].fflags |= (1 << f);
					treestack[ntreestack+1].gridpt[3].fflags |= (1 << f);
				}
			}
			gpt.y = r2d_d.y*(curnode.jmin + curnode.joff);
			treestack[ntreestack].gridpt[1].fflags = 0x00;
			treestack[ntreestack+1].gridpt[2].fflags = 0x00;
			for(f = 0; f < nfaces; ++f) {
				gor = faces[f].d + dot(gpt, faces[f].n);
				treestack[ntreestack].gridpt[1].fdist[f] = gor;
				treestack[ntreestack+1].gridpt[2].fdist[f] = gor;
				if(gor > 0.0) {
					treestack[ntreestack].gridpt[1].fflags |= (1 << f);
					treestack[ntreestack+1].gridpt[2].fflags |= (1 << f);
				}
			}
			ntreestack += 2;
			continue;
		}
		else { // j > i

			// LEFT NODE
			treestack[ntreestack].imin = curnode.imin;
			treestack[ntreestack].ioff = curnode.ioff;
			treestack[ntreestack].jmin = curnode.jmin;
			treestack[ntreestack].joff = j;
			treestack[ntreestack].gridpt[0] = curnode.gridpt[0];
			treestack[ntreestack].gridpt[3] = curnode.gridpt[3];

			// RIGHT NODE
			treestack[ntreestack+1].imin = curnode.imin;
			treestack[ntreestack+1].ioff = curnode.ioff;
			treestack[ntreestack+1].jmin = curnode.jmin + j;
			treestack[ntreestack+1].joff = curnode.joff - j;
			treestack[ntreestack+1].gridpt[1] = curnode.gridpt[1];
			treestack[ntreestack+1].gridpt[2] = curnode.gridpt[2];

			// FILL IN COMMON POINTS
			gpt.x = r2d_d.x*curnode.imin;
			gpt.y = r2d_d.y*(curnode.jmin + j);
			treestack[ntreestack].gridpt[2].fflags = 0x00;
			treestack[ntreestack+1].gridpt[3].fflags = 0x00;
			for(f = 0; f < nfaces; ++f) {
				gor = faces[f].d + dot(gpt, faces[f].n);
				treestack[ntreestack].gridpt[2].fdist[f] = gor;
				treestack[ntreestack+1].gridpt[3].fdist[f] = gor;
				if(gor > 0.0) {
					treestack[ntreestack].gridpt[2].fflags |= (1 << f);
					treestack[ntreestack+1].gridpt[3].fflags |= (1 << f);
				}
			}
			gpt.x = r2d_d.x*(curnode.imin + curnode.ioff);
			treestack[ntreestack].gridpt[1].fflags = 0x00;
			treestack[ntreestack+1].gridpt[0].fflags = 0x00;
			for(f = 0; f < nfaces; ++f) {
				gor = faces[f].d + dot(gpt, faces[f].n);
				treestack[ntreestack].gridpt[1].fdist[f] = gor;
				treestack[ntreestack+1].gridpt[0].fdist[f] = gor;
				if(gor > 0.0) {
					treestack[ntreestack].gridpt[1].fflags |= (1 << f);
					treestack[ntreestack+1].gridpt[0].fflags |= (1 << f);
				}
			}
			ntreestack += 2;
			continue;
		}
	}

#else // !USE_TREE

	// make sure the grid buffer is large enough
	if(r2d_gblen < (nx+1)*(ny+1)) {
		r2d_gblen = (nx+1)*(ny+1)*2; // extra factor of two to minimize reallocations 
		r2d_gridpoint* ngb = (r2d_gridpoint*) realloc((void*) r2d_gridbuffer, r2d_gblen*sizeof(r2d_gridpoint));
		if(ngb == r2d_gridbuffer) {
			current_info.errmsg = "Bad reallocation of grid buffer";
			current_info.good = 0;
			return current_info;
		}
		r2d_gridbuffer = ngb;//(r2d_gridpoint*) realloc((void*) r2d_gridbuffer, r2d_gblen*sizeof(r2d_gridpoint));
	}

	// check all grid vertices in the patch against each tet face
	for(i = 0; i <= nx; ++i)
	for(j = 0; j <= ny; ++j) {

		gpt.x = i*r2d_d.x; gpt.y = j*r2d_d.y;
		vv0 = gind(i, j);
		r2d_gridbuffer[vv0].fflags = 0x00;

		// flag the vertex for each face it lies inside
		// also save its distance from each face
		for(f = 0; f < nfaces; ++f) {
			gor = faces[f].d + dot(gpt, faces[f].n);
			if(gor > 0.0) r2d_gridbuffer[vv0].fflags |= (1 << f);
			r2d_gridbuffer[vv0].fdist[f] = gor;
		}
	}

	// iterate over all voxels in the patch
	for(i = 0; i < nx; ++i)
	for(j = 0; j < ny; ++j) {

		// precompute flattened grid indices
		vv[0] = gind(i+1, j);
		vv[1] = gind(i+1, j+1);
		vv[2] = gind(i, j+1);
		vv[3] = gind(i, j);
	
		// check inclusion of each voxel within the tet
		orcmp = 0x00;
        andcmp = 0x0f; // longer
		for(v = 0; v < 4; ++v) {
			orcmp |= r2d_gridbuffer[vv[v]].fflags; 
			andcmp &= r2d_gridbuffer[vv[v]].fflags; 
		}

		if(andcmp == 0x0f) { // longer
			// the voxel is entirely inside the quad

#ifndef NO_REDUCTION
#ifndef NO_SHIFTING
			locvol = r2d_d.x*r2d_d.y;
#else
			/*
			locvol = 0.0;
			locvol -= (i*r2d_d.x)*(j*r2d_d.y)*(k*r2d_d.z);
			locvol += ((i+1)*r2d_d.x)*(j*r2d_d.y)*(k*r2d_d.z);
			locvol -= ((i+1)*r2d_d.x)*((j+1)*r2d_d.y)*(k*r2d_d.z);
			locvol += (i*r2d_d.x)*((j+1)*r2d_d.y)*(k*r2d_d.z);
			locvol += (i*r2d_d.x)*(j*r2d_d.y)*((k+1)*r2d_d.z);
			locvol -= ((i+1)*r2d_d.x)*(j*r2d_d.y)*((k+1)*r2d_d.z);
			locvol += ((i+1)*r2d_d.x)*((j+1)*r2d_d.y)*((k+1)*r2d_d.z);
			locvol -= (i*r2d_d.x)*((j+1)*r2d_d.y)*((k+1)*r2d_d.z);
			*/
#endif // NO_SHIFTING
			r2d_dest_grid[vind(i, j)] += C*locvol;

			current_info.vtot += locvol;
			if(locvol < current_info.vox_min) current_info.vox_min = locvol;
			if(locvol > current_info.vox_max) current_info.vox_max = locvol;

#endif // NO_REDUCTION
		}	
		else if(orcmp == 0x0f) { // longer
			// the voxel crosses the boundary of the tet
			// need to process further
			
			// initialize the unit cube connectivity 
			// if shifting is turned off, use absolute coordinates
#ifdef NO_SHIFTING
			rbounds = {
				{i*r2d_d.x, j*r2d_d.y}, {(i+1)*r2d_d.x, (j+1)*r2d_d.y}
			};
#endif // NO_SHIFTING
			r2du_init_box(vertbuffer, &nverts, rbounds);
			for(v = 0; v < nverts; ++v) {
				vertbuffer[v].fflags = r2d_gridbuffer[vv[v]].fflags;
				for(f = 0; f < 4; ++f)
					vertbuffer[v].fdist[f] = r2d_gridbuffer[vv[v]].fdist[f];
			}

			// Clipping
#ifndef NO_CLIPPING
		/*r2du_clip_quad(vertbuffer, &nverts, andcmp);*/
		r2du_clip_n(vertbuffer, &nverts, faces, 4);
#endif // NO_CLIPPING

#ifndef NO_REDUCTION
		r2du_reduce(vertbuffer, &nverts, C, vind(i, j));
#endif // NO_REDUCTION
		}
	}
#endif // USE_TREE

	return current_info;
}




inline void r2du_clip_n_andcmp(r2d_vertex64* vertbuffer, r2d_int* nverts, r2d_ulong andcmp, r2d_int nfaces) {

	r2d_int v, f, ff;
	r2d_int vstart, vnext, vcur, searching;
	r2d_ulong fmask, ffmask;

	// loop over faces
	for(f = 0; f < nfaces; ++f) {

		fmask = (1 << f);
		if(andcmp & fmask) continue;

		// find the first vertex lying outside of the face
		// only need to find one (taking advantage of convexity)
		// TODO: add a smarter graph traversal here!
		vcur = 41663;
		for(v = 0; vcur >= 41663 && v < *nverts; ++v) 
			if(!(vertbuffer[v].fflags & (CLIP_MASK_64 | fmask))) vcur = v;
		if(vcur >= 41663) continue;
				
/*
		// traverse to the vertex lying oustide of the clip plane
		while(r2d_vertbuffer[vcur].fflags & fmask) {
			if(r2d_vertbuffer[vcur].fdist)
			if(r2d_vertbuffer[r2d_vertbuffer[vcur].pnbrs[0]].fdist[f] < r2d_vertbuffer[r2d_vertbuffer[vcur].pnbrs[1]].fdist[f]) {
				vcur = r2d_vertbuffer[vcur].pnbrs[0];
			}
			else {
				vcur = r2d_vertbuffer[vcur].pnbrs[1];
			}

		}
*/

		// traverse to the right
		vstart = vcur;
		searching = 1;
		while(searching) { 
			vnext = vertbuffer[vcur].pnbrs[1];
			if(fmask & vertbuffer[vnext].fflags) {
				// vnext is inside the face

				// compute the intersection point using a weighted
				// average of perpendicular distances to the plane
				wav(vertbuffer[vcur].pos, vertbuffer[vnext].fdist[f],
					vertbuffer[vnext].pos, -vertbuffer[vcur].fdist[f],
					vertbuffer[*nverts].pos);
				
				vertbuffer[*nverts].pnbrs[1] = vnext;

				// reciprocal connetivity
				vertbuffer[vnext].pnbrs[0] = *nverts;

				// do face intersections and flags
				vertbuffer[*nverts].fflags = 0x00;
				for(ff = f + 1; ff < nfaces; ++ff) {

					// skip if all initial verts are inside ff
					ffmask = (1 << ff); 
					if(andcmp & ffmask) continue;

					// weighted average keeps us in a relative coordinate system
					
					vertbuffer[*nverts].fdist[ff] = (vertbuffer[vcur].fdist[ff]*vertbuffer[vnext].fdist[f] 
							- vertbuffer[vnext].fdist[ff]*vertbuffer[vcur].fdist[f])
							/(vertbuffer[vnext].fdist[f] - vertbuffer[vcur].fdist[f]);
					if(vertbuffer[*nverts].fdist[ff] > 0.0) vertbuffer[*nverts].fflags |= ffmask;
				}
				++(*nverts);
				searching = 0;
			}
			vertbuffer[vcur].fflags |= CLIP_MASK_64;
			vcur = vnext;
		}

		// traverse to the left
		vcur = vstart;
		searching = 1;
		while(searching) { 
			vnext = vertbuffer[vcur].pnbrs[0];
			if(fmask & vertbuffer[vnext].fflags) {
				// vnext is inside the face

				// compute the intersection point using a weighted
				// average of perpendicular distances to the plane
				wav(vertbuffer[vcur].pos, vertbuffer[vnext].fdist[f],
					vertbuffer[vnext].pos, -vertbuffer[vcur].fdist[f],
					vertbuffer[*nverts].pos);
				
				vertbuffer[*nverts].pnbrs[0] = vnext;

				// reciprocal connetivity
				vertbuffer[vnext].pnbrs[1] = *nverts;

				// do face intersections and flags
				vertbuffer[*nverts].fflags = 0x00;
				for(ff = f + 1; ff < nfaces; ++ff) {

					// skip if all verts are inside ff
					ffmask = (1 << ff); 
					if(andcmp & ffmask) continue;

					// weighted average keeps us in a relative coordinate system
					vertbuffer[*nverts].fdist[ff] = (vertbuffer[vcur].fdist[ff]*vertbuffer[vnext].fdist[f] 
							- vertbuffer[vnext].fdist[ff]*vertbuffer[vcur].fdist[f])
							/(vertbuffer[vnext].fdist[f] - vertbuffer[vcur].fdist[f]);
					if(vertbuffer[*nverts].fdist[ff] > 0.0) vertbuffer[*nverts].fflags |= ffmask;

				}
				++(*nverts);
				searching = 0;
			}
			vertbuffer[vcur].fflags |= CLIP_MASK_64;
			vcur = vnext;
		}
		vertbuffer[*nverts-2].pnbrs[0] = *nverts-1; 
		vertbuffer[*nverts-1].pnbrs[1] = *nverts-2;
	}
}


inline void r2du_clip_n(r2d_vertex* vertbuffer, r2d_int* nverts,
							r2d_plane* faces, r2d_int nfaces) {

	unsigned char v, f;
	unsigned char vstart, vnext, vcur, searching;
	unsigned char sdir, skipface;
	r2d_real fcur, ftmp, fmin;

	// start at the first unclipped vertex
	vcur = 127;
	for(v = 0; vcur >= 127 && v < *nverts; ++v) 
		if(!(vertbuffer[v].fflags & CLIP_MASK)) vcur = v;
	
	// exit if no unclipped vertex exists
	if(vcur >= 127) return;

	// loop over faces
	for(f = 0; f < nfaces; ++f) {
		
		// if the vertex is inside the clip plane,	
		// traverse to a vertex lying outside of the clip plane
		fcur = faces[f].d + dot(vertbuffer[vcur].pos, faces[f].n);
		skipface = 0;
		if(fcur > 0.0) {

			fmin = fcur;
			sdir = 127;

			// search in the direction of steepest descent
			ftmp = faces[f].d + dot(vertbuffer[vertbuffer[vcur].pnbrs[0]].pos, faces[f].n);
			if(ftmp < fmin) {
				sdir = 0;
				fmin = ftmp;
			}
			ftmp = faces[f].d + dot(vertbuffer[vertbuffer[vcur].pnbrs[1]].pos, faces[f].n);
			if(ftmp < fmin) {
				sdir = 1;
				fmin = ftmp;
			}
			
			if(sdir >= 127) {
				// skip the face, all verts lie inside
				continue;
			}

			fcur = fmin;
			vcur = vertbuffer[vcur].pnbrs[sdir];
	
			while(fcur > 0.0) {
	
				// check neighbor
				fmin = faces[f].d + dot(vertbuffer[vertbuffer[vcur].pnbrs[sdir]].pos, faces[f].n);
				if(fmin < fcur) {
					fcur = fmin;
					vcur = vertbuffer[vcur].pnbrs[sdir];
				}
				else {
					// vcur is the closest to the face, but still not inside
					// skip the face
					skipface = 1;
					break;
				}
			}
		}
		// continue if no verts lie outside of the face
		if(skipface) continue;

		// vcur is now guaranteed to lie outside of the clip plane
		//
		// TODO: what if all verts need to be clipped?? Check left == right?

		// traverse to the right
		vstart = vcur;
		searching = 1;
		while(vcur != vstart && searching) { 
			vnext = vertbuffer[vcur].pnbrs[1];
			if(faces[f].d + dot(vertbuffer[vnext].pos, faces[f].n) > 0.0) {
				// vnext is inside the face

				// compute the intersection point using a weighted
				// average of perpendicular distances to the plane
				wav(vertbuffer[vcur].pos, vertbuffer[vnext].fdist[f],
					vertbuffer[vnext].pos, -vertbuffer[vcur].fdist[f],
					vertbuffer[*nverts].pos);
				
				vertbuffer[*nverts].pnbrs[1] = vnext;

				// reciprocal connetivity
				vertbuffer[vnext].pnbrs[0] = *nverts;

				// add flags in case this vert is clipped later on
				vertbuffer[*nverts].fflags = 0x00;

				++(*nverts);
				searching = 0;
			}
			vertbuffer[vcur].fflags |= CLIP_MASK;
			vcur = vnext;
		}
		if(vcur == vstart) {
			// All verts are behind the clip plane!
			for(v = 0; vcur >= 127 && v < *nverts; ++v) 
				vertbuffer[v].fflags |= CLIP_MASK;
			return;
		}

		// traverse to the left
		vcur = vstart;
		searching = 1;
		while(searching) { 
			vnext = vertbuffer[vcur].pnbrs[0];
			if(faces[f].d + dot(vertbuffer[vnext].pos, faces[f].n) > 0.0) {
				// vnext is inside the face

				// compute the intersection point using a weighted
				// average of perpendicular distances to the plane
				wav(vertbuffer[vcur].pos, vertbuffer[vnext].fdist[f],
					vertbuffer[vnext].pos, -vertbuffer[vcur].fdist[f],
					vertbuffer[*nverts].pos);
				
				vertbuffer[*nverts].pnbrs[0] = vnext;

				// reciprocal connetivity
				vertbuffer[vnext].pnbrs[1] = *nverts;

				// do face intersections and flags
				vertbuffer[*nverts].fflags = 0x00;

				++(*nverts);
				searching = 0;
			}
			vertbuffer[vcur].fflags |= CLIP_MASK;
			vcur = vnext;
		}
		vertbuffer[*nverts-2].pnbrs[0] = *nverts-1; 
		vertbuffer[*nverts-1].pnbrs[1] = *nverts-2;
	}
}
#endif


