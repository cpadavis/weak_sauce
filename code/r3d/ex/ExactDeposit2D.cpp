
extern "C" {
#include "r2d.h"
}

// ..... all this stuff goes inside a function or whatever .......

		// things that I assume you initialize on your own
		r2d_int Nx, Ny; // image size in pixels
		r2d_rvec2 image_size; // dimensions of physical image space (e.g. (0.0 to 1.0)^2 )
		r2d_real* dest_grid; // must be Nx*Ny big

		// initialize all the r2d stuff
		r2d_info info;
		info = r2d_init(32*32*32);

		// set the destination grid parameters
		r2d_dvec2 dest_dims = {Nx, Ny};
		r2d_rvec2 dest_window = {image_size.x, image_size.y};
		info = r2d_set_dest_grid(dest_grid, dest_dims, dest_window);

		// temporary storage for each quadrilateral to be rasterized
		r2d_rvec2 verts[4];
		r2d_plane faces[4];

		// for passing working grid indices
		r2d_dvec2 ibounds[2];

		for(quad in all_quads) {

				// this function lives elsewhere in my code
				// but is pretty easy to write on your own :)
				Box bbox = getBoundingBox(quad);

				// copy vertices of the current quad
				for(int v = 0; v < 4; ++v) {
					verts[v].x = quad.verts[v].x;
					verts[v].y = quad.verts[v].y;
				}

				// calculate the area (by splitting into two triangles)
				quadvol = 0.0;
				quadvol -= 0.5*((verts[3].x - verts[0].x)*(verts[1].y - verts[0].y) - (verts[1].x - verts[0].x)*(verts[3].y - verts[0].y)); 
				quadvol += 0.5*((verts[3].x - verts[2].x)*(verts[1].y - verts[2].y) - (verts[1].x - verts[2].x)*(verts[3].y - verts[2].y)); 

				// get the grid indices corresponding
				// to the minimal bounding box
				r2du_get_box_index_bounds(&ibounds[0], bbox.xmin, bbox.ymin, bbox.xmax, bbox.ymax);


				// get the quad in a face description
				// THE QUAD HAD DAMN WELL BETTER BE CONVEX
				r2du_faces_from_verts(&verts[0], 4, &faces[0]);

				// figure out the flux density for this quad
				r2d_real rhoav = quad.flux/quadvol;

				// conservatively rasterize the quad
				info = r2d_rasterize_quad(rhoav, faces, ibounds);

		}

		// free the internal buffer
		r2d_finalize();
		
// ..... done!
