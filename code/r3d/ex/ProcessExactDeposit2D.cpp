/*************************************************************
 *
 *		ProcessExactDeposit2D.cpp
 *		Processing plugin for exact mass deposit using raster2d
 *
 *		This file is part of PSI - the Phase Sheet Intersector.
 *
 *		Devon Powell
 *		5 January 2015
 *
 *************************************************************/

#include "Input.h"
#include "Process.h"
#include "Common.h"

extern "C" {
#include "r2d.h"
}

#define __PLUGIN__ "ExactDeposit2D"

using namespace std;

class ProcessExactDeposit2D : public Process {

private:
	bool zoom;
	bool dorho, donstr, domom;

public:

	ProcessExactDeposit2D(Config* config) {

		dorho = false;
		donstr = false;
		domom = false;
		zoom = false;
		Nx = 0; Ny = 0; Nz = 0;

		// parse options passed from Config
		vector<vector<string>>& optionslist = config->getOptionsList(__PLUGIN__);
		vector<vector<string>>::iterator it;
		for(it = optionslist.begin(); it != optionslist.end(); ++it) {
			vector<string>& thisoption = (*it);
			vector<string>::iterator ittok = thisoption.begin();
			string optname = *ittok;
			int nargs = thisoption.size() - 1;

			// actual parsing of option string tokens
			if(optname.compare("Fields") == 0) {
				myError(nargs >= 1, "\'Fields\' requires at least 1 argument.", __FILE__, __LINE__);
				++ittok;
				for(; ittok != thisoption.end(); ++ittok) {
					string field = *ittok;
					if(field.compare("rho") == 0) dorho = true;
					if(field.compare("nstr") == 0) donstr = true;
					if(field.compare("mom") == 0) domom = true;
				}
			}
			if(optname.compare("Grid") == 0) {
				myError(nargs >= 3, "\'Grid\' requires 3 arguments.", __FILE__, __LINE__);
				++ittok;
				try { Nx = stoi(*ittok, NULL); }
				catch(exception e) {
					myError(false, "\'Grid\' requires integer arguments.", __FILE__, __LINE__);
				}
				++ittok;
				try { Ny = stoi(*ittok, NULL); }
				catch(exception e) {
					myError(false, "\'Grid\' requires integer arguments.", __FILE__, __LINE__);
				}
				++ittok;
				try { Nz = stoi(*ittok, NULL); }
				catch(exception e) {
					myError(false, "\'Grid\' requires integer arguments.", __FILE__, __LINE__);
				}
			}
			if(optname.compare("Window") == 0) {
				myError(nargs >= 6, "\'Window\' requires 6 arguments.", __FILE__, __LINE__);
				++ittok;
				try { window.xmin = stod(*ittok, NULL);	}
				catch(exception e) {
					myError(false, "\'Window\' requires numerical arguments.", __FILE__, __LINE__);
				}
				++ittok;
				try { window.ymin = stod(*ittok, NULL); }
				catch(exception e) {
					myError(false, "\'Window\' requires numerical arguments.", __FILE__, __LINE__);
				}
				++ittok;
				try { window.zmin = stod(*ittok, NULL); }
				catch(exception e) {
					myError(false, "\'Window\' requires numerical arguments.", __FILE__, __LINE__);
				}
				++ittok;
				try { window.xmax = stod(*ittok, NULL); }
				catch(exception e) {
					myError(false, "\'Window\' requires numerical arguments.", __FILE__, __LINE__);
				}
				++ittok;
				try { window.ymax = stod(*ittok, NULL); }
				catch(exception e) {
					myError(false, "\'Window\' requires numerical arguments.", __FILE__, __LINE__);
				}
				++ittok;
				try { window.zmax = stod(*ittok, NULL); }
				catch(exception e) {
					myError(false, "\'Window\' requires numerical arguments.", __FILE__, __LINE__);
				}
				zoom = true;
			}
		}

		if(dorho) fields3d["RHO"] = vector<psi_real>(Nx*Ny*Nz);
		if(donstr) fields3d["NSTR"] = vector<psi_real>(Nx*Ny*Nz);
		if(domom) {
			fields3d["PX"] = vector<psi_real>(Nx*Ny*Nz);
			fields3d["PY"] = vector<psi_real>(Nx*Ny*Nz);
			fields3d["PZ"] = vector<psi_real>(Nx*Ny*Nz);
		}

		printf("Grid: %d %d %d\n", Nx, Ny, Nz);
		if(zoom) printf("Window: %f %f %f %f %f %f\n", window.xmin, window.ymin, window.zmin, window.xmax, window.ymax, window.zmax);

		// make sure we have all required parameters
		myError(dorho || donstr || domom, "Must specify at least one Field parameter.", __FILE__, __LINE__);
		myError(Nx > 0 && Ny > 0 && Nz > 0, "Must specify grid dimensions > 0.", __FILE__, __LINE__);

	}

	void process(Input* input) {

		setbuf(stdout, NULL);

		// figure out basic grid parameters
		if(!zoom) {
			window.xmin = 0.0;
			window.ymin = 0.0;
			window.zmin = 0.0;
			window.xmax = input->getBoxSize();
			window.ymax = input->getBoxSize();
			window.zmax = input->getBoxSize();
		}
		psi_real dx = (window.xmax - window.xmin)/Nx;
		psi_real dy = (window.ymax - window.ymin)/Ny;
		psi_real dz = (window.zmax - window.zmin)/Nz;
		psi_real invdx = 1.0/dx;
		psi_real invdy = 1.0/dy;
		psi_real invdz = 1.0/dz;
		psi_real dv = dx*dy*dz;
		psi_real invdv = 1.0/dv;

		// get references to grids by name
		vector<psi_real>& rho = fields3d["RHO"];
		vector<psi_real>& nstreams = fields3d["NSTR"];
		vector<psi_real>& px = fields3d["PX"];
		vector<psi_real>& py = fields3d["PY"];
		vector<psi_real>& pz = fields3d["PZ"];
		if(!dorho) fields3d.erase("RHO");
		if(!donstr) fields3d.erase("NSTR");
		if(!domom) {
			fields3d.erase("PX");
			fields3d.erase("PY");
			fields3d.erase("PZ");
		}

		r2d_info info;
		info = r2d_init(32*32*32);

//		printf("INFO, r2d_init: good = %d, vtot = %f, vox_min = %f, vox_max = %f\n", info.good, info.vtot, info.vox_min, info.vox_max);

		r2d_dvec2 dest_dims = {Nx, Ny};
		r2d_rvec2 dest_window = {(r2d_real) (window.xmax - window.xmin), 
			(r2d_real) (window.ymax - window.ymin)};
		info = r2d_set_dest_grid(rho.data(), dest_dims, dest_window);


//		printf("INFO, r2d_set_dest_grid: good = %d, vtot = %f, vox_min = %f, vox_max = %f\n", info.good, info.vtot, info.vox_min, info.vox_max);

		r2d_rvec2 verts[4];
		r2d_plane faces[4];
		r2d_dvec2 ibounds[2];

		// Loop through all tets and deposit them to the grid
		psi_tet tet;
		psi_box bbox;
		r2d_real tetvol;
		int64_t i, tot;
		vector<psi_tet>& tetbuffer = input->getTetBuffer();
		while(input->fillTetBuffer(i, tot)) {
			for(vector<psi_tet>::iterator it = tetbuffer.begin(); it != tetbuffer.end(); ++it) {

				tet = *it;

#if 0

				// Transform the mesh
				double r0 = 0.333333333333;
				double sfac = 0.1;
				double rfac = 0.00;
				double xp, xpp;	
                double yp, ypp;
                double r, s;
				for(int v = 0; v < 4; ++v) {

#define S(x) (1+sfac*16*(x)*(x)*((x)-1)*((x)-1))
#define S(x) (1+sfac*2*(x)*((x)-1)*((x)-0.5))

					xp = tet.verts[v].x-0.5;
					yp = tet.verts[v].y-0.5;
					r = sqrt(xp*xp+yp*yp);
					if(r < r0) {


						s = S(r/r0);
						//if(r/r0 < 0.25) s = 4*r/r0;
						//else if(r/r0 < 0.75) s = -4*(r/r0-0.5);
						//else if(r/r0 < 1.0) s = 4*(r/r0-1.0);
						//s *= sfac;
						//s += 1.0;



						if(r/r0 < 0.4) printf("r/r0 = %f, S(r/r0) = %f\n", r/r0, s);

						// scale
						xp *= s;
						yp *= s;

						// rotate
						xpp = xp*cos(s*rfac) - yp*sin(s*rfac);
						ypp = xp*sin(s*rfac) + yp*cos(s*rfac);


						// add box center back in
						tet.verts[v].x = xpp + 0.5;
						tet.verts[v].y = ypp + 0.5;
					}
				}
#endif


				bbox = getBoundingBox(tet);

				//printf("Got tet:\n");
				//for(int ii = 0; ii < 4; ++ii) 
					//printf("  ( %f %f )\n", tet.verts[ii].x, tet.verts[ii].y);
					//
				for(int v = 0; v < 4; ++v) {
					verts[v].x = tet.verts[v].x;
					verts[v].y = tet.verts[v].y;
				}

				tetvol = 0.0;
				tetvol -= 0.5*((verts[3].x - verts[0].x)*(verts[1].y - verts[0].y) - (verts[1].x - verts[0].x)*(verts[3].y - verts[0].y)); 
				tetvol += 0.5*((verts[3].x - verts[2].x)*(verts[1].y - verts[2].y) - (verts[1].x - verts[2].x)*(verts[3].y - verts[2].y)); 

				// get the grid indices corresponding
				// to the minimal bounding box
				r2du_get_box_index_bounds(&ibounds[0], bbox.xmin - window.xmin,
													bbox.ymin - window.ymin,
													bbox.xmax - window.xmin,
													bbox.ymax - window.ymin);


				// get the quad in a face description
				r2du_faces_from_verts(&verts[0], 4, &faces[0]);


				psi_real rhoav = (1.0/1024)/tetvol;//tet.rhos[0];
				info = r2d_rasterize_quad(rhoav*invdv, faces, ibounds);

				//printf("Deposited quad:\n");
				//for(int ii = 0; ii < 4; ++ii) {
					//printf("  ( %f %f )\n", verts[ii].x, verts[ii].y);
				//}


				//printf("INFO, r2d_voxelize_constant: good = %d, vtot = %f, vox_min = %f, vox_max = %f\n", info.good, info.vtot, info.vox_min, info.vox_max);
				//printf("%d %.20e %.20e %.20e %.20e %.20e\n", info.good, tetvol, info.vtot, (tetvol - info.vtot)/tetvol, info.vox_min/(dx*dy*dz), info.vox_max/(dz*dy*dz));


#if 0
				// before reducing the patch to the main grid, check the index bounds
				int pimin = imin, pjmin = jmin, pkmin = kmin;
				if(imin < 0) imin = 0; if(imax >= Nx) imax = Nx;
				if(jmin < 0) jmin = 0; if(jmax >= Ny) jmax = Ny;
				if(kmin < 0) kmin = 0; if(kmax >= Nz) kmax = Nz; 
				
				// ---------------- Density accumulation ----------------------------------
				if(dorho)	
					for(int i = imin; i < imax; ++i)
						for(int j = jmin; j < jmax; ++j)
							for(int k = kmin; k < kmax; ++k)
								rho[Ny*Nz*i + Nz*j + k] += rhoav*invdv*i_c[nj*nk*(i - pimin) + nk*(j - pjmin) + (k - pkmin)];
			
				// ---------------- Momentum density accumulation ----------------------------------
				if(domom) {
					
					BCTransform bct = getBCTransform(tet);
			
					vec3 gradvx;
					psi_real intvx;
					psi_real vertvx[4] = {tet.vels[0].x, tet.vels[1].x, tet.vels[2].x, tet.vels[3].x};
					getLinearField(tet, bct, &vertvx[0], gradvx, intvx);
				
					vec3 gradvy;
					psi_real intvy;
					psi_real vertvy[4] = {tet.vels[0].y, tet.vels[1].y, tet.vels[2].y, tet.vels[3].y};
					getLinearField(tet, bct, &vertvy[0], gradvy, intvy);
				
					vec3 gradvz;
					psi_real intvz;
					psi_real vertvz[4] = {tet.vels[0].z, tet.vels[1].z, tet.vels[2].z, tet.vels[3].z};
					getLinearField(tet, bct, &vertvz[0], gradvz, intvz);
				
						for(int i = imin; i < imax; ++i) 
							for(int j = jmin; j < jmax; ++j) 
								for(int k = kmin; k < kmax; ++k) {
									
									// quadratic mass-weighted velocity
									px[Ny*Nz*i + Nz*j + k] += rhoav*intvx	*invdv*i_c[nj*nk*(i - pimin) + nk*(j - pjmin) + (k - pkmin)];
									px[Ny*Nz*i + Nz*j + k] += rhoav*gradvx.x*invdv*i_x[nj*nk*(i - pimin) + nk*(j - pjmin) + (k - pkmin)];
									px[Ny*Nz*i + Nz*j + k] += rhoav*gradvx.y*invdv*i_y[nj*nk*(i - pimin) + nk*(j - pjmin) + (k - pkmin)];
									px[Ny*Nz*i + Nz*j + k] += rhoav*gradvx.z*invdv*i_z[nj*nk*(i - pimin) + nk*(j - pjmin) + (k - pkmin)];
									py[Ny*Nz*i + Nz*j + k] += rhoav*intvy	*invdv*i_c[nj*nk*(i - pimin) + nk*(j - pjmin) + (k - pkmin)];
									py[Ny*Nz*i + Nz*j + k] += rhoav*gradvy.x*invdv*i_x[nj*nk*(i - pimin) + nk*(j - pjmin) + (k - pkmin)];
									py[Ny*Nz*i + Nz*j + k] += rhoav*gradvy.y*invdv*i_y[nj*nk*(i - pimin) + nk*(j - pjmin) + (k - pkmin)];
									py[Ny*Nz*i + Nz*j + k] += rhoav*gradvy.z*invdv*i_z[nj*nk*(i - pimin) + nk*(j - pjmin) + (k - pkmin)];
									pz[Ny*Nz*i + Nz*j + k] += rhoav*intvz	*invdv*i_c[nj*nk*(i - pimin) + nk*(j - pjmin) + (k - pkmin)];
									pz[Ny*Nz*i + Nz*j + k] += rhoav*gradvz.x*invdv*i_x[nj*nk*(i - pimin) + nk*(j - pjmin) + (k - pkmin)];
									pz[Ny*Nz*i + Nz*j + k] += rhoav*gradvz.y*invdv*i_y[nj*nk*(i - pimin) + nk*(j - pjmin) + (k - pkmin)];
									pz[Ny*Nz*i + Nz*j + k] += rhoav*gradvz.z*invdv*i_z[nj*nk*(i - pimin) + nk*(j - pjmin) + (k - pkmin)];
								}
				}
	
				//------------------- Number of streams ----------------------------------------------------
				if(donstr)	
					for(int i = imin; i < imax; ++i) 
						for(int j = jmin; j < jmax; ++j) 
							for(int k = kmin; k < kmax; ++k) 
								nstreams[Ny*Nz*i + Nz*j + k] += invdv*i_c[nj*nk*(i - pimin) + nk*(j - pjmin) + (k - pkmin)];
	
				if(patch) free(patch);
	
				/*
				printf("verts: \n");
				for(int v = 0; v < 4; ++v) {
					printf(" ( %f %f %f )", tet.verts[v].x, tet.verts[v].y, tet.verts[v].z);
				}
				printf("\n");
				printf("vels: \n");
				for(int v = 0; v < 4; ++v) {
					printf(" ( %f %f %f )", tet.vels[v].x, tet.vels[v].y, tet.vels[v].z);
				}
				printf("\n");
				*/
#endif
			}
			printf("Tet %ld of %ld, %.02f%%\r", i, tot, 100.0*i/tot);
		}
		printf("\n");
		r2d_finalize();
		
	}

	string getPluginName() { 
		return __PLUGIN__;
	}

};

// register this plugin with the factory
namespace {
	ProcessFactoryDerived<ProcessExactDeposit2D> factory(__PLUGIN__);
}

#undef __PLUGIN__
