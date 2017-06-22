#!/u/ki/mbaumer/anaconda/bin/python

import sys, os
sys.path.append('/u/ki/mbaumer/random_pixel_size/weak_sauce/code')
import pandas as pd
import copy
import numpy as np
from weak_sauce.grid import MoveableGrid
from weak_sauce.sources import Source
from weak_sauce.movers import FixedIlluminationMover, SimpleVerticesMover, UniformIlluminationMover
import galsim
import time

def batchScienceAna(*args):
   
    # e, fwhm, theta, flux 
    params = [float(args[0][1]), float(args[0][2]), float(args[0][3]), float(args[0][4])]

    #for LSST (50x50)
    #saved_mg = MoveableGrid('../data/test100k_iter.pkl')
    #flat=np.load('7th_order_LSST_50x50.npy')
    
    #for LSST (full amp)
    saved_mg = MoveableGrid('../data/best_lsst_amp3_mg.pkl')
    flat = np.load('../data/big_lsst_flat.npy')

    #for LSST (10x worse PRNU)
    #saved_mg = MoveableGrid('../data/lsst_model_10x_worsePRNU.pkl')
    #flat = np.load('../data/big_lsst_flat.npy')
    #flat = 10*flat-9

    #for DES tree rings
    #saved_mg = MoveableGrid('/nfs/slac/g/ki/ki19/des/mbaumer/ccd_mg_model_fits/des_chip04_maxit2500_step0.5_decay0.0001/mg.pkl')
    #flat = saved_mg.source.fluxes #TODO: is this a bug? doesn't seem like it would matter, since the model gets so close...

    cutout_side_length = 64 #should be even
    
    res_df = pd.DataFrame()
    
    psf = galsim.Moffat(4.765,fwhm=params[3])
    
    #obj = galsim.Gaussian(fwhm=params[1],flux=1)
    #obj = galsim.Sersic(4,half_light_radius=params[1]) #deV profile        
    #obj = galsim.Convolve(obj,psf)
    obj = psf
    obj = obj.shift(cutout_side_length/2,cutout_side_length/2)
    #obj = obj.shear(e=params[0],beta=galsim.Angle(params[2],galsim.degrees))
    
    #oversampled galaxy to plop down
    stationary_source = Source(num_x=5*(cutout_side_length)+1, flux_func=obj)
    illuminator = FixedIlluminationMover(stationary_source)
    areaFinder = UniformIlluminationMover()
    
    fake_img = np.zeros(flat.shape)

    ideal_grid = Source(num_x=cutout_side_length+1)
    ideal_mg = MoveableGrid(ideal_grid,illuminator)
    ideal_mg.step()
    #ideal_mg.source.plot_pixel_grid()
    ideal_res = ideal_mg.evaluate_psf()
    #ideal_sex_res = ideal_mg.evaluate_sex()['FLUX_AUTO']
    #if len(ideal_sex_res == 0): ideal_res['FLUX_AUTO'] = ideal_sex_res[0]
    #else: ideal_res['FLUX_AUTO'] = 0
    #print ideal_res['FLUX_AUTO']
    #ideal_mg.plot_pixel_grid()
    
    #for xctr in np.floor(np.random.uniform(cutout_side_length/2+1,flat.shape[0]-cutout_side_length/2-1,5)): #num=25 normally
    #    for yctr in np.floor(np.random.uniform(cutout_side_length/2+1,flat.shape[1]-cutout_side_length/2-1,5)): #num=170
    for i in [1]: #for back compatibility
	for j in np.arange(10000):
	    xctr = np.floor(np.random.uniform(cutout_side_length/2+1,flat.shape[0]-cutout_side_length/2-1,1))[0]
            yctr = np.floor(np.random.uniform(cutout_side_length/2+1,flat.shape[1]-cutout_side_length/2-1,1))[0]
            temp = Source(num_x=cutout_side_length+2)
            #print temp.vertices
            print xctr, yctr
            temp.vertices = saved_mg.source.vertices[xctr-cutout_side_length/2:xctr+cutout_side_length/2+2,\
                                                     yctr-cutout_side_length/2:yctr+cutout_side_length/2+2,:].copy()
            #print temp.vertices
            cutoutMover = SimpleVerticesMover(mu_x=(cutout_side_length/2-xctr),mu_y=(cutout_side_length/2-yctr))
            shiftingGrid = MoveableGrid(temp,cutoutMover)
	    shiftingGrid.step()
	    temp.update_centroids()
            areaFindingGrid = MoveableGrid(temp,areaFinder)
            areaFindingGrid.step()
            #temp.fluxes = saved_mg.source.fluxes[yctr-cutout_side_length/2:yctr+cutout_side_length/2, \
            #                                     xctr-cutout_side_length/2:xctr+cutout_side_length/2]
            temp.fluxes = np.zeros_like(temp.fluxes) #flush the CCD!
            #print temp.fluxes.shape
          
            #temp.plot_pixel_grid()
            mg = MoveableGrid(temp, illuminator)
            mg.step()
            mg.source.fluxes /= flat[xctr-cutout_side_length/2:xctr+cutout_side_length/2+1,\
                                     yctr-cutout_side_length/2:yctr+cutout_side_length/2+1]
            #mg.plot_pixel_grid()
            fitted_res = mg.evaluate_psf()
            #sex_res = mg.evaluate_sex()['FLUX_AUTO']
            #if len(sex_res == 0): fitted_res['FLUX_AUTO'] = sex_res[0]
            #else: ideal_res['FLUX_AUTO'] = 0
            #print fitted_res

	    #fake_img[xctr-cutout_side_length/2:xctr+cutout_side_length/2+1,\
            #              yctr-cutout_side_length/2:yctr+cutout_side_length/2+1] += mg.source.fluxes
	    residual = fitted_res-ideal_res
            residual['xctr'] = xctr
            residual['yctr'] = yctr
            res_df = res_df.append(residual)
    res_df["inputE"] = params[0]
    res_df["inputS"] = params[1]
    res_df["inputTheta"] = params[2]
    res_df["inputPSF"] = params[3]
    name = 'lsst_10k_' + str(time.time())
    #res_df.to_pickle('./'+name+'.pkl')
    res_df.to_pickle(os.path.expandvars('$LSST_DATA/scienceImpactDFs/for_paper/'+name+'.pkl'))
    #res_df.to_pickle(os.path.expandvars('$LSST_DATA/scienceImpactDFs/lsst_10xworse/deV_e'+args[0][1]+'_s'+args[0][2]+'_t'+args[0][3]+'_f'+args[0][4]+'.pkl'))
    #np.save('fake_img',fake_img)
    return

if __name__ == "__main__":
    batchScienceAna(sys.argv)
