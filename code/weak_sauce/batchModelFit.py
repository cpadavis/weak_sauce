#!/u/ki/mbaumer/anaconda/bin/python
from __future__ import division
import sys, os
sys.path.append('/u/ki/mbaumer/random_pixel_size/weak_sauce/code')
import weak_sauce as ws
import weak_sauce.data_tools
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

from weak_sauce.movers import UniformGaussianMover
from weak_sauce.grid import MoveableGrid
from weak_sauce.sources import Source
from weak_sauce.movers import UniformIlluminationMover, FixedIlluminationMover
from weak_sauce.fit_flat import FlatFitter

def fitModel(*args):
	fitType = 'LSST'
	maxit_arg = int(args[0][1])
	step_arg = float(args[0][2])
	decay_arg = float(args[0][3])
	perturb_arg = float(args[0][4])
	
	if fitType == 'DES':
		foldername = '/nfs/slac/g/ki/ki19/des/mbaumer/ccd_mg_model_fits/des_chip04_maxit'+str(maxit_arg)+'_step'+str(step_arg)+'_decay'+str(decay_arg)
	elif fitType == 'LSST':
		foldername = '/nfs/slac/g/ki/ki19/lsst/mbaumer/ccd_mg_model_fits/lsst_amp03_maxit'+str(maxit_arg)+'_step'+str(step_arg)+'_decay'+str(decay_arg)
	elif fitType == 'validation':
		foldername = '/nfs/slac/g/ki/ki19/lsst/mbaumer/ccd_mg_model_fits/validation_gaussian_synthetic4_lsst_amp03_maxit'+str(maxit_arg)+'_step'+str(step_arg)+'_decay'+str(decay_arg)+'_perturb'+str(perturb_arg)
	else:
		raise IOError('select valid fitType (DES,LSST,validation) to find input data')
	if os.path.exists(foldername):
		raise IOError('folder for this config already exists! delete if you want to do again!')
	
	if fitType == 'DES':
		full_amp_img = fits.getdata('/nfs/slac/g/ki/ki19/des/mbaumer/DES_flatcor_supercal/coadds/coadd_r_04.fits')
		#full_amp_img = fits.getdata('/u/ki/mbaumer/random_pixel_size/weak_sauce/data/coadd_r_04.fits')
		full_amp_img = full_amp_img[100:-100,200:924]
	elif fitType == 'LSST':
		full_amp_img = np.load('/u/ki/mbaumer/random_pixel_size/weak_sauce/data/lsst_ultraflat_75ke_amp3.npy')
		full_amp_img = full_amp_img[100:-100,100:-100]
	elif fitType == 'validation':
		#for perturbation validation
		#full_amp_img = np.load('/u/ki/mbaumer/random_pixel_size/weak_sauce/data/lsst_ultraflat_75ke_amp3.npy')
		#for LSST match validation
		#data_mg = MoveableGrid('/nfs/slac/g/ki/ki19/lsst/mbaumer/ccd_mg_model_fits/lsst_amp03_maxit100000_step0.1_decay0.0/best_mg.pkl') 
		#full_amp_img = data_mg.source.fluxes
		#for gaussian model validation
		input_mg = MoveableGrid('./gaussian_test4') 
		full_amp_img = input_mg.source.fluxes
	else:
		raise IOError('select valid fitType (DES,LSST,validation) to find input data')
	
	fitted = ws.data_tools.fitIlluminationVariation(full_amp_img)
	data_rel_flux_map = (full_amp_img-fitted)/fitted+1
	sigma = np.std(data_rel_flux_map)
	small_img = data_rel_flux_map
	sig = 5
	mask = np.logical_not(np.logical_and(small_img < np.mean(small_img)+sig*sigma,small_img > np.mean(small_img)-sig*sigma))
	locs = np.where(mask)
	small_img[locs] = 1
	small_img = small_img + (1 - np.mean(small_img))
	print 'found ' + str(len(locs[0])) + ' bad pixels.'
    
    #unperturbed
	data_like_source = Source(num_x=small_img.shape[1]+1,num_y=small_img.shape[0]+1) 
	data_like_source.fluxes += 1 #fit to flat field
    
    #perturbed
	#mover = UniformGaussianMover(mu_x=0, mu_y=0, sigma_xx=perturb_arg, sigma_yy=perturb_arg, sigma_xy=0.0)
	#data_like_source = Source(num_x=small_img.shape[1]+1,num_y=small_img.shape[0]+1) 
	#data_like_source.fluxes+=1
	#perturbed_mg = MoveableGrid(data_like_source,mover)
	#perturbed_mg.step() #perturb the initial conditions
    
	small_img = small_img.transpose()
	data_mg = FlatFitter(data_like_source,small_img)
	data_mg.fit(maxiter=maxit_arg,verbose=False, step_size=step_arg,learning_rate_decay=decay_arg)
	
	finalLL = data_mg.lnlike()
	std_resid = np.std(small_img-data_mg.source.fluxes)
	std_small_img = np.std(small_img)
	os.makedirs(foldername)
	os.makedirs(foldername+'/finalLL'+str(finalLL))
	os.makedirs(foldername+'/residFrac'+str(std_resid/std_small_img))	
	data_mg.saveto(foldername+'/mg.pkl')

if __name__ == "__main__":
	fitModel(sys.argv)

