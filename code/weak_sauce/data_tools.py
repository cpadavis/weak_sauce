"""
data_tools.py: some methods for working with laboratory data in weak_sauce
"""
from __future__ import division
import numpy as np
import scipy as sp
from scipy import ndimage, stats
import matplotlib.pyplot as plt
import weak_sauce
from weak_sauce.movers import UniformIlluminationMover
from weak_sauce.shifted_cmap import shiftedColorMap
from weak_sauce.movers import UniformIlluminationMover
from weak_sauce.sources import Source

def mycrop(image,border):
    return image[border:image.shape[0]-border,border:image.shape[1]-border]

def detrend(image):
    import warnings
    warnings.warn('use fitIlluminationVariation(); this is an outdated smoothing method', Warning)

    kernel = (1/273)*np.array([[1,4,7,4,1],
                [4,16,26,16,4],
                [7,26,41,26,7],
                [4, 16,26,16,4],
                [1,4,7,4,1]])
    smoothed = sp.ndimage.convolve(image,kernel)
    return image-smoothed

def fitIlluminationVariation(img, order=7):
    """
        returns 7th order polynomial fit ("predicted")
        detrended = img - predicted
    """
    from sklearn.linear_model import LinearRegression
    x = np.arange(img.shape[1])
    y = np.arange(img.shape[0])
    xx, yy = np.meshgrid(x,y)
    xx = xx.flatten().astype('g') #astype np.longdoubles to avoid overflow
    yy = yy.flatten().astype('g')

    truth = img.flatten()

    max_order=order
    datamat = np.array([xx, yy],dtype='g').transpose()
    current_order=1
    while current_order < max_order:
        datamat = np.hstack([datamat,xx[:,np.newaxis]**current_order, yy[:,np.newaxis]**current_order])
        current_order += 1
    clf = LinearRegression()
    clf.fit(datamat,truth)
    predicted = clf.predict(datamat).reshape((img.shape[0], img.shape[1]))
    return predicted.astype(float)

def makeCorr(img_to_use, rescale_cmap=True, N=5):
    img_to_use = img_to_use - np.mean(img_to_use)
    corr_arr = np.zeros((N * 2 + 1, N * 2 + 1))

    for y in np.arange(N * 2 + 1):
        for x in np.arange(N * 2 + 1):
            corr_arr[y,x] = np.sum(np.roll(np.roll(img_to_use,y-N,axis=0),x-N,axis=1)*img_to_use)/np.sum(img_to_use*img_to_use)

    fig,ax = plt.subplots()
    b = np.max(corr_arr)
    a = np.min(corr_arr)
    c = 0
    midpoint = (c - a) / (b - a)
    if rescale_cmap:
        cmap = shiftedColorMap(plt.cm.RdBu_r,midpoint=midpoint)
        vmin=a
        vmax=b
    # elif np.all(corr_arr >= 0):
    #     cmap = plt.cm.Reds
    #     vmin=None
    #     vmax=None
    # elif np.all(corr_arr <= 0):
    #     cmap = plt.cm.Bu_r
    #     vmin=None
    #     vmax=None
    else:
        cmap = plt.cm.RdBu_r
        vmin=None
        vmax=None
    plt.imshow(corr_arr,interpolation='None',cmap=cmap,extent=(-N,N,-N,N),
               vmin=vmin, vmax=vmax)
    plt.title('Pixel-Neighbor Correlations')
    #eventually print text in pixels...
    #for x_val, y_val in zip(np.arange(11), np.arange(11)):
    #    c = np.round(100*corr_arr[y_val,x_val],0)
    #    ax.text(y_val, x_val, c, va='center', ha='center')
    cbar = plt.colorbar()
    print np.round(100*corr_arr[N-2:N+3,N-2:N+3],1)
    #plt.figure()
    #tmp = plt.hist(img_to_use.flatten(),bins=20)
    return corr_arr

def neighborScatter(img):
    crop = np.zeros([img.shape[0]-2,img.shape[1]-2])
    vertSum = np.zeros([img.shape[0]-2,img.shape[1]-2])
    horizSum = np.zeros([img.shape[0]-2,img.shape[1]-2])
    diagSum = np.zeros([img.shape[0]-2,img.shape[1]-2])
    for i in np.arange(img.shape[0]-2):
        for j in np.arange(img.shape[1]-2):
            crop[i,j] = img[i,j]
            vertSum[i,j] = img[i-1,j] + img[i+1,j]
            horizSum[i,j] = img[i,j-1] + img[i,j+1]
            diagSum[i,j] = img[i-1,j-1] + img[i-1,j+1] + img[i+1,j+1] + img[i+1,j-1]
    plt.figure()
    plt.title('Vertical Neighbors')
    plt.hist2d(crop.flatten(),vertSum.flatten(),bins=50)
    plt.xlabel('Pixel value')
    plt.ylabel('Sum of vertical neighbors')
    print 'vert corr: ' + str(sp.stats.pearsonr(crop.flatten(),vertSum.flatten()))
    plt.figure()
    plt.title('Horizontal Neighbors')
    plt.hist2d(crop.flatten(),horizSum.flatten(),bins=50)
    plt.xlabel('Pixel value')
    plt.ylabel('Sum of horizontal neighbors')
    print 'horiz corr: ' + str(sp.stats.pearsonr(crop.flatten(),horizSum.flatten()))
    plt.figure()
    plt.title('Diagonal Neighbors')
    plt.hist2d(crop.flatten(),diagSum.flatten(),bins=50)
    plt.xlabel('Pixel value')
    plt.ylabel('Sum of diagonal neighbors')
    print 'diag corr: ' + str(sp.stats.pearsonr(crop.flatten(),diagSum.flatten()))

def moverFlatTestPlot(mover, title=None):
    source = Source(num_x=15, min_x=-40, max_x=40,
                    num_y=15)
    # calculate uniform illumination flux
    UniformIlluminationMover()(source)
    old_mean_flux = source.fluxes.mean()
    source.fluxes -= old_mean_flux  # so that we have 0 flux initially
    mover(source)

    # make the fluxes relative to old mean flux
    source.fluxes = (source.fluxes - old_mean_flux) / old_mean_flux

    fig, ax = source.plot_vertices()
    ax.set_title('Vertices {0}'.format(title))

    fig, ax = source.plot_real_grid()
    ax.set_title('Real Coordinates {0}'.format(title))

    fig, ax = source.plot_pixel_grid()
    ax.set_title('Pixel Coordinates {0}'.format(title))
