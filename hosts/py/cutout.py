import copy
import pickle
import numpy as np
import glob

# matplotlib
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

# astropy
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.utils.data import download_file
from astropy.coordinates import SkyCoord
from astropy import wcs
import astropy.units as u

# scipy
from scipy.optimize import fsolve
import math

# warnings
import warnings 
warnings.simplefilter("ignore")


def cutout(filename,position, size, cutout_filename = 'example_cutout.fits'):
    
    # Download the image
    # commented out since I already have file, could be useful just change filename argument to url
    # filename = download_file(url)


    # Load the image and the WCS
    hdu = fits.open(filename)[0]
    w = wcs.WCS(hdu.header)
    
    # w.sip = None temporarily removes sip coefficients from a wcs obj; image distortion
    # https://docs.astropy.org/en/stable/wcs/note_sip.html
    w.sip = None

    # Make the cutout, including the WCS
    cut = Cutout2D(hdu.data, position=position, size=size, wcs=w)

    # Put the cutout image in the FITS HDU 
    hdu.data = cut.data

    # Update FITS header with the cutout WCS
    hdu.header.update(cut.wcs.to_header())

    # Write the cutout to a new FITS file
    hdu.writeto(cutout_filename, overwrite=True)


def cxx(p):
    a,b,theta = p
    return (np.cos(theta)/a)**2 + (np.sin(theta)/b)**2
def cyy(p):
    a,b,theta = p
    return (np.sin(theta)/a)**2 + (np.cos(theta)/b)**2
def cxy(p):
    a,b,theta=p
    return 2*np.sin(theta)*np.cos(theta)*(1/a**2 - 1/b**2)

def equations(p):
    cxx,cyy,cxy=1,2,3
    #cxx,cyy,cxy,a,b,theta = p
    a,b,theta = p
    return (cxx - np.cos(theta)/a)**2 - (np.sin(theta)/b)**2, cyy - (np.sin(theta)/a)**2 - (np.cos(theta)/b)**2, cxy - 2*np.sin(theta)*np.cos(theta)*(1/a**2 - 1/b**2)

def b_theta_equations(p,cxx,cyy,cxy):
    #cxx,cyy,cxy=cij
    #cxx,cyy,cxy=4,5,7
    b,theta=p
    return (2*np.tan(theta)*cxx-cxy-2*np.tan(theta)*np.sin(theta)**2/b**2 - 2*np.cos(theta)*np.sin(theta)/b**2, 2*np.tan(theta)**-1*cyy-cxy-2*np.tan(theta)**-1*np.cos(theta)**2/b**2-2*np.cos(theta)*np.sin(theta)/b**2)

def a_theta_equations(p,b,cxx,cyy,cxy):
    #cxx,cyy,cxy=cij
    #cxx,cyy,cxy=4,5,7
    a,theta=p
    return cxx-np.cos(theta)**2/a**2-np.sin(theta)**2/b**2,cyy-np.cos(theta)**2/b**2-np.sin(theta)**2/a**2

"""
# tst_case to remind of steps on how to use above to get from cij to a,b,theta
cij = (3,6,4)
cxx,cyy,cxy=cij

init_p=(2,2.5) # b,theta initial guess
test_p=fsolve(b_theta_equations,init_p,args=(cxx,cyy,cxy)) # b,theta best solutions
print(test_p) 
b,theta=test_p
print(b_theta_equations(test_p,cxx,cyy,cxy)) # if not ~ 0 soln didnt do well

init_p=(1,2) # a init guess
test_p=fsolve(a_theta_equations,init_p,args=(cxx,cyy,cxy)) # a best soln
print(test_p)
a,theta=test_p
print(a_theta_equations(test_p,cxx,cyy,cxy))

print(cxx((a,b,theta)),cyy((a,b,theta)),cxy((a,b,theta)))
"""

def ellipse(file,a,b,theta,gal_position,sn_position,conv=None,sip=True,title='', save = False,saveas='saved.pdf',show=False,
    logscale=True,
    diverging=True,
    color_idx = 0,
    val_min = 0.015,
    val_max=1.15):
    
    # access the file
    fits_image = fits.open(file)
    print(file)
    hdu = fits_image[0]
    if conv:
        # the deg/pixel for the fits image along ra and dec
        # the ellipse (add_patch) needs to know pixel dimensions of a,b to plot on the image 
        # a,b in units of pixels from the drz hst wfc3ir images ~ .065''/pixel
        std_res = .065 # ''/pixel 
        ra_conv = hdu.header.get('CD1_1')*u.deg
        dec_conv = hdu.header.get('CD1_1')*u.deg
        ra_conv,dec_conv=ra_conv.to(u.arcsec),dec_conv.to(u.arcsec)
        ra_conv,dec_conv =ra_conv.value,dec_conv.value # arcsec
        # the ra,dec convs ought to be the same use either one
        a,b=np.abs(a*std_res/ra_conv),np.abs(b*std_res/ra_conv)    
    # the wcs obj from header of the image;
    w = wcs.WCS(hdu.header)
    # the image is drizzled ~ distortion corrected, and should be  w.sip = None already
    # just doing it again for consistency/peace of mind
    if sip:
        warnings.simplefilter('ignore')
        w.sip = None

    # initialize the figure
    f= plt.figure()
    #print('w~{}'.format(w))
    ax=f.gca(projection=w)
    # choosing how to normalize the image 
    
    # the image min,max values 
    #img_min = np.abs(np.amin(hdu.data))
    #img_max = np.amax(hdu.data)+2*np.abs(np.amin(hdu.data))
    
    if logscale == True: # is the default
        norm_scale = matplotlib.colors.LogNorm(vmin=val_min,vmax=val_max)
    else: # will go to linear normalization, may be better go to fits and zscale to determine good vmin,vmax 
        norm_scale = matplotlib.colors.Normalize(vmin=val_min,vmax=val_max)
    
    # choosing which color map 
    if diverging == True: # is the default
        color_map = ['PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
            'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']
    else: # will go to sequential here 
        color_map =       ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
            'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
            'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']
    

    ax.imshow(hdu.data,norm=norm_scale,cmap=color_map[color_idx])

    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel(r'$\delta$')
    ax.set_title(title)
    
    sn_pix_coords = wcs.utils.skycoord_to_pixel(sn_position, w,origin=1,mode='all')
    #TODO the pixel physical vs image x,y on clash cutouts is sometimes very different
    #for some reason this is causing an issue when using SN skycoord to pixel
    #baffled as to how the same exact fcn is working for the host though 
    #print('got sn pixel coords ~ {}'.format(sn_pix_coords))
    #ax.scatter([int(sn_pix_coords[0])],[int(sn_pix_coords[1])],marker='X',c='black',edgecolor='white',label='SN')
    pix_coords = []
    
    # now into the possible hosts to get major,minor,position angle and plot ellipse at loc 
    # get the pixel coords for the center of the host using wcs
    # origin 0 or 1?
    gal_pixel_coords = wcs.utils.skycoord_to_pixel(gal_position, w, origin=1, mode='all')
    #print('got gal pixel coords ~ {}'.format(gal_pixel_coords))
    x = int(gal_pixel_coords[0])
    y = int(gal_pixel_coords[1])
    # the major,minor are given in pixels...
    # a/b,theta [pixel,deg ccw from x]
    colors = ['red','blue','green','cyan','violet','turquoise','purple','yellow','orange','lime','navy','pink','saddlebrown','silver','gold','navajowhite']
    c = colors[0]
    #print('a,b',a,b)
    ax.add_patch(Ellipse((x,y),width = 5*a, height = 5*b, angle = theta,color='red',facecolor=None,fill=False))
    #print('made ellipse')
    plt.legend(loc='upper right', bbox_to_anchor=(1.4, 1.0))
    if save == True:
        plt.savefig(saveas)
    if show == True:
        plt.show()
    return f