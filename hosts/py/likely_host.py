import astropy
from astropy.io import ascii,fits
from astropy import wcs
from astropy.table import vstack, Table
import numpy as np
import glob
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
import warnings
import copy
warnings.filterwarnings("ignore")

cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Ob0=0.05)
snwcsunits,catwcsunits = (u.hourangle,u.deg), u.deg
scale = .06
                            
def stack(path):
    # to stack catalogs on top of each other
    # will only work if the format is the same i.e. identical colnames
    tabs = []
    cats = glob.glob(path)
    for i in cats:
        tab = ascii.read(i)
        tabs.append(tab)
    stacked_cats = vstack(tabs)
    return stacked_cats

def set(path,types=None):
    # to get list of all catalogs on a path
    # uses meta to give names to each catalog in list using the filname
    # optionally specify which types of files you want

    # types ~ default None will grab all formats in the path folder
    if types==None:
        catset = glob.glob(path) 
        catname = copy.copy(catset) # list of the globbed files on path
        tmp=  []
        j = 0
        for i in catset:
            try:
                tmp.append(ascii.read(i))
                tmp[j].meta['catname_'+str(j)] = catname[j] # give meta property catname to each catalog
                # TODO ignore the path and just take filename for meta, ie final slice string after final / 
            except:
                # maybe read-me don't bother
                continue
            j+=1
        catset = tmp
        return catset
    
    # optionally give types as list of desired formats; eg ['.txt','.cat']
    elif types:
        tmp = [path+'/*'+i for i in types] 
        #catset = [glob.glob(e)[0] for e in [path+'/*.RF', path+'/*.cat']] # grabs the .cat and .RF catalogs 
        catset = [glob.glob(e)[0] for e in tmp] # grabs the .cat and .RF catalogs 
        # TODO the [0] here grabs just the first of each type; fine for my application
        # however should really use another loop to keep this general 
        catname = copy.copy(catset)
        tmp=  []
        j = 0
        for i in catset:
            try:
                tmp.append(ascii.read(i))
                tmp[j].meta['catname_'+str(j)] = catname[j]
            except:
                # maybe read-me don't bother
                continue
            j+=1
        catset = tmp
        return catset

# star flags, 1=star
def prune(cat,prune_colnames=['class_star','stel','pointsource'],verbose=None):
    # to take out stars from catalog or list of catalogs
    
    if type(cat) == list:
        tmps,n = [],0 # will loop through each cat in list and attempt prune then stick back into list  
        for i in cat:
            colnames = i.colnames
            #colnames = [name.lower() for name in colnames]
            for j in colnames:
                if j.lower() in prune_colnames:
                    tmp = i[i[j] < 0.1]
                    n+=1 # counting the successful prunes of catalogs in list
                    break
                else:
                    tmp = i
                    continue
            tmps.append(tmp)
        if verbose:
            print('pruned {}/{} of the catalogs in list'.format(n,len(cat)))
        return tmps
    else: # cat ~ single table; not a list 
        for i in cat.colnames:
            if i.lower() in prune_colnames:
                tmp = cat[cat[i] < 0.1]
                if verbose:
                    print('pruned catalog from {} --> {}'.format(len(cat),len(tmp)))
            else:
                continue
        return tmp

def candels_fields(field='COS',sn=ascii.read('../../sn/candels.csv'),name='name'):
    tmp = []
    for i in sn:
        if field in i[name]:
            tmp.append(i)

    tmp = vstack(tmp)
    return tmp

# would maybe be useful to raise a warning if neighbors were > X arcsec away
def neighbors(sn,catfile,N=3,snradec=('ra','dec'),catradec=('ra','dec'),
    snwcsunits = (u.hourangle,u.deg),catwcsunits = u.deg,do_diff = True,verbose=True,
    ID = 'id'):
    # find N nearby galaxies using just ra dec (astropy match to catalog sky fcn)
    # returns table with these N neighbors 
    # eventually needs a single catalog ~ astropy table object to look through
    # however can enter str of path to file or a list of catalogs 
    # will try to read the catalog on path and/or stack the list (won't succeed if different colnames)

    snra,sndec = snradec[0],snradec[1]
    catra,catdec = catradec[0],catradec[1]
    sn_coord = SkyCoord(sn[snra],sn[sndec],unit=snwcsunits)
    if type(catfile)==str:
        cat = ascii.read(catfile)
        cat_coords = SkyCoord(cat[catra],cat[catdec],unit=catwcsunits)

    elif type(catfile)==list:
        # set ~ list of tables from a source
        cat = catfile 
        for i in cat:
            try:
                cat_coords = SkyCoord(i[catra],i[catdec],unit=catwcsunits)
                cat = i
                break
            except:
                continue

    else: # is single table
        cat = catfile
        cat_coords = SkyCoord(cat[catra],cat[catdec],unit=catwcsunits)
    
    #print(type(cat))
    """
    if type(cat) == list:
        try: # if give list of catalog tables try and stack them all
            # relevant to ff,relics which were seperated into the clusters
            cat = vstack([i for i in cat])
        except:
            print('no can stack this list')        
    """


    tmp = []
    for i in range(N):
        idx = sn_coord.match_to_catalog_sky(cat_coords,nthneighbor=i+1)
        match = cat[idx[0]]
        tmp.append(match)
    
    tmp = vstack(tmp) # table of neighbors 
    neighbor_coords = SkyCoord(tmp[catra],tmp[catdec],unit=catwcsunits)
    if type(catfile)==list:
        # if gave it a list of tables use match id to return set of neighbors
        tmp = [match_id(i,catfile,ID) for i in tmp] # [[set neighbor 1],...,[set neighbor N]]
        # optional todo vstack? so we have like [[N neighbors from cat 1],...[N neighbors from len(set)]]
    
    if do_diff:
        dra, ddec = sn_coord.spherical_offsets_to(neighbor_coords)
        sep = sn_coord.separation(neighbor_coords)
        if verbose:
            print('{} neighbors, w seps ~ {} arcsec'.format(len(tmp),sep.arcsec))
        return tmp
    else:
        return tmp

def match_id(obj,files,ID,**kwargs):
    # to get all the catalogs from set of files for a given obj
    # uses id in obj to match w other catalogs
    # if multiple id colnames, enter as a list w the obj id as zeroth initial entry

    # read in the globbed files
    if type(files[0]) == str:
        cats = []
        for i in files:
            cat = ascii.read(i)
            cats.append(cat)
    else: # don't want to read each time use in for loop
        cats = files # list of tables

    matches = []
    if type(ID) == str:
        obj_id = obj[ID]
        for cat in cats:
            match = cat[cat[ID] == obj_id]
            matches.append(match)
    else: # list of possible ids to try
        obj_id = obj[ID[0]] # put the obj id first in list
        for cat in cats:
            for i in ID:
                try:
                    match = cat[cat[i] == obj_id]
                    matches.append(match)
                except:
                    continue   
    return matches


def abt(params,deg=True):
    # get cij from a,b,theta
    a,b,theta = params
    if deg:
        theta*=np.pi/180
    cxx = (np.cos(theta)/a)**2 + (np.sin(theta)/b)**2
    cyy = (np.sin(theta)/a)**2 + (np.cos(theta)/b)**2
    cxy = 2*np.sin(theta)*np.cos(theta)*(1/a**2 - 1/b**2)
    return (cxx,cyy,cxy)
def eAt(params,deg=True):
    # get cij from ell,area,theta
    ell,area,theta = params
    if deg:
        theta*=np.pi/180 
    # ell ~ (a-b)/a, area ~ pi*a*b (pixels)  
    if deg:
        theta*=np.pi/180
    a = (area/(np.pi*(1-ell)))**.5
    b = area/(np.pi*a)
    # sanity checks to look against table values
    #check_area = np.pi*a*b
    #check_ell = (a-b)/a
    cxx = (np.cos(theta)/a)**2 + (np.sin(theta)/b)**2
    cyy = (np.sin(theta)/a)**2 + (np.cos(theta)/b)**2
    cxy = 2*np.sin(theta)*np.cos(theta)*(1/a**2 - 1/b**2)
    return (cxx,cyy,cxy)
def eA(params):
    # TODO don't have theta assigning zero 
    # get cij from ell,area
    theta = 0
    ell,area = params
    a = (area/(np.pi*(1-ell)))**.5
    b = area/(np.pi*a)
    # sanity checks
    #check_area = np.pi*a*b
    #check_ell = (a-b)/a
    cxx = (np.cos(theta)/a)**2 + (np.sin(theta)/b)**2
    cyy = (np.sin(theta)/a)**2 + (np.cos(theta)/b)**2
    cxy = 2*np.sin(theta)*np.cos(theta)*(1/a**2 - 1/b**2)
    return (cxx,cyy,cxy)


def eff_off(sn,obj,snradec,catradec,objparams=None,method=None,verbose=True,**kwargs):
    # scale ~ .06''/pix for the archive catalogs where most of the parameters come from
    snra,sndec=snradec[0],snradec[1]
    catra,catdec=catradec[0],catradec[1]
    sn_coord = SkyCoord(sn[snra],sn[sndec],unit=snwcsunits)
    if type(obj) != list:
        # we have a single table obj w the wcs and sextractor params
        obj_coord = SkyCoord(obj[catra],obj[catdec],unit=catwcsunits)
        sex_params = [obj[i] for i in objparams]
    else:
        # set ~ we have a list of table objects for the neighbor
        # this expected as norm, will be feeding in sets to the host in ipynb 
        for table in obj:
            try:
                obj_coord = SkyCoord(table[catra],table[catdec],unit=catwcsunits)
                break
            except:
                continue
        for table in obj:
            try:
                sex_params = [table[i] for i in objparams]
                break
            except:
                continue

    sn_ra,sn_dec = sn_coord.ra, sn_coord.dec
    obj_ra,obj_dec = obj_coord.ra, obj_coord.dec
    dra, ddec = sn_ra - obj_ra, sn_dec - obj_dec
    x,y = dra.arcsec*scale**-1,ddec.arcsec*scale**-1 # pixels; scale conv ~ ''/pixel for the catalog
    x,y = abs(x),abs(y) # TODO is this true or is sign important? 
    
    # convert sex_parms to cij in /pixel^2
    # these methods are all for sextractor params (i.e. provided in units of pixels not physical) 
    methods = {'abt':abt(sex_params),'eAt':eAt(sex_params),'cij':sex_params}
    # TODO 'eA':eA(sex_params) 
    cxx,cyy,cxy = methods[method] # /pixel^2

    eo = np.sqrt(cxx*(x**2) + cyy*(y**2) + cxy*(x*y))
    if type(eo) == astropy.table.column.Column:
        eo = eo[0]
    
    if verbose:
        try:
            catname = obj.meta['catname']
            print(catname+' eff offset ~ {}'.format(eo))
        except:
            print('eff offset ~ {}'.format(eo))
    
    return eo

def host(sn,neighbors,snradec,catradec,objparams=None,method=None,verbose=True,no_host=False,**kwargs):
    offsets = []
    # neighbors either table with N neighbors; or list length N neighbors that have set of tables in lists 
    for i in range(len(neighbors)): 
        neighbor = neighbors[i]
        tmp = eff_off(sn,neighbor,snradec,catradec,objparams=objparams,method=method,**kwargs)
        offsets.append(tmp)

    li = np.argmin(offsets)
    if verbose:
        eo = offsets[li]
        print('host eo ~ {}'.format(eo))
    
    if no_host:
        # optional argument to id those wo a likely host
        # will return none if offset > 5
        # default false since easy enough to get w one liner from 'best' host (min eo) if skip
        eo = offsets[li]
        if eo > 5:
            likely_host = None
        else:
            likely_host = [neighbors[li],eo]
    else:
        likely_host = [neighbors[li],eo]
    
    return likely_host

# loop_pickle here doesn't work atm; rather is achieved in jupyter nb
# if want to revisit later: *args is what needs work and see host.py 
def loop_pickle(name=None,sne=None,verbose=None,*args):
    print('loop_pickle')
    for arg in args:
        print(arg)

    j=0
    neighbors_list,likely_hosts = [],[]
    if verbose:
        print(name)
        print('length of SNe ~ {}'.format(len(sne)))
        print('')

    if name=='cos' or name=='egs' or name=='goodsn' or name=='goodss' or name=='uds':
        for sn in sne:
            if verbose:
                print('working on {}'.format(j))
            j+=1
            if verbose:
                print('3dhst')
            neighbors_3dhst = neighbors(sn,catset_3dhst,snradec=snradec,catradec=catradec_3dhst,ID=ID_3dhst)
            likely_host_3dhst = host(sn,neighbors_3dhst,snradec=snradec,catradec=catradec_3dhst,
                                                 objparams=objparams_3dhst,method=method_3dhst)
            if verbose:
                print('archive')
            neighbors_archive = neighbors(sn,catset_archive,snradec=snradec,catradec=catradec_archive,ID=ID_archive)
            likely_host_archive = host(sn,neighbors_archive,snradec=snradec,catradec=catradec_archive,
                                                  objparams=objparams_archive,method=method_archive)
            neighbors_list.append([sn,neighbors_archive,neighbors_3dhst])
            likely_hosts.append([sn,likely_host_archive,likely_host_3dhst]) 
            if verbose:
                print('')

    elif name=='clash':
        for sn in sne:
            if verbose:
                print('working on {}'.format(j))
            j+=1
            if verbose:
                print('fullfield')
            neighbors_fullfield = neighbors(sn,catfile=catset_fullfield,snradec=snradec,catradec=catradec_fullfield,ID=ID_fullfield)
            likely_host_fullfield = host(sn,neighbors_fullfield,snradec=snradec,catradec=catradec_fullfield,
                                                    objparams=objparams_fullfield,method=method_fullfield)
            if verbose:
                print('molino')
            neighbors_molino = neighbors(sn,catfile=catset_molino,snradec=snradec,catradec=catradec_molino,ID=ID_molino)
            likely_host_molino = host(sn,neighbors_molino,snradec=snradec,catradec=catradec_molino,
                                                 objparams=objparams_molino,method=method_molino)
            if verbose:
                print('archive (likely host none dont have any way for theta)')
            neighbors_archive = neighbors(sn,catfile=catset_archive,snradec=snradec,
                                                      catradec=catradec_archive,ID=ID_archive)
            likely_host_archive = None
            """
            # have an area and ellipticity but no way for theta cant do an eo
            likely_host_archive = host(sn,neighbors_archive,snradec=snradec,catradec=catradec_archive,
                                                  objparams=objparams_archive,method=method_archive)
            """
            neighbors_list.append([sn,neighbors_archive,neighbors_molino,neighbors_fullfield])
            likely_hosts.append([sn,likely_host_archive,likely_host_molino,likely_host_fullfield])
            if verbose:
                print('')

    elif name=='ff':
        for sn in sne:
            if verbose:
                print('working on {}'.format(j))
            j+=1
            if verbose:
                print('archive')
            neighbors_archive = neighbors(sn,catfile=catset_archive,snradec=snradec,catradec=catradec_archive,ID=ID_archive)
            if verbose:
                print('I dont have parameters to do eo for these host is using nearest neighbor')
            likely_host_archive = neighbors(sn,catfile=catset_archive,snradec=snradec,catradec=catradec_archive,N=1,ID=ID_archive)
            neighbors_list.append([sn,neighbors_archive])
            likely_hosts.append([sn,likely_host_archive])

    elif name=='relics':    
        for sn in sne:
            if verbose:
                print('working on {}'.format(j))
            j+=1
            if verbose:
                print('archive')
            neighbors_archive = neighbors(sn,catfile=catset_archive,snradec=snradec,catradec=catradec_archive,ID=ID_archive)
            likely_host_archive = host(sn,neighbors_archive,snradec=snradec,catradec=catradec_archive,
                                                  objparams=objparams_archive,method=method_archive)

            neighbors_list.append([sn,neighbors_archive])
            likely_hosts.append([sn,likely_host_archive])
            
    pkl,pickleas = True, name+'_hosts.pkl'
    if pkl:
        pickle.dump([likely_hosts,neighbors_list],open(pickleas,'wb'))
        if verbose:
            print('pickle dumped ~ {}'.format(pickleas))