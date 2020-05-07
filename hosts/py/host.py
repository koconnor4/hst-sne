from likely_host import *

# adjust the truth of these to decide which you want to do rn
cos,egs,goodsn,goodss,uds = False,False,False,False,False
candels = False # will do all the candels
clash = False
relics = False
ff = True
_all = False # will do all 
if candels:
    cos,egs,goodsn,goodss,uds=True,True,True,True,True

if _all:
    cos,egs,goodsn,goodss,uds,clash,relics,ff=True,True,True,True,True,True,True,True

snpath = '../../sn/'
path = '../catalogs/'

if ff:
    name = 'ff' 
    sne = ascii.read(snpath+'ffsn.csv')
    snradec = ('ra','dec')
    catset_archive = set(path+'ff_catalogs/archive_stsci.cats/*')
    catset_archive = prune(catset_archive)
    print(name,': have archive catalogs')
    ID_archive = ('id','col1')
    catradec_archive = ('RA','Dec')
    objparams_archive = ('ell','area') # I don't have theta, cant do an eo
    method_archive = 'eA'
    print('I dont have theta for archive cant do an eo')
    loop_pickle(name=name,sne=sne,verbose=True)

if relics:
    name = 'relics' 
    sne = ascii.read(snpath+'relics.csv')
    snradec = ('ra','dec')
    catset_archive = set(path+'relics_catalogs/archive_stsci.cats/*')
    catset_archive = prune(catset_archive)
    print(name,': have archive catalogs')
    ID_archive = ('id','col1')
    catradec_archive = ('RA','Dec')
    objparams_archive = ('ell','area','theta')
    method_archive = 'eAt'
    loop_pickle(name=name,sne=sne,verbose=True)

    
if clash:
    name = 'clash' # for pickle
    sne = ascii.read(snpath+'clash.csv')
    snradec = ('ra','dec')
    catset_fullfield = stack(path+'clash_catalogs/fullfield.cats/*.cat')
    catset_molino = stack(path+'clash_catalogs/archive_stsci.cats/molino/*')
    catset_archive = set(path+'clash_catalogs/archive_stsci.cats/hst/*')
    catset_molino = prune(catset_molino)
    catset_fullfield = prune(catset_fullfield)
    catset_archive = prune(catset_archive)
    print(name,': have archive, molino, and fullfield catalogs')
    ID_molino = 'CLASHID'
    ID_archive = 'id'
    ID_fullfield = 'id'
    catradec_molino=('RA','Dec') 
    catradec_archive = ('RA','Dec')
    catradec_fullfield = ('ALPHA_J2000','DELTA_J2000')
    objparams_molino=('a','b','theta')
    objparams_archive = ('ell','area') # no theta
    objparams_fullfield = ('CXX_IMAGE','CYY_IMAGE','CXY_IMAGE')
    method_molino = 'abt'
    method_archive = 'eA'
    method_fullfield = 'cij'
    print('archive has no theta cant do an eo using that catalog')
    loop_pickle(name=name,sne=sne,verbose=True)

    
if cos:
    name = 'cos' # for pkl
    sne = candels_fields(field='COS')
    snradec = ('RA','DEC')
    catset_3dhst = set(path+'cos_catalogs/*3d*/*',types=['.RF','.cat','.fout'])
    catset_archive = set(path+'cos_catalogs/archive_stsci.cats/*')
    catset_archive = prune(catset_archive)
    catset_3dhst = prune(catset_3dhst)
    print(name,': have 3dhst and archive catalogs')
    ID_archive = 'ID'
    ID_3dhst = 'id'
    catradec_archive=('RA','DEC') 
    catradec_3dhst = ('ra','dec')
    objparams_archive = ('CXX_IMAGE','CYY_IMAGE','CXY_IMAGE')
    objparams_3dhst = ('a_image','b_image','theta_J2000')
    method_archive = 'cij'
    method_3dhst = 'abt'
    loop_pickle(name=name,sne=sne,verbose=True)

    
if egs:
    name = 'egs' # to pkl
    sne = candels_fields(field='EG')
    snradec = ('RA','DEC')
    catset_3dhst = set(path+'egs_catalogs/*3d*/*',types=['.RF','.cat','.fout']) 
    catset_archive = set(path+'egs_catalogs/archive_stsci.cats/*')
    catset_archive = prune(catset_archive)
    catset_3dhst = prune(catset_3dhst)
    print(name,': have 3dhst and archive catalogs')
    ID_archive = ('ID','col1')
    ID_3dhst = 'id'
    catradec_archive=('RA','DEC') 
    catradec_3dhst = ('ra','dec')
    objparams_archive = ('CXX_IMAGE','CYY_IMAGE','CXY_IMAGE')
    objparams_3dhst = ('a_image','b_image','theta_J2000')
    method_archive = 'cij'
    method_3dhst = 'abt'
    loop_pickle(name=name,sne=sne,verbose=True)

if goodsn: 
    name = 'goodsn' # to pkl
    sne = candels_fields(field='GN')
    snradec = ('RA','DEC')
    catset_3dhst = set(path+'goodsn_catalogs/*3d*/*',types=['.RF','.cat','.fout']) 
    catset_archive = set(path+'goodsn_catalogs/archive_stsci.cats/*')
    catset_3dhst = prune(catset_3dhst)
    catset_archive = prune(catset_archive)
    print(name,': have 3dhst and archive catalogs')
    ID_archive = ('ID','col1')
    ID_3dhst = 'id'
    catradec_archive=('RA','DEC') 
    catradec_3dhst = ('ra','dec')
    objparams_archive = ('CXX_IMAGE','CYY_IMAGE','CXY_IMAGE')
    objparams_3dhst = ('a_image','b_image','theta_J2000')
    method_archive = 'cij'
    method_3dhst = 'abt'
    loop_pickle(name=name,sne=sne,verbose=True)

if goodss: 
    name = 'goodss' # for pkl
    sne_D = candels_fields(field='GSD')
    sne_W = candels_fields(field='GSW')
    sne = vstack([sne_D,sne_W])
    snradec = ('RA','DEC')
    catset_3dhst = set(path+'goodss_catalogs/*3d*/*',types=['.RF','.cat','.fout']) 
    catset_archive = set(path+'goodss_catalogs/archive_stsci.cats/*')
    catset_3dhst = prune(catset_3dhst)
    catset_archive = prune(catset_archive)
    print(name,': have 3dhst and archive catalogs')
    ID_archive = ('ID','Seq')
    ID_3dhst = 'id'
    catradec_archive=('RA','DEC') 
    catradec_3dhst = ('ra','dec')
    objparams_archive = ('CXX_IMAGE','CYY_IMAGE','CXY_IMAGE')
    objparams_3dhst = ('a_image','b_image','theta_J2000')
    method_archive = 'cij'
    method_3dhst = 'abt'
    loop_pickle(name=name,sne=sne,verbose=True)

if uds:
    name = 'uds' # for pickling
    sne = candels_fields(field='UDS')
    snradec = ('RA','DEC')
    catset_3dhst = set(path+'uds_catalogs/*3d*/*',types=['.RF','.cat','.fout']) 
    catset_archive = set(path+'uds_catalogs/archive_stsci.cats/*')
    catset_3dhst = prune(catset_3dhst)
    catset_archive = prune(catset_archive)
    print(name,': have 3dhst and archive catalogs')
    ID_archive = ('ID','Seq')
    ID_3dhst = 'id'
    catradec_archive=('RA','Dec') 
    catradec_3dhst = ('ra','dec')
    objparams_archive = ('cxx_image','cyy_image','cxy_image')
    objparams_3dhst = ('a_image','b_image','theta_J2000')
    method_archive = 'cij'
    method_3dhst = 'abt'
    loop_pickle(name=name,sne=sne,verbose=True)

