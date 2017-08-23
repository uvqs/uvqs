# Module for doing things with the UVQS and other samples

# Imports
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import glob, os, sys
import pdb

from pkg_resources import resource_filename

from astropy.io import ascii, fits
from astropy import units as u
from astropy.table import Table, Column, QTable
from astropy import coordinates as coords
from astropy.coordinates import SkyCoord
from astropy.coordinates.angles import Angle

from linetools import utils as ltu


####
def get_candidates(cand_fil=None, verbose=True):
    """
    Read in QSO candidates and prioritize according to FUV
    """
    if cand_fil is None:
        cand_fil = resource_filename('uvqs', 'data/candidates/z1qso_ALL_FN0.3_W0.6_candidates.fits.gz')
    if verbose:
        print('z1qso_analy: Reading z1qso candidates from {:s}'.format(cand_fil))

    # Load candidates
    cand = Table(fits.open(cand_fil)[1].data)
    ncand = len(cand)

    # Prioritize (by FUV mag)
    priority = np.zeros(ncand, dtype='int')

    gd1 = np.where(cand['FUV'] < 17.5)[0] 
    priority[gd1] = 1
    gd2 = np.where( (cand['FUV'] < 18.0) & priority == 0)[0] 
    priority[gd2] = 2
    gd3 = np.where( cand['FUV'] > 18.0 )[0]
    priority[gd3] = 3

    cand.add_column( Column( priority, name='PRI') )
    
    return cand

####
def load_casbah(clobber=False):
    ''' 
    Read in QSOs from CASBAH
    '''
    import glob

    casbah_fil = os.getenv('z1Q')+'/Analysis/'+'casbah_qsos.fits'

    # Check for file and return it as desired
    files = glob.glob(casbah_fil)
    if (len(files) > 0) & (not clobber): # Read file
        casbah = Table( fits.open(casbah_fil)[1].data )
        return casbah
    
    # Generate the file
    
    # List for convenience
    lcasbah = [ ['FBQS 0751+2919', '07:51:12.3062','+29:19:38.271'],
               ['PG1148+549', '11:51:20.4625', '+54:37:33.078'], 
               ['PG1206+459', '12:08:58.0129', '+45:40:35.487'], 
               ['PG1338+416', '13:41:00.7747', '+41:23:14.072'], 
               ['PG1407+265', '14:09:23.8979', '+26:18:21.043'], 
               ['PG1522+101', '15:24:24.5538', '+09:58:29.774'],
               ['PG1630+377', '16:32:01.1077', '+37:37:49.984'], 
               ['PHL1377',    '02:35:07.34',   '-04:02:05.4'], 
               ['LBQS 1435-0134', '14:37:48.2871','-01:47:10.781'] 
               ]
    ncasbah = len(lcasbah)
    RA = Column(np.zeros(ncasbah), name='RA')
    DEC = Column(np.zeros(ncasbah), name='DEC')

    names = []
    # Make a Table
    for ii,qso in enumerate(lcasbah):
        names.append(qso[0])
        coord = ltu.radec_to_coord((qso[1],qso[2]))
        RA[ii], DEC[ii] = coord.ra.value, coord.dec.value

    name = Column(names, name='QSO')
    casbah = Table([name,RA,DEC])

    # Cross match to milliq
    milliq = load_milliq()

    c_casbah = SkyCoord(ra=casbah['RA']*u.degree, dec=casbah['DEC']*u.degree)
    c_mq = SkyCoord(ra=milliq['RA']*u.degree, dec=milliq['DEC']*u.degree)

    idx, d2d, d3d = coords.match_coordinates_sky(c_casbah, c_mq, nthneighbor=1)

    # Check for no match
    bad = np.where(d2d.to('arcsec').value > 5.)[0] 
    if len(bad) > 0:
        raise ValueError('Missed a CASBAH QSO')

    # Fill up z, FUV, NUV
    new_keys = ['Z', 'FUV', 'NUV']
    for key in new_keys:
        val = milliq[key][idx]
        casbah.add_column( Column( val, name=key ))
    
    
    # Generate a file
    print('z1qso_analy.casbah: Writing CASBAH table -- {:s}'.format(casbah_fil))
    casbah.write(casbah_fil, overwrite=True)
    #xxf.table_to_fits( casbah, casbah_fil)

    # Return
    return casbah

    

####
def pair_up(targ, qsos, sep, z_targ=None, silent=False):
    """
    Given a list of targets, find all quasars within an angle or physical distance.

    Parameters:
    ----------
    targ: array of SkyCoord
    qsos: array of SkyCoord for quasars
    sep:  Angle or Quantity
       Separation to search for
    z_targ: array (None)
      Array of redshifts of the targets.  Required if sep=Distance
      Calculations are done in physical

    Returns:
    ----------
    idxt:    Indices of targets that were hits
    idxqso:  Indices of qsos that were hits
    d2d:     Angular or physical separation between qso-targ pairs
    """


    # Set angle
    if not z_targ is None: 
        from astropy.cosmology import WMAP9 as cosmo
        medz = np.median(z_targ)
        pkpc_amin = cosmo.kpc_comoving_per_arcmin( medz )/(1+medz) # physical kpc per arcmin
        ang_sep = 3.*(sep / pkpc_amin).to('arcsec')  # Physical
        if not silent:
            print('z1qa.pair_up: Searching to {:g}'.format(ang_sep))
    else:
        ang_sep = sep

    # Search about
    idxt, idxq, d2d, d3d = qsos.search_around_sky(targ, ang_sep)

    # Cut on physical distance?
    if not z_targ is None:
        pkpc_amin = cosmo.kpc_comoving_per_arcmin( z_targ[idxt] )/(1+z_targ[idxt]) # physical kpc per arcmin
        phys_sep = d2d.to('arcmin') * pkpc_amin
        gdp = np.where(phys_sep < sep)[0]
        # Cut
        if len(gdp) > 0:
            idxt = idxt[gdp]
            idxq = idxq[gdp]
            d2d = phys_sep[gdp]

    if not silent:
        print('z1qa.pair_up: Found {:d} matches!'.format(len(gdp)))

    # Return
    return idxt, idxq, d2d


####
def mk_z1qso_table(milliq=None, outfil=None, only_new=True):
    """
    Generate a binary FITS table of the new sources
    """
    # Outfil
    if outfil is None:
        outfil = 'z1qso_new.fits'


    # Load objects and candidates
    z1qso, cand = load_z1qso(good_z=True)

    # New?
    if only_new:
        # Load observed targets with JXP cuts
        if milliq is None:
            milliq = load_milliq()
        # Match to known quasars
        print('Cross matching')
        mq_rad = SkyCoord(ra=milliq['RA']*u.degree, dec=milliq['DEC']*u.degree)
        zq_rad = SkyCoord(ra=z1qso['RA']*u.degree, dec=z1qso['DEC']*u.degree)
        idx, d2d, d3d = coords.match_coordinates_sky(zq_rad, mq_rad, nthneighbor=1)

        # New ones
        newq = d2d.to('arcsec').value > 5. 
        z1qso = z1qso[newq]

    print('z1qa.mk_z1qso_table: Cleaning up..')

    # Rename columns, units, QTable
    z1qso = QTable(z1qso)
    #z1qso.rename_column('RAD', 'RA')
    #z1qso.rename_column('DECD', 'DEC')
    z1qso['RA'].unit = u.degree
    z1qso['DEC'].unit = u.degree

    # Clean up the names!
    for jj,iz1qso in enumerate(z1qso):
        pdb.set_trace()
        #z1qso['NAME'][jj] = xrad.dtos1( (iz1qso['RA'], iz1qso['DEC']), fmt=1)

    # Write
    print('z1qso_analy: Writing {:d} quasars to {:s} and compressing...'.format(len(z1qso),outfil))
    xxf.table_to_fits( z1qso, outfil, compress=True, comment='z1QSO New QSOs')


#### ########################## #########################
#### ########################## #########################
#### ########################## #########################

def main(flg):

    if flg == 'all':
        flg = np.sum( np.array( [2**ii for ii in range(1)] ))
    else:
        flg = int(flg)

    # Generate new QSO list
    if (flg % 2**1) >= 2**0:
        mk_z1qso_table()

    # Test target search
    if (flg % 2**2) >= 2**1:
        # BOSS LRGs
        boss = Table.read('Candidates/LRGs/wisconsin_pca_bc03-v5_7_0.fits.gz')
        uni_boss = Table.read('Candidates/LRGs/uniq_wisconsin_pca_bc03-v5_7_0.fits.gz')
        gdb = np.where(uni_boss['CHOOSE'] == 1)[0]
        boss = boss[gdb]
        c_boss = SkyCoord(ra=boss['RA']*u.degree, dec=boss['DEC']*u.degree)
        # z1QSO
        z1qso = new_z1qso()
        c_zq = SkyCoord(ra=z1qso['RAD']*u.degree, dec=z1qso['DECD']*u.degree)
        # Call
        sep = Angle(30.*u.arcsec)
        idx,idxq,d2d = pair_up(c_boss, c_zq, sep)
        #
        print('Matches: ',idx)


    
# Command line execution
if __name__ == '__main__':

    if len(sys.argv) == 1: #
        flg = 0 
        flg += 1  # Generate the "New" list
        #flg += 2**1  # Test target search with LRGs
    else:
        flg = sys.argv[1]

    main(flg)
