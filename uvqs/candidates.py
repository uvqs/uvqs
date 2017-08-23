# Module for Candidate code for UVQ project
# Imports
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import glob, os, sys

from pkg_resources import resource_filename

from astropy import units as u
from astropy.table import Table
from astropy import coordinates as coords
from astropy.coordinates import SkyCoord


def load_phot(use_NUV=False, verbose=True):
    """ Load the WISE and GALEX photometry
    Uses MIS data where applicable

    Parameters
    ----------
    use_NUV: bool, optional

    Returns
    -------
    z1_galex: Table
      GALEX photometry
    z1_wise: Table
      WISE-matched photometry
    """
    # WISE and GALEX AIS
    if use_NUV:
        galex_wise_fil = resource_filename('uvqs','/data/photom/GALEX_NODUP_NUV_WISE_2015apr06.fits')
    else:
        galex_wise_fil = resource_filename('uvqs','/data/photom/GALEX_NODUP_WISE_2014apr16.fits')
    z1_galex = Table.read(galex_wise_fil)
    z1_wise = Table.read(galex_wise_fil, exten=2)

    # Update for GALEX MIS photometry (might need to get NUV MIS)
    galex_MIS = resource_filename('uvqs','/data/photom/GALEX_MIS_xavier_2014apr18.fit.gz')
    mis = Table.read(galex_MIS)
    nmis = len(mis)

    z1_rad = SkyCoord(ra=z1_galex['RA']*u.degree, dec=z1_galex['DEC']*u.degree)
    mis_rad = SkyCoord(ra=mis['ra']*u.degree, dec=mis['dec']*u.degree)
    idx, d2d, d3d = coords.match_coordinates_sky(z1_rad, mis_rad, nthneighbor=1)

    mt = d2d.to('arcsec').value < 5. 
    z1_galex[mt]['NUV'] = mis[idx[mt]]['nuv']
    z1_galex[mt]['FUV'] = mis[idx[mt]]['fuv']

    if verbose:
        print('Updated {:d} sources with MIS'.format(len(mt)))

    # Return
    return z1_galex, z1_wise
