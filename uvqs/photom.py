# Module for Photometry code related to the UVQ project
# Imports
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import glob, os, sys
import pdb

from pkg_resources import resource_filename

from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord, match_coordinates_sky


def galex_photom(coords):
    """ Match an input set of coordiantes to the GALEX photometry

    Parameters
    ----------
    coords : SkyCoord

    Returns
    -------
    fuv: ndarray of FUV values;  -1 for non-matched sources
    nuv: ndarray of NUV values;  -1 for non-matched sources
    """
    # Init
    fuv = -1 * np.ones(len(coords))
    nuv = -1 * np.ones(len(coords))

    # GALEX files
    galex_ais_fil = resource_filename('uvqs', 'data/UVQS/GALEX_NODUP_xavier_2014apr16.fits.gz')
    galex_mis_fil = resource_filename('uvqs', 'data/UVQS/GALEX_MIS_xavier_2014apr18.fit.gz')
    galex_nuv_fil = resource_filename('uvqs', 'data/UVQS/GALEX_NODUP_NUV_xavier_2015apr06.fits.gz')
    galex_files = [galex_ais_fil, galex_mis_fil, galex_nuv_fil]

    # Loop on em
    for gfile in galex_files:
        # Read
        try:
            galex_tbl = Table.read(gfile)
        except FileNotFoundError:
            raise FileNotFoundError("You probably need to run uvqs_get_photom!")
        # Coordinates
        try:
            galex_coord = SkyCoord(ra=galex_tbl['RA'], dec=galex_tbl['DEC'], unit='deg')
        except KeyError:
            galex_coord = SkyCoord(ra=galex_tbl['ra'], dec=galex_tbl['dec'], unit='deg')
        # Match
        idx, d2d, d3d = match_coordinates_sky(coords, galex_coord, nthneighbor=1)
        # Fill
        mtch = d2d.to('arcsec') < 5.*u.arcsec
        try:
            fuv[mtch] = galex_tbl['FUV'][idx[mtch]]
        except KeyError:
            fuv[mtch] = galex_tbl['fuv'][idx[mtch]]
            nuv[mtch] = galex_tbl['nuv'][idx[mtch]]
        else:
            nuv[mtch] = galex_tbl['NUV'][idx[mtch]]
    # Return
    return fuv, nuv
