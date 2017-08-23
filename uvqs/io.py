""" I/O routines reltaed to UVQS
"""

from pkg_resources import resource_filename

from astropy.table import Table

def load_milliq(milliq_fil=None, good_z=False, good_FUV=False, good_NUV=False,
                verbose=True):
    """
    Read in QSOs from Flesch+15 (with JXP enhancements)

    Parameters
    ----------
    milliq_fil
    good_z: Bool (False)
      Require a good z [Should already have been done by JXP in his file]
    good_FUV: Bool (False)
      Require a good FUV value (>0.)
    good_NUV: Bool (False)
      Require a good NUV value (>0.)

    Returns
    -------
    milliq: Table

    """
    if milliq_fil is None:
        milliq_fil = resource_filename('uvqs', 'data/MilliQuas/xmq_GALEX_cutz.fits.gz')
    if verbose:
        print('z1qso_analy: Reading known quasars from {:s}'.format(milliq_fil))

    # Load
    milliq = Table.read(milliq_fil)
    msk = milliq['Z'] == milliq['Z']

    # Restrict to solid redshifts
    if good_z:
        for jj,item in enumerate(milliq['TYPE']):
            msk[jj] = ('Q' in item) | ('A' in item)

    # Restrict to sources with FUV magnitudes
    if good_FUV:
        bad = np.where(np.abs(milliq['FUV']) < 1.e-3)[0]
        if verbose:
            print('load_milliq: Cutting {:d} sources without a FUV measurement.'.format(len(bad)))
        msk[bad] = False
    # Restrict to sources with NUV magnitudes
    if good_NUV:
        bad = np.where(np.abs(milliq['NUV']) < 1.e-3)[0]
        if verbose:
            print('load_milliq: Cutting {:d} sources without a NUV measurement.'.format(len(bad)))
        msk[bad] = False

    # Final
    gd = np.where(msk)[0]
    if verbose:
        print('load_milliq: Return {:d} good sources.'.format(len(gd)))
    return milliq[msk]


def load_z1qso(z1qso_fil=None, cand_fil=None, good_z=False, new=None, NOSTOP=False):
    """
    DEPRECATED!?!
    Read in sources observed and update for priority
    good_z: Bool (False)
      If True, Require Z_QUAL > 1
    new: Table (None)
      Table of sources to match against.  Parse to the new ones
    """

    #if not NOSTOP:
    #    pdb.set_trace()  # You probably want the UVQ DR1

    if z1qso_fil is None:
        z1qso_fil = resource_filename('uvqs', 'data/UVQS/uvqs_dr1_sources.fits')

    print('z1qso_analy: Reading z1qso data from {:s}'.format(z1qso_fil))

    # Load candidates
    cand = get_candidates(cand_fil)

    # Load objects
    z1qso = Table( fits.open(z1qso_fil)[1].data )
    nz1 = len(z1qso)

    priority = np.zeros(nz1, dtype='int')  # 0 = Not a candidate
    priority[:] = -1

    # Fill up
    cqso = SkyCoord(ra=cand['RA']*u.degree, dec=cand['DEC']*u.degree)
    oqso = SkyCoord(ra=z1qso['RA']*u.degree, dec=z1qso['DEC']*u.degree)
    idx, d2d, d3d = coords.match_coordinates_sky(oqso, cqso, nthneighbor=1)

    # Close
    mtch = np.where( d2d.to('arcsec').value < 5. )[0]
    priority[mtch] = cand['PRI'][idx[mtch]]

    # Quick sanity check
    hi = np.where(priority == 1)[0]
    if np.max(z1qso['FUV'][hi]) > 17.5:
        raise ValueError('Oops..')

    #if np.max(d2d)*3600. > 1.*u.Unit('arcsec'):

    z1qso.add_column( Column( priority, name='PRI') )

    # Good z?
    if good_z:
        gdz = np.where(z1qso['Z_QUAL'] > 1)[0]
        z1qso = z1qso[gdz]

    # New?
    if not new is None:
        nq_rad = SkyCoord(ra=new['RA']*u.degree, dec=new['DEC']*u.degree)
        zq_rad = SkyCoord(ra=z1qso['RA']*u.degree, dec=z1qso['DEC']*u.degree)
        idx, d2d, d3d = coords.match_coordinates_sky(zq_rad, nq_rad, nthneighbor=1)

        newq = d2d.to('arcsec').value > 5.
        z1qso = z1qso[newq]
        print('load_z1qso: Keeping only the {:d} new ones.'.format(len(z1qso)))

    return z1qso, cand

####
def load_new_z1qso(z1qso_fil=None):
    """ Read in new z1QSO with good z
    """

    if z1qso_fil is None:
        z1qso_fil = resource_filename('uvqs', 'data/UVQS/z1qso_new.fits.gz')

    print('z1qso_analy: Reading z1qso data from {:s}'.format(z1qso_fil))
    z1qso = Table.read(z1qso_fil)

    return z1qso

def load_uvq_dr1(verbose=True):
    """
    Load the UVQ DR1 dataset

    Parameters:
    ----------
    verbose: bool, optional

    Returns
    -------
    uvq: Table
    """
    uvq_fil=resource_filename('uvqs', 'data/UVQS/uvq_dr1.fits.gz')
    if verbose:
        print('Reading {:s}'.format(uvq_fil))
    uvq = Table.read(uvq_fil)

    #Return
    return uvq
