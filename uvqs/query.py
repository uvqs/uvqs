# Module for UVQS queries

# Imports
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import time, subprocess
import glob, os, sys
import argparse

from astropy import units as u
from astropy.table import Table, Column, QTable, vstack
from astropy.coordinates import SkyCoord

from linetools import utils as ltu

from uvqs import database as uvqs_db
from uvqs import io as uvqs_io

#from xastropy.xutils import xdebug as xdb
#from xastropy.obs import radec as xrad
#from xastropy.obs import finder as xfinder

# def qck_srch

#sys.path.append(os.path.abspath(os.getenv('z1Q')+"/Analysis/py"))
#import z1qso_analy as z1qa
#sys.path.append(os.path.abspath(os.getenv('z1Q')+"/Analysis/Candidates/py"))
#import z1qso_cand as z1qc
#sys.path.append(os.path.abspath(os.getenv('z1Q')+"/Database/py"))
#import z1qso_database as z1qd


#
def radial_srch(icoord, radius=10*u.arcmin, uvq_only=False, allq=None, silent=False, max_src=None):
    """Perform a simple search of the UVQ database at an input coordinates within a given radius

    Parameters
    ----------
    icoord: tuple or coord or str
      Coordinates of the source
    radius: Quantity, optional
      Angular radius to search
    uvq_only: bool, optional
      Exclude Million quasars?  [False]
    allq: Table
      Table of all sources to search
    max_src: int, optional
      Maximum number of sources to return

    Returns
    -------
    cands: Table
       Candidates matching the criteria
    allq: Table
       Full table searched
    """

    # Convert to SkyCoord
    coord = ltu.radec_to_coord(icoord)

    # Tables
    if allq is None:
        # Load UVQ
        uvq1 = uvqs_io.load_uvq_dr1()

        # Milliq?
        if not uvq_only:
            milliq = uvqs_io.load_milliq(good_z=True)
            allq = vstack([Table(uvq1),Table(milliq)])
        else:
            allq = uvq1

    # Coords
    q_coord = SkyCoord(ra=allq['RA'], dec=allq['DEC'])

    # Search 
    sep = coord.separation(q_coord)
    gd = np.where( sep.to('arcmin') < radius)[0]

    # Return
    if len(gd) > 0:
        # Calculate separation
        if max_src is not None:
            nkeep = min(max_src,len(gd))
            isrt = np.argsort(sep[gd])
            gd = gd[isrt[np.arange(nkeep)]]
        # Generate
        tmp = allq[gd]
        # Add Separation
        sepc = Column(sep[gd].to('arcmin'), name='SEP')
        tmp.add_column(sepc)
        if not silent:
            print('Here are the quasars within your radius {:g}'.format(radius))
            print(tmp[['SEP', 'NAME','RA','DEC','FUV','NUV','Z']])
        return tmp, allq
    else:
        if not silent:
            print('I searched UVQ and MilliQ at {:s}'.format(coord))
            print('There are no quasars in your search radius {:g}. Sorry'.format(radius))

#
def uvq_inspect(icoord, show_finder=True, show_spec_pdf=True, silent=False): 
    """Show data for an input UVQ source
    NEEDS TO LINK TO SPECDB

    Parameters
    ----------
    icoord: tuple or coord or str
      Coordinates of the source

    Returns
    -------
    src : Table
      sources matching
    finder : str
      path+file of the finder chart in the DR1 database
    spec_fil : str
      path+file of the 1D spectrum
    """
    # Get source
    src = uvq_srch(icoord,radius=1.*u.arcsec, silent=silent)
    if src is None:
        return
    row = src[0]

    # Refine coord
    coord = SkyCoord(ra=row['RA'], dec=row['DEC'])

    # Open Finder
    finder_path = z1qd.uvq_path('finders')
    finder= finder_path+row['FINDER']
    if show_finder:
        if os.path.isfile(finder):
            subprocess.call( ['open', finder] )
        else:
            print('PDF finder: {:s} not found'.format(finder))

    # Spectrum 
    spec_path=z1qd.uvq_path('1dspec')
    spec_fil = spec_path+row['SPEC_FILE']
    if show_spec_pdf:
        # Search for PDF
        pdf_fil = spec_fil[0:spec_fil.rfind('.fits')]+'.pdf'
        if os.path.isfile(pdf_fil):
            subprocess.call( ['open', pdf_fil] )
        else:
            print('PDF for spec file: {:s} not found'.format(pdf_fil))

    # Return
    return src, finder, spec_fil


def uvq_srch(icoord, radius=1*u.arcsec, silent=True):
    '''Perform a quick search of the UVQ database at an input coordinates within a given radius
    Returns only the closest source within the search radius

    Parameters:
    -----------
    icoord: tuple or coord or str
      Coordinates of the source
    radius: Quantity, optional
      Angular radius to search
    silent: bool, optional
      Keep quiet?

    Returns:
    -----------
    '''

    # Use radial search  
    tmp,_ = radial_srch(icoord, max_src=1, radius=radius, uvq_only=True, silent=silent)
    if tmp is None:
        print('No UVQ source found within {:g}'.format(radius))
        return

    # Print
    print('======= Here is your source in UVQ =======')
    keys = tmp.keys()
    keys.sort()
    # Put NAME, RA, DEC first
    for reord in ['DEC', 'RA', 'NAME']:
        keys.remove(reord)
        keys = [reord] + keys
    for key in keys:
        print('{:s}: '.format(key), tmp[key][0])

    #if not silent:
    #    print(tmp[['NAME','Z','FUV','NUV','SCI_FILE']])

    # Return
    return tmp


#### ########################## #########################
#### ########################## #########################
#### ########################## #########################

def main():
    # MAKE THIS A SCRIPT, AS DESIRED

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('flag', type = str, help = '1=Radial, 2=UVQ Database')
    parser.add_argument('coord', type = str, help = 'JXXXXXX.XX+XXXXXX.X')

    args = parser.parse_args()

    # Option
    if args.flag=='1': # Radial search
        tmp, allq = radial_srch(args.coord)
    elif args.flag=='2': # UVQ search
        tmp = uvq_srch(args.coord)


# Command line execution
if __name__ == '__main__':

    if len(sys.argv) == 1: # Testing
        #main(['dum', '1', 'J010018.69-741815.9']) # Radial search
        #main(['dum', '2', 'J010018.69-741815.9'])  # UVQ search
        main(['dum', '2', 'J195151.72-054816.6'])
    else:
        main()

