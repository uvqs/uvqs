#!/usr/bin/env python
"""
Grab files for Photometry from Dropbox
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import os
import pdb

try:  # Python 3
    ustr = unicode
except NameError:
    ustr = str

def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(description='Grab the Photometry for UVQS')
    #parser.add_argument("-v", "--version", default='v02', help="DB version to generate")
    #parser.add_argument("-llist", default='ISM', action='store_true', help="Name of LineList:  ISM, HI, H2, CO, etc.")

    if options is None:
        pargs = parser.parse_args()
    else:
        pargs = parser.parse_args(options)
    return pargs


def main(pargs):
    """ Run
    Parameters
    ----------
    pargs

    Returns
    -------

    """
    import subprocess
    from pkg_resources import resource_filename

    # URLs
    urls = []
    newfiles = []
    # MIS
    urls.append('https://www.dropbox.com/s/zq52f606afcgmlk/GALEX_MIS_xavier_2014apr18.fit.gz?dl=0')
    newfiles.append(resource_filename('uvqs', 'data/UVQS/GALEX_MIS_xavier_2014apr18.fit.gz'))
    # Original WISE (FUV)
    urls.append('https://www.dropbox.com/s/y9ejblkf3zwgiw6/GALEX_NODUP_WISE_2014apr16.fits.gz?dl=0')
    newfiles.append(resource_filename('uvqs', 'data/UVQS/GALEX_NODUP_WISE_2014apr16.fits.gz'))
    # NUV WISE
    urls.append('https://www.dropbox.com/s/1m38btzj1aole2i/GALEX_NODUP_NUV_WISE_2015apr06.fits.gz?dl=0')
    newfiles.append(resource_filename('uvqs', 'data/UVQS/GALEX_NODUP_NUV_WISE_2015apr06.fits.gz'))
    # NUV GALEX
    urls.append('https://www.dropbox.com/s/3j3th9fslz633bu/GALEX_NODUP_NUV_xavier_2015apr06.fits.gz?dl=0')
    newfiles.append(resource_filename('uvqs', 'data/UVQS/GALEX_NODUP_NUV_xavier_2015apr06.fits.gz'))
    # Original GALEX (FUV)
    urls.append('https://www.dropbox.com/s/quybab7xzo57lcv/GALEX_NODUP_xavier_2014apr16.fits.gz?dl=0')
    newfiles.append(resource_filename('uvqs', 'data/UVQS/GALEX_NODUP_xavier_2014apr16.fits.gz'))

    # wget command
    for ss,url in enumerate(urls):
        subprocess.call(['wget', '--continue', '--timestamping', url, '-O', newfiles[ss]])

