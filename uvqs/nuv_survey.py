# Module for generating NUV z1QSO candidates

# Imports
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import glob, os, sys

from astropy.io import ascii, fits
from astropy import units as u
from astropy.table import Table, Column, QTable
from astropy import coordinates as coords
from astropy.coordinates import SkyCoord

#from xastropy.xutils import xdebug as xdb
#from xastropy.obs import radec as xrad
#from xastropy.xutils import fits as xxf

#sys.path.append(os.path.abspath(os.getenv('z1Q')+"/Analysis/py"))
#import z1qso_analy as z1qa
#sys.path.append(os.path.abspath(os.getenv('z1Q')+"/Analysis/Candidates/py"))
#import z1qso_cand as z1qc

from uvqs import candidates as uvqs_cand
from uvqs import io as uvqs_io

#
def get_nuv_cand(outfil=None, known_fil=None):
    '''
    Generate a list of NUV AGN candidates from WISE + GALEX

    Paramaeters:
    -----------
    known_fil: string (None)
      Filename for a Table of NUV candidates with known redshifts
    '''
    #  Read photometry
    z1_galex_nuv, z1_wise_nuv = uvqs_cand.load_phot(use_NUV=True)

    # COLOR CUTS
    gdc = color_cut(z1_galex_nuv, z1_wise_nuv)
    srt = np.argsort(z1_wise_nuv[gdc]['RA'])
    gdc = gdc[srt]
    w1mw2 = z1_wise_nuv[gdc]['W1MPRO']-z1_wise_nuv[gdc]['W2MPRO']

    # Generate the QTable
    z1_cand_nuv = QTable( [Column( z1_wise_nuv[gdc]['RA'], name='RA', unit=u.degree)] )
    z1_cand_nuv.add_column( Column( z1_wise_nuv[gdc]['DEC'], name='DEC', unit=u.degree) )
    z1_cand_nuv.add_column( Column( z1_wise_nuv[gdc]['W1MPRO'], name='W1') )
    z1_cand_nuv.add_column( Column( z1_wise_nuv[gdc]['W2MPRO'], name='W2') )
    z1_cand_nuv.add_column( Column( z1_wise_nuv[gdc]['W1MPRO']-z1_wise_nuv[gdc]['W2MPRO'], name='W1-W2') )
    z1_cand_nuv.add_column( Column( z1_galex_nuv[gdc]['FUV'], name='FUV') )
    z1_cand_nuv.add_column( Column( z1_galex_nuv[gdc]['NUV'], name='NUV') )
    z1_cand_nuv.add_column( Column( z1_galex_nuv[gdc]['FUV']-z1_galex_nuv[gdc]['NUV'], name='FUV-NUV') )

    msk = z1_cand_nuv['RA'] == z1_cand_nuv['RA'] 
    c_z1_cand = SkyCoord(ra=z1_cand_nuv['RA'], dec=z1_cand_nuv['DEC'])

    # Check against z1QSO 
    cut_z1qso = uvqs_io.new_z1qso() # z1QSO with good z, not in milliq
    c_zq= SkyCoord(ra=cut_z1qso['RA'], dec=cut_z1qso['DEC'])
    idx, d2d, d3d = coords.match_coordinates_sky(c_z1_cand, c_zq, nthneighbor=1)
    bad = np.where( d2d.to('arcsec') < (5.*u.arcsec))[0]
    print('z1q_nuv: Rejecting {:d} quasars from z1QSO'.format(len(bad)))
    msk[bad] = False

    # Check against Flesch+15
    milliq = uvqs_io.load_milliq(good_z=True, good_FUV=True)
    c_mq = SkyCoord(ra=milliq['RA']*u.degree, dec=milliq['DEC']*u.degree)
    idx, d2d, d3d = coords.match_coordinates_sky(c_z1_cand, c_mq, nthneighbor=1)
    # Fill in z
    z = np.zeros(len(c_z1_cand))
    mt = np.where( (d2d.to('arcsec') < (5.*u.arcsec)) )[0]
    z[mt] = milliq['Z'][idx][mt]
    z1_cand_nuv.add_column( Column( z, name='Z'))

    # Write?
    if known_fil is not None:
        import pdb; pdb.set_trace()
        xxf.table_to_fits(Table(z1_cand_nuv[mt]), known_fil, compress=True, comment='Known NUV Candidates')
        print('z1qso_nuv.get_cand: Writing {:d} candidates to {:s}'.format(len(z1_cand_nuv), outfil))
        
    # Parse out low-z
    bad = np.where( (d2d.to('arcsec') < (5.*u.arcsec)) & (milliq['Z'][idx] < 0.5))[0]
    print('z1q_nuv: Rejecting {:d} quasars from Flesch+15'.format(len(bad)))
    msk[bad] = False

    # Add priority flags
    z1_cand_nuv = z1_cand_nuv[msk]
    priority = np.zeros(len(z1_cand_nuv),dtype='int')
    
    gd1 = np.where(z1_cand_nuv['NUV'] < 16.5)[0] 
    priority[gd1] = 1
    gd2 = np.where( (z1_cand_nuv['NUV'] < 17.0) & priority == 0)[0] 
    priority[gd2] = 2
    gd3 = np.where( z1_cand_nuv['NUV'] > 17.0 )[0]
    priority[gd3] = 3

    z1_cand_nuv.add_column( Column( priority, name='PRI') )

    # Fill in z1qso observed
    z1qso, cand = z1qa.load_z1qso() # z1QSO with good z, not in milliq
    c_zn= SkyCoord(ra=z1_cand_nuv['RA'], dec=z1_cand_nuv['DEC'])
    c_zq= SkyCoord(ra=z1qso['RA']*u.degree, dec=z1qso['DEC']*u.degree)
    idx, d2d, d3d = coords.match_coordinates_sky(c_zn, c_zq, nthneighbor=1)
    mt = np.where( (d2d.to('arcsec') < (5.*u.arcsec)) )[0]
    spec_qual = np.zeros(len(z1_cand_nuv),dtype='int')
    z_qual = np.zeros(len(z1_cand_nuv),dtype='int')

    spec_qual[mt] = z1qso['SPEC_QUAL'][idx[mt]]
    z_qual[mt] = z1qso['Z_QUAL'][idx[mt]]

    z1_cand_nuv.add_column( Column( spec_qual, name='SPEC_QUAL') )
    z1_cand_nuv.add_column( Column( z_qual, name='Z_QUAL') )

    # Names
    names = []
    for jj,z1c in enumerate(z1_cand_nuv):
        name = 'z1c_J{:s}{:s}'.format(c_zn[jj].ra.to_string(unit=u.hour,pad=True,sep='',precision=2), c_zn[jj].dec.to_string(pad=True,alwayssign=True,sep='',precision=2))
        names.append(str(name))
    z1_cand_nuv.add_column( Column( names, name='NAME') )

    # Sort
    srt = np.argsort(z1_cand_nuv['RA'])
    z1_cand_nuv = z1_cand_nuv[srt]

    # Write?
    if not outfil is None:
        xxf.table_to_fits(Table(z1_cand_nuv), outfil, compress=True, comment='z1QSO NUV Candidates')
        print('z1qso_nuv.get_cand: Writing {:d} candidates to {:s}'.format(len(z1_cand_nuv), outfil))

    # Return
    return z1_cand_nuv


####
#  Plot the standard candidates + modify those that have been observed
def sky_plot(i_cand, outfil=None, cand_only=False, cut_z1q=True):
    ''' 
    Plots a mollweide projection of the z1QSO NUV 
    '''

    # Imports
    import matplotlib as mpl
    mpl.rcParams['font.family'] = 'stixgeneral'
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib import pyplot as plt
    import matplotlib.gridspec as gridspec


    # Load observed targets
    #z1qso = Table( fits.open('z1qso_final_spec.fits.gz')[1].data )

    if cut_z1q:
        cut = np.where( (i_cand['Z_QUAL'] < 2) & (i_cand['SPEC_QUAL'] < 4))[0]
        NUV_cand = i_cand[cut]
    else:
        NUV_cand = i_cand

    # Init
    if outfil == None:
        outfil='fig_z1QSO_NUV_skyplot.pdf'
    i_z1qso = {'CUTS': [17.5, 18.0, 18.5], 'CLRS': ['blue', 'darkgreen', 'lightgrey'],
               'LBLS': ['NUV<16.5', '16.5<NUV<17', 'NUV>17'], 'PRI': [1,2,3]
               }


    # Start the plot
    if outfil != None:
        pp = PdfPages(outfil)

    plt.figure(figsize=(8, 5))
    plt.clf()

    ax = plt.axes(projection='mollweide')
    ax = plt.axes()
    ax.set_xlabel('RA')
    ax.set_ylabel('DEC')
    ax.set_xticklabels(np.arange(30,331,30))
    ax.grid(True)

    # Cuts
    allpri = i_z1qso['PRI']
    allc = i_z1qso['CLRS']
    alabel = i_z1qso['LBLS']
    psz = 3. # Scatter point size
    plw = 0.5 # Scatter point line-width

    # Candidates
    for ii,pri in enumerate(allpri):
        pidx  = np.where(NUV_cand['PRI'] == pri)[0]
        plt.scatter((NUV_cand['RA'][pidx].value-180.)*np.pi/180., NUV_cand['DEC'][pidx].value*np.pi/180.,
                    s=psz, lw=plw, edgecolors=allc[ii], facecolors='none', label=alabel[ii])
    
    '''
    # Objects
    if cand_only is False:
        alabel = i_z1qso['LBLS']
        # Obj with good spectrum
        for ii,pri in enumerate(allpri):
            pidx  = np.where((z1qso['PRI'] == pri) & (z1qso['SPEC_QUAL'] > 2))[0]
            plt.scatter((z1qso['RAD'][pidx]-180.)*np.pi/180., z1qso['DECD'][pidx]*np.pi/180.,
                        marker='+', facecolors=allc[ii], s=psz, lw=0.2)
    
        # Obj with good z
        for ii,pri in enumerate(allpri):
            pidx  = np.where((z1qso['PRI'] == pri) & (z1qso['Z_QUAL'] > 1))[0]
            plt.scatter((z1qso['RAD'][pidx]-180.)*np.pi/180., z1qso['DECD'][pidx]*np.pi/180.,
                        s=psz, lw=plw, edgecolors=allc[ii], facecolors=allc[ii], label=alabel[ii])
    '''
    
    # Legend
    legend = plt.legend(loc='upper right', scatterpoints=1, borderpad=0.2,
                        handletextpad=0.1, fontsize='small')

    # End
    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
    if outfil != None:
        pp.savefig()
        pp.close()
        print('Writing {:s}'.format(outfil))
    else: 
        plt.show()


#
def color_cut(wise, W1_W2=0.6):
    '''
    Performs the color cuts for the NUV selection

    Parameters:
    --------
    W1_W2: float (0.6)
      WISE color cut for AGNs (e.g. Stern+XX)
    '''

    #if not keyword_set(FUV_NUV) then  FUV_NUV = 0.3
    #if not keyword_set(MAX_FUV_NUV) then  MAX_FUV_NUV = 99.
    gd = np.where( (wise['RA'] > 0.) & (wise['W1SNR'] > 5.) &
                   (wise['W2SNR'] > 5.) & ( (wise['W1MPRO']-wise['W2MPRO']) > W1_W2))[0]
    ncand = len(gd)
    if ncand == 0:
        raise ValueError('Uh oh')
    # 
    print('color_cut: Found {:d} sources matching the cuts.'.format(len(gd)))
    return gd
