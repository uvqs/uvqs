
# Module for creating the database of z1QSO

# Imports
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import time, subprocess
import glob, os, sys

from pkg_resources import resource_filename

from scipy.interpolate import splev, splrep # bspline

from astropy.io import ascii, fits
from astropy import units as u
from astropy.table import Table, Column
from astropy import coordinates as coords
from astropy.coordinates import SkyCoord
from astropy.coordinates.angles import Angle
from astropy.nddata import StdDevUncertainty

from linetools.spectra import io as lsio

#from xastropy.xutils import xdebug as xdb
#from xastropy.obs import radec as xrad
#from xastropy.obs import finder as xfinder
#from xastropy.xutils import fits as xxf

# def load_uvq_dr1
# def uvq_name
# def uvq_path
# def scan_runs
# def parse_zass
# def mk_images
# def build_1dspec
# def pdf_1dspec
# def bspline_ymnx
# def get_1dspec
# def targ_dict
# def flux
# def build_full_database

#sys.path.append(os.path.abspath(os.getenv('z1Q')+"/Analysis/py"))
#import z1qso_analy as z1qa
#sys.path.append(os.path.abspath(os.getenv('z1Q')+"/Analysis/Candidates/py"))
#import z1qso_cand as z1qc



####
def uvq_name(coord, latex=False):
    """
    Generate the UVQ name

    Parameters:
    ----------
    coord: SkyCoord
    latex: bool, optional
      Format for LateX?  (False)
    """
    if latex:
        name = 'UVQSJ{:s}${:s}$'.format(coord.ra.to_string(unit=u.hour,pad=True,sep='',precision=2),
                                            coord.dec.to_string(pad=True,alwayssign=True,sep='',precision=1))
    else:
        name = 'UVQSJ{:s}{:s}'.format(coord.ra.to_string(unit=u.hour,pad=True,sep='',precision=2),
                                            coord.dec.to_string(pad=True,alwayssign=True,sep='',precision=1))
    #Return
    return name

####
def uvq_path(ptype):
    '''Generate the path
    '''
    if ptype == 'finders':
        db_path=os.getenv('UVQ_DATA')+'/Images/'
    elif ptype == '2dspec':
        db_path=os.getenv('UVQ_DATA')+'/2D_Spectra/'
    elif ptype == '1dspec':
        db_path=os.getenv('UVQ_DATA')+'/1D_Spectra/'
    else:
        raise ValueError('Not a valid path')
    # Return
    return db_path

####
def scan_runs():
    '''
    Load the runs for z1Q
    '''
    z1q_runs = ascii.read(os.getenv('z1Q')+'/Database/z1qso_obs_runs.dat',comment='#')

    return z1q_runs

####
def parse_zass(z1_galex=None, z1_wise=None, outfil=None, verbose=False):
    '''
    Loops through the 2D zassess files and generate a UVQ table
    '''
    # Load Photometry
    if (z1_galex is None) | (z1_wise is None):
        z1_galex, z1_wise = z1qc.load_phot()#use_NUV=True)
    c_wise = SkyCoord(ra=z1_wise['RA']*u.degree, dec=z1_wise['DEC']*u.degree)

    # Load runs
    z1q_runs = scan_runs()

    # Read other targets (not UVQ)
    other = QTable.read(os.getenv('z1Q')+'/Database/other_targ.lst',comment='#', format='ascii')
    xrad.stod_table(other)
    c_other = SkyCoord( ra=other['RA'], dec=other['DEC'])

    # Start the Total dict
    adict = [targ_dict()]
    all_qso = QTable(rows=adict)
    #xdb.set_trace()

    for irun in z1q_runs:
        if verbose:
            print('Working on run at {:s} file from {:s} on {:s}'.format(irun['zass_fil'], irun['obs'], irun['date']))
        #if irun['obs'] != 'LCO':
        #    continue

        # Read zass file
        zfil = os.getenv('DROPBOX_DIR')+'z1QSO/database/2D_Spectra/'+irun['dbx_pth']+'/'+irun['zass_fil']
        zass = xxf.bintab_to_table(zfil)
        #xdb.set_trace()

        # Loop on targets
        for targ in zass:
            # sci_file?
            if len(targ['SCI_FILE'].strip()) == 0:
                continue

            # Check RA/DEC against complete list
            mt = np.where( (np.abs(targ['RAD']-all_qso['RA']) < 1e-4) &
                           (np.abs(targ['DECD']-all_qso['DEC']) < 1e-4) )[0]
            if len(mt) == 0:
                pass
            elif len(mt) == 1:
                '''
                #Check night/date [NOT SURE WHY I DID THIS PREVIOUSLY]
                if ((all_qso['OBSV'][mt] != irun['obs']) | (all_qso['OBS_DATE'][mt] != irun['date'])):
                '''
                # Differing: Take new one if  higher z_qual
                # or equal but better spectral quality
                #if (np.abs(targ['RAD']-217.007833) < 1e-4):
                #    xdb.set_trace() 
                if ((targ['Z_QUAL'] > all_qso[mt]['Z_QUAL']) |
                    ((targ['Z_QUAL'] == all_qso[mt]['Z_QUAL']) &
                    (targ['SPEC_QUAL'] > all_qso[mt]['SPEC_QUAL']))):
                    print('Replacing {:s} at {:s} on {:s}'.format(all_qso[mt[0]]['NAME'],
                                                                  all_qso[mt[0]]['OBSV'],
                                                                  all_qso[mt[0]]['OBS_DATE']))
                    # Fill in new values
                    for key in targ.dtype.names:
                        if key in all_qso.keys():
                            all_qso[mt[0]][key] = targ[key]
                    all_qso[mt[0]]['RA'] = targ['RAD']
                    all_qso[mt[0]]['DEC'] = targ['DECD']
                    all_qso[mt[0]]['OBSV'] = irun['obs']
                    all_qso[mt[0]]['OBS_DATE'] = irun['date']
                    all_qso[mt[0]]['instr'] = irun['instr']
                    all_qso[mt[0]]['dbx_path'] = irun['dbx_pth']
                    all_qso[mt[0]]['ZASS_FIL'] = irun['zass_fil']
                                #xdb.set_trace()
            else:
                xdb.set_trace()
                raise ValueError('Uh oh')
            if len(mt) == 1:
                continue

            #  Search for a match in WISE + GALEX
            coord = SkyCoord(ra=targ['RAD']*u.degree, dec=targ['DECD']*u.degree)
            sep = coord.separation(c_wise)
            wise_mt = np.where( sep.to('arcsec').value < 2.)[0]
            noth = 0
            if len(wise_mt) != 1:
                # Other target??
                osep = SkyCoord(ra=targ['RAD']*u.degree, dec=targ['DECD']*u.degree).separation(c_other)
                #xdb.set_trace()
                if np.min(osep) > (1.*u.arcsec):
                    print(targ)
                    raise ValueError('Might be a new NUV candidate.  Need the other WISE file!')
                oth_mt = np.argmin(osep)
                noth=1
            if noth == 1:
                if verbose:
                    print('Skipping other target {:s}'.format(other[oth_mt]['Name']))
                continue

            galex_mt = wise_mt #; These 2 files are synced

            #; Fill it all up
            tmp = targ_dict()
            for key in targ.dtype.names:
                if key in tmp.keys():
                    #xdb.set_trace()
                    tmp[key] = targ[key]
            tmp['RA'] = targ['RAD']
            tmp['DEC'] = targ['DECD']
            tmp['FUV'] = z1_galex[galex_mt]['FUV']
            tmp['NUV'] = z1_galex[galex_mt]['NUV']
            tmp['W1'] = z1_wise[wise_mt]['W1MPRO']
            tmp['W2'] = z1_wise[wise_mt]['W2MPRO']
            tmp['NAME'] = uvq_name(coord)

            # Run
            tmp['OBSV'] = irun['obs']
            tmp['OBS_DATE'] = irun['date']
            tmp['instr'] = irun['instr']
            tmp['dbx_path'] = irun['dbx_pth']
            tmp['ZASS_FIL'] = irun['zass_fil']

            #; Save
            #xdb.set_trace()
            all_qso.add_row(tmp)

    # Sort
    all_qso.sort('RA')
    all_qso = all_qso[1:]

    # Units
    all_qso['RA'].unit = 'degree'
    all_qso['DEC'].unit = 'degree'

    # Primary candidate?
    flg_c = np.ones(len(all_qso),dtype='int')+1
    gdc = np.where( all_qso['FUV']-all_qso['NUV'] > 0.3)[0]
    flg_c[gdc] = 1 # Primary
    all_qso.add_column(Column(flg_c, name='CAND'))

    return all_qso

def write_uvq(uvq, outfil=None):

    # Write
    if outfil is None:
        outfil = 'uvq_dr1_sources.fits'
    date = 'Date: '+time.strftime('%d-%b-%Y')
    xxf.table_to_fits(uvq, outfil, comment=date)
    print('-------------------------------')
    print('Wrote table to {:s}'.format(outfil))

####
def mk_images(uvq, clobber=False, USE_PDF=False):
    '''
    Generate the Finder images for the sample

    Parameters:
    -----------
    uvq: QTable
      Summary file of the UVQ sources
    '''
    reload(xfinder)

    # Outpath
    db_path=uvq_path('finders') 

    # Loop to loop
    finders = []
    if USE_PDF:
        exten = '.pdf'
    else:
        exten = '.png'
    for row in uvq:
        # Outfil
        coord = SkyCoord(ra=row['RA'], dec=row['DEC'])

        name = uvq_name(coord)
        outroot = name+'_'
        fil = glob.glob(db_path+outroot+'*'+exten)
        if (len(fil) > 0) and (not clobber):
            print('File already generated: {:s}'.format(fil[0]))
            i0 = fil[0].rfind('/')
            finders.append(fil[0][i0+1:])
            continue

        # Generate
        title = name 
        if USE_PDF:
            ftype = '.pdf'
            out_type = 'PDF'
        else:
            ftype = '.png'
            out_type = 'PNG'
        fil = name+ftype
        #xdb.set_trace()
        oBW = xfinder.main([title, coord.ra.to_string(unit=u.hour,pad=True,sep=':',precision=2),
                            coord.dec.to_string(pad=True,alwayssign=True,sep=':',precision=1)], 
                            radec=1, imsize=1.0*u.arcmin, OUT_TYPE=out_type)

        # Move to Folder
        if oBW == 1:
            suff = 'DSS'
        else:
            suff = 'SDSS'
        finders.append(outroot+suff+ftype)
        subprocess.call( ['mv', fil, db_path+outroot+suff+ftype] )

    # Return
    print('Finished with the finders!')
    return finders

########
def build_1d_spec(uvq, clobber=False, USE_PDF=False):
    '''
    Generate a 1D spec file

    Parameters:
    -----------
    uvq: QTable
      Summary file of the UVQ sources
    '''

    datapath=uvq_path('2dspec') 
    db_path=uvq_path('1dspec') 

    # Geneate flux dict
    flux_dict ={}
    runs = ['Lick_22Mar2014', 'Lick_22Mar2014r', 
        'Lick_30Mar2014', 'Lick_30Mar2014r', 
        'MMT_Apr2014', 
        'Lick_May2014', 'Lick_May2014r', 
        'Lick_Aug2014', 'Lick_Aug2014r', 
        'LCO_Aug2014', 'LCO_Feb2015', 'MMT_Jan2015',
        'Lick_Jan2015', 'Lick_Jan2015r']
    for run in runs:
        flux_dict[run] = {}
    flux_fils = [
        'Lick/March22_2014/Flux/BD+33_b600.fits', 
        'Lick/March22_2014/Flux/BD+33_r831.fits', 
        'Lick/March30_2014/Flux/Feige34_b600.fits', 
        'Lick/March30_2014/Flux/Feige34_r600.fits', 
        'MMT/BCS_Apr2014/g191b2b_300sens.fits', 
        'Lick/May_2014/Flux/HZ44_b600.fits', 
        'Lick/May_2014/Flux/HZ44_r600.fits', 
        'Lick/May_2014/Flux/HZ44_b600.fits', 
        'Lick/May_2014/Flux/HZ44_r600.fits',  # Aug2014
        'LCO/August2014/04aug2014/Flux/feige110_sens.fits', 
        'LCO/August2014/04aug2014/Flux/feige110_sens.fits',
        'MMT/BCS_Apr2014/g191b2b_300sens.fits', 
        'Lick/Jan_2015/Flux/G191b2b_b600.fits', 
        'Lick/Jan_2015/Flux/G191b2b_r600.fits' 
        ]

    for kk,flux_fil in enumerate(flux_fils):
        # Read Flux file(s)
        flux_dict[runs[kk]]['FLUX_FIL'] = flux_fil
        mag_set = xxf.bintab_to_table(datapath+flux_fil)
        flux_dict[runs[kk]]['MAG_SET'] = mag_set

        # bspline
        knots = np.array(mag_set['FULLBKPT']).flatten()
        coeff = np.array(mag_set['COEFF']).flatten()
        nord = mag_set['NORD'][0]-1  # Offset of 1 in the defintion
        tck = (knots, coeff, nord)

        # Save
        flux_dict[runs[kk]]['TCK'] = tck

    # Loop on uvq
    spec_files = []
    for row in uvq:
        if row['OBSV'] not in ['Lick','LCO', 'CAHA', 'Magellan', 'Keck', 'MMT']:
            raise ValueError('Not ready for this Observatory')

        # Outfil
        coord = SkyCoord(ra=row['RA'], dec=row['DEC'])
        name = uvq_name(coord)
        outroot = name+'_'
        outfil = db_path+outroot+row['OBSV']+'.fits'
        spec_files.append(outroot+row['OBSV']+'.fits')

        # Clobber?
        if os.path.isfile(outfil) and (not clobber):
            print('Not over-writing 1D Spec file: {:s}'.format(outfil))
            continue

        # ###############
        # Read unfluxed data
        sci_file = datapath + row['dbx_path'] + '/' + row['SCI_FILE']
        if row['instr'].strip() in ['Kast','BCS']:
            # 2D output
            try:
                hdu = fits.open(sci_file.strip())
            except IOError:
                xdb.set_trace()
                raise IOError('File does not exist: {:s}'.format(sci_file))
            # Grab spectrum
            obj = row['OBJ_IDX']
            try:
                bint = Table( Table(hdu[5].data)[obj] ) # Yes, this is crazy looking
            except IndexError:
                print('Obj = {:d}'.format(obj))
                xdb.set_trace()
            if 'JUNK' in bint.colnames: # There is at least one JUNK file on blue side
                spec = None
            else:
                spec = lsio.readspec(bint)
            # 2nd camera??
            if row['instr'].strip() in ['Kast']:
                scnd_file = datapath + row['dbx_path'] + '/' + row['SCND_FILE']
                obj2 = row['OBJ_IDX2']
                # Read
                hdu2 = fits.open(scnd_file.strip())
                rint = Table( Table(hdu2[5].data)[obj2] ) 
                #
                rspec = lsio.readspec(rint)
                rspec.head = hdu2[0].header
            # Header
            try:
                spec.head = hdu[0].header
            except AttributeError:
                pass
        elif row['instr'].strip() in ['CAFOS', 'MagE', 'MBC', 'ESI']:
            spec = lsio.readspec(sci_file)

        #print(row)
        # #######
        # Flux 
        flux_key = row['OBSV']+'_'+row['OBS_DATE']
        try:
            mag_set = flux_dict[flux_key]['MAG_SET']
        except KeyError:
            if flux_key not in ['CAHA_Aug2014', 'Magellan_Jul2014', 'Magellan_Aug2014', 'Keck_Apr2014']:
                xdb.set_trace()
        else:
            #
            tck = flux_dict[flux_key]['TCK']
            if spec is not None:
                spec = flux(spec,mag_set,tck)
 
            # Red camera?
            if row['instr'].strip() in ['Kast']:
                flux_keyr = row['OBSV']+'_'+row['OBS_DATE']+'r'
                mag_set = flux_dict[flux_keyr]['MAG_SET']
                tck = flux_dict[flux_keyr]['TCK']
                # Flux
                rspec = flux(rspec,mag_set,tck)
                # Combine the spectra
                if spec is None:
                    spec = rspec
                else:    
                    spec = spec.splice(rspec)
                    spec.head = hdu[0].header
            # Header
            spec.head.set('FLUX_FIL', 
                value=str(flux_dict[flux_key]['FLUX_FIL']))

        # History
        spec.head.add_history('UVQ Survey DR1 (v1.0)')
        spec.head.add_history('Relative fluxing only.')  
        spec.head.add_history('Vacuum, heliocentric wavelengths')

        # Write
        spec.write_to_fits(outfil)

        # ####
        # Generate the PDF
        if USE_PDF:
            ftype = '.pdf'
        else:
            ftype = '.png'
        pdf_fil = db_path+outroot+row['OBSV']+ftype
        title = name 
        pdf_1dspec(row, spec, pdf_fil, title=title)

    # Finish
    print('Finished with the 1D spectra!')
    return spec_files


########
def pdf_1dspec(uvq, spec, pdf_fil, title=None, ymnx=None):
    '''
    Generate a 1D spec file

    Parameters:
    -----------
    uvq: QTable row
      From UVQ summary
    spec: Spectrum1D
      Spectrum to plot
    pdf_fil: string
      Name of the PDF file
    title: string
      Title for the plot
    '''
    # Import
    import matplotlib as mpl
    mpl.rcParams['font.family'] = 'stixgeneral'
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib import pyplot as plt
    import matplotlib.gridspec as gridspec
    from xastropy.plotting import utils as xputils

    # b-spline for limits
    gd = np.where( (spec.sig > 0.) & (spec.flux > 0.))[0]
    try:
        xmnx = [np.min(spec.wavelength.value[gd]), np.max(spec.wavelength.value[gd])]
    except ValueError:
        if uvq['SCI_FILE'] in ['August2014/CAHA_none.fits']: # None file for CAHA stars
            ymnx = np.array([-0.05,1.])
    # Special
    if uvq['OBSV'] == 'CAHA':
        xmnx = (3800., 9000.)
        gd = np.where( (spec.sig > 0.) & (spec.flux > 0.) & 
            (spec.wavelength > xmnx[0]*u.AA))[0]
    elif uvq['instr'] == 'MagE':
        xmnx = (3300., 9000.)
        gd = np.where( (spec.sig > 0.) & (spec.flux > 0.) & 
            (spec.wavelength > xmnx[0]*u.AA))[0]
    elif uvq['instr'] == 'Kast':
        xmnx[1] = 7500. 
    elif uvq['instr'] == 'ESI':
        xmnx = (4200., 8000.)
        gd = np.where( (spec.sig > 0.) & (spec.flux > 0.) & 
            (spec.wavelength > xmnx[0]*u.AA))[0]
    if ymnx is None:
        ymnx = bspline_ymnx(spec,gd)

    yscl = 10.**(int(np.log10(ymnx[1])))

    # Start the plot
    #pp = PdfPages(pdf_fil)

    plt.figure(figsize=(8, 5))
    plt.clf()
    gs = gridspec.GridSpec(1, 1)
    ax = plt.subplot(gs[0,0])

    # Axes
    #ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
    #ax.xaxis.set_major_locator(plt.MultipleLocator(1.))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    #ax.yaxis.set_major_locator(plt.MultipleLocator(1.))
    ax.set_xlim(xmnx)
    ax.set_ylim(ymnx/yscl)

    ax.set_xlabel(r'Wavelength ($\AA$)')
    ax.set_ylabel('Relative Flux')

    # Plot
    ax.plot(spec.wavelength, spec.flux/yscl, 'k-',drawstyle='steps-mid', linewidth=0.9)
    ax.plot(spec.wavelength, spec.sig/yscl, 'r:')

    # Title
    y0 = 0.92
    if title is not None:
        ax.text(0.5, y0, title, transform=ax.transAxes, size=18, color='blue',
               ha='center', bbox={'facecolor':'white', 'edgecolor':'white'})
    # Other data
    x1_lbl = 0.94
    yoff = 0.05
    ax.text(x1_lbl, y0, uvq['OBSV'], transform=ax.transAxes, size=14, color='black', ha='right')
    ax.text(x1_lbl, y0-yoff, uvq['instr'], transform=ax.transAxes, size=14, color='black', ha='right')
    ax.text(x1_lbl, y0-2*yoff, uvq['OBS_DATE'], transform=ax.transAxes, size=14, color='black',
        ha='right')
    ax.text(x1_lbl, y0-3*yoff, 'z={:0.3f}'.format(uvq['Z']), transform=ax.transAxes, size=14,
        color='black', ha='right')
    ax.text(x1_lbl, y0-4*yoff, r'$z_Q=$'+'{:d}'.format(uvq['Z_QUAL']), transform=ax.transAxes,
        size=14, color='black', ha='right')
    ax.text(x1_lbl, y0-5*yoff, r'$S_Q=$'+'{:d}'.format(uvq['SPEC_QUAL']), transform=ax.transAxes,
        size=14, color='black', ha='right')

    # Template
    #xdb.set_trace()
    zfil = os.getenv('UVQ_DATA')+'2D_Spectra/'+uvq['dbx_path']+'/'+uvq['ZASS_FIL']
    try:
        new_eigen = fits.open(zfil)[2].data
    except:
        pass
    else:
        rebin_wv2 = fits.open(zfil)[3].data
        sz_eigen = new_eigen.shape
        model_flux = np.dot(new_eigen.transpose(), uvq['THETA'][0:sz_eigen[0]])
        ax.plot(rebin_wv2*(1+uvq['Z']), model_flux/yscl, 'g-', linewidth=0.5, alpha=0.6)
        #xdb.set_trace()

    # Label a few lines 
    line_lbls = {1903.91: 'CIII]', 4341.69: r'H$\gamma$', 4862.70: r'H$\beta$', 
        5008.24: '[OIII]', 2799.40: 'MgII', 6564.63: r'H$\alpha$'}
    for key in line_lbls.keys():
        xplt = (1+uvq['Z'])*key
        if (xplt>xmnx[0]) & (xplt<xmnx[1]):
            ax.text((1+uvq['Z'])*key, 0.8*ymnx[1]/yscl, line_lbls[key], size=11, color='green', ha='center')

    # Font size
    xputils.set_fontsize(ax,15.)

    # Layout and save
    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
    plt.savefig(pdf_fil)
    #pp.savefig(bbox_inches='tight')
    #pp.close()
    plt.close()
    print('Wrote spectrum to {:s}'.format(pdf_fil))

########
def bspline_ymnx(spec,gd=None, frac=1.3):
    '''Figure out ymnx using a b-spline approach.

    Parameters:
    -----------
    spec: Spectrum1D
    frac: float
      Fraction above the max to plot
    '''
    if gd is None:
        gd = np.where( (spec.sig > 0.) & (spec.flux > 0.))[0]
    #
    ngd = len(gd)
    idx_knots = np.arange(10, ngd-10, 20) # A knot every good 20 pixels
    knots = spec.wavelength.value[gd[idx_knots]]
    try:
        tck = splrep( spec.wavelength.value[gd], spec.flux[gd], w=1./spec.sig[gd], k=3, t=knots)
    except ValueError:
        # Could be Nebular emission (not well-modelded by a b-spline)
        isrt = np.argsort(spec.flux[gd])
        ymnx = np.array([0.0, spec.flux[gd[isrt[-4]]]*frac]) # Hope for the best..
        # 
        return ymnx
        #
        #print('pdf_1dspec: Probably bad knot spacing') # https://github.com/scipy/scipy/issues/2343
        #xdb.set_trace()

    # Continuing with the b-spline
    bs = splev(spec.wavelength, tck,ext=1)

    bmax = np.max(bs)
    ymnx = np.array([0.0*bmax, bmax*frac])
    # Return
    return ymnx

########
def get_1dspec(name):
    ''' Grab the 1D spectrum from the UVQ database

    Parameters:
    -----------
    name: string 
      UVQ name
    '''
    db_path = uvq_path('1dspec') 
    fils = glob.glob(db_path+name+'*.fits*')
    # More than one?
    if len(fils) > 1:
        raise ValueError('Multiple files: Not sure what to do here')
    elif len(fils) == 0:
        xdb.set_trace()
        raise ValueError('No UVQ file with the name: {:s}'.format(name))
    else:
        return lsio.readspec(fils[0])

################################ ################################
def targ_dict():
    '''
    Generate the dict for a UVQ target
    '''
    dums = str(' '*80)
    tmp =  {'NAME': dums,  # JXXXXXX.XX+XXXXXX.X
            'RA': 0.,
            'DEC': 0.,
            'FUV': 0., #  GALEX photometry
            'NUV': 0., #
            'W1': 0., #  ;: WISE photometry
            'W2': 0., #
            'OBSV': dums, # ; Observatory
            'OBS_DATE': dums, # ; Date of Observation
            'instr': dums, # ; Instrument
            'dbx_path': dums, #
            'SCI_FILE': dums, #
            'SCND_FILE': dums, #
            'comment': dums, #  ;; User comment
            'ZASS_FIL': dums, #
            'THETA': np.zeros(20), #
            'SPEC_QUAL': 0, #  ;; 1=Poorest, 5=Best
            'OBJ': 0, # Object Flag
            'OBJ_IDX': 0, # Object Flag
            'OBJ_IDX2': 0, # Object Flag
            'Z_QUAL': 0, #  ;; 1=Poorest, 5=Best
            'Z': 0., #
            'Z_SIG': 0. #  ;; Usually from SDSS z code
            }
    return tmp

def flux(spec,mag_set,tck):
    '''Flux
    '''
    npix = len(spec.wavelength)
    wave_min = mag_set['WAVE_MIN'] * u.AA
    wave_max = mag_set['WAVE_MAX'] * u.AA

    ext = np.ones(npix)
    MAX_EXTRAP = 0. * u.AA
    inds = np.where( (spec.wavelength >= (wave_min - MAX_EXTRAP)) &
                 (spec.wavelength <= (wave_max + MAX_EXTRAP)))[0]
    # Bspline
    try:
        mag_func = splev(spec.wavelength[inds], tck)
    except ValueError:
        xdb.set_trace()
    sens = 10.0**(0.4*mag_func)
    scale = np.zeros(npix)
    scale[inds] = sens*ext[inds]

    # Apply
    spec.flux = spec.flux * scale
    spec.uncertainty = StdDevUncertainty(spec.sig*scale)

    # Return
    return spec


################################ ################################
def build_full_database(skip_uvq=False,USE_PDF=False):
    '''
    Script to generate the full UVQ database (DR1)
    '''
    # Generate the Table
    if skip_uvq: # Mainly for debugging
        uvq = QTable.read('uvq_dr1_sources.fits')
        uvq.remove_column('SPEC_FILE')
        uvq.remove_column('FINDER')
    else:
        uvq = parse_zass()

    # Finder Images
    finders = mk_images(uvq, USE_PDF=USE_PDF)

    # 1D Spectra
    spec_files = build_1d_spec(uvq, USE_PDF=USE_PDF)

    # Add to uvq
    uvq.add_column(Column(spec_files, name='SPEC_FILE'))
    uvq.add_column(Column(finders, name='FINDER'))

    # Write the Table
    write_uvq(uvq)
    print('Remember to push the output file to the Database!!!!!')


#### ########################## #########################
#### ########################## #########################
#### ########################## #########################

def main(flg):

    if flg == 'all':
        flg = np.sum( np.array( [2**ii for ii in range(1)] ))
    else:
        flg = int(flg)

    # Test Runs
    if (flg % 2**1) >= 2**0:
        uvq_runs = scan_runs()
        print(uvq_runs)

    # Make the Binary table 
    if (flg % 2**2) >= 2**1:
        xdb.set_trace() # Should not run this any longer. Only *full* database
        uvq=parse_zass()

    # Make the Images
    if (flg % 2**3) >= 2**2:
        xdb.set_trace() # Should not run this any longer. Only *full* database
        uvq = QTable.read('uvq_dr1_sources.fits')
        mk_images(uvq)

    # Test 1Dspec
    if (flg % 2**4) >= 2**3:
        uvq = QTable.read('uvq_dr1_sources.fits')
        #idx = np.where(uvq['SCI_FILE'] == '14feb2015/Science/sci-ccd0033.fits.gz')[0]
        #idx = np.where(uvq['SCI_FILE'] == 'Jan_2015/Science/sci-b333.fits.gz')[0]
        #idx = np.where(uvq['SCI_FILE'] == 'MagE_Aug2014/J023507-040205_F.fits')[0]
        #idx = np.where(uvq['SCI_FILE'] == 'August2014/J180437.45+542540.7.fits')[0]
        idx = np.where(uvq['SCI_FILE'] == '10aug2014/Science/sci-ccd1023.fits.gz')[0]
        #idx = np.where(uvq['SCI_FILE'] == 'ESI_Apr2014/J180650.67+69a_F.fits')[0]
        #idx = np.where(uvq['SCI_FILE'] == 'March30_2014/Science/sci-b228.fits.gz')[0]
        #idx = np.where(uvq['SCI_FILE'] == 'March22_2014/Science/sci-b2111.fits.gz')[0]
        #idx = np.where(uvq['SCI_FILE'] == 'May_2014/Science/sci-b1071.fits.gz')[0]
        #idx = np.where(uvq['SCI_FILE'] == 'BCS_Apr2014/J1353+6345_mmt300.fits.gz')[0]
        #idx = np.where(uvq['SCI_FILE'] == '09aug2014/Science/sci-ccd0856.fits.gz')[0] # special case
        if len(idx) == 0:
            xdb.set_trace()
        spec_files = build_1d_spec(uvq[idx], clobber=True)

    # Build 1Dspec
    if (flg % 2**5) >= 2**4:
        xdb.set_trace() # Should not run this any longer. Only *full* database
        uvq = QTable.read('uvq_dr1_sources.fits')
        spec_files = build_1d_spec(uvq, clobber=False)

    # Build full database
    if (flg % 2**11) >= 2**10: # 1024
        build_full_database(skip_uvq=True)


# Command line execution
if __name__ == '__main__':

    if len(sys.argv) == 1: #
        flg = 0
        #flg += 2**0  # Test Observing Runs
        #flg += 2**1  # Generate the UVQ Table
        #flg += 2**2  # Make the images
        #flg += 2**3  # Test 1D spec
        flg += 2**4  # Make 1D spec
        #flg += 2**10  # Build full database (1024)
    else:
        flg = sys.argv[1]

    main(flg)
