#!/usr/bin/env python

import sys, pyfits, os.path, pylab, string, re, time
from glob import glob
from copy import copy
import scipy
from scipy import interpolate, optimize

c = 299792458e10 

''' add SeqNr column to FITS table '''
def add_SeqNr(file,extension=1):
    p = pyfits.open(file)
    ext_str = True 
    try: 
        extension = int(extension)
        ext_str = False
    except: pass

    cols = []
    for col in p[extension].columns:
        cols.append(col)

    cols.append(pyfits.Column(name='SeqNr',format='L',array=range(len(p[extension].data))))

    hdu = pyfits.PrimaryHDU()
    hdulist = pyfits.HDUList([hdu])
    tbhu = pyfits.new_table(cols)
    hdulist.append(tbhu)
    if ext_str:
        hdulist[1].header.update('EXTNAME',extension)
    outcat = file.replace('.fits','.seqnr.fits')
    os.system('rm ' + outcat)
    hdulist.writeto(outcat)
    print 'WRITTEN TO ', outcat
    

def readtxtfile(file):
    f = open(file,'r').readlines()
    file_out = []
    for l in f:
        res = re.split('\s+',l)
        if l[0] != '#':
            if res[0] == '': res = res[1:]           
            if res[-1] == '': res = res[:-1]
            file_out.append([(x) for x in res])
    return file_out


def synth(p,spectra,filters,show=False):

    mags = {} 

    for filt in filters:
        specall = scipy.zeros(len(spectra[0][0][:,1]))
        val = 0 
        for coeff,specfull  in [[p[0],spectra[0]]]: #,[p[1],spectra[1]],[1.-p[0]-p[1],spectra[2]]]: 
            spec = specfull[0]
            specStep = spec[1:,0] - spec[0:-1,0] # wavelength increment                   
            #print specStep[400:600], 'specStep'
            resampFilter = filt['spline'](spec[:,0]) # define an interpolating function
            #print resampFilter
            #print filt_name

            if False: #string.find(filt_name,'SDSS') != -1:
                pylab.plot(spec[:,0],resampFilter) 
                pylab.show()
            ''' need to multiply by polynomial '''
            val += abs(coeff)*sum(specStep * resampFilter[:-1] * spec[:-1,0] * spec[:-1,1]) # photon counting!!

        logEff = scipy.log10(val)                                        
        logNorm = scipy.log10(sum(resampFilter[:-1]*c*specStep/spec[:-1,0]))
        mag = 2.5*(logNorm - logEff) # to calculated an AB magnitude
    
        mags[filt['name']]=mag

    if show:
        pylab.plot(spec[:,0], specall)
        pylab.show()
    return mags

def cas_locus(fits=True):

    if fits:
        import pyfits
        locus_list_mag = pyfits.open(os.environ['BIGMACS'] + '/lociCAS.fits')['STDTAB']
    else:
        import pickle
        f = open(os.environ['BIGMACS'] + 'lociCAS','r')
        m = pickle.Unpickler(f)
        locus_list_mag = m.load()

    #print locus_list_mag


    return locus_list_mag


def synthesize_expected_locus_for_observations(filters):
    c_locus = cas_locus()    

    ''' add SDSS filters '''        
    SDSS_filters = [{'name':'USDSS','filter':'SDSS-u.res'},{'name':'GSDSS','filter':'SDSS-g.res'},{'name':'RSDSS','filter':'SDSS-r.res'},{'name':'ISDSS','filter':'SDSS-i.res'},{'name':'ZSDSS','filter':'SDSS-z.res'}]
    filter_info = get_filters([[a['name'],a['filter']] for a in SDSS_filters])
    for i in range(len(filter_info)):
        SDSS_filters[i].update(filter_info[i])

    loci = [a[:-1] for a in open(os.environ['BIGMACS'] + '/LOCUS_SPECTRA/spliced_spectra','r').readlines()]

    locus = []

    print 'STARTING SAMPLING LOCUS'

    for i in range(len(loci[:])):
        locus_point = loci[i]
        locus_index = int(locus_point.replace('.dat',''))
        print 'CONVOLVING RESPONSE FUNCTIONS WITH SPECTRUM ' + str(locus_point)
        stitchSpec = scipy.genfromtxt(os.environ['BIGMACS'] + '/LOCUS_SPECTRA/' + locus_point)

        ''' do not synthesize 2MASS filters '''
        mags = synth([1.,0,0,0],[[stitchSpec]],filter(lambda x: string.find(x['filter'],'2MASS') == -1, filters + SDSS_filters)) 

        ''' not synthesizing the 2MASS magnitudes '''
        for filt in filters:
            if filt['filter'] == 'J2MASS.res':

                mags[filt['mag']] = (mags['ZSDSS'] - c_locus.data.field('ZSDSS_JTMASS')[locus_index]) #+ (mags['ISDSS'] - c_locus['ZSDSS_JTMASS'][2*i]))/2.


                #print mags[filt['mag']], mags['ZSDSS'], c_locus.data.field('ZSDSS_JTMASS')
                #raw_input()

        locus.append(mags)

    print 'FINISHED SAMPLING LOCUS'
    return locus


def parse_columns(columns_description, fitSDSS=False, noHoldExceptSDSS=False, noHoldExcept2MASS=False):
    f = filter(lambda x: x[0] != '#', filter(lambda x: len(x) > 0, readtxtfile(columns_description)))

    input_info = [] 
    for l in f:
        dict = {'mag_err': l[1], 'filter': l[2]}

        if len(l[0].split('#')) > 1:                
            dict['mag'] = l[0].split('#')[0]
            dict['plotName'] = l[0].split('#')[1]
        else:
            dict['mag'] = dict['plotName'] = l[0]

        ''' do not hold any filter fixed if noHold is True '''
        if noHoldExceptSDSS and (dict['filter'][:4] != 'SDSS' and dict['filter'][:4] != 'sdss'):
            dict['HOLD_VARY'] = 'VARY'
        elif noHoldExcept2MASS and (dict['filter'][:4] != '2MASS'):
            dict['HOLD_VARY'] = 'VARY'
        else:
            dict['HOLD_VARY'] = l[3] 
            if dict['HOLD_VARY'] == 'HOLD':
                dict['ZP'] = float(l[4])
        input_info.append(dict)        

    return input_info 
    

def get_filters(flist = [['USDSS','SDSS-u.res'],['GSDSS','SDSS-g.res'],['RSDSS','SDSS-r.res'],['ISDSS','SDSS-i.res'],['ZSDSS','SDSS-z.res']]):

    filt_dir = os.environ['BIGMACS'] + '/FILTERS/'

    #flist = [{'mag':'USDSS','filter':'SDSS-u.res'},{'mag':'GSDSS','filter':'SDSS-g.res'},{'mag':'RSDSS','filter':'SDSS-r.res'},{'mag':'ISDSS','filter':'SDSS-i.res'},{'mag':'ZSDSS','filter':'SDSS-z.res'}]

    filters = []
    for filt_name, filt_file in flist:
        file = filt_dir + filt_file
        filt = scipy.loadtxt(file)
        step = filt[1,0] - filt[0,0]
        if filt[0,0] > filt[-1,0]:
            filt_list = filt.tolist()
            filt_list.reverse()
            filt = scipy.array(filt_list)

        filterSpline = scipy.interpolate.interp1d(filt[:,0], filt[:,1], 
                                       bounds_error = False, 
                                       fill_value = 0.)
        filters.append({'wavelength':filt[:,0],'response':filt[:,1],'spline':copy(filterSpline),'step':copy(step),'name':copy(filt_name),'center wavelength': scipy.average(filt[:,0],weights=filt[:,1])})

    return filters


def odonnell(input,wavelength=True):
    if wavelength:
        x = (input/ 10**4.)**-1. # convert from angstroms to micrometers
    else:
        x = input
    y = (x-1.82)
    a = 1. + 0.104*y - 0.609*y**2. + 0.701*y**3. + 1.137*y**4. - 1.718*y**5. - 0.827*y**6. + 1.647*y**7. - 0.505*y**8.
    b = 1.952*y + 2.908*y**2. - 3.989*y**3. - 7.985*y**4. + 11.102*y**5. + 5.491*y**6. - 10.805*y**7. + 3.347*y**8.
    R_v = 3.1
    A = a + b/R_v
    return A

''' compute filter extinction coefficients from O'Donnell extinction law '''
def compute_ext(filt, N=0.78):

    print '''COMPUTING EXTINCTION COEFFICIENT USING FITZPATRICK99 EXTINCTION LAW'''
    odonnell_ext_1_um = 1.32 # A(1 um) / E(B-V) where R_v = 3.1
    import scipy, math
    sed = scipy.loadtxt(os.environ['BIGMACS'] + '/munari.sed')
    ''' now stitch together with blackbody spectrum '''    
    longwave = sed[-1,0]    
    flux = sed[-1,1]

    wavelength = sed[:,0].tolist()
    source = sed[:,1].tolist()

    a = flux / longwave**-3.

    for wave in scipy.arange(11500,20000,25):
        wavelength.append(wave)
        source.append(a*wave**-3.)

    import scipy
    from scipy import interpolate
    sedSpline = interpolate.interp1d(wavelength, source, 
                                   bounds_error = True, 
                                   )


    #s_od = N*odonnell_ext_1_um*odonnell(scipy.arange(3000,20000))
    #s_od = fitzpatrick(scipy.arange(3000,20000))
    #import pylab
    #pylab.clf()
    #pylab.plot(scipy.arange(3000,20000),s_od)
    #pylab.xlim([3000,20000])
    #pylab.show()

    ''' source flux is ergs / s / Ang '''
    filt_wavelength = filt['wavelength']
    filt_response = filt['response']
    throw_out = scipy.zeros(len(filt_wavelength))

    ''' trim off zero-valued tails of response function'''
    for i in range(len(filt_wavelength)):
        if filt_response[i] == 0:
            throw_out[i]= 1.
        else: break

    for i in range(len(filt_wavelength)):
        if filt_response[len(filt_response)-1-i] == 0:
            throw_out[len(filt_response)-1-i]= 1.
        else: break

    filt_wavelength = filt_wavelength[throw_out==0.] 
    filt_response = filt_response[throw_out==0.] 

    #print scipy.array([(filt_wavelength[i]) for i in range(len(filt_wavelength[:-1]))])
    #print scipy.array([fitzpatrick(filt_wavelength[i]) for i in range(len(filt_wavelength[:-1]))])
    numerator = scipy.array([10.**(fitzpatrick(filt_wavelength[i])/-2.5)*sedSpline(filt_wavelength[i])*filt_wavelength[i]*(filt_response[i])*(filt_wavelength[i+1]-filt_wavelength[i]) for i in range(len(filt_wavelength[:-1]))])
    denom = scipy.array([source[i]*filt_wavelength[i]*(filt_response[i])*(filt_wavelength[i+1]-filt_wavelength[i]) for i in range(len(filt_wavelength[:-1]))])

    coeff = -2.5*math.log10(numerator.sum()/denom.sum()) 

    print filt['name'], coeff, 'coeff'
    return coeff


''' returns A(lambda)/A(1 um) '''
def fitzpatrick(input,wavelength=True,plot=False):
    if wavelength:
        x = (input/ 10**4.)**-1. # convert from angstroms to micrometers
    else:
        x = input

    ''' for R_v = 3.1 '''
    wavelength = [0.000,0.377,0.820,1.667,1.828,2.141,2.433,3.704,3.846]
    ratio = [0.000,0.265,0.829,2.688,3.055,3.806,4.315,6.265,6.591]

    #wavelength = [1.667,1.828,2.141,2.433]        
    #ratio = [2.688,3.055,3.806,4.315]

    import scipy
    from scipy import interpolate
    fitzSpline = scipy.interpolate.interp1d(wavelength,ratio,
                                   kind='cubic') #, bounds_error=False)
    if plot:
        import pylab, scipy
        pylab.clf()
        x_range = scipy.arange(wavelength[0],wavelength[-1],0.01)
        pylab.plot(x_range,fitzSpline(x_range))
        pylab.scatter(wavelength,ratio)
        pylab.box()
        pylab.xlim([0,4])
        pylab.savefig('/Users/pkelly/Dropbox/spline.png')


    ''' normalized so that A_1 um = 1 mag '''
    A = fitzSpline(x) / fitzSpline(1.)

    return  A 

