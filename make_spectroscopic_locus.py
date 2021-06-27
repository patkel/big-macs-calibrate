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
    

''' search SDSS online database for spectra with similar synthetic and photometric colors to locus point '''
def get_sdss_spectra(umg, imz, gmr,rmi,number=4,tol=0.01,S_N=5):   

    def mag(c):
        return '(-2.5*LOG10(abs(s.spectrosynflux_' + c + ')) + 22.5)'

    dict_keys = ['s.plate', 's.MJD', 's.fiberID', 's.ra', 's.dec','s.survey'] +  [mag(c)  + ' as mag_' + c  for c in ['g','r','i']]

    dict_names = ['plate', 'MJD', 'fiberID', 'ra', 'dec','survey'] +  ['mag_' + c for c in ['g','r','i']]

    import sqlcl
    if rmi < 0.7: pattern = 's.elodiesptype like "%V%" and elodieSpType not like "%var%"'
    #elif 0.7 < rmi < 1.0: pattern = '(zbelodiesptype like "%G%v%" or zbelodiesptype like "%K%v%" or zbelodiesptype like "%M%v%")'
    else: pattern = 's.elodiesptype like "%M%V%"' # and s.elodiesptype not like "%var%"'

    ''' do not use stars strongly affected by the u-band red leak '''
    ''' http://www.sdss3.org/dr8/imaging/caveats.php '''
    ''' not synthetic magnitude!!! also affected by red leak '''    
    if rmi < 1.5:
        u_band_selection = ' and abs(' + str(umg) + ' - (p.psfMag_u - p.psfMag_g)) < 0.15'
    else:
        u_band_selection = ''

    query = 'select top ' + str(number) + ' ' + reduce(lambda x,y: x + ',' + y, ['' + x for x in dict_keys]) + ' from specobjall as s join specphoto as p on s.specobjid = p.specobjid join sppParams sp on sp.specobjid = s.specobjid where extinction_g < 0.05 and  ' + pattern +  ' and abs(' + mag('g') + ' - ' + mag('r') + '  - ' + str(gmr) + ') < ' + str(tol) + ' and abs(' + mag('r') + ' - ' + mag('i') + ' - ' + str(rmi) + ') < ' + str(tol) + ' and abs(' + mag('g') + ' - ' + mag('i') + ' - ' + str(gmr + rmi) + ') < ' + str(tol) + ' and s.snMedian > ' + str(S_N) + ' and abs(' + mag('g') + '  - ' + mag('r') + ' - (p.fibermag_g - p.fibermag_r)) < 0.1 and abs(' + mag('r') + ' - ' + mag('i') + ' - (p.fibermag_r - p.fibermag_i)) < 0.1 and abs(' + str(imz) + ' - (p.psfMag_i - p.psfMag_z)) < 0.5 ' + u_band_selection + ' \
order by -1.*s.snMedian'

    query = 'select top ' + str(number) + ' ' + reduce(lambda x,y: x + ',' + y, ['' + x for x in dict_keys]) + ' from specobjall as s join specphoto as p on s.specobjid = p.specobjid join sppParams sp on sp.specobjid = s.specobjid where extinction_g < 0.1 and s.spectrosynflux_g > 0 and  s.spectrosynflux_r > 0 and  s.spectrosynflux_i > 0 and \
 ' + pattern +  ' and abs(' + mag('g') + ' - ' + mag('r') + '  - ' + str(gmr) + ') < ' + str(tol) + ' and abs(' + mag('r') + ' - ' + mag('i') + ' - ' + str(rmi) + ') < ' + str(tol) + ' and abs(' + mag('g') + ' - ' + mag('i') + ' - ' + str(gmr + rmi) + ') < ' + str(tol) + ' and s.snMedian > ' + str(S_N) + ' \
and abs(' + str(imz) + ' - (p.psfMag_i - p.psfMag_z)) < 0.05 \
and abs(' + str(umg) + ' - (p.psfMag_u - p.psfMag_g)) < 0.05 \
order by -1.*s.snMedian'

    ''' cannot query SDSS database more than once per second '''
    time.sleep(1.5)

    print 'DOWNLOADING SDSS SPECTRA W/ SIMILAR COLORS'
    print query
    lines = sqlcl.query(query).readlines()
    print len(lines) - 1, 'STAR(S) FOUND'
    print lines

    dicts = []                                                                                   

    if lines[0] != 'N':

        for line in lines[1:]:
            dict = {}
            line = line.replace('\n','')
            res = re.split(',',line)
            for i in range(len(res)): 
                if dict_names[i] == 'fiberID' or dict_names[i] == 'plate' or dict_names[i] == 'MJD':
                    dict[dict_names[i]] = int(res[i])
                else:
                    dict[dict_names[i]] = (res[i])
            dicts.append(dict)
                                                                                                     

    return dicts

def download_sdss_spectrum(dict,plot=False): 
    dict['gmr'] = float(dict['mag_g']) - float(dict['mag_r'])
    dict['rmi'] = float(dict['mag_r']) - float(dict['mag_i'])
    file = "http://das.sdss.org/spectro/1d_26/%(plate)04d/1d/spSpec-%(MJD)d-%(plate)04d-%(fiberID)03d.fit" % dict      

    file = "http://data.sdss3.org/sas/dr8/spectro/1d_26/%(plate)04d/1d/spSpec-%(MJD)d-%(plate)04d-%(fiberID)03d.fit" % dict      
    #file = 'http://data.sdss3.org/returnSpec/fits?plateid=%(plate)04d&mjd=%(MJD)d&fiber=%(fiberID)03d' % dict

    #file = 'http://skyserver.sdss3.org/dr9/en/get/specByPF.asp?P=%(plate)04d&F=%(fiberID)03d' % dict

    bases = {'boss': ['data.sdss3.org/sas/dr9/boss/spectro/redux/v5_4_45/spectra/'],
            'sdss': ['data.sdss3.org/sas/dr9/sdss/spectro/redux/26/spectra/'],
            #'segue1': 'data.sdss3.org/sas/dr9/sdss/spectro/redux/26/spectra/',
            'segue1': ['data.sdss3.org/sas/dr9/sdss/spectro/redux/103/spectra/','data.sdss3.org/sas/dr9/sdss/spectro/redux/26/spectra/'],
            'segue2': ['data.sdss3.org/sas/dr9/sdss/spectro/redux/104/spectra/']}

    for root in bases[dict['survey']]:
        file = 'http://' + root + '%(plate)04d/spec-%(plate)04d-%(MJD)d-%(fiberID)04d.fits' % dict
        print dict['survey'], bases[dict['survey']]
        print file
    
        try: 
            p = pyfits.open(file)
            break
        except:
            print 'keep trying'


    print p, p[1].data, 'here'
    mask = p[1].data.field('and_mask')
    flux = p[1].data.field('flux')            
    print flux, mask
    indices = scipy.array(range(len(flux)))
    print indices

    COEFF0 = float(p[0].header['COEFF0'])
    COEFF1 = float(p[0].header['COEFF1'])
    print COEFF0, COEFF1
    wavelength = 10.**(COEFF0 + COEFF1*indices)
    spectrum = []
    for i in range(len(indices)):
        spectrum.append([wavelength[i],flux[i]])
    spectrum = scipy.array(spectrum)
    if plot:
        pylab.plot(spectrum[:,0], spectrum[:,1])
        pylab.xlabel('angstroms')
        pylab.ylabel('flux')
        pylab.show()
    return spectrum 

def build_new_spectrum(locus_index,plot_live=False, plot_directory=None, spectra_complete=None):

    filters = get_filters()
    locus_list = cas_locus()

    comp_list =  filter(lambda x: string.find(x.replace('SDSS_',''),'SDSS')!=-1 and string.find(x,'SDSS_')!=-1, locus_list.keys())
    
    gmr_all = locus_list['GSDSS_RSDSS'][:]
    rmi_all = locus_list['RSDSS_ISDSS'][:]
    umg_all = locus_list['USDSS_GSDSS'][:]
    imz_all = locus_list['ISDSS_ZSDSS'][:]

    print 'locus_index', locus_index
    gmr = locus_list['GSDSS_RSDSS'][locus_index]
    rmi = locus_list['RSDSS_ISDSS'][locus_index]
    umg = locus_list['USDSS_GSDSS'][locus_index]
    imz = locus_list['ISDSS_ZSDSS'][locus_index]

    #print gmr, rmi  

    if plot:
        pylab.clf()                         
        pylab.scatter(gmr_all,rmi_all,color='blue')
        pylab.scatter(gmr,rmi,color='red')
        pylab.title(str(locus_index))
        pylab.draw()


    good = False
    gmr_off = 0
    rmi_off = 0
    trys = 0
    tol = 0.01
    while not good:
        trys += 1
        #print gmr, rmi                                            
        dicts = get_sdss_spectra(umg,imz,gmr-gmr_off,rmi-rmi_off,tol=tol)

        if len(dicts):
            #print dicts                                                
            gmr_diffs = []
            rmi_diffs = []
            for dict in dicts:
                if 1: #try:
                    spectrum = download_sdss_spectrum(dict,plot=False)    
                    mags = synth([1.],[[spectrum]],filters,show=False)
                    #print mags
                    gmr_diffs.append(mags['GSDSS'] - mags['RSDSS'] - gmr)
                    rmi_diffs.append(mags['RSDSS'] - mags['ISDSS'] - rmi)
                #except: 
                #    print 'DOWNLOADING PROBLEM/POSSIBLY CORRUPTED FILE'
                #print mags['GSDSS'] - mags['RSDSS'], gmr
                #print float(dict['mag_0']) - float(dict['mag_1'])
                #print mags['RSDSS'] - mags['ISDSS'], rmi
                #print float(dict['mag_1']) - float(dict['mag_2'])

            gmr_diffs.sort()
            rmi_diffs.sort()
                                                                      
            median_gmr = gmr_diffs[int(len(gmr_diffs)/2)]
            median_rmi = rmi_diffs[int(len(rmi_diffs)/2)]
                                                                       
            if abs(median_gmr) > tol or abs(median_rmi) > tol:
                gmr_off += median_gmr
                rmi_off += median_rmi
            else: good = True            

            ''' increase tolerance if no star found, to break loop '''
            tol += 0.005
           
            #print gmr_off, rmi_off                                                           
            #print gmr_diffs, rmi_diffs
            #print median_gmr, median_rmi
            #print gmr, rmi
        else: tol += 0.01

        if not good: print 'NO SPECTRUM MEETS CRITERIA, RELAXING MATCHING TOLERANCE\nRECENTERING SEARCH'              
   
    #print spectrum 
    #print comp_list

    if False: #plot:
        max = spectrum[:,1].max()               
        pylab.plot(spectrum[:,0],spectrum[:,1]/max)
        pylab.xlim(3000,11000)
        pylab.show()

    sdssSpec, pickleSpec = find_most_similar_pickles_spectrum(spectrum, spectra_complete)

    stitchSpec = splice_spectra_together(sdssSpec,pickleSpec, locus_index,plot_live=plot_live, plot_directory=plot_directory)

    return stitchSpec 

''' assemble a new locus '''
def make_spectroscopic_model_locus(start_over=True):
    
    locus_list = cas_locus()

    keys = locus_list.keys()

    #filters = get_filters(sdss=False)

    spectra_dir = os.environ['BIGMACSMAKE'] + '/LOCUS_SPECTRA/'
    spectra_list_file = spectra_dir + 'spliced_spectra'

    if start_over:
        print  spectra_dir
        os.system('rm -rf ' + spectra_dir)
        os.system('mkdir -p ' + spectra_dir)

    spectra_complete = load_pickles_spectra() 

    #print len(locus_list['GSDSS_ISDSS'])

    for i in 2*scipy.array(range(len(locus_list['GSDSS_ISDSS'])/2)):

        if glob(spectra_list_file):
            spectra_files = [a[:-1] for a in open(spectra_list_file,'r').readlines()]
        else: spectra_files = []

        if not filter(lambda x: x==str(i)+'.dat', spectra_files): 

            stitchSpec = build_new_spectrum(i,plot_live=True,plot_directory=spectra_dir,spectra_complete=spectra_complete)                                  
            scipy.savetxt(spectra_dir + str(i) + '.dat',stitchSpec)                
            open(spectra_list_file,'a').write(str(i) + '.dat\n')


            
def splice_spectra_together(specSDSS,pickleSpec,locus_index,plot_live=False, plot_directory=None):
    filters = get_filters()
    locus_list = cas_locus()

    print 'SPLICING TOGETHER SDSS AND PICKLES SPECTRA'

    comp_list =  filter(lambda x: string.find(x.replace('SDSS_',''),'SDSS')!=-1 and string.find(x,'SDSS_')!=-1, locus_list.keys())
    #print comp_list
    grdiff = (locus_list['GSDSS_RSDSS'][locus_index])
    
    sdssSpline = scipy.interpolate.interp1d(specSDSS[:,0], specSDSS[:,1], 
                    bounds_error = False, 
                    fill_value = 0.)

    sdssLimits = [4200,8500] #specSDSS[-1,0]]
    zOverLap = [8000,9000]
    uOverLap = [4100,4600]

    #print sdssLimits
    specSDSS_new = []
    for l in specSDSS:
        if sdssLimits[0] < l[0] < sdssLimits[1]:
            specSDSS_new.append(l)

    specSDSS = scipy.array(specSDSS_new)
    
    uSpec = []
    zSpec = []
    zOverLapData = []
    uOverLapData = []
    for l in pickleSpec:
        if l[0] < sdssLimits[0]:
            uSpec.append(l)
        if l[0] > sdssLimits[1]:
            zSpec.append(l)
        if zOverLap[0] < l[0] < zOverLap[1]:
            zOverLapData.append(l)
        if uOverLap[0] < l[0] < uOverLap[1]:
            uOverLapData.append(l)

    uOverLapData = scipy.array(uOverLapData)        
    zOverLapData = scipy.array(zOverLapData)        

    uSpec = scipy.array(uSpec)        
    zSpec = scipy.array(zSpec)        

    zRescale = scipy.median(sdssSpline(zOverLapData[:,0])/ zOverLapData[:,1])
    uRescale = scipy.median(sdssSpline(uOverLapData[:,0])/ uOverLapData[:,1])

    uSpec = scipy.array(zip(uSpec[:,0],uRescale*uSpec[:,1]))
    zSpec = scipy.array(zip(zSpec[:,0],zRescale*zSpec[:,1]))

    fit_list = ['USDSS_GSDSS','GSDSS_ZSDSS','ISDSS_ZSDSS']        

    def errfunc(p,plot_it=False, getSpec=False):

        zWarp = scipy.interpolate.interp1d(zOverLap + [11000], [1.,1.,abs(p[0])], 
                    bounds_error = False, 
                    fill_value = 1.)

        ''' if locus is point is sufficiently red that we can fit for u-band '''
        if len(p) == 2: 
            uWarp = scipy.interpolate.interp1d([2500]+uOverLap, [abs(p[1]),1.,1.],  
                        bounds_error = False, 
                        fill_value = 1.)
        else: 
            ''' just array filled with ones, when not fitting for u-band slope  '''
            uWarp = scipy.interpolate.interp1d([2500]+uOverLap, [1.,1.,1.],  
                        bounds_error = False, 
                        fill_value = 1.)


        specStitch_0 = (uSpec[:,0].tolist() + specSDSS[:,0].tolist() + zSpec[:,0].tolist())
        specStitch_1 = (uSpec[:,1].tolist() + specSDSS[:,1].tolist() + zSpec[:,1].tolist())

        specStitch = scipy.array(zip(specStitch_0,specStitch_1*uWarp(specStitch_0)*zWarp(specStitch_0)))
        mags = synth([1.,0,0,0],[[specStitch]],filters) 

        if plot_it:

            pylab.clf()                              
            pylab.plot(pickleSpec[:,0],uWarp(pickleSpec[:,0])*zWarp(pickleSpec[:,0]),color='red')
            pylab.axvline(sdssLimits[0],c='black')
            pylab.axvline(sdssLimits[1],c='black')
            pylab.xlim([3000,10000])
            [x0,x1,y0,y1] = pylab.axis()
            pylab.ylim([y0-0.1,y1+0.1])
            pylab.text((sdssLimits[0] + 3000.)/2.,0.2*(y1-y0) + y0,'Pickles',color='black',horizontalalignment='center',fontsize=16)
            pylab.text((sdssLimits[0] + sdssLimits[1])/2.,0.2*(y1-y0) + y0,'SDSS',color='blue',horizontalalignment='center',fontsize=16)
            pylab.text((sdssLimits[1] + 10000.)/2.,0.2*(y1-y0) + y0,'Pickles',color='black',horizontalalignment='center',fontsize=16)
            pylab.xlabel('Angstroms',fontsize='large')            
            pylab.ylabel('Normalization',fontsize='large')
            pylab.savefig(plot_directory + str(locus_index) + '_norm.png')
            #pylab.show()

            pylab.clf()                              
            pylab.plot(specStitch[:,0],specStitch[:,1])
            pylab.axvline(sdssLimits[0],c='black')
            pylab.axvline(sdssLimits[1],c='black')
            pylab.xlim([3000,10000])
            [x0,x1,y0,y1] = pylab.axis()
            pylab.text((sdssLimits[0] + 3000.)/2.,0.8*(y1-y0) + y0,'Pickles',color='black',horizontalalignment='center',fontsize=16)
            pylab.text((sdssLimits[0] + sdssLimits[1])/2.,0.8*(y1-y0) + y0,'SDSS',color='blue',horizontalalignment='center',fontsize=16)
            pylab.text((sdssLimits[1] + 10000.)/2.,0.8*(y1-y0) + y0,'Pickles',color='black',horizontalalignment='center',fontsize=16)
            pylab.xlabel('Angstroms',fontsize='large')            
            pylab.ylabel('Flux',fontsize='large')
            pylab.savefig(plot_directory + str(locus_index) + '_spec.png')
            #pylab.show()

        ugdiff = (mags['USDSS'] - mags['GSDSS'] - locus_list['USDSS_GSDSS'][locus_index])
        gzdiff = (mags['GSDSS'] - mags['ZSDSS'] - locus_list['GSDSS_ZSDSS'][locus_index])
        izdiff = (mags['ISDSS'] - mags['ZSDSS'] - locus_list['ISDSS_ZSDSS'][locus_index])
        ridiff = (mags['RSDSS'] - mags['ISDSS'] - locus_list['RSDSS_ISDSS'][locus_index])
        if len(p) == 2:
            stat = ( ugdiff**2. + gzdiff**2.  + izdiff**2. + ridiff**2.)
        else:
            stat = ( gzdiff**2.  + izdiff**2. + ridiff**2.)
       
        if getSpec: return specStitch
        else:
            return stat

    ''' DISABLED if locus point is sufficiently blue that we can fit for u-band '''
    if True: #locus_list['RSDSS_ISDSS'][locus_index] < 0.45:
        pinit = [1.,1.]                                   
        out = scipy.optimize.fmin(errfunc,pinit,args=()) 
    else:
        pinit = [1.]                                   
        out = scipy.optimize.fmin(errfunc,pinit,args=()) 


    stitchSpec = errfunc(out,plot_it=True,getSpec=True)

    mags = synth([1.,0,0,0],[[stitchSpec]],filters) 


    return stitchSpec


def find_most_similar_pickles_spectrum(input, spectra_complete=None):

    print 'FINDING PICKLES SPECTRUM WITH MOST SIMILAR SHAPE TO SDSS SPECTRUM'

    sdssSpectrum = copy(input)
    sdssSpectrum[:,1] = sdssSpectrum[:,1] / (scipy.ones(len(sdssSpectrum[:,1]))*scipy.median(sdssSpectrum[:,1]))


    diffs = []

    for i in range(len(spectra_complete)):
        sp = spectra_complete[i] 
        spectrum = sp[0]
        picklesSpline = scipy.interpolate.interp1d(spectrum[:,0], spectrum[:,1], 
                                   bounds_error = False, 
                                   fill_value = 0.)


        specInterp = picklesSpline(sdssSpectrum[:,0])
        
        specInterp = specInterp / (scipy.ones(len(sdssSpectrum[:,1]))*scipy.median(specInterp))

        diff = specInterp - sdssSpectrum[:,1]
        diff = diff - scipy.ones(len(diff))*scipy.median(diff)
        stat = abs(diff).sum()
        diffs.append([stat,i])

    diffs.sort()


    sp = spectra_complete[diffs[0][1]] 
    spectrum = sp[0]
    picklesSpline = scipy.interpolate.interp1d(spectrum[:,0], spectrum[:,1], 
                               bounds_error = False, 
                               fill_value = 0.)
                                                                                          
    specInterp = picklesSpline(sdssSpectrum[:,0])
    
    specInterp = specInterp / (scipy.ones(len(sdssSpectrum[:,1]))*scipy.median(specInterp))
                                                                                          
    diff = specInterp - sdssSpectrum[:,1]
    diff = diff - scipy.ones(len(diff))*scipy.median(diff)

    specAll = scipy.array(zip(spectrum[:,0], spectrum[:,1] / (scipy.ones(len(spectrum[:,1]))*scipy.median(specInterp))))
    if False: 
        pylab.clf()                                     
        pylab.plot(specAll[:,0],specAll[:,1])
        pylab.plot(sdssSpectrum[:,0],sdssSpectrum[:,1])
        pylab.plot(sdssSpectrum[:,0],diff)
                                                        
        pylab.xlim(3000,11000)
        pylab.show()

    ''' need to fit spectral ends to reproduce locus color '''        
    return sdssSpectrum, specAll

    

#def load_spectra():
#    f = open('picklespectra','r')
#    m = pickle.Unpickler(f)
#    spectra = m.load()
#    
#    return spectra

def cas_locus():
    import pickle
    f = open(os.environ['BIGMACSMAKE'] + 'lociCAS','r')
    m = pickle.Unpickler(f)
    locus_list_mag = m.load()

    #print locus_list_mag


    return locus_list_mag

    #print len(locus_list_mag)

    #locus = []
    #for i in range(len(locus_list_mag['GSDSS_RSDSS'])):
    #    mag = {}
    #    mag['psfmag_g'] = 0 
    #    mag['psfmag_u'] = locus_list_mag['USDSS_GSDSS'][i] + mag['psfmag_g']
    #    mag['psfmag_r'] = mag['psfmag_g'] - locus_list_mag['GSDSS_RSDSS'][i]
    #    mag['psfmag_i'] = mag['psfmag_r'] - locus_list_mag['RSDSS_ISDSS'][i]
    #    mag['psfmag_z'] = mag['psfmag_i'] - locus_list_mag['ISDSS_ZSDSS'][i]
    #    mag['psfmag_J'] = mag['psfmag_z'] - locus_list_mag['ZSDSS_JTMASS'][i]
#
#        mag['USDSS'] = mag['psfmag_u']
#        mag['GSDSS'] = mag['psfmag_g']
#        mag['RSDSS'] = mag['psfmag_r']
#        mag['ISDSS'] = mag['psfmag_i']
#        mag['ZSDSS'] = mag['psfmag_z']
#        if i % 6 == 0:
#            locus.append(mag)
#
#    return locus 




def covey_locus():
    f = open(os.environ['BIGMACSMAKE'] + '/locus.txt','r').readlines()
    id = -1
    rows = {}
    colors = {}
    for i in range(len(f)):
        l = f[i]
        if l[0] != ' ':
            rows[i] = l[:-1]
        else: 
            id += 1 
            colors[rows[id]] = [float(x) for x in re.split('\s+',l[:-1])[1:]]

    return colors


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

''' load pickles spectra from text file '''
def load_pickles_spectra():
    pickles_location = os.environ['BIGMACSMAKE'] + '/PICKLES/'
    print 'LOADING PICKLES SPECTRA FROM ' + pickles_location 
    spectrafiles = glob(pickles_location + '/*v.dat')[:] # ONLY DWARF STARS
    #print spectrafiles
    if not spectrafiles: 
        print 'need to add pickles library and put in the ./pickles/ subdirectory'
        print 'available at: http://www.ifa.hawaii.edu/users/pickles/AJP/hilib.html'
        raise Exception
    spectra = [[scipy.loadtxt(s),s] for s in spectrafiles]

    return spectra


''' calculate synthetic SDSS magnitudes from Pickles spectra '''
def pickles_synthetic_magnitudes():

    spectra = load_pickles_spectra()

    filters = get_filters()
    nspectra = len(spectra)
        
    ''' interpolate only on the filter '''

    spec_mags = []

    for spec,name in spectra:
        star = {'name':name} 
        for filterSpline, step, filt_name in filters:
            specStep = spec[1,0] - spec[0,0] # wavelength increment                   
            resampFilter = filterSpline(spec[:,0]) # define an interpolating function
            val =    sum(specStep * resampFilter * spec[:,1])
            logEff = scipy.log10(val)
            logNorm = scipy.log10(sum(resampFilter*c*specStep/spec[:,0]**2))
            mag = 2.5*(logNorm - logEff) # to calculated an AB magnitude
            star[filt_name] = mag
        spec_mags.append(star)

    #f = open('picklelocus_MACS','w')
    #m = pickle.Pickler(f)
    #pickle.dump(spec_mags,m)
    #f.close()
        
    return spec_mags

def synth(p,spectra,filters,show=False):

    mags ={} 

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


def plot():
   
    spectra_complete = load_pickles_spectra() 

    filters = get_filters(False)

    stars = pickles_synthetic_magnitudes() 

    locus_list = cas_locus()
    comp_list =  filter(lambda x: string.find(x.replace('SDSS_',''),'SDSS')!=-1 and string.find(x,'SDSS_')!=-1, locus_list.keys())

    #print locus_list.keys()
    close_locus = []

    if 0:
        fit_mags = []                                                                                             
                                                                                                                  
        for i in 5 * scipy.array(range(len(locus_list[comp_list[0]])/5)): 
            star_stats = []
            for s in range(len(stars)):
                stat = 0                 
                for combo in comp_list:
                    res = re.split('\_',combo)
                    f1 = res[0]
                    f2 = res[1]
                    stat += ((stars[s][f1]-stars[s][f2]) - locus_list[combo][i])**2.
                star_stats.append([stat,copy(s)])
                                                                                                                  
            star_stats.sort()
                                                                                                                  
            spectra_sub = [spectra_complete[x[1]] for x in star_stats[:4]]
                                                                                                                  
            if True:
                mags = synth([1,0,0,0],spectra_sub,filters)
                for combo in comp_list:                                                                        
                    res = re.split('\_',combo)
                    f1 = res[0]
                    f2 = res[1]
                    #print mags[f1] - mags[f2], locus_list[combo][star_stats[0][1]], f1, f2
                                                                                                                  
                                                                                                                  
            close_locus.append(star_stats[0][1])
            close_locus.append(star_stats[1][1])

            #print spectra_sub
            pinit = [1,0,0,0] 
            locus_index = i
            out = scipy.optimize.fmin(errfunc,pinit,xtol=0.005,ftol=0.001,args=(spectra_sub,locus_list,locus_index,comp_list,filters)) 
            mags = synth(out,spectra_sub,filters,show=False) 
            fit_mags.append([mags,out,spectra_sub,copy(i)])
                                     
        f = open('maglocus','w')
        m = pickle.Pickler(f)
        pickle.dump(fit_mags,m)

    f = open('maglocus_SYNTH','r')
    m = pickle.Unpickler(f)
    fit_mags = m.load()


    pylab.clf()

    c1 = []
    c2 = []
    for i in range(len(stars)):
        c1.append(stars[i]['GSDSS']-stars[i]['RSDSS'])
        c2.append(stars[i]['RSDSS']-stars[i]['ISDSS'])

    pylab.scatter(c1,c2,color='green')
    
    pylab.scatter(locus_list['GSDSS_RSDSS'], locus_list['RSDSS_ISDSS'])

    c1 = []
    c2 = []
    for i in close_locus:
        #print len(stars)
        c1.append(stars[i]['GSDSS']-stars[i]['RSDSS'])
        c2.append(stars[i]['RSDSS']-stars[i]['ISDSS'])
    #print c1, c2

    c1 = []
    c2 = []
    for i in range(len(fit_mags)):
        c1.append(fit_mags[i]['GSDSS']-fit_mags[i]['RSDSS'])
        c2.append(fit_mags[i]['RSDSS']-fit_mags[i]['ISDSS'])

        #c1.append(fit_mags[i]['RSDSS']-fit_mags[i]['CAPAKIS'])
        #c2.append(fit_mags[i]['RSDSS']-fit_mags[i]['WSISUBARU'])


    pylab.scatter(c1,c2,color='red')

    pylab.show()

def synthesize_expected_locus_for_observations(filters):
    c_locus = cas_locus()    

    ''' add SDSS filters '''        
    SDSS_filters = [{'name':'USDSS','filter':'SDSS-u.res'},{'name':'GSDSS','filter':'SDSS-g.res'},{'name':'RSDSS','filter':'SDSS-r.res'},{'name':'ISDSS','filter':'SDSS-i.res'},{'name':'ZSDSS','filter':'SDSS-z.res'}]
    filter_info = get_filters([[a['name'],a['filter']] for a in SDSS_filters])
    for i in range(len(filter_info)):
        SDSS_filters[i].update(filter_info[i])

    loci = [a[:-1] for a in open(os.environ['BIGMACSMAKE'] + '/LOCUS_SPECTRA/spliced_spectra','r').readlines()]

    locus = []

    print 'STARTING SAMPLING LOCUS'

    for i in range(len(loci[:])):
        locus_point = loci[i]
        print 'CONVOLVING RESPONSE FUNCTIONS WITH SPECTRUM ' + str(locus_point)
        stitchSpec = scipy.genfromtxt(os.environ['BIGMACSMAKE'] + '/LOCUS_SPECTRA/' + locus_point)
        mags = synth([1.,0,0,0],[[stitchSpec]],filter(lambda x: string.find(x['filter'],'TMASS') == -1, filters + SDSS_filters)) 

        for filt in filters:
            if filt['filter'] == 'JTMASS.res':
                mags[filt['mag']] = (mags['ZSDSS'] - c_locus['ZSDSS_JTMASS'][2*i]) #+ (mags['ISDSS'] - c_locus['ZSDSS_JTMASS'][2*i]))/2.

        locus.append(mags)

    print 'FINISHED SAMPLING LOCUS'
    return locus

def parse_columns(columns_description, fitSDSS=False, noHoldExceptSDSS=False):
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
        else:
            dict['HOLD_VARY'] = l[3] 
            if dict['HOLD_VARY'] == 'HOLD':
                dict['ZP'] = float(l[4])
        input_info.append(dict)        

    return input_info 
    

def get_filters(flist = [['USDSS','SDSS-u.res'],['GSDSS','SDSS-g.res'],['RSDSS','SDSS-r.res'],['ISDSS','SDSS-i.res'],['ZSDSS','SDSS-z.res']]):

    filt_dir = os.environ['BIGMACSMAKE'] + '/FILTERS/'

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
    sed = scipy.loadtxt(os.environ['BIGMACSMAKE'] + '/LIBRARIES/munari.sed')
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

