import utilities, pyfits, os

import pickle
f = open(os.environ['BIGMACSMAKE'] + 'lociCAS','r')
m = pickle.Unpickler(f)
locus_list_mag = m.load()




c_locus = locus_list_mag #utilities.cas_locus()

outcat = 'lociCAS.fits'
extension = 'STDTAB'

cols = []
for key in c_locus.keys():
    cols.append(pyfits.Column(name=key,format='D',array=c_locus[key]))
                                                                                        
hdu = pyfits.PrimaryHDU()
hdulist = pyfits.HDUList([hdu])
tbhu = pyfits.new_table(cols)
hdulist.append(tbhu)
hdulist[1].header.update('EXTNAME',extension)
os.system('rm ' + outcat)
hdulist.writeto(outcat)
print 'WRITTEN TO ', outcat


    
