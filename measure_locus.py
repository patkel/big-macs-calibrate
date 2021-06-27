
def convert_to_pogson(p):

    import pyfits, scipy

    cols = []
    for col in p.columns:
        if col.name[0:7] == 'psfFlux':  
            array = -2.5*scipy.log10(p.data.field(col.name)) + 22.5
            cols.append(pyfits.Column(name=col.name.replace('psfFlux','psfPog'), format=col.format, array=array)) 
        cols.append(col)

    hdu = pyfits.PrimaryHDU()                                                  
    hdulist = pyfits.HDUList([hdu])
    print cols
    tbhu = pyfits.new_table(cols)
    #hdulist.append(tbhu)
    #hdulist[1].header.update('EXTNAME','STDTAB')
    #outcat = '/tmp/test' #path + 'PHOTOMETRY/' + type + '.cat'                
    #os.system('rm ' + f + '.tab')
    #hdulist.writeto(f + '.tab')
    return tbhu
        



import pyfits, pylab, scipy

p = pyfits.open('SLR_2MASS_pkelly50.fit')[1].data#[0:20000]

#p = convert_to_pogson(p).data

#p = p[p.field('extinction_u') < 0.3]

p = p[:200000]

#p2 = pyfits.open('locus_pkelly50.fit')[1]

pylab.close()
pylab.clf()

import make_spectroscopic_locus
cl = make_spectroscopic_locus.covey_locus()

#for c1,c2,c1_lim,c2_lim in [[['g','i'],['g','r'],[-0.5,3.0],[-0.5,2.0]]]:

#for c1,c2,c1_lim,c2_lim in [[['g','r'],['u','g'],[-0.5,2.0],[-1,5]]]:

loci = {}

c1 = ['g','i']
c1_lim = [0.3,3.5]


for c2 in [['u','g'],['r','i'],['g','r'],['i','z'],['g','i'],['z','J'],['g','z']]:
    #pylab.scatter(p.field('psfMag_' + c1[0])-p.field('psfMag_' + c1[1]),p.field('psfMag_' + c2[0])- p.field('psfMag_' + c2[1]),s=0.001)

    points = []
    incre = 0.02
    for x in scipy.arange(c1_lim[0],c1_lim[1],incre):        
        mask = (x < p.field('psfMag_'+c1[0]) - p.field('psfMag_' + c1[1])) * (x + incre > p.field('psfMag_'+c1[0]) - p.field('psfMag_'+c1[1])) #* (19 > p.field('psfMag_r')) # - p.field('psfMag_r')) * (-99 < p.field('psfMag_g') - p.field('psfMag_r'))

        ''' if including J band, require that 2MASS magnitude is accurate '''
        if c2[1] == 'J' or c2[0] == 'J':
            mask *= p.field('Jgood') == 1

        filt = p[mask]
        colors = filt.field('psfMag_' + c2[0]) - filt.field('psfMag_' + c2[1])
        #print colors
        points.append([x+incre/2., scipy.median(colors)] )                    
        
    points = scipy.array(points)
    print points


    if c2[1] != 'J':
        loci[c2[0].upper() + 'SDSS_' + c2[1].upper() + 'SDSS'] = points[:,1]
    else:
        loci[c2[0].upper() + 'SDSS_' + c2[1].upper() + 'TMASS'] = points[:,1]

    if True: #record:
        import pickle           
        f = open('lociCAS','w')
        m = pickle.Pickler(f)
        pickle.dump(loci,m) 
        f.close()

    else:
        import pickle 
        f = open('lociCAS','r')
        m = pickle.Unpickler(f)
        locus_list_mag = m.load()
        print locus_list_mag.keys()

        pylab.scatter(locus_list_mag['g_r'][:,1],locus_list_mag['r_i'][:,1],color='green')

    pylab.scatter(points[:,0], points[:,1],color='red')

    #c1_points = cl[c1[0].upper() + 'SDSS_' + c1[1].upper() + 'SDSS']
    #c2_points = cl[c2[0].upper() + 'SDSS_' + c2[1].upper() + 'SDSS']
    #pylab.scatter(c1_points[:], c2_points[:],color='green')

    #pylab.scatter(p2.field('psfMag_' + c1[0])-p2.field('psfMag_' + c1[1]),p2.field('psfMag_' + c2[0])- p2.field('psfMag_' + c2[1]),s=0.001, color='red')
    pylab.xlim(c1_lim)
    #pylab.ylim(c2_lim)
    pylab.xlabel(c1[0] + '-' + c1[1])
    pylab.ylabel(c2[0] + '-' + c2[1])
    
    #p = pyfits.open('locusext2_pkelly100.fit')[1]
    #pylab.scatter(p.field('Column2')-p.field('Column3'),p.field('Column1')- p.field('Column2'),s=0.05,color='red')
    
    #pylab.show()
