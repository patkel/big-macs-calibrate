detector = 't2kb'

for detector in ['t2ka','t2kb']:

    f = open(detector + '.response').readlines()
    
    wave_0,trans_0,x_0,y_0 = [float(x) for x in f[0][:-1].split(' ')]
    wave_1,trans_1,x_1,y_1 = [float(x) for x in f[1][:-1].split(' ')]
    
    waves = []
    resps = []
    
    for l in f[2:]:
        x,y = [float(p) for p in l[:-1].split(' ')]
        waves.append( 10. * ((x-x_0)/(x_1-x_0)*(wave_1-wave_0) + wave_0) )
        resps.append( 0.01 * (trans_1 + (y_1-y)/(y_0-y_1)*(trans_1-trans_0)) )
    
    file = open(detector + '.txt','w')
    for w,r in zip(waves,resps):
        file.write(str(w) + ' ' + str(r) + '\n')        
    
    file.close()    
        
    
    import pylab
    pylab.plot(waves,resps)
pylab.show()
        
