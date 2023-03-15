#!/usr/bin/env python
import math
import numpy as np

wv = 50
am = 20

# lam in micron
# T in Kelvin
# output in photon/sec/m^2/nm/arcsec^2
def bb_radiation(lam, T):
    #h = 6.626e-27 # erg sec (planck constant)
    #kb = 1.38e-16 # erg/K (boltzman constant) 
    #c = 2.998e+10 # cm/sec (speed of light)

    # Black Body radiation in [photon/sec/cm^2/A/sr]
    bb = 5.995e+18*(1.0/lam**4)*(1.0/(math.exp(14387.8/(lam*T))-1.0))*(1.0e+5/4.25e+10)
    
    return bb

ft = open('mktrans_zm_'+str(wv)+'_'+str(am)+'.txt')
skt_lam = []
skt_tr = []
for line in ft:
    param = line.strip().split()
    skt_lam.append(float(param[0]))
    skt_tr.append(float(param[1]))
ft.close()

skt_lam = np.array(skt_lam)
skt_tr = np.array(skt_tr)

fin = open('mk_skybg_zm_'+str(wv)+'_'+str(am)+'_ph.txt')
fout = open('mk_skybg_zm_'+str(wv)+'_'+str(am)+'_ph_wo_skybb.txt', "w")

for line in fin:
    param = line.strip().split()

    lam = float(param[0])
    bg = float(param[1])

    sky_e = 1.0 - np.interp(lam/1000.0, skt_lam, skt_tr)
    sky_cont = sky_e * bb_radiation(lam/1000.0, 273.0)
    bg = bg - sky_cont
    #print lam
    fout.write('%f %f %f\n' % (lam, bg, sky_cont))

fin.close()
fout.close()


    

