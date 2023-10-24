#!/usr/bin/env python
import math
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import numpy as np
import pandas as pd

print('WFI Sensitivity Calculator')
print('Input desired values or press return to use default values')
print('**************************')

instrument = input('Filter options are \'moircs\', \'swims\' or \'custom\' (default is moircs): ') or 'moircs'

if instrument == 'moircs':
    filter = input('MOIRCS filter options are \'bb\' or \'nb\' (default is bb): ') or 'bb'
    filename = instrument+'_'+filter+'.asc'
elif instrument == 'swims':
    filter = input('SWIMS filter options are \'bb\', \'mb\' or \'nb\' (default is bb): ') or 'bb'
    filename = instrument+'_'+filter+'.asc'
else:
    filename = 'custom.asc'
    filter = 'custom'

ao = input('AO options are \'glao\' or \'noao\' (default is glao): ') or 'glao'

hrs = float(input('Exposure time in hours (default is 1 hours): ') or '1.0')
snr = float(input('Desired SNR (default is 5): ') or '5')
#NOAO = 1.2/1.9/3.0
#GLAO = 0.65/1.15/2.25
Tr_inst_filter = float(input('Desired Instrument+Filter Throughput (default is 0.45): ') or '0.45')

wv_og = float(input('Water Vapour, choose from 1.0, 1.6, 3.0, 5.0 mm (default is 1.0): ') or '1.0')
am_og = float(input('Airmass, choose from 1.0, 1.5, 2.0 (default is 1.0): ') or '1.0')
wv = int(str(wv_og).replace('.',''))
am = int(str(am_og).replace('.',''))

width = 3.487
height = width / 1.618

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth']=1.5
plt.rcParams['xtick.major.size']=10
plt.rcParams['xtick.major.width']=1.5
plt.rcParams['ytick.major.size']=10
plt.rcParams['ytick.major.width']=1.5

# Telescope mirror reflectivity 
R_Tel_M1 = 0.95 
R_Tel_M2 = 0.98
R_Tel_CC = 0.98

# Telescope mirror emissivity 
e_Tel_M1 = 1.0 - R_Tel_M1
e_Tel_M2 = 1.0 - R_Tel_M2
e_Tel_CC = 1.0 - R_Tel_CC

# Primary Over illumination factor
rho = 1.0

# Telescope diameter and Central Obscuration
Dp = 808.19 #Primary mirror diameter in [cm]
vc = 0.277  #Central obscuration due to the center cone 

# Telescope Focal length
fl = 10000.0 # Effective focal length in [cm]
fr = fl / Dp # Focal ratio

# Telescope environment temperature 
Temp_Env = 273.0 # Environment temperature [K]

# physical constants
h = 6.626e-27 # erg sec (planck constant)
c = 2.998e+8 # m/sec (speed of light)

# lam in micron
# T in Kelvin
# output in photon/sec/cm^2/A/sr
def bb_radiation(lam, T):
    #h = 6.626e-27 # erg sec (planck constant)
    #kb = 1.38e-16 # erg/K (boltzman constant) 
    #c = 2.998e+10 # cm/sec (speed of light)

    # Black Body radiation in [photon/sec/cm^2/A/sr]
    bb = 5.995e+18*(1.0/lam**4)*(1.0/(math.exp(14387.8/(lam*T))-1.0))
    
    return bb

# return photon count at the detector [e-/sec] for m0 magnitude star in AB
# m0 : AB magnitude
# Tr_sky: Sky throughput
# Tr_tel: Telescope mirror throughput
# Tr_inst: Instrument throughput (optics, detector QE)
# Tr_filter: filter throughput
# l0 : Filter central wavelength [micron]
# delta: Filter bandwidth [micron]
def f0(m0, Tr_sky, Tr_tel, Tr_inst, Tr_filter, l0, delta):

    # fulx density that corresponds to m0 in AB [erg/sec/cm^2/A]
    f0_cal = 10.0**(-0.4*(m0 + 5.0*math.log10(l0*1.0e+4) + 2.408))

    # detector count [e-/sec]
    single_photon_energy = h * c * 1.0e+6 / l0 # [erg]
    ph = f0_cal * delta*1.0e+4 * (math.pi*0.25*Dp**2*(1.0 - vc**2)) * Tr_sky * Tr_tel * Tr_inst * Tr_filter / single_photon_energy

    return ph # [e-/sec]

# Sky background from Gemini
# http://www.gemini.edu/sciops/telescopes-and-sites/observing-condition-constraints/ir-background-spectra
# Return [photn/sec/A/arcsec^2]
# Tsky: Sky temperature
# f_sky: scale factor to the sky emission line
def sky_bg_calc(l0, delta, f_sky, Tsky):

#    fp = open("MaunaKea_Sky/backup/mk_skybg_zm_16_15_ph_wo_skybb.dat")
    fp = open('MaunaKea_Sky/ultimate/mk_skybg_zm_'+str(wv)+'_'+str(am)+'_ph_wo_skybb.txt')

    sky_lam_arr = []
    sky_bg_arr = []
    
    for line in fp:
        param = line.strip().split()
        lam_skybg = float(param[0])
        bg_skybg = float(param[1])
    
        if lam_skybg >= (l0-0.5*delta)*1.0e+3 and lam_skybg <= (l0+0.5*delta)*1.0e+3:

            sky_throughput = np.interp(lam_skybg/1000.0, skt_lam, skt_tr)
            sky_emissivity = 1.0 - sky_throughput
            bg_skybg = f_sky*bg_skybg + sky_emissivity * bb_radiation(lam_skybg/1000.0, Tsky)*(1.0e+5/4.25e+10)

            sky_lam_arr.append(lam_skybg/1.0e+3) # micron
            sky_bg_arr.append(bg_skybg*1.0e-5*(math.pi*0.25*Dp**2)) # ph/sec/A/arcsec^2
            
    fp.close()
    return np.array(sky_lam_arr), np.array(sky_bg_arr)

# Sky transmittance from Gemini
# http://www.gemini.edu/sciops/telescopes-and-sites/observing-condition-constraints/ir-transmission-spectra
# Maunakea sky transmittance
ft = open('MaunaKea_Sky/ultimate/mktrans_zm_'+str(wv)+'_'+str(am)+'.txt')
#ft = open('MaunaKea_Sky/mktrans_zm_const.txt')
skt_lam = []
skt_tr = []
for line in ft:
    param = line.strip().split()
    skt_lam.append(float(param[0]))
    skt_tr.append(float(param[1]))
ft.close()

skt_lam = np.array(skt_lam) # wavelength in micron
skt_tr = np.array(skt_tr)   # transmittance

# Return Transmittance
def sky_trans(l0, delta):
    
    tr_sum = 0
    n_sum = 0
    for i in range(len(skt_lam)):
        if skt_lam[i] > l0 - delta*0.5 and skt_lam[i] < l0 + delta*0.5:
            tr_sum += skt_tr[i]
            n_sum += 1
    
    return tr_sum / n_sum

# l0: Filter central wavelength in micron
# delta: Filter bandwidth in micron
# Tsky: Sky temperature in Kelvin
# fsky: Sky emission scale factor
def bgcalc(l0, delta, fsky, Tsky, Tr_filter, Tr_inst_filter):
    
    Tr_sky = sky_trans(l0, delta) # average sky transmittance
    Tr_tel = R_Tel_M1 * R_Tel_M2 # Primary and Secondary
#    Tr_inst = 0.863 * 0.957 * 0.85 # FourStar, Persson et al. 2013
#    Tr_filter = 0.96 # MOIRCS Ks-band filter transmittance
#    Tr_inst = eta/(Tr_sky*Tr_tel*Tr_filter)
    Tr_inst = Tr_inst_filter/Tr_filter
    
    # calculate sky background level [photon/sec/arcsec^2]
    sky_lam, sky_bg = sky_bg_calc(l0, delta, fsky, Tsky)
    
    # average wavelength step width in Angstrom
    dw_ave = 0
    dw_num = 0
    for i in range(len(sky_lam)-1):
        dw_ave += (sky_lam[i+1] - sky_lam[i])*1.0e+4
        dw_num += 1
    dw_ave = dw_ave / dw_num

    # Total detector count from sky background [e-/sec/arcsec^2]
    f_skybg = np.sum(sky_bg*dw_ave*Tr_tel*Tr_inst*Tr_filter)

    m0 = 15.0
    #print m0 - 2.5*math.log10(f_skybg/f0(m0, Tr_sky, Tr_tel, Tr_inst, Tr_filter, l0, delta))

    A_arcsec2 = (fl * (1.0/3600.0) * (math.pi/180.0))**2 # [cm^2]
    Omega_t = (math.pi/(4.0*fr**2))*((e_Tel_M2+e_Tel_M1*(1.0-e_Tel_M2))*(1.0-vc**2) + e_Tel_CC*vc**2 + rho**2 - 1.0)
    
    f_telbg = 0
    for i in range(len(sky_lam)):
        f_telbg += bb_radiation(sky_lam[i], Temp_Env) * Tr_inst * Tr_filter * dw_ave * A_arcsec2 * Omega_t # [e-/sec/arcsec^2]

    mbg = m0 - 2.5*math.log10((f_skybg+f_telbg)/f0(m0, Tr_sky, Tr_tel, Tr_inst, Tr_filter, l0, delta))
    
    return mbg, f_skybg, f_telbg, Tr_tel, Tr_sky, Tr_inst

def filter_info(file_name):
    fbg = pd.read_csv(file_name, header = None, delimiter=" |\t", engine='python')
    filter_name = fbg.iloc[:, 0].to_list()
    filter_wc = fbg.iloc[:, 1].to_list()
    filter_dw = fbg.iloc[:, 2].to_list()
    Tr_filter = fbg.iloc[:, 3].to_list()
    Tr_inst = fbg.iloc[:, 4].to_list()
    return filter_name, filter_wc, filter_dw, Tr_filter, Tr_inst

def filter_info_custom(file_name):
    fbg = pd.read_csv(file_name, delimiter=" |\t", engine='python')
    filter_name = fbg['filter_name'].to_list()
    filter_wc = fbg['filter_central'].to_list()
    filter_dw = fbg['filter_width'].to_list()
    Tr_filter = fbg['filter_tr'].to_list()
#    filter_wc_pre = [1.03,1.25,1.29,1.64,2.14,1.5,1.61,1.73,2.02,2.17,2.32]
#    Tr_inst_pre = [0.4,0.4,0.4,0.47,0.44,0.48,0.55,0.47,0.48,0.49,0.45]
#    Tr_inst_spline = UnivariateSpline(filter_wc_pre, Tr_inst_pre, s=1)
#    Tr_inst = Tr_inst_spline(filter_wc)
#    Tr_inst = fbg.iloc[:, 4].to_list()
    Tr_inst = 1
    return filter_name, filter_wc, filter_dw, Tr_filter, Tr_inst

def wfi_cs_limiting_mag(filter_name, filter_wc, filter_dw, Tr_filter, instrument, hrs, snr):

    wfi_cs_limiting_mag = []

    for i in range(len(filter_name)):
        
        df = pd.read_csv('D50.csv')
        wavenc = df['Wav'].to_numpy()
        
        glao_good = df['D50_glao_good'].to_numpy()
        glao_median = df['D50_glao_median'].to_numpy()
        glao_bad = df['D50_glao_bad'].to_numpy()
        noao_good = df['D50_noao_good'].to_numpy()
        noao_median = df['D50_noao_median'].to_numpy()
        noao_bad = df['D50_noao_bad'].to_numpy()
        
        sglao_good = UnivariateSpline(wavenc, glao_good, k=2)
        sglao_median = UnivariateSpline(wavenc, glao_median, k=2)
        sglao_bad = UnivariateSpline(wavenc, glao_bad, k=2)
        snoao_good = UnivariateSpline(wavenc, noao_good, k=2)
        snoao_median = UnivariateSpline(wavenc, noao_median, k=2)
        snoao_bad = UnivariateSpline(wavenc, noao_bad, k=2)
        
        sl_glao_good = sglao_good(filter_wc[i])
        sl_glao_median = sglao_median(filter_wc[i])
        sl_glao_bad = sglao_bad(filter_wc[i])
        sl_noao_good = snoao_good(filter_wc[i])
        sl_noao_median = snoao_median(filter_wc[i])
        sl_noao_bad = snoao_bad(filter_wc[i])
        
        fsky_min = 1.0
        Tsky_min = 273.0

        if filter_wc[i] > 1.4 and filter_wc[i] < 1.8:
            fsky_min = 0.58

        if filter_wc[i] < 1.3:
            fsky_min = 1.0

        bgres = bgcalc(filter_wc[i], filter_dw[i], fsky_min, Tsky_min, Tr_filter[i], Tr_inst_filter)

        bg_total = bgres[1] + bgres[2] # e/sec/arcsec^2
        mbg = bgres[0]
        f_skybg = bgres[1]
        f_telbg = bgres[2]
        Tr_tel = bgres[3]
        Tr_sky = bgres[4]
        Tr_inst = bgres[5]
        
        texp = 3600.0 * hrs # sec
        snr = snr
        A = math.pi * Dp**2 * (1 - vc**2) / 4.0

        # eta & D50
        eta = Tr_filter[i] * Tr_inst * Tr_sky * Tr_tel
        
        # Single photon energy
        p =  h * c * 1.0e+6 / filter_wc[i] # [erg]
        
        flam1 = 2.0 * snr * p * math.sqrt(bg_total * (math.pi*sl_glao_good**2/4.0) * texp) / (A * texp * eta * filter_dw[i]*1.0e+4)
        flam2 = 2.0 * snr * p * math.sqrt(bg_total * (math.pi*sl_glao_median**2/4.0) * texp) / (A * texp * eta * filter_dw[i]*1.0e+4)
        flam3 = 2.0 * snr * p * math.sqrt(bg_total * (math.pi*sl_glao_bad**2/4.0) * texp) / (A * texp * eta * filter_dw[i]*1.0e+4)
        flam4 = 2.0 * snr * p * math.sqrt(bg_total * (math.pi*sl_noao_good**2/4.0) * texp) / (A * texp * eta * filter_dw[i]*1.0e+4)
        flam5 = 2.0 * snr * p * math.sqrt(bg_total * (math.pi*sl_noao_median**2/4.0) * texp) / (A * texp * eta * filter_dw[i]*1.0e+4)
        flam6 = 2.0 * snr * p * math.sqrt(bg_total * (math.pi*sl_noao_bad**2/4.0) * texp) / (A * texp * eta * filter_dw[i]*1.0e+4)
        
        m1 = -2.5*math.log10(flam1) - 5.0*math.log10(filter_wc[i]*1.0e+4) - 2.408
        m2 = -2.5*math.log10(flam2) - 5.0*math.log10(filter_wc[i]*1.0e+4) - 2.408
        m3 = -2.5*math.log10(flam3) - 5.0*math.log10(filter_wc[i]*1.0e+4) - 2.408
        m4 = -2.5*math.log10(flam4) - 5.0*math.log10(filter_wc[i]*1.0e+4) - 2.408
        m5 = -2.5*math.log10(flam5) - 5.0*math.log10(filter_wc[i]*1.0e+4) - 2.408
        m6 = -2.5*math.log10(flam6) - 5.0*math.log10(filter_wc[i]*1.0e+4) - 2.408
        
        wfi_cs_limiting_mag_filter = [filter_name[i], filter_wc[i], filter_dw[i], Tr_filter[i], Tr_inst, Tr_tel, Tr_sky, eta, m1, m2, m3, m4, m5, m6, f_skybg, f_telbg, bg_total]
        wfi_cs_limiting_mag.append(wfi_cs_limiting_mag_filter)
    
    lm = pd.DataFrame(wfi_cs_limiting_mag, columns = ['name', 'wav_central', 'filter_width', 'filter_tr', 'inst_tr', 'tel_tr', 'sky_tr', 'eta', 'mag_glao25_seeing_good', 'mag_glao50_seeing_moderate', 'mag_glao75_seeing_bad', 'mag_noao25_seeing_good', 'mag_noao50_seeing_moderate', 'mag_noao75_seeing_bad', 'counts_skybg', 'counts_telbg', 'counts_bg_total'])
    
    lm = lm.round(3)
    
    lm.to_csv('wfi_'+instrument+'_'+filter+'_am='+str(am)+'_wv='+str(wv)+'_tr='+str(Tr_inst_filter)+'.csv', index=False)
    return lm

def plot_sens(lm, filters, ao, instrument, snr, hrs):
    if filters == 'nb':
#        lm = lm[(lm.Wave_w < 0.05)]
        lgd_dim = [0.03, 0.125]
    elif filters == 'mb':
#        lm = lm[(lm.Wave_w > 0.09) & (lm.Wave_w < 0.15)]
        lgd_dim = [0.02, 0.125]
    elif filters == 'bb':
#        lm = lm[(lm.Wave_w > 0.15)]
        lgd_dim = [0.01, 0.1]
    else:
        lgd_dim = [0, 0.125]
    
    name = lm['name'].to_numpy()
    wc_arr = lm['wav_central'].to_numpy()
    dw_arr = lm['filter_width'].to_numpy()*0.5
    ft_arr = lm['filter_tr'].to_numpy()
    g25_arr = lm['mag_glao25_seeing_good'].to_numpy()
    g50_arr = lm['mag_glao50_seeing_moderate'].to_numpy()
    g75_arr = lm['mag_glao75_seeing_bad'].to_numpy()
    n25_arr = lm['mag_noao25_seeing_good'].to_numpy()
    n50_arr = lm['mag_noao50_seeing_moderate'].to_numpy()
    n75_arr = lm['mag_noao75_seeing_bad'].to_numpy()
    
    if ao == 'glao':
        y = [g25_arr, g50_arr, g75_arr]
        title = 'WFI at Cassegrain (With GLAO)'+' Airmass='+str(am_og)+' Water Vapor='+str(wv_og)+'mm'
    else:
        y = [n25_arr, n50_arr, n75_arr]
        title = 'WFI at Cassegrain (No GLAO)'+' Airmass='+str(am_og)+' Water Vapor='+str(wv_og)+'mm'
    
    plt.errorbar(wc_arr, y[0], xerr=dw_arr, fmt='o', ecolor='blue', capsize=5, markersize=10,  color = "blue", label='Good seeing')
    plt.errorbar(wc_arr, y[1], xerr=dw_arr, fmt='o', ecolor='red', capsize=5, markersize=10,  color = "red", label='Moderate seeing')
    plt.errorbar(wc_arr, y[2], xerr=dw_arr, fmt='o', ecolor='green', capsize=5, markersize=10,  color = "green", label='Bad seeing')
    
    for i, txt in enumerate(name):
        plt.annotate(txt, (wc_arr[i], y[2][i]), xytext=(wc_arr[i]-lgd_dim[0], y[2][i]-lgd_dim[1]))

    plt.xlabel('Wavelength [micron]')
    plt.ylabel(str(snr)+'-sigma limiting AB mag (exp time = '+str(hrs)+' hrs)')
#    plt.xlim(1.1,2.4)
#    plt.ylim(25.0,26.45)
    plt.grid()
    plt.legend(loc='lower left')
    plt.title(title)
    
    plt.savefig('wfi_'+instrument+'_'+filters+'_am='+str(am_og)+'_wv='+str(wv_og)+'_tr='+str(Tr_inst_filter)+'.png', format='png', dpi=400, width=width, height=height, bbox_inches='tight')
    return

try:
    filter_info =  filter_info(filename)
except IndexError:
    filter_info = filter_info_custom(filename)

lm = wfi_cs_limiting_mag(filter_info[0], filter_info[1], filter_info[2], filter_info[3], instrument=instrument, hrs=hrs, snr=snr)
plot_sens(lm, filters=filter, ao=ao, instrument=instrument, snr=snr, hrs=hrs)
