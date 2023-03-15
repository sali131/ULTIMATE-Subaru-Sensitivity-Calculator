#!/usr/bin/env python
import math
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import numpy as np
import pandas as pd
import sys
sys.path.append('/Users/sali/python')
from figure import figure

wv = 16
am = 15

fbg_s = pd.read_csv('mk_skybg_zm_'+str(wv)+'_'+str(am)+'_ph_wo_skybb.txt', sep=' ', header=None)
fbg_s = fbg_s[(fbg_s.iloc[:, 0] > 900) & (fbg_s.iloc[:, 0] < 2501)]
fbg_s = fbg_s[(fbg_s.iloc[:, 1] > 0.1) & (fbg_s.iloc[:, 1] < 1.3)]

wav_s = fbg_s.iloc[:, 0].to_list()
int_s = fbg_s.iloc[:, 1].to_list()

wav = np.arange (900, 2501, 0.02)

skybg_spline = UnivariateSpline(wav_s, int_s, s=10000)
skybg = skybg_spline(wav)

fbg = pd.DataFrame({'wav': np.round(wav, 2),
                   'int': skybg})
fbg.to_csv('mk_skybg_zm_'+str(wv)+'_'+str(am)+'_ph_const_wo_skybb.txt', sep='\t', index=False)

plt.ylim([0,1.5])
plt.xlim([900, 2500])
plt.plot(wav_s, int_s, color='red', linewidth=0.1, label='Gemini Model')
plt.plot(wav, skybg, color='blue', linewidth=1, label='Spline Fit')

figure('Sky_bg_model', 'Wavelength [nm]', 'Intensity [Photons/s/arcsec^2/nm/m^2')
