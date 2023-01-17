import sys
import numpy as np
import pandas as pd
from scipy.stats import linregress
from scipy.stats import t
from scipy.stats import norm
from sklearn.utils import resample as bootstrap
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq, ifft, rfft, irfft

sys.path.append('../../processing')

import processing.stations as stns
import processing.utils as ut

filepath = '../../../megafires_data/QN/UmbralesdeOlasdeCalor(Diurna)_2023-01-07_20_52.csv'
df_umb = pd.read_csv(filepath, sep=';')

df = stns.get_QN_daily_tmax_1911_2022()

df['clim'] = df.loc[:,'1991':'2020'].mean(axis=1)

n = 30+31

temp_novdic = df.loc[:,'2022'].iloc[-n:]

p90_nov = df_umb.loc[:,'Nov'].iloc[:-1]
p90_dic = df_umb.loc[:,'Dic']


temp = temp_novdic.values
p90 = np.concatenate((p90_nov.values, p90_dic.values))
x = np.arange(n)


# Number of sample points
N = df.loc[:,'1991':'2020'].size
# sample spacing
T = 1/365 
# x = np.linspace(0.0, N*T, N, endpoint=False)
y = np.ravel(df.loc[:,'1991':'2020'].values, order='F')-np.mean(np.ravel(df.loc[:,'1991':'2020'].values, order='F'))
yf = rfft(y)
xf = fftfreq(N, T)#[:N//2]
newyf = yf*0
newyf[15] = yf[15]
newyf[30] = yf[30]
newy = irfft(newyf) + np.mean(np.ravel(df.loc[:,'1991':'2020'].values, order='F'))

newclim = newy[-n:]


plt.rcParams["font.family"] = 'Arial'
plt.rcParams["font.size"] = '10'

fig = plt.figure(figsize=(12,4))
plt.fill_between(np.arange(30,61), 19, 38, color='lightgrey')
plt.fill_between(x, newclim, temp, where=temp>newclim, interpolate=True, color='gold')
plt.fill_between(x, temp, newclim, where=temp<=newclim, interpolate=True, color='royalblue')
plt.plot(x, temp, color='k', lw=0.8)
plt.plot(x, newclim, color='grey', lw=0.8)
plt.plot(x, p90, color='r', lw=1.2, ls='--', label='90-percentile')
plt.xlim([-0,15+31+14])
plt.ylabel('Daily maximum temperature (ºC)')
plt.ylim([19,38])
plt.legend()
#plt.axvline(15+20, color='grey')
#plt.axvline(15+26, color='grey')
plt.xticks([0, 14, 30, 30+14, 30+14+16], ['2022-11-01', '2022-11-15', '2022-12-01', '2022-12-15', '2022-12-31'])
#plt.gca().tick_params(direction="in")
plt.savefig('../../../megafires_data/png/heat_wave_dec_2022.png', dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.show()



