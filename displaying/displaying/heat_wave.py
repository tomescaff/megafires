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

temp_dic = df.loc[:,'2016'].iloc[-15:]
temp_enefeb = df.loc[:,'2017'].iloc[:46]

clim_dic = df.loc[:,'clim'].iloc[-15:]
clim_enefeb = df.loc[:,'clim'].iloc[:46]

p90_dic = df_umb.loc[:,'Dic'].iloc[-15:]
p90_ene = df_umb.loc[:,'Ene']
p90_feb = df_umb.loc[:,'Feb'].iloc[:15]

temp = np.concatenate((temp_dic.values, temp_enefeb.values))
clim = np.concatenate((clim_dic.values, clim_enefeb.values))
p90 = np.concatenate((p90_dic.values, p90_ene.values, p90_feb.values))
x = np.arange(15+31+15)


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

newclim = np.concatenate((newy[-15:], newy[:46]))


plt.rcParams["font.family"] = 'Arial'
plt.rcParams["font.size"] = '10'

fig = plt.figure(figsize=(12,4))
plt.fill_between(np.arange(15,30+18), 23, 38, color='lightgrey')
plt.fill_between(x, newclim, temp, where=temp>newclim, interpolate=True, color='gold')
plt.fill_between(x, temp, newclim, where=temp<=newclim, interpolate=True, color='royalblue')
plt.plot(x, temp, color='k', lw=0.8)
plt.plot(x, newclim, color='grey', lw=0.8)
plt.plot(x, p90, color='r', lw=1.2, ls='--', label='90-percentile')
plt.xlim([-0,15+31+14])
plt.ylabel('Daily maximum temperature (ºC)')
plt.ylim([23,38])
plt.legend()
#plt.axvline(15+20, color='grey')
#plt.axvline(15+26, color='grey')
plt.xticks([0, 15, 30, 30+17, 15+31+14], ['2016-12-17', '2017-01-01', '2017-01-15', '2017-02-01', '2017-02-15'])
#plt.gca().tick_params(direction="in")
# plt.savefig('../../../megafires_data/png/heat_wave.png', dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.show()



