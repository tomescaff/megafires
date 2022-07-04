import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'Arial'
import sys

sys.path.append('../../processing')

import processing.stations as stns
from scipy.stats import norm
from scipy.stats import genextreme
from scipy.stats import gumbel_r, gumbel_l
from scipy import stats
from scipy.stats import kstest
import numpy as np

# get Quinta Normal time series
jan = stns.get_QN_tmax_jan().sel(time=slice('1950','2021')).values
feb = stns.get_QN_tmax_feb().sel(time=slice('1950','2021')).values
janfeb = np.append(jan, feb)

# fit normal dist
normfit_jan = norm.fit(jan)
normfit_feb = norm.fit(feb)
normfit_janfeb = norm.fit(janfeb)

# # test fit
stat, p_jan = kstest(jan, 'norm', normfit_jan)
stat, p_feb = kstest(feb, 'norm', normfit_feb)
stat, p_janfeb = kstest(janfeb, 'norm', normfit_janfeb)

xmin, xmax = 26, 34
ymin, ymax = 0, 0.5
# create figure
fig, axs = plt.subplots(2, 2, figsize=(9,7.5), constrained_layout=True)

plt.sca(axs[0,0])
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
x = np.linspace(xmin, xmax, 100)
plt.plot(x, norm.pdf(x, *normfit_jan), 'k', linewidth=2, label = 'Jan normfit (p = {:0.2f})'.format(p_jan))
plt.plot(x, norm.pdf(x, *normfit_feb), 'g', linewidth=2, label = 'Feb normfit (p = {:0.2f})'.format(p_feb))
plt.plot(x, norm.pdf(x, *normfit_janfeb), 'r', linewidth=2, label = 'Jan-Feb normfit (p = {:0.2f})'.format(p_janfeb))
plt.grid(lw=0.2, ls='--', color='grey')
plt.legend(prop={'size': 8})
plt.xlabel('Monthly Tmax (ºC)')
plt.ylabel('PDF')

plt.sca(axs[0,1])
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
hist, bins = np.histogram(jan, bins=np.arange(xmin,xmax+0.5, 0.5), density=True)
width = 0.9 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
plt.bar(center, hist, align='center', width=width, edgecolor='k', facecolor='grey', color='k', alpha = 0.75, label = 'QN Jan Tmax')
plt.grid(lw=0.2, ls='--', color='grey')
plt.legend(prop={'size': 8})
plt.xlabel('Monthly Tmax (ºC)')
plt.ylabel('PDF')

plt.sca(axs[1,1])
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
hist, bins = np.histogram(feb, bins=np.arange(xmin,xmax+0.5, 0.5), density=True)
width = 0.9 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
plt.bar(center, hist, align='center', width=width, edgecolor='k', facecolor='green', color='k', alpha = 0.75, label = 'QN Feb Tmax')
plt.grid(lw=0.2, ls='--', color='grey')
plt.legend(prop={'size': 8})
plt.xlabel('Monthly Tmax (ºC)')
plt.ylabel('PDF')

plt.sca(axs[1,0])
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
hist, bins = np.histogram(janfeb, bins=np.arange(xmin,xmax+0.5, 0.5), density=True)
width = 0.9 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
plt.bar(center, hist, align='center', width=width, edgecolor='k', facecolor='red', color='k', alpha = 0.75, label = 'QN Jan-Feb Tmax')
plt.grid(lw=0.2, ls='--', color='grey')
plt.legend(prop={'size': 8})
plt.xlabel('Monthly Tmax (ºC)')
plt.ylabel('PDF')

for ax in axs.ravel():
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(9) 
        tick.label.set_weight('light') 
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(9) 
        tick.label.set_weight('light')

plt.savefig('../../../megafires_data/png/QN_pdf_jan_feb.png', dpi=300)
plt.show()
    

