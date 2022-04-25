import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'Arial'
import sys

sys.path.append('../../processing')

import processing.utils as ut
from scipy.stats import norm
from scipy.stats import genextreme
from scipy.stats import gumbel_r, gumbel_l
from scipy import stats
from scipy.stats import kstest
import numpy as np

# get Quinta Normal time series
da = ut.get_QN_series_detrended()

# fit normal dist
normfit = norm.fit(da.values)

# fit gev
gevfit = genextreme.fit(da.values)

# fit gumbel_r
gumrfit = gumbel_r.fit(da.values)

# fit gumbel_l
gumlfit = gumbel_l.fit(da.values)

# test fit
stat, p_norm = kstest(da.values, 'norm', normfit)
stat, p_gev = kstest(da.values, 'genextreme', gevfit)
stat, p_gumr = kstest(da.values, 'gumbel_r', gumrfit)


# compute histogram
hist, bins = np.histogram(da.values, bins=np.arange(26,34+2, 0.5), density=True)

# create figure
fig = plt.figure(figsize=(8,6))

# plot the histogram
width = 0.9 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
plt.bar(center, hist, align='center', width=width, edgecolor='orange', facecolor='navajowhite', color='orange', alpha = 1, label = 'QN Jan Tmax')

# plot the PDF
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)

plt.plot(x, norm.pdf(x, *normfit), 'k', linewidth=2, label = 'Normal fit (p = {:g})'.format(p_norm))
plt.plot(x, genextreme.pdf(x, *gevfit), 'g', linewidth=2, label = 'GEV fit (p = {:g})'.format(p_gev))
plt.plot(x, gumbel_r.pdf(x, *gumrfit), 'r', linewidth=2, label = 'GUM_r fit (p = {:g})'.format(p_gumr))

# set grid
plt.grid(lw=0.2, ls='--', color='grey')

# set legend
plt.legend()

# set xlim
plt.xlim([26,34])

# set title and labels
plt.xlabel('January Tmax (ºC)')
plt.ylabel('PDF')
plt.title('January Tmax distribution at Quinta Normal')
plt.savefig('../../../megafires_data/png/QN_pdf_detrended.png', dpi=300)
plt.show()
    

