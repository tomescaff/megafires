import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import sys
from scipy.stats import norm

sys.path.append('../../processing')

import processing.lens as lens
import processing.stations as stns

lens1 = lens.get_LENS_jan_tmax_QNWE()
lens2 = lens.get_LENS2_jan_tmax_QNWE()

# get Quinta Normal time series
qn = stns.get_QN_tmax_jan()
qn_mean_1950_1970 = qn.sel(time=slice('1950', '1970')).mean('time')

lens1_mean_1950_1970 = lens1.sel(time=slice('1950', '1970')).mean('time')
lens2_mean_1950_1970 = lens2.sel(time=slice('1950', '1970')).mean('time')

lens1_anom = lens1 - lens1_mean_1950_1970
lens2_anom = lens2 - lens2_mean_1950_1970

lens1_fixed_1950_1970 = lens1_anom + qn_mean_1950_1970
lens2_fixed_1950_1970 = lens2_anom + qn_mean_1950_1970

lens1_fixed_1950_1970_ensmean = lens1_fixed_1950_1970.mean('run')
lens2_fixed_1950_1970_ensmean = lens2_fixed_1950_1970.mean('run')

# 1925-1960
qn_mean_1925_1960 = qn.sel(time=slice('1925', '1960')).mean('time')

lens1_mean_1925_1960 = lens1.sel(time=slice('1925', '1960')).mean('time')
lens2_mean_1925_1960 = lens2.sel(time=slice('1925', '1960')).mean('time')

lens1_anom = lens1 - lens1_mean_1925_1960
lens2_anom = lens2 - lens2_mean_1925_1960

lens1_fixed_1925_1960 = lens1_anom + qn_mean_1925_1960
lens2_fixed_1925_1960 = lens2_anom + qn_mean_1925_1960

lens1_fixed_1925_1960_ensmean = lens1_fixed_1925_1960.mean('run')
lens2_fixed_1925_1960_ensmean = lens2_fixed_1925_1960.mean('run')

init = '1925'
end = '1960'
normfit_qn = norm.fit(np.ravel(qn.sel(time=slice(init, end)).values))
normfit_lens1_1950_1970 = norm.fit(np.ravel(lens1_fixed_1950_1970.sel(time=slice(init, end)).values))
normfit_lens1_1925_1960 = norm.fit(np.ravel(lens1_fixed_1925_1960.sel(time=slice(init, end)).values))
normfit_lens2_1950_1970 = norm.fit(np.ravel(lens2_fixed_1950_1970.sel(time=slice(init, end)).values))
normfit_lens2_1925_1960 = norm.fit(np.ravel(lens2_fixed_1925_1960.sel(time=slice(init, end)).values))

temp_ee = qn.sel(time='2017')
xmin = 20
xmax = 36
x = np.linspace(xmin, xmax, 100)

# create figure
fig = plt.figure(figsize=(12,7))
plt.plot(x, norm.pdf(x, *normfit_qn), 'green', linewidth=2, label = 'QN')
plt.plot(x, norm.pdf(x, *normfit_lens1_1950_1970), 'b', linewidth=1.2, ls='--',label = 'LENS1 fixed 1950-1970')
plt.plot(x, norm.pdf(x, *normfit_lens1_1925_1960), 'b', linewidth=1.2, label = 'LENS1 fixed 1925-1960')
plt.plot(x, norm.pdf(x, *normfit_lens2_1950_1970), 'r', linewidth=1.2, ls='--',label = 'LENS2 fixed 1950-1970')
plt.plot(x, norm.pdf(x, *normfit_lens2_1925_1960), 'r', linewidth=1.2, label = 'LENS2 fixed 1925-1960')
plt.axvline(temp_ee, color='grey', linewidth=1.2, label = 'EE Temperature')
plt.xlim([xmin,xmax])
plt.ylim([0,0.5])
plt.grid(lw=0.2, ls='--', color='grey')
# set legend
ax = plt.gca()
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.01,
                 box.width, box.height * 0.99])

# Put a legend below current axis
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
          fancybox=True, shadow=True, ncol=6)
# set title and labels
plt.title('Tmax distribution for period '+init+'-'+end)
plt.xlabel('January Tmax (ºC)')
plt.ylabel('PDF')
plt.tight_layout()
plt.savefig('../../../megafires_data/png/LENS_LENS2_bias_correction_all_runs_QNWE_pdf_'+init+'_'+end+'.png', dpi=300)
plt.show()