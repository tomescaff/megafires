import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from os.path import join

relpath = '../../../megafires_data/output/'

filename = 'MLE_tasmax_jan_QN_GMST_1000_normal_validation.nc'
qn = xr.open_dataset(join(relpath, filename))

filename = 'MLE_tasmax_jan_LENS1_GMST_100_normal_validation_QN_NN.nc'
lens1 = xr.open_dataset(join(relpath, filename))

filename = 'MLE_tasmax_jan_LENS2_GMST_100_normal_validation_QN_NN.nc'
lens2 = xr.open_dataset(join(relpath, filename))

model_names = ['Observation', 'CESM1-LENS', 'CESM2-LENS', 'EC_Earth3', 'CMIP6']
models = [qn, lens1, lens2, lens2, lens2]

fig, axs = plt.subplots(3,1,figsize=(8,7))
plt.rcParams["font.family"] = 'Arial'
plt.rcParams["font.size"] = 10

# mu
varname = 'mu0'
center = [np.quantile(m[varname].values, 0.5, axis=0) for m in models]
lower  = [np.quantile(m[varname].values, 0.025, axis=0) for m in models]
upper  = [np.quantile(m[varname].values, 0.975, axis=0) for m in models]
width  = np.array(upper) - np.array(lower)

plt.sca(axs[0])
y_pos = np.arange(len(model_names))
ax = plt.gca()
barlist = ax.barh(y_pos, width=width, left=lower, height=0.4, align='center')
colors=['blue', 'blue', 'blue', 'none', 'none']
for bar, color in zip(barlist, colors):
    bar.set_color(color)
plt.scatter(center, y_pos, s=200, marker='|', color=['k', 'k', 'k', 'none', 'none'], zorder=4)
plt.yticks(y_pos, model_names)
plt.xlim([28,30])
ax.set_axisbelow(True)
plt.grid( lw=0.4, ls='--', color='grey', zorder=-4)
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel(varname)
ax.spines.right.set_visible(False)
ax.tick_params(direction="in")
ax.spines.top.set_visible(False)

# sigma
varname = 'sigma0'
center = [np.quantile(m[varname].values, 0.5, axis=0) for m in models]
lower  = [np.quantile(m[varname].values, 0.025, axis=0) for m in models]
upper  = [np.quantile(m[varname].values, 0.975, axis=0) for m in models]
width  = np.array(upper) - np.array(lower)

plt.sca(axs[1])
y_pos = np.arange(len(model_names))
ax = plt.gca()
barlist = ax.barh(y_pos, width=width, left=lower, height=0.4, align='center')
colors=['blue', 'blue', 'blue', 'none', 'none']
for bar, color in zip(barlist, colors):
    bar.set_color(color)
plt.scatter(center, y_pos, s=200, marker='|', color=['k', 'k', 'k', 'none', 'none'], zorder=4)
plt.yticks(y_pos, model_names)
plt.xlim([0.0,2.0])
ax.set_axisbelow(True)
plt.grid( lw=0.4, ls='--', color='grey', zorder=-4)
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel(varname)
ax.spines.right.set_visible(False)
ax.tick_params(direction="in")
ax.spines.top.set_visible(False)

# alpha
varname = 'alpha'
center = [np.quantile(m[varname].values, 0.5, axis=0) for m in models]
lower  = [np.quantile(m[varname].values, 0.025, axis=0) for m in models]
upper  = [np.quantile(m[varname].values, 0.975, axis=0) for m in models]
width  = np.array(upper) - np.array(lower)

plt.sca(axs[2])
y_pos = np.arange(len(model_names))
ax = plt.gca()
barlist = ax.barh(y_pos, width=width, left=lower, height=0.4, align='center')
colors=['blue', 'blue', 'blue', 'none', 'none']
for bar, color in zip(barlist, colors):
    bar.set_color(color)
plt.scatter(center, y_pos, s=200, marker='|', color=['k', 'k', 'k', 'none', 'none'], zorder=4)
plt.yticks(y_pos, model_names)
plt.xlim([0.5, 2.5])
ax.set_axisbelow(True)
ax.invert_yaxis()  # labels read top-to-bottom
plt.grid( lw=0.4, ls='--', color='grey', zorder=-4)
ax.set_xlabel(varname)
ax.spines.right.set_visible(False)
ax.tick_params(direction="in")
ax.spines.top.set_visible(False)

plt.tight_layout()
plt.savefig('../../../megafires_data/png/QN_all_models_mle_validation_NN.png', dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.show()


