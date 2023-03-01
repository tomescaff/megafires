import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from os.path import join

relpath = '../../../megafires_data/output/'
obs_mle_1880 = pd.read_csv(join(relpath, 'metrics_RI_MLE_2017_2070_normfit_1960_2021_nboot_10.csv'), index_col=0)
lens1_mle_cr = pd.read_csv(join(relpath, 'metrics_LENS1_MLE_2017_2070_normfit_2010_2100_by_return_period_RI_nboot_10.csv'), index_col=0)
lens2_mle_1880 = pd.read_csv(join(relpath, 'metrics_LENS2_MLE_2017_2070_normfit_2010_2100_by_return_period_RI_nboot_10.csv'), index_col=0)

models = [ 
            obs_mle_1880, 
            lens1_mle_cr,
            lens2_mle_1880, 
         ]
model_names = [
            'CR2MET',
            'CESM1-LENS',
            'CESM2-LENS']

center_rr = np.array([x.loc['rr a-f', 'raw'] for x in models])
lower_rr = np.array([x.loc['rr a-f', '95ci lower'] for x in models])
upper_rr = np.array([x.loc['rr a-f', '95ci upper'] for x in models])
width_rr = upper_rr - lower_rr

center_ac = np.array([x.loc['tau ac', 'raw'] for x in models])
lower_ac = np.array([x.loc['tau ac', '95ci lower'] for x in models])
upper_ac = np.array([x.loc['tau ac', '95ci upper'] for x in models])
width_ac = upper_ac - lower_ac

center_fu = np.array([x.loc['tau fu', 'raw'] for x in models])
lower_fu = np.array([x.loc['tau fu', '95ci lower'] for x in models])
upper_fu = np.array([x.loc['tau fu', '95ci upper'] for x in models])
width_fu = upper_fu - lower_fu

center_ma = -1*np.array([x.loc['delta a-f', 'raw'] for x in models])
lower_ma = -1*np.array([x.loc['delta a-f', '95ci upper'] for x in models])
upper_ma = -1*np.array([x.loc['delta a-f', '95ci lower'] for x in models])
width_ma = upper_ma - lower_ma

fig, axs = plt.subplots(4,1, figsize=(12,7))
plt.rcParams["font.family"] = 'Arial'
plt.rcParams["font.size"] = 10
y_pos = np.arange(len(model_names))

ax=axs[0]
plt.sca(ax)
plt.grid(color='grey', lw=0.8, ls='--', axis='x', zorder=0)
barlist = ax.barh(y_pos, width=width_ac, left=lower_ac, height=0.4, align='center', zorder=3)
colors=['blue', 'red', 'fuchsia']
for bar, color in zip(barlist, colors):
    bar.set_color(color)
plt.scatter(center_ac, y_pos, s=200, marker='|', color='k', zorder=4)
plt.yticks(np.append(y_pos, [3, 4]), model_names + ['EC_Earth3', 'CMIP6'])
plt.ylim([-0.5, 4.5])
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel('Factual return period (yr)')
ax.set_xscale('log')
ax.spines.right.set_visible(False)
ax.spines.left.set_visible(False)
ax.spines.top.set_visible(False)
ax.tick_params(direction="in")
plt.xlim([0.1, 1e9])
plt.xticks([0.1, 1.0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9])
# plt.axvline(1, color='k', lw=1.0)


ax = axs[1]
plt.sca(ax)
plt.grid(color='grey', lw=0.8, ls='--', axis='x', zorder=0)
barlist = ax.barh(y_pos, width=width_fu, left=lower_fu, height=0.4, align='center', zorder=3)
for bar, color in zip(barlist, colors):
    bar.set_color(color)
plt.scatter(center_fu, y_pos, s=200, marker='|', color='k', zorder=4)
plt.yticks(np.append(y_pos, [3, 4]), model_names + ['EC_Earth3', 'CMIP6'])
plt.ylim([-0.5, 4.5])
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel('Future return period (yr)')
ax.set_xscale('log')
ax.spines.right.set_visible(False)
ax.spines.left.set_visible(False)
ax.spines.top.set_visible(False)
plt.xlim([0.1, 1e9])
plt.xticks([0.1, 1.0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9])
# plt.axvline(1, color='k', lw=1.0)


ax = axs[2]
plt.sca(ax)
plt.grid(color='grey', lw=0.8, ls='--', axis='x', zorder=0)
barlist = ax.barh(y_pos, width=width_rr, left=lower_rr, height=0.4, align='center', zorder=3)
for bar, color in zip(barlist, colors):
    bar.set_color(color)
plt.scatter(center_rr, y_pos, s=200, marker='|', color='k', zorder=4)
plt.yticks(np.append(y_pos, [3, 4]), model_names + ['EC_Earth3', 'CMIP6'])
plt.ylim([-0.5, 4.5])
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel('Probability ratio')
ax.set_xscale('log')
ax.spines.right.set_visible(False)
ax.spines.left.set_visible(False)
ax.spines.top.set_visible(False)
plt.xlim([0.1, 1e9])
plt.xticks([0.1, 1.0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9])
plt.axvline(1, color='k', lw=1.0)

ax = axs[3]
plt.sca(ax)
plt.grid(color='grey', lw=0.8, ls='--', axis='x', zorder=0)
barlist = ax.barh(y_pos, width=width_ma, left=lower_ma, height=0.4, align='center', zorder=3)
for bar, color in zip(barlist, colors):
    bar.set_color(color)
plt.scatter(center_ma, y_pos, s=200, marker='|', color='k', zorder=4)
plt.yticks(np.append(y_pos, [3, 4]), model_names + ['EC_Earth3', 'CMIP6'])
plt.ylim([-0.5, 4.5])
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel('Intensity change (ºC)')
ax.spines.right.set_visible(False)
ax.spines.left.set_visible(False)
ax.spines.top.set_visible(False)
plt.xlim([-4.0, 4.0])
plt.xticks([-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0])
plt.axvline(0, color='k', lw=1.0)

plt.tight_layout()
plt.savefig('../../../megafires_data/png/synthesis_plot_metrics_normfit_MLE_RI_future.png', dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.show()




