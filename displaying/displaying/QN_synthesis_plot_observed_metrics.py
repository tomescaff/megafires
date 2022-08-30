import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from os.path import join

relpath = '../../../megafires_data/output/'

moa_30y = pd.read_csv(join(relpath, 'metrics_of_attr_QN_30y.csv'), index_col=0)
moa_50y = pd.read_csv(join(relpath, 'metrics_of_attr_QN_50y.csv'), index_col=0)
moa_mle_1943 = pd.read_csv(join(relpath, 'metrics_of_attr_QN_MLE_1943_2017.csv'), index_col=0)
moa_30y_gum = pd.read_csv(join(relpath, 'metrics_of_attr_QN_30y_gumbel.csv'), index_col=0)
moa_50y_gum = pd.read_csv(join(relpath, 'metrics_of_attr_QN_50y_gumbel.csv'), index_col=0)
moa_mle_1943_gum = pd.read_csv(join(relpath, 'metrics_of_attr_QN_MLE_gumbel_1943_2017.csv'), index_col=0)
moa_mle_1880 = pd.read_csv(join(relpath, 'metrics_of_attr_QN_MLE_1880_2017.csv'), index_col=0)
moa_mle_1880_gum = pd.read_csv(join(relpath, 'metrics_of_attr_QN_MLE_gumbel_1880_2017.csv'), index_col=0)

models = [  moa_30y, 
            moa_50y, 
            moa_mle_1943, 
            moa_30y_gum, 
            moa_50y_gum, 
            moa_mle_1943_gum, 
            moa_mle_1880, 
            moa_mle_1880_gum]
model_names = (  
            '30yr - Normal (1943 vs. 2017)', 
            '50yr - Normal (1943 vs. 2017)', 
            'MLE - Normal (1943 vs. 2017)',
            '30yr - Gumbel (1943 vs. 2017)',
            '50yr - Gumbel (1943 vs. 2017)',
            'MLE - Gumbel (1943 vs. 2017)',
            'MLE - Normal (1880 vs. 2017)',
            'MLE - Gumbel (1880 vs. 2017)' )

center = np.array([x.loc['rr c-a', 'raw'] for x in models])
lower = np.array([x.loc['rr c-a', '95ci lower'] for x in models])
upper = np.array([x.loc['rr c-a', '95ci upper'] for x in models])
width = upper - lower

# rr
fig = plt.figure(figsize=(12,4))
plt.rcParams["font.family"] = 'Arial'
plt.rcParams["font.size"] = 10
y_pos = np.arange(len(model_names))
ax = plt.gca()
barlist = ax.barh(y_pos, width=width, left=lower, height=0.4, align='center')
colors=['blue', 'blue', 'blue', 'red', 'red', 'red', 'green', 'green']
for bar, color in zip(barlist, colors):
    bar.set_color(color)
plt.scatter(center, y_pos, s=200, marker='|', color='k', zorder=4)
plt.yticks(y_pos, model_names)
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel('Probability ratio')
ax.set_xscale('log')
ax.spines.right.set_visible(False)
ax.spines.left.set_visible(False)
ax.spines.top.set_visible(False)
plt.xlim([0.001, 100000])
plt.tight_layout()
plt.axvline(1, color='k', lw=0.4, ls='--')
plt.savefig('../../../megafires_data/png/QN_sysnthesis_plot_observed_metrics_rr.png', dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.show()

center = np.array([x.loc['tau ac', 'raw'] for x in models])
lower = np.array([x.loc['tau ac', '95ci lower'] for x in models])
upper = np.array([x.loc['tau ac', '95ci upper'] for x in models])
width = upper - lower

# return periods
fig = plt.figure(figsize=(12,4))
plt.rcParams["font.family"] = 'Arial'
plt.rcParams["font.size"] = 10
y_pos = np.arange(len(model_names))
ax = plt.gca()
barlist = ax.barh(y_pos, width=width, left=lower, height=0.4, align='center')
colors=['blue', 'blue', 'blue', 'red', 'red', 'red', 'green', 'green']
for bar, color in zip(barlist, colors):
    bar.set_color(color)
plt.scatter(center, y_pos, s=200, marker='|', color='k', zorder=4)
plt.yticks(y_pos, model_names)
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel('Return period (yr)')
ax.set_xscale('log')
ax.spines.right.set_visible(False)
ax.spines.left.set_visible(False)
ax.spines.top.set_visible(False)
plt.xlim([0.001, 100000])
plt.tight_layout()
plt.axvline(1, color='k', lw=0.4, ls='--')
plt.savefig('../../../megafires_data/png/QN_sysnthesis_plot_observed_metrics_tau.png', dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.show()