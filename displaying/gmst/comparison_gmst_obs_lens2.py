import sys
import matplotlib.pyplot as plt
from os.path import join, abspath, dirname

currentdir = dirname(abspath(__file__))
sys.path.append(join(currentdir, '../../processing'))

from processing import gmst

lens2 = gmst.get_gmst_annual_lens2_ensmean()
smt = gmst.get_gmst_annual_5year_smooth()

fig = plt.figure(figsize=(10,6))
plt.rcParams["font.family"] = 'Arial'
plt.rcParams["font.size"] = '10'
plt.plot(smt.time.dt.year, smt.values, linewidth=2.0, color='b')
plt.plot(lens2.time.dt.year, lens2.values, linewidth=2.0, color='r')
plt.xlabel('Year')
plt.ylabel('GMST (ºC)')
plt.legend(['NASA GMST - Lowess (5)', 'LENS2 GMST - ensmean'])
plt.grid(color='grey', linestyle='--', linewidth=0.4)
ax = plt.gca()
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
ax.tick_params(direction="in")
filepath = join(currentdir, '../../../megafires_data/png/comparison_GMST_obs_lens2.png')
plt.savefig(filepath, dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.show()