import sys
import matplotlib.pyplot as plt
from os.path import join, abspath, dirname

currentdir = dirname(abspath(__file__))
sys.path.append(join(currentdir, '../../processing'))

from processing import gmst

jan = gmst.get_gmst_jan()
yea = gmst.get_gmst_annual()
smt = gmst.get_gmst_annual_5year_smooth()

fig = plt.figure(figsize=(10,6))
plt.rcParams["font.family"] = 'Arial'
plt.rcParams["font.size"] = '10'
plt.plot(smt.time.dt.year, smt.values, linewidth=2.0, color='k')
plt.plot(jan.time.dt.year, jan.values, linewidth=1.0, color='r')
plt.plot(yea.time.dt.year, yea.values, linewidth=1.0, color='b')
plt.xlabel('Year')
plt.ylabel('GMST (ºC)')
plt.legend(['Lowess (5)', 'January', 'Annual (J-D)'])
plt.grid(color='grey', linestyle='--', linewidth=0.4)
ax = plt.gca()
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
ax.tick_params(direction="in")
filepath = join(currentdir, '../../../megafires_data/png/comparison_GMST.png')
plt.savefig(filepath, dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.show()