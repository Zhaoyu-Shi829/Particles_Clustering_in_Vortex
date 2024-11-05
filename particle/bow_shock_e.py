import sys
import os
import time
import numpy as np
from scipy.io import FortranFile
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as patches
import seaborn as sns
import h5py
from mgtools import MGreadH5, ScalarField, VectorField, h5Viewer, vtkViewer
from scipy.interpolate import griddata
from matplotlib import rcParams

start = time.time()

dir_name = 'txt_cont_1'
base_name_0 = 'St3/ntot_d1'
base_name_1 = 'St3/pnts_d1_xy'
base_name_2 = 'St3/poly_pos_area'
base_name_3 = 'St3/poly_pos_area_up'
base_name_4 = 'St3/poly_pos_area_normal'

# PLOT AREA VS X (5 TIME STEPS SUPERPOSED)
poly_front_base = []
with open('txt_1/St3/poly_pos_area_normal.txt') as fu:
    poly_base = fu.readlines()
    newlist = []
    for line in poly_base:
        newlist.append(line.split())
    poly_front_base = [x for x in np.array(newlist) if -2.5 < float(x[0]) < -0.51 and -0.2 < float(x[1]) < 0.2]
print(len(poly_front_base))

# read and append poly_front from continuous runs
N = 10          # number of continuous files
dir_name_lst = ['txt_cont_{}'.format(i+1) for i in range(N)]
fu = ['fu{}'.format(i+1) for i in range(N)]
poly_cont = ['poly_cont_{}'.format(i+1) for i in range(N)]
poly_front = ['poly_front_{}'.format(i+1) for i in range(N)]
newlist = ['newlist_{}'.format(i+1) for i in range(N)]
for i in range(N):
    with open(os.path.join(dir_name_lst[i], base_name_4 + '.' + 'txt'), 'r') as fu[i]:
        poly_cont[i] = fu[i].readlines()
    newlist[i] = []
    for line in poly_cont[i]:
        newlist[i].append(line.split())
    poly_front[i] = [x for x in np.array(newlist[i]) if -2.5 < float(x[0]) < -0.51 and -0.1 < float(x[1]) < 0.1]
    poly_front_base.extend(poly_front[i])
poly_front_base = np.array(poly_front_base).astype(np.float)
print(poly_front_base.shape)
print('\n'.join(str(len(x)) for x in poly_front))
# descending x to be able to plot
id = np.lexsort((poly_front_base[:, 2], poly_front_base[:, 1], -1*poly_front_base[:, 0]))
poly_front_upx = poly_front_base[id]
# np.set_printoptions(threshold=np.inf)
# print(poly_front_upx)
# np.savetxt('./poly_front_upx.txt', poly_front_upx, fmt='%10.7e')

# Normal scatter plot(log or linear)
# plt.scatter(poly_front_upx[:, 0], poly_front_upx[:, 2], s=3, c='aqua')
# plt.yscale('log')
# plt.ylim(pow(10, -6), pow(10, -1))
# plt.show()


# Jointplot(voro area or normalised voro area)
# calculate slope and intercept ahead and behind the shock wave
poly_reg_seg_1 = []
poly_reg_seg_2 = []
poly_reg_seg_3 = []
for m in poly_front_upx:
    if -2.5 <= m[0] < -1.26:
        poly_reg_seg_1.append(m)
    if -1.26 <= m[0] < -1.17:
        poly_reg_seg_2.append(m)
    if -1.17 <= m[0] < -0.51:
        poly_reg_seg_3.append(m)
poly_reg_seg_1 = np.array(poly_reg_seg_1)
poly_reg_seg_2 = np.array(poly_reg_seg_2)
poly_reg_seg_3 = np.array(poly_reg_seg_3)
slope, intercept, r_value, p_value, str_err = stats.linregress(poly_reg_seg_3[:, 0], poly_reg_seg_3[:, 2])
print(slope, intercept)

# non-normalised
sns.set(font_scale=1.5, style='ticks', rc={'ytick.direction': 'in', 'xtick.direction': 'in', 'grid.linestyle': '--', 'grid.linewidth': 1})
g = sns.JointGrid(x=poly_front_upx[:, 0], y=poly_front_upx[:, 2], space=0)
g.plot_joint(plt.scatter, color='royalblue', alpha=0.7)
g.ax_joint.set_yscale('log')
g.set_axis_labels('$x(m)$', '$Voronoi\ area (m^2)$', fontsize=30, fontweight='bold')
plt.grid(True, which='both', linewidth=1.5)
# normalised (can not see the sharp shock interface)
# sns.set(font_scale=1.5, style='ticks', rc={'ytick.direction': 'in', 'xtick.direction': 'in', 'grid.linestyle': '--', 'grid.linewidth': 1})
# g = sns.JointGrid(x=poly_front_upx[:, 0], y=poly_front_upx[:, 2]/0.002, space=0)
# g.plot_joint(plt.scatter, color='royalblue', alpha=0.7)
# g.set_axis_labels('$x(m)$', '$V(m^2)$', fontsize=30, fontweight='bold')
# plt.grid(True, which='both')

g.plot_marginals(sns.distplot, bins=30, kde_kws={'color': 'red', 'lw': 3, 'label': 'kernel density'}, hist_kws={'color': 'royalblue'})
g.ax_marg_x.yaxis.set_ticks([0.5, 1, 1.5])
g.ax_marg_x.yaxis.set_ticks_position('left')
g.ax_marg_x.yaxis.set_ticklabels(('0.5', '1', '1.5'))
# g.ax_marg_x.set_ylim([0, 400])
g.ax_marg_x.grid('on', linewidth=1.5)
g.ax_marg_y.set_visible(False)

g.ax_joint.plot(poly_reg_seg_1[:, 0], poly_reg_seg_1[:, 0]*1.7082e-4 + 1.7278036e-3, 'red', linewidth=3.5)
g.ax_joint.plot(poly_reg_seg_3[:, 0], poly_reg_seg_3[:, 0]*2.4637e-4 + 1.1026555e-3, 'red', linewidth=3.5)
g.ax_joint.annotate('$y=1.708e^{-4}x+1.728e^{-3}$', xy=(-2.0, 2e-3), xytext=(-2.0, 5e-3), fontSize=20, arrowprops=dict(facecolor='black'),
                    horizontalalignment='center', verticalalignment='center')
g.ax_joint.annotate('$y=2.464e^{-4}x+1.103e^{-3}$', xy=(-0.8, 1.5e-3), xytext=(-0.8, 4.5e-3), fontSize=20, arrowprops=dict(facecolor='black'),
                    horizontalalignment='center', verticalalignment='center')
g.ax_joint.add_patch(patches.Rectangle((-1.92, 4.2e-3), 0.22, 0.002, linewidth=2.5, linestyle='--', edgecolor='red', facecolor='none'))
g.ax_joint.add_patch(patches.Rectangle((-0.72, 3.8e-3), 0.22, 0.002, linewidth=2.5, linestyle='--', edgecolor='red', facecolor='none'))
plt.ylim([pow(10, -5), pow(10, -2)])
plt.legend()
plt.show()


end = time.time()
print(end - start)