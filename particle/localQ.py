import sys
import os
import time
import numpy as np
from scipy.io import FortranFile
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1 import make_axes_locatable
import h5py
import ast
from mgtools import MGreadH5, ScalarField, VectorField, h5Viewer, vtkViewer
from scipy.interpolate import RegularGridInterpolator
from matplotlib import rcParams
import matplotlib.tri as tri

start = time.time()
dir_name = 'txt_cont_0/St1'
base_name_0 = 'Sk1/ntot_d1'
base_name_1 = 'Sk1/pnts_d1_xy'
base_name_2 = 'Sk1/poly_pos_area'
base_name_3 = 'Sk1/poly_pos_area_up'
base_name_4 = 'Sk1/poly_pos_area_normal'
'''
# Read flow Q at the first time step in the intermediate file generated from saveh5.py
Nstart = 0
Nstop = 439

x = np.linspace(-16.32, 16.32, 1020 + 1)
y = np.linspace(-8.0, 8.0, 500 + 1) 
z = np.linspace(10.0, 14.032, 126 + 1)
qvortex = ScalarField(x=x, y=y, z=z)

# extract flow_Q and xyz
f1 = '/media/zhaoyus/Zhaoyu/Trial_case/Par_codes_trial/Re100/IJMF_DATA/Re100_refine/Sticky/long_0.016/St1/par_10k/par_10k_St1.h5'
viewer = h5Viewer(f1)
for grid in range(Nstart, Nstop):
    tic = time.perf_counter()
    viewer.set_group('/grid-{}'.format(grid))
    qvortexGrid = viewer.read("Q")
    qvortex.merge(qvortexGrid, strip=1)

viewer.close()

# save flow_Q into new h5 file
f2 = '/media/zhaoyus/Zhaoyu/Trial_case/Par_codes_trial/Re100/IJMF_DATA/Re100_refine/Sticky/long_0.016/St1/par_10k'
fq = h5py.File(os.path.join(f2, 'flow_Q.h5'), 'w')
fq.create_dataset('flow_Q', data=qvortex.data, dtype='<f4')
fq.create_dataset('x', data=qvortex.x, dtype='<f4')
fq.create_dataset('y', data=qvortex.y, dtype='<f4')
fq.create_dataset('z', data=qvortex.z, dtype='<f4')
fq.close()
'''

# read projected particle position (2D)
with open(os.path.join(dir_name, base_name_1 + '.' + 'txt'), 'r') as fpos2:
    pxy = fpos2.readlines()
    ppos_2d = []
    for line in pxy:
        line = eval(line)
        ppos_2d.append(line)

ppos_2d_x = []
ppos_2d_y = []
for item in ppos_2d:
    item = np.asarray(item)
    item = ' '.join(item)
    item = np.fromstring(item, dtype=np.float32, sep=' ')
    ppos_2d_x.append(item[0])
    ppos_2d_y.append(item[1])
ppos_2d_x = np.asarray(ppos_2d_x)
ppos_2d_y = np.asarray(ppos_2d_y)

# read flow_Q of the whole 3D domain
# for laminar flow, just take one slice at z=const=N
f3 = '/media/zhaoyus/Zhaoyu/Trial_case/Par_codes_trial/Re100/IJMF_DATA/Re100_refine/Sticky/long_0.016/St1/par_10k/flow_Q.h5'
with h5py.File(f3, 'r') as ffQ:
    flowQ = ffQ.get('flow_Q')[()]
    x = ffQ.get('x')[()]
    y = ffQ.get('y')[()]
    z = ffQ.get('z')[()]

# reshape flowQ(values) and ppos(points) into griddata type
N = 12
# flowQ_arr = []
# for i in range(len(x)):
#     for j in range(len(y)):
#         flowQ_arr.append(flowQ[i, j, N])
# print(max(flowQ_arr), min(flowQ_arr))
# flowQ_arr = np.asarray(flowQ_arr, dtype=np.float32)

# points = []
# for i in range(len(x)):
#     grid_xy = np.zeros((len(y), 2))
#     for j in range(len(y)):
#         grid_xy[j, 0] = x[i]
#         grid_xy[j, 1] = y[j]
#     points.extend(list(grid_xy))
#     # print(*points, sep='\n')
# points_arr = np.asarray(points, dtype=np.float32)

'''
# plot flow_Q and particle distribution
ax = plt.gca()
plt.xlabel('x(m)', fontsize=22)
plt.ylabel('y(m)', fontsize=22)
plt.tick_params(labelsize=20)
plt.imshow(flowQ[:, :, N].T, extent=(-16.32, 16.32, -8, 8), origin='lower', cmap='RdBu_r')
ax.scatter(ppos_2d_x, ppos_2d_y, marker='.', s=1, color='black')
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='1.5%', pad=0.05)
cbar = plt.colorbar(extend='both', cax=cax)
cbar.set_label(label='Flow_Q', size=17, weight='bold')
cbar.ax.tick_params(labelsize=17)
plt.clim(-0.05, 0.05)

figure = plt.gcf()
figure.set_size_inches(18, 12)
plt.show()
'''

# interpolate flowQ into local flowQ
intp_fun = RegularGridInterpolator((x, y), flowQ[:, :, N], method='linear')

pts = list(zip(ppos_2d_x, ppos_2d_y))
pts_f = []
for line in pts:
    arr = np.asarray(line)
    if (-2 < arr[0] < 15 and -3 < arr[1] < 3):
        pts_f.append(arr)
pts_f = np.array(pts_f)
x = np.array(pts_f[:, 0], dtype=np.float16)
y = np.array(pts_f[:, 1], dtype=np.float16)
localQ = intp_fun(pts_f)
localQ = localQ.astype('float16')

plt.tricontour(x, y, localQ, levels=np.linspace(-1, 1, 50))
plt.tricontourf(x, y, localQ, levels=np.linspace(-1, 1, 50))

circle = plt.Circle((0, 0), 0.5, fill=True, color='gray')
plt.gcf().gca().add_artist(circle)
plt.axis('equal')
plt.colorbar()

plt.show()

# np.set_printoptions(threshold=np.inf)
end = time.time()
print(end-start)