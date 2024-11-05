import sys
import os
import numpy as np
from scipy.io import FortranFile
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d, ConvexHull, convex_hull_plot_2d
from shapely.geometry import Polygon
from ast import literal_eval
import ast
import time
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.tri as tri
import matplotlib.colors as colors
from matplotlib import rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns


filename = 'fort.22'
itype='<i4'
rtype='<f4'
htype='<u4'

part_type = [('XPART_X', '<f4'),
             ('XPART_Y', '<f4'),
             ('XPART_Z', '<f4'),
             ('START_VEC_X', '<f4'),
             ('START_VEC_Y', '<f4'),
             ('START_VEC_Z', '<f4'),
             ('DS', '<f4'),
             ('PART_IDX', '<i4'),
             ('DIAMCLASS', '<f4'),
             ('UPART_X', '<f4'),
             ('UPART_Y', '<f4'),
             ('UPART_Z', '<f4'),
             ('OLD_UPART_X', '<f4'),
             ('OLD_UPART_Y', '<f4'),
             ('OLD_UPART_Z', '<f4')]

fp = open(filename, 'rb')
f = FortranFile(fp, header_dtype=htype)

[idtag] = f.read_ints(itype)
print("IDTAG:         ", idtag)

[version] = f.read_ints(itype)
print("VERSION:       ", version)

[npart] = f.read_ints(itype)
print("NPART:         ", npart)

[rand_walk_max] = f.read_ints(itype)
print("RAND_WALK_MAX: ", rand_walk_max)

[ntotgrids] = f.read_ints(itype)
print("NTOTGRIDS:     ", ntotgrids)

[nsections] = f.read_ints(itype)
print("NSECTIONS:     ", nsections)

sectiontype = f.read_record('|S16')
print("SECTIONTYPE:     ", sectiontype)

sectionname = f.read_record('|S16')
print("SECTIONNAME:     ", sectionname)

npart = np.zeros(ntotgrids, dtype=itype)
totnpart = 0

for i in range(ntotgrids):
    [igrid, npart[i]] = f.read_ints(itype)
    print("IGRID, NPART:    ", igrid, npart[i])
    totnpart = totnpart + npart[i]
print('totnpart:', totnpart)

[numvars] = f.read_ints(itype)
print("NUMVARS:       ", numvars)

for i in range(numvars):
    datadesc = f.read_record('|S16')
    print("DATADESC:        ", datadesc)

# READ TOTAL PARTICLE NUMBERS FOR EACH ST
# change file sequence txt_cont_ / St() / Sk()
dir_name = 'txt_cont_0/St2'
base_name_1 = 'Sk12/pnts_d2_xy'
base_name_2 = 'Sk12/poly_pos_area'
base_name_3 = 'Sk12/poly_pos_area_up'
base_name_4 = 'Sk12/poly_pos_area_normal'

with open(os.path.join(dir_name, 'npdiam' + '.' + 'txt'), 'r') as fd:
    npdiam_tot = fd.readlines()
    for n in npdiam_tot:
        npdiam_tot = ast.literal_eval(n)
    ndiam = npdiam_tot
    # print(ndiam[1])

# VORONOI 2D-PLOT PROJECTED BY ALL PARTICLES
with open(os.path.join(dir_name, base_name_1 + '.' + 'txt'), 'r') as fo1:
    pnts_d1_xy = fo1.readlines()
    pnts_d1_vor2D = []
    for line in pnts_d1_xy:
        pnts_d1_xy = np.array(literal_eval(line), dtype=rtype)
        pnts_d1_vor2D.append(pnts_d1_xy)
# for line in pnts_d1_vor2D:
#     print(line)
# print(len(pnts_d1_vor2D))

points2D_1 = np.array(pnts_d1_vor2D)
# print(points2D_1.shape, points2D_1[0].shape)

# CUT BOUNDARY POINTS WITHIN DELTA
# difference between tmp_2D and vor2D_1/2/3 decreases when delta is smaller
# still the frond boundary and four corners have dashed lines
# delta = 0.9 gives least dashed lines
lxf = -5
lxb = 16.384
lyt = 8.192
lyb = -8.192
delta = 0.95
tmp_2D = []
for p in pnts_d1_vor2D:
    if p[0] < lxf + delta:
        tmp_2D.append([p[0], p[1]])
    if p[0] > lxb - delta:
        tmp_2D.append([p[0], p[1]])
    if p[1] < lyb + delta:
        tmp_2D.append([p[0], p[1]])
    if p[1] > lyt - delta:
        tmp_2D.append([p[0], p[1]])
print(len(pnts_d1_vor2D) - len(tmp_2D))

pnts_d1_vor2D = list(map(tuple, pnts_d1_vor2D))
tmp_2D = list(map(tuple, tmp_2D))
print(len(pnts_d1_vor2D), len(tmp_2D))
vor2D_1 = list(set(pnts_d1_vor2D) - set(tmp_2D))
print(len(vor2D_1))
# np.savetxt('txt_1/St3/vor2D_1.txt', vor2D_1, fmt="%10.7e")

vor2D_pnts_1 = Voronoi(np.array(vor2D_1))
fig1 = voronoi_plot_2d(vor2D_pnts_1, show_points=False, show_vertices=False, line_colors='black', line_width=1, line_alpha=0.4, point_size=2)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()
'''
# COMPUTE VORONOI 2D AREA AND OUTPUT CORRESPONDING POINTS
# NOTE: ConvexHull -> circumscribed sphere volume(v=4pi*rv^3/3) and area(s=4pi*rs^2)(voronoi cell = 1/(rv^2/rs^3) )
# NOTE: 2D voronoi -> polygon.area
start = time.time()
dim = 2
vor = Voronoi(np.array(vor2D_1))
poly_pos_area = []
for i in range(len(vor2D_1)):
    vi = vor.point_region[i]            # corresponding region index to a given input point
    v = vor.regions[i]                  # Voronoi vertices indices for each Voronoi region. -1 indicates infinite vertex
    l = len(v)                          # number of side ridges(vertices)
    if vi >= 0 and l > dim and min(v) >= 0:
        poly_pos_area_i = np.zeros(3)
        polygon = Polygon(vor.vertices[vor.regions[vi]])
        # print(vor.vertices[vor.regions[vi]])
        poly_pos_area_i[0] = vor.points[i, 0]
        poly_pos_area_i[1] = vor.points[i, 1]
        poly_pos_area_i[2] = polygon.area
        poly_pos_area.append(poly_pos_area_i)
# print(len(poly_pos_area))
# np.savetxt(os.path.join(dir_name, base_name_2 + '.' + 'txt'), poly_pos_area, fmt='%10.7e')

# ascending sort by volume value and find out the big vor_vol value between 10^1 - 10^7
idices = np.lexsort(np.array(poly_pos_area).T)
poly_pos_area_up = np.array(poly_pos_area)[idices, :]
# print(np.amax(poly_pos_area_up))
# poly_pos_area_1 = np.array([i for i in poly_pos_area_up if 10e-6 < i[2] < 2.5*10e-5])
# poly_pos_area_2 = np.array([i for i in poly_pos_area_up if 10e-5 < i[2] < 10e-4])
# poly_pos_area_3 = np.array([i for i in poly_pos_area_up if 10e-4 < i[2] < 10e-3])
# poly_pos_area_4 = np.array([i for i in poly_pos_area_up if 10e-3 < i[2] < 10e-2])
# poly_pos_area_5 = np.array([i for i in poly_pos_area_up if 10e-2 < i[2] < 10e-1])
poly_pos_area_ill = np.array([i for i in poly_pos_area_up if i[2] > 1])     # ill points lying on boundaries
# np.set_printoptions(threshold=np.inf)
# print(len(poly_pos_area_5), len(poly_pos_area))
# np.savetxt(os.path.join(dir_name, base_name_3 + '.' + 'txt'), poly_pos_area_up, fmt='%10.7e')
'''
'''
# PLOT VORONOI CELL AND CORRESPONDING POINTS
voronoi_plot_2d(vor, show_points=False, show_vertices=False, point_size=1, line_alpha=0.3)
# colors = np.random.rand(len(vor_sort_big))
# scatter_1 = plt.scatter(poly_pos_area_1[:, 0], poly_pos_area_1[:, 1], s=25, c='red', alpha=1.0)
# scatter_2 = plt.scatter(poly_pos_area_2[:, 0], poly_pos_area_2[:, 1], s=25, c='blue', alpha=1.0)
# scatter_3 = plt.scatter(poly_pos_area_3[:, 0], poly_pos_area_3[:, 1], s=25, c='aqua', alpha=1)
# scatter_4 = plt.scatter(poly_pos_area_4[:, 0], poly_pos_area_4[:, 1], s=40, c='blue', alpha=1, label='$O(10^{-3})\sim O(10^{-2})$')
# scatter_5 = plt.scatter(poly_pos_area_5[:, 0], poly_pos_area_5[:, 1], s=40, c='red', alpha=1, label='$O(10^{-2})\sim O(10^{-1})$')
scatter = plt.scatter(poly_pos_area_ill[:, 0], poly_pos_area_ill[:, 1], s=80, c='lime', alpha=1, label='$>O(1)$')
circle = plt.Circle((0, 0), 0.5, fill=True, color='gray')
# plt.gcf().gca().add_artist(scatter_1)
# plt.gcf().gca().add_artist(scatter_2)
# plt.gcf().gca().add_artist(scatter_3)
# plt.gcf().gca().add_artist(scatter_4)
# plt.gcf().gca().add_artist(scatter_5)
plt.gcf().gca().add_artist(scatter)
plt.gcf().gca().add_artist(circle)
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel('X/D', fontsize=30, fontweight='bold')
plt.ylabel('Y/D', fontsize=30, fontweight='bold')
plt.tick_params(axis='both', labelsize=30)
plt.gca().set_aspect('equal', adjustable='box')
plt.legend(markerscale=2, loc=2, bbox_to_anchor=(0.08, 0.92), fontsize=22)
# figure = plt.gcf()
# figure.set_size_inches(18, 14)
plt.show()
'''
'''
# PLOT VORONOI AREA CONTOUR
poly_pos_area_normal = np.array([i for i in poly_pos_area_up if i[2] <= 1])
# print(len(poly_pos_area_normal))
# np.savetxt(os.path.join(dir_name, base_name_4 + '.' + 'txt'), poly_pos_area_normal, fmt='%10.7e')

# minima = poly_pos_area_normal[0, 2]
# maxima = poly_pos_area_normal[-1, 2]
minima = 10e-5
maxima = 10e-2
# norm = mpl.colors.Normalize(vmin=minima, vmax=maxima, clip=True)   # linear normalization
norm = colors.LogNorm(vmin=minima, vmax=maxima, clip=True)           # log normalization
mapper = plt.cm.ScalarMappable(norm=norm, cmap=cm.RdBu_r)
mapper.set_array([])
voronoi_plot_2d(vor, show_points=False, show_vertices=False, point_size=0.5, line_alpha=0.37)
vor_area_sort = []
polygon_sort = []
for i in range(len(vor.point_region)):
    v = vor.regions[i]
    l = len(v)
    vi = vor.point_region[i]
    region = vor.regions[vi]
    polygon = Polygon(vor.vertices[vor.regions[vi]])
    vor_area = polygon.area
    if vi >= 0 and l > dim and min(v) >= 0 and vor_area < 1:
        vor_area_sort.append(vor_area)
        polygon_sort.append([vor.vertices[i] for i in region])
for j in range(len(vor_area_sort)):
    # print(j, vor_area_sort[j])
    plt.fill(*zip(*polygon_sort[j]), color=mapper.to_rgba(vor_area_sort[j]))

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.xlabel('X/D', fontsize=30, fontweight='bold')
plt.ylabel('Y/D', fontsize=30, fontweight='bold')
plt.tick_params(labelsize=30)
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes('right', size='1.4%', pad=0.12)
cbar = plt.colorbar(mapper, extend='both', cax=cax)
cbar.set_label(label='Voro/$\mathbf{D^2}$', size=28, weight='bold')
cbar.ax.tick_params(labelsize=24)

plt.show()
'''

