import numpy as np
from scipy.io import FortranFile

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
ngrid = np.zeros(ntotgrids, dtype=itype)
totnpart = 0
for i in range(ntotgrids):
    [igrid, npart[i]] = f.read_ints(itype)
    print("IGRID, NPART:    ", igrid, npart[i])
    ngrid[i] = igrid
    totnpart = totnpart + npart[i]
print('totnpart:', totnpart)

[numvars] = f.read_ints(itype)
print("NUMVARS:       ", numvars)

for i in range(numvars):
    datadesc = f.read_record('|S16')
    print("DATADESC:        ", datadesc)

ppos_tmp = []
pvel_tmp = []
pvelf_tmp = []
pdiam_tmp = []
pidx_tmp = []
# grid_rm = list(np.arange(1, 9)) + list(np.arange(29, 61))
# grids = [x for x in list(np.arange(1, ntotgrids+1)) if x not in grid_rm]
for i in range(ntotgrids):
    ppos = np.zeros([npart[i], 3])
    pvel = np.zeros([npart[i], 3])
    pvelf = np.zeros([npart[i], 3])
    pidx = np.zeros(npart[i], dtype=np.int32)
    pdiam = np.zeros(npart[i])
    for j in range(npart[i]):
        data = f.read_record(part_type)
        if data['DIAMCLASS'] == 1:
            ppos[j, 0] = data['XPART_X']
            ppos[j, 1] = data['XPART_Y']
            ppos[j, 2] = data['XPART_Z']
            pvel[j, 0] = data['UPART_X']
            pvel[j, 1] = data['UPART_Y']
            pvel[j, 2] = data['UPART_Z']
            pvelf[j, 0] = data['OLD_UPART_X']
            pvelf[j, 1] = data['OLD_UPART_Y']
            pvelf[j, 2] = data['OLD_UPART_Z']
            pidx[j] = data['PART_IDX']
            pdiam[j] = data['DIAMCLASS']

            ppos_tmp.append(ppos[j, :])
            pvel_tmp.append(pvel[j, :])
            pvelf_tmp.append(pvelf[j, :])
            pdiam_tmp.append(pdiam[j])
            pidx_tmp.append(pidx[j])
            # print(i + 1, j, npart[i], ppos[j, 0])
    if (npart[i] < 1):
        continue

ppos_tot = []
pvel_tot = []
pvelf_tot = []
pdiam_tot = []
pidx_tot = []
ind = [i for i, elem in enumerate(ppos_tmp) if -2 < elem[0] < 20]
for i in ind:
    ppos_tot.append(ppos_tmp[i])
    pvel_tot.append(pvel_tmp[i])
    pvelf_tot.append(pvelf_tmp[i])
    pdiam_tot.append(pdiam_tmp[i])
    pidx_tot.append(pidx_tmp[i])
totnpart = len(ppos_tot)
ppos_tot = np.array(ppos_tot)
pvel_tot = np.array(pvel_tot)
pvelf_tot = np.array(pvelf_tot)
print('useful pp:', ppos_tot.shape)

fo = open('ppgrid-Sk1.vtk', 'wb')
fo.write(b'# vtk DataFile Version 3.0\n')
fo.write(b'MGLET particles\n')
fo.write(b'ASCII\n')
fo.write(b'\n')
fo.write(b'DATASET POLYDATA\n')

fo.write('POINTS {} float\n'.format(totnpart).encode())
for i in range(len(ppos_tot)):
    fo.write('{} {} {}\n'.format(ppos_tot[i, 0], ppos_tot[i, 1], ppos_tot[i, 2]).encode())

fo.write('VERTICES 1 {}\n'.format(totnpart+1).encode())
fo.write('{}\n'.format(totnpart).encode())
for i in range(totnpart):
    fo.write('{}\n'.format(i).encode())

fo.write('POINT_DATA {}\n'.format(totnpart).encode())
fo.write(b'SCALARS idx int 1\n')
fo.write(b'LOOKUP_TABLE default\n')
for i in range(len(pidx_tot)):
    fo.write('{}\n'.format(pidx_tot[i]).encode())

fo.write(b'VECTORS velocity float\n')
for i in range(len(pvel_tot)):
    fo.write('{} {} {}\n'.format(pvelf_tot[i, 0], pvelf_tot[i, 1], pvelf_tot[i, 2]).encode())

fo.write(b'SCALARS diameter float 1\n')
fo.write(b'LOOKUP_TABLE default\n')
for i in range(len(pdiam_tot)):
    fo.write('{}\n'.format(pdiam_tot[i]).encode())

fo.close()
