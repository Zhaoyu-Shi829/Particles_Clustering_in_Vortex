import numpy as np
from scipy.io import FortranFile

nlogs = 5000

for i in range(nlogs):
    filename = 'PARTLOG/part-0016{:06d}.dat'.format(i+1)
    print(i, filename)
    log = np.loadtxt(filename, delimiter='\n', dtype=str)
    log_arr = np.reshape(log, (int(len(log)/4), 4))

    ittot, diam = [], []
    x, y, z = [], [], []
    u, v, w = [], [], []
    for j in range(len(log_arr)):
        log_1 = log_arr[:, 0][j].split()
        ittot.append(log_1[0])
        diam.append(log_1[1])
        x.append(log_1[2])
        y.append(log_1[3])
        z.append(log_1[4])

        log_2 = log_arr[:, 1][j].split()
        u.append(log_2[0])
        v.append(log_2[1])
        w.append(log_2[2])

    ntimes = len(ittot)
    print(ntimes)

    fo = open('PARTLOGVTK/part-0016{}.vtk'.format(i+1), 'wb')
    fo.write(b'# vtk DataFile Version 3.0\n')
    fo.write(b'MGLET particles\n')
    fo.write(b'ASCII\n')
    fo.write(b'\n')
    fo.write(b'DATASET POLYDATA\n')

    fo.write('POINTS {} float\n'.format(ntimes).encode())
    for j in range(ntimes):
        fo.write('{} {} {}\n'.format(x[j], y[j], z[j]).encode())
    
    fo.write('VERTICES 1 {}\n'.format(ntimes+1).encode())
    fo.write('{}\n'.format(ntimes).encode())
    for j in range(ntimes):
        fo.write('{}\n'.format(j).encode())
    
    fo.write('LINES 1 {}\n'.format(3*(ntimes-1)).encode())
    for j in range(ntimes-1):
        fo.write('2 {} {}\n'.format(j, j+1).encode())

    fo.write('POINT_DATA {}\n'.format(ntimes).encode())    
    fo.write(b'VECTORS velocity float\n')
    for j in range(ntimes):
        fo.write('{} {} {}\n'.format(u[j], v[j], w[j]).encode())

    fo.write(b'SCALARS diameter float 1\n')
    fo.write(b'LOOKUP_TABLE default\n')
    for j in range(ntimes):
        fo.write('{}\n'.format(diamclass[j]).encode())

    fo.close()


