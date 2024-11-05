#!/usr/bin/env python3

import time
import numpy as np
from multiprocessing import Process, SimpleQueue, Lock
from mgtools import MGreadH5, ScalarField, VectorField, TensorField, h5Viewer, vtkViewer

# import logging
# logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')


def worker(q, N, l):
    reader = MGreadH5('./mglet_fieldgrid.h5')

    while not q.empty():
        startGrid = q.get()
    
        for grid in range(startGrid, startGrid + N):
            tic = time.perf_counter()
            
            U = reader[grid].read('U')
            V = reader[grid].read('V')
            W = reader[grid].read('W')
          
            velocity = VectorField.velocity(U, V, W)
            J = TensorField.gradient(U, V, W)
            lambda2 = J.lambda2()
            vorticity = J.vorticity()

            l.acquire()
            viewer = h5Viewer('cyl_34_1.h5:/grid-{}'.format(grid))
            viewer.write(velocity, "velocity")
            viewer.write(lambda2, "lambda2")
            viewer.write(vorticity, "vorticity")
            viewer.close()
            l.release()

            toc = time.perf_counter()
            print("Processed grid {0} in {1:.3f} sec.".format(grid, toc-tic))

    reader.close()

def main():
    
    workers = 10

    queue = SimpleQueue()
    lock = Lock()

    # Add the four first grids to the queue
    N = 1
    for i in range(0, 1131, N):
        queue.put(i)

    # Start the workers
    processes = []
    for w in range(workers):
        p = Process(target=worker, args=(queue, N, lock))
        p.start()
        processes.append(p)

    # Wait for workers to finish
    for p in processes:
        p.join()

if __name__ == '__main__':
    main()

