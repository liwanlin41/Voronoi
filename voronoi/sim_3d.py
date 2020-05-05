# run 3D simulations for large numbers of points without any drawing included

import time
import random
from point3d import Point3D
from sphere import VoronoiSphere

points = set()

point_nums = [3 * 2**i for i in range(10)]

for n in point_nums:
    # generate random points on unit sphere
    for i in range(n):
        point_selected = False
        while not point_selected:
            raw_x = random.uniform(-1, 1)
            raw_y = random.uniform(-1, 1)
            raw_z = random.uniform(-1, 1)
            magnitude = raw_x **2 + raw_y **2 + raw_z **2
            # scale to sphere
            x = raw_x / magnitude
            y = raw_y / magnitude
            z = raw_z / magnitude
            if (x, y, z) == (0, 0, 1): # avoid north pole
                continue
            point_selected = True
        points.add(Point3D(x, y, z))

    assert len(points) == n

    # track time
    start_time = time.process_time()
    voronoi = VoronoiSphere(points, False) # no printing
    while not voronoi.done():
        voronoi.step()
    output_info = voronoi.output() # extra processing occurs during output

    # computation complete
    elapsed_time = time.process_time() - start_time

    print("Took %f seconds to compute on %d points" %(elapsed_time, n))
    
    # setup for next iteration
    points.clear()
