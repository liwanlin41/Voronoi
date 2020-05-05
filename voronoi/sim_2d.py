# run 2D simulations for large numbers of points, no drawing included

import time
import random
from point import Point
from plane import Voronoi

points = set()

point_nums = [3 * 2**i for i in range(10)]

rand_min = -100
rand_max = 100

for n in point_nums:
    # genreate random points
    for i in range(n):
        raw_x = random.uniform(rand_min, rand_max)
        raw_y = random.uniform(rand_min, rand_max)
        # keep only one decimal precision
        x = round(raw_x * 10) / 10
        y = round(raw_y * 10) / 10
        points.add(Point(x,y))

    assert len(points) == n

    # track time
    start_time = time.process_time()
    voronoi = Voronoi(points, False) # no printed output
    while not voronoi.done():
        voronoi.step()
    output_tup = voronoi.output()

    # computation complete
    elapsed_time = time.process_time() - start_time

    print("Took %f seconds to compute on %d points" %(elapsed_time, n))
    
    # setup for next iteration
    points.clear()

