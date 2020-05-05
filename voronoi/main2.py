import matplotlib.pyplot as plt
import numpy as np
import re
from matplotlib.widgets import Button
from point import Point
from plane import Voronoi

import time

def onclick(event): # point clicked in plane
    if select_allowed and event.inaxes == ax:
        data_string = ax.format_coord(event.xdata, event.ydata)
        coord_list_parsed = re.split(r'[xy=,\s]\s*', data_string)
        coord_list = np.array([float(s) for s in coord_list_parsed if len(s) > 0])
        # keep only one decimal precision
        x = round(coord_list[0] * 10) / 10
        y = round(coord_list[1] * 10) / 10
        points.add(Point(x,y))
        print("x = %f, y = %f" %(x, y))
        ax.scatter(x, y, c='r') # show points in red
        fig.canvas.draw()

def button_click(event):
    if event.inaxes == clear_ax:
        points.clear()
        ax.clear()
        ax.set_xlim(min_coord, max_coord)
        ax.set_ylim(min_coord, max_coord)
        ax.set_aspect('equal')
        # re-allow point listening
        select_allowed = True
        print() # separate for new point inputs
    elif event.inaxes == button_ax:
        select_allowed = False # stop listening for new points

        start_time = time.process_time()
        voronoi = Voronoi(points, verbose)
        while(not voronoi.done()):
            voronoi.step()
#            retry = str(input())
        edge_dict = voronoi.output()[0]

        # computation is complete
        elapsed_time = time.process_time() - start_time

        for edge in edge_dict:
            point_list = list(edge_dict[edge]) # hold the points to draw
            if len(point_list) == 1: # extend to infinity
                print("defective")
                # these edges should go to +infinity
                print(edge)
                site1, site2 = edge.get_sites()
                # make sure site1 < site2 for convenience
                if site1 > site2: site1, site2 = site2, site1
                midpoint_x = (site1.get_x() + site2.get_x())/2
                midpoint_y = (site1.get_y() + site2.get_y())/2
                midpoint = Point(midpoint_x, midpoint_y)
                # rotate site2 90 about midpoint
                rotated_x = midpoint_x - (site2.get_y() - midpoint_y)
                rotated_y = midpoint_y + site2.get_x() - midpoint_x
                rotated = Point(rotated_x, rotated_y)
                existing_point = point_list[0]
                x_diff = rotated_x - midpoint_x
                y_diff = rotated_y - midpoint_y
                if rotated == existing_point:
                    # get a point along the segment
                    extended = Point(2 * rotated_x - midpoint_x, 2 * rotated_y - midpoint_y)
                else: 
                    extended = Point(existing_point.x + x_diff, existing_point.y + y_diff)
                try:
                    vertex1, vertex2 = extend_ray(existing_point, extended)
                    draw_segment(vertex1, vertex2)
                except ValueError:
                    continue
            elif len(point_list) > 1: # crop to box
                for i in range(len(point_list)-1):
                    try:
                        draw_segment(point_list[i], point_list[i+1])
                    except ValueError:
                        print("shouldn't happen")
                        continue
        fig.canvas.draw()

        print("Took %f seconds to compute on %d points" %(elapsed_time, len(points)))
                

def draw_segment(p1, p2):
    ''' draw line segment between Points p1, p2 '''
    xs = np.array([p1.get_x(), p2.get_x()])
    ys = np.array([p1.get_y(), p2.get_y()])
    ax.plot(xs, ys)
    fig.canvas.draw()
#    wait = str(input())


# helper functions
#def has_box_intersection(p1, p2):
#    ''' given two points p1, p2, determine if p1 p2 intersects the bounding box '''
#    # determine if intersection exists by computing sign of bounding box coords
#    # with respect to line
#    x1 = p1.get_x()
#    y1 = p1.get_y()
#    x2 = p2.get_x()
#    y2 = p2.get_y()
#    some_sign = None
#    for x in [min_coord, max_coord]:
#        for y in [min_coord, max_coord]:
#            # get equation in the form ax + by + c = 0
#            a = y2 - y1
#            b = x1 - x2
#            c = -a * x1 - b * y1
#            sign = a * x + b * y + c
#            if some_sign is None:
#                some_sign = sign
#            if some_sign != sign:
#                return True
#    return False
#
#def is_between(p, p1, p2): # return true if p lies strictly between p1 and p2
#    return p1 <= p <= p2 or p2 <= p <= p1
#
#
def extend_ray(a, b):
    ''' given two points a, b, return a tuple of points representing the intersection
    of ray ab with the figure box 
    raise ValueError if there is no such intersection or if input values are bad'''
    a_x = a.get_x()
    a_y = a.get_y()
    b_x = b.get_x()
    b_y = b.get_y()
    # deal with bad cases
    if np.isnan(a_x) or np.isnan(b_x):
#or not has_box_intersection(a, b):
        print("not a number")
        raise ValueError

    # intersection exists
    try:
        slope = (b_y-a_y)/(b_x-a_x)
        intercept = a_y - slope * a_x
        horizontal_lim = max_coord if b_x > a_x else min_coord
        vertical_lim = max_coord if b_y > a_y else min_coord
        if slope == 0:
            extended = Point(horizontal_lim, a_y)
        else: # otherwise this line is not axis-parallel
            horizontal_intersect_y = slope * horizontal_lim + intercept
            if min_coord <= horizontal_intersect_y <= max_coord:
                extended = Point(horizontal_lim, horizontal_intersect_y)
            else:
                vertical_intersect_x = (vertical_lim - intercept) / slope
                extended = Point(vertical_intersect_x, vertical_lim)
    except (RuntimeWarning, ZeroDivisionError) as e: # vertical line
#        print("divide by zero")
        y_dir = max_coord if b_y > a_y else min_coord
        extended = Point(a_x, y_dir)

#    cast_a = Point(a_x, a_y) # no infinities
    if not inside_box(extended): # doesn't actually intersesct
#        print("no intersection")
        raise ValueError
    if inside_box(a):
        return (a, extended)
    return extend_ray(extended, a) # this will crop a to fit the box

def inside_box(p):
    ''' determine if a point p is inside the grid box '''
    return (min_coord <= p.get_x() <= max_coord) and (min_coord <= p.get_y() <= max_coord)

#def intersect_box(a, b):
#    ''' given two points a, b, return a tuple of points representing the same
#    line that fits inside the figure box 
#    raise ValueError if no intersection'''
#    print("intersect box with %s, %s" %(a, b))
#    if inside_box(a):
#        if inside_box(b): return (a, b)
#        bounded_points = extend_ray(a,b)
##        return extend_ray(a, b)
#    else:
#        bounded_points = extend_ray(b, a) # crop a to fit box, also crop b if necessary
#    # check to see if the order of these points along the line is correct
#    print("Computed points: %s, %s" %(bounded_points[0], bounded_points[1]))
#    if not is_between(bounded_points[0], a, b) or not is_between(bounded_points[1], a, b):
#        print("no intersection here")
#        raise ValueError
#    return bounded_points


if __name__ == '__main__':
    fig = plt.figure(figsize=(16,10))
    ax = fig.add_subplot(111)
    min_coord = -10
    max_coord = 10
    ax.set_xlim(min_coord,max_coord)
    ax.set_ylim(min_coord,max_coord)
    ax.set_aspect('equal')

    points = set()

    ### manual testing ###

#    points.add(Point(x = 6.500000, y = 5.500000))
#    points.add(Point(x = 1.100000, y = 3.500000))
#    points.add(Point(x = 3.800000, y = 4.500000))
#    points.add(Point(x = 0.000000, y = 1.300000))
#
#    for point in points:
#        ax.scatter(point.get_x(), point.get_y(), c='r')

    verbose = False # set this for printed output
    select_allowed = True

    # create button
    button_ax = plt.axes([0.4, 0, 0.1, 0.075])
    button = Button(button_ax, "START")
    button.on_clicked(button_click)

    clear_ax = plt.axes([0.5, 0, 0.1, 0.075])
    clear_button = Button(clear_ax, "CLEAR")
    clear_button.on_clicked(button_click)

    cid = fig.canvas.mpl_connect('button_release_event', onclick)

    plt.show()
    plt.draw()
