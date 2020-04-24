import matplotlib.pyplot as plt
import numpy as np
import re
from matplotlib.widgets import Button
from point import Point
from plane import Voronoi

def onclick(event): # point clicked in plane
    if select_allowed and event.inaxes == ax:
        data_string = ax.format_coord(event.xdata, event.ydata)
        coord_list_parsed = re.split(r'[xy=,\s]\s*', data_string)
        coord_list = np.array([float(s) for s in coord_list_parsed if len(s) > 0])
        x, y = coord_list
#        points.add(Point(x,y))
        print(data_string)
#        ax.scatter(x, y, c='r') # show points in red
        fig.canvas.draw()

def button_click(event):
    if event.inaxes == clear_ax:
        points.clear()
        ax.clear()
        ax.set_xlim(min_coord, max_coord)
        ax.set_ylim(min_coord, max_coord)
        # re-allow point listening
        select_allowed = True
    if event.inaxes == button_ax:
        for point in points:
            ax.scatter(point.get_x(), point.get_y(), c='r')
        select_allowed = False # stop listening for new points
        voronoi = Voronoi(points)
        while(not voronoi.done()):
            voronoi.step()
        edge_dict = voronoi.output()
        for edge in edge_dict:
            point_list = list(edge_dict[edge]) # hold the points to draw
            if len(point_list) == 1: # extend to infinity
                site1, site2 = edge.get_sites()
                # make sure site1 < site2 for convenience
                if site1 > site2: site1, site2 = site2, site1
                midpoint_x = (site1.get_x() + site2.get_x())/2
                midpoint_y = (site1.get_y() + site2.get_y())/2
                midpoint = Point(midpoint_x, midpoint_y)
                if midpoint == point_list[0]: # rotate site2 90 about midpoint
                    bisector_point_x = midpoint_x - (site2.get_y() - midpoint_y)
                    bisector_point_y = midpoint_y + site2.get_x() - midpoint_x
                    bisector = Point(bisector_point_x, bisector_point_y)
                else: bisector = midpoint
                try:
                    print("extend_ray calling for %s, %s" %(site1, site2))
                    vertex1, vertex2 = extend_ray(point_list[0], bisector)
                    draw_segment(vertex1, vertex2)
                except ValueError:
                    continue
            elif len(point_list) == 2: # crop to box
                try:
                    print("intersect_box calling with %s, %s" %(point_list[0], point_list[1]))
                    crop1, crop2 = intersect_box(point_list[0], point_list[1])
                    draw_segment(crop1, crop2)
                except ValueError:
                    print("something went wrong")
                    continue
        fig.canvas.draw()
                

def draw_segment(p1, p2):
    ''' draw line segment between Points p1, p2 '''
    xs = np.array([p1.get_x(), p2.get_x()])
    ys = np.array([p1.get_y(), p2.get_y()])
    ax.plot(xs, ys)


# helper function
def extend_ray(a, b):
    ''' given two points a, b, return a tuple of points representing the intersection
    of ray ab with the figure box '''
    a_x = a.get_x()
    a_y = a.get_y()
    b_x = b.get_x()
    b_y = b.get_y()
    # deal with bad cases
    if np.isnan(a_x) or np.isnan(b_x):
        raise ValueError
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
    except ZeroDivisionError: # vertical line
        y_dir = max_coord if b_y > a_y else min_coord
        extended = Point(a_x, y_dir)

    cast_a = Point(a_x, a_y) # no infinities
    if inside_box(cast_a):
        return (cast_a, extended)
    return extend_ray(extended, cast_a) # this will crop a to fit the box

def inside_box(p):
    ''' determine if a point p is inside the grid box '''
    return (min_coord <= p.get_x() <= max_coord) and (min_coord <= p.get_y() <= max_coord)

def intersect_box(a, b):
    ''' given two points a, b, return a tuple of points representing the same
    line that fits inside the figure box '''
    if inside_box(a):
        if inside_box(b): return (a, b)
        return extend_ray(a, b)
    return extend_ray(b, a) # crop a to fit the box, also crop b if necessary


if __name__ == '__main__':
    fig = plt.figure(figsize=(16,10))
    ax = fig.add_subplot(111)
    min_coord = -10
    max_coord = 10
    ax.set_xlim(min_coord,max_coord)
    ax.set_ylim(min_coord,max_coord)

    points = {Point(0,2), Point(0,0), Point(0,-2), Point(2,0), Point(-2,0)}
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
