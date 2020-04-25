from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Button
import matplotlib.pyplot as plt
import numpy as np
import re
from point3d import Point3D
from sphere import VoronoiSphere

default_elev = 30
default_azim = 300

def disable_vert_rot(event):
    # from https://stackoverflow.com/questions/37457680/how-disable-vertical-camera-rotation-for-3d-plot-in-matplotlib
    azim = ax.azim
    ax.view_init(elev=default_elev, azim = azim)

def onclick(event):
    # from https://stackoverflow.com/questions/6748184/matplotlib-plot-surface-get-the-x-y-z-values-written-in-the-bottom-right-cor/9673338#9673338
    if select_allowed and event.inaxes == ax and event.button == 1:
        try:
            data_string = ax.format_coord(event.xdata, event.ydata)
        except:
            return
        coord_list_parsed = re.split(r'[xyz=,\s]\s*', data_string) # contains empty strings
        coord_list = np.array([float(s) for s in coord_list_parsed if len(s) > 0])
        coord_list[0] = round(coord_list[0] * 10) / 10
        coord_list[1] = round(coord_list[1] * 10) / 10
        mag = np.linalg.norm(coord_list)
        scaled_point = coord_list/mag
        x, y, z = scaled_point[0], scaled_point[1], scaled_point[2]
#        scaled_point.resize(1,3) # for extraction purposes
        print("x = %f, y = %f, z = %f" %(x, y, z))
        ax.scatter(x, y, z, c='r')
#        ax.view_init(elev=default_elev, azim = default_azim) # default values
        points.add(Point3D(x, y, z))
        fig.canvas.draw()

def button_click(event):
    if event.inaxes == clear_ax:
        points.clear()
        ax.clear()
        ax.set_zlim(-1,1)
        ax.set_aspect('equal')
        # re-allow point listening; this is currently broken
        select_allowed = True
        draw_sphere()
        print() # separate for new point inputs
    elif event.inaxes == start_ax:
        select_allowed = False
        voronoi = VoronoiSphere(points, verbose)
        while(not voronoi.done()):
            voronoi.step()
        edge_dict = voronoi.output()
        for edge in edge_dict:
            point_list = list(edge_dict[edge])
            if len(point_list) == 1:
                print("defective")
            elif len(point_list) > 1:
                for i in range(len(point_list) - 1):
                    draw_segment(point_list[i], point_list[i+1])
        fig.canvas.draw()

def draw_segment(p1, p2):
    ''' draw line segment between points p1, p2 '''
    xs = np.array([p1.x, p2.x])
    ys = np.array([p1.y, p2.y])
    zs = np.array([p1.z, p2.z])
    ax.plot(xs, ys, zs)
    fig.canvas.draw()

def draw_sphere():
    # sphere coordinates
    num_sample = 20j
    # convert to 2D arrays for ease of multiplication
    theta, phi = np.mgrid[0:2*np.pi:2*num_sample, -0.5*np.pi:0.5*np.pi:num_sample]

    # convert to spherical coordinates, get num_sample^2 points as 2D array
    x = np.cos(theta) * np.cos(phi)
    y = np.sin(theta) * np.cos(phi)
    z = np.sin(phi)

    ax.plot_surface(x,y,z, alpha=.3)
    fig.canvas.draw()


if __name__ == '__main__':
    fig = plt.figure(figsize=(16,10))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect("equal")
    ax.set_zlim(-1, 1)

    # sphere coordinates
    num_sample = 20j
    # convert to 2D arrays for ease of multiplication
    theta, phi = np.mgrid[0:2*np.pi:2*num_sample, -0.5*np.pi:0.5*np.pi:num_sample]

    # convert to spherical coordinates, get num_sample^2 points as 2D array
    x = np.cos(theta) * np.cos(phi)
    y = np.sin(theta) * np.cos(phi)
    z = np.sin(phi)

    ax.plot_surface(x,y,z, alpha=.3)
    
    # setup for events
    points = set()
    verbose = False # set to True for printed output
    select_allowed = True

    # create button
    start_ax = plt.axes([0.4, 0, 0.1, 0.075])
    start = Button(start_ax, "START")
    start.on_clicked(button_click)

    clear_ax = plt.axes([0.5, 0, 0.1, 0.075])
    clear_button = Button(clear_ax, "CLEAR")
    clear_button.on_clicked(button_click)

    fig.canvas.mpl_connect('motion_notify_event', disable_vert_rot)
    cid = fig.canvas.mpl_connect('button_release_event', onclick)

    plt.show()
    plt.draw()
