from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Button
import matplotlib.pyplot as plt
import numpy as np
import re

default_elev = 30
default_azim = 300

def disable_vert_rot(event):
    # from https://stackoverflow.com/questions/37457680/how-disable-vertical-camera-rotation-for-3d-plot-in-matplotlib
    azim = ax.azim
    ax.view_init(elev=default_elev, azim = azim)

def onclick(event):
    # from https://stackoverflow.com/questions/6748184/matplotlib-plot-surface-get-the-x-y-z-values-written-in-the-bottom-right-cor/9673338#9673338
    try:
        data_string = ax.format_coord(event.xdata, event.ydata)
    except:
        return
    coord_list_parsed = re.split(r'[xyz=,\s]\s*', data_string) # contains empty strings
    coord_list = np.array([float(s) for s in coord_list_parsed if len(s) > 0])
    mag = np.linalg.norm(coord_list)
    scaled_point = coord_list/mag
#    scaled_point.resize(1,3) # for extraction purposes
    print(data_string)
    ax.scatter(scaled_point[0], scaled_point[1], scaled_point[2], c='r')
#    ax.view_init(elev=default_elev, azim = default_azim) # default values
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

    fig.canvas.mpl_connect('motion_notify_event', disable_vert_rot)
    cid = fig.canvas.mpl_connect('button_release_event', onclick)

    plt.show()
    plt.draw()
