from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Button
import matplotlib.pyplot as plt
import numpy as np
import re
from point3d import Point3D
from sphere import VoronoiSphere

default_elev = 30
default_azim = 300
scale = 0.92

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
        ax.scatter(x, y, z, c='r', zorder = 5, linewidth = 5.0)
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
        edge_dict_far, edge_dict_near = voronoi.output()
        for edge in edge_dict_near:
            midpoint = get_midpoint(edge)
            point_list = list(edge_dict_near[edge])
            if len(point_list) == 1:
                print("defective")
            elif len(point_list) == 2:
                draw_arc(point_list[0], point_list[1], midpoint)
            else:
                print("this is a weird number")
                for i in range(len(point_list) - 1):
                    draw_arc(point_list[i], point_list[i+1], midpoint)
#        input()
#        for edge in edge_dict_far:
#            point_list = list(edge_dict_far[edge])
#            if len(point_list) == 1:
#                print("defective")
#            elif len(point_list) == 2:
#                midpoint = get_midpoint(edge)
#                draw_arc(point_list[0], point_list[1], midpoint)
#            elif len(point_list) > 2:
#                print("why are there so many")
##                for i in range(len(point_list) - 1):
##                    draw_arc(point_list[i], point_list[i+1], edge.get_sites())
        fig.canvas.draw()

def get_midpoint(edge):
    ''' find the midpoint of the 3D site points of the edge '''
    site1, site2 = edge.get_sites()
    print("edge between %s, %s" %(site1, site2))
    midpoint_arr = np.array([(site1.x + site2.x)/2, (site1.y + site2.y)/2, (site1.z + site2.z)/2])
    midpoint = midpoint_arr / np.linalg.norm(midpoint_arr)
    return midpoint


def draw_segment(p1, p2):
    ''' draw line segment between points p1, p2 '''
    xs = np.array([p1.x, p2.x])
    ys = np.array([p1.y, p2.y])
    zs = np.array([p1.z, p2.z])
    ax.plot(xs, ys, zs)
    fig.canvas.draw()

def draw_arc(p1, p2, midpoint):
    ''' draw the great circle arc from p1 to p2 containing midpoint '''
    u = np.array([p1.x, p1.y, p1.z])
    v = np.array([p2.x, p2.y, p2.z])
    uv = np.cross(u,v) # get a direction vector
    w = np.cross(uv, u) # u, w parametrize the great circle
    w = w / np.linalg.norm(w) # scale to unit vector
    # w and v should be in the same direction relative to u
    theta_max = np.arccos(np.dot(u,v)) # angle between u and v
    # determine if midpoint is on arc from u to v or from v to u
    mid_dir = np.cross(u, midpoint)
    if np.dot(mid_dir, uv) >= 0: # same direction
        print("same direction")
        theta = np.mgrid[0:theta_max:num_sample]
    else:
        print("opposite direction")
        theta = np.mgrid[theta_max:2*np.pi:num_sample]
    x = np.cos(theta) * u[0] + np.sin(theta) * w[0]
    y = np.cos(theta) * u[1] + np.sin(theta) * w[1]
    z = np.cos(theta) * u[2] + np.sin(theta) * w[2]
    ax.plot(x,y,z,linewidth=7.0, zorder = 10)
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

    # draw a solid internal sphere to aid visualization
#    ax.plot_surface(scale*x, scale*y, scale*z, color='b', zorder = 0)

    ax.plot_surface(x,y,z, alpha=.5, zorder = 1)
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

#    ax.plot_surface(scale*x, scale*y, scale*z, color='w', zorder = 0)
    ax.plot_surface(x,y,z, alpha=.5, zorder = 1)
    
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
