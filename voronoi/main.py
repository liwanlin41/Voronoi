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
        if (x, y, z) == (0, 0, 1):
            print("Please don't pick the north pole")
            return
#        scaled_point.resize(1,3) # for extraction purposes
        print("x = %f, y = %f, z = %f" %(x, y, z))
#        ax.view_init(elev=default_elev, azim = default_azim) # default values
        new_point = Point3D(x,y,z)
        points.add(new_point)
        ax.scatter(x, y, z, c='r', linewidth = 3.0)
        fig.canvas.draw()

def button_click(event):
    global select_allowed # toggle point selection enabling
    if event.inaxes == clear_ax:
        points.clear()
        ax.clear()
        ax.set_zlim(-1,1)
        ax.set_aspect('equal')
        # re-allow point listening
        select_allowed = True
        draw_sphere()
        print() # separate for new point inputs
    elif event.inaxes == start_ax:
        select_allowed = False
        voronoi = VoronoiSphere(points, verbose)
        while(not voronoi.done()):
            voronoi.step()
#        eta = voronoi.eta
#        ax.scatter(eta.x, eta.y, eta.z, color = 'g', linewidth = 10)
        edge_dict_far, edge_dict_near = voronoi.output()
        for edge in edge_dict_near:
            vertex_set = edge_dict_near[edge]
            point_list = list(vertex_set)
            if len(point_list) == 1:
                print("defective")
            elif len(point_list) == 2:
#                site1, site2 = edge.get_sites()
#                draw_segment(site1, site2)
                draw_arc(point_list[0], point_list[1])
            else:
                print("this is a weird number")
                for i in range(len(point_list) - 1):
                    draw_arc(point_list[i], point_list[i+1])
#        input() # pause between parts
        for edge in edge_dict_far:
#            print(edge)
            point_list = list(edge_dict_far[edge])
            if len(point_list) == 1:
                print("defective")
#                print(edge)
#                print(point_list[0])
            elif len(point_list) == 2:
                draw_arc(point_list[0], point_list[1])
            else:
                print("why are there so many")
                for i in range(len(point_list) - 1):
                    draw_arc(point_list[i], point_list[i+1])
        fig.canvas.draw()

def draw_segment(p1, p2):
    ''' draw line segment between points p1, p2 '''
    xs = np.array([p1.x, p2.x])
    ys = np.array([p1.y, p2.y])
    zs = np.array([p1.z, p2.z])
    ax.plot(xs, ys, zs) 
    fig.canvas.draw()

def draw_arc(p1, p2):
    ''' draw the minor arc of the great circle from p1 to p2 '''
    u = np.array([p1.x, p1.y, p1.z])
    v = np.array([p2.x, p2.y, p2.z])
    uv = np.cross(u,v) # get a direction vector
    w = np.cross(uv, u) # u, w parametrize the great circle 
    w = w / np.linalg.norm(w) # scale to unit vector
    # w and v should be in the same direction relative to u
    theta_max = np.arccos(np.dot(u,v)) # angle between u and v
    theta = np.mgrid[0:theta_max:num_sample]
    x = np.cos(theta) * u[0] + np.sin(theta) * w[0]
    y = np.cos(theta) * u[1] + np.sin(theta) * w[1]
    z = np.cos(theta) * u[2] + np.sin(theta) * w[2]
    ax.plot(x,y,z,linewidth=7.0)
    fig.canvas.draw()
#    input() # to draw one at a time

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


### TESTING ###
#    point1 = Point3D(0,0,1)
#    point2 = Point3D(0,1,0)
#    midpoint = np.array([0,0.6,0.8])
#    draw_arc(point2, point1, midpoint, False)
#    point1 = Point3D(0,1,0)
#    point2 = Point3D(1,0,0)
#    point3 = Point3D(0,0,-1)
#    point4 = Point3D(-0.5, 0.5, 0.5**0.5)
#    point1 = Point3D(1,0,0)
#    point2 = Point3D(0,0.6,-0.8)
#    point3 = Point3D(-0.5,0.5, 0.5**0.5)
#    point4 = Point3D(0.36, -0.48, 0.8)
#    ax.scatter(point1.x, point1.y, point1.z) 
#    ax.scatter(point2.x, point2.y, point2.z)
#    ax.scatter(point3.x, point3.y, point3.z)
#    ax.scatter(point4.x, point4.y, point4.z)
#    points = {point1, point2, point3, point4}
#    print(point1)
#    print(point2)
#    print(point3)
#    print(point4)

    plt.show()
    plt.draw()
