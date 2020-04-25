import numpy as np
from mayavi import mlab
from point3d import Point3D
from sphere import VoronoiSphere

def pick_point(picker_obj):
    picked = picker_obj.actors
    if select_allowed and (sphere.actor.actor._vtk_obj in [o._vtk_obj for o in picked]):
        r, c = np.lib.index_tricks.unravel_index(picker_obj.point_id, shape)
        x = x_arr[r][c]
        y = y_arr[r][c]
        z = z_arr[r][c]
        print("x = %f, y = %f, z = %f" %(x, y, z))
        mlab.points3d(x,y,z, 0.03, color=(1,0,0),scale_factor=1)
        points.add(Point3D(x,y,z))

def draw_arc(p1,p2, midpoint_arr):
    ''' draw the great circle arc from p1 to p2 containing midpoint,
    represented as an np array'''
    u = np.array([p1.x, p1.y, p1.z])
    v = np.array([p2.x, p2.y, p2.z])
    uv = np.cross(u,v) # get a direction vector
    w = np.cross(uv, u) # u, w parametrize the great circle
    w = w / np.linalg.norm(w)
    theta_v = np.arccos(np.dot(u,v)) # angle between u and v
    # determine if midpoint is on arc from u to v or from v to u
    mid_dir = np.cross(u, midpoint_arr)
    if np.dot(mid_dir, uv) >= 0: # same direction
        theta = np.mgrid[0:theta_v:num_sample]
    else:
        theta = np.mgrid[theta_v:2*np.pi:num_sample]
    x = np.cos(theta) * u[0] + np.sin(theta) * w[0]
    y = np.cos(theta) * u[1] + np.sin(theta) * w[1]
    z = np.cos(theta) * u[2] + np.sin(theta) * w[2]
    mlab.plot3d(x,y,z)


if __name__ == '__main__':
    fig = mlab.figure(1)
    mlab.title("Voronoi Sphere")
    # sphere coordinates
    num_sample = 180j
    # convert to 2D arrays for ease of multiplication
    theta, phi = np.mgrid[0:2*np.pi:2*num_sample, -0.5*np.pi:0.5*np.pi:num_sample]

    # convert to spherical coordinates, get num_sample^2 points as 2D array
    x_arr = np.cos(theta) * np.cos(phi)
    y_arr = np.sin(theta) * np.cos(phi)
    z_arr = np.sin(phi)
    shape = x_arr.shape
    # draw sphere
    sphere = mlab.mesh(x_arr, y_arr, z_arr, color=(0,0,1))
    
    # setup for events
    points = set()
    verbose = False # set to True for printed output
    select_allowed = True

    fig.on_mouse_pick(pick_point)
    point1 = Point3D(0,0,1)
    point2 = Point3D(1,0,0)
    midpoint = np.array([0.6,0,0.8])
    draw_arc(point1, point2, -1*midpoint)

    mlab.show()

