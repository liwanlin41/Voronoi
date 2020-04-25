import numpy as np
from mayavi import mlab

def pick_point(picker_obj):
    picked = picker_obj.actors
    if sphere.actor.actor._vtk_obj in [o._vtk_obj for o in picked]:
        r, c = np.lib.index_tricks.unravel_index(picker_obj.point_id, shape)
        x = x_arr[r][c]
        y = y_arr[r][c]
        z = z_arr[r][c]
        print("x = %f, y = %f, z = %f" %(x, y, z))
        mlab.points3d(x,y,z, 0.03, color=(1,0,0),scale_factor=1)


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
    cursor3d = mlab.points3d(0,0,0, mode='axes', color=(0,0,0), scale_factor=0.5)

    
    # setup for events
    points = set()
    verbose = False # set to True for printed output
    select_allowed = True

    fig.on_mouse_pick(pick_point)

    mlab.show()

