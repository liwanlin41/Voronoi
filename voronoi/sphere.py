from plane import Voronoi
from point import Point
from point3d import Point3D
from event import are_collinear
import numpy as np

class VoronoiSphere:
    ''' represent a voronoi diagram on a sphere '''

    def __init__(self, point_set, verbose = True):
        ''' create a new VoronoiSphere object from given site points
        point_set is a set of Point3D objects representing the sites 
        verbose indicates printing verbosity'''
        self.verbose = verbose
        self.points = point_set
        # maps for inversion from sphere to plane and vice versa
        self.sphere_to_plane = {}
        self.plane_to_sphere = {}
        for point in self.points:
            inverse = point.invert()
            self.sphere_to_plane[point] = inverse
            self.plane_to_sphere[inverse] = point
        # also want to pick a second center of inversion inside the convex hull
        planar_iter = iter(self.plane_to_sphere)
        triangle = [next(planar_iter)] * 3
        while(are_collinear(triangle[0], triangle[1], triangle[2])):
            triangle[2] = next(planar_iter)
        # triangle will now be a valid non-degenerate triangle
        weight = 0
        while(True):
            eta_inv_x = 0.5 * triangle[0].get_x() + (0.25-weight) * triangle[1].get_x() + (0.25+weight) * triangle[2].get_x()
            eta_inv_y = 0.5 * triangle[0].get_y() + (0.25-weight) * triangle[1].get_y() + (0.25+weight) * triangle[2].get_y()
            eta_inv = Point(eta_inv_x, eta_inv_y)
            if eta_inv in self.plane_to_sphere or eta_inv == Point(0,0): # this point exists
                weight = 1/8 if weight == 0 else weight/2
            else:
                break # good eta found
       a, b, c = eta_inv.project_to_sphere()
       self.eta = Point3D(a, b, c) # second center of inversion

       # maps for inversion to second plane under coordinate transform,
       # where inverting is harder
       self.sphere_to_plane2 = {}
       self.plane2_to_sphere = {}
       for point in self.points:
           inverted = np.array(point.invert_through(self.eta)).reshape((3,1))
           # if x' are new coordinates, matrix_transform @ x' are old coords
           # new coordinates have z = -1 always
           matrix_transform = np.array([[b, -a*c, a],[-a,-b*c,b],[0,a**2+b**2,c]])
           new_coords = np.linalg.inv(matrix_transform) @ inverted
           new_coords.reshape((3,))
           inverted_point = Point(new_coords[0], new_coords[1])
           self.sphere_to_plane2[point] = inverted_point
           self.plane2_to_sphere[inverted_point] = point

        # construct the voronoi diagrams
        self.voronoi1 = Voronoi(set(self.plane_to_sphere.keys()), verbose)
        self.voronoi2 = Voronoi(set(self.plane2_to_sphere.keys()), verbose)

    def step(self):
        ''' handle the next events concurrently '''
        if not self.voronoi1.done():
            self.voronoi1.step()
        if not self.voronoi2.done():
            self.voronoi2.step()

    def output(self):
        ''' combine the two Voronoi diagrams to get the output '''
        raise NotImplementedError

    def done(self):
        return self.voronoi1.done() and self.voronoi2.done()

### helper functions ###

def compute_center(p1, p2, p3, invert = Point3D(0,0,1)):
    ''' given three Point3Ds, return Point3D on sphere that is equidistant
    from all of them, on the side not containing invert '''
    # compute normal to the plane using cross product
    vec1 = np.array([p2.get_x() - p1.get_x(), p2.get_y() - p1.get_y(), p2.get_z() - p1.get_z()])
    vec2 = np.array([p3.get_x() - p1.get_x(), p3.get_y() - p1.get_y(), p3.get_z() - p1.get_z()])
    normal = np.cross(vec1, vec2).reshape((1,3))
    print("NORMAL VECTOR")
    print(normal)
    d = normal @ np.array([[p1.get_x(), p1.get_y(), p1.get_z()]]).T
    # d represents ax + by + cz = d
    invert_point = np.array([[invert.get_x(), invert.get_y(), invert.get_z()]])
    mag = np.linalg.norm(normal)
    if np.sign(normal @ invert_point.T - d) == np.sign(-d): 
        # invert on the same side as the origin
        coords = normal/mag
    else:
        coords = -normal/mag
    return Point3D(coords[0][0], coords[0][1], coords[0][2])

if __name__ == '__main__':
    point_set = {Point3D(0,1,0), Point3D(0,-1,0), Point3D(1,0,0), Point3D(-1,0,0)}
    VoronoiSphere(point_set)


