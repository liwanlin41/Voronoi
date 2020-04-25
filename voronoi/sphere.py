from plane import Voronoi
from point import Point
from point3d import Point3D
from event import are_collinear

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
        for point in point_set:
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
            if eta_inv in self.plane_to_sphere: # this point exists
                weight = 1/8 if weight == 0 else weight/2
            else:
                break # good eta found
       eta_x, eta_y, eta_z = eta_inv.project_to_sphere()
       self.eta = Point3D(eta_x, eta_y, eta_z) # second center of inversion


if __name__ == '__main__':
    point_set = {Point3D(0,1,0), Point3D(0,-1,0), Point3D(1,0,0), Point3D(-1,0,0)}
    VoronoiSphere(point_set)


