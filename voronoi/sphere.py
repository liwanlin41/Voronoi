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
       # if x' are new coordinates, matrix_transform @ x' are old coords
       # new coordinatnes have z = -1 always
       matrix_transform = np.array([[b, -a*c, a],[-a,-b*c,b],[0,a**2+b**2,c]])
       inverse_transform = np.linalg.inv(matrix_transform)
       # track the image of (0,0,1)
       q = np.array(Point3d(0,0,1).invert_through(self.eta)).reshape((3,1))
       inverted_q = (inverse_transform @ q).reshape((3,))
       # image of (0,0,1) under second inversion
       self.q = Point(inverted_q[0], inverted_q[1]) 
       for point in self.points: # these are 3d points
           inverted = np.array(point.invert_through(self.eta)).reshape((3,1))
           new_coord = inverse_transform @ inverted
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

    def find_far_section(self):
        ''' find the Voronoi region containing the image of (0,0,1)
        under inversion about eta '''
        # first find site point of this Voronoi region
        site_point = None
        cur_distance = None
        for point in self.plane2_to_sphere: # point is a 2D point
            distance = point.distance(self.q)
            if site_point is None:
                site_point = point
                cur_distance = distance
            elif distance < cur_distance:
                cur_distance = distance
                site_point = point
        # site_point now contains the correct site point
        for edge in self.voronoi2.voronoi_vertices:
            plane_sites = edge.get_sites()
            if site_point in plane_sites:
                # get the points on the sphere
                sphere_site1 = self.plane2_to_sphere[plane_sites[0]]
                sphere_site2 = self.plane2_to_sphere[plane_sites[1]]
                sphere_edge = Edge(sphere_site1, sphere_site2)
                for circle_set in self.voronoi2.voronoi_vertices[edge]:
                    set_iterator = iter(circle_set)
                    circle_list = []
                    for i in range(3):
                        plane_focus = next(set_iterator)
                        circle_list.append(self.plane2_to_sphere[plane_focus])
                    # circle_list contains the three points on the sphere determining
                    # the Voronoi vertex on the sphere
                    # convert directly to coordinates
                    vertex = compute_center(circle_list[0], circle_list[1], circle_list[2], self.eta)
                    if sphere_edge in self.voronoi_edges:
                        self.voronoi_edges[sphere_edge].add(vertex)
                    else:
                        self.voronoi_edges[sphere_edge] = {vertex}
        # now self.voronoi_edges contains the preimage of the entire section surrounding q = (0,0,1)
                    

    def find_near_section(self):
        ''' lift the z=-1 inverted image back onto the sphere '''
        for edge in self.voronoi1.voronoi_vertices:
            plane_sites = edge.get_sites()
            sphere_site1 = self.plane_to_sphere[plane_sites[0]]
            sphere_site2 = self.plane_to_sphere[plane_sites[1]]
            sphere_edge = Edge(sphere_site1, sphere_site2)
            for circle_set in self.voronoi1.voronoi_vertices[edge]:
                set_iterator = iter(circle_set)
                # extract just a triangle
                circle_list = []
                for i in range(3):
                    plane_focus = next(set_iterator)
                    circle_list.append(self.plane_to_sphere[plane_focus])
                vertex = compute_center(circle_list[0], circle_list[1], circle_list[2])
                if sphere_edge in self.voronoi_edges:
                    self.voronoi_edges[sphere_edge].add(vertex)
                else:
                    self.voronoi_edges[sphere_edge] = {vertex}

    def output(self):
        ''' combine the two Voronoi diagrams to get the output '''
        # store as dictionary of edge: circumcenters
        self.voronoi_edges = {}
        self.find_far_section()
        self.find_near_section()

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


