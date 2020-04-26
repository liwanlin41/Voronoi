from plane import Voronoi
from point import Point
from point3d import Point3D
from event import are_collinear, compute_circumcenter
import numpy as np
from edge import Edge

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
        triangle = [next(planar_iter) for i in range(3)]
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
        self.inverse_transform = inverse_transform # store for later
        # track the image of (0,0,1)
        q = np.array(Point3D(0,0,1).invert_through(self.eta)).reshape((3,1))
        inverted_q = (inverse_transform @ q).reshape((3,))
        # image of (0,0,1) under second inversion
        self.q = Point(inverted_q[0], inverted_q[1]) 
        for point in self.points: # these are 3d points
            inverted = np.array(point.invert_through(self.eta)).reshape((3,1))
            new_coords = inverse_transform @ inverted
            new_coords.reshape((3,))
            inverted_point = Point(new_coords[0][0], new_coords[1][0])
            self.sphere_to_plane2[point] = inverted_point
            self.plane2_to_sphere[inverted_point] = point

        # construct the voronoi diagrams
        self.voronoi1 = Voronoi(set(self.plane_to_sphere.keys()), verbose)
        self.voronoi2 = Voronoi(set(self.plane2_to_sphere.keys()), verbose)
#        self.voronoi2 = Voronoi(set(self.plane2_to_sphere.keys()).add(self.q), verbose)

    def step(self):
        ''' handle the next events concurrently '''
        if not self.voronoi1.done():
            self.voronoi1.step()
        if not self.voronoi2.done():
            self.voronoi2.step()

    def find_far_section(self):
        ''' find the Voronoi region containing the image of (0,0,1)
        under inversion about eta '''
        voronoi_edges = {}
        # voronoi_edges will be a dictionary of sphere edge: {vertices}, contains_midpoint
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
                    vertex = compute_center(circle_list[0], circle_list[1], circle_list[2], self.eta, self.inverse_transform, self.sphere_to_plane2)
                    if sphere_edge in voronoi_edges:
                        voronoi_edges[sphere_edge].add(vertex)
                    else:
                        voronoi_edges[sphere_edge] = {vertex}
        return voronoi_edges
        # now voronoi_edges contains the preimage of the entire section surrounding q = (0,0,1)

    def find_near_section(self):
        ''' lift the z=-1 inverted image back onto the sphere '''
        voronoi_edges = {}
        # voronoi_edges will be a dictionary of sphere edge: {vertices}, contains_midpoint
        for edge in self.voronoi1.voronoi_vertices:
            plane_sites = edge.get_sites()
            sphere_site1 = self.plane_to_sphere[plane_sites[0]]
            sphere_site2 = self.plane_to_sphere[plane_sites[1]]
            sphere_edge = Edge(sphere_site1, sphere_site2)
            # iterate over circumcenters as sets of 2d points and convert to coords
            for circle_set in self.voronoi1.voronoi_vertices[edge]:
                set_iterator = iter(circle_set)
                # extract just a triangle
                circle_list = []
#                print("CIRCLE SET HERE")
#                print(circle_set)
                try:
                    for i in range(3):
                        plane_focus = next(set_iterator)
                        circle_list.append(self.plane_to_sphere[plane_focus])
                except StopIteration: # hacky fix for now
                    print("NOT ENOUGH POINTS") # this is not supposed to happen?
                    continue
                vertex = compute_center(circle_list[0], circle_list[1], circle_list[2])
                if sphere_edge in voronoi_edges:
                    voronoi_edges[sphere_edge].add(vertex)
                else:
                    voronoi_edges[sphere_edge] = {vertex}
        # endpoints at infinity to join with the second voronoi structure
        for edge in voronoi_edges:
            if len(voronoi_edges[edge]) == 1:
                site1, site2 = edge.get_sites()
                vertex = compute_infinite_center(site1, site2)
                voronoi_edges[edge].add(vertex)
        return voronoi_edges
#        join_points = {}
#        for edge in voronoi_edges: # handle things at infinity
#            if len(voronoi_edges[edge]) == 1:
#                site1, site2 = edge.get_sites()
#                vertex = compute_center(site1, site2, Point3D(0,0,1))
#                join_points[edge] = vertex
#        return voronoi_edges, join_points

    def output(self):
        ''' return a coordinate representation of the two Voronoi diagrams to get the output '''
        # store as dictionary of edge: circumcenters
#        near_edges, join_points = self.find_near_section()
        near_edges = self.find_near_section()
        far_edges = self.find_far_section()
        return (far_edges, near_edges)

    def done(self):
        return self.voronoi1.done() and self.voronoi2.done()

### helper functions ###

def is_in_circle(p, p1, p2, p3):
    ''' given 4 2d points, determine if point p is inside the circle
    determined by p1, p2, p3'''
    center = compute_circumcenter(p1, p2, p3)
    return center.distance(p) <= center.distance(p1)

def compute_center(p1, p2, p3, invert = None, inverse_transform = None, sphere_to_plane = None):
    ''' given three Point3Ds, return Point3D on sphere that is equidistant
    from all of them and whose image is inside the circumcircle
    of the inverted triangle 
    invert through (0,0,1) by default'''
    # compute normal to the plane using cross product to get equidistant line
    vec1 = np.array([p2.x - p1.x, p2.y - p1.y, p2.z - p1.z])
    vec2 = np.array([p3.x - p1.x, p3.y - p1.y, p3.z - p1.z])
#    normal = np.cross(vec1, vec2).reshape((1,3))
    normal = np.cross(vec1, vec2)
    mag = np.linalg.norm(normal)
    coords = normal / mag # one of the circumline points on the sphere
#    print("NORMAL VECTOR")
#    print(normal)
    centerpoint = Point3D(coords[0], coords[1], coords[2])
    if invert:
        inverted_center_arr = np.array(centerpoint.invert_through(invert)).reshape((3,1))
        p_coords = (inverse_transform @ inverted_center_arr).reshape((3,))
        p = Point(p_coords[0], p_coords[1])
        im_p1 = sphere_to_plane[p1]
        im_p2 = sphere_to_plane[p2]
        im_p3 = sphere_to_plane[p3]
    else: # do the same for normal inversion
        p = centerpoint.invert()
        im_p1 = p1.invert()
        im_p2 = p2.invert()
        im_p3 = p3.invert()
    if is_in_circle(p, im_p1, im_p2, im_p3):
        return centerpoint
    return Point3D(-coords[0], -coords[1], -coords[2])
#    d = normal @ np.array([[p1.x, p1.y, p1.z]]).T
    # d represents ax + by + cz = d
#    invert_point = np.array([[invert.x, invert.y, invert.z]])
#    mag = np.linalg.norm(normal)
#    invert_sign = np.sign(normal @ invert_point.T - d)
#    if invert_sign == 0 or invert_sign == np.sign(-d):
#        # invert on the same side as the origin or on the plane
#        coords = normal/mag
#    else:
#        coords = -normal/mag
#    return Point3D(coords[0][0], coords[0][1], coords[0][2])

def compute_infinite_center(p1, p2):
    ''' given two 3d points p1, p2, return the circumcenter with (0,0,1)
    representing the endpoint of this edge toward infinity 
    not guaranteed to be correct if p1, p2, (0,0,1), (0,0,0) are coplanar'''
    # mostly copied from above
    vec1 = np.array([p2.x - p1.x, p2.y - p1.y, p2.z - p1.z])
    vec2 = np.array([- p1.x, - p1.y, 1 - p1.z])
    normal = np.cross(vec1, vec2)
    mag = np.linalg.norm(normal)
    coords = normal / mag 
    if coords[2] >= 0:
#        print("center point %f, %f, %f" %(coords[0], coords[1], coords[2]))
        return Point3D(coords[0], coords[1], coords[2])
#    print("center point %f, %f, %f" %(-coords[0], -coords[1], -coords[2]))
    return Point3D(-coords[0], -coords[1], -coords[2])


if __name__ == '__main__':
#    point_set = {Point3D(0,1,0), Point3D(0,-1,0), Point3D(1,0,0), Point3D(-1,0,0)}
#    VoronoiSphere(point_set)

#    p = Point(0,1.5)
#    p1 = Point(1,0)
#    p2 = Point(0,1)
#    p3 = Point(-1,0)
#    print(is_in_circle(p, p1, p2, p3))

#    p1 = Point3D(0,1,0)
#    p2 = Point3D(-0.8,0.6,0)
#    p3 = Point3D(0.64,0.48,0.6)
#    print(compute_center(p1, p2, p3))
    p1 = Point3D(-0.726667, -0.121111, -0.676230)
    p2 = Point3D(-0.210467, -0.701558, 0.680823)
    print(compute_infinite_center(p1, p2))



