from point import Point
from BST import BST
from edge import Edge
from queue import PriorityQueue
from event import Event, EventType

class Voronoi:
    ''' a class for computing Voronoi diagrams in the plane 
    represented by site points, event queue, beachline'''

    def __init__(self, point_set, verbose = True):
        ''' create a new Voronoi object from given site points 
        point_set is a set of Point objects representing the sites '''
        # create own copy
        self.points = set()
        self.verbose = verbose
        for point in point_set:
            self.points.add(point)

        self.voronoi_edges = {}
        self.voronoi_vertices = {}
        self.event_queue = PriorityQueue()
        self.beachline = BST() # BST representing beachline
        self.time = float('-inf') # current timing

        # initialize event queue with site events, sorted in
        # decreasing order of y coordinate
        for point in self.points:
            event = Event(point, EventType.SITE)
            self.event_queue.push(event)

    def output(self):
        ''' output the Voronoi diagram corresponding to these site points,
        represented as a dictionary mapping edges to point lists 
        and a dictionary mapping edges to circumcenters;
        different options depending on plane or spherical case'''
        return (self.voronoi_edges, self.voronoi_vertices)

    def get_next_event(self):
        ''' return the next event to be processed, popping the
        event off the queue
        also updates the current time'''
        # recall self.event_queue contains tuples of (timing, event)
        next_event = self.event_queue.pop()
        while next_event.get_timing() < self.time:
            next_event = self.event_queue.pop()
        self.time = next_event.get_timing()
        return next_event

    def step(self):
        ''' handle the next event '''
        next_event = self.get_next_event()
        if self.verbose:
            print("------------------STEP---------------------")
            print(next_event)
            print("Timing: %f" %(next_event.get_timing(),))
        new_events = next_event.handle(self.beachline, self.voronoi_edges, self.voronoi_vertices, self.verbose)
        if self.verbose:
            print("NEW EVENTS:")
            for event in new_events:
                print(event)
            print() # add extra line to aid readability
        for event in new_events:
            self.event_queue.push(event)
#        self.beachline.assert_invariant()
        if self.done(): # add something for edges going to infinity
            if self.verbose:
                print("FINAL BEACHLINE")
                print(self.beachline)
                print("ADD FINAL ENDPOINTS")
            for node in self.beachline.traverse():
                if node.is_breakpoint():
                    edge = Edge(node.pointer[0], node.pointer[1])
                    # compute a far out point for the sweep line
                    far_timing = next_event.get_timing()
                    if far_timing > 0: far_timing = -10 * far_timing - 10
                    else:
                        far_timing = (far_timing - 1) * 10
#                    cur_point = node.keyfunc(next_event.get_timing())[0]
                    cur_point = node.keyfunc(far_timing)[0]
                    if edge in self.voronoi_edges:
                        self.voronoi_edges[edge].add(cur_point)
                    else:
                        assert len(set(edge.get_sites())) == 2
                        self.voronoi_edges[edge] = {cur_point}

    def done(self):
        ''' return boolean for whether the algorithm has completed '''
        return (len(self.event_queue) == 0)

### helper functions ###

def are_collinear(p1, p2, p3):
    ''' determine if points p1, p2, p3 are collinear '''
    # first compute a good order of the points to avoid division by zero
    sorted_points = [p1, p2, p3]
    sorted_points.sort(key = lambda p: p.get_y())
    b = sorted_points[1]
    if sorted_points[0].get_y() == sorted_points[2].get_y(): # horizontally collinear
        return True
    if sorted_points[0].get_y() != sorted_points[1].get_y():
        a = sorted_points[0]
        c = sorted_points[2]
    else:
        a = sorted_points[2]
        c = sorted_points[0]
    # compute circumcenter
    slope_ab = (b.get_x() - a.get_x())/(a.get_y() - b.get_y())
    slope_ca = (c.get_x() - a.get_x())/(a.get_y() - c.get_y())
    return slope_ab == slope_ca

def compute_circumcenter(p1, p2, p3):
    ''' given three non-collinear Points, return the Point representing their circumcenter '''
    # first compute a good order of the points to avoid division by zero
    sorted_points = [p1, p2, p3]
    sorted_points.sort(key = lambda p: p.get_y())
    b = sorted_points[1]
    if sorted_points[0].get_y() != sorted_points[1].get_y():
        a = sorted_points[0]
        c = sorted_points[2]
    else:
        a = sorted_points[2]
        c = sorted_points[0]
    # compute circumcenter
    slope_ab = (b.get_x() - a.get_x())/(a.get_y() - b.get_y())
    intercept_ab = (a.get_y() + b.get_y())/2 - (a.get_x() + b.get_x())/2 * slope_ab
    slope_ca = (c.get_x() - a.get_x())/(a.get_y() - c.get_y())
    intercept_ca = (a.get_y() + c.get_y())/2 - (a.get_x() + c.get_x())/2 * slope_ca
    x_intersection = (intercept_ca - intercept_ab)/(slope_ab - slope_ca)
    y_intersection = slope_ab * x_intersection + intercept_ab
    return Point(x_intersection, y_intersection)

### actually execute ###

if __name__ == '__main__':
    # create new Voronoi object
    print("Enter list of site points in general position:")
#    test_string = "(0,0), (0.000000001,2), (-0.00000001,-2), (2,0.00000001), (-2,-0.00000001)"
#    test_string = "(0,0), (0,2), (0,-2), (2,0), (-2,0)"
#    voronoi = Voronoi(str(input()))
#    point_set = {Point(0,0), Point(0,2), Point(0,-2), Point(2,0), Point(-2,0)}
    point_set = {Point(-1.5,7.1), Point(-4.1, 4.1), Point(1.8,2), Point(-5.9,-2.3)}
    voronoi = Voronoi(point_set)
    while not voronoi.done(): # events left to handle
#        print("------------------STEP---------------------")
        voronoi.step()
        stop = input()
    edge_dict = voronoi.output()
    for edge in edge_dict:
        output_string = edge.__str__() + ": "
        for point in edge_dict[edge]:
            output_string += str(point) + " "
        print(output_string)

