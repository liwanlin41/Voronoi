import re
import heapq
from enum import Enum

from point import Point
from BST import BST, Arc, BreakPoint # from https://rosettacode.org/wiki/AVL_tree#Python
from parabola import Parabola
from edge import Edge

### helper classes ###

class Graph:
    ''' represent a mutable graph on a collection of points with given vertices
    and edges'''

    def __init__(self, points, edges):
        ''' initialize a Graph object with a given set of vertices
        and a dictionary of edges represented as {point: adjacent set}'''
        self.point_set = points
        self.edge_dict = edges

    def insert_point(self, point):
        ''' insert a new isolated Point point into the graph'''
        if point not in self.point_set:
            self.point_set.add(point)
            self.edge_dict[point] = set()

    def insert_edge(self, point_a, point_b):
        ''' insert a new undirected edge (point_a, point_b) into the graph'''
        if point_a not in self.point_set:
            self.insert_point(point_a)
        if point_b not in self.point_set:
            self.insert_point(point_b)
        self.edge_dict[point_a].add(point_b)
        self.edge_dict[point_b].add(point_a)

    def get_neighbors(self, point):
        ''' return the set of points that are neighbors to point
        raise KeyError if point is not in this graph'''
        return self.edge_dict[point]

    def get_edges(self):
        ''' return the set of all edges in this graph, represented
        as a set of tuples'''
        all_edges = set()
        for point in self.point_set:
            for neighbor in self.edge_dict[point]:
                if point < neighbor: # self-loops not allowed
                    all_edges.add((point, neighbor))
        return all_edges

### main classes ###

class EventType(Enum):
    ''' represent the different kind of events that can occur
    during the algorithm '''
    SITE = 1
    CIRCLE = 2

class Event:
    ''' class representing events that occur as the sweep line descends'''

    def __init__(self, point, event_type, radius = 0):
        ''' create an event at location point and of type event_type 
        point is a Point and event_type is an EventType'''
        self.location = point
        self.event_type = event_type
        self.radius = radius

    def get_timing(self):
        ''' return a representation of when this event occurs during
        the sweep line algorithm '''
        # return the negative of the y coordinate because the sweep line
        # moves from top to bottom
        return -1 * (self.location.get_y() - self.radius)

    def __lt__(self, other):
        return self.get_timing() < other.get_timing()

    def handle(self, beachline, voronoi_edges):
        ''' handle this event, mutating the beachline as necessary
        adds edges and/or vertices to the dictionary voronoi_edges
        return list of new events '''
        if self.event_type == EventType.CIRCLE:
            return self.circle_handle(beachline, voronoi_edges)
        else:
            return self.site_handle(beachline, voronoi_edges)

    def circle_handle(self, beachline, voronoi_edges):
        ''' handle a circle event by modifying the beachline and graph
        remove the corresponding arc from the beachline
        return list of new events
        '''
        print(self.location)
        print(self.radius)
        raise NotImplementedError

    def site_handle(self, beachline, voronoi_edges):
        ''' handle a site event by modifying the beachline and graph
        add an arc to the beachline, possibly splitting arcs
        return list of new events
        '''
        ### beachline maintenance ###

        # currently the sweep line goes through self.location
        directrix = self.location.get_y()
        parabola = Parabola(self.location) # new parabola with this focus
        # remove and replace the beachline objects above this site
        nodes = beachline.delete(self.location.get_x(), self.location.get_y())

        if len(nodes) == 0: # no objects yet, first insert; handle this by itself
            arc = Arc(None, parabola, Parabola.NEG_INF(), Parabola.INF())
            beachline.insert(arc, directrix)
            return []
        
        # nodes is nicely sorted from left to right, may contain 1 or 3 objects
        foci = [] # list of foci of parabolas from left to right (by parabola order)
#        print("PRINTING NODES HERE")
#        for node in nodes:
#            print(node.pointer[1])
        if len(nodes) == 3: # this site lies below a breakpoint
            left_pointer = nodes[0].pointer
            right_pointer = nodes[2].pointer
            foci = left_pointer[:2] + [self.location] + right_pointer[1:]
        elif len(nodes) == 1:
            pointer = nodes[0].pointer # nodes[0] is a single arc
            foci = pointer[0:2] + [self.location] + pointer[1:]
        else:
            raise Exception("should not get here")
        new_nodes = [] # new arcs and breakpoints to be inserted; should be 3/2
        for i in range(1,len(foci) - 1): # len(foci) should be 4
            arc = Arc(None, Parabola(foci[i]), Parabola(foci[i-1]), Parabola(foci[i+1]))
#            print("New breakpoint added with left, right")
#            print(foci[i])
#            print(foci[i+1])
            breakpoint = BreakPoint(None, Parabola(foci[i]), Parabola(foci[i+1]))
#            print(breakpoint)
            new_nodes.append(arc)
            new_nodes.append(breakpoint)
        # this gives an extra duplicated breakpoint at the end of the list
        # NOTE: make sure to insert in left to right order to handle degeneracies
        for i in range(len(new_nodes)-1):
            beachline.insert(new_nodes[i], directrix)

        ### graph maintenance using voronoi_edges ###
        # first handle the common case
        if foci[1] == foci[3]: # site point breaks up a single arc, tracing new edge
            voronoi_edges[Edge(foci[1], foci[2])] = set()
        else: # site event coincides with circle event
            new_vertex = nodes[1].keyfunc(directrix)[0] # nodes[1] is a BreakPoint representing the collision point, coord is repeated
            voronoi_edges[Edge(foci[1], foci[3])].add(new_vertex)
            voronoi_edges[Edge(foci[1], foci[2])] = {new_vertex}
            voronoi_edges[Edge(foci[2], foci[3])] = {new_vertex}

        ### compute new events ###
        # the only potential circle events added by this vertex so far are
        # ones in which the arcs to the left and right disappear
        # circle event occurs when sweep line hits the bottom of the circle
        # recall foci has length 5 and contains the 5 foci associated with
        # the relevant parabolas here; degeneracy when this site event
        # is also a circle event is handled by the site event
        new_event_list = []
#        print("FOCUS HERE")
#        for focus in foci:
#            print(focus)
        # check for extreme cases
        if foci[0] != Point.NEG_INF():
            left_circle = compute_circumcenter(foci[0], foci[1], foci[2])
            # left_circle is where the arc to the left of this site disappears
            left_radius = foci[2].distance(left_circle)
            left_event = Event(left_circle, EventType.CIRCLE, left_radius)
            # only add if this event is new
            if not left_circle.get_y() - left_radius == directrix:
                new_event_list.append(left_event)
        if foci[4] != Point.INF():
            right_circle = compute_circumcenter(foci[2], foci[3], foci[4])
            right_radius = foci[2].distance(right_circle)
            right_event = Event(right_circle, EventType.CIRCLE, right_radius)
            if not right_circle.get_y() - right_radius == directrix:
                new_event_list.append(right_event)
        return new_event_list


class Voronoi:
    ''' a class for computing Voronoi diagrams in the plane 
    represented by site points, event queue, beachline'''

    def __init__(self, point_string):
        ''' create a new Voronoi object with given site points
        point_string is a string representation of the site points,
        taking the form "(a,b), (c,d) ..." with possible whitespace'''
        # strip parens and whitespace, leaving just the raw values
        # with possible empty spaces
        point_string_parse = re.split(r'[(,)\s]\s*', point_string)
        # extract values into list [x0, y0, x1, y1, ...]
        point_string_values = []
        for string in point_string_parse:
            if len(string) > 0:
                point_string_values.append(float(string))

        self.points = set()
        # turn raw values into point objects
        for i in range(len(point_string_values)//2):
            point = Point(point_string_values[2*i], point_string_values[2*i+1])
            self.points.add(point)

#        self.graph = Graph(set(),{}) # create empty graph
        self.voronoi_edges = {}
        self.event_queue = []
        self.beachline = BST() # BST representing beachline

        # initialize event queue with site events, sorted in
        # decreasing order of y coordinate
        for point in self.points:
            event = Event(point, EventType.SITE)
            heapq.heappush(self.event_queue, event)
#            heapq.heappush(self.event_queue, (-1 * point.get_y(), event))

    def output(self):
        ''' output the Voronoi diagram corresponding to these site points,
        represented as a Graph'''
        # TODO: see if this is the desired representation,
        # maybe directly draw the result instead
        edges = []
        for edge in self.voronoi_edges:
            endpoints = self.voronoi_edges[edge]
            edges.append(list(endpoints))
        return edges

    def get_next_event(self):
        ''' return the next event to be processed, popping the
        event off the queue'''
        # recall self.event_queue contains tuples of (timing, event)
        return heapq.heappop(self.event_queue)
#        return heapq.heappop(self.event_queue)[1]

    def step(self):
        ''' handle the next event '''
        next_event = self.get_next_event()
        new_events = next_event.handle(self.beachline, self.voronoi_edges)
        for event in new_events:
            heapq.heappush(self.event_queue, event)
#            heapq.heappush(self.event_queue, (event.get_timing(), event))

    def done(self):
        ''' return boolean for whether the algorithm has completed '''
        return (len(self.event_queue) == 0)

### helper functions ###

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
    voronoi = Voronoi(str(input()))
    while not voronoi.done(): # events left to handle
        voronoi.step()
    for tup in voronoi.output():
        for point in tup:
            print(point)

