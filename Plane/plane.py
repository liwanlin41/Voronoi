import re
from enum import Enum

from point import Point
from BST import BST, Arc, BreakPoint # from https://rosettacode.org/wiki/AVL_tree#Python
from parabola import Parabola
from edge import Edge
from queue import PriorityQueue

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
    CIRCLE = 0
    SITE = 1

class Event:
    ''' class representing events that occur as the sweep line descends'''

    def __init__(self, point, event_type, radius = 0.0):
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
        # process left to right
        # in case of collision, process circle event first
        return (self.get_timing(), self.location, self.event_type.value) \
                < (other.get_timing(), other.location, other.event_type.value)

    def __eq__(self, other):
        return self.__class__ == other.__class__ and self.location == other.location and self.radius == other.radius

    def __hash__(self):
        return self.location.__hash__() + self.radius.__hash__()

    def __str__(self):
        return ("Event at %s of type %s and radius %f" %(self.location, self.event_type, self.radius))

    def handle(self, beachline, voronoi_edges):
        ''' handle this event, mutating the beachline as necessary
        adds edges and/or vertices to the dictionary voronoi_edges
        return list of new events '''
#        print("BEACHLINE HERE")
        print(beachline)
        if self.event_type == EventType.CIRCLE:
            return self.circle_handle(beachline, voronoi_edges)
        else:
            return self.site_handle(beachline, voronoi_edges)

    def circle_handle(self, beachline, voronoi_edges):
        ''' handle a circle event by modifying the beachline and graph
        remove the corresponding arc from the beachline
        return list of new events
        '''
        # get current sweep line location
        directrix = self.location.get_y() - self.radius

#        print("CIRCLE EVENT HERE:")
#        print(self)
        # first check if this event actually occurs, i.e. parabolas do intersect here
        if len(beachline.find_exact(self.location.get_x(), directrix)) == 0:
            print("THIS EVENT DOES NOT OCCUR")
            return [] # do nothing

        ### maintain beachline ###
        # remove and replace beachline objects above this site
        nodes = beachline.delete(self.location.get_x(), directrix)
        # if there are multiple intersections here, the only ones that
        # survive are the leftmost and rightmost arcs
        # nodes should take the form arc, breakpoint, arc, breakpoint, arc, etc.
        # where the arcs in the middle have all degenerated
#        print("CIRCLE EVENT NODES DELETED:")
#        for node in nodes:
#            print(node._str__())
#        print("BEACHLINE AFTER DELETION")
#        print(beachline)
        # get the list of foci associated with these beachline elements
        foci = [nodes[0].pointer[0]] # determines the parabola for the left breakpoint of the first arc
        for i in range(len(nodes)//2): # nodes has odd length
            breakpoint = nodes[2*i+1] # extract the breakpoints to get the foci
            foci.append(breakpoint.pointer[0]) # add only the left parabolas
        foci = foci + nodes[-1].pointer[1:] # add rightmost parabola + determinator of its right breakpoint
        # remember middle arcs all disappear, leaving only two arcs
        left_arc = Arc(None, Parabola(foci[1]), Parabola(foci[0]), Parabola(foci[-2]))
        breakpoint = BreakPoint(None, Parabola(foci[1]), Parabola(foci[-2]))
        right_arc = Arc(None, Parabola(foci[-2]), Parabola(foci[1]), Parabola(foci[-1]))
#        print("INSERTED NODES")
#        print(left_arc._str__())
#        print(breakpoint._str__())
#        print(right_arc._str__())
        beachline.insert(left_arc, directrix)
        beachline.insert(breakpoint, directrix)
        beachline.insert(right_arc, directrix)
#        print("NEW BEACHLINE")
#        print(beachline)

        ### graph maintenance using voronoi_edges ###
        # all arcs except the leftmost and rightmost disappear
        for i in range(1,len(foci)-2): # len(foci)-2 is the rightmost arc
            edge = Edge(foci[i], foci[i+1])
            if edge in voronoi_edges:
                voronoi_edges[edge].add(self.location)
            else:
                voronoi_edges[edge] = {self.location}
        edge = Edge(foci[1], foci[-2])
        if edge in voronoi_edges:
            voronoi_edges[edge].add(self.location)
        else:
            voronoi_edges[edge] = {self.location}

        ### event creation ###
        # new adjacent arcs are foci[0], foci[1], foci[-2] and foci[1], foci[-2], foci[-1]
        new_event_list = []
        left_point_set = {foci[0], foci[1], foci[-2]}
        right_point_set = {foci[1], foci[-2], foci[-1]}
#        print("ALL FOCI POTENTIALLY IN CIRCLE EVENT")
#        for focus in foci:
#            print(focus)
        if foci[0] != Point.NEG_INF() and len(left_point_set) == 3:
            if not are_collinear(foci[0], foci[1], foci[-2]):
                left_circle = compute_circumcenter(foci[0], foci[1], foci[-2])
                left_radius = left_circle.distance(foci[1])
                left_event = Event(left_circle, EventType.CIRCLE, left_radius)
                if self < left_event:
                    new_event_list.append(left_event)
        if foci[-1] != Point.INF() and len(right_point_set) == 3:
            if not are_collinear(foci[1], foci[-2], foci[-1]):
                right_circle = compute_circumcenter(foci[1], foci[-2], foci[-1])
                right_radius = right_circle.distance(foci[-2])
                right_event = Event(right_circle, EventType.CIRCLE, right_radius)
#                print("CIRCLE EVENT CREATED %s" %(right_event,))
                if self < right_event:
                    new_event_list.append(right_event)
        return new_event_list


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
#            print("NO NODES FOUND")
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
#            print("shouldn't get here!!!")
#            print(self)
#            for node in nodes:
#                print(node._str__())
            raise Exception("should not get here")
        new_nodes = [] # new arcs and breakpoints to be inserted; should be 3/2
        for i in range(1,len(foci) - 1): # len(foci) should be 4
            arc = Arc(None, Parabola(foci[i]), Parabola(foci[i-1]), Parabola(foci[i+1]))
#            print("New breakpoint added")
            breakpoint = BreakPoint(None, Parabola(foci[i]), Parabola(foci[i+1]))
#            print(breakpoint)
            new_nodes.append(arc)
            new_nodes.append(breakpoint)
        # this gives an extra duplicated breakpoint at the end of the list
        # NOTE: make sure to insert in left to right order to handle degeneracies
        for i in range(len(new_nodes)-1):
#            print("Inserting")
#            print(new_nodes[i])
            beachline.insert(new_nodes[i], directrix)
#            print("resulting beachline:")
#            print(beachline)
#            print("end beachline")

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
#        print("FOCUS FOR EVENTS FROM SITE HERE")
#        for focus in foci:
#            print(focus)
        # check for extreme cases
        if foci[0] != Point.NEG_INF():
            if not are_collinear(foci[0], foci[1], foci[2]):
                left_circle = compute_circumcenter(foci[0], foci[1], foci[2])
                # left_circle is where the arc to the left of this site disappears
                left_radius = foci[2].distance(left_circle)
                left_event = Event(left_circle, EventType.CIRCLE, left_radius)
                # only add if this event is new
#                if not left_circle.get_y() - left_radius == directrix:
                if self < left_event:
#                    print("CIRCLE EVENT CREATED %s" %(left_event))
                    new_event_list.append(left_event)
    #                # associate event with left arc, new_nodes[0]
    #                new_nodes[0].events.add(left_event)
        if foci[4] != Point.INF():
            if not are_collinear(foci[2], foci[3], foci[4]):
                right_circle = compute_circumcenter(foci[2], foci[3], foci[4])
                right_radius = foci[2].distance(right_circle)
                right_event = Event(right_circle, EventType.CIRCLE, right_radius)
#                if not right_circle.get_y() - right_radius == directrix:
                if self < right_event:
#                    print("CIRCLE EVENT CREATED %s" %(right_event,))
                    new_event_list.append(right_event)
    #                # recall new_nodes has an extra breakpoint at the end
    #                new_nodes[-2].events.add(right_event)
        return new_event_list


class Voronoi:
    ''' a class for computing Voronoi diagrams in the plane 
    represented by site points, event queue, beachline'''

#    def __init__(self, point_string):
#        ''' create a new Voronoi object with given site points
#        point_string is a string representation of the site points,
#        taking the form "(a,b), (c,d) ..." with possible whitespace'''
#        # strip parens and whitespace, leaving just the raw values
#        # with possible empty spaces
#        point_string_parse = re.split(r'[(,)\s]\s*', point_string)
#        # extract values into list [x0, y0, x1, y1, ...]
#        point_string_values = []
#        for string in point_string_parse:
#            if len(string) > 0:
#                point_string_values.append(float(string))
#
#        self.points = set()
#        # turn raw values into point objects
#        for i in range(len(point_string_values)//2):
#            point = Point(point_string_values[2*i], point_string_values[2*i+1])
#            self.points.add(point)

    def __init__(self, point_set):
        ''' create a new Voronoi object from given site points 
        point_set is a set of Point objects representing the sites '''
        # create own copy
        self.points = set()
        for point in point_set:
            self.points.add(point)

        self.graph = Graph(set(),{}) # create empty graph
        self.voronoi_edges = {}
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
        represented as a dictionary mapping edges to point lists '''
        # TODO: see if this is the desired representation,
        # maybe directly draw the result instead
#        edges = []
##        print("EDGE LIST")
#        for edge in self.voronoi_edges:
##            print(edge)
#            endpoints = self.voronoi_edges[edge]
#            edges.append(list(endpoints))
#        return edges
        return self.voronoi_edges

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
        print(next_event)
        print("Timing: %f" %(next_event.get_timing(),))
        new_events = next_event.handle(self.beachline, self.voronoi_edges)
        print("NEW EVENTS:")
        for event in new_events:
            print(event)
            self.event_queue.push(event)
#            heapq.heappush(self.event_queue, (event.get_timing(), event))
        if self.done(): # add something for edges going to infinity
            print("ADD FINAL ENDPOINTS")
            for node in self.beachline.traverse():
                if node.is_breakpoint():
                    edge = Edge(node.pointer[0], node.pointer[1])
                    cur_point = node.keyfunc(next_event.get_timing())[0]
                    if edge in self.voronoi_edges:
                        self.voronoi_edges[edge].add(cur_point)
                    else:
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
    point_set = {Point(0,0), Point(0,2), Point(0,-2), Point(2,0), Point(-2,0)}
    voronoi = Voronoi(point_set)
    while not voronoi.done(): # events left to handle
        print("------------------STEP---------------------")
        voronoi.step()
        stop = input()
    edge_dict = voronoi.output()
    for edge in edge_dict:
        output_string = edge.__str__() + ": "
        for point in edge_dict[edge]:
            output_string += str(point) + " "
        print(output_string)

