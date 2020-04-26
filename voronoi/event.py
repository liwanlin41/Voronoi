from enum import Enum

from BST import BST, Arc, BreakPoint
from parabola import Parabola
from edge import Edge
from point import Point

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

    def handle(self, beachline, voronoi_edges, voronoi_vertices, verbose):
        ''' handle this event, mutating the beachline as necessary
        adds edges and/or vertices to the dictionary voronoi_edges

        adds to the dictionary voronoi_vertices (which is basically
        a copy of voronoi_edges but with segments marked by circumcenters
        instead of the actual coordinates; used for the spherical case)
        values of voronoi_vertices are (vertex_set, contains_midpoint)
        where vertex_set contains the 2D site point determining the
        voronoi vertex and contains_midpoint determines whether the voronoi edge
        contains the midpoint of the key edge

        return list of new events '''
        if verbose:
            print("OLD BEACHLINE HERE")
            print(beachline)
        if self.event_type == EventType.CIRCLE:
            return self.circle_handle(beachline, voronoi_edges, voronoi_vertices, verbose)
        else:
            return self.site_handle(beachline, voronoi_edges, voronoi_vertices, verbose)

    def circle_handle(self, beachline, voronoi_edges, voronoi_vertices, verbose):
        ''' handle a circle event by modifying the beachline and graph
        remove the corresponding arc from the beachline
        return list of new events
        '''
        # get current sweep line location
        directrix = self.location.get_y() - self.radius

#        print("CIRCLE EVENT HERE:")
#        print(self)
        # first check if this event actually occurs, i.e. parabolas do intersect here
        if len(beachline.find_exact(self.location.get_x(), directrix)) <= 3:
            if verbose:
                print(self.location)
                print("THIS EVENT DOES NOT OCCUR OR HAS ALREADY BEEN PROCESSED")
            return [] # do nothing
        # this line is for testing only
        exact_find = beachline.find_exact(self.location.get_x(), directrix)
        if len(exact_find) == 3:
            print("ALREADY PROCESSED")

        ### maintain beachline ###
        # remove and replace beachline objects above this site
        nodes = beachline.delete(self.location.get_x(), directrix)
        # if there are multiple intersections here, the only ones that
        # survive are the leftmost and rightmost arcs
        # nodes should take the form arc, breakpoint, arc, breakpoint, arc, etc.
        # where the arcs in the middle have all degenerated
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
        beachline.insert(left_arc, directrix)
        beachline.insert(breakpoint, directrix)
        beachline.insert(right_arc, directrix)

        ### graph maintenance using voronoi_edges and voronoi_vertices ###
        collisions = frozenset(foci[1:-1]) # all sites colliding here
#        if len(collisions) != 3:
#            print("voronoi vertex determined by")
#            for point in collisions:
#                print(point)
#            print("current location")
#            print(self.location)
#            print("exact find")
#            for node in exact_find:
#                print(node._str__())
#                print("Break points")
#                for break_point in node.keyfunc(directrix):
#                    print(break_point)
#            print("deleted nodes")
#            for node in nodes:
#                print(node._str__())
#            print("beachline")
#            print(beachline)
#            raise StopIteration
        # all arcs except the leftmost and rightmost disappear
        for i in range(1,len(foci)-2): # len(foci)-2 is the rightmost arc
            edge = Edge(foci[i], foci[i+1])
            other = foci[i-1] if i > 1 else foci[-2] # some vertex for comparison
            # foci[i], foci[i+1] are adjacent
            contains_midpoint = same_side(other, foci[i], foci[i+1], self.location)
            if edge in voronoi_edges:
                voronoi_edges[edge].add(self.location)
                voronoi_vertices[edge].add((collisions, contains_midpoint))
            elif foci[i] != foci[i+1]:
                voronoi_edges[edge] = {self.location}
                voronoi_vertices[edge] = {(collisions, contains_midpoint)}
        edge = Edge(foci[1], foci[-2])
        contains_midpoint = same_side(foci[2], foci[1], foci[-2], self.location)
        if edge in voronoi_edges:
            voronoi_edges[edge].add(self.location)
            voronoi_vertices[edge].add((collisions, contains_midpoint))
        elif foci[1] != foci[-2]:
            voronoi_edges[edge] = {self.location}
            voronoi_vertices[edge] = {(collisions, contains_midpoint)}

        ### event creation ###
        # new adjacent arcs are foci[0], foci[1], foci[-2] and foci[1], foci[-2], foci[-1]
        new_event_list = []
        left_point_set = {foci[0], foci[1], foci[-2]}
        right_point_set = {foci[1], foci[-2], foci[-1]}
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


    def site_handle(self, beachline, voronoi_edges, voronoi_vertices, verbose):
        ''' handle a site event by modifying the beachline and graph
        add an arc to the beachline, possibly splitting arcs
        return list of new events
        '''
        ### beachline maintenance ###

        # currently the sweep line goes through self.location
        directrix = self.location.get_y()
        parabola = Parabola(self.location) # new parabola with this focus
        # remove and replace the beachline objects above this site
        nodes = beachline.delete(self.location.get_x(), directrix)

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
        new_nodes = [] # new arcs and breakpoints to be inserted; should be 3/2 unless two points have the same y coord
        if self.location.get_y() == foci[1].get_y(): # same y coordinates, parabolas will only have one intersection
            # events are processed left to right, the rightmost arc will also have to change
            far_left = Parabola(foci[0])
            center_left = Parabola(foci[1]) # which is equal to foci[3]
            center_right = Parabola(foci[2])
            far_right = Parabola(foci[4])
            arc1 = Arc(None, center_left, far_left, center_right)
            breakpoint1 = BreakPoint(None, center_left, center_right)
            arc2 = Arc(None, center_right, center_left, far_right)
            breakpoint2 = BreakPoint(None, center_right, far_right)
            new_nodes.append(arc1)
            new_nodes.append(breakpoint1)
            new_nodes.append(arc2)
            new_nodes.append(breakpoint2)
            # the arc of far_right has not been removed, but its left parabola needs to be changed
            try:
                far_right_arc = nodes[-1].successor().successor()
                if verbose:
                    print("This is the far right arc")
                    print(far_right_arc)
                # manually reset for far_right_arc
                far_right_arc.pointer[0] = self.location
                def f(directrix):
                    left_breakpoint = Parabola.ordered_intersect(center_right, far_right)(directrix)
                    right_breakpoint = Parabola.ordered_intersect(far_right, Parabola(far_right_arc.pointer[-1]))(directrix)
                    return (left_breakpoint, right_breakpoint)
                # reset keyfunc as well
                far_right_arc.keyfunc = f
                assert far_right_arc.pointer[1] == foci[4]
                new_nodes.append(None) # because iteration will not see the last item
            except: # this point lies below the rightmost arc
                if verbose:
                    print("beneath rightmost arc")

        else: # non-degenerate case
            for i in range(1,len(foci) - 1): # len(foci) should be 5
                arc = Arc(None, Parabola(foci[i]), Parabola(foci[i-1]), Parabola(foci[i+1]))
#               print("New breakpoint added")
                breakpoint = BreakPoint(None, Parabola(foci[i]), Parabola(foci[i+1]))
#               print(breakpoint)
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
            assert foci[1] != foci[2]
            voronoi_edges[Edge(foci[1], foci[2])] = set()
            voronoi_vertices[Edge(foci[1], foci[2])] = set()
        else: # site event coincides with circle event
            new_vertex = nodes[1].keyfunc(directrix)[0] # nodes[1] is a BreakPoint representing the collision point, coord is repeated
#            print("ADDING EDGES AMONG:")
#            print(foci[1], foci[2], foci[3])
            collision = frozenset({foci[1], foci[2], foci[3]})
            voronoi_edges[Edge(foci[1], foci[3])].add(new_vertex)
            voronoi_vertices[Edge(foci[1], foci[3])].add(collision)
            voronoi_edges[Edge(foci[1], foci[2])] = {new_vertex}
            voronoi_edges[Edge(foci[2], foci[3])] = {new_vertex}
            voronoi_vertices[Edge(foci[1], foci[2])] = {collision}
            voronoi_vertices[Edge(foci[2], foci[3])] = {collision}

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
        if foci[4] != Point.INF():
            if not are_collinear(foci[2], foci[3], foci[4]):
                right_circle = compute_circumcenter(foci[2], foci[3], foci[4])
                right_radius = foci[2].distance(right_circle)
                right_event = Event(right_circle, EventType.CIRCLE, right_radius)
#                if not right_circle.get_y() - right_radius == directrix:
                if self < right_event:
#                    print("CIRCLE EVENT CREATED %s" %(right_event,))
                    new_event_list.append(right_event)
                    # recall new_nodes has an extra breakpoint at the end
        return new_event_list

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
    ''' given three non-collinear Points, return the Point representing their
    circumcenter'''
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

def same_side(q, p1, p2, r):
    ''' return whether q is on the same side of p1 p2 as r 
    q, p1, p2 must be non-collinear'''
    # write line in form ax + by + c = 0
    if len({q, p1, p2, r}) != 4:
        raise RuntimeError("expected 4 distinct points")
    a = p2.get_y() - p1.get_y()
    b = p1.get_x() - p2.get_x()
    c = -(a * p1.get_x() + b * p1.get_y())
    loc_r = a * r.get_x() + b * r.get_y() + c
    if loc_r == 0:
        return True
    sign_q = (a * q.get_x() + b * q.get_y() + c) > 0
    sign_r = loc_r > 0
    return sign_q == sign_r

