# modified from https://rosettacode.org/wiki/AVL_tree#Python

from enum import Enum
from parabola import Parabola

class BSTNode:
    ''' node in the BST, may be an arc or a breakpoint '''

    def __init__(self, parent, f):
        ''' create a new node
        f is a function with input y (current sweep line location)
        and output the current range of this node as a tuple;
        for arcs this output is (left breakpoint, right breakpoint)
        and for breakpoints this is (breakpoint, breakpoint)
        '''
        self.parent = parent
        self.left = None
        self.right = None
        self.keyfunc = f
        self.height = 0

    def is_breakpoint(self):
        ''' return true if this is a breakpoint '''
        # breakpoints have left and right breakpoints the same everywhere
        # it is enough to check at two distinct sweep line locations
        points0 = self.keyfunc(0)
        points1 = self.keyfunc(1)
        return points0[0] == points0[1] and points1[0] == points1[1]

    def successor(self):
        ''' find the successor of this node '''
        if self.right is not None:
            successor = self.right
            while successor.left is not None:
                successor = successor.left
            return successor
        # otherwise find the first node for which this is to the left
        uppernode = self
        # keep going up to the left
        while uppernode.parent is not None and uppernode is uppernode.parent.right:
            uppernode = uppernode.parent
        return uppernode.parent # may be None

    def predecessor(self):
        ''' find the predecessor of this node '''
        if self.left is not None:
            predecessor = self.left
            while predecessor.right is not None:
                predecessor = predecessor.right
            return predecessor
        uppernode = self
        while uppernode.parent is not None and uppernode is uppernode.parent.left:
            uppernode = uppernode.parent
        return uppernode.parent

    def find(self, x, d):
        ''' when the sweep line is at y=d, find and return the 
        node/nodes (up to 3, as a list) whose x-coord intervals contain x'''
        # normally x will be inside an arc
        # the issue is if x is a breakpoint; then we also want the predecessor
        # and successor
        breakpoints = self.keyfunc(d)
#        print("NODE AT (parabola or right_parabola)")
#        print(self.pointer[1])
#        print("INSIDE FIND with x: %d, d: %d" %(x, d))
#        print("breakpoints:")
#        for point in breakpoints:
#            print(point)
        nodes = []
        # first handle the case where x is found here
        if x == breakpoints[0].get_x(): # the predecessor will also have x
            # keep finding predecessors until we have everything above this point
            predecessor = self.predecessor()
            nodes.append(predecessor) # insert in backwards order
            while x == predecessor.keyfunc(d)[0].get_x():
                predecessor = predecessor.predecessor()
                nodes.append[predecessor]
            nodes.reverse()
        if breakpoints[0].get_x() <= x <= breakpoints[1].get_x():
            nodes.append(self)
        if x == breakpoints[1].get_x():
            # same as above
            successor = self.successor()
            nodes.append(successor)
            while x == successor.keyfunc(d)[1].get_x():
                successor = successor.successor()
                nodes.append(successor)
        if len(nodes) > 0:
            return nodes

        if breakpoints[0].get_x() < x:
            if self.left is None:
                return []
            return self.left.find(x, d)
        if self.right is None:
            return []
        return self.right.find(x, d)

    def find_exact(self, x, d):
        ''' as opposed to find, find and return only the nodes with a breakpoint
        with x-coordinate x at sweep line location d'''
        breakpoints = self.keyfunc(d)
        nodes = []
        # first handle the case where x is found here
        if x == breakpoints[0].get_x(): # the predecessor will also have x
            # keep finding predecessors until we have everything above this point
            predecessor = self.predecessor()
            nodes.append(predecessor) # insert in backwards order
            while x == predecessor.keyfunc(d)[0].get_x():
                predecessor = predecessor.predecessor()
                nodes.append[predecessor]
            nodes.reverse()
        if breakpoints[0].get_x() == x or breakpoints[1].get_x() == x:
            nodes.append(self)
        if x == breakpoints[1].get_x():
            # same as above
            successor = self.successor()
            nodes.append(successor)
            while x == successor.keyfunc(d)[1].get_x():
                successor = successor.successor()
                nodes.append(successor)
        if len(nodes) > 0:
            return nodes

        if breakpoints[0].get_x() < x:
            if self.left is None:
                return []
            return self.left.find(x, d)
        if self.right is None:
            return []
        return self.right.find(x, d)

    def insert(self, node, d):
        ''' insert node into subtree rooted here when sweep line is at y=d
        any arc breakings have already happened, so all (left, right) breakpoint
        pairs will be distinct and ordered '''
        self_breakpoints = self.keyfunc(d)
        node_breakpoints = node.keyfunc(d)
#        print("inserting node %s" %(node,))
#        print("self breakpoints %s, %s" %(self_breakpoints[0], self_breakpoints[1]))
#        print("node breakpoints %s, %s" %(node_breakpoints[0], node_breakpoints[1]))
        if node_breakpoints[0].get_x() < self_breakpoints[0].get_x():
            if self.left is None:
                node.parent = self
                self.left = node
#                print("inserted to left of %s" %(self,))
            else:
                self.left.insert(node, d)
        else:
            if self.right is None:
                node.parent = self
                self.right = node
#                print("inserted to right of %s" %(self,))
            else:
                self.right.insert(node, d)

    def delete(self):
        ''' delete this node from tree, return it for rebalancing '''
        if self.left is None or self.right is None:
            if self is self.parent.left:
                self.parent.left = self.left or self.right
                if self.parent.left is not None:
                    self.parent.left.parent = self.parent
            else:
                self.parent.right = self.left or self.right
                if self.parent.right is not None:
                    self.parent.right.parent = self.parent
            return self
        else:
            s = self.successor()
            self.keyfunc, s.keyfunc = s.keyfunc, self.keyfunc
            self.pointer, s.pointer = s.pointer, self.pointer
#            self.events, s.events = s.events, self.events
            return s.delete()

### helper functions for balancing the tree ###

def height(node):
    if node is None:
        return -1
    return node.height

def update_height(node):
    node.height = max(height(node.left), height(node.right)) + 1

class Arc(BSTNode):
    ''' BSTNode representing an arc '''

    def __init__(self, parent, parabola, left_parabola, right_parabola):
        ''' portion of arc parabola defined by left and right parabolas '''
        def f(directrix): # return (left_breakpoint, right_breakpoint) as function of directrix
            left_breakpoint = Parabola.ordered_intersect(left_parabola, parabola)(directrix)
            right_breakpoint = Parabola.ordered_intersect(parabola, right_parabola)(directrix)
            return (left_breakpoint, right_breakpoint)
        BSTNode.__init__(self, parent, f)
        # points to the site associated with this arc, as tuple
        # for convenience
        self.pointer = [left_parabola.get_focus(), parabola.get_focus(), right_parabola.get_focus()]
#        # set of associated circle events
#        self.events = set()

    def __str__(self):
        return "[%s, %s, %s]" %(self.pointer[0], self.pointer[1], self.pointer[2])

class BreakPoint(BSTNode):
    ''' BSTNode representing a breakpoint'''

    def __init__(self, parent, left_parabola, right_parabola):
        ''' create breakpoint at the left/right (lr) intersection 
        of left and right parabolas '''
        def f(directrix):
            breakpoint = Parabola.ordered_intersect(left_parabola, right_parabola)(directrix)
            return (breakpoint, breakpoint)
        BSTNode.__init__(self, parent, f)
        # points to the points for which this point traces a Voronoi edge
        self.pointer = [left_parabola.get_focus(), right_parabola.get_focus()]
#        self.events = set() # this will always be empty

    def __str__(self):
        return "[%s, %s]" %(self.pointer[0], self.pointer[1])

class BST:
    ''' a binary search tree representing the beachline as a collection of
    arcs and breakpoints '''

    def __init__(self):
        ''' create an empty BST tree '''
        self.root = None

    def find(self, x, d):
        ''' see BSTNode find'''
        if self.root is None:
            return []
        return self.root.find(x, d)

    def find_exact(self, x, d):
        ''' see BSTNode find_exact'''
        if self.root is None:
            return []
        return self.root.find_exact(x,d)

    def left_rotate(self, x):
        y = x.right
        y.parent = x.parent
        if x.parent is None:
            self.root = y
        else:
            if y.parent.left is x:
                y.parent.left = y
            elif y.parent.right is x:
                y.parent.right = y
        x.right = y.left
        if y.left is not None:
            x.right.parent = x
        y.left = x
        x.parent = y
        update_height(x)
        update_height(y)

    def right_rotate(self, x):
        y = x.left
        y.parent = x.parent
        if x.parent is None:
            self.root = y
        else:
            if x.parent.left is x:
                x.parent.left = y
            elif x.parent.right is x:
                x.parent.right = y
        x.left = y.right
        if y.right is not None:
            y.right.parent = x
        y.right = x
        x.parent = y
        update_height(x)
        update_height(y)

    def rebalance(self, node):
        while node is not None:
            update_height(node)
            if height(node.left) >= 2 + height(node.right): # left heavy
                if height(node.left.left) >= height(node.left.right):
                    self.right_rotate(node)
                else:
                    self.left_rotate(node.left)
                    self.right_rotate(node)
            elif height(node.right) >= 2 + height(node.left):
                if height(node.right.right) >= height(node.right.left):
                    self.left_rotate(node)
                else:
                    self.right_rotate(node.right)
                    self.left_rotate(node)
            node = node.parent # propagate up the tree

    def insert(self, node, d):
        ''' insert new node into tree at sweep line y=d'''
        if self.root is None: # empty tree
            self.root = node
        else:
            self.root.insert(node, d)

    def delete(self, x, d):
        ''' delete the nodes located above x coord x when sweep line is at y = d
        return in sorted order'''
        nodes = self.find(x, d)
        for node in nodes:
            if node is None:
                continue
            if node is self.root:
                pseudoroot = BSTNode(None, 0) # placeholder for structure
                pseudoroot.left = self.root
                self.root.parent = pseudoroot
                deleted = self.root.delete()
                self.root = pseudoroot.left
                if self.root is not None:
                    self.root.parent = None
                self.rebalance(self.root)
            else:
                deleted = node.delete()
                self.rebalance(deleted.parent)
        return nodes
