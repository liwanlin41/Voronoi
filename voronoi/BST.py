# modified from https://rosettacode.org/wiki/AVL_tree#Python

from enum import Enum
from parabola import Parabola

# error value
#eps = 3.2e-6 # chosen because discriminant error has precision of 1e-11
eps = 1e-5

### helper function ###
def approx_equal(val1, val2, error = eps):
    ''' determine if val1 and val2 are equal up to a certain amount of error'''
    return val1 - error <= val2 <= val1 + error

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

    def __str__(self):
        return "Beachline components:\n" + self._string_() + "end"
    
    def _string_(self):
        output = ""
        if self.left:
            output += self.left._string_()
            output += "\n"
        output += self._str__()
        output += "\n"
        if self.right:
            output += self.right._string_()
        return output
    

    def is_breakpoint(self):
        ''' return true if this is a breakpoint '''
        return len(self.pointer) == 2

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
#        print("INSIDE FIND with x: %d, d: %d" %(x, d))
#        print("SEARCHING INSIDE %s" %(self._str__()))
#        print("breakpoints:")
#        for point in breakpoints:
#            print(point)
        nodes = []
        # first handle the case where x is found here
        if approx_equal(x, breakpoints[0].get_x()): # the predecessor will also have x
            # keep finding predecessors until we have everything above this point
            predecessor = self.predecessor()
            nodes.append(predecessor) # insert in backwards order
            while approx_equal(x, predecessor.keyfunc(d)[0].get_x()):
                predecessor = predecessor.predecessor()
                nodes.append(predecessor)
            nodes.reverse()
        # allow error tolerance
        if breakpoints[0].get_x()-eps <= x <= breakpoints[1].get_x()+eps:
            nodes.append(self)
        if approx_equal(x, breakpoints[1].get_x()):
            # same as above
            successor = self.successor()
            nodes.append(successor)
            while approx_equal(x, successor.keyfunc(d)[1].get_x()):
                successor = successor.successor()
                nodes.append(successor)
        if len(nodes) > 0:
            return nodes

        if breakpoints[0].get_x() > x:
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
#        print("FIND EXACT FOR COORDINATE %f" %(x,))
#        print(self._str__())
#        print("CURRENT BREAKPOINTS:")
#        print(breakpoints[0])
#        print(breakpoints[1])
        # first handle the case where x is found here
        if approx_equal(x, breakpoints[0].get_x()): # the predecessor will also have x
            # keep finding predecessors until we have everything above this point
            predecessor = self.predecessor()
            nodes.append(predecessor) # insert in backwards order
            while approx_equal(x, predecessor.keyfunc(d)[0].get_x()):
                predecessor = predecessor.predecessor()
                nodes.append(predecessor)
            nodes.reverse()
        if approx_equal(breakpoints[0].get_x(), x) or approx_equal(breakpoints[1].get_x(), x):
            nodes.append(self)
        if approx_equal(x, breakpoints[1].get_x()):
            # same as above
            successor = self.successor()
            nodes.append(successor)
            while approx_equal(x, successor.keyfunc(d)[1].get_x()):
                successor = successor.successor()
                nodes.append(successor)
        if len(nodes) > 0:
            return nodes

        if breakpoints[0].get_x() > x + eps:
            if self.left is None:
                return []
            return self.left.find_exact(x, d)
        if self.right is None:
            return []
        return self.right.find_exact(x, d)

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
        else: # swap with successor and delete
            s = self.successor()
            # first handle the case where s is the right child of self,
            # in which case the code below causes a loop; otherwise everything is fine
            if self is s.parent: # s cannot have any left children
                # left pointer relation
                self.left.parent = s
                s.left = self.left
                # parent relation
                if self is self.parent.left:
                    self.parent.left = s
                else:
                    self.parent.right = s
                s.parent = self.parent
                # nothing needs to be done on the right
                return self

            # change locations of pointers to self and s
            if self is self.parent.left:
                self.parent.left = s
            else:
                self.parent.right = s
            if s is s.parent.left:
                s.parent.left = self
            else:
               s.parent.right = self

            if self.left:
                self.left.parent = s
            if self.right:
                self.right.parent = s
            if s.left:
                s.left.parent = self
            if s.right:
                s.right.parent = self

            # switch pointers for self and s
            self.parent, s.parent = s.parent, self.parent
            self.left, s.left = s.left, self.left
            self.right, s.right = s.right, self.right
            return self.delete()

    def traverse(self):
        ''' iterate over the nodes of this subtree '''
        if self.left:
            for node in self.left.traverse():
                yield node
        yield self
        if self.right:
            for node in self.right.traverse():
                yield node


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

    def _str__(self):
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

    def _str__(self):
        return "[%s, %s]" %(self.pointer[0], self.pointer[1])


### helper functions for balancing the tree ###

def height(node):
    if node is None:
        return -1
    return node.height

def update_height(node):
    node.height = max(height(node.left), height(node.right)) + 1


class BST:
    ''' a binary search tree representing the beachline as a collection of
    arcs and breakpoints '''

    def __init__(self):
        ''' create an empty BST tree '''
        self.root = None

    def __str__(self):
        if self.root is None: return 'empty tree'
#        print("self.root.right here")
#        print(self.root.right)
        return str(self.root)

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
        current = node
        while current is not None:
#            print("REBALANCING AT CURRENT %s" %(current._str__()))
            update_height(current)
            if height(current.left) >= 2 + height(current.right): # left heavy
                if height(current.left.left) >= height(current.left.right):
                    self.right_rotate(current)
                else:
                    self.left_rotate(current.left)
                    self.right_rotate(current)
            elif height(current.right) >= 2 + height(current.left):
                if height(current.right.right) >= height(current.right.left):
                    self.left_rotate(current)
                else:
                    self.right_rotate(current.right)
                    self.left_rotate(current)
            if current is current.parent:
                raise Exception("self loop")
            current = current.parent # propagate up the tree

    def insert(self, node, d):
        ''' insert new node into tree at sweep line y=d'''
        if self.root is None: # empty tree
            self.root = node
        else:
            self.root.insert(node, d)
#        print("current beachline after insertion")
#        print(self)

    def delete(self, x, d):
        ''' delete the nodes located above x coord x when sweep line is at y = d
        return in sorted order'''
        nodes = self.find(x, d)
        for node in nodes:
            if node is None:
                continue
#            print("DELETING %s" %(node._str__()))
            if node is self.root:
                pseudoroot = BSTNode(None, 0) # placeholder for structure to track root
                pseudoroot.left = self.root
                self.root.parent = pseudoroot
                deleted = self.root.delete()
                self.root = pseudoroot.left
                if self.root is not None:
                    self.root.parent = None
#                print("REBALANCE AT ROOT")
                self.rebalance(self.root)
#                print("DELETED: %s" %(deleted._str__()))
            else:
                deleted = node.delete()
#                print("REBALANCE NOT AT ROOT")
                self.rebalance(deleted.parent)
#            print("Beachline after deletion of %s" %(node._str__()))
#            print(self)
        return nodes

    def traverse(self):
        ''' in-order traversal of tree '''
        for node in self.root.traverse():
            yield node

    def assert_invariant(self):
        ordered_components = []
        for node in self.traverse():
            ordered_components.append(node)
        for i in range(len(ordered_components)-1):
            assert len(ordered_components[i].pointer) == (3 if i%2 == 0 else 2)
            assert ordered_components[i].pointer[-2:] == ordered_components[i+1].pointer[:2]
            assert ordered_components[i].pointer[-2] != ordered_components[i].pointer[-1]
