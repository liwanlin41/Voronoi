class Point:
    ''' immutable class representing a point (x,y) in the plane'''

    def __init__(self, x, y):
        ''' create a new Point object with coordinates (x, y)'''
        self.x = x
        self.y = y

    def NEG_INF():
        ''' return Point object infinitely far to the left '''
        return Point(float('-inf'), 0)

    def INF():
        ''' return Point object infinitely far to the right '''
        return Point(float('inf'), 0)

    def get_x(self):
        return self.x

    def get_y(self):
        return self.y

    def distance(self, other):
        ''' compute the Euclidean distance between self and another Point other'''
        return ((self.x - other.get_x())**2 + (self.y - other.get_y())**2)**0.5

    def project_to_sphere(self):
        ''' project onto unit sphere through (0,0,1), assuming self
        is the point (x, y, -1) '''
        # point must take the form r(x, y, -1) + s(0, 0, 1) for r + s = 1
        r = 4 / (4 + self.x**2 + self.y**2)
        return (r * self.x, r * self.y, 1 - 2 * r)

    def __eq__(self, other):
        return self.__class__ == other.__class__ and self.x == other.get_x() and self.y == other.get_y()

    def __hash__(self):
        return hash((self.x, self.y))

    def __str__(self):
        return str((self.x, self.y))

    def __lt__(self, other):
        return (self.x, self.y) < (other.get_x(), other.get_y())

    def __le__(self, other):
        return (self.x, self.y) <= (other.get_x(), other.get_y())

