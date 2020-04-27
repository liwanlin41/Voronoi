from point import Point

# helper for error tolerance
def approx_equal(a, b, error = 1e-12):
    return b - error <= a <= b + error


class Point3D:
    ''' immutable class representing point (x, y, z) on the sphere '''

    def __init__(self, x, y, z):
        ''' create new Point object with coordinates (x, y, z) '''
        self.x = x
        self.y = y
        self.z = z
        # tolerate some error
#        assert 1-1e-10 <= x**2 + y**2 + z**2 <= 1+1e-10
    
    def invert(self):
        ''' invert through (0,0,1) and return the resulting 2DPoint 
        in the plane z = -1; this is projecting from (0,0,1)
        down onto the plane '''
        delta_z = -1 - self.z
        delta_y = delta_z * self.y / (self.z - 1)
        delta_x = delta_z * self.x / (self.z - 1)
        return Point(self.x + delta_x, self.y + delta_y)

    def invert_through(self, point):
        ''' invert through the 3D point point onto the plane tangent
        to the antipode of point, returning the (x, y, z) coordinates '''
        a, b, c = point.x, point.y, point.z
        # tangent plane has equation ax + by + cz = -1
        scale_factor = (a * self.x + b * self.y + c * self.z - 1)/-2
        x = (self.x - a) / scale_factor + a
        y = (self.y - b) / scale_factor + b
        z = (self.z - c) / scale_factor + c
        return [x, y, z]

    def distance(self, other):
        ''' return Euclidean distance between self and other '''
        return ((self.x - other.x)**2 + (self.y - other.y)**2 + (self.z - other.z)**2)**0.5

    def on_sphere(self):
        ''' return whether self is on the unit sphere '''
        # allow some error
        return 1-1e-10 <= self.distance(Point.origin()) <= 1+1e-10

    def origin():
        ''' return origin '''
        return Point3D(0,0,0)

    def __eq__(self, other):
        return approx_equal(self.x, other.x) and approx_equal(self.y, other.y) and approx_equal(self.z, other.z)

    def __hash__(self):
        round_x = int(self.x * 1e10) / 1e10
        round_y = int(self.y * 1e10) / 1e10
        round_z = int(self.z * 1e10) / 1e10
        return hash((round_x, round_y, round_z))

    def __str__(self):
        return str((self.x, self.y, self.z))

if __name__ == '__main__':
    p1 = Point3D(-0.567659576000496, 0.5040595140029139, 0.6509121385548771)
    p2 = Point3D(-0.5676595760004961, 0.5040595140029139, 0.6509121385548771)
    print(p1 == p2)
#    south = Point3D(0,0,-1)
#    print(south.invert())
#    print(south.invert().project_to_sphere())
#    point = Point3D(0.5**0.5, -0.5, 0.5)
#    print(point.invert())
#    print(point.invert().project_to_sphere())
#    print(south.invert_through(Point3D(0,1,0)))

