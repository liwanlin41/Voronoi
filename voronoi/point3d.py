from point import Point

class Point3D:
    ''' immutable class representing point (x, y, z) on the sphere '''

    def __init__(self, x, y, z):
        ''' create new Point object with coordinates (x, y, z) '''
        self.x = x
        self.y = y
        self.z = z
        assert x**2 + y**2 + z**2 == 1

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
        return (x, y, z)

    def __eq__(self, other):
        return (self.x, self.y, self.z) == (other.x, other.y, other.z)

    def __hash__(self):
        return hash((self.x, self.y, self.z))

    def __str__(self):
        return str((self.x, self.y, self.z))

if __name__ == '__main__':
    south = Point3D(0,0,-1)
    print(south.invert())
    print(south.invert().project_to_sphere())
    point = Point3D(0.5**0.5, -0.5, 0.5)
    print(point.invert())
    print(point.invert().project_to_sphere())
    print(south.invert_through(Point3D(0,1,0)))

