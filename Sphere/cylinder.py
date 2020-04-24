import numpy as np

class Point:
    ''' immutable class representing a point (theta, z) on the unit sphere in cylindrical coordinates'''

    def __init__(self, theta, z):
        ''' create new Point object with coordinates (theta, z) 
        where theta is between 0 and 2pi and z is between -1, 1'''
        self.theta = theta
        self.z = z

    def convert_to_rectangular(self):
        ''' return tuple of Cartesian coordinates representing this point'''
        radius = (1-self.z**2)**0.5
        x = np.cos(self.theta) * radius
        y = np.sin(self.theta) * radius
        return (x, y, self.z)

    def get_theta(self):
        return self.theta

    def get_z(self):
        return self.z

    def euclidean_distance(self, other):
        ''' compute Euclidean distance between self and other '''
        x, y, z = self.convert_to_rectangular()
        a, b, c = other.convert_to_rectangular()
        return ((x-a)**2 + (y-b)**2 + (c-z)**2)**0.5

    def distance(self, other):
        ''' compute great circle distance between self and other '''
        return 2 * np.arcsin(self.euclidean_distance(other)/2)

    def __eq__(self, other):
        return self.__class__ == other.__class__ and (self.theta, self.z) == (other.get_theta(), other.get_z())

    def __hash__(self):
        return hash((self.theta, self.z))

    def __str__(self):
        return str((self.theta, self.z))

    def __lt__(self, other):
        return (self.theta, self.z) < (other.get_theta(), other.get_z())

    def __le__(self, other):
        return (self.theta, self.z) <= (other.get_theta(), other.get_z())

if __name__ == '__main__':
    point1 = Point(0, 0)
    point2 = Point(0, 1)
    # manually check that these are approximately equal
    print(point1.distance(point2))
    print(np.pi/2)

