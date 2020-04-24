import numpy as np

class Point:
    ''' immutable class representing a point (theta, phi) on the unit sphere '''

    def __init__(self, theta, phi):
        ''' create new Point object with coordinates (theta, phi) 
        where theta is between 0 and 2pi and phi is between -pi/2, pi/2'''
        self.theta = theta
        self.phi = phi

    def convert_to_rectangular(self):
        ''' return tuple of Cartesian coordinates representing this point'''
        x = np.cos(self.theta) * np.cos(self.phi)
        y = np.sin(self.theta) * np.cos(self.phi)
        z = np.sin(self.phi)
        return (x, y, z)

    def get_theta(self):
        return self.theta

    def get_phi(self):
        return self.phi

    def distance(self, other):
        ''' compute great circle distance between self and other '''
        d_theta = self.get_theta() - other.get_theta()
        d_phi = self.get_phi() - other.get_phi()
        return 2 * np.arcsin(((1 - np.cos(d_theta) * np.cos(d_phi))/2)**0.5)

    def __eq__(self, other):
        return self.__class__ == other.__class__ and (self.theta, self.phi) == (other.get_theta(), other.get_phi())

    def __hash__(self):
        return hash((self.theta, self.phi))

    def __str__(self):
        return str((self.theta, self.phi))

    def __lt__(self, other):
        return (self.theta, self.phi) < (other.get_theta(), other.get_phi())

    def __le__(self, other):
        return (self.theta, self.phi) <= (other.get_theta(), other.get_phi())

if __name__ == '__main__':
    point1 = Point(0, 0)
    point2 = Point(0, np.pi/2)
    # manually check that these are approximately equal
    print(point1.distance(point2))
    print(np.pi/2)

