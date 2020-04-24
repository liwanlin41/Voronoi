#import numpy as np
from cylinder import Point
from sympy import sin, cos, pi, Symbol, solve, Eq, evalf

class Locus:
    ''' represent locus of points equidistant from point and sweep plane '''

    def __init__(self, focus):
        ''' create a new locus focused at a given site '''
        self.focus = focus

    def get_focus(self):
        return self.focus

    def get_locus(self, z):
        ''' get the equation of the locus when sweep plane is at z 
        when evaluated at a site event z0 = z, this is not necessarily
        accurate and just returns z '''
        def f(theta):
            # handle degenerate north/south pole cases
            x0, y0, z0 = self.focus.convert_to_rectangular()
            if z0 == z: return z
            var_z = Symbol('z', real = True)
            eq = Eq(((1-var_z**2)**0.5 * cos(theta) - x0) ** 2 + \
                    ((1-var_z**2)**0.5 * sin(theta) - y0) ** 2 + \
                    (z0 - var_z) ** 2 + 2 * (1-var_z**2)**0.5 * (1-z**2)**0.5 \
                    - 2 * var_z * z, 2)
            print("about to solve")
            values = solve(eq, var_z)
            print("solved")
            for expr in values:
                num = float(expr)
                if -1 <= num <= 1: return num
            raise ValueError("no solution")
        return f

    def intersect(self, other):
        ''' return a function representing the intersections of self and other,
        where the function takes in a sweep plane value phi and outputs
        the intersection point '''
        theta0 = self.focus.get_theta()
        phi0 = self.focus.get_phi()
        theta1 = other.focus.get_theta()
        phi1 = other.focus.get_phi()
        def f(phi):
            v_theta = Symbol('theta', real = True)
            v_phi = Symbol('phi', real = True)
            locus0 = Eq(cos(v_theta - theta0) * cos(v_phi - phi0) - cos(v_phi - phi),0)
            locus1 = Eq(cos(v_theta - theta1) * cos(v_phi = phi1) - cos(v_phi - phi),0)
            values = solve([locus0, locus1], v_theta, v_phi)
            print(values)
            return values
        return f



if __name__ == '__main__':
    point = Point(0, 0.5**0.5)
    cap = Locus(point)
    f = cap.get_locus(0)
    for ind in range(6):
        print(f(pi/3 * ind))
