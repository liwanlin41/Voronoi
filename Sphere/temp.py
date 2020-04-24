#import numpy as np
from point import Point
from sympy import sin, cos, pi, Symbol, solve, Eq, evalf

class Locus:
    ''' represent locus of points equidistant from point and sweep plane '''

    def __init__(self, focus):
        ''' create a new locus focused at a given site '''
        self.focus = focus

    def get_focus(self):
        return self.focus

    def get_locus(self, phi):
        ''' get the equation of the locus when sweep plane is at phi 
        when evaluated at a site event f_phi = phi, this is not necessarily
        accurate and just returns phi '''
        def f(theta):
            # handle degenerate north/south pole cases
            d_theta = self.focus.get_theta() - theta if self.focus.get_phi() != pi/2 else 0
            f_phi = self.focus.get_phi()
            if f_phi == phi: return phi
            phi_val = Symbol('phi', real = True)
            print(d_theta)
            print(f_phi)
            print(phi)
            eq = Eq(cos(d_theta) * cos(phi_val - f_phi) - cos(phi_val - phi),0)
            print("about to solve")
            values = solve(eq, phi_val)
            print("solved")
            for expr in values:
                num = float(expr)
                if -pi/2 <= num <= pi/2: return num
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
    north_pole = Point(0, pi/4)
    cap = Locus(north_pole)
    f = cap.get_locus(0)
    for ind in range(6):
        print(f(pi/3 * ind))
