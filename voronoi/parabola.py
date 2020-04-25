from point import Point
import warnings
warnings.filterwarnings("error")

class Parabola:
    ''' represent an entire parabolic locus'''

    def __init__(self, focus):
        ''' create a new parabola with a given focus '''
        self.focus = focus

    def NEG_INF():
        ''' return pseudo parabola object all the way to the left '''
        return Parabola(Point.NEG_INF())

    def INF():
        ''' return pseudo parabola object all the way to the right '''
        return Parabola(Point.INF())

    def get_focus(self):
        return self.focus

    # don't evaluate the locus at a site event, it won't be accurate - just returns 0
    def get_locus(self, directrix):
        ''' get the equation of the parabola at a directrix line '''
        def f(x):
#            print("Evaluate locus at %f, directrix = %f" %(x,directrix))
            fx = self.focus.get_x()
            fy = self.focus.get_y()
            if fy == directrix:
                return 0
            # locus is given by y = 1/(2fy-2d) (x^2-2x fx+fx^2+fy^2-d^2)
            return ((x-fx)**2 + fy**2 - directrix**2)/(2 * (fy - directrix))
        return f

    def intersect(self, other):
        ''' return a pair of functions representing the intersections
        of self and other,
        where the function takes in a directrix line y=d and outputs 
        the intersection point when the sweep line is at this directrix '''
        fx0 = self.focus.get_x()
        fy0 = self.focus.get_y()
        fx1 = other.focus.get_x()
        fy1 = other.focus.get_y()

        # functions for left and right
        def f_left(d):
            # represent as ax^2 + bx + c = 0;
            a = fy1 - fy0
            b = 2 * (fx1 * (fy0-d) - fx0 * (fy1-d))
            c = (fx0**2 + fy0**2 - d**2) * (fy1-d) - (fx1**2 + fy1**2 - d**2) * (fy0 - d)
            discr_square = b**2 - 4 * a * c
            if discr_square < 0:
#                print("complex discriminant")
#                assert -1e-11 < discr_square
                discr_square = 0
#            try:
            discr = discr_square**0.5
#            except RuntimeWarning: # account for rounding error
#                assert -1e-11 < discr_square < 1e-11
#                except:
#                    print("Assertion error")
#                    print(discr_square)
#                    input()
#                discr = 0
            # account for plus or minus error to ensure the left point is returned
            diff = discr if a > 0 else -1 * discr
            intersect_x = (fx0 + fx1)/2 if (a == 0) else -0.5 * (b + diff) / a
                # compute y coord from intersect_x
            intersect_y = self.get_locus(d)(intersect_x)
            return Point(intersect_x, intersect_y)
        def f_right(d):
            a = fy1 - fy0
            b = 2 * (fx1 * (fy0-d) - fx0 * (fy1-d))
            c = (fx0**2 + fy0**2 - d**2) * (fy1-d) - (fx1**2 + fy1**2 - d**2) * (fy0 - d)
            discr_square = b**2 - 4 * a * c
            if discr_square < 0:
#                print("negative discriminant")
#                try:
#                assert -1e-11 < discr_square < 1e-11
#                except AssertionError:
#                    print("Assertion error")
#                    input()
#                    print(discr_square)
                discr_square = 0
#            try:
            discr = discr_square**0.5
#            except RuntimeWarning: # account for rounding error
#                assert -1e-11 < discr_square < 1e-11
#                discr = 0
            diff = discr if a > 0 else -1 * discr
            intersect_x = (fx0 + fx1)/2 if (a == 0) else -0.5 * (b - diff) / a
            intersect_y = self.get_locus(d)(intersect_x)
            return Point(intersect_x, intersect_y)
        return (f_left, f_right)

    def ordered_intersect(left_parabola, right_parabola):
        ''' return only one of the functions from intersect '''
        # check degenerate cases
        if left_parabola.get_focus() == Point.NEG_INF():
            return lambda d: Point.NEG_INF()
        if right_parabola.get_focus() == Point.INF():
            return lambda d: Point.INF()
        # do something normal
        intersections = left_parabola.intersect(right_parabola)
        # if a left parabola has focus above the right one, it will be wider
        # and there are two intersection points to the right; we want
        # the closer one
        ind = 0 if left_parabola.get_focus().get_y() > right_parabola.get_focus().get_y() else 1
        return intersections[ind]


# for testing
if __name__ == '__main__':
    parabola1 = Parabola(Point(0,0))
    parabola2 = Parabola(Point(2,0))
    parabola3 = Parabola(Point(2,1))
    print("Locus of points on parabola 1 with focus (0,0) and directrix -1")
    for i in range(10):
        print(parabola1.get_locus(-1)(i))
    print("Right intersection of focus (0,0) and (2,0) with moving directrix")
    for i in range(-1, -10, -1):
        print(parabola1.intersect(parabola2)[1](i))
    print("Left intersection of focus (0,0) and (2,1)")
    for i in range(-1, -10, -1):
        print(parabola1.intersect(parabola3)[0](i))
    print("Right intersection of focus (2,0) and (2,1)")
    for i in range(-1, -10, -1):
        print(parabola2.intersect(parabola3)[1](i))
    print("Left is (0,0) and right is (2,1)")
    for i in range(-1, -10, -1):
        print(Parabola.ordered_intersect(parabola1, parabola3)(i))
