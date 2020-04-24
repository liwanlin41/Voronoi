from point import Point

class Edge:
    ''' represent a Voronoi edge on the sphere '''

    def __init__(self, site1, site2):
        ''' create a new Edge object representing the Voronoi edge between
        Point site1 and Point site2 '''
        # ensure site1 < site2
        self.site1, self.site2 = site1, site2 if site1 < site2 else site2, site1
        self.endpoints = set()

    def get_sites(self):
        ''' return list of sites associated with this Edge '''
        return [self.site1, self.site2]

    def __eq__(self, other):
        return self.__class__ == other.__class__ and self.get_sites() == other.get_sites()

    def __hash__(self):
        return self.site1.__hash__() + self.site2.__hash__()

    def __str__(self):
        return "%s, %s" %(self.site1, self.site2)

if __name__ == '__main__':
    site1 = Point(0,0)
    site2 = Point(0,0)
    assert site1 == site2
