class Edge:
    ''' represent a Voronoi edge; mutable endpoints '''

    def __init__(self, site1, site2):
        ''' create a new Edge object representing the Voronoi edge
        between site1 and site2 '''
        self.site1 = site1
        self.site2 = site2

    def get_sites(self):
        ''' return list of sites associated with this Edge '''
        return [self.site1, self.site2]

    def __eq__(self, other):
        return self.__class__ == other.__class__ and set(self.get_sites()) == set(other.get_sites())

    def __hash__(self):
        return self.site1.__hash__() + self.site2.__hash__()

    def __str__(self):
        return "%s, %s" %(self.site1, self.site2)

# testing
if __name__ == '__main__':
    edge_set = set()
    site1 = Point(0,0)
    site2 = Point(1,2)
    edge1 = Edge(site1, site2)
    edge2 = Edge(site2, site1)
    edge_set.add(edge1)
    edge_set.add(edge2)
    edge_set.add(Edge(site1, site1))
    for edge in edge_set:
        print(edge)


