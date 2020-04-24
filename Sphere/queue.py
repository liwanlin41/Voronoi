import heapq

class PriorityQueue:
    ''' create a priority queue with unique objects'''
    def __init__(self):
        self.queue_set = set()
        self.queue_heap = []

    def push(self, elt):
        if not elt in self.queue_set:
            heapq.heappush(self.queue_heap, elt)
            self.queue_set.add(elt)

    def pop(self):
        elt = heapq.heappop(self.queue_heap)
        self.queue_set.remove(elt)
        return elt

    def __len__(self):
        return len(self.queue_heap)
