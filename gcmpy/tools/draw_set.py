import random


class DrawSet(object):

    def __init__(self):
        self._edge_hashmap: dict = {}
        self._edges: list = []

    def __contains__(self, e: tuple) -> bool:
        return e in self._edge_hashmap

    def __iter__(self):
        return iter(self._edges)

    def __len__(self):
        return len(self._edges)

    def add(self, e: tuple) -> None:
        if e in self._edge_hashmap:
            return
        self._edges.append(e)
        self._edge_hashmap[e] = len(self._edges)-1

    def remove(self, e: tuple) -> None:
        position = self._edge_hashmap.pop(e)
        last_item = self._edges.pop()
        if position != len(self._edges):
            self._edges[position] = last_item
            self._edge_hashmap[last_item] = position

    def draw(self) -> tuple:
        return random.choice(self._edges)
