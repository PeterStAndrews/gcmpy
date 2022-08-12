import networkx as nx


class MessagePassingMixin():

    def __init__(self, cover_type: str, G: nx.Graph):
        '''
        Mixing class to pull the motif cover labelling from a network model.
        The cover labelling depends on the motifs in the cover, whereas the
        graph label is constant.

        :param cover_type: str of motifs in cover
        :param G: networkx graph with edge labels.
        '''
        self._CoverType: str = cover_type
        self._G: nx.Graph = G

    def get_edge_cover_label(i: int, j: int) -> str:
        '''
        Interrogates the graph `G' for edge <i,j>'s cover label. Each label
        will depend on the cover that is being used.

        :param i: vertex id
        :param j: vertex id

        :returns string: the cover label
        '''
        pass

    def get_motif_topology(label: str) -> str:
        '''
        Parses the cover label to get the topology of the edge.

        :param label: the cover label
        :returns string: topology
        '''
        pass

    def get_motif_ID(label: str) -> str:
        '''
        Parses the cover label to return the unique ID of
        the motif.

        :param str: cover label
        :returns str: unique motif ID
        '''
        pass

    def get_vertices_in_motif(label: str) -> list:
        '''
        Parses the cover label to return the vertices in
        the motif as a list of integers.

        :param str: cover label
        :returns list: vertex IDs in motif.
        '''
        pass
