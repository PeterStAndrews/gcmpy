import ast
import networkx as nx


class MessagePassingMixin():

    def __init__(self, cover_type: str, G: nx.Graph):
        '''
        Mixin class to pull the motif cover labelling from a network model.
        The cover labelling depends on the motifs in the cover, whereas the
        graph label is constant.
        
        Cover label is assumed to be of the form 
        
            f"{<key>}-{[<vertices>]}-{[<edges>]}-{UID}".
        
        For instance for a triangle we have 
        
            f"{3}-{[n1,n2,n3]}-{[(n1,n2),(n1,n3),(n2,n3)]}-{ID}"

        The key is a way of identifying the topology of the motif and should
        represent each site of each motif uniquely.

        :param cover_type: str of motifs in cover
        :param G: networkx graph with edge labels.
        '''
        self._CoverType: str = cover_type
        self._G: nx.Graph = G

    def get_edge_cover_label(self, i: int, j: int) -> str:
        '''
        Interrogates the graph `G' for edge <i,j>'s cover label. Each label
        will depend on the cover that is being used.

        :param i: vertex id
        :param j: vertex id

        :returns string: the cover label
        '''
        return self._G.edges[i,j]['CoverLabel']

    def get_motif_topology(self, label: str) -> int:
        '''
        Parses the cover label to get the topology of the edge.

        :param label: the cover label
        :returns int: topology
        '''
        return int(label.split('-')[0])

    def get_motif_ID(self, label: str) -> int:
        '''
        Parses the cover label to return the unique ID of
        the motif.

        :param str: cover label
        :returns str: unique motif ID
        '''
        return int(label.split('-')[-1])

    def get_vertices_in_motif(self, label: str) -> list:
        '''
        Parses the cover label to return the vertices in
        the motif as a list of integers.

        :param str: cover label
        :returns list: vertex IDs in motif.
        '''
        return ast.literal_eval(label.split('-')[1])
    
    def get_edges_in_motif(self, label: str) -> list:
        '''
        Parses the cover label to return the edges in
        the motif as a list of tuples of integers.

        :param str: cover label
        :returns list: edges in the motif.
        '''
        return ast.literal_eval(label.split('-')[2])