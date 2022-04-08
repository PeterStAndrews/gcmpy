import random
from collections import Counter
from abc import ABC, abstractmethod


class JointDegree(ABC):

    """
    Abstract class whose purpose is to generate a joint degree
    distribution (jdd) dict through a variety of methods. A jdd
    can then be sampled at a later point during the GCM algorithm
    to create a joint degree sequence (jds).

    The abstract method `create_jdd' should be defined in a subclass
    and is the primary method that will be called to create an object
    to be sampled.
    """

    _type: str = ""

    def __init__(self):
        self._jdd: dict = None
        self._motif_sizes: list = None

    def __new__(cls, *args, **kwargs):
        if cls is JointDegree:
            raise TypeError(
                "The JointDegree class is abstract and may not be instantiated"
            )
        return object.__new__(cls)

    @abstractmethod
    def create_jdd(self) -> None:
        raise NotImplementedError(
            "Error attempting to call virtual method on JointDegree: create_joint_degree"
        )

    def handshaking_lemma(self, jds: list) -> list:
        """
        Samples a joint degree sequence from the joint degree distribution.
        Then ensures the handshaking lemma is satisfied by adding another
        subgraph to the network.
        """
        # ensure that the sum of the numbers is divisible by the motif size
        ntops = list(map(sum, zip(*jds)))
        for i, ntop in enumerate(ntops):
            if ntop % self._motif_sizes[i] != 0:
                # if not, round up to add one more motif to the network
                for j in range(self._motif_sizes[i] - ntop % self._motif_sizes[i]):
                    j = random.randrange(0, len(jds))
                    t = list(jds[j])
                    t[i] += 1
                    jds[j] = t
        return jds

    def sample_jds_from_jdd(self, N: int) -> list:
        """
        Chooses a list of `N' keys from `self._jdd' with replacement
        according to their weighting in the joint degree distribution.
        Checks the handshaking lemma is satisfied.
        :param N int: number of keys to choose
        :returns list: list of key choices
        """
        keys = list(self._jdd.keys())
        weights = list(self._jdd.values())
        jds = random.choices(population=keys, weights=weights, k=N)
        return self.handshaking_lemma(jds)

    def normalise_jdd(self) -> None:
        summation: float = sum(self._jdd.values())
        for key in self._jdd:
            self._jdd[key] /= summation

    def convert_jds_to_jdd(self, jds: list) -> None:
        n_samples: int = len(jds)
        self._jdd = {}
        d = Counter(jds)
        for k, v in d.items():
            self._jdd[k] = v / n_samples

    @property
    def jdd(self) -> dict:
        return self._jdd

    @jdd.setter
    def jdd(self, value: dict) -> None:
        self._jdd = value

    @property
    def motif_sizes(self) -> list:
        return self._motif_sizes

    @motif_sizes.setter
    def motif_sizes(self, value: list) -> None:
        self._motif_sizes = value
