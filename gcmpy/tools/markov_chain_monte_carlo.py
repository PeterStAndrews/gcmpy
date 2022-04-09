from functools import wraps


class MarkovChainMonteCarlo:

    """
    Base class for an MCMC experiment. Defines some decorators
    that yield the metadata of the Markov chain.
    """

    _proposal_count: int = 0
    _proposals_accepted: int = 0

    @classmethod
    def proposal_efficiency(self, func) -> callable:
        """
        Decorator to count the number of times a function is called
        and when it returns True. This is used to measure the efficiency
        of the Markov chain acceptance vs rejection.
        """

        @wraps(func)
        def wrapper(*args, **kwargs):
            self._proposal_count += 1
            result: bool = func(*args, **kwargs)
            if result:
                self._proposals_accepted += 1
            return result

        return wrapper
