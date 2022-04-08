class JointExcessJointDegreeKeysView:

    """
    Utility class to access correct view of the joint
    excess joint degree keys.
    :param keys: ordered list of keys [u0, u1, v0, v1]
    """

    def __init__(self, keys: list):
        self._keys: list = keys

    def get_u0u1(self):
        return self._keys[0] + self._keys[1]

    def get_u1u0(self):
        return self._keys[1] + self._keys[0]

    def get_v0v1(self):
        return self._keys[2] + self._keys[3]

    def get_v1v0(self):
        return self._keys[3] + self._keys[2]

    def get_u0v1(self):
        return self._keys[0] + self._keys[3]

    def get_v0u1(self):
        return self._keys[2] + self._keys[1]
