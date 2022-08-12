from gcmpy.tools.joint_excess_joint_degree_matrices import (
    JointExcessJointDegreeMatrices,
)


class JointExcessFromEjk:
    """
    Extracts the joint excess degree distributions from a
    list of joint excess joint degree distributions. This is
    done by summing the rows of the ejk matrices.
    """

    @staticmethod
    def get_excess_joint_distributions(
        ejks: JointExcessJointDegreeMatrices,
    ) -> dict:

        if len(ejks._ejks) != len(ejks.excess_degree_keys):
            raise ("Error in JointExcessDDFromEjk: incorrect ejk or keys provided")

        qks: dict = {}
        for key in ejks._ejks:
            ejk = ejks._ejks[key]
            keys = ejks._excess_degree_keys[key]
            q = {}
            for left_key in keys:
                for right_key in keys:
                    q[left_key] = q.get(left_key, 0.0) + ejk[left_key + right_key]
            qks[key] = q
        return qks
