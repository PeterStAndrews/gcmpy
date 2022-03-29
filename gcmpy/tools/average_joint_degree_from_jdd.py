class AverageJointDegreeFromJDD:
    @staticmethod
    def get_average_joint_degrees(jdd: dict) -> list:
        joint_degrees = list(jdd.keys())
        num_topologies = len(joint_degrees[0])
        average_degrees = [0.0] * num_topologies
        for joint_degree in joint_degrees:
            for index in range(num_topologies):
                average_degrees[index] += joint_degree[index] * jdd[joint_degree]
        return average_degrees
