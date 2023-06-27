from sage.homology.homology_group import HomologyGroup

class MomentAngleComplex:
    def __init__(self, delta):
        self.delta = delta
        self.z_delta = None
        self.z_delta_topological = None

    def construct(self):
        self.z_delta = self.delta.moment_angle_complex()
        self.z_delta_topological = CellComplex(self.z_delta).topological_space()

    def vertices(self):
        return self.z_delta.vertices()

    def simplices(self):
        return self.z_delta.simplices()

    def homology(self):
        return self.z_delta_topological.homology()
