def test_add_vertex_boundary_320():
    G = Graph(320)
    assert G.add_vertex() == 320
    assert G.add_vertex() == 321

