try:
    from graph_tool.all import Graph, load_graph
    from graph_tool import topology
except:
    from sys import stderr
    stderr.write('Please install graph-tool')
    raise

N = 10

def write_dummy_graph(n=N, cyclic=True):
    graph_file = 'my_graph.gt'

    g = Graph(directed=False)
    vertices = [g.add_vertex() for n in range(0, N)]
    vertices_pairs = zip(vertices, vertices[1:] + (vertices[0:1] if cyclic else []))
    [g.add_edge(v_1, v_2) for (v_1, v_2) in vertices_pairs]

    vertex_types = g.new_vertex_property("string")
    g.vertex_properties['type'] = vertex_types
    for (i, v) in enumerate(g.vertices()):
        vertex_types[v] = 'C' + str(i)

    g.save(graph_file)
    return (graph_file, vertex_types)

if __name__ == '__main__':
    graph_file, vertex_types = write_dummy_graph()
    g = load_graph(graph_file)

    g2 = Graph(g)


    print topology.subgraph_isomorphism(
        g,
        g2,
        vertex_label=(
            g.vertex_properties['type'],
            g2.vertex_properties['type'],
        ),
    )
