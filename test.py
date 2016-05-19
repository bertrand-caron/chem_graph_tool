from itertools import cycle, izip

try:
    from graph_tool.all import Graph, load_graph, graph_draw
    from graph_tool import topology
except:
    from sys import stderr
    stderr.write('Please install graph-tool')
    raise

from atb_helpers.pdb import is_pdb_atom_line, is_pdb_connect_line, pdb_fields

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

def graph_from_pdb(pdb_str):
    g = Graph(directed=False)

    vertex_types = g.new_vertex_property("string")
    g.vertex_properties['type'] = vertex_types

    vertices = []

    for line in pdb_str.splitlines():
        if is_pdb_atom_line(line):
            fields = pdb_fields(line)
            v = g.add_vertex()
            vertex_types[v] = fields[11]
            vertices.append(v)
        elif is_pdb_connect_line(line):
            connect_ids = [int(str_id) for str_id in line.split()[1:]]
            for (i, j) in map(lambda (i, j): (i - 1, j - 1), izip(cycle(connect_ids[0:1]), connect_ids[1:])):
                if i > j:
                    continue
                g.add_edge(vertices[i], vertices[j])
        else:
            pass

    return g

def draw_graph(graph, vertex_text=None):
    graph_draw(
        graph,
        vertex_text=(graph.vertex_index if vertex_text is None else vertex_text),
        vertex_font_size=18,
        output_size=(200, 200),
        output="graph.png",
)

TEST_PDB = 'data/test.pdb'

if __name__ == '__main__':
    with open(TEST_PDB) as fh:
        molecule_graph = graph_from_pdb(fh.read())

    draw_graph(
        molecule_graph,
        vertex_text=molecule_graph.vertex_properties['type'],
    )
    molecule_graph.save(TEST_PDB.replace('.pdb', '.gt'))
    exit()

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

