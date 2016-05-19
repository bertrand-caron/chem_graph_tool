from itertools import cycle, izip, product
from os.path import exists, join

try:
    from graph_tool.all import Graph, load_graph, graph_draw
    from graph_tool import topology
except:
    from sys import stderr
    stderr.write('Please install graph-tool')
    raise

from atb_helpers.pdb import is_pdb_atom_line, is_pdb_connect_line, pdb_fields

N = 10

PATTERNS = {
    'alcohol I': (
        ('J', 'H', 'H', 'C', 'O', 'H'),
        ((0, 3), (1, 3), (2, 3), (3, 4), (4, 5)),
    ),
    'alcohol II': (
        ('C', 'C', 'H', 'C', 'O', 'H'),
        ((0, 3), (1, 3), (2, 3), (3, 4), (4, 5)),
    ),
    'alcohol III': (
        ('C', 'C', 'C', 'C', 'O', 'H'),
        ((0, 3), (1, 3), (2, 3), (3, 4), (4, 5)),
    ),
    'benzene': (
        ('C', 'C', 'C', 'C', 'C', 'C', 'J', 'J', 'J', 'J', 'J', 'J'),
        ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (0, 6), (1, 7), (2, 8), (3, 9), (4, 10), (5, 11)),
    ),
}

ATOM_CLASSES = {
    'J': ('C', 'H'),
    'X': ('F', 'BR', 'CL', 'I'),
}

def atoms_for_class(atom_class):
    if atom_class in ATOM_CLASSES:
        return ATOM_CLASSES[atom_class]
    else:
        return (atom_class,)

def atom_classes(types):
    return sum([1 for a_type in types if a_type in ATOM_CLASSES.keys()])

def pattern_graph_for_pattern(pattern):
    vertices_types, edges = pattern

    graph = Graph(directed=False)

    vertex_types = graph.new_vertex_property("string")
    graph.vertex_properties['type'] = vertex_types

    vertices = []

    for (i, vertex_type) in enumerate(vertices_types):
        v = graph.add_vertex()
        vertex_types[v] = vertex_type
        vertices.append(v)

    for (i, j) in edges:
        graph.add_edge(vertices[i], vertices[j])

    draw_graph(
        graph,
    )

    return graph

def graphs_for_pattern_graph(pattern_graph):
    get_vertex_type = lambda v: pattern_graph.vp.type[v]

    type_permutations = product(*[atoms_for_class(get_vertex_type(v)) for v in pattern_graph.vertices()])

    graphs = []
    for type_permutation in type_permutations:
        new_graph = Graph(pattern_graph)
        for (v, vertex_type) in zip(new_graph.vertices(), type_permutation):
            new_graph.vp.type[v] = vertex_type

        graphs.append(new_graph)
    return graphs

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
            vertex_types[v] = fields[11].strip().upper()
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

def draw_graph(graph, fnme='graph'):
    try:
        vertex_text=graph.vertex_properties['type']
    except:
        vertex_text=graph.vertex_index


    if not '.png' in fnme:
        fnme += '.png'

    if exists(fnme):
        return

    graph_draw(
        graph,
        vertex_text=vertex_text,
        vertex_font_size=18,
        output_size=(200, 200),
        output=fnme,
)

pattern_graphs = [pattern_graph_for_pattern(pattern) for (moiety, pattern) in PATTERNS.items()]
interpreted_pattern_graphs = [(moiety, graphs_for_pattern_graph(pattern_graph)) for (moiety, pattern_graph) in zip(PATTERNS.keys(), pattern_graphs)]
for (moiety, graph_list) in interpreted_pattern_graphs:
    [
        draw_graph(
            graph,
            fnme=join('patterns', moiety.replace(' ', '_') + '_' + str(i)),
        )
        for (i, graph) in enumerate(graph_list)
    ]

def moieties_in_graph(super_graph):
    def match(moiety, graph_list):
        return any([
            topology.subgraph_isomorphism(
                pattern_graph,
                super_graph,
                vertex_label=(
                    pattern_graph.vertex_properties['type'],
                    super_graph.vertex_properties['type'],
                ),
            )
        for (i, pattern_graph) in enumerate(graph_list)
        ])

    return [
        moiety
        for (moiety, graph_list) in interpreted_pattern_graphs
        if match(moiety, graph_list)
    ]

from glob import glob

TEST_PDBS = glob('data/*.pdb')

if __name__ == '__main__':
    for test_pdb in TEST_PDBS:
        with open(test_pdb) as fh:
            molecule_graph = graph_from_pdb(fh.read())

        draw_graph(
            molecule_graph,
            fnme=test_pdb.replace('.pdb', '.png')
        )
        molecule_graph.save(test_pdb.replace('.pdb', '.gt'))

        print test_pdb
        print moieties_in_graph(molecule_graph)
