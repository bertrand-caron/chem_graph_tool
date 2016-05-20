from itertools import cycle, izip, product
from os.path import exists, join
from re import search

try:
    from graph_tool.all import Graph, load_graph, graph_draw
    from graph_tool import topology
except:
    from sys import stderr
    stderr.write('Please install graph-tool')
    raise

from atb_helpers.pdb import is_pdb_atom_line, is_pdb_connect_line, pdb_fields

DRAW_PATTERN_GRAPHS = True

PATTERNS = {
    'alcohol I': (
        ('J', 'H', 'H', 'C4', 'O2', 'H',),
        ((0, 3), (1, 3), (2, 3), (3, 4), (4, 5),),
    ),
    'alcohol II': (
        ('C', 'C', 'H', 'C4', 'O2', 'H',),
        ((0, 3), (1, 3), (2, 3), (3, 4), (4, 5),),
    ),
    'alcohol III': (
        ('C', 'C', 'C', 'C4', 'O2', 'H',),
        ((0, 3), (1, 3), (2, 3), (3, 4), (4, 5),),
    ),
    'benzene': (
        ('C3', 'C3', 'C3', 'C3', 'C3', 'C3',),
        ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),)
    ),
    'cyclohexane': (
        ('C4', 'C4', 'C4', 'C4', 'C4', 'C4',),
        ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),),
    ),
    'amine I': (
        ('J{3,4}', 'N3', 'H', 'H'),
        ((0,1), (1, 2), (1, 3),),
    ),
    'amine II': (
        ('C{3,4}', 'N3', 'C{3,4}', 'H'),
        ((0,1), (1, 2), (1, 3),),
    ),
    'amine III': (
        ('C{3,4}', 'N3', 'C{3,4}', 'C{3,4}'),
        ((0,1), (1, 2), (1, 3),),
    ),
}

MONOVALENT = (1,)
HALOGEN = MONOVALENT
CHALCOGEN = (1, 2, 3),

DEFAULT_VALENCES = {
    'H': MONOVALENT,
    'F': HALOGEN,
    'BR': HALOGEN,
    'CL': HALOGEN,
    'I': HALOGEN,
    'C': (2, 3, 4,),
    'O': CHALCOGEN,
    'N': (1, 2, 3, 4,),
    'S': CHALCOGEN,
}

ATOM_CLASSES = {
    'J': ('C', 'H',),
    'X': ('F', 'BR', 'CL', 'I',),
}

# Construct DEFAULT_VALENCES of atom_classes automatically
for (atom_class_name, atom_class_atoms) in ATOM_CLASSES.items():
    DEFAULT_VALENCES[atom_class_name] = tuple(
        set(
            reduce(
                lambda acc, e: acc + e,
                [DEFAULT_VALENCES[atom_type] for atom_type in atom_class_atoms],
                (),
            )
        )
    )

def parse_atom_class(atom_class):
    m =  search('^([A-Z]+){?([0-9]?),?([0-9]?)}?', atom_class)
    assert m
    return m

def atoms_for_class(atom_class):
    m = parse_atom_class(atom_class)

    type_class = m.group(1)

    if type_class in ATOM_CLASSES:
        return ATOM_CLASSES[type_class]
    else:
        return (type_class,)

def valences_for_class(atom_class):
    m = parse_atom_class(atom_class)

    if not m.group(2) or not m.group(3):
        start = [int(group) for group in (m.group(2), m.group(3)) if group]
        if len(start) == 0:
            return DEFAULT_VALENCES[atom_class]
        else:
            start = start[0]
        end = start + 1
    else:
        start, end = int(m.group(2)), int(m.group(3)) + 1

    return tuple(range(start, end))

def types_and_valences_for_class(atom_class):
    return map(
        lambda (atom_class, valence): '{0}{1}'.format(atom_class, valence),
        product(
            atoms_for_class(atom_class),
            valences_for_class(atom_class),
        ),
    )

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

    return graph

def graphs_for_pattern_graph(pattern_graph):
    get_vertex_type = lambda v: pattern_graph.vp.type[v]

    type_permutations = product(*[types_and_valences_for_class(get_vertex_type(v)) for v in pattern_graph.vertices()])

    graphs = []
    for type_permutation in type_permutations:
        new_graph = Graph(pattern_graph)
        for (v, vertex_type) in zip(new_graph.vertices(), type_permutation):
            new_graph.vp.type[v] = vertex_type

        graphs.append(new_graph)
    return graphs

def write_dummy_graph(n=10, cyclic=True):
    graph_file = 'data/dummy.gt'

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

    def get_connect_list():
        connect_ids_list = [
            [int(str_id) for str_id in line.split()[1:]]
            for line in pdb_str.splitlines()
            if is_pdb_connect_line(line)
        ]

        return reduce(
            lambda acc, e: acc + e,
            [
                [
                    (i, j)
                    for (i, j) in map(
                        lambda (i, j): (i - 1, j - 1),
                        izip(cycle(connect_ids[0:1]), connect_ids[1:]),
                    )
                    if i < j
                ]
                for connect_ids in connect_ids_list
            ],
            [],
        )

    connects = get_connect_list()

    def get_valence(atom_id):
        return sum([1 for connect in connects if atom_id in connect])

    atom_lines = [
        line for line in pdb_str.splitlines()
        if is_pdb_atom_line(line)
    ]

    for (atom_id, line) in enumerate(atom_lines):
        fields = pdb_fields(line)
        v = g.add_vertex()
        vertex_types[v] = '{0}{1}'.format(
            fields[11].strip().upper(),
            get_valence(atom_id),
        )
        vertices.append(v)

    for (i, j) in connects:
        g.add_edge(vertices[i], vertices[j])

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

if DRAW_PATTERN_GRAPHS:
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

def test_atom_class_parsing():
    for atom_class in ('C4', 'C{4,7}', 'C', 'J{3,4}'):
        print list(types_and_valences_for_class(atom_class))

def moieties_in_pdb_file(pdb_file, should_draw_graph=True, should_dump_graph=False):
    with open(pdb_file) as fh:
        molecule_graph = graph_from_pdb(fh.read())

    if should_draw_graph:
        draw_graph(
            molecule_graph,
            fnme=pdb_file.replace('.pdb', '.png')
        )

    if should_dump_graph:
        molecule_graph.save(pdb_file.replace('.pdb', '.gt'))

    return moieties_in_graph(molecule_graph)

if __name__ == '__main__':
    if False:
        test_atom_class_parsing()

    for test_pdb in TEST_PDBS:
        print test_pdb
        print moieties_in_pdb_file(test_pdb)
