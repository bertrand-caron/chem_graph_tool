from itertools import product, combinations
from os.path import exists, join, dirname, abspath
from re import search, sub
from functools import reduce
from typing import List, Tuple, Any, Optional
from glob import glob

from chem_graph_tool.pdb import Graph, load_graph, graph_draw, topology, graph_from_pdb, type_identifier_for, sfdp_layout, PropertyMap

DRAW_PATTERN_GRAPHS = True

DISABLE_PATTERNS = False

Bond = Tuple[int, int]

Atom_Pattern = str

Moiety = str

Element = str

Graph_Pattern = Tuple[List[Atom_Pattern], List[Bond]]

ROOT_DIR = dirname(abspath(__file__))

def OR(*x: List[str]) -> str:
 return '|'.join(x)

H = 'H1'
R_NO_H = 'C4'
R = ALKYL = OR(R_NO_H, H)
R2 = OR(R, 'C3')

PHENYL_CORE = (
    ['C3', 'C3', 'C3', 'C3', 'C3', 'C3'],
    [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)],
)

MONO_HALOGENO = lambda X: (
        ['C4', X, R, R, R2],
        [(0, 1), (0, 2), (0, 3), (0, 4)],
    )

DI_HALOGENO = lambda X: (
        ['C4', X, X, R, R2],
        [(0, 1), (0, 2), (0, 3), (0, 4)],
    )

TRI_HALOGENO = lambda X: (
        ['C4', X, X, X, R2],
        [(0, 1), (0, 2), (0, 3), (0, 4)],
    )

TETRA_HALOGENO = lambda X: (
        ['C4', X, X, X, X],
        [(0, 1), (0, 2), (0, 3), (0, 4)],
    )

PATTERNS = ({
    'alkane': (
        ['J', R, R, 'C4', 'C4', H, H, H],
        [(0, 3), (1, 3), (2, 3), (3, 4), (4, 5), (4, 6), (4, 7)],
    ),
    'alkene': (
        ['J', R, 'C3', 'C3', R, R],
        [(0, 2), (1, 2), (2, 3), (3, 4), (3, 5)],
    ),
    'alkyne': (
        [R, 'C2', 'C2', R],
        [(0, 1), (1, 2), (2, 3)],
    ),
    'alcohol I': (
        ['J', H, H, 'C4', 'O2', H],
        [(0, 3), (1, 3), (2, 3), (3, 4), (4, 5)],
    ),
    'alcohol II': (
        ['C', 'C', H, 'C4', 'O2', H],
        [(0, 3), (1, 3), (2, 3), (3, 4), (4, 5)],
    ),
    'alcohol III': (
        ['C', 'C', 'C', 'C4', 'O2', H],
        [(0, 3), (1, 3), (2, 3), (3, 4), (4, 5)],
    ),
    'thiol': (
        [R, 'S2', H],
        [(0, 1), (1, 2)],
    ),
    'thioether': (
        ['C4', 'S2', 'C4'],
        [(0, 1), (1, 2)],
    ),
    'thiophene': (
        ['S', 'C3', 'C3', 'C3', 'C3'],
        [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)],
    ),
    'disulfide': (
        [R, 'S2', 'S2', R],
        [(0, 1), (1, 2), (2, 3)],
    ),
    'phenyl': (
        PHENYL_CORE[0] + [R, R, R, R, R, 'H1|C3|C4'],
        PHENYL_CORE[1] + [(0, 6), (1, 7), (2, 8), (3, 9), (4, 10), (5, 11)],
    ),
    'phenol': (
        PHENYL_CORE[0] + ['O2', H],
        PHENYL_CORE[1] + [(5, 6), (6, 7)],
    ),
    'phenoxy': (
        PHENYL_CORE[0] + ['O2', R_NO_H],
        PHENYL_CORE[1] + [(5, 6), (6, 7)],
    ),
    'aniline': (
        PHENYL_CORE[0] + ['N3', H, H],
        PHENYL_CORE[1] + [(5, 6), (6, 7), (6, 8)],
    ),
    'fluorophenyl': (
        PHENYL_CORE[0] + ['F'],
        PHENYL_CORE[1] + [(5, 6)],
    ),
    'chlorophenyl': (
        PHENYL_CORE[0] + ['CL'],
        PHENYL_CORE[1] + [(5, 6)],
    ),
    'bromophenyl': (
        PHENYL_CORE[0] + ['BR'],
        PHENYL_CORE[1] + [(5, 6)],
    ),
    'iodophenyl': (
        PHENYL_CORE[0] + ['I'],
        PHENYL_CORE[1] + [(5, 6)],
    ),
    'pyridine': (
        ['N2', 'C3', 'C3', 'C3', 'C3', 'C3'],
        [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)],
    ),
    'cyclohexane': (
        ['C{3,4}', 'C{3,4}', 'C4', 'C4', 'C4', 'C4'],
        [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)],
    ),
    'amine I': (
        [R, 'N3', H, H],
        [(0, 1), (1, 2), (1, 3)],
    ),
    'amine II': (
        ['C4', 'N3', 'C4', H],
        [(0, 1), (1, 2), (1, 3)],
    ),
    'amine III': (
        ['C4', 'N3', 'C4', 'C4'],
        [(0, 1), (1, 2), (1, 3)],
    ),
    'ammonium': (
        ['C4', 'N4', 'C4', 'C4', 'C4'],
        [(0, 1), (1, 2), (1, 3), (1, 4)],
    ),
    'enamine': (
        ['C3', 'C3', 'N3', 'C4', 'C4', R],
        [(0, 1), (1, 2), (2, 3), (2, 4), (1, 5)],
    ),
    'ester': (
        [R2, 'C3', 'O1', 'O2', 'C{3,4}'],
        [(0, 1), (1, 2), (1, 3), (3, 4)],
    ),
    'carboxylic acid': (
        [R2, 'C3', 'O1', 'O2', H],
        [(0, 1), (1, 2), (1, 3), (3, 4)],
    ),
    'carboxylate': (
        [R2, 'C3', 'O1', 'O1'],
        [(0, 1), (1, 2), (1, 3)],
    ),
    'ether': (
        ['C4', 'O2', 'C4'],
        [(0, 1), (1, 2)],
    ),
    'peroxide': (
        [R, 'O2', 'O2', R],
        [(0, 1), (1, 2), (2, 3)],
    ),
    'amide': (
        [R2, 'C3', 'O1', 'N3', R, R],
        [(0, 1), (1, 2), (1, 3), (3, 4), (3, 5)],
    ),
    'hemiketal': (
        [R, R, 'C4', 'O2', H, 'O2', R],
        [(0, 2), (1, 2), (2, 3), (3, 4), (2, 5), (5, 6)],
    ),
    'ketone': (
        ['C{3,4}', 'C3', 'O1', 'C{3,4}'],
        [(0, 1), (1, 2), (1, 3)],
    ),
    'aldehyde': (
        [H, 'C3', 'O1', 'C{3,4}'],
        [(0, 1), (1, 2), (1, 3)],
    ),
    'nitrile': (
        [R2, 'C2', 'N1'],
        [(0, 1), (1, 2)],
    ),
    'nitro': (
        [R, 'N3', 'O1', 'O1'],
        [(0, 1), (1, 2), (1, 3)],
    ),
    'nitrophenyl': (
        PHENYL_CORE[0] + ['N3', 'O1', 'O1'],
        PHENYL_CORE[1] + [(5, 6), (6, 7), (6, 8)],
    ),
    'monochloro': MONO_HALOGENO('CL'),
    'dichloro': DI_HALOGENO('CL'),
    'trichloro': TRI_HALOGENO('CL'),
    'tetrachloro': TETRA_HALOGENO('CL'),
    'monofluoro': MONO_HALOGENO('F'),
    'difluoro': DI_HALOGENO('F'),
    'trifluoro': TRI_HALOGENO('F'),
    'tetrafluoro': TETRA_HALOGENO('F'),
    'monobromo': MONO_HALOGENO('BR'),
    'dibromo': DI_HALOGENO('BR'),
    'tribromo': TRI_HALOGENO('BR'),
    'tetrabromo': TETRA_HALOGENO('BR'),
    'monoiodo': MONO_HALOGENO('I'),
    'diiodo': DI_HALOGENO('I'),
    'triiodo': TRI_HALOGENO('I'),
    'tetraiodo': TETRA_HALOGENO('I'),
    'sulfonyl': (
        [R, 'S4', 'O1', 'O1', R],
        [(0, 1), (1, 2), (1, 3), (1, 4)],
    ),
    'sulfinyl': (
        [R, 'S3', 'O1', R],
        [(0, 1), (1, 2), (1, 3)],
    ),
    'nitrate': (
        [R, 'O2', 'N', 'O1', 'O1'],
        [(0, 1), (1, 2), (2, 3), (2, 4)],
    ),
    'hydrazine': (
        [R, R, 'N3', 'N3', R, R],
        [(0, 2), (1, 2), (2, 3), (3, 4), (3, 5)],
    ),
    'cyclopropane': (
        ['C4', 'C4', 'C4'],
        [(0, 1), (1, 2), (2, 0)],
    ),
    'cyclobutane': (
        ['C4', 'C4', 'C4', 'C4'],
        [(0, 1), (1, 2), (2, 3), (3, 0)],
    ),
    'urea': (
        [R2, R2, 'N3', 'C3', 'O1', 'N3', R2, R2],
        [(0, 2), (1, 2), (2, 3), (3, 4), (3, 5), (5, 6), (5, 7)],
    ),
    'phosphate': (
        [R2, 'O2', 'P4', 'O1', 'O{1,2}', 'O{1,2}'],
        [(0, 1), (1, 2), (2, 3), (2, 4), (2, 5)],
    ),
    'imidazole ': (
        ['N2', 'C3', 'C3', 'N3', 'C3'],
        [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)],
    ),
    'aromatic amine II': (
        ['N{2,3}', 'C3', 'C3', 'C3', 'C3'],
        [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)],
    ),
    'silane': (
        ['SI', R, R, R, R],
        [(0, 1), (0, 2), (0, 3), (0, 4)],
    ),
} if not DISABLE_PATTERNS else {})

MONOVALENT = [1]
HALOGEN = MONOVALENT
CHALCOGEN = [1, 2, 3]

DEFAULT_VALENCES = {
    'H': MONOVALENT,
    'F': HALOGEN,
    'BR': HALOGEN,
    'CL': HALOGEN,
    'I': HALOGEN,
    'C': [2, 3, 4],
    'O': CHALCOGEN,
    'N': [1, 2, 3, 4],
    'S': CHALCOGEN,
    'SI': [2, 3, 4],
}

ATOM_CLASSES = {
    'J': ['C', 'H'],
    'X': ['F', 'BR', 'CL', 'I'],
}

def parse_atom_class(atom_pattern: Atom_Pattern) -> Any:
    m =  search('^([A-Z]+){?([0-9]?),?([0-9]?)}?', atom_pattern)
    assert m
    return m

def types_and_valences_for_class(atom_pattern: Atom_Pattern) -> List[str]:
    def atoms_for_class(atom_pattern: Atom_Pattern) -> List[Element]:
        m = parse_atom_class(atom_pattern)

        type_class = m.group(1)

        if type_class in ATOM_CLASSES:
            return ATOM_CLASSES[type_class]
        else:
            return [type_class]

    def valences_for_class(atom_pattern: Atom_Pattern) -> List[int]:
        m = parse_atom_class(atom_pattern)

        if not m.group(2) or not m.group(3):
            start = [int(group) for group in (m.group(2), m.group(3)) if group]
            if len(start) == 0:
                if atom_pattern in DEFAULT_VALENCES:
                    return DEFAULT_VALENCES[atom_pattern]
                else:
                    raise Exception(atom_pattern)
            else:
                start = start[0]
            end = start + 1
        else:
            start, end = int(m.group(2)), int(m.group(3)) + 1

        return list(range(start, end))

    if '|' in atom_pattern:
        type_valence_list = atom_pattern.split('|')
    else:
        type_valence_list = list(map(
            lambda atom_class_valence: type_identifier_for(atom_class_valence[0], atom_class_valence[1]),
             reduce(
                lambda acc, e: acc + e,
                [
                    list(
                        product(
                            (atom, ),
                            valences_for_class(
                                sub(
                                    '^[A-Z]+',
                                    atom,
                                    atom_pattern,
                                ),
                            ),
                        ),
                    )
                    for atom in atoms_for_class(atom_pattern)
                ],
                [],
            ),
        ))
    return type_valence_list

def atom_classes(types: List[str]) -> int:
    return sum([1 for a_type in types if a_type in list(ATOM_CLASSES.keys())])

def pattern_graph_for_graph_pattern(graph_pattern: Graph_Pattern) -> Graph:
    vertices_types, edges = graph_pattern

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

MAX_NUMBER_PERMUTATIONS = 100

def graphs_for_pattern_graph(pattern_graph: Graph, pattern_identifier: str = '') -> List[Graph]:
    get_vertex_type = lambda v: pattern_graph.vp.type[v]

    type_permutations = list(
        product(
            *[
                types_and_valences_for_class(get_vertex_type(v))
                for v in pattern_graph.vertices()
            ],
        ),
    )

    assert len(type_permutations) <= MAX_NUMBER_PERMUTATIONS, '''Error: Unreasonably large number ({0}) of graphs for pattern identifier '{1}'. Aborting ...'''.format(
        len(type_permutations),
        pattern_identifier,
    )

    graphs = []
    for type_permutation in type_permutations:
        new_graph = Graph(pattern_graph)
        for (v, vertex_type) in zip(new_graph.vertices(), type_permutation):
            new_graph.vp.type[v] = vertex_type

        graphs.append(new_graph)
    return graphs

def write_dummy_graph(N: int = 10, cyclic: bool = True) -> Tuple[str, Any]:
    graph_file = 'data/dummy.gt'

    g = Graph(directed=False)
    vertices = [g.add_vertex() for n in range(0, N)]
    vertices_pairs = list(zip(vertices, vertices[1:] + (vertices[0:1] if cyclic else [])))
    [g.add_edge(v_1, v_2) for (v_1, v_2) in vertices_pairs]

    vertex_types = g.new_vertex_property("string")
    g.vertex_properties['type'] = vertex_types
    for (i, v) in enumerate(g.vertices()):
        vertex_types[v] = 'C' + str(i)

    g.save(graph_file)
    return (graph_file, vertex_types)

def draw_graph(
    graph: Graph,
    pos: Optional[PropertyMap] = None,
    fnme: str = 'graph',
    force_regen: bool = False,
    output_size: Tuple[float, float] = (400, 400),
    default_extension: str = '.pdf',
) -> None:
    try:
        vertex_text = graph.vertex_properties['type']
    except:
        vertex_text = graph.vertex_index

    try:
        vertex_color = graph.vertex_properties['color']
    except:
        vertex_color = graph.vertex_index

    try:
        edge_text = graph.edge_properties['type']
    except:
        edge_text = ''

    if not fnme.endswith(default_extension):
        fnme += default_extension

    if exists(fnme) and not force_regen:
        return
    else:
        graph_draw(
            graph,
            pos=pos,
            vertex_text=vertex_text,
            vertex_font_size=18,
            vertex_fill_color=vertex_color,
            edge_text=edge_text,
            edge_font_size=18,
            output_size=output_size,
            output=fnme,
        )

def get_interpreted_pattern_graphs() -> List[Tuple[Moiety, List[Graph]]]:
    pattern_graphs = [pattern_graph_for_graph_pattern(pattern) for (moiety, pattern) in list(PATTERNS.items())]
    interpreted_pattern_graphs = [
        (moiety, graphs_for_pattern_graph(pattern_graph, pattern_identifier=moiety))
        for (moiety, pattern_graph) in zip(list(PATTERNS.keys()), pattern_graphs)
    ]

    if DRAW_PATTERN_GRAPHS:
        for (moiety, graph_list) in interpreted_pattern_graphs:
            [
                draw_graph(
                    graph,
                    fnme=join(ROOT_DIR, 'patterns', moiety.replace(' ', '_') + '_' + str(i)),
                )
                for (i, graph) in enumerate(graph_list)
            ]

    return interpreted_pattern_graphs

def does_graph_match_graph(super_graph: Graph, sub_graph: Graph) -> bool:
    return topology.subgraph_isomorphism(
        sub_graph,
        super_graph,
        vertex_label=(
            sub_graph.vertex_properties['type'],
            super_graph.vertex_properties['type'],
        ),
        generator=False,
    )

def moieties_in_graph(super_graph: Graph, interpreted_pattern_graphs) -> List[Moiety]:
    def match(moiety: Moiety, graph_list: List[Graph]) -> bool:
        return any(
            does_graph_match_graph(super_graph, pattern_graph)
            for (i, pattern_graph) in enumerate(graph_list)
        )

    return [
        moiety
        for (moiety, graph_list) in interpreted_pattern_graphs
        if match(moiety, graph_list)
    ]

TEST_PDBS = glob('data/*.pdb')

def test_atom_class_parsing() -> None:
    TEST_DATA = (
        ('C', ['C2', 'C3', 'C4',]),
        ('C3', ['C3',]),
        ('C{4,5}', ['C4', 'C5',]),
        ('C4|H1', ['C4', 'H1',]),
        ('CL', ['CL1']),
        ('CL{4,5}', ['CL4', 'CL5']),
    )

    for test_class, test_result in TEST_DATA:
        r = list(types_and_valences_for_class(test_class))
        assert r == test_result, 'Error: pattern "{2}": {0} != {1}'.format(r, test_result, test_class)
    exit()

def moieties_in_pdb(pdb_str: str, should_draw_graph: bool = True, should_dump_graph: bool = False, interpreted_pattern_graphs: Optional[Any] = None, pdb_file: Optional[str] = None) -> List[Moiety]:
    if interpreted_pattern_graphs is None:
        interpreted_pattern_graphs = get_interpreted_pattern_graphs()

    molecule_graph = graph_from_pdb(pdb_str)

    if should_draw_graph and pdb_file is not None:
        draw_graph(
            molecule_graph,
            fnme=pdb_file.replace('.pdb', '.pdf')
        )

    if should_dump_graph and pdb_file is not None:
        molecule_graph.save(pdb_file.replace('.pdb', '.gt'))

    return moieties_in_graph(molecule_graph, interpreted_pattern_graphs)

def moieties_in_pdb_file(pdb_file: str, should_draw_graph: bool = True, should_dump_graph: bool = False, interpreted_pattern_graphs: Optional[Any] = None) -> List[Moiety]:
    with open(pdb_file) as fh:
        pdb_str = fh.read()

    return moieties_in_pdb(
        pdb_str,
        should_draw_graph=should_draw_graph,
        should_dump_graph=should_dump_graph,
        interpreted_pattern_graphs=interpreted_pattern_graphs,
        pdb_file=pdb_file,
    )

def enforce_disjoint_patterns(interpreted_pattern_graphs: List[Tuple[Moiety, List[Graph]]]) -> None:
    print('''INFO: Will make sure not two moieties are embedded in each others. NB: Takes several minutes ...''')
    for ((moiety_1, graphs_1), (moiety_2, graphs_2)) in combinations(interpreted_pattern_graphs, r=2):
        if any(does_graph_match_graph(graph_1, graph_2) or does_graph_match_graph(graph_2, graph_1) for graph_1 in graphs_1 for graph_2 in graphs_2):
            print('ERROR: Found conflicting moieties: {0}'.format([moiety_1, moiety_2]))
    print('''INFO: Success''')

if __name__ == '__main__':
    interpreted_pattern_graphs = get_interpreted_pattern_graphs()

    enforce_disjoint_patterns(interpreted_pattern_graphs)

    if True:
        test_atom_class_parsing()

    for test_pdb in TEST_PDBS:
        print(test_pdb)
        print(moieties_in_pdb_file(test_pdb, interpreted_pattern_graphs=interpreted_pattern_graphs))
