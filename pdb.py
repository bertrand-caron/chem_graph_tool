from itertools import cycle, product
from functools import reduce
try:
    from graph_tool.all import Graph, load_graph, graph_draw
    from graph_tool import topology
except:
    from sys import stderr
    stderr.write('Please install graph-tool')
    raise

from atb_helpers.pdb import is_pdb_atom_line, is_pdb_connect_line, pdb_fields

def type_identifier_for(atom_type, valence):
    return '{0}{1}'.format(
        atom_type,
        valence,
    )

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
                        lambda i_j: (i_j[0] - 1, i_j[1] - 1),
                        zip(cycle(connect_ids[0:1]), connect_ids[1:]),
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
        vertex_types[v] = type_identifier_for(
            fields[11].strip().upper(),
            get_valence(atom_id),
        )
        vertices.append(v)

    for (i, j) in connects:
        g.add_edge(vertices[i], vertices[j])

    return g
