# -*- coding: utf-8 -*-

# Some network analysis utilities.

__author__           = "Dilawar Singh"
__copyright__        = "Copyright 2018-19, NCBS Bangalore"
__maintainer__       = "Dilawar Singh"
__email__            = "dilawars@ncbs.res.in"

import sys
import hashlib
import moose._moose as _moose
import re

import logging
logger_ = logging.getLogger('moose.utils.graph')

pathPat_ = re.compile(r'.+?\[\d+\]$')

def morphologyToGraphviz(filename=None, pat='/##[TYPE=Compartment]'):
    '''Write Electrical network to a dot graph.

    Params:

    :filename: Default None. Write graphviz file to this path. If None, write to
        stdout.
    :pat: Compartment path. By default, search for all moose.Compartment.
    '''

    def fix(path):
        '''Fix a given path so it can be written to a graphviz file'''
        # If no [0] is at end of the path then append it.
        global pathPat_
        if not pathPat_.match(path):
            path = path + '[0]'
        return path

    compList = _moose.wildcardFind(pat)
    if not compList:
        logger_.warning("No compartment found")

    dot = []
    dot.append("digraph G {")
    dot.append("\tconcentrate=true;")
    for c in compList:
        if c.neighbors['raxial']:
            for n in c.neighbors['raxial']:
                lhs = fix(c.path)
                rhs = fix(n.path)
                dot.append('\t"{}" -> "{}";'.format(lhs, rhs))
        elif c.neighbors['axial']:
            for n in c.neighbors['axial']:
                lhs = fix(c.path)
                rhs = fix(n.path)
                dot.append('\t"{}" -> "{}" [dir=back];'.format(lhs, rhs))
        else:
            p = fix(c.path)
            dot.append('\t"{}"'.format(p))

    dot.append('}')
    dot = '\n'.join(dot)
    if not filename:
        print(dot)
    else:
        with open(filename, 'w') as graphviz:
            logger_.info("Writing compartment topology to file {}".format(filename))
            graphviz.write(dot)
    return True


# ---------------------------------------------------------------------------
# Compartment connectivity traversal
# ---------------------------------------------------------------------------

# Shared-message fields linking a compartment to its parent (proximal side).
# Asymmetric Compartment uses `axial`; SymCompartment uses `proximal`.
_PARENT_LINK_FIELDS = ('axial', 'proximal')

# All shared-message fields linking a compartment to any neighbour, used when
# walking connectivity from a seed compartment.  Asymmetric compartments use
# axial/raxial; SymCompartment adds proximal/distal/sibling.
_ALL_LINK_FIELDS = ('axial', 'raxial', 'proximal', 'distal', 'sibling')

# className -> set of shared-message field names it actually defines.  Cached
# so we probe only fields that exist (querying an absent neighbour field emits
# a C++ warning rather than raising).
_sharedFieldCache = {}


def _sharedFields(comp):
    cls = comp.className
    if cls not in _sharedFieldCache:
        _sharedFieldCache[cls] = set(_moose.getFieldNames(comp, 'sharedFinfo'))
    return _sharedFieldCache[cls]


def getCompartmentNeighbors(comp, fields=_ALL_LINK_FIELDS):
    """Return compartments directly linked to `comp` via any of `fields`.

    Topology comes from the compartment's connection messages, not the element
    tree, so this works whether compartments are nested or laid out flat under
    a container.  Fields the compartment's class does not define (e.g.
    `proximal` on an asymmetric Compartment) are skipped, so it works uniformly
    for Compartment and SymCompartment.

    :param comp: moose Compartment (or SymCompartment) element.
    :param fields: iterable of shared-message field names to follow.
        Defaults to all axial/raxial/proximal/distal/sibling links; pass
        ``_PARENT_LINK_FIELDS`` to get only parent-side neighbours.
    :returns: list of unique neighbouring compartment elements.
    """
    valid = _sharedFields(comp)
    out, seen = [], set()
    for field in fields:
        if field not in valid:
            continue
        for n in comp.neighbors[field]:
            el = _moose.element(n)
            if el.path not in seen:
                seen.add(el.path)
                out.append(el)
    return out


def getConnectedCompartments(seed):
    """Return every compartment connected to `seed`, via a breadth-first walk
    of axial/raxial connectivity.

    Passing any single compartment (typically the soma) therefore returns the
    whole connected cell, in breadth-first order from the seed.

    :param seed: moose Compartment (or SymCompartment) element.
    :returns: list of compartment elements including `seed`.
    """
    from collections import deque
    order, visited = [], {seed.path}
    queue = deque([seed])
    while queue:
        comp = queue.popleft()
        order.append(comp)
        for nbr in getCompartmentNeighbors(comp):
            if nbr.path not in visited:
                visited.add(nbr.path)
                queue.append(nbr)
    return order


def getCompartments(root):
    """Resolve `root` to the list of compartments it refers to.

    * A single Compartment / SymCompartment -> the whole connected cell,
      found by walking axial/raxial connectivity (see
      :func:`getConnectedCompartments`).
    * Any other element (Neuron, Neutral, ...) -> every compartment in its
      subtree.

    :param root: path string or moose element.
    :returns: list of compartment elements.
    """
    root = _moose.element(root)
    if _moose.wildcardFind('{}[ISA=Compartment]'.format(root.path)):
        return getConnectedCompartments(root)
    return list(_moose.wildcardFind('{}/##[ISA=Compartment]'.format(root.path)))


def chemicalReactionNetworkToGraphviz(compt, path=None):
    """Write chemical reaction network to a graphviz file.

    :param compt: Given compartment.
    :param filepath: Save to this filepath. If None, write to stdout.
    """
    dot = _crn(compt)
    if path is None:
        print(dot, file=sys.stdout)
        return
    with open(path, 'w') as f:
        f.write(dot)

# aliases
crnToDot = chemicalReactionNetworkToGraphviz
crnToGraphviz = chemicalReactionNetworkToGraphviz

def _fixLabel(name):
    name = name.replace('*', 'star')
    name = name.replace('.', '_')
    return "\"{}\"".format(name)

def _addNode(n, nodes, dot):
    node = hashlib.sha224(n.path.encode()).hexdigest()
    nodeType = 'pool'
    if isinstance(n, _moose.Reac) or isinstance(n, _moose.ZombieReac):
        node = 'r'+node
        nodeType = 'reac'
    else:
        node = 'p'+node
    if node in nodes:
        return node

    nLabel = n.name
    if nodeType == 'reac':
        nLabel = "kf=%g kb=%g"%(n.Kf, n.Kb)
        dot.append('\t%s [label="%s", kf=%g, kb=%g, shape=rect]' % (
            node, nLabel, n.Kf, n.Kb))
    else:
        dot.append('\t%s [label="%s", concInit=%g]' % (
            node, nLabel, n.concInit))
    return node

def _crn(compt):
    nodes = {}
    reacs = _moose.wildcardFind(compt.path+'/##[TYPE=Reac]')
    reacs += _moose.wildcardFind(compt.path+'/##[TYPE=ZombieReac]')
    dot = ['digraph %s {\n\t overlap=false' % compt.name ]
    for r in reacs:
        rNode = _addNode(r, nodes, dot)
        for s in r.neighbors['sub']:
            sNode = _addNode(s, nodes, dot)
            dot.append('\t%s -> %s' % (sNode, rNode))
        for p in r.neighbors['prd']:
            pNode = _addNode(p, nodes, dot)
            dot.append('\t%s -> %s' % (rNode, pNode))
    dot.append('}')
    return '\n'.join(dot)
