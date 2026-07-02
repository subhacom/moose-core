# -*- coding: utf-8 -*-

__author__           = "Dilawar Singh"
__copyright__        = "Copyright 2013, NCBS Bangalore"
__credits__          = ["NCBS Bangalore", "Bhalla Lab"]
__license__          = "GPL"
__version__          = "1.0.0"
__maintainer__       = "Dilawar Singh"
__email__            = "dilawars@ncbs.res.in"
__status__           = "Development"

import numpy as np
import matplotlib.pyplot as plt
import moose
import moose.print_utils as pu
import re
from collections import OrderedDict

def plotAscii(yvec, xvec = None, file=None):
    """Plot two list-like object in terminal using gnuplot.
    If file is given then save data to file as well.
    """
    if xvec is None:
        plotInTerminal(yvec, list(range( len(yvec))), file=file)
    else:
        plotInTerminal(yvec, xvec, file=file)

def plotInTerminal(yvec, xvec = None, file=None):
    '''
    Plot given vectors in terminal using gnuplot.

    If file is not None then write the data to a file.
    '''
    import subprocess
    g = subprocess.Popen(["gnuplot"], stdin=subprocess.PIPE)
    g.stdin.write("set term dumb 100 25\n")
    g.stdin.write("plot '-' using 1:2 title '{}' with linespoints\n".format(file))
    if file:
        saveAsGnuplot(yvec, xvec, file=file)
    for i,j in zip(xvec, yvec):
        g.stdin.write("%f %f\n" % (i, j))
    g.stdin.write("\n")
    g.stdin.flush()

def xyToString( yvec, xvec, sepby = ' '):
    """ Given two list-like objects, returns a text string.
    """
    textLines = []
    for y, x in zip( yvec, xvec ):
        textLines.append("{}{}{}".format(y, sepby, x))
    return "\n".join(textLines)


def saveNumpyVec( yvec, xvec, file):
    """save the numpy vectors to a data-file

    """
    if file is None:
        return
    print(("[INFO] Saving plot data to file {}".format(file)))
    textLines = xyToString( yvec, xvec)
    with open(file, "w") as dataF:
        dataF.write(textLines)

def saveAsGnuplot( yvec, xvec, file):
    ''' Save the plot as stand-alone gnuplot script '''
    if file is None:
        return
    print(("[INFO] Saving plot data to a gnuplot-script: {}".format(file)))
    dataText = xyToString( yvec, xvec )
    text = []
    text.append("#!/bin/bash")
    text.append("gnuplot << EOF")
    text.append("set term post eps")
    text.append("set output \"{0}.eps\"".format(file))
    text.append("plot '-' using 0:1 title '{0}'".format(file))
    text.append(dataText)
    text.append("EOF")
    with open(file+".gnuplot","w") as gnuplotF:
        gnuplotF.write("\n".join(text))

def scaleVector(vec, scaleF):
    """ Scale a vector by a factor """
    if scaleF == 1.0 or scaleF is None:
        return vec
    else:
        return [ x*scaleF for x in vec ]

def scaleAxis(xvec, yvec, scaleX, scaleY):
    """ Multiply each elements by factor """
    xvec = scaleVector( xvec, scaleX )
    yvec = scaleVector( yvec, scaleY )
    return xvec, yvec

def reformatTable(table, kwargs):
    """ Given a table return x and y vectors with proper scaling """
    clock = moose.Clock('/clock')
    if isinstance(table, moose.Table):
        vecY = table.vector
        vecX = np.arange(0, clock.currentTime, len(vecY))
    elif isinstance(table, tuple):
        vecX, vecY = table
    return (vecX, vecY)

def plotTable(table, **kwargs):
    """Plot a given table. It plots table.vector

    This function can scale the x-axis. By default, y-axis and x-axis scaling is
    done by a factor of 1.

    Pass 'xscale' and/or 'yscale' argument to function to modify scales.

    """
    if not isinstance(table, moose.Table):
        msg = "Expected moose.Table, got {}".format( type(table) )
        raise TypeError(msg)

    vecX, vecY = reformatTable(table, kwargs)
    plt.plot(vecX, vecY, label = kwargs.get('label', ""))
    # This may not be available on older version of matplotlib.
    try:
        plt.legend(loc='best', framealpha=0.4)
    except TypeError:
        plt.legend(loc='best')

def plotTables(tables, outfile=None, **kwargs):
    """Plot a dict of tables.

    Each axis will be labeled with dict keys. The dict values must be
    moose.Table/moose.Table2.

    :param outfile: Default `None`. If given, plot will be saved to this filepath.
    :param grid: A tuple with (cols, rows), default is (len(tables)//2+1, 2)
    :param figsize: Size of figure (W, H) in inches. Default (10, 1.5*len(tables))
    """
    if not isinstance(tables, dict): raise TypeError("Expected a dict of moose.Table")
    plt.figure(figsize=kwargs.get('figsize', (10, 1.5*len(tables))))
    subplot = kwargs.get('subplot', True)
    gridSize = kwargs.get('grid', (len(tables)//2+1, 2))
    for i, tname in enumerate(tables):
        if subplot:
            if gridSize[0] > 9: raise ValueError("gridSize rows must be <= 9, got %d" % gridSize[0])
            plt.subplot(100*gridSize[0]+10*gridSize[1]+(i+1))
        yvec = tables[tname].vector
        xvec = np.linspace(0, moose.Clock('/clock').currentTime, len(yvec))
        plt.plot(xvec, yvec, label=tname)

        # This may not be available on older version of matplotlib.
        try:
            plt.legend(loc='best', framealpha=0.4)
        except TypeError:
            plt.legend(loc='best')

    plt.tight_layout()
    if outfile:
        pu.dump("PLOT", "Saving plots to file {}".format(outfile))
        try:
            plt.savefig(outfile, transparent=True)
        except Exception as e:
            pu.dump("WARN", "Failed to save figure. Errror %s"%e)
            plt.show()
    else:
        plt.show()

def plotVector(vec, xvec = None, **options):
    """plotVector: Plot a given vector. On x-axis, plot the time.

    :param vec: Given vector.
    :param **kwargs: Optional to pass to maplotlib.
    """

    if not isinstance(vec, np.ndarray): raise TypeError("Expected np.ndarray, got %s" % type(vec))
    legend = options.get('legend', True)

    if xvec is None:
        clock = moose.Clock('/clock')
        xx = np.linspace(0, clock.currentTime, len(vec))
    else:
        xx = xvec[:]

    if len(xx) != len(vec): raise ValueError("Expecting %d points, got %d" % (len(vec), len(xx)))

    plt.plot(xx, vec, label=options.get('label', ''))
    if legend:
        # This may not be available on older version of matplotlib.
        try:
            plt.legend(loc='best', framealpha=0.4)
        except TypeError:
            plt.legend(loc='best')

    if xvec is None:
        plt.xlabel('Time (sec)')
    else:
        plt.xlabel(options.get('xlabel', ''))

    plt.ylabel = options.get('ylabel', '')
    plt.title(options.get('title', ''))

    if(options.get('legend', True)):
        try:
            plt.legend(loc='best', framealpha=0.4, prop={'size' : 9})
        except TypeError:
            plt.legend(loc='best', prop={'size' : 9})


def saveRecords(records, xvec = None, **kwargs):
    """saveRecords
    Given a dictionary of data with (key, numpy array) pair, it saves them to a
    file 'outfile'

    :param outfile
    :param dataDict:
    :param **kwargs:
        comment: Adds comments below the header.
    """
    if len(records) == 0:
        pu.warn("No data in dictionary to save.")
        return False

    outfile = kwargs.get('outfile', 'data.moose')
    clock = moose.Clock('/clock')
    if clock.currentTime <= 0: raise RuntimeError("Simulation has not been run (currentTime=%g)" % clock.currentTime)
    yvecs = [ ]
    text = "time," + ",".join([ str(x) for x in records ])
    for k in records:
        try:
            yvec = records[k].vector
        except AttributeError as e:
            yevc = records[k]
        yvecs.append(yvec)
    xvec = np.linspace(0, clock.currentTime, len(yvecs[0]))
    yvecs = [ xvec ] + yvecs
    if kwargs.get('comment', ''):
        text += ("\n"  + kwargs['comment'] )
    np.savetxt(outfile, np.array(yvecs).T, delimiter=',' , header = text)
    pu.info("Done writing data to %s" % outfile)

def plotRecords(records, xvec = None, **kwargs):
    """Wrapper
    """
    dataDict = OrderedDict( )
    try:
        for k in sorted(records.keys(), key=str.lower):
            dataDict[k] = records[k]
    except Exception:
        dataDict = records

    outfile = kwargs.get('outfile', None)
    subplot = kwargs.get('subplot', False)
    filters = [ x.lower() for x in kwargs.get('filter', [])]

    plt.figure(figsize=(10, 1.5*len(dataDict)))
    #plt.rcParams.update( { 'font-size' : 10 } )
    for i, k in enumerate(dataDict):
        pu.info("+ Plotting for %s" % k)
        plotThis = False
        if not filters: plotThis = True
        for accept in filters:
            if accept in k.lower():
                plotThis = True
                break

        if plotThis:
            if not subplot:
                yvec = dataDict[k].vector
                plotVector(yvec, xvec, label=k, **kwargs)
            else:
                plt.subplot(len(dataDict), 1, i+1)
                yvec = dataDict[k].vector
                plotVector(yvec, xvec, label=k, **kwargs)

    # title in Image.
    if 'title' in kwargs:
        plt.title(kwargs['title'])

    if subplot:
        try:
            plt.tight_layout()
        except: pass

    if outfile:
        pu.info("Writing plot to %s" % outfile)
        plt.savefig("%s" % outfile, transparent=True)
    else:
        plt.show()
    plt.close( )


def plotTablesByRegex( regex = '.*', **kwargs ):
    """plotTables Plot all moose.Table/moose.Table2 matching given regex. By
    default plot all tables. Table names must be unique. Table name are used as
    legend.

    :param regex: Python regular expression to be matched.
    :param **kwargs:
        - subplot = True/False; if True, each Table is plotted in a subplot.
        - outfile = filepath; If given, plot will be saved to this path.
    """

    tables = moose.wildcardFind( '/##[TYPE=Table]' )
    tables += moose.wildcardFind( '/##[TYPE=Table2]' )
    toPlot = OrderedDict( )
    for t in sorted(tables, key = lambda x: x.name):
        if re.search( regex, t.name ):
            toPlot[ t.name ] = t
    return plotRecords( toPlot, None, **kwargs )


# ---------------------------------------------------------------------------
# Morphology display
# ---------------------------------------------------------------------------

def _segmentColor(comp, color):
    """Resolve the drawing colour for one compartment.

    `color` may be a plain matplotlib colour, a callable ``comp -> colour``,
    or a dict looked up by compartment name (with 'k' as the fallback).
    """
    if callable(color):
        return color(comp)
    if isinstance(color, dict):
        return color.get(comp.name, 'k')
    return color


# Colours for the SWC structure types soma/axon/basal/apical, keyed by the
# leading token of the compartment name.  Used by ``color='type'``.
MORPH_TYPE_COLORS = {
    'soma':   'k',   # soma
    'axon':   'g',   # axon
    'apical': 'r',   # apical dendrite
    'basal':  'b',   # basal dendrite
    'dend':   'b',   # generic dendrite
    None:     '0.5', # anything else
}


def compartmentTypeColor(comp, colors=MORPH_TYPE_COLORS):
    """Return a colour for `comp` inferred from its name.

    Matches the leading alphabetic token of the compartment name against
    common conventions (soma, axon, apical/apic/prim, basal/dend); unmatched
    names get the ``None`` entry of `colors`.  Pass as ``color=`` to
    :func:`plotMorphology` (or via the ``color='type'`` shorthand).
    """
    name = re.split(r'[\d_\[]', comp.name, 1)[0].lower()
    if name.startswith(('apic', 'prim', 'glom')):
        name = 'apical'
    elif name.startswith('ax'):
        name = 'axon'
    elif name.startswith('bas'):
        name = 'basal'
    return colors.get(name, colors.get(None, 'k'))


def _paddedLimits(mins, maxs, margin=0.05):
    """Return (lo, hi) axis limits with a margin, never zero-width.

    A flat (single-line) morphology gives ``min == max`` on some axes, which
    makes matplotlib warn about a singular transform.  Degenerate axes are
    padded using the largest non-degenerate span (or 1.0 if all are flat) so
    every axis has a finite range.
    """
    mins = np.asarray(mins, dtype=float)
    maxs = np.asarray(maxs, dtype=float)
    spans = maxs - mins
    ref = spans.max()
    if ref <= 0:
        ref = max(float(np.abs(maxs).max()), 1.0)
    pad = np.where(spans > 0, spans * margin, ref * margin)
    return mins - pad, maxs + pad


def plotMorphology(root, ax=None, projection='3d', color='type',
                   diam_scale=None, linewidth=1.0, autoscale=True, **kwargs):
    """Display a neuron's morphology by traversing axial/raxial connectivity.

    Each compartment is drawn as a line from its parent's distal point to its
    own distal point ``(x, y, z)``.  Topology is taken from the compartment
    connection messages (see :func:`moose.network_utils.getCompartments`), not
    the element tree, so it works regardless of how the compartments are
    arranged under their container.

    :param root: path string or moose element of either

        * a container (e.g. a ``Neuron`` returned by :func:`moose.loadSwc`)
          holding the compartments, or
        * a single ``Compartment`` / ``SymCompartment`` (e.g. the soma), in
          which case the whole connected cell is traversed.
    :param ax: matplotlib Axes to draw into; created if omitted (a 3-D axes
        when ``projection='3d'``).
    :param projection: ``'3d'`` for a 3-D view, or ``'xy'`` / ``'xz'`` /
        ``'yz'`` for a 2-D projection onto the named plane.
    :param color: a single matplotlib colour, a ``{compartment_name: colour}``
        dict, or a callable ``comp -> colour``.  The shorthand ``'type'``
        colours each compartment by its inferred type (see
        :func:`compartmentTypeColor`).
    :param diam_scale: if given, each segment's line width is
        ``diameter * diam_scale`` (diameter is in metres, so e.g.
        ``diam_scale=1e6`` gives width in microns).  If ``None`` (default),
        ``linewidth`` is used for every segment.
    :param linewidth: constant line width used when ``diam_scale`` is ``None``.
    :param autoscale: if True, set data limits / equal aspect from the extent.
    :returns: the matplotlib Axes the morphology was drawn into.
    """
    from moose.network_utils import getCompartments, getCompartmentNeighbors, \
        _PARENT_LINK_FIELDS

    axisNames = ('x', 'y', 'z')
    if color == 'type':
        color = compartmentTypeColor
    comps = getCompartments(root)
    if not comps:
        raise ValueError('No compartments found under {!r}'.format(root))

    segments, colors, widths = [], [], []
    for comp in comps:
        parents = getCompartmentNeighbors(comp, _PARENT_LINK_FIELDS)
        distal = (comp.x, comp.y, comp.z)
        if parents:
            p = parents[0]
            proximal = (p.x, p.y, p.z)
        else:
            # Root compartment: fall back to its own proximal end (x0,y0,z0).
            proximal = (comp.x0, comp.y0, comp.z0)
        segments.append((proximal, distal))
        colors.append(_segmentColor(comp, color))
        widths.append(comp.diameter * diam_scale if diam_scale else linewidth)

    if projection == '3d':
        from mpl_toolkits.mplot3d.art3d import Line3DCollection
        if ax is None:
            ax = plt.figure().add_subplot(projection='3d')
        ax.add_collection3d(
            Line3DCollection(segments, colors=colors, linewidths=widths,
                             **kwargs))
        ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
        if autoscale:
            pts = np.array([pt for seg in segments for pt in seg])
            los, his = _paddedLimits(pts.min(axis=0), pts.max(axis=0))
            ax.set_xlim(los[0], his[0])
            ax.set_ylim(los[1], his[1])
            ax.set_zlim(los[2], his[2])
            try:
                ax.set_box_aspect(his - los)
            except (AttributeError, ValueError):
                pass
    else:
        from matplotlib.collections import LineCollection
        try:
            h, v = axisNames.index(projection[0]), axisNames.index(projection[1])
        except (ValueError, IndexError):
            raise ValueError(
                "projection must be '3d', 'xy', 'xz' or 'yz', got {!r}".format(
                    projection))
        if ax is None:
            ax = plt.figure().add_subplot()
        segs2d = [((s[0][h], s[0][v]), (s[1][h], s[1][v])) for s in segments]
        ax.add_collection(
            LineCollection(segs2d, colors=colors, linewidths=widths, **kwargs))
        ax.set_xlabel(projection[0]); ax.set_ylabel(projection[1])
        if autoscale:
            pts = np.array([pt for seg in segs2d for pt in seg])
            los, his = _paddedLimits(pts.min(axis=0), pts.max(axis=0))
            ax.set_xlim(los[0], his[0])
            ax.set_ylim(los[1], his[1])
            ax.set_aspect('equal')
    return ax


def _springLayout(n, edges, dim=2, iterations=50, k=None, seed=0):
    """Fruchterman-Reingold force-directed layout in pure numpy.

    Distributes `n` nodes so that connected nodes are pulled together and all
    nodes repel each other, giving a readable graph even when the original
    coordinates are collinear.  Implemented here (rather than via networkx) to
    avoid adding a dependency.

    :param n: number of nodes.
    :param edges: iterable of ``(i, j)`` index pairs (undirected).
    :param dim: 2 or 3 output dimensions.
    :param iterations: number of relaxation steps.
    :param k: optimal edge length; defaults to ``1/sqrt(n)``.
    :param seed: RNG seed for the initial random placement.
    :returns: ``(n, dim)`` array of node positions in roughly [0, 1].
    """
    rng = np.random.default_rng(seed)
    pos = rng.random((n, dim))
    if n < 2:
        return pos
    if k is None:
        k = 1.0 / np.sqrt(n)
    adj = np.zeros((n, n))
    for i, j in edges:
        adj[i, j] = adj[j, i] = 1.0
    temp = 0.1                       # max displacement per step (cools to 0)
    cooling = temp / (iterations + 1)
    for _ in range(iterations):
        delta = pos[:, None, :] - pos[None, :, :]          # (n, n, dim)
        dist = np.sqrt((delta * delta).sum(axis=-1))
        # Diagonal delta is zero (contributes no force); set a finite non-zero
        # distance there to avoid 0/0, and floor coincident nodes.
        np.fill_diagonal(dist, 1.0)
        dist = np.maximum(dist, 1e-9)
        # Repulsion on every pair (k^2/d), attraction along edges (-d^2/k).
        force = (k * k) / dist - adj * (dist * dist) / k
        disp = (delta / dist[..., None] * force[..., None]).sum(axis=1)
        length = np.sqrt((disp * disp).sum(axis=-1))
        length = np.where(length < 1e-9, 1.0, length)
        pos += disp / length[:, None] * np.minimum(length, temp)[:, None]
        temp -= cooling
    return pos


def plotMorphologyGraph(root, ax=None, dim=2, color='type', node_size=20,
                        edge_color='0.6', with_labels=False,
                        iterations=50, k=None, seed=0, **kwargs):
    """Draw a neuron as a force-directed graph (schematic topology).

    Unlike :func:`plotMorphology`, this ignores the physical ``(x, y, z)``
    coordinates and lays the compartments out with a force-directed algorithm
    (see :func:`_springLayout`).  It is meant for abstract/collinear models
    (e.g. ``traub91_CA3``) that collapse to a line when drawn to scale, letting
    you see the branching structure instead.

    Topology is taken from axial/raxial connectivity via
    :func:`moose.network_utils.getCompartments`.

    :param root: path string or moose element -- a container (Neuron) or a
        single Compartment (the whole connected cell is used).
    :param ax: matplotlib Axes to draw into; created if omitted (3-D when
        ``dim=3``).
    :param dim: 2 for a 2-D graph, 3 for a 3-D graph.
    :param color: node colour -- a matplotlib colour, dict, callable, or the
        ``'type'`` shorthand (colour by inferred compartment type).
    :param node_size: marker size for the nodes (``0`` to hide them).
    :param edge_color: colour of the connecting edges.
    :param with_labels: if True, annotate each node with its compartment name.
    :param iterations, k, seed: forwarded to :func:`_springLayout`.
    :returns: the matplotlib Axes.
    """
    from moose.network_utils import getCompartments, getCompartmentNeighbors

    if color == 'type':
        color = compartmentTypeColor
    comps = getCompartments(root)
    if not comps:
        raise ValueError('No compartments found under {!r}'.format(root))

    index = {c.path: i for i, c in enumerate(comps)}
    edges = set()
    for c in comps:
        i = index[c.path]
        for nbr in getCompartmentNeighbors(c):
            j = index.get(nbr.path)
            if j is not None and j != i:
                edges.add((min(i, j), max(i, j)))
    edges = list(edges)

    pos = _springLayout(len(comps), edges, dim=dim,
                        iterations=iterations, k=k, seed=seed)
    node_colors = [_segmentColor(c, color) for c in comps]

    if dim == 3:
        from mpl_toolkits.mplot3d.art3d import Line3DCollection
        if ax is None:
            ax = plt.figure().add_subplot(projection='3d')
        segs = [(pos[i], pos[j]) for i, j in edges]
        ax.add_collection3d(Line3DCollection(segs, colors=edge_color))
        ax.scatter(pos[:, 0], pos[:, 1], pos[:, 2], s=node_size,
                   c=node_colors, depthshade=False, **kwargs)
        if with_labels:
            for c, i in ((c, index[c.path]) for c in comps):
                ax.text(pos[i, 0], pos[i, 1], pos[i, 2], c.name)
        los, his = _paddedLimits(pos.min(axis=0), pos.max(axis=0))
        ax.set_xlim(los[0], his[0]); ax.set_ylim(los[1], his[1])
        ax.set_zlim(los[2], his[2])
    else:
        from matplotlib.collections import LineCollection
        if ax is None:
            ax = plt.figure().add_subplot()
        segs = [(pos[i], pos[j]) for i, j in edges]
        ax.add_collection(LineCollection(segs, colors=edge_color, zorder=1))
        ax.scatter(pos[:, 0], pos[:, 1], s=node_size, c=node_colors,
                   zorder=2, **kwargs)
        if with_labels:
            for c in comps:
                i = index[c.path]
                ax.annotate(c.name, (pos[i, 0], pos[i, 1]))
        los, his = _paddedLimits(pos.min(axis=0), pos.max(axis=0))
        ax.set_xlim(los[0], his[0]); ax.set_ylim(los[1], his[1])
        ax.set_aspect('equal')
        ax.axis('off')
    return ax
