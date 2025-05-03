# test_neighbors.py ---
#
# Filename: test_neighbors.py
# Description:
# Author: Subhasis Ray
# Created: Thu Apr 24 22:04:16 2025 (+0530)
#

# Code:


import moose
import math
import pytest


@pytest.fixture
def container():
    """
    Setup:
        Create a Neutral object as model container and ce to it before
        yielding.

    Teardown:
        ce back to parent and delete the model
        container.

    """
    ret = moose.Neutral('/test')
    moose.ce(ret)
    yield ret
    moose.ce('..')
    moose.delete(ret)


def test_all_neighbors(container):
    neighbors = moose.neighbors('/')
    names = {node.path for node in neighbors}
    assert names == {'/Msgs', '/clock', '/classes', '/postmaster', '/test'}


def test_incoming(container):
    neighbors = moose.neighbors('/', direction=moose.INMSG)
    assert neighbors == []


def test_in_neighbors(container):
    comp = moose.Compartment(f'{container.path}/comp')
    tab = moose.Table(f'{container.path}/tab')
    moose.connect(tab, 'requestOut', comp, 'getVm')
    neighbors = moose.neighbors(comp, direction=moose.INMSG)
    assert set(neighbors) == {
        moose.element('/clock'),
        moose.element(container),
        tab,
    }


def test_infield_neighbors(container):
    comp = moose.Compartment(f'{container.path}/comp')
    tab = moose.Table(f'{container.path}/tab')
    moose.connect(tab, 'requestOut', comp, 'getVm')
    neighbors = moose.neighbors(comp, field='getVm', direction=moose.INMSG)
    assert set(neighbors) == {tab}


def test_out_neighbors(container):
    comp = moose.Compartment(f'{container.path}/comp')
    tab = moose.Table(f'{container.path}/tab')
    moose.connect(tab, 'requestOut', comp, 'getVm')
    neighbors = moose.neighbors(tab, direction=moose.OUTMSG)
    assert set(neighbors) == {comp}


def test_outfield_neighbors(container):
    comp = moose.Compartment(f'{container.path}/comp')
    tab = moose.Table(f'{container.path}/tab')
    moose.connect(tab, 'requestOut', comp, 'getVm')
    neighbors = moose.neighbors(
        tab, field='requestOut', direction=moose.OUTMSG
    )
    assert set(neighbors) == {comp}


#
# test_neighbors.py ends here
