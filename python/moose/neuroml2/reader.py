# -*- coding: utf-8 -*-
from __future__ import print_function, division, absolute_import

# Description: NeuroML2 reader.
#     Implementation of reader for NeuroML 2 models.
#     TODO: handle morphologies of more than one segment...
# Author: Subhasis Ray, Padraig Gleeson
# Maintainer:  Dilawar Singh <dilawars@ncbs.res.in>, Padraig Gleeson
# Created: Wed Jul 24 15:55:54 2013 (+0530)
# Notes:
#    For update/log, please see git-blame documentation or browse the github
#    repo https://github.com/BhallaLab/moose-core

import os
import math
import logging
import numpy as np
import moose


try:
    import neuroml as nml
    import pyneuroml.pynml as pynml
except ImportError as error:
    raise ImportError(
        'Could not import neuroml/pyneuroml.'
        ' Please make sure you have pyneuroml installed'
        ' (`pip install pyneuroml`)'
    ) from error


from moose.neuroml2.hhfit import exponential2, sigmoid2, linoid2
from moose.neuroml2.units import SI


logger_ = logging.getLogger('moose.nml2')

RATE_FN_MAP = {
    'HHExpRate': exponential2,
    'HHSigmoidRate': sigmoid2,
    'HHSigmoidVariable': sigmoid2,
    'HHExpLinearRate': linoid2,
}

# Usually [Ca2+] is 50-100 nM and rise up to 100x during electrical
# activity. That is 1e4 nM = 1e-5 M, keeping 1e-4 as a safe
# bet. Reference: Grienberger & Konnerth, 2012. Neuron
CACONC_MAX = 1e-4

def _write_flattened_nml(doc, outfile):
    """_write_flattened_nml
    Concat all NML2 read by moose and generate one flattened NML file.
    Only useful when debugging.

    :param doc: NML document (nml.doc)
    :param outfile: Name of the output file.
    """
    import neuroml.writers

    neuroml.writers.NeuroMLWriter.write(doc, outfile)
    logger_.debug("Wrote flattened NML model to %s" % outfile)


def _gates_sorted(all_gates):
    """_gates_sorted

    Parameters
    ----------
    all_gates (list)
        List of all moose.HHChannel.gates

    Notes
    -----
    If the id of gates are subset of 'x', 'y' or 'z' then sort them so they load in
    X, Y or Z gate respectively. Otherwise do not touch them i.e. first gate
    will be loaded into X, second into Y and so on.
    """
    allMooseGates = ['x', 'y', 'z']
    allGatesDict = {g.id: g for g in all_gates}
    gateNames = [g.id.lower() for g in all_gates]
    if set(gateNames).issubset(set(allMooseGates)):
        sortedGates = []
        for gid in allMooseGates:
            sortedGates.append(allGatesDict.get(gid))
        return sortedGates
    return all_gates


def _unique(sequence):
    """Return a list with unique items from `sequence`"""
    return list(dict.fromkeys(sequence))


def _isConcDep(ct):
    """_isConcDep
    Check if componet is dependent on concentration. Most HHGates are
    dependent on voltage.

    :param ct: ComponentType
    :type ct: nml.ComponentType

    :return: True if Component is depenant on conc, False otherwise.
    """
    # logger_.debug(f"{'#' * 10} EXTENDS {ct.extends}")
    if 'ConcDep' in ct.extends:
        return True
    return False


def _findCaConcVariableName():
    """_findCaConcVariableName
    Find a suitable CaConc for computing HHGate tables.
    This is a hack, though it is likely to work in most cases.
    """
    caConcs = moose.wildcardFind('/library/##[TYPE=CaConc]')
    assert (
        len(caConcs) >= 1
    ), "No moose.CaConc found. Currently moose \
            supports HHChannel which depends only on moose.CaConc ."
    return caConcs[0].name


def sarea(comp):
    """sarea
    Return the surface area (2ϖrL) of compartment from length and diameter.

    :param comp: Compartment instance.
    :type comp: str
    :return: surface area of `comp`.
    :rtype: float
    """
    if comp.length > 0:
        return math.pi * comp.diameter * comp.length
    else:
        return math.pi * comp.diameter * comp.diameter


def xarea(compt):
    """xarea
    Return the cross sectional area (𝜋r²) from the diameter of the compartment.

    Note:
    ----
    How to do it for spherical compartment?

    :param compt: Compartment in moose.
    :type compt: moose.Compartment
    :return: cross sectional area.
    :rtype: float
    """
    return math.pi * (compt.diameter / 2.0) ** 2.0


def setRa(comp, resistivity):
    """Calculate total raxial from specific value `resistivity`"""
    if comp.length > 0:
        comp.Ra = resistivity * comp.length / xarea(comp)
    else:
        comp.Ra = resistivity * 4.0 / (comp.diameter * np.pi)


def setRm(comp, condDensity):
    """Set membrane resistance"""
    comp.Rm = 1 / (condDensity * sarea(comp))


def setEk(comp, erev):
    """Set reversal potential"""
    comp.setEm(erev)


def getSegments(nmlcell, component, sg_to_segments):
    """Get the list of segments the `component` is applied to"""
    sg = component.segment_groups
    if sg is None:
        segments = nmlcell.morphology.segments
    elif sg == 'all':
        segments = [seg for seglist in sg_to_segments.values() for seg in seglist]
    else:
        segments = sg_to_segments[sg]
    return _unique(segments)


class NML2Reader(object):
    """Reads NeuroML2 and creates MOOSE model.

    NML2Reader.read(filename) reads an NML2 model under `/library`
    with the toplevel name defined in the NML2 file.

    Example:

    >>> import moose
    >>> reader = moose.NML2Reader()
    >>> reader.read('moose/neuroml2/test_files/Purk2M9s.nml')

    creates a passive neuronal morphology `/library/Purk2M9s`.
    """

    def __init__(self, verbose=False):
        global logger_
        self.lunit = 1e-6  # micron is the default length unit
        self.verbose = verbose
        if self.verbose:
            logger_.setLevel(logging.DEBUG)
        self.doc = None
        self.filename = None
        self.nml_cells_to_moose = {}  # NeuroML object to MOOSE object
        self.nml_segs_to_moose = {}  # NeuroML object to MOOSE object
        self.nml_chans_to_moose = {}  # NeuroML object to MOOSE object
        self.nml_conc_to_moose = {}  # NeuroML object to MOOSE object
        self.moose_to_nml = {}  # Moose object to NeuroML object
        self.proto_cells = {}  # map id to prototype cell in moose
        self.proto_chans = {}  # map id to prototype channels in moose
        self.proto_pools = {}  # map id to prototype pools (Ca2+, Mg2+)
        self.includes = {}  # Included files mapped to other readers
        self.lib = (
            moose.element('/library')
            if moose.exists('/library')
            else moose.Neutral('/library')
        )
        self.id_to_ionChannel = {}
        self._cell_to_sg = (
            {}
        )  # nml cell to dict - the dict maps segment groups to segments
        self._variables = {}

        self.cells_in_populations = {}
        self.pop_to_cell_type = {}
        self.seg_id_to_comp_name = {}
        self.paths_to_chan_elements = {}
        self.network = None

    def read(self, filename, symmetric=True):
        filename = os.path.realpath(filename)
        self.doc = nml.loaders.read_neuroml2_file(
            filename, include_includes=True, verbose=self.verbose
        )

        self.filename = filename

        logger_.info('Parsed the NeuroML2 file: %s' % filename)
        if self.verbose:
            _write_flattened_nml(self.doc, '%s__flattened.xml' % self.filename)

        if len(self.doc.networks) >= 1:
            self.network = self.doc.networks[0]
            moose.celsius = self._getTemperature()

        self.importConcentrationModels(self.doc)
        self.importIonChannels(self.doc)
        self.importInputs(self.doc)

        for cell in self.doc.cells:
            # logger_.debug(f"{'%' * 10} Creating cell prototype {cell}")
            self.createCellPrototype(cell, symmetric=symmetric)

        if len(self.doc.networks) >= 1:
            self.createPopulations()
            self.createInputs()

    def _getTemperature(self):
        if self.network is not None:
            if self.network.type == "networkWithTemperature":
                return SI(self.network.temperature)
            else:
                # Why not, if there's no temp dependence in nml..?
                return 0
        return SI('25')

    def getCellInPopulation(self, pop_id, index):
        return self.cells_in_populations[pop_id][index]

    def getComp(self, pop_id, cellIndex, segId):
        return moose.element(
            '%s/%s/%s/%s'
            % (
                self.lib.path,
                pop_id,
                cellIndex,
                self.seg_id_to_comp_name[self.pop_to_cell_type[pop_id]][segId],
            )
        )

    def createPopulations(self):
        for pop in self.network.populations:
            ep = '%s/%s' % (self.lib.path, pop.id)
            mpop = moose.element(ep) if moose.exists(ep) else moose.Neutral(ep)
            self.cells_in_populations[pop.id] = {}
            for i in range(pop.size):
                logger_.info(
                    "Creating %s/%s instances of cell: %s under %s"
                    % (i, pop.size, pop.component, mpop)
                )
                self.pop_to_cell_type[pop.id] = pop.component
                chid = moose.copy(self.proto_cells[pop.component], mpop, '%s' % (i))
                self.cells_in_populations[pop.id][i] = chid

    def getInput(self, input_id):
        return moose.element('%s/inputs/%s' % (self.lib.path, input_id))

    def createInputs(self):
        for el in self.network.explicit_inputs:
            pop_id = el.target.split('[')[0]
            i = el.target.split('[')[1].split(']')[0]
            seg_id = 0
            if '/' in el.target:
                seg_id = el.target.split('/')[1]
            input = self.getInput(el.input)
            moose.connect(input, 'output', self.getComp(pop_id, i, seg_id), 'injectMsg')

        for il in self.network.input_lists:
            for ii in il.input:
                input = self.getInput(il.component)
                moose.connect(
                    input,
                    'output',
                    self.getComp(
                        il.populations, ii.get_target_cell_id(), ii.get_segment_id()
                    ),
                    'injectMsg',
                )

    def createCellPrototype(self, cell, symmetric=True):
        """To be completed - create the morphology, channels in prototype"""
        ep = '%s/%s' % (self.lib.path, cell.id)
        nrn = moose.element(ep) if moose.exists(ep) else moose.Neuron(ep)
        self.proto_cells[cell.id] = nrn
        self.nml_cells_to_moose[cell.id] = nrn
        self.moose_to_nml[nrn] = cell
        self.createMorphology(cell, nrn, symmetric=symmetric)
        self.importBiophysics(cell, nrn)
        return cell, nrn

    def createMorphology(self, nmlcell, moosecell, symmetric=True):
        """Create the MOOSE compartmental morphology in `moosecell` using the
        segments in NeuroML2 cell `nmlcell`. Create symmetric
        compartments if `symmetric` is True.

        """
        morphology = nmlcell.morphology
        segments = morphology.segments
        id_to_segment = dict([(seg.id, seg) for seg in segments])
        if symmetric:
            compclass = moose.SymCompartment
        else:
            compclass = moose.Compartment
        # segment names are used as compartment names - assuming
        # naming convention does not clash with that in MOOSE
        cellpath = moosecell.path
        id_to_comp = {}
        self.seg_id_to_comp_name[nmlcell.id] = {}
        for seg in segments:
            if seg.name is not None:
                ep = '%s/%s' % (cellpath, seg.name)
                id_to_comp[seg.id] = (
                    moose.element(ep) if moose.exists(ep) else compclass(ep)
                )
                self.seg_id_to_comp_name[nmlcell.id][seg.id] = seg.name
            else:
                name = 'comp_%s' % seg.id
                ep = '%s/%s' % (cellpath, name)
                id_to_comp[seg.id] = (
                    moose.element(ep) if moose.exists(ep) else compclass(ep)
                )
                self.seg_id_to_comp_name[nmlcell.id][seg.id] = name
        # Now assign the positions and connect up axial resistance
        if not symmetric:
            src, dst = 'axial', 'raxial'
        else:
            src, dst = 'proximal', 'distal'
        for segid, comp in id_to_comp.items():
            segment = id_to_segment[segid]
            try:
                parent = id_to_segment[segment.parent.segments]
            except AttributeError:
                parent = None
            self.moose_to_nml[comp] = segment
            self.nml_segs_to_moose[segment.id] = comp
            p0 = segment.proximal
            if p0 is None:
                if parent:
                    p0 = parent.distal
                else:
                    raise Exception(
                        'No proximal point and no parent segment for segment:'
                        + 'name=%s, id=%s' % (segment.name, segment.id)
                    )
            comp.x0, comp.y0, comp.z0 = (
                x * self.lunit for x in map(float, (p0.x, p0.y, p0.z))
            )
            p1 = segment.distal
            comp.x, comp.y, comp.z = (
                x * self.lunit for x in map(float, (p1.x, p1.y, p1.z))
            )
            comp.length = np.sqrt(
                (comp.x - comp.x0) ** 2
                + (comp.y - comp.y0) ** 2
                + (comp.z - comp.z0) ** 2
            )

            # This can pose problem with moose where both ends of
            # compartment have same diameter. We are averaging the two
            # - may be splitting the compartment into two is better?
            comp.diameter = (float(p0.diameter) + float(p1.diameter)) * self.lunit / 2
            if parent:
                pcomp = id_to_comp[parent.id]
                moose.connect(comp, src, pcomp, dst)
        sg_to_segments = {}
        for sg in morphology.segment_groups:
            sg_to_segments[sg.id] = [id_to_segment[m.segments] for m in sg.members]
        for sg in morphology.segment_groups:
            if not sg.id in sg_to_segments:
                sg_to_segments[sg.id] = []
            for inc in sg.includes:
                for cseg in sg_to_segments[inc.segment_groups]:
                    sg_to_segments[sg.id].append(cseg)

        if not 'all' in sg_to_segments:
            sg_to_segments['all'] = [s for s in segments]

        self._cell_to_sg[nmlcell.id] = sg_to_segments
        return id_to_comp, id_to_segment, sg_to_segments

    def importBiophysics(self, nmlcell, moosecell):
        """Create the biophysical components in moose Neuron `moosecell`
        according to NeuroML2 cell `nmlcell`."""
        bp = nmlcell.biophysical_properties
        if bp is None:
            logger_.info(
                'Warning: %s in %s has no biophysical properties'
                % (nmlcell.id, self.filename)
            )
            return
        self.importMembraneProperties(nmlcell, moosecell, bp.membrane_properties)
        self.importIntracellularProperties(
            nmlcell, moosecell, bp.intracellular_properties
        )

    def importMembraneProperties(self, nmlcell, moosecell, mp):
        """Create the membrane properties from nmlcell in moosecell"""
        if self.verbose:
            logger_.info('Importing membrane properties')
        self.importCapacitances(nmlcell, moosecell, mp.specific_capacitances)
        self.importChannelsToCell(nmlcell, moosecell, mp)
        self.importInitMembPotential(nmlcell, moosecell, mp)

    def importCapacitances(self, nmlcell, moosecell, specificCapacitances):
        sg_to_segments = self._cell_to_sg[nmlcell.id]
        for specific_cm in specificCapacitances:
            cm = SI(specific_cm.value)
            for seg in sg_to_segments[specific_cm.segment_groups]:
                comp = self.nml_segs_to_moose[seg.id]
                comp.Cm = sarea(comp) * cm

    def importInitMembPotential(self, nmlcell, moosecell, membraneProperties):
        sg_to_segments = self._cell_to_sg[nmlcell.id]
        for imp in membraneProperties.init_memb_potentials:
            initv = SI(imp.value)
            for seg in sg_to_segments[imp.segment_groups]:
                comp = self.nml_segs_to_moose[seg.id]
                comp.initVm = initv

    def importIntracellularProperties(self, nmlcell, moosecell, properties):
        self.importAxialResistance(nmlcell, properties)
        self.importSpecies(nmlcell, properties)

    def importSpecies(self, nmlcell, properties):
        sg_to_segments = self._cell_to_sg[nmlcell.id]
        for species in properties.species:
            # Developer note: Not sure if species.concentration_model should be
            # a nml element of just plain string. I was getting plain text from
            # nml file here.
            concModel = species.concentration_model
            if (concModel is not None) and (concModel not in self.proto_pools):
                logger_.warn("No concentrationModel '%s' found." % concModel)
                continue
            segments = getSegments(nmlcell, species, sg_to_segments)
            for seg in segments:
                comp = self.nml_segs_to_moose[seg.id]
                self.copySpecies(species, comp)

    def copySpecies(self, species, compartment):
        """Copy the prototype pool `species` to compartment. Currently only
        decaying pool of Ca2+ supported"""
        proto_pool = None
        concModel = species.concentration_model
        if concModel in self.proto_pools:
            proto_pool = self.proto_pools[concModel]
        else:
            for innerReader in self.includes.values():
                if concModel in innerReader.proto_pools:
                    proto_pool = innerReader.proto_pools[concModel]
                    break
        if not proto_pool:
            msg = 'No prototype pool for %s referred to by %s' % (concModel, species.id)
            logger_.error(msg)
            raise RuntimeError(msg)
        pool_id = moose.copy(proto_pool, compartment, species.id)
        pool = moose.element(pool_id)
        if compartment.length <= 0:
            vol = (
                4
                * np.pi
                * (
                    0.5 * compartment.diameter**3
                    - (0.5 * compartment.diameter - pool.thick) ** 3
                )
                / 3
            )
        else:
            vol = (
                np.pi
                * compartment.length
                * (0.5 * compartment.diameter + pool.thick)
                * (0.5 * compartment.diameter - pool.thick)
            )
        pool.B = pool.B / vol
        return pool

    def importAxialResistance(self, nmlcell, intracellularProperties):
        sg_to_segments = self._cell_to_sg[nmlcell.id]
        for r in intracellularProperties.resistivities:
            segments = getSegments(nmlcell, r, sg_to_segments)
            for seg in segments:
                comp = self.nml_segs_to_moose[seg.id]
                setRa(comp, SI(r.value))

    def isPassiveChan(self, chan):
        if chan.type == 'ionChannelPassive':
            return True
        if hasattr(chan, 'gates'):
            return len(chan.gate_hh_rates) + len(chan.gates) == 0
        return False

    def evaluate_moose_component(self, ct, variables):
        print("[INFO ] Not implemented.")
        return False

    def loopEvalHHRate(self, ct, varname, arrayvars, req_vars=None):
        """Evaluate the parameter `component` for extended HHGate.

        This is used when alpha and beta rates have already been
        calculated using standard HH formulation, but then some
        tweaking is done for calculating tau and/or inf (steady-state)
        value. The formula for tweaking is evaluated using the
        `pynml.evaluate_component` function.

        - When tau is tweaked, the return parameter is `t`. alpha and beta are already calculated and should be passed in `req_vars`.
        - When inf is tweaked, the return parameter is ?
        - When alpha or beta are tweaked, the return parameter is `r`.

        :ct: neuroml2 component type
        :varname: name of the return variable computed by pynml
        :arrayvars: dict array variables. These should have identical
                    length and the same name as in the neuroml expression.
        :req_vars: dict of required variables

        """
        # Evaluate for each value of v, alpha, beta and accumulate the value
        if req_vars is None:
            req_vars = {}
        rate = []
        for k, v in arrayvars.items():
            # extract the common length of array vars
            length = len(v)
            break
        # eval neuroml formula for each index in the array vars
        for ii in range(length):
            req_vars.update({k: v[ii] for k, v in arrayvars.items()})
            val_dict = pynml.evaluate_component(ct, req_variables=req_vars)
            rate.append(val_dict[varname])
        logger_.debug(f'.... rate: {rate}')
        return np.r_[rate]

    def calculateExtendedHHRateFn(self, component, arrayvars, req_vars=None):
        """Compute rate functions for extensions of HHRates."""
        logger_.info(f'AAAAAAA Calculating rate for {component} AAAAAAA')
        for ct in self.doc.ComponentType:
            if ct.name != component.type:
                continue
            # we could do ct.extends.endswith('Time') but explicit is better than implicit: subha
            if ct.extends == 'baseVoltageDepTime':
                varname = 't'
            elif ct.extends == 'baseVoltageDepRate':
                varname = 'r'
            elif ct.extends == 'baseVoltageDepVariable':
                varname = 'x'
            elif ct.extends == 'baseVoltageConcDepRate':
                varname = 'r'
                # Calcium concentration dependent - create a lookup
                # range.  Assuming there is always a vtab that can be
                # referenced for the number of divisions.
                if 'caConc' not in arrayvars:
                    arrayvars['caConc'] = np.linspace(np.finfo(float).smallest_normal, CACONC_MAX, len(arrayvars['v']))
            elif ct.extends == 'baseVoltageConcDepVariable':
                varname = 'x'
                if 'caConc' not in arrayvars:
                    arrayvars['caConc'] = np.linspace(np.finfo(float).smallest_normal, CACONC_MAX, len(arrayvars['v']))
            elif ct.extends == 'baseVoltageConcDepTime':
                varname = 't'
                if 'caConc' not in arrayvars:
                    arrayvars['caConc'] = np.linspace(np.finfo(float).smallest_normal, CACONC_MAX, len(arrayvars['v']))
                                                    
            return self.loopEvalHHRate(ct, varname, arrayvars, req_vars)

    def calculateRateFn(self, ratefn, vtab, vshift='0mV'):
        """Returns A / B table from ngate. If the gate type could not
        be recognized, it returns `None`.

        """
        if ratefn.type in RATE_FN_MAP:
            midpoint, rate, scale = map(
                SI, (ratefn.midpoint, ratefn.rate, ratefn.scale)
            )
            logger_.info(f'{ratefn.type}: Standard HH channel rate function')
            return RATE_FN_MAP[ratefn.type](vtab, rate, scale, midpoint)
        
        req_vars = {'vShift': vshift, 'temperature': self._getTemperature()}
        return self.calculateExtendedHHRateFn(ratefn, {'v': vtab}, req_vars)

    def importChannelsToCell(self, nmlcell, moosecell, membrane_properties):
        sg_to_segments = self._cell_to_sg[nmlcell.id]
        for chdens in (
            membrane_properties.channel_densities
            + membrane_properties.channel_density_v_shifts
        ):
            segments = getSegments(nmlcell, chdens, sg_to_segments)
            condDensity = SI(chdens.cond_density)
            erev = SI(chdens.erev)
            try:
                ionChannel = self.id_to_ionChannel[chdens.ion_channel]
            except KeyError:
                logger_.info('No channel with id: %s' % chdens.ion_channel)
                continue

            if self.verbose:
                logger_.info(
                    'Setting density of channel %s in %s to %s; erev=%s (passive: %s)'
                    % (
                        chdens.id,
                        segments,
                        condDensity,
                        erev,
                        self.isPassiveChan(ionChannel),
                    )
                )

            if self.isPassiveChan(ionChannel):
                for seg in segments:
                    comp = self.nml_segs_to_moose[seg.id]
                    setRm(comp, condDensity)
                    setEk(comp, erev)
            else:
                for seg in segments:
                    self.copyChannel(
                        chdens, self.nml_segs_to_moose[seg.id], condDensity, erev
                    )

    def copyChannel(self, chdens, comp, condDensity, erev):
        """Copy moose prototype for `chdens` condutcance density to `comp`
        compartment.

        """
        proto_chan = None
        if chdens.ion_channel in self.proto_chans:
            proto_chan = self.proto_chans[chdens.ion_channel]
        else:
            for innerReader in self.includes.values():
                if chdens.ionChannel in innerReader.proto_chans:
                    proto_chan = innerReader.proto_chans[chdens.ion_channel]
                    break
        if not proto_chan:
            raise Exception(
                'No prototype channel for %s referred to by %s'
                % (chdens.ion_channel, chdens.id)
            )

        if self.verbose:
            logger_.info(
                'Copying the %s to %s, %s; erev=%s'
                % (chdens.id, comp, condDensity, erev)
            )
        orig = chdens.id
        chid = moose.copy(proto_chan, comp, chdens.id)
        chan = moose.element(chid)
        els = list(self.paths_to_chan_elements.keys())
        for p in els:
            pp = p.replace('%s/' % chdens.ion_channel, '%s/' % orig)
            self.paths_to_chan_elements[pp] = self.paths_to_chan_elements[p].replace(
                '%s/' % chdens.ion_channel, '%s/' % orig
            )
        # logger_.info(self.paths_to_chan_elements)
        chan.Gbar = sarea(comp) * condDensity
        chan.Ek = erev
        moose.connect(chan, 'channel', comp, 'channel')
        return chan

    def getQ10Scale(self, ngate):
        """Compute q10 scale factor for NeuroML gate object `ngate`."""
        q10_scale = 1
        if ngate.q10_settings is None:
            return q10_scale
        
        if ngate.q10_settings.type == 'q10Fixed':
            return float(ngate.q10_settings.fixed_q10)
        elif ngate.q10_settings.type == 'q10ExpTemp':
            return math.pow(
                float(ngate.q10_settings.q10_factor),
                (self._getTemperature() - SI(ngate.q10_settings.experimental_temp))
                / 10,
            )
        else:
            raise Exception(
                'Unknown Q10 scaling type %s: %s'
                % (ngate.q10_settings.type, ngate.q10_settings)
            )

    def createHHChannel(self, chan, vmin=-150e-3, vmax=100e-3, vdivs=5000):
        mchan = moose.HHChannel('%s/%s' % (self.lib.path, chan.id))
        mgates = [moose.element(x) for x in [mchan.gateX, mchan.gateY, mchan.gateZ]]
        assert len(chan.gate_hh_rates) <= 3, "We handle only up to 3 gates in HHCHannel"

        if self.verbose:
            logger_.info(
                '== Creating channel: %s (%s) -> %s (%s)'
                % (chan.id, chan.gate_hh_rates, mchan, mgates)
            )

        all_gates = chan.gates + chan.gate_hh_rates
        # Sort all_gates such that they come in x, y, z order.
        all_gates = _gates_sorted(all_gates)
        # ngate - neuroml gate, mgate - moose gate
        for ngate, mgate in zip(all_gates, mgates):
            logger_.info(f'CCCCCCCCC {chan.id}\n======\n{ngate}\n---------- {mgate}\nCCCCCCCCCC')
            if ngate is None:
                continue
            if mgate.name.endswith('X'):
                mchan.Xpower = ngate.instances
            elif mgate.name.endswith('Y'):
                mchan.Ypower = ngate.instances
            elif mgate.name.endswith('Z'):
                mchan.Zpower = ngate.instances
            mgate.min = vmin
            mgate.max = vmax
            mgate.divs = vdivs

            # I saw only examples of GateHHRates in
            # HH-channels, the meaning of forwardRate and
            # reverseRate and steadyState are not clear in the
            # classes GateHHRatesInf, GateHHRatesTau and in
            # FateHHTauInf the meaning of timeCourse and
            # steady state is not obvious. Is the last one
            # refering to tau_inf and m_inf??
            fwd = ngate.forward_rate
            rev = ngate.reverse_rate
            self.paths_to_chan_elements['%s/%s' % (chan.id, ngate.id)] = '%s/%s' % (
                chan.id,
                mgate.name,
            )

            # Start with standard HHGate of type `gateHHRates`
            # These get tweaked in extended versions
            vtab = np.linspace(vmin, vmax, vdivs + 1)
            alpha = None
            beta = None
            tau = None
            inf = None
            if ngate.forward_rate is not None:
                alpha = self.calculateRateFn(ngate.forward_rate, vtab)
            if ngate.reverse_rate is not None:
                beta = self.calculateRateFn(rev, vtab)
                
            if (alpha is not None) and (beta is not None):
                tau = 1 / (alpha + beta)
                inf = alpha / (alpha + beta)
                
            # Non-standard HH-type gates have a type attribute.
            # The rates/time constants/state variables are tweaked with custom rules.
            if hasattr(ngate, 'type'):
                req_vars = {'temperature': self._getTemperature()}
                if ngate.type == 'gateHHratesTau':
                    assert hasattr(
                        ngate, 'time_course'
                    ), f'{ngate.type} - expected `time_course` attribute'
                    tau = self.calculateExtendedHHRateFn(
                        ngate.time_course,
                        {'v': vtab, 'alpha': alpha, 'beta': beta},
                        req_vars,
                    )
                elif ngate.type == 'gateHHratesInf':
                    assert hasattr(
                        ngate, 'steady_state'
                    ), f'{ngate.type} - expected `steady_state` attribute'
                    inf = self.calculateExtendedHHRateFn(
                        ngate.steady_state,
                        {'v': vtab, 'alpha': alpha, 'beta': beta},
                        req_vars,
                    )
                elif ngate.type == 'gateHHratesTauInf':
                    assert hasattr(ngate, 'steady_state') and hasattr(
                        ngate, 'time_course'
                    ), f'{ngate.type} - expected `steady_state` and `time_course` attributes'                    
                    inf = self.calculateExtendedHHRateFn(
                        ngate.steady_state,
                        {'v': vtab, 'alpha': alpha, 'beta': beta},
                        req_vars,
                    )
                    tau = self.calculateExtendedHHRateFn(
                        ngate.time_course,
                        {'v': vtab, 'alpha': alpha, 'beta': beta},
                        req_vars,
                    )

                # logger_.info(f'inf="{inf}", tau="{tau}"')
            q10_scale = self.getQ10Scale(ngate)
            tau /= q10_scale
            mgate.tableA = inf / tau
            mgate.tableB = 1 / tau

            # DEBUG: BEGIN
        #     import matplotlib.pyplot as plt

        #     fig, axes = plt.subplots(nrows=2)
        #     fig.suptitle(f'{chan.id} {ngate.id} {ngate.type}')
        #     axes[0].plot(vtab, inf)
        #     axes[1].plot(vtab, tau)
        #     plt.show()

        # logger_.info(f'{"$" * 20} {self.filename}: Created {mchan.path} for {chan.id}')

        # DEBUG: END
        return mchan

    def createPassiveChannel(self, chan):
        epath = '%s/%s' % (self.lib.path, chan.id)
        if moose.exists(epath):
            mchan = moose.element(epath)
        else:
            mchan = moose.Leakage(epath)
        logger_.info('%s: Created %s for %s' % (self.filename, mchan.path, chan.id))
        return mchan

    def importInputs(self, doc):
        epath = '%s/inputs' % (self.lib.path)
        if moose.exists(epath):
            minputs = moose.element(epath)
        else:
            minputs = moose.Neutral(epath)

        for pg_nml in doc.pulse_generators:
            epath = '%s/%s' % (minputs.path, pg_nml.id)
            pg = moose.element(epath) if moose.exists(epath) else moose.PulseGen(epath)
            pg.firstDelay = SI(pg_nml.delay)
            pg.firstWidth = SI(pg_nml.duration)
            pg.firstLevel = SI(pg_nml.amplitude)
            pg.secondDelay = 1e9
            logger_.debug(f'{"$" * 10} Created input {epath}')

    def importIonChannels(self, doc, vmin=-150e-3, vmax=100e-3, vdivs=5000):
        logger_.info('%s: Importing the ion channels' % self.filename)

        for chan in doc.ion_channel + doc.ion_channel_hhs:
            if chan.type == 'ionChannelHH':
                mchan = self.createHHChannel(chan)
            elif self.isPassiveChan(chan):
                mchan = self.createPassiveChannel(chan)
            else:
                mchan = self.createHHChannel(chan)

            self.id_to_ionChannel[chan.id] = chan
            self.nml_chans_to_moose[chan.id] = mchan
            self.proto_chans[chan.id] = mchan
            logger_.info(
                self.filename
                + ': Created ion channel %s for %s %s'
                % (mchan.path, chan.type, chan.id)
            )

    def importConcentrationModels(self, doc):
        for concModel in doc.decaying_pool_concentration_models:
            self.createDecayingPoolConcentrationModel(concModel)

    def createDecayingPoolConcentrationModel(self, concModel):
        """Create prototype for concentration model"""
        if hasattr(concModel, 'name') and concModel.name is not None:
            name = concModel.name
        else:
            name = concModel.id

        ca = moose.CaConc('%s/%s' % (self.lib.path, name))
        ca.CaBasal = SI(concModel.resting_conc)
        ca.tau = SI(concModel.decay_constant)
        ca.thick = SI(concModel.shell_thickness)
        ca.B = 5.2e-6  # B = 5.2e-6/(Ad) where A is the area of the
        # shell and d is thickness - must divide by
        # shell volume when copying
        self.proto_pools[concModel.id] = ca
        self.nml_conc_to_moose[concModel.id] = ca
        self.moose_to_nml[ca] = concModel
        logger_.debug(
            'Created moose element: %s for nml conc %s' % (ca.path, concModel.id)
        )
