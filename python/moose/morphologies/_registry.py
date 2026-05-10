"""
moose.morphologies._registry
=============================
Bundled SWC morphology files shipped with pymoose.

Each entry dict may contain:
  name        str   Short identifier used in moose.morphologies.load('name', ...)
  filename    str   SWC filename inside moose/morphologies/data/swc/
  species     str   Species (e.g. 'rat', 'mouse', 'human')
  cell_type   str   Cell type (e.g. 'CA1 pyramidal', 'L5 pyramidal')
  region      str   Brain region
  source      str   Original database / publication
  source_url  str   URL to original data record
  license     str   License or terms of use
  description str   One-line description
"""

from moose._registry_base import Registry

_registry = Registry('Morphology', 'moose.morphologies.load()')

_registry.add([
    {
        'name':        'traub91_CA1',
        'filename':    'swc/traub91_CA1.swc',
        'species':     'rat',
        'cell_type':   'CA1 pyramidal',
        'region':      'hippocampus CA1',
        'source':      'Traub et al. 1991, J Neurophysiol 66:635-650',
        'description': '19-compartment CA1 pyramidal (Traub et al. 1991)',
    },
    {
        'name':        'traub91_CA3',
        'filename':    'swc/traub91_CA3.swc',
        'species':     'rat',
        'cell_type':   'CA3 pyramidal',
        'region':      'hippocampus CA3',
        'source':      'Traub et al. 1991, J Neurophysiol 66:635-650',
        'description': '19-compartment CA3 pyramidal (Traub et al. 1991)',
    },
    {
        'name':        'purk_eds1994_full',
        'filename':    'swc/purk_eds1994_full.swc',
        'species':     'guinea pig',
        'cell_type':   'Purkinje',
        'region':      'cerebellum',
        'source':      'De Schutter & Bower 1994, J Neurophysiol 71:375-400',
        'description': '1600-compartment Purkinje cell (De Schutter & Bower 1994)',
    },
    {
        'name':        'gran_bhalla1991_ob',
        'filename':    'swc/gran_bhalla1991_ob.swc',
        'species':     'rat',
        'cell_type':   'granule',
        'region':      'olfactory bulb',
        'source':      'Bhalla & Bower 1993, J Neurophysiol 69:1948-1965',
        'description': '112-compartment olfactory bulb granule cell (Bhalla 1991)',
    },
    {
        'name':        'mit_bhalla1991',
        'filename':    'swc/mit_bhalla1991.swc',
        'species':     'rat',
        'cell_type':   'mitral',
        'region':      'olfactory bulb',
        'source':      'Bhalla & Bower 1993, J Neurophysiol 69:1948-1965',
        'description': '286-compartment detailed mitral cell (Bhalla 1991)',
    },
    {
        'name':        'mit_davison_reduced',
        'filename':    'swc/mit_davison_reduced.swc',
        'species':     'rat',
        'cell_type':   'mitral',
        'region':      'olfactory bulb',
        'source':      'Davison et al. 2003, J Comput Neurosci 15:375-384',
        'description': '7-compartment reduced mitral cell (Davison/Bhalla/Bower)',
    },
    {
        'name':        'gran_migliore_olfactory',
        'filename':    'swc/gran_migliore_olfactory.swc',
        'species':     'rat',
        'cell_type':   'granule',
        'region':      'olfactory bulb',
        'source':      'Migliore & Shepherd 2008, PLoS Comput Biol 4:e1000011',
        'description': '2-compartment olfactory bulb granule cell (Migliore & Shepherd 2008)',
    },

    # ── Allen Cell Types Database morphologies ────────────────────────────────
    # License: Allen Institute Terms of Use (non-commercial research use only)
    # https://alleninstitute.org/legal/terms-use/
    {
        'name':        'allen_mouse_VISp_L5_485909730',
        'filename':    'swc/allen_mouse_VISp_L5_485909730.swc',
        'species':     'mouse',
        'cell_type':   'L5 pyramidal',
        'region':      'primary visual cortex (VISp) layer 5',
        'source':      'Allen Institute for Brain Science, Cell Types Database, specimen 485909730',
        'source_url':  'https://celltypes.brain-map.org/experiment/morphology/485909730',
        'license':     'Allen Institute Terms of Use (non-commercial research use): https://alleninstitute.org/legal/terms-use/',
        'description': 'Full reconstruction of mouse VISp L5 pyramidal neuron (Allen Cell Types DB #485909730)',
    },
    {
        'name':        'allen_mouse_VISp_L5_515249852',
        'filename':    'swc/allen_mouse_VISp_L5_515249852.swc',
        'species':     'mouse',
        'cell_type':   'L5 pyramidal',
        'region':      'primary visual cortex (VISp) layer 5',
        'source':      'Allen Institute for Brain Science, Cell Types Database, specimen 515249852',
        'source_url':  'https://celltypes.brain-map.org/experiment/morphology/515249852',
        'license':     'Allen Institute Terms of Use (non-commercial research use): https://alleninstitute.org/legal/terms-use/',
        'description': 'Full reconstruction of mouse VISp L5 pyramidal neuron (Allen Cell Types DB #515249852)',
    },
    {
        'name':        'allen_mouse_VISp_L6a_580007431',
        'filename':    'swc/allen_mouse_VISp_L6a_580007431.swc',
        'species':     'mouse',
        'cell_type':   'L6a pyramidal',
        'region':      'primary visual cortex (VISp) layer 6a',
        'source':      'Allen Institute for Brain Science, Cell Types Database, specimen 580007431',
        'source_url':  'https://celltypes.brain-map.org/experiment/morphology/580007431',
        'license':     'Allen Institute Terms of Use (non-commercial research use): https://alleninstitute.org/legal/terms-use/',
        'description': 'Full reconstruction of mouse VISp L6a pyramidal neuron (Allen Cell Types DB #580007431)',
    },
    {
        'name':        'allen_mouse_VISp_L6b_589128331',
        'filename':    'swc/allen_mouse_VISp_L6b_589128331.swc',
        'species':     'mouse',
        'cell_type':   'L6b pyramidal',
        'region':      'primary visual cortex (VISp) layer 6b',
        'source':      'Allen Institute for Brain Science, Cell Types Database, specimen 589128331',
        'source_url':  'https://celltypes.brain-map.org/experiment/morphology/589128331',
        'license':     'Allen Institute Terms of Use (non-commercial research use): https://alleninstitute.org/legal/terms-use/',
        'description': 'Full reconstruction of mouse VISp L6b pyramidal neuron (Allen Cell Types DB #589128331)',
    },
    {
        'name':        'allen_human_MTG_L2_616647103',
        'filename':    'swc/allen_human_MTG_L2_616647103.swc',
        'species':     'human',
        'cell_type':   'L2 pyramidal',
        'region':      'middle temporal gyrus (MTG) layer 2',
        'source':      'Allen Institute for Brain Science, Cell Types Database, specimen 616647103',
        'source_url':  'https://celltypes.brain-map.org/experiment/morphology/616647103',
        'license':     'Allen Institute Terms of Use (non-commercial research use): https://alleninstitute.org/legal/terms-use/',
        'description': 'Full reconstruction of human MTG L2 pyramidal neuron (Allen Cell Types DB #616647103)',
    },
    {
        'name':        'allen_human_MTG_L6_614635228',
        'filename':    'swc/allen_human_MTG_L6_614635228.swc',
        'species':     'human',
        'cell_type':   'L6 pyramidal',
        'region':      'middle temporal gyrus (MTG) layer 6',
        'source':      'Allen Institute for Brain Science, Cell Types Database, specimen 614635228',
        'source_url':  'https://celltypes.brain-map.org/experiment/morphology/614635228',
        'license':     'Allen Institute Terms of Use (non-commercial research use): https://alleninstitute.org/legal/terms-use/',
        'description': 'Full reconstruction of human MTG L6 pyramidal neuron (Allen Cell Types DB #614635228)',
    },
])

get         = _registry.get
all_entries = _registry.all_entries
