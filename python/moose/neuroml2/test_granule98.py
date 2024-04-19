# test_granule98.py --- 
# 
# Filename: test_granule98.py
# Description: 
# Author: Subhasis Ray
# Created: Mon Apr  8 21:41:22 2024 (+0530)
# Last-Updated: Fri Apr 19 12:01:39 2024 (+0530)
#           By: Subhasis Ray
# 

# Code:
"""Test code for the Granule cell model

"""
import os
import sys
import numpy as np
# import unittest
import logging
import argparse

LOGLEVEL = os.environ.get('LOGLEVEL', 'INFO')
logging.basicConfig(level=LOGLEVEL)


import moose
from moose.neuroml2.reader import NML2Reader


def run(modeldir, nogui=True, refdir=None):
    reader = NML2Reader()
    filename = os.path.join(modeldir, 'GranuleCell.net.nml')
    reader.read(filename)
    soma = reader.getComp(reader.doc.networks[0].populations[0].id, 0, 0)
    data = moose.Neutral('/data')
    pg = reader.getInput('Gran_10pA')
    inj = moose.Table(f'{data.path}/pulse')
    moose.connect(inj, 'requestOut', pg, 'getOutputValue')
    vm = moose.Table(f'{data.path}/Vm')
    moose.connect(vm, 'requestOut', soma, 'getVm')
    naf_m = moose.Table(f'{data.path}/NaF_m')
    moose.connect(naf_m, 'requestOut', f'{soma.path}/Gran_NaF_98_all', 'getX')
    naf_h = moose.Table(f'{data.path}/NaF_h')
    moose.connect(naf_h, 'requestOut', f'{soma.path}/Gran_NaF_98_all', 'getY')
    # kdr_h = moose.Table(f'{data.path}/KDr_h')
    # moose.connect(kdr_h, 'requestOut', f'{soma.path}/Gran_KDr_98_all', 'getX')
    # kca_m = moose.Table(f'{data.path}/KCa_m')
    # moose.connect(kca_m, 'requestOut', f'{soma.path}/Gran_KCa_98_all', 'getX')
    # ka_m = moose.Table(f'{data.path}/KA_m')
    # moose.connect(ka_m, 'requestOut', f'{soma.path}/Gran_KA_98_all', 'getX')
    # ka_h = moose.Table(f'{data.path}/KA_h')
    # moose.connect(ka_h, 'requestOut', f'{soma.path}/Gran_KA_98_all', 'getY')
    # h_n = moose.Table(f'{data.path}/H_n')
    # moose.connect(h_n, 'requestOut', f'{soma.path}/Gran_H_98_all', 'getX')
    # cahva_m = moose.Table(f'{data.path}/CaHVA_m')
    # moose.connect(cahva_m, 'requestOut', f'{soma.path}/Gran_CaHVA_98_all', 'getX')
    # cahva_h = moose.Table(f'{data.path}/CaHVA_h')
    # moose.connect(cahva_h, 'requestOut', f'{soma.path}/Gran_CaHVA_98_all', 'getY')
    # caconc = moose.Table(f'{data.path}/CaPool')
    # moose.connect(caconc, 'requestOut', f'{soma.path}/ca', 'getCa')
    
    
    
    data_files = {'Gran_0': vm,
                  'Gran_0.Gran_NaF_98_m': naf_m,
                  'Gran_0.Gran_NaF_98_h': naf_h,
                  # 'Gran_0.Gran_KDr_98_h': kdr_h,
                  # 'Gran_0.Gran_KCa_98_m': kca_m,
                  # 'Gran_0.Gran_KA_98_m': ka_m,
                  # 'Gran_0.Gran_KA_98_h': ka_h,
                  # 'Gran_0.Gran_H_98_n': h_n,
                  # 'Gran_0.Gran_CaHVA_98_m': cahva_m,
                  # 'Gran_0.Gran_CaHVA_98_h': cahva_h,
                  # 'Gran_0.Gran_CaPool_98_CONC_ca': caconc
                  }
                  
                  

                  
    simtime = pg.delay[0] + pg.width[0] + 0.1
    moose.reinit()
    moose.start(simtime)

    for fname, tab in data_files.items():
        t = np.arange(len(tab.vector)) * tab.dt
        results = np.c_[t, tab.vector]
        np.savetxt(f'{fname}.dat', results)
        print(f'Saved data in {fname}.dat')
    if not nogui:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(nrows=len(data_files), sharex='all')
        for ax, fname in zip(axes, data_files):
            tab = data_files[fname]
            t = np.arange(len(tab.vector)) * tab.dt
            ax.plot(t, tab.vector, label=f'{fname} moose')
            if refdir is not None:
                refdata = np.loadtxt(os.path.join(refdir, f'{fname}.dat'))
                ax.plot(refdata[:, 0], refdata[:, 1], '--', label='reference')
            ax.legend()
        fig, axes = plt.subplots(nrows=1, ncols=2)
        naf_x = moose.element(f'{soma.path}/Gran_NaF_98_all/gateX')
        naf_y = moose.element(f'{soma.path}/Gran_NaF_98_all/gateY')
        v = np.linspace(naf_x.min, naf_x.max, naf_x.divs + 1)
        minf = naf_x.tableA / naf_x.tableB
        mtau = 1 / naf_x.tableB
        hinf = naf_y.tableA / naf_y.tableB
        htau = 1 / naf_y.tableB
        print(f'taum min: {min(mtau)} max {max(mtau)}')
        print(f'tauh min: {min(htau)} max {max(htau)}')
        print(f'minf min: {min(minf)} max {max(minf)}')
        print(f'hinf min: {min(hinf)} max {max(hinf)}')
        np.savetxt('Gran_NaF_98.minf.txt', np.c_[v, minf])
        np.savetxt('Gran_NaF_98.hinf.txt', np.c_[v, hinf])
        np.savetxt('Gran_NaF_98.mtau.txt', np.c_[v, mtau])
        np.savetxt('Gran_NaF_98.htau.txt', np.c_[v, htau])
        axes[0].plot(v, minf, label='moose minf')
        axes[0].plot(v, hinf, label='moose hinf')
        axes[1].plot(v, mtau, label='moose mtau')
        axes[1].plot(v, htau, label='moose htau')
        # if refdir is not None:
        #     minf_lems = np.loadtxt(os.path.join(refdir, 'Gran_NaF_98.m.inf.lems.dat'))
        #     hinf_lems = np.loadtxt(os.path.join(refdir, 'Gran_NaF_98.h.inf.lems.dat'))
        #     mtau_lems = np.loadtxt(os.path.join(refdir, 'Gran_NaF_98.m.tau.lems.dat'))
        #     htau_lems = np.loadtxt(os.path.join(refdir, 'Gran_NaF_98.h.tau.lems.dat'))
        #     axes[0].plot(minf_lems[:, 0], minf_lems[:, 1], label='lems minf')
        #     axes[0].plot(hinf_lems[:, 0], hinf_lems[:, 1], label='lems hinf')
        #     axes[1].plot(mtau_lems[:, 0], mtau_lems[:, 1], label='lems mtau')
        #     axes[1].plot(htau_lems[:, 0], htau_lems[:, 1], label='lems htau')
        axes[0].legend()
        axes[1].legend()
        plt.show()
        return reader
        


        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run Maex and De Schutter 98 Granule cell model')
    parser.add_argument('modeldir', type=str, default='.', help='directory that contains the NeuroML2 model')
    parser.add_argument('-r', '--refdir', type=str, default='', help='directory that contains reference data files')
    parser.add_argument('-n', '--nogui', help='disable gui', action='store_true')
    
    args = parser.parse_args(sys.argv[1:])
    if args.refdir == '':
        args.refdir = args.modeldir
        
    reader = run(args.modeldir, nogui=args.nogui, refdir=args.refdir)


# 
# test_granule98.py ends here
