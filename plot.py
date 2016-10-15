#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np
import matplotlib.pyplot as plt

LX = 50
LY = 50
MAX_TIMESTEP = 100

INPUT_PATH = 'data/'
OUTPUT_PATH = 'plot/'

origin = 'lower'
#origin = 'upper'


def input_filename(prefix, ts):
    return INPUT_PATH + prefix + '%06d' % ts + '.txt'

def output_filename(prefix, ts):
    return OUTPUT_PATH + prefix + '%06d' % ts + '.png'

def setup_data(filename, index):
    return np.loadtxt(filename, delimiter=' \t ')[:, index].reshape(LX, LY)


def main():
    argvs = sys.argv
    argc = len(argvs)

    if (argc != 3):
        quit()

    prefix = argvs[1]
    index  = int(argvs[2])
    
    for ts in range(MAX_TIMESTEP + 1):
        print input_filename(prefix, ts)
        print index

        v = setup_data(input_filename(prefix, ts), index)
        vmax = 0.05#np.max(v)
        vmin = -vmax

        plt.imshow(v, origin = 'lower', interpolation = 'none', cmap=plt.cm.jet,vmin=vmin,vmax=vmax)
        plt.colorbar()

        plt.savefig(output_filename(prefix, ts))

        plt.clf()
    
if __name__ == '__main__':
    main()
