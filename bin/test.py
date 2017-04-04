# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 15:43:02 2017

@author: Bstone
"""

import antitile
import numpy as np

def thing(co, n, freq):
    a, b = freq
    cco = antitile.xmath.float2d_to_complex(co)
    rot_cco = cco * np.exp(1j*np.pi*n/2)
    if n == 0:
        offset = 0
    elif n == 1:
        offset = 1
    elif n == 2:
        offset = 1+1j
    elif n == 3:
        offset = 1j
    shift_cco = rot_cco + offset    
    cxy = shift_cco * (a + b*1j) + b        
    return antitile.xmath.complex_to_float2d(cxy)

def rot_4(coord, n, freq):
    a, b = freq
    coord = coord.astype(float)
    cco = antitile.xmath.float2d_to_complex(coord).flatten()
    rot_cco = cco * np.exp(1j*np.pi*n/2)
    if n == 0:
        offset = 0
    elif n == 1:
        offset = 1
    elif n == 2:
        offset = 1+1j
    elif n == 3:
        offset = 1j
    shift_cco = rot_cco + offset
    cxy = shift_cco * (a + b*1j) + b
    lindex = antitile.xmath.complex_to_float2d(cxy)
    return lindex.astype(int)

freq = 4, 2
bkdn = antitile.breakdown.Breakdown(*freq,4)
cx = bkdn.coord

cx0 = thing(cx, 0, freq)
cx1 = thing(cx, 1, freq)
cx2 = thing(cx, 2, freq)
cx3 = thing(cx, 3, freq)

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(2,2)
    fig.set_size_inches(6, 6)
    plt.axis('equal')
    ax0.scatter(cx0[..., 0], cx0[..., 1])
    ax1.scatter(cx1[..., 0], cx1[..., 1])
    ax2.scatter(cx2[..., 0], cx2[..., 1])
    ax3.scatter(cx3[..., 0], cx3[..., 1])
    