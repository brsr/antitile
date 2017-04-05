# -*- coding: utf-8 -*-
"""
Creates SVG files for each subdivision square
"""

import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from math import gcd
import antitile
from antitile import xmath


def square_explicit(frame, line_pts_n=61, degenerate=False):
    if ~degenerate:
        frame = frame.reshape((-1,2,2))
        norm = np.linalg.norm(frame[:,0]-frame[:,1], axis=-1)
        bad = np.isclose(norm, 0)
        frame = frame[~bad]
    t = np.linspace(0, 1, line_pts_n)[:, np.newaxis]
    frame_0 = frame[..., 0, np.newaxis, :]
    frame_1 = frame[..., 1, np.newaxis, :]
    return xmath.lerp(frame_0, frame_1, t)


def plot_square_circle(n=4, m=2, inverse=True, steplimit=90):
    if inverse:
        limits = [-2, 2]
    else:
        limits = [-1.05, 1.05]
    frame = antitile.breakdown.frame_square(n, m)
    euc_lines = square_explicit(frame)
    lines = antitile.projection.square_to_circle(euc_lines)
    fig, ax = plt.subplots()
    fig.set_size_inches(8, 8)
    fig.set_tight_layout(True)
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.axis('off')#why do we need both? it is a mystery
    if inverse:
        norm = np.linalg.norm(lines, axis=-1, keepdims=True)
        inverse = lines / norm**2
        inverse[..., 0] = -inverse[..., 0]
        dx = np.linalg.norm(inverse[:, 1:] - inverse[:, :-1], axis=-1)
        if steplimit:
            index = dx > steplimit
            pad = np.zeros((index.shape[0], 1), dtype=bool)
            index = np.concatenate((index, pad), axis=-1)
            print(index.shape, inverse.shape)
            inverse = np.where(index[..., np.newaxis], np.nan, inverse)
        lc = LineCollection(inverse, color='grey', zorder=0)
        ax.add_collection(lc)
    lc = LineCollection(lines, color='k', zorder=1)
    ax.add_collection(lc)
    if m > 0:
        color = 'blue'
    else:
        color = 'black'
    circle = plt.Circle((0, 0), 1, color=color, fill=False, zorder=3)
    ax.add_artist(circle)
    #plot_m1_pts(ax, base_pts, method, n, m)
    return fig, ax

def lines(n, m, digits, poly, style='stroke="black"'):
    shape = poly['shape']
    frame_fun = poly['frame']
    statements = []
    lineitem = ('<line x1="{x1:g}" y1="{y1:g}" x2="{x2:g}" y2="{y2:g}" ' +
        style + '/>')
    frame = frame_fun(n, m)
    frame = frame.reshape((-1,2,2))
    if m > 0:#remove degenerate lines
        norm = np.linalg.norm(frame[:,0]-frame[:,1], axis=-1)
        bad = np.isclose(norm, 0)
    else: #remove lines that overlap the black triangle (and degenerates)
        norm = np.linalg.norm(frame[..., np.newaxis, :] - shape, axis=-1)
        bad = np.any(np.isclose(norm, 0), axis=(1,2))
    frame = np.round(frame[~bad], digits)
    for line in frame:
        pt1 = line[0]
        pt2 = line[1]
        if np.all(np.isclose(pt1, pt2)):
            pass
        else:
            thisline = lineitem.format(x1=pt1[0], x2=pt2[0],
                                       y1=pt1[1], y2=pt2[1])
            statements.append(thisline)
    return statements

def svg_breakdowns(shapeinfo, n, m):
    g = gcd(n, m)
    digits = 1
    size = shapeinfo['size']
    shape = shapeinfo['shape']
    header1 = '<?xml version="1.0" encoding="utf-8" standalone="no"?>'
    header2 = ('<svg height="{h}" width="{w}" '.format(w=size[0], h=size[1]) +
        'xmlns="http://www.w3.org/2000/svg" ' +
        'xmlns:xlink="http://www.w3.org/1999/xlink">')
    footer1 = ('<polygon points="' +
            ' '.join([format(pt[0], 'g') + ',' + format(pt[1], 'g')
                for pt in np.round(shape, digits)]) +
            '" fill="none" stroke=')
    if m > 0:
        footer1 += '"blue"/>'
    else:
        footer1 += '"black"/>'
    footer2 = '</svg>'
    statements = [header1, header2]
    if g > 1 and m > 0:
        style = 'stroke="red" stroke-dasharray="3"'
        statements += lines(g, 0, digits, shapeinfo, style)
    statements +=  lines(n, m, digits, shapeinfo)
    statements += [footer1, footer2]
    return '\n'.join(statements)

square = {'name': 'square',
          'shape': np.array([[1,  1],
                             [181, 1],
                             [181, 181],
                             [1,  181]]),
          'size': [182,]*2,
          'frame': lambda n, m: antitile.breakdown.frame_square(n, m)*180 + 1,
          'diagonal': 'Diagonal subdivision'}

tri_shape = np.array([[1, 157],
                      [181, 157],
                      [91, 1.1154273188010393]]
        )
triangle = {'name': 'triangle',
            'shape': tri_shape,
            'size': [182, 158],
            'frame': lambda n, m: antitile.breakdown.frame_triangle(
                                        tri_shape, n, m, interp=xmath.lerp),
            'diagonal': 'Ortho subdivision'}

def factors(number):
    """Cheezy function to get the factors of a small number
    other than 1 and itself"""
    for i in range(2,number):
        if number % i == 0:
            yield i

def create_svgs(shapeinfo, end = 16):
    """Produces a bunch of SVG files and a csv file for use with pattypan"""
    shapename = shapeinfo['name']
    desc_end = (". See [[:Category:Subdivision triangles; transparent]] for " +
               "details. Please do not edit this file without reading the " +
               "category page.")
    name_template = 'Subdivided ' +shapename+' {:02d} {:02d}'
    wikifile = 'File:' + name_template + '.svg'

    with open(shapeinfo['name'] + '.csv', mode='w', newline='') as c:
        wr = csv.writer(c)
        wr.writerow(['path','name','description','other_versions','category'])
        for i in range(1, end+1):
            for j in range(end+1):
                g = gcd(i, j)
                svg = svg_breakdowns(shapeinfo, i, j)
                name = name_template.format(i, j)
                filename = (name +'.svg')

                with open(filename, 'w') as f:
                    f.write(svg)

                gallery = ''
                if j == 0:
                    category = 'Parallel subdivision '+shapename+'s; transparent'
                elif j == i:
                    category = shapeinfo['diagonal']+' '+shapename+'s; transparent'
                else:
                    category = 'Skew subdivision '+shapename+'s; transparent'
                    gallery += (wikifile.format(j, i) +
                                      '|Mirror image ')
                if g > 1 and j > 0:
                    gallery += (wikifile.format(g, 0) +
                                      '|Red lines indicate a composition of this...')
                    gallery += (wikifile.format(i//g, j//g) +
                                      '|...with this.')
                elif j == 0:
                    for n in factors(i):
                        gallery += (wikifile.format(n, 0) +
                                           '|Contained subdivision')
                if len(gallery) > 0:
                    other_versions = '<gallery>' + gallery + '</gallery>'
                else:
                    other_versions = ''
                desc = name + desc_end
                wr.writerow([filename, filename, desc,
                             other_versions, category])

def create_all(end=16):
    for shapeinfo in [triangle, square]:
        create_svgs(shapeinfo, end)
