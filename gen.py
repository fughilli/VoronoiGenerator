#!/usr/bin/python

import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
import argparse
import sys
from math import sin,cos,pi,sqrt,log
from PIL import Image, ImageDraw

# Defaults
geom = (2160, 2160)
numcenters = int(np.random.random((1,))[0] * 4000) + 100

# Argparse
parser = argparse.ArgumentParser(description="Generate a fractalized gradient.")
parser.add_argument('-g', '--geometry', type=int, nargs=2, help='A pair of integers specifying the dimensions of the output image.')

# Parse arguments
args = parser.parse_args(sys.argv[1:])

if args.geometry:
    geom = tuple(args.geometry)

def rand_in_rect(aspect, count):
    ret = []
    while(len(ret) < count):
        candidates = np.random.random((count,2))
        ret.extend(filter(lambda x : x[1] < aspect, candidates))
    return ret[:count]

colorf_coeffs = np.random.random((3,5)) - 0.1

def colorf(vertex):
    xyc = np.concatenate((vertex, (0.5,cos(vertex[0]),sin(vertex[1]))), 0)
    return np.clip(
            (
                np.dot(colorf_coeffs[0], xyc)/3,
                np.dot(colorf_coeffs[1], xyc)/3,
                np.dot(colorf_coeffs[2], xyc)/3
            ), 0, 1) + (np.random.random((1,))[0] * 0.05 - 0.025)

aspect = float(geom[1])/geom[0]

vertices = rand_in_rect(aspect, numcenters)
vor = Voronoi(vertices)
colors = map(colorf, vertices)
colors_and_regions = zip(colors, map(lambda i : vor.regions[i], vor.point_region))
complete_colors_and_regions = filter(lambda cr : -1 not in cr[1], colors_and_regions)
colors_and_vertices = map(lambda cr : (cr[0], map(lambda i : vor.vertices[i], cr[1])), complete_colors_and_regions)
z = np.random.random(len(colors_and_vertices))

im = Image.new("RGB", geom, "white")
draw = ImageDraw.Draw(im)

vor_verts = reduce(lambda a,b : a+b, map(lambda cv : cv[1], colors_and_vertices))

vor_xs = map(lambda xy : xy[0], vor_verts)
vor_ys = map(lambda xy : xy[1], vor_verts)

max_left = 0.2
min_right = 0.8
max_bottom = 0.2
min_top = 0.8

def linmap(v, im, ix, om, ox):
    return ((ox-om)*(v-im)/(ix - im)) + om

# make a list of polygons
polys = []
fill_colors = []
outline_colors = []
for color,vertices in colors_and_vertices:
    scaled_color = tuple(int(255 * x) for x in color)
    darkened_scaled_color = tuple(int(x * 0.9) for x in scaled_color)
    scaled_vertices = map(lambda v : (im.size[0] * linmap(v[0], max_left, min_right, 0.0, 1.0), im.size[1] * linmap(v[1]/aspect, max_bottom, min_top, 0.0, 1.0)), vertices)
    polys.append(scaled_vertices)
    fill_colors.append(scaled_color)
    outline_colors.append(darkened_scaled_color)

for poly,fill in zip(polys, fill_colors):
    draw.polygon(poly, fill=fill)

for poly,outline in zip(polys, outline_colors):
    draw.polygon(poly, outline=outline)

del draw

im.save("output.png")
