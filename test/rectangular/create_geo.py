#/!usr/bin/env python3

import numpy as np

alpha = 0.25*np.pi # contact angle
width = 500e-6 # channel width
height = 500e-6 # triple contact point
sample_points = 100

x = np.linspace(-width*0.5, width*0.5, num=sample_points, endpoint=True)
y = height + 0.5*width*np.tan(alpha) - np.sqrt(width**2/(4*np.cos(alpha)**2) - x**2)

point_index = 0
line_index = 0
with open("channel.geo", "w") as geo:
    geo.write(f"cl_ = {0.01*width};\n")
    point_index += 1
    geo.write(f"Point({point_index}) = {{{0.5*width}, 0.0, 0, cl_}};\n")
    point_index += 1
    geo.write(f"Point({point_index}) = {{{-0.5*width}, 0.0, 0, cl_}};\n")
    line_index += 1
    geo.write(f"Line({line_index}) = {{1, 2}};\n")

    for xx, yy in zip(x, y):
        point_index += 1
        geo.write(f"Point({point_index}) = {{{xx}, {yy}, 0, cl_}};\n")
        line_index += 1
        geo.write(f"Line({line_index}) = {{{point_index-1}, {point_index}}};\n")

    line_index += 1
    geo.write(f"Line({line_index}) = {{{point_index}, 1}};\n")

    line_list = ",".join([str(i) for i in range(1, line_index+1)]).rstrip(",")
    geo.write(f"Curve Loop(1) = {{{line_list}}};\n")
    geo.write(f"Plane Surface(1) = {{1}};\n")
