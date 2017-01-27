#!/bin/sh

DEST=../../Euler

# riemann-states.py makes riemann_comp.pdf
python3 riemann-states.py
cp -f riemann_comp.pdf ${DEST}


# ppm.py makes piecewise-constant.pdf piecewise-linear.pdf piecewise-parabolic.pdf
python3 ppm.py
cp -f piecewise-constant.pdf piecewise-linear.pdf piecewise-parabolic.pdf ${DEST}


# states.py makes states.pdf
python3 states.py
cp -f states.pdf ${DEST}


