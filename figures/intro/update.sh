#!/bin/sh

# integrals.py creates rectange.pdf trapezoid.pdf simpsons.pdf
python3 integrals.py 
cp -f rectangle.pdf trapezoid.pdf simpsons.pdf ../../intro/

# rk4_plot.py creates rk4_k[1-4].pdf, rk4_final.pdf
python3 rk4_plot.py
cp -f rk4_k[1-4].pdf rk4_final.pdf ../../intro

