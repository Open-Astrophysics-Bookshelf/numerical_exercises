#!/bin/sh

# derivatives.py creates derivs.pdf
python3 derivatives.py
cp -f derivs.pdf ../../intro/

# deriv_error.py creates deriv_error.pdf
python3 deriv_error.py
cp -f deriv_error.pdf ../../intro/

# integrals.py creates rectange.pdf trapezoid.pdf simpsons.pdf
python3 integrals.py 
cp -f rectangle.pdf trapezoid.pdf simpsons.pdf ../../intro/

# roots_plot.py creates newton_0[0-3].pdf
python3 roots_plot.py
cp -f newton_0[0-3].pdf ../../intro

# rk4_plot.py creates rk4_k[1-4].pdf, rk4_final.pdf
python3 rk4_plot.py
cp -f rk4_k[1-4].pdf rk4_final.pdf ../../intro

# fft_simple_examples.py creates fft-sine-phase.pdf
python3 fft_simple_examples.py
cp -f fft-sine-phase.pdf ../../intro

