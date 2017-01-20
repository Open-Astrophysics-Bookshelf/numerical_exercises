#!/bin/sh

DEST=../../burgers

# characteristics.py makes burgers-characteristics.pdf
python3 characteristics.py
cp -f burgers-characteristics.pdf ${DEST}

# rh.py makes rh.pdf
python3 rh.py
cp -f rh.pdf ${DEST}

# hydro_examples: burgers.py makes fv-burger-rarefaction.pdf fv-burger-sine.pdf 

