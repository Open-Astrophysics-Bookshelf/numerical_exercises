#!/bin/sh

# fd.py creates fd_grid.pdf
python3 fd.py
cp -f fd_grid.pdf ../../finite-volume/

# ccfd.py creates ccfd_grid.pdf
python3 ccfd.py
cp -f ccfd_grid.pdf ../../finite-volume/

# fv.py creates fv_grid.pdf
python3 fv.py
cp -f fv_grid.pdf ../../finite-volume/

# simplegrid_gc.py creates simplegrid_gc.pdf
python3 simplegrid_gc.py
cp -f simplegrid_gc.pdf ../../finite-volume/

# domain.py creates domain.pdf
python3 domain.py
cp -r domain.pdf ../../finite-volume/
