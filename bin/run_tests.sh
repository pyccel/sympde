#!/bin/bash

python3 -m pytest sympde/core/tests/test_geometry.py
python3 -m pytest sympde/core/tests/test_derivatives.py
python3 -m pytest sympde/core/tests/test_space.py
python3 -m pytest sympde/core/tests/test_mapping.py
python3 -m pytest sympde/core/tests/test_model_1d.py

#python3 -m pytest sympde/core/tests/test_expr_1d.py
python3 -m pytest sympde/core/tests/test_expr_2d.py
#python3 -m pytest sympde/core/tests/test_expr_3d.py

python3 -m pytest sympde/printing/tests/test_latex.py

python3 -m pytest sympde/gallery/tests/test_poisson.py
python3 -m pytest sympde/gallery/tests/test_stokes.py
python3 -m pytest sympde/gallery/tests/test_wave.py

