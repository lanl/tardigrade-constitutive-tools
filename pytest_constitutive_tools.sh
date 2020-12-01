set -Eeuxo pipefail
cd src/python
# Build the cython tests
python setup.py build_ext --inplace
MPLBACKEND=svg  pytest --cov-report=xml:$(pwd)/cobertura.xml --cov=$(pwd)/. --cov-report=html:$(pwd)/htmlcov
