"""Project root for all scripts in this repository.

GMO_ROOT is taken from the environment variable GMO_ROOT; if it is not set, the
root of this repository (the directory containing this file) is used.
Large inputs that are not shipped with the repository (genome assemblies,
proteomes, reads, NCBI nr) are expected under $GMO_ROOT/data/ (see README.md).
"""
import os

GMO_ROOT = os.environ.get("GMO_ROOT") or os.path.dirname(os.path.abspath(__file__))
