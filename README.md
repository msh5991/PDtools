# PDtools
PDtools is a Python utility that provides a workflow for constructing and visualizing phase diagrams from ab-initio structural optimization calculations.

This utility contains the following functions:
- Build a phase diagram from the output files of vc-relax calculation performed by Quantum Espresso. (The pymatgen module are used internally.)
- Convert PhaseDiagram object (defined in pymatgen module) to the Igor Text format, which can be opened with Igor Pro software.
- Export the formation enthalpy and enthalpy above convex hull for each composition in the phase diagram to the Igor Text format, which can be opened with Igor Pro software.

## Usage
You can try a tutorial from ['tutorial/tutorial.ipynb'](tutorial/tutorial.ipynb).

## Required Packages
The author has tested the operation in the following environment.
- Python 3.8.3
- pymatgen 2020.12.31
- Igor Pro 6.37

## Licence
PDtools is distributed under the MIT License.
