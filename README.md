# PDtools
PDtools is a pymatgen-based Python utility to visualize phase diagrams with Igor Pro from ab-initio structural optimization calculations under pressure performed by Quantum Espresso.

This utility contains the following functions:
- Build a phase diagram from the output files of vc-relax calculation performed by Quantum Espresso. (The [pymatgen module](https://github.com/materialsproject/pymatgen) is used internally.)
- Convert PhaseDiagram object (defined in pymatgen module) to the Igor Text format, which can be opened with Igor Pro software.
- Export the formation enthalpy and enthalpy above convex hull for each composition in the phase diagram to the Igor Text format, which can be opened with Igor Pro software.

## Usage
You can try a tutorial from ['tutorial/tutorial.ipynb'](tutorial/tutorial.ipynb).

## Required Packages
The author has tested the operation in the following environment.
- Python 3.8.3
- pymatgen 2020.12.31
- Igor Pro 6.37

## License
PDtools is distributed under the MIT License.
