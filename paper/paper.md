---
title: 'GrapeMR: A Julia package for Gradient Ascent Pulse Engineering for Magnetic Ressonance Applications.'
tags:
  - Julia
  - magnetic ressonance
  - nmr
  - mri
  - optimal control
authors:
  - name: Amanda Nicotina
    orcid: 0009-0004-8524-6436
    equal-contrib: false
    affiliation: 1
  - name: David Dodel
    affiliation: 2
affiliations:
 - name: Steffen J. Glaser, Professor, TUM, Germany
   index: 1
   orcid: 0000-0003-4099-3177
 # TODO(ddodel): check if this is fine
 - name: Independent Researcher, Germany
   index: 2
date: 11 June 2025
bibliography: paper.bib

---

# Summary
In NMR and MRI, radiofrequency (RF) pulses are the means to manipulate spin systems. Designing optimal RF pulse sequences is critical for applications ranging from robust excitation and inversion to saturation contrast and spatial encoding. However, these pulses must account for hardware constraints, field inhomogeneities, relaxation effects, and multi-spin interactions. This makes optimal control a natural framework for RF pulse design.

GrapeMR.jl is an open-source Julia package for simulating and optimizing radiofrequency (RF) pulses in Nuclear Magnetic Resonance (NMR) and Magnetic Resonance Imaging (MRI) using Gradient Ascent Pulse Engineering (GRAPE). The package provides a high-performance, modular, and extensible framework for applying optimal control theory to spin dynamics governed by the Bloch equations.

# Statement of need

GrapeMR.jl is a Julia-based package designed for optimizing RF pulses using the GRAPE algorithm [@khaneja2005optimal], with applications in NMR/MRI. 
It uses Julia’s high-performance computing capabilities to provide an open-source and eﬀicient alternative to previous GRAPE implementations [4, 28], developed in MATLAB.Moreover, the package is designed to be modular and extensible for future updates and new functionalities.
While building on previous algorithmic foundations, GrapeMR.jl focuses on significant performance gains, as demonstrated in the tests presented in Appendix ??.
Although there is room for further optimization, this implementation has outperformed its MATLAB counterparts in many tests.

The performance gains introduced by GrapeMR.jl allow for usage of hyperparameter tuning commonly used in Machine Learning.
Additionally, GrapeMR.jl offers a standalone application via `PackageCompiler.jl` for non-developer users. 
The standalone application can be configured via a simple human-readable `TOML` file and its results exported to scanners commonly used in the area, e.g. Bruker's TopSpin format.

# Acknowledgements

We acknowledge contributions from David Dodel, Guillaume Dalle, and Leo van Damme.
Additionally, we are grateful for the support from the Julia and specifically the JuliaHealth community.

# References
