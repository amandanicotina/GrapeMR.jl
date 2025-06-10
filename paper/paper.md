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
    # TODO(anicotina): Register for an ORCID
    orcid: 0000-0000-0000-0000
    equal-contrib: true
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
date: 9 June 2025
bibliography: paper.bib

---

# Summary

TODO

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
