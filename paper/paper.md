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
  #- name: David Dodel
   # affiliation: 2
affiliations:
 - name: Technical University of Munich, TUM School of Natural Sciences, Germany
   index: 1
   orcid: 0000-0003-4099-3177
 # TODO(ddodel): check if this is fine
 #- name: Independent Researcher, Germany
  # index: 2
date: 11 June 2025
bibliography: paper.bib

---

# Summary

In Nuclear Magnetic Resonance (NMR) and Magnetic Resonance Imaging (MRI), radiofrequency (RF) pulses are used to manipulate spin systems for tasks such as robust excitation, inversion, and saturation contrast. Designing these RF pulses requires consideration of hardware limitations, field inhomogeneities, and relaxation effects. This makes Optimal Control Theory (OCT) a powerful framework for pulse optimization.

GrapeMR.jl is an open-source Julia package for simulating and optimizing RF pulses using the GRadient Ascent Pulse Engineering (GRAPE) algorithm [@khaneja2005optimal]. GRAPE discretizes the control pulse into time steps and iteratively adjusts the control fields to minimize (or maximize) a cost function that reflects the desired spin dynamics. At each iteration, the gradient of the cost function concerning the control fields is computed, enabling efficient updates toward optimal control. GRAPE is well-suited to high-dimensional and constrained optimization problems in NMR and MRI contexts.

# Statement of need

**GrapeMR.jl** was designed to be extensible, efficient, easy to use, and accessible even for users who are not deeply familiar with programming. Optimal control in magnetic resonance applications can be powerful, but it remains significantly underutilized in practice, primarily due to the complexity and inaccessibility of existing tools and their performance limitations.

Several GRAPE implementations exist, such as the MATLAB-based toolkit by Van Reeth et al. for MRI contrast optimization [@van2017optimal]. However, they often rely on proprietary software, are not extensible for different applications, and lack modern features like automated hyperparameter tuning. GrapeMR.jl addresses these limitations by leveraging Julia's features like type stability, memory allocation, and static arrays. It consistently outperforms its MATLAB counterparts, achieving speedups of up to 20× in representative test cases. In addition, GrapeMR.jl automatically computes the gradient of the cost function, eliminating the need for users to manually derive and implement these expressions when testing new optimization goals.

To further enhance accessibility, a standalone application built via PackageCompiler.jl is provided. This implementation allows users to run optimizations without interacting directly with Julia code, making the package suitable for non-developers. Optimization parameters can be defined through a Julia script or a human-readable `TOML` configuration file. Results are automatically exported in formats compatible with widely used NMR software, such as Bruker TopSpin. The export system is modular and easily extensible, supporting integration with other scanner formats.

The main goal of **GrapeMR.jl** lowers the barrier to entry to promote broader adoption of optimal control-designed RF pulses in research and clinical settings.

# How the Package Works

The **GrapeMR.jl** package allows users to perform GRAPE-based optimization by defining the spin system, inhomogeneity distributions, and optimization parameters. Users specify the number of time steps, cost function, control fields to optimize, and the initial guess for the RF pulse, commonly generated via a cubic spline interpolation using `spline_RF()` (or alternatives like `hard_RF()`, `sinc_RF()`, and `gaussian_RF()`).

Optimizations can be run in two main ways:

- **Script-based workflow:** Users define all relevant variables in a Julia script and call the main functions directly.
- **TOML-based workflow:** Users define a configuration file containing all parameters and execute the optimization from the REPL or a compiled binary without writing code.

The core function `grape()` executes the optimization and returns optimized fields, cost function history, and magnetization trajectories. Visualization functions such as `plot_magnetization_control_field()` and `plot_cost_values()` assist in interpreting the results. If no file path is specified, results are automatically saved in a default folder within the package directory.

Thanks to Julia's multiple dispatch,**GrapeMR.jl** supports seamless integration with optimization algorithms from Optim.jl [@mogensen2018optim]. This allows users to switch between gradient descent, BFGS, L-BFGS, and others, depending on their problem's characteristics and constraints.

Hyperparameter search is handled via a hyper_opt flag: if set to `true`, the package automatically performs a hyperparameter search; otherwise, it uses the provided user values.

# Acknowledgements

We acknowledge contributions from David Dodel, Guillaume Dalle, and Leo van Damme.
Additionally, we are grateful for the support from the Julia and specifically the JuliaHealth community.

# References
