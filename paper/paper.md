---
title: 'GrapeMR: A Julia package for Gradient Ascent Pulse Engineering for Magnetic Ressonance Applications.'
tags:
  - Julia
  - magnetic ressonance
  - nmr
  - mri
  - optimal control
authors:
  - name: Amanda Nicotina Pereira
    orcid: 0000-0000-0000-0000
    equal-contrib: true
    affiliation: 1
  - name: David Dodel
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 2
affiliations:
 - name: Steffen J. Glaser, Professor, TUM, Germany
   index: 1
   ror: 00hx57361
 - name: Independent Researcher, Germany
   index: 2
date: 9 June 2025
bibliography: paper.bib

---

# Summary

The forces on stars, galaxies, and dark matter under external gravitational
fields lead to the dynamical evolution of structures in the universe. The orbits
of these bodies are therefore key to understanding the formation, history, and
future state of galaxies. The field of "galactic dynamics," which aims to model
the gravitating components of galaxies to study their structure and evolution,
is now well-established, commonly taught, and frequently used in astronomy.
Aside from toy problems and demonstrations, the majority of problems require
efficient numerical tools, many of which require the same base code (e.g., for
performing numerical orbit integration).

# Statement of need

GrapeMR.jl is a Julia-based package designed for optimizing RF pulses using the GRAPE algorithm [@khaneja:2005], with applications in NMR/MRI. 
It uses Julia’s high-performance computing capabilities to provide an open-source and eﬀicient alternative to previous GRAPE implementations [4, 28], developed in MATLAB.Moreover, the package is designed to be modular and extensible for future updates and new functionalities.
While building on previous algorithmic foundations, GrapeMR.jl focuses on significant performance gains, as demonstrated in the tests presented in Appendix ??.
Although there is room for further optimization, this implementation has outperformed its MATLAB counterparts in many tests.

The performance gains introduced by GrapeMR.jl allow for usage of hyperparameter tuning commonly used in Machine Learning.
Additionally, GrapeMR.jl offers a standalone application via `PackageCompiler.jl` for non-developer users. 
The standalone application can be configured via a simple human-readable `TOML` file and its results exported to scanners commonly used in the area, e.g. Bruker's TopSpin format.

<!-- # Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% } -->

# Acknowledgements

We acknowledge contributions from David Dodel, Guillaume Dalle, and Leo van Damme.
Additionally, we are grateful for the support from the Julia and specifically the JuliaHealth community.

# References