![Grape Logo](./docs/images/grape_logo.png)
# GrapeMR
 GrapeMR.jl is a Julia package that numerically designes pulse sequences for NMR/MRI applications. By defining a cost funciton, the algorithm calculates step-by-step what is the best pulse sequence for a specific problem.


## Getting Started

```bash
$ julia
julia>] activate .
pkg>
julia>exit()
$ cd docs && julia --project make.jl
$ open build/index.html
```

## Test Coverage

In your global julia installation, ensure you have `TestTools` installed.

```julia
using Pkg
Pkg.add("TestTools")
using TestTools; TestTools.install()
```

then add `export PATH="$PATH:/Users/daviddodel/.julia/bin"` to your `~/.zshrc` or equivalent shell configuration file.

Note: You might need to restart the julia language server and VSCode to pick up the changes in the shell environment.

### Local Coverage Report

```julia
using Pkg; Pkg.test("GrapeMR"; coverage=true)
;jlcoverage
# Optional: Clean up all .cov files
using Coverage; Coverage.clean_folder(".");
```

### VSCode

You might want to consider installing the [Coverage Gutters](https://marketplace.visualstudio.com/items?itemName=ryanluker.vscode-coverage-gutters) extension to enable visual cues inside the code editor.

To this end, you'll need to generate an `lcov.info` in your root directory:

```julia
using Pkg; Pkg.test("GrapeMR"; coverage=true)
using Coverage; LCOV.writefile("lcov.info", process_folder())
```

Go to the Command Palette: `Coverage Gutters: Watch`
