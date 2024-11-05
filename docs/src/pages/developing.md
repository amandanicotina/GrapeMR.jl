# Developing

In case you are interested in locally developing on GrapeMR, please follow the steps below.
Consider opening an issue if you run into any problems or if you have any suggestions for improvements.
In any case, please open an issue to discuss implementation ideas before starting to code or opening a pull request (PR) directly.

## Getting Started

```bash
julia
julia>] activate .
(GrapeMR) pkg> instantiate
julia>using GrapeMR
```

## Documentation

```bash
julia -e 'using Pkg; Pkg.add("LiveServer")'
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs -e 'using GrapeMR, LiveServer; ENV["DEV"]=true; servedocs(skip_dirs=["docs/src/generated"])'
```

## Testing

To run the tests, run the command below in your shell:

```bash
julia --project=. -e 'using Pkg; Pkg.test("GrapeMR"; coverage=true)'
```

### Test Coverage

In your global julia installation, ensure you have `TestTools` installed.

```julia
using Pkg
Pkg.add("TestTools")
using TestTools; TestTools.install()
```

then add `export PATH="$PATH:/Users/daviddodel/.julia/bin"` to your `~/.zshrc` or equivalent shell configuration file.

Note: You might need to restart the julia language server and VSCode to pick up the changes in the shell environment.

#### Local Coverage Report

```julia
using Pkg; Pkg.test("GrapeMR"; coverage=true)
# ';' is not a typo and will get you a shell within the current julia session
;jlcoverage
# Optional: Clean up all .cov files
using Coverage; Coverage.clean_folder(".");
```

#### VSCode

You might want to consider installing the [Coverage Gutters](https://marketplace.visualstudio.com/items?itemName=ryanluker.vscode-coverage-gutters) extension to enable visual cues inside the code editor.

To this end, you'll need to generate an `lcov.info` in your root directory:

```julia
using Pkg; Pkg.test("GrapeMR"; coverage=true)
using Coverage; LCOV.writefile("lcov.info", process_folder())
```

Go to the Command Palette: `Coverage Gutters: Watch`
