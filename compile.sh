julia -e 'using PackageCompiler; create_app(\".\", \"GrapeMRCompiled\", force=true; precompile_execution_file=\"src/compilation_main.jl\")'
