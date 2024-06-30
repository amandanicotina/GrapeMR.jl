using Distributed
const num_processes::Int = 6;
addprocs(num_processes)
println(nprocs())
println(workers())

@everywhere begin
    using GrapeMR
    using Wandb
    using Plots
    using JLD
    using Hyperopt

    M0 = [0.0, 0.0, 1.0]
    T1 = [0.3, 1.0]
    T2 = [0.08, 0.6]
    label  = ["T1=300ms", "T1=1000ms"]  
    target = ["max", "min"]
    B0 = 30.0
    offset = collect(-B0/2:1:B0/2) 
    ΔB1 = [1.0]
    spins = GrapeMR.Spin(M0, T1, T2, offset, ΔB1, target, label)

    # Grape Parameters 
    grape_params = GrapeParams(1000, :saturation_contrast, [true true false])
end

const wandb_project::String = "GrapeMR"
logger::WandbLogger = Wandb.WandbLogger(; project = wandb_project, name = nothing)
# bohb = @time hyperoptimization(spins, grape_params, LinRange(0.05, 1.0, 100), 5000, i=100)
# rand_hopt = random_sample(spins, grape_params, LinRange(0.05, 1.0, 20), range(1000, stop=5000, step=100), i=20)
rand_hopt = random_sample(spins, grape_params, LinRange(0.05, 1.0, 20), range(1000, stop=5000, step=100), i=20)

# close(logger)
# @everywhere
function bohb_params(bohb)
    cost  = []
    t_c   = []
    start = []
    deg   = []
    max   = []
    for i ∈ eachindex(bohb.results)
        push!(cost, bohb.results[i][1])
        push!(t_c, bohb.history[i][1])
        push!(start, bohb.history[i][2])
        push!(deg, bohb.history[i][3])
        push!(max, bohb.history[i][4])
    end
    function get_eps(start, deg, max_iter)
        return start / (1 - (max_iter/2 - 1) / max_iter)^deg
    end
    eps = [get_eps(start, deg, it) for (start, deg, it) in zip(start, deg, max)]
    t_min = bohb.minimizer[1]
    st_min = bohb.minimizer[2]
    d_min = bohb.minimizer[3]
    m_min = bohb.minimizer[4]
    c_min = round(bohb.minimum[1], digits=4)
    order = collect(1:length(t_c))
    m_iter = 1000:100:5000
    p_cost =  scatter(t_c, max, zcolor=cost, 
        markerstrokecolor = :auto, label = false,
        xlabel = "Control time [s]", ylabel = "Iter", colorbar_title="Cost Value",
        title = "Optimization",
        color = :viridis)
        scatter!([t_min], [m_min], label = "Minimum = $c_min", 
        marker = :star5, markersize = 8, color = :red)
    p_order =  scatter(cost, m_iter, zcolor=order, 
        markerstrokecolor = :auto, label = false, 
        xlabel = "Cost Function", ylabel = "Resources - log scale", colorbar_title="Iterations",
        title = "Hyperparameter Tuning - Random Sampler", 
        color = :viridis)
        scatter!([c_min], [m_min], label = "Minimum = $c_min", 
        marker = :star5, markersize = 8, color = :red)

    return p_cost, p_order
end
# end
# p1, p2 = bohb_params(rand_hopt)
# p1
# p2

# rand_hopt

# using JLD2 
# JLD2.@save "/Users/daviddodel/git-repos/GrapeMR.jl/rand_hopt.jld2" rand_hopt
# JLD2.@save "/Users/daviddodel/Documents/amanda/rand_hopt.jl" rand_hopt
# JLD.@save "/Users/daviddodel/Documents/amanda/rand_hopt.jld" rand_hopt
# JLD.@save "/Users/daviddodel/Documents/amanda/rand_hopt.jl" rand_hopt



