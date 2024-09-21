

# function outside_optimization_dynamics(iso::Vector{Isochromat}, cf::ControlField)
#     # Dynamics
#     plot_magnetization_2D(iso)

#     # Cost?
#     c = [cost[:,1][i][1] for i in 1:length(spins_range)]
#     g = [cost[:,1][i][2] for i in 1:length(spins_range)]

#     all_costs = Matrix{Float64}(undef, length(T1), length(T2));
#     all_costs = reshape(c, (length(T1), length(T2)))
#     h1 = contourf(T1, T2, all_costs*1e3, color=:viridis, framestyle=:box,
#             xlabel="T1 [s]", ylabel="T2 [s]", title="Cost Function Map")
# end


function topspin_dynamics_comparison(file_path::String, iso::Vector{GrapeMR.Isochromat}; s::String)
    # Call the function and store the result
    data = read_nmrsim_data(file_path)
    if s == "Transverse"
        p = plot_tmagnetization_simul_vs_topspin(data, iso, rf_time)
    else 
        p = plot_magnetization_simul_vs_topspin(data, iso, rf_time) 
    end

    return p
end

function read_nmrsim_data(file_path::String)
    # Empty arrays to hold the data
    Mx = Float64[]
    My = Float64[]
    Mz = Float64[]

    # Open file
    open(file_path, "r") do file
        for line in eachline(file)
            if startswith(line, ";")
                continue
            end
            # Split the line into components by whitespace and convert to Float64
            values = split(strip(line), r"\s+")
            if length(values) == 3
                push!(Mx, parse(Float64, values[1]))
                push!(My, parse(Float64, values[2]))
                push!(Mz, parse(Float64, values[3]))
            end
        end
    end

    # Return the data as a DataFrame
    return DataFrame(Mx=Mx, My=My, Mz=Mz)
end


function plot_magnetization_simul_vs_topspin(df::DataFrame, iso::GrapeMR.Isochromat, t::Float64)
        # TopSpin nmrsim
        Mx_ts = df.Mx
        My_ts = df.My
        Mz_ts = df.Mz
    
        # Grape
        m  = iso.magnetization.dynamics
        s = iso.spin
        time = range(0.0, t, length = length(m[1,:]))

        pMx = plot(time, m[2,:], label = false, lw = 2, color = cgrad(:viridis)[128],
            ylabel = "Mx",
            title  = "Magnetization Dynamics - Sample = $(s.label), Target = $(s.target)",
            xlims = [0.0, 0.5],
            titlefontsize = 12)
            scatter!(time, Mx_ts, label = false, markersize = 0.8, color = "black")
    
        pMy = plot(time, -m[3,:], label = false, lw = 2, color = cgrad(:viridis)[128],
            xlims = [0.0, 0.5],
            ylabel = "My")
            scatter!(time, My_ts, label = false, markersize = 0.8, color = "black")

    
        pMz =   plot(time, m[4,:], label = "GrapeMR", lw = 2, color = cgrad(:viridis)[128], 
                ylabel = "Mz",
                xlabel = "t [s]",
                xlims = [0.0, 0.5],
                legend = :bottomleft)
                scatter!(time, Mz_ts, label = "TopSpin", markersize = 0.8, color = "black")
        pMag = plot(pMx, pMy, pMz, layout = (3,1))
    
        return pMag
end


function plot_tmagnetization_simul_vs_topspin(df::DataFrame, iso::GrapeMR.Isochromat, t::Float64)
    # TopSpin nmrsim
    Mz_ts = df.Mz
    Mt_ts = sqrt.(df.Mx.^2 + df.My.^2)

    # Grape
    m  = iso.magnetization.dynamics
    s = iso.spin
    time = range(0.0, t, length = length(m[1,:]))
    mtrans = sqrt.(m[2,:].^2 + m[3,:].^2)

    pTrans = plot(time, mtrans, label = false, lw = 2, color = cgrad(:viridis)[128],
                ylabel = "Transverse",                
                xlabel = "t [s]",
                xlims = [0.9, time[end]],
                legend = :bottomleft)
                scatter!(time, Mt_ts, label = false, markersize = 1, color = "black")
    pMz =   plot(time, m[4,:], label = "GrapeMR", lw = 2, color = cgrad(:viridis)[128], 
                xlims = [0.9, time[end]],
                ylabel = "Longitudinal",
                title  = "Magnetization Dynamics - Sample = $(s.label)", #, Target = $(s.target)",
                titlefontsize = 12)

                scatter!(time, Mz_ts, label = "TopSpin", markersize = 1, color = "black")
    pMag = plot(pMz, pTrans, layout = (2,1))

    return pMag
end

 




