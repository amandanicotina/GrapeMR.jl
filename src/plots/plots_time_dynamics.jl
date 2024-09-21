
function plot_magnetization_time(iso::Isochromat, t::Float64)
    mag  = iso.magnetization.dynamics
    spin = iso.spin
    time = range(0.0, t, length = length(mag[1,:]))

    # Create a plot of the magnetization
    p = plot(time, mag[2:end,:]', label = ["Mx" "My" "Mz"], lw = 2,
        xlabel = "t [sec]",
        ylabel = "Magnitude",
        title  = "Magnetization Dynamics - Sample = $(spin.label), Target = $(spin.target)",
        titlefontsize = 12,
        )
    return p
end
 
# function plot_magnetization_simul_vs_topspin(df::DataFrame, iso::Isochromat, t::Float64)
#     # TopSpin nmrsim
#     Mx_ts = df.Mx
#     My_ts = df.My
#     Mz_ts = df.Mz
#     Mt_ts = sqrt.(Mx_ts.^2 + My_ts.^2)

#     # Grape
#     m  = iso.magnetization.dynamics
#     s = iso.spin
#     time = range(0.0, t, length = length(m[1,:]))

#     # pMx = plot(time, m[2,:], label = false, lw = 2, color = cgrad(:viridis)[128],
#     #     ylabel = "Mx",
#     #     title  = "Magnetization Dynamics - Sample = $(s.label), Target = $(s.target)",
#     #     titlefontsize = 12)
#     #     plot!(time, Mx_ts, label = false, lw = 2, color = "black")

#     # pMy = plot(time, m[3,:], label = false, lw = 2, color = cgrad(:viridis)[128],
#     #     ylabel = "My")
#     #     plot!(time, My_ts, label = false, lw = 2, color = "black")
#     mtrans = sqrt.(m[2,:].^2 + m[3,:].^2)

#     pTrans = plot(time, mtrans, label = false, lw = 2, color = cgrad(:viridis)[128],
#             ylabel = "Transverse",
#             title  = "Magnetization Dynamics - Sample = $(s.label)", #, Target = $(s.target)",
#             titlefontsize = 12)
#             plot!(time, Mt_ts, label = false, lw = 2, color = "black")

#     pMz =   plot(time, m[4,:], label = "GrapeMR", lw = 2, color = cgrad(:viridis)[128],
#             ylabel = "Longitudinal",
#             xlabel = "t [s]")
#             plot!(time, Mz_ts, label = "TopSpin", lw = 2, color = "black")
    
#     pMag = plot(pTrans, pMz, layout = (2,1))

#     return pMag
# end

