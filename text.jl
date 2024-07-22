t = Template(;
           user="amandanicotina",
           authors=["Amanda Nicotina", "David Dodel"],
           plugins=[
               License(name="MIT"),
               Git(),
               GitHubActions(),
           ],
       )




pTrans = plot()
isos = grape_output3.isochromats
for iso ∈ isos
    m = iso.magnetization.dynamics
    Mt =  m[2,:] + im*m[3,:]
    Mx = real.(Mt)
    My = imag.(Mt)
    plot!(pTrans, Mx, My, color = "black", label = false, 
    title = "Transverse Magnetization",
    titlefontsize = 12,
        xlabel="Mxy", ylabel="Mz")
    scatter!(pTrans, [Mx[end]], [My[end]], color = "black", label = false, framestyle=:box)
end
pMag = plot_magnetization_2D(grape_output3.isochromats)
cf = grape_output3.control_field
time = range(0.0, cf.t_control, length = length(cf.B1x))*1e3
Bx = cf.B1x
By = cf.B1y
B1 = vec(Bx + im * By)

p_Bx = plot(time, abs.(B1), linewidth=2, grid = false, label=false, color=:viridis, 
        framestyle=:box, ylabel="|B1| [Hz]", title="Control Fields", titlefontsize=12)
p_By = plot(time, angle.(B1), linewidth=2, grid = false, label=false, color=:viridis, 
        framestyle=:box, ylabel="ϕ [rad]", xlabel="t [ms]", 
        yticks=([-π, -π/2, 0, π/2, π], ["-π", "-π/2", "0", "π/2", "π"]))

l = @layout [[a{0.7h}; b{0.3h}] c]
pMagCf = plot(p_Bx, p_By, pMag, layout = l)
display(pMagCf)

