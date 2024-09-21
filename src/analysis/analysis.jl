"""
   run_analysis()

Run magnetization analysis 

"""
## Gycloproteins
T1_gycons = collect(range(650, stop = 1100, length=100).*1e-3)
T2_gycons = collect(range(10, stop = 100, length=100).*1e-3)
label_gycons  = fill("Glycoproteins", length(T1_gycons)) 
target_gycons = fill("max", length(T1_gycons))
## Water
T1_water = collect(range(1150, stop = 3000, length=100).*1e-3)
T2_water = collect(range(1500, stop = 3000, length=100).*1e-3)
label_water  = fill("Water", length(T1_water)) 
target_water = fill("min", length(T1_water))
## Outside Optimization range
T1_outside = collect(range(1150, stop = 3000, length=100).*1e-3)
T2_outside = collect(range(1500, stop = 3000, length=100).*1e-3)
label_outside  = fill("Outside", length(T1_water)) 
target_outside = fill("-", length(T1_water))
spins_outside = GrapeMR.Spin(M0, [T1_gycons; T1_water; T1_outside], [T2_gycons; T2_water; T1_outside], [0.0], ΔB1, [target_gycons; target_water; target_outside], [label_gycons; label_water; label_outside])

# No Relaxation
cf = res_grape.control_field
rf_time = res_grape.control_field.t_control
spin_noRelax = GrapeMR.Spin([0.0, 0.0, 1.0], [1e10], [1e10], [0.0], [1.0], ["-"], ["noRelax"])
mag_noRelax = forward_propagation(cf, spin_noRelax[1])
dyn_noRelax = GrapeMR.Magnetization(mag_noRelax)
iso_noRelax = Isochromat(dyn_noRelax, spin_noRelax[1])
plot_magnetization_time(iso_noRelax, rf_time)
file_path = "/Users/amandanicotina/Documents/Documents/ProgressReports/ResultsGrapeMR/Metabolomics/TopSpin/oc_grape_gycloproteins_1.txt"

function run_analysis()
    
end

"""
    run_analysis()

Run RF pulse analysis 

"""
function run_analysis()
   
"""
    run_analysis()

Run cost function analysis 

"""
end
function run_analysis()
    
end