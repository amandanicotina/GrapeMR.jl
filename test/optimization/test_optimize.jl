using GrapeMR

# Spin System
M0 = [0.0, 0.0, 1.0]
ΔB1 = [1.0]
B0 = 0.0
offsets = collect(-B0:1:B0)

# Water
T1_water = 0.5
T2_water = 0.1
label_water = "S1"
target_water = "[0.0, 1.0, 0.0]"

spins = GrapeMR.Spin(M0, [T1_water], [T2_water], offsets, ΔB1, [target_water], [label_water])

# Grape Parameters 
grape_params = GrapeParams(1000, GrapeMR.spin_target, Dict("B1x" => true, "B1y" => true, "Bz" => false))

# Optimization Parameters
random_opt = random_hyperopt(spins, grape_params, LinRange(0.01, 1.0, 15), range(10, 100, step = 10))
bohb_opt = bohb_hyperopt(spins, grape_params, LinRange(0.01, 1.0, 15), 100)
hband_opt = hband_hyperopt(spins, grape_params, LinRange(0.01, 1.0, 15), 100)

_ = run_grape_optimization(joinpath(@__DIR__, "../../src/default_config.toml"))