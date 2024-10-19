using PythonCall
using GrapeMR

# Spin System
## Glycoproteins
# T1_gycons = collect(range(650, stop = 1100, length=2).*1e-3)
# T2_gycons = collect(range(10, stop = 100, length=2).*1e-3)
# label_gycons  = fill("Gycloproteins", length(T1_gycons)) 
# target_gycons = fill("max", length(T1_gycons))
## Water
# T1_water = collect(range(1150, stop = 3000, length=2).*1e-3)
# T2_water = collect(range(1500, stop = 3000, length=2).*1e-3)
# label_water  = fill("Water", length(T1 _water)) 
# target_water = fill("min", length(T1_water))

M0 = [0.0, 0.0, 1.0]
ΔB1 = [1.0]
B0 = 5.0
offsets = collect(-B0:1:B0)

# Water
T1_water = 0.5
T2_water = 0.1
label_water = "S1"
target_water = "[0.0, 1.0, 0.0]"

# Glycerol
# T1_glycerol = collect(range(0.3, 0.4, 3))
# T2_glycerol = collect(range(0.07, 0.11, 3))
# label_glycerol= fill("Glycerol", length(T1_glycerol)) 
# target_glycerol = fill("max", length(T1_glycerol))

# T1_glycerol = 0.35
# T2_glycerol = 0.1
# label_glycerol= "Water"
# target_glycerol = "min"

# spins = GrapeMR.Spin(M0, [T1_water; T1_glycerol], [T2_water; T2_glycerol],
#                     offset, ΔB1, [target_water; target_glycerol], [label_water; label_glycerol])

spins = GrapeMR.Spin(M0, [T1_water], [T2_water], offsets, ΔB1, [target_water], [label_water])

# Grape Parameters 
grape_params = GrapeParams(1500, :spin_target, [true true false])

# Optimization Parameters
# bohb = @time hyperoptimization(spins, grape_params, LinRange(0.01, 1.0, 15), 1500)
# plot_bohb(bohb)
# spline_bohb = @time hyperoptimization(spins, grape_params, LinRange(0.01, 1.0, 15), 1500)
# Tc, poly_start, poly_degree, max_iter = spline_bohb.minimizer
Tc = 0.5
poly_start = 0.1
poly_degree = 1.0
max_iter = 10.0
# max_iter = 1000.0
opt_params   = OptimizationParams(poly_start, poly_degree, Int(ceil(max_iter)))

# Parameters 
params = Parameters(grape_params, opt_params)

# Initial RF Pulse
B1ref = 1.0
control_field = spline_RF(grape_params.N, Tc, B1ref) 

# Run Optimization
grape_output = @time grape(params, control_field, spins); 
no_threads_grape_output = @time no_threads_grape(params, control_field, spins);
threads_grape_output = @time threads_grape(params, control_field, spins);

spin = grape_output.isochromats[1].spin
@time run_cost_analysis(grape_output.control_field, spin, 50.0, 50, grape_params.cost_function)

# # # Save data
# folder_path = "/Users/amandanicotina/Documents/PhD/Thesis/SimulationResults/"
# experiment_folder = save_grape_data(grape_output; folder_path)
# go = load_grape_data(experiment_folder)

# # Export data
# export_bruker(grape_output)

# Plots
plot_cost_values(threads_grape_output.cost_values, grape_params)
plot_magnetization_2D(grape_output.isochromats)
plot_magnetization_control_field(threads_grape_output.control_field, grape_output.isochromats)
plot_control_fields(grape_output.control_field; unit="Hz")
plot_magnetization_time(grape_output.isochromats[11], grape_output.control_field.t_control)


function __init__()
    @pyexec """
    import time

    import ConfigSpace as CS
    from hpbandster.core.worker import Worker

    class GrapeWorker(Worker):
        import juliacall
        from juliacall import Pkg as jlPkg
        JL = juliacall.newmodule("PyGrapeMR")

        jlPkg.activate(".")
        JL.seval("using GrapeMR")

        def __init__(self, *args, sleep_interval=0, **kwargs):
            import juliacall
            print("kwargs: ", kwargs)
            super().__init__(*args, **kwargs)
            M0 = [0.0, 0.0, 1.0]
            delta_B1 = [1.0] 
            B0 = 5.0

            # Water
            T1_water = 0.5
            T2_water = 0.1
            label_water = "S1"
            target_water = "[0.0, 1.0, 0.0]"

            # self.spins = self.JL.GrapeMR.Spin(M0, [T1_water], [T2_water], range(-B0, B0+1, 1), delta_B1, [target_water], [label_water])
            self.gp = self.JL.GrapeParams(1500, "spin_target")


        def compute(self, config, budget, **kwargs):
            print(config, budget)
            # B1ref = 1.0
            # control_field = GrapeWorker.JL.spline_RF(self.gp.N, config['Tc'], B1ref) 
            
            # # Optimize
            # opt_params = GrapeWorker.JL.OptimizationParams(config['poly_start'], config['poly_degree'], int(budget))
            # params     = GrapeWorker.JL.Parameters(self.gp, opt_params)
            # res = GrapeWorker.JL.grape(params, control_field, self.spins)

            # loss = res.cost_values[-1]

            return({
                        'loss': float(10),  # this is the a mandatory field to run hyperband
                        'info': {}  # can be used for any user-defined information - also mandatory
                    })
        
        @staticmethod
        def get_configspace():
            import ConfigSpace as CS
            config_space = CS.ConfigurationSpace()
            config_space.add_hyperparameter(CS.UniformFloatHyperparameter('Tc', lower=0.05, upper=1))
            config_space.add_hyperparameter(CS.UniformFloatHyperparameter('poly_start', lower=1e-2, upper=1e-1))
            config_space.add_hyperparameter(CS.UniformIntegerHyperparameter('poly_degree', lower=1, upper=3))
            return(config_space)

    """ => GrapeWorker
end

# python imports
hpns = pyimport("hpbandster.core.nameserver")
optimizers = pyimport("hpbandster.optimizers")

NS = hpns.NameServer(run_id="example1", host="127.0.0.1", port=nothing)
NS.start()

GrapeWorker = __init__()
w = GrapeWorker(nameserver="127.0.0.1", run_id="example1")
w.run(background=true)
bohb = optimizers.BOHB(  configspace = w.get_configspace(),
                         run_id = "example1", nameserver="127.0.0.1",
                         min_budget=0.01, max_budget=1.0
                      )
res = bohb.run(n_iterations=10)

bohb.shutdown(shutdown_workers=true)
NS.shutdown()
# Step 5: Analysis
# Each optimizer returns a hpbandster.core.result.Result object.
# It holds informations about the optimization run like the incumbent (=best) configuration.
# For further details about the Result object, see its documentation.
# Here we simply print out the best config and some statistics about the performed runs.
id2config = res.get_id2config_mapping()
incumbent = res.get_incumbent_id()

println("Best found configuration: ", id2config[incumbent]["config"])
println("A total of $(length(keys(id2config))) unique configurations were sampled.")
println("A total of $(length(res.get_all_runs())) runs were executed.")
println("Total budget corresponds to $(sum([r.budget for r in res.get_all_runs()])/1.0) full function evaluations.")