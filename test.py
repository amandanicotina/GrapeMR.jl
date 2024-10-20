import ConfigSpace as CS
from hpbandster.core.worker import Worker
import hpbandster.core.nameserver as hpns
import hpbandster.optimizers as optimizers
import threading
import logging

from juliacall import Main as jl
from juliacall import Pkg as jlPkg

jlPkg.activate(".")
jlPkg.instantiate()
jl.seval("using GrapeMR")

logging.basicConfig(level=logging.DEBUG)
jl_yield = getattr(jl, "yield")


class GrapeWorker(Worker):

    def __init__(self, *args, **kwargs):
        # jlPkg.activate(".")
        # jlPkg.instantiate()
        # jl.seval("using GrapeMR")
        print("kwargs: ", kwargs)
        super().__init__(*args, **kwargs)
        M0 = [0.0, 0.0, 1.0]
        delta_B1 = [1.0] 
        B0 = 5

        # Water
        T1_water = 0.5
        T2_water = 0.1
        label_water = "S1"
        target_water = "[0.0, 1.0, 0.0]"

        self.spins = jl.GrapeMR.Spin(M0, [T1_water], [T2_water], range(-B0, B0+1, 1), delta_B1, [target_water], [label_water])
        self.gp = jl.GrapeParams(1500, "spin_target")

        self._timer = threading.Timer(1e-2, jl_yield)
        self._timer.start()

    def compute(self, config_id, config, budget, working_directory):
        B1ref = 1.0
        test = jl.Array([0, 1, 2, 3, 10])
        test._jl_display()
        control_field = jl.spline_RF(self.gp.N, config['Tc'], B1ref)
        control_field._jl_display()
        # # Optimize
        # opt_params = jl.OptimizationParams(config['poly_start'], config['poly_degree'], int(budget))
        # # opt_params._jl_display()
        # params     = jl.Parameters(self.gp, opt_params)
        # # params._jl_display()
        
        # res = jl.grape(params, control_field, self.spins)
        # res._jl_display()

        # loss = res.cost_values[-1]
        loss = test[0]

        return ({
            'loss': float(loss),  # this is the a mandatory field to run hyperband
            'info': {}  # can be used for any user-defined information - also mandatory
        })
    
    @staticmethod
    def get_configspace():
        config_space = CS.ConfigurationSpace()
        config_space.add_hyperparameter(CS.UniformFloatHyperparameter('Tc', lower=0.05, upper=1))
        config_space.add_hyperparameter(CS.UniformFloatHyperparameter('poly_start', lower=1e-2, upper=1e-1))
        config_space.add_hyperparameter(CS.UniformIntegerHyperparameter('poly_degree', lower=1, upper=3))
        return(config_space)


if __name__ == "__main__":
    jl.seval("ccall(:jl_enter_threaded_region, Cvoid, ())")
    
    print("Starting")

    NS = hpns.NameServer(run_id="example1", host="127.0.0.1", port=None)
    NS.start()

    w = GrapeWorker(nameserver="127.0.0.1", run_id="example1")
    w.run(background=True)
    
    bohb = optimizers.BOHB(  configspace = GrapeWorker.get_configspace(),
                            run_id = "example1", nameserver="127.0.0.1",
                            min_budget=1, max_budget=1
                        )
    res = bohb.run(n_iterations=10)

    bohb.shutdown(shutdown_workers=True)
    NS.shutdown()
    # Step 5: Analysis
    # Each optimizer returns a hpbandster.core.result.Result object.
    # It holds informations about the optimization run like the incumbent (=best) configuration.
    # For further details about the Result object, see its documentation.
    # Here we simply print out the best config and some statistics about the performed runs.
    id2config = res.get_id2config_mapping()
    incumbent = res.get_incumbent_id()

    print('Best found configuration:', id2config[incumbent]['config'])
    print('A total of %i unique configurations where sampled.' % len(id2config.keys()))
    print('A total of %i runs where executed.' % len(res.get_all_runs()))
    print('Total budget corresponds to %.1f full function evaluations.'%(sum([r.budget for r in res.get_all_runs()])/1.0))
    jl.seval("ccall(:jl_exit_threaded_region, Cvoid, ())")
