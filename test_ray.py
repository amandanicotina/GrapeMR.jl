from bohb import BOHB
import bohb.configspace as cs

from juliacall import Main as jl
from juliacall import Pkg as jlPkg

jl_yield = getattr(jl, "yield")


class GrapeWorker:

    def __init__(self, *args, **kwargs):
        jlPkg.activate(".")
        jl.seval("using GrapeMR")
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


    def evaluate(self, config, budget):
        B1ref = 1.0
        # jl_yield()
        control_field = jl.GrapeMR.spline_RF(self.gp.N, config['Tc'], B1ref)
        # jl_yield()
        # Optimize
        opt_params = jl.GrapeMR.OptimizationParams(config['poly_start'], config['poly_degree'], int(budget))
        opt_params._jl_display()
        params     = jl.GrapeMR.Parameters(self.gp, opt_params)
        params._jl_display()
        
        res = jl.GrapeMR.grape(params, control_field, self.spins)

        loss = res.cost_values[-1]
        # loss = 10

        return float(loss)


if __name__ == "__main__":
    # Optional: Pass the parameter space yourself
    tc = cs.UniformHyperparameter('Tc', lower=0.05, upper=1)
    ps = cs.CategoricalHyperparameter('poly_start', choices=[1e-2, 1e-1])
    pd = cs.IntegerUniformHyperparameter('poly_degree', lower=1, upper=3)
    config_space = cs.ConfigurationSpace([tc, ps, pd])
    
    worker = GrapeWorker()
    opt = BOHB(configspace=config_space, evaluate=worker.evaluate, max_budget=10, min_budget=1)
    logs = opt.optimize()
    
    print(logs)
