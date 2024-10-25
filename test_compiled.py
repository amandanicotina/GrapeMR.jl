from datetime import datetime
import ConfigSpace as CS
from hpbandster.core.worker import Worker
import hpbandster.core.nameserver as hpns
import hpbandster.optimizers as optimizers
import logging
from subprocess import Popen, PIPE
import re
from pathlib import Path
from argparse import ArgumentParser

LOSS_REGEX = re.compile(r".*Final Cost Function Value = (\d+\.\d+)")
logging.basicConfig(level=logging.DEBUG)


def run_grape(args: list[str]) -> float:
    p = Popen(["./GrapeMRCompiled/bin/GrapeMR"] + args, stdout=PIPE, stderr=PIPE, text=True)
    stdout, _ = p.communicate()
    return float(LOSS_REGEX.search(stdout).group(1))


class GrapeWorker(Worker):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def compute(self, config_id, config, budget, working_directory):
        loss = run_grape(list(map(str, [config["Tc"], config["poly_start"], config["poly_degree"], budget])))

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
    parser = ArgumentParser()
    parser.add_argument("--worker", action="store_true", default=False)
    args = parser.parse_args()
    
    print("Starting")
    
    if args.worker:
        w = GrapeWorker(nameserver="127.0.0.1", run_id="example1")
        w.run(background=False)
        exit()
    
    else:
        NS = hpns.NameServer(run_id="example1", host="127.0.0.1", port=None)
        NS.start()

        w = GrapeWorker(nameserver="127.0.0.1", run_id="example1")
        w.run(background=True)
    
        import hpbandster.core.result as hpres
        results_path = Path("results") / datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        results_path.mkdir(parents=True, exist_ok=True)
        result_logger = hpres.json_result_logger(directory=results_path, overwrite=True)
        
        bohb = optimizers.BOHB(  configspace = GrapeWorker.get_configspace(),
                                run_id = "example1", nameserver="127.0.0.1",
                                result_logger=result_logger,
                                min_budget=500, max_budget=2000
                            )
        res = bohb.run(n_iterations=500)

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
        
        import matplotlib.pyplot as plt
        import hpbandster.core.result as hpres
        import hpbandster.visualization as hpvis



        # load the example run from the log files
        result = hpres.logged_results_to_HBS_result(str(results_path))

        # get all executed runs
        all_runs = result.get_all_runs()

        # get the 'dict' that translates config ids to the actual configurations
        id2conf = result.get_id2config_mapping()


        # Here is how you get he incumbent (best configuration)
        inc_id = result.get_incumbent_id()

        # let's grab the run on the highest budget
        inc_runs = result.get_runs_by_id(inc_id)
        inc_run = inc_runs[-1]


        # We have access to all information: the config, the loss observed during
        #optimization, and all the additional information
        inc_loss = inc_run.loss
        inc_config = id2conf[inc_id]['config']
        # inc_test_loss = inc_run.info['test accuracy']

        print('Best found configuration:')
        print(inc_config)
        print(inc_loss)


        # Let's plot the observed losses grouped by budget,
        hpvis.losses_over_time(all_runs)

        # the number of concurent runs,
        hpvis.concurrent_runs_over_time(all_runs)

        # and the number of finished runs.
        hpvis.finished_runs_over_time(all_runs)

        # This one visualizes the spearman rank correlation coefficients of the losses
        # between different budgets.
        hpvis.correlation_across_budgets(result)

        # For model based optimizers, one might wonder how much the model actually helped.
        # The next plot compares the performance of configs picked by the model vs. random ones
        # hpvis.performance_histogram_model_vs_random(all_runs, id2conf)

        plt.show()
