import getpass
import logging

from seiir_model_pipeline.core.task import RegressionTask, ForecastTask, DiagnosticTask
from jobmon.client.swarm.workflow.workflow import Workflow

log = logging.getLogger(__name__)

PROJECT = 'proj_covid'


class SEIIRWorkFlow(Workflow):
    def __init__(self, directories):
        """
        Create a workflow for the SEIIR pipeline.

        :param directories: (Directories)
        :return (Workflow)
        """
        self.directories = directories

        user = getpass.getuser()
        working_dir = f'/ihme/homes/{user}'

        workflow_args = f'seiir-model-'
        if directories.regression_version is not None:
            workflow_args += f'{directories.regression_version} '
        if directories.forecast_version is not None:
            workflow_args += f'{directories.forecast_version} '

        super().__init__(
            workflow_args=workflow_args,
            project=PROJECT,
            stderr=f'/ihme/scratch/users/mnorwood/covid/logs',
            stdout=f'/ihme/scratch/users/mnorwood/covid/logs',
            working_dir=working_dir,
            seconds_until_timeout=60*60*24,
            resume=True
        )

    def attach_regression_tasks(self, n_draws, **kwargs):
        """
        Attach n_draws DrawTasks.

        :param n_draws: (int)
        **kwargs: keyword arguments to DrawTask
        :return: self
        """
        tasks = [RegressionTask(draw_id=i, **kwargs) for i in range(n_draws)]
        self.add_tasks(tasks)
        return tasks

    def attach_forecast_tasks(self, location_ids, add_splicer, **kwargs):
        splicer_tasks = []
        tasks = [
            ForecastTask(location_id=loc, **kwargs)
            for loc in location_ids
        ]
        if add_splicer:
            splicer_tasks = [
                task.add_splicer_task() for task in tasks
            ]
        self.add_tasks(tasks)
        self.add_tasks(splicer_tasks)
        return tasks + splicer_tasks

    def attach_diagnostic_tasks(self, regression_version, forecast_version, **kwargs):
        task = DiagnosticTask(
            regression_version=regression_version,
            forecast_version=forecast_version, **kwargs
        )
        self.add_task(task)