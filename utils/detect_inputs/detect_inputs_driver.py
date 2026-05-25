from hamilton import driver
from hamilton.execution import executors
from typing import Dict, Any
from utils.connect_s3.download_lookup_tables import fetch_lookup_tables

import utils.detect_inputs.detect_inputs_subdag as detect_inputs_subdag


def instantiate_driver(config: Dict[str, Any]):
    """
    Downloads lookup tables to tmp/
    Builds a Hamilton driver, passes lookup table paths and config to driver.
    """
    remote_executor = executors.MultiThreadingExecutor(max_tasks=config.get("cpus", 1))

    lookup_tables = fetch_lookup_tables()

    dr = (
        driver.Builder()
        .enable_dynamic_execution(allow_experimental_mode=True)
        .with_remote_executor(remote_executor)
        .with_modules(detect_inputs_subdag)
        .build()
    )

    result = dr.execute(
        ['return_collected_results'], 
        inputs={'config': config, 'lookup_tables': lookup_tables}
    )['return_collected_results']

    return result

