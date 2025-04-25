import json
import os
import subprocess
import time
from typing import Union

from heutopia.config import RunnerConfig


def standalone_run(
    solver_command: str,
    input_fullpath: str,
    runner_config: RunnerConfig,
    show_detail: bool,
) -> dict[str, Union[str, float, int]]:
    comm = (
        runner_config.run_command.replace("{INPUT_FILE}", input_fullpath)
        .replace("{SOLVER_CMD}", solver_command)
        .replace("{SOLVER_OUTPUT}", "/dev/null")
    )

    begin_ns = time.perf_counter_ns()
    process = subprocess.run(
        comm, shell=True, check=False, capture_output=True, text=True
    )
    end_ns = time.perf_counter_ns()
    spent_ms = (end_ns - begin_ns) // 1000000

    input_filename = os.path.basename(input_fullpath)

    if process.returncode:
        print("Error: {}".format(input_filename))

    ret: dict[str, Union[str, float, int]] = {"input_filename": input_filename}
    ret |= json.loads(process.stdout)
    ret["msec"] = spent_ms

    if show_detail:
        print("Run {}: {}".format(input_filename, ret))

    return ret
