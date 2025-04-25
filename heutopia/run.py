import datetime
import multiprocessing
import os
import pathlib
import random
import shutil
import string
import subprocess
import tempfile
import time
from multiprocessing.pool import AsyncResult
from typing import Any

import pandas as pd  # type: ignore
import yaml
from pydantic import RootModel

from heutopia.config import HeutopiaConfig
from heutopia.cpphelper import expandcpp
from heutopia.exporter import (
    AbstractExporter,
    BatchResult,
    GSheetsExporter,
    LocalFileExporter,
    StdoutExporter,
)
from heutopia.runner import standalone_run

if __name__ == "__main__":
    with open("config.yaml", "r") as f:
        config_dict = yaml.safe_load(f)

    config = HeutopiaConfig(**config_dict)

    print("Config:")
    print(RootModel(config).model_dump_json(indent=4, exclude_none=True))
    print("")

    start_time = time.perf_counter()
    pool = multiprocessing.get_context("fork").Pool(processes=config.num_process)
    num_case_limit = config.num_case_limit

    run_id = "{}{}".format(
        datetime.datetime.now().strftime("%y%m%d%H%M%S"),
        "".join(random.choice(string.ascii_lowercase) for _ in range(2)),
    )

    exporters: list[AbstractExporter] = [
        StdoutExporter(config.print_columns),
        LocalFileExporter(output_dirpath=config.result_dir.format(RUN_ID=run_id)),
    ]

    if config.google_sheets:
        exporters.append(
            GSheetsExporter(
                run_id=run_id,
                gcp_sa_path=config.google_sheets.gcp_sa_path,
                sheet_key=config.google_sheets.sheet_key,
            )
        )

    print("\nEnabled exporters:")
    for exporter in exporters:
        print("  - {}".format(exporter.summary()))
    print("")

    cpp_src = expandcpp(os.path.abspath("./main.cpp"))
    for exporter in exporters:
        exporter.dump_src(cpp_src)

    input_files: list[str] = list()
    for root, _, files in os.walk(config.dataset_dir):
        input_files.extend([str((pathlib.Path(root) / f).resolve()) for f in files])
    input_files = sorted(input_files)[: config.num_case_limit]
    print("Input: {} ({} cases)".format(input_files, len(input_files)))

    result_list: list[Any] = list()

    with tempfile.TemporaryDirectory() as tmpdir:
        subprocess.check_output("make", shell=True)
        solver_path = str((pathlib.Path(tmpdir) / "solver.out").resolve())
        shutil.copy2("./solver.out", solver_path)

        if config.parallel:
            process_list: list[AsyncResult] = [
                pool.apply_async(
                    standalone_run,
                    args=(solver_path, input_filepath, config.runner, True),
                )
                for input_filepath in input_files
            ]
            pool.close()
            pool.join()

            result_list = [x.get() for x in process_list if x.successful()]
        else:
            for input_filepath in input_files:
                result_list.append(
                    standalone_run(solver_path, input_filepath, config.runner, True)
                )

    end_time = time.perf_counter()

    result = BatchResult(
        df=pd.DataFrame(result_list)
        .sort_values("input_filename", ascending=True)
        .reset_index(drop=True)
    )

    for exporter in exporters:
        exporter.dump(result)

    print("Elapsed time: {} sec.".format(end_time - start_time))
