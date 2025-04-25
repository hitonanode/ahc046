import os
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))


def expandcpp(main_cpp_abs_path: str) -> str:
    return subprocess.check_output(
        'LOG_LEVEL={} python3 scripts/expander.py "{}" -c'.format(
            "WARNING",  # Suppress log info such as `[INFO] include: common.hpp`
            main_cpp_abs_path,
        ),
        shell=True,
        text=True,
        cwd=SCRIPT_DIR + "/..",
    )
