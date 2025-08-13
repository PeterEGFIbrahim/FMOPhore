from ..modules.FMOPhore_utility import EnvironmentGuard
# EnvironmentGuard().enforce()
def execute_DFTB(input_file, DFTB_file="DFTB"):
    import os
    import shutil
    import re
    import subprocess

    base_name = os.path.splitext(input_file)[0]
    os.makedirs(base_name, exist_ok=True)
    shutil.move(input_file, os.path.join(base_name, input_file))
    shutil.copy(DFTB_file, base_name)
    os.chdir(base_name)
    with open(DFTB_file, "r") as file:
        DFTB_content = file.read()
    DFTB_content = re.sub(r"-N ", f"-N x{base_name}", DFTB_content)
    DFTB_content = DFTB_content.replace("input.inp", input_file)
    DFTB_content = DFTB_content.replace("output.log", f"{base_name}.log")
    with open(DFTB_file, "w") as file:
        file.write(DFTB_content)
    subprocess.run(["qsub", "-jc", "short", DFTB_file])
    os.chdir("..")