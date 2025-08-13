from ..modules.FMOPhore_utility import EnvironmentGuard
# EnvironmentGuard().enforce()
def execute_MP2(input_file, MP2_file="MP2"):
    import os
    import shutil
    import re
    import subprocess

    base_name = os.path.splitext(input_file)[0]
    os.makedirs(base_name, exist_ok=True)
    shutil.move(input_file, os.path.join(base_name, input_file))
    shutil.copy(MP2_file, base_name)
    os.chdir(base_name)
    with open(MP2_file, "r") as file:
        MP2_content = file.read()
    MP2_content = re.sub(r"-N ", f"-N x{base_name}", MP2_content)
    MP2_content = MP2_content.replace("input.inp", input_file)
    MP2_content = MP2_content.replace("output.log", f"{base_name}.log")
    with open(MP2_file, "w") as file:
        file.write(MP2_content)
    subprocess.run(["qsub", "-jc", "long", MP2_file])
    os.chdir("..")