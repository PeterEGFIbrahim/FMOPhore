from ..modules.FMOPhore_utility import EnvironmentGuard
# EnvironmentGuard().enforce()
def execute_MP2_personalized(input_file, MP2_file="MP2_personalized"):
    import os
    import shutil
    import re
    import subprocess
    try:
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        os.makedirs(base_name, exist_ok=True)
        shutil.move(input_file, os.path.join(base_name, os.path.basename(input_file)))
        shutil.copy(MP2_file, base_name)
        os.chdir(base_name)
        with open(MP2_file, "r") as file:
            MP2_content = file.read()
        MP2_content = re.sub(r"-N ", f"-N x{re.escape(base_name)}", MP2_content)
        MP2_content = MP2_content.replace("input.inp", os.path.basename(input_file))
        MP2_content = MP2_content.replace("output.log", f"{base_name}.log")
        with open(MP2_file, "w") as file:
            file.write(MP2_content)
        os.chmod(MP2_file, 0o755)
        subprocess.Popen(f"nohup ./{MP2_file} &", shell=True)
    except Exception as e:
        print(f"An error occurred: {e}")
    finally:
        os.chdir("..")
