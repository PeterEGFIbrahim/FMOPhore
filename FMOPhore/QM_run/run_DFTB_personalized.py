from ..modules.FMOPhore_utility import EnvironmentGuard
# EnvironmentGuard().enforce()

def execute_DFTB_personalized(input_file, DFTB_file="DFTB_personalized"):
    import os
    import shutil
    import re
    import subprocess
    try:
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        os.makedirs(base_name, exist_ok=True)
        shutil.move(input_file, os.path.join(base_name, os.path.basename(input_file)))
        shutil.copy(DFTB_file, base_name)
        os.chdir(base_name)
        with open(DFTB_file, "r") as file:
            DFTB_content = file.read()
        DFTB_content = re.sub(r"-N ", f"-N x{re.escape(base_name)}", DFTB_content)
        DFTB_content = DFTB_content.replace("input.inp", os.path.basename(input_file))
        DFTB_content = DFTB_content.replace("output.log", f"{base_name}.log")
        with open(DFTB_file, "w") as file:
            file.write(DFTB_content)
        os.chmod(DFTB_file, 0o755)
        subprocess.Popen(f"nohup ./{DFTB_file} &", shell=True)
    except Exception as e:
        print(f"An error occurred: {e}")
    finally:
        os.chdir("..")
