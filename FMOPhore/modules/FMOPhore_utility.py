import os
import socket
import sys

class EnvironmentGuard:
    _logged = False  # Class-level flag

    def __init__(self, log_file="FMOPhore.log"):
        self.log_file = log_file
        self.allowed_hosts = [
            "hpc-ddu-gpu-001.compute.dundee.ac.uk",
            "hpc-ddu-gpu-002.compute.dundee.ac.uk",
            "hpc-ddu-gpu-003.compute.dundee.ac.uk",
            "hpc-ddu-gpu-004.compute.dundee.ac.uk",
            "hpc-ddu-gpu-005.compute.dundee.ac.uk",
            "hpc-ddu-gpu-006.compute.dundee.ac.uk",
            "hpc-ddu-gpu-007.compute.dundee.ac.uk",
            "hpc-ddu-gpu-008.compute.dundee.ac.uk",
            "gpu-25.compute.dundee.ac.uk",
            "gpu-26.compute.dundee.ac.uk",
            "gpu-27.compute.dundee.ac.uk",
            "gpu-28.compute.dundee.ac.uk",
            "gpu-29.compute.dundee.ac.uk",
            "gpu-30.compute.dundee.ac.uk",
            "gpu-31.cluster.lifesci.dundee.ac.uk",
            "gpu-28.cluster.lifesci.dundee.ac.uk"
        ]

    def log_message(self, message):
        with open(self.log_file, "a") as log_file:
            log_file.write(f"{message}\n")

    def check_qsub_environment(self):
        required_env_vars = ["JOB_ID", "QUEUE"]
        if not all(var in os.environ for var in required_env_vars):
            self.log_message("🚫 FMOPhore must be run within a valid qsub job. Aborting!")
            sys.exit(1)

    def check_host_authorization(self):
        actual_hostname = socket.getfqdn()
        if not EnvironmentGuard._logged:
            self.log_message(f"Current hostname: {actual_hostname}")
            EnvironmentGuard._logged = True

        if actual_hostname not in self.allowed_hosts:
            self.log_message("🚫 FMOPhore is restricted to authorized compute nodes. Aborting!")
            sys.exit(1)

    def enforce(self):
        self.check_qsub_environment()
        self.check_host_authorization()
