import yaml
from subprocess import run
from functools import partial

shell = partial(run, capture_output=True, shell=True)

def load_config(file: str):
    with open(file, "r") as stream:
        return yaml.safe_load(stream)

def submit(file: str, condition:str = None, previous = None):
    cmd = f"sbatch {file}" if condition is None else f"sbatch --dependency={condition}:{previous} {file}"
    exe = shell(cmd)

    if exe.returncode:
        raise Exception(f"Error while submitting job {file}...\n\n{exe.stderr}")
    
    jobid = exe.stdout.decode("UTF-8").strip().split()[-1]
    print(f"Submitted job {file} with jobid: {jobid}")
    
    return jobid

def bind_dependency(jobid: str, dependency: str, condition: str = "afterok"):
    cmd = f"scontrol update jobid={jobid} dependency={condition}:{dependency}"

    exe = shell(cmd)

    if exe.returncode:
        raise Exception(f"Error while updating job {jobid}...\n\n{exe.stderr}")

def main(args):
    config = load_config(args.config)

    jobs = config["Jobs"]
    condition = config["Condition"]

    previous_jobid = submit(jobs[0])

    for job in jobs[1:]:
        previous_jobid = submit(job, condition, previous_jobid)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("--config", help="Config file")

    args = parser.parse_args()

    main(args)