"""Reproduce the exact computations for the integrated manuscript."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory


HERE = Path(__file__).resolve().parent


def require(condition, description):
    if not condition:
        raise ValueError(description)


def data(path):
    return json.loads(path.read_text(encoding="utf-8"))


def mathematical_content(value):
    if isinstance(value, dict):
        return {k: mathematical_content(v) for k, v in value.items()
                if k != "seconds"}
    if isinstance(value, list):
        return [mathematical_content(v) for v in value]
    return value


def run(tree, group, job, output):
    name = job["id"]
    directory = tree / group
    expected_paths = job.get("outputs", [])
    for relative in expected_paths:
        (directory / relative).unlink(missing_ok=True)
    env = os.environ.copy()
    for key in ("PYTHONOPTIMIZE", "PYTHONPATH", "PYTHONHOME"):
        env.pop(key, None)
    env["PYTHONINTMAXSTRDIGITS"] = "0"
    command = ([sys.executable, "-E", "-B"] if job["runtime"] == "python" else ["node"])
    command += [job["script"], *job.get("args", [])]
    log = output / f"{name}.log"
    try:
        with log.open("wb") as stream:
            subprocess.run(command, cwd=directory, stdout=stream,
                           stderr=subprocess.STDOUT, env=env, check=True)
        for relative in expected_paths:
            generated = directory / relative
            reference = HERE / group / relative
            actual, expected = data(generated), data(reference)
            require(mathematical_content(actual) == mathematical_content(expected),
                    f"Certificate differs in {relative}")
        if job.get("compare_stdout"):
            reference = HERE / group / job["compare_stdout"]
            same = log.read_text(encoding="utf-8").splitlines()
            require(same == reference.read_text(encoding="utf-8").splitlines(),
                    f"Transcript differs for {name}")
        print(f"PASS {name}", flush=True)
        return {"id": name, "status": "passed", "log": log.name}
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        print(f"FAIL {name} ({error})", flush=True)
        return {"id": name, "status": "failed", "error": str(error), "log": log.name}


def execute(tree, output, mode, jobs):
    assignments = data(HERE / "checks.json")
    certificates = assignments["certificates"]
    subresonant = assignments["subresonant"]
    additional = assignments["additional"]
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {}

        def launch(group, job):
            futures[job["id"]] = pool.submit(run, tree, group, job, output)

        if mode == "subresonant":
            for job in subresonant:
                launch("subresonant", job)
        if mode == "quick":
            for job in certificates:
                if job.get("quick"):
                    launch("certificates", job)
            for job in subresonant:
                if job.get("quick"):
                    launch("subresonant", job)
            for job in additional:
                launch("certificates", job)
        if mode == "full":
            pending = {job["id"]: ("certificates", job)
                       for job in certificates + additional}
            pending.update({job["id"]: ("subresonant", job) for job in subresonant})
            while pending:
                ready = [name for name, (_, job) in pending.items()
                         if all(dep in futures for dep in job.get("after", []))]
                require(ready, "Dependency cycle in the checker catalogue")
                for name in ready:
                    group, job = pending.pop(name)

                    def run_after(group, job):
                        if all(futures[dep].result()["status"] == "passed"
                               for dep in job.get("after", [])):
                            return run(tree, group, job, output)
                        return {"id": job["id"], "status": "skipped",
                                "reason": "A prerequisite checker failed"}

                    futures[name] = pool.submit(run_after, group, job)
        return [future.result() for future in futures.values()]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("quick", "full", "subresonant"))
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    require(not sys.flags.optimize, "Run the checks with Python assertions enabled")
    with TemporaryDirectory(prefix="proof-check-", dir=HERE.parent) as temporary:
        fresh = Path(temporary)
        output = args.output_dir.resolve() if args.output_dir else fresh / "logs"
        output.mkdir(parents=True, exist_ok=False)
        for name in ("certificates", "subresonant"):
            shutil.copytree(HERE / name, fresh / name)
        results = execute(fresh, output, args.mode, args.jobs)
        passed = all(result["status"] == "passed" for result in results)
        report = {"mode": args.mode, "status": "passed" if passed else "failed",
                  "commands": results}
        (output / "report.json").write_text(json.dumps(report, indent=2) + "\n",
                                             encoding="utf-8")
        print(json.dumps({"status": report["status"], "commands": len(results)},
                         indent=2), flush=True)
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
