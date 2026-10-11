"""Step 05: run every modelling job in parallel single-threaded subprocesses.

Usage: python scripts/05_run_all.py [--workers 8] [--max-pancan 5] [--task tumor_status|cancer_type]
Jobs: cancer_type x 9 views (PANCAN) and tumor_status x 9 views x each in-scope cancer.
Longest jobs start first. A job with job_done.json is skipped; a fold with done.json is
resumed. Per-job logs: logs/jobs/<task>__<view>__<cohort>.log; status: logs/jobs/_status.tsv
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time

import pandas as pd
from _common import REPO, log

from ctsm.config import PAN, load_config


def job_list(cfg, only_task: str | None) -> list[tuple[str, str, str, float]]:
    """(task, view, cohort, cost estimate) sorted by decreasing estimated cost."""
    scope = pd.read_csv(cfg.results_dir / "data" / "task_scope.csv")
    cohorts = scope.loc[scope["tumor_status_in_scope"], ["project", "tumour_cases", "normal_cases"]]
    jobs = []
    for view in cfg.views:
        if only_task in (None, "cancer_type"):
            jobs.append(("cancer_type", view, PAN, 1e9))
        if only_task in (None, "tumor_status"):
            for _, r in cohorts.iterrows():
                jobs.append(("tumor_status", view, r["project"], float(r["tumour_cases"] + r["normal_cases"])))
    return sorted(jobs, key=lambda j: -j[3])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--max-pancan", type=int, default=5, help="memory cap on concurrent PANCAN jobs")
    ap.add_argument("--task", choices=["tumor_status", "cancer_type"], default=None)
    args = ap.parse_args()
    cfg = load_config()
    py = sys.executable
    logdir = REPO / "logs" / "jobs"
    logdir.mkdir(parents=True, exist_ok=True)
    pending = [j for j in job_list(cfg, args.task) if not (cfg.run_dir(*j[:3]) / "job_done.json").exists()]
    log(f"{len(pending)} jobs to run with {args.workers} workers")
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               VECLIB_MAXIMUM_THREADS="1", PYTHONUNBUFFERED="1")
    running: dict[subprocess.Popen, tuple] = {}
    status = []
    while pending or running:
        n_pan = sum(1 for j in running.values() if j[2] == PAN)
        while pending and len(running) < args.workers:
            pick = next((i for i, j in enumerate(pending) if j[2] != PAN or n_pan < args.max_pancan), None)
            if pick is None:
                break
            job = pending.pop(pick)
            fh = open(logdir / f"{job[0]}__{job[1]}__{job[2]}.log", "a")
            p = subprocess.Popen([py, str(REPO / "scripts" / "04_run_job.py"), *job[:3]], stdout=fh,
                                 stderr=subprocess.STDOUT, env=env, cwd=REPO)
            running[p] = (*job[:3], time.time(), fh)
            n_pan += job[2] == PAN
            log(f"start {'/'.join(job[:3])} ({len(running)} running, {len(pending)} pending)")
        time.sleep(5)
        for p in [p for p in running if p.poll() is not None]:
            task, view, cohort, t0, fh = running.pop(p)
            fh.close()
            mins = (time.time() - t0) / 60
            state = "ok" if p.returncode == 0 else f"FAILED ({p.returncode})"
            status.append({"task": task, "view": view, "cohort": cohort, "state": state, "minutes": round(mins, 1)})
            pd.DataFrame(status).to_csv(logdir / "_status.tsv", sep="\t", index=False)
            log(f"{state:>6} {task}/{view}/{cohort} in {mins:.1f} min")
    failed = [s for s in status if s["state"] != "ok"]
    log(f"all done; {len(failed)} failed")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
