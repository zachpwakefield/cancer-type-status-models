"""Step 33: run the subtype jobs in parallel single-threaded subprocesses.

Usage: python scripts/33_subtype_run_all.py [--workers 6] [--cohorts BRCA SARC ...] [--views gex se ...]

Jobs: every view x every cohort in scope (results/subtype/data/cohort_scope.csv). Run order is
fixed in advance (D-36): the histology cohorts, then the genomic/IHC cohorts, then the
expression-defined ones; within a group, the largest cohort first. A job with job_done.json is
skipped, a fold with done.json is resumed. Per-job logs: logs/jobs_subtype/<view>__<cohort>.log;
status: logs/jobs_subtype/_status.tsv. The machine is shared: 6 workers at most by default.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time

import pandas as pd
from _common import REPO, log

from ctsm.subtype import TASK, subtype_config

GROUP = {"histology": 0, "genomic": 1, "IHC": 1, "expression-defined": 2}


def job_list(scfg, cohorts=None, views=None) -> list[tuple[str, str]]:
    scope = pd.read_csv(scfg.results_dir / TASK / "data" / "cohort_scope.csv")
    scope = scope[scope["in_scope"]]
    if cohorts:
        scope = scope[scope["cohort"].isin(cohorts)]
    scope = scope.assign(g=scope["label_type"].map(GROUP)).sort_values(["g", "labelled"], ascending=[True, False])
    return [(v, c) for c in scope["cohort"] for v in scfg.views if not views or v in views]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--cohorts", nargs="*")
    ap.add_argument("--views", nargs="*")
    args = ap.parse_args()
    scfg = subtype_config()
    logdir = REPO / "logs" / "jobs_subtype"
    logdir.mkdir(parents=True, exist_ok=True)
    pending = [j for j in job_list(scfg, args.cohorts, args.views)
               if not (scfg.run_dir(TASK, *j) / "job_done.json").exists()]
    log(f"{len(pending)} jobs to run with {args.workers} workers")
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               VECLIB_MAXIMUM_THREADS="1", PYTHONUNBUFFERED="1")
    running: dict[subprocess.Popen, tuple] = {}
    status_path = logdir / "_status.tsv"
    status = pd.read_csv(status_path, sep="\t").to_dict("records") if status_path.exists() else []
    n_before = len(status)
    while pending or running:
        while pending and len(running) < args.workers:
            view, cohort = pending.pop(0)
            fh = open(logdir / f"{view}__{cohort}.log", "a")
            p = subprocess.Popen([sys.executable, str(REPO / "scripts" / "32_subtype_run_job.py"), view, cohort],
                                 stdout=fh, stderr=subprocess.STDOUT, env=env, cwd=REPO)
            running[p] = (view, cohort, time.time(), fh)
            log(f"start {view}/{cohort} ({len(running)} running, {len(pending)} pending)")
        time.sleep(10)
        for p in [p for p in running if p.poll() is not None]:
            view, cohort, t0, fh = running.pop(p)
            fh.close()
            mins = (time.time() - t0) / 60
            state = "ok" if p.returncode == 0 else f"FAILED ({p.returncode})"
            status.append({"view": view, "cohort": cohort, "state": state, "minutes": round(mins, 1),
                           "finished": time.strftime("%Y-%m-%d %H:%M")})
            pd.DataFrame(status).to_csv(status_path, sep="\t", index=False)
            log(f"{state:>6} {view}/{cohort} in {mins:.1f} min")
    failed = [s for s in status[n_before:] if s["state"] != "ok"]
    log(f"all done; {len(failed)} failed")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
