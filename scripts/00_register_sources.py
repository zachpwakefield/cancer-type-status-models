"""Step 00: record size and SHA-256 of every raw source (results/data/source_manifest.json).

A later change to any source is then detectable. Dataless (iCloud-evicted) files fail loudly.
"""
from __future__ import annotations

import os
from pathlib import Path

from _common import log

from ctsm.config import load_config, portable
from ctsm.provenance import sha256_file, write_json


def walk(node, prefix=""):
    if isinstance(node, dict):
        for k, v in node.items():
            yield from walk(v, f"{prefix}{k}.")
    else:
        yield prefix.rstrip("."), node


def main() -> None:
    cfg = load_config()
    paths = {f"tcga.{k}": cfg.tcga_path(*k.split(".")) for k, _ in walk(cfg.sources["tcga"])}
    paths["gencode_v45_gtf"] = cfg.gencode_gtf
    # subtype label sources (docs/SUBTYPE_METHODS.md, D-30): absolute paths
    paths.update({f"subtype_labels.{k}": Path(v) for k, v in cfg.sources.get("subtype_labels", {}).items()})
    out = {}
    for key, p in paths.items():
        st = os.stat(p)
        if getattr(st, "st_flags", 0) & 0x40000000:          # SF_DATALESS: evicted by iCloud
            raise SystemExit(f"{p} is an iCloud placeholder; download it first (brctl download)")
        log(f"hashing {key}: {p.name} ({st.st_size / 1e9:.2f} GB)")
        out[key] = {"path": portable(p), "bytes": st.st_size, "sha256": sha256_file(p)}   # ${CTSM_*}, not machine paths
    write_json(cfg.results_dir / "data" / "source_manifest.json", out)
    log(f"{len(out)} sources registered")


if __name__ == "__main__":
    main()
