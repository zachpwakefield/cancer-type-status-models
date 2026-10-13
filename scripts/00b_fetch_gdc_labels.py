"""Step 00b: fetch the open-access GDC tables the subtype tasks and their follow-ups use, and verify them.

All come from GDC publication pages (no token, no controlled data):
  TCGASubtype.20170308.tsv                    PanImmune page; per-patient molecular subtype calls
                                              collected for the PanCanAtlas (TCGAbiolinks table)
  clinical_PANCAN_patient_with_followup.tsv   PanCanAtlas page; per-patient clinical fields
                                              (BRCA ER/PR/HER2 IHC, HNSC HPV, PRAD Gleason, grade)
  added 2026-10-12 for the BRCA follow-up (D-57..), with the data user's approval:
  PANCAN_ArmCallsAndAneuploidyScore_092817.txt  arm-level gain/loss calls (Taylor et al. 2018)
  TCGA_mastercalls.abs_tables_JSedit.fixed.txt  ABSOLUTE purity and ploidy
  ABSOLUTE_scores.tsv                           aneuploidy and LOH scores
  TCGA_mastercalls.abs_segtabs.fixed.txt        ABSOLUTE allele-specific segments (LOH per segment)
  mc3.v0.2.8.PUBLIC.maf.gz                      MC3 public somatic mutation calls
Each file is pinned by UUID, byte count and MD5. A file already on disk is re-hashed and kept if
it matches; a mismatch stops the step instead of overwriting. Downloads refuse anything larger
than the pinned size. Written to <data_dir>/external/gdc_pancanatlas/ (git-ignored: the tables
are per-patient), with a manifest of what was fetched and when.
"""
from __future__ import annotations

import hashlib
import json
import time
import urllib.request

from _common import log

from ctsm.config import load_config

API = "https://api.gdc.cancer.gov/data/"
FILES = {
    "TCGASubtype.20170308.tsv": {
        "uuid": "0f31b768-7f67-4fc4-abc3-06ac5bd90bf0", "bytes": 603_836,
        "md5": "d0a6054b799a8a61488636ad400c3af9",
        "page": "https://gdc.cancer.gov/about-data/publications/panimmune"},
    "clinical_PANCAN_patient_with_followup.tsv": {
        "uuid": "0fc78496-818b-4896-bd83-52db1f533c5c", "bytes": 18_633_685,
        "md5": "ffcb35edda305dd8d615497f9214eb92",
        "page": "https://gdc.cancer.gov/about-data/publications/pancanatlas"},
    "PANCAN_ArmCallsAndAneuploidyScore_092817.txt": {
        "uuid": "4c35f34f-b0f3-4891-8794-4840dd748aad", "bytes": 1_079_849,
        "md5": "a610bbd79b25e09e3cda61b99706e0cf",
        "page": "https://gdc.cancer.gov/about-data/publications/PanCan-CellOfOrigin"},
    "TCGA_mastercalls.abs_tables_JSedit.fixed.txt": {
        "uuid": "4f277128-f793-4354-a13d-30cc7fe9f6b5", "bytes": 901_812,
        "md5": "8ea2ca92c8ae58350538999dfa1174da",
        "page": "https://gdc.cancer.gov/about-data/publications/pancanatlas"},
    "ABSOLUTE_scores.tsv": {
        "uuid": "0e8831f4-dd7e-4673-8624-b4519c2e0d65", "bytes": 456_006,
        "md5": "bef23328ca806e81f3591973319e3605",
        "page": "https://gdc.cancer.gov/about-data/publications/panimmune"},
    "TCGA_mastercalls.abs_segtabs.fixed.txt": {
        "uuid": "0f4f5701-7b61-41ae-bda9-2805d1ca9781", "bytes": 253_061_161,
        "md5": "585c8793730f294d7bf0144566bb37fa",
        "page": "https://gdc.cancer.gov/about-data/publications/pancanatlas"},
    "mc3.v0.2.8.PUBLIC.maf.gz": {
        "uuid": "1c8cfe5f-e52d-41ba-94da-f15ea1337efc", "bytes": 753_339_089,
        "md5": "639ad8f8386e98dacc22e439188aa8fa",
        "page": "https://gdc.cancer.gov/about-data/publications/pancanatlas"},
}


def md5_file(path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest, max_bytes: int) -> None:
    tmp = dest.with_suffix(dest.suffix + ".part")
    n = 0
    with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as fh:
        while chunk := r.read(1 << 20):
            n += len(chunk)
            if n > max_bytes:
                tmp.unlink()
                raise SystemExit(f"{dest.name}: more than the pinned {max_bytes:,} bytes; stopped")
            fh.write(chunk)
    tmp.rename(dest)


def main() -> None:
    cfg = load_config()
    out = cfg.data_dir / "external" / "gdc_pancanatlas"
    out.mkdir(parents=True, exist_ok=True)
    manifest = {}
    for name, f in FILES.items():
        dest = out / name
        if not dest.exists():
            log(f"fetching {name} ({f['bytes'] / 1e6:.1f} MB)")
            fetch(API + f["uuid"], dest, f["bytes"])
        got = md5_file(dest)
        if dest.stat().st_size != f["bytes"] or got != f["md5"]:
            raise SystemExit(f"{name}: {dest.stat().st_size:,} bytes, md5 {got}; "
                             f"expected {f['bytes']:,} bytes, md5 {f['md5']}")
        log(f"{name} verified")
        manifest[name] = {**f, "url": API + f["uuid"], "verified": time.strftime("%Y-%m-%d")}
    (out / "MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
