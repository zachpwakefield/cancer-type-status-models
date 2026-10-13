"""Step 00c (D-69): fetch the TCGA PanCancer Atlas MSI scores from cBioPortal, and verify them.

Approved by the data user on 2026-10-13 (about 5 MB of JSON, 64 requests, no token and no login).
Source: the public cBioPortal REST API, one request per study and attribute:
  https://www.cbioportal.org/api/studies/<study>_tcga_pan_can_atlas_2018/clinical-data
      ?clinicalDataType=SAMPLE&attributeId=<attribute>&projection=SUMMARY
  MSI_SCORE_MANTIS   MANTIS score (Bonneville et al. 2017)
  MSI_SENSOR_SCORE   MSIsensor score
Only sample ID, study and the two scores are kept, as one TSV sorted by sample, written to
<data_dir>/external/cbioportal_msi/ (git-ignored: per-sample values, D-12) with a manifest. The
TSV is pinned by byte count and MD5. A TSV already on disk is re-hashed and kept if it matches; a
mismatch stops the step. Each response is refused above 2 MB, the whole fetch above 10 MB.
"""
from __future__ import annotations

import hashlib
import io
import json
import time
import urllib.request

import pandas as pd
from _common import log

from ctsm.config import load_config

API = "https://www.cbioportal.org/api/studies/{study}_tcga_pan_can_atlas_2018/clinical-data" \
      "?clinicalDataType=SAMPLE&attributeId={attr}&projection=SUMMARY"
STUDIES = ["acc", "blca", "brca", "cesc", "chol", "coadread", "dlbc", "esca", "gbm", "hnsc", "kich", "kirc",
           "kirp", "laml", "lgg", "lihc", "luad", "lusc", "meso", "ov", "paad", "pcpg", "prad", "sarc", "skcm",
           "stad", "tgct", "thca", "thym", "ucec", "ucs", "uvm"]
ATTRS = {"MSI_SCORE_MANTIS": "mantis", "MSI_SENSOR_SCORE": "msisensor"}
NAME = "tcga_pancan_msi_scores.tsv"
PIN = {"bytes": 320_388, "md5": "81b5e01fa0547b18eadd71ecfdefcdb2"}
MAX_ONE, MAX_ALL = 2_000_000, 10_000_000


def get(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        body = r.read(MAX_ONE + 1)
    if len(body) > MAX_ONE:
        raise SystemExit(f"{url}: response above {MAX_ONE:,} bytes; stopped")
    return body


def main() -> None:
    out = load_config().data_dir / "external" / "cbioportal_msi"
    out.mkdir(parents=True, exist_ok=True)
    dest = out / NAME
    if not dest.exists():
        rows, total = [], 0
        for study in STUDIES:
            for attr, col in ATTRS.items():
                body = get(API.format(study=study, attr=attr))
                total += len(body)
                if total > MAX_ALL:
                    raise SystemExit(f"fetch above {MAX_ALL:,} bytes; stopped")
                for rec in json.loads(body):
                    rows.append({"sample_id": str(rec["sampleId"]), "study": study, "attr": col,
                                 "value": pd.to_numeric(rec.get("value"), errors="coerce")})
                time.sleep(0.2)
            log(f"{study}: {total / 1e6:.2f} MB so far")
        t = (pd.DataFrame(rows).pivot_table(index=["sample_id", "study"], columns="attr", values="value", aggfunc="first")
             .reindex(columns=list(ATTRS.values())).reset_index().sort_values(["sample_id", "study"]))
        buf = io.StringIO()
        t.to_csv(buf, sep="\t", index=False, float_format="%.6g")
        dest.write_text(buf.getvalue())
        log(f"fetched {len(rows):,} records ({total / 1e6:.2f} MB of JSON) -> {len(t):,} samples")
        manifest = {"source": API, "studies": STUDIES, "attributes": ATTRS, "records": len(rows),
                    "json_bytes": total, "fetched": time.strftime("%Y-%m-%d")}
        (out / "MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    data = dest.read_bytes()
    got = {"bytes": len(data), "md5": hashlib.md5(data).hexdigest()}
    if PIN["md5"] is None:
        log(f"{NAME}: {got['bytes']:,} bytes, md5 {got['md5']} (not yet pinned)")
    elif got != PIN:
        raise SystemExit(f"{NAME}: {got}; expected {PIN}")
    else:
        log(f"{NAME} verified")


if __name__ == "__main__":
    main()
