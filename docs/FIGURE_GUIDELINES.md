# Figure guidelines

Every figure in this repo follows the lab's house style. Its source is the figure guidelines
of the sister project, `survival_prediction/docs/figure_guidelines/`: the README with sha256
prefix `60f6fa0294a1`, together with its three example figures and style code. Read that
README and open its example figures before you change a figure. They are not copied here.

## Where the style lives

[`src/ctsm/figstyle.py`](../src/ctsm/figstyle.py) carries the house style. Every figure script
imports it.

| From | Taken |
|---|---|
| `figstyle.py` (`caad6e7c5f50`) | rcParams: Arial at 6.5–8 pt; 0.6 pt axes; no top or right spines; no grid; no legend frames; `pdf.fonttype` 42 and `svg.fonttype` `none`; 400-dpi PNG |
| `fig_style_snhg.py` (`512ac8891be9`) | Ink greys, tumour `#2a78d6` and normal `#eb6834`, the blue tiers, the diverging map, and the `fp()` p-value formatter |

**Colour roles.**
- In tumour–normal figures, tumour is blue and normal is orange.
- Elsewhere, blue `#2a78d6` is the point of the figure, such as splicing or the AS ensemble.
- Dark blue `#184f95` is a second accent, such as fusion or "observed".
- Grey or ink marks the reference, such as GEX, single views, chance or background.
- The diverging map marks signed changes. The sequential tier blues mark magnitudes.

## Rules applied here

- **Titles.** No in-plot titles. A small heading is allowed only where it names the subject of
  a panel, such as a cancer, a class pair or an event.
- **Captions.** Explanations go in `<stem>.caption.md` next to the figure. A caption opens with
  the question the figure answers, in bold, and then gives the method, n, the test and where
  the values are.
- **Exports.** Final figures are written as SVG, PDF and a 400-dpi PNG, plus `<stem>.csv` with
  the exact plotted values. A separate `<stem>_tests.csv` holds the tests when they are not
  part of the plotted table.
- **Individual-level data (D-12).** When a figure plots individual samples or patients, its
  CSV holds the aggregates. The per-sample values go to the git-ignored
  `<data_dir>/figure_data/<stem>.csv`.
- **Statistics.** Figures show raw p as values or stars. The tiers are * < 0.05, ** < 0.01 and
  *** < 0.001, and each caption repeats them. q appears only in the CSV, with its FDR family
  named in a `q_family` column. A missing test is shown as `–`. A measured zero is printed as
  `0`, because a zero-length bar is invisible and would otherwise read as a missing value.
- **Gene labels.** Genes appear as HGNC symbols. A gene without one is shown by its locus
  (`chr6:41.78 Mb`, from `_figbase.locus_label`), never as an Ensembl ID, which stays in the
  CSV. Two rows of one figure never carry the same label: where one gene contributes two
  events, the label adds what tells them apart, such as the exon's length.
- **Fewer plots.** Each figure answers one question. Diagnostics are kept as tables, not
  figures; examples are the error-overlap Jaccard and the tumour-status fusion.
- **Colour.** Blue `#2a78d6` is tumour and orange `#eb6834` is normal, on every figure where
  the contrast is tumour against normal tissue, including the colour scale of a tumour − normal
  difference (`figstyle.TN_DIVERGE`). `figstyle.DIVERGE` (blue below, red above) is for signed
  quantities that are not that contrast, such as an odds ratio or a hazard ratio.
- **Where they live.** The figures of the interpretation layer are grouped by narrative, not by
  the analysis that produced them: `results/figures/tumour_status/S1`–`S10`
  ([TUMOUR_STATUS.md](../results/TUMOUR_STATUS.md)) and `results/figures/cancer_type/T1`–`T9`
  ([CANCER_TYPE.md](../results/CANCER_TYPE.md)), drawn by `scripts/20` and `scripts/21`. The
  tables they read stay under `results/biology/<NN>_<analysis>/`. A caption refers only to
  figures of its own narrative.
- **Per-job panels.** These are `results/<task>/<view>/<cohort>/fig_*.png`, one per job, 153
  in all. They are PNG only, to keep the repository small, and each has a caption file. Their
  values are the job's SHAP and metrics tables.

## Before you call a figure done

1. Look at the PNG. Check for overlapping or clipped labels, text that would be unreadable at
   print size, and legends covering data.
2. Check that the CSV reproduces exactly what is plotted.
