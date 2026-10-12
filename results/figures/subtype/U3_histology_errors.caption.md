**When a histology model is wrong, is the error molecularly coherent?**

Each panel is one histology cohort's expression (GEX) model. Rows: true histology → the model's out-of-fold call,
with patients in brackets; columns: the published molecular class of the same patients (PanCanAtlas, D-43); cells:
patients, coloured by the share of the row. Only patients with both labels are counted. Counts only, no test (D-41).
The AS ensemble's table is in the CSV (`plotted` = False); it shows the same pattern.

- **LGG.** 126 of 128 codeleted oligodendrogliomas are called oligodendroglioma, but 41 of the 58 without a
  codeletion are called astrocytoma; 6 of 6 codeleted astrocytomas are called oligodendroglioma.
- **UCEC.** 9 of the 12 endometrioid tumours called serous are copy-number high, against 28 of all 100 endometrioid
  tumours.
- **STAD.** 17 of 17 genomically stable diffuse tumours are called diffuse, and so are 7 of the 12 genomically
  stable intestinal ones.
- **PCPG.** 11 of the 11 adrenal tumours called extra-adrenal are in the pseudohypoxia cluster, and 22 of 22 tumours
  with cortical admixture are called adrenal.

The labels here were assigned by other groups over the full TCGA cohort; that is a property of the label, not a
leak, and no molecular label entered any histology model.
