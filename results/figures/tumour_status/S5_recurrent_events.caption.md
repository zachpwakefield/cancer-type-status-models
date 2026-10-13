**Figure S5. Different cancers rely on the same events, and those events move the same way across cancers.** Across
the 8 AS views, 404 events rank in the top 100 (mean |SHAP|) of at least 4 of the 15 tumour-vs-normal models. 3.1
would be expected if the rankings were independent (Poisson-binomial over each event's measured cohorts;
`recurrence_vs_chance.csv`).

The figure shows the 40 most recurrent events. Rows give the host gene (HGNC), view and coordinates (hg38, 1-based).
Colour is the median paired tumour − normal change in every cancer, including cancers where the event did not rank
(so the pattern across cancers is not selected on): blue where the event is higher in the tumour, orange where it is
higher in the normal tissue, clipped at ±30. Dots mark the cancers whose top 100 contain the event.

Recurrent events tend to move the same way across cancers. For the median event, 73% of the cancers where it is
measured share its majority direction; 46% of events reach 80%.

Several genes appear as switch pairs, with two exons of the same gene moving in opposite directions (EXOSC7 last
exons, NDE1 first exons, TPM1 mutually exclusive exons).

Recurrent RI and A5SS events mostly rise in tumours (73% and 77%), matching the global shift in S2. The full list
with direction consistency is in `recurrent_events.csv`.
