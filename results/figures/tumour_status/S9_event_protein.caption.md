**Figure S9. The most shared events exchange protein segments through mutually exclusive and terminal exons.** The 9
recurrent events of S5 whose exon belongs to a protein-coding transcript (of ten, by number of models;
`11_event_cards/tumour_status_cards.csv` holds all ten with the rest of what is known about them). Each bar is one
protein, scaled to its own length (printed at the right), from its N terminus to its C terminus: grey blocks are
Pfam/CDD domains of the chosen transcript, the line underneath marks disordered regions (MobiDB-lite), and the
coloured block is the segment the event adds or removes -- blue where the segment is higher in tumours, orange where
it is lower (D-17, GENCODE v45). Three of the transcripts have no InterPro annotation at all, which is said so in
the margin: a bare bar there means the reference is silent, not that the segment misses every domain.

The annotation is what the swap does to the reading frame. For a mutually exclusive pair that is the difference
between the two exons, not either length, because the swap replaces one with the other. A first or last exon
replaces one end of the protein rather than shifting a frame, so it is annotated by the length of the segment it
contributes instead. A retained intron has no protein coordinates and is not drawn.

These are annotations of one transcript per gene, not isoform quantification: the models see the event, not the
protein, and no protein was measured here.
