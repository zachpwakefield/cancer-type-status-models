**Can a single splicing event name a lineage?**

One event per class, picked in two steps (D-26): of the five events with the largest mean |SHAP| in that class's
own explanation -- across all 8 splicing views, so this is the models' own ranking -- the one whose class median is
furthest from the nearest other class. The ten classes with the clearest such event are drawn, most separated
first. The second step chooses which examples to show and nothing else; a feature can rank first for a class and
still not separate it alone, because the model reads it together with the others.

Each row is one event: the filled dot is the median over the samples of the class it marks, the open dot the
nearest other class, and the small dots the remaining 32, spread vertically so that classes sharing a value can
be counted. The gap between the two large dots is the margin a single event gives on its own, with no model and
no other feature.

Medians over all samples of each class, out-of-fold membership not required (this is a description of the data,
not of a prediction). The CSV holds every drawn dot, one row per event × class. Per-sample values stay local
(D-12); `11_event_cards/cancer_type_cards.csv` adds each event's geometry, protein consequence and regulon
membership.
