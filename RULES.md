# World-building rules

## Rule 1. Scientific data is the source of truth

Before any decision about the map is made (climate, continent outlines, resources or anything else), it is checked
whether the decision can be justified scientifically.

## Rule 2. Every change to canon is recorded

Every time the map departs from canon, an entry is added to `CANON_DEVIATIONS.md` that explains why this particular
decision was made.

## Rule 3. Requests that contradict the rules

If any request from the project owner contradicts the accepted rules, the work stops and the owner is told about it
immediately. If we decide to break a rule this time, an entry is added to `RULES_DEVIATIONS.md`.

## Rule 4. Layers are generated, not stored

The map is built from layers, but the layers themselves are not kept in the repository. Each layer is produced by a
script that generates it. Every layer has its own document that describes the generation process in detail. Every run
of the procedure with the same parameters must produce the same result.

## Rule 5. What is kept for each layer

For each layer the repository keeps:

1. the research data;
2. every experiment, with its input data, results and conclusions. The results of an experiment are an assessment,
   not generated layers: whether the experiment reached its goal, whether its results satisfy us and, if not, why.
   Intermediate layers are not kept (see Rule 4);
3. once the procedure that generates the layer is ready, a document that briefly summarises the results of the work.
