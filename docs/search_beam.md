# mxi search-beam: the beam centre from the spots

`mxi search-beam` finds where the beam meets the detector from the strong spots
alone, as `dials.search_beam_position` does by default -- LABELIT's method
(Sauter, Grosse-Kunstleve & Adams 2004, J. Appl. Cryst. 37, 399-409). For when
the header's centre is too far out for indexing to find the lattice: indexing
refines the centre once it has a lattice, and not before.

```sh
mxi find imported.rflx
mxi search-beam strong.rflx          # writes optimised.rflx
mxi index optimised.rflx
```

It prints the centre as given and as found, in pixels and millimetres, and the
shift, and writes the experiments with the detector moved, the spots unchanged
(`-o` for another name, `.expt` and `.refl` for the pair). Options as DIALS's:
`--scope` (4 mm either way), `--binning` (the grid's step, 2 pixels), `--d-min`,
`--max-cell`, `--max-reflections` (10000).

## How it works

With the centre right, the reciprocal lattice's planes pass through the origin,
so along each of the lattice's strong periodicities t every spot's reciprocal
lattice point r projects to an integer, r.t, and the phase of S = mean exp(2 pi
i r.t) is zero. A wrong centre displaces the lattice from the origin, and the
phases move off zero. So the strongest periodicities are found once, from the
centre as given; the detector is moved across the beam over a grid, the spots
mapped again for each trial and scored; and the best refined.

**The periodicities** are the three-dimensional FFT's candidate vectors, as `mxi
index` finds them (DIALS uses LABELIT's one-dimensional search for the same
thing), sharing indexing's scale -- the resolution, the largest cell, the grid.
But only the coherent: on 30 degrees of cubic insulin 4 of the strongest 20
candidates were the lattice's, |S| 0.95 at the true centre, the rest 0.01, and
two short ones at 0.07 and 0.13 pointing the wrong way, which cancelled half the
signal. So those with |S| at least 40 per cent of the best's, with their sums
and differences -- lattice vectors too, more directions -- up to 20, none within
5 degrees of a better.

**The score**, at each trial: for each periodicity, its vector refitted to the
spots -- three steps of Gauss-Newton flattening the phases about their mean --
and Re S, |S| times the cosine of its phase; summed. Refitted because the FFT's
vectors are only as fine as its grid, and a direction slightly off scatters the
phases, which a small move of the detector partly mends: even on an exact
lattice the best score lay 0.1 mm off the truth until the vectors were refitted.
|S| as well as the phase -- DIALS sums the cosine alone -- because of aliases,
below.

**The search**: a grid over the scope, its step the pixel times the binning;
then every local maximum within 90 per cent of the best, at most five, refined
by a pattern search to a thousandth of a millimetre, and the highest refined
score taken -- the nearest the start only among those within half a per cent of
it.

## Aliases

A detector moved by D lambda |g| shifts the reciprocal lattice by one of its own
vectors g, and every spot lies on a lattice point again: for insulin's 78 A cell
at 170 mm that is about 3 mm, for a 200 A cell some 1.2 mm. Only the Ewald
sphere's curvature, which makes the shift uneven, tells such an alias from the
truth. DIALS takes, of the scores within 90 per cent of the best, the nearest
the centre as given -- on insulin started 1.8 mm out that was an alias. Here the
refined vectors and |S| in the score make the truth score clearly higher: on
insulin, centres given 1.8, 2.4 and 4.1 mm out all came back within 3 pixels of
the refined centre, the last beyond half the alias spacing. An alias may still
win where the data are few or the cell long; `mxi index` will then fail or index
badly, and the scope can be narrowed (`--scope`) about a better guess.

## Tested

`tests/test_beam_search.cc`: on an exact lattice, a cubic 78 A crystal with
insulin's geometry, its reflections predicted and taken as the spots -- from the
true centre the detector moves by no more than 0.05 mm; displaced by 1.0 and 0.7
mm across the beam, it comes back to within 0.05 mm. On real data (30 degrees of
insulin, by hand): the four starts above, each within 3 pixels; on one core some
30 to 55 s, the 3000 trials of the grid most of it.
