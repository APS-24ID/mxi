# The algorithms, stage by stage

What each program does, by what method, from whom, and where it lives.
Where mxi follows DIALS the method is DIALS's and the reference is DIALS's
source as much as the paper; where it departs, the departure is named and the
reason is in the document given. `docs/review.md` says the same for a
reviewer, as differences from DIALS; this is for someone who wants to know how
it works.

## Import -- `mxi_import` (`docs/import.md`)

**Reading NXmx** (HDF5): the beam, the detector's panels with their
transformations and parallax, the goniometer's axes and the scan, as
dials.import reads them, into one experiment a sweep; the detector's own pixel
mask with it. `--split N` cuts a sweep into N for testing.

## Spot finding -- `mxi_find` (`docs/spots.md`)

**The extended dispersion threshold**, DIALS's default (`src/dext*`): a
pixel is signal where the local index of dispersion and the local mean say so,
over windows of 7 x 7, 5 x 5 and 11 x 11, on the CPU, CUDA or Metal; then
**connected components** over images, the spots' centroids and variances.
Frames are read straight from their compressed chunks (bitshuffle-LZ4) by
pread, outside HDF5's lock.

## Indexing -- `mxi_index`

**The three-dimensional FFT** (Bricogne 1986), as DIALS's fft3d: the spots'
reciprocal-space positions sampled onto a grid a step of 1 / (2.5 max_cell),
transformed, the largest peaks taken as candidate basis vectors, the largest
cell estimated as DIALS's find_max_cell does. The basis chosen as dials.index's
filtering ranker chooses; a **Buerger reduction** of the cell; a macrocycle of
assignment and refinement. Several lattices (`--max-lattices`) found in turn,
each spot assigned again to the lattice it lies nearest
(`docs/multi-lattice.md`).

## Refinement -- `mxi_refine`

**Weighted least squares on the spot centroids**, x and y on the detector and
the rotation angle, by **Gauss-Newton with Levenberg damping**
(`src/refine.hh`), with **analytical derivatives** for crystal, detector and
beam, checked against finite differences. **Scan-varying** crystal as a cubic
B-spline of control points, one every 10 degrees by default -- DIALS uses a
Gaussian smoother; both smooth and local, not the same model. The Jacobian is
sparse and threaded.

**Prediction** (`src/predict.hh`): a reciprocal lattice point carried round
the rotation axis crosses the Ewald sphere where |r0|^2 + 2 s0 . R(phi) r0 = 0,
which about the axis is B cos(phi) + C sin(phi) + A0 = 0: closed form, one atan2
and one acos, two roots, no iteration.

## Integration -- `mxi_integrate` (`docs/integration.md`)

* **The profile model** after Kabsch (2010): the spot's spread in the
  diffracted beam's direction, sigma_b, from each refined spot's pixels about
  its predicted s1; the reflecting range, sigma_m, by maximum likelihood over
  the images holding foreground -- as DIALS's defaults (`basic`; or weighted by
  counts, `--sigma-m-estimator extended`). The foreground spans 3 sigma in
  Kabsch's (eps1, eps2, eps3) frame.
* **Shoeboxes and masks**: each reflection's box, its foreground, and every
  neighbouring reflection's foreground -- any lattice's -- masked out of its
  background, as DIALS's 3d_threaded integrator does.
* **The background**: the **robust Poisson GLM** of Parkhurst et al. (2016),
  Huber's weights, as DIALS's; fitted over a histogram of the background's
  counts where they are photon counts.
* **Summation** after Leslie (1999), the Poisson term from the counts
  expected.
* **Reference profiles** after Kabsch (2010): strong spots transformed into
  Kabsch's frame and summed on a grid, 13 a side, by detector region and scan
  block of 10 degrees.
* **Profile fitting** against the shoebox's pixels -- the reference profile
  interpolated between neighbouring cells and carried onto the pixels, not on
  the grid -- by iteratively reweighted least squares, on the CPU, CUDA or
  Metal; one pass over the images.

## Symmetry -- `mxi_symmetry` (`docs/symmetry.md`)

* **Centred input** taken to its primitive cell first, by the standard change
  of basis for its centring.
* **Several crystals to one setting** by **cosym** (Gildea & Winter 2018): the
  correlation of every pair of data sets under every rotation of the lattice,
  the matrix's leading eigenvectors (Householder and QL), each data set's
  operation the one nearest the reference.
* **The lattice's symmetry** by **Le Page's method** (1982), through gemmi's
  implementation after Zwart et al. (2006).
* **Normalisation**: an anisotropic Wilson model with a light atom's
  scattering.
* **Each symmetry element and each subgroup scored** after Evans (2011),
  likelihoods from the correlations of reflections the element relates, to a
  resolution limit from CC half.
* **The space group by its absences**, then the data reindexed into it.

## Scaling -- `mxi_scale` (`docs/scaling.md`)

* **The physical model** of Beilsten-Edmands et al. (2020), as dials.scale's:
  a scale C(phi) and a decay B(phi), each a cubic B-spline, and absorption by
  spherical harmonics in the crystal's frame.
* **The fit**: Levenberg-Marquardt with the **variable-projection Jacobian**,
  the merged intensities eliminated; the normal equations a group of
  equivalents at a time; the scales' covariance from one factorisation.
* **Partials**: a reflection split by integration's blocks combined first;
  partiality divides the summed intensity only.
* **Profile, summation or a mix**, chosen by Rmeas as dials.scale chooses.
* **Outlier rejection**: an observation against the mean of its group without
  it, above 6 sigma (Evans 2006), repeated.
* **The error model** sigma'^2 = a^2 [sigma^2 + (b I)^2]: by regression across
  intensity bins, split into a for the variance with no signal and a for the
  rest where integration records it; or as dials.scale's default, a from the
  normal probability plot and b from the bins in turn (`--error-model-fit
  individual`). The scales' uncertainty carried into each variance before it.
* **The resolution limit** as dials.estimate_resolution's, CC half fitted.
* **The cell's e.s.d.s from two-theta**, as dials.two_theta_refine.

## Clustering -- `mxi_correlation` (`docs/correlation.md`)

As dials.correlation_matrix (Thompson et al. 2025):

* **The correlation of every pair of data sets**, sigma-weighted as DIALS's
  weighted_cchalf, each pair weighted by **Kish's effective sample size**.
* **Cosym's coordinates** minimising sum w_ij (r_ij - x_i . x_j)^2 by L-BFGS
  with the curvatures, from DIALS's own start; **the dimension** at the elbow of
  the functional over 1 to 50, as distl's resolution method (Zhang et al.
  2006); the coordinates' principal components.
* **OPTICS** (Ankerst et al. 1999) with **xi** cluster extraction, ported from
  scikit-learn; min_samples chosen by the **silhouette** with a noise penalty
  over every candidate, where DIALS samples them by optuna; DIALS's
  single-cluster heuristic.
* **Ward's hierarchical clustering** for the report's dendrograms, by the
  Lance-Williams update, as scipy's.

## Export -- `mxi_export` (`docs/export.md`)

**An unmerged MTZ** as dials.export writes it: the intensities and their
corrections, the columns, a batch an image with its header.

## Reading the papers

* Ankerst, Breunig, Kriegel & Sander (1999). OPTICS. *SIGMOD '99*.
* Beilsten-Edmands et al. (2020). *Acta Cryst.* D76, 385-399. dials.scale.
* Bricogne (1986). The FFT in indexing.
* Evans (2006). *Acta Cryst.* D62, 72-82. Scaling and outliers.
* Evans (2011). *Acta Cryst.* D67, 282-292. Symmetry scoring.
* Gildea & Winter (2018). *Acta Cryst.* D74, 405-410. cosym.
* Kabsch (2010). *Acta Cryst.* D66, 133-144. Profile model and integration.
* Le Page (1982). *J. Appl. Cryst.* 15, 255-259. Lattice symmetry.
* Leslie (1999). *Acta Cryst.* D55, 1696-1702. Summation.
* Parkhurst et al. (2016). *J. Appl. Cryst.* 49, 1912-1921. The background
  GLM.
* Thompson et al. (2025). *Acta Cryst.* D81, 278-290. Intensity clustering.
* Zhang et al. (2006). *J. Appl. Cryst.* 39, 112-119. distl.
* Zwart, Grosse-Kunstleve & Adams (2006). Lattice symmetry, as gemmi has it.
