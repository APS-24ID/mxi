# mxi_import

STATUS: first version, 2 October 2026; several masters, 5 October. Right on
masters from three writers -- DECTRIS's for insulin (the 4M cut of ins10_1),
redhorn-archive's for thaumatin (thau_1_3), Diamond I03's for threonine
(thr_7_1, thr_7_3, two sweeps) -- its goniometer, axes, sensor, mu, trusted
range and scan identical to dials.import's on insulin, and on threonine every
model of both experiments but the trusted range and the exposure time. To be made right on the population of
data sets one at a time, `mxeq compare-expt` saying what differs.

`mxi_import master.nxs` writes `imported.expt`, as `dials.import master.nxs`
does, so that nothing of DIALS is needed for the chain. `mxi_import a.nxs
b.nxs` writes one experiment a master, as dials.import writes several sweeps:
each with a beam, detector, goniometer, scan and image set of its own, the
experiments' indices numbered past those before, each its own identifier. The
overrides apply to every master alike.

## What it reads

* **The beam:** `incident_wavelength` (instrument's or sample's NXbeam), or
  `incident_energy`; a wavelength that varies through the scan is its first
  value, said so. Direction along +z in DIALS's frame (the beam travelling
  along -z), and dials.import's polarisation defaults: fraction 0.999, normal
  +y.
* **The detector:** a panel for each NXdetector_module, named by its path.
  The module's fast and slow pixel directions give the pixel size (their
  values) and directions (their vectors); everything they depend on places the
  module. A chain is followed from its dataset through each `depends_on`, a
  relative path resolved from the dataset's group, to `"."`; each step is
  `R p + offset` for a rotation and `p + offset + value vector` for a
  translation, innermost first, at the first image's value. An offset without
  `offset_units` is in the transformation's own `units`, as nxmx reads it:
  Diamond's Eiger masters write offsets so, in metres -- taken as millimetres
  they put the detector a thousandth of the way out. A rotation's offset with
  none is taken as millimetres, said so. Then from NeXus's
  McStas frame to DIALS's imgCIF, a half turn about y: x and z negated.
  `data_size` and `data_origin`, slow then fast, are the image size and raw
  offset, fast then slow.
* **The sensor:** `sensor_material` (Silicon, CdTe, GaAs, Ge to their symbols),
  `sensor_thickness`, and mu -- not in the file -- from Hubbell and Seltzer's NIST
  mass attenuation coefficients, interpolated log-log in energy, times the
  density: on insulin's 0.953738 A, 3.663092965478474 per mm, dials.import's to
  the last digit. Only silicon is tabulated; another material has mu 0, no
  parallax correction, and says so -- `--mu`.
* **The trusted range:** 0 to `saturation_value`, or the DECTRIS
  `countrate_correction_count_cutoff`, or failing both just below the bad-pixel
  marker, said so -- `--trusted-max`.
* **The goniometer:** the rotations of `/entry/sample/depends_on`'s chain,
  innermost first, written as the file gives their vectors (not normalised, as
  dials.import writes them); the one whose values change the scan axis at angle
  0, the others at their values. More than one moving, or none, is a failure.
* **The scan:** one image a value of the scan axis, images 1 to n; the
  oscillation those values; exposure time `count_time`; epochs 0.
* **The image set:** an ImageSequence of the master, by absolute path, and a new
  identifier.

Every unit is read from its attribute -- lengths to mm, angles to degrees -- and
one missing is taken as mm or degrees and said.

## Overrides

`--wavelength A`; `--distance MM`, moving the detector along its normal;
`--beam-centre X,Y`, pixels fast and slow, moving it in its plane so that the
beam meets it there; `--mu`; `--trusted-max`; `--image-range A,B`. Each override
used is said.

## Against dials.import

    dials.import master.nxs output.experiments=dials.expt
    mxi_import master.nxs
    mxeq compare-expt imported.expt dials.expt

prints each model's parts, "same" or the difference: origins in mm, axes in
degrees, numbers relatively. On insulin, against dials.import of the full 16M
master, the differences were the origin and image size -- a different file --
and the exposure time, 0.0026 s here where dials.import wrote 0; to settle on a
matched pair.

On threonine, two sweeps from Diamond's I03, against dials.import of both
masters, `mxeq compare-expt` comparing experiment by experiment: every model
the same in each but two. The exposure time again, 0.01 s against 0. And the
trusted range's top, 65534 here against dials.import's 133201: the master has no
saturation value, from which dxtbx's NeXus reader would take 0x7FFFFFFF, so the
133201 came from elsewhere -- most likely the detector's `_meta.h5` beside the
master, which dxtbx's Diamond format finds. Open until that file is seen; give
`--trusted-max` meanwhile, as data of 32 bits are wrongly distrusted at 65534.

## Untested

A detector of several modules (each becomes a panel; nothing downstream has
been run on one); a moving detector; a beam direction or polarisation from the
file; scans about any axis but the outermost's; materials other than silicon.
