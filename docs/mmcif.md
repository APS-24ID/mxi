# Export: PDBx/mmCIF, merged and unmerged

STATUS: plan, 6 November 2026; Graeme's decisions below. `docs/outstanding.md`,
28.

`mxi export scaled.rflx -o scaled.cif` to write a PDBx/mmCIF file that a wwPDB
deposition takes: the scaled observations, unmerged, and the data merged,
with what mxi knows of how they were measured and processed -- and nothing it
does not know.

## Decisions (Graeme, 6 November)

1. **The new dialect**: `_pdbx_diffrn_unmerged_refln` and its categories, not
   the older `_diffrn_refln`.
2. **A merged block** beside the unmerged one.
3. **Scaled data only**: no export of integrated, unscaled data as mmCIF.
4. **The beamline's metadata** -- source, beamline, detector type, how they
   are found or given -- deferred, Graeme to decide later. Until then, what
   the models say and nothing more.

## The dialect, as it stands

The new categories are in the wwPDB's development dictionary,
`mmcif_pdbx_v5_next`, not in its release (V5), and the dictionary marks them
used in no PDB entry yet: a target that may move. They form a whole, "data
sections" for merged and unmerged data:

* `_pdbx_diffrn_data_section` -- each section's id, scattering, whether
  merged, whether scaled; and its lists, `_pdbx_diffrn_merge_crystal_list`,
  `_pdbx_diffrn_merge_image_list`, `_pdbx_diffrn_merge_rejected_list`.
* `_pdbx_diffrn_unmerged_refln` -- an observation a row: indices (which must
  be in the cell of `_pdbx_diffrn_unmerged_cell`), scan, first and last image,
  the scaled intensity and its sigma, the scale applied, the profile-fitted and
  summed intensities, partiality, scan angle.
* `_pdbx_diffrn_unmerged_cell`, `_pdbx_diffrn_scan`, `_pdbx_diffrn_symmetry`.
* `_pdbx_diffrn_merged_refln`, `_pdbx_diffrn_merged_cell` -- the merged data.
* `_pdbx_diffrn_merge_stat`, `_pdbx_diffrn_merge_stat_shell` -- the merging
  statistics.

`dials.export format=mmcif` writes this dialect as its `pdb_version=v5_next`,
the unmerged observations and scans and cells, with the statistics in the
established `_reflns` and `_reflns_shell` (`dials/util/export_mmcif.py`, read
for this): the working reference, as dials.export's MTZ writer was for
`mxi export`'s MTZ. The established `_reflns` and `_reflns_shell` are written
here too, beside `_pdbx_diffrn_merge_stat`: OneDep takes a structure's
data-collection table from them.

## What goes in, and where mxi has it

| Content | Categories | From |
| --- | --- | --- |
| each observation | `_pdbx_diffrn_unmerged_refln` | scaled.rflx's table: indices, experiment's scan and images, intensity.scale and inverse_scale_factor, intensity.prf, intensity.sum, partiality, xyzcal's angle |
| the merged data | `_pdbx_diffrn_merged_refln`, `_pdbx_diffrn_merged_cell` | merged as scaling merges for its statistics |
| the statistics | `_pdbx_diffrn_merge_stat(_shell)`, `_reflns`, `_reflns_shell` | scaling's summary record, as scaling found them -- not recomputed |
| Wilson B | `_reflns.B_iso_Wilson_estimate` | the ML fit mxi symmetry and mxi triage use (laue.hh) |
| scans, images, angles | `_pdbx_diffrn_scan`, `_diffrn_scan(_axis)` | each sweep's scan |
| crystals and cells | `_pdbx_diffrn_unmerged_cell`, `_cell`, `_exptl_crystal` | each crystal's cell; the median for `_cell` |
| space group | `_pdbx_diffrn_symmetry`, `_symmetry` | the crystals' group |
| wavelengths | `_diffrn_radiation_wavelength`, `_diffrn_source.pdbx_wavelength_list` | the beams |
| the processing | `_software`, `_audit`, `_pdbx_diffrn_data_section.details` | mxi's version, the summary's steps |

Rejected observations -- scaling's outliers -- in
`_pdbx_diffrn_merge_rejected_list`, not left out without trace.

## How

In C++, `mxi export`'s `-o NAME.cif` (`.cif.gz` compressed), through gemmi's
CIF document and writer (`third_party/gemmi`: cifdoc.hpp, to_cif.hpp,
header-only): not through gemmi's MTZ-to-mmCIF converter, which would keep only
what an MTZ keeps. Some 100 MB uncompressed for a million observations: the
writing streamed, the loops not held twice.

## Phases

1. **Unmerged**: the unmerged data section -- observations, scans, cells,
   symmetry, wavelengths, rejected observations -- the statistics, Wilson B,
   software.
2. **Merged**: the merged data section beside it, intensities. Amplitudes
   (French and Wilson) not here: mxi has none yet.
3. **The beamline's metadata**: when Graeme decides (above).

## How it will be judged

* **Against the dictionary**: every item in `mmcif_pdbx_v5_next`, its type and
  its mandatory items, by gemmi's validator.
* **Read back**: by gemmi and by cctbx's iotbx.cif; the merging statistics
  recomputed from the file's observations against scaling's.
* **Beside dials.export** (`format=mmcif pdb_version=v5_next`) on the same
  data: category by category, every difference explained.
* **A deposition** test, OneDep's, when Graeme has a structure to try it with.

## Open

* The dictionary's development version may change before release: which
  version a file was written against goes in `_audit_conform`.
* Which of the merged statistics OneDep takes from the new categories, and
  which from `_reflns`: both written until that is known.
