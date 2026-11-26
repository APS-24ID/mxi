# mxi

```
░▒▓ █▀▄▀█ ▀▄▀ █
░▒▓ █ ▀ █ ▄▀▄ █
```

mxi is a data processing package for NXmx format macromolecular crystallography
(MX) data: from the diffraction images to scaled intensities ready for structure
solution -- spot finding, indexing, refinement, integration, symmetry, scaling,
merging, the clustering of data from many crystals, and what the scaled data say
of twinning. Small molecules too.

**Almost instant data processing.** mxi aims to process data as DIALS does, with
results of similar quality, an order of magnitude more quickly -- fast enough to
keep pace with modern detectors, and to give results while the crystal is still
mounted. Spot finding and integration can run on a GPU, through CUDA or Apple's
Metal.

**A new, standalone package.** mxi is written from scratch in C++. It contains no
DIALS code and does not need DIALS installed, but it follows DIALS' methods,
reads and writes DIALS' files, and any step can be exchanged for DIALS' own.

**NXmx data.** [NXmx](https://journals.iucr.org/m/issues/2020/05/00/ti5018/index.html)
is the community standard format for X-ray diffraction images, written directly
by modern detectors such as DECTRIS's EIGER; older formats can be converted, CBF
for example by [cbf2nxmx](https://github.com/graeme-winter/cbf2nxmx).

## Installing

```sh
pip install mxi          # on the CPU; Metal on Apple silicon
pip install mxi-cuda     # Linux with an NVIDIA GPU
```

Install **one or the other** in an environment, not both: each brings the same
`mxi` command. Wheels are built for Linux on x86_64 and aarch64 (glibc 2.28 or
later: RHEL 8 and 9 and their clones, Ubuntu 20.04 on, Debian 11 on), macOS on
Apple silicon (13 or later) and Windows on x86_64 -- not macOS on Intel. Python
3.10 or later. `mxi-cuda` needs the NVIDIA driver for CUDA 12; which GPUs each
build supports is in the
[README](https://github.com/graeme-winter/mxi/blob/main/README.md#installing).

## Using it

```sh
mxi-ui                                    # the interface, in a directory of data
mxi process master.nxs -j 8 --d-min-auto  # every step in turn, for routine data
```

Or step by step:

```sh
mxi import    master.nxs        # the experiment list            -> imported.rflx
mxi find      imported.rflx     # strong spots on every image    -> strong.rflx
mxi index     strong.rflx       # the lattice, by 3D FFT         -> indexed.rflx
mxi refine    indexed.rflx      # beam, detector, crystal        -> refined.rflx
mxi integrate refined.rflx      # summation and profile fitting  -> integrated.rflx
mxi symmetry  integrated.rflx   # the Laue and space group       -> symmetrized.rflx
mxi scale     symmetrized.rflx  # scaled, with merging statistics -> scaled.rflx
mxi export    scaled.rflx       # an unmerged MTZ                -> scaled.mtz
mxi merge     scaled.rflx       # a merged MTZ, with amplitudes  -> merged.mtz
```

`--gpu` puts spot finding and integration on the GPU. Every command says what it
takes with `--help`; all of them are in the
[command reference](https://github.com/graeme-winter/mxi/blob/main/docs/commands.md).

## What the package holds

- `mxi`, the programs, compiled for each platform: the steps above and the tools
  around them -- `mxi show`, `mxi report`, `mxi convert` to and from DIALS'
  files, `mxi combine` and `mxi split`, `mxi correlation` for data from many
  crystals, `mxi triage`, and more.
- `mxi-ui`, the interface: each step as a form, building and running the real
  command, with its log, summary and plots, and each step's images in a viewer.
- `mxi-view`, the image viewer and mask editor.
- `mxeq`, the referee: comparisons between runs of a pipeline, reflection by
  reflection, mxi's or DIALS'.

## Licence

mxi's source is under the BSD 3-Clause licence. The programs in these wheels are
built with FFTW, which is under the GPL (version 2 or later), so the programs as
distributed here are subject to the GPL's terms as a whole; the wheels carry
every notice for what they contain, gemmi's (MPL-2.0), HDF5's, bitshuffle's,
LZ4's, zlib's and mimalloc's among them, and in `mxi-cuda` NVIDIA's terms for
the CUDA runtime.

## More

- [README](https://github.com/graeme-winter/mxi/blob/main/README.md): installing, building, using it
- [How it works](https://github.com/graeme-winter/mxi/blob/main/docs/algorithms.md): each stage's methods and sources
- [How it differs from DIALS, and why](https://github.com/graeme-winter/mxi/blob/main/docs/review.md)
- [A tutorial](https://github.com/graeme-winter/mxi/blob/main/docs/tutorials/cows-pigs-people.md): many crystals, and which belong together
- [Release notes](https://github.com/graeme-winter/mxi/blob/main/RELEASE_NOTES.md)
- [The source](https://github.com/graeme-winter/mxi)
