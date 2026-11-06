// mxi triage: what the merged intensities say of twinning, as
// phenix.xtriage (Zwart, Grosse-Kunstleve & Adams, CCP4 newsletters 42 and
// 43, 2005) says it, a start -- the statistics that do not depend on twin
// laws (the moments, N(Z), the L test) and those that do (the twin laws the
// lattice allows, each one's R and H-test fraction), and a summary that, as
// xtriage's, calls the data twinned only where both kinds say so (docs/
// triage.md). Translational NCS, the Britton plot, Wilson and the data's
// health are to come.
//
// On the data as they are given: merged, in their space group, Friedel
// mates merged; the cell regularised to the group.

#ifndef MXI_TRIAGE_HH
#define MXI_TRIAGE_HH

#include <string>
#include <vector>

#include "geometry.hh"
#include "symmetry.hh"

namespace mxi {

//: Merged intensities in the asymmetric unit, one each.
struct TriageData {
  std::vector<Miller> hkl;
  std::vector<double> intensity, sigma;
};

struct TriageOptions {
  double d_max = 10.0;        //: low resolution of the twin tests, as xtriage's
  double d_min = 0.0;         //: high resolution; 0 for xtriage's choice
  double isigi_cut = 3.0;     //: in that choice: I/sigma over this ...
  double completeness = 0.85; //: ... in over this fraction of a shell
  double resolution_at_least = 3.5; //: the choice never coarser than this
  int strength_bins = 15;           //: shells of the choice, equal volume
  double max_delta = 3.0;           //: obliquity of twin laws (degrees)
  unsigned seed = 1;                //: the L test's choice of neighbours
};

//: The moments of the normalised intensities, E^2 = I / (eps <I/eps>), and
//: their values untwinned and perfectly twinned (xtriage's table).
struct Moments {
  std::size_t n = 0;
  double i2_over_i_squared = 0; //: <I^2>/<I>^2: 2 acentric, 3 centric
  double f_squared_over_f2 = 0; //: <F>^2/<F^2>: 0.785, 0.637
  double e2_minus_one = 0;      //: <|E^2 - 1|>: 0.736, 0.968
};

//: The cumulative distribution of E^2, N(z), z = 0.1 .. 1.0, observed and
//: Wilson's -- 1 - exp(-z) acentric, erf(sqrt(z / 2)) centric.
struct CumulativeIntensity {
  std::vector<double> z, acentric, acentric_theory, centric, centric_theory;
  double max_deviation_acentric = 0, max_deviation_centric = 0;
};

//: Padilla & Yeates' L test (2003), as xtriage's: L = (E1^2 - E2^2) /
//: (E1^2 + E2^2) for each positive acentric reflection and one neighbour
//: drawn for it, h + (2h', 2k', 2l'), |dh|+|dk|+|dl| 2, 4 or 6 -- none where
//: that neighbour was not measured. Untwinned <|L|> = 1/2, <L^2> = 1/3;
//: perfectly twinned 3/8, 1/5. The multivariate Z score is xtriage's: the
//: Mahalanobis distance of (<|L|>, <L^2>) from those of some 1200 PDB data
//: sets, (0.4851, 0.3164), by their covariance -- how far from typical good
//: data, whatever their number, not a test of the ideal values.
struct LTest {
  std::size_t pairs = 0;
  double mean_abs_l = 0, mean_l2 = 0, z_score = 0;
};

//: A twin law: a rotation of the lattice outside the point group, one of its
//: coset; merohedral where it keeps the cell's metric exactly. R_obs =
//: sum |I1 - I2| / sum (I1 + I2) over twin-related pairs; the H test's twin
//: fraction (Yeates 1997) from the strongest half of the pairs, as xtriage:
//: H = (I1 - I2) / (I1 + I2), <|H|> = 1/2 - alpha. Both on the intensities
//: as merged, xtriage's: R over every reflection with a mate, centric too,
//: sum |I1 - I2| / sum |I1 + I2|; H over the acentric.
struct TwinLaw {
  Rotation op{};
  std::string text; //: as xtriage writes it, "-h,k,-l"
  bool merohedral = false;
  double r_metric = 0; //: sum (G' - G)^2 / sum G^2, the cell's metric
  std::size_t pairs = 0;
  double r_obs = 0, h_alpha = 0;
};

struct Triage {
  double d_max = 0, d_min = 0; //: the twin tests' range
  //: As xtriage prepares the data before the twin tests: put on an absolute
  //: scale with the fall-off of an anisotropic Wilson model removed (laue.hh,
  //: its f2(s) left in), over d*^2 0.008 to 0.690; and outliers removed.
  bool sharpened = false;
  std::size_t outliers = 0;
  std::size_t used = 0, acentric = 0, centric = 0;
  Moments acentric_moments, centric_moments;
  CumulativeIntensity nz;
  LTest l_test;
  std::vector<TwinLaw> twin_laws;
  //: xtriage's abnormal intensity statistics: Z at least 3.5 and <|L|> below
  //: 1/2, as twinning moves it.
  bool abnormal = false;
  //: Z at least 3.5 the other way, <|L|> 1/2 or more: a spread wider than
  //: good data's, which twinning never makes -- translational NCS can.
  bool wide = false;
  std::string verdict; //: xtriage's summary, in words
};

//: epsilon: how many of the point group's rotations leave h as it is.
int epsilon(const SpaceGroup &group, const Miller &hkl);

//: The twin laws `cell`'s lattice allows within max_delta, beyond `group`.
std::vector<TwinLaw> twin_laws(const UnitCell &cell, const SpaceGroup &group,
                               double max_delta);

//: xtriage's high resolution limit for the twin tests (data_statistics.py,
//: i_sigi_completeness_stats): the data's range in 15 shells of equal volume,
//: the last left out; of those finer than resolution_at_least, the finest in
//: which more than `completeness` of the reflections possible have I/sigma
//: over isigi_cut -- or resolution_at_least where none is, or it is coarser.
double twin_test_d_min(const TriageData &data, const UnitCell &cell,
                       const SpaceGroup &group, const TriageOptions &o);

Triage triage(const TriageData &data, const UnitCell &cell,
              const SpaceGroup &group, const TriageOptions &options = {});

//: A rotation on Miller indices, h' = R^T h, as xtriage writes a twin law.
std::string twin_law_text(const Rotation &r);

} // namespace mxi

#endif // MXI_TRIAGE_HH
