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
  double isigi_cut = 3.0;     //: in that choice: I/sigma at least this ...
  double completeness = 0.85; //: ... in at least this fraction of a shell
  double coarsest_d_min = 3.5; //: the choice never coarser than this
  double max_delta = 3.0;      //: obliquity of twin laws (degrees)
  unsigned seed = 1;           //: the L test's choice of neighbours
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

//: Padilla & Yeates' L test (2003): L = (I1 - I2) / (I1 + I2) for acentric
//: reflections and a neighbour h + (2h', 2k', 2l'), 2 <= |dh|+|dk|+|dl| <= 8.
//: Untwinned <|L|> = 1/2, <L^2> = 1/3; perfectly twinned 3/8, 1/5. The
//: multivariate Z score: the Mahalanobis distance of (<|L|>, <L^2>) from the
//: untwinned values, L uniform on [-1, 1] -- over 3.5, abnormal, as xtriage.
struct LTest {
  std::size_t pairs = 0;
  double mean_abs_l = 0, mean_l2 = 0, z_score = 0;
};

//: A twin law: a rotation of the lattice outside the point group, one of its
//: coset; merohedral where it keeps the cell's metric exactly. R_obs =
//: sum |I1 - I2| / sum (I1 + I2) over twin-related pairs; the H test's twin
//: fraction (Yeates 1997) from the strongest half of the pairs, as xtriage:
//: H = (I1 - I2) / (I1 + I2), <|H|> = 1/2 - alpha.
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
  std::size_t used = 0, acentric = 0, centric = 0;
  Moments acentric_moments, centric_moments;
  CumulativeIntensity nz;
  LTest l_test;
  std::vector<TwinLaw> twin_laws;
  bool abnormal = false; //: the L test's Z over 3.5, either way
  //: Abnormal as twinning makes it: <|L|> below 1/2. Above, the spread is
  //: wider than untwinned data's -- not twinning: the zones a higher
  //: symmetry makes centric, taken as acentric, or translational NCS. xtriage's
  //: Z is two-sided; its own example of a symmetry missed has <|L|> 0.512.
  bool twin_like = false;
  std::string verdict; //: xtriage's summary, in words
};

//: epsilon: how many of the point group's rotations leave h as it is.
int epsilon(const SpaceGroup &group, const Miller &hkl);

//: The twin laws `cell`'s lattice allows within max_delta, beyond `group`.
std::vector<TwinLaw> twin_laws(const UnitCell &cell, const SpaceGroup &group,
                               double max_delta);

//: xtriage's high resolution limit for the twin tests: the finest d, from
//: d_max down, of shells in which at least `completeness` of the reflections
//: possible have I/sigma at least isigi_cut -- never coarser than
//: coarsest_d_min, never finer than the data.
double twin_test_d_min(const TriageData &data, const UnitCell &cell,
                       const SpaceGroup &group, const TriageOptions &o);

Triage triage(const TriageData &data, const UnitCell &cell,
              const SpaceGroup &group, const TriageOptions &options = {});

//: A rotation on Miller indices, h' = R^T h, as xtriage writes a twin law.
std::string twin_law_text(const Rotation &r);

} // namespace mxi

#endif // MXI_TRIAGE_HH
