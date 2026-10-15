#pragma once

// Scaling: grouping observations by symmetry, the merged intensities, and the
// least-squares fit of the scaling model (Beilsten-Edmands et al. 2020, eqns
// 2 and 3).

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "expt.hh"
#include "refl.hh"
#include "scale_model.hh"
#include "symmetry.hh"

namespace mxi {

enum class IntensityChoice { profile, summation };

struct ScaleDataOptions {
  IntensityChoice use = IntensityChoice::profile;
  //: dials.scale's default: a partial below 0.4 is too uncertain to scale up.
  double partiality_cutoff = 0.4;
  double d_min = 0.0;
  double d_max = 0.0; //: reflections with d above it left out; 0, none
  //: Friedel mates kept apart: for an acentric reflection I(+) and I(-) are
  //: separate groups, so scaling, the error model and outlier rejection do not
  //: take a real anomalous difference for error, as dials.scale's
  //: anomalous=True. Centric reflections stay one group, their mates equal by
  //: symmetry.
  bool anomalous = false;
};

//: The observations scaling uses, as parallel arrays.
struct ScaleData {
  //: Corrected as dials.scale corrects them before scaling: a profile-fitted
  //: intensity I lp / qe, a summed one I lp / (qe partiality), the variance
  //: by the square of the same factor; a reflection split into pieces by
  //: integration made one observation first (build_scale_data).
  std::vector<double> intensity, variance;
  //: Pieces of split reflections merged into another, as dials.scale's
  //: "Combined N partial reflections with other partial reflections".
  std::size_t partials_combined = 0;
  std::vector<ScaleObservation> observation;
  std::vector<std::size_t> group; //: index into `unique`
  std::vector<Miller> unique;     //: the symmetry-unique index of each group
  std::vector<std::size_t> row;   //: in the reflection table
  std::vector<double> d;
  //: A byte each, not std::vector<bool>: outlier rejection writes these from
  //: several threads at once, and neighbouring bits of one word are one write.
  std::vector<std::uint8_t> outlier;
  //: Both estimates, corrected likewise, for combining them; has_sum is false
  //: where summation failed, as it does across a module gap.
  std::vector<double> prf, prf_variance, sum, sum_variance;
  //: Each observation's variance with no signal, profile and summed, as
  //: integration wrote them (corrected as the variances are), and as
  //: combined; its correction factor; and, once a model exists, its intensity
  //: under the model, g I_h. The variance is then built on the model, not on
  //: the observation's own noise (model_variances).
  std::vector<double> prf_variance_bg, sum_variance_bg, variance_background;
  std::vector<double> factor, model_intensity;
  bool has_background = false;
  //: The intensity the error model's b term and bins take: the model's once
  //: there is one, else the observation's.
  double b_intensity(std::size_t i) const {
    return model_intensity.empty() ? intensity[i] : model_intensity[i];
  }
  std::vector<bool> has_sum;
  //: Which half of its Friedel pair each observation is (true for I+ and
  //: for centric reflections), and whether each group is centric.
  std::vector<bool> plus;
  std::vector<bool> centric;
  //: The variances as integration gave them, corrected, before an error model.
  std::vector<double> variance_before;
  //: I^2 var(g) / g^2, the uncertainty of the scale carried into the
  //: observation's own units: var(I / g) = (sigma^2 + this) / g^2. Zero until
  //: propagate_scale_variances sets it; the error model is refined and applied
  //: to sigma^2 plus this, so that a and b correct what remains after it.
  std::vector<double> scale_term;
  //: The term for observation i, zero where none was propagated -- including
  //: where the vector is shorter, as it is in data built by hand: reading it
  //: unguarded, the error model's tests ran off the end of an empty one.
  double scale_term_at(std::size_t i) const {
    return i < scale_term.size() ? scale_term[i] : 0.0;
  }
  // With anomalous, the Friedel pairs as well: each observation's pair,
  // numbered as the groups are without anomalous, so that the merging
  // statistics and the resolution estimate -- which split a pair into I(+) and
  // I(-) themselves -- report the same groups either way. Empty otherwise.
  std::vector<std::size_t> pair;
  std::vector<Miller> pair_unique;
  std::vector<bool> pair_centric;
  //: The grouping statistics are reported by: Friedel pairs together, whether
  //: or not scaling keeps them apart.
  const std::vector<std::size_t> &stats_group() const {
    return pair.empty() ? group : pair;
  }
  const std::vector<Miller> &stats_unique() const {
    return pair.empty() ? unique : pair_unique;
  }
  const std::vector<bool> &stats_centric() const {
    return pair.empty() ? centric : pair_centric;
  }
  std::size_t size() const { return intensity.size(); }
};

//: The observations of one sweep fit to scale: integrated by the chosen
//: method, with a positive variance, partiality at the cutoff or above, a
//: finite d beyond d_min, and not absent in the space group. Each carries its
//: position in the scan, 1/2d^2, and -- if the shape has an absorption term --
//: its harmonics at s1 and the reverse incident beam in the crystal frame.
ScaleData build_scale_data(const ExperimentList &experiments,
                           const Table &reflections, const SpaceGroup &group,
                           const ScaleModelShape &shape,
                           const ScaleDataOptions &options = {});

//: g for every observation.
std::vector<double> inverse_scales(const ScaleModel &model,
                                   const ScaleData &data);

//: <I_h> = sum w g I / sum w g^2 over each group's observations that are not
//: outliers, w = 1 / variance (the paper's eqn 3). A group with none is zero.
std::vector<double> merged_intensities(const ScaleData &data,
                                       const std::vector<double> &g);

//: A random subset of whole groups for fitting, as the paper's section 4.4:
//: groups are drawn until there are at least min_groups of them and
//: min_observations observations among them, or all are used.
std::vector<std::size_t>
select_for_fitting(const ScaleData &data, std::size_t min_groups = 2000,
                   std::size_t min_observations = 50000, unsigned seed = 1);

struct ScaleFitOptions {
  int max_iterations = 50;
  //: Weights of the restraints toward zero, sum B_i^2 and sum P_lm^2:
  //: dials.scale's decay_restraint and, for the absorption surface, its
  //: surface_weight at the default absorption level, low -- 5e5. It was 1000,
  //: 500 times weaker, for no recorded reason.
  double decay_restraint = 0.1;
  double absorption_restraint = 5e5;
};

struct ScaleFitResult {
  int iterations = 0;
  bool converged = false;
  double target_start = 0.0, target_end = 0.0;
  std::size_t observations = 0;
};

//: Fit the model by Levenberg-Marquardt over `use` (observation indices; all
//: that are not outliers if empty). The residuals are sqrt(w)(I - g <I_h>),
//: and are differentiated with <I_h> held: at the least-squares <I_h> the
//: target's derivative with respect to it is zero, so the gradient held is the
//: whole gradient and the fixed point is the true minimum. After each step the
//: scale is normalised, which a common factor with <I_h> leaves unseen.
ScaleFitResult fit_scale_model(ScaleModel &model, const ScaleData &data,
                               const ScaleFitOptions &options = {},
                               const std::vector<std::size_t> &use = {});

//: The covariance of the model's parameters at a fit, from the normal
//: equations of the variable-projection Jacobian and the restraints. Two
//: directions are fixed by the model's normalisation -- the scale's mean at one
//: and the relative B's at zero -- and the covariance is that of the parameters
//: under those constraints, Z (Z^T N Z)^-1 Z^T for Z spanning them: the first
//: of those directions is exactly null in N, and without the constraints
//: nothing can be inverted. Scaled by the goodness of fit, chi^2 over the
//: degrees of freedom, which count every merged <I> as a parameter fitted.
struct ParameterCovariance {
  std::vector<double> matrix; //: n x n, row major
  double goodness_of_fit = 1.0;
  std::size_t degrees_of_freedom = 0;
  bool ok = false;
};
ParameterCovariance parameter_covariance(const ScaleModel &model,
                                         const ScaleData &data,
                                         const ScaleFitOptions &options,
                                         const std::vector<std::size_t> &use);
//: var(g) for every observation: grad_g^T C grad_g.
std::vector<double>
inverse_scale_variances(const ScaleModel &model, const ScaleData &data,
                        const ParameterCovariance &covariance);
//: Set each observation's scale_term, I^2 var(g) / g^2.
void propagate_scale_variances(ScaleData &data, const std::vector<double> &g,
                               const std::vector<double> &g_variance);

//: Observations whose normalised deviation from their group's weighted mean,
//: excluding themselves, exceeds zmax (Evans 2006), flagged in data.outlier.
//: Every observation is retested, earlier outliers too; within a group the
//: worst is removed and the rest retested, so one cannot hide another. Returns
//: how many are flagged.
std::size_t reject_outliers(ScaleData &data, const std::vector<double> &g,
                            double zmax = 6.0);

struct ErrorModel {
  double a = 1.0, b = 0.02;
  //: Split (docs/scaling.md, "The error model, split"): a_bg for the variance
  //: with no signal, a for the rest -- the signal's and the scale's -- and b
  //: with a as before. Not split, a_bg is a.
  bool split = false;
  double a_bg = 1.0;
  std::size_t used = 0; //: observations it was determined from
};
//: sigma'^2 = a^2 (sigma^2 + (b I)^2), fitted as dials.scale's default: by
//: regression across bins of merged intensity, on groups with <I> above 25 and
//: <I / sigma^2> above 0.85. Then, where integration wrote each variance with
//: no signal, split: a_bg^2 V_bg + a^2 (V_signal + the scale's term) + a^2 b^2
//: I^2, fitted the same way over every group, the weak ones deciding a_bg --
//: kept if all three come out positive, else the model above.
//: `individual`: as dials.scale's default, a and b fitted in turn (a from the
//: normal probability plot, b from intensity bins); otherwise by regression
//: across intensity bins, with the split where integration wrote each
//: variance with no signal (docs/scaling.md).
ErrorModel refine_error_model(const ScaleData &data,
                              const std::vector<double> &g,
                              bool individual = false);
//: The error model applied to every observation's variance, from the variance
//: before any: applying it twice does not compound.
void apply_error_model(ScaleData &data, const ErrorModel &model);

//: Each observation's variance from the model rather than its own noise: its
//: variance with no signal, plus the signal's Poisson variance at the model's
//: intensity g I_h, corrected as the observation is. An observation's own
//: variance moves with its own intensity, so an inverse-variance weighted
//: mean leans toward the low ones: on bare background at a hundredth of a
//: count a pixel, 210 observations' at -1.9 sigma after integration's
//: variance from expected counts, and profile fitting's further
//: (docs/integration.md, "The fitted variance"). False, and nothing done, if
//: integration wrote no variances with no signal.
bool model_variances(ScaleData &data, const std::vector<double> &g);

//: How well the sigmas describe the scatter, bin by bin: each observation's
//: deviation from the mean of its equivalents without it, over the sigma of
//: the difference -- 1 where the sigmas are right -- as a mean square and
//: robustly, (1.4826 MAD)^2, so a few wild observations show as the two
//: apart. With the variances before the error model ("raw": integration's,
//: or the model's, and the scale's term) and after it ("final").
struct CalibrationRow {
  double lo = 0.0, hi = 0.0; //: the bin's range: d in A, or <I>
  std::size_t n = 0;
  double raw_ms = 0.0, raw_robust = 0.0, final_ms = 0.0, final_robust = 0.0;
};
std::vector<CalibrationRow> calibration(const ScaleData &data,
                                        const std::vector<double> &g,
                                        bool by_resolution, int bins);

//: I = w I_prf + (1 - w) I_sum with w = 1 / (1 + (I_sum / I_mid)^3), eqns 13
//: and 14; the variance as for fully correlated estimates, which they are,
//: being from the same pixels. I_mid of zero is profile fitting alone, and of
//: infinity summation alone.
void combine_intensities(ScaleData &data, double i_mid);
//: Rmeas with the scales given, over observations that are not outliers.
double rmeas(const ScaleData &data, const std::vector<double> &g);
//: The I_mid of lowest Rmeas, among profile alone, summation alone, and powers
//: of ten between; leaves data combined with it.
double choose_intensity_combination(ScaleData &data,
                                    const std::vector<double> &g);

//: One resolution shell of merging statistics, or all of them.
struct MergingShell {
  double d_max = 0.0, d_min = 0.0;
  std::size_t observations = 0, unique = 0, possible = 0, possible_acentric = 0;
  double multiplicity = 0.0, completeness = 0.0; //: completeness as a fraction
  double mean_i = 0.0;       //: mean over merged reflections of <I>
  double i_over_sigma = 0.0; //: mean over merged reflections of <I>/sigma(<I>)
  //: With Friedel mates merged, and with them apart ("I+/-").
  double rmerge = 0.0, rmeas = 0.0, rpim = 0.0;
  double rmerge_anom = 0.0, rmeas_anom = 0.0, rpim_anom = 0.0;
  double cc_half = 0.0;
  //: As iotbx.merging_statistics, which dials.scale reports: acentric
  //: reflections with both mates measured over the acentric ones possible; the
  //: observations over the groups with mates apart; the correlation of
  //: I+ - I- between random halves; the slope of the normal probability plot
  //: of dI / sigma(dI) over |x| < 0.9; sqrt(2 <(F+ - F-)^2> / <F+^2 + F-^2>);
  //: and mean |dI| over mean sigma(dI).
  double anom_completeness = 0.0, anom_multiplicity = 0.0, cc_anom = 0.0;
  double anom_slope = 0.0, df_over_f = 0.0, di_over_sig_di = 0.0;
  //: The control for the slope: the same slope over random halves of each
  //: pair's observations, I+ and I- pooled and dealt out again in the same
  //: numbers -- no anomalous signal, so near 1 where the sigmas are right.
  //: The slope above it, the signal; both above 1, the sigmas too small.
  double anom_slope_random = 0.0;
  std::size_t anomalous_pairs = 0;
};

//: Everything for the groups given except what needs the cell -- possible
//: reflections, completeness, the resolution limits -- over their observations
//: that are not outliers. The random halves are drawn with a fixed seed.
MergingShell merge_groups(const ScaleData &data, const std::vector<double> &g,
                          const std::vector<std::size_t> &groups);

//: merge_groups for each resolution shell -- equal volume in reciprocal
//: space, highest resolution LAST -- and overall, with completeness: every
//: reflection the crystal's cell allows to d_min, in the asymmetric unit and
//: not absent, and the acentric ones among them.
std::vector<MergingShell> merging_statistics(const ScaleData &data,
                                             const std::vector<double> &g,
                                             const SpaceGroup &group,
                                             const Crystal &crystal, int shells,
                                             MergingShell *overall);

struct ScaleRunOptions {
  //: Each variance from the model's intensity once there is one
  //: (model_variances); false keeps the integrated variances, as before.
  bool model_variances = true;
  bool anomalous = false; //: Friedel mates apart, as ScaleDataOptions
  bool combine = true;    //: choose profile, summation or a mix by Rmeas
  bool absorption = true; //: when the sweep is wide enough for it
  //: The absorption surface's degree: -1 for the sweep's default -- 4 from 60
  //: degrees, none below -- or 0 for none; lmax (lmax + 2) terms, 24 at 4 and
  //: 48 at 6.
  int lmax = -1;
  //: The degree a sweep wide enough for absorption takes when lmax is not
  //: given: 4 at dials.scale's absorption level low, 6 at medium and high.
  int level_lmax = 4;
  double d_min = 0.0;
  double d_max = 0.0; //: low resolution limit; 0, none
  ScaleFitOptions fit;
  //: The error model fitted as dials.scale's default does, a and b in turn,
  //: not by regression (refine_error_model).
  bool error_model_individual = false;
};

struct ScaleRun {
  ScaleData data;
  ScaleModel model{ScaleModelShape{}};
  std::vector<double> g;
  ErrorModel error_model;
  double i_mid = 0.0; //: 0 profile alone, infinity summation alone
  std::vector<ScaleFitResult> fits;
  std::vector<double> g_variance; //: var(g) for every observation
  ParameterCovariance covariance;
  std::size_t outliers = 0;
  //: Observations in the subset the model was fitted to.
  std::size_t fitted_on = 0;
  //: Wall seconds of each step, in order, for --timing.
  std::vector<std::pair<std::string, double>> timing;
};

//: Scale one sweep as the paper's figure 2: outliers on the unscaled data; fit;
//: outliers; the intensity combination; outliers; error model; fit; outliers;
//: error model; a final fit, outlier rejection and error model.
ScaleRun scale_sweep(const ExperimentList &experiments,
                     const Table &reflections, const SpaceGroup &group,
                     const ScaleRunOptions &options = {});

namespace flag {
//: dials.scale's, from dials/array_family/reflection_table.h, which is what
//: dials.merge selects by.
constexpr std::int64_t kOutlierInScaling = 1 << 23;
constexpr std::int64_t kExcludedForScaling = 1 << 24;
constexpr std::int64_t kScaled = 1 << 26;
} // namespace flag

//: The scaling written into the table as dials.scale writes it:
//: inverse_scale_factor and its variance, intensity.scale.value and variance
//: (corrected and combined, the error model's variance), and the flags --
//: scaled, an outlier in scaling, or excluded -- ADDED to those already there.
//: Rows not scaled get an inverse scale of 1 and a variance of 0. The
//: intensity's variance already carries the scale's, as dials.scale's does, so
//: nothing downstream should add inverse_scale_factor_variance to it again.
void write_scaling(Table &reflections, const ScaleData &data,
                   const std::vector<double> &g,
                   const std::vector<double> &g_variance);

} // namespace mxi
