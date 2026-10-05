#pragma once

// The scaling model: g = C(r) . exp(B(t) / 2d^2) . S(s0, s1), the physical
// model of Beilsten-Edmands et al. (2020), parameterised with cubic B-splines
// for the scale C and relative B factor B, where DIALS averages the nearest
// parameters under a Gaussian, and the paper's spherical harmonics (eqn 7) for
// the absorption surface S.

#include <cstddef>
#include <vector>

#include "linalg.hh"

namespace mxi {

//: Real spherical harmonics of degree 1 to lmax at a unit vector, orthonormal
//: on the sphere, in the order l = 1..lmax and m = -l..l: lmax (lmax + 2) of
//: them. Without the Condon-Shortley sign, which would only flip the sign of a
//: fitted coefficient: degree 1 is sqrt(3 / 4 pi) (y, z, x) for m = -1, 0, 1.
void real_spherical_harmonics(int lmax, const Vec3 &unit,
                              std::vector<double> *out);
std::size_t harmonic_count(int lmax);

//: What the model needs of one observation.
struct ScaleObservation {
  double rotation = 0.0; //: position in the scan's rotation, 0 to 1, for C
  double time = 0.0;     //: position in the scan's time, 0 to 1, for B
  double inv_2d2 = 0.0;  //: 1 / (2 d^2)
  std::size_t sweep = 0; //: which experiment it was measured in, its id
  //: [Y(s1) + Y(s0)] / 2 in the crystal frame, harmonic_count(lmax) of them;
  //: empty without absorption.
  std::vector<double> absorption;
};

struct ScaleModelShape {
  std::size_t scale_points = 1;
  std::size_t decay_points = 0; //: none: no decay term
  int lmax = 0;                 //: 0: no absorption term
};

//: The number of control points DIALS' physical model uses for a sweep, in
//: effect: 6 and 5 for scale and decay on 30 degrees; for a sweep of 90 degrees
//: or more, one per 15 and 20 degrees and two more. Absorption from 60
//: degrees, as dials.scale's automatic choice.
ScaleModelShape default_shape(double degrees);

//: The physical model, for one sweep or several. Several sweeps' parameters
//: lie one sweep's block after another in `parameters`, each block a sweep's
//: scale points, decay points and harmonics; an observation's `sweep` chooses
//: its block. The accessors without a sweep are the first's, which is all of
//: it for one sweep.
class ScaleModel {
public:
  explicit ScaleModel(const ScaleModelShape &shape);
  //: One block a sweep, each its own shape -- its own width of rotation.
  explicit ScaleModel(const std::vector<ScaleModelShape> &shapes);

  std::size_t sweeps() const { return shapes_.size(); }
  const ScaleModelShape &shape(std::size_t sweep = 0) const {
    return shapes_[sweep];
  }
  std::size_t size() const { return parameters.size(); }
  std::vector<double> parameters;
  //: Where a sweep's block, its decay points and its harmonics begin.
  std::size_t first_scale(std::size_t sweep = 0) const {
    return offsets_[sweep];
  }
  std::size_t first_decay(std::size_t sweep = 0) const {
    return offsets_[sweep] + shapes_[sweep].scale_points;
  }
  std::size_t first_absorption(std::size_t sweep = 0) const {
    return first_decay(sweep) + shapes_[sweep].decay_points;
  }

  double inverse_scale(
      const ScaleObservation &o,
      std::vector<std::pair<std::size_t, double>> *gradient = nullptr) const;

  //: The scale's mean over every sweep's points made one, and the decay's
  //: mean over every sweep's taken away: the one overall scale and the one
  //: overall B the merged intensities can absorb. Not each sweep's apart --
  //: that would erase the sweeps' relative scale, which is what scaling them
  //: together is for.
  void normalise();

private:
  std::vector<ScaleModelShape> shapes_;
  std::vector<std::size_t> offsets_;
};

} // namespace mxi
