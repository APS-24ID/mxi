// The unit cell from two-theta, as dials.two_theta_refine refines it: the
// cell alone, against the scattering angle of every integrated reflection's
// observed centroid -- which depends on the cell, the wavelength and where
// the spot fell on the detector, not on the crystal's orientation -- the
// cell constrained to its crystal system, its e.s.d.s from the fit's
// covariance. Cheap and, Graeme's word, more truthful than refinement's
// covariance, which fits the orientation, the detector and the cell
// together on the strong spots (docs/scaling.md, "The cell from two-theta").

#ifndef MXI_TWO_THETA_HH
#define MXI_TWO_THETA_HH

#include <array>
#include <string>
#include <vector>

#include "expt.hh"
#include "geometry.hh"
#include "refl.hh"
#include "symmetry.hh"

namespace mxi {

enum class CrystalSystem {
  triclinic,
  monoclinic,
  orthorhombic,
  tetragonal,
  trigonal, //: hexagonal axes, or rhombohedral, as the cell is given
  hexagonal,
  cubic
};

//: From the space group's number, 1 to 230.
CrystalSystem crystal_system_of(int space_group_number);
const char *name_of(CrystalSystem system);

struct TwoThetaObservation {
  Miller hkl{};
  double two_theta = 0.0; //: radians
  double variance = 0.0;  //: of two_theta, radians squared
  double wavelength = 0.0;
};

struct TwoThetaResult {
  bool valid = false;
  std::string why; //: when not valid
  UnitCell cell;
  //: a, b, c in A, alpha, beta, gamma in degrees; a parameter the system
  //: fixes or ties to another has its parent's or none.
  std::array<double, 6> esd{};
  std::size_t reflections = 0;
  int parameters = 0;
  double rmsd_deg = 0.0; //: of two-theta, unweighted
  CrystalSystem system = CrystalSystem::triclinic;
};

//: Refine `start` against the observations, constrained to `system`:
//: Gauss-Newton on sum w (2theta_obs - 2theta_calc)^2, w the inverse
//: variances; the e.s.d.s from (J^T W J)^-1 chi^2 / (n - p).
TwoThetaResult refine_two_theta(const std::vector<TwoThetaObservation> &obs,
                                const UnitCell &start, CrystalSystem system);

//: The observations from a scaled table -- each scaled reflection, not an
//: outlier, of I/sigma `min_isigma` or more, whose centroid has a variance --
//: two-theta from its centroid's place on its panel, the variance from the
//: centroid's, at least a pixel's own (size^2 / 12) on each axis; then one
//: cell refined over every crystal, starting from the
//: first's, in the space group's crystal system.
TwoThetaResult two_theta_refine(const ExperimentList &experiments,
                                const Table &reflections,
                                const SpaceGroup &group,
                                double min_isigma = 3.0);

} // namespace mxi

#endif // MXI_TWO_THETA_HH
