#ifndef EICD_HELPERS_HH
#define EICD_HELPERS_HH


#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <string>
#include <vector>

#include <Math/Vector4D.h>

#include "eicd/TrackParametersCollection.h"
#include "eicd/ReconstructedParticleCollection.h"
#include "eicd/ReconstructedParticleData.h"

namespace eicd::helpers {

  /** Four momentum from track and mass.
   * Get a vector of 4-momenta from raw tracking info, using an externally
   * provided particle mass assumption. 
   */
  inline auto momenta_from_tracking(const std::vector<eic::TrackParametersData>& tracks,
                                    const double                                 mass)
  {
    std::vector<ROOT::Math::PxPyPzMVector> momenta{tracks.size()};
    // transform our raw tracker info into proper 4-momenta
    std::transform(tracks.begin(), tracks.end(), momenta.begin(), [mass](const auto& track) {
      // make sure we don't divide by zero
      if (fabs(track.qOverP) < 1e-9) {
        return ROOT::Math::PxPyPzMVector{};
      }
      const double p  = fabs(1. / track.qOverP);
      const double px = p * cos(track.direction.phi) * sin(track.direction.theta);
      const double py = p * sin(track.direction.phi) * sin(track.direction.theta);
      const double pz = p * cos(track.direction.theta);
      return ROOT::Math::PxPyPzMVector{px, py, pz, mass};
    });
    return momenta;
  }
} // namespace eicd::helpers
#endif
