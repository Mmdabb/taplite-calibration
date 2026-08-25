#ifndef TAPLITE_AUTO_CALIBRATION_FURTHER_DEVELOPMENT_H
#define TAPLITE_AUTO_CALIBRATION_FURTHER_DEVELOPMENT_H

#include "AutoCalibrationRefinement.h"

namespace taplite {

// Experimental copy of the refined fixed-volume inverse.  Unlike the
// refinement implementation, this model does not impose
// alpha = theta * Qcp * Qcd^s.  It estimates alpha independently so that the
// period-average speed can be fitted without disturbing the exact duration
// and trough-speed inverses.
FixedVolumeOracleResult SolveIndependentAlphaOracle(
    const FixedVolumeOracleInput& input,
    const FixedVolumeOracleConfig& config,
    double preferred_plf);

}  // namespace taplite

#endif  // TAPLITE_AUTO_CALIBRATION_FURTHER_DEVELOPMENT_H
