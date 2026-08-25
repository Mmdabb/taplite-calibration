#include "AutoCalibrationRefinement.h"
#include "AutoCalibrationFurtherDevelopment.h"

#include <cassert>
#include <cmath>
#include <iostream>

namespace {

bool Near(double left, double right, double tolerance = 1e-7) {
    return std::fabs(left - right) <= tolerance;
}

taplite::FixedVolumeOracleInput ExactCase() {
    taplite::FixedVolumeOracleInput input;
    input.volume = 2400.0;
    input.lanes = 3.0;
    input.period_hours = 3.0;
    input.capacity_vphpl = 2000.0;
    input.free_speed_mph = 65.0;
    input.cutoff_speed_mph = 45.0;
    input.qn = 1.0;
    input.qs = 4.0;
    input.observed_duration_hour = 2.0;
    input.observed_trough_speed_mph = 30.0;
    const double theta = 8.0 / 15.0;
    const double queue_factor = 1.0 + theta * (45.0 / 30.0 - 1.0);
    const double reference_speed = 55.0;
    const double queue_speed = reference_speed / queue_factor;
    input.observed_average_speed_mph =
        (2.0 / 3.0) * queue_speed + (1.0 / 3.0) * (reference_speed + 65.0) / 2.0;
    return input;
}

}  // namespace

int main() {
    taplite::FixedVolumeOracleConfig config;
    config.use_bounds = false;
    const taplite::FixedVolumeOracleInput input = ExactCase();
    taplite::FixedVolumeOracleResult exact =
        taplite::SolveFixedVolumeOracle(input, config);
    assert(exact.status == taplite::OracleStatus::ExactFeasible);
    assert(exact.exact_feasible);
    assert(Near(exact.doc_raw, 0.5));
    assert(Near(exact.raw_prediction.duration_hour, input.observed_duration_hour));
    assert(Near(exact.raw_prediction.trough_speed_mph,
                input.observed_trough_speed_mph));
    assert(Near(exact.raw_prediction.average_speed_mph,
                input.observed_average_speed_mph));

    config.use_bounds = true;
    config.maximum_plf = 0.20;
    config.maximum_qcd = 3.0;
    taplite::FixedVolumeOracleResult bounded =
        taplite::SolveFixedVolumeOracle(input, config);
    assert(bounded.status == taplite::OracleStatus::MultiBoundLimited);
    assert(bounded.bound_count == 2);
    assert(!Near(bounded.applied_prediction.duration_hour,
                 input.observed_duration_hour) ||
           !Near(bounded.applied_prediction.trough_speed_mph,
                 input.observed_trough_speed_mph) ||
           !Near(bounded.applied_prediction.average_speed_mph,
                 input.observed_average_speed_mph));

    taplite::FixedVolumeOracleInput incompatible = input;
    incompatible.observed_average_speed_mph = 5.0;
    taplite::FixedVolumeOracleResult failed =
        taplite::SolveFixedVolumeOracle(incompatible, config);
    assert(failed.status == taplite::OracleStatus::AverageSpeedIncompatible);

    incompatible = input;
    incompatible.observed_trough_speed_mph = 50.0;
    failed = taplite::SolveFixedVolumeOracle(incompatible, config);
    assert(failed.status == taplite::OracleStatus::QcpSignIncompatible);

    // The further-development copy releases alpha from theta*Qcp*Qcd^s.
    // This target is not tied to the refinement's queue factor, but all three
    // observables are exactly reconstructable with an independent alpha.
    taplite::FixedVolumeOracleInput independent = input;
    independent.observed_average_speed_mph = 35.0;
    config.use_bounds = false;
    config.minimum_plf = 0.01;
    config.maximum_plf = 10.0;
    taplite::FixedVolumeOracleResult further =
        taplite::SolveIndependentAlphaOracle(independent, config, 0.40);
    assert(further.status == taplite::OracleStatus::ExactFeasible);
    assert(further.exact_feasible);
    assert(Near(further.raw_prediction.duration_hour,
                independent.observed_duration_hour));
    assert(Near(further.raw_prediction.trough_speed_mph,
                independent.observed_trough_speed_mph));
    assert(Near(further.raw_prediction.average_speed_mph,
                independent.observed_average_speed_mph));

    // The further-development inverse must not encode a non-identifiable
    // average-speed target as an arbitrarily large routing-cost coefficient.
    config.use_bounds = true;
    config.minimum_plf = 0.001;
    config.maximum_plf = 100.0;
    config.minimum_qcd = 1e-9;
    config.maximum_qcd = 1e9;
    config.minimum_qcp = 1e-12;
    config.maximum_qcp = 1e9;
    config.minimum_alpha = 1e-6;
    config.maximum_alpha = 0.01;
    taplite::FixedVolumeOracleResult alpha_bounded =
        taplite::SolveIndependentAlphaOracle(independent, config, 0.40);
    assert(alpha_bounded.status == taplite::OracleStatus::AlphaBoundLimited);
    assert(alpha_bounded.bound_count == 1);
    assert(Near(alpha_bounded.alpha_applied, config.maximum_alpha));
    assert(!alpha_bounded.exact_feasible);

    std::cout << "auto-calibration refinement oracle tests passed\n";
    return 0;
}
