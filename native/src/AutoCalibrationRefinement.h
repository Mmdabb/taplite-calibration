#ifndef TAPLITE_AUTO_CALIBRATION_REFINEMENT_H
#define TAPLITE_AUTO_CALIBRATION_REFINEMENT_H

#include <string>

namespace taplite {

enum class OracleStatus {
    NotApplicable = 0,
    ExactFeasible,
    InvalidTarget,
    QcpSignIncompatible,
    AverageSpeedIncompatible,
    SaturatedNonidentifiable,
    PlfBoundLimited,
    QcdBoundLimited,
    QcpBoundLimited,
    AlphaBoundLimited,
    MultiBoundLimited,
    SemanticReview,
};

struct QVDFPrediction {
    double doc;
    double duration_hour;
    double trough_speed_mph;
    double reference_speed_mph;
    double queue_speed_mph;
    double average_speed_mph;

    QVDFPrediction();
};

struct FixedVolumeOracleInput {
    double volume;
    double lanes;
    double period_hours;
    double capacity_vphpl;
    double free_speed_mph;
    double cutoff_speed_mph;
    double qn;
    double qs;
    double observed_duration_hour;
    double observed_trough_speed_mph;
    double observed_average_speed_mph;

    FixedVolumeOracleInput();
};

struct FixedVolumeOracleConfig {
    double theta;
    bool use_bounds;
    double minimum_plf;
    double maximum_plf;
    double minimum_qcd;
    double maximum_qcd;
    double minimum_qcp;
    double maximum_qcp;
    double minimum_alpha;
    double maximum_alpha;
    double residual_tolerance;
    double average_speed_feasibility_mph;
    double saturated_speed_tolerance_mph;

    FixedVolumeOracleConfig();
};

struct FixedVolumeOracleResult {
    OracleStatus status;
    std::string detail;
    bool exact_feasible;
    bool duration_semantic_review;
    int bound_count;

    double target_reference_speed_mph;
    double doc_raw;
    double plf_raw;
    double qcd_raw;
    double qcp_raw;
    double alpha_raw;
    double beta_raw;
    double plf_bound_ratio;
    double qcd_bound_ratio;
    double qcp_bound_ratio;
    double alpha_bound_ratio;

    double plf_applied;
    double qcd_applied;
    double qcp_applied;
    double alpha_applied;
    double beta_applied;

    QVDFPrediction raw_prediction;
    QVDFPrediction applied_prediction;
    double raw_duration_residual;
    double raw_trough_residual;
    double raw_average_residual;
    double applied_duration_residual;
    double applied_trough_residual;
    double applied_average_residual;

    FixedVolumeOracleResult();
};

QVDFPrediction EvaluateRefinedQVDF(
    const FixedVolumeOracleInput& input,
    double plf,
    double qcd,
    double qcp,
    double alpha,
    double beta);

FixedVolumeOracleResult SolveFixedVolumeOracle(
    const FixedVolumeOracleInput& input,
    const FixedVolumeOracleConfig& config);

const char* OracleStatusName(OracleStatus status);

}  // namespace taplite

#endif  // TAPLITE_AUTO_CALIBRATION_REFINEMENT_H
