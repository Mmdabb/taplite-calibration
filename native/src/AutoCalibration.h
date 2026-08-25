#ifndef TAPLITE_AUTO_CALIBRATION_H
#define TAPLITE_AUTO_CALIBRATION_H

#include <map>
#include <string>
#include <vector>

namespace taplite {

// Analytical fixed-volume inverse used by the production refinement.  These
// declarations live with the engine so consumers have one calibration API and
// one native implementation to build.
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

enum class ObservationClass {
    Unobserved = 0,
    Episode = 1,
    NoEpisode = 2,
};

struct AutoCalibrationConfig {
    bool enabled;
    int max_outer_iterations;
    int minimum_outer_iterations;
    int maximum_retries;
    int workers;
    double theta;
    double period_start_hour;
    double period_end_hour;
    double departure_bin_minutes;
    double arrival_tail_minutes;
    double minimum_column_od_volume;
    int maximum_column_paths_per_od;
    double minimum_column_path_share;
    double plf_damping;
    double q_damping;
    double rejected_step_reduction;
    double departure_plf_weight;
    double mode1_plf_weight;
    double group_plf_weight;
    double minimum_plf;
    double maximum_plf;
    double minimum_qcd;
    double maximum_qcd;
    double minimum_qcp;
    double maximum_qcp;
    double objective_relative_tolerance;
    double parameter_relative_tolerance;
    double plf_relative_tolerance;
    double route_policy_tolerance;
    double unobserved_volume_tolerance;
    double system_vmt_tolerance;
    double system_vht_tolerance;
    double inner_gap_tolerance_pct;
    double speed_scale_mph;
    double no_episode_min_duration_hour;
    double huber_delta;
    double weight_duration;
    double weight_trough_speed;
    double weight_average_speed;
    double weight_s3;
    double weight_volume_envelope;
    double weight_count;
    double weight_no_episode;
    double weight_q_prior;
    double weight_plf_prior;
    std::string calibration_fit_mode;
    std::string candidate_direction_mode;
    std::string acceptance_policy;
    std::string class_n_role;
    std::string duration_target_policy;
    double maximum_duration_ratio;
    bool oracle_use_bounds;
    bool include_oracle_plf_in_target;
    bool write_oracle_audit;
    double oracle_residual_tolerance;
    double average_speed_feasibility_mph;
    double saturated_speed_tolerance_mph;
    double oracle_plf_weight;
    double max_e_fit_degradation;
    // Observed N links retain their accepted PLF as the base of the next
    // controller step. E links remain anchored by the analytical P/vt2/average
    // speed oracle, while N-link speed corrections can accumulate across outers.
    bool observed_plf_accumulation;
    double speed_control_gain;
    // Optional N-link override; a negative value reuses speed_control_gain.
    double non_episode_speed_control_gain;
    // Scales speed/S3/envelope/screenline PLF control after the analytical
    // oracle on E links. Zero gives the oracle exclusive control of E links.
    double episode_external_control_scale;
    // Holds accepted N-link parameters fixed during an episode-equilibrium
    // closure pass after the N-link speed polish.
    bool freeze_non_episode_parameters;
    double s3_control_gain;
    double volume_envelope_control_gain;
    // Asymmetric refinement control: under-assigned links below both the CUBE
    // and S3 volumes receive a stronger upward-volume PLF pressure. Links
    // above the envelope continue to use volume_envelope_control_gain.
    double volume_below_envelope_control_gain;
    double count_control_gain;
    double maximum_log_plf_control;
    std::string departure_profile_file;
    std::string volume_constraint_file;
    // Optional final-only export written from the accepted DTAC-v2 route pool.
    // Empty keeps the historical in-memory-only behavior.
    std::string arrival_profile_output_file;
    std::string history_output_file;
    std::string audit_output_file;
    std::string volume_constraint_audit_output_file;
    std::string summary_output_file;
    std::string oracle_audit_output_file;

    AutoCalibrationConfig();
    static AutoCalibrationConfig Load(
        const std::string& path,
        double period_start_hour,
        double period_end_hour);
};

struct CalibrationLink {
    int internal_index;
    int external_link_id;
    int from_node_id;
    int to_node_id;
    int vdf_code;
    std::string corridor;
    std::string facility_class;
    std::string target_tmc;
    std::string calibration_exclusion_reason;
    bool calibration_eligible;
    ObservationClass observation_class;

    double length_miles;
    double lanes;
    double capacity_vphpl;
    double free_speed_mph;
    double cutoff_speed_mph;
    double volume;
    double vehicle_volume;
    double travel_time_minutes;

    double plf;
    double qcd;
    double qcp;
    double qn;
    double qs;
    double alpha;
    double beta;

    double mode1_plf;
    double mode1_qcd;
    double mode1_qcp;
    double observed_duration_hour;
    double observed_trough_speed_mph;
    double observed_average_speed_mph;
    double s3_volume;
    double cube_vehicle_volume;
    double observation_quality;

    CalibrationLink();
};

struct CalibrationRoute {
    int mode_index;
    int origin_zone_id;
    int destination_zone_id;
    std::string mode_name;
    double od_demand;
    double pce;
    double share;
    std::vector<int> link_indices;
    std::vector<double> link_travel_times_minutes;
};

struct DepartureProfilePoint {
    double minute_of_day;
    double share;
};

// Compact regional handoff. The C++ bridge derives this directly from the
// existing DTAC-v2 pool instead of duplicating every path and link sequence.
struct CalibrationRouteSummary {
    std::vector<double> departure_plf;
    std::vector<double> policy_weights;
};

struct CalibrationEvaluation {
    double objective;
    double episode_fit_loss;
    double duration_loss;
    double trough_speed_loss;
    double average_speed_loss;
    double s3_loss;
    double volume_envelope_loss;
    double count_loss;
    double no_episode_loss;
    double prior_loss;
    double route_policy_distance;
    double unobserved_volume_deviation;
    double vmt_change_fraction;
    double vht_change_fraction;
    double relative_gap_pct;
    bool guardrails_pass;

    CalibrationEvaluation();
};

struct CalibrationProposal {
    std::vector<double> departure_plf;
    std::vector<double> plf;
    std::vector<double> qcd;
    std::vector<double> qcp;
    std::vector<double> alpha;
    std::vector<double> beta;
    std::vector<FixedVolumeOracleResult> oracle;
    double maximum_parameter_change;
    double maximum_plf_change;

    CalibrationProposal();
};

struct CalibrationIteration {
    int outer_iteration;
    int retry;
    bool accepted;
    double step_scale;
    CalibrationEvaluation evaluation;
    double maximum_parameter_change;
    double maximum_plf_change;
};

struct CalibrationVolumeConstraint {
    std::string constraint_id;
    std::string constraint_type;
    double target_vehicle_volume;
    double weight;
    std::vector<int> link_indices;
    std::vector<double> coefficients;

    CalibrationVolumeConstraint();
};

class AutoCalibrationEngine {
public:
    explicit AutoCalibrationEngine(const AutoCalibrationConfig& config);

    bool LoadDepartureProfiles(std::string* error_message);
    bool LoadVolumeConstraints(
        const std::vector<CalibrationLink>& links,
        std::string* error_message);
    void Initialize(
        const std::vector<CalibrationLink>& links,
        const std::vector<CalibrationRoute>& routes,
        double relative_gap_pct);
    void Initialize(
        const std::vector<CalibrationLink>& links,
        const CalibrationRouteSummary& route_summary,
        double relative_gap_pct);
    CalibrationProposal Propose(
        const std::vector<CalibrationLink>& links,
        const std::vector<CalibrationRoute>& routes,
        double step_scale) const;
    CalibrationProposal Propose(
        const std::vector<CalibrationLink>& links,
        const CalibrationRouteSummary& route_summary,
        double step_scale) const;
    CalibrationEvaluation Evaluate(
        const std::vector<CalibrationLink>& links,
        const std::vector<CalibrationRoute>& routes,
        double relative_gap_pct) const;
    CalibrationEvaluation Evaluate(
        const std::vector<CalibrationLink>& links,
        const CalibrationRouteSummary& route_summary,
        double relative_gap_pct) const;
    bool ShouldAccept(const CalibrationEvaluation& evaluation) const;
    bool Converged(
        const CalibrationEvaluation& evaluation,
        const CalibrationProposal& proposal,
        int accepted_outer_iterations) const;
    void Record(
        int outer_iteration,
        int retry,
        bool accepted,
        double step_scale,
        const CalibrationEvaluation& evaluation,
        const CalibrationProposal& proposal);
    void Accept(const CalibrationEvaluation& evaluation);
    bool WriteOutputs(
        const std::vector<CalibrationLink>& links,
        const CalibrationProposal& final_proposal,
        std::string* error_message) const;

    const AutoCalibrationConfig& config() const { return config_; }
    const CalibrationEvaluation& accepted_evaluation() const {
        return accepted_evaluation_;
    }
    const std::vector<CalibrationIteration>& history() const { return history_; }
    bool uses_origin_specific_profiles() const { return uses_origin_profiles_; }
    bool uses_mode_fallback_profiles() const { return uses_mode_profiles_; }
    const std::vector<DepartureProfilePoint>* DepartureProfileFor(
        int origin_zone_id,
        const std::string& mode_name) const;

private:
    typedef std::pair<int, std::string> ProfileKey;

    const std::vector<DepartureProfilePoint>* FindProfile(
        int origin_zone_id,
        const std::string& mode_name) const;
    std::vector<double> BuildDeparturePLF(
        const std::vector<CalibrationLink>& links,
        const std::vector<CalibrationRoute>& routes) const;
    CalibrationRouteSummary SummarizeRoutes(
        const std::vector<CalibrationLink>& links,
        const std::vector<CalibrationRoute>& routes) const;
    double RoutePolicyDistance(const CalibrationRouteSummary& summary) const;

    AutoCalibrationConfig config_;
    std::map<ProfileKey, std::vector<DepartureProfilePoint> > profiles_;
    std::vector<CalibrationVolumeConstraint> volume_constraints_;
    std::vector<CalibrationLink> baseline_links_;
    std::vector<double> baseline_policy_weights_;
    CalibrationEvaluation accepted_evaluation_;
    std::vector<CalibrationIteration> history_;
    bool initialized_;
    bool uses_origin_profiles_;
    bool uses_mode_profiles_;
};

const char* ObservationClassName(ObservationClass value);
ObservationClass ParseObservationClass(const std::string& value);

}  // namespace taplite

#endif  // TAPLITE_AUTO_CALIBRATION_H
