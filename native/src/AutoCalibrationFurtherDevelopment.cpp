#include "AutoCalibrationFurtherDevelopment.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace taplite {
namespace {

const double kEpsilon = 1e-12;

bool Positive(double value) {
    return std::isfinite(value) && value > 0.0;
}

double NaN() {
    return std::numeric_limits<double>::quiet_NaN();
}

double Clamp(double value, double lower, double upper) {
    return std::max(lower, std::min(upper, value));
}

double ReferenceSpeed(const FixedVolumeOracleInput& input, double doc) {
    return doc < 1.0
        ? (1.0 - doc) * input.free_speed_mph + doc * input.cutoff_speed_mph
        : input.cutoff_speed_mph;
}

struct Candidate {
    bool valid;
    bool exact_average;
    double doc;
    double alpha;
    double score;
    double average_residual;

    Candidate()
        : valid(false), exact_average(false), doc(NaN()), alpha(NaN()),
          score(std::numeric_limits<double>::infinity()),
          average_residual(NaN()) {}
};

Candidate BuildCandidate(
    const FixedVolumeOracleInput& input,
    double doc,
    double preferred_doc,
    double beta) {
    Candidate result;
    if (!Positive(doc) || !Positive(beta)) {
        return result;
    }
    const double reference = ReferenceSpeed(input, doc);
    const double congested_fraction = Clamp(
        input.observed_duration_hour / input.period_hours, 0.0, 1.0);
    if (!(congested_fraction > 0.0)) {
        return result;
    }
    const double uncongested_component = (1.0 - congested_fraction) *
        (reference + input.free_speed_mph) / 2.0;
    const double desired_queue_speed =
        (input.observed_average_speed_mph - uncongested_component) /
        congested_fraction;
    const double queue_speed = Clamp(desired_queue_speed, kEpsilon, reference);
    const double alpha = std::max(
        kEpsilon,
        (reference / std::max(queue_speed, kEpsilon) - 1.0) /
            std::pow(doc, beta));
    const double modeled_average = congested_fraction *
        (reference / (1.0 + alpha * std::pow(doc, beta))) +
        uncongested_component;
    result.valid = Positive(alpha) && std::isfinite(modeled_average);
    result.exact_average = desired_queue_speed > 0.0 &&
        desired_queue_speed <= reference + 1e-9;
    result.doc = doc;
    result.alpha = alpha;
    result.average_residual = modeled_average - input.observed_average_speed_mph;
    const double plf_distance = std::fabs(std::log(
        doc / std::max(preferred_doc, kEpsilon)));
    result.score = result.exact_average
        ? plf_distance
        : 1e6 + std::fabs(result.average_residual) + 1e-6 * plf_distance;
    return result;
}

void SetResiduals(
    const FixedVolumeOracleInput& input,
    const QVDFPrediction& prediction,
    double* duration,
    double* trough,
    double* average) {
    *duration = prediction.duration_hour - input.observed_duration_hour;
    *trough = prediction.trough_speed_mph - input.observed_trough_speed_mph;
    *average = prediction.average_speed_mph - input.observed_average_speed_mph;
}

}  // namespace

FixedVolumeOracleResult SolveIndependentAlphaOracle(
    const FixedVolumeOracleInput& input,
    const FixedVolumeOracleConfig& config,
    double preferred_plf) {
    FixedVolumeOracleResult result;
    if (!Positive(input.volume) || !Positive(input.lanes) ||
        !Positive(input.period_hours) || !Positive(input.capacity_vphpl) ||
        !Positive(input.free_speed_mph) || !Positive(input.cutoff_speed_mph) ||
        input.cutoff_speed_mph >= input.free_speed_mph || !Positive(input.qn) ||
        !Positive(input.qs) || !Positive(input.observed_duration_hour) ||
        !Positive(input.observed_trough_speed_mph) ||
        !Positive(input.observed_average_speed_mph) || !Positive(preferred_plf)) {
        result.status = OracleStatus::InvalidTarget;
        result.detail = "invalid independent-alpha target or preferred PLF";
        return result;
    }
    const double severity = input.cutoff_speed_mph /
        input.observed_trough_speed_mph - 1.0;
    if (!(severity > 0.0)) {
        result.status = OracleStatus::QcpSignIncompatible;
        result.detail = "observed trough is at or above cutoff speed";
        return result;
    }

    const double volume_factor = input.volume /
        (input.lanes * input.period_hours * input.capacity_vphpl);
    const double preferred_doc = volume_factor / preferred_plf;
    const double lower_plf = std::max(config.minimum_plf, 1e-4);
    const double upper_plf = std::max(lower_plf, config.maximum_plf);
    double lower_doc = std::max(1e-5, volume_factor / upper_plf);
    double upper_doc = std::min(100.0, volume_factor / lower_plf);
    if (!(upper_doc >= lower_doc)) {
        result.status = OracleStatus::InvalidTarget;
        result.detail = "PLF range produces an empty DOC search interval";
        return result;
    }
    const double beta = input.qn * input.qs;
    std::vector<double> docs;
    docs.reserve(405);
    docs.push_back(Clamp(preferred_doc, lower_doc, upper_doc));
    docs.push_back(Clamp(1.0, lower_doc, upper_doc));
    const double log_lower = std::log(lower_doc);
    const double log_upper = std::log(upper_doc);
    for (int i = 0; i <= 400; ++i) {
        docs.push_back(std::exp(
            log_lower + (log_upper - log_lower) * i / 400.0));
    }

    Candidate best;
    for (std::size_t i = 0; i < docs.size(); ++i) {
        const Candidate candidate = BuildCandidate(
            input, docs[i], preferred_doc, beta);
        if (candidate.valid && candidate.score < best.score) {
            best = candidate;
        }
    }
    if (!best.valid) {
        result.status = OracleStatus::AverageSpeedIncompatible;
        result.detail = "no valid independent-alpha DOC candidate";
        return result;
    }

    result.doc_raw = best.doc;
    result.plf_raw = volume_factor / best.doc;
    result.qcd_raw = input.observed_duration_hour /
        std::pow(best.doc, input.qn);
    result.qcp_raw = severity /
        std::pow(input.observed_duration_hour, input.qs);
    result.alpha_raw = best.alpha;
    result.beta_raw = beta;
    result.target_reference_speed_mph = ReferenceSpeed(input, best.doc);
    result.plf_applied = result.plf_raw;
    result.qcd_applied = config.use_bounds
        ? Clamp(result.qcd_raw, config.minimum_qcd, config.maximum_qcd)
        : result.qcd_raw;
    result.qcp_applied = config.use_bounds
        ? Clamp(result.qcp_raw, config.minimum_qcp, config.maximum_qcp)
        : result.qcp_raw;
    result.alpha_applied = config.use_bounds
        ? Clamp(result.alpha_raw, config.minimum_alpha, config.maximum_alpha)
        : result.alpha_raw;
    result.beta_applied = result.beta_raw;
    result.plf_bound_ratio = 1.0;
    result.qcd_bound_ratio = result.qcd_applied == result.qcd_raw ? 1.0 :
        result.qcd_raw / std::max(result.qcd_applied, kEpsilon);
    result.qcp_bound_ratio = result.qcp_applied == result.qcp_raw ? 1.0 :
        result.qcp_raw / std::max(result.qcp_applied, kEpsilon);
    result.alpha_bound_ratio = result.alpha_applied == result.alpha_raw ? 1.0 :
        result.alpha_raw / std::max(result.alpha_applied, kEpsilon);
    result.bound_count = static_cast<int>(result.qcd_bound_ratio != 1.0) +
        static_cast<int>(result.qcp_bound_ratio != 1.0) +
        static_cast<int>(result.alpha_bound_ratio != 1.0);
    result.raw_prediction = EvaluateRefinedQVDF(
        input, result.plf_raw, result.qcd_raw, result.qcp_raw,
        result.alpha_raw, result.beta_raw);
    result.applied_prediction = EvaluateRefinedQVDF(
        input, result.plf_applied, result.qcd_applied, result.qcp_applied,
        result.alpha_applied, result.beta_applied);
    SetResiduals(input, result.raw_prediction,
        &result.raw_duration_residual, &result.raw_trough_residual,
        &result.raw_average_residual);
    SetResiduals(input, result.applied_prediction,
        &result.applied_duration_residual, &result.applied_trough_residual,
        &result.applied_average_residual);
    const double speed_tolerance = std::max(
        config.average_speed_feasibility_mph, config.residual_tolerance);
    result.exact_feasible = best.exact_average && result.bound_count == 0 &&
        std::fabs(result.raw_duration_residual) <= config.residual_tolerance &&
        std::fabs(result.raw_trough_residual) <= config.residual_tolerance &&
        std::fabs(result.raw_average_residual) <= speed_tolerance;
    if (result.exact_feasible) {
        result.status = OracleStatus::ExactFeasible;
        result.detail = "independent alpha exactly reconstructs P, vt2, and vbar";
    } else if (result.bound_count > 1) {
        result.status = OracleStatus::MultiBoundLimited;
        result.detail = "independent-alpha inverse is limited by multiple parameter bounds";
    } else if (result.bound_count == 1 && result.qcd_bound_ratio != 1.0) {
        result.status = OracleStatus::QcdBoundLimited;
        result.detail = "independent-alpha inverse is limited by Qcd bounds";
    } else if (result.bound_count == 1 && result.qcp_bound_ratio != 1.0) {
        result.status = OracleStatus::QcpBoundLimited;
        result.detail = "independent-alpha inverse is limited by Qcp bounds";
    } else if (result.bound_count == 1) {
        result.status = OracleStatus::AlphaBoundLimited;
        result.detail = "independent-alpha inverse is limited by alpha bounds";
    } else {
        result.status = OracleStatus::AverageSpeedIncompatible;
        result.detail = "average speed lies outside the QVDF range at bounded PLF";
    }
    return result;
}

}  // namespace taplite
