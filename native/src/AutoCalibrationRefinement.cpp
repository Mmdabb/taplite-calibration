#include "AutoCalibrationRefinement.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace taplite {
namespace {

const double kEpsilon = 1e-12;

double NaN() {
    return std::numeric_limits<double>::quiet_NaN();
}

bool Positive(double value) {
    return std::isfinite(value) && value > 0.0;
}

double Clamp(double value, double lower, double upper) {
    return std::max(lower, std::min(upper, value));
}

bool Changed(double raw, double applied) {
    return std::isfinite(raw) && std::isfinite(applied) &&
        std::fabs(raw - applied) > 1e-12 * std::max(1.0, std::fabs(raw));
}

double BoundRatio(double raw, double lower, double upper) {
    if (!Positive(raw)) {
        return NaN();
    }
    if (raw < lower) {
        return raw / std::max(lower, kEpsilon);
    }
    if (raw > upper) {
        return raw / std::max(upper, kEpsilon);
    }
    return 1.0;
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

QVDFPrediction::QVDFPrediction()
    : doc(NaN()),
      duration_hour(NaN()),
      trough_speed_mph(NaN()),
      reference_speed_mph(NaN()),
      queue_speed_mph(NaN()),
      average_speed_mph(NaN()) {}

FixedVolumeOracleInput::FixedVolumeOracleInput()
    : volume(NaN()),
      lanes(NaN()),
      period_hours(NaN()),
      capacity_vphpl(NaN()),
      free_speed_mph(NaN()),
      cutoff_speed_mph(NaN()),
      qn(NaN()),
      qs(NaN()),
      observed_duration_hour(NaN()),
      observed_trough_speed_mph(NaN()),
      observed_average_speed_mph(NaN()) {}

FixedVolumeOracleConfig::FixedVolumeOracleConfig()
    : theta(8.0 / 15.0),
      use_bounds(false),
      minimum_plf(0.10),
      maximum_plf(1.25),
      minimum_qcd(0.01),
      maximum_qcd(20.0),
      minimum_qcp(0.001),
      maximum_qcp(20.0),
      minimum_alpha(1e-6),
      maximum_alpha(10.0),
      residual_tolerance(1e-7),
      average_speed_feasibility_mph(0.10),
      saturated_speed_tolerance_mph(0.25) {}

FixedVolumeOracleResult::FixedVolumeOracleResult()
    : status(OracleStatus::NotApplicable),
      exact_feasible(false),
      duration_semantic_review(false),
      bound_count(0),
      target_reference_speed_mph(NaN()),
      doc_raw(NaN()),
      plf_raw(NaN()),
      qcd_raw(NaN()),
      qcp_raw(NaN()),
      alpha_raw(NaN()),
      beta_raw(NaN()),
      plf_bound_ratio(NaN()),
      qcd_bound_ratio(NaN()),
      qcp_bound_ratio(NaN()),
      alpha_bound_ratio(NaN()),
      plf_applied(NaN()),
      qcd_applied(NaN()),
      qcp_applied(NaN()),
      alpha_applied(NaN()),
      beta_applied(NaN()),
      raw_duration_residual(NaN()),
      raw_trough_residual(NaN()),
      raw_average_residual(NaN()),
      applied_duration_residual(NaN()),
      applied_trough_residual(NaN()),
      applied_average_residual(NaN()) {}

QVDFPrediction EvaluateRefinedQVDF(
    const FixedVolumeOracleInput& input,
    double plf,
    double qcd,
    double qcp,
    double alpha,
    double beta) {
    QVDFPrediction result;
    const double denominator = input.lanes * input.period_hours * plf *
        input.capacity_vphpl;
    if (!Positive(denominator) || !Positive(qcd) || !Positive(qcp) ||
        !Positive(input.qn) || !Positive(input.qs) || !Positive(alpha) ||
        !Positive(beta) || !Positive(input.free_speed_mph) ||
        !Positive(input.cutoff_speed_mph)) {
        return result;
    }
    result.doc = input.volume / denominator;
    if (!std::isfinite(result.doc) || result.doc < 0.0) {
        return QVDFPrediction();
    }
    result.duration_hour = qcd * std::pow(result.doc, input.qn);
    result.trough_speed_mph = input.cutoff_speed_mph /
        std::max(kEpsilon, 1.0 + qcp *
            std::pow(result.duration_hour, input.qs));
    result.reference_speed_mph = result.doc < 1.0
        ? (1.0 - result.doc) * input.free_speed_mph +
            result.doc * input.cutoff_speed_mph
        : input.cutoff_speed_mph;
    result.queue_speed_mph = result.reference_speed_mph /
        std::max(kEpsilon, 1.0 + alpha * std::pow(result.doc, beta));
    result.average_speed_mph = result.duration_hour > input.period_hours
        ? result.queue_speed_mph
        : (result.duration_hour / input.period_hours) * result.queue_speed_mph +
            (1.0 - result.duration_hour / input.period_hours) *
                (result.reference_speed_mph + input.free_speed_mph) / 2.0;
    return result;
}

FixedVolumeOracleResult SolveFixedVolumeOracle(
    const FixedVolumeOracleInput& input,
    const FixedVolumeOracleConfig& config) {
    FixedVolumeOracleResult result;
    if (!Positive(input.volume) || !Positive(input.lanes) ||
        !Positive(input.period_hours) || !Positive(input.capacity_vphpl) ||
        !Positive(input.free_speed_mph) || !Positive(input.cutoff_speed_mph) ||
        input.cutoff_speed_mph >= input.free_speed_mph || !Positive(input.qn) ||
        !Positive(input.qs) || !Positive(input.observed_duration_hour) ||
        !Positive(input.observed_trough_speed_mph) ||
        !Positive(input.observed_average_speed_mph) || !Positive(config.theta)) {
        result.status = OracleStatus::InvalidTarget;
        result.detail = "nonpositive/nonfinite input or cutoff not below free speed";
        return result;
    }
    result.duration_semantic_review = input.observed_duration_hour > 24.0;
    const double severity = input.cutoff_speed_mph /
        input.observed_trough_speed_mph - 1.0;
    if (!(severity > 0.0)) {
        result.status = OracleStatus::QcpSignIncompatible;
        result.detail = "observed trough is at or above cutoff speed";
        return result;
    }
    const double queue_factor = 1.0 + config.theta * severity;
    if (!Positive(queue_factor)) {
        result.status = OracleStatus::InvalidTarget;
        result.detail = "invalid analytical queue factor";
        return result;
    }

    const double p = input.observed_duration_hour;
    if (p <= input.period_hours) {
        const double a = p / input.period_hours;
        const double coefficient = a / queue_factor + (1.0 - a) / 2.0;
        const double free_coefficient = (1.0 - a) / 2.0;
        if (std::fabs(coefficient) <= kEpsilon) {
            result.status = OracleStatus::InvalidTarget;
            result.detail = "average-speed inverse has a zero coefficient";
            return result;
        }
        result.target_reference_speed_mph =
            (input.observed_average_speed_mph -
             free_coefficient * input.free_speed_mph) / coefficient;
    } else {
        result.target_reference_speed_mph =
            queue_factor * input.observed_average_speed_mph;
    }

    const double speed_tolerance = std::max(
        config.average_speed_feasibility_mph, 0.0);
    if (std::fabs(result.target_reference_speed_mph -
                  input.cutoff_speed_mph) <=
        std::max(config.saturated_speed_tolerance_mph, 0.0)) {
        result.status = OracleStatus::SaturatedNonidentifiable;
        result.detail = "average speed identifies the saturated branch but not DOC";
        return result;
    }
    if (result.target_reference_speed_mph <=
            input.cutoff_speed_mph + speed_tolerance ||
        result.target_reference_speed_mph >
            input.free_speed_mph + speed_tolerance) {
        result.status = OracleStatus::AverageSpeedIncompatible;
        result.detail = "required reference speed is outside (cutoff, free flow]";
        return result;
    }
    result.target_reference_speed_mph = std::min(
        result.target_reference_speed_mph, input.free_speed_mph);
    result.doc_raw = (input.free_speed_mph - result.target_reference_speed_mph) /
        (input.free_speed_mph - input.cutoff_speed_mph);
    if (!(result.doc_raw > 0.0 && result.doc_raw < 1.0)) {
        result.status = OracleStatus::AverageSpeedIncompatible;
        result.detail = "inverse reference speed does not produce DOC in (0,1)";
        return result;
    }

    result.plf_raw = input.volume /
        (input.lanes * input.period_hours * input.capacity_vphpl * result.doc_raw);
    result.qcd_raw = p / std::pow(result.doc_raw, input.qn);
    result.qcp_raw = severity / std::pow(p, input.qs);
    result.alpha_raw = config.theta * result.qcp_raw *
        std::pow(result.qcd_raw, input.qs);
    result.beta_raw = input.qn * input.qs;
    if (!Positive(result.plf_raw) || !Positive(result.qcd_raw) ||
        !Positive(result.qcp_raw) || !Positive(result.alpha_raw) ||
        !Positive(result.beta_raw)) {
        result.status = OracleStatus::InvalidTarget;
        result.detail = "inverse solution produced a nonpositive parameter";
        return result;
    }

    result.plf_bound_ratio = BoundRatio(
        result.plf_raw, config.minimum_plf, config.maximum_plf);
    result.qcd_bound_ratio = BoundRatio(
        result.qcd_raw, config.minimum_qcd, config.maximum_qcd);
    result.qcp_bound_ratio = BoundRatio(
        result.qcp_raw, config.minimum_qcp, config.maximum_qcp);
    result.plf_applied = config.use_bounds
        ? Clamp(result.plf_raw, config.minimum_plf, config.maximum_plf)
        : result.plf_raw;
    result.qcd_applied = config.use_bounds
        ? Clamp(result.qcd_raw, config.minimum_qcd, config.maximum_qcd)
        : result.qcd_raw;
    result.qcp_applied = config.use_bounds
        ? Clamp(result.qcp_raw, config.minimum_qcp, config.maximum_qcp)
        : result.qcp_raw;
    result.alpha_applied = config.theta * result.qcp_applied *
        std::pow(result.qcd_applied, input.qs);
    result.beta_applied = result.beta_raw;
    const bool plf_limited = Changed(result.plf_raw, result.plf_applied);
    const bool qcd_limited = Changed(result.qcd_raw, result.qcd_applied);
    const bool qcp_limited = Changed(result.qcp_raw, result.qcp_applied);
    result.bound_count = static_cast<int>(plf_limited) +
        static_cast<int>(qcd_limited) + static_cast<int>(qcp_limited);
    if (result.bound_count > 1) {
        result.status = OracleStatus::MultiBoundLimited;
        result.detail = "multiple inverse parameters are outside configured bounds";
    } else if (plf_limited) {
        result.status = OracleStatus::PlfBoundLimited;
        result.detail = "inverse PLF is outside configured bounds";
    } else if (qcd_limited) {
        result.status = OracleStatus::QcdBoundLimited;
        result.detail = "inverse Qcd is outside configured bounds";
    } else if (qcp_limited) {
        result.status = OracleStatus::QcpBoundLimited;
        result.detail = "inverse Qcp is outside configured bounds";
    } else {
        result.status = OracleStatus::ExactFeasible;
        result.detail = result.duration_semantic_review
            ? "exact algebraic fit; duration exceeds 24 hours and requires semantic review"
            : "exact algebraic fit";
    }

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
    const double tolerance = std::max(config.residual_tolerance, 0.0);
    result.exact_feasible =
        std::fabs(result.raw_duration_residual) <= tolerance &&
        std::fabs(result.raw_trough_residual) <= tolerance &&
        std::fabs(result.raw_average_residual) <=
            std::max(tolerance, config.average_speed_feasibility_mph);
    if (!result.exact_feasible && result.status == OracleStatus::ExactFeasible) {
        result.status = OracleStatus::SemanticReview;
        result.detail = "inverse parameters did not reconstruct targets within tolerance";
    }
    return result;
}

const char* OracleStatusName(OracleStatus status) {
    switch (status) {
        case OracleStatus::ExactFeasible: return "EXACT_FEASIBLE";
        case OracleStatus::InvalidTarget: return "INVALID_TARGET";
        case OracleStatus::QcpSignIncompatible: return "QCP_SIGN_INCOMPATIBLE";
        case OracleStatus::AverageSpeedIncompatible: return "AVG_SPEED_INCOMPATIBLE";
        case OracleStatus::SaturatedNonidentifiable: return "SATURATED_NONIDENTIFIABLE";
        case OracleStatus::PlfBoundLimited: return "PLF_BOUND_LIMITED";
        case OracleStatus::QcdBoundLimited: return "QCD_BOUND_LIMITED";
        case OracleStatus::QcpBoundLimited: return "QCP_BOUND_LIMITED";
        case OracleStatus::AlphaBoundLimited: return "ALPHA_BOUND_LIMITED";
        case OracleStatus::MultiBoundLimited: return "MULTI_BOUND_LIMITED";
        case OracleStatus::SemanticReview: return "SEMANTIC_REVIEW";
        default: return "NOT_APPLICABLE";
    }
}

}  // namespace taplite
