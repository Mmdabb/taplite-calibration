#include "AutoCalibration.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace taplite {
namespace {

const double kEpsilon = 1e-9;

double Clamp(double value, double minimum, double maximum) {
    return std::max(minimum, std::min(maximum, value));
}

double SafeLog(double value) {
    return std::log(std::max(value, kEpsilon));
}

bool IsFinitePositive(double value) {
    return std::isfinite(value) && value > 0.0;
}

std::string Trim(const std::string& value) {
    std::string::size_type first = 0;
    while (first < value.size() &&
           std::isspace(static_cast<unsigned char>(value[first]))) {
        ++first;
    }
    std::string::size_type last = value.size();
    while (last > first &&
           std::isspace(static_cast<unsigned char>(value[last - 1]))) {
        --last;
    }
    return value.substr(first, last - first);
}

std::string Lower(const std::string& value) {
    std::string result = Trim(value);
    std::transform(result.begin(), result.end(), result.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return result;
}

std::vector<std::string> SplitCSV(const std::string& line) {
    std::vector<std::string> fields;
    std::string field;
    bool quoted = false;
    for (std::string::size_type i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (c == '"') {
            if (quoted && i + 1 < line.size() && line[i + 1] == '"') {
                field.push_back('"');
                ++i;
            } else {
                quoted = !quoted;
            }
        } else if (c == ',' && !quoted) {
            fields.push_back(Trim(field));
            field.clear();
        } else {
            field.push_back(c);
        }
    }
    fields.push_back(Trim(field));
    return fields;
}

bool ParseDouble(const std::string& text, double* value) {
    if (value == NULL) {
        return false;
    }
    const std::string trimmed = Trim(text);
    if (trimmed.empty()) {
        return false;
    }
    char* end = NULL;
    const double parsed = std::strtod(trimmed.c_str(), &end);
    if (end == trimmed.c_str() || *end != '\0' || !std::isfinite(parsed)) {
        return false;
    }
    *value = parsed;
    return true;
}

bool ParseInt(const std::string& text, int* value) {
    double parsed = 0.0;
    if (!ParseDouble(text, &parsed) || std::floor(parsed) != parsed) {
        return false;
    }
    *value = static_cast<int>(parsed);
    return true;
}

bool ParseBool(const std::string& text, bool fallback) {
    const std::string value = Lower(text);
    if (value == "1" || value == "true" || value == "yes" || value == "on") {
        return true;
    }
    if (value == "0" || value == "false" || value == "no" || value == "off") {
        return false;
    }
    return fallback;
}

double ParseTimeMinute(const std::string& text) {
    const std::string value = Trim(text);
    const std::string::size_type colon = value.find(':');
    if (colon == std::string::npos) {
        double numeric = 0.0;
        return ParseDouble(value, &numeric) ? numeric :
            std::numeric_limits<double>::quiet_NaN();
    }
    double hour = 0.0;
    double minute = 0.0;
    if (!ParseDouble(value.substr(0, colon), &hour) ||
        !ParseDouble(value.substr(colon + 1), &minute)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return hour * 60.0 + minute;
}

double PeriodDurationHours(double start_hour, double end_hour) {
    double duration = end_hour - start_hour;
    if (duration <= 0.0) {
        duration += 24.0;
    }
    return std::max(duration, 1.0 / 60.0);
}

double UnwrapMinute(double minute, double start_minute, double end_minute) {
    double result = minute;
    if (end_minute <= start_minute && result < start_minute) {
        result += 1440.0;
    }
    return result;
}

double Huber(double value, double delta) {
    const double absolute = std::fabs(value);
    if (absolute <= delta) {
        return 0.5 * value * value;
    }
    return delta * (absolute - 0.5 * delta);
}

double Median(std::vector<double> values, double fallback) {
    values.erase(
        std::remove_if(values.begin(), values.end(),
                       [](double value) { return !IsFinitePositive(value); }),
        values.end());
    if (values.empty()) {
        return fallback;
    }
    const std::size_t middle = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + middle, values.end());
    double result = values[middle];
    if (values.size() % 2 == 0) {
        std::nth_element(values.begin(), values.begin() + middle - 1,
                         values.begin() + middle);
        result = 0.5 * (result + values[middle - 1]);
    }
    return result;
}

std::string CSVText(const std::string& value) {
    std::string escaped = value;
    std::string::size_type position = 0;
    while ((position = escaped.find('"', position)) != std::string::npos) {
        escaped.insert(position, 1, '"');
        position += 2;
    }
    return std::string("\"") + escaped + "\"";
}

void ApplyConfigValue(AutoCalibrationConfig* config,
                      const std::string& raw_key,
                      const std::string& raw_value) {
    const std::string key = Lower(raw_key);
    const std::string value = Trim(raw_value);
    double number = 0.0;
    int integer = 0;
    if (key == "enabled" || key == "auto_calibration") {
        config->enabled = ParseBool(value, config->enabled);
    } else if (key == "max_outer_iterations" && ParseInt(value, &integer)) {
        config->max_outer_iterations = std::max(0, integer);
    } else if (key == "minimum_outer_iterations" && ParseInt(value, &integer)) {
        config->minimum_outer_iterations = std::max(0, integer);
    } else if (key == "maximum_retries" && ParseInt(value, &integer)) {
        config->maximum_retries = std::max(1, integer);
    } else if (key == "workers" && ParseInt(value, &integer)) {
        config->workers = std::max(1, std::min(20, integer));
    } else if (key == "theta" && ParseDouble(value, &number)) {
        config->theta = number;
    } else if (key == "departure_bin_minutes" && ParseDouble(value, &number)) {
        config->departure_bin_minutes = number;
    } else if (key == "arrival_tail_minutes" && ParseDouble(value, &number)) {
        config->arrival_tail_minutes = number;
    } else if (key == "minimum_column_od_volume" && ParseDouble(value, &number)) {
        config->minimum_column_od_volume = number;
    } else if (key == "maximum_column_paths_per_od" && ParseInt(value, &integer)) {
        config->maximum_column_paths_per_od = integer;
    } else if (key == "minimum_column_path_share" && ParseDouble(value, &number)) {
        config->minimum_column_path_share = number;
    } else if (key == "plf_damping" && ParseDouble(value, &number)) {
        config->plf_damping = number;
    } else if (key == "q_damping" && ParseDouble(value, &number)) {
        config->q_damping = number;
    } else if (key == "rejected_step_reduction" && ParseDouble(value, &number)) {
        config->rejected_step_reduction = number;
    } else if (key == "departure_plf_weight" && ParseDouble(value, &number)) {
        config->departure_plf_weight = number;
    } else if (key == "mode1_plf_weight" && ParseDouble(value, &number)) {
        config->mode1_plf_weight = number;
    } else if (key == "group_plf_weight" && ParseDouble(value, &number)) {
        config->group_plf_weight = number;
    } else if (key == "minimum_plf" && ParseDouble(value, &number)) {
        config->minimum_plf = number;
    } else if (key == "maximum_plf" && ParseDouble(value, &number)) {
        config->maximum_plf = number;
    } else if (key == "minimum_qcd" && ParseDouble(value, &number)) {
        config->minimum_qcd = number;
    } else if (key == "maximum_qcd" && ParseDouble(value, &number)) {
        config->maximum_qcd = number;
    } else if (key == "minimum_qcp" && ParseDouble(value, &number)) {
        config->minimum_qcp = number;
    } else if (key == "maximum_qcp" && ParseDouble(value, &number)) {
        config->maximum_qcp = number;
    } else if (key == "objective_relative_tolerance" && ParseDouble(value, &number)) {
        config->objective_relative_tolerance = number;
    } else if (key == "parameter_relative_tolerance" && ParseDouble(value, &number)) {
        config->parameter_relative_tolerance = number;
    } else if (key == "plf_relative_tolerance" && ParseDouble(value, &number)) {
        config->plf_relative_tolerance = number;
    } else if (key == "route_policy_tolerance" && ParseDouble(value, &number)) {
        config->route_policy_tolerance = number;
    } else if (key == "unobserved_volume_tolerance" && ParseDouble(value, &number)) {
        config->unobserved_volume_tolerance = number;
    } else if (key == "system_vmt_tolerance" && ParseDouble(value, &number)) {
        config->system_vmt_tolerance = number;
    } else if (key == "system_vht_tolerance" && ParseDouble(value, &number)) {
        config->system_vht_tolerance = number;
    } else if (key == "inner_gap_tolerance_pct" && ParseDouble(value, &number)) {
        config->inner_gap_tolerance_pct = number;
    } else if (key == "speed_scale_mph" && ParseDouble(value, &number)) {
        config->speed_scale_mph = number;
    } else if (key == "no_episode_min_duration_hour" && ParseDouble(value, &number)) {
        config->no_episode_min_duration_hour = number;
    } else if (key == "huber_delta" && ParseDouble(value, &number)) {
        config->huber_delta = number;
    } else if (key == "weight_duration" && ParseDouble(value, &number)) {
        config->weight_duration = number;
    } else if (key == "weight_trough_speed" && ParseDouble(value, &number)) {
        config->weight_trough_speed = number;
    } else if (key == "weight_average_speed" && ParseDouble(value, &number)) {
        config->weight_average_speed = number;
    } else if (key == "weight_s3" && ParseDouble(value, &number)) {
        config->weight_s3 = number;
    } else if (key == "weight_volume_envelope" && ParseDouble(value, &number)) {
        config->weight_volume_envelope = number;
    } else if (key == "weight_count" && ParseDouble(value, &number)) {
        config->weight_count = number;
    } else if (key == "weight_no_episode" && ParseDouble(value, &number)) {
        config->weight_no_episode = number;
    } else if (key == "weight_q_prior" && ParseDouble(value, &number)) {
        config->weight_q_prior = number;
    } else if (key == "weight_plf_prior" && ParseDouble(value, &number)) {
        config->weight_plf_prior = number;
    } else if (key == "calibration_fit_mode") {
        config->calibration_fit_mode = Lower(value);
    } else if (key == "candidate_direction_mode") {
        config->candidate_direction_mode = Lower(value);
    } else if (key == "acceptance_policy") {
        config->acceptance_policy = Lower(value);
    } else if (key == "class_n_role") {
        config->class_n_role = Lower(value);
    } else if (key == "duration_target_policy") {
        config->duration_target_policy = Lower(value);
    } else if (key == "maximum_duration_ratio" && ParseDouble(value, &number)) {
        config->maximum_duration_ratio = number;
    } else if (key == "oracle_use_bounds") {
        config->oracle_use_bounds = ParseBool(value, config->oracle_use_bounds);
    } else if (key == "include_oracle_plf_in_target") {
        config->include_oracle_plf_in_target = ParseBool(
            value, config->include_oracle_plf_in_target);
    } else if (key == "write_oracle_audit") {
        config->write_oracle_audit = ParseBool(value, config->write_oracle_audit);
    } else if (key == "oracle_residual_tolerance" && ParseDouble(value, &number)) {
        config->oracle_residual_tolerance = number;
    } else if (key == "average_speed_feasibility_mph" && ParseDouble(value, &number)) {
        config->average_speed_feasibility_mph = number;
    } else if (key == "saturated_speed_tolerance_mph" && ParseDouble(value, &number)) {
        config->saturated_speed_tolerance_mph = number;
    } else if (key == "oracle_plf_weight" && ParseDouble(value, &number)) {
        config->oracle_plf_weight = number;
    } else if (key == "max_e_fit_degradation" && ParseDouble(value, &number)) {
        config->max_e_fit_degradation = number;
    } else if (key == "observed_plf_accumulation") {
        config->observed_plf_accumulation = ParseBool(
            value, config->observed_plf_accumulation);
    } else if (key == "speed_control_gain" && ParseDouble(value, &number)) {
        config->speed_control_gain = number;
    } else if (key == "non_episode_speed_control_gain" &&
               ParseDouble(value, &number)) {
        config->non_episode_speed_control_gain = number;
    } else if (key == "episode_external_control_scale" &&
               ParseDouble(value, &number)) {
        config->episode_external_control_scale = number;
    } else if (key == "freeze_non_episode_parameters") {
        config->freeze_non_episode_parameters = ParseBool(
            value, config->freeze_non_episode_parameters);
    } else if (key == "s3_control_gain" && ParseDouble(value, &number)) {
        config->s3_control_gain = number;
    } else if (key == "volume_envelope_control_gain" && ParseDouble(value, &number)) {
        config->volume_envelope_control_gain = number;
    } else if (key == "volume_below_envelope_control_gain" &&
               ParseDouble(value, &number)) {
        config->volume_below_envelope_control_gain = number;
    } else if (key == "count_control_gain" && ParseDouble(value, &number)) {
        config->count_control_gain = number;
    } else if (key == "maximum_log_plf_control" && ParseDouble(value, &number)) {
        config->maximum_log_plf_control = number;
    } else if (key == "departure_profile_file") {
        config->departure_profile_file = value;
    } else if (key == "volume_constraint_file") {
        config->volume_constraint_file = value;
    } else if (key == "arrival_profile_output_file") {
        config->arrival_profile_output_file = value;
    } else if (key == "history_output_file") {
        config->history_output_file = value;
    } else if (key == "audit_output_file") {
        config->audit_output_file = value;
    } else if (key == "volume_constraint_audit_output_file") {
        config->volume_constraint_audit_output_file = value;
    } else if (key == "summary_output_file") {
        config->summary_output_file = value;
    } else if (key == "oracle_audit_output_file") {
        config->oracle_audit_output_file = value;
    }
}

FixedVolumeOracleConfig BuildOracleConfig(
    const AutoCalibrationConfig& config) {
    FixedVolumeOracleConfig oracle;
    oracle.theta = config.theta;
    oracle.use_bounds = config.oracle_use_bounds;
    oracle.minimum_plf = config.minimum_plf;
    oracle.maximum_plf = config.maximum_plf;
    oracle.minimum_qcd = config.minimum_qcd;
    oracle.maximum_qcd = config.maximum_qcd;
    oracle.minimum_qcp = config.minimum_qcp;
    oracle.maximum_qcp = config.maximum_qcp;
    oracle.residual_tolerance = config.oracle_residual_tolerance;
    oracle.average_speed_feasibility_mph =
        config.average_speed_feasibility_mph;
    oracle.saturated_speed_tolerance_mph =
        config.saturated_speed_tolerance_mph;
    return oracle;
}

double EffectiveDurationTarget(
    const CalibrationLink& link,
    double period_hours,
    const AutoCalibrationConfig& config) {
    const double observed = link.observed_duration_hour;
    if (!IsFinitePositive(observed) ||
        config.calibration_fit_mode == "diagnostic_exact" ||
        config.duration_target_policy != "censor_to_period") {
        return observed;
    }
    const double ceiling = period_hours *
        std::max(config.maximum_duration_ratio, kEpsilon);
    return std::min(observed, ceiling);
}

FixedVolumeOracleInput BuildOracleInput(
    const CalibrationLink& link,
    double period_hours,
    const AutoCalibrationConfig& config) {
    FixedVolumeOracleInput input;
    input.volume = link.volume;
    input.lanes = link.lanes;
    input.period_hours = period_hours;
    input.capacity_vphpl = link.capacity_vphpl;
    input.free_speed_mph = link.free_speed_mph;
    input.cutoff_speed_mph = link.cutoff_speed_mph;
    input.qn = link.qn;
    input.qs = link.qs;
    input.observed_duration_hour = EffectiveDurationTarget(
        link, period_hours, config);
    input.observed_trough_speed_mph = link.observed_trough_speed_mph;
    input.observed_average_speed_mph = link.observed_average_speed_mph;
    return input;
}

bool IsRefinedMode(const AutoCalibrationConfig& config) {
    return config.calibration_fit_mode == "diagnostic_exact" ||
        config.calibration_fit_mode == "equilibrium_regularized" ||
        config.calibration_fit_mode == "refined_fixed_point";
}

double ModeledAverageSpeed(const CalibrationLink& link) {
    return link.travel_time_minutes > kEpsilon && link.length_miles > 0.0
        ? link.length_miles / (link.travel_time_minutes / 60.0)
        : link.free_speed_mph;
}

// Returns log(current/target).  A positive residual means the link carries too
// much traffic and therefore needs a downward PLF pressure (higher DOC/cost).
// When CUBE and S3 agree on direction, the nearest edge of their envelope is
// the strong target.  Inside the envelope only the weaker S3 controller stays
// active, as requested by the refinement design.
enum class VolumeEnvelopePosition {
    InsideOrSingle = 0,
    Below = 1,
    Above = 2,
};

const char* VolumeEnvelopePositionName(VolumeEnvelopePosition position) {
    switch (position) {
    case VolumeEnvelopePosition::Below:
        return "below_both";
    case VolumeEnvelopePosition::Above:
        return "above_both";
    default:
        return "inside_or_single";
    }
}

double VolumeDirectionResidual(
    const CalibrationLink& link,
    VolumeEnvelopePosition* position) {
    if (position != NULL) {
        *position = VolumeEnvelopePosition::InsideOrSingle;
    }
    const bool has_s3 = IsFinitePositive(link.s3_volume);
    const bool has_cube = IsFinitePositive(link.cube_vehicle_volume);
    const double current = std::max(0.0, link.vehicle_volume);
    double target = std::numeric_limits<double>::quiet_NaN();
    if (has_s3 && has_cube) {
        const double lower = std::min(link.s3_volume, link.cube_vehicle_volume);
        const double upper = std::max(link.s3_volume, link.cube_vehicle_volume);
        if (current < lower) {
            target = lower;
            if (position != NULL) {
                *position = VolumeEnvelopePosition::Below;
            }
        } else if (current > upper) {
            target = upper;
            if (position != NULL) {
                *position = VolumeEnvelopePosition::Above;
            }
        } else {
            target = link.s3_volume;
        }
    } else if (has_s3) {
        target = link.s3_volume;
    } else if (has_cube) {
        target = link.cube_vehicle_volume;
    }
    return IsFinitePositive(target)
        ? SafeLog((current + 1.0) / (target + 1.0)) : 0.0;
}

}  // namespace

AutoCalibrationConfig::AutoCalibrationConfig()
    : enabled(false),
      max_outer_iterations(6),
      minimum_outer_iterations(2),
      maximum_retries(6),
      workers(20),
      theta(8.0 / 15.0),
      period_start_hour(7.0),
      period_end_hour(8.0),
      departure_bin_minutes(15.0),
      arrival_tail_minutes(720.0),
      minimum_column_od_volume(0.05),
      maximum_column_paths_per_od(2),
      minimum_column_path_share(0.01),
      plf_damping(0.35),
      q_damping(0.30),
      rejected_step_reduction(0.5),
      departure_plf_weight(4.0),
      mode1_plf_weight(1.0),
      group_plf_weight(1.0),
      minimum_plf(0.10),
      maximum_plf(1.25),
      minimum_qcd(0.01),
      maximum_qcd(20.0),
      minimum_qcp(0.001),
      maximum_qcp(20.0),
      objective_relative_tolerance(0.002),
      parameter_relative_tolerance(0.005),
      plf_relative_tolerance(0.005),
      route_policy_tolerance(0.35),
      unobserved_volume_tolerance(0.35),
      system_vmt_tolerance(0.20),
      system_vht_tolerance(0.35),
      inner_gap_tolerance_pct(0.01),
      speed_scale_mph(5.0),
      no_episode_min_duration_hour(0.25),
      huber_delta(1.5),
      weight_duration(1.0),
      weight_trough_speed(1.0),
      weight_average_speed(0.75),
      weight_s3(0.15),
      weight_volume_envelope(0.75),
      weight_count(0.35),
      weight_no_episode(0.75),
      weight_q_prior(0.02),
      weight_plf_prior(0.01),
      // The packaged implementation promotes the stable refinement as its
      // default. The historical direct-update mode remains available only for
      // compatibility and focused regression tests.
      calibration_fit_mode("refined_fixed_point"),
      candidate_direction_mode("single"),
      acceptance_policy("objective_first"),
      class_n_role("objective"),
      duration_target_policy("raw"),
      maximum_duration_ratio(1.0),
      oracle_use_bounds(false),
      include_oracle_plf_in_target(true),
      write_oracle_audit(true),
      oracle_residual_tolerance(1e-7),
      average_speed_feasibility_mph(0.10),
      saturated_speed_tolerance_mph(0.25),
      oracle_plf_weight(2.0),
      max_e_fit_degradation(0.02),
      observed_plf_accumulation(false),
      speed_control_gain(0.20),
      non_episode_speed_control_gain(-1.0),
      episode_external_control_scale(1.0),
      freeze_non_episode_parameters(false),
      s3_control_gain(0.12),
      volume_envelope_control_gain(0.50),
      volume_below_envelope_control_gain(0.80),
      count_control_gain(0.20),
      maximum_log_plf_control(0.40),
      departure_profile_file("departure_profiles.csv"),
      volume_constraint_file(""),
      arrival_profile_output_file(""),
      history_output_file("auto_calibration_history.csv"),
      audit_output_file("auto_calibration_link_audit.csv"),
      volume_constraint_audit_output_file(
          "auto_calibration_volume_constraint_audit.csv"),
      summary_output_file("auto_calibration_summary.json"),
      oracle_audit_output_file("auto_calibration_oracle_audit.csv") {}

AutoCalibrationConfig AutoCalibrationConfig::Load(
    const std::string& path,
    double period_start,
    double period_end) {
    AutoCalibrationConfig config;
    config.period_start_hour = period_start;
    config.period_end_hour = period_end;
    std::ifstream stream(path.c_str());
    if (!stream.is_open()) {
        return config;
    }
    std::string line;
    if (!std::getline(stream, line)) {
        return config;
    }
    const std::vector<std::string> header = SplitCSV(line);
    if (header.size() >= 2 && Lower(header[0]) == "key" &&
        Lower(header[1]) == "value") {
        while (std::getline(stream, line)) {
            const std::vector<std::string> row = SplitCSV(line);
            if (row.size() >= 2) {
                ApplyConfigValue(&config, row[0], row[1]);
            }
        }
    } else if (std::getline(stream, line)) {
        const std::vector<std::string> row = SplitCSV(line);
        for (std::size_t i = 0; i < header.size() && i < row.size(); ++i) {
            ApplyConfigValue(&config, header[i], row[i]);
        }
    }
    config.workers = std::max(1, std::min(20, config.workers));
    config.arrival_tail_minutes = std::max(
        config.departure_bin_minutes, config.arrival_tail_minutes);
    config.minimum_column_od_volume = std::max(
        0.0, config.minimum_column_od_volume);
    config.maximum_column_paths_per_od = std::max(
        1, config.maximum_column_paths_per_od);
    config.minimum_column_path_share = Clamp(
        config.minimum_column_path_share, 0.0, 0.49);
    config.plf_damping = Clamp(config.plf_damping, 0.0, 1.0);
    config.q_damping = Clamp(config.q_damping, 0.0, 1.0);
    config.rejected_step_reduction = Clamp(
        config.rejected_step_reduction, 0.05, 0.95);
    if (!IsRefinedMode(config) &&
        config.calibration_fit_mode != "legacy_regularized") {
        throw std::invalid_argument(
            "Unsupported calibration_fit_mode '" +
            config.calibration_fit_mode +
            "'; use refined_fixed_point");
    }
    return config;
}

CalibrationLink::CalibrationLink()
    : internal_index(-1),
      external_link_id(0),
      from_node_id(0),
      to_node_id(0),
      vdf_code(0),
      calibration_eligible(true),
      observation_class(ObservationClass::Unobserved),
      length_miles(0.0),
      lanes(1.0),
      capacity_vphpl(0.0),
      free_speed_mph(0.0),
      cutoff_speed_mph(0.0),
      volume(0.0),
      vehicle_volume(0.0),
      travel_time_minutes(0.0),
      plf(1.0),
      qcd(1.0),
      qcp(0.28125),
      qn(1.0),
      qs(4.0),
      alpha(0.15),
      beta(4.0),
      mode1_plf(1.0),
      mode1_qcd(1.0),
      mode1_qcp(0.28125),
      observed_duration_hour(std::numeric_limits<double>::quiet_NaN()),
      observed_trough_speed_mph(std::numeric_limits<double>::quiet_NaN()),
      observed_average_speed_mph(std::numeric_limits<double>::quiet_NaN()),
      s3_volume(std::numeric_limits<double>::quiet_NaN()),
      cube_vehicle_volume(std::numeric_limits<double>::quiet_NaN()),
      observation_quality(0.5) {}

CalibrationEvaluation::CalibrationEvaluation()
    : objective(std::numeric_limits<double>::infinity()),
      episode_fit_loss(0.0),
      duration_loss(0.0),
      trough_speed_loss(0.0),
      average_speed_loss(0.0),
      s3_loss(0.0),
      volume_envelope_loss(0.0),
      count_loss(0.0),
      no_episode_loss(0.0),
      prior_loss(0.0),
      route_policy_distance(0.0),
      unobserved_volume_deviation(0.0),
      vmt_change_fraction(0.0),
      vht_change_fraction(0.0),
      relative_gap_pct(std::numeric_limits<double>::infinity()),
      guardrails_pass(false) {}

CalibrationProposal::CalibrationProposal()
    : maximum_parameter_change(0.0), maximum_plf_change(0.0) {}

CalibrationVolumeConstraint::CalibrationVolumeConstraint()
    : target_vehicle_volume(std::numeric_limits<double>::quiet_NaN()),
      weight(1.0) {}

AutoCalibrationEngine::AutoCalibrationEngine(
    const AutoCalibrationConfig& config)
    : config_(config),
      initialized_(false),
      uses_origin_profiles_(false),
      uses_mode_profiles_(false) {}

const char* ObservationClassName(ObservationClass value) {
    switch (value) {
        case ObservationClass::Episode:
            return "E";
        case ObservationClass::NoEpisode:
            return "N";
        default:
            return "U";
    }
}

ObservationClass ParseObservationClass(const std::string& raw_value) {
    const std::string value = Lower(raw_value);
    if (value == "e" || value == "episode" || value == "detected_episode") {
        return ObservationClass::Episode;
    }
    if (value == "n" || value == "no_episode" ||
        value == "no_detected_episode") {
        return ObservationClass::NoEpisode;
    }
    return ObservationClass::Unobserved;
}

bool AutoCalibrationEngine::LoadDepartureProfiles(std::string* error_message) {
    profiles_.clear();
    uses_origin_profiles_ = false;
    uses_mode_profiles_ = false;
    std::ifstream stream(config_.departure_profile_file.c_str());
    if (!stream.is_open()) {
        if (error_message != NULL) {
            *error_message = "Cannot open departure profile file: " +
                config_.departure_profile_file;
        }
        return false;
    }
    std::string line;
    if (!std::getline(stream, line)) {
        if (error_message != NULL) {
            *error_message = "Departure profile file is empty: " +
                config_.departure_profile_file;
        }
        return false;
    }
    const std::vector<std::string> header = SplitCSV(line);
    std::map<std::string, int> index;
    for (std::size_t i = 0; i < header.size(); ++i) {
        index[Lower(header[i])] = static_cast<int>(i);
    }
    const bool long_format =
        (index.count("share") || index.count("probability")) &&
        (index.count("time_min") || index.count("time") ||
         index.count("time_of_day"));
    if (long_format) {
        const int origin_column = index.count("origin_zone_id")
            ? index["origin_zone_id"] : -1;
        const int mode_column = index.count("mode") ? index["mode"] :
            (index.count("mode_type") ? index["mode_type"] : -1);
        const int time_column = index.count("time_min") ? index["time_min"] :
            (index.count("time") ? index["time"] : index["time_of_day"]);
        const int share_column = index.count("share") ? index["share"] :
            index["probability"];
        while (std::getline(stream, line)) {
            const std::vector<std::string> row = SplitCSV(line);
            int origin = 0;
            double minute = 0.0;
            double share = 0.0;
            if (origin_column >= 0 && origin_column < static_cast<int>(row.size())) {
                ParseInt(row[origin_column], &origin);
            }
            const std::string mode = mode_column >= 0 &&
                mode_column < static_cast<int>(row.size())
                ? Lower(row[mode_column]) : "*";
            if (time_column >= static_cast<int>(row.size()) ||
                share_column >= static_cast<int>(row.size())) {
                continue;
            }
            minute = ParseTimeMinute(row[time_column]);
            if (!ParseDouble(row[share_column], &share) ||
                !std::isfinite(minute) || share < 0.0) {
                continue;
            }
            profiles_[ProfileKey(origin, mode)].push_back(
                DepartureProfilePoint{minute, share});
            uses_origin_profiles_ = uses_origin_profiles_ || origin != 0;
            uses_mode_profiles_ = uses_mode_profiles_ || origin == 0;
        }
    } else {
        int origin_column = -1;
        int mode_column = 0;
        if (index.count("origin_zone_id") && index.count("mode")) {
            origin_column = index["origin_zone_id"];
            mode_column = index["mode"];
        } else if (index.count("time_series_label")) {
            mode_column = index["time_series_label"];
        } else if (index.count("mode_type")) {
            mode_column = index["mode_type"];
        }
        std::vector<double> header_minutes(header.size(),
            std::numeric_limits<double>::quiet_NaN());
        for (std::size_t i = 0; i < header.size(); ++i) {
            header_minutes[i] = ParseTimeMinute(header[i]);
        }
        while (std::getline(stream, line)) {
            const std::vector<std::string> row = SplitCSV(line);
            if (mode_column >= static_cast<int>(row.size())) {
                continue;
            }
            int origin = 0;
            if (origin_column >= 0 && origin_column < static_cast<int>(row.size())) {
                ParseInt(row[origin_column], &origin);
            }
            const std::string mode = Lower(row[mode_column]);
            std::vector<DepartureProfilePoint>& profile =
                profiles_[ProfileKey(origin, mode)];
            for (std::size_t i = 0; i < row.size() && i < header.size(); ++i) {
                double share = 0.0;
                if (std::isfinite(header_minutes[i]) &&
                    ParseDouble(row[i], &share) && share >= 0.0) {
                    profile.push_back(
                        DepartureProfilePoint{header_minutes[i], share});
                }
            }
            uses_origin_profiles_ = uses_origin_profiles_ || origin != 0;
            uses_mode_profiles_ = uses_mode_profiles_ || origin == 0;
        }
    }
    for (std::map<ProfileKey, std::vector<DepartureProfilePoint> >::iterator it =
             profiles_.begin(); it != profiles_.end();) {
        std::vector<DepartureProfilePoint>& profile = it->second;
        std::sort(profile.begin(), profile.end(),
                  [](const DepartureProfilePoint& left,
                     const DepartureProfilePoint& right) {
                      return left.minute_of_day < right.minute_of_day;
                  });
        double sum = 0.0;
        for (std::size_t i = 0; i < profile.size(); ++i) {
            sum += profile[i].share;
        }
        if (profile.empty() || sum <= 0.0) {
            profiles_.erase(it++);
        } else {
            ++it;
        }
    }
    if (profiles_.empty()) {
        if (error_message != NULL) {
            *error_message = "No usable departure profiles were found in " +
                config_.departure_profile_file;
        }
        return false;
    }
    return true;
}

bool AutoCalibrationEngine::LoadVolumeConstraints(
    const std::vector<CalibrationLink>& links,
    std::string* error_message) {
    volume_constraints_.clear();
    if (Trim(config_.volume_constraint_file).empty()) {
        return true;
    }
    std::ifstream stream(config_.volume_constraint_file.c_str());
    if (!stream.is_open()) {
        if (error_message != NULL) {
            *error_message = "Cannot open volume constraint file: " +
                config_.volume_constraint_file;
        }
        return false;
    }
    std::string line;
    if (!std::getline(stream, line)) {
        if (error_message != NULL) {
            *error_message = "Volume constraint file is empty: " +
                config_.volume_constraint_file;
        }
        return false;
    }
    const std::vector<std::string> header = SplitCSV(line);
    std::map<std::string, int> column;
    for (std::size_t i = 0; i < header.size(); ++i) {
        column[Lower(header[i])] = static_cast<int>(i);
    }
    const char* required[] = {
        "constraint_id", "target_vehicle_volume", "from_node_id", "to_node_id"};
    for (std::size_t i = 0; i < sizeof(required) / sizeof(required[0]); ++i) {
        if (!column.count(required[i])) {
            if (error_message != NULL) {
                *error_message = std::string("Volume constraint file is missing column '") +
                    required[i] + "': " + config_.volume_constraint_file;
            }
            return false;
        }
    }
    std::map<std::pair<int, int>, std::vector<int> > node_pair_links;
    std::map<int, int> external_link_indices;
    for (std::size_t i = 0; i < links.size(); ++i) {
        node_pair_links[std::make_pair(
            links[i].from_node_id, links[i].to_node_id)].push_back(
                static_cast<int>(i));
        external_link_indices[links[i].external_link_id] = static_cast<int>(i);
    }
    std::map<std::string, int> constraint_indices;
    int row_number = 1;
    while (std::getline(stream, line)) {
        ++row_number;
        if (Trim(line).empty()) {
            continue;
        }
        const std::vector<std::string> row = SplitCSV(line);
        const int id_column = column["constraint_id"];
        if (id_column >= static_cast<int>(row.size())) {
            continue;
        }
        const std::string constraint_id = Trim(row[id_column]);
        int from_node = 0;
        int to_node = 0;
        double target = 0.0;
        if (constraint_id.empty() ||
            column["from_node_id"] >= static_cast<int>(row.size()) ||
            column["to_node_id"] >= static_cast<int>(row.size()) ||
            column["target_vehicle_volume"] >= static_cast<int>(row.size()) ||
            !ParseInt(row[column["from_node_id"]], &from_node) ||
            !ParseInt(row[column["to_node_id"]], &to_node) ||
            !ParseDouble(row[column["target_vehicle_volume"]], &target) ||
            target < 0.0) {
            if (error_message != NULL) {
                std::ostringstream message;
                message << "Invalid volume constraint row " << row_number
                        << " in " << config_.volume_constraint_file;
                *error_message = message.str();
            }
            return false;
        }
        int link_index = -1;
        if (column.count("link_id") &&
            column["link_id"] < static_cast<int>(row.size())) {
            int external_link_id = 0;
            if (ParseInt(row[column["link_id"]], &external_link_id)) {
                std::map<int, int>::const_iterator found =
                    external_link_indices.find(external_link_id);
                if (found != external_link_indices.end()) {
                    link_index = found->second;
                }
            }
        }
        if (link_index < 0) {
            const std::map<std::pair<int, int>, std::vector<int> >::const_iterator
                found = node_pair_links.find(std::make_pair(from_node, to_node));
            if (found != node_pair_links.end() && found->second.size() == 1) {
                link_index = found->second[0];
            }
        }
        if (link_index < 0) {
            if (error_message != NULL) {
                std::ostringstream message;
                message << "Volume constraint row " << row_number
                        << " does not identify one network link ("
                        << from_node << "->" << to_node << ")";
                *error_message = message.str();
            }
            return false;
        }
        int constraint_index = -1;
        std::map<std::string, int>::const_iterator existing =
            constraint_indices.find(constraint_id);
        if (existing == constraint_indices.end()) {
            CalibrationVolumeConstraint item;
            item.constraint_id = constraint_id;
            item.target_vehicle_volume = target;
            if (column.count("constraint_type") &&
                column["constraint_type"] < static_cast<int>(row.size())) {
                item.constraint_type = Trim(row[column["constraint_type"]]);
            }
            if (column.count("weight") &&
                column["weight"] < static_cast<int>(row.size())) {
                double parsed_weight = 0.0;
                if (ParseDouble(row[column["weight"]], &parsed_weight) &&
                    parsed_weight > 0.0) {
                    item.weight = parsed_weight;
                }
            }
            volume_constraints_.push_back(item);
            constraint_index = static_cast<int>(volume_constraints_.size()) - 1;
            constraint_indices[constraint_id] = constraint_index;
        } else {
            constraint_index = existing->second;
            const double previous =
                volume_constraints_[constraint_index].target_vehicle_volume;
            if (std::fabs(previous - target) >
                1e-8 * std::max(1.0, std::fabs(previous))) {
                if (error_message != NULL) {
                    *error_message = "Conflicting targets for volume constraint '" +
                        constraint_id + "'";
                }
                return false;
            }
        }
        double coefficient = 1.0;
        if (column.count("coefficient") &&
            column["coefficient"] < static_cast<int>(row.size())) {
            ParseDouble(row[column["coefficient"]], &coefficient);
        }
        if (!std::isfinite(coefficient) || coefficient == 0.0) {
            if (error_message != NULL) {
                *error_message = "Invalid zero/nonfinite coefficient for volume constraint '" +
                    constraint_id + "'";
            }
            return false;
        }
        volume_constraints_[constraint_index].link_indices.push_back(link_index);
        volume_constraints_[constraint_index].coefficients.push_back(coefficient);
    }
    return true;
}

const std::vector<DepartureProfilePoint>*
AutoCalibrationEngine::FindProfile(
    int origin_zone_id,
    const std::string& raw_mode_name) const {
    const std::string mode_name = Lower(raw_mode_name);
    const ProfileKey keys[] = {
        ProfileKey(origin_zone_id, mode_name),
        ProfileKey(origin_zone_id, "*"),
        ProfileKey(0, mode_name),
        ProfileKey(0, "*")};
    for (std::size_t i = 0; i < sizeof(keys) / sizeof(keys[0]); ++i) {
        std::map<ProfileKey,
                 std::vector<DepartureProfilePoint> >::const_iterator found =
            profiles_.find(keys[i]);
        if (found != profiles_.end()) {
            return &found->second;
        }
    }
    return NULL;
}

const std::vector<DepartureProfilePoint>*
AutoCalibrationEngine::DepartureProfileFor(
    int origin_zone_id,
    const std::string& mode_name) const {
    return FindProfile(origin_zone_id, mode_name);
}

void AutoCalibrationEngine::Initialize(
    const std::vector<CalibrationLink>& links,
    const std::vector<CalibrationRoute>& routes,
    double relative_gap_pct) {
    Initialize(links, SummarizeRoutes(links, routes), relative_gap_pct);
}

void AutoCalibrationEngine::Initialize(
    const std::vector<CalibrationLink>& links,
    const CalibrationRouteSummary& route_summary,
    double relative_gap_pct) {
    baseline_links_ = links;
    baseline_policy_weights_ = route_summary.policy_weights;
    initialized_ = true;
    accepted_evaluation_ = Evaluate(links, route_summary, relative_gap_pct);
    history_.clear();
    CalibrationProposal baseline_proposal;
    Record(0, 0, true, 0.0, accepted_evaluation_, baseline_proposal);
}

std::vector<double> AutoCalibrationEngine::BuildDeparturePLF(
    const std::vector<CalibrationLink>& links,
    const std::vector<CalibrationRoute>& routes) const {
    std::vector<double> result(links.size(), 1.0);
    if (links.empty() || routes.empty() || profiles_.empty()) {
        for (std::size_t i = 0; i < links.size(); ++i) {
            result[i] = Clamp(IsFinitePositive(links[i].plf) ? links[i].plf : 1.0,
                              config_.minimum_plf, config_.maximum_plf);
        }
        return result;
    }
    const double start_minute = config_.period_start_hour * 60.0;
    double end_minute = config_.period_end_hour * 60.0;
    if (end_minute <= start_minute) {
        end_minute += 1440.0;
    }
    const double duration_minutes = end_minute - start_minute;
    double maximum_route_minutes = 0.0;
    for (std::size_t i = 0; i < routes.size(); ++i) {
        double route_minutes = 0.0;
        for (std::size_t j = 0; j < routes[i].link_travel_times_minutes.size(); ++j) {
            route_minutes += std::max(0.0, routes[i].link_travel_times_minutes[j]);
        }
        maximum_route_minutes = std::max(maximum_route_minutes, route_minutes);
    }
    const double bin_minutes = std::max(1.0, config_.departure_bin_minutes);
    const int horizon_bins = std::max(
        1, static_cast<int>(std::ceil(
               (duration_minutes + maximum_route_minutes + bin_minutes) /
               bin_minutes)));
    std::vector<double> arrival(links.size() * static_cast<std::size_t>(horizon_bins),
                                0.0);

#ifdef _OPENMP
#pragma omp parallel for num_threads(config_.workers) schedule(dynamic, 16)
#endif
    for (int route_index = 0; route_index < static_cast<int>(routes.size());
         ++route_index) {
        const CalibrationRoute& route = routes[route_index];
        const std::vector<DepartureProfilePoint>* profile =
            FindProfile(route.origin_zone_id, route.mode_name);
        if (profile == NULL || route.link_indices.empty()) {
            continue;
        }
        double period_share = 0.0;
        for (std::size_t p = 0; p < profile->size(); ++p) {
            const double minute = UnwrapMinute(
                (*profile)[p].minute_of_day, start_minute, end_minute);
            if (minute >= start_minute && minute < end_minute) {
                period_share += (*profile)[p].share;
            }
        }
        if (period_share <= 0.0) {
            continue;
        }
        const double route_flow = std::max(0.0, route.od_demand) *
            std::max(0.0, route.pce) * std::max(0.0, route.share);
        if (route_flow <= 0.0) {
            continue;
        }
        // The supplied mode profiles contain a full day (96 bins). Filter and
        // normalize the active period once per route rather than rescanning the
        // full profile for every downstream link on that route.
        std::vector<double> period_departure_minutes;
        std::vector<double> period_departure_flows;
        period_departure_minutes.reserve(profile->size());
        period_departure_flows.reserve(profile->size());
        for (std::size_t p = 0; p < profile->size(); ++p) {
            const double departure_minute = UnwrapMinute(
                (*profile)[p].minute_of_day, start_minute, end_minute);
            if (departure_minute < start_minute ||
                departure_minute >= end_minute) {
                continue;
            }
            period_departure_minutes.push_back(departure_minute);
            period_departure_flows.push_back(
                route_flow * (*profile)[p].share / period_share);
        }
        double cumulative_minutes = 0.0;
        for (std::size_t position = 0; position < route.link_indices.size();
             ++position) {
            const int link_index = route.link_indices[position];
            if (link_index >= 0 && link_index < static_cast<int>(links.size())) {
                for (std::size_t p = 0;
                     p < period_departure_minutes.size(); ++p) {
                    const int bin = static_cast<int>(std::floor(
                        (period_departure_minutes[p] + cumulative_minutes -
                         start_minute) /
                        bin_minutes));
                    if (bin < 0 || bin >= horizon_bins) {
                        continue;
                    }
                    const std::size_t offset =
                        static_cast<std::size_t>(link_index) * horizon_bins + bin;
#ifdef _OPENMP
#pragma omp atomic
#endif
                    arrival[offset] += period_departure_flows[p];
                }
            }
            if (position < route.link_travel_times_minutes.size()) {
                cumulative_minutes += std::max(
                    0.0, route.link_travel_times_minutes[position]);
            }
        }
    }

    const double period_hours = PeriodDurationHours(
        config_.period_start_hour, config_.period_end_hour);
    const double bin_hours = bin_minutes / 60.0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(config_.workers) schedule(static)
#endif
    for (int link_index = 0; link_index < static_cast<int>(links.size());
         ++link_index) {
        double total = 0.0;
        double peak = 0.0;
        const std::size_t offset =
            static_cast<std::size_t>(link_index) * horizon_bins;
        for (int bin = 0; bin < horizon_bins; ++bin) {
            const double value = arrival[offset + bin];
            total += value;
            peak = std::max(peak, value);
        }
        double plf = IsFinitePositive(links[link_index].plf)
            ? links[link_index].plf : 1.0;
        if (total > kEpsilon && peak > kEpsilon) {
            const double peak_shape_share = peak / total;
            plf = bin_hours / (period_hours * peak_shape_share);
        }
        result[link_index] = Clamp(
            plf, config_.minimum_plf, config_.maximum_plf);
    }
    return result;
}

CalibrationProposal AutoCalibrationEngine::Propose(
    const std::vector<CalibrationLink>& links,
    const std::vector<CalibrationRoute>& routes,
    double step_scale) const {
    return Propose(links, SummarizeRoutes(links, routes), step_scale);
}

CalibrationProposal AutoCalibrationEngine::Propose(
    const std::vector<CalibrationLink>& links,
    const CalibrationRouteSummary& route_summary,
    double step_scale) const {
    CalibrationProposal proposal;
    proposal.departure_plf = route_summary.departure_plf;
    if (proposal.departure_plf.size() != links.size()) {
        proposal.departure_plf.resize(links.size());
        for (std::size_t i = 0; i < links.size(); ++i) {
            proposal.departure_plf[i] = IsFinitePositive(links[i].plf)
                ? links[i].plf : 1.0;
        }
    }
    proposal.plf.resize(links.size());
    proposal.qcd.resize(links.size());
    proposal.qcp.resize(links.size());
    proposal.alpha.resize(links.size());
    proposal.beta.resize(links.size());
    proposal.oracle.resize(links.size());
    std::map<int, std::vector<double> > group_plf_values;
    for (std::size_t i = 0; i < links.size(); ++i) {
        const double prior = IsFinitePositive(links[i].mode1_plf)
            ? links[i].mode1_plf : links[i].plf;
        group_plf_values[links[i].vdf_code].push_back(prior);
    }
    std::map<int, double> group_plf;
    for (std::map<int, std::vector<double> >::const_iterator it =
             group_plf_values.begin(); it != group_plf_values.end(); ++it) {
        group_plf[it->first] = Median(it->second, 1.0);
    }
    const double period_hours = PeriodDurationHours(
        config_.period_start_hour, config_.period_end_hour);
    const double plf_step = Clamp(config_.plf_damping * step_scale, 0.0, 1.0);
    const double q_step = Clamp(config_.q_damping * step_scale, 0.0, 1.0);
    const bool refined_mode = IsRefinedMode(config_);
    const FixedVolumeOracleConfig oracle_config = BuildOracleConfig(config_);
    std::vector<double> count_log_control(links.size(), 0.0);
    std::vector<double> count_control_weight(links.size(), 0.0);
    for (std::size_t c = 0; c < volume_constraints_.size(); ++c) {
        const CalibrationVolumeConstraint& constraint = volume_constraints_[c];
        double modeled = 0.0;
        for (std::size_t p = 0; p < constraint.link_indices.size(); ++p) {
            const int link_index = constraint.link_indices[p];
            if (link_index >= 0 && link_index < static_cast<int>(links.size())) {
                modeled += constraint.coefficients[p] *
                    std::max(0.0, links[link_index].vehicle_volume);
            }
        }
        if (!IsFinitePositive(constraint.target_vehicle_volume)) {
            continue;
        }
        const double residual = SafeLog(
            (std::max(0.0, modeled) + 1.0) /
            (constraint.target_vehicle_volume + 1.0));
        const double weight = std::max(0.0, constraint.weight);
        const double shift = -std::max(0.0, config_.count_control_gain) *
            Clamp(residual, -3.0, 3.0);
        for (std::size_t p = 0; p < constraint.link_indices.size(); ++p) {
            const int link_index = constraint.link_indices[p];
            if (link_index >= 0 && link_index < static_cast<int>(links.size())) {
                count_log_control[link_index] += weight * shift;
                count_control_weight[link_index] += weight;
            }
        }
    }
    for (std::size_t i = 0; i < count_log_control.size(); ++i) {
        if (count_control_weight[i] > kEpsilon) {
            count_log_control[i] /= count_control_weight[i];
        }
    }
#ifdef _OPENMP
#pragma omp parallel for num_threads(config_.workers) schedule(static)
#endif
    for (int i = 0; i < static_cast<int>(links.size()); ++i) {
        const CalibrationLink& link = links[i];
        const double current_plf = IsFinitePositive(link.plf) ? link.plf : 1.0;
        const double mode1_plf = IsFinitePositive(link.mode1_plf)
            ? link.mode1_plf : current_plf;
        const double departure_plf = IsFinitePositive(proposal.departure_plf[i])
            ? proposal.departure_plf[i] : current_plf;
        const double grouped_plf = group_plf[link.vdf_code];
        const double total_weight = std::max(
            kEpsilon, config_.departure_plf_weight +
                config_.mode1_plf_weight + config_.group_plf_weight);
        const double target_log_plf =
            (config_.departure_plf_weight * SafeLog(departure_plf) +
             config_.mode1_plf_weight * SafeLog(mode1_plf) +
             config_.group_plf_weight * SafeLog(grouped_plf)) /
            total_weight;
        double next_plf = Clamp(
            std::exp((1.0 - plf_step) * SafeLog(current_plf) +
                     plf_step * target_log_plf),
            config_.minimum_plf, config_.maximum_plf);
        if (config_.observed_plf_accumulation &&
            link.calibration_eligible &&
            link.observation_class == ObservationClass::NoEpisode &&
            IsFinitePositive(link.observed_average_speed_mph)) {
            // Keep the last accepted state as the controller base. In the
            // original refinement loop, every N-link proposal was pulled back
            // toward the fixed departure/mode-1 prior before applying its
            // speed residual, preventing the correction from accumulating.
            next_plf = current_plf;
        }

        double next_qcd = link.qcd;
        double next_qcp = link.qcp;
        FixedVolumeOracleResult oracle;
        if (refined_mode && link.calibration_eligible &&
            link.observation_class == ObservationClass::Episode) {
            const FixedVolumeOracleInput oracle_input = BuildOracleInput(
                link, period_hours, config_);
            oracle = SolveFixedVolumeOracle(oracle_input, oracle_config);
            const bool usable_oracle = IsFinitePositive(oracle.plf_applied) &&
                IsFinitePositive(oracle.qcd_applied) &&
                IsFinitePositive(oracle.qcp_applied);
            if (usable_oracle) {
                double oracle_target_plf = oracle.plf_applied;
                if (config_.calibration_fit_mode == "equilibrium_regularized" &&
                    config_.candidate_direction_mode != "exact_e") {
                    const double legacy_weight = std::max(
                        kEpsilon, config_.departure_plf_weight +
                            config_.mode1_plf_weight + config_.group_plf_weight);
                    oracle_target_plf = std::exp(
                        (legacy_weight * target_log_plf +
                         std::max(0.0, config_.oracle_plf_weight) *
                             SafeLog(oracle.plf_applied)) /
                        std::max(kEpsilon, legacy_weight +
                            std::max(0.0, config_.oracle_plf_weight)));
                }
                if (config_.include_oracle_plf_in_target) {
                    next_plf = std::exp(
                        (1.0 - plf_step) * SafeLog(current_plf) +
                        plf_step * SafeLog(oracle_target_plf));
                    if (config_.oracle_use_bounds) {
                        next_plf = Clamp(next_plf,
                            config_.minimum_plf, config_.maximum_plf);
                    }
                }
                next_qcd = std::exp(
                    (1.0 - q_step) * SafeLog(
                        IsFinitePositive(link.qcd) ? link.qcd : oracle.qcd_applied) +
                    q_step * SafeLog(oracle.qcd_applied));
                next_qcp = std::exp(
                    (1.0 - q_step) * SafeLog(
                        IsFinitePositive(link.qcp) ? link.qcp : oracle.qcp_applied) +
                    q_step * SafeLog(oracle.qcp_applied));
            }
        }
        if (link.calibration_eligible) {
            double log_control = count_log_control[i];
            const double modeled_average = ModeledAverageSpeed(link);
            if (IsFinitePositive(link.observed_average_speed_mph) &&
                IsFinitePositive(modeled_average)) {
                const double speed_residual =
                    (modeled_average - link.observed_average_speed_mph) /
                    std::max(config_.speed_scale_mph, 0.1);
                // Faster-than-observed links receive a lower PLF and therefore
                // a higher DOC/cost; slower links receive the opposite pressure.
                const double configured_speed_gain =
                    link.observation_class == ObservationClass::NoEpisode &&
                            config_.non_episode_speed_control_gain >= 0.0
                        ? config_.non_episode_speed_control_gain
                        : config_.speed_control_gain;
                log_control -= std::max(0.0, configured_speed_gain) *
                    Clamp(speed_residual, -3.0, 3.0);
            }
            VolumeEnvelopePosition envelope_position =
                VolumeEnvelopePosition::InsideOrSingle;
            const double volume_residual = VolumeDirectionResidual(
                link, &envelope_position);
            const double volume_gain =
                envelope_position == VolumeEnvelopePosition::Below
                ? std::max(0.0, config_.volume_below_envelope_control_gain)
                : envelope_position == VolumeEnvelopePosition::Above
                ? std::max(0.0, config_.volume_envelope_control_gain)
                : std::max(0.0, config_.s3_control_gain);
            log_control -= volume_gain * Clamp(volume_residual, -3.0, 3.0);
            const double maximum_control = std::max(
                0.0, config_.maximum_log_plf_control);
            if (link.observation_class == ObservationClass::Episode) {
                log_control *= std::max(
                    0.0, config_.episode_external_control_scale);
            }
            log_control = Clamp(log_control, -maximum_control, maximum_control);
            next_plf = Clamp(
                next_plf * std::exp(plf_step * log_control),
                config_.minimum_plf, config_.maximum_plf);
        }
        if (link.observation_class == ObservationClass::Episode &&
            (!refined_mode || !IsFinitePositive(oracle.qcd_applied)) &&
            IsFinitePositive(EffectiveDurationTarget(link, period_hours, config_)) &&
            IsFinitePositive(link.observed_trough_speed_mph) &&
            IsFinitePositive(link.capacity_vphpl) &&
            IsFinitePositive(link.lanes) && IsFinitePositive(link.volume) &&
            IsFinitePositive(link.qn) && IsFinitePositive(link.qs)) {
            const double doc = link.volume /
                (link.lanes * period_hours * next_plf * link.capacity_vphpl);
            if (IsFinitePositive(doc)) {
                const double effective_duration = EffectiveDurationTarget(
                    link, period_hours, config_);
                const double raw_qcd = effective_duration /
                    std::pow(doc, link.qn);
                const double severity = link.cutoff_speed_mph /
                    link.observed_trough_speed_mph - 1.0;
                const double raw_qcp = severity > 0.0
                    ? severity /
                        std::pow(effective_duration, link.qs)
                    : std::numeric_limits<double>::quiet_NaN();
                const double quality = Clamp(
                    std::isfinite(link.observation_quality)
                        ? link.observation_quality : 0.5,
                    0.0, 1.0);
                if (IsFinitePositive(raw_qcd)) {
                    const double prior = IsFinitePositive(link.mode1_qcd)
                        ? link.mode1_qcd : link.qcd;
                    const double target = Clamp(std::exp(
                        (1.0 - quality) * SafeLog(prior) +
                        quality * SafeLog(raw_qcd)),
                        config_.minimum_qcd, config_.maximum_qcd);
                    // Some legacy networks legitimately start outside the
                    // configured calibration bounds.  Interpolate from that
                    // accepted state toward the bounded target instead of
                    // snapping to the bound: otherwise every trust-region
                    // retry reports the same large relative Q change.
                    const double current = IsFinitePositive(link.qcd)
                        ? link.qcd : config_.minimum_qcd;
                    next_qcd = std::exp(
                        (1.0 - q_step) * SafeLog(current) +
                        q_step * SafeLog(target));
                    next_qcd = Clamp(next_qcd,
                        std::min(current, config_.minimum_qcd),
                        std::max(current, config_.maximum_qcd));
                }
                if (IsFinitePositive(raw_qcp)) {
                    const double prior = IsFinitePositive(link.mode1_qcp)
                        ? link.mode1_qcp : link.qcp;
                    const double target = Clamp(std::exp(
                        (1.0 - quality) * SafeLog(prior) +
                        quality * SafeLog(raw_qcp)),
                        config_.minimum_qcp, config_.maximum_qcp);
                    const double current = IsFinitePositive(link.qcp)
                        ? link.qcp : config_.minimum_qcp;
                    next_qcp = std::exp(
                        (1.0 - q_step) * SafeLog(current) +
                        q_step * SafeLog(target));
                    next_qcp = Clamp(next_qcp,
                        std::min(current, config_.minimum_qcp),
                        std::max(current, config_.maximum_qcp));
                }
            }
        }
        const bool freeze_non_episode =
            config_.freeze_non_episode_parameters &&
            link.observation_class == ObservationClass::NoEpisode;
        proposal.plf[i] = freeze_non_episode ? link.plf : next_plf;
        proposal.qcd[i] = freeze_non_episode ? link.qcd : next_qcd;
        proposal.qcp[i] = freeze_non_episode ? link.qcp : next_qcp;
        proposal.alpha[i] = freeze_non_episode
            ? link.alpha
            : config_.theta * next_qcp * std::pow(next_qcd, link.qs);
        proposal.beta[i] = freeze_non_episode
            ? link.beta
            : link.qn * link.qs;
        proposal.oracle[i] = oracle;
    }
    for (std::size_t i = 0; i < links.size(); ++i) {
        proposal.maximum_plf_change = std::max(
            proposal.maximum_plf_change,
            std::fabs(std::exp(SafeLog(proposal.plf[i]) -
                                SafeLog(links[i].plf)) - 1.0));
        proposal.maximum_parameter_change = std::max(
            proposal.maximum_parameter_change,
            std::max(
                std::fabs(std::exp(SafeLog(proposal.qcd[i]) -
                                    SafeLog(links[i].qcd)) - 1.0),
                std::max(
                    std::fabs(std::exp(SafeLog(proposal.qcp[i]) -
                                        SafeLog(links[i].qcp)) - 1.0),
                    std::fabs(std::exp(SafeLog(proposal.alpha[i]) -
                                        SafeLog(links[i].alpha)) - 1.0))));
    }
    return proposal;
}

CalibrationRouteSummary AutoCalibrationEngine::SummarizeRoutes(
    const std::vector<CalibrationLink>& links,
    const std::vector<CalibrationRoute>& routes) const {
    CalibrationRouteSummary summary;
    summary.departure_plf = BuildDeparturePLF(links, routes);
    int maximum_mode = 0;
    for (std::size_t i = 0; i < routes.size(); ++i) {
        maximum_mode = std::max(maximum_mode, routes[i].mode_index);
    }
    summary.policy_weights.assign(
        static_cast<std::size_t>(maximum_mode + 1) * links.size(), 0.0);
    for (std::size_t i = 0; i < routes.size(); ++i) {
        const CalibrationRoute& route = routes[i];
        const double flow = std::max(0.0, route.od_demand) *
            std::max(0.0, route.pce) * std::max(0.0, route.share);
        if (flow <= 0.0 || route.mode_index < 0) {
            continue;
        }
        for (std::size_t p = 0; p < route.link_indices.size(); ++p) {
            const int link_index = route.link_indices[p];
            if (link_index >= 0 && link_index < static_cast<int>(links.size())) {
                summary.policy_weights[
                    static_cast<std::size_t>(route.mode_index) * links.size() +
                    static_cast<std::size_t>(link_index)] += flow;
            }
        }
    }
    return summary;
}

double AutoCalibrationEngine::RoutePolicyDistance(
    const CalibrationRouteSummary& summary) const {
    if (baseline_policy_weights_.empty() ||
        baseline_policy_weights_.size() != summary.policy_weights.size()) {
        return 0.0;
    }
    double baseline_total = 0.0;
    double current_total = 0.0;
    double overlap = 0.0;
    for (std::size_t i = 0; i < baseline_policy_weights_.size(); ++i) {
        const double baseline = std::max(0.0, baseline_policy_weights_[i]);
        const double current = std::max(0.0, summary.policy_weights[i]);
        baseline_total += baseline;
        current_total += current;
        overlap += std::min(baseline, current);
    }
    const double denominator = std::max(baseline_total, current_total);
    return denominator > kEpsilon ? Clamp(1.0 - overlap / denominator, 0.0, 1.0)
                                  : 0.0;
}

CalibrationEvaluation AutoCalibrationEngine::Evaluate(
    const std::vector<CalibrationLink>& links,
    const std::vector<CalibrationRoute>& routes,
    double relative_gap_pct) const {
    return Evaluate(links, SummarizeRoutes(links, routes), relative_gap_pct);
}

CalibrationEvaluation AutoCalibrationEngine::Evaluate(
    const std::vector<CalibrationLink>& links,
    const CalibrationRouteSummary& route_summary,
    double relative_gap_pct) const {
    CalibrationEvaluation evaluation;
    evaluation.relative_gap_pct = relative_gap_pct;
    const double period_hours = PeriodDurationHours(
        config_.period_start_hour, config_.period_end_hour);
    double duration_sum = 0.0;
    double trough_sum = 0.0;
    double episode_average_sum = 0.0;
    double average_sum = 0.0;
    double s3_sum = 0.0;
    double volume_envelope_sum = 0.0;
    double no_episode_sum = 0.0;
    double prior_sum = 0.0;
    int duration_count = 0;
    int trough_count = 0;
    int episode_average_count = 0;
    int average_count = 0;
    int s3_count = 0;
    int volume_envelope_count = 0;
    int no_episode_count = 0;
    int prior_count = 0;
    double unobserved_deviation = 0.0;
    int unobserved_count = 0;
    double baseline_vmt = 0.0;
    double current_vmt = 0.0;
    double baseline_vht = 0.0;
    double current_vht = 0.0;

    for (std::size_t i = 0; i < links.size(); ++i) {
        const CalibrationLink& link = links[i];
        const double doc_denominator = link.lanes * period_hours * link.plf *
            link.capacity_vphpl;
        const double doc = doc_denominator > kEpsilon
            ? link.volume / doc_denominator : 0.0;
        const double modeled_duration = IsFinitePositive(doc) &&
            IsFinitePositive(link.qcd) && IsFinitePositive(link.qn)
            ? link.qcd * std::pow(doc, link.qn) : 0.0;
        const double modeled_trough = IsFinitePositive(link.cutoff_speed_mph)
            ? link.cutoff_speed_mph /
                (1.0 + link.qcp * std::pow(std::max(0.0, modeled_duration), link.qs))
            : 0.0;
        const double modeled_average = ModeledAverageSpeed(link);

        if (link.calibration_eligible &&
            link.observation_class == ObservationClass::Episode) {
            const double effective_duration = EffectiveDurationTarget(
                link, period_hours, config_);
            if (IsFinitePositive(effective_duration) &&
                IsFinitePositive(modeled_duration)) {
                duration_sum += Huber(
                    SafeLog(modeled_duration / effective_duration),
                    config_.huber_delta);
                ++duration_count;
            }
            if (IsFinitePositive(link.observed_trough_speed_mph) &&
                IsFinitePositive(modeled_trough)) {
                trough_sum += Huber(
                    (modeled_trough - link.observed_trough_speed_mph) /
                        std::max(config_.speed_scale_mph, 0.1),
                    config_.huber_delta);
                ++trough_count;
            }
            if (IsFinitePositive(link.observed_average_speed_mph) &&
                IsFinitePositive(modeled_average)) {
                episode_average_sum += Huber(
                    (modeled_average - link.observed_average_speed_mph) /
                        std::max(config_.speed_scale_mph, 0.1),
                    config_.huber_delta);
                ++episode_average_count;
            }
        }
        const bool class_n_in_objective =
            config_.class_n_role != "report_only";
        if (link.calibration_eligible &&
            IsFinitePositive(link.observed_average_speed_mph) &&
            IsFinitePositive(modeled_average)) {
            average_sum += Huber(
                (modeled_average - link.observed_average_speed_mph) /
                    std::max(config_.speed_scale_mph, 0.1),
                config_.huber_delta);
            ++average_count;
        }
        if (link.calibration_eligible &&
            IsFinitePositive(link.s3_volume)) {
            s3_sum += Huber(
                SafeLog((link.vehicle_volume + 1.0) /
                        (link.s3_volume + 1.0)),
                config_.huber_delta);
            ++s3_count;
        }
        if (link.calibration_eligible &&
            IsFinitePositive(link.s3_volume) &&
            IsFinitePositive(link.cube_vehicle_volume)) {
            VolumeEnvelopePosition envelope_position =
                VolumeEnvelopePosition::InsideOrSingle;
            const double residual = VolumeDirectionResidual(
                link, &envelope_position);
            if (envelope_position != VolumeEnvelopePosition::InsideOrSingle) {
                volume_envelope_sum += Huber(
                    residual, config_.huber_delta);
                ++volume_envelope_count;
            }
        }
        if (link.calibration_eligible && class_n_in_objective &&
            link.observation_class == ObservationClass::NoEpisode) {
            const bool modeled_episode =
                modeled_duration >= config_.no_episode_min_duration_hour &&
                modeled_trough < link.cutoff_speed_mph;
            if (modeled_episode) {
                const double duration_ratio = modeled_duration /
                    std::max(config_.no_episode_min_duration_hour, 0.01);
                const double speed_ratio =
                    (link.cutoff_speed_mph - modeled_trough) /
                    std::max(config_.speed_scale_mph, 0.1);
                no_episode_sum += Huber(duration_ratio, config_.huber_delta) +
                    Huber(speed_ratio, config_.huber_delta);
            }
            ++no_episode_count;
        }
        if (IsFinitePositive(link.mode1_qcd) && IsFinitePositive(link.qcd)) {
            const double delta = SafeLog(link.qcd / link.mode1_qcd);
            prior_sum += delta * delta;
            ++prior_count;
        }
        if (IsFinitePositive(link.mode1_qcp) && IsFinitePositive(link.qcp)) {
            const double delta = SafeLog(link.qcp / link.mode1_qcp);
            prior_sum += delta * delta;
            ++prior_count;
        }
        if (IsFinitePositive(link.mode1_plf) && IsFinitePositive(link.plf)) {
            const double delta = SafeLog(link.plf / link.mode1_plf);
            prior_sum += config_.weight_plf_prior * delta * delta;
            ++prior_count;
        }
        if (i < baseline_links_.size()) {
            const CalibrationLink& baseline = baseline_links_[i];
            if (link.observation_class == ObservationClass::Unobserved &&
                link.volume >= 0.0 && baseline.volume >= 0.0) {
                unobserved_deviation += std::fabs(
                    SafeLog((link.volume + 1.0) / (baseline.volume + 1.0)));
                ++unobserved_count;
            }
            baseline_vmt += baseline.volume * baseline.length_miles;
            current_vmt += link.volume * link.length_miles;
            baseline_vht += baseline.volume * baseline.travel_time_minutes / 60.0;
            current_vht += link.volume * link.travel_time_minutes / 60.0;
        }
    }
    evaluation.duration_loss = duration_count > 0
        ? duration_sum / duration_count : 0.0;
    evaluation.trough_speed_loss = trough_count > 0
        ? trough_sum / trough_count : 0.0;
    evaluation.average_speed_loss = average_count > 0
        ? average_sum / average_count : 0.0;
    evaluation.s3_loss = s3_count > 0 ? s3_sum / s3_count : 0.0;
    evaluation.volume_envelope_loss = volume_envelope_count > 0
        ? volume_envelope_sum / volume_envelope_count : 0.0;
    const int episode_fit_count = duration_count + trough_count +
        episode_average_count;
    evaluation.episode_fit_loss = episode_fit_count > 0
        ? (duration_sum + trough_sum + episode_average_sum) /
            episode_fit_count : 0.0;
    evaluation.no_episode_loss = no_episode_count > 0
        ? no_episode_sum / no_episode_count : 0.0;
    evaluation.prior_loss = prior_count > 0 ? prior_sum / prior_count : 0.0;
    double count_sum = 0.0;
    double count_weight_sum = 0.0;
    for (std::size_t c = 0; c < volume_constraints_.size(); ++c) {
        const CalibrationVolumeConstraint& constraint = volume_constraints_[c];
        double modeled = 0.0;
        for (std::size_t p = 0; p < constraint.link_indices.size(); ++p) {
            const int link_index = constraint.link_indices[p];
            if (link_index >= 0 && link_index < static_cast<int>(links.size())) {
                modeled += constraint.coefficients[p] *
                    std::max(0.0, links[link_index].vehicle_volume);
            }
        }
        if (IsFinitePositive(constraint.target_vehicle_volume)) {
            const double weight = std::max(0.0, constraint.weight);
            count_sum += weight * Huber(
                SafeLog((std::max(0.0, modeled) + 1.0) /
                    (constraint.target_vehicle_volume + 1.0)),
                config_.huber_delta);
            count_weight_sum += weight;
        }
    }
    evaluation.count_loss = count_weight_sum > kEpsilon
        ? count_sum / count_weight_sum : 0.0;
    evaluation.objective =
        config_.weight_duration * evaluation.duration_loss +
        config_.weight_trough_speed * evaluation.trough_speed_loss +
        config_.weight_average_speed * evaluation.average_speed_loss +
        config_.weight_s3 * evaluation.s3_loss +
        config_.weight_volume_envelope * evaluation.volume_envelope_loss +
        config_.weight_count * evaluation.count_loss +
        config_.weight_no_episode * evaluation.no_episode_loss +
        config_.weight_q_prior * evaluation.prior_loss;
    evaluation.route_policy_distance = RoutePolicyDistance(route_summary);
    evaluation.unobserved_volume_deviation = unobserved_count > 0
        ? unobserved_deviation / unobserved_count : 0.0;
    evaluation.vmt_change_fraction = baseline_vmt > kEpsilon
        ? std::fabs(current_vmt / baseline_vmt - 1.0) : 0.0;
    evaluation.vht_change_fraction = baseline_vht > kEpsilon
        ? std::fabs(current_vht / baseline_vht - 1.0) : 0.0;
    evaluation.guardrails_pass =
        evaluation.route_policy_distance <= config_.route_policy_tolerance &&
        evaluation.unobserved_volume_deviation <=
            config_.unobserved_volume_tolerance &&
        evaluation.vmt_change_fraction <= config_.system_vmt_tolerance &&
        evaluation.vht_change_fraction <= config_.system_vht_tolerance &&
        (!std::isfinite(config_.inner_gap_tolerance_pct) ||
         config_.inner_gap_tolerance_pct <= 0.0 ||
         (std::isfinite(relative_gap_pct) && relative_gap_pct >= 0.0 &&
          relative_gap_pct <= config_.inner_gap_tolerance_pct));
    return evaluation;
}

bool AutoCalibrationEngine::ShouldAccept(
    const CalibrationEvaluation& evaluation) const {
    if (!evaluation.guardrails_pass || !std::isfinite(evaluation.objective)) {
        return false;
    }
    if (!std::isfinite(accepted_evaluation_.objective)) {
        return true;
    }
    if (config_.acceptance_policy == "lexicographic_e_fit") {
        const double primary_scale = std::max(
            std::fabs(accepted_evaluation_.episode_fit_loss), 1e-9);
        const double primary_noise = std::max(
            0.0, config_.objective_relative_tolerance) * primary_scale;
        const double maximum_primary = accepted_evaluation_.episode_fit_loss +
            std::max(primary_noise,
                std::max(0.0, config_.max_e_fit_degradation) * primary_scale);
        if (!std::isfinite(evaluation.episode_fit_loss) ||
            evaluation.episode_fit_loss > maximum_primary) {
            return false;
        }
        // A clear improvement in P/vt2/period-speed fit wins.  Inside the
        // primary numerical band, use the all-GP speed/volume/count objective
        // as the second lexicographic level.
        if (evaluation.episode_fit_loss <
            accepted_evaluation_.episode_fit_loss - primary_noise) {
            return true;
        }
    }
    // Permit only the configured numerical-noise band.  This avoids six full
    // trust-region retries for sub-roundoff equilibrium differences while the
    // convergence test still prevents the band from accumulating indefinitely.
    const double noise_band = std::max(0.0, config_.objective_relative_tolerance) *
        std::max(std::fabs(accepted_evaluation_.objective), 1e-9);
    return evaluation.objective <= accepted_evaluation_.objective + noise_band;
}

bool AutoCalibrationEngine::Converged(
    const CalibrationEvaluation& evaluation,
    const CalibrationProposal& proposal,
    int accepted_outer_iterations) const {
    if (accepted_outer_iterations < config_.minimum_outer_iterations ||
        !evaluation.guardrails_pass) {
        return false;
    }
    const double denominator = std::max(
        std::fabs(accepted_evaluation_.objective), 1e-9);
    const double objective_change = std::fabs(
        evaluation.objective - accepted_evaluation_.objective) / denominator;
    return objective_change <= config_.objective_relative_tolerance &&
        proposal.maximum_parameter_change <=
            config_.parameter_relative_tolerance &&
        proposal.maximum_plf_change <= config_.plf_relative_tolerance;
}

void AutoCalibrationEngine::Record(
    int outer_iteration,
    int retry,
    bool accepted,
    double step_scale,
    const CalibrationEvaluation& evaluation,
    const CalibrationProposal& proposal) {
    CalibrationIteration item;
    item.outer_iteration = outer_iteration;
    item.retry = retry;
    item.accepted = accepted;
    item.step_scale = step_scale;
    item.evaluation = evaluation;
    item.maximum_parameter_change = proposal.maximum_parameter_change;
    item.maximum_plf_change = proposal.maximum_plf_change;
    history_.push_back(item);
}

void AutoCalibrationEngine::Accept(const CalibrationEvaluation& evaluation) {
    accepted_evaluation_ = evaluation;
}

bool AutoCalibrationEngine::WriteOutputs(
    const std::vector<CalibrationLink>& links,
    const CalibrationProposal& final_proposal,
    std::string* error_message) const {
    std::ofstream history(config_.history_output_file.c_str());
    if (!history.is_open()) {
        if (error_message != NULL) {
            *error_message = "Cannot write " + config_.history_output_file;
        }
        return false;
    }
    history << "outer_iteration,retry,accepted,step_scale,objective,episode_fit_loss,"
            << "duration_loss,trough_speed_loss,average_speed_loss,s3_loss,"
            << "volume_envelope_loss,count_loss,no_episode_loss,"
            << "prior_loss,route_policy_distance,unobserved_volume_deviation,"
            << "vmt_change_fraction,vht_change_fraction,relative_gap_pct,"
            << "guardrails_pass,max_parameter_change,max_plf_change\n";
    history << std::setprecision(12);
    for (std::size_t i = 0; i < history_.size(); ++i) {
        const CalibrationIteration& item = history_[i];
        const CalibrationEvaluation& value = item.evaluation;
        history << item.outer_iteration << ',' << item.retry << ','
                << (item.accepted ? 1 : 0) << ',' << item.step_scale << ','
                << value.objective << ',' << value.episode_fit_loss << ','
                << value.duration_loss << ','
                << value.trough_speed_loss << ',' << value.average_speed_loss
                << ',' << value.s3_loss << ',' << value.volume_envelope_loss
                << ',' << value.count_loss << ',' << value.no_episode_loss << ','
                << value.prior_loss << ',' << value.route_policy_distance << ','
                << value.unobserved_volume_deviation << ','
                << value.vmt_change_fraction << ',' << value.vht_change_fraction
                << ',' << value.relative_gap_pct << ','
                << (value.guardrails_pass ? 1 : 0) << ','
                << item.maximum_parameter_change << ','
                << item.maximum_plf_change << '\n';
    }
    history.close();

    std::ofstream audit(config_.audit_output_file.c_str());
    if (!audit.is_open()) {
        if (error_message != NULL) {
            *error_message = "Cannot write " + config_.audit_output_file;
        }
        return false;
    }
    audit << "link_id,from_node_id,to_node_id,vdf_code,corridor,facility_class,"
          << "calibration_eligible,calibration_exclusion_reason,observation_class,"
          << "baseline_pce_volume,baseline_vehicle_volume,pce_volume,vehicle_volume,baseline_doc,doc,"
          << "travel_time_min,mode1_plf,departure_plf,"
          << "final_plf,mode1_qcd,"
          << "final_qcd,mode1_qcp,final_qcp,qn,qs,alpha,beta,"
          << "observed_p_raw_hr,observed_p_effective_hr,duration_censored,"
          << "modeled_p_hr,observed_vt2_mph,modeled_vt2_mph,"
          << "observed_avg_speed_mph,modeled_avg_speed_mph,s3_volume,"
          << "cube_vehicle_volume,volume_envelope_status,volume_direction_log_residual,"
          << "volume_direction_control_gain,"
          << "observation_quality\n";
    audit << std::setprecision(12);
    const double period_hours = PeriodDurationHours(
        config_.period_start_hour, config_.period_end_hour);
    for (std::size_t i = 0; i < links.size(); ++i) {
        const CalibrationLink& link = links[i];
        const double denominator = link.lanes * period_hours * link.plf *
            link.capacity_vphpl;
        const double doc = denominator > kEpsilon ? link.volume / denominator : 0.0;
        const CalibrationLink& baseline = i < baseline_links_.size()
            ? baseline_links_[i] : link;
        const double baseline_denominator = baseline.lanes * period_hours *
            baseline.plf * baseline.capacity_vphpl;
        const double baseline_doc = baseline_denominator > kEpsilon
            ? baseline.volume / baseline_denominator : 0.0;
        const double modeled_duration = link.qcd *
            std::pow(std::max(0.0, doc), link.qn);
        const double modeled_trough = IsFinitePositive(link.cutoff_speed_mph)
            ? link.cutoff_speed_mph /
                (1.0 + link.qcp *
                 std::pow(std::max(0.0, modeled_duration), link.qs))
            : 0.0;
        const double modeled_average = ModeledAverageSpeed(link);
        VolumeEnvelopePosition envelope_position =
            VolumeEnvelopePosition::InsideOrSingle;
        const double volume_direction_residual = VolumeDirectionResidual(
            link, &envelope_position);
        const double volume_direction_control_gain =
            envelope_position == VolumeEnvelopePosition::Below
            ? std::max(0.0, config_.volume_below_envelope_control_gain)
            : envelope_position == VolumeEnvelopePosition::Above
            ? std::max(0.0, config_.volume_envelope_control_gain)
            : std::max(0.0, config_.s3_control_gain);
        const double departure_plf = i < final_proposal.departure_plf.size()
            ? final_proposal.departure_plf[i] : link.plf;
        const double effective_duration = EffectiveDurationTarget(
            link, period_hours, config_);
        const bool duration_censored = IsFinitePositive(link.observed_duration_hour) &&
            IsFinitePositive(effective_duration) &&
            effective_duration + 1e-12 < link.observed_duration_hour;
        audit << link.external_link_id << ',' << link.from_node_id << ','
              << link.to_node_id << ',' << link.vdf_code << ','
              << CSVText(link.corridor) << ','
              << CSVText(link.facility_class) << ','
              << (link.calibration_eligible ? 1 : 0) << ','
              << CSVText(link.calibration_exclusion_reason) << ','
              << ObservationClassName(link.observation_class) << ','
              << baseline.volume << ',' << baseline.vehicle_volume << ','
              << link.volume << ','
              << link.vehicle_volume << ',' << baseline_doc << ',' << doc << ','
              << link.travel_time_minutes << ','
              << link.mode1_plf << ',' << departure_plf << ',' << link.plf << ','
              << link.mode1_qcd << ',' << link.qcd << ',' << link.mode1_qcp
              << ',' << link.qcp << ',' << link.qn << ',' << link.qs << ','
              << link.alpha << ',' << link.beta << ','
              << link.observed_duration_hour << ',' << effective_duration << ','
              << (duration_censored ? 1 : 0) << ',' << modeled_duration << ','
              << link.observed_trough_speed_mph << ',' << modeled_trough << ','
              << link.observed_average_speed_mph << ',' << modeled_average << ','
              << link.s3_volume << ',' << link.cube_vehicle_volume << ','
              << VolumeEnvelopePositionName(envelope_position) << ','
              << volume_direction_residual << ','
              << volume_direction_control_gain << ','
              << link.observation_quality << '\n';
    }
    audit.close();

    std::ofstream constraint_audit(
        config_.volume_constraint_audit_output_file.c_str());
    if (!constraint_audit.is_open()) {
        if (error_message != NULL) {
            *error_message = "Cannot write " +
                config_.volume_constraint_audit_output_file;
        }
        return false;
    }
    constraint_audit << "constraint_id,constraint_type,target_vehicle_volume,"
        << "modeled_vehicle_volume,log_residual,weight,member_count\n";
    constraint_audit << std::setprecision(12);
    for (std::size_t c = 0; c < volume_constraints_.size(); ++c) {
        const CalibrationVolumeConstraint& constraint = volume_constraints_[c];
        double modeled = 0.0;
        for (std::size_t p = 0; p < constraint.link_indices.size(); ++p) {
            const int link_index = constraint.link_indices[p];
            if (link_index >= 0 && link_index < static_cast<int>(links.size())) {
                modeled += constraint.coefficients[p] *
                    std::max(0.0, links[link_index].vehicle_volume);
            }
        }
        const double residual = IsFinitePositive(
            constraint.target_vehicle_volume)
            ? SafeLog((std::max(0.0, modeled) + 1.0) /
                (constraint.target_vehicle_volume + 1.0)) : 0.0;
        constraint_audit << CSVText(constraint.constraint_id) << ','
            << CSVText(constraint.constraint_type) << ','
            << constraint.target_vehicle_volume << ',' << modeled << ','
            << residual << ',' << constraint.weight << ','
            << constraint.link_indices.size() << '\n';
    }
    constraint_audit.close();

    if (config_.write_oracle_audit) {
        std::ofstream oracle_stream(config_.oracle_audit_output_file.c_str());
        if (!oracle_stream.is_open()) {
            if (error_message != NULL) {
                *error_message = "Cannot write " + config_.oracle_audit_output_file;
            }
            return false;
        }
        oracle_stream << "link_id,from_node_id,to_node_id,corridor,facility_class,"
            << "calibration_eligible,observation_class,oracle_status,oracle_detail,"
            << "duration_semantic_review,exact_feasible,bound_count,volume_before,"
            << "volume_after,doc_before,doc_after,target_reference_speed_mph,"
            << "oracle_doc_raw,oracle_plf_raw,oracle_qcd_raw,oracle_qcp_raw,"
            << "oracle_alpha_raw,oracle_beta_raw,plf_bound_ratio,qcd_bound_ratio,"
            << "qcp_bound_ratio,oracle_plf_applied,oracle_qcd_applied,"
            << "oracle_qcp_applied,oracle_alpha_applied,oracle_beta_applied,"
            << "observed_p_raw_hr,observed_p_effective_hr,duration_censored,"
            << "oracle_p_raw,oracle_p_applied,modeled_p_posteq,"
            << "oracle_p_raw_residual,oracle_p_applied_residual,posteq_p_residual,"
            << "observed_vt2_mph,oracle_vt2_raw,oracle_vt2_applied,modeled_vt2_posteq,"
            << "oracle_vt2_raw_residual,oracle_vt2_applied_residual,posteq_vt2_residual,"
            << "observed_avg_speed_mph,oracle_avg_speed_raw,oracle_avg_speed_applied,"
            << "modeled_avg_speed_posteq,oracle_avg_raw_residual,"
            << "oracle_avg_applied_residual,posteq_avg_residual,mode1_plf,"
            << "departure_plf,final_plf,mode1_qcd,final_qcd,mode1_qcp,final_qcp\n";
        oracle_stream << std::setprecision(12);
        const FixedVolumeOracleConfig oracle_config = BuildOracleConfig(config_);
        for (std::size_t i = 0; i < links.size(); ++i) {
            const CalibrationLink& link = links[i];
            const CalibrationLink& baseline = i < baseline_links_.size()
                ? baseline_links_[i] : link;
            const double baseline_denominator = baseline.lanes * period_hours *
                baseline.plf * baseline.capacity_vphpl;
            const double baseline_doc = baseline_denominator > kEpsilon
                ? baseline.volume / baseline_denominator : 0.0;
            const double denominator = link.lanes * period_hours * link.plf *
                link.capacity_vphpl;
            const double doc = denominator > kEpsilon
                ? link.volume / denominator : 0.0;
            FixedVolumeOracleResult oracle;
            if (link.calibration_eligible &&
                link.observation_class == ObservationClass::Episode) {
                const FixedVolumeOracleInput oracle_input = BuildOracleInput(
                    link, period_hours, config_);
                oracle = SolveFixedVolumeOracle(oracle_input, oracle_config);
            }
            const double modeled_duration = link.qcd *
                std::pow(std::max(0.0, doc), link.qn);
            const double modeled_trough = IsFinitePositive(link.cutoff_speed_mph)
                ? link.cutoff_speed_mph /
                    (1.0 + link.qcp * std::pow(
                        std::max(0.0, modeled_duration), link.qs))
                : 0.0;
            const double modeled_average = ModeledAverageSpeed(link);
            const double departure_plf = i < final_proposal.departure_plf.size()
                ? final_proposal.departure_plf[i] : link.plf;
            const double effective_duration = EffectiveDurationTarget(
                link, period_hours, config_);
            const bool duration_censored =
                IsFinitePositive(link.observed_duration_hour) &&
                IsFinitePositive(effective_duration) &&
                effective_duration + 1e-12 < link.observed_duration_hour;
            oracle_stream << link.external_link_id << ',' << link.from_node_id << ','
                << link.to_node_id << ',' << CSVText(link.corridor) << ','
                << CSVText(link.facility_class) << ','
                << (link.calibration_eligible ? 1 : 0) << ','
                << ObservationClassName(link.observation_class) << ','
                << OracleStatusName(oracle.status) << ',' << CSVText(oracle.detail) << ','
                << (oracle.duration_semantic_review ? 1 : 0) << ','
                << (oracle.exact_feasible ? 1 : 0) << ',' << oracle.bound_count << ','
                << baseline.volume << ',' << link.volume << ',' << baseline_doc << ','
                << doc << ',' << oracle.target_reference_speed_mph << ','
                << oracle.doc_raw << ',' << oracle.plf_raw << ',' << oracle.qcd_raw
                << ',' << oracle.qcp_raw << ',' << oracle.alpha_raw << ','
                << oracle.beta_raw << ',' << oracle.plf_bound_ratio << ','
                << oracle.qcd_bound_ratio << ',' << oracle.qcp_bound_ratio << ','
                << oracle.plf_applied << ',' << oracle.qcd_applied << ','
                << oracle.qcp_applied << ',' << oracle.alpha_applied << ','
                << oracle.beta_applied << ',' << link.observed_duration_hour << ','
                << effective_duration << ',' << (duration_censored ? 1 : 0) << ','
                << oracle.raw_prediction.duration_hour << ','
                << oracle.applied_prediction.duration_hour << ',' << modeled_duration
                << ',' << oracle.raw_duration_residual << ','
                << oracle.applied_duration_residual << ','
                << modeled_duration - effective_duration << ','
                << link.observed_trough_speed_mph << ','
                << oracle.raw_prediction.trough_speed_mph << ','
                << oracle.applied_prediction.trough_speed_mph << ','
                << modeled_trough << ',' << oracle.raw_trough_residual << ','
                << oracle.applied_trough_residual << ','
                << modeled_trough - link.observed_trough_speed_mph << ','
                << link.observed_average_speed_mph << ','
                << oracle.raw_prediction.average_speed_mph << ','
                << oracle.applied_prediction.average_speed_mph << ','
                << modeled_average << ',' << oracle.raw_average_residual << ','
                << oracle.applied_average_residual << ','
                << modeled_average - link.observed_average_speed_mph << ','
                << link.mode1_plf << ',' << departure_plf << ',' << link.plf << ','
                << link.mode1_qcd << ',' << link.qcd << ',' << link.mode1_qcp
                << ',' << link.qcp << '\n';
        }
    }

    std::ofstream summary(config_.summary_output_file.c_str());
    if (!summary.is_open()) {
        if (error_message != NULL) {
            *error_message = "Cannot write " + config_.summary_output_file;
        }
        return false;
    }
    int episode_links = 0;
    int no_episode_links = 0;
    int unobserved_links = 0;
    for (std::size_t i = 0; i < links.size(); ++i) {
        if (links[i].observation_class == ObservationClass::Episode) {
            ++episode_links;
        } else if (links[i].observation_class == ObservationClass::NoEpisode) {
            ++no_episode_links;
        } else {
            ++unobserved_links;
        }
    }
    summary << std::setprecision(12)
            << "{\n"
            << "  \"status\": \"complete\",\n"
            << "  \"outer_attempts\": "
            << (history_.empty() ? 0 : history_.size() - 1) << ",\n"
            << "  \"final_objective\": " << accepted_evaluation_.objective << ",\n"
            << "  \"final_episode_fit_loss\": "
            << accepted_evaluation_.episode_fit_loss << ",\n"
            << "  \"final_average_speed_loss\": "
            << accepted_evaluation_.average_speed_loss << ",\n"
            << "  \"final_s3_loss\": " << accepted_evaluation_.s3_loss << ",\n"
            << "  \"final_volume_envelope_loss\": "
            << accepted_evaluation_.volume_envelope_loss << ",\n"
            << "  \"final_count_loss\": " << accepted_evaluation_.count_loss << ",\n"
            << "  \"final_relative_gap_pct\": "
            << accepted_evaluation_.relative_gap_pct << ",\n"
            << "  \"route_policy_distance\": "
            << accepted_evaluation_.route_policy_distance << ",\n"
            << "  \"guardrails_pass\": "
            << (accepted_evaluation_.guardrails_pass ? "true" : "false") << ",\n"
            << "  \"uses_origin_specific_profiles\": "
            << (uses_origin_profiles_ ? "true" : "false") << ",\n"
            << "  \"uses_mode_fallback_profiles\": "
            << (uses_mode_profiles_ ? "true" : "false") << ",\n"
            << "  \"episode_links\": " << episode_links << ",\n"
            << "  \"no_episode_links\": " << no_episode_links << ",\n"
            << "  \"unobserved_links\": " << unobserved_links << ",\n"
            << "  \"volume_constraints\": "
            << volume_constraints_.size() << "\n"
            << "}\n";
    return true;
}

}  // namespace taplite

// ---------------------------------------------------------------------------
// Production fixed-volume refinement oracle
// ---------------------------------------------------------------------------
// Kept in this translation unit so AutoCalibration.cpp/.h are the single
// native implementation surface for calibration.  The equations and bounds
// are unchanged from the validated refinement implementation.

namespace taplite {
namespace {

const double kOracleEpsilon = 1e-12;

double OracleNaN() {
    return std::numeric_limits<double>::quiet_NaN();
}

bool OraclePositive(double value) {
    return std::isfinite(value) && value > 0.0;
}

double OracleClamp(double value, double lower, double upper) {
    return std::max(lower, std::min(upper, value));
}

bool OracleChanged(double raw, double applied) {
    return std::isfinite(raw) && std::isfinite(applied) &&
        std::fabs(raw - applied) >
            1e-12 * std::max(1.0, std::fabs(raw));
}

double OracleBoundRatio(double raw, double lower, double upper) {
    if (!OraclePositive(raw)) {
        return OracleNaN();
    }
    if (raw < lower) {
        return raw / std::max(lower, kOracleEpsilon);
    }
    if (raw > upper) {
        return raw / std::max(upper, kOracleEpsilon);
    }
    return 1.0;
}

void SetOracleResiduals(
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
    : doc(OracleNaN()),
      duration_hour(OracleNaN()),
      trough_speed_mph(OracleNaN()),
      reference_speed_mph(OracleNaN()),
      queue_speed_mph(OracleNaN()),
      average_speed_mph(OracleNaN()) {}

FixedVolumeOracleInput::FixedVolumeOracleInput()
    : volume(OracleNaN()),
      lanes(OracleNaN()),
      period_hours(OracleNaN()),
      capacity_vphpl(OracleNaN()),
      free_speed_mph(OracleNaN()),
      cutoff_speed_mph(OracleNaN()),
      qn(OracleNaN()),
      qs(OracleNaN()),
      observed_duration_hour(OracleNaN()),
      observed_trough_speed_mph(OracleNaN()),
      observed_average_speed_mph(OracleNaN()) {}

FixedVolumeOracleConfig::FixedVolumeOracleConfig()
    : theta(8.0 / 15.0),
      use_bounds(false),
      minimum_plf(0.10),
      maximum_plf(1.25),
      minimum_qcd(0.01),
      maximum_qcd(20.0),
      minimum_qcp(0.001),
      maximum_qcp(20.0),
      residual_tolerance(1e-7),
      average_speed_feasibility_mph(0.10),
      saturated_speed_tolerance_mph(0.25) {}

FixedVolumeOracleResult::FixedVolumeOracleResult()
    : status(OracleStatus::NotApplicable),
      exact_feasible(false),
      duration_semantic_review(false),
      bound_count(0),
      target_reference_speed_mph(OracleNaN()),
      doc_raw(OracleNaN()),
      plf_raw(OracleNaN()),
      qcd_raw(OracleNaN()),
      qcp_raw(OracleNaN()),
      alpha_raw(OracleNaN()),
      beta_raw(OracleNaN()),
      plf_bound_ratio(OracleNaN()),
      qcd_bound_ratio(OracleNaN()),
      qcp_bound_ratio(OracleNaN()),
      plf_applied(OracleNaN()),
      qcd_applied(OracleNaN()),
      qcp_applied(OracleNaN()),
      alpha_applied(OracleNaN()),
      beta_applied(OracleNaN()),
      raw_duration_residual(OracleNaN()),
      raw_trough_residual(OracleNaN()),
      raw_average_residual(OracleNaN()),
      applied_duration_residual(OracleNaN()),
      applied_trough_residual(OracleNaN()),
      applied_average_residual(OracleNaN()) {}

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
    if (!OraclePositive(denominator) || !OraclePositive(qcd) ||
        !OraclePositive(qcp) || !OraclePositive(input.qn) ||
        !OraclePositive(input.qs) || !OraclePositive(alpha) ||
        !OraclePositive(beta) || !OraclePositive(input.free_speed_mph) ||
        !OraclePositive(input.cutoff_speed_mph)) {
        return result;
    }
    result.doc = input.volume / denominator;
    if (!std::isfinite(result.doc) || result.doc < 0.0) {
        return QVDFPrediction();
    }
    result.duration_hour = qcd * std::pow(result.doc, input.qn);
    result.trough_speed_mph = input.cutoff_speed_mph /
        std::max(kOracleEpsilon, 1.0 + qcp *
            std::pow(result.duration_hour, input.qs));
    result.reference_speed_mph = result.doc < 1.0
        ? (1.0 - result.doc) * input.free_speed_mph +
            result.doc * input.cutoff_speed_mph
        : input.cutoff_speed_mph;
    result.queue_speed_mph = result.reference_speed_mph /
        std::max(kOracleEpsilon, 1.0 + alpha *
            std::pow(result.doc, beta));
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
    if (!OraclePositive(input.volume) || !OraclePositive(input.lanes) ||
        !OraclePositive(input.period_hours) ||
        !OraclePositive(input.capacity_vphpl) ||
        !OraclePositive(input.free_speed_mph) ||
        !OraclePositive(input.cutoff_speed_mph) ||
        input.cutoff_speed_mph >= input.free_speed_mph ||
        !OraclePositive(input.qn) || !OraclePositive(input.qs) ||
        !OraclePositive(input.observed_duration_hour) ||
        !OraclePositive(input.observed_trough_speed_mph) ||
        !OraclePositive(input.observed_average_speed_mph) ||
        !OraclePositive(config.theta)) {
        result.status = OracleStatus::InvalidTarget;
        result.detail =
            "nonpositive/nonfinite input or cutoff not below free speed";
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
    if (!OraclePositive(queue_factor)) {
        result.status = OracleStatus::InvalidTarget;
        result.detail = "invalid analytical queue factor";
        return result;
    }

    const double p = input.observed_duration_hour;
    if (p <= input.period_hours) {
        const double a = p / input.period_hours;
        const double coefficient =
            a / queue_factor + (1.0 - a) / 2.0;
        const double free_coefficient = (1.0 - a) / 2.0;
        if (std::fabs(coefficient) <= kOracleEpsilon) {
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
        result.detail =
            "average speed identifies the saturated branch but not DOC";
        return result;
    }
    if (result.target_reference_speed_mph <=
            input.cutoff_speed_mph + speed_tolerance ||
        result.target_reference_speed_mph >
            input.free_speed_mph + speed_tolerance) {
        result.status = OracleStatus::AverageSpeedIncompatible;
        result.detail =
            "required reference speed is outside (cutoff, free flow]";
        return result;
    }
    result.target_reference_speed_mph = std::min(
        result.target_reference_speed_mph, input.free_speed_mph);
    result.doc_raw =
        (input.free_speed_mph - result.target_reference_speed_mph) /
        (input.free_speed_mph - input.cutoff_speed_mph);
    if (!(result.doc_raw > 0.0 && result.doc_raw < 1.0)) {
        result.status = OracleStatus::AverageSpeedIncompatible;
        result.detail =
            "inverse reference speed does not produce DOC in (0,1)";
        return result;
    }

    result.plf_raw = input.volume /
        (input.lanes * input.period_hours * input.capacity_vphpl *
         result.doc_raw);
    result.qcd_raw = p / std::pow(result.doc_raw, input.qn);
    result.qcp_raw = severity / std::pow(p, input.qs);
    result.alpha_raw = config.theta * result.qcp_raw *
        std::pow(result.qcd_raw, input.qs);
    result.beta_raw = input.qn * input.qs;
    if (!OraclePositive(result.plf_raw) ||
        !OraclePositive(result.qcd_raw) ||
        !OraclePositive(result.qcp_raw) ||
        !OraclePositive(result.alpha_raw) ||
        !OraclePositive(result.beta_raw)) {
        result.status = OracleStatus::InvalidTarget;
        result.detail = "inverse solution produced a nonpositive parameter";
        return result;
    }

    result.plf_bound_ratio = OracleBoundRatio(
        result.plf_raw, config.minimum_plf, config.maximum_plf);
    result.qcd_bound_ratio = OracleBoundRatio(
        result.qcd_raw, config.minimum_qcd, config.maximum_qcd);
    result.qcp_bound_ratio = OracleBoundRatio(
        result.qcp_raw, config.minimum_qcp, config.maximum_qcp);
    result.plf_applied = config.use_bounds
        ? OracleClamp(result.plf_raw, config.minimum_plf, config.maximum_plf)
        : result.plf_raw;
    result.qcd_applied = config.use_bounds
        ? OracleClamp(result.qcd_raw, config.minimum_qcd, config.maximum_qcd)
        : result.qcd_raw;
    result.qcp_applied = config.use_bounds
        ? OracleClamp(result.qcp_raw, config.minimum_qcp, config.maximum_qcp)
        : result.qcp_raw;
    result.alpha_applied = config.theta * result.qcp_applied *
        std::pow(result.qcd_applied, input.qs);
    result.beta_applied = result.beta_raw;
    const bool plf_limited =
        OracleChanged(result.plf_raw, result.plf_applied);
    const bool qcd_limited =
        OracleChanged(result.qcd_raw, result.qcd_applied);
    const bool qcp_limited =
        OracleChanged(result.qcp_raw, result.qcp_applied);
    result.bound_count = static_cast<int>(plf_limited) +
        static_cast<int>(qcd_limited) + static_cast<int>(qcp_limited);
    if (result.bound_count > 1) {
        result.status = OracleStatus::MultiBoundLimited;
        result.detail =
            "multiple inverse parameters are outside configured bounds";
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
        input, result.plf_applied, result.qcd_applied,
        result.qcp_applied, result.alpha_applied, result.beta_applied);
    SetOracleResiduals(input, result.raw_prediction,
        &result.raw_duration_residual, &result.raw_trough_residual,
        &result.raw_average_residual);
    SetOracleResiduals(input, result.applied_prediction,
        &result.applied_duration_residual, &result.applied_trough_residual,
        &result.applied_average_residual);
    const double tolerance = std::max(config.residual_tolerance, 0.0);
    result.exact_feasible =
        std::fabs(result.raw_duration_residual) <= tolerance &&
        std::fabs(result.raw_trough_residual) <= tolerance &&
        std::fabs(result.raw_average_residual) <=
            std::max(tolerance, config.average_speed_feasibility_mph);
    if (!result.exact_feasible &&
        result.status == OracleStatus::ExactFeasible) {
        result.status = OracleStatus::SemanticReview;
        result.detail =
            "inverse parameters did not reconstruct targets within tolerance";
    }
    return result;
}

const char* OracleStatusName(OracleStatus status) {
    switch (status) {
        case OracleStatus::ExactFeasible: return "EXACT_FEASIBLE";
        case OracleStatus::InvalidTarget: return "INVALID_TARGET";
        case OracleStatus::QcpSignIncompatible:
            return "QCP_SIGN_INCOMPATIBLE";
        case OracleStatus::AverageSpeedIncompatible:
            return "AVG_SPEED_INCOMPATIBLE";
        case OracleStatus::SaturatedNonidentifiable:
            return "SATURATED_NONIDENTIFIABLE";
        case OracleStatus::PlfBoundLimited: return "PLF_BOUND_LIMITED";
        case OracleStatus::QcdBoundLimited: return "QCD_BOUND_LIMITED";
        case OracleStatus::QcpBoundLimited: return "QCP_BOUND_LIMITED";
        case OracleStatus::MultiBoundLimited: return "MULTI_BOUND_LIMITED";
        case OracleStatus::SemanticReview: return "SEMANTIC_REVIEW";
        default: return "NOT_APPLICABLE";
    }
}

}  // namespace taplite
