#include "AutoCalibration.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <vector>

namespace {

bool Near(double left, double right, double tolerance = 1e-7) {
    return std::fabs(left - right) <= tolerance;
}

taplite::FixedVolumeOracleInput ExactOracleCase() {
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
        (2.0 / 3.0) * queue_speed +
        (1.0 / 3.0) * (reference_speed + 65.0) / 2.0;
    return input;
}

void TestFixedVolumeOracle() {
    taplite::FixedVolumeOracleConfig config;
    config.use_bounds = false;
    const taplite::FixedVolumeOracleInput input = ExactOracleCase();
    taplite::FixedVolumeOracleResult exact =
        taplite::SolveFixedVolumeOracle(input, config);
    assert(exact.status == taplite::OracleStatus::ExactFeasible);
    assert(exact.exact_feasible);
    assert(Near(exact.doc_raw, 0.5));
    assert(Near(exact.raw_prediction.duration_hour,
                input.observed_duration_hour));
    assert(Near(exact.raw_prediction.trough_speed_mph,
                input.observed_trough_speed_mph));
    assert(Near(exact.raw_prediction.average_speed_mph,
                input.observed_average_speed_mph));

    config.use_bounds = true;
    config.maximum_plf = 0.20;
    config.maximum_qcd = 3.0;
    const taplite::FixedVolumeOracleResult bounded =
        taplite::SolveFixedVolumeOracle(input, config);
    assert(bounded.status == taplite::OracleStatus::MultiBoundLimited);
    assert(bounded.bound_count == 2);

    taplite::FixedVolumeOracleInput incompatible = input;
    incompatible.observed_average_speed_mph = 5.0;
    taplite::FixedVolumeOracleResult failed =
        taplite::SolveFixedVolumeOracle(incompatible, config);
    assert(failed.status == taplite::OracleStatus::AverageSpeedIncompatible);

    incompatible = input;
    incompatible.observed_trough_speed_mph = 50.0;
    failed = taplite::SolveFixedVolumeOracle(incompatible, config);
    assert(failed.status == taplite::OracleStatus::QcpSignIncompatible);
}

}  // namespace

int main()
{
    TestFixedVolumeOracle();
    const char* profile_path = "auto_calibration_test_departure.csv";
    {
        std::ofstream profile(profile_path);
        profile << "origin_zone_id,mode,time_min,share\n"
                << "1,sov,360,0.10\n"
                << "1,sov,375,0.20\n"
                << "1,sov,390,0.50\n"
                << "1,sov,405,0.20\n";
    }

    taplite::AutoCalibrationConfig config;
    assert(config.calibration_fit_mode == "refined_fixed_point");
    // The proposal checks below exercise the historical direct-update branch;
    // refinement/oracle behavior has its own focused test executable.
    config.calibration_fit_mode = "legacy_regularized";
    config.enabled = true;
    config.period_start_hour = 6.0;
    config.period_end_hour = 9.0;
    config.departure_profile_file = profile_path;
    config.workers = 2;
    config.inner_gap_tolerance_pct = 0.01;
    taplite::AutoCalibrationEngine engine(config);
    std::string error;
    assert(engine.LoadDepartureProfiles(&error));
    assert(engine.uses_origin_specific_profiles());

    std::vector<taplite::CalibrationLink> links(2);
    for (int i = 0; i < 2; ++i)
    {
        links[i].internal_index = i;
        links[i].external_link_id = 100 + i;
        links[i].from_node_id = i + 1;
        links[i].to_node_id = i + 2;
        links[i].vdf_code = 101;
        links[i].length_miles = 1.0;
        links[i].lanes = 2.0;
        links[i].capacity_vphpl = 1000.0;
        links[i].free_speed_mph = 60.0;
        links[i].cutoff_speed_mph = 42.0;
        links[i].volume = 1600.0;
        links[i].vehicle_volume = 1600.0;
        links[i].travel_time_minutes = 1.2;
        links[i].plf = links[i].mode1_plf = 1.0;
        links[i].qcd = links[i].mode1_qcd = 1.0;
        links[i].qcp = links[i].mode1_qcp = 0.25;
        links[i].qn = 1.0;
        links[i].qs = 4.0;
        links[i].observed_average_speed_mph = 50.0;
        links[i].s3_volume = 1500.0;
        links[i].observation_quality = 0.9;
    }
    links[0].observation_class = taplite::ObservationClass::Episode;
    links[0].observed_duration_hour = 0.8;
    links[0].observed_trough_speed_mph = 25.0;
    links[1].observation_class = taplite::ObservationClass::NoEpisode;

    taplite::CalibrationRoute route;
    route.mode_index = 1;
    route.origin_zone_id = 1;
    route.destination_zone_id = 2;
    route.mode_name = "sov";
    route.od_demand = 1600.0;
    route.pce = 1.0;
    route.share = 1.0;
    route.link_indices.push_back(0);
    route.link_indices.push_back(1);
    route.link_travel_times_minutes.push_back(1.2);
    route.link_travel_times_minutes.push_back(1.2);
    std::vector<taplite::CalibrationRoute> routes(1, route);

    engine.Initialize(links, routes, 0.001);
    assert(engine.history().size() == 1);
    assert(engine.history()[0].outer_iteration == 0);
    assert(engine.history()[0].accepted);
    taplite::CalibrationEvaluation worse = engine.accepted_evaluation();
    worse.objective = engine.accepted_evaluation().objective +
        std::max(1e-9, std::fabs(engine.accepted_evaluation().objective)) *
            (config.objective_relative_tolerance + 0.01);
    worse.guardrails_pass = true;
    assert(!engine.ShouldAccept(worse));
    taplite::CalibrationProposal proposal = engine.Propose(links, routes, 1.0);
    assert(proposal.plf.size() == links.size());
    assert(proposal.plf[0] >= config.minimum_plf);
    assert(proposal.plf[0] <= config.maximum_plf);
    assert(proposal.qcd[0] != links[0].qcd);
    assert(proposal.qcd[1] == links[1].qcd);
    assert(std::fabs(proposal.beta[0] - 4.0) < 1e-12);
    assert(std::fabs(
        proposal.alpha[0] - config.theta * proposal.qcp[0] *
            std::pow(proposal.qcd[0], links[0].qs)) < 1e-12);

    // A legacy value below the calibration bound must approach the bound
    // continuously as the retry step shrinks, rather than snapping to it.
    links[0].qcd = links[0].mode1_qcd = 0.00001;
    config.minimum_qcd = 0.01;
    taplite::AutoCalibrationEngine bounded_engine(config);
    assert(bounded_engine.LoadDepartureProfiles(&error));
    bounded_engine.Initialize(links, routes, 0.001);
    const taplite::CalibrationProposal bounded_proposal =
        bounded_engine.Propose(links, routes, 0.25);
    assert(bounded_proposal.qcd[0] > links[0].qcd);
    assert(bounded_proposal.qcd[0] < config.minimum_qcd);
    assert(bounded_proposal.maximum_parameter_change < 10.0);

    std::remove(profile_path);
    return 0;
}
