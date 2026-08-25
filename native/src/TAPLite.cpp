// TAPLite.cpp : This file contains the 'main' function. Program execution begins and ends here.

// This code is built based on the implementation available at:
// http://www.bgu.ac.il/~bargera/tntp/FW.zip
// tntp from Dr. Hillel Bar-Gera. It includes one of the most efficient Deque implementations of 
// the modified label correcting (MLC) algorithm in C. For more details, see mincostroutes.cpp. 
// The enhanced C++ counterpart has served as the path engine for Path4GMNS and TransOMS.
// https://github.com/jdlph/TAP101
#define M_PI 3.14159265358979323846
#define FLOAT_ACCURACY 1.0E-15
#define NO_COSTPARAMETERS 4
#define IVTT 0
#define OVTT 1
#define MONETARY 2
#define DIST 3
#define INVALID -1
#define MAX_MODE_TYPES  10

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define SQR(x) ((x) * (x))
#define FABS(x) ((x) >= 0 ? (x) : (-x))
#define VALID(x) ((x) != -1)

#include "TAPLite.h"
#include "AutoCalibration.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <memory>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <cstring>
#include <algorithm>
#include <limits>
#include <new>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <chrono>  // for high_resolution_clock
#include <deque>
#include <queue>
#include <utility>
#include <iomanip> // For std::setw and std::setfill
#include <unordered_map>  // For hash table (unordered_map)

using namespace std;

typedef double cost_vector[NO_COSTPARAMETERS];

#ifndef _WIN32
void fopen_s(FILE** file, const char* fileName, const char* mode)
{
   *file = fopen(fileName, mode);
}
#endif

struct link_record {
	int internal_from_node_id;
	int internal_to_node_id;
	int link_id;            // internal sequential link id (row index k)
	int external_link_id;   // original link_id as read from link.csv
	int link_type;
	int external_from_node_id;
	int external_to_node_id;

	double Lane_Capacity;
	double Lane_SaturationFlowRate; 
	double Link_Capacity;
	double lanes;
	double FreeTravelTime;
	double free_speed;
	double Cutoff_Speed;


	double VDF_Alpha;
	double VDF_A;        // ARC modified-BPR linear term: fftt*(1 + A*(V/C) + alpha*(V/C)^beta); 0 = standard BPR
	double green_ratio;  // SANDAG signal green ratio g/C (vdf_type=6 intersection delay); default 0.45
	double BoverC;
	double VDF_Beta;
	double VDF_plf;
	int VDF_type;        // 0 = BPR (default), 1 = Spiess conic
	double Conic_a;      // conic alpha (per functional type)
	double Conic_b;      // conic beta  (= (2a-1)/(2a-2) per Spiess)
	double Q_cd, Q_n, Q_cp, Q_s;
	// Mode 2 auto-calibration observations and immutable Mode 1 priors. These
	// fields are populated once from link.csv and then consumed in memory by
	// AutoCalibration.cpp; no outer iteration re-reads the converted network.
	std::string Calibration_observation_class;
	std::string Calibration_corridor;
	std::string Calibration_facility_class;
	std::string Calibration_target_tmc;
	std::string Calibration_exclusion_reason;
	bool Calibration_eligible;
	double Calibration_observed_P;
	double Calibration_observed_vt2;
	double Calibration_observed_average_speed;
	double Calibration_S3_volume;
	double Calibration_cube_vehicle_volume;
	double Calibration_observation_quality;
	double Mode1_plf, Mode1_Q_cd, Mode1_Q_cp;
	// Optional reporting-profile selector from link.csv. -1 preserves the
	// historical link_type/observed-t2 activation rule.
	int QVDF_profile_mode;
	double QVDF_t0_observed; // observed episode start, decimal hour; NaN = unavailable
	double QVDF_t2;     // observed queue-profile trough time, decimal hour; NaN = midpoint fallback
	double QVDF_t3_observed; // observed episode end, decimal hour; NaN = unavailable
	double QVDF_start_speed_observed; // observed first profile sample, mph; NaN = modeled fallback
	double QVDF_end_speed_observed; // observed last profile sample, mph; NaN = modeled fallback

	double length;
	double Speed;
	std::string allowed_uses;
	int mode_allowed_use[MAX_MODE_TYPES];
	double mode_MainVolume[MAX_MODE_TYPES];
	double mode_Base_demand_volume[MAX_MODE_TYPES];
	double mode_SubVolume[MAX_MODE_TYPES];
	double mode_SDVolume[MAX_MODE_TYPES];

	double mode_Toll[MAX_MODE_TYPES];
	double mode_AdditionalCost[MAX_MODE_TYPES];
	int non_uturn_flag; 
	double Travel_time;  // final travel time used in assignment 
	double BPR_TT;  // BPR_based travel time  for the entire assignment period  
	double QVDF_TT;  // QVDF_based travel time for the entire assignment period  

	bool b_withmovement_restrictions;
	double GenCost;
	double GenCostDer;
	double Ref_volume;
	double Base_demand_volume;
	double background_volume; 
	double Obs_volume[MAX_MODE_TYPES];
	double reassignment_volume;
	std::string geometry;

	int timing_arc_flag, cycle_length, start_green_time, end_green_time;


	link_record()
	{
		// Zero-initialize all scalar members first. Previously these records were
		// allocated with calloc (zero-filled); now they are allocated with new[],
		// which runs this constructor but leaves any member we do not assign
		// indeterminate. Reproduce the old zero default here, then override the
		// fields with non-zero defaults below.
		internal_from_node_id = internal_to_node_id = 0;
		link_id = external_link_id = 0;
		external_from_node_id = external_to_node_id = 0;
		Lane_Capacity = Link_Capacity = lanes = 0;
		FreeTravelTime = free_speed = Cutoff_Speed = 0;
		BoverC = 0;
		length = Speed = 0;
		GenCost = GenCostDer = 0;
		reassignment_volume = 0;
		for (int m = 0; m < MAX_MODE_TYPES; m++)
		{
			mode_allowed_use[m] = 0;
			mode_MainVolume[m] = 0;
			mode_Base_demand_volume[m] = 0;
			mode_SubVolume[m] = 0;
			mode_SDVolume[m] = 0;
			mode_Toll[m] = 0;
			mode_AdditionalCost[m] = 0;
			Obs_volume[m] = 0;
		}

		link_type = 1;

		VDF_Alpha = 0.15;
		VDF_A = 0;
		green_ratio = 0.45;
		VDF_Beta = 4;
		VDF_plf = 1;
		VDF_type = 0;
		Conic_a = 0.0;   // 0 => fall back to VDF_Alpha as conic a (staged convention)
		Conic_b = 0.0;   // 0 => fall back to VDF_Beta  as conic b

		Q_cd = 1.0;
		Q_n = 1.0;
		Q_cp = 0.28125 /*0.15*15/8*/;
		Q_s = 4;
		Calibration_observation_class = "U";
		Calibration_corridor = "";
		Calibration_facility_class = "";
		Calibration_target_tmc = "";
		Calibration_exclusion_reason = "";
		Calibration_eligible = true;
		Calibration_observed_P = std::numeric_limits<double>::quiet_NaN();
		Calibration_observed_vt2 = std::numeric_limits<double>::quiet_NaN();
		Calibration_observed_average_speed = std::numeric_limits<double>::quiet_NaN();
		Calibration_S3_volume = std::numeric_limits<double>::quiet_NaN();
		Calibration_cube_vehicle_volume = std::numeric_limits<double>::quiet_NaN();
		Calibration_observation_quality = 0.5;
		Mode1_plf = Mode1_Q_cd = Mode1_Q_cp = std::numeric_limits<double>::quiet_NaN();
		QVDF_profile_mode = -1;
		QVDF_t0_observed = std::numeric_limits<double>::quiet_NaN();
		QVDF_t2 = std::numeric_limits<double>::quiet_NaN();
		QVDF_t3_observed = std::numeric_limits<double>::quiet_NaN();
		QVDF_start_speed_observed = std::numeric_limits<double>::quiet_NaN();
		QVDF_end_speed_observed = std::numeric_limits<double>::quiet_NaN();

		Travel_time = 0;
		BPR_TT = 0;
		QVDF_TT = 0;
		Ref_volume = 0;
		Base_demand_volume = 0;
		for(int m = 1; m< MAX_MODE_TYPES-1; m++)
		{
		Obs_volume[m] = -1;
		}
		Lane_SaturationFlowRate = 1800;

		timing_arc_flag = 0;
		cycle_length = 60;
		start_green_time = 0;
		end_green_time = 30; 
		background_volume = 0; 
		non_uturn_flag = 0;
		b_withmovement_restrictions = false;

	}
	void setup(int num_of_modes)
	{
		link_type = 1;
		VDF_Alpha = 0.15;
		VDF_A = 0;
		green_ratio = 0.45;
		VDF_Beta = 4;
		VDF_plf = 1;
		VDF_type = 0;
		Conic_a = 0.0;   // 0 => fall back to VDF_Alpha as conic a (staged convention)
		Conic_b = 0.0;   // 0 => fall back to VDF_Beta  as conic b

		Q_cd = 1.0;
		Q_n = 1.0;
		Q_cp = 0.28125 /*0.15*15/8*/;
		Q_s = 4;
		Calibration_observation_class = "U";
		Calibration_corridor = "";
		Calibration_facility_class = "";
		Calibration_target_tmc = "";
		Calibration_exclusion_reason = "";
		Calibration_eligible = true;
		Calibration_observed_P = std::numeric_limits<double>::quiet_NaN();
		Calibration_observed_vt2 = std::numeric_limits<double>::quiet_NaN();
		Calibration_observed_average_speed = std::numeric_limits<double>::quiet_NaN();
		Calibration_S3_volume = std::numeric_limits<double>::quiet_NaN();
		Calibration_cube_vehicle_volume = std::numeric_limits<double>::quiet_NaN();
		Calibration_observation_quality = 0.5;
		Mode1_plf = Mode1_Q_cd = Mode1_Q_cp = std::numeric_limits<double>::quiet_NaN();
		QVDF_profile_mode = -1;
		QVDF_t0_observed = std::numeric_limits<double>::quiet_NaN();
		QVDF_t2 = std::numeric_limits<double>::quiet_NaN();
		QVDF_t3_observed = std::numeric_limits<double>::quiet_NaN();
		QVDF_start_speed_observed = std::numeric_limits<double>::quiet_NaN();
		QVDF_end_speed_observed = std::numeric_limits<double>::quiet_NaN();

		Travel_time = 0;
		BPR_TT = 0;
		QVDF_TT = 0;
		Ref_volume = 0;
		Base_demand_volume = 0;

		for (int m = 1; m <= num_of_modes; m++)
		{
			mode_Base_demand_volume[m] = 0;
			Obs_volume[m] = -1;
		}
	}
};


// Define a 2D map: ib_link_id → ob_link_id → restricted (bool)
std::map<int, std::map<int, bool>> global_movement_restrictions;

// Function to insert a movement restriction
void InsertMovementRestriction(int ib_link_id, int ob_link_id, bool is_restricted)
{
	global_movement_restrictions[ib_link_id][ob_link_id] = is_restricted;
}

// Function to check if a movement restriction exists
bool IsMovementRestricted(int ib_link_id, int ob_link_id)
{
	if (global_movement_restrictions.find(ib_link_id) != global_movement_restrictions.end())
	{
		if (global_movement_restrictions[ib_link_id].find(ob_link_id) != global_movement_restrictions[ib_link_id].end())
		{
			return global_movement_restrictions[ib_link_id][ob_link_id];
		}
	}
	return false; // Default: no restriction
}
struct mode_type {
	std::string mode_type ;
	float vot;
	float op_cost = 0;   // ARC generalized cost: $/mile operating cost (0 = off)
	float pce;
	float occ;
	int dedicated_shortest_path;
	std::string demand_file;

		
};


mode_type g_mode_type_vector[MAX_MODE_TYPES];
int g_metric_system_flag = 1; 
void StatusMessage(const char* group, const char* format, ...);

struct link_record* Link;
// Sparse agency link ids are translated once at the input boundary. All engine
// arrays use dense 1..number_of_links indices; outputs use Link[k].external_link_id.
static std::unordered_map<int, int> g_map_external_link_id_2_link_seq_no;

static inline int LinkSeqFromExternal(int external_link_id)
{
	std::unordered_map<int, int>::const_iterator it =
		g_map_external_link_id_2_link_seq_no.find(external_link_id);
	return (it == g_map_external_link_id_2_link_seq_no.end()) ? 0 : it->second;
}

int* FirstLinkFrom;  // used in shortet path algorithm
int* LastLinkFrom;

// Turn restrictions (movement.csv + non_uturn_flag). When none are present the
// kernel uses the fast node-based Minpath() unchanged; when any exist it switches
// to the exact link-state Minpath_TR(). g_MinPathPredLink_link[m][Orig][k] is the
// link traversed immediately before link k on the shortest path (link-state
// predecessor), used only by the turn-restricted back-trace.
int g_has_movement_restrictions = 0;
int*** g_MinPathPredLink_link = NULL;


double Link_Travel_Time(int k, double* Volume);
double Link_GenCost(int k, double* Volume);
double LinkCost_Integral(int k, double* Volume);
double Link_GenCostDer(int k, double* Volume);

void* Alloc_1D(int dim1, size_t size);
void** Alloc_2D(int dim1, int dim2, size_t size);
void Free_2D(void** Array, int dim1, int dim2);
void*** Alloc_3D(int dim1, int dim2, int dim3, size_t size);
void Free_3D(void*** Array, int dim1, int dim2, int dim3);

/* LINK VECTOR FUNCTIONS */
void ClearVolume(double* VolumeArray);
void VolumeDifference(double* Volume1, double* Volume2, double* Difference);
void UpdateVolume(double* MainVolume, double* SubVolume, double Lambda);

/* LINK COST FUNCTIONS */
void UpdateLinkAdditionalCost(void);
double UpdateLinkCost(double* Volume);
void UpdateLinkCostDer(double* Volume);
void GetLinkTravelTimes(double* Volume, double* TravelTime);
double TotalLinkCost(double* Volume);

/*  LINK OBJECTIVE FUNCTION */
double OF_Links(double* MainVolume);
double OF_LinksDirectionalDerivative(double* MainVolume, double* SDVolume, double Lambda);

void InitLinks(void);
void CloseLinks(void);

int Read_ODflow(double* TotalODflow, int* number_of_modes, int* no_zones);

int Minpath(int m, int Orig, int* PredLink, double* Cost_to);
int Minpath_TR(int mode, int Orig, int* PredLink, int* PredLinkOfLink, double* CostTo);
double FindMinCostRoutes(int** MinPathPredLink);
void All_or_Nothing_Assign(double** ODflow, int** MinPathPredLink, double* Volume);


double LinksSDLineSearch(double* MainVolume, double* SDVolume);
void ApplyConjugateFW(double* MainVolume, double* SubVolume, double* SDVolume, int iteration_no);


/* Gloabal variables */

int no_zones, number_of_modes, no_nodes, number_of_links, FirstThruNode;
int number_of_internal_zones = 1; 
int TotalAssignIterations = 20;
double demand_period_starting_hours = 7;
double	demand_period_ending_hours = 8;

static double DemandPeriodMidpointHours()
{
	return (demand_period_starting_hours + demand_period_ending_hours) / 2.0;
}

static bool IsHourWithinDemandPeriod(double hour)
{
	if (!std::isfinite(hour) || hour < 0.0 || hour > 24.0)
		return false;
	if (demand_period_ending_hours >= demand_period_starting_hours)
		return hour >= demand_period_starting_hours && hour <= demand_period_ending_hours;
	// Overnight periods wrap across midnight.
	return hour >= demand_period_starting_hours || hour <= demand_period_ending_hours;
}

int first_through_node_id_input = -1; 
int g_tap_log_file = 0;
int g_base_demand_mode = 1;
// MAG-style volume delay function adds an additive per-length delay to the BPR
// congested time: travel_time = T_c + g_added_delay_per_mile * length(miles).
// 0 = off (pure BPR/conic/QVDF). Set via settings.csv 'added_delay_per_mile'
// (MAG uses 1.4 min/mile on all facility types: "1.4*L + Tc").
double g_added_delay_per_mile = 0.0;
// SCAG piecewise BPR (VDF_type==7): the BPR exponent below capacity (x/C <= 1).
// SCAG's Validation Report fixes it at 4.0 globally; the calibrated per-link
// VDF_Beta (5/6/8) applies only above capacity. Constant, not a setting.
const double SCAG_UNCONGESTED_BETA = 4.0;
// Relative-gap stopping criterion: stop Frank-Wolfe once the gap (%) falls below
// this. 0 = off (run the full number_of_iterations). Set via settings.csv
// 'convergence_gap_pct'.
double g_convergence_gap_pct = 0.0;
// link_performance output detail (settings 'link_output'): 2 = full CSV (default:
// per-link QVDF queue model + per-5-min speed profile + all columns), 1 = COMPACT
// CSV (link_id/from/to/volume/doc/travel_time/speed/VMT), 3 = BINARY (same compact
// fields, packed .bin -- fastest), 0 = none.
int g_link_performance_output = 2;
// Skip the OD-accessibility post-processing (OutputODPerformance +
// GenerateAggregatedPerformanceAndAccessibility) -- the dominant cost on large
// networks. settings 'accessibility_output' (default 1 = on; 0 = fast verify runs).
int g_accessibility_output = 1;
// In FULL output mode, run the expensive per-link Link_QueueVDF speed profile ONLY
// for otherwise-eligible links carrying positive volume at or above this threshold.
// Skipped links retain a flat numeric compatibility profile. settings
// 'qvdf_volume_threshold' (default 0 = all positive-volume eligible links).
double g_qvdf_volume_threshold = 0.0;
// Shortest-path algorithm (settings 'sp_algorithm'): 1 = binary-heap Dijkstra
// label-setting (DEFAULT -- equal-or-faster everywhere, ~20-50x faster on large
// networks; correct since all VDFs give non-negative link costs); 0 = D'Esopo-Pape
// deque label-correcting (legacy). Same FirstThruNode/mode/turn rules either way.
int g_sp_algorithm = 1;
// Number of CONSECUTIVE iterations the gap must stay below convergence_gap_pct
// before stopping (ARC: gap < 1e-4 for 3 successive iters). Default 1 = stop on
// the first iteration below the gap (previous behavior).
int g_convergence_consecutive = 1;
// Relative-gap denominator: 0 (default) = normalize by the all-or-nothing total
// (system_least, legacy TAPLite); 1 = normalize by the CURRENT total (system_wide),
// the AequilibraE/agency-standard relative gap so a 1e-4 target means the same
// across agencies. Set via settings.csv 'relative_gap_standard'.
int g_relative_gap_standard = 0;
// (Bi)conjugate Frank-Wolfe (Mitradjieva-Lindberg 2013): 0 = plain FW (default),
// 1 = conjugate FW (CFW, one prior direction), 2 = bi-conjugate FW (BFW, two prior
// directions). Faster gap closure to the same UE. Set via settings 'assignment_method'.
int g_assignment_method = 0;
double* g_cfw_s1 = NULL;   // previous conjugate auxiliary point (volume) s_{k-1}
double* g_cfw_s2 = NULL;   // s_{k-2}
double* g_cfw_aux = NULL;  // scratch: this iteration's conjugate auxiliary
double* g_cfw_H = NULL;    // diagonal Hessian dCost/dV
int g_cfw_nstored = 0;     // how many prior auxiliaries are valid (0,1,2)
// Demand input format: 0 = CSV (default), 1 = binary (.bin produced by
// `python -m dtalite_qa demand-bin`). Binary skips text parsing -> much faster
// OD load on large regional models. Set via settings.csv 'demand_format'.
int g_demand_format = 0;
int g_ODME_mode = 0;
int g_ODME_target_od = -1;
double g_ODME_obs_VMT = -1;
double g_System_VMT = 0;

double g_ODME_link_volume_penalty = 0.1;  // relative weight on volume , convert the deviation of link volume to the travel time 
double g_ODME_VMT_penalty = 1e-6;  // relative weight on VMT , convert the deviation of VMT to the travel time 
double g_ODME_od_penalty = 0.01;  // relative weight on od , 
FILE* summary_log_file;


double** ODflow, TotalODflow;
double* TotalOFlow;
int* zone_outbound_link_size;

double*** seed_MDODflow;
double*** old_MDODflow;
double*** candidate_MDODflow;
double*** gradient_MDODflow;


double*** MDODflow;
double*** MDDiffODflow;  // D^c - D^b
double*** targetMDODflow; // for ODME

double*** MDRouteCost;

// ---- P0 efficiency trio (see github_taplite/docs/EFFICIENCY_STUDY.md) ----
// L1 warm start: settings.csv 'warm_start_times' = path to a prior run's
// link_performance CSV (compact or full) or link_performance.bin (DTLP,
// link_output=3). Congested Travel_time is preloaded by EXTERNAL link id right
// before the FIRST FindMinCostRoutes, so iteration-0 all-or-nothing routes on
// congested instead of free-flow times. Only the STARTING POINT changes --
// every later iteration recomputes costs from volumes normally, so the
// equilibrium is unchanged. Default empty = off (bit-identical cold start).
std::string g_warm_start_times_file;
// ---- P1/P2 increments (see EFFICIENCY_STUDY.md P1/P2) ----
// L2 warm start from FLOWS (settings 'warm_start_flows', default "" = off):
// restore MainVolume + per-mode volumes from a DTLR link_flows.bin snapshot of
// a prior run on the SAME demand and skip the iteration-0 all-or-nothing, so
// Frank-Wolfe continues from the prior run's point. Guarded by a 16-byte
// demand fingerprint stored in the snapshot: on mismatch the kernel warns
// LOUDLY and demotes to L1 behavior (the snapshot flows only seed congested
// TIMES, then a normal cold AoN start) -- never a silent restart from
// incompatible flows.
std::string g_warm_start_flows_file;
// settings 'flow_snapshot' (default 0): 1 = write link_flows.bin (DTLR) at the
// end of the run. Also written automatically when link_output=3.
int g_flow_snapshot = 0;
// settings 'column_output' (default 0) -- LEVELED:
//   0 = off (default; bit-identical to the pre-feature kernel).
//   1 = DTAC v1: the LAST-iteration shortest path per (mode,o,d) only --
//       route_assignment parity, light (no in-memory pool, no per-iteration
//       accumulation cost).
//   2 = DTAC v2: the FULL theta-share column set (the L3 store, RECOMMENDED
//       for the warm-start/rerun recipe). The in-memory pool accumulates the
//       path SET across Frank-Wolfe iterations -- each iteration's
//       all-or-nothing paths enter with the FW step weight (the same lambda
//       cascade computeTheta() applies to route_assignment.csv), duplicate
//       paths per (mode,o,d) merged by summing shares; per-path float theta
//       (shares per OD sum to 1) + demand fingerprint on disk. A PARALLEL
//       mechanism to the 5D linkIndices route store (route_output), untouched.
//   3+ = reserved (future per-iteration/full trace); warns and behaves as 2.
int g_column_output = 0;
// Keep column_output as the in-memory policy level while allowing large
// calibration runs to suppress the optional final DTAC file.  Inner and outer
// calibration logic still sees level 2; only the end-of-run serialization is
// skipped when column_file_output=0.
int g_column_file_output = 1;
// Mode 2 equilibrium-coupled auto calibration. The public auto-calibration API
// can force this on without requiring callers to mutate settings.csv. The
// implementation still reads network/demand/profile inputs only during the
// normal AssignmentAPI initialization and emits calibration files only once.
int g_auto_calibration = 0;
static bool g_force_auto_calibration = false;
// The explicit one-shot API runs in a disposable worker process.  It leaves
// the enormous nested DTAC-v2 pool for the OS to reclaim in bulk after all
// result files are closed, while the ordinary in-process API still performs
// full reusable cleanup.
static bool g_auto_calibration_one_shot = false;
std::string g_auto_calibration_config_file = "auto_calibration_settings.csv";
std::string g_departure_profile_file;
// Regional DTAC-v2 columns are bounded only for auto calibration. Assignment
// volumes still use every OD row. Retained high-volume multipath columns
// supply the arrival-time shape, which is rescaled to exact per-mode link
// volume during route summarization; omitted tiny ODs never remove volume.
static double g_calibration_minimum_column_od_volume = 0.0;
static int g_calibration_maximum_column_paths_per_od =
	std::numeric_limits<int>::max();
static double g_calibration_minimum_column_path_share = 0.0;
// L3 warm start from COLUMNS (settings 'warm_start_columns', default "" = off):
// load a DTAC column store, scatter theta x CURRENT OD volume to
// MainVolume/mode_MainVolume (demand may DIFFER from the stored fingerprint --
// that is the point: theta over a path set is a demand-invariant routing
// policy), drop paths crossing deleted/now-banned links (theta renormalized,
// counted + warned), fall back to a one-shot all-or-nothing on warm times for
// OD pairs with no stored column, then continue normal Frank-Wolfe from
// iteration 1 (self-healing). With number_of_iterations=0 this is pure
// freeze/replay. A fresh-shortest-path TRUE relative gap is ALWAYS reported so
// replay output is never mistaken for equilibrium.
std::string g_warm_start_columns_file;
// settings 'column_adjust_sweeps' (default 0): after loading columns, N
// fixed-policy gradient-projection sweeps that shift share between the STORED
// paths of each (mode, origin, destination) toward the cheaper stored path
// under current costs (costs recomputed between sweeps), WITHOUT any new
// shortest-path call. The per-sweep gap is the RESTRICTED gap (stored paths
// only); the mandatory fresh-SP TRUE gap is printed afterwards.
int g_column_adjust_sweeps = 0;
// SP-equivalence auto-detect (always on, exact): g_rep_mode[m] is the
// representative mode whose shortest-path tree mode m consumes. Two modes are
// SP-equivalent iff their mode_allowed_use masks AND per-link
// mode_AdditionalCost vectors are identical on EVERY link (PCE/occupancy never
// enter routing, so they never force a dedicated tree). rep==m computes its
// own tree; members reuse it -- Minpath on identical inputs is deterministic,
// so results are bit-identical, redundant solves are just skipped. A value of
// 0 (never grouped) is treated as "own tree" for safety.
int g_rep_mode[MAX_MODE_TYPES] = { 0 };
int g_sp_scratch_row[MAX_MODE_TYPES] = { 0 };  // CostTo scratch row per mode within a processor block
int g_sp_scratch_rows = 1;                     // rows per processor: 1 shared + 1 per shared-tree representative
// Zero-demand-origin skip (per mode): Minpath is not called for a (mode,
// origin) whose SP group carries no demand out of the origin, and the
// MDRouteCost row is BIGM-filled directly (exactly what the full destination
// loop would have produced with zero demand).
char** g_mode_origin_active = NULL;   // [mode][orig] 1 iff any MDODflow[m][O][D] > 0.000001
char** g_tree_origin_active = NULL;   // [rep mode][orig] 1 iff any member mode of the SP group is active

/* Local Declarations */
/* void FW(void); Should there be a function for fw, or should it be included in main? */
static void Init(int input_number_of_modes, int input_no_zones);
static void Close();
static void InitODflow(int input_number_of_modes, int input_no_zones);
static void CloseODflow(void);
/* End of local declarations. */

FILE* logfile;
int shortest_path_log_flag = 0;
int vehicle_log_flag = 0;
int baseODDemand_loaded_flag = 0;
int baselinkvolume_loaded_flag = 0;

FILE* link_performance_file;

#define INVALID -1      /* Represents an invalid value. */
#define BIGM 9999999      /* Represents an invalid value. */
#define WAS_IN_QUEUE -7 /* Shows that the node was in the queue before. (7 is for luck.) */

bool rewriteFile(const std::string& filename, const std::string& content) {
	// Optional: Check if content is empty and inform the user
	if (content.empty()) {
		std::cout << "The provided content is empty. The file will be cleared." << std::endl;
	}

	std::ofstream outputFile(filename);  // Open the file (truncates by default)

	if (!outputFile.is_open()) {
		std::cerr << "Error opening file: " << filename << std::endl;
		return false;
	}

	outputFile << content;  // Write the (empty) content to the file
	outputFile.close();     // Close the file
	return true;
}



double Link_QueueVDF(int k, double Volume, double& IncomingDemand, double& DOC, double& P, double& t0, double& t2, double& t3, double& vt2, double& Q_mu, double& Q_gamma, double& congestion_ref_speed,
	double& avg_queue_speed, double& avg_QVDF_period_speed, double& Severe_Congestion_P, double model_speed[300]);

class Node {
public:
	int node_id;
	int zone_id; 
	int internal_zone_no; 
	double x, y; 
	vector<int> m_incoming_link_seq_no_vector;
	vector<int> m_outgoing_link_seq_no_vector;

	Node() : node_id(-1) { zone_id = -1; internal_zone_no = -1; }
};

vector<Node> g_node_vector;

// Pre-allocated per-thread SSSP scratch, reused across ALL Minpath calls so no
// malloc/free or fresh allocation happens per shortest-path computation. Sized once
// (InitSPScratch) after the network is read; indexed by omp_get_thread_num().
static std::vector<std::vector<int>> g_sp_QueueNext;                    // Pape deque
static std::vector<std::vector<int>> g_sp_PrevLink;                     // arriving link
static std::vector<std::vector<std::pair<double, int>>> g_sp_heap;      // Dijkstra heap
static int g_sp_nnodes = 0;

void InitSPScratch()
{
	// scratch is indexed by omp_get_thread_num(), so size by the OpenMP thread count
	// (>= any concurrent thread id, independent of g_number_of_processors).
	int nthreads = 1;
#ifdef _OPENMP
	nthreads = omp_get_max_threads();
#endif
	if (nthreads < 1)
		nthreads = 1;
	g_sp_QueueNext.assign(nthreads, std::vector<int>(no_nodes + 1, INVALID));
	g_sp_PrevLink.assign(nthreads, std::vector<int>(no_nodes + 1, INVALID));
	g_sp_heap.assign(nthreads, std::vector<std::pair<double, int>>());
	for (int t = 0; t < nthreads; t++)
		g_sp_heap[t].reserve(2048);
	g_sp_nnodes = no_nodes;
}

static int GetSPScratchThreadIndex()
{
	int tid = 0;
#ifdef _OPENMP
	tid = omp_get_thread_num();
#endif
	if (tid < 0 || static_cast<std::size_t>(tid) >= g_sp_QueueNext.size() ||
		static_cast<std::size_t>(tid) >= g_sp_PrevLink.size() ||
		static_cast<std::size_t>(tid) >= g_sp_heap.size())
	{
		fprintf(stderr,
			"ERROR: OpenMP thread %d is outside shortest-path scratch range [0, %zu). "
			"Configure OpenMP before InitSPScratch().\n",
			tid, g_sp_QueueNext.size());
		if (summary_log_file != NULL)
		{
			fprintf(summary_log_file,
				"ERROR: OpenMP thread %d is outside shortest-path scratch range [0, %zu).\n",
				tid, g_sp_QueueNext.size());
			fflush(summary_log_file);
		}
		fflush(stderr);
		abort();
	}
	return tid;
}

// Min-heap comparator for the reusable Dijkstra heap (smallest cost on top).
static inline bool g_sp_heap_gt(const std::pair<double, int>& a, const std::pair<double, int>& b)
{
	return a.first > b.first;
}

// Label-setting SSSP: binary-heap Dijkstra (lazy deletion) on pre-allocated scratch.
// Same FirstThruNode gating / mode_allowed_use / movement-restriction rules as the
// deque Minpath (settings sp_algorithm=1). Each node settled once; correct for
// non-negative link costs (BPR/conic/QVDF all >= 0).
// External node id (agency data) -> dense node seq 1..N.
// Returns 0 when the id is not present in node.csv.
static inline int NodeSeqFromExternal(int external_node_id)
{
	std::unordered_map<int, int>::const_iterator it =
		g_map_external_node_id_2_node_seq_no.find(external_node_id);
	return (it == g_map_external_node_id_2_node_seq_no.end()) ? 0 : it->second;
}

// External zone id (agency data) -> dense zone seq 1..Z (internal renumbering).
// Returns 0 when the id is not a zone in node.csv.
static inline int ZoneSeqFromExternal(int external_zone_id)
{
	std::unordered_map<int, int>::const_iterator it =
		g_map_old_zone_id_2_zone_seq.find(external_zone_id);
	return (it == g_map_old_zone_id_2_zone_seq.end()) ? 0 : it->second;
}

int Minpath_Dijkstra(int mode, int Orig, int* PredLink, double* CostTo)
{
	int tid = GetSPScratchThreadIndex();
	int* PrevLink = g_sp_PrevLink[tid].data();
	std::vector<std::pair<double, int>>& heap = g_sp_heap[tid];
	heap.clear();

	for (int node = 1; node <= no_nodes; node++)
	{
		CostTo[node] = BIGM;
		PredLink[node] = INVALID;
		PrevLink[node] = INVALID;
	}
	int origin = Orig;   // dense zone seq == centroid node seq (internal renumbering)
	CostTo[origin] = 0.0;
	heap.push_back(std::make_pair(0.0, origin));

	while (!heap.empty())
	{
		std::pop_heap(heap.begin(), heap.end(), g_sp_heap_gt);
		std::pair<double, int> top = heap.back();
		heap.pop_back();
		int now = top.second;
		if (top.first > CostTo[now])
			continue;                         // stale heap entry
		if (!(now >= FirstThruNode || now == origin))
			continue;                         // do not scan through a centroid/zone
		for (int k = FirstLinkFrom[now]; k <= LastLinkFrom[now]; k++)
		{
			if (Link[k].mode_allowed_use[mode] == 0)
				continue;
			if (PrevLink[now] != INVALID && Link[PrevLink[now]].b_withmovement_restrictions == true)
			{
				if (IsMovementRestricted(PrevLink[now], k))
					continue;
			}
			int NewNode = Link[k].internal_to_node_id;
			double NewCost = CostTo[now] + Link[k].Travel_time + Link[k].mode_AdditionalCost[mode];
			if (CostTo[NewNode] > NewCost)
			{
				CostTo[NewNode] = NewCost;
				PredLink[NewNode] = k;
				PrevLink[NewNode] = k;
				heap.push_back(std::make_pair(NewCost, NewNode));
				std::push_heap(heap.begin(), heap.end(), g_sp_heap_gt);
			}
		}
	}
	return 0;
}

int Minpath(int mode, int Orig, int* PredLink, double* CostTo)
{
	if (g_sp_algorithm == 1)
		return Minpath_Dijkstra(mode, Orig, PredLink, CostTo);
	int node, now, NewNode, k, Return2Q_Count = 0;
	double NewCost;
	int QueueFirst, QueueLast;

	int tid = GetSPScratchThreadIndex();
	int* QueueNext = g_sp_QueueNext[tid].data();   // pre-allocated, reused (no malloc)
	int* PrevLink = g_sp_PrevLink[tid].data();

	for (node = 1; node <= no_nodes; node++)
	{
		QueueNext[node] = INVALID;
		CostTo[node] = BIGM;
		PredLink[node] = INVALID;
		PrevLink[node] = INVALID;
	}

	now = Orig;   // dense zone seq == centroid node seq (internal renumbering)
	int internal_node_id_for_origin_zone = now;
	QueueNext[now] = WAS_IN_QUEUE;
	PredLink[now] = INVALID;
	PrevLink[now] = INVALID;
	CostTo[now] = 0.0;

	QueueFirst = QueueLast = INVALID;

	while ((now != INVALID) && (now != WAS_IN_QUEUE))
	{
		if (now >= FirstThruNode || now == internal_node_id_for_origin_zone)
		{
			for (k = FirstLinkFrom[now]; k <= LastLinkFrom[now]; k++)
			{
				if (Link[k].mode_allowed_use[mode] == 0)
					continue;

				NewNode = Link[k].internal_to_node_id;
				NewCost = CostTo[now] + Link[k].Travel_time + Link[k].mode_AdditionalCost[mode];

				// Check if the previous link (PrevLink[now]) has restrictions on this link (k)
				if (PrevLink[now] != INVALID && Link[PrevLink[now]].b_withmovement_restrictions == true)
				{
					if(IsMovementRestricted(PrevLink[now],k))
						continue;  // Skip restricted movement
				}

				if (CostTo[NewNode] > NewCost)
				{
					CostTo[NewNode] = NewCost;
					PredLink[NewNode] = k;
					PrevLink[NewNode] = k;

					if (QueueNext[NewNode] == WAS_IN_QUEUE)
					{
						QueueNext[NewNode] = QueueFirst;
						QueueFirst = NewNode;
						if (QueueLast == INVALID)
							QueueLast = NewNode;
						Return2Q_Count++;
					}
					else if (QueueNext[NewNode] == INVALID && NewNode != QueueLast)
					{
						if (QueueLast != INVALID)
						{
							QueueNext[QueueLast] = NewNode;
							QueueLast = NewNode;
						}
						else
						{
							QueueFirst = QueueLast = NewNode;
							QueueNext[QueueLast] = INVALID;
						}
					}
				}
			}
		}

		now = QueueFirst;
		if ((now == INVALID) || (now == WAS_IN_QUEUE))
			break;

		QueueFirst = QueueNext[now];
		QueueNext[now] = WAS_IN_QUEUE;
		if (QueueLast == now)
			QueueLast = INVALID;
	}

	return (Return2Q_Count);
}

/* Exact turn-restricted shortest path (link-state label correcting).
   Only called when g_has_movement_restrictions != 0. The search state is a
   directed link ("arrived at head(k) via link k"), so the allowed next moves
   depend on the incoming link -- which is exactly what turn restrictions need
   and what the one-label-per-node Minpath() cannot represent.

   Outputs (per origin):
     CostTo[node]        - least cost to reach node (min over incoming links)
     PredLink[node]      - best incoming link into node (start of back-trace)
     PredLinkOfLink[k]   - link traversed immediately before link k (the chain)
   Same per-link cost terms (Travel_time + mode_AdditionalCost) and the same
   FirstThruNode gating and mode_allowed_use filtering as Minpath(). Complexity
   is O(links + movements) per origin; arrays are O(links). */
int Minpath_TR(int mode, int Orig, int* PredLink, int* PredLinkOfLink, double* CostTo)
{
	int origin_node = Orig;   // dense zone seq == centroid node seq (internal renumbering)

	double* link_cost = (double*)Alloc_1D(number_of_links, sizeof(double));
	int* QueueNextLink = (int*)Alloc_1D(number_of_links, sizeof(int));

	for (int node = 1; node <= no_nodes; node++)
	{
		CostTo[node] = BIGM;
		PredLink[node] = INVALID;
	}
	for (int k = 1; k <= number_of_links; k++)
	{
		link_cost[k] = BIGM;
		PredLinkOfLink[k] = INVALID;
		QueueNextLink[k] = INVALID;
	}
	CostTo[origin_node] = 0.0;

	int QueueFirst = INVALID, QueueLast = INVALID;

	// Seed the queue with the origin node's outgoing links.
	for (int k = FirstLinkFrom[origin_node]; k <= LastLinkFrom[origin_node]; k++)
	{
		if (Link[k].mode_allowed_use[mode] == 0)
			continue;

		double c = Link[k].Travel_time + Link[k].mode_AdditionalCost[mode];
		link_cost[k] = c;
		PredLinkOfLink[k] = INVALID;

		int hn = Link[k].internal_to_node_id;
		if (CostTo[hn] > c) { CostTo[hn] = c; PredLink[hn] = k; }

		// enqueue k at the back
		if (QueueLast != INVALID) { QueueNextLink[QueueLast] = k; QueueLast = k; }
		else { QueueFirst = QueueLast = k; QueueNextLink[k] = INVALID; }
	}

	// Label-correcting (deque / Pape) scan over link states.
	while ((QueueFirst != INVALID) && (QueueFirst != WAS_IN_QUEUE))
	{
		int k = QueueFirst;                 // current state: just traversed link k
		QueueFirst = QueueNextLink[k];
		QueueNextLink[k] = WAS_IN_QUEUE;
		if (QueueLast == k) QueueLast = INVALID;

		int v = Link[k].internal_to_node_id;   // we are at node v
		// FirstThruNode gating: do not let through-traffic pass a centroid/zone
		// connector node (same rule as Minpath()), except at the origin.
		if (!(v >= FirstThruNode || v == origin_node))
			continue;

		for (int k2 = FirstLinkFrom[v]; k2 <= LastLinkFrom[v]; k2++)
		{
			if (Link[k2].mode_allowed_use[mode] == 0)
				continue;

			// Turn restriction on the movement (incoming link k -> outgoing k2).
			if (Link[k].b_withmovement_restrictions && IsMovementRestricted(k, k2))
				continue;

			double nc = link_cost[k] + Link[k2].Travel_time + Link[k2].mode_AdditionalCost[mode];
			if (link_cost[k2] > nc)
			{
				link_cost[k2] = nc;
				PredLinkOfLink[k2] = k;

				int hn = Link[k2].internal_to_node_id;
				if (CostTo[hn] > nc) { CostTo[hn] = nc; PredLink[hn] = k2; }

				if (QueueNextLink[k2] == WAS_IN_QUEUE)        // re-add to front
				{
					QueueNextLink[k2] = QueueFirst;
					QueueFirst = k2;
					if (QueueLast == INVALID) QueueLast = k2;
				}
				else if (QueueNextLink[k2] == INVALID && k2 != QueueLast)  // add to back
				{
					if (QueueLast != INVALID) { QueueNextLink[QueueLast] = k2; QueueLast = k2; }
					else { QueueFirst = QueueLast = k2; QueueNextLink[k2] = INVALID; }
				}
			}
		}
	}

	free(link_cost);
	free(QueueNextLink);
	return 0;
}


/* Find minimum cost routes .
Input: 	None
Output:	RouteCost - route generalized cost, by origin and destination
		MinPathSuccLink - trees of minimum cost routes, by destination and node. */
std::vector<std::vector<int>> Processor_origin_zones;

int g_number_of_processors = 4;
static constexpr int MAX_SAFE_PROCESSOR_BUCKETS = 128;

int ProcessorCountValidationStatus(int requested_processors)
{
	return requested_processors >= 1 &&
		requested_processors <= MAX_SAFE_PROCESSOR_BUCKETS ? 0 : 2;
}

static bool ValidateProcessorCount()
{
	if (ProcessorCountValidationStatus(g_number_of_processors) == 0)
		return true;

	fprintf(stderr,
		"ERROR: number_of_processors=%d is outside the accepted range [1, %d]. "
		"This safety ceiling prevents unreasonable OpenMP thread requests and "
		"processor-dependent memory allocations.\n",
		g_number_of_processors, MAX_SAFE_PROCESSOR_BUCKETS);
	if (summary_log_file != NULL)
	{
		fprintf(summary_log_file,
			"ERROR: number_of_processors=%d is outside the accepted range [1, %d]. "
			"This safety ceiling prevents unreasonable OpenMP thread requests and "
			"processor-dependent memory allocations.\n",
			g_number_of_processors, MAX_SAFE_PROCESSOR_BUCKETS);
		fflush(summary_log_file);
	}
	return false;
}

static bool ConfigureOpenMPRuntime()
{
	int compiled = 0;
	int openmp_version = 0;
	int max_threads = 1;
	int num_procs = 1;
	int probe_team_size = 1;
	int dynamic = 0;

#ifdef _OPENMP
	compiled = 1;
	openmp_version = _OPENMP;
	omp_set_dynamic(0);
	omp_set_num_threads(g_number_of_processors);
	max_threads = omp_get_max_threads();
	num_procs = omp_get_num_procs();
	dynamic = omp_get_dynamic();
#pragma omp parallel
	{
#pragma omp single
		probe_team_size = omp_get_num_threads();
	}
#endif

	printf("OpenMP status: compiled=%d, _OPENMP=%d, requested_processors=%d, "
		"max_threads=%d, num_procs=%d, probe_team_size=%d, dynamic=%d\n",
		compiled, openmp_version, g_number_of_processors, max_threads,
		num_procs, probe_team_size, dynamic);
	if (summary_log_file != NULL)
	{
		fprintf(summary_log_file,
			"OpenMP status: compiled=%d, _OPENMP=%d, requested_processors=%d, "
			"max_threads=%d, num_procs=%d, probe_team_size=%d, dynamic=%d\n",
			compiled, openmp_version, g_number_of_processors, max_threads,
			num_procs, probe_team_size, dynamic);
		fflush(summary_log_file);
	}

	if (compiled && g_number_of_processors > num_procs)
	{
		fprintf(stderr,
			"WARNING: number_of_processors=%d exceeds the %d processors reported by "
			"the OpenMP runtime; oversubscription may reduce performance.\n",
			g_number_of_processors, num_procs);
		if (summary_log_file != NULL)
		{
			fprintf(summary_log_file,
				"WARNING: number_of_processors=%d exceeds the %d processors reported by "
				"the OpenMP runtime; oversubscription may reduce performance.\n",
				g_number_of_processors, num_procs);
			fflush(summary_log_file);
		}
	}

	if (g_number_of_processors > 1 && probe_team_size == 1)
	{
		fprintf(stderr,
			"WARNING: %d processors were requested, but the OpenMP probe created only one thread. "
			"Check the native build and external OpenMP runtime limits.\n",
			g_number_of_processors);
		if (summary_log_file != NULL)
		{
			fprintf(summary_log_file,
				"WARNING: %d processors were requested, but the OpenMP probe created only one thread.\n",
				g_number_of_processors);
			fflush(summary_log_file);
		}
	}

	return true;
}

static bool InitializeProcessorOriginZones()
{
	try
	{
		Processor_origin_zones.assign(
			static_cast<std::size_t>(g_number_of_processors),
			std::vector<int>());
	}
	catch (const std::bad_alloc&)
	{
		fprintf(stderr,
			"ERROR: unable to allocate %d processor-origin buckets. "
			"Reduce number_of_processors or free memory and try again.\n",
			g_number_of_processors);
		if (summary_log_file != NULL)
		{
			fprintf(summary_log_file,
				"ERROR: unable to allocate %d processor-origin buckets. "
				"Reduce number_of_processors or free memory and try again.\n",
				g_number_of_processors);
			fflush(summary_log_file);
		}
		return false;
	}
	return true;
}

static void WarnIfProcessorsExceedOrigins(int assignable_origin_zones)
{
	if (g_number_of_processors <= assignable_origin_zones)
		return;

	fprintf(stderr,
		"WARNING: number_of_processors=%d exceeds the %d assignable origin zones; "
		"some processor buckets will be empty.\n",
		g_number_of_processors, assignable_origin_zones);
	if (summary_log_file != NULL)
	{
		fprintf(summary_log_file,
			"WARNING: number_of_processors=%d exceeds the %d assignable origin zones; "
			"some processor buckets will be empty.\n",
			g_number_of_processors, assignable_origin_zones);
		fflush(summary_log_file);
	}
}

double FindMinCostRoutes(int*** MinPathPredLink)
{
	double** CostTo;

	// CostTo is per-node shortest-path scratch that is consumed immediately after
	// each Minpath() call (read into MDRouteCost in the same loop iteration), so it
	// only needs one row PER PROCESSOR, not per origin zone. Sizing it by
	// g_number_of_processors instead of no_zones removes the O(zones x nodes)
	// growth with the OD count (the legacy accessibility code uses the same
	// processor-indexed scratch trick). Each parallel iteration p owns row p.
	// With SP-equivalence groups (g_rep_mode) a representative's cost row must
	// additionally stay valid until its member modes have consumed it within the
	// same origin iteration, so each processor owns g_sp_scratch_rows rows:
	// row 0 is the shared immediate-consumption scratch, plus one retained row
	// per representative that shares its tree (g_sp_scratch_row[m]).
	CostTo = (double**)Alloc_2D(g_number_of_processors * g_sp_scratch_rows, no_nodes, sizeof(double));
	StatusMessage("Minpath", "Starting the minpath calculations.");
	double* system_least_travel_time_org_zone = (double*)Alloc_1D(no_zones, sizeof(double));

#pragma omp parallel for
	for (int p = 0; p < g_number_of_processors; p++)
	{

		for (int i = 0; i < Processor_origin_zones[p].size(); i++)
		{

		int Orig = Processor_origin_zones[p][i];  // get origin zone id


		if(TotalOFlow[Orig] < 0.00001)  // only work on positive zone flow 
			continue;

		system_least_travel_time_org_zone[Orig] = 0;  // reset it before mode based computing 

		for (int m = 1; m <= number_of_modes; m++)
		{
			if (g_mode_type_vector[m].dedicated_shortest_path == 0)  // skip the shortest path computing
				continue;

			int rep = (g_rep_mode[m] >= 1) ? g_rep_mode[m] : m;   // SP-equivalence group representative
			double* cost_row = CostTo[p * g_sp_scratch_rows + g_sp_scratch_row[m]];

			if (rep == m)   // representative: computes the tree (SP-equivalent members reuse it)
			{
				// zero-demand-origin skip: no SP when neither this mode nor any
				// member of its SP group has demand out of this origin. Exact:
				// the destination loop below then BIGM-fills the row, which is
				// precisely what the full loop produces with zero demand.
				if (g_tree_origin_active == NULL || g_tree_origin_active[m][Orig])
				{
					if (g_has_movement_restrictions)
						Minpath_TR(m, Orig, MinPathPredLink[m][Orig], g_MinPathPredLink_link[m][Orig], cost_row);
					else
						Minpath(m, Orig, MinPathPredLink[m][Orig], cost_row);
				}
			}
			// else: SP-equivalent member -- rep (< m) already filled cost_row and
			// MinPathPredLink[rep][Orig] for this origin; identical inputs make
			// Minpath deterministic, so reuse is bit-identical.

			if (MDRouteCost != NULL)
			{
				if (g_mode_origin_active != NULL && !g_mode_origin_active[m][Orig])
				{
					// no demand for (m, Orig): the loop below would leave every
					// destination at BIGM (no MDODflow gate passes), so write
					// that directly without consulting cost_row (possibly stale).
					for (int Dest = 1; Dest <= no_zones; Dest++)
						MDRouteCost[m][Orig][Dest] = BIGM;
				}
				else
				for (int Dest = 1; Dest <= no_zones; Dest++)
				{
					MDRouteCost[m][Orig][Dest] = BIGM; // initialization

					if (MDODflow[m][Orig][Dest] > 0.000001)
					{
						if (cost_row[Dest] < 0.0)
							ExitMessage("Negative cost %lg from Origin %d to Destination %d.",
								(double)cost_row[Dest], Orig, Dest);

						// CostTo is coded as internal node id
						int  internal_node_id_for_destination_zone  = Dest;   // dense zone seq == node seq

						if (cost_row[internal_node_id_for_destination_zone] <= BIGM - 1)  // feasible cost
						{
							MDRouteCost[m][Orig][Dest] = cost_row[internal_node_id_for_destination_zone];

							system_least_travel_time_org_zone[Orig] += MDRouteCost[m][Orig][Dest] * MDODflow[m][Orig][Dest] * g_mode_type_vector[m].pce;

						}
						else
						{
							int debug_flag = 1;
						}
					}
				}
			}
		}
	}

	}

	double system_least_travel_time = 0;
	for (int Orig = 1; Orig <= no_zones; Orig++)
	{
		system_least_travel_time += system_least_travel_time_org_zone[Orig];
	}


	// free(CostTo);
	Free_2D((void**)CostTo, g_number_of_processors * g_sp_scratch_rows, no_nodes);
	free(system_least_travel_time_org_zone);
	StatusMessage("Minpath", "Found all minpath.");
	return system_least_travel_time;
}

// Structure to hold path information for each OD pair
struct ODPathInfo {
	int mode = 0;
	double distance = -1;
	double freeFlowTime = 0.0;
	double congestionTime = 0.0;
	double cost = -1.0;
};




//#include <iomanip>

void WriteOutputFiles(const char* filename, ODPathInfo** odPathInfoMatrix)
{


	auto start = std::chrono::high_resolution_clock::now();
	std::ofstream outputFile(filename);
	if (!outputFile.is_open()) {
		std::cerr << "Error: Could not open output od_performance.csv file." << std::endl;
		return;
	}

	// Write headers (if needed)
	outputFile << "o_zone_id,d_zone_id,total_distance_mile,total_distance_km,total_free_flow_travel_time\n";

	// Global statistics initialization
	double grand_totalDistance = 0.0, grand_totalFreeFlowTravelTime = 0.0, grand_totalTravelTime = 0.0, grand_total_count = 0;

	// Set fixed precision to 2 decimal places to reduce file size
	outputFile << std::fixed << std::setprecision(3);

	for (int Orig = 0; Orig < number_of_internal_zones; ++Orig) {
		for (int Dest = 0; Dest < number_of_internal_zones; ++Dest) {
			if (Orig == Dest)
				continue;

			const ODPathInfo& pathInfo = odPathInfoMatrix[Orig][Dest];
			if (pathInfo.cost < 0.0)
				continue;

			float volume = 1;

			grand_totalDistance += pathInfo.distance * volume;
			grand_totalFreeFlowTravelTime += pathInfo.freeFlowTime * volume;
			grand_totalTravelTime += pathInfo.congestionTime * volume;
			grand_total_count += volume;

			outputFile << g_zone_seq_2_old_zone_id[Orig + 1] << ","
				<< g_zone_seq_2_old_zone_id[Dest + 1] << ","
				<< pathInfo.distance << "," << (pathInfo.distance * 1.609) << ","
				<< pathInfo.freeFlowTime << "\n";
		}
	}

	outputFile.close();
	std::cout << "Output written to " << filename << std::endl;

	if (grand_total_count > 0) {
		std::cout << "---------- Summary Statistics ----------\n"
			<< "Average path distance: " << grand_totalDistance / grand_total_count << " miles\n"
			<< "Average free flow travel time: " << grand_totalFreeFlowTravelTime / grand_total_count << " minutes\n";
	}

	// Calculate elapsed time
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = end - start;
	auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
	auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration % std::chrono::hours(1));
	auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration % std::chrono::minutes(1));
	auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(duration % std::chrono::seconds(1));

	printf("Printing out OD accessibility: %lld hours %lld minutes %lld seconds %lld ms\n",
		hours.count(), minutes.count(), seconds.count(), milliseconds.count());
}


// Aggregated statistics for a zone
struct ZoneAggregates {
	double totalDistance = 0.0;
	double totalFreeFlowTime = 0.0;
	double totalCongestionTime = 0.0;
	int count = 0;
};


void WriteZoneAccessibilityCSV(const char* filename, ODPathInfo** odPathInfoMatrix)
{
	auto start = std::chrono::high_resolution_clock::now();

	std::ofstream outputFile(filename);
	if (!outputFile.is_open()) {
		std::cerr << "Error: Could not open output zone_accessibility.csv file." << std::endl;
		return;
	}

	// Write header line.
	// Columns:
	// zone_id,
	// origin_count, origin_avg_distance_mile, origin_avg_distance_km, origin_avg_free_flow, origin_avg_congestion,
	// destination_count, destination_avg_distance_mile, destination_avg_distance_km, destination_avg_free_flow, destination_avg_congestion
	outputFile << "zone_id,origin_count,origin_avg_distance_mile,origin_avg_distance_km,origin_avg_free_flow,origin_avg_congestion,"
		<< "destination_count,destination_avg_distance_mile,destination_avg_distance_km,destination_avg_free_flow,destination_avg_congestion\n";

	// Create aggregation vectors (one element per zone).
	std::vector<ZoneAggregates> originAgg(number_of_internal_zones);
	std::vector<ZoneAggregates> destAgg(number_of_internal_zones);

	// Process every OD pair.
	for (int Orig = 0; Orig < number_of_internal_zones; ++Orig) {
		for (int Dest = 0; Dest < number_of_internal_zones; ++Dest) {
			if (Orig == Dest)
				continue;  // Skip self trips.
			const ODPathInfo& pathInfo = odPathInfoMatrix[Orig][Dest];
			if (pathInfo.cost < 0.0)
				continue;  // Skip unreachable pairs.

			// Update Origin aggregates for zone Orig.
			originAgg[Orig].totalDistance += pathInfo.distance;
			originAgg[Orig].totalFreeFlowTime += pathInfo.freeFlowTime;
			originAgg[Orig].totalCongestionTime += pathInfo.congestionTime;
			originAgg[Orig].count++;

			// Update Destination aggregates for zone Dest.
			destAgg[Dest].totalDistance += pathInfo.distance;
			destAgg[Dest].totalFreeFlowTime += pathInfo.freeFlowTime;
			destAgg[Dest].totalCongestionTime += pathInfo.congestionTime;
			destAgg[Dest].count++;
		}
	}

	// Set fixed precision for output numbers.
	outputFile << std::fixed << std::setprecision(2);

	// Write out the aggregated values for each zone.
	for (int zoneIdx = 0; zoneIdx < number_of_internal_zones; ++zoneIdx) {
		// Map the internal zone index to a node (zone) id.
		int zone_id = g_zone_seq_2_old_zone_id[zoneIdx + 1];

		// Origin measures.
		int origCount = originAgg[zoneIdx].count;
		double origAvgDistance = (origCount > 0) ? (originAgg[zoneIdx].totalDistance / origCount) : 0.0;
		double origAvgDistanceKM = origAvgDistance * 1.609;
		double origAvgFreeFlow = (origCount > 0) ? (originAgg[zoneIdx].totalFreeFlowTime / origCount) : 0.0;
		double origAvgCongestion = (origCount > 0) ? (originAgg[zoneIdx].totalCongestionTime / origCount) : 0.0;

		// Destination measures.
		int destCount = destAgg[zoneIdx].count;
		double destAvgDistance = (destCount > 0) ? (destAgg[zoneIdx].totalDistance / destCount) : 0.0;
		double destAvgDistanceKM = destAvgDistance * 1.609;
		double destAvgFreeFlow = (destCount > 0) ? (destAgg[zoneIdx].totalFreeFlowTime / destCount) : 0.0;
		double destAvgCongestion = (destCount > 0) ? (destAgg[zoneIdx].totalCongestionTime / destCount) : 0.0;

		// Write CSV row.
		outputFile << zone_id << ","
			<< origCount << "," << origAvgDistance << "," << origAvgDistanceKM << "," << origAvgFreeFlow << "," << origAvgCongestion << ","
			<< destCount << "," << destAvgDistance << "," << destAvgDistanceKM << "," << destAvgFreeFlow << "," << destAvgCongestion << "\n";
	}

	outputFile.close();
	std::cout << "Zone-based accessibility output written to " << filename << std::endl;

	auto end = std::chrono::high_resolution_clock::now();
	auto duration = end - start;
	auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
	auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration % std::chrono::hours(1));
	auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration % std::chrono::minutes(1));
	auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(duration % std::chrono::seconds(1));

	printf("Elapsed time: %lld hours %lld minutes %lld seconds %lld ms\n",
		hours.count(), minutes.count(), seconds.count(), milliseconds.count());
}
// Function to compute accessibility and OD costs
int ComputeAccessibilityAndODCosts_v1(const char* filename)
{
	// Start timing
	auto start0 = std::chrono::high_resolution_clock::now();

	// Allocate memory for cost matrix
	double** CostTo = (double**)Alloc_2D(g_number_of_processors, no_nodes+1, sizeof(double));



	for (int k = 1; k <= number_of_links; k++)
		Link[k].Travel_time = Link[k].FreeTravelTime;

	// Store the link sequences for each OD pair
	// Simplified to 2D vector: first dimension is origin zone, second is destination zone

	ODPathInfo** odPathInfoMatrix = (ODPathInfo**)Alloc_2D(number_of_internal_zones, number_of_internal_zones, sizeof(ODPathInfo));

	// Allocate memory for predecessor links
	int** PredLink = (int**)Alloc_2D(g_number_of_processors, no_nodes + 1, sizeof(int));

	cout << " Memory allocation completes. Starting the minpath calculations." << endl; 

	// Parallel processing of origin zones
#pragma omp parallel for
	for (int p = 0; p < g_number_of_processors; p++)
	{
		for (int i = 0; i < Processor_origin_zones[p].size(); i++)
		{
			int Orig = Processor_origin_zones[p][i];  // get origin zone id

			if (p == 0)
				cout << "Accessibility computing for zone " << g_zone_seq_2_old_zone_id[Orig] << endl;

			// Process each mode
			for (int m = 1; m <= number_of_modes; m++)
			{
				// Compute shortest path from origin to all destinations
				Minpath(m, Orig, PredLink[p], CostTo[p]);

				// Store path information for each destination
				for (int j = 1; j <= no_nodes; j++)
				{
					if(g_node_vector[j].zone_id >0)
					{
						int Dest = g_node_vector[j].zone_id;
						if (CostTo[p][j] < BIGM - 1)
						{
							// Calculate path metrics
							double pathDistance = 0.0;
							double freeFlowTime = 0.0;
							double congestionTime = 0.0;

							int currentNode = Dest;   // dense zone seq == node seq

							while (currentNode != Orig && currentNode != 0)
							{
								int linkId = PredLink[p][currentNode];
								if (linkId == -1)
									break;
								pathDistance += Link[linkId].length;
								freeFlowTime += Link[linkId].FreeTravelTime;
								congestionTime += Link[linkId].Travel_time;
								currentNode = Link[linkId].internal_from_node_id;
							}

							int internal_node_id = Orig;   // dense zone seq == node seq
							int internal_o_zone_no = g_node_vector[internal_node_id].internal_zone_no;
							int internal_d_zone_no = g_node_vector[j].internal_zone_no;
							// Store in our 2D matrix
							ODPathInfo& pathInfo = odPathInfoMatrix[internal_o_zone_no][internal_d_zone_no];
							pathInfo.mode = m;
							pathInfo.distance = pathDistance;
							pathInfo.freeFlowTime = freeFlowTime;
							pathInfo.congestionTime = congestionTime;
							pathInfo.cost = CostTo[p][Dest];
						}
						else
						{
							int internal_node_id = Orig;   // dense zone seq == node seq
							int internal_o_zone_no = g_node_vector[internal_node_id].internal_zone_no;
							int internal_d_zone_no = g_node_vector[j].internal_zone_no;
							// Store in our 2D matrix
							ODPathInfo& pathInfo = odPathInfoMatrix[internal_o_zone_no][internal_d_zone_no];
							pathInfo.mode = m;
							pathInfo.distance = -1;
							pathInfo.freeFlowTime = -1;
							pathInfo.congestionTime = -1;
							pathInfo.cost = -1;
						}
					}
				}
			}
		}
	}

	// Free the memory for costs and predecessor links
	Free_2D((void**)CostTo, g_number_of_processors, no_nodes+1);
	Free_2D((void**)PredLink, g_number_of_processors, no_nodes+1);

	// Calculate elapsed time
	auto end0 = std::chrono::high_resolution_clock::now();
	auto duration = end0 - start0;

	// Convert to hours, minutes, seconds
	auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
	auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration % std::chrono::hours(1));
	auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration % std::chrono::minutes(1));
	auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(duration % std::chrono::seconds(1));

	printf("All OD accessibility computing: %lld hours %lld minutes %lld seconds %lld ms\n",
		hours.count(), minutes.count(), seconds.count(), milliseconds.count());

	// Output the results to files

	if(number_of_internal_zones<10000)
		WriteOutputFiles(filename, odPathInfoMatrix);

	WriteZoneAccessibilityCSV("zone_accessibility.csv", odPathInfoMatrix);

	Free_2D((void**)odPathInfoMatrix, number_of_internal_zones, number_of_internal_zones);

	return 1;
}

// Function to compute accessibility and OD costs
int ComputeAccessibilityAndODCosts_v2(const char* filename)
{
	auto start0 = std::chrono::high_resolution_clock::now();

	// Allocate cost and predecessor arrays by the configured processor count.
	double** CostTo = (double**)Alloc_2D(g_number_of_processors, no_nodes + 1, sizeof(double));
	int** PredLink = (int**)Alloc_2D(g_number_of_processors, no_nodes + 1, sizeof(int));



	// Set each link's travel time to its free travel time.
	for (int k = 1; k <= number_of_links; k++) {
		Link[k].Travel_time = Link[k].FreeTravelTime;
	}

	// Prepare aggregation vectors for zone-based statistics.
	std::vector<ZoneAggregates> originAgg(number_of_internal_zones);
	std::vector<ZoneAggregates> destAgg(number_of_internal_zones);

	// Open the OD performance output file and write header.
	std::ofstream outputFile(filename);
	if (!outputFile.is_open()) {
		std::cerr << "Error: Could not open output file " << filename << std::endl;
		return -1;
	}
	outputFile << "o_zone_id,d_zone_id,total_distance_mile,total_distance_km,total_free_flow_travel_time\n";
	outputFile << std::fixed << std::setprecision(3);

	// Define the batch size. For example, process 100 origin zones at a time.
	const int batch_size = 100;
	// Loop over origin zones in batches.
	for (int batchStart = 0; batchStart < number_of_internal_zones; batchStart += batch_size) {

		cout << "outputing accessibility for zone id: \n" << batchStart + 1 << endl;

		int currentBatchSize = std::min(batch_size, number_of_internal_zones - batchStart);

		// Allocate a temporary OD matrix for this batch:
		// Each row corresponds to an origin in [batchStart, batchStart + currentBatchSize)
		// and each column corresponds to a destination (internal zone index).
		ODPathInfo** batchMatrix = (ODPathInfo**)Alloc_2D(currentBatchSize, number_of_internal_zones, sizeof(ODPathInfo));

		// Initialize the batchMatrix.
		for (int i = 0; i < currentBatchSize; ++i) {
			for (int j = 0; j < number_of_internal_zones; ++j) {
				batchMatrix[i][j].mode = 0;
				batchMatrix[i][j].distance = -1;
				batchMatrix[i][j].freeFlowTime = 0.0;
				batchMatrix[i][j].congestionTime = 0.0;
				batchMatrix[i][j].cost = -1.0;
			}
		}

		// Process each origin zone in the current batch.
		for (int localOrigin = 0; localOrigin < currentBatchSize; ++localOrigin) {
			int internal_o_zone_no = batchStart + localOrigin;
			int Orig = internal_o_zone_no + 1; // dense zone seq == centroid node seq

			// Process each transportation mode.
			for (int m = 1; m <= number_of_modes; m++) {
				// Compute shortest path from origin Orig to all destinations.
				// (For simplicity, we use processor index 0. Adapt if you are using parallel processing.)
				Minpath(m, Orig, PredLink[0], CostTo[0]);

				// Process each destination node.
				for (int j = 1; j <= no_nodes; j++) {
					if (g_node_vector[j].zone_id > 0) {
						int Dest = g_node_vector[j].zone_id;
						if (CostTo[0][j] < BIGM - 1) {
							double pathDistance = 0.0;
							double freeFlowTime = 0.0;
							double congestionTime = 0.0;

							int currentNode = Dest;   // dense zone seq == node seq
							// Backtrack using predecessor links to accumulate path metrics.
							while (currentNode != Orig && currentNode != 0) {
								int linkId = PredLink[0][currentNode];
								if (linkId == -1)
									break;
								pathDistance += Link[linkId].length;
								freeFlowTime += Link[linkId].FreeTravelTime;
								congestionTime += Link[linkId].Travel_time;
								currentNode = Link[linkId].internal_from_node_id;
							}

							int internal_d_zone_no = g_node_vector[j].internal_zone_no;

							// Update batchMatrix for this origin-destination pair.
							ODPathInfo& pathInfo = batchMatrix[localOrigin][internal_d_zone_no];
							pathInfo.mode = m;
							pathInfo.distance = pathDistance;
							pathInfo.freeFlowTime = freeFlowTime;
							pathInfo.congestionTime = congestionTime;
							pathInfo.cost = CostTo[0][j];

							// Update aggregation if the path is reachable.
							if (pathInfo.cost >= 0.0) {
								originAgg[internal_o_zone_no].totalDistance += pathDistance;
								originAgg[internal_o_zone_no].totalFreeFlowTime += freeFlowTime;
								originAgg[internal_o_zone_no].totalCongestionTime += congestionTime;
								originAgg[internal_o_zone_no].count++;

								destAgg[internal_d_zone_no].totalDistance += pathDistance;
								destAgg[internal_d_zone_no].totalFreeFlowTime += freeFlowTime;
								destAgg[internal_d_zone_no].totalCongestionTime += congestionTime;
								destAgg[internal_d_zone_no].count++;
							}
						}
						else {
							// Unreachable destination.
							int internal_d_zone_no = g_node_vector[j].internal_zone_no;
							ODPathInfo& pathInfo = batchMatrix[localOrigin][internal_d_zone_no];
							pathInfo.mode = m;
							pathInfo.distance = -1;
							pathInfo.freeFlowTime = -1;
							pathInfo.congestionTime = -1;
							pathInfo.cost = -1;
						}
					}
				} // End for each destination
			} // End for each mode

			// After processing all modes for the current origin, write its OD data.
			for (int Dest = 0; Dest < number_of_internal_zones; ++Dest) {
				ODPathInfo& info = batchMatrix[localOrigin][Dest];
				// Skip self-trips and unreachable pairs.
				if (info.cost < 0 || internal_o_zone_no == Dest)
					continue;
				outputFile << g_zone_seq_2_old_zone_id[internal_o_zone_no + 1] << ","
					<< g_zone_seq_2_old_zone_id[Dest + 1] << ","
					<< info.distance << "," << (info.distance * 1.609) << ","
					<< info.freeFlowTime << "\n";
			}
		} // End for each origin in batch

		// Free the temporary batch matrix.
		Free_2D((void**)batchMatrix, currentBatchSize, number_of_internal_zones);
	} // End batch loop

	// Free the cost and predecessor arrays.
	Free_2D((void**)CostTo, g_number_of_processors, no_nodes + 1);
	Free_2D((void**)PredLink, g_number_of_processors, no_nodes + 1);

	// Write the zone accessibility aggregated statistics to a separate CSV file.
	std::ofstream zoneFile("zone_accessibility.csv");
	if (!zoneFile.is_open()) {
		std::cerr << "Error: Could not open zone_accessibility.csv file." << std::endl;
		return -1;
	}


		cout << "outputing zone_accessibility.csv file." << std::endl;

		zoneFile << "zone_id,origin_count,origin_avg_distance_mile,origin_avg_distance_km,origin_avg_free_flow,origin_avg_congestion,"
		<< "destination_count,destination_avg_distance_mile,destination_avg_distance_km,destination_avg_free_flow,destination_avg_congestion\n";
	zoneFile << std::fixed << std::setprecision(2);
	for (int zoneIdx = 0; zoneIdx < number_of_internal_zones; zoneIdx++) {
		int zone_id = g_zone_seq_2_old_zone_id[zoneIdx + 1];

		int origCount = originAgg[zoneIdx].count;
		double origAvgDistance = (origCount > 0) ? (originAgg[zoneIdx].totalDistance / origCount) : 0.0;
		double origAvgDistanceKM = origAvgDistance * 1.609;
		double origAvgFreeFlow = (origCount > 0) ? (originAgg[zoneIdx].totalFreeFlowTime / origCount) : 0.0;
		double origAvgCongestion = (origCount > 0) ? (originAgg[zoneIdx].totalCongestionTime / origCount) : 0.0;

		int destCount = destAgg[zoneIdx].count;
		double destAvgDistance = (destCount > 0) ? (destAgg[zoneIdx].totalDistance / destCount) : 0.0;
		double destAvgDistanceKM = destAvgDistance * 1.609;
		double destAvgFreeFlow = (destCount > 0) ? (destAgg[zoneIdx].totalFreeFlowTime / destCount) : 0.0;
		double destAvgCongestion = (destCount > 0) ? (destAgg[zoneIdx].totalCongestionTime / destCount) : 0.0;

		zoneFile << zone_id << ","
			<< origCount << "," << origAvgDistance << "," << origAvgDistanceKM << "," << origAvgFreeFlow << "," << origAvgCongestion << ","
			<< destCount << "," << destAvgDistance << "," << destAvgDistanceKM << "," << destAvgFreeFlow << "," << destAvgCongestion << "\n";
	}
	zoneFile.close();

	auto end0 = std::chrono::high_resolution_clock::now();
	auto duration = end0 - start0;
	auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
	auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration % std::chrono::hours(1));
	auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration % std::chrono::minutes(1));
	auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(duration % std::chrono::seconds(1));

	printf("All OD accessibility computing: %lld hours %lld minutes %lld seconds %lld ms\n",
		hours.count(), minutes.count(), seconds.count(), milliseconds.count());

	return 1;
}

/* Assign OD flows to links according to the routes in MinPathPredLink. */

// Define a global 3D vector to store link indices for each OD pair

// Global 5D vector for storing link sequences
std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> linkIndices;

void AddLinkSequence(int m, int Orig, int Dest, int route_id, const std::vector<int>& linkIDs)
{
	if (linkIndices.size() == 0)
		return; 

	// Ensure we are within bounds before adding the link sequence
	if (Orig > 0 && Orig < linkIndices[m].size() &&
		Dest > 0 && Dest < linkIndices[m][Orig].size() &&
		route_id >= 0 && route_id < linkIndices[m][Orig][Dest].size())
	{
		linkIndices[m][Orig][Dest][route_id] = linkIDs;  // Add link sequence to the 5D vector
	}
	else
	{
		std::cerr << "Error: Invalid indices for adding link sequence." << std::endl;
	}
}





void All_or_Nothing_Assign(int Assignment_iteration_no, double*** ODflow, int*** MinPathPredLink, double* Volume)
{
//	printf("All or nothing assignment\n");

	auto start0 = std::chrono::high_resolution_clock::now();
	double** ProcessorVolume;
	double*** ProcessorModeVolume;


	ProcessorVolume = (double**)Alloc_2D(number_of_links, g_number_of_processors, sizeof(double));

	if (ProcessorVolume == nullptr) {
		std::cerr << "Error: Memory allocation for ProcessorVolume failed." << std::endl;
		// Handle the error, e.g., exit the program, clean up resources, etc.
		exit(EXIT_FAILURE);
	}

	ProcessorModeVolume = (double***)Alloc_3D(number_of_links, number_of_modes, g_number_of_processors, sizeof(double));
	if (ProcessorModeVolume == nullptr) {
		std::cerr << "Error: Memory allocation for ProcessorModeVolume failed." << std::endl;
		// Handle the error, e.g., exit the program, clean up resources, etc.
		exit(EXIT_FAILURE);
	}

#pragma omp parallel for
	for (int k = 1; k <= number_of_links; k++)
	{


		for (int p = 0; p < g_number_of_processors; p++)
		{
			ProcessorVolume[k][p] = 0;

			for (int m = 1; m <= number_of_modes; m++)
				ProcessorModeVolume[k][m][p] = 0;
		}
	}
	//StatusMessage("Assign", "Starting assign.");

	if(Assignment_iteration_no == 0)
	{
		printf("The list of zero-volume zones:");
		for (int Orig = 1; Orig <= no_zones; Orig++)
		{
			
			if (TotalOFlow[Orig] < 0.00001)  // only work on positive zone flow 
			{ 
				printf("%d,", g_zone_seq_2_old_zone_id[Orig]);
			}

		}

		printf("\n");

		printf("The list of zones without outbound connecting links:");
		for (int Orig = 1; Orig <= no_zones; Orig++)
		{

			if (zone_outbound_link_size[Orig] == 0)  // there is no outbound link from the origin 
			{
				printf("%d,", g_zone_seq_2_old_zone_id[Orig]);
			}
		}
		printf("\n");
	}
#pragma omp parallel for
	for (int p = 0; p < g_number_of_processors; p++)
	{

		for (int i = 0; i < Processor_origin_zones[p].size(); i++)
		{

			int Orig = Processor_origin_zones[p][i];  // get origin zone id
			
			if (TotalOFlow[Orig] < 0.00001)  // only work on positive zone flow 
				continue;

			if (zone_outbound_link_size[Orig] == 0)  // there is no outbound link from the origin 
				continue; 


			int Dest, k;
			int CurrentNode;
			double RouteFlow;
			std::vector<int> currentLinkSequence; // Temporary vector to store link indices


			for (int m = 1; m <= number_of_modes; m++)
			{


				//	printf("Assign", "Assigning origin %6d.", Orig);
				for (Dest = 1; Dest <= no_zones; Dest++)  // Dest zone 
				{
					if (Dest == Orig)
						continue;

					if (shortest_path_log_flag || Assignment_iteration_no == 0)
						currentLinkSequence.clear();

					RouteFlow = ODflow[m][Orig][Dest];  //test
					if (RouteFlow == 0)
						continue;

					if (g_mode_type_vector[m].dedicated_shortest_path == 0)  // skip the shortest path computing
					{
						if (MDRouteCost[1][Orig][Dest] >= BIGM - 1)  // if feasible 
							continue;
					}
					else
					{
						if (MDRouteCost[m][Orig][Dest] >= BIGM - 1)  // if feasible 
							continue;
					}

					CurrentNode = Dest;
					CurrentNode = Dest;   // dense zone seq == node seq (internal renumbering)
					int internal_node_for_origin_node = Orig;   // dense zone seq == node seq (internal renumbering)
					// MinPathPredLink is coded as internal node id 
					// 
					//double total_travel_time = 0;
					//double total_length = 0;
					//double total_FFTT = 0;

					int _trace_hops = 0;
					// mode slot holding this OD's predecessor tree (legacy
					// dedicated_shortest_path=0 modes share mode 1's tree;
					// SP-equivalent modes share their group representative's tree,
					// g_rep_mode). prev_k drives the link-state back-trace when
					// turn restrictions are active.
					int mpred = (g_mode_type_vector[m].dedicated_shortest_path == 0)
						? 1 : ((g_rep_mode[m] >= 1) ? g_rep_mode[m] : m);
					int prev_k = INVALID;
					while (CurrentNode != internal_node_for_origin_node)
					{
						// Robustness guard: a malformed SP tree (e.g. a pred cycle
						// from cordon/First-Through-Node edge cases) would make this
						// back-trace never reach the origin -> infinite loop / hang.
						// Cap hops at the node count and skip the OD if exceeded.
						// A turn-restricted shortest path is a simple path in LINK-state
						// space (each directed link used at most once), so bound hops by
						// links there; otherwise by nodes.
						int _hop_cap = g_has_movement_restrictions ? (number_of_links + 2) : (no_nodes + 2);
						if (++_trace_hops > _hop_cap)
						{
							printf("A problem in All_or_Nothing_Assign() path-trace exceeded %d hops (cycle?) for node seq no %d Orig zone= %d (OD skipped)\n", _hop_cap, CurrentNode, Orig);
#ifdef DUMP_CYCLE
							{
								int cn = Dest;   // dense zone seq == node seq
								printf("  CYCLE Orig=%d Dest=%d chain(node<-predlink->fromnode,cost):", Orig, Dest);
								for (int hh = 0; hh < 50 && cn != internal_node_for_origin_node; hh++)
								{
									int kk = MinPathPredLink[mpred][Orig][cn];
									printf(" n%d<-L%d", cn, kk);
									if (kk <= 0 || kk > number_of_links) { printf("(bad)"); break; }
									printf("(tt=%.3f,from=%d)", Link[kk].Travel_time, Link[kk].internal_from_node_id);
									cn = Link[kk].internal_from_node_id;
								}
								printf("\n");
							}
#endif
							break;
						}
						if (g_has_movement_restrictions)
						{
							// Link-state back-trace: first hop is the best entry link into
							// the destination node; thereafter follow the link-state
							// predecessor chain (NOT the node tree, which is not exact under
							// turn restrictions).
							if (prev_k == INVALID)
								k = MinPathPredLink[mpred][Orig][CurrentNode];
							else
								k = g_MinPathPredLink_link[mpred][Orig][prev_k];
						}
						else
						{
							k = MinPathPredLink[mpred][Orig][CurrentNode];
						}

						if (k <= 0 || k > number_of_links || k == INVALID)
						{
							printf("A problem in All_or_Nothing_Assign() Invalid pred for node seq no %d Orig zone= %d \n\n", CurrentNode, Orig);
							break;
						}
						ProcessorVolume[k][p] += RouteFlow * g_mode_type_vector[m].pce;
						ProcessorModeVolume[k][m][p] += RouteFlow;  // pure volume

						CurrentNode = Link[k].internal_from_node_id;
						prev_k = k;

						if (CurrentNode <= 0 || CurrentNode > no_nodes )
						{
							printf("A problem in All_or_Nothing_Assign() Invalid node seq no %d Orig zone = %d \n\n", CurrentNode, Orig);
							break;
						}

						if(linkIndices.size() >0)
						{
						if (shortest_path_log_flag || Assignment_iteration_no == 0)
						{
#pragma omp critical
							{
								currentLinkSequence.push_back(k); // Store the link index
							}

						}
					}

							if (linkIndices.size() > 0)
							{
								if (shortest_path_log_flag || Assignment_iteration_no == 0)
								{
									AddLinkSequence(m, Orig, Dest, Assignment_iteration_no, currentLinkSequence);
									// Store the link sequence for this OD pair

								}
						}
						}

				}
			}

		}


	}

	if (Assignment_iteration_no == 0)  // with base demand
	{
		for (int k = 1; k <= number_of_links; k++)
		{
			Volume[k] = Link[k].Base_demand_volume;

			for (int p = 0; p < g_number_of_processors; p++)
			{
				Volume[k] += ProcessorVolume[k][p];
			}

			for (int m = 1; m <= number_of_modes; m++)
			{
				Link[k].mode_MainVolume[m] = Link[k].mode_Base_demand_volume[m];

				for (int p = 0; p < g_number_of_processors; p++)
				{
					Link[k].mode_MainVolume[m] += ProcessorModeVolume[k][m][p];
				}
			}

		}

	}
	else
	{  // all or nothing 
#pragma omp parallel for
		for (int k = 1; k <= number_of_links; k++)
		{
			Volume[k] = 0;

			for (int p = 0; p < g_number_of_processors; p++)
			{
				Volume[k] += ProcessorVolume[k][p];
			}

			for (int m = 1; m <= number_of_modes; m++)
			{
				Link[k].mode_SubVolume[m] = 0.0;

				for (int p = 0; p < g_number_of_processors; p++)
				{
					Link[k].mode_SubVolume[m] += ProcessorModeVolume[k][m][p];
				}
			}

		}
	}


	Free_2D((void**)ProcessorVolume, number_of_links, g_number_of_processors);
	Free_3D((void***)ProcessorModeVolume, number_of_links, number_of_modes, g_number_of_processors);

	auto end0 = std::chrono::high_resolution_clock::now();

	// Calculate the duration in milliseconds

	// Calculate the duration in seconds
	auto duration = end0 - start0;

	// Convert to hours, minutes, seconds
	auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
	auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration % std::chrono::hours(1));
	auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration % std::chrono::minutes(1));
	auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(duration % std::chrono::seconds(1));

	printf("All or nothing assignment: %lld hours %lld minutes %lld seconds %lld ms\n", hours.count(), minutes.count(), seconds.count(), milliseconds.count());

}


//---------------------------------------------------------------------------
// ODME via Gradient Descent
//---------------------------------------------------------------------------
std::vector<double> computeTheta(const std::vector<double>& lambda) {
	// Start with the initial path which has a proportion of 1.
	std::vector<double> theta;
	theta.push_back(1.0);

	// For each iteration, update the current theta values and add the new path's proportion.
	for (size_t k = 0; k < lambda.size(); k++) {
		// Scale all current proportions by (1 - lambda[k]).
		for (size_t i = 0; i < theta.size(); i++) {
			theta[i] *= (1.0 - lambda[k]);
		}
		// Append the new path's proportion, which is lambda[k].
		theta.push_back(lambda[k]);
	}

	return theta;
}



void performODME(std::vector<double> theta, double* MainVolume, struct link_record* Link)
{
    // Open log file (this will overwrite any existing file with the same name)
    std::ofstream logFile("ODME_log.txt");
    if (!logFile.is_open()) {
        std::cerr << "Error: Could not open ODME_log.txt for writing." << std::endl;
        return;
    }
    logFile << "Starting ODME process.\n";

    // Define target and weights
    double VMT_target = g_ODME_obs_VMT;
    double w_link = g_ODME_link_volume_penalty;
    double w_od = g_ODME_od_penalty;
    double w_vmt = g_ODME_VMT_penalty;

    // Gradient descent and line search parameters
    double step_size = 0.2;   // initial step size guess
    int maxIter = 400;
    double tol = 1;
    const double armijo_c = 0.05; // Armijo condition constant

    // Get the number of modes and zones (note: index 0 is unused)
    int numModes = number_of_modes;
    int numZones = no_zones;

    // --------------------------------------------------
    // Helper lambda to compute the overall objective function.
    // The objective is the sum of:
    //   - OD demand deviations: w_od * (MDODflow - targetMDODflow)^2,
    //   - Link flow deviations: w_link * (linkFlow - Link[l].Obs_volume)^2,
    //   - VMT deviations: w_vmt * (VMT - VMT_target)^2.
	auto computeObjectiveForMDOD = [&](double*** mdod, double& r2) -> double {
		// Compute link flows from the given OD matrix "mdod"
		std::vector<double> linkFlows(number_of_links + 1, 0.0);
		for (int m = 1; m <= numModes; ++m) {
			for (int Orig = 1; Orig <= numZones; ++Orig) {
				for (int Dest = 1; Dest <= numZones; ++Dest) {
					if (Orig == Dest)
						continue;
					float od_flow = mdod[m][Orig][Dest];
					if (od_flow <= 0)
						continue;
					// Use available routes (assume equal split)
					int nRoutes = std::min(theta.size(), linkIndices[m][Orig][Dest].size());
					for (int route_id = 0; route_id < nRoutes; route_id++) {
						float routeFlow = od_flow * theta[route_id];
						if (linkIndices[m][Orig][Dest][route_id].empty())
							continue;
						for (int link_id : linkIndices[m][Orig][Dest][route_id]) {
							if (link_id >= 0 && link_id < static_cast<int>(linkFlows.size()))
								linkFlows[link_id] += routeFlow;
						}
					}
				}
			}
		}

		// Compute total VMT
		double VMT = 0.0;
		for (size_t l = 1; l < linkFlows.size(); ++l) {
			VMT += linkFlows[l] * Link[l].length;
		}

		// Sum up the objective terms
		double obj = 0.0;

		// OD demand term:
		for (int m = 1; m <= numModes; ++m) {
			for (int Orig = 1; Orig <= numZones; ++Orig) {
				for (int Dest = 1; Dest <= numZones; ++Dest) {
					if (Orig == Dest)
						continue;
					double diff = mdod[m][Orig][Dest] - targetMDODflow[m][Orig][Dest];
					obj += w_od * diff * diff;
				}
			}
		}

		// For R² calculation, first compute the mean of observed link volumes (for valid links)
		double sumObs = 0.0;
		int validCount = 0;
		for (size_t l = 1; l < linkFlows.size(); ++l) {
			if (Link[l].Obs_volume[1] > 1) {
				sumObs += Link[l].Obs_volume[1];
				validCount++;
			}
		}
		double meanObs = (validCount > 0) ? (sumObs / validCount) : 0.0;

		// Compute SSE and SStot for the link flow term
		double SSE_link = 0.0;
		double SStot = 0.0;
		for (size_t l = 1; l < linkFlows.size(); ++l) {
			if (Link[l].Obs_volume[1] > 1) {
				double diff = linkFlows[l] - Link[l].Obs_volume[1];
				SSE_link += w_link * diff * diff;
				double diffMean = Link[l].Obs_volume[1] - meanObs;
				SStot += diffMean * diffMean;
				// Also add to the overall objective function:
				obj += w_link * diff * diff;
			}
		}
		// Calculate R² (coefficient of determination) for link flows
		if (SStot > 0)
			r2 = 1.0 - SSE_link / SStot;
		else
			r2 = 1.0; // if SStot is 0, then the data is constant or perfect

		// VMT term:
		if (VMT_target > 1) {
			double diff = VMT - VMT_target;
			obj += w_vmt * diff * diff;
		}

		return obj;
		};

    // --------------------------------------------------
    // Main gradient descent iterations with line search
	double F_current_best = 0; 


	int count_Armijo_condition_NOT_met = 0; 
	double r2 = 0; 
	for (int iter = 0; iter < maxIter; ++iter) {

		// Compute the current objective value
		double F_current = computeObjectiveForMDOD(MDODflow, r2);
		logFile << "--------------------------------------------------\n";
//		cout << "Iteration " << iter << " R2 = " << r2 << "\n";

		if (iter == 0)
			F_current_best = F_current;
		else
		{
			if (F_current < F_current_best)
				F_current_best = F_current;

		}

		logFile << "Current objective = " << F_current << "\n";
		std::vector<double> linkFlows(number_of_links + 1, 0.0);
		// --- Log simple deviation measures ---

			// OD deviation: compute sum of absolute differences (and average)

		double totalODDev = 0.0;
		int odCount = 0;
		if (g_ODME_target_od >= 1)
		{
			for (int m = 1; m <= numModes; ++m) {
				for (int Orig = 1; Orig <= numZones; ++Orig) {
					for (int Dest = 1; Dest <= numZones; ++Dest) {
						if (Orig == Dest)
							continue;
						double diff = MDODflow[m][Orig][Dest] - targetMDODflow[m][Orig][Dest];
						totalODDev += diff;
						odCount++;
					}
				}
			}
		}
		double avgODDev = (odCount > 0 ? totalODDev / odCount : 0.0);

		// Link deviation: compute link flows and compare with observed volume

		for (int m = 1; m <= numModes; ++m) {
			for (int Orig = 1; Orig <= numZones; ++Orig) {
				for (int Dest = 1; Dest <= numZones; ++Dest) {
					if (Orig == Dest)
						continue;
					float od_flow = MDODflow[m][Orig][Dest];
					if (od_flow <= 0)
						continue;
					int nRoutes = std::min(theta.size(), linkIndices[m][Orig][Dest].size());
					for (int route_id = 0; route_id < nRoutes; route_id++) {
						float routeFlow = od_flow * theta[route_id];
						if (linkIndices[m][Orig][Dest][route_id].empty())
							continue;
						for (int link_id : linkIndices[m][Orig][Dest][route_id]) {
							if (link_id >= 0 && link_id < static_cast<int>(linkFlows.size()))
								linkFlows[link_id] += routeFlow;
						}
					}
				}
			}
		}
		double totalLinkDev = 0.0;
		int linkCount = 0;
		for (size_t l = 0; l < linkFlows.size(); ++l) {
			if (Link[l].Obs_volume[1] > 1) {
				totalLinkDev += linkFlows[l] - Link[l].Obs_volume[1];
				linkCount++;
			}
		}
		double avgLinkDev = (linkCount > 0 ? totalLinkDev / linkCount : 0.0);

		// VMT deviation: compute total VMT and its deviation from target
		double VMT = 0.0;
		for (size_t l = 0; l < linkFlows.size(); ++l) {
			VMT += linkFlows[l] * Link[l].length;
		}
		double vmtDev = (VMT_target > 0 ? (VMT - VMT_target) : 0.0);

		cout << "iter = " << iter << " [Deviation Log] OD avg deviation = " << avgODDev
			<< ", Link avg deviation = " << avgLinkDev
			<< ", VMT deviation = " << vmtDev << " ";

		logFile << "iter = " << iter << " [Deviation Log] OD avg deviation = " << avgODDev
			<< ", Link avg deviation = " << avgLinkDev
			<< ", VMT deviation = " << vmtDev << "\n";



		double sumGradSq = 0.0;
		// --- Compute gradients for each OD pair ---
		for (int m = 1; m <= numModes; ++m) {
			for (int Orig = 1; Orig <= numZones; ++Orig) {
				for (int Dest = 1; Dest <= numZones; ++Dest) {
					if (Orig == Dest)
						continue;


					float od_flow = MDODflow[m][Orig][Dest];
					if (od_flow <= 0)
						continue;

					gradient_MDODflow[m][Orig][Dest] = 0;

					double grad_link = 0.0;
					double grad_vmt = 0.0;
					// For link contributions, count usage of each link for this OD pair.
					int routeCount = 0;
					std::unordered_map<int, float> linkUsage;
					for (int route_id = 0; route_id < static_cast<int>(linkIndices[m][Orig][Dest].size()); ++route_id) {
						if (linkIndices[m][Orig][Dest][route_id].empty())
							continue;
						routeCount++;
						for (int link_id : linkIndices[m][Orig][Dest][route_id]) {

							if (Link[link_id].Obs_volume[1] > 1) {
								linkUsage[link_id] += 1;
							}
						}
					}
					if (routeCount == 0)
						continue;

					// In this example the link and VMT gradient terms are approximated.
					// (Depending on your formulation, you might want to recompute the sensitivities using the updated link flows.)
					for (auto& entry : linkUsage) {
						int link_id = entry.first;
						int count_usage = entry.second;
						double sensitivity = static_cast<double>(count_usage) / routeCount;
						if (Link[link_id].Obs_volume[1] > 1) {
							grad_link += 2.0 * w_link * (linkFlows[link_id] - Link[link_id].Obs_volume[1]) * sensitivity;
						}
						if (VMT_target > 1) {
							grad_vmt += 2.0 * w_vmt * (Link[link_id].length * sensitivity);
						}
					}
					double grad_od = 0;
					
					if(targetMDODflow[m][Orig][Dest]>=0.0)
						grad_od = 2.0 * w_od * (od_flow - targetMDODflow[m][Orig][Dest]);

					double total_grad = grad_link + grad_vmt + grad_od;
					gradient_MDODflow[m][Orig][Dest] = total_grad;
					sumGradSq += total_grad * total_grad;
				}
			}
		}

		// Check for convergence based on the norm of the gradient.
		if (std::sqrt(sumGradSq) < tol) {
			logFile << "Convergence reached after " << iter << " iterations.\n";
			std::cout << "Convergence reached after " << iter << " iterations.\n";
			break;
		}

		// --- Backtracking Line Search ---
		double alpha = step_size;
		// Save a copy of the current MDODflow.
		for (int m = 1; m <= numModes; ++m)
			for (int Orig = 1; Orig <= numZones; ++Orig)
				for (int Dest = 1; Dest <= numZones; ++Dest)
				{
					old_MDODflow[m][Orig][Dest] = MDODflow[m][Orig][Dest];
					candidate_MDODflow[m][Orig][Dest] = MDODflow[m][Orig][Dest];

				}
		// Prepare a candidate OD matrix.

		int line_search_iteration = 0; 

		while (true) {

		//	cout << "Line search iteration " << line_search_iteration << ": Step size chosen = " << alpha << endl; 
			line_search_iteration++; 
			// Update candidate MDODflow: candidate = oldMDODflow - alpha * gradient.
			for (int m = 1; m <= numModes; ++m) 
				for (int Orig = 1; Orig <= numZones; ++Orig) 
					for (int Dest = 1; Dest <= numZones; ++Dest) 
					{
						if (Orig == Dest)
							continue;

						if (old_MDODflow[m][Orig][Dest] > 0.001)
						{
							candidate_MDODflow[m][Orig][Dest] = old_MDODflow[m][Orig][Dest] - alpha * gradient_MDODflow[m][Orig][Dest];

							double larger_value = fmax(seed_MDODflow[m][Orig][Dest], targetMDODflow[m][Orig][Dest]);

							double lesser_value = seed_MDODflow[m][Orig][Dest];
								if(targetMDODflow[m][Orig][Dest] >=0.0)
									lesser_value = fmin(seed_MDODflow[m][Orig][Dest], targetMDODflow[m][Orig][Dest]);

							if (candidate_MDODflow[m][Orig][Dest] < lesser_value *0.5)
								candidate_MDODflow[m][Orig][Dest] = lesser_value * 0.5; // ensure non-negativity and resonableness

							if (candidate_MDODflow[m][Orig][Dest] > larger_value* 1.5)
								candidate_MDODflow[m][Orig][Dest] = larger_value * 1.5; // ensure uppper bound resonableness

						}
					}
				




				double F_candidate = computeObjectiveForMDOD(candidate_MDODflow, r2);
				logFile << "[Line Search] alpha = " << alpha
					<< " | F_candidate = " << F_candidate
					<< " | F_current = " << F_current
					<< " | Threshold = " << (F_current - armijo_c * alpha * sumGradSq) << "\n";
				// Check the Armijo condition:
				if (F_candidate <= F_current - armijo_c * alpha * sumGradSq) {
					//cout << "[Line Search] Armijo condition met: F_candidate (" << F_candidate
					//	<< ") <= F_current (" << F_current << ") - "
					//	<< armijo_c << " * alpha (" << alpha << ") * sumGradSq (" << sumGradSq << ")\n";

					logFile << "[Line Search] Armijo condition met: F_candidate (" << F_candidate
						<< ") <= F_current (" << F_current << ") - "
						<< armijo_c << " * alpha (" << alpha << ") * sumGradSq (" << sumGradSq << ")\n";
					count_Armijo_condition_NOT_met = 0;
					break;
				}
				else {
					logFile << "[Line Search] Armijo condition NOT met. Reducing alpha from " << alpha;
					alpha *= 0.5; // reduce step size
					logFile << " to " << alpha << "\n";
					if (alpha < 1e-8) {  // minimal threshold to avoid endless looping
						count_Armijo_condition_NOT_met++;
						cout << "[Line Search] Alpha dropped below threshold (1e-8). Exiting line search loop.\n";
						logFile << "[Line Search] Alpha dropped below threshold (1e-8). Exiting line search loop.\n";						break;
					}
				}
			
			logFile << "Iteration " << iter << ": Step size chosen = " << alpha
				<< ", Objective: " << F_current << " -> " << computeObjectiveForMDOD(candidate_MDODflow, r2) << "\n";
			cout << "Iteration " << iter << ": Step size chosen = " << alpha
				<< ", Objective: " << F_current << " -> " << computeObjectiveForMDOD(candidate_MDODflow, r2) << "\n";
			// Update MDODflow with the candidate from the line search.

			if (count_Armijo_condition_NOT_met >= 3)
			{
				cout << "count_Armijo_condition_NOT_met >= 3. Reaching final estimate\n";
				break;
			}
		}
			for (int m = 1; m <= numModes; ++m)
				for (int Orig = 1; Orig <= numZones; ++Orig)
					for (int Dest = 1; Dest <= numZones; ++Dest)
					{
						if (Orig == Dest)
							continue;
						else
						{
							MDODflow[m][Orig][Dest] = candidate_MDODflow[m][Orig][Dest];
						}
					}

			if (iter >= 5)
			{
				if (F_current > F_current_best + 1 || F_current > F_current_best * 1.1)
				{

					logFile << "Convergence reached as F_current > F_current_best\n";
					std::cout << "Convergence reached as F_current > F_current_best\n";
					break;
				}
			}
			// end of gradient descent loop

		
	}
    // --------------------------------------------------
    // Final synchronization of path flow and link volume 
    std::vector<double> linkFlows(number_of_links + 1, 0.0);
    std::vector<std::vector<double>> mode_linkFlows(numModes + 1, std::vector<double>(number_of_links + 1, 0.0));
    for (int m = 1; m <= numModes; ++m) {
        for (int Orig = 1; Orig <= numZones; ++Orig) {
            for (int Dest = 1; Dest <= numZones; ++Dest) {
                if (Orig == Dest)
                    continue;
                float od_flow = MDODflow[m][Orig][Dest];
                if (od_flow <= 0)
                    continue;
                int nRoutes = std::min(theta.size(), linkIndices[m][Orig][Dest].size());
                for (int route_id = 0; route_id < nRoutes; route_id++) {
                    float routeFlow = od_flow * theta[route_id];
                    if (linkIndices[m][Orig][Dest][route_id].empty())
                        continue;
                    for (int link_id : linkIndices[m][Orig][Dest][route_id]) {
                        if (link_id >= 0 && link_id < static_cast<int>(linkFlows.size())) {
                            linkFlows[link_id] += routeFlow;
                            mode_linkFlows[m][link_id] += routeFlow;
                        }
                    }
                }
            }
        }
        for (size_t l = 0; l < linkFlows.size(); ++l) 
		{
            MainVolume[l] = linkFlows[l];
            Link[l].mode_MainVolume[m] = mode_linkFlows[m][l];
        }



		// Start the R² calculation process for link flows.
		logFile << "Starting R² calculation for link flows.\n";

		// ---------------------------------------------------------------------
		// First pass: compute mean of observed volumes (for valid links)
		// ---------------------------------------------------------------------
		double sumObs = 0.0;
		int validCount = 0;
		for (size_t l = 1; l < linkFlows.size(); ++l) {
			if (Link[l].Obs_volume[1] > 1) {
				logFile << "Link " << l << " is valid. Observed volume: " << Link[l].Obs_volume[1] << "\n";
				sumObs += Link[l].Obs_volume[1];
				validCount++;
			}
			else {
				logFile << "Link " << l << " is invalid. Observed volume: " << Link[l].Obs_volume[1] << "\n";
			}
		}
		double meanObs = (validCount > 0) ? (sumObs / validCount) : 0.0;
		logFile << "Total valid links: " << validCount
			<< ", Sum of observed volumes: " << sumObs
			<< ", Mean observed volume: " << meanObs << "\n";

		// ---------------------------------------------------------------------
		// Second pass: Compute SSE and SStot for regression R² and accumulate
		// values for Pearson correlation computation.
		// ---------------------------------------------------------------------
		double SSE = 0.0;
		double SStot = 0.0;

		// Variables for Pearson correlation
		double sumObs_corr = 0.0;
		double sumPred_corr = 0.0;
		double sumObsPred = 0.0;
		double sumObsSq = 0.0;
		double sumPredSq = 0.0;
		int count = 0;

		for (size_t l = 1; l < linkFlows.size(); ++l) {
			if (Link[l].Obs_volume[1] > 1) {
				double obs = Link[l].Obs_volume[1];
				double pred = linkFlows[l];
				double diff = pred - obs;
				double squaredDiff = diff * diff;
				SSE += squaredDiff;
				double diffMean = obs - meanObs;
				double squaredDiffMean = diffMean * diffMean;
				SStot += squaredDiffMean;

				// Accumulate for correlation computation.
				count++;
				sumObs_corr += obs;
				sumPred_corr += pred;
				sumObsPred += obs * pred;
				sumObsSq += obs * obs;
				sumPredSq += pred * pred;

				logFile << "Link " << l << ": Predicted volume: " << pred
					<< ", Observed volume: " << obs
					<< ", Difference: " << diff
					<< ", Squared diff: " << squaredDiff
					<< ", Squared diff from mean: " << squaredDiffMean << "\n";
			}
		}
		logFile << "Total SSE (link): " << SSE << ", Total SStot: " << SStot << "\n";

		// ---------------------------------------------------------------------
		// Compute R² (regression) = 1 - SSE/SStot
		// ---------------------------------------------------------------------
		double r2_regression = 1.0;
		if (SStot > 0) {
			r2_regression = 1.0 - SSE / SStot;
			logFile << "Computed R2 (regression): " << r2_regression << " (1 - SSE/SStot).\n";
		}
		else {
			logFile << "SStot is 0. Data is constant or perfect. Set R2 (regression) to 1.0.\n";
		}

		// ---------------------------------------------------------------------
		// Compute Pearson correlation coefficient and then R2 (from correlation)
		// ---------------------------------------------------------------------
		double numerator = count * sumObsPred - sumObs_corr * sumPred_corr;
		double denominator_term1 = count * sumObsSq - sumObs_corr * sumObs_corr;
		double denominator_term2 = count * sumPredSq - sumPred_corr * sumPred_corr;
		double r_corr = 0.0;
		if (denominator_term1 > 0 && denominator_term2 > 0) {
			double denominator = std::sqrt(denominator_term1 * denominator_term2);
			if (denominator != 0) {
				r_corr = numerator / denominator;
			}
		}
		double r2_correlation = r_corr * r_corr;
		logFile << "Computed Pearson correlation: " << r_corr
			<< ", R2 (from correlation): " << r2_correlation << "\n";

		logFile << "ODME process finished.\n";

		// ---------------------------------------------------------------------
		// Output the results
		// ---------------------------------------------------------------------
		std::cout << "R2 (regression) = " << r2_regression << std::endl;
		std::cout << "R2 (from correlation) = " << r2_correlation << std::endl;

		logFile << "ODME process finished.\n";
    }


    logFile << "ODME process finished.\n";

    logFile.close();
}

// Define a structure to hold all relevant route data.
struct RouteData {
	int firstRouteID;                // The first encountered route_id (optional)
	int unique_route_id;             // A sequential unique id for output
	std::string nodeIDsStr;          // Node IDs (string)
	std::string linkIDsStr;          // Link IDs (string)
	std::string linkLengthStr;          // Link Length (string)
	double totalDistance;
	double totalFreeFlowTravelTime;
	double totalTravelTime;
	double accumulatedTheta = 0;         // Sum of theta for all routes with this key
	std::string routeKey;            // Unique key (e.g., based on node and link sums)
};

void OutputRouteDetails(const std::string& filename, std::vector<double> theta)
{
	std::ofstream outputFile(filename);  // Open the file for writing

	if (linkIndices.size() == 0)
		return; 
	// Write the CSV header in lowercase
	outputFile << "mode,route_id,o_zone_id,d_zone_id,unique_route_id,prob,node_ids,link_ids,distance_mile,total_distance_km,total_free_flow_travel_time,total_travel_time,route_key,seed_od_volume,target_od_volume,final_est_od_volume,volume,";

	outputFile << "\n";


	for (int m = 1; m < linkIndices.size(); ++m)
	{
		for (int Orig = 1; Orig < linkIndices[m].size(); ++Orig)
		{
			for (int Dest = 1; Dest < linkIndices[m][Orig].size(); ++Dest)
			{
				std::unordered_map<std::string, RouteData> uniqueRoutes;

				// Reset unique route counter.
				int unique_route_id = 1;

				// Loop over all candidate routes.
				for (int route_id = 0; route_id < linkIndices[m][Orig][Dest].size(); ++route_id)
				{
					if (!linkIndices[m][Orig][Dest][route_id].empty())
					{
						double totalDistance = 0.0;
						double totalFreeFlowTravelTime = 0.0;
						double totalTravelTime = 0.0;
						std::string nodeIDsStr;
						std::string linkIDsStr;
						std::string linkLengthStr;
						
						int nodeSum = 0;  // Sum of node IDs (for uniqueness key)
						int linkSum = 0;  // Sum of link IDs (for uniqueness key)

						// Process the route in reverse order (as in your code).
						for (int i = linkIndices[m][Orig][Dest][route_id].size() - 1; i >= 0; --i)
						{
							int k = linkIndices[m][Orig][Dest][route_id][i];

							// Append the from_node_id and accumulate its value.
							int fromNodeID = Link[k].external_from_node_id;
							nodeIDsStr += std::to_string(fromNodeID) + ";";
							nodeSum += fromNodeID;

							// Append the external (link.csv) link id; keep the internal
							// index k only for the uniqueness checksum below.
							linkIDsStr += std::to_string(Link[k].external_link_id) + ";";
							//double length = Link[k].length;
							//linkLengthStr += std::to_string(length) + ";";
							linkSum += k;

							// Sum the distance and travel times.
							totalDistance += Link[k].length;
							totalFreeFlowTravelTime += Link[k].FreeTravelTime;
							totalTravelTime += Link[k].Travel_time;

							// For the first link in the (reversed) order, also add the to_node_id.
							if (i == 0)
							{
								int toNodeID = Link[k].external_to_node_id;
								nodeIDsStr += std::to_string(toNodeID);
								nodeSum += toNodeID;
							}
						}

						// Create a unique key for the route.
						std::string routeKey = std::to_string(nodeSum) + "_" + std::to_string(linkSum);

						// Check if this unique route already exists.
						auto it = uniqueRoutes.find(routeKey);
						if (it == uniqueRoutes.end())
						{
							// Not found: create a new RouteData entry.
							RouteData rd;
							rd.firstRouteID = route_id;      // Save the first encountered candidate's route_id.
							rd.unique_route_id = unique_route_id;  // Unique id for output.
							rd.nodeIDsStr = nodeIDsStr;
							rd.linkIDsStr = linkIDsStr;
							//rd.linkLengthStr = linkLengthStr; 
							rd.totalDistance = totalDistance;
							rd.totalFreeFlowTravelTime = totalFreeFlowTravelTime;
							rd.totalTravelTime = totalTravelTime;

							if(route_id < theta.size() )
							{
							rd.accumulatedTheta = theta[route_id];  // Initialize with current theta.
							}
							else
							{
								rd.accumulatedTheta = unique_route_id *1.0/uniqueRoutes.size();

							}
							rd.routeKey = routeKey;

							// Insert the new unique route.
							uniqueRoutes[routeKey] = rd;
							unique_route_id++;
						}
						else
						{
							// Duplicate route: accumulate theta.
							it->second.accumulatedTheta += theta[route_id];
							// Optionally, you could verify that other route details (e.g., distances) match.
						}
					}
				}

				// Now output the unique routes with their accumulated theta values.
				for (const auto& pair : uniqueRoutes)
				{
					const RouteData& rd = pair.second;
					float od_volume = 1;
					float seed_od_volume = 1;
					float target_od_volume = 0;

					float route_volume = 1;

					double accumulatedTheta = rd.accumulatedTheta;
					if (TotalAssignIterations <= 1)
						accumulatedTheta = 1.0;

					if(MDODflow!=NULL)
					{
						seed_od_volume = seed_MDODflow[m][Orig][Dest];
						od_volume = MDODflow[m][Orig][Dest];
						target_od_volume = targetMDODflow[m][Orig][Dest];
						route_volume = od_volume * accumulatedTheta;
					} 



					// Route-output volume floor: on large-zone networks (>=1000)
					// only OD pairs with volume >= floor are written (file-size
					// guard). Overridable via env TAPLITE_ROUTE_VOL_MIN so
					// matrix-operator extraction (OD->path->link) can request
					// full coverage. Default 1.0 preserves prior behavior.
					static double route_vol_min = -1.0;
					if (route_vol_min < 0.0)
					{
						const char* ev = getenv("TAPLITE_ROUTE_VOL_MIN");
						route_vol_min = (ev != NULL) ? atof(ev) : 1.0;
					}
					if(no_zones <1000 || (no_zones>=1000 && od_volume >= route_vol_min))
					{
					// (Optional) Remove trailing semicolon from linkIDsStr if needed.
					std::string cleanedLinkIDsStr = rd.linkIDsStr;
					if (!cleanedLinkIDsStr.empty() && cleanedLinkIDsStr.back() == ';')
						cleanedLinkIDsStr.pop_back();

					outputFile << g_mode_type_vector[m].mode_type.c_str() << ","
						<< rd.firstRouteID << ","  // or route_id (the first candidate id)
						<< g_zone_seq_2_old_zone_id[Orig] << ","
						<< g_zone_seq_2_old_zone_id[Dest] << ","
						<< rd.unique_route_id << ","
						<< accumulatedTheta << ","
						<< rd.nodeIDsStr << ","
						<< cleanedLinkIDsStr << ","
						//<< rd.linkLengthStr << ","
						<< rd.totalDistance << ","
						<< rd.totalDistance * 1.609 << ","
						<< rd.totalFreeFlowTravelTime << ","
						<< rd.totalTravelTime << ","
						<< rd.routeKey << ","
						<< seed_od_volume << "," 
						<< target_od_volume << "," 
						<<  od_volume << ","

						<< route_volume << "\n";
					}
					else
					{

						// Process the route in reverse order (as in your code).
						for (int i = linkIndices[m][Orig][Dest][rd.firstRouteID].size() - 1; i >= 0; --i)
						{
							int k = linkIndices[m][Orig][Dest][rd.firstRouteID][i];
							Link[k].background_volume += route_volume;
						}
					}
				}

			}
		}
	}

	// Close the file after writing
	outputFile.close();


	std::cout << "Output written to " << filename << std::endl;

}




constexpr auto LCG_a = 17364;
constexpr auto LCG_c = 0;
constexpr auto LCG_M = 65521; // 16-bit memory-saving LCG
unsigned int RandomSeed = 12345; // Initial random seed

// Function to generate random numbers using LCG
float generate_random() {
	RandomSeed = (LCG_a * RandomSeed + LCG_c) % LCG_M;
	return float(RandomSeed) / LCG_M; // Normalize to [0, 1)
}


// Helper function to convert minutes to hh:mm:ss format
std::string convert_minutes_to_hhmmss(double minutes) {
	int total_seconds = static_cast<int>(minutes * 60); // Convert minutes to total seconds
	int hours = total_seconds / 3600;                   // Extract hours
	int remaining_seconds = total_seconds % 3600;       // Remaining seconds after hours
	int mins = remaining_seconds / 60;                  // Extract minutes
	int secs = remaining_seconds % 60;                  // Extract seconds

	// Format as hh:mm:ss
	std::ostringstream oss;
	oss << std::setw(2) << std::setfill('0') << hours << ":"
		<< std::setw(2) << std::setfill('0') << mins << ":"
		<< std::setw(2) << std::setfill('0') << secs;

	return oss.str();
}


// Additional functions (find_links_to, init_link_pointers, init_links) would follow
// with similar C++ implementations...
// Calculate departure time function
std::vector<double> calculate_departure_time(double T0, double T1, double T3, double beta, int route_volume) {
	std::vector<double> departure_times;
	double time_span = T3 - T0;

	for (int i = 0; i < route_volume; i++) {
		double normalized_time = static_cast<double>(i) / route_volume;
		double time = T0 + (beta * pow(normalized_time, 2) + normalized_time) * time_span;
		departure_times.push_back(time);
	}

	return departure_times;
}
void OutputVehicleDetails(const std::string& filename, std::vector<double> theta)
{
	// Record the start time
	auto start = std::chrono::high_resolution_clock::now();

	std::ofstream outputFile(filename);  // Open the file for writing

	if (linkIndices.size() == 0)
		return;
	// Write the CSV header in lowercase
	outputFile << "agent_id,departure_time,departure_time_hhmmss,mode,route_id,o_zone_id,d_zone_id,unique_route_id,node_ids,link_ids,total_distance_mile,total_distance_km,total_free_flow_travel_time,total_travel_time,route_key,route_volume,\n";

	long agent_id = 1; 

	for (int m = 1; m < linkIndices.size(); ++m)
	{
		for (int Orig = 1; Orig < linkIndices[m].size(); ++Orig)
		{
			for (int Dest = 1; Dest < linkIndices[m][Orig].size(); ++Dest)
			{
				std::unordered_map<std::string, RouteData> uniqueRoutes;

				// Reset unique route counter.
				int unique_route_id = 1;

				// Loop over all candidate routes.
				for (int route_id = 0; route_id < linkIndices[m][Orig][Dest].size(); ++route_id)
				{
					if (!linkIndices[m][Orig][Dest][route_id].empty())
					{
						double totalDistance = 0.0;
						double totalFreeFlowTravelTime = 0.0;
						double totalTravelTime = 0.0;
						std::string nodeIDsStr;
						std::string linkIDsStr;

						int nodeSum = 0;  // Sum of node IDs (for uniqueness key)
						int linkSum = 0;  // Sum of link IDs (for uniqueness key)

						// Process the route in reverse order (as in your code).
						for (int i = linkIndices[m][Orig][Dest][route_id].size() - 1; i >= 0; --i)
						{
							int k = linkIndices[m][Orig][Dest][route_id][i];

							// Append the from_node_id and accumulate its value.
							int fromNodeID = Link[k].external_from_node_id;
							nodeIDsStr += std::to_string(fromNodeID) + ";";
							nodeSum += fromNodeID;

							// Append the external (link.csv) link id; keep the internal
							// index k only for the uniqueness checksum below.
							linkIDsStr += std::to_string(Link[k].external_link_id) + ";";
							linkSum += k;

							// Sum the distance and travel times.
							totalDistance += Link[k].length;
							totalFreeFlowTravelTime += Link[k].FreeTravelTime;
							totalTravelTime += Link[k].Travel_time;

							// For the first link in the (reversed) order, also add the to_node_id.
							if (i == 0)
							{
								int toNodeID = Link[k].external_to_node_id;
								nodeIDsStr += std::to_string(toNodeID);
								nodeSum += toNodeID;
							}
						}

						// Create a unique key for the route.
						std::string routeKey = std::to_string(nodeSum) + "_" + std::to_string(linkSum);

						// Check if this unique route already exists.
						auto it = uniqueRoutes.find(routeKey);
						if (it == uniqueRoutes.end())
						{
							// Not found: create a new RouteData entry.
							RouteData rd;
							rd.firstRouteID = route_id;      // Save the first encountered candidate's route_id.
							rd.unique_route_id = unique_route_id;  // Unique id for output.
							rd.nodeIDsStr = nodeIDsStr;
							rd.linkIDsStr = linkIDsStr;
							rd.totalDistance = totalDistance;
							rd.totalFreeFlowTravelTime = totalFreeFlowTravelTime;
							rd.totalTravelTime = totalTravelTime;
							rd.accumulatedTheta = theta[route_id];  // Initialize with current theta.
							rd.routeKey = routeKey;

							// Insert the new unique route.
							uniqueRoutes[routeKey] = rd;
							unique_route_id++;
						}
						else
						{
							// Duplicate route: accumulate theta.
							it->second.accumulatedTheta += theta[route_id];
							// Optionally, you could verify that other route details (e.g., distances) match.
						}
					}
				}

				// Now output the unique routes with their accumulated theta values.
				for (const auto& pair : uniqueRoutes)
				{
					const RouteData& rd = pair.second;
					float od_volume = MDODflow[m][Orig][Dest];
					float route_volume = od_volume * rd.accumulatedTheta;

					if (!linkIndices[m][Orig][Dest][rd.firstRouteID].empty())
					{
						double totalDistance = 0.0;
						double totalFreeFlowTravelTime = 0.0;
						double totalTravelTime = 0.0;
						std::string nodeIDsStr;
						std::string linkIDsStr;

						std::vector<double> linkTravelTimes; // Store travel times for each link
						// Record the travel time for this link
							for (int i = linkIndices[m][Orig][Dest][rd.firstRouteID].size() - 1; i >= 0; --i)
						{
							int k = linkIndices[m][Orig][Dest][rd.firstRouteID][i];
							linkTravelTimes.push_back(Link[k].Travel_time);
						}

						// Remove trailing space from the link IDs string
						if (!linkIDsStr.empty())
							linkIDsStr.pop_back();

						// Use the route_volume computed above (od_volume * accumulatedTheta).
						// Previously this block re-declared od_volume/route_volume and reset
						// route_volume to 0, so integer_volume was always 0 and no vehicles
						// were ever generated.
						int integer_volume = static_cast<int>(route_volume);
						float residual = route_volume - integer_volume;

						// Add residual probability
						if (generate_random() < residual) {
							++integer_volume;
						}
						int total_time_window_min = (demand_period_ending_hours - demand_period_starting_hours) * 60;

						float departure_time = demand_period_starting_hours * 60;


						if (integer_volume > 0)
						{


							float route_time_increment = static_cast<float>(total_time_window_min) / integer_volume;


							for (int v = 0; v < integer_volume; ++v)
							{
								float time_increment_in_min = route_time_increment;
								if (integer_volume > 10)
								{
									departure_time = demand_period_starting_hours * 60 + v * route_time_increment;
									time_increment_in_min = route_time_increment;
								}
								else
								{
									departure_time = demand_period_starting_hours * 60 + agent_id % total_time_window_min;
									time_increment_in_min = 1;
								}

								std::vector<double> nodeArrivalTimes; // Store arrival times for each node
								double cumulativeTime = departure_time;
								// Calculate arrival times for each node
								for (double travelTime : linkTravelTimes) {
									cumulativeTime += travelTime;
									nodeArrivalTimes.push_back(cumulativeTime);
								}

								std::string timeSequenceStr;
								// Convert arrival times to a semicolon-separated string
								timeSequenceStr.clear();
								for (double time : nodeArrivalTimes) {
									timeSequenceStr += convert_minutes_to_hhmmss(time) + ";";
								}
								if (!timeSequenceStr.empty()) timeSequenceStr.pop_back();

								//departure_time += std::rand() * time_increment_in_min;

								std::string departure_time_hhmmss = convert_minutes_to_hhmmss(departure_time);

								// Write the data for this OD pair and route to the CSV file
								// Use the route geometry/metrics already stored on the unique
								// route (rd). The local nodeIDsStr/linkIDsStr/totalDistance/
								// totalFreeFlowTravelTime/totalTravelTime in this block were
								// never populated, so these columns came out empty/zero.
								outputFile << agent_id << "," << departure_time << ","
									<< departure_time_hhmmss << ","
									<< g_mode_type_vector[m].mode_type.c_str() << ","
									<< rd.firstRouteID << "," << g_zone_seq_2_old_zone_id[Orig] << "," << g_zone_seq_2_old_zone_id[Dest] << "," << unique_route_id << ","
									<< rd.nodeIDsStr << "," << rd.linkIDsStr << ","
									<< rd.totalDistance << "," << rd.totalDistance * 1.609 << "," << rd.totalFreeFlowTravelTime << ","
									<< rd.totalTravelTime << "," << rd.firstRouteID << "," << route_volume << "\n";

								agent_id++;
							}

							unique_route_id++;
						}
						//else
						//{
						//    // Duplicate path found, skipping output
						//    //std::cout << "Duplicate route skipped for Origin: " << Orig << ", Destination: " << Dest << "\n";
						//}
					}
				}

			}

		}
	}
	// Close the file after writing
	outputFile.close();


	std::cout << "Output written to " << filename << std::endl;
	// Record the end time

}


void OutputODPerformance(const std::string& filename)
{
	std::ofstream outputFile(filename);
	std::ofstream googleMapsFile("google_maps_od_distance.csv");

	if (!googleMapsFile.is_open())
	{
		std::cerr << "Error: Could not open Google Maps OD distance file." << std::endl;
		return;
	}

	googleMapsFile << "route_key,mode,o_zone_id,d_zone_id,volume,total_distance_mile,total_distance_km,straight_line_distance_mile,straight_line_distance_km,distance_ratio,total_free_flow_travel_time,total_congestion_travel_time,google_maps_http_link,WKT_geometry\n";

	outputFile << "route_key,mode,o_zone_id,d_zone_id,o_x_coord,o_y_coord,d_x_coord,d_y_coord,total_distance_mile,total_distance_km,straight_line_distance_mile,straight_line_distance_km,distance_ratio,total_free_flow_travel_time,total_congestion_travel_time,volume,WKT_geometry\n";

	double grand_totalDistance = 0.0;
	double grand_totalFreeFlowTravelTime = 0.0;
	double grand_totalTravelTime = 0.0;
	double grand_total_count = 0;
	double grand_total_straight_line_distance = 0.0;
	double grand_total_distance_ratio = 0.0;

	for (int m = 1; m < linkIndices.size(); ++m)
	{
		for (int Orig = 1; Orig < linkIndices[m].size(); ++Orig)
		{
			for (int Dest = 1; Dest < linkIndices[m][Orig].size(); ++Dest)
			{
				std::unordered_map<std::string, bool> uniqueRoutes;
				for (int route_id = 0; route_id < linkIndices[m][Orig][Dest].size(); ++route_id)
				{
					if (!linkIndices[m][Orig][Dest][route_id].empty())
					{
						double totalDistance = 0.0;
						double totalFreeFlowTravelTime = 0.0;
						double totalTravelTime = 0.0;

						std::string nodeIDsStr;
						std::string linkIDsStr;

						long nodeSum = 0;
						long linkSum = 0;


						for (int i = linkIndices[m][Orig][Dest][route_id].size() - 1; i >= 0; --i)
						{
							long k = linkIndices[m][Orig][Dest][route_id][i];

							int fromNodeID = Link[k].external_from_node_id;
							int internal_fromNodeID = Link[k].internal_from_node_id;
							double from_x = g_node_vector[internal_fromNodeID].x;
							double from_y = g_node_vector[internal_fromNodeID].y;

							nodeIDsStr += std::to_string(fromNodeID) + ";";
							nodeSum += fromNodeID;

							linkIDsStr += std::to_string(k) + ";";
							linkSum += k;

							totalDistance += Link[k].length;
							totalFreeFlowTravelTime += Link[k].FreeTravelTime;
							totalTravelTime += Link[k].Travel_time;


							if (i == 0)
							{
								int toNodeID = Link[k].external_to_node_id;
								int internal_toNodeID = Link[k].internal_to_node_id;
								double to_x = g_node_vector[internal_toNodeID].x;
								double to_y = g_node_vector[internal_toNodeID].y;

								nodeIDsStr += std::to_string(toNodeID);
								nodeSum += toNodeID;
							}
						}

						std::string WKT_geometry = "LINESTRING(";
						if (no_zones < 300 || (no_zones >= 300 && MDODflow[m][Orig][Dest] >= 10.0))
						{

							// Construct WKT format from node coordinates
							for (int i = linkIndices[m][Orig][Dest][route_id].size() - 1; i >= 0; --i)
							{
								long k = linkIndices[m][Orig][Dest][route_id][i];

								int fromNodeID = Link[k].external_from_node_id;
								int internal_fromNodeID = Link[k].internal_from_node_id;
								double from_x = g_node_vector[internal_fromNodeID].x;
								double from_y = g_node_vector[internal_fromNodeID].y;



								WKT_geometry += std::to_string(from_x) + " " + std::to_string(from_y) + ",";

								if (i == 0)
								{
									int toNodeID = Link[k].external_to_node_id;
									int internal_toNodeID = Link[k].internal_to_node_id;
									double to_x = g_node_vector[internal_toNodeID].x;
									double to_y = g_node_vector[internal_toNodeID].y;

									WKT_geometry += std::to_string(to_x) + " " + std::to_string(to_y);
								}
							}


						}
						WKT_geometry += ")"; // Close the LINESTRING

						std::string routeKey = std::to_string(Orig) + "_" + std::to_string(Dest) + "_" + g_mode_type_vector[m].mode_type;

						if (uniqueRoutes.find(routeKey) == uniqueRoutes.end())
						{
							uniqueRoutes[routeKey] = true;
							float volume = MDODflow[m][Orig][Dest];

							grand_totalDistance += totalDistance * volume;
							grand_totalFreeFlowTravelTime += totalFreeFlowTravelTime * volume;
							grand_totalTravelTime += totalTravelTime * volume;
							grand_total_count += volume;

							int internal_o_node_id = Orig;   // dense zone seq == node seq
							int internal_d_node_id = Dest;   // dense zone seq == node seq
							double o_x_coord = g_node_vector[internal_o_node_id].x;
							double o_y_coord = g_node_vector[internal_o_node_id].y;
							double d_x_coord = g_node_vector[internal_d_node_id].x;
							double d_y_coord = g_node_vector[internal_d_node_id].y;

							double straight_line_distance_mile = 3958.8 * 2 * atan2(sqrt(pow(sin((d_y_coord - o_y_coord) * M_PI / 360.0), 2) + cos(o_y_coord * M_PI / 180.0) * cos(d_y_coord * M_PI / 180.0) * pow(sin((d_x_coord - o_x_coord) * M_PI / 360.0), 2)), sqrt(1 - pow(sin((d_y_coord - o_y_coord) * M_PI / 360.0), 2) - cos(o_y_coord * M_PI / 180.0) * cos(d_y_coord * M_PI / 180.0) * pow(sin((d_x_coord - o_x_coord) * M_PI / 360.0), 2)));
							double straight_line_distance_km = straight_line_distance_mile * 1.609;

							double distance_ratio = (straight_line_distance_mile > 0.001) ? totalDistance / straight_line_distance_mile : 1.0;

							grand_total_straight_line_distance += straight_line_distance_mile * volume;
							grand_total_distance_ratio += distance_ratio * volume;

							outputFile << routeKey << "," << g_mode_type_vector[m].mode_type << "," << g_zone_seq_2_old_zone_id[Orig] << "," << g_zone_seq_2_old_zone_id[Dest] << ","
								<< o_x_coord << "," << o_y_coord << "," << d_x_coord << "," << d_y_coord << ","
								<< totalDistance << "," << totalDistance * 1.609 << "," << straight_line_distance_mile << "," << straight_line_distance_km << ","
								<< distance_ratio << "," << totalFreeFlowTravelTime << "," << totalTravelTime << "," << volume << ",\"" << WKT_geometry << "\"\n";


							if (volume > 10.0)
							{
								std::string googleMapsLink = "https://www.google.com/maps/dir/" + std::to_string(o_y_coord) + "," + std::to_string(o_x_coord) + "/" + std::to_string(d_y_coord) + "," + std::to_string(d_x_coord) + "/";

								googleMapsFile << routeKey << "," << g_mode_type_vector[m].mode_type << "," << g_zone_seq_2_old_zone_id[Orig] << "," << g_zone_seq_2_old_zone_id[Dest] << "," << volume << ","
									<< totalDistance << "," << totalDistance * 1.609 << "," << straight_line_distance_mile << "," << straight_line_distance_km << ","
									<< distance_ratio << "," << totalFreeFlowTravelTime << "," << totalTravelTime << ",\"" << googleMapsLink << "\",\"" << WKT_geometry << "\"\n";
							}
						}
					}
				}
			}
		}
	}
	outputFile.close();
	googleMapsFile.close();
	std::cout << "Output written to " << filename << std::endl;
	std::cout << "Google Maps links saved to google_maps_od_distance.csv" << std::endl;
}

// ===== Data Structures for Aggregation =====

// Aggregation at the mode (system) level.
struct AggregationData {
	double total_volume = 0.0;
	double total_VMT = 0.0;      // VMT in miles * volume
	double total_VHT = 0.0;      // VHT in hours * volume
	double total_FF_VHT = 0.0;   // Free-flow VHT in hours * volume
	double total_delay = 0.0;    // (congested - free-flow time) in hours * volume
	int count_OD = 0;
};

// Aggregation for accessibility at origins.
struct OriginStats {
	double total_volume = 0.0;
	double total_distance = 0.0;     // Sum(distance * volume) in miles
	double total_travel_time = 0.0;  // Sum(congested travel time * volume) in minutes
	int count_dest = 0;              // Number of destination flows (nonzero OD pairs)
};

// Aggregation for accessibility at destinations.
struct DestinationStats {
	double total_volume = 0.0;
	double total_distance = 0.0;
	double total_travel_time = 0.0;
	int count_orig = 0;              // Number of origin flows (nonzero OD pairs)
};

// ===== New Function for Aggregated Performance and Accessibility =====
void GenerateAggregatedPerformanceAndAccessibility()
{
	// Mode-level aggregator.
	std::map<std::string, AggregationData> modeAgg;

	// Aggregators for origin and destination accessibility.
	std::vector<OriginStats> originStats(no_zones + 1);         // 1-indexed zones.
	std::vector<DestinationStats> destinationStats(no_zones + 1);

	// Process accessible OD pairs.
	// (Assumes linkIndices and MDODflow are 1-indexed for modes and zones.)
	for (int m = 1; m < linkIndices.size(); ++m)
	{
		std::string modeStr = g_mode_type_vector[m].mode_type;
		for (int Orig = 1; Orig < linkIndices[m].size(); ++Orig)
		{
			for (int Dest = 1; Dest < linkIndices[m][Orig].size(); ++Dest)
			{
				// Process only the first non-empty route per OD pair.
				for (int route_id = 0; route_id < linkIndices[m][Orig][Dest].size(); ++route_id)
				{
					if (linkIndices[m][Orig][Dest][route_id].empty())
						continue;

					double totalDistance = 0.0;           // in miles
					double totalFreeFlowTravelTime = 0.0;   // in minutes
					double totalTravelTime = 0.0;           // in minutes
					for (int i = linkIndices[m][Orig][Dest][route_id].size() - 1; i >= 0; --i)
					{
						long k = linkIndices[m][Orig][Dest][route_id][i];
						totalDistance += Link[k].length;
						totalFreeFlowTravelTime += Link[k].FreeTravelTime;
						totalTravelTime += Link[k].Travel_time;
					}

					float volume = MDODflow[m][Orig][Dest];
					if (volume <= 0)
						break;  // No flow for this OD pair; skip further processing.

					// Update mode-level aggregator.
					AggregationData& agg = modeAgg[modeStr];
					agg.total_volume += volume;
					agg.total_VMT += totalDistance * volume;
					// Convert travel times from minutes to hours.
					agg.total_VHT += (totalTravelTime / 60.0) * volume;
					agg.total_FF_VHT += (totalFreeFlowTravelTime / 60.0) * volume;
					agg.total_delay += ((totalTravelTime - totalFreeFlowTravelTime) / 60.0) * volume;
					agg.count_OD += 1;

					// Update origin accessibility aggregator.
					originStats[Orig].total_volume += volume;
					originStats[Orig].total_distance += totalDistance * volume;
					originStats[Orig].total_travel_time += totalTravelTime * volume;
					originStats[Orig].count_dest += 1;

					// Update destination accessibility aggregator.
					destinationStats[Dest].total_volume += volume;
					destinationStats[Dest].total_distance += totalDistance * volume;
					destinationStats[Dest].total_travel_time += totalTravelTime * volume;
					destinationStats[Dest].count_orig += 1;

					// Process only the first non-empty route for this OD pair.
					break;
				}
			}
		}
	}

	// ----- Write system_performance.csv -----
	std::ofstream systemPerfFile("system_performance.csv");
	if (!systemPerfFile.is_open())
	{
		std::cerr << "Error: Could not open system_performance.csv for writing." << std::endl;
		return;
	}
	// Header with extra columns.
	systemPerfFile << "mode_type,total_volume,PMT (VMT in miles),PHT (VHT in hours),Delay (hours),TTI,"
		<< "avg_distance_mile,avg_distance_km,avg_travel_time_min,avg_travel_delay_in_min,"
		<< "avg_speed_mph,avg_speed_kmph\n";
	for (const auto& kv : modeAgg)
	{
		const std::string& mode = kv.first;
		const AggregationData& data = kv.second;
		double TTI = (data.total_FF_VHT > 0) ? (data.total_VHT / data.total_FF_VHT) : 1.0;
		double avg_distance_mile = (data.total_volume > 0) ? (data.total_VMT / data.total_volume) : 0.0;
		double avg_distance_km = avg_distance_mile * 1.609;
		double avg_travel_time_min = (data.total_volume > 0) ? ((data.total_VHT / data.total_volume) * 60.0) : 0.0;
		double avg_travel_delay_min = (data.total_volume > 0) ? ((data.total_delay / data.total_volume) * 60.0) : 0.0;
		double avg_speed_mph = (avg_travel_time_min > 0) ? (avg_distance_mile / (avg_travel_time_min / 60.0)) : 0.0;
		double avg_speed_kmph = (avg_travel_time_min > 0) ? (avg_distance_km / (avg_travel_time_min / 60.0)) : 0.0;

		systemPerfFile << mode << ","
			<< data.total_volume << ","
			<< data.total_VMT << ","
			<< data.total_VHT << ","
			<< data.total_delay << ","
			<< TTI << ","
			<< avg_distance_mile << ","
			<< avg_distance_km << ","
			<< avg_travel_time_min << ","
			<< avg_travel_delay_min << ","
			<< avg_speed_mph << ","
			<< avg_speed_kmph << "\n";
	}
	systemPerfFile.close();

	// ----- Write origin_accessibility.csv -----
	// This file lists every origin zone (from 1 to no_zones) and includes a Google Maps search link.
	std::ofstream originFile("origin_accessibility.csv");
	if (!originFile.is_open())
	{
		std::cerr << "Error: Could not open origin_accessibility.csv for writing." << std::endl;
		return;
	}
	originFile << "origin_zone_id,total_volume,number_of_destinations,avg_distance_mile,avg_travel_time_min,google_maps_link\n";
	for (int o = 1; o <= no_zones; ++o)
	{
		double vol = originStats[o].total_volume;
		double avg_distance = (vol > 0) ? (originStats[o].total_distance / vol) : 0.0;
		double avg_travel_time = (vol > 0) ? (originStats[o].total_travel_time / vol) : 0.0;
		std::string googleLink = "";
		if (o < (int)g_node_vector.size())   // dense zone seq == centroid node seq
		{
			double x = g_node_vector[o].x;
			double y = g_node_vector[o].y;
			googleLink = "https://www.google.com/maps/search/?api=1&query=" + std::to_string(y) + "," + std::to_string(x);
		}
		// Enclose the google link in double quotes.
		originFile << g_zone_seq_2_old_zone_id[o] << ","
			<< vol << ","
			<< originStats[o].count_dest << ","
			<< avg_distance << ","
			<< avg_travel_time << ","
			<< "\"" << googleLink << "\"" << "\n";
	}
	originFile.close();

	// ----- Write destination_accessibility.csv -----
	// This file lists every destination zone (from 1 to no_zones) with a Google Maps search link.
	std::ofstream destFile("destination_accessibility.csv");
	if (!destFile.is_open())
	{
		std::cerr << "Error: Could not open destination_accessibility.csv for writing." << std::endl;
		return;
	}
	destFile << "destination_zone_id,total_volume,number_of_origins,avg_distance_mile,avg_travel_time_min,google_maps_link\n";
	for (int d = 1; d <= no_zones; ++d)
	{
		double vol = destinationStats[d].total_volume;
		double avg_distance = (vol > 0) ? (destinationStats[d].total_distance / vol) : 0.0;
		double avg_travel_time = (vol > 0) ? (destinationStats[d].total_travel_time / vol) : 0.0;
		std::string googleLink = "";
		if (d < (int)g_node_vector.size())   // dense zone seq == centroid node seq
		{
			double x = g_node_vector[d].x;
			double y = g_node_vector[d].y;
			googleLink = "https://www.google.com/maps/search/?api=1&query=" + std::to_string(y) + "," + std::to_string(x);
		}
		destFile << g_zone_seq_2_old_zone_id[d] << ","
			<< vol << ","
			<< destinationStats[d].count_orig << ","
			<< avg_distance << ","
			<< avg_travel_time << ","
			<< "\"" << googleLink << "\"" << "\n";
	}
	destFile.close();

	// ----- Write inaccessible_od.csv -----
	// For each mode, origin, and destination, check for a non-empty route.
	// If no non-empty route exists and the OD flow (volume) is > 0, then output this OD pair.
	std::ofstream inacFile("inaccessible_od.csv");
	if (!inacFile.is_open())
	{
		std::cerr << "Error: Could not open inaccessible_od.csv for writing." << std::endl;
		return;
	}
	inacFile << "mode_type,origin_zone_id,destination_zone_id,google_maps_http_link\n";
	for (int m = 1; m < linkIndices.size(); ++m)
	{
		std::string modeStr = g_mode_type_vector[m].mode_type;
		for (int Orig = 1; Orig < linkIndices[m].size(); ++Orig)
		{
			for (int Dest = 1; Dest < linkIndices[m][Orig].size(); ++Dest)
			{
				// Check if any route exists for this OD.
				bool routeExists = false;
				for (int route_id = 0; route_id < linkIndices[m][Orig][Dest].size(); ++route_id)
				{
					if (!linkIndices[m][Orig][Dest][route_id].empty())
					{
						routeExists = true;
						break;
					}
				}
				// If no route exists but the flow (volume) is > 0, then mark as inaccessible.
				if (!routeExists && MDODflow[m][Orig][Dest] > 0)
				{
					std::string googleLink = "";
					if (Orig < (int)g_node_vector.size() &&
						Dest < (int)g_node_vector.size())   // dense zone seq == node seq
					{
						int originNodeIndex = Orig;
						int destNodeIndex = Dest;
						double o_x = g_node_vector[originNodeIndex].x;
						double o_y = g_node_vector[originNodeIndex].y;
						double d_x = g_node_vector[destNodeIndex].x;
						double d_y = g_node_vector[destNodeIndex].y;
						// Construct a Google Maps directions link.
						googleLink = "https://www.google.com/maps/dir/" + std::to_string(o_y) + "," + std::to_string(o_x) + "/" +
							std::to_string(d_y) + "," + std::to_string(d_x) + "/";
					}
					inacFile << modeStr << ","
						<< g_zone_seq_2_old_zone_id[Orig] << ","
						<< g_zone_seq_2_old_zone_id[Dest] << ","
						<< "\"" << googleLink << "\"" << "\n";
				}
			}
		}
	}
	inacFile.close();

	std::cout << "System performance written to system_performance.csv" << std::endl;
	std::cout << "Origin accessibility written to origin_accessibility.csv" << std::endl;
	std::cout << "Destination accessibility written to destination_accessibility.csv" << std::endl;
	std::cout << "Inaccessible OD pairs written to inaccessible_od.csv" << std::endl;

}
int get_number_of_nodes_from_node_file(int& number_of_zones, int& l_FirstThruNode)
{
	number_of_zones = 0;
	CDTACSVParser parser_node;
	l_FirstThruNode = 1;

	// INTERNAL RENUMBERING. The kernel historically required dense ids with
	// zone_id == node_id (1..Z centroids first); sparse agency ids made every
	// zone-dimensioned array scale with the LARGEST zone id (GBs of memory and
	// crashes on small subareas). We now renumber internally: centroids take
	// node seq 1..Z in file order (dense zone id == centroid node seq), other
	// nodes follow. Every output still reports the ORIGINAL ids via
	// g_node_vector[].node_id / g_map_node_seq_no_2_external_node_id /
	// g_zone_seq_2_old_zone_id; zone-id inputs (demand, stored routes) are
	// translated through g_map_old_zone_id_2_zone_seq at read time.
	struct NodeRow { int node_id; int zone_id; double x; double y; };
	std::vector<NodeRow> rows;
	int max_zone_id = 0;

	if (parser_node.OpenCSVFile("node.csv", true))
	{
		while (parser_node.ReadRecord())
		{
			NodeRow r; r.node_id = 0; r.zone_id = 0; r.x = 0.0; r.y = 0.0;
			parser_node.GetValueByFieldName("node_id", r.node_id);
			parser_node.GetValueByFieldName("zone_id", r.zone_id);
			parser_node.GetValueByFieldName("x_coord", r.x);
			parser_node.GetValueByFieldName("y_coord", r.y);
			if (r.zone_id >= 1 && r.zone_id > max_zone_id)
				max_zone_id = r.zone_id;
			rows.push_back(r);
		}
		parser_node.CloseCSVFile();
	}
	int number_of_nodes = (int)rows.size();

	std::vector<size_t> ordered;
	ordered.reserve(rows.size());
	for (size_t i = 0; i < rows.size(); i++)
		if (rows[i].zone_id >= 1)
			ordered.push_back(i);
	int Z = (int)ordered.size();
	for (size_t i = 0; i < rows.size(); i++)
		if (rows[i].zone_id < 1)
			ordered.push_back(i);

	if (Z > 0 && max_zone_id != Z)
		printf("NOTE: sparse zone ids detected (%d zones, max zone id %d). "
			"Renumbered internally to 1..%d; all outputs keep the original ids.\n",
			Z, max_zone_id, Z);

	g_map_external_node_id_2_node_seq_no.clear();
	g_map_external_node_id_2_node_seq_no.reserve((size_t)number_of_nodes * 2 + 1);
	g_map_node_seq_no_2_external_node_id.assign(number_of_nodes + 1, 0);
	g_map_old_zone_id_2_zone_seq.clear();
	g_map_old_zone_id_2_zone_seq.reserve((size_t)Z * 2 + 1);
	g_node_vector.clear();
	g_node_vector.resize(number_of_nodes + 1);
	g_zone_seq_2_old_zone_id.assign(Z + 1, 0);

	for (int seq = 1; seq <= number_of_nodes; seq++)
	{
		const NodeRow& r = rows[ordered[seq - 1]];
		if (g_map_external_node_id_2_node_seq_no.find(r.node_id) !=
			g_map_external_node_id_2_node_seq_no.end())
			printf("WARNING: duplicate node_id %d in node.csv; keeping the first.\n",
				r.node_id);
		else
			g_map_external_node_id_2_node_seq_no[r.node_id] = seq;
		g_map_node_seq_no_2_external_node_id[seq] = r.node_id;
		g_node_vector[seq].x = r.x;
		g_node_vector[seq].y = r.y;
		g_node_vector[seq].node_id = r.node_id;        // ORIGINAL id, for outputs
		if (r.zone_id >= 1)                            // seq <= Z by construction
		{
			g_node_vector[seq].zone_id = seq;           // dense internal zone id
			g_node_vector[seq].internal_zone_no = seq - 1;
			if (g_map_old_zone_id_2_zone_seq.find(r.zone_id) !=
				g_map_old_zone_id_2_zone_seq.end())
				printf("WARNING: zone_id %d appears on more than one node; "
					"demand maps to the first.\n", r.zone_id);
			else
				g_map_old_zone_id_2_zone_seq[r.zone_id] = seq;
			g_zone_seq_2_old_zone_id[seq] = r.zone_id;
		}
		if (g_tap_log_file == 1)
			fprintf(logfile, "node_id = %d, node_seq_no = %d\n", r.node_id, seq);
	}

	number_of_zones = Z;
	if (first_through_node_id_input == -1)
	{
		// auto: first non-centroid seq (zones-first ordering). When EVERY node
		// is a zone (small teaching networks), keep the legacy behavior of
		// FirstThruNode = 1 so all nodes remain through-traversable.
		l_FirstThruNode = (Z < number_of_nodes) ? (Z + 1) : 1;
	}

	return number_of_nodes;
}

int get_number_of_links_from_link_file()
{
	int number_of_links = 0;
	std::ifstream link_file("link.csv");
	std::string line;
	if (link_file.is_open())
	{
		std::getline(link_file, line); // header
		while (std::getline(link_file, line))
		{
			if (line.find_first_not_of(" \t\r") != std::string::npos)
				number_of_links++; // allocation upper bound; ReadLinks sets accepted count
		}
	}

	return number_of_links;
}

void createSettingsFile(const std::string& fileName) {
	std::ofstream file(fileName);
	if (!file.is_open()) {
		std::cerr << "Could not create the file: " << fileName << std::endl;
		return;
	}

	// Writing the headers of the CSV
	file << "number_of_iterations,number_of_processors,demand_period_starting_hours,demand_period_ending_hours,first_through_node_id,base_demand_mode,route_output,vehicle_output,log_file,odme_mode,odme_vmt\n";

	// Writing the sample data (from your provided file)
	file << "20,8,7,8,-1,0,1,0,0,0,0\n";

	file.close();
	std::cout << "sample_settings.csv file created successfully!" << std::endl;
}

void read_settings_file()
{

	createSettingsFile("sample_settings.csv");

	CDTACSVParser parser_settings;

	if (parser_settings.OpenCSVFile("settings.csv", true))
	{
		while (parser_settings.ReadRecord())  // if this line contains [] mark, then we will also read
			// field headers.
		{
			g_number_of_processors = 4;

			parser_settings.GetValueByFieldName("number_of_iterations", TotalAssignIterations);

			if (TotalAssignIterations == 0)
				g_accessibility_only_mode = 1;
			else
				g_accessibility_only_mode = 0; 

			parser_settings.GetValueByFieldName("number_of_processors", g_number_of_processors);
			parser_settings.GetValueByFieldName("demand_period_starting_hours", demand_period_starting_hours);
			parser_settings.GetValueByFieldName("demand_period_ending_hours", demand_period_ending_hours);
			parser_settings.GetValueByFieldName("first_through_node_id", first_through_node_id_input);

			parser_settings.GetValueByFieldName("log_file", g_tap_log_file);
			parser_settings.GetValueByFieldName("base_demand_mode", g_base_demand_mode);
			parser_settings.GetValueByFieldName("demand_format", g_demand_format);
			parser_settings.GetValueByFieldName("added_delay_per_mile", g_added_delay_per_mile);
			parser_settings.GetValueByFieldName("convergence_gap_pct", g_convergence_gap_pct);
			parser_settings.GetValueByFieldName("link_output", g_link_performance_output, false, false);
			parser_settings.GetValueByFieldName("accessibility_output", g_accessibility_output, false, false);
			parser_settings.GetValueByFieldName("qvdf_volume_threshold", g_qvdf_volume_threshold, false, false);
			parser_settings.GetValueByFieldName("sp_algorithm", g_sp_algorithm, false, false);
			parser_settings.GetValueByFieldName("convergence_consecutive", g_convergence_consecutive);
			parser_settings.GetValueByFieldName("relative_gap_standard", g_relative_gap_standard);
			parser_settings.GetValueByFieldName("assignment_method", g_assignment_method);
			parser_settings.GetValueByFieldName("odme_mode", g_ODME_mode);
			parser_settings.GetValueByFieldName("odme_vmt", g_ODME_obs_VMT);
			parser_settings.GetValueByFieldName("route_output", shortest_path_log_flag);
			parser_settings.GetValueByFieldName("vehicle_output", vehicle_log_flag);
			parser_settings.GetValueByFieldName("warm_start_times", g_warm_start_times_file, false);
			parser_settings.GetValueByFieldName("warm_start_flows", g_warm_start_flows_file, false);
			parser_settings.GetValueByFieldName("flow_snapshot", g_flow_snapshot, false, false);
			parser_settings.GetValueByFieldName("column_output", g_column_output, false, false);
			parser_settings.GetValueByFieldName("column_file_output", g_column_file_output, false, false);
			parser_settings.GetValueByFieldName("auto_calibration", g_auto_calibration, false, false);
			parser_settings.GetValueByFieldName("auto_calibration_config", g_auto_calibration_config_file, false);
			parser_settings.GetValueByFieldName("departure_profile_file", g_departure_profile_file, false);
			if (g_force_auto_calibration)
				g_auto_calibration = 1;
			// Calibration needs the accumulated in-memory route-policy columns for
			// origin-profile time shifting. Level 2 is mandatory for every inner
			// equilibrium, but the DTAC file is still written only after the final
			// accepted outer iteration.
			if (g_auto_calibration != 0 && g_column_output != 2)
			{
				printf("auto_calibration: forcing column_output=2 for in-memory link-arrival profiles (requested %d).\n",
					g_column_output);
				g_column_output = 2;
			}
			if (g_auto_calibration != 0 && g_number_of_processors > 0)
			{
				taplite::AutoCalibrationConfig thread_config =
					taplite::AutoCalibrationConfig::Load(
						g_auto_calibration_config_file,
						demand_period_starting_hours,
						demand_period_ending_hours);
				const int calibrated_workers = std::min(
					g_number_of_processors, thread_config.workers);
				g_calibration_minimum_column_od_volume =
					thread_config.minimum_column_od_volume;
				g_calibration_maximum_column_paths_per_od =
					thread_config.maximum_column_paths_per_od;
				g_calibration_minimum_column_path_share =
					thread_config.minimum_column_path_share;
				printf("auto_calibration: bounded DTAC-v2 columns at OD volume >= %.6f, "
					"max %d paths/OD, minimum retained share %.6f.\n",
					g_calibration_minimum_column_od_volume,
					g_calibration_maximum_column_paths_per_od,
					g_calibration_minimum_column_path_share);
				if (calibrated_workers != g_number_of_processors)
				{
					printf("auto_calibration: limiting assignment workers from %d to %d.\n",
						g_number_of_processors, calibrated_workers);
					g_number_of_processors = calibrated_workers;
				}
			}
			// Leveled: 0 off, 1 = DTAC v1 last-iteration paths, 2 = DTAC v2
			// theta-share columns (L3), 3+ reserved -> behaves as 2 with a warning.
			if (g_column_output >= 3)
			{
				printf("WARNING: column_output=%d is reserved for a future per-iteration trace level; treating it as 2 (DTAC v2 theta-share columns).\n",
					g_column_output);
				g_column_output = 2;
			}
			parser_settings.GetValueByFieldName("warm_start_columns", g_warm_start_columns_file, false);
			parser_settings.GetValueByFieldName("column_adjust_sweeps", g_column_adjust_sweeps, false, false);

			// L3 freeze/replay: number_of_iterations=0 normally means
			// accessibility-only mode, but with warm_start_columns it means PURE
			// REPLAY -- scatter the stored columns, run the (optional) adjustment
			// sweeps, report the TRUE gap, write all normal outputs -- so keep
			// the full assignment path alive.
			if (TotalAssignIterations == 0 && !g_warm_start_columns_file.empty())
			{
				g_accessibility_only_mode = 0;
				printf("warm_start_columns + number_of_iterations=0: FREEZE/REPLAY mode (no FW iterations; scatter + adjust sweeps + TRUE-gap report + normal outputs).\n");
			}

			break;
		}

		parser_settings.CloseCSVFile();
	}

	return;
}

void createModeTypeFile(const std::string& fileName) {
	std::ofstream file(fileName);
	if (!file.is_open()) {
		std::cerr << "Could not create the file: " << fileName << std::endl;
		return;
	}

	// Writing the headers of the CSV
	file << "mode_type_id,mode_type,name,vot,pce,occ,demand_file\n";

	// Writing the sample data (from your provided file)
	file << "1,sov,DRIVE, 10, 1, 1,demand.csv\n";
	file << "2,hov,HOV, 10, 1, 2,demand_hov.csv\n";
	file << "3,trk,truck, 10, 2, 1,demand_trk.csv\n";

	std::cout << "sample_mode_type.csv file created successfully!" << std::endl;
}


void read_mode_type_file()
{
	createModeTypeFile("sample_mode_type.csv");
	CDTACSVParser parser_mode_type;

	if (parser_mode_type.OpenCSVFile("mode_type.csv", true))
	{

		number_of_modes = 0;
		while (parser_mode_type.ReadRecord())  // if this line contains [] mark, then we will also read
			// field headers.
		{
			number_of_modes += 1;
			if (number_of_modes < MAX_MODE_TYPES)
			{
				g_mode_type_vector[number_of_modes].mode_type = "auto";
				g_mode_type_vector[number_of_modes].vot = 10;
				g_mode_type_vector[number_of_modes].pce = 1;
				g_mode_type_vector[number_of_modes].occ = 1; 
				g_mode_type_vector[number_of_modes].demand_file = "demand.csv";

				parser_mode_type.GetValueByFieldName("mode_type", g_mode_type_vector[number_of_modes].mode_type);
				parser_mode_type.GetValueByFieldName("vot", g_mode_type_vector[number_of_modes].vot);
				parser_mode_type.GetValueByFieldName("pce", g_mode_type_vector[number_of_modes].pce);
				parser_mode_type.GetValueByFieldName("occ", g_mode_type_vector[number_of_modes].occ);
				parser_mode_type.GetValueByFieldName("operating_cost", g_mode_type_vector[number_of_modes].op_cost);
				parser_mode_type.GetValueByFieldName("demand_file", g_mode_type_vector[number_of_modes].demand_file);

				g_mode_type_vector[number_of_modes].dedicated_shortest_path = 1;
				parser_mode_type.GetValueByFieldName("dedicated_shortest_path", g_mode_type_vector[number_of_modes].dedicated_shortest_path);

				if (number_of_modes == 1)
					g_mode_type_vector[1].dedicated_shortest_path = 1;  // reset; 

			}

		}

		parser_mode_type.CloseCSVFile();
	}

	if (number_of_modes == 0)  //default 
	{
		g_mode_type_vector[1].demand_file = "demand.csv";
		g_mode_type_vector[1].mode_type = "auto";
		g_mode_type_vector[1].vot = 10;
		g_mode_type_vector[1].pce = 1;
		g_mode_type_vector[1].occ = 1;
		g_mode_type_vector[1].dedicated_shortest_path = 1;  // reset; 
		number_of_modes = 1;
	}

	printf("number_of_modes = %d\n", number_of_modes);
	return;
}
// Initialize the 5D vector
void InitializeLinkIndices(int num_modes, int no_zones, int max_routes)
{


	auto start = std::chrono::high_resolution_clock::now();  // Start timing



	// MEMORY/SPEED OPT (multi-mode, volume-only Frank-Wolfe): the 5D
	// route-policy store scales as modes x zones x zones x routes and is the
	// dominant memory user (billions of std::vector headers for large zone
	// counts). It is ONLY needed to store/output path results. Link volumes
	// (mode_MainVolume[m]) are accumulated by walking MinPathPredLink and do
	// NOT need it. The route-storage path also runs under #pragma omp critical
	// (huge serialization cost on big networks).
	//
	// So: when route output is OFF (shortest_path_log_flag == 0) we leave
	// linkIndices EMPTY (size 0). Every populate/consume site is guarded by
	// `if (linkIndices.size() > 0)` (and AddLinkSequence re-checks per slice),
	// so an empty store cleanly disables all path storage -> fast FW, low mem,
	// correct per-mode link volumes. When route output IS requested, allocate
	// the full inner structure only for dedicated (main) modes.
	linkIndices.clear();
	if (shortest_path_log_flag == 1)
	{
		linkIndices = std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>>(
			num_modes + 1);
		for (int m = 1; m <= num_modes; ++m)
		{
			if (g_mode_type_vector[m].dedicated_shortest_path == 1)
			{
				linkIndices[m] = std::vector<std::vector<std::vector<std::vector<int>>>>(
					no_zones + 1,
					std::vector<std::vector<std::vector<int>>>(
						no_zones + 1,
						std::vector<std::vector<int>>(
							max_routes + 1
						)
					)
				);
			}
			// non-dedicated modes: linkIndices[m] stays empty -> memory saving
		}
	}
	// route output OFF: linkIndices stays empty -> no path storage at all

	// Record the end time
	auto end = std::chrono::high_resolution_clock::now();

	// Calculate the duration in milliseconds

	// Calculate the duration in seconds
	auto duration = end - start;

	// Convert to hours, minutes, seconds
	auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
	auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration % std::chrono::hours(1));
	auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration % std::chrono::minutes(1));
	auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(duration % std::chrono::seconds(1));

	printf("Memmory creation time for 5D link path matrix: %lld hours %lld minutes %lld seconds %lld ms\n", hours.count(), minutes.count(), seconds.count(), milliseconds.count());
}

// SP-equivalence auto-detect (exact, always on): after all inputs are read,
// group modes whose routing problems are IDENTICAL -- same mode_allowed_use
// permission and same mode_AdditionalCost on every link. Members reuse the
// representative's shortest-path tree in FindMinCostRoutes /
// All_or_Nothing_Assign instead of recomputing it (bit-identical results,
// redundant Minpath calls skipped). PCE/occupancy do not enter routing and
// never split a group. Legacy manual dedicated_shortest_path=0 (hardcoded
// "reuse mode 1") is left untouched and excluded from auto-grouping.
static void BuildSPModeGroups()
{
	for (int m = 1; m <= number_of_modes; m++)
	{
		g_rep_mode[m] = m;
		g_sp_scratch_row[m] = 0;
	}

	for (int m = 2; m <= number_of_modes; m++)
	{
		if (g_mode_type_vector[m].dedicated_shortest_path == 0)
			continue;  // legacy manual reuse of mode 1; keep its existing code path

		for (int r = 1; r < m; r++)
		{
			if (g_mode_type_vector[r].dedicated_shortest_path == 0 || g_rep_mode[r] != r)
				continue;  // representative candidates only

			bool same = true;
			for (int k = 1; k <= number_of_links && same; k++)
			{
				if (Link[k].mode_allowed_use[r] != Link[k].mode_allowed_use[m] ||
					Link[k].mode_AdditionalCost[r] != Link[k].mode_AdditionalCost[m])
					same = false;
			}

			if (same)
			{
				g_rep_mode[m] = r;
				break;
			}
		}
	}

	// CostTo scratch rows: row 0 is the shared immediate-consumption scratch;
	// each representative whose tree is shared by members gets a private
	// retained row so the member modes can consume it later in the mode loop.
	g_sp_scratch_rows = 1;
	for (int r = 1; r <= number_of_modes; r++)
	{
		if (g_rep_mode[r] != r)
			continue;
		bool has_member = false;
		for (int m = r + 1; m <= number_of_modes; m++)
			if (g_rep_mode[m] == r)
				has_member = true;
		if (has_member)
		{
			g_sp_scratch_row[r] = g_sp_scratch_rows++;
			for (int m = r + 1; m <= number_of_modes; m++)
				if (g_rep_mode[m] == r)
					g_sp_scratch_row[m] = g_sp_scratch_row[r];
		}
	}

	// always log the grouping -- silent cleverness loses trust
	if (number_of_modes > 1)
	{
		std::string msg = "SP tree groups:";
		for (int r = 1; r <= number_of_modes; r++)
		{
			if (g_mode_type_vector[r].dedicated_shortest_path == 0 || g_rep_mode[r] != r)
				continue;
			msg += " {" + g_mode_type_vector[r].mode_type;
			for (int m = r + 1; m <= number_of_modes; m++)
				if (g_rep_mode[m] == r)
					msg += "," + g_mode_type_vector[m].mode_type;
			msg += "} -> " + g_mode_type_vector[r].mode_type + ";";
		}
		for (int m = 1; m <= number_of_modes; m++)
			if (g_mode_type_vector[m].dedicated_shortest_path == 0)
				msg += " {" + g_mode_type_vector[m].mode_type + "} -> " + g_mode_type_vector[1].mode_type + " (manual dedicated_shortest_path=0);";
		printf("%s\n", msg.c_str());
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "%s\n", msg.c_str());
	}
}

// Zero-demand-origin skip: precompute, per (mode, origin), whether any
// destination carries demand (> 0.000001, the same gate FindMinCostRoutes
// uses), and per SP-group representative whether any member mode is active.
// FindMinCostRoutes then skips the Minpath call for inactive trees and
// BIGM-fills the MDRouteCost row for inactive modes -- exactly what the full
// computation produces when the demand row is zero.
static void BuildOriginActivity()
{
	if (MDODflow == NULL)
		return;

	g_mode_origin_active = (char**)Alloc_2D(number_of_modes, no_zones, sizeof(char));
	g_tree_origin_active = (char**)Alloc_2D(number_of_modes, no_zones, sizeof(char));

	long long active_pairs = 0, total_pairs = 0;
	for (int m = 1; m <= number_of_modes; m++)
	{
		for (int Orig = 1; Orig <= no_zones; Orig++)
		{
			char active = 0;
			for (int Dest = 1; Dest <= no_zones; Dest++)
			{
				if (MDODflow[m][Orig][Dest] > 0.000001)
				{
					active = 1;
					break;
				}
			}
			g_mode_origin_active[m][Orig] = active;
			total_pairs++;
			if (active)
				active_pairs++;
		}
	}

	for (int m = 1; m <= number_of_modes; m++)
	{
		int rep = (g_rep_mode[m] >= 1) ? g_rep_mode[m] : m;
		for (int Orig = 1; Orig <= no_zones; Orig++)
			if (g_mode_origin_active[m][Orig])
				g_tree_origin_active[rep][Orig] = 1;
	}

	printf("zero-demand-origin skip: %lld of %lld (mode, origin) pairs active\n",
		active_pairs, total_pairs);
	if (summary_log_file != NULL)
		fprintf(summary_log_file, "zero-demand-origin skip: %lld of %lld (mode, origin) pairs active\n",
			active_pairs, total_pairs);
}

// L1 warm start (settings 'warm_start_times', default empty = off): preload
// congested Travel_time from a prior run's link_performance file -- compact or
// full CSV (link_id + travel_time columns) or the DTLP binary written by
// link_output=3 -- matching by EXTERNAL link id, never row order. GenCost is
// recomposed exactly as Link_GenCost does (mode_AdditionalCost[1] +
// Travel_time); Minpath itself reads Travel_time + mode_AdditionalCost[mode],
// which is static, so nothing else needs refreshing. Called after the initial
// UpdateLinkCost and before the first FindMinCostRoutes: the first
// all-or-nothing routes on congested times, then every later iteration
// recomputes costs from volumes normally -- only the starting point changes,
// the equilibrium does not. Missing/unreadable file = loud warning, cold start.
static void ApplyWarmStartTimes()
{
	if (g_warm_start_times_file.empty())
		return;

	int matched = 0, records = 0;

	// binary DTLP? (header magic 0x504C5444, int32 version, int32 n_links;
	// per link: int32 link_id, from, to; float volume, doc, travel_time, speed, vmt)
	FILE* f = NULL;
	fopen_s(&f, g_warm_start_times_file.c_str(), "rb");
	if (f == NULL)
	{
		printf("WARNING: warm_start_times file '%s' cannot be opened; starting cold.\n",
			g_warm_start_times_file.c_str());
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "WARNING: warm_start_times file '%s' cannot be opened; starting cold.\n",
				g_warm_start_times_file.c_str());
		return;
	}

	int magic = 0;
	bool is_dtlp = (fread(&magic, sizeof(int), 1, f) == 1 && magic == 0x504C5444);
	if (is_dtlp)
	{
		int version = 0, n_links_in = 0;
		if (fread(&version, sizeof(int), 1, f) == 1 && fread(&n_links_in, sizeof(int), 1, f) == 1)
		{
			for (int i = 0; i < n_links_in; i++)
			{
				int ids[3];
				float vals[5];
				if (fread(ids, sizeof(int), 3, f) != 3 || fread(vals, sizeof(float), 5, f) != 5)
					break;
				records++;
				int k = LinkSeqFromExternal(ids[0]);
				if (k != 0 && vals[2] > 0.0f)
				{
					Link[k].Travel_time = vals[2];
					Link[k].GenCost = Link[k].mode_AdditionalCost[1] + Link[k].Travel_time;
					matched++;
				}
			}
		}
		fclose(f);
	}
	else
	{
		fclose(f);
		// CSV: both compact and full link_performance formats carry
		// link_id + travel_time columns; match by header name.
		CDTACSVParser parser_warm;
		if (!parser_warm.OpenCSVFile(g_warm_start_times_file, false))
		{
			printf("WARNING: warm_start_times file '%s' is not a DTLP .bin and cannot be parsed as CSV; starting cold.\n",
				g_warm_start_times_file.c_str());
			if (summary_log_file != NULL)
				fprintf(summary_log_file, "WARNING: warm_start_times file '%s' unreadable; starting cold.\n",
					g_warm_start_times_file.c_str());
			return;
		}
		while (parser_warm.ReadRecord())
		{
			int link_id = -1;
			double tt = -1;
			if (!parser_warm.GetValueByFieldName("link_id", link_id, false) ||
				!parser_warm.GetValueByFieldName("travel_time", tt, false))
				continue;
			records++;
			int k = LinkSeqFromExternal(link_id);
			if (k != 0 && tt > 0.0)
			{
				Link[k].Travel_time = tt;
				Link[k].GenCost = Link[k].mode_AdditionalCost[1] + Link[k].Travel_time;
				matched++;
			}
		}
		parser_warm.CloseCSVFile();
	}

	printf("warm_start_times: preloaded congested times on %d of %d links (%d records read from '%s', %s)\n",
		matched, number_of_links, records, g_warm_start_times_file.c_str(), is_dtlp ? "DTLP binary" : "CSV");
	if (summary_log_file != NULL)
		fprintf(summary_log_file, "warm_start_times: preloaded congested times on %d of %d links (%d records read from '%s', %s)\n",
			matched, number_of_links, records, g_warm_start_times_file.c_str(), is_dtlp ? "DTLP binary" : "CSV");
	if (matched == 0)
	{
		printf("WARNING: warm_start_times matched 0 links (id mismatch?); effectively a cold start.\n");
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "WARNING: warm_start_times matched 0 links (id mismatch?); effectively a cold start.\n");
	}
}

// ---- P1: L2 warm start from FLOWS + DTLR flow snapshot ----
// 16-byte demand fingerprint over the CURRENT OD table (MDODflow): {total
// demand, positive-cell count, min zone id, max zone id touched by demand}
// folded into 4 floats. Deterministic: fixed mode/origin/destination loop
// order with a double accumulator, so identical demand inputs always fold to
// identical bytes. Written into every DTLR snapshot; a warm_start_flows file
// whose fingerprint differs from the current OD table triggers the loud L1
// demotion in ApplyWarmStartFlows.
static void ComputeDemandFingerprint(float fp[4])
{
	double sum = 0.0, cnt = 0.0;
	int min_zone = 0, max_zone = 0;
	if (MDODflow != NULL)
	{
		for (int m = 1; m <= number_of_modes; m++)
			for (int Orig = 1; Orig <= no_zones; Orig++)
				for (int Dest = 1; Dest <= no_zones; Dest++)
				{
					double v = MDODflow[m][Orig][Dest];
					if (v > 0.0)
					{
						sum += v;
						cnt += 1.0;
						int lo = (Orig < Dest) ? Orig : Dest;
						int hi = (Orig > Dest) ? Orig : Dest;
						if (min_zone == 0 || lo < min_zone) min_zone = lo;
						if (hi > max_zone) max_zone = hi;
					}
				}
	}
	fp[0] = (float)sum;
	fp[1] = (float)cnt;
	fp[2] = (float)min_zone;
	fp[3] = (float)max_zone;
}

// DTLR flow snapshot writer (settings 'flow_snapshot'=1, or automatic with
// link_output=3): link_flows.bin. Header: int32 magic 0x524C5444 ('DTLR'),
// int32 version=1, int32 n_links, int32 n_modes, then the 16-byte demand
// fingerprint (4 floats). Per link: int32 external_link_id; double MainVolume;
// number_of_modes doubles of per-mode (pce-free) volume. Input format for
// settings 'warm_start_flows'.
static void WriteFlowSnapshotDTLR(const char* filename, double* MainVolume)
{
	FILE* f = NULL;
	fopen_s(&f, filename, "wb");
	if (f == NULL)
	{
		printf("WARNING: flow_snapshot: cannot open '%s' for writing; snapshot skipped.\n", filename);
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "WARNING: flow_snapshot: cannot open '%s' for writing; snapshot skipped.\n", filename);
		return;
	}
	int hdr[4] = { 0x524C5444 /*'DTLR'*/, 1, number_of_links, number_of_modes };
	fwrite(hdr, sizeof(int), 4, f);
	float fp[4];
	ComputeDemandFingerprint(fp);
	fwrite(fp, sizeof(float), 4, f);
	for (int k = 1; k <= number_of_links; k++)
	{
		fwrite(&Link[k].external_link_id, sizeof(int), 1, f);
		fwrite(&MainVolume[k], sizeof(double), 1, f);
		for (int m = 1; m <= number_of_modes; m++)
			fwrite(&Link[k].mode_MainVolume[m], sizeof(double), 1, f);
	}
	fclose(f);
	printf("flow_snapshot: wrote %s (DTLR v1, %d links, %d modes, demand fingerprint [%.6g, %.0f, %.0f, %.0f])\n",
		filename, number_of_links, number_of_modes, fp[0], fp[1], fp[2], fp[3]);
	if (summary_log_file != NULL)
		fprintf(summary_log_file, "flow_snapshot: wrote %s (DTLR v1, %d links, %d modes, demand fingerprint [%.6g, %.0f, %.0f, %.0f])\n",
			filename, number_of_links, number_of_modes, fp[0], fp[1], fp[2], fp[3]);
}

// L2 warm start from FLOWS (settings 'warm_start_flows', default "" = off).
// Reads a DTLR snapshot and, when its demand fingerprint matches the CURRENT
// OD table exactly, restores MainVolume + Link.mode_MainVolume by EXTERNAL
// link id and returns true: the caller then skips the iteration-0
// all-or-nothing entirely and Frank-Wolfe continues from the restored point
// (a convex combination of AoN assignments of the same demand, so a valid FW
// iterate). On ANY incompatibility -- unreadable file, wrong magic, mode-count
// mismatch, or demand-fingerprint mismatch -- it warns LOUDLY and demotes to
// L1 behavior: the snapshot flows are used ONLY to derive congested times
// (Link_Travel_Time at the snapshot volumes, GenCost recomposed as in
// ApplyWarmStartTimes), MainVolume is left at its cold initialization, and it
// returns false so the normal cold AoN start runs. Never fails the run.
static bool ApplyWarmStartFlows(double* MainVolume)
{
	if (g_warm_start_flows_file.empty())
		return false;

	FILE* f = NULL;
	fopen_s(&f, g_warm_start_flows_file.c_str(), "rb");
	if (f == NULL)
	{
		printf("WARNING: warm_start_flows file '%s' cannot be opened; starting cold.\n",
			g_warm_start_flows_file.c_str());
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "WARNING: warm_start_flows file '%s' cannot be opened; starting cold.\n",
				g_warm_start_flows_file.c_str());
		return false;
	}

	int hdr[4] = { 0, 0, 0, 0 };
	float fp_in[4] = { 0, 0, 0, 0 };
	if (fread(hdr, sizeof(int), 4, f) != 4 || hdr[0] != 0x524C5444 ||
		fread(fp_in, sizeof(float), 4, f) != 4)
	{
		printf("WARNING: warm_start_flows file '%s' is not a DTLR flow snapshot (magic mismatch); starting cold.\n",
			g_warm_start_flows_file.c_str());
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "WARNING: warm_start_flows file '%s' is not a DTLR flow snapshot; starting cold.\n",
				g_warm_start_flows_file.c_str());
		fclose(f);
		return false;
	}
	int n_links_in = hdr[2];
	int n_modes_in = hdr[3];

	float fp_now[4];
	ComputeDemandFingerprint(fp_now);
	bool fp_match = (fp_now[0] == fp_in[0] && fp_now[1] == fp_in[1] &&
		fp_now[2] == fp_in[2] && fp_now[3] == fp_in[3]);
	bool modes_match = (n_modes_in == number_of_modes);
	bool restore = fp_match && modes_match;

	if (!fp_match)
	{
		printf("WARNING: warm_start_flows '%s' demand fingerprint [%.6g, %.0f, %.0f, %.0f] does NOT match the current OD table [%.6g, %.0f, %.0f, %.0f].\n"
			"         FALLING BACK to L1: snapshot flows seed congested TIMES only; a normal cold all-or-nothing start follows.\n",
			g_warm_start_flows_file.c_str(), fp_in[0], fp_in[1], fp_in[2], fp_in[3],
			fp_now[0], fp_now[1], fp_now[2], fp_now[3]);
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "WARNING: warm_start_flows demand fingerprint mismatch (file [%.6g, %.0f, %.0f, %.0f] vs current [%.6g, %.0f, %.0f, %.0f]); demoted to L1 times-only warm start.\n",
				fp_in[0], fp_in[1], fp_in[2], fp_in[3], fp_now[0], fp_now[1], fp_now[2], fp_now[3]);
	}
	if (!modes_match)
	{
		printf("WARNING: warm_start_flows '%s' has %d modes but the current run has %d; demoted to L1 times-only warm start.\n",
			g_warm_start_flows_file.c_str(), n_modes_in, number_of_modes);
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "WARNING: warm_start_flows mode-count mismatch (%d vs %d); demoted to L1 times-only warm start.\n",
				n_modes_in, number_of_modes);
	}

	// L1-demotion path accumulates snapshot volumes here (cold MainVolume kept).
	std::vector<double> snap_vol;
	if (!restore)
	{
		snap_vol.resize(number_of_links + 1);
		for (int k = 1; k <= number_of_links; k++)
			snap_vol[k] = MainVolume[k];
	}

	int matched = 0, records = 0, unmatched_file = 0;
	std::vector<double> mode_vol(n_modes_in + 1, 0.0);
	std::vector<char> link_matched(number_of_links + 1, 0);
	for (int i = 0; i < n_links_in; i++)
	{
		int ext_id = 0;
		double vol = 0.0;
		if (fread(&ext_id, sizeof(int), 1, f) != 1 || fread(&vol, sizeof(double), 1, f) != 1)
			break;
		bool ok = true;
		for (int m = 1; m <= n_modes_in; m++)
			if (fread(&mode_vol[m], sizeof(double), 1, f) != 1) { ok = false; break; }
		if (!ok)
			break;
		records++;
		int k = LinkSeqFromExternal(ext_id);
		if (k == 0)
		{
			if (vol > 0.0)
				unmatched_file++;
			continue;
		}
		matched++;
		link_matched[k] = 1;
		if (restore)
		{
			MainVolume[k] = vol;
			for (int m = 1; m <= number_of_modes; m++)
				Link[k].mode_MainVolume[m] = mode_vol[m];
		}
		else
			snap_vol[k] = vol;
	}
	fclose(f);

	if (matched == 0)
	{
		printf("WARNING: warm_start_flows matched 0 of %d snapshot links (id mismatch?); starting cold.\n", records);
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "WARNING: warm_start_flows matched 0 of %d snapshot links; starting cold.\n", records);
		return false;
	}
	if (unmatched_file > 0)
	{
		printf("WARNING: warm_start_flows: %d snapshot links carrying volume have no matching link in this network (network edited?); their flow is dropped.\n",
			unmatched_file);
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "WARNING: warm_start_flows: %d snapshot links with volume unmatched in this network; flow dropped.\n",
				unmatched_file);
	}

	if (restore)
	{
		int unmatched_net = 0;
		for (int k = 1; k <= number_of_links; k++)
			if (!link_matched[k])
				unmatched_net++;
		if (unmatched_net > 0)
			printf("WARNING: warm_start_flows: %d network links absent from the snapshot keep their cold initialization.\n", unmatched_net);
		printf("warm_start_flows: restored flows on %d of %d links from '%s' (%d records, fingerprint match); iteration-0 all-or-nothing will be SKIPPED.\n",
			matched, number_of_links, g_warm_start_flows_file.c_str(), records);
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "warm_start_flows: restored flows on %d of %d links from '%s' (%d records, fingerprint match); iteration-0 AoN skipped.\n",
				matched, number_of_links, g_warm_start_flows_file.c_str(), records);
		if (shortest_path_log_flag)
		{
			printf("WARNING: warm_start_flows + route_output: iteration-0 paths are not stored, so route_assignment.csv prob shares cover post-warm-start iterations only.\n");
			if (summary_log_file != NULL)
				fprintf(summary_log_file, "WARNING: warm_start_flows + route_output: route_assignment prob shares cover post-warm-start iterations only.\n");
		}
		return true;
	}

	// L1 demotion: derive congested times at the snapshot volumes; volumes untouched.
	for (int k = 1; k <= number_of_links; k++)
	{
		Link_Travel_Time(k, snap_vol.data());  // writes Link[k].Travel_time
		Link[k].GenCost = Link[k].mode_AdditionalCost[1] + Link[k].Travel_time;
	}
	printf("warm_start_flows (L1 fallback): congested times derived from %d matched snapshot links; cold all-or-nothing start.\n", matched);
	if (summary_log_file != NULL)
		fprintf(summary_log_file, "warm_start_flows (L1 fallback): congested times derived from %d matched snapshot links; cold AoN start.\n", matched);
	return false;
}

// ---- P2 increment 1: DTAC binary column store (settings 'column_output'=1) ----
// route_columns.bin stores, per (mode, origin, destination) with positive
// demand and a feasible route, the LAST-ITERATION shortest path -- i.e. the
// FINAL all-or-nothing direction whose Frank-Wolfe step produced the reported
// volumes -- as an external-link-id sequence, back-traced from MinPathPredLink
// exactly like All_or_Nothing_Assign / route_assignment.csv (same mpred /
// prev_k walk, same hop-cap guards). ROUTE_ASSIGNMENT PARITY ONLY: increment 1
// does NOT store the full multi-path column set with theta shares (that is
// increment 2 / the L3 warm start), so replaying these single paths against
// the OD table reproduces the final AoN direction, not the FW-blended volumes.
// Layout: header {int32 magic 0x43415444 'DTAC', int32 version=1,
// int32 n_modes, int32 n_zones}; then for each mode m=1..n_modes and origin
// O=1..n_zones (row-major, every (m,O) present) one sparse-CSR block:
// {int32 n_dest; int32 dest_zone_id[n_dest]; int32 offsets[n_dest+1];
// int32 external_link_ids[offsets[n_dest]]} with links in origin->destination
// order. Reader/verifier: dtalite_qa/columns.py.
static void WriteColumnsDTAC(const char* filename, int*** MinPathPredLink)
{
	FILE* f = NULL;
	fopen_s(&f, filename, "wb");
	if (f == NULL)
	{
		printf("WARNING: column_output=1 but '%s' cannot be opened for writing; column store skipped.\n", filename);
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "WARNING: column_output=1 but '%s' cannot be opened; column store skipped.\n", filename);
		return;
	}
	int hdr[4] = { 0x43415444 /*'DTAC'*/, 1, number_of_modes, no_zones };
	fwrite(hdr, sizeof(int), 4, f);

	long long total_paths = 0, total_links = 0, skipped = 0;
	std::vector<int> dests, offsets, links_out, path;
	for (int m = 1; m <= number_of_modes; m++)
	{
		// mode slot holding this mode's predecessor tree: legacy
		// dedicated_shortest_path=0 modes share mode 1's tree; SP-equivalent
		// modes share their group representative's (g_rep_mode) -- identical
		// to the All_or_Nothing_Assign back-trace.
		int mpred = (g_mode_type_vector[m].dedicated_shortest_path == 0)
			? 1 : ((g_rep_mode[m] >= 1) ? g_rep_mode[m] : m);
		int cost_mode = (g_mode_type_vector[m].dedicated_shortest_path == 0) ? 1 : m;
		for (int Orig = 1; Orig <= no_zones; Orig++)
		{
			dests.clear();
			offsets.clear();
			links_out.clear();
			offsets.push_back(0);
			bool orig_ok = (TotalOFlow[Orig] >= 0.00001 && zone_outbound_link_size[Orig] != 0);
			if (orig_ok)
			{
				int internal_origin = Orig;   // dense zone seq == node seq
				for (int Dest = 1; Dest <= no_zones; Dest++)
				{
					if (Dest == Orig)
						continue;
					double od_demand = MDODflow[m][Orig][Dest];
					if (od_demand <= 0.0 ||
						(g_auto_calibration != 0 &&
						 od_demand < g_calibration_minimum_column_od_volume))
						continue;
					if (MDRouteCost[cost_mode][Orig][Dest] >= BIGM - 1)  // infeasible
						continue;

					path.clear();
					int CurrentNode = Dest;   // dense zone seq == node seq
					int prev_k = INVALID;
					int _trace_hops = 0;
					int _hop_cap = g_has_movement_restrictions ? (number_of_links + 2) : (no_nodes + 2);
					bool ok = true;
					while (CurrentNode != internal_origin)
					{
						if (++_trace_hops > _hop_cap) { ok = false; break; }
						int k;
						if (g_has_movement_restrictions)
							k = (prev_k == INVALID) ? MinPathPredLink[mpred][Orig][CurrentNode]
							: g_MinPathPredLink_link[mpred][Orig][prev_k];
						else
							k = MinPathPredLink[mpred][Orig][CurrentNode];
						if (k <= 0 || k > number_of_links) { ok = false; break; }
						path.push_back(k);
						CurrentNode = Link[k].internal_from_node_id;
						prev_k = k;
						if (CurrentNode <= 0 || CurrentNode > no_nodes) { ok = false; break; }
					}
					if (!ok || path.empty())
					{
						skipped++;
						continue;
					}
					dests.push_back(g_zone_seq_2_old_zone_id[Dest]);   // file stores EXTERNAL zone ids
					for (int i = (int)path.size() - 1; i >= 0; --i)  // reverse: origin -> destination
						links_out.push_back(Link[path[i]].external_link_id);
					offsets.push_back((int)links_out.size());
				}
			}
			int n_dest = (int)dests.size();
			fwrite(&n_dest, sizeof(int), 1, f);
			if (n_dest > 0)
			{
				fwrite(dests.data(), sizeof(int), n_dest, f);
				fwrite(offsets.data(), sizeof(int), n_dest + 1, f);
				fwrite(links_out.data(), sizeof(int), links_out.size(), f);
			}
			total_paths += n_dest;
			total_links += (long long)links_out.size();
		}
	}
	fclose(f);
	printf("column_output: wrote %s (DTAC v1, increment 1: last-iteration AoN path per OD) -- %lld paths, %lld link entries, %lld back-traces skipped\n",
		filename, total_paths, total_links, skipped);
	if (summary_log_file != NULL)
		fprintf(summary_log_file, "column_output: wrote %s (DTAC v1, increment 1) -- %lld paths, %lld link entries, %lld skipped\n",
			filename, total_paths, total_links, skipped);
}

// ---- P2 increment 2 (L3 full): theta-share column pool ----------------------
// In-memory sparse-CSR-style pool of the path SET per (mode, origin,
// destination). This is a PARALLEL mechanism to the 5D linkIndices route store
// (route_output) -- that store is untouched. The pool is active when
// column_output=1 (accumulation across FW iterations) and/or warm_start_columns
// is set (loaded columns for scatter / adjustment sweeps).
//
// Theta bookkeeping uses the exact cascade computeTheta() applies to
// route_assignment.csv: after a FW step of size lambda every stored share is
// scaled by (1-lambda) and the new all-or-nothing path enters with share
// lambda. The scaling is done in O(1) through the global multiplier g_col_mult
// (effective theta = raw_theta * g_col_mult); duplicate paths merge by summing.
struct ColPath
{
	std::vector<int> links;   // INTERNAL link ids, origin -> destination order
	double raw_theta;
};
struct ColOD
{
	int dest;
	std::vector<ColPath> paths;
};
static std::vector<std::vector<ColOD>> g_col_pool;  // row = (m-1)*no_zones + (Orig-1), ColODs sorted by dest
static bool g_col_pool_active = false;
static double g_col_mult = 1.0;

static inline std::vector<ColOD>& ColPoolRow(int m, int Orig)
{
	return g_col_pool[(size_t)(m - 1) * no_zones + (Orig - 1)];
}

// find (or insert, keeping the row sorted by dest) the ColOD of a destination
static ColOD* ColPoolFind(std::vector<ColOD>& row, int Dest, bool insert)
{
	std::vector<ColOD>::iterator it = std::lower_bound(row.begin(), row.end(), Dest,
		[](const ColOD& a, int d) { return a.dest < d; });
	if (it != row.end() && it->dest == Dest)
		return &(*it);
	if (!insert)
		return NULL;
	it = row.insert(it, ColOD());
	it->dest = Dest;
	return &(*it);
}

static void InitColumnPool()
{
	g_col_pool.assign((size_t)number_of_modes * no_zones, std::vector<ColOD>());
	g_col_mult = 1.0;
	g_col_pool_active = true;
}

static void FreeColumnPool()
{
	// Each (mode, origin) row owns independent path allocations.  Regional
	// calibration can contain more than 100 million paths, so release rows in
	// parallel instead of serially walking every small vector on the API thread.
#ifdef _OPENMP
#pragma omp parallel for num_threads(g_number_of_processors) schedule(dynamic, 1)
#endif
	for (int row = 0; row < static_cast<int>(g_col_pool.size()); row++)
		std::vector<ColOD>().swap(g_col_pool[row]);
	std::vector<std::vector<ColOD>>().swap(g_col_pool);
	g_col_pool_active = false;
	g_col_mult = 1.0;
}

// fold g_col_mult into the raw thetas so raw == effective (used before the
// adjustment sweeps, which edit thetas in place, and on underflow guard)
static void ColPoolFlatten()
{
	if (g_col_mult == 1.0)
		return;
	for (size_t r = 0; r < g_col_pool.size(); r++)
		for (size_t d = 0; d < g_col_pool[r].size(); d++)
			for (size_t j = 0; j < g_col_pool[r][d].paths.size(); j++)
				g_col_pool[r][d].paths[j].raw_theta *= g_col_mult;
	g_col_mult = 1.0;
}

// Back-trace the current shortest path (m-tree slot 'mpred', Orig -> Dest) from
// MinPathPredLink into 'path' (internal link ids, origin->destination order).
// Same mpred / g_rep_mode sharing, same link-state walk under turn restrictions
// and the same hop-cap guards as All_or_Nothing_Assign. false = malformed tree.
static bool TraceMinPath(int mpred, int Orig, int Dest, int*** MinPathPredLink, std::vector<int>& path)
{
	path.clear();
	int internal_origin = Orig;   // dense zone seq == node seq (internal renumbering)
	int CurrentNode = Dest;
	int prev_k = INVALID;
	int hops = 0;
	int hop_cap = g_has_movement_restrictions ? (number_of_links + 2) : (no_nodes + 2);
	while (CurrentNode != internal_origin)
	{
		if (++hops > hop_cap)
			return false;
		int k;
		if (g_has_movement_restrictions)
			k = (prev_k == INVALID) ? MinPathPredLink[mpred][Orig][CurrentNode]
			: g_MinPathPredLink_link[mpred][Orig][prev_k];
		else
			k = MinPathPredLink[mpred][Orig][CurrentNode];
		if (k <= 0 || k > number_of_links)
			return false;
		path.push_back(k);
		CurrentNode = Link[k].internal_from_node_id;
		prev_k = k;
		if (CurrentNode <= 0 || CurrentNode > no_nodes)
			return false;
	}
	if (path.empty())
		return false;
	std::reverse(path.begin(), path.end());
	return true;
}

// Fold the CURRENT shortest-path trees into the column pool with FW step
// weight 'lambda'. Called right after the line search of every iteration
// (lambda = the step actually taken) and after the iteration-0 all-or-nothing
// (lambda = 1.0). No-op unless column accumulation is on (column_output=1).
static void ColumnPoolUpdate(double lambda, int*** MinPathPredLink)
{
	if (!g_col_pool_active || g_column_output < 2)
		return;
	if (lambda <= 0.0)
		return;   // zero step: shares unchanged, the new AoN direction carries no mass

	if (lambda >= 1.0 - 1e-12)
	{
		// full step: every previous share vanishes exactly
		for (size_t r = 0; r < g_col_pool.size(); r++)
			g_col_pool[r].clear();
		g_col_mult = 1.0;
		lambda = 1.0;
	}
	else
	{
		g_col_mult *= (1.0 - lambda);
		if (g_col_mult < 1e-200)
			ColPoolFlatten();  // underflow guard (rare)
	}
	double new_raw = lambda / g_col_mult;

#pragma omp parallel for
	for (int p = 0; p < g_number_of_processors; p++)
	{
		std::vector<int> path;
		for (size_t i = 0; i < Processor_origin_zones[p].size(); i++)
		{
			int Orig = Processor_origin_zones[p][i];
			if (TotalOFlow[Orig] < 0.00001 || zone_outbound_link_size[Orig] == 0)
				continue;
			for (int m = 1; m <= number_of_modes; m++)
			{
				int mpred = (g_mode_type_vector[m].dedicated_shortest_path == 0)
					? 1 : ((g_rep_mode[m] >= 1) ? g_rep_mode[m] : m);
				int cost_mode = (g_mode_type_vector[m].dedicated_shortest_path == 0) ? 1 : m;
				std::vector<ColOD>& row = ColPoolRow(m, Orig);
				for (int Dest = 1; Dest <= no_zones; Dest++)
				{
					if (Dest == Orig)
						continue;
					if (MDODflow[m][Orig][Dest] <= 0.0)
						continue;
					if (MDRouteCost[cost_mode][Orig][Dest] >= BIGM - 1)  // infeasible
						continue;
					if (!TraceMinPath(mpred, Orig, Dest, MinPathPredLink, path))
						continue;
					ColOD* od = ColPoolFind(row, Dest, true);
					bool merged = false;
					for (size_t j = 0; j < od->paths.size(); j++)
						if (od->paths[j].links == path)
						{
							od->paths[j].raw_theta += new_raw;
							merged = true;
							break;
						}
					if (!merged)
					{
						od->paths.push_back(ColPath());
						od->paths.back().links = path;
						od->paths.back().raw_theta = new_raw;
					}
					if (g_auto_calibration != 0 && !od->paths.empty())
					{
						std::stable_sort(od->paths.begin(), od->paths.end(),
							[](const ColPath& left, const ColPath& right) {
								return left.raw_theta > right.raw_theta;
							});
						double raw_sum = 0.0;
						for (size_t j = 0; j < od->paths.size(); j++)
							raw_sum += fmax(0.0, od->paths[j].raw_theta);
						if (raw_sum > 0.0)
						{
							size_t keep = 0;
							while (keep < od->paths.size() &&
								keep < static_cast<size_t>(
									g_calibration_maximum_column_paths_per_od) &&
								(keep == 0 || od->paths[keep].raw_theta / raw_sum >=
									g_calibration_minimum_column_path_share))
								keep++;
							od->paths.resize(std::max<size_t>(1, keep));
							double kept_raw_sum = 0.0;
							for (size_t j = 0; j < od->paths.size(); j++)
								kept_raw_sum += fmax(0.0, od->paths[j].raw_theta);
							const double target_raw_sum = 1.0 / g_col_mult;
							const double scale = target_raw_sum /
								fmax(kept_raw_sum, 1e-300);
							for (size_t j = 0; j < od->paths.size(); j++)
								od->paths[j].raw_theta *= scale;
						}
					}
				}
			}
		}
	}
}

// DTAC v2 writer (settings 'column_output'=1, increment 2 / L3 full): the FULL
// column set with theta shares. Layout: header {int32 magic 0x43415444 'DTAC',
// int32 version=2, int32 n_modes, int32 n_zones}, then the 16-byte demand
// fingerprint (4 floats, ComputeDemandFingerprint); then for each mode
// m=1..n_modes and origin O=1..n_zones (row-major, every (m,O) present) one
// block: {int32 n_dest; int32 dest_zone_id[n_dest]; int32 path_offsets[n_dest+1]
// (path_offsets[0]==0, n_paths = path_offsets[n_dest]); float theta[n_paths];
// int32 link_offsets[n_paths+1]; int32 external_link_ids[link_offsets[n_paths]]
// origin->destination order}. Thetas per (mode,O,D) sum to 1: they are
// renormalized at write (the pre-normalization drift is reported -- it is fp
// noise for a cold full run, and a REAL, loudly-warned share deficit when the
// run was warm-started from flows so pre-warm-start iterations have no paths).
static void WriteColumnsDTACv2(const char* filename)
{
	FILE* f = NULL;
	fopen_s(&f, filename, "wb");
	if (f == NULL)
	{
		printf("WARNING: column_output=1 but '%s' cannot be opened for writing; column store skipped.\n", filename);
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "WARNING: column_output=1 but '%s' cannot be opened; column store skipped.\n", filename);
		return;
	}
	int hdr[4] = { 0x43415444 /*'DTAC'*/, 2, number_of_modes, no_zones };
	fwrite(hdr, sizeof(int), 4, f);
	float fp[4];
	ComputeDemandFingerprint(fp);
	fwrite(fp, sizeof(float), 4, f);

	long long total_ods = 0, total_paths = 0, total_links = 0, dropped_tiny = 0;
	double max_sum_dev = 0.0;   // max |sum(theta)-1| before renormalization
	std::vector<int> dests, path_offsets, link_offsets, links_out;
	std::vector<float> thetas;
	for (int m = 1; m <= number_of_modes; m++)
	{
		for (int Orig = 1; Orig <= no_zones; Orig++)
		{
			dests.clear(); path_offsets.clear(); link_offsets.clear();
			links_out.clear(); thetas.clear();
			path_offsets.push_back(0);
			link_offsets.push_back(0);
			std::vector<ColOD>& row = ColPoolRow(m, Orig);
			for (size_t d = 0; d < row.size(); d++)
			{
				ColOD& od = row[d];
				double sum = 0.0;
				for (size_t j = 0; j < od.paths.size(); j++)
					sum += od.paths[j].raw_theta * g_col_mult;
				if (sum <= 0.0 || od.paths.empty())
					continue;
				double dev = fabs(sum - 1.0);
				if (dev > max_sum_dev)
					max_sum_dev = dev;
				int n_written = 0;
				for (size_t j = 0; j < od.paths.size(); j++)
				{
					double eff = od.paths[j].raw_theta * g_col_mult / sum;  // renormalized
					if (eff < 1e-7)   // negligible share: drop (mass goes to the kept paths on the reader side via sum-to-1)
					{
						dropped_tiny++;
						continue;
					}
					thetas.push_back((float)eff);
					for (size_t i = 0; i < od.paths[j].links.size(); i++)
						links_out.push_back(Link[od.paths[j].links[i]].external_link_id);
					link_offsets.push_back((int)links_out.size());
					n_written++;
				}
				if (n_written == 0)
				{
					// all shares tiny (cannot happen with sum>0 unless fp weirdness): keep the largest
					size_t jbest = 0;
					for (size_t j = 1; j < od.paths.size(); j++)
						if (od.paths[j].raw_theta > od.paths[jbest].raw_theta)
							jbest = j;
					thetas.push_back(1.0f);
					for (size_t i = 0; i < od.paths[jbest].links.size(); i++)
						links_out.push_back(Link[od.paths[jbest].links[i]].external_link_id);
					link_offsets.push_back((int)links_out.size());
					n_written = 1;
				}
				dests.push_back(g_zone_seq_2_old_zone_id[od.dest]);   // file stores EXTERNAL zone ids
				path_offsets.push_back((int)thetas.size());
			}
			int n_dest = (int)dests.size();
			fwrite(&n_dest, sizeof(int), 1, f);
			if (n_dest > 0)
			{
				fwrite(dests.data(), sizeof(int), n_dest, f);
				fwrite(path_offsets.data(), sizeof(int), n_dest + 1, f);
				fwrite(thetas.data(), sizeof(float), thetas.size(), f);
				fwrite(link_offsets.data(), sizeof(int), thetas.size() + 1, f);
				fwrite(links_out.data(), sizeof(int), links_out.size(), f);
			}
			total_ods += n_dest;
			total_paths += (long long)thetas.size();
			total_links += (long long)links_out.size();
		}
	}
	fclose(f);
	printf("column_output: wrote %s (DTAC v2, theta-share column set) -- %lld OD pairs, %lld paths, %lld link entries; %lld negligible paths (<1e-7) dropped; max |sum(theta)-1| before renormalization = %.3g; demand fingerprint [%.6g, %.0f, %.0f, %.0f]\n",
		filename, total_ods, total_paths, total_links, dropped_tiny, max_sum_dev, fp[0], fp[1], fp[2], fp[3]);
	if (summary_log_file != NULL)
		fprintf(summary_log_file, "column_output: wrote %s (DTAC v2) -- %lld OD pairs, %lld paths, %lld link entries; %lld tiny dropped; max |sum(theta)-1| pre-norm = %.3g\n",
			filename, total_ods, total_paths, total_links, dropped_tiny, max_sum_dev);
	if (max_sum_dev > 1e-3)
	{
		printf("WARNING: column_output: theta shares deviated from 1 by up to %.3g before renormalization -- expected when the run was warm-started (pre-warm iterations carry no stored paths); shares cover the stored iterations only.\n", max_sum_dev);
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "WARNING: column theta shares deviated from 1 by up to %.3g pre-normalization (warm start?).\n", max_sum_dev);
	}
}

// L3 warm start from COLUMNS (settings 'warm_start_columns', default "" = off).
// Loads a DTAC column store (v1 single-path or v2 theta-share), keeps it in the
// in-memory pool, and scatters theta x CURRENT OD volume onto
// MainVolume/mode_MainVolume. The CURRENT demand may DIFFER from the stored
// fingerprint -- that is the point (rescale mode): theta over a path set is a
// demand-invariant routing policy. Both fingerprints and the coverage % are
// printed. Paths crossing links that no longer exist, are banned for the mode,
// break the node chain, or violate a movement restriction are DROPPED with
// theta renormalized over the survivors (counted + warned). OD pairs with
// positive demand but no usable stored column fall back to a one-shot
// all-or-nothing on the warm times implied by the scattered flows (those paths
// join the pool with theta=1). Returns true when the scatter happened: the
// caller then skips the iteration-0 all-or-nothing and Frank-Wolfe continues
// from iteration 1 (self-healing). On unreadable/incompatible input it warns
// LOUDLY and returns false (cold start) -- never fails the run.
static bool ApplyWarmStartColumns(double* MainVolume, int*** MDMinPathPredLink)
{
	if (g_warm_start_columns_file.empty())
		return false;

	auto col_load_t0 = std::chrono::high_resolution_clock::now();  // SETUP: DTAC file read + pool build

	FILE* f = NULL;
	fopen_s(&f, g_warm_start_columns_file.c_str(), "rb");
	if (f == NULL)
	{
		printf("WARNING: warm_start_columns file '%s' cannot be opened; starting cold.\n",
			g_warm_start_columns_file.c_str());
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "WARNING: warm_start_columns file '%s' cannot be opened; starting cold.\n",
				g_warm_start_columns_file.c_str());
		return false;
	}
	int hdr[4] = { 0, 0, 0, 0 };
	if (fread(hdr, sizeof(int), 4, f) != 4 || hdr[0] != 0x43415444 /*'DTAC'*/)
	{
		printf("WARNING: warm_start_columns file '%s' is not a DTAC column store (magic mismatch); starting cold.\n",
			g_warm_start_columns_file.c_str());
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "WARNING: warm_start_columns file '%s' not DTAC; starting cold.\n",
				g_warm_start_columns_file.c_str());
		fclose(f);
		return false;
	}
	int version = hdr[1], n_modes_in = hdr[2], n_zones_in = hdr[3];
	if (version != 1 && version != 2)
	{
		printf("WARNING: warm_start_columns '%s' has unknown DTAC version %d (this kernel reads v1/v2); starting cold.\n",
			g_warm_start_columns_file.c_str(), version);
		fclose(f);
		return false;
	}
	if (n_modes_in != number_of_modes)
	{
		printf("WARNING: warm_start_columns '%s' stores %d modes but the current run has %d -- mode semantics unknown; starting cold.\n",
			g_warm_start_columns_file.c_str(), n_modes_in, number_of_modes);
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "WARNING: warm_start_columns mode-count mismatch (%d vs %d); starting cold.\n",
				n_modes_in, number_of_modes);
		fclose(f);
		return false;
	}
	float fp_in[4] = { 0, 0, 0, 0 };
	bool has_fp = (version >= 2);
	if (has_fp && fread(fp_in, sizeof(float), 4, f) != 4)
	{
		printf("WARNING: warm_start_columns '%s' truncated in the header; starting cold.\n",
			g_warm_start_columns_file.c_str());
		fclose(f);
		return false;
	}
	float fp_now[4];
	ComputeDemandFingerprint(fp_now);
	if (has_fp)
	{
		bool fp_match = (fp_in[0] == fp_now[0] && fp_in[1] == fp_now[1] &&
			fp_in[2] == fp_now[2] && fp_in[3] == fp_now[3]);
		printf("warm_start_columns: stored demand fingerprint [%.6g, %.0f, %.0f, %.0f], CURRENT [%.6g, %.0f, %.0f, %.0f] -- %s (columns are rescaled to the CURRENT demand: theta is a demand-invariant policy, so a difference is expected and allowed).\n",
			fp_in[0], fp_in[1], fp_in[2], fp_in[3], fp_now[0], fp_now[1], fp_now[2], fp_now[3],
			fp_match ? "MATCH" : "DIFFER");
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "warm_start_columns: stored fp [%.6g, %.0f, %.0f, %.0f] vs current [%.6g, %.0f, %.0f, %.0f] (%s; rescale mode).\n",
				fp_in[0], fp_in[1], fp_in[2], fp_in[3], fp_now[0], fp_now[1], fp_now[2], fp_now[3],
				fp_match ? "match" : "differ");
	}
	else
		printf("warm_start_columns: DTAC v1 file (no stored demand fingerprint, single path per OD, theta=1); current demand fingerprint [%.6g, %.0f, %.0f, %.0f].\n",
			fp_now[0], fp_now[1], fp_now[2], fp_now[3]);
	if (n_zones_in != no_zones)
		printf("WARNING: warm_start_columns: stored n_zones %d != current %d; matching by zone id, out-of-range zones skipped.\n",
			n_zones_in, no_zones);

	if (!g_col_pool_active)
		InitColumnPool();
	ColPoolFlatten();

	long long paths_read = 0, paths_kept = 0, paths_dropped = 0, ods_renormalized = 0;
	long long stored_ods_read = 0, stored_ods_no_demand = 0;
	bool truncated = false;
	std::vector<int> dest_buf, poff_buf, loff_buf, link_buf;
	std::vector<float> theta_buf;
	std::vector<int> path_internal;
	for (int m = 1; m <= n_modes_in && !truncated; m++)
	{
		for (int Orig = 1; Orig <= n_zones_in && !truncated; Orig++)
		{
			int n_dest = 0;
			if (fread(&n_dest, sizeof(int), 1, f) != 1) { truncated = true; break; }
			if (n_dest < 0 || n_dest > n_zones_in) { truncated = true; break; }
			if (n_dest == 0)
				continue;
			dest_buf.resize(n_dest);
			if (fread(dest_buf.data(), sizeof(int), n_dest, f) != (size_t)n_dest) { truncated = true; break; }
			int n_paths_blk = 0;
			if (version == 1)
			{
				// v1 block: offsets[n_dest+1] then links; one path per dest, theta=1
				poff_buf.resize(n_dest + 1);
				if (fread(poff_buf.data(), sizeof(int), n_dest + 1, f) != (size_t)(n_dest + 1)) { truncated = true; break; }
				n_paths_blk = n_dest;
				theta_buf.assign(n_dest, 1.0f);
				loff_buf = poff_buf;                 // link offsets == v1 offsets
				for (int j = 0; j <= n_dest; j++)
					poff_buf[j] = j;                 // path offsets: 1 path per dest
			}
			else
			{
				poff_buf.resize(n_dest + 1);
				if (fread(poff_buf.data(), sizeof(int), n_dest + 1, f) != (size_t)(n_dest + 1)) { truncated = true; break; }
				n_paths_blk = poff_buf[n_dest];
				if (n_paths_blk < 0) { truncated = true; break; }
				theta_buf.resize(n_paths_blk);
				loff_buf.resize(n_paths_blk + 1);
				if (n_paths_blk > 0 && fread(theta_buf.data(), sizeof(float), n_paths_blk, f) != (size_t)n_paths_blk) { truncated = true; break; }
				if (fread(loff_buf.data(), sizeof(int), n_paths_blk + 1, f) != (size_t)(n_paths_blk + 1)) { truncated = true; break; }
			}
			int n_links_blk = loff_buf.empty() ? 0 : loff_buf[n_paths_blk];
			if (n_links_blk < 0) { truncated = true; break; }
			link_buf.resize(n_links_blk);
			if (n_links_blk > 0 && fread(link_buf.data(), sizeof(int), n_links_blk, f) != (size_t)n_links_blk) { truncated = true; break; }

			if (Orig < 1 || Orig > no_zones)
				continue;   // origin block outside this network's zone range
			int internal_origin = Orig;   // positional origin == dense zone seq
			std::vector<ColOD>& row = ColPoolRow(m, Orig);
			for (int d = 0; d < n_dest; d++)
			{
				int Dest = ZoneSeqFromExternal(dest_buf[d]);   // file stores EXTERNAL zone ids
				stored_ods_read++;
				if (Dest == 0 || Dest == Orig)
					continue;
				if (MDODflow[m][Orig][Dest] <= 0.0)
				{
					stored_ods_no_demand++;
					continue;   // stored policy for an OD the current demand does not use
				}
				int internal_dest = Dest;   // dense zone seq == node seq
				ColOD* od = NULL;
				int kept_here = 0, dropped_here = 0;
				for (int j = poff_buf[d]; j < poff_buf[d + 1]; j++)
				{
					paths_read++;
					path_internal.clear();
					bool ok = (loff_buf[j + 1] > loff_buf[j]);
					int prev_k = INVALID, chain_node = internal_origin;
					for (int i = loff_buf[j]; ok && i < loff_buf[j + 1]; i++)
					{
						int k = LinkSeqFromExternal(link_buf[i]);
						if (k == 0) { ok = false; break; }                            // link no longer exists
						if (Link[k].mode_allowed_use[m] == 0) { ok = false; break; }   // now banned for the mode
						if (Link[k].internal_from_node_id != chain_node) { ok = false; break; }  // broken chain
						if (g_has_movement_restrictions && prev_k != INVALID &&
							Link[prev_k].b_withmovement_restrictions && IsMovementRestricted(prev_k, k))
						{
							ok = false; break;                                          // now-banned movement
						}
						path_internal.push_back(k);
						chain_node = Link[k].internal_to_node_id;
						prev_k = k;
					}
					if (ok && chain_node != internal_dest)
						ok = false;
					if (!ok)
					{
						paths_dropped++;
						dropped_here++;
						continue;
					}
					if (od == NULL)
						od = ColPoolFind(row, Dest, true);
					bool merged = false;
					for (size_t jj = 0; jj < od->paths.size(); jj++)
						if (od->paths[jj].links == path_internal)
						{
							od->paths[jj].raw_theta += (double)theta_buf[j];
							merged = true;
							break;
						}
					if (!merged)
					{
						od->paths.push_back(ColPath());
						od->paths.back().links = path_internal;
						od->paths.back().raw_theta = (double)theta_buf[j];
					}
					paths_kept++;
					kept_here++;
				}
				if (kept_here > 0 && dropped_here > 0)
					ods_renormalized++;   // survivors renormalized at scatter time (theta/sum)
			}
		}
	}
	fclose(f);
	if (truncated)
		printf("WARNING: warm_start_columns '%s' is truncated/corrupt; using the %lld paths read so far.\n",
			g_warm_start_columns_file.c_str(), paths_kept);
	if (paths_dropped > 0)
	{
		printf("WARNING: warm_start_columns: DROPPED %lld of %lld stored paths (deleted links / mode bans / broken chains / banned movements); theta renormalized over the survivors for %lld OD pairs.\n",
			paths_dropped, paths_read, ods_renormalized);
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "WARNING: warm_start_columns dropped %lld of %lld stored paths; theta renormalized for %lld ODs.\n",
				paths_dropped, paths_read, ods_renormalized);
	}
	if (paths_kept == 0)
	{
		printf("WARNING: warm_start_columns: NO usable stored path matches this network/demand; starting cold.\n");
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "WARNING: warm_start_columns: 0 usable paths; starting cold.\n");
		return false;
	}

	{
		std::chrono::duration<double> load_dur = std::chrono::high_resolution_clock::now() - col_load_t0;
		printf("phase timing: DTAC column load + pool build (SETUP) = %.2f s\n", load_dur.count());
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "phase timing: DTAC column load + pool build (SETUP) = %.2f s\n", load_dur.count());
	}
	auto col_scatter_t0 = std::chrono::high_resolution_clock::now();  // COMPUTE: scatter + uncovered AoN

	// ---- scatter theta x CURRENT OD volume + find uncovered OD pairs ----
	for (int k = 1; k <= number_of_links; k++)
		for (int m = 1; m <= number_of_modes; m++)
			Link[k].mode_MainVolume[m] = Link[k].mode_Base_demand_volume[m];
	// (MainVolume[] already holds Base_demand_volume from the cold initialization)

	long long covered_pairs = 0, uncovered_pairs = 0, unroutable_pairs = 0;
	double covered_demand = 0.0, uncovered_demand = 0.0, unroutable_demand = 0.0;
	std::vector<int> unc_m, unc_o, unc_d;    // uncovered OD list for the one-shot AoN
	for (int m = 1; m <= number_of_modes; m++)
	{
		for (int Orig = 1; Orig <= no_zones; Orig++)
		{
			bool orig_ok = (TotalOFlow[Orig] >= 0.00001 && zone_outbound_link_size[Orig] != 0);
			std::vector<ColOD>& row = ColPoolRow(m, Orig);
			for (int Dest = 1; Dest <= no_zones; Dest++)
			{
				if (Dest == Orig)
					continue;
				double T = MDODflow[m][Orig][Dest];
				if (T <= 0.0)
					continue;
				if (!orig_ok)
				{
					unroutable_pairs++;
					unroutable_demand += T;
					continue;   // same ODs All_or_Nothing_Assign drops
				}
				ColOD* od = ColPoolFind(row, Dest, false);
				double sum = 0.0;
				if (od != NULL)
					for (size_t j = 0; j < od->paths.size(); j++)
						sum += od->paths[j].raw_theta;
				if (od == NULL || sum <= 0.0)
				{
					uncovered_pairs++;
					uncovered_demand += T;
					unc_m.push_back(m); unc_o.push_back(Orig); unc_d.push_back(Dest);
					continue;
				}
				covered_pairs++;
				covered_demand += T;
				double pce = g_mode_type_vector[m].pce;
				for (size_t j = 0; j < od->paths.size(); j++)
				{
					double share = od->paths[j].raw_theta / sum;
					od->paths[j].raw_theta = share;      // store renormalized theta
					double flow = share * T;
					if (flow <= 0.0)
						continue;
					for (size_t i = 0; i < od->paths[j].links.size(); i++)
					{
						int k = od->paths[j].links[i];
						MainVolume[k] += flow * pce;
						Link[k].mode_MainVolume[m] += flow;
					}
				}
			}
		}
	}

	double total_demand = covered_demand + uncovered_demand + unroutable_demand;
	printf("warm_start_columns: coverage %lld of %lld routable OD pairs (%.2f%%), %.1f of %.1f trips (%.2f%%); %lld stored ODs carried no current demand (skipped); %lld OD pairs (%.1f trips) have no stored column.\n",
		covered_pairs, covered_pairs + uncovered_pairs,
		100.0 * covered_pairs / fmax(1.0, (double)(covered_pairs + uncovered_pairs)),
		covered_demand, total_demand,
		100.0 * covered_demand / fmax(1e-9, total_demand),
		stored_ods_no_demand, uncovered_pairs, uncovered_demand);
	if (summary_log_file != NULL)
		fprintf(summary_log_file, "warm_start_columns: coverage %lld/%lld OD pairs (%.2f%% of trips); %lld uncovered pairs -> one-shot AoN.\n",
			covered_pairs, covered_pairs + uncovered_pairs,
			100.0 * covered_demand / fmax(1e-9, total_demand), uncovered_pairs);
	if (unroutable_pairs > 0)
		printf("warm_start_columns: %lld OD pairs (%.1f trips) from zones with no outbound link are dropped (same as all-or-nothing).\n",
			unroutable_pairs, unroutable_demand);

	// ---- one-shot all-or-nothing on warm times for the uncovered OD pairs ----
	if (!unc_m.empty())
	{
		if (covered_demand > 0.0)
			UpdateLinkCost(MainVolume);   // warm times implied by the scattered flows
		printf("warm_start_columns: routing %lld uncovered OD pairs with a one-shot all-or-nothing on warm times...\n",
			(long long)unc_m.size());
		FindMinCostRoutes(MDMinPathPredLink);
		long long infeasible = 0;
		std::vector<int> path;
		for (size_t u = 0; u < unc_m.size(); u++)
		{
			int m = unc_m[u], Orig = unc_o[u], Dest = unc_d[u];
			int mpred = (g_mode_type_vector[m].dedicated_shortest_path == 0)
				? 1 : ((g_rep_mode[m] >= 1) ? g_rep_mode[m] : m);
			int cost_mode = (g_mode_type_vector[m].dedicated_shortest_path == 0) ? 1 : m;
			if (MDRouteCost[cost_mode][Orig][Dest] >= BIGM - 1 ||
				!TraceMinPath(mpred, Orig, Dest, MDMinPathPredLink, path))
			{
				infeasible++;
				continue;
			}
			double T = MDODflow[m][Orig][Dest];
			double pce = g_mode_type_vector[m].pce;
			for (size_t i = 0; i < path.size(); i++)
			{
				MainVolume[path[i]] += T * pce;
				Link[path[i]].mode_MainVolume[m] += T;
			}
			ColOD* od = ColPoolFind(ColPoolRow(m, Orig), Dest, true);
			od->paths.push_back(ColPath());
			od->paths.back().links = path;
			od->paths.back().raw_theta = 1.0;
		}
		if (infeasible > 0)
			printf("WARNING: warm_start_columns: %lld uncovered OD pairs are infeasible on the current network (demand dropped, same as all-or-nothing).\n",
				infeasible);
	}

	{
		std::chrono::duration<double> scatter_dur = std::chrono::high_resolution_clock::now() - col_scatter_t0;
		printf("phase timing: column scatter + uncovered one-shot AoN (COMPUTE) = %.2f s\n", scatter_dur.count());
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "phase timing: column scatter + uncovered one-shot AoN (COMPUTE) = %.2f s\n", scatter_dur.count());
	}
	printf("warm_start_columns: restored %lld paths from '%s' (DTAC v%d); iteration-0 all-or-nothing will be SKIPPED; Frank-Wolfe continues from the scattered flows.\n",
		paths_kept, g_warm_start_columns_file.c_str(), version);
	if (summary_log_file != NULL)
		fprintf(summary_log_file, "warm_start_columns: restored %lld paths (DTAC v%d); iteration-0 AoN skipped.\n",
			paths_kept, version);
	return true;
}

// Exact link-cost derivative dTravel_time/dVolume at the current volumes via a
// forward difference on the ACTUAL VDF (any vdf_type) -- same pattern as
// ApplyConjugateFW. Link_Travel_Time mutates Link[k].Travel_time, so restore it.
static double ColLinkTTDer(int k, double* MainVolume)
{
	double v = MainVolume[k];
	double keep_tt = Link[k].Travel_time;
	double t0 = Link_Travel_Time(k, MainVolume);
	double dv = fmax(1.0, fabs(v) * 1e-4);
	MainVolume[k] = v + dv;
	double t1 = Link_Travel_Time(k, MainVolume);
	MainVolume[k] = v;
	Link[k].Travel_time = keep_tt;
	return fmax(0.0, (t1 - t0) / dv);
}

// restricted gap over the STORED paths at the current Link.Travel_time (read-only)
static double ColRestrictedGap(void)
{
	double total_cost = 0.0, least_cost = 0.0;
	for (int m = 1; m <= number_of_modes; m++)
	{
		double pce = g_mode_type_vector[m].pce;
		for (int Orig = 1; Orig <= no_zones; Orig++)
		{
			std::vector<ColOD>& row = ColPoolRow(m, Orig);
			for (size_t d = 0; d < row.size(); d++)
			{
				ColOD& od = row[d];
				double T = MDODflow[m][Orig][od.dest];
				if (T <= 0.0 || od.paths.empty())
					continue;
				double c_best = 1e30;
				for (size_t j = 0; j < od.paths.size(); j++)
				{
					double cj = 0.0;
					for (size_t i = 0; i < od.paths[j].links.size(); i++)
					{
						int k = od.paths[j].links[i];
						cj += Link[k].Travel_time + Link[k].mode_AdditionalCost[m];
					}
					total_cost += od.paths[j].raw_theta * T * pce * cj;
					if (cj < c_best)
						c_best = cj;
				}
				least_cost += T * pce * c_best;
			}
		}
	}
	return (total_cost - least_cost) /
		fmax(0.1, g_relative_gap_standard ? total_cost : least_cost) * 100.0;
}

// Fixed-policy gradient-projection adjustment (settings 'column_adjust_sweeps',
// default 0): N sweeps over the loaded column pool that shift share from
// costlier stored paths to the CHEAPEST stored path of each (mode, origin,
// destination) under current costs, WITHOUT any new shortest-path call. The
// flow shift per path is the classic GP (Newton) step min(movable flow,
// cost difference / sum of exact link-cost derivatives over the two paths'
// non-shared links). Sweeps are GAUSS-SEIDEL: after each OD's shift the
// travel times of the touched links are re-evaluated from the actual VDF, so
// later ODs see fresh costs (a frozen-cost simultaneous sweep overshoots and
// diverges -- verified). The printed per-sweep gap is the RESTRICTED gap
// (best of the STORED paths as the benchmark) -- it is NOT the true
// equilibrium gap; the caller always prints the fresh-shortest-path TRUE gap
// afterwards. Sweeps are serial (panel recommendation: shared-link write
// hazards). Returns the final system-wide travel time (UpdateLinkCost applied
// after the last sweep).
static double RunColumnAdjustSweeps(double* MainVolume)
{
	ColPoolFlatten();
	std::vector<int> stamp(number_of_links + 1, 0);
	int stamp_cnt = 0;
	double system_wide = 0.0;
	for (int sweep = 1; sweep <= g_column_adjust_sweeps; sweep++)
	{
		auto sweep_t0 = std::chrono::high_resolution_clock::now();
		double rgap_before = ColRestrictedGap();
		long long shifts = 0;
		double moved_flow = 0.0;
		for (int m = 1; m <= number_of_modes; m++)
		{
			double pce = g_mode_type_vector[m].pce;
			for (int Orig = 1; Orig <= no_zones; Orig++)
			{
				std::vector<ColOD>& row = ColPoolRow(m, Orig);
				for (size_t d = 0; d < row.size(); d++)
				{
					ColOD& od = row[d];
					double T = MDODflow[m][Orig][od.dest];
					int nP = (int)od.paths.size();
					if (T <= 0.0 || nP < 2)
						continue;
					// path costs under the CURRENT (Gauss-Seidel-fresh) link costs
					double c_best = 1e30;
					int jbest = 0;
					std::vector<double> c(nP);
					for (int j = 0; j < nP; j++)
					{
						double cj = 0.0;
						for (size_t i = 0; i < od.paths[j].links.size(); i++)
						{
							int k = od.paths[j].links[i];
							cj += Link[k].Travel_time + Link[k].mode_AdditionalCost[m];
						}
						c[j] = cj;
						if (cj < c_best) { c_best = cj; jbest = j; }
					}
					ColPath& best = od.paths[jbest];
					stamp_cnt++;
					for (size_t i = 0; i < best.links.size(); i++)
						stamp[best.links[i]] = stamp_cnt;
					bool touched = false;
					for (int j = 0; j < nP; j++)
					{
						if (j == jbest)
							continue;
						double theta_j = od.paths[j].raw_theta;
						double diff = c[j] - c_best;
						if (theta_j <= 0.0 || diff <= 1e-12)
							continue;
						// exact derivative sum over the two paths' NON-SHARED links
						double s = 0.0;
						for (size_t i = 0; i < od.paths[j].links.size(); i++)
							if (stamp[od.paths[j].links[i]] != stamp_cnt)
								s += ColLinkTTDer(od.paths[j].links[i], MainVolume);
						for (size_t i = 0; i < best.links.size(); i++)
						{
							int k = best.links[i];
							bool on_j = false;
							for (size_t ii = 0; ii < od.paths[j].links.size(); ii++)
								if (od.paths[j].links[ii] == k) { on_j = true; break; }
							if (!on_j)
								s += ColLinkTTDer(k, MainVolume);
						}
						double f_j = theta_j * T;
						double dflow = (s * pce > 1e-12) ? fmin(f_j, diff / (s * pce)) : f_j;
						if (dflow <= 0.0)
							continue;
						od.paths[j].raw_theta -= dflow / T;
						best.raw_theta += dflow / T;
						for (size_t i = 0; i < od.paths[j].links.size(); i++)
						{
							MainVolume[od.paths[j].links[i]] -= dflow * pce;
							Link[od.paths[j].links[i]].mode_MainVolume[m] -= dflow;
						}
						for (size_t i = 0; i < best.links.size(); i++)
						{
							MainVolume[best.links[i]] += dflow * pce;
							Link[best.links[i]].mode_MainVolume[m] += dflow;
						}
						shifts++;
						moved_flow += dflow;
						touched = true;
					}
					if (touched)
					{
						// Gauss-Seidel freshness: re-evaluate travel time on every
						// link whose volume just changed (Link_Travel_Time writes
						// Link[k].Travel_time), so subsequent ODs price correctly.
						for (size_t j2 = 0; j2 < od.paths.size(); j2++)
							for (size_t i = 0; i < od.paths[j2].links.size(); i++)
								Link_Travel_Time((int)od.paths[j2].links[i], MainVolume);
					}
				}
			}
		}
		system_wide = UpdateLinkCost(MainVolume);   // full cost re-sync between sweeps
		double rgap_after = ColRestrictedGap();
		std::chrono::duration<double> sweep_dur = std::chrono::high_resolution_clock::now() - sweep_t0;
		printf("column_adjust sweep %d: RESTRICTED gap (STORED paths only, NOT the true equilibrium gap) %f %% -> %f %%, %lld share shifts, %.1f trips moved, %.2f s\n",
			sweep, rgap_before, rgap_after, shifts, moved_flow, sweep_dur.count());
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "column_adjust sweep %d: RESTRICTED gap (stored paths only) %f %% -> %f %%, %lld shifts, %.1f trips moved, %.2f s\n",
				sweep, rgap_before, rgap_after, shifts, moved_flow, sweep_dur.count());
	}
	return system_wide;
}

struct AutoCalibrationParameterSnapshot
{
	std::vector<double> plf;
	std::vector<double> qcd;
	std::vector<double> qcp;
	std::vector<double> alpha;
	std::vector<double> beta;
};

static AutoCalibrationParameterSnapshot CaptureAutoCalibrationParameters()
{
	AutoCalibrationParameterSnapshot snapshot;
	snapshot.plf.resize(number_of_links + 1);
	snapshot.qcd.resize(number_of_links + 1);
	snapshot.qcp.resize(number_of_links + 1);
	snapshot.alpha.resize(number_of_links + 1);
	snapshot.beta.resize(number_of_links + 1);
	for (int k = 1; k <= number_of_links; k++)
	{
		snapshot.plf[k] = Link[k].VDF_plf;
		snapshot.qcd[k] = Link[k].Q_cd;
		snapshot.qcp[k] = Link[k].Q_cp;
		snapshot.alpha[k] = Link[k].VDF_Alpha;
		snapshot.beta[k] = Link[k].VDF_Beta;
	}
	return snapshot;
}

static void RestoreAutoCalibrationParameters(
	const AutoCalibrationParameterSnapshot& snapshot)
{
	for (int k = 1; k <= number_of_links; k++)
	{
		Link[k].VDF_plf = snapshot.plf[k];
		Link[k].Q_cd = snapshot.qcd[k];
		Link[k].Q_cp = snapshot.qcp[k];
		Link[k].VDF_Alpha = snapshot.alpha[k];
		Link[k].VDF_Beta = snapshot.beta[k];
		Link[k].BoverC = Link[k].Link_Capacity > 0.0
			? Link[k].VDF_Alpha /
				pow(Link[k].Link_Capacity, Link[k].VDF_Beta)
			: 0.0;
	}
}

static std::vector<taplite::CalibrationLink> BuildAutoCalibrationLinks(
	const double* MainVolume)
{
	std::vector<taplite::CalibrationLink> links(number_of_links);
	for (int k = 1; k <= number_of_links; k++)
	{
		taplite::CalibrationLink& value = links[k - 1];
		value.internal_index = k - 1;
		value.external_link_id = Link[k].external_link_id;
		value.from_node_id = Link[k].external_from_node_id;
		value.to_node_id = Link[k].external_to_node_id;
		value.vdf_code = Link[k].link_type;
		value.corridor = Link[k].Calibration_corridor;
		value.facility_class = Link[k].Calibration_facility_class;
		value.target_tmc = Link[k].Calibration_target_tmc;
		value.calibration_exclusion_reason =
			Link[k].Calibration_exclusion_reason;
		value.calibration_eligible = Link[k].Calibration_eligible;
		value.observation_class = taplite::ParseObservationClass(
			Link[k].Calibration_observation_class);
		value.length_miles = Link[k].length;
		value.lanes = Link[k].lanes;
		value.capacity_vphpl = Link[k].Lane_Capacity;
		value.free_speed_mph = Link[k].free_speed;
		value.cutoff_speed_mph = Link[k].Cutoff_Speed;
		value.volume = MainVolume[k];
		value.vehicle_volume = 0.0;
		for (int m = 1; m <= number_of_modes; m++)
			value.vehicle_volume += fmax(0.0, Link[k].mode_MainVolume[m]);
		if (value.vehicle_volume <= 0.0 && value.volume > 0.0)
			value.vehicle_volume = value.volume;
		value.travel_time_minutes = Link[k].Travel_time;
		value.plf = Link[k].VDF_plf;
		value.qcd = Link[k].Q_cd;
		value.qcp = Link[k].Q_cp;
		value.qn = Link[k].Q_n;
		value.qs = Link[k].Q_s;
		value.alpha = Link[k].VDF_Alpha;
		value.beta = Link[k].VDF_Beta;
		value.mode1_plf = Link[k].Mode1_plf;
		value.mode1_qcd = Link[k].Mode1_Q_cd;
		value.mode1_qcp = Link[k].Mode1_Q_cp;
		value.observed_duration_hour = Link[k].Calibration_observed_P;
		value.observed_trough_speed_mph = Link[k].Calibration_observed_vt2;
		value.observed_average_speed_mph =
			Link[k].Calibration_observed_average_speed;
		value.s3_volume = Link[k].Calibration_S3_volume;
		value.cube_vehicle_volume = Link[k].Calibration_cube_vehicle_volume;
		value.observation_quality = Link[k].Calibration_observation_quality;
	}
	return links;
}

static std::vector<taplite::CalibrationRoute> BuildAutoCalibrationRoutes()
{
	std::vector<taplite::CalibrationRoute> routes;
	if (!g_col_pool_active || g_column_output < 2)
		return routes;
	for (int m = 1; m <= number_of_modes; m++)
	{
		for (int Orig = 1; Orig <= no_zones; Orig++)
		{
			std::vector<ColOD>& row = ColPoolRow(m, Orig);
			for (size_t d = 0; d < row.size(); d++)
			{
				ColOD& od = row[d];
				double sum = 0.0;
				for (size_t j = 0; j < od.paths.size(); j++)
					sum += od.paths[j].raw_theta * g_col_mult;
				if (sum <= 0.0 || od.dest < 1 || od.dest > no_zones)
					continue;
				double demand = MDODflow[m][Orig][od.dest];
				if (demand <= 0.0)
					continue;
				for (size_t j = 0; j < od.paths.size(); j++)
				{
					double share = od.paths[j].raw_theta * g_col_mult / sum;
					if (share <= 1e-9 || od.paths[j].links.empty())
						continue;
					taplite::CalibrationRoute route;
					route.mode_index = m;
					route.origin_zone_id = g_zone_seq_2_old_zone_id[Orig];
					route.destination_zone_id = g_zone_seq_2_old_zone_id[od.dest];
					route.mode_name = g_mode_type_vector[m].mode_type;
					route.od_demand = demand;
					route.pce = g_mode_type_vector[m].pce;
					route.share = share;
					route.link_indices.reserve(od.paths[j].links.size());
					route.link_travel_times_minutes.reserve(
						od.paths[j].links.size());
					for (size_t p = 0; p < od.paths[j].links.size(); p++)
					{
						int k = od.paths[j].links[p];
						if (k >= 1 && k <= number_of_links)
						{
							route.link_indices.push_back(k - 1);
							route.link_travel_times_minutes.push_back(
								Link[k].Travel_time);
						}
					}
					if (!route.link_indices.empty())
						routes.push_back(route);
				}
			}
		}
	}
	return routes;
}

static double AutoCalibrationUnwrapMinute(
	double minute, double start_minute, double end_minute)
{
	while (minute < start_minute)
		minute += 1440.0;
	while (minute >= start_minute + 1440.0 && minute >= end_minute)
		minute -= 1440.0;
	return minute;
}

static std::vector<std::pair<int, double> > AutoCalibrationDepartureBins(
	const taplite::AutoCalibrationEngine& engine,
	int origin_zone_id,
	const std::string& mode_name,
	double start_minute,
	double end_minute,
	double bin_minutes)
{
	std::vector<std::pair<int, double> > bins;
	const std::vector<taplite::DepartureProfilePoint>* profile =
		engine.DepartureProfileFor(origin_zone_id, mode_name);
	if (profile == NULL)
		return bins;
	double period_share = 0.0;
	for (size_t i = 0; i < profile->size(); i++)
	{
		double minute = AutoCalibrationUnwrapMinute(
			(*profile)[i].minute_of_day, start_minute, end_minute);
		if (minute >= start_minute && minute < end_minute)
			period_share += fmax(0.0, (*profile)[i].share);
	}
	if (period_share <= 0.0)
		return bins;
	for (size_t i = 0; i < profile->size(); i++)
	{
		double minute = AutoCalibrationUnwrapMinute(
			(*profile)[i].minute_of_day, start_minute, end_minute);
		if (minute < start_minute || minute >= end_minute)
			continue;
		int bin = static_cast<int>(floor(
			(minute - start_minute) / bin_minutes));
		bins.push_back(std::make_pair(
			bin, fmax(0.0, (*profile)[i].share) / period_share));
	}
	return bins;
}

// Stream DTAC-v2 paths into bounded link-delay histograms.  The former bridge
// copied every path (and every link travel time) into CalibrationRoute objects,
// then created string route keys.  On the full NVTA OD matrix that exceeded
// 120 GiB before the first outer iteration.  This version keeps column_output=2
// as the authoritative route policy but adds only O(workers*links*bins) scratch.
static taplite::CalibrationRouteSummary BuildAutoCalibrationRouteSummary(
	const taplite::AutoCalibrationEngine& engine,
	const std::vector<taplite::CalibrationLink>& links,
	bool include_departure_plf,
	const std::string& arrival_profile_output_file = "")
{
	taplite::CalibrationRouteSummary summary;
	const size_t link_count = links.size();
	summary.policy_weights.assign(
		static_cast<size_t>(number_of_modes + 1) * link_count, 0.0);
	for (int m = 1; m <= number_of_modes; m++)
	{
		double pce = fmax(0.0, g_mode_type_vector[m].pce);
		for (int k = 1; k <= number_of_links; k++)
			summary.policy_weights[static_cast<size_t>(m) * link_count + k - 1] =
				fmax(0.0, Link[k].mode_MainVolume[m]) * pce;
	}

	if (!include_departure_plf || !g_col_pool_active || g_column_output < 2)
		return summary;

	summary.departure_plf.resize(link_count, 1.0);
	if (link_count == 0)
		return summary;
	const taplite::AutoCalibrationConfig& config = engine.config();
	const double start_minute = config.period_start_hour * 60.0;
	double end_minute = config.period_end_hour * 60.0;
	if (end_minute <= start_minute)
		end_minute += 1440.0;
	const double bin_minutes = fmax(1.0, config.departure_bin_minutes);
	const double duration_minutes = end_minute - start_minute;
	// A configurable post-period tail bounds the dense arrival histogram. Trips
	// beyond it are skipped rather than forcing a route-wide pre-scan.
	const int horizon_bins = std::max(1, static_cast<int>(ceil(
		(duration_minutes + config.arrival_tail_minutes + bin_minutes) /
		bin_minutes)));
	const int workers = std::max(1, std::min(config.workers,
		g_number_of_processors));
	std::vector<double> arrival(link_count * horizon_bins, 0.0);

	for (int m = 1; m <= number_of_modes; m++)
	{
		const std::string mode_name = g_mode_type_vector[m].mode_type;
		const bool shared_mode_profile = !engine.uses_origin_specific_profiles();
		std::vector<std::pair<int, double> > mode_departures =
			AutoCalibrationDepartureBins(engine, 0, mode_name,
				start_minute, end_minute, bin_minutes);
		const size_t values_per_thread = link_count * horizon_bins;
		std::vector<double> thread_bins(
			static_cast<size_t>(workers) * values_per_thread, 0.0);
		long long path_count = 0;
		long long link_entry_count = 0;

#ifdef _OPENMP
#pragma omp parallel for num_threads(workers) schedule(dynamic, 1) reduction(+:path_count,link_entry_count)
#endif
		for (int Orig = 1; Orig <= no_zones; Orig++)
		{
			int thread_id = 0;
#ifdef _OPENMP
			thread_id = omp_get_thread_num();
#endif
			double* local = &thread_bins[
				static_cast<size_t>(thread_id) * values_per_thread];
			std::vector<std::pair<int, double> > departures = shared_mode_profile
				? mode_departures
				: AutoCalibrationDepartureBins(
					engine, g_zone_seq_2_old_zone_id[Orig], mode_name,
					start_minute, end_minute, bin_minutes);
			if (departures.empty())
				continue;
			std::vector<ColOD>& row = ColPoolRow(m, Orig);
			for (size_t d = 0; d < row.size(); d++)
			{
				ColOD& od = row[d];
				if (od.dest < 1 || od.dest > no_zones)
					continue;
				double demand = MDODflow[m][Orig][od.dest];
				if (demand <= 0.0)
					continue;
				double theta_sum = 0.0;
				for (size_t j = 0; j < od.paths.size(); j++)
					theta_sum += od.paths[j].raw_theta * g_col_mult;
				if (theta_sum <= 0.0)
					continue;
				for (size_t j = 0; j < od.paths.size(); j++)
				{
					double share = od.paths[j].raw_theta * g_col_mult / theta_sum;
					double route_flow = demand *
						fmax(0.0, g_mode_type_vector[m].pce) * share;
					if (route_flow <= 1e-9 || od.paths[j].links.empty())
						continue;
					path_count++;
					double cumulative_minutes = 0.0;
					for (size_t p = 0; p < od.paths[j].links.size(); p++)
					{
						int k = od.paths[j].links[p];
						if (k < 1 || k > number_of_links)
							continue;
						link_entry_count++;
						if (shared_mode_profile)
						{
							int delay_bin = static_cast<int>(floor(
								cumulative_minutes / bin_minutes));
							if (delay_bin >= 0 && delay_bin < horizon_bins)
								local[static_cast<size_t>(k - 1) * horizon_bins +
									delay_bin] += route_flow;
						}
						else
						{
							for (size_t q = 0; q < departures.size(); q++)
							{
								int arrival_bin = departures[q].first +
									static_cast<int>(floor(
										cumulative_minutes / bin_minutes));
								if (arrival_bin >= 0 && arrival_bin < horizon_bins)
									local[static_cast<size_t>(k - 1) * horizon_bins +
										arrival_bin] += route_flow * departures[q].second;
							}
						}
						cumulative_minutes += fmax(0.0, Link[k].Travel_time);
					}
				}
			}
		}

#ifdef _OPENMP
#pragma omp parallel for num_threads(workers) schedule(static)
#endif
		for (int link_index = 0; link_index < static_cast<int>(link_count);
			link_index++)
		{
			const size_t base = static_cast<size_t>(link_index) * horizon_bins;
			double routed_total = 0.0;
			for (int bin = 0; bin < horizon_bins; bin++)
				for (int thread = 0; thread < workers; thread++)
					routed_total += thread_bins[
						static_cast<size_t>(thread) * values_per_thread + base + bin];
			const double exact_mode_link_volume =
				fmax(0.0, Link[link_index + 1].mode_MainVolume[m]) *
				fmax(0.0, g_mode_type_vector[m].pce);
			if (exact_mode_link_volume <= 0.0)
				continue;
			const double route_scale = routed_total > 1e-12
				? exact_mode_link_volume / routed_total : 0.0;
			for (int bin = 0; bin < horizon_bins; bin++)
			{
				double value = 0.0;
				for (int thread = 0; thread < workers; thread++)
					value += thread_bins[
						static_cast<size_t>(thread) * values_per_thread + base + bin];
				value *= route_scale;
				if (value <= 0.0)
					continue;
				if (shared_mode_profile)
				{
					for (size_t q = 0; q < mode_departures.size(); q++)
					{
						int target = mode_departures[q].first + bin;
						if (target >= 0 && target < horizon_bins)
							arrival[base + target] +=
								value * mode_departures[q].second;
					}
				}
				else
					arrival[base + bin] += value;
			}
			// If no retained high-volume OD path reaches this link, preserve its
			// exact assigned volume with the mode departure shape at zero delay.
			if (routed_total <= 1e-12)
				for (size_t q = 0; q < mode_departures.size(); q++)
				{
					const int target = mode_departures[q].first;
					if (target >= 0 && target < horizon_bins)
						arrival[base + target] +=
							exact_mode_link_volume * mode_departures[q].second;
				}
		}
		printf("auto_calibration route summary: mode %s streamed %lld paths and %lld link entries.\n",
			mode_name.c_str(), path_count, link_entry_count);
		fflush(stdout);
	}

	const double period_hours = fmax(
		(end_minute - start_minute) / 60.0, 1e-6);
	const double bin_hours = bin_minutes / 60.0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(workers) schedule(static)
#endif
	for (int link_index = 0; link_index < static_cast<int>(link_count);
		link_index++)
	{
		const size_t base = static_cast<size_t>(link_index) * horizon_bins;
		double total = 0.0;
		double peak = 0.0;
		for (int bin = 0; bin < horizon_bins; bin++)
		{
			double value = arrival[base + bin];
			total += value;
			peak = fmax(peak, value);
		}
		double plf = links[link_index].plf > 0.0
			? links[link_index].plf : 1.0;
		if (total > 1e-9 && peak > 1e-9)
			plf = bin_hours / (period_hours * (peak / total));
		summary.departure_plf[link_index] = fmin(
			config.maximum_plf, fmax(config.minimum_plf, plf));
	}

	// The regional route pool is intentionally not serialized during calibration.
	// When requested, stream only canonical links (nonempty selected TMC) from the
	// accepted final pool.  Every row is a 15-minute PCE arrival bin propagated
	// along the final equilibrium paths using the final link travel times.
	if (!arrival_profile_output_file.empty())
	{
		std::ofstream output(arrival_profile_output_file.c_str());
		if (!output.is_open())
		{
			printf("ERROR: cannot write canonical link arrival profile '%s'.\n",
				arrival_profile_output_file.c_str());
		}
		else
		{
			auto csv_text = [](const std::string& value) {
				std::string escaped = value;
				size_t position = 0;
				while ((position = escaped.find('"', position)) != std::string::npos)
				{
					escaped.insert(position, 1, '"');
					position += 2;
				}
				return std::string("\"") + escaped + "\"";
			};
			auto hhmm = [](double minute) {
				int rounded = static_cast<int>(floor(minute + 0.5));
				int day_minute = ((rounded % 1440) + 1440) % 1440;
				char buffer[6];
				std::snprintf(buffer, sizeof(buffer), "%02d:%02d",
					day_minute / 60, day_minute % 60);
				return std::string(buffer);
			};
			output << "link_id,from_node_id,to_node_id,target_tmc,corridor,"
				<< "facility_class,bin_index,minutes_from_period_start,day_offset,"
				<< "time_start,time_end,within_assignment_period,arrival_pce,"
				<< "arrival_rate_pce_per_hour,arrival_share_of_link_total\n";
			output << std::setprecision(12);
			long long rows_written = 0;
			int canonical_links = 0;
			const int period_bins = std::max(1, static_cast<int>(ceil(
				duration_minutes / bin_minutes)));
			for (int link_index = 0; link_index < static_cast<int>(link_count);
				link_index++)
			{
				const taplite::CalibrationLink& link = links[link_index];
				if (link.target_tmc.empty())
					continue;
				canonical_links++;
				const size_t base = static_cast<size_t>(link_index) * horizon_bins;
				double total = 0.0;
				int last_bin = period_bins - 1;
				for (int bin = 0; bin < horizon_bins; bin++)
				{
					total += arrival[base + bin];
					if (arrival[base + bin] > 1e-12)
						last_bin = bin;
				}
				last_bin = std::min(last_bin, horizon_bins - 1);
				for (int bin = 0; bin <= last_bin; bin++)
				{
					const double minute_start = start_minute + bin * bin_minutes;
					const double minute_end = minute_start + bin_minutes;
					const double value = arrival[base + bin];
					const int day_offset = static_cast<int>(floor(minute_start / 1440.0));
					output << link.external_link_id << ',' << link.from_node_id << ','
						<< link.to_node_id << ',' << csv_text(link.target_tmc) << ','
						<< csv_text(link.corridor) << ',' << csv_text(link.facility_class)
						<< ',' << bin << ',' << bin * bin_minutes << ',' << day_offset
						<< ',' << hhmm(minute_start) << ',' << hhmm(minute_end) << ','
						<< (bin < period_bins ? 1 : 0) << ',' << value << ','
						<< value / bin_hours << ','
						<< (total > 1e-12 ? value / total : 0.0) << '\n';
					rows_written++;
				}
			}
			output.close();
			printf("auto_calibration: wrote %lld arrival bins for %d canonical links to %s.\n",
				rows_written, canonical_links, arrival_profile_output_file.c_str());
		}
	}
	return summary;
}

static void ApplyAutoCalibrationProposal(
	const taplite::CalibrationProposal& proposal)
{
	for (int k = 1; k <= number_of_links; k++)
	{
		size_t i = static_cast<size_t>(k - 1);
		Link[k].VDF_plf = proposal.plf[i];
		Link[k].Q_cd = proposal.qcd[i];
		Link[k].Q_cp = proposal.qcp[i];
		Link[k].VDF_Alpha = proposal.alpha[i];
		Link[k].VDF_Beta = proposal.beta[i];
		Link[k].BoverC = Link[k].Link_Capacity > 0.0
			? Link[k].VDF_Alpha /
				pow(Link[k].Link_Capacity, Link[k].VDF_Beta)
			: 0.0;
	}
}

static double RunAutoCalibrationInnerEquilibrium(
	double* MainVolume,
	double* SubVolume,
	double* SDVolume,
	int*** MDMinPathPredLink,
	double& system_wide_travel_time,
	int outer_iteration,
	double calibration_gap_tolerance_pct)
{
	system_wide_travel_time = UpdateLinkCost(MainVolume);
	double system_least_travel_time = system_wide_travel_time;
	double gap = std::numeric_limits<double>::infinity();
	int gap_below_count = 0;
	double Lambda = 1.0;
	for (int inner_iteration = 1;
		inner_iteration < TotalAssignIterations; inner_iteration++)
	{
		system_least_travel_time = FindMinCostRoutes(MDMinPathPredLink);
		gap = (system_wide_travel_time - system_least_travel_time) /
			fmax(0.1, g_relative_gap_standard
				? system_wide_travel_time : system_least_travel_time) * 100.0;
		printf("auto_calibration outer %d inner %d: sys. TT %.1f, least TT %.1f, gap %.6f %%\n",
			outer_iteration, inner_iteration, system_wide_travel_time,
			system_least_travel_time, gap);
		fflush(stdout);
		const double stopping_gap = calibration_gap_tolerance_pct > 0.0
			? calibration_gap_tolerance_pct : g_convergence_gap_pct;
		if (stopping_gap > 0.0 && gap >= 0.0 && gap < stopping_gap)
		{
			gap_below_count++;
			if (gap_below_count >= g_convergence_consecutive)
				break;
		}
		else
			gap_below_count = 0;

		All_or_Nothing_Assign(inner_iteration, MDODflow,
			MDMinPathPredLink, SubVolume);
		VolumeDifference(SubVolume, MainVolume, SDVolume);
		ApplyConjugateFW(MainVolume, SubVolume, SDVolume, inner_iteration);
		Lambda = LinksSDLineSearch(MainVolume, SDVolume);
		ColumnPoolUpdate(Lambda, MDMinPathPredLink);
		UpdateVolume(MainVolume, SDVolume, Lambda);
		system_wide_travel_time = UpdateLinkCost(MainVolume);
	}
	// Every outer candidate is evaluated using a fresh full-route-space gap,
	// never the restricted column gap or the pre-update gap from the last FW step.
	system_least_travel_time = FindMinCostRoutes(MDMinPathPredLink);
	gap = (system_wide_travel_time - system_least_travel_time) /
		fmax(0.1, g_relative_gap_standard
			? system_wide_travel_time : system_least_travel_time) * 100.0;
	return gap;
}

static int RunAutoCalibration(
	double* MainVolume,
	double* SubVolume,
	double* SDVolume,
	int*** MDMinPathPredLink,
	double& system_wide_travel_time)
{
	if (g_auto_calibration == 0)
		return 0;

	taplite::AutoCalibrationConfig config =
		taplite::AutoCalibrationConfig::Load(
			g_auto_calibration_config_file,
			demand_period_starting_hours,
			demand_period_ending_hours);
	config.enabled = true;
	if (!g_departure_profile_file.empty())
		config.departure_profile_file = g_departure_profile_file;
	taplite::AutoCalibrationEngine engine(config);
	std::string error_message;
	if (!engine.LoadDepartureProfiles(&error_message))
	{
		printf("ERROR: auto_calibration: %s\n", error_message.c_str());
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "ERROR: auto_calibration: %s\n",
				error_message.c_str());
		return 4;
	}

	double initial_least = FindMinCostRoutes(MDMinPathPredLink);
	double initial_gap = (system_wide_travel_time - initial_least) /
		fmax(0.1, g_relative_gap_standard
			? system_wide_travel_time : initial_least) * 100.0;
	std::vector<taplite::CalibrationLink> links =
		BuildAutoCalibrationLinks(MainVolume);
	if (!engine.LoadVolumeConstraints(links, &error_message))
	{
		printf("ERROR: auto_calibration: %s\n", error_message.c_str());
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "ERROR: auto_calibration: %s\n",
				error_message.c_str());
		return 4;
	}
	if (!g_col_pool_active || g_column_output < 2)
	{
		printf("ERROR: auto_calibration: column_output=2 produced no in-memory routes.\n");
		return 5;
	}
	taplite::CalibrationRouteSummary route_summary =
		BuildAutoCalibrationRouteSummary(engine, links, false);
	engine.Initialize(links, route_summary, initial_gap);
	printf("auto_calibration: initialized %d links and %zu compact mode-link policy cells; "
		"origin profiles=%s, mode fallback profiles=%s, initial objective=%.8f, gap=%.6f%%.\n",
		number_of_links, route_summary.policy_weights.size(),
		engine.uses_origin_specific_profiles() ? "yes" : "no",
		engine.uses_mode_fallback_profiles() ? "yes" : "no",
		engine.accepted_evaluation().objective, initial_gap);

	int accepted_iterations = 0;
	bool stop = false;
	for (int outer = 1; outer <= config.max_outer_iterations && !stop; outer++)
	{
		double step_scale = 1.0;
		bool accepted_this_outer = false;
		for (int retry = 0; retry < config.maximum_retries; retry++)
		{
			links = BuildAutoCalibrationLinks(MainVolume);
			route_summary = BuildAutoCalibrationRouteSummary(
				engine, links, true);
			taplite::CalibrationProposal proposal =
				engine.Propose(links, route_summary, step_scale);
			AutoCalibrationParameterSnapshot snapshot =
				CaptureAutoCalibrationParameters();
			ApplyAutoCalibrationProposal(proposal);
				double candidate_gap = RunAutoCalibrationInnerEquilibrium(
				MainVolume, SubVolume, SDVolume, MDMinPathPredLink,
				system_wide_travel_time, outer,
				config.inner_gap_tolerance_pct);
			std::vector<taplite::CalibrationLink> candidate_links =
				BuildAutoCalibrationLinks(MainVolume);
			taplite::CalibrationRouteSummary candidate_summary =
				BuildAutoCalibrationRouteSummary(
					engine, candidate_links, false);
			taplite::CalibrationEvaluation evaluation = engine.Evaluate(
				candidate_links, candidate_summary, candidate_gap);
			bool accepted = engine.ShouldAccept(evaluation);
			bool converged = accepted && engine.Converged(
				evaluation, proposal, accepted_iterations + 1);
			engine.Record(outer, retry, accepted, step_scale,
				evaluation, proposal);
			printf("auto_calibration outer %d retry %d: objective %.8f, E-fit %.8f, "
				"all-GP speed %.8f, S3 %.8f, envelope %.8f, screenline %.8f, "
				"route distance %.5f, unobserved volume deviation %.5f, "
				"max Q change %.5f, max PLF change %.5f, %s.\n",
				outer, retry, evaluation.objective,
				evaluation.episode_fit_loss,
				evaluation.average_speed_loss,
				evaluation.s3_loss,
				evaluation.volume_envelope_loss,
				evaluation.count_loss,
				evaluation.route_policy_distance,
				evaluation.unobserved_volume_deviation,
				proposal.maximum_parameter_change,
				proposal.maximum_plf_change,
				accepted ? "ACCEPTED" : "REJECTED");
			if (accepted)
			{
				engine.Accept(evaluation);
				accepted_iterations++;
				accepted_this_outer = true;
				stop = converged;
				break;
			}

			// Restore the accepted parameter state, then restore equilibrium in
			// memory before retrying with a smaller trust-region step.
			RestoreAutoCalibrationParameters(snapshot);
			RunAutoCalibrationInnerEquilibrium(
				MainVolume, SubVolume, SDVolume, MDMinPathPredLink,
				system_wide_travel_time, outer,
				config.inner_gap_tolerance_pct);
			step_scale *= config.rejected_step_reduction;
		}
		if (!accepted_this_outer)
		{
			printf("auto_calibration: no acceptable trust-region step at outer %d; stopping at the last accepted state.\n",
				outer);
			break;
		}
	}

	links = BuildAutoCalibrationLinks(MainVolume);
	route_summary = BuildAutoCalibrationRouteSummary(
		engine, links, true, config.arrival_profile_output_file);
	// step_scale=0 preserves final parameters while refreshing departure PLF
	// diagnostics for the accepted final route policy.
	taplite::CalibrationProposal final_proposal =
		engine.Propose(links, route_summary, 0.0);
	if (!engine.WriteOutputs(links, final_proposal, &error_message))
	{
		printf("ERROR: auto_calibration: %s\n", error_message.c_str());
		return 6;
	}
	printf("auto_calibration: final objective %.8f; outputs written once to %s, %s, and %s.\n",
		engine.accepted_evaluation().objective,
		config.history_output_file.c_str(), config.audit_output_file.c_str(),
		config.summary_output_file.c_str());
	return 0;
}

int AssignmentAPI()
{
	auto api_setup_t0 = std::chrono::high_resolution_clock::now();  // SETUP timer (parse + reads + allocations)

	fopen_s(&summary_log_file, "summary_log_file.txt", "w");

	double* MainVolume, * SubVolume, * SDVolume, Lambda;
	int*** MDMinPathPredLink;

	read_settings_file();
	if (!ValidateProcessorCount())
	{
		if (summary_log_file != NULL)
			fclose(summary_log_file);
		summary_log_file = NULL;
		return 2;
	}
	ConfigureOpenMPRuntime();
	if (!InitializeProcessorOriginZones())
	{
		if (summary_log_file != NULL)
			fclose(summary_log_file);
		summary_log_file = NULL;
		return 3;
	}
	read_mode_type_file(); 


	fopen_s(&logfile, "TAP_log.csv", "w");  // Open the log file for writing.
	no_nodes = get_number_of_nodes_from_node_file(no_zones, FirstThruNode);


	if (first_through_node_id_input >= 1)
	{ 
		int first_through_seq = NodeSeqFromExternal(first_through_node_id_input);
		if (first_through_seq != 0)
		{
			FirstThruNode = first_through_seq;
		}
		else
		{
			printf("cannot find node id matching first_through_node_id_input %d", first_through_node_id_input);
		}
	}else if (first_through_node_id_input == 0)
		FirstThruNode = -1; // 


	number_of_links = get_number_of_links_from_link_file();

	int first_through_external_id =
		(FirstThruNode >= 1 && FirstThruNode <= no_nodes)
		? g_node_vector[FirstThruNode].node_id : first_through_node_id_input;
	printf("# of nodes= %d, # of zones = %d, First Through Node ID = %d, number of links = %d\n",
		no_nodes, no_zones, first_through_external_id, number_of_links);

	fprintf(summary_log_file, "no_nodes= %d, no_zones = %d, FirstThruNode Node ID = %d, number_of_links = %d\n", no_nodes, no_zones,
		first_through_external_id, number_of_links);

	fopen_s(&link_performance_file, "link_performance.csv", "w");
	if (link_performance_file == NULL)
	{
		printf("Error opening file!\n");
		return 1;
	}
	fclose(link_performance_file);

	fopen_s(&link_performance_file, "link_performance.csv", "a+");


	double system_wide_travel_time = 0;
	double system_least_travel_time = 0;

	//fprintf(logfile,
	//    "iteration_no,link_id,internal_from_node_id,internal_to_node_id,volume,capacity,voc,"
	//    "fftt,travel_time,delay\n");



	Init(number_of_modes, no_zones);

	InitSPScratch();   // pre-allocate reusable per-thread shortest-path scratch buffers

	if(g_accessibility_only_mode ==0)
		 InitializeLinkIndices(number_of_modes, no_zones, TotalAssignIterations);

	// Origins are the dense zone seq 1..no_zones (zones-first renumbering makes
	// zone seq == centroid node seq, so downstream Minpath uses it directly).
	int assignable_origin_zones = no_zones;
	for (int z = 1; z <= no_zones; z++)
		Processor_origin_zones[z % g_number_of_processors].push_back(z);
	WarnIfProcessorsExceedOrigins(assignable_origin_zones);



	if (g_accessibility_only_mode == 1)
	{
		number_of_internal_zones = 0;
		// Store path information for each destination
		for (int j = 1; j <= no_nodes; j++)
		{
			if (g_node_vector[j].internal_zone_no >= 0)
				number_of_internal_zones++;
		}

		if (number_of_internal_zones < 10000)
			ComputeAccessibilityAndODCosts_v1("od_performance.csv");
		else
		{
			rewriteFile("od_performance.csv", "");
			ComputeAccessibilityAndODCosts_v2("od_performance.csv");
		}

		// The in-process Python/shared-library path returns to a caller that
		// immediately reads the generated CSVs. Close the Windows file handles
		// and release link storage just as the normal assignment path does.
		CloseLinks();
		fclose(link_performance_file);
		fclose(logfile);
		fclose(summary_log_file);
		return 0;
		}

	BuildSPModeGroups();    // exact SP-equivalence mode grouping (allowed_use + AdditionalCost)
	BuildOriginActivity();  // per-(mode, origin) zero-demand skip tables

	int iteration_no = 0;
	int gap_below_count = 0;   // consecutive iterations with gap < convergence_gap_pct
	MainVolume = (double*)Alloc_1D(number_of_links, sizeof(double));
	SDVolume = (double*)Alloc_1D(
		number_of_links, sizeof(double)); /* Compute search direction and sub-volume in the same place. */
	SubVolume = (double*)Alloc_1D(
		number_of_links, sizeof(double)); /* Compute search direction and sub-volume in the same place. */
	MDMinPathPredLink = (int***)Alloc_3D(number_of_modes, no_zones, no_nodes, sizeof(int));

	// Link-state predecessor tree, only needed (and only allocated) when turn
	// restrictions are present; sized by links rather than nodes.
	if (g_has_movement_restrictions)
		g_MinPathPredLink_link = (int***)Alloc_3D(number_of_modes, no_zones, number_of_links, sizeof(int));

	// Record the start time
	auto start = std::chrono::high_resolution_clock::now();

	// SETUP vs COMPUTE accounting (print-only): everything from AssignmentAPI
	// entry to here is input parse + demand read + allocations. The DTAC column
	// file load (warm_start_columns) is also SETUP and prints its own timer
	// inside ApplyWarmStartColumns.
	{
		std::chrono::duration<double> setup_dur = start - api_setup_t0;
		printf("phase timing: SETUP (network+demand parse, allocations) = %.2f s\n", setup_dur.count());
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "phase timing: SETUP (network+demand parse, allocations) = %.2f s\n", setup_dur.count());
	}

	for (int k = 1; k <= number_of_links; k++)
	{
		MainVolume[k] = Link[k].Base_demand_volume;  // assign the base volume  to main volume
	}

	system_wide_travel_time = UpdateLinkCost(MainVolume);  // set up the cost first using FFTT

	ApplyWarmStartTimes();  // L1 warm start (no-op when settings 'warm_start_times' is empty)

	// P2 increment 2: activate the column pool when column_output>=2 (DTAC v2
	// theta-share level) so every FW iteration's AoN paths are folded in with
	// the theta cascade. Level 1 (DTAC v1, last-iteration paths) needs no pool.
	if (g_column_output >= 2)
	{
		InitColumnPool();
		if (g_assignment_method > 0)
			printf("WARNING: column_output with assignment_method=%d (CFW/BFW): theta shares use the plain-FW lambda cascade (same approximation route_assignment.csv uses); exact only for plain FW.\n",
				g_assignment_method);
		if (g_base_demand_mode == 1 && baselinkvolume_loaded_flag == 1)
			printf("WARNING: column_output with base_demand_mode base volumes: theta shares cover the OD-table part of the flow only; base link volumes are not path-decomposed.\n");
	}

	// P2 L3 warm start from COLUMNS: takes precedence over warm_start_flows.
	bool columns_restored = false;
	if (!g_warm_start_columns_file.empty())
	{
		if (!g_warm_start_flows_file.empty())
			printf("WARNING: both warm_start_columns and warm_start_flows are set; columns take precedence (warm_start_flows is used only if the column load fails).\n");
		if (g_base_demand_mode == 1 && baselinkvolume_loaded_flag == 1)
			printf("WARNING: warm_start_columns with base_demand_mode base volumes: columns scatter the FULL OD table on top of the base link volumes -- the diff-assignment contract is NOT applied.\n");
		columns_restored = ApplyWarmStartColumns(MainVolume, MDMinPathPredLink);
	}

	// P1 L2 warm start: restore flows from a DTLR snapshot (no-op/false when
	// 'warm_start_flows' is empty; false with L1-style times seeding when the
	// snapshot's demand fingerprint does not match the current OD table).
	bool flows_restored = columns_restored ? false : ApplyWarmStartFlows(MainVolume);
	bool sp_trees_computed = false;  // MDMinPathPredLink validity (for column_output)

	if (g_link_performance_output == 1)
		fprintf(link_performance_file,
			"iteration_no,link_id,from_node_id,to_node_id,volume,doc,travel_time,speed_mph,VMT\n");
	else
	fprintf(link_performance_file,
		"iteration_no,link_id,from_node_id,to_node_id,volume,ref_volume,base_demand_volume,obs_volume,background_volume,"
		"link_capacity,lane_capacity,D,doc,vdf_fftt,travel_time,vdf_alpha,vdf_beta,vdf_plf,speed_mph,speed_kmph,VMT,VHT,PMT,PHT,VHT_QVDF,PHT_QVDF,geometry,");

	fprintf(logfile, "iteration_no,link_id,from_node_id,to_node_id,volume,ref_volume,obs_volume,capacity,doc,fftt,travel_time,delay,");

	for (int m = 1; m <= number_of_modes; m++)
		fprintf(logfile, "mod_vol_%s,", g_mode_type_vector[m].mode_type.c_str());

	fprintf(logfile, "sub_main_vol,");

	for (int m = 1; m <= number_of_modes; m++)
		fprintf(logfile, "sub_mod_vol_%s,", g_mode_type_vector[m].mode_type.c_str());

	fprintf(logfile, "\n");


	if (g_link_performance_output != 1) {
	for (int m = 1; m <= number_of_modes; m++)
		fprintf(link_performance_file, "mod_vol_%s,", g_mode_type_vector[m].mode_type.c_str());

	for (int m = 1; m <= number_of_modes; m++)
		fprintf(link_performance_file, "obs_vol_%s,", g_mode_type_vector[m].mode_type.c_str());

	fprintf(link_performance_file, "qvdf_profile_status,P,t0,t2,t3,vt2_mph,vt2_kmph,mu,Q_gamma,free_speed_mph,cutoff_speed_mph,free_speed_kmph,cutoff_speed_kmph,congestion_ref_speed_mph,avg_queue_speed_mph,avg_QVDF_period_speed_mph,congestion_ref_speed_kmph,avg_queue_speed_kmph,avg_QVDF_period_speed_kmph,avg_QVDF_period_travel_time,Severe_Congestion_P,");

	for (int t = demand_period_starting_hours * 60; t < demand_period_ending_hours * 60; t += 5)
	{
		int hour = t / 60;
		int minute = t - hour * 60;

		fprintf(link_performance_file, "spd_mph_%02d:%02d,", hour, minute);
	}

	fprintf(link_performance_file, "\n");
	}

	double gap;
	if (!flows_restored && !columns_restored)
	{
	system_least_travel_time = FindMinCostRoutes(MDMinPathPredLink);
	sp_trees_computed = true;
	if(system_wide_travel_time>0)
	{
		gap = (system_wide_travel_time - system_least_travel_time) /
		fmax(0.1, g_relative_gap_standard ? system_wide_travel_time : system_least_travel_time) * 100;

	if (gap < 0)
	{
		int ii = 0;
	}

	{
		std::chrono::duration<double> el = std::chrono::high_resolution_clock::now() - start;
		printf("iter No = %d, sys. TT =  %lf, least TT =  %lf, gap = %f %%, elapsed = %.2f s\n", iteration_no,
			system_wide_travel_time, system_least_travel_time, gap, el.count());
		fflush(stdout);   // per-iteration progress must reach pipes/notebooks live
		fprintf(summary_log_file, "iter No = %d, sys. TT =  %lf, least TT =  %lf, gap = %f %%, elapsed = %.2f s\n", iteration_no,
			system_wide_travel_time, system_least_travel_time, gap, el.count());
	}
	}



	All_or_Nothing_Assign(iteration_no, MDDiffODflow, MDMinPathPredLink, MainVolume);  // here we use MDDiffODflow as our OD search direction of D^c - D^b,

	ColumnPoolUpdate(1.0, MDMinPathPredLink);  // iteration-0 paths enter the pool with share 1 (no-op unless column_output=1)

	}
	else if (flows_restored)
	{
		// P1 L2: MainVolume + per-mode volumes restored from the DTLR snapshot;
		// the iteration-0 all-or-nothing is skipped and FW continues from here.
		printf("warm_start_flows: iteration-0 all-or-nothing skipped; Frank-Wolfe continues from the restored flows.\n");
		fprintf(summary_log_file, "warm_start_flows: iteration-0 AoN skipped; FW continues from restored flows.\n");
		if (g_column_output >= 2)
			printf("WARNING: warm_start_flows + column_output: the restored flows carry no path information, so column theta shares cover post-warm-start iterations only (renormalized at write, warned there too).\n");
	}


	system_wide_travel_time = UpdateLinkCost(MainVolume);

	// P2 L3: fixed-policy adjustment sweeps + the MANDATORY fresh-shortest-path
	// TRUE-gap report -- never let freeze/replay output masquerade as
	// equilibrium. (With number_of_iterations=0 this is the whole run.)
	if (columns_restored)
	{
		if (g_column_adjust_sweeps > 0)
			system_wide_travel_time = RunColumnAdjustSweeps(MainVolume);
		// Both gaps, clearly labeled, ALWAYS (also in freeze/replay with 0 sweeps):
		// restricted = equilibration within the STORED path sets only; true = fresh
		// shortest paths over the full route space. Never let the restricted number
		// masquerade as equilibrium.
		double restricted_gap = ColRestrictedGap();
		system_least_travel_time = FindMinCostRoutes(MDMinPathPredLink);
		sp_trees_computed = true;
		double true_gap = (system_wide_travel_time - system_least_travel_time) /
			fmax(0.1, g_relative_gap_standard ? system_wide_travel_time : system_least_travel_time) * 100;
		printf("warm_start_columns: RESTRICTED gap (STORED paths only, NOT an equilibrium claim) = %f %%; TRUE relative gap (fresh shortest paths, full route space) = %f %% (sys. TT = %.1f, least TT = %.1f).\n",
			restricted_gap, true_gap, system_wide_travel_time, system_least_travel_time);
		fprintf(summary_log_file, "warm_start_columns: RESTRICTED gap (stored paths only) = %f %%; TRUE relative gap (fresh SP) = %f %% (sys TT %.1f, least TT %.1f).\n",
			restricted_gap, true_gap, system_wide_travel_time, system_least_travel_time);
		if (g_ODME_mode == 1)
			printf("WARNING: odme_mode=1 with warm_start_columns: ODME walks the route_output store (linkIndices), which covers post-warm-start iterations only in this run; ODME proceeds on that basis.\n");
		if (shortest_path_log_flag)
			printf("WARNING: warm_start_columns + route_output: route_assignment.csv prob shares cover post-warm-start iterations only (the loaded columns live in the DTAC pool, not the route store).\n");
	}
	else if (g_column_adjust_sweeps > 0)
	{
		printf("WARNING: column_adjust_sweeps=%d but no columns were loaded (warm_start_columns empty or failed); sweeps skipped.\n",
			g_column_adjust_sweeps);
		if (summary_log_file != NULL)
			fprintf(summary_log_file, "WARNING: column_adjust_sweeps=%d ignored (no columns loaded).\n", g_column_adjust_sweeps);
	}

	if (g_tap_log_file == 1)
	{
		for (int k = 1; k <= number_of_links; k++)
		{
			fprintf(logfile, "%d,%d,%d,%d,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,",
				iteration_no, k, Link[k].external_from_node_id, Link[k].external_to_node_id,
				MainVolume[k], Link[k].Ref_volume, Link[k].Obs_volume[1], Link[k].background_volume, Link[k].Link_Capacity,
				MainVolume[k] / fmax(0.01, Link[k].Link_Capacity), Link[k].FreeTravelTime,
				Link[k].Travel_time, Link[k].Travel_time - Link[k].FreeTravelTime);

			for (int m = 1; m <= number_of_modes; m++)
				fprintf(logfile, "%2lf,", Link[k].mode_MainVolume[m]);

			for (int m = 1; m <= number_of_modes; m++)
				fprintf(logfile, "%2lf,", Link[k].Obs_volume[m]);

			fprintf(logfile, "%2lf,", MainVolume[k]);

			for (int m = 1; m <= number_of_modes; m++)
				fprintf(logfile, "%2lf,", Link[k].mode_SubVolume[m]);

			fprintf(logfile, "\n");

		}
	}
	std::chrono::duration<double, std::milli> duration_FindMinCostRoutes;
	std::chrono::duration<double, std::milli> duration_All_or_Nothing_Assign;
	std::chrono::duration<double, std::milli> duration_LinksSDLineSearch;
	auto start0 = std::chrono::high_resolution_clock::now();  // Start timing


	std::vector<double> m_lambda;

	Lambda = 1;
	for (iteration_no = 1; iteration_no < TotalAssignIterations; iteration_no++)
	{
		system_least_travel_time = FindMinCostRoutes(MDMinPathPredLink);  // the one right before the assignment iteration
		sp_trees_computed = true;

		gap = (system_wide_travel_time - system_least_travel_time) /
			fmax(0.1, g_relative_gap_standard ? system_wide_travel_time : system_least_travel_time) * 100;


	

		if (gap < 0)
		{
			int ii = 0;
		}

		{
			std::chrono::duration<double> el = std::chrono::high_resolution_clock::now() - start;
			printf("iter No = %d, Lambda = %f, g_System_VMT = %.1f, sys. TT =  %.1f, least TT =  %.1f, gap = %f %%, elapsed = %.2f s\n",
				iteration_no, Lambda, g_System_VMT, system_wide_travel_time, system_least_travel_time, gap, el.count());
			fflush(stdout);   // per-iteration progress must reach pipes/notebooks live
			fprintf(summary_log_file, "iter No = %d, Lambda = %f, g_System_VMT = %f, sys. TT =  %lf, least TT =  %lf, gap = %f %%, elapsed = %.2f s \n",
				iteration_no, Lambda, g_System_VMT, system_wide_travel_time, system_least_travel_time, gap, el.count());
		}

		// Relative-gap stopping criterion: stop once the gap stays below the target
		// for g_convergence_consecutive successive iterations (ARC: gap<1e-4 x3).
		if (g_convergence_gap_pct > 0.0 && gap >= 0.0 && gap < g_convergence_gap_pct)
		{
			gap_below_count++;
			if (gap_below_count >= g_convergence_consecutive)
			{
				printf("Converged: gap %.4f%% < %.4f%% for %d consecutive iters at iter %d; stopping.\n",
					gap, g_convergence_gap_pct, gap_below_count, iteration_no);
				fprintf(summary_log_file, "Converged: gap %.4f%% < %.4f%% x%d at iter %d\n",
					gap, g_convergence_gap_pct, gap_below_count, iteration_no);
				break;
			}
		}
		else
		{
			gap_below_count = 0;
		}

		g_System_VMT = 0;

		for (int k = 1; k <= number_of_links; k++)
		{
			if (Link[k].link_type >= 1)  // only include physical links
			{
				g_System_VMT += MainVolume[k] * Link[k].length;
			}
		}

		
		All_or_Nothing_Assign(iteration_no, MDODflow, MDMinPathPredLink, SubVolume); // this uees D^c, subvolume is Y, 
		auto end1 = std::chrono::high_resolution_clock::now();    // End timing
		auto start2 = std::chrono::high_resolution_clock::now();  // Start timing

		VolumeDifference(SubVolume, MainVolume, SDVolume); /* Which yields the search direction. SDVolume = y-X */

		// (Bi)conjugate Frank-Wolfe: rewrite SDVolume with a conjugate direction when
		// assignment_method = 1 (CFW) or 2 (BFW). No-op for the default plain FW.
		ApplyConjugateFW(MainVolume, SubVolume, SDVolume, iteration_no);

		Lambda = LinksSDLineSearch(MainVolume, SDVolume);

		m_lambda.push_back(Lambda);

		// P2 increment 2: fold this iteration's AoN paths (trees of the
		// FindMinCostRoutes above) into the column pool with the FW step weight
		// -- the exact theta cascade. No-op unless column_output=1.
		ColumnPoolUpdate(Lambda, MDMinPathPredLink);
		// MSA options
	 //   Lambda = 1.0 / (iteration_no + 1);



		UpdateVolume(MainVolume, SDVolume, Lambda);  // x(k+1) = x(k) +lambda * (y-x) MainVolume is MainVolume
		system_wide_travel_time = UpdateLinkCost(MainVolume);

		//system_least_travel_time = FindMinCostRoutes(MinPathPredLink);  // the one right after the updated link cost 


		if (g_tap_log_file == 1)
		{
			for (int k = 1; k <= number_of_links; k++)
			{
				fprintf(logfile, "%d,%d,%d,%d,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,",
					iteration_no, k, Link[k].external_from_node_id, Link[k].external_to_node_id,
					MainVolume[k], Link[k].Ref_volume, Link[k].Lane_Capacity, Link[k].Link_Capacity,
					MainVolume[k] / fmax(0.01, Link[k].Link_Capacity), Link[k].FreeTravelTime,
					Link[k].Travel_time, Link[k].Travel_time - Link[k].FreeTravelTime);

				for (int m = 1; m <= number_of_modes; m++)
					fprintf(logfile, "%2lf,", Link[k].mode_MainVolume[m]);

				fprintf(logfile, "%2lf,", SDVolume[k]);

				for (int m = 1; m <= number_of_modes; m++)
					fprintf(logfile, "%2lf,", Link[k].mode_SDVolume[m]);

				fprintf(logfile, "\n");

			}
		}
		auto end2 = std::chrono::high_resolution_clock::now();    // End timing

	}

	// Mode 2 outer loop: all inner equilibria reuse the initialized network,
	// demand, route-column pool, and link arrays. No intermediate link/column
	// output is written and no candidate state is read back from disk.
	int auto_calibration_status = RunAutoCalibration(
		MainVolume, SubVolume, SDVolume, MDMinPathPredLink,
		system_wide_travel_time);

	// Record the end time
	auto end = std::chrono::high_resolution_clock::now();

	// Calculate the duration in milliseconds

	// Calculate the duration in seconds
	auto duration = end - start;

	// Convert to hours, minutes, seconds
	auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
	auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration % std::chrono::hours(1));
	auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration % std::chrono::minutes(1));
	auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(duration % std::chrono::seconds(1));

	printf("CPU running time: %lld hours %lld minutes %lld seconds %lld ms\n",
		hours.count(), minutes.count(), seconds.count(), milliseconds.count());

	// Log the result with hours, minutes, seconds, and milliseconds
	fprintf(summary_log_file, "CPU running time: %lld hours %lld minutes %lld seconds %lld ms\n",
		hours.count(), minutes.count(), seconds.count(), milliseconds.count());

	//for (int k = 1; k <= number_of_links; k++)
	//{
	//    fprintf(logfile, "%d,%d,%d,%d,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf\n", iteration_no, k,
	//        Link[k].external_from_node_id, Link[k].external_to_node_id, MainVolume[k],
	//        Link[k].Capacity, MainVolume[k] / fmax(0.01, Link[k].Capacity),
	//        Link[k].FreeTravelTime, Link[k].Travel_time,
	//        Link[k].Travel_time - Link[k].FreeTravelTime);
	//}



	// Compute the flow proportions (theta).
	std::vector<double> m_theta;
	if (m_lambda.size() == 0)
		m_lambda.push_back(1.0);

	m_theta = computeTheta(m_lambda);

	if (g_ODME_mode == 1)
		performODME(m_theta, MainVolume, Link);

	if (g_accessibility_output)
	{
		OutputODPerformance("od_performance.csv");

		GenerateAggregatedPerformanceAndAccessibility();
	}

	if (shortest_path_log_flag)
	{
		OutputRouteDetails("route_assignment.csv", m_theta);
		if (vehicle_log_flag)
			OutputVehicleDetails("vehicle.csv", m_theta);
	}

	// Leveled DTAC column store:
	//   column_output=1 -> DTAC v1, the last-iteration path per OD (light,
	//   route_assignment parity; back-traced from the final SP trees).
	//   column_output>=2 -> DTAC v2, the accumulated theta-share path set
	//   (plus any warm-started columns); falls back to the v1 writer only if
	//   the pool was somehow never activated (defensive).
	if (g_column_file_output != 0 &&
		(g_column_output == 1 || (g_column_output >= 2 && !g_col_pool_active)))
	{
		if (!sp_trees_computed)
		{
			printf("column_output: shortest-path trees never built this run (warm start + <=1 iteration); computing them now.\n");
			FindMinCostRoutes(MDMinPathPredLink);
		}
		WriteColumnsDTAC("route_columns.bin", MDMinPathPredLink);
	}
	else if (g_column_file_output != 0 && g_column_output >= 2)
		WriteColumnsDTACv2("route_columns.bin");
	else if (g_column_output > 0)
		printf("column_output: retained level %d in memory; final route_columns.bin skipped because column_file_output=0.\n",
			g_column_output);
	if (g_auto_calibration_one_shot)
		printf("column_output: one-shot worker defers column-pool reclamation to process exit.\n");
	else
		FreeColumnPool();

	// output link_performance.csv

	// Binary compact output (link_output=3): link_performance.bin. Header "DTLP",
	// int32 version, int32 n_links; then per link: int32 link_id,from,to; float
	// volume,doc,travel_time,speed,vmt. Same fields as compact CSV, fastest to write.
	FILE* link_bin_file = NULL;
	if (g_link_performance_output == 3)
	{
		fopen_s(&link_bin_file, "link_performance.bin", "wb");
		if (link_bin_file)
		{
			int hdr[3] = { 0x504C5444 /*'DTLP'*/, 1, number_of_links };
			fwrite(hdr, sizeof(int), 3, link_bin_file);
		}
	}

	for (int k = 1; k <= number_of_links; k++)
	{
		// Compact output (CSV=1 or binary=3): skip the per-link Link_QueueVDF queue
		// simulation + speed profile (dominant cost on large nets); essentials only.
		if (g_link_performance_output == 1 || g_link_performance_output == 3)
		{
			double Hh = fmax(demand_period_ending_hours - demand_period_starting_hours, 1e-6);
			double lanes_eff = Link[k].Lane_Capacity > 0 ? Link[k].Link_Capacity / Link[k].Lane_Capacity : 1.0;
			double plf_c = fmax(Link[k].VDF_plf, 1e-6);
			double doc_c = Link[k].Lane_Capacity > 0
				? (MainVolume[k] / fmax(lanes_eff, 1.0) / Hh / plf_c) / Link[k].Lane_Capacity : 0.0;
			double spd_c = Link[k].length / fmax(Link[k].Travel_time / 60.0, 0.001);
			double vmt_c = MainVolume[k] * Link[k].length;
			if (g_link_performance_output == 3 && link_bin_file)
			{
				int ids[3] = { Link[k].external_link_id, Link[k].external_from_node_id, Link[k].external_to_node_id };
				float vals[5] = { (float)MainVolume[k], (float)doc_c, (float)Link[k].Travel_time, (float)spd_c, (float)vmt_c };
				fwrite(ids, sizeof(int), 3, link_bin_file);
				fwrite(vals, sizeof(float), 5, link_bin_file);
			}
			else
				fprintf(link_performance_file, "%d,%d,%d,%d,%.4lf,%.4lf,%.4lf,%.3lf,%.4lf\n",
					iteration_no, Link[k].external_link_id, Link[k].external_from_node_id,
					Link[k].external_to_node_id, MainVolume[k], doc_c, Link[k].Travel_time, spd_c, vmt_c);
			continue;
		}

		double P = 0;
		double vt2 = Link[k].Cutoff_Speed;
		double mu = Link[k].Lane_Capacity;
		double Severe_Congestion_P;
		double model_speed[300];
		double t0 = 0;
		double t2 = 0;
		double t3 = 0;
		double Q_gamma = 0;
		double congestion_ref_speed = 0;
		double avg_queue_speed = 0;
		double avg_QVDF_period_speed = 0;
		double IncomingDemand = 0;
		double DOC = 0;
		// Select the reporting profile independently from the assignment VDF. A
		// missing qvdf_profile_mode (-1) preserves the historical behavior exactly:
		// freeway link_type encodings or an observed t2 request a QVDF profile.
		// Explicit mode 0 disables it, mode 1 selects it independent of link_type,
		// and mode 2 gates it on an observed t2. Positive assigned volume and the
		// volume threshold remain hard computational guards for every mode that is
		// otherwise eligible.
		bool is_freeway_link_type =
			(Link[k].link_type == 1) ||
			(Link[k].link_type >= 100 && Link[k].link_type % 100 == 1);
		bool has_observed_qvdf_t2 = std::isfinite(Link[k].QVDF_t2);
		bool has_observed_start_speed =
			std::isfinite(Link[k].QVDF_start_speed_observed);
		bool has_observed_end_speed =
			std::isfinite(Link[k].QVDF_end_speed_observed);
		bool qvdf_eligible = false;
		const char* qvdf_profile_status = "flat_legacy_not_selected";
		switch (Link[k].QVDF_profile_mode)
		{
		case 0:
			qvdf_profile_status = "flat_disabled";
			break;
		case 1:
			qvdf_eligible = true;
			qvdf_profile_status = "generated_model";
			break;
		case 2:
			qvdf_eligible = has_observed_qvdf_t2;
			qvdf_profile_status = qvdf_eligible
				? "generated_observed" : "flat_missing_observation";
			break;
		default:
			qvdf_eligible = is_freeway_link_type || has_observed_qvdf_t2;
			if (qvdf_eligible)
				qvdf_profile_status = is_freeway_link_type
					? "generated_legacy_link_type"
					: "generated_legacy_observed_t2";
			break;
		}
		bool has_positive_qvdf_volume =
			std::isfinite(MainVolume[k]) && MainVolume[k] > 0.0;
		bool do_qvdf = qvdf_eligible && has_positive_qvdf_volume &&
			(MainVolume[k] >= g_qvdf_volume_threshold);
		if (!has_positive_qvdf_volume)
			qvdf_profile_status = "flat_zero_volume";
		else if (qvdf_eligible && !do_qvdf)
			qvdf_profile_status = "flat_below_volume_threshold";
		// Boundary observations are also useful when QVDF generation is skipped.
		// Except for explicit mode 0, use any valid side as an observed-only
		// fallback and retain the assigned period-average speed on a missing side.
		// Preserve the skip reason in the status so the connector is not mistaken
		// for an analytically generated QVDF profile.
		bool use_smoothed_boundary_fallback =
			has_positive_qvdf_volume && !do_qvdf &&
			Link[k].QVDF_profile_mode != 0 &&
			(has_observed_start_speed || has_observed_end_speed);
		if (use_smoothed_boundary_fallback)
		{
			if (qvdf_eligible)
				qvdf_profile_status = "smoothed_boundary_below_volume_threshold";
			else if (Link[k].QVDF_profile_mode == 2)
				qvdf_profile_status = "smoothed_boundary_missing_observation";
			else
				qvdf_profile_status = "smoothed_boundary_legacy_not_selected";
		}
		if (do_qvdf)
		{
			t2 = std::isfinite(Link[k].QVDF_t2)
				? Link[k].QVDF_t2 : DemandPeriodMidpointHours();
			Link_QueueVDF(k, MainVolume[k], IncomingDemand, DOC, P, t0, t2, t3, vt2, mu, Q_gamma, congestion_ref_speed, avg_queue_speed, avg_QVDF_period_speed, Severe_Congestion_P, model_speed);

			// A valid observed anchor at or below the Hermite margin uses the
			// direct anchor-to-vt2 connector in Link_QueueVDF. Report that
			// observation/model conflict explicitly instead of presenting it as
			// an ordinary generated profile.
			double generated_boundary_speed =
				fmax(congestion_ref_speed, avg_queue_speed);
			double generated_anchor_margin = fmax(
				2.0, 0.10 * fmax(0.0, generated_boundary_speed - vt2));
			bool used_low_anchor_connector = has_observed_qvdf_t2 &&
				((has_observed_start_speed &&
					Link[k].QVDF_start_speed_observed < generated_boundary_speed &&
					Link[k].QVDF_start_speed_observed <=
						vt2 + generated_anchor_margin) ||
				 (has_observed_end_speed &&
					Link[k].QVDF_end_speed_observed < generated_boundary_speed &&
					Link[k].QVDF_end_speed_observed <=
						vt2 + generated_anchor_margin));
			if (used_low_anchor_connector)
				qvdf_profile_status = "generated_low_anchor_connector";
		}
		else
		{
			double Hh = fmax(demand_period_ending_hours - demand_period_starting_hours, 1e-6);
			double lanes_eff = Link[k].Lane_Capacity > 0 ? Link[k].Link_Capacity / Link[k].Lane_Capacity : 1.0;
			double plf_e = fmax(Link[k].VDF_plf, 1e-6);
			IncomingDemand = MainVolume[k] / fmax(lanes_eff, 1.0) / Hh / plf_e;
			DOC = Link[k].Lane_Capacity > 0 ? IncomingDemand / Link[k].Lane_Capacity : 0.0;
			Severe_Congestion_P = 0.0;
			// Report a usable, floor-free period-average fallback instead of an
			// apparently failed QVDF calculation. Only invalid/zero travel times use
			// the link free speed as a defensive fallback. These scalar fields remain
			// assignment-consistent even when boundary observations shape the emitted
			// reporting samples below.
			double travel_time_hours = Link[k].Travel_time / 60.0;
			double spd_e =
				std::isfinite(Link[k].length) && Link[k].length >= 0.0 &&
				std::isfinite(travel_time_hours) && travel_time_hours > 0.0
				? Link[k].length / travel_time_hours
				: Link[k].free_speed;
			vt2 = spd_e;
			congestion_ref_speed = spd_e;
			avg_queue_speed = spd_e;
			avg_QVDF_period_speed = spd_e;
			Link[k].QVDF_TT = Link[k].Travel_time;
			for (int i = 0; i < 300; i++) model_speed[i] = spd_e;

			if (use_smoothed_boundary_fallback)
			{
				int profile_start_min =
					static_cast<int>(demand_period_starting_hours * 60);
				int period_end_min =
					static_cast<int>(demand_period_ending_hours * 60);
				int profile_last_min =
					std::max(profile_start_min, period_end_min - 5);
				double start_speed = has_observed_start_speed
					? Link[k].QVDF_start_speed_observed : spd_e;
				double end_speed = has_observed_end_speed
					? Link[k].QVDF_end_speed_observed : spd_e;
				double profile_span_min =
					static_cast<double>(profile_last_min - profile_start_min);
				auto smoothstep01 = [](double value)
					{
						double position = fmin(1.0, fmax(0.0, value));
						return position * position * (3.0 - 2.0 * position);
					};

				for (int t_in_min = profile_start_min;
					t_in_min <= profile_last_min; t_in_min += 5)
				{
					double factor = profile_span_min > 1e-9
						? (t_in_min - profile_start_min) / profile_span_min
						: 0.0;
					double smooth_factor = smoothstep01(factor);
					model_speed[t_in_min / 5] =
						(1.0 - smooth_factor) * start_speed
						+ smooth_factor * end_speed;
				}

				// The end observation wins when an unusually short period emits only
				// one sample, matching the generated-profile anchor precedence.
				model_speed[profile_start_min / 5] = start_speed;
				if (profile_last_min > profile_start_min || has_observed_end_speed)
					model_speed[profile_last_min / 5] = end_speed;

				Severe_Congestion_P = 0.0;
				for (int t_in_min = profile_start_min;
					t_in_min <= profile_last_min; t_in_min += 5)
				{
					if (model_speed[t_in_min / 5] < Link[k].free_speed * 0.5)
						Severe_Congestion_P += 5.0 / 60;
				}
			}
		}


		// MainVolume is the total PCE-weighted link volume, so vehicle-equivalent
		// measures are link totals and must be calculated exactly once. Only the
		// occupancy-weighted person measures belong inside the mode loop.
		double VMT = MainVolume[k] * Link[k].length;
		double VHT = MainVolume[k] * Link[k].Travel_time / 60.0;
		double VHT_QVDF = MainVolume[k] * Link[k].QVDF_TT / 60.0;
		double PMT = 0;
		double PHT = 0;
		double PHT_QVDF = 0;
		for (int m = 1; m <= number_of_modes; m++)
		{
			PMT += Link[k].mode_MainVolume[m] * g_mode_type_vector[m].occ * Link[k].length;
			PHT += Link[k].mode_MainVolume[m] * g_mode_type_vector[m].occ * Link[k].Travel_time / 60.0;

			PHT_QVDF += Link[k].mode_MainVolume[m] * g_mode_type_vector[m].occ * Link[k].QVDF_TT / 60.0;


		}

		fprintf(link_performance_file, "%d,%d,%d,%d,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,",
			iteration_no, Link[k].external_link_id, Link[k].external_from_node_id, Link[k].external_to_node_id,
			MainVolume[k], Link[k].Ref_volume, Link[k].Base_demand_volume, Link[k].Obs_volume[1], Link[k].background_volume, Link[k].Link_Capacity, Link[k].Lane_Capacity, IncomingDemand, DOC, Link[k].FreeTravelTime,
			Link[k].Travel_time, Link[k].VDF_Alpha, Link[k].VDF_Beta, Link[k].VDF_plf, Link[k].length / fmax(Link[k].Travel_time / 60.0, 0.001), Link[k].length / fmax(Link[k].Travel_time / 60.0, 0.001) * 1.609, Link[k].Travel_time - Link[k].FreeTravelTime);

		fprintf(link_performance_file, "%2lf,%2lf,%2lf,%2lf,%2lf,%2lf,", VMT, VHT, PMT, PHT, VHT_QVDF, PHT_QVDF);

		fprintf(link_performance_file, "\"%s\",",
			Link[k].geometry.c_str());

		for (int m = 1; m <= number_of_modes; m++)
			fprintf(link_performance_file, "%2lf,", Link[k].mode_MainVolume[m]);
	
		for (int m = 1; m <= number_of_modes; m++)
			fprintf(link_performance_file, "%2lf,", Link[k].Obs_volume[m]);

		fprintf(link_performance_file, "%s,", qvdf_profile_status);
		fprintf(link_performance_file, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,", P, t0, t2, t3, vt2, vt2 * 1.609, mu, Q_gamma,
			Link[k].free_speed, Link[k].Cutoff_Speed, Link[k].free_speed * 1.609, Link[k].Cutoff_Speed * 1.609,
			congestion_ref_speed, avg_queue_speed, avg_QVDF_period_speed,
			congestion_ref_speed * 1.609, avg_queue_speed * 1.609, avg_QVDF_period_speed * 1.609,
			Link[k].QVDF_TT, Severe_Congestion_P);
		for (int t = demand_period_starting_hours * 60; t < demand_period_ending_hours * 60; t += 5)
		{
			int t_interval = t / 5;
			double speed = model_speed[t_interval];
			fprintf(link_performance_file, "%.3f,", speed);
		}



		fprintf(link_performance_file, "\n");
	}

	if (link_bin_file)
		fclose(link_bin_file);

	// P1: DTLR flow snapshot -- written at the end of EVERY run when
	// link_output=3 or flow_snapshot=1; input for warm_start_flows.
	if (g_flow_snapshot == 1 || g_link_performance_output == 3)
		WriteFlowSnapshotDTLR("link_flows.bin", MainVolume);

	if (!g_auto_calibration_one_shot)
	{
		free(MainVolume);
		free(SubVolume);
		free(SDVolume);

		Free_3D((void***)MDMinPathPredLink, number_of_modes, no_zones, no_nodes);

		if (g_MinPathPredLink_link != NULL)
		{
			Free_3D((void***)g_MinPathPredLink_link, number_of_modes, no_zones, number_of_links);
			g_MinPathPredLink_link = NULL;
		}

		Close();
	}

	fclose(link_performance_file);
	fclose(logfile);  // Close the log file when you're done with it.
	fclose(summary_log_file);
	return auto_calibration_status;
}

static void Init(int int_number_of_modes, int input_no_zones)
{
	// tuiInit(tuiFileName);
	InitLinks();

	if (g_accessibility_only_mode == 0)
	{
		baseODDemand_loaded_flag = Read_ODflow(&TotalODflow, &int_number_of_modes, &input_no_zones);
	}
	if (baseODDemand_loaded_flag == 0)
	{
		// reset 
		for (int k = 1; k <= number_of_links; k++)
		{

			Link[k].Base_demand_volume = 0;

			for (int m = 1; m <= int_number_of_modes; m++)
			{
				Link[k].mode_Base_demand_volume[m] = 0;
			}
		}
	}
}

static void Close()
{
	StatusMessage("General", "Closing all modules");
	// tuiClose(tuiFileName);
	CloseLinks();
	CloseODflow();

}


static void CloseODflow(void)
{

	free(TotalOFlow);

	Free_3D((void***)seed_MDODflow, number_of_modes, no_zones, no_zones);

	Free_3D((void***)MDODflow, number_of_modes, no_zones, no_zones);
	Free_3D((void***)targetMDODflow, number_of_modes, no_zones, no_zones);

	if (g_ODME_mode == 1)
	{
		Free_3D((void***)old_MDODflow, number_of_modes, no_zones, no_zones);
		Free_3D((void***)candidate_MDODflow, number_of_modes, no_zones, no_zones);
		Free_3D((void***)gradient_MDODflow, number_of_modes, no_zones, no_zones);
	}

	Free_3D((void***)MDDiffODflow, number_of_modes, no_zones, no_zones);
	Free_3D((void***)MDRouteCost, number_of_modes, no_zones, no_zones);

	if (g_mode_origin_active != NULL)
	{
		Free_2D((void**)g_mode_origin_active, number_of_modes, no_zones);
		g_mode_origin_active = NULL;
	}
	if (g_tree_origin_active != NULL)
	{
		Free_2D((void**)g_tree_origin_active, number_of_modes, no_zones);
		g_tree_origin_active = NULL;
	}
}






void ReadLinks()
{

	CDTACSVParser parser_link;
	std::vector<CLink> links;

	std::string line;

	if (parser_link.OpenCSVFile("link.csv", true))
	{
		float total_base_link_volume = 0;
		int line_no = 0;

		int k = 1;  // link start from 1

		while (parser_link.ReadRecord())  // if this line contains [] mark, then we will also read
			// field headers.
		{
			// CLink link;
			int lanes = 1;
			float capacity = 0;
			float free_speed = 10;

			parser_link.GetValueByFieldName("lanes", lanes);
			parser_link.GetValueByFieldName("capacity", capacity);
			parser_link.GetValueByFieldName("free_speed", free_speed);

			if (g_metric_system_flag == 1)
				free_speed = free_speed / 1.609;

			parser_link.GetValueByFieldName("vdf_free_speed_mph", free_speed);

			if (lanes <= 0 || capacity < 0.0001 || free_speed < 0.0001) {
				std::cerr << "Error: Invalid link attributes detected. "
					<< "lanes = " << lanes
					<< ", capacity = " << capacity
					<< ", free_speed = " << free_speed
					<< ". Skipping this link." << std::endl;
				continue;
			}


			Link[k].setup(number_of_modes);
			std::string value;



			// Read link_id
			parser_link.GetValueByFieldName("from_node_id", Link[k].external_from_node_id);
			parser_link.GetValueByFieldName("to_node_id", Link[k].external_to_node_id);
			parser_link.GetValueByFieldName("link_id", Link[k].external_link_id);

			Link[k].link_id = k; // rewrite using internal link id; external id preserved in external_link_id
			parser_link.GetValueByFieldName("link_type", Link[k].link_type);



			Link[k].internal_from_node_id = Link[k].external_from_node_id;
			Link[k].internal_to_node_id = Link[k].external_to_node_id;

			int internal_from_node_id = NodeSeqFromExternal(Link[k].external_from_node_id);
			if (internal_from_node_id != 0)
			{
				Link[k].internal_from_node_id = internal_from_node_id;
			}
			else
			{
				printf("Error in from_node_id =%d for link_id = %d\n", Link[k].external_from_node_id, Link[k].link_id);

				continue; 
			}

			if (Link[k].internal_from_node_id == 0)
			{
				printf("Error in Link[k].internal_from_node_id\n"); 
			}
			int internal_to_node_id = NodeSeqFromExternal(Link[k].external_to_node_id);
			if (internal_to_node_id != 0)
			{
				Link[k].internal_to_node_id = internal_to_node_id;
			}
			else
			{
				printf("Error in to_node_id =%d for link_id = %d\n", Link[k].external_to_node_id, Link[k].link_id);
				continue;
			}

			parser_link.GetValueByFieldName("length", Link[k].length);  // meter
			Link[k].length = Link[k].length / 1609.0;  // miles

			double vdf_length_mi = -1; 
			parser_link.GetValueByFieldName("vdf_length_mi", vdf_length_mi);

			if (vdf_length_mi >= 0)  // valid. 
				Link[k].length = vdf_length_mi;

			parser_link.GetValueByFieldName("ref_volume", Link[k].Ref_volume);

			parser_link.GetValueByFieldName("non_uturn_flag", Link[k].non_uturn_flag);
			

			for (int m = 1; m <= number_of_modes; m++)
			{
				char CSV_field_name[50];

				if (number_of_modes == 1)
				{
					sprintf(CSV_field_name, "obs_volume");
				}
				else
				{
					sprintf(CSV_field_name, "obs_volume_%s", g_mode_type_vector[m].mode_type.c_str());
				}
				parser_link.GetValueByFieldName(CSV_field_name, Link[k].Obs_volume[m], false, false);

			}

			if (g_base_demand_mode == 1)
			{

				parser_link.GetValueByFieldName("base_demand_volume", Link[k].Base_demand_volume);
				total_base_link_volume += Link[k].Base_demand_volume;

				if (number_of_modes == 1)  // single mode 
					Link[k].mode_Base_demand_volume[1] = Link[k].Base_demand_volume;

				for (int m = 1; m <= number_of_modes; m++)
				{
					std::string field_name = "base_vol_" + std::string(g_mode_type_vector[m].mode_type);
					parser_link.GetValueByFieldName(field_name.c_str(), Link[k].mode_Base_demand_volume[m]);

				}
			}



			Link[k].lanes = lanes;
			Link[k].Lane_Capacity = capacity;
			Link[k].Link_Capacity = lanes * capacity;


			parser_link.GetValueByFieldName("allowed_use", Link[k].allowed_uses);

			if (g_tap_log_file == 1)
			{
				fprintf(logfile, "link %d->%d, node_seq_no %d->%d\n", Link[k].external_from_node_id, Link[k].external_to_node_id, Link[k].internal_from_node_id, Link[k].internal_to_node_id);

			}


			for (int m = 1; m <= number_of_modes; m++)
			{
				Link[k].mode_allowed_use[m] = 1;
				Link[k].mode_MainVolume[m] = 0;
				Link[k].mode_SubVolume[m] = 0;
				Link[k].mode_Toll[m] = 0;
				Link[k].mode_AdditionalCost[m] = 0;

			}

			if (Link[k].allowed_uses.size() > 0 && Link[k].allowed_uses != "all")
			{
				for (int m = 1; m <= number_of_modes; m++)
				{
					if (Link[k].allowed_uses.find(g_mode_type_vector[m].mode_type) != std::string::npos)  // otherwise, only an agent type is listed in this "allowed_uses", then this agent type is allowed to travel on this link
					{
						Link[k].mode_allowed_use[m] = 1;  // found 
					}
					else
					{
						Link[k].mode_allowed_use[m] = 0;
					}
				}


			}



			// Read internal_from_node_id

			// Read length

			// Read capacity

			Link[k].FreeTravelTime = Link[k].length / fmax(0.001, free_speed) * 60.0;

			parser_link.GetValueByFieldName("vdf_fftt", Link[k].FreeTravelTime,true);

			parser_link.GetValueByFieldName("vdf_alpha", Link[k].VDF_Alpha, true);
			parser_link.GetValueByFieldName("vdf_beta", Link[k].VDF_Beta, true);
			parser_link.GetValueByFieldName("vdf_A", Link[k].VDF_A, false, false);  // ARC modified-BPR linear term (default 0)
			parser_link.GetValueByFieldName("green_ratio", Link[k].green_ratio, false, false);  // SANDAG vdf_type=6 (default 0.45)
			parser_link.GetValueByFieldName("vdf_plf", Link[k].VDF_plf, true);

			// Conic (Spiess) VDF support: vdf_type=1 selects conic using per-link
			// conic_a / conic_b (per functional type). Optional columns; default
			// stays BPR (vdf_type=0). If conic_b absent, derive Spiess beta below.
			parser_link.GetValueByFieldName("vdf_type", Link[k].VDF_type);
			parser_link.GetValueByFieldName("conic_a", Link[k].Conic_a);
			parser_link.GetValueByFieldName("conic_b", Link[k].Conic_b);
			if (Link[k].VDF_type == 1 && Link[k].Conic_b <= 0.0 && Link[k].Conic_a > 1.0)
				Link[k].Conic_b = (2.0 * Link[k].Conic_a - 1.0) / (2.0 * Link[k].Conic_a - 2.0);


			
			// Link-level single-class toll fallback: the designed `vdf_toll` column
			// is used when a per-class `toll_<mode>` field is absent.
			double vdf_toll_val = 0;
			parser_link.GetValueByFieldName("vdf_toll", vdf_toll_val, false, false);
			for (int m = 1; m <= number_of_modes; m++)
			{
				char CSV_field_name[50];

				if(number_of_modes == 1)
				{
					sprintf(CSV_field_name, "toll");
				}
				else
				{
				sprintf(CSV_field_name, "toll_%s", g_mode_type_vector[m].mode_type.c_str());
				}
				if (!parser_link.GetValueByFieldName(CSV_field_name, Link[k].mode_Toll[m], false, false))
					Link[k].mode_Toll[m] = vdf_toll_val;

				// ARC generalized cost (in minutes): time + (toll + distance*op_cost)/VOT*60.
				// op_cost defaults to 0 -> distance term off; toll term unchanged.
				Link[k].mode_AdditionalCost[m] =
					(Link[k].mode_Toll[m] + Link[k].length * g_mode_type_vector[m].op_cost)
					/ g_mode_type_vector[m].vot * 60.0;
			}



			if (capacity > 0)
				Link[k].BoverC = Link[k].VDF_Alpha / pow(capacity * lanes, Link[k].VDF_Beta);
			else
				Link[k].BoverC = 0;


			parser_link.GetValueByFieldName("vdf_cp", Link[k].Q_cp);
			parser_link.GetValueByFieldName("vdf_cd", Link[k].Q_cd);
			parser_link.GetValueByFieldName("vdf_n", Link[k].Q_n);
			parser_link.GetValueByFieldName("vdf_s", Link[k].Q_s);

			// Auto-calibration handoff. Mode 1 priors default to the converted
			// QVDF values and may be overridden explicitly. No-episode links keep
			// P/vt2 as NaN; a zero P is never manufactured here.
			Link[k].Mode1_plf = Link[k].VDF_plf;
			Link[k].Mode1_Q_cd = Link[k].Q_cd;
			Link[k].Mode1_Q_cp = Link[k].Q_cp;
			parser_link.GetValueByFieldName("mode1_plf", Link[k].Mode1_plf, false, false);
			parser_link.GetValueByFieldName("mode1_qcd", Link[k].Mode1_Q_cd, false, false);
			parser_link.GetValueByFieldName("mode1_qcp", Link[k].Mode1_Q_cp, false, false);
			if (!parser_link.GetValueByFieldName(
				"calibration_observation_class", Link[k].Calibration_observation_class, false))
				parser_link.GetValueByFieldName(
					"tmc_observation_class", Link[k].Calibration_observation_class, false);
			parser_link.GetValueByFieldName("corridor", Link[k].Calibration_corridor, false);
			parser_link.GetValueByFieldName(
				"facility_class", Link[k].Calibration_facility_class, false);
			parser_link.GetValueByFieldName(
				"target_tmc", Link[k].Calibration_target_tmc, false);
			parser_link.GetValueByFieldName(
				"calibration_exclusion_reason",
				Link[k].Calibration_exclusion_reason, false);
			std::string calibration_facility_class =
				Link[k].Calibration_facility_class;
			std::transform(calibration_facility_class.begin(),
				calibration_facility_class.end(),
				calibration_facility_class.begin(),
				[](unsigned char c) { return static_cast<char>(std::tolower(c)); });
			Link[k].Calibration_eligible =
				calibration_facility_class.empty() ||
				calibration_facility_class == "gp";
			if (!Link[k].Calibration_eligible)
				Link[k].Calibration_observation_class = "U";
			parser_link.GetValueByFieldName("observed_p_hr", Link[k].Calibration_observed_P, false, false);
			parser_link.GetValueByFieldName("observed_vt2_mph", Link[k].Calibration_observed_vt2, false, false);
			parser_link.GetValueByFieldName("observed_avg_speed_mph", Link[k].Calibration_observed_average_speed, false, false);
			parser_link.GetValueByFieldName("s3_volume", Link[k].Calibration_S3_volume, false, false);
			parser_link.GetValueByFieldName("cube_vehicle_volume", Link[k].Calibration_cube_vehicle_volume, false, false);
			parser_link.GetValueByFieldName("observation_quality", Link[k].Calibration_observation_quality, false, false);
			if (Link[k].Calibration_eligible &&
				taplite::ParseObservationClass(Link[k].Calibration_observation_class) ==
				taplite::ObservationClass::Unobserved)
			{
				if (std::isfinite(Link[k].Calibration_observed_P) &&
					Link[k].Calibration_observed_P > 0.0 &&
					std::isfinite(Link[k].Calibration_observed_vt2) &&
					Link[k].Calibration_observed_vt2 > 0.0)
					Link[k].Calibration_observation_class = "E";
				else if (std::isfinite(Link[k].Calibration_observed_average_speed) ||
					std::isfinite(Link[k].Calibration_S3_volume))
					Link[k].Calibration_observation_class = "N";
			}
			if (!Link[k].Calibration_eligible)
				Link[k].Calibration_observation_class = "U";

			// Optional QVDF reporting-profile mode: blank/missing = legacy auto,
			// 0 = disabled, 1 = model-generated, 2 = observed-t2 gated. Parse as
			// text so blank cells remain distinguishable from an explicit zero.
			std::string qvdf_profile_mode_text;
			if (parser_link.GetValueByFieldName(
				"qvdf_profile_mode", qvdf_profile_mode_text, false))
			{
				TrimString(qvdf_profile_mode_text);
				if (!qvdf_profile_mode_text.empty())
				{
					char* endp = NULL;
					long parsed_mode = strtol(
						qvdf_profile_mode_text.c_str(), &endp, 10);
					if (endp != qvdf_profile_mode_text.c_str() &&
						*endp == '\0' && parsed_mode >= 0 && parsed_mode <= 2)
						Link[k].QVDF_profile_mode = static_cast<int>(parsed_mode);
					else
					{
						fprintf(stderr,
							"WARNING: link_id=%d has invalid qvdf_profile_mode='%s'; expected blank, 0, 1, or 2. Using legacy auto.\n",
							Link[k].external_link_id, qvdf_profile_mode_text.c_str());
						if (summary_log_file != NULL)
							fprintf(summary_log_file,
								"WARNING: link_id=%d has invalid qvdf_profile_mode='%s'; expected blank, 0, 1, or 2. Using legacy auto.\n",
								Link[k].external_link_id, qvdf_profile_mode_text.c_str());
					}
				}
			}

			// Optional observed episode endpoints. They may extend outside this
			// assignment period because they are used only to recover the
			// before/after-t2 episode proportion. Missing, partial, malformed, or
			// unordered endpoint data silently retains the symmetric fallback.
			auto read_optional_episode_hour =
				[&parser_link](const char* column_name, double& target)
				{
					std::string text;
					if (!parser_link.GetValueByFieldName(column_name, text, false))
						return;
					TrimString(text);
					if (text.empty())
						return;
					char* endp = NULL;
					double value = strtod(text.c_str(), &endp);
					if (endp != text.c_str() && *endp == '\0' &&
						std::isfinite(value) && value >= 0.0 && value <= 24.0)
						target = value;
				};
			read_optional_episode_hour("t0_hour", Link[k].QVDF_t0_observed);
			read_optional_episode_hour("t3_hour", Link[k].QVDF_t3_observed);

			// Optional observed speeds for the first and last emitted five-minute
			// profile samples. Each endpoint falls back independently to the modeled
			// profile when its column is missing, blank, malformed, or non-positive.
			auto read_optional_boundary_speed =
				[&](const char* column_name, double& target)
				{
					std::string text;
					if (!parser_link.GetValueByFieldName(column_name, text, false))
						return;
					TrimString(text);
					if (text.empty())
						return;
					char* endp = NULL;
					double value = strtod(text.c_str(), &endp);
					if (endp != text.c_str() && *endp == '\0' &&
						std::isfinite(value) && value > 0.0)
					{
						target = value;
						return;
					}

					fprintf(stderr,
						"WARNING: link_id=%d has invalid %s='%s'; expected a finite positive speed in mph. Using modeled boundary fallback.\n",
						Link[k].external_link_id, column_name, text.c_str());
					if (summary_log_file != NULL)
						fprintf(summary_log_file,
							"WARNING: link_id=%d has invalid %s='%s'; expected a finite positive speed in mph. Using modeled boundary fallback.\n",
							Link[k].external_link_id, column_name, text.c_str());
				};
			read_optional_boundary_speed(
				"qvdf_start_speed_mph", Link[k].QVDF_start_speed_observed);
			read_optional_boundary_speed(
				"qvdf_end_speed_mph", Link[k].QVDF_end_speed_observed);

			// Optional observed QVDF trough time. The CBI handoff uses t2_hour;
			// accept t2 as a compatibility alias. Empty/missing cells retain the
			// historical assignment-period midpoint. Invalid supplied values are
			// rejected with a link-specific diagnostic rather than clamped.
			std::string qvdf_t2_text;
			const char* qvdf_t2_column = "t2_hour";
			bool qvdf_t2_supplied =
				parser_link.GetValueByFieldName(qvdf_t2_column, qvdf_t2_text, false);
			if (!qvdf_t2_supplied)
			{
				qvdf_t2_column = "t2";
				qvdf_t2_supplied =
					parser_link.GetValueByFieldName(qvdf_t2_column, qvdf_t2_text, false);
			}
			if (qvdf_t2_supplied)
			{
				TrimString(qvdf_t2_text);
				char* endp = NULL;
				double observed_t2 = strtod(qvdf_t2_text.c_str(), &endp);
				bool parsed = endp != qvdf_t2_text.c_str() && *endp == '\0';
				if (parsed && IsHourWithinDemandPeriod(observed_t2))
					Link[k].QVDF_t2 = observed_t2;
				else
				{
					double fallback_t2 = DemandPeriodMidpointHours();
					fprintf(stderr,
						"WARNING: link_id=%d has invalid %s='%s'; expected a finite decimal hour within the configured demand period %.3f-%.3f. Using midpoint fallback %.3f.\n",
						Link[k].external_link_id, qvdf_t2_column, qvdf_t2_text.c_str(),
						demand_period_starting_hours, demand_period_ending_hours, fallback_t2);
					if (summary_log_file != NULL)
						fprintf(summary_log_file,
							"WARNING: link_id=%d has invalid %s='%s'; expected a finite decimal hour within the configured demand period %.3f-%.3f. Using midpoint fallback %.3f.\n",
							Link[k].external_link_id, qvdf_t2_column, qvdf_t2_text.c_str(),
							demand_period_starting_hours, demand_period_ending_hours, fallback_t2);
				}
			}

			Link[k].free_speed = free_speed;
			// Default cutoff (speed at capacity) = 0.75 * free_speed when not given.
			// If link.csv provides an explicit cutoff_speed column (mph, same unit
			// convention as vdf_free_speed_mph), use it. The QVDF queue speed model
			// is sensitive to this value, so an explicit calibrated cutoff_speed
			// (e.g. 49 mph in the I-10 corridor) must not be silently overridden.
			Link[k].Cutoff_Speed = free_speed * 0.75;
			parser_link.GetValueByFieldName("cutoff_speed", Link[k].Cutoff_Speed, false);
			parser_link.GetValueByFieldName("geometry", Link[k].geometry, false);
			k++;
		}

		const int input_link_rows = number_of_links;
		const int valid_links = k - 1;

		// The zones-first node renumbering can change the ordering of from-node
		// ids even when link.csv was sorted by the external ids. The shortest
		// path engine uses contiguous outgoing-link ranges, so group links by
		// dense from-node id with a stable counting sort. This is O(E + N), not
		// O(E log E), and gives every accepted link a dense 1..E internal id.
		std::vector<int> link_count_by_from(no_nodes + 1, 0);
		for (int old_k = 1; old_k <= valid_links; old_k++)
			link_count_by_from[Link[old_k].internal_from_node_id]++;

		std::vector<int> next_link_position(no_nodes + 1, 0);
		int next_position = 1;
		for (int node_seq = 1; node_seq <= no_nodes; node_seq++)
		{
			next_link_position[node_seq] = next_position;
			next_position += link_count_by_from[node_seq];
		}

		link_record* ordered_links = new link_record[valid_links + 1];
		bool reordered = false;
		for (int old_k = 1; old_k <= valid_links; old_k++)
		{
			int new_k = next_link_position[Link[old_k].internal_from_node_id]++;
			reordered = reordered || (new_k != old_k);
			ordered_links[new_k] = std::move(Link[old_k]);
			ordered_links[new_k].link_id = new_k;
		}
		delete[] Link;
		Link = ordered_links;
		number_of_links = valid_links;

		g_map_external_link_id_2_link_seq_no.clear();
		g_map_external_link_id_2_link_seq_no.reserve((size_t)number_of_links * 2 + 1);
		for (int node_seq = 1; node_seq <= no_nodes; node_seq++)
		{
			g_node_vector[node_seq].m_incoming_link_seq_no_vector.clear();
			g_node_vector[node_seq].m_outgoing_link_seq_no_vector.clear();
		}
		for (int internal_link_id = 1; internal_link_id <= number_of_links; internal_link_id++)
		{
			int external_link_id = Link[internal_link_id].external_link_id;
			if (g_map_external_link_id_2_link_seq_no.find(external_link_id) !=
				g_map_external_link_id_2_link_seq_no.end())
			{
				printf("WARNING: duplicate link_id %d in link.csv; external-id lookups keep the first.\n",
					external_link_id);
			}
			else
			{
				g_map_external_link_id_2_link_seq_no[external_link_id] = internal_link_id;
			}
			g_node_vector[Link[internal_link_id].internal_to_node_id].
				m_incoming_link_seq_no_vector.push_back(internal_link_id);
			g_node_vector[Link[internal_link_id].internal_from_node_id].
				m_outgoing_link_seq_no_vector.push_back(internal_link_id);
		}

		if (input_link_rows != valid_links)
			printf("WARNING: accepted %d of %d link.csv rows; engine arrays use the %d valid links only.\n",
				valid_links, input_link_rows, valid_links);
		if (reordered)
			printf("NOTE: links grouped by dense internal from-node id in O(E+N); external link ids are preserved in outputs.\n");

		printf("total_base_link_volume = %f\n", total_base_link_volume);

		if (total_base_link_volume > 0)
			baselinkvolume_loaded_flag = 1;
		parser_link.CloseCSVFile();
	}



}


void Load_Movement_Restrictions(const std::string& filename)
{
	// movement.csv ib_link_id / ob_link_id are EXTERNAL link ids (link.csv link_id).
	// The shortest path keys restrictions by INTERNAL link index, so translate
	// external -> internal once here. (Previously the external ids were used as
	// internal indices, which only worked when they happened to coincide.)
	int n_restrictions = 0;

	CDTACSVParser parser_settings;

	if (parser_settings.OpenCSVFile(filename, true)) // Open CSV file
	{
		while (parser_settings.ReadRecord()) // Read each record (row)
		{
			int mvmt_id = 0, node_id = 0, ib_ext = 0, ob_ext = 0;
			float penalty = 0;

			parser_settings.GetValueByFieldName("mvmt_id", mvmt_id);
			parser_settings.GetValueByFieldName("node_id", node_id);
			parser_settings.GetValueByFieldName("ib_link_id", ib_ext);
			parser_settings.GetValueByFieldName("ob_link_id", ob_ext);
			parser_settings.GetValueByFieldName("penalty", penalty);

			// penalty >= 10 marks a forbidden (hard-banned) movement.
			if (penalty >= 10)
			{
				int ib = LinkSeqFromExternal(ib_ext);
				int ob = LinkSeqFromExternal(ob_ext);
				if (ib != 0 && ob != 0)
				{
					InsertMovementRestriction(ib, ob, true);
					Link[ib].b_withmovement_restrictions = true;
					n_restrictions++;
				}
				else
				{
					printf("Warning: movement.csv references unknown link id (ib=%d, ob=%d); skipped.\n", ib_ext, ob_ext);
				}
			}
		}

		parser_settings.CloseCSVFile(); // Close file after reading
	}
	else
	{
		cout << "Warning: Could not open " << filename << std::endl;
	}

	// non_uturn_flag (link.csv): forbid the U-turn that immediately returns along
	// the reverse link (head->tail back to tail->head). Folded into the same
	// restriction mechanism so it is handled exactly by the link-state search.
	int n_uturn = 0;
	for (int k = 1; k <= number_of_links; k++)
	{
		if (Link[k].non_uturn_flag == 0)
			continue;
		int a = Link[k].internal_from_node_id;
		int b = Link[k].internal_to_node_id;   // arrived at b via link k
		// the U-turn is any outgoing link from b that returns to a
		for (int k2 = FirstLinkFrom[b]; k2 <= LastLinkFrom[b]; k2++)
		{
			if (Link[k2].internal_to_node_id == a)
			{
				InsertMovementRestriction(k, k2, true);
				Link[k].b_withmovement_restrictions = true;
				n_uturn++;
			}
		}
	}

	if (n_restrictions > 0 || n_uturn > 0)
	{
		g_has_movement_restrictions = 1;
		printf("Turn restrictions active: %d from movement.csv, %d U-turn bans from non_uturn_flag.\n", n_restrictions, n_uturn);
		fprintf(summary_log_file, "Turn restrictions active: %d from movement.csv, %d U-turn bans from non_uturn_flag.\n", n_restrictions, n_uturn);
	}
}
static void InitLinkPointers(char* LinksFileName)
{
	(void)LinksFileName;
	FirstLinkFrom = (int*)Alloc_1D(no_nodes, sizeof(int));
	LastLinkFrom = (int*)Alloc_1D(no_nodes, sizeof(int));

	// ReadLinks has already grouped links by dense internal from-node id.
	// Build the contiguous adjacency ranges in one linear scan and represent a
	// node with no outgoing links as the empty range [0, -1].
	int k = 1;
	for (int Node = 1; Node <= no_nodes; Node++)
	{
		int first = k;
		while (k <= number_of_links && Link[k].internal_from_node_id == Node)
			k++;
		if (k > first)
		{
			FirstLinkFrom[Node] = first;
			LastLinkFrom[Node] = k - 1;
		}
		else
		{
			FirstLinkFrom[Node] = 0;
			LastLinkFrom[Node] = -1;
		}
	}

	if (g_tap_log_file == 1)
	{
		for (int Node = 1; Node <= no_nodes; Node++)
			fprintf(logfile, "node_id = %d, FirstLinkFrom = %d, LastLinkFrom = %d \n",
				g_map_node_seq_no_2_external_node_id[Node],
				FirstLinkFrom[Node], LastLinkFrom[Node]);
	}
}


void InitLinks()
{
	char LinksFileName[100] = "link.csv";
	char* InterFlowFileName;
	FILE* LinksFile;
	int k;

	// Use new[] (not calloc) so each link_record is properly constructed.
	// link_record contains std::string members (allowed_uses, geometry) and a
	// non-trivial constructor; calloc would skip construction and assigning to
	// those uninitialized std::string objects later is undefined behavior.
	Link = new link_record[number_of_links + 1];
	ReadLinks();

	InitLinkPointers(LinksFileName);
	UpdateLinkAdditionalCost();
	Load_Movement_Restrictions("movement.csv");
}



void StatusMessage(const char* group, const char* format, ...)
{
	double new_time;
}

// Binary demand format (little-endian, packed):
//   header: char magic[4]="DTAB", int32 version=1, int64 n_records
//   record: int32 o_zone, int32 d_zone, double volume   (16 bytes each)
// Written by `python -m dtalite_qa demand-bin <scenario>`. Reading it avoids
// per-record text parsing -> large OD loads (millions of pairs) go from seconds
// of fscanf to a few bulk freads. Returns false (and leaves tables untouched) if
// the file is missing or not a valid .bin, so the caller can fall back to CSV.
static bool ReadBinaryDemandFile(const std::string& path, int m,
	double*** ODtable, double*** DiffODtable, double*** Seed_ODtable,
	int no_zones, double& total_volume, long long& n_loaded, bool write_seed)
{
	FILE* f = NULL;
	fopen_s(&f, path.c_str(), "rb");
	if (f == NULL)
		return false;

	char magic[4] = { 0,0,0,0 };
	int version = 0;
	long long n = 0;
	if (fread(magic, 1, 4, f) != 4 || memcmp(magic, "DTAB", 4) != 0 ||
		fread(&version, sizeof(int), 1, f) != 1 ||
		fread(&n, sizeof(long long), 1, f) != 1 || version != 1 || n < 0)
	{
		fclose(f);
		return false;
	}

	const size_t REC = 16;                 // int32 + int32 + double, packed
	const long long CHUNK = 1 << 16;       // records per bulk read
	std::vector<unsigned char> buf((size_t)CHUNK * REC);
	total_volume = 0;
	n_loaded = 0;
	long long remaining = n;
	while (remaining > 0)
	{
		long long take = remaining < CHUNK ? remaining : CHUNK;
		size_t got = fread(buf.data(), REC, (size_t)take, f);
		if (got == 0)
			break;
		for (size_t i = 0; i < got; i++)
		{
			unsigned char* p = buf.data() + i * REC;
			int o, d; double vol;
			memcpy(&o, p, 4);
			memcpy(&d, p + 4, 4);
			memcpy(&vol, p + 8, 8);
			// external agency zone ids -> dense internal zone seq (renumbering)
			int o_seq = ZoneSeqFromExternal(o);
			int d_seq = ZoneSeqFromExternal(d);
			if (o_seq >= 1 && d_seq >= 1)
			{
				ODtable[m][o_seq][d_seq] = vol;
				DiffODtable[m][o_seq][d_seq] = vol;
				if (write_seed)
					Seed_ODtable[m][o_seq][d_seq] = vol;   // Seed only needed for ODME / route output
				total_volume += vol;
				n_loaded++;
			}
		}
		remaining -= got;
	}
	fclose(f);
	return true;
}

// Read one mode's demand (binary .bin if demand_format=1, else CSV) into the OD
// tables. Extracted from Read_ODtable so per-mode reads can run in parallel (each
// mode writes a distinct ODtable[m] slice -> no races). write_seed=false skips the
// Seed copy (only ODME / route output consume seed_MDODflow).
static void ReadOneModeDemand(int m, double*** ODtable, double*** DiffODtable,
	double*** Seed_ODtable, int no_zones, bool write_seed)
{
	if (g_demand_format == 1)
	{
		std::string bin = g_mode_type_vector[m].demand_file;
		size_t pos = bin.rfind(".csv");
		bin = (pos != std::string::npos) ? bin.substr(0, pos) + ".bin" : bin + ".bin";

		double bin_total = 0; long long bin_loaded = 0;
		if (ReadBinaryDemandFile(bin, m, ODtable, DiffODtable, Seed_ODtable,
			no_zones, bin_total, bin_loaded, write_seed))
		{
			printf("read demand file %s (binary): %lld pairs, total_volume = %f\n",
				bin.c_str(), bin_loaded, bin_total);
#pragma omp critical
			fprintf(summary_log_file, " mode type = %s, total_volume = %f (binary)\n",
				g_mode_type_vector[m].mode_type.c_str(), bin_total);
			return;
		}
		printf("binary demand %s not found/invalid; falling back to CSV\n", bin.c_str());
	}

	FILE* file = NULL;
	fopen_s(&file, g_mode_type_vector[m].demand_file.c_str(), "r");
	printf("read demand file %s\n", g_mode_type_vector[m].demand_file.c_str());

	if (file == NULL)
	{
		if (g_mode_type_vector[m].demand_file.length() > 0)
		{
			printf("Failed to open demand file %s\n", g_mode_type_vector[m].demand_file.c_str());
			if (number_of_modes == 1)
			{
				for (int o = 1; o <= no_zones; o++)      // dense zone seq: all valid
					for (int d = 1; d <= no_zones; d++)
					{
						ODtable[m][o][d] = 1.0;
						DiffODtable[m][o][d] = 1.0;
						if (write_seed) Seed_ODtable[m][o][d] = 1.0;
					}
			}
		}
		return;   // skip this mode; its OD stays zero
	}

	int o_zone_id, d_zone_id;
	double volume;
	char header[100];
	if (fgets(header, sizeof(header), file) == NULL)
	{
		printf("Failed to read header for %s\n", g_mode_type_vector[m].demand_file.c_str());
		fclose(file);
		return;
	}

	double total_volume = 0;
	long long unknown_zone_rows = 0;
	int result;
	while ((result = fscanf(file, "%d,%d,%lf", &o_zone_id, &d_zone_id, &volume)) != EOF)
	{
		if (result != 3)
		{
			printf("Error reading %s (o=%d)\n", g_mode_type_vector[m].demand_file.c_str(), o_zone_id);
			break;
		}
		// external agency zone ids -> dense internal zone seq (renumbering)
		int o_seq = ZoneSeqFromExternal(o_zone_id);
		int d_seq = ZoneSeqFromExternal(d_zone_id);
		if (o_seq == 0 || d_seq == 0)
		{
			if (unknown_zone_rows < 3)
				printf("WARNING: demand row (o=%d, d=%d) references a zone id not in "
					"node.csv (%s); row skipped.\n",
					o_zone_id, d_zone_id, g_mode_type_vector[m].demand_file.c_str());
			unknown_zone_rows++;
			continue;
		}
		ODtable[m][o_seq][d_seq] = volume;
		DiffODtable[m][o_seq][d_seq] = volume;
		if (write_seed) Seed_ODtable[m][o_seq][d_seq] = volume;
		total_volume += volume;
	}
	if (unknown_zone_rows > 0)
		printf("WARNING: %lld demand rows skipped in %s (zone ids not in node.csv).\n",
			unknown_zone_rows, g_mode_type_vector[m].demand_file.c_str());

	printf(" mode type = %s, total_volume = %f\n", g_mode_type_vector[m].mode_type.c_str(), total_volume);
#pragma omp critical
	fprintf(summary_log_file, " mode type = %s, total_volume = %f\n", g_mode_type_vector[m].mode_type.c_str(), total_volume);
	fclose(file);
}

int Read_ODtable(double*** ODtable, double*** DiffODtable, double*** Seed_ODtable, double*** target_ODtable, int no_zones)
{
	char ch, Coloumn[2], Semicoloumn[2]; /* Reserve room for the '\0' in the fscanf_s. */
	int Orig, Dest, NewOrig, NewDest;
	double Value;

	// current OD demand. Seed_ODtable (seed_MDODflow) is only consumed by ODME and
	// route output, so skip populating it otherwise -> a third fewer OD writes.
	bool write_seed = (g_ODME_mode != 0) || (g_ODME_obs_VMT > 0) || (shortest_path_log_flag != 0);

	// Read modes sequentially. (Reading them in parallel was measured to REGRESS the
	// CSV path ~1.8x on NVTA -- 6 large files + ~5 GB of concurrent dense-array writes
	// are I/O / memory-bandwidth bound, so one-at-a-time has better locality. The
	// binary path is already I/O-light; use the binary format + skip_seed for speed.)
	for (int m = 1; m <= number_of_modes; m++)
	{
		ReadOneModeDemand(m, ODtable, DiffODtable, Seed_ODtable, no_zones, write_seed);
	}


	// baseline OD demand --. diff OD  = current OD - baseline ODdemand 

	if (g_base_demand_mode == 1 && baselinkvolume_loaded_flag == 1)
	{


		for (int m = 1; m <= number_of_modes; m++)
		{
			FILE* file;
			std::string original_filename = g_mode_type_vector[m].demand_file;
			std::string modified_filename;

			// Check if the original filename ends with ".csv" and replace it with "_base.csv"
			size_t pos = original_filename.find(".csv");
			if (pos != std::string::npos) {
				modified_filename = original_filename.substr(0, pos) + "_base.csv";
			}
			else {
				// If the file does not have ".csv", append "_base"
				modified_filename = original_filename + "_base";
			}

			fopen_s(&file, modified_filename.c_str(), "r");
			printf("read base demand file %s\n", modified_filename.c_str());

			if (file == NULL)
			{
				// by default, we can skip this requirement, but if we load baseline link volume we should have base OD demand for consistency 
				break;
			}

			int o_zone_id, d_zone_id;
			double volume;

			// Skip the header line
			char header[100];
			if (fgets(header, sizeof(header), file) == NULL)
			{
				printf("Failed to read header\n");
				return 0;
			}

			int line_count = 0;
			double total_volume_diff = 0;
			// Read the data
			int result;
			while ((result = fscanf(file, "%d,%d,%lf", &o_zone_id, &d_zone_id, &volume)) != EOF)
			{
				if (result == 3)  // we have read all the 3 values correctly
				{
					if (line_count <= 3)
					{
						printf("o_zone_id: %d, d_zone_id: %d, volume: %.4lf\n", o_zone_id, d_zone_id,
							volume);

					}
					int o_seq = ZoneSeqFromExternal(o_zone_id);
					int d_seq = ZoneSeqFromExternal(d_zone_id);
					if (o_seq >= 1 && d_seq >= 1)
					{
					DiffODtable[m][o_seq][d_seq] = ODtable[m][o_seq][d_seq] - volume;  // diff OD demand  = current demand - base demand 
					total_volume_diff += ODtable[m][o_seq][d_seq] - volume;
					}
					line_count++;
				}
				else
				{
					printf("Error reading line %d\n", line_count);
					break;
				}
			}

			printf(" mode type = %s, total_volume_diff = %f\n", g_mode_type_vector[m].mode_type.c_str(), total_volume_diff);

			fclose(file);

		}
	}

	// ODME mode

	if (g_ODME_mode == 1)
	{


		for (int m = 1; m <= number_of_modes; m++)
		{
			FILE* file;
			std::string original_filename = g_mode_type_vector[m].demand_file;
			std::string modified_filename;

			// Check if the original filename ends with ".csv" and replace it with "_base.csv"
			size_t pos = original_filename.find(".csv");
			if (pos != std::string::npos) {
				modified_filename = original_filename.substr(0, pos) + "_target.csv";
			}
			else {
				// If the file does not have ".csv", append "_target"
				modified_filename = original_filename + "_target";
			}

			fopen_s(&file, modified_filename.c_str(), "r");
			printf("read target demand file %s\n", modified_filename.c_str());

			if (file == NULL)
			{

				g_ODME_target_od = -1; 
				std::cerr << "Error:  ODME requires the target demand files " << modified_filename << std::endl;
				std::exit(1);
				// by default, we can skip this requirement, but if we load baseline link volume we should have base OD demand for consistency 

			}

			g_ODME_target_od = 1;

			int o_zone_id, d_zone_id;
			double volume;

			// Skip the header line
			char header[100];
			if (fgets(header, sizeof(header), file) == NULL)
			{
				printf("Failed to read header\n");
				return 0;
			}

			int line_count = 0;
			double total_volume_diff = 0;
			// Read the data
			int result;
			while ((result = fscanf(file, "%d,%d,%lf", &o_zone_id, &d_zone_id, &volume)) != EOF)
			{
				{
					int o_seq = ZoneSeqFromExternal(o_zone_id);
					int d_seq = ZoneSeqFromExternal(d_zone_id);
					if (o_seq >= 1 && d_seq >= 1)
						target_ODtable[m][o_seq][d_seq] = volume;  // target
				}

				if (result == 3)  // we have read all the 3 values correctly
				{
					if (line_count <= 3)
					{
						printf("o_zone_id: %d, d_zone_id: %d, volume: %.4lf\n", o_zone_id, d_zone_id,
							volume);

					}


					line_count++;
				}
				else
				{
					printf("Error reading line %d\n", line_count);
					break;
				}
			}

			printf(" mode type = %s, total_volume_diff = %f\n", g_mode_type_vector[m].mode_type.c_str(), total_volume_diff);

			fclose(file);

		}
	}
	return 1;
}

double Link_Travel_Time(int k, double* Volume)
{
	double IncomingDemand = Volume[k] / fmax(0.01, Link[k].lanes) / fmax(0.001, demand_period_ending_hours - demand_period_starting_hours) / fmax(0.0001, Link[k].VDF_plf);

	if (Link[k].VDF_type == 1)
	{
		// Spiess conic congestion function (per-lane V/C):
		//   t = t0 * ( 2 + sqrt(a^2 (1-x)^2 + b^2) - a(1-x) - b ),  x = V/C
		// a (Conic_a) per functional type; b (Conic_b) per Spiess. Asymptotically
		// linear (unlike BPR power). Uses Lane_Capacity for a true per-lane V/C.
		double x = IncomingDemand / fmax(0.1, Link[k].Lane_Capacity);
		// Conic a/b: prefer explicit conic_a/conic_b columns; otherwise the
		// staged convention stores conic a in vdf_alpha and b in vdf_beta.
		double a = (Link[k].Conic_a > 0.0) ? Link[k].Conic_a : Link[k].VDF_Alpha;
		double b = (Link[k].Conic_b > 0.0) ? Link[k].Conic_b : Link[k].VDF_Beta;
		double om = 1.0 - x;
		Link[k].Travel_time =
			Link[k].FreeTravelTime * (2.0 + sqrt(a * a * om * om + b * b) - a * om - b);
	}
	else if (Link[k].VDF_type == 2)
	{
		// QVDF (Queue-VDF) as the assignment cost: period-average congested
		// travel time from the queue speed model (calibrated Q_cd/Q_n/Q_cp/Q_s +
		// VDF_Alpha/Beta + Cutoff_Speed). Same period-average speed the kernel
		// reports; monotonic in D/C, so valid as a link cost. (Cheap part only;
		// the time-sliced queue profile stays in Link_QueueVDF for reporting.)
		double DOC = IncomingDemand / fmax(0.1, Link[k].Lane_Capacity);
		double cong_ref = (DOC < 1.0)
			? ((1.0 - DOC) * Link[k].free_speed + DOC * Link[k].Cutoff_Speed)
			: Link[k].Cutoff_Speed;
		double avg_queue_speed = cong_ref / (1.0 + Link[k].VDF_Alpha * pow(DOC, Link[k].VDF_Beta));
		double P = Link[k].Q_cd * pow(DOC, Link[k].Q_n);
		double H = fmax(0.001, demand_period_ending_hours - demand_period_starting_hours);
		double avg_period_speed = (P > H)
			? avg_queue_speed
			: (P / H) * avg_queue_speed + (1.0 - P / H) * (cong_ref + Link[k].free_speed) / 2.0;
		Link[k].Travel_time = Link[k].length / fmax(0.1, avg_period_speed) * 60.0;
	}
	else if (Link[k].VDF_type == 3)
	{
		// BPR2 (AequilibraE): the BPR exponent DOUBLES above capacity, giving a
		// steeper over-saturation penalty. x = per-lane V/C.
		double x = IncomingDemand / fmax(0.1, Link[k].Lane_Capacity);
		double e = (x <= 1.0) ? Link[k].VDF_Beta : 2.0 * Link[k].VDF_Beta;
		Link[k].Travel_time = Link[k].FreeTravelTime * (1.0 + Link[k].VDF_Alpha * pow(x, e));
	}
	else if (Link[k].VDF_type == 4)
	{
		// INRETS (AequilibraE). vdf_alpha is the uncongested ratio (~0.9-1.0).
		//   x<=1: t0*(1.1 - a*x)/(1.1 - x);   x>1: t0*((1.1 - a)/0.1)*x^2.
		double x = IncomingDemand / fmax(0.1, Link[k].Lane_Capacity);
		double a = Link[k].VDF_Alpha;
		Link[k].Travel_time = (x <= 1.0)
			? Link[k].FreeTravelTime * (1.1 - a * x) / fmax(0.05, 1.1 - x)
			: Link[k].FreeTravelTime * ((1.1 - a) / 0.1) * x * x;
	}
	else if (Link[k].VDF_type == 5)
	{
		// Akcelik time-dependent (VDOT-allowed, AequilibraE):
		//   t = t0 + alpha*( z + sqrt(z^2 + beta*x) ),  z = x - 1   (alpha in minutes).
		double x = IncomingDemand / fmax(0.1, Link[k].Lane_Capacity);
		double z = x - 1.0;
		Link[k].Travel_time = Link[k].FreeTravelTime
			+ Link[k].VDF_Alpha * (z + sqrt(z * z + Link[k].VDF_Beta * x));
	}
	else if (Link[k].VDF_type == 6)
	{
		// SANDAG-style: BPR running time + Webster uniform signal delay.
		//   d = 0.5*C*(1 - g/C)^2 / (1 - min(1,x)*g/C)   [C=cycle_length sec -> min].
		double x = IncomingDemand / fmax(0.1, Link[k].Lane_Capacity);
		double bpr = Link[k].FreeTravelTime * (1.0 + Link[k].VDF_Alpha * pow(x, Link[k].VDF_Beta));
		double gc = fmin(0.95, fmax(0.05, Link[k].green_ratio));
		double C = fmax(0.0, Link[k].cycle_length);
		double d_min = 0.5 * C * (1.0 - gc) * (1.0 - gc) / fmax(0.05, 1.0 - fmin(1.0, x) * gc) / 60.0;
		Link[k].Travel_time = bpr + d_min;
	}
	else if (Link[k].VDF_type == 7)
	{
		// SCAG piecewise BPR (Validation Report Table 16-2). The exponent is set to
		// SCAG_UNCONGESTED_BETA (=4.0) below capacity and to the calibrated per-link
		// VDF_Beta (5/6/8, by facility x posted speed x area type) at/above capacity.
		// The two branches meet at x=1 (both give t0*(1+alpha)), so t is continuous.
		// VDF_Alpha is the table alpha (1.00 freeway, 0.80 others).
		double x = IncomingDemand / fmax(0.1, Link[k].Lane_Capacity);
		double e = (x <= 1.0) ? SCAG_UNCONGESTED_BETA : Link[k].VDF_Beta;
		Link[k].Travel_time = Link[k].FreeTravelTime * (1.0 + Link[k].VDF_Alpha * pow(x, e));
	}
	else if (Link[k].VDF_type == 8)
	{
		// SCAG freeway on-ramp meter delay (Validation Report, facility 82/84):
		//   t = L/FFS + [ (PLPHx/120) * 5.0 * (1 + x)^8 ] / 60   (hours -> stored in min)
		// PLPHx = per-lane-per-hour flow = IncomingDemand; x = per-lane V/C. The bracket
		// is a metered-queue delay that grows with the 8th power of the demand ratio.
		double x = IncomingDemand / fmax(0.1, Link[k].Lane_Capacity);
		double ramp_delay_hr = (IncomingDemand / 120.0) * 5.0 * pow(1.0 + x, 8.0) / 60.0;
		Link[k].Travel_time = Link[k].FreeTravelTime + ramp_delay_hr * 60.0;
	}
	else
	{
		// BPR. IncomingDemand above is a PER-LANE hourly rate, so the V/C ratio
		// must divide by the per-lane capacity (Lane_Capacity), not the total
		// link capacity (Link_Capacity = lanes * Lane_Capacity). Using
		// Link_Capacity understated D/C by a factor of #lanes. This now matches
		// the conic (VDF_type==1) and QVDF (VDF_type==2) branches above.
		double x_voc = IncomingDemand / fmax(0.1, Link[k].Lane_Capacity);
		// ARC modified BPR: fftt*(1 + A*(V/C) + alpha*(V/C)^beta). VDF_A defaults to
		// 0, which recovers the standard BPR exactly (other networks unaffected).
		Link[k].Travel_time =
			Link[k].FreeTravelTime * (1.0 + Link[k].VDF_A * x_voc + Link[k].VDF_Alpha * pow(x_voc, Link[k].VDF_Beta));
	}

	// MAG volume delay function: add the per-length term (1.4*L) to the congested
	// time T_c, for ALL facility types. length is in miles -> minutes. Off when
	// g_added_delay_per_mile == 0 (other networks unaffected).
	if (g_added_delay_per_mile != 0.0)
		Link[k].Travel_time += g_added_delay_per_mile * Link[k].length;

	if (Link[k].Travel_time < 0)
		Link[k].Travel_time = 0;
	Link[k].BPR_TT = Link[k].Travel_time;

	return (Link[k].Travel_time);
}


double Link_QueueVDF(int k, double Volume, double& IncomingDemand, double& DOC, double& P, double& t0, double& t2, double& t3, double& vt2, double& Q_mu, double& Q_gamma, double& congestion_ref_speed,
	double& avg_queue_speed, double& avg_QVDF_period_speed, double& Severe_Congestion_P, double model_speed[300])
{

	IncomingDemand = Volume / fmax(0.01, Link[k].lanes) / fmax(0.001, demand_period_ending_hours - demand_period_starting_hours) / fmax(0.0001, Link[k].VDF_plf);
	DOC = IncomingDemand / fmax(0.1, Link[k].Lane_Capacity);

	double Travel_time =
		Link[k].FreeTravelTime * (1.0 + Link[k].VDF_Alpha * (pow(DOC, Link[k].VDF_Beta)));

	congestion_ref_speed = Link[k].Cutoff_Speed;
	if (DOC < 1)
		congestion_ref_speed = (1 - DOC) * Link[k].free_speed + DOC * Link[k].Cutoff_Speed;


	//step 3.2 calculate speed from VDF based on D/C ratio
	avg_queue_speed = congestion_ref_speed / (1.0 + Link[k].VDF_Alpha * pow(DOC, Link[k].VDF_Beta));


	P = Link[k].Q_cd * pow(DOC, Link[k].Q_n);  // applifed for both uncongested and congested conditions

	double H = demand_period_ending_hours - demand_period_starting_hours;

	if (P > H)
		avg_QVDF_period_speed = avg_queue_speed;
	else
		avg_QVDF_period_speed = P / H * avg_queue_speed + (1.0 - P / H) * (congestion_ref_speed + Link[k].free_speed) / 2.0;


	Link[k].QVDF_TT = Link[k].length / fmax(0.1, avg_QVDF_period_speed) * 60.0;






	double base = Link[k].Q_cp * pow(P, Link[k].Q_s) + 1.0;
	vt2 = Link[k].Cutoff_Speed / fmax(0.001, base);

	// Observed t0/t3 determine only how the analytical duration P is divided
	// around t2. No additional user-selectable QVDF parameter is introduced.
	// Missing or unusable endpoint data retains the historical symmetric split.
	double observed_left_fraction = 0.5;
	if (std::isfinite(Link[k].QVDF_t0_observed) &&
		std::isfinite(Link[k].QVDF_t2) &&
		std::isfinite(Link[k].QVDF_t3_observed) &&
		Link[k].QVDF_t0_observed < Link[k].QVDF_t2 &&
		Link[k].QVDF_t2 < Link[k].QVDF_t3_observed)
	{
		double observed_duration =
			Link[k].QVDF_t3_observed - Link[k].QVDF_t0_observed;
		if (observed_duration > 1e-6)
		{
			observed_left_fraction =
				(Link[k].QVDF_t2 - Link[k].QVDF_t0_observed) /
				fmax(1e-6, observed_duration);
			observed_left_fraction =
				fmin(0.95, fmax(0.05, observed_left_fraction));
		}
	}

	// Keep the reported/profiled congestion window inside the configured
	// assignment-period band. P remains the analytical QVDF duration; after
	// clamping, the visible t0--t3 span can be shorter than P.
	t0 = fmax(
		demand_period_starting_hours,
		t2 - observed_left_fraction * P);
	t3 = fmin(
		demand_period_ending_hours,
		t2 + (1.0 - observed_left_fraction) * P);

	Q_mu = std::min(Link[k].Lane_Capacity, IncomingDemand / std::max(0.01, P));

	//use  as the lower speed compared to 8/15 values for the congested states
	double RTT = Link[k].length / fmax(0.01, congestion_ref_speed);
	double wt2 = Link[k].length / vt2 - RTT; // in hour


	//step 5 compute gamma parameter is controlled by the maximum queue
	Q_gamma = wt2 * 64 * Q_mu / pow(P, 4);  // because q_tw = w*mu =1/4 * gamma (P/2)^4, => 1/vt2 * mu = 1/4 * gamma  * (P/2)^4


	double td_w = 0;
	//step scan the entire analysis period
	Severe_Congestion_P = 0;
	double boundary_speed = fmax(congestion_ref_speed, avg_queue_speed);
	auto smoothstep01 = [](double value)
		{
			double position = fmin(1.0, fmax(0.0, value));
			return position * position * (3.0 - 2.0 * position);
		};


	for (int t_in_min = demand_period_starting_hours * 60; t_in_min <= demand_period_ending_hours * 60; t_in_min += 5)  // 5 min interval
	{
		int t_interval = t_in_min / 5;
		double t = t_in_min / 60.0;  // t in hour
		double td_queue = 0;
		double td_speed = 0;
		model_speed[t_interval] = Link[k].free_speed;

		if (t0 <= t && t <= t3)  // within the band-limited congestion window
		{
			// Normalize an asymmetric beta-shaped queue curve so its maximum
			// remains exactly at the observed t2 after either window edge is
			// clamped. With no clamping, peak_fraction=0.5 and this reduces to
			// the historical quartic 1/4*gamma*(t-t0)^2*(t-t3)^2.
			double window_span = t3 - t0;
			double queue_shape = 0.0;
			if (window_span > 1e-9)
			{
				double x = fmin(1.0, fmax(0.0, (t - t0) / window_span));
				double peak_fraction =
					fmin(1.0, fmax(0.0, (t2 - t0) / window_span));
				if (peak_fraction <= 1e-9)
					queue_shape = pow(1.0 - x, 4.0);
				else if (peak_fraction >= 1.0 - 1e-9)
					queue_shape = pow(x, 4.0);
				else
				{
					double left_exponent = 4.0 * peak_fraction;
					double right_exponent = 4.0 * (1.0 - peak_fraction);
					queue_shape =
						pow(x / peak_fraction, left_exponent) *
						pow((1.0 - x) / (1.0 - peak_fraction), right_exponent);
				}
			}
			else if (fabs(t - t2) <= 1e-9)
				queue_shape = 1.0;

			queue_shape = fmin(1.0, fmax(0.0, queue_shape));
			double peak_queue = wt2 * Q_mu;
			td_queue = peak_queue * queue_shape;
			td_w = td_queue / fmax(0.001, Q_mu);
			//L/[(w(t)+RTT_in_hour]
			td_speed = Link[k].length / (td_w + RTT);
		}
		else if (t < t0) //outside
		{
			td_queue = 0;
			double factor = (t - demand_period_starting_hours) / fmax(0.001, t0 - demand_period_starting_hours);
			double smooth_factor = smoothstep01(factor);
			td_speed = (1.0 - smooth_factor) * Link[k].free_speed
				+ smooth_factor * boundary_speed;
		}
		else  // t> t3
		{
			td_queue = 0;
			double factor = (t - t3) / fmax(0.001, demand_period_ending_hours - t3);
			double smooth_factor = smoothstep01(factor);
			td_speed = (1.0 - smooth_factor) * boundary_speed
				+ smooth_factor * Link[k].free_speed;
		}

		// dtalog.output() << "td_queue t" << t << " =  " << td_queue << ", speed =" << td_speed << '\n';
		// g_DTA_log_file << "td_queue t" << t << " =  " << td_queue << ", speed =" << td_speed << '\n';



		if (t_in_min <= 410)
		{
			int idebug = 1;
		}
		double td_flow = 0; // default: get_volume_from_speed(td_speed, vf, k_critical, s3_m);
		model_speed[t_interval] = td_speed;

		if (td_speed < Link[k].free_speed * 0.5)
			Severe_Congestion_P += 5.0 / 60;  // 5 min interval
	}

	// Anchor optional observed speeds to the first and last emitted samples.
	// The CSV profile is half-open [period_start, period_end), so its last sample
	// is period_end - 5 minutes. With an observed t2, a low anchor connects
	// directly to vt2 with cubic smoothstep so the raw shoulder cannot create an
	// artificial reversal. An anchor clearly between vt2 and the modeled boundary
	// instead uses a monotone cubic Hermite splice: zero slope at the observation
	// and the raw QVDF slope at the inward join. The historical blend remains the
	// fallback, and every path preserves the analytical scalar values.
	bool has_start_speed = std::isfinite(Link[k].QVDF_start_speed_observed);
	bool has_end_speed = std::isfinite(Link[k].QVDF_end_speed_observed);
	if (has_start_speed || has_end_speed)
	{
		int profile_start_min = static_cast<int>(demand_period_starting_hours * 60);
		int period_end_min = static_cast<int>(demand_period_ending_hours * 60);
		int profile_last_min = std::max(profile_start_min, period_end_min - 5);
		double profile_start_hour = profile_start_min / 60.0;
		double profile_last_hour = profile_last_min / 60.0;
		double anchor_pivot = fmin(
			profile_last_hour, fmax(profile_start_hour, t2));
		bool has_observed_t2 = std::isfinite(Link[k].QVDF_t2);
		double anchor_margin = fmax(
			2.0, 0.10 * fmax(0.0, boundary_speed - vt2));
		auto use_monotone_hermite = [&](double anchor_speed)
			{
				return has_observed_t2 &&
					anchor_speed > vt2 + anchor_margin &&
					anchor_speed < boundary_speed;
			};
		auto use_low_anchor_connector = [&](double anchor_speed)
			{
				return has_observed_t2 &&
					anchor_speed < boundary_speed &&
					anchor_speed <= vt2 + anchor_margin;
			};
		auto monotone_hermite_value = [](
			double start_speed, double end_speed,
			double start_slope, double end_slope,
			double span_hours, double factor)
			{
				double position = fmin(1.0, fmax(0.0, factor));
				double position2 = position * position;
				double position3 = position2 * position;
				double h00 = 2.0 * position3 - 3.0 * position2 + 1.0;
				double h10 = position3 - 2.0 * position2 + position;
				double h01 = -2.0 * position3 + 3.0 * position2;
				double h11 = position3 - position2;
				double value =
					h00 * start_speed + h10 * span_hours * start_slope +
					h01 * end_speed + h11 * span_hours * end_slope;
				return fmin(fmax(start_speed, end_speed),
					fmax(fmin(start_speed, end_speed), value));
			};

		if (has_start_speed)
		{
			double start_speed = Link[k].QVDF_start_speed_observed;
			bool used_low_anchor_connector = false;
			if (use_low_anchor_connector(start_speed))
			{
				double left_span = anchor_pivot - profile_start_hour;
				if (left_span > 1e-9)
				{
					for (int t_in_min = profile_start_min;
						t_in_min <= profile_last_min &&
						t_in_min / 60.0 <= anchor_pivot;
						t_in_min += 5)
					{
						double t = t_in_min / 60.0;
						double smooth_factor = smoothstep01(
							(t - profile_start_hour) / left_span);
						model_speed[t_in_min / 5] =
							(1.0 - smooth_factor) * start_speed +
							smooth_factor * vt2;
					}
				}
				used_low_anchor_connector = true;
			}
			bool used_monotone_hermite = false;
			if (!used_low_anchor_connector &&
				use_monotone_hermite(start_speed))
			{
				for (int join_min = profile_start_min + 5;
					join_min <= profile_last_min &&
					join_min / 60.0 <= anchor_pivot;
					join_min += 5)
				{
					double join_speed = model_speed[join_min / 5];
					if (join_speed >= start_speed)
						continue;
					double span_hours = (join_min - profile_start_min) / 60.0;
					double secant_slope =
						(join_speed - start_speed) / span_hours;
					int next_min = std::min(profile_last_min, join_min + 5);
					double slope_span_hours = (next_min - join_min) / 60.0;
					double join_slope = slope_span_hours > 1e-9
						? (model_speed[next_min / 5] - join_speed) /
							slope_span_hours
						: 0.0;
					if (join_slope > 0.0 ||
						join_slope / secant_slope > 3.0)
						continue;

					for (int t_in_min = profile_start_min;
						t_in_min <= join_min; t_in_min += 5)
					{
						double factor = static_cast<double>(
							t_in_min - profile_start_min) /
							(join_min - profile_start_min);
						model_speed[t_in_min / 5] = monotone_hermite_value(
							start_speed, join_speed, 0.0, join_slope,
							span_hours, factor);
					}
					used_monotone_hermite = true;
					break;
				}
			}

			if (!used_low_anchor_connector && !used_monotone_hermite)
			{
				double left_span = anchor_pivot - profile_start_hour;
				for (int t_in_min = profile_start_min;
					t_in_min <= profile_last_min &&
					t_in_min / 60.0 <= anchor_pivot;
					t_in_min += 5)
				{
					int t_interval = t_in_min / 5;
					double t = t_in_min / 60.0;
					double observed_weight = left_span > 1e-9
						? 1.0 - smoothstep01(
							(t - profile_start_hour) / left_span)
						: (t_in_min == profile_start_min ? 1.0 : 0.0);
					model_speed[t_interval] =
						observed_weight * start_speed +
						(1.0 - observed_weight) * model_speed[t_interval];
				}
			}
			model_speed[profile_start_min / 5] =
				start_speed;
		}

		if (has_end_speed)
		{
			double end_speed = Link[k].QVDF_end_speed_observed;
			bool used_low_anchor_connector = false;
			if (use_low_anchor_connector(end_speed))
			{
				double right_span = profile_last_hour - anchor_pivot;
				if (right_span > 1e-9)
				{
					for (int t_in_min = profile_start_min;
						t_in_min <= profile_last_min; t_in_min += 5)
					{
						double t = t_in_min / 60.0;
						if (t < anchor_pivot)
							continue;
						double smooth_factor = smoothstep01(
							(t - anchor_pivot) / right_span);
						model_speed[t_in_min / 5] =
							(1.0 - smooth_factor) * vt2 +
							smooth_factor * end_speed;
					}
				}
				used_low_anchor_connector = true;
			}
			bool used_monotone_hermite = false;
			if (!used_low_anchor_connector &&
				use_monotone_hermite(end_speed))
			{
				for (int join_min = profile_last_min - 5;
					join_min >= profile_start_min &&
					join_min / 60.0 >= anchor_pivot;
					join_min -= 5)
				{
					double join_speed = model_speed[join_min / 5];
					if (join_speed >= end_speed)
						continue;
					double span_hours = (profile_last_min - join_min) / 60.0;
					double secant_slope =
						(end_speed - join_speed) / span_hours;
					int adjacent_min = join_min > anchor_pivot * 60.0 + 1e-9
						? std::max(profile_start_min, join_min - 5)
						: std::min(profile_last_min, join_min + 5);
					double slope_span_hours = fabs(adjacent_min - join_min) / 60.0;
					double join_slope = 0.0;
					if (slope_span_hours > 1e-9)
					{
						join_slope = adjacent_min < join_min
							? (join_speed - model_speed[adjacent_min / 5]) /
								slope_span_hours
							: (model_speed[adjacent_min / 5] - join_speed) /
								slope_span_hours;
					}
					if (join_slope < 0.0 ||
						join_slope / secant_slope > 3.0)
						continue;

					for (int t_in_min = join_min;
						t_in_min <= profile_last_min; t_in_min += 5)
					{
						double factor = static_cast<double>(
							t_in_min - join_min) /
							(profile_last_min - join_min);
						model_speed[t_in_min / 5] = monotone_hermite_value(
							join_speed, end_speed, join_slope, 0.0,
							span_hours, factor);
					}
					used_monotone_hermite = true;
					break;
				}
			}

			if (!used_low_anchor_connector && !used_monotone_hermite)
			{
				double right_span = profile_last_hour - anchor_pivot;
				for (int t_in_min = profile_start_min;
					t_in_min <= profile_last_min; t_in_min += 5)
				{
					double t = t_in_min / 60.0;
					if (t < anchor_pivot)
						continue;
					int t_interval = t_in_min / 5;
					double observed_weight = right_span > 1e-9
						? smoothstep01((t - anchor_pivot) / right_span)
						: (t_in_min == profile_last_min ? 1.0 : 0.0);
					model_speed[t_interval] =
						(1.0 - observed_weight) * model_speed[t_interval] +
						observed_weight * end_speed;
				}
			}
			model_speed[profile_last_min / 5] =
				end_speed;
		}

		// Severe congestion is profile-derived, so recompute it from the anchored,
		// emitted samples whenever either boundary observation changes the profile.
		Severe_Congestion_P = 0.0;
		for (int t_in_min = profile_start_min;
			t_in_min <= profile_last_min; t_in_min += 5)
		{
			if (model_speed[t_in_min / 5] < Link[k].free_speed * 0.5)
				Severe_Congestion_P += 5.0 / 60;
		}
	}

	return P;
}

double Link_Travel_Time_Integral(int k, double* Volume)
{
	double IncomingDemand = Volume[k] / fmax(0.01, Link[k].lanes) / fmax(0.001, demand_period_ending_hours - demand_period_starting_hours) / fmax(0.0001, Link[k].VDF_plf);
	double integral = 0;
	if (Link[k].VDF_Beta >= 0.0)
		integral += IncomingDemand + (Volume[k] * Link[k].FreeTravelTime *
			(1.0 + (Link[k].BoverC / (Link[k].VDF_Beta + 1)) * pow(IncomingDemand, Link[k].VDF_Beta + 1)));

	return integral;
}

double Link_Travel_Time_Der(int k, double* Volume)
{
	if (Link[k].VDF_Beta == 0.0)
		return 0.0;
	else
		return (Link[k].FreeTravelTime * Link[k].BoverC * Link[k].VDF_Beta *
			pow(Volume[k], (Link[k].VDF_Beta - 1)));
}

double AdditionalCost(int k, int m)
{
	double AddCost = 0;

	AddCost = (Link[k].mode_Toll[m] + Link[k].length * g_mode_type_vector[m].op_cost) / g_mode_type_vector[m].vot * 60.0;

	return AddCost;
}

double Link_GenCost(int k, double* Volume)
{
	return (Link[k].mode_AdditionalCost[1] + Link_Travel_Time(k, Volume));
}

double LinkCost_Integral(int k, double* Volume)
{
	return (Link[k].mode_AdditionalCost[1] * Volume[k] + Link_Travel_Time_Integral(k, Volume));
}

double Link_GenCostDer(int k, double* Volume)
{
	return (Link_Travel_Time_Der(k, Volume));
}

/* External functions */

void ClearVolume(double* VolumeArray)
{
	int k;
	for (k = 1; k <= number_of_links; k++)
		VolumeArray[k] = 0.0;
}

void VolumeDifference(double* Volume1, double* Volume2, double* Difference)
{  // SubVolume, MainVolume
	int k;
	for (k = 1; k <= number_of_links; k++)
	{
		Difference[k] = Volume1[k] - Volume2[k];

		for (int m = 1; m <= number_of_modes; m++)
		{
			Link[k].mode_SDVolume[m] = Link[k].mode_SubVolume[m] - Link[k].mode_MainVolume[m];

			if (fabs(Difference[k] - Link[k].mode_SDVolume[m]) > 0.01)
			{
				printf("");
			}
		}
	}
}

void UpdateVolume(double* MainVolume, double* SDVolume, double Lambda)
{
	int k;
	for (k = 1; k <= number_of_links; k++)
	{
		MainVolume[k] += Lambda * SDVolume[k];
	}


	for (int k = 1; k <= number_of_links; k++)
	{
		for (int m = 1; m <= number_of_modes; m++)
		{
			Link[k].mode_MainVolume[m] += Lambda * Link[k].mode_SDVolume[m];
		}
	}
}

void UpdateLinkAdditionalCost(void)
{
	int k;

	for (k = 1; k <= number_of_links; k++)
		for (int m = 1; m <= number_of_modes; m++)
			AdditionalCost(k, m);
}

// (Bi)conjugate Frank-Wolfe: rewrite the FW search direction SDVolume (= y - x) with a
// conjugate direction. The new direction is built from a CONVEX-combination auxiliary
//   CFW: aux = a*s1 + (1-a)*y        (one prior auxiliary s1)
//   BFW: aux = b0*y + b1*s1 + b2*s2  (two prior auxiliaries)
// with coefficients chosen so (aux - x) is H-conjugate to the previous direction(s),
// H = diag(dCost/dV). Because s1,s2,y are all feasible (convex combos of all-or-nothing
// solutions), x + lambda*(aux - x) stays feasible for lambda in [0,1] -- so the existing
// bisection line search is safe. Falls back to plain FW when degenerate / non-descent.
void ApplyConjugateFW(double* MainVolume, double* SubVolume, double* SDVolume, int iteration_no)
{
	if (g_assignment_method <= 0) return;
	int n = number_of_links;
	if (g_cfw_s1 == NULL) {
		g_cfw_s1 = (double*)Alloc_1D(n + 1, sizeof(double));
		g_cfw_s2 = (double*)Alloc_1D(n + 1, sizeof(double));
		g_cfw_aux = (double*)Alloc_1D(n + 1, sizeof(double));
		g_cfw_H = (double*)Alloc_1D(n + 1, sizeof(double));
	}
	if (iteration_no <= 1) g_cfw_nstored = 0;   // reset at the start of an assignment

	// diagonal Hessian via a forward difference on the link cost (valid for any VDF).
	for (int k = 1; k <= n; k++) {
		double v = MainVolume[k];
		double t0 = Link_Travel_Time(k, MainVolume);
		double dv = fmax(1.0, fabs(v) * 1e-4);
		MainVolume[k] = v + dv;
		double t1 = Link_Travel_Time(k, MainVolume);
		MainVolume[k] = v;
		Link[k].Travel_time = t0;               // restore (Link_Travel_Time mutates it)
		g_cfw_H[k] = fmax(0.0, (t1 - t0) / dv);
	}

	bool built = false;
	if (g_cfw_nstored >= 2 && g_assignment_method == 2) {
		// BFW. direction = b0*d + b1*e1 + b2*e2 with d=y-x, e1=s1-x, e2=s2-x, b0=1-b1-b2.
		// => direction = d + b1*(e1-d) + b2*(e2-d); conjugate to e1 and e2.
		double A11 = 0, A12 = 0, A21 = 0, A22 = 0, r1 = 0, r2 = 0;
		for (int k = 1; k <= n; k++) {
			double H = g_cfw_H[k], d = SDVolume[k];
			double e1 = g_cfw_s1[k] - MainVolume[k], e2 = g_cfw_s2[k] - MainVolume[k];
			double p1 = e1 - d, p2 = e2 - d;
			A11 += H * e1 * p1; A12 += H * e1 * p2;
			A21 += H * e2 * p1; A22 += H * e2 * p2;
			r1 -= H * e1 * d;   r2 -= H * e2 * d;
		}
		double det = A11 * A22 - A12 * A21;
		if (fabs(det) > 1e-10) {
			double b1 = (r1 * A22 - A12 * r2) / det;
			double b2 = (A11 * r2 - r1 * A21) / det;
			double b0 = 1.0 - b1 - b2;
			if (b0 >= 0 && b1 >= 0 && b2 >= 0) {     // convex -> feasible
				for (int k = 1; k <= n; k++)
					g_cfw_aux[k] = b0 * SubVolume[k] + b1 * g_cfw_s1[k] + b2 * g_cfw_s2[k];
				built = true;
			}
		}
	}
	if (!built && g_cfw_nstored >= 1) {
		// CFW. aux = a*s1 + (1-a)*y; conjugacy of (aux-x) to e1=s1-x gives
		//   a = (e1.H.d) / (e1.H.d - e1.H.e1),  d = y-x.
		double e1Hd = 0, e1He1 = 0;
		for (int k = 1; k <= n; k++) {
			double H = g_cfw_H[k], d = SDVolume[k], e1 = g_cfw_s1[k] - MainVolume[k];
			e1Hd += H * e1 * d; e1He1 += H * e1 * e1;
		}
		double denom = e1Hd - e1He1;
		double a = (fabs(denom) > 1e-10) ? e1Hd / denom : 0.0;
		if (a < 0) a = 0; if (a > 1.0 - 1e-4) a = 1.0 - 1e-4;
		for (int k = 1; k <= n; k++)
			g_cfw_aux[k] = a * g_cfw_s1[k] + (1.0 - a) * SubVolume[k];
		built = true;
	}
	if (!built)
		for (int k = 1; k <= n; k++) g_cfw_aux[k] = SubVolume[k];   // first iter: aux = y (FW)

	// descent guard: require Sum GenCost*(aux - x) < 0; else fall back to FW (aux = y).
	double dd = 0;
	for (int k = 1; k <= n; k++) dd += Link[k].GenCost * (g_cfw_aux[k] - MainVolume[k]);
	if (dd >= 0)
		for (int k = 1; k <= n; k++) g_cfw_aux[k] = SubVolume[k];

	// shift the auxiliary history and set the new search direction = aux - x.
	for (int k = 1; k <= n; k++) {
		g_cfw_s2[k] = g_cfw_s1[k];
		g_cfw_s1[k] = g_cfw_aux[k];
		SDVolume[k] = g_cfw_aux[k] - MainVolume[k];
	}
	if (g_cfw_nstored < 2) g_cfw_nstored++;
}

double UpdateLinkCost(double* MainVolume)
{
	int k;
	double system_wide_travel_time = 0;

	for (k = 1; k <= number_of_links; k++)
	{
		Link[k].Travel_time = Link_Travel_Time(k, MainVolume);

		Link[k].GenCost = Link_GenCost(k, MainVolume);
		system_wide_travel_time += (MainVolume[k] * Link[k].Travel_time);

		// The system total must use the SAME generalized cost (travel time +
		// additional/toll cost) that the shortest-path benchmark uses in
		// FindMinCostRoutes (route cost = travel time + mode_AdditionalCost,
		// flow weighted by pce). Without the additional-cost term here, the
		// relative gap becomes negative whenever tolls/additional costs are
		// present, because the benchmark includes tolls but the system total
		// did not. mode_MainVolume is the pure per-mode flow; weight by pce to
		// match. For toll-free networks mode_AdditionalCost == 0, so this term
		// vanishes and the gap is byte-for-byte unchanged.
		for (int m = 1; m <= number_of_modes; m++)
			system_wide_travel_time +=
				Link[k].mode_MainVolume[m] * g_mode_type_vector[m].pce * Link[k].mode_AdditionalCost[m];
	}

	return system_wide_travel_time;
}

void UpdateLinkCostDer(double* MainVolume)
{
	int k;

	for (k = 1; k <= number_of_links; k++)
	{
		Link[k].GenCostDer = Link_GenCostDer(k, MainVolume);
	}
}

void GetLinkTravelTimes(double* Volume, double* TravelTime)
{
	int k;

	for (k = 1; k <= number_of_links; k++)
	{
		TravelTime[k] = Link_Travel_Time(k, Volume);
	}
}

double TotalLinkCost(double* Volume)
{
	int k;
	double Sum = 0;

	for (k = 1; k <= number_of_links; k++)
		Sum += Link[k].GenCost * Volume[k];
	return Sum;
}

double OFscale = 1.0;

double OF_Links(double* MainVolume)
{
	int k;
	double Sum = 0;

	for (k = 1; k <= number_of_links; k++)
		Sum += LinkCost_Integral(k, MainVolume);

	return Sum / OFscale;
}

double OF_LinksDirectionalDerivative(double* MainVolume, double* SDVolume, double Lambda)
{
	//
	//purpose:
	//    this function calculates the directional derivative of the objective function with respect to the step size lambda.in optimization, the directional derivative indicates how the objective function changes as you move in a specific direction.
	//
	//        parameters :
	//        mainvolume : current flow volumes on each link.
	//        sdvolume : search direction volumes(difference between aon assignment and current flows).
	//        lambda : step size parameter indicating how far to move along the search direction.

		//Return Normalized Directional Derivative :
		//Normalize the sum by OFscale to get the directional derivative :
		//DirectionalDerivative
		//    =
		//    LinkCostSum
		//    OFscale
		//    DirectionalDerivative =
		//    OFscale
		//    LinkCostSum
		//    ?

		//    Mathematical Interpretation :
		//The directional derivative represents how the total system cost changes as you adjust the flow along the search direction by Lambda.In the context of traffic assignment, it helps determine whether increasing or decreasing Lambda will reduce the overall congestion and travel time.

	int k;
	double* Volume;
	double LinkCostSum = 0;

	Volume = (double*)Alloc_1D(number_of_links, sizeof(double));

	for (k = 1; k <= number_of_links; k++)
	{
		Volume[k] = MainVolume[k] + Lambda * SDVolume[k];
	}
	for (k = 1; k <= number_of_links; k++)
	{
		LinkCostSum += Link_GenCost(k, Volume) * SDVolume[k];

		//   LinkCostSum += Link_Travel_Time_Integral(k, Volume) * SDVolume[k];

	}

	free(Volume);
	return LinkCostSum / OFscale;
}


double Sum_ODtable(double*** ODtable, double* total_o_table, int no_zones)
{
	int Orig, Dest;
	double sum = 0.0;

	for (Orig = 1; Orig <= no_zones; Orig++)
		total_o_table[Orig] = 0;

	for (int m = 1; m <= number_of_modes; m++)
		for (Orig = 1; Orig <= no_zones; Orig++)
			for (Dest = 1; Dest <= no_zones; Dest++)
			{
				total_o_table[Orig] += ODtable[m][Orig][Dest];
				sum += ODtable[m][Orig][Dest];
			}


	return (sum);
}

int Read_ODflow(double* TotalODflow, int* number_of_modes, int* no_zones)
{
	FILE* ODflowFile;
	double RealTotal, InputTotal;

	MDODflow = (double***)Alloc_3D(*number_of_modes, *no_zones, *no_zones, sizeof(double));
	seed_MDODflow = (double***)Alloc_3D(*number_of_modes, *no_zones, *no_zones, sizeof(double));

	if (g_ODME_mode == 1)
	{
		old_MDODflow = (double***)Alloc_3D(*number_of_modes, *no_zones, *no_zones, sizeof(double));
		candidate_MDODflow = (double***)Alloc_3D(*number_of_modes, *no_zones, *no_zones, sizeof(double));
		gradient_MDODflow = (double***)Alloc_3D(*number_of_modes, *no_zones, *no_zones, sizeof(double));
	}

	targetMDODflow = (double***)Alloc_3D(*number_of_modes, *no_zones, *no_zones, sizeof(double));
		 
	TotalOFlow = (double*)Alloc_1D(*no_zones, sizeof(double));

	MDDiffODflow = (double***)Alloc_3D(*number_of_modes, *no_zones, *no_zones, sizeof(double));

	MDRouteCost = (double***)Alloc_3D(*number_of_modes, *no_zones, *no_zones, sizeof(double));

	int with_basedemand = Read_ODtable(MDODflow, MDDiffODflow, seed_MDODflow, targetMDODflow, *no_zones);

	RealTotal = (double)Sum_ODtable(MDODflow, TotalOFlow, *no_zones);


	zone_outbound_link_size = (int*)Alloc_1D(*no_zones, sizeof(int));

	for (int n = 1; n <= *no_zones; n++)
	{
		zone_outbound_link_size[n] = 0;
	}
	for (int k = 1; k <= number_of_links; k++)
	{
		// centroids occupy node seq 1..Z after the internal renumbering, so the
		// INTERNAL from-node seq doubles as the dense zone seq.
		int from_node_seq = Link[k].internal_from_node_id;
		if (from_node_seq >= 1 && from_node_seq <= *no_zones)
		{
			zone_outbound_link_size[from_node_seq] += 1;
		}
	}

	

	// checking 
	int total_infeasible_outbound_zones = 0; 
	float total_infeasible_outbound_zone_demand = 0;
	float total_zone_demand = 0;
	for (int z = 1; z <= *no_zones; z++)
	{
		if (zone_outbound_link_size[z] == 0 && TotalOFlow[z]>0.01)
		{
			printf("Error: There is no outbound link from zone %d with positive demand %f\n",
				g_zone_seq_2_old_zone_id[z], TotalOFlow[z]);
			total_infeasible_outbound_zones++; 
			total_infeasible_outbound_zone_demand += TotalOFlow[z];
		}
		total_zone_demand += TotalOFlow[z];
	}

	if(total_infeasible_outbound_zones)
	{
	printf("Error: %d zones have no outbound link with positive demand %f, = %f percentage of total demand\n", total_infeasible_outbound_zones,
		total_infeasible_outbound_zone_demand, total_infeasible_outbound_zone_demand*100/fmax(0.01, total_zone_demand));
	}

	*TotalODflow = (double)RealTotal;

	return with_basedemand;

}

void CloseLinks(void)
{
	int Node;
	free(zone_outbound_link_size);
	delete[] Link;   // matches new link_record[] allocation in InitLinks()
	free(FirstLinkFrom);
	free(LastLinkFrom);

}

double LinksSDLineSearch(double* MainVolume, double* SDVolume) {
	int n;
	double lambdaleft = 0, lambdaright = 1, lambda = 0.5;
	double grad;

	// Initial check at lambda = 0
	grad = OF_LinksDirectionalDerivative(MainVolume, SDVolume, 0.0);
	StatusMessage("Line search step", "0.0");
	StatusMessage("Line search grad", "%lf", grad);
	if (grad >= 0) {
		LastLambda = 0.0;
		return 0.0;
	}

	// Check at lambda = 1
	grad = OF_LinksDirectionalDerivative(MainVolume, SDVolume, 1.0);
	StatusMessage("Line search step", "1.0");
	StatusMessage("Line search grad", "%lf", grad);
	if (grad <= 0) {
		LastLambda = 1.0;
		return 1.0;
	}

	// Bisection method for line search within [0, 1]
	for (n = 1; n <= MinLineSearchIterations; n++) {
		grad = OF_LinksDirectionalDerivative(MainVolume, SDVolume, lambda);

		if (grad <= 0.0) {
			lambdaleft = lambda;
		}
		else {
			lambdaright = lambda;
		}

		lambda = 0.5 * (lambdaleft + lambdaright);
	}

	// Additional iterations to ensure convergence, if necessary
	for (; lambdaleft == 0 && n <= MAX_NO_BISECTITERATION; n++) {
		grad = OF_LinksDirectionalDerivative(MainVolume, SDVolume, lambda);

		if (grad <= 0.0) {
			lambdaleft = lambda;
		}
		else {
			lambdaright = lambda;
		}

		lambda = 0.5 * (lambdaleft + lambdaright);
	}

	ActualIterations = n - 1;
	LastLambda = lambdaleft;
	return lambdaleft;
}

/// <summary>
///  below is simulation codes
/// </summary>
/// 
/// 
/// 
/// 
// Add at the beginning of your code
enum LogLevel {
	LOG_ERROR = 0,    // Critical errors
	LOG_WARN = 1,     // Warnings
	LOG_INFO = 2,     // Important simulation events
	LOG_DEBUG = 3,    // Detailed debugging info
	LOG_TRACE = 4     // Most detailed level
};

class SimulationLogger {
private:
	ofstream debug_file;
	ofstream info_file;
	ofstream error_file;
	LogLevel current_level;

	static const char* getLevelString(LogLevel level) {
		switch (level) {
		case LOG_ERROR: return "ERROR";
		case LOG_WARN:  return "WARN";
		case LOG_INFO:  return "INFO";
		case LOG_DEBUG: return "DEBUG";
		case LOG_TRACE: return "TRACE";
		default: return "UNKNOWN";
		}
	}

public:
	SimulationLogger(LogLevel level = LOG_INFO) : current_level(level) {
		debug_file.open("sim_debug.csv");
		info_file.open("sim_info.csv");
		error_file.open("sim_error.log");

		// Write headers
		debug_file << "timestamp,level,event,agent_id,link_id,node_id,details\n";
		info_file << "timestamp,level,event,count,details\n";
	}

	void log(LogLevel level, const string& event, int agent_id, int link_id, int node_id, const string& details) {
		if (level > current_level) return;

		string timestamp = to_string(time(nullptr));

		// Always log errors to error file
		if (level == LOG_ERROR) {
			error_file << timestamp << " [ERROR] " << event << " Agent:" << agent_id
				<< " Link:" << link_id << " Node:" << node_id << " - " << details << endl;
		}

		// Detailed debugging info
		if (level <= LOG_DEBUG) {
			debug_file << timestamp << "," << getLevelString(level) << "," << event << ","
				<< agent_id << "," << link_id << "," << node_id << "," << details << endl;
		}

		// Summary information
		if (level <= LOG_INFO) {
			info_file << timestamp << "," << getLevelString(level) << "," << event << ","
				<< "," << details << endl;
		}
	}

	void logAgentMovement(int t, int agent_id, int from_link, int to_link, const string& status) {
		if (LOG_TRACE <= current_level) {
			debug_file << t << ",TRACE,MOVE," << agent_id << "," << from_link << "->"
				<< to_link << "," << status << endl;
		}
	}

	void logQueueStatus(int t, int link_id, int entrance_size, int exit_size) {
		if (LOG_DEBUG <= current_level) {
			debug_file << t << ",DEBUG,QUEUE," << link_id << "," << entrance_size
				<< "," << exit_size << endl;
		}
	}

	void logSimulationStatus(int t, int loaded, int transfers, int completed) {
		if (LOG_INFO <= current_level) {
			info_file << t << ",INFO,STATUS," << loaded << "," << transfers
				<< "," << completed << endl;
		}
	}

	~SimulationLogger() {
		debug_file.close();
		info_file.close();
		error_file.close();
	}
};

//// Add to ParallelQueueSimulator class:
//SimulationLogger logger;
//
//// Then in your simulation methods, add logging:
//logger.logEvent(t, "agent_load", agent->agent_uid, l,
//	link.entrance_queue.size(), "Agent loaded to link");
class CAgent_Simu {
public:
	int agent_uid;
	int o_zone_id, d_zone_id;
	double distannce; 
	double travel_time; 
	int agent_seq_no; 
	double departure_time_in_min;
	int departure_time_in_simu_interval; 
	vector<int> path_link_sequence;
	int current_link_seq_no; 
	bool loaded;
	double desired_free_travel_time_ratio;

	vector<double> link_arrival_times;  // in simu time interval
	vector<double> link_departure_times; // in simu time interval

	CAgent_Simu() :
		agent_uid(-1),
		agent_seq_no (-1),
		distannce (0),
		travel_time(0),
		departure_time_in_min(0),
		departure_time_in_simu_interval (0),
		current_link_seq_no(0),
		loaded(false),
		desired_free_travel_time_ratio(1.0) {}

	void initializeTimes(int path_size) {
		link_arrival_times.resize(path_size);
		link_departure_times.resize(path_size);
	}
};

int g_ParserIntSequence(std::string str, std::vector<int>& vect)
{
	std::stringstream ss(str);
	int i;
	while (ss >> i)
	{
		vect.push_back(i);
		if (ss.peek() == ';')
			ss.ignore();
	}

	return vect.size();
}
int read_vehicle_file(vector<shared_ptr<CAgent_Simu>>& agents) {
	CDTACSVParser parser_vehicle;
	int count = 0;

	if (!parser_vehicle.OpenCSVFile("vehicle.csv", true)) {
		cerr << "Error: Cannot open vehicle.csv" << endl;
		return 0;
	}

	while (parser_vehicle.ReadRecord()) {
		int agent_id;
		double departure_time;  // in min 
		string link_ids;

		if (!parser_vehicle.GetValueByFieldName("agent_id", agent_id) ||
			!parser_vehicle.GetValueByFieldName("departure_time", departure_time) ||
			!parser_vehicle.GetValueByFieldName("link_ids", link_ids)) {
			continue;
		}

		auto agent = make_shared<CAgent_Simu>();
		agent->agent_uid = agent_id;
		agent->departure_time_in_min = departure_time;
	 
		vector<int> external_path_sequence;
		if (g_ParserIntSequence(link_ids, external_path_sequence) > 0) {
			vector<int> internal_path_sequence;
			internal_path_sequence.reserve(external_path_sequence.size());
			bool valid_path = true;
			for (size_t i = 0; i < external_path_sequence.size(); i++)
			{
				int internal_link_id = LinkSeqFromExternal(external_path_sequence[i]);
				if (internal_link_id == 0)
				{
					printf("WARNING: vehicle %d references unknown external link_id %d; skipped.\n",
						agent_id, external_path_sequence[i]);
					valid_path = false;
					break;
				}
				internal_path_sequence.push_back(internal_link_id);
			}
			if (!valid_path)
				continue;
			agent->path_link_sequence = internal_path_sequence;
			agent->initializeTimes(internal_path_sequence.size());
			agent->agent_seq_no = agents.size();
			parser_vehicle.GetValueByFieldName("o_zone_id", agent->o_zone_id);
			parser_vehicle.GetValueByFieldName("d_zone_id", agent->d_zone_id);
			parser_vehicle.GetValueByFieldName("total_distance_mile", agent->distannce);
			agents.push_back(agent);
			count++;
		}
	}

	// Sorting agents by departure time in ascending order
	sort(agents.begin(), agents.end(), [](const shared_ptr<CAgent_Simu>& a, const shared_ptr<CAgent_Simu>& b) {
		return a->departure_time_in_min < b->departure_time_in_min;
		});



	parser_vehicle.CloseCSVFile();

	//// Print the number of vehicles loaded
	//cout << "Number of vehicles loaded: " << count << endl;
	return count;
}

// Global state management
struct SimulationState {
	// Per-second state variables
	vector<vector<double>> linkOutFlowCapacity;    // State per second
	vector<vector<double>> linkOutFlowState;       // State per second

	// Per-minute state variables 
	vector<vector<double>> linkCumulativeArrival;  // Aggregated per minute
	vector<vector<double>> linkCumulativeDeparture;// Aggregated per minute
	vector<vector<double>> linkWaitingTime;        // Tracked per minute

	void initialize(int numLinks, int demandPeriodStartHours, int demandPeriodEndHours) {
		const int BUFFER_HOURS = 2;  // 120 minutes buffer
		int intervalsPerMin = (demandPeriodEndHours - demandPeriodStartHours + BUFFER_HOURS) * 60;
		int intervalsPerSec = intervalsPerMin * 60;

		// Initialize per-second arrays
		linkOutFlowCapacity.resize(numLinks, vector<double>(intervalsPerSec));
		linkOutFlowState.resize(numLinks, vector<double>(intervalsPerSec));

		// Initialize per-minute arrays
		linkCumulativeArrival.resize(numLinks, vector<double>(intervalsPerMin));
		linkCumulativeDeparture.resize(numLinks, vector<double>(intervalsPerMin));
		linkWaitingTime.resize(numLinks, vector<double>(intervalsPerMin));
	}

	void cleanup() {
		linkOutFlowCapacity.clear();
		linkOutFlowState.clear();
		linkCumulativeArrival.clear();
		linkCumulativeDeparture.clear();
		linkWaitingTime.clear();
	}
};

class TrafficSimulator {
private:
	SimulationState state;
	int numLinks;
	int numIntervals;
	const double SECONDS_PER_INTERVAL = 0.1;  // 6 seconds

	void initializeMemory() {
		state.initialize(numLinks, demand_period_starting_hours, demand_period_ending_hours);
	}

	void processSignalizedIntersections() {
		for (int link = 0; link < numLinks; ++link) {
			if (!Link[link].timing_arc_flag) continue;

			int cycleLength = Link[link].cycle_length;
			int startGreen = Link[link].start_green_time;
			int endGreen = Link[link].end_green_time;

			if (endGreen < startGreen) {
				endGreen += cycleLength;
			}

			processSignalTimings(link, cycleLength, startGreen, endGreen);
		}
	}

	void processSignalTimings(int link, int cycleLength, int startGreen, int endGreen) {
		int numCycles = (demand_period_ending_hours - demand_period_starting_hours) * 60 / fmax(1.0f, cycleLength);

		for (int cycle = 0; cycle < numCycles; ++cycle) {
			for (int t = cycle * cycleLength + startGreen;
				t <= cycle * cycleLength + endGreen; t++) {

				float saturationFlow = Link[link].Lane_SaturationFlowRate *
					Link[link].lanes / 3600.0;

				state.linkOutFlowCapacity[link][t] = (t % 2 == 0) ?
					Link[link].lanes : 0;
				state.linkOutFlowState[link][t] = 1;
			}
		}
	}

public:
	TrafficSimulator(int links, int intervals) : numLinks(links), numIntervals(intervals) {
		initializeMemory();
	}

	~TrafficSimulator() {
		state.cleanup();
	}

};
int number_of_seconds_per_interval = 6; 
class LinkQueue {
public:
	deque<int> loading_entrance_queue;       // Queue of vehicle IDs waiting to enter at loading stage
	deque<int> entrance_queue;       // Queue of vehicle IDs waiting to enter
	deque<int> exit_queue;          // Queue of vehicle IDs waiting to exit
	double capacity_per_time_step;   // Maximum number of vehicles that can enter/exit per time step
	double free_flow_time;          // Free flow travel time in minutes
	int number_of_lanes;            // Number of lanes
	double length;                  // Link length
	double spatial_capacity;        // Maximum number of vehicles that can be on the link

	LinkQueue() :
		capacity_per_time_step(0),
		free_flow_time(0),
		number_of_lanes(1),
		length(0),
		spatial_capacity(0) {}

	void reset() {
		entrance_queue.clear();
		exit_queue.clear();
	}

	void setup(link_record link_record_element)
	{
		free_flow_time = link_record_element.FreeTravelTime;
		number_of_lanes = link_record_element.lanes;
		length = link_record_element.length;   // internal unit: miles
		capacity_per_time_step = link_record_element.Link_Capacity / 3600;
		double kjam = 300; 
		spatial_capacity = link_record_element.length * kjam;
	}
};

class AgentLogger {
private:
	ofstream agent_file;

	string convertTimeToHHMMSS(double time_in_min) {
		int hours = (int)(time_in_min / 60);
		int minutes = (int)(time_in_min) % 60;
		int seconds = (int)((time_in_min - (int)time_in_min) * 60);

		char buffer[9];
		sprintf(buffer, "%02d:%02d:%02d", hours, minutes, seconds);
		return string(buffer);
	}

	string linkPathToExternalString(const vector<int>& path, const string& delimiter = ";") {
		stringstream ss;
		for (size_t i = 0; i < path.size(); ++i) {
			if (i > 0) ss << delimiter;
			int internal_link_id = path[i];
			if (internal_link_id >= 1 && internal_link_id <= number_of_links)
				ss << Link[internal_link_id].external_link_id;
		}
		return ss.str();
	}

	string formatTimesVector(const vector<double>& times, int l_demand_period_starting_hours, const string& delimiter = ";") {
		stringstream ss;
		for (size_t i = 0; i < times.size(); ++i) {
			if (i > 0) ss << delimiter;
			if (times[i] >= 0) {
				ss << convertTimeToHHMMSS(l_demand_period_starting_hours*60.0+ times[i] * number_of_seconds_per_interval / 60.0);
			}
			else {
				ss << "NA";
			}
		}
		return ss.str();
	}

public:
	AgentLogger() {
		agent_file.open("trajectory.csv");
		// Write header
		agent_file << "agent_id,departure_time,departure_time_hhmmss,loaded_status,o_zone_id,d_zone_id,distance,travel_time,"
			<< "current_link_seq_no,link_ids,arrival_times,departure_times,geometry\n";
	}


	string getTrajectoryGeometry(const vector<int>& path_link_sequence) {
		stringstream ss;

		ss << fixed << setprecision(7);  // Fixed notation with 7 decimal places
		ss << "\"LINESTRING (";

		bool first_point = true;
		for (int link_id : path_link_sequence) {
			// Get from node of link
			if (!first_point) {
				ss << ", ";
			}

			int from_node = Link[link_id].internal_from_node_id;
			ss << g_node_vector[from_node].x << " " << g_node_vector[from_node].y;

			// Add to node only for the last link
			if (link_id == path_link_sequence.back()) {
				int to_node = Link[link_id].internal_to_node_id;
				ss << ", " << g_node_vector[to_node].x << " " << g_node_vector[to_node].y;
			}

			first_point = false;
		}
		ss << ")\"";
		return ss.str();
	}
	void logAgent(const shared_ptr<CAgent_Simu>& agent, int l_demand_period_starting_hours) {



		double travel_time = (agent->link_departure_times[agent->current_link_seq_no] -
			agent->link_arrival_times[0]) * number_of_seconds_per_interval / 60.0; 

		agent_file << agent->agent_uid << ","
			<< agent->departure_time_in_min << ","
			<< convertTimeToHHMMSS(agent->departure_time_in_min) << ","
			<< (agent->loaded ? "1" : "0") << ","
			<< agent->o_zone_id << "," 
			<< agent->d_zone_id << ","
			<< agent->distannce << ","
			<< agent->current_link_seq_no << ","
			<< linkPathToExternalString(agent->path_link_sequence) << ","
			<< formatTimesVector(agent->link_arrival_times, l_demand_period_starting_hours) << ","
			<< formatTimesVector(agent->link_departure_times, l_demand_period_starting_hours) << ","
			<< getTrajectoryGeometry(agent->path_link_sequence) << "\n";
	}

	void logAllAgents(const vector<shared_ptr<CAgent_Simu>>& agents, int l_demand_period_starting_hours) {
		for (const auto& agent : agents) {
			logAgent(agent,l_demand_period_starting_hours);
		}
		agent_file.flush();
	}

	~AgentLogger() {
		agent_file.close();
	}
};

class ParallelQueueSimulator {
public: 
	SimulationLogger logger;
	// Add debugging counters
	int total_loaded_agents = 0;
	int total_link_transfers = 0;
	int total_node_transfers = 0;
	int total_completed_trips = 0;
	ofstream debug_file;


	SimulationState state;
	vector<LinkQueue> link_queues;
	vector<shared_ptr<CAgent_Simu>> agents;
	int num_links;
	int num_nodes;

	void setupLink_queue()
	{
		for (int l = 1; l <= num_links; l++) {
					link_queues[l].setup(Link[l]);
					}
	}


	double calculateOutflowCapacity(int link_id, int time_step, int total_intervals) {
		const LinkQueue& link = link_queues[link_id];
		double base_capacity;

		// 1. Basic capacity calculation based on lanes and flow rate
		float discharge_rate_per_sec = link.capacity_per_time_step;
		base_capacity = discharge_rate_per_sec * number_of_seconds_per_interval;
		int BUFFER_HOURS = 2;  // 120 minutes buffer
		// 4. Special conditions handling
		if (total_intervals - time_step <= BUFFER_HOURS * 60 * number_of_seconds_per_interval) {
			// After loading period, increase capacity to clear network
			base_capacity *= 10;  // Multiplier for quick network clearance
		}
		//// 2. Signal timing handling
		//if (g_link_vector[link_id].timing_arc_flag) {
		//	int cycle_length = g_link_vector[link_id].cycle_length;
		//	int start_green = g_link_vector[link_id].start_green_time;
		//	int end_green = g_link_vector[link_id].end_green_time;

		//	if (end_green < start_green) {
		//		end_green += cycle_length;  // Handle wraparound
		//	}

		//	int time_in_cycle = time_step % cycle_length;
		//	if (time_in_cycle >= start_green && time_in_cycle <= end_green) {
		//		// Green phase - use saturation flow rate
		//		float saturation_flow_rate = g_link_vector[link_id].saturation_flow_rate *
		//			link.number_of_lanes / 3600.0;
		//		base_capacity = saturation_flow_rate * number_of_seconds_per_interval;
		//	}
		//	else {
		//		// Red phase
		//		return 0;
		//	}
		//}

		//// 3. Capacity reduction handling (incidents, work zones, etc.)
		//if (!g_link_vector[link_id].capacity_reduction_map.empty()) {
		//	int current_time_in_min = time_step * number_of_seconds_per_interval / 60;
		//	auto it = g_link_vector[link_id].capacity_reduction_map.find(current_time_in_min);

		//	if (it != g_link_vector[link_id].capacity_reduction_map.end()) {
		//		base_capacity *= it->second;  // Apply reduction factor
		//	}
		//}



		//// 5. Queue spillback check
		//if (link.exit_queue.size() >= link.spatial_capacity) {
		//	// If link is at spatial capacity, reduce outflow
		//	base_capacity = 0;
		//}

		// 6. Randomization for more realistic flow
		unsigned int RandomSeed = 101 + time_step;  // Simple seed based on time
		float residual = base_capacity - floor(base_capacity);
		RandomSeed = (LCG_a * RandomSeed + LCG_c) % LCG_M;  // Linear Congruential Generator
		float random_ratio = float(RandomSeed) / LCG_M;

		if (random_ratio < residual) {
			base_capacity = ceil(base_capacity);
		}
		else {
			base_capacity = floor(base_capacity);
		}

		// 7. Final bounds check
		return max(0.0, base_capacity);  // Maximum 2 vehicles per lane per interval
	}
	void estimateLinkCapacities(int t,  int total_intervals) {

		for (int l = 1; l <= num_links; l++) {
			double outflow_capacity = calculateOutflowCapacity(l, t , total_intervals);
			state.linkOutFlowCapacity[l][t] = outflow_capacity;
		}
	}

	void synchronizeBottlenecksFromNodetoLink(int t) {
#pragma omp parallel for
		for (int i = 1; i <= num_nodes; i++) {
			auto& node = g_node_vector[i];
			double total_demand = 0;

			// Calculate total demand
			for (int link_id : node.m_incoming_link_seq_no_vector) {
				total_demand += link_queues[link_id].exit_queue.size();
			}

			//if (total_demand >= 1)
			//{
			//	// Distribute capacity if bottleneck exists
			//	for (int link_id : node.m_incoming_link_seq_no_vector) {
			//		if (total_demand > state.linkOutFlowCapacity[link_id][t]) {
			//			double capacity_ratio = state.linkOutFlowCapacity[link_id][t] / total_demand;
			//			state.linkOutFlowCapacity[link_id][t] *= capacity_ratio;
			//		}
			//	}
			//}
		}
	}

	void logSimulationStatus(int t) {
		int t_in_min = t * number_of_seconds_per_interval / 60;
		if (t % 600 == 0) {  // Log every minute
			fprintf(summary_log_file, "\n=== Simulation Time: %d minutes ===\n", t_in_min);
			fprintf(summary_log_file, "Total agents loaded: %d\n", total_loaded_agents);
			fprintf(summary_log_file, "Total link transfers: %d\n", total_link_transfers);
			fprintf(summary_log_file, "Total node transfers: %d\n", total_node_transfers);
			fprintf(summary_log_file, "Total completed trips: %d\n", total_completed_trips);
			debug_file << t << ",status," << total_loaded_agents << ","
				<< total_link_transfers << "," << total_node_transfers
				<< ",minute_summary\n";

			printf("\n=== Simulation Time: %d minutes ===\n", t_in_min);
			printf( "Total agents loaded: %d\n", total_loaded_agents);
			printf( "Total link transfers: %d\n", total_link_transfers);
			printf("Total node transfers: %d\n", total_node_transfers);
			printf( "Total completed trips: %d\n", total_completed_trips);



		}
	}

	void loadNewAgentsPerSimuInterval(int t) {
		int agents_loaded_this_step = 0;
		for (int l = 1; l <= num_links; l++) {
			LinkQueue& link = link_queues[l];
			while (!link.loading_entrance_queue.empty()) {
				int agent_seq_no = link.loading_entrance_queue.front();
				auto agent = agents[agent_seq_no];

				if (agent->departure_time_in_simu_interval == t) {
					link.loading_entrance_queue.pop_front();
					link.entrance_queue.push_back(agent->agent_seq_no);
					agent->link_arrival_times[agent->current_link_seq_no] = t;

					total_loaded_agents++;
					agents_loaded_this_step++;


				}
				else {
					break;
				}
			}
		}

		
	}

	void processNodeTransfers(int t) {
		int t_in_min = t * number_of_seconds_per_interval / 60;

#pragma omp parallel for
		for (int i = 1; i <= num_nodes; i++) {
			auto& node = g_node_vector[i];
			int node_transfers = 0;

			for (int link_id : node.m_incoming_link_seq_no_vector) {
				auto& link = link_queues[link_id];

				while (!link.exit_queue.empty()) {
					int agent_seq_no = link.exit_queue.front();
					auto agent = agents[agent_seq_no];

					if (agent->link_departure_times[agent->current_link_seq_no] > t) {
						break;
					}

					if (agent->current_link_seq_no >= agent->path_link_sequence.size() - 1) {
						// Agent completed trip
						link.exit_queue.pop_front();
#pragma omp atomic
						total_completed_trips++;
						debug_file << t << ",complete," << agent->agent_uid << ","
							<< link_id << ",0,trip_completed\n";
						continue;
					}

					int next_link_id = agent->path_link_sequence[agent->current_link_seq_no + 1];
					auto& next_link = link_queues[next_link_id];

					if (state.linkOutFlowCapacity[link_id][t] >= 1 &&
						next_link.entrance_queue.size() < next_link.capacity_per_time_step) {

						link.exit_queue.pop_front();
						next_link.entrance_queue.push_back(agent_seq_no);

						agent->current_link_seq_no++;
						agent->link_arrival_times[agent->current_link_seq_no] = t;

						node_transfers++;
#pragma omp atomic
						total_node_transfers++;

						debug_file << t << ",transfer," << agent->agent_uid << ","
							<< link_id << "->" << next_link_id << ","
							<< next_link.entrance_queue.size() << ",node_transfer\n";
					}
					else {
						// Log capacity constraint
						debug_file << t << ",blocked," << agent->agent_uid << ","
							<< link_id << "," << next_link.entrance_queue.size()
							<< ",capacity_constraint\n";
						break;
					}
				}
			}

			if (node_transfers > 0) {
#pragma omp critical
				{
					fprintf(summary_log_file, "Time %d: Node %d processed %d transfers\n",
						t, i, node_transfers);
				}
			}
		}
	}


	void processLinkTransfer(int t) {

		int t_in_min = t * number_of_seconds_per_interval / 60; 
#pragma omp parallel for
		for (int l = 1; l <= num_links; l++) {
			auto& link = link_queues[l];

			while (!link.entrance_queue.empty()) {
				int agent_seq_no = link.entrance_queue.front();
				auto agent = agents[agent_seq_no];

				double arrival_time_in_simu_interval = t;
				double free_flow_time = link.free_flow_time * agent->desired_free_travel_time_ratio;
				agent->link_departure_times[agent->current_link_seq_no] =
					t + (int)(free_flow_time * 60 / number_of_seconds_per_interval);

				link.entrance_queue.pop_front();
				link.exit_queue.push_back(agent_seq_no);
				//state.linkCumulativeArrival[l][t_in_min]++;
			}
		}
	}

public:
	ParallelQueueSimulator(int links_count, int nodes_count)
		: num_links(links_count), num_nodes(nodes_count) {
		link_queues.resize(num_links + 1);
		state.initialize(num_links + 1, demand_period_starting_hours, demand_period_ending_hours);

		// Initialize debug file
		debug_file.open("simulation_debug.csv");
		debug_file << "time,event_type,agent_id,link_id,queue_size,details\n";
	}

	~ParallelQueueSimulator() {
		debug_file.close();
	}

	void simulate(int total_intervals, int l_demand_period_starting_hours) {
		AgentLogger agent_logger;
		for (int t = 0; t < total_intervals; t++) 
		{
			loadNewAgentsPerSimuInterval(t);
			estimateLinkCapacities(t, total_intervals);
			synchronizeBottlenecksFromNodetoLink(t);
			processNodeTransfers(t);
			processLinkTransfer(t);
			logSimulationStatus(t);

			if (t>=600 && total_completed_trips > 0 && total_completed_trips == total_loaded_agents)
				break;
		}

		agent_logger.logAllAgents(agents, l_demand_period_starting_hours);
	}
	void loadAgents(const vector<shared_ptr<CAgent_Simu>>& agent_list, int demand_starting_time_in_simu_inteval) {
		fprintf(summary_log_file, "\n=== Loading Agents ===\n");
		fprintf(summary_log_file, "Total agents to load: %zu\n", agent_list.size());
		fprintf(summary_log_file, "Demand starting time in simulation interval: %d\n", demand_starting_time_in_simu_inteval);

		agents = agent_list;
		int agents_with_empty_path = 0;
		map<int, int> link_loading_counts;  // Track how many agents are assigned to each first link

		for (auto& agent : agents) {
			agent->loaded = false;
			agent->current_link_seq_no = 0;

			// Log agent initial setup
			debug_file << "0,init," << agent->agent_uid << ",0,0,departure_time="
				<< agent->departure_time_in_min << "\n";

			if (!agent->path_link_sequence.empty()) {
				int first_link = agent->path_link_sequence[0];
				agent->departure_time_in_simu_interval =
					(int)(agent->departure_time_in_min * 60 / number_of_seconds_per_interval) - demand_starting_time_in_simu_inteval;

				link_queues[first_link].loading_entrance_queue.push_back(agent->agent_seq_no);
				link_loading_counts[first_link]++;

				// Log successful agent assignment
				debug_file << "0,assign," << agent->agent_uid << "," << first_link << ","
					<< link_queues[first_link].loading_entrance_queue.size()
					<< ",departure_interval=" << agent->departure_time_in_simu_interval << "\n";
			}
			else {
				agents_with_empty_path++;
				// Log error for agents with empty paths
				debug_file << "0,error," << agent->agent_uid << ",0,0,empty_path\n";
			}
		}


		//// Log link loading distribution
		//fprintf(summary_log_file, "\nInitial link loading distribution:\n");
		//for (const auto& [link_id, count] : link_loading_counts) {
		//	fprintf(summary_log_file, "Link %d: %d agents\n", link_id, count);
		//	debug_file << "0,link_summary," << link_id << ",0," << count
		//		<< ",initial_loading\n";
		}

}; 
int SimulationAPI()
{

	fopen_s(&summary_log_file, "summary_log_file.txt", "w");

	double* MainVolume, * SubVolume, * SDVolume, Lambda;
	int*** MDMinPathPredLink;

	read_settings_file();
	if (!ValidateProcessorCount())
	{
		if (summary_log_file != NULL)
			fclose(summary_log_file);
		summary_log_file = NULL;
		return 2;
	}
	ConfigureOpenMPRuntime();
	read_mode_type_file();
	fopen_s(&logfile, "TAP_log.csv", "w");  // Open the log file for writing.
	no_nodes = get_number_of_nodes_from_node_file(no_zones, FirstThruNode);
	number_of_links = get_number_of_links_from_link_file();

	printf("# of nodes= %d, # of zones = %d, First Through Node (Seq No) = %d, number of links = %d\n", no_nodes, no_zones,
		FirstThruNode, number_of_links);

	fprintf(summary_log_file, "no_nodes= %d, no_zones = %d, FirstThruNode (seq No) = %d, number_of_links = %d\n", no_nodes, no_zones,
		FirstThruNode, number_of_links);


	double system_wide_travel_time = 0;
	double system_least_travel_time = 0;

	//fprintf(logfile,
	//    "iteration_no,link_id,internal_from_node_id,internal_to_node_id,volume,capacity,voc,"
	//    "fftt,travel_time,delay\n");



	Init(number_of_modes, no_zones);


	InitializeLinkIndices(number_of_modes, no_zones, TotalAssignIterations);
	int total_simu_time_intervals = (demand_period_ending_hours - demand_period_starting_hours) * 3600 / number_of_seconds_per_interval;

	vector<shared_ptr<CAgent_Simu>> agents;
	read_vehicle_file(agents);
	ParallelQueueSimulator simulator(number_of_links, no_nodes);
	simulator.setupLink_queue();
	simulator.loadAgents(agents, demand_period_starting_hours*3600/number_of_seconds_per_interval);
	simulator.simulate(total_simu_time_intervals, demand_period_starting_hours);

	fclose(summary_log_file);
	return 0;
}


//// map_matching.cpp
//#include <cmath>
//#include <algorithm>
//#include <vector>
//#include <map>
//#include <string>
//#include <iostream>
//#include <cstdio>
//#include <cstdlib>

namespace MapMatching {

	// Constants
	constexpr double PI = 3.141592653589793;
	constexpr double EARTH_RADIUS = 6371000; // in meters
	constexpr double DEFAULT_GRID_RESOLUTION = 0.01; // adjust as needed

	// ---------------------------------------------------------
	// Geometry Utilities
	// ---------------------------------------------------------
	inline double toRadians(double degrees) {
		return degrees * PI / 180.0;
	}

	// Haversine distance (meters) given two lat/lon pairs.
	double haversineDistance(double lat1, double lon1, double lat2, double lon2) {
		lat1 = toRadians(lat1); lon1 = toRadians(lon1);
		lat2 = toRadians(lat2); lon2 = toRadians(lon2);
		double dlat = lat2 - lat1;
		double dlon = lon2 - lon1;
		double a = std::pow(std::sin(dlat / 2), 2) +
			std::cos(lat1) * std::cos(lat2) * std::pow(std::sin(dlon / 2), 2);
		double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));
		return EARTH_RADIUS * c;
	}

	// A simple point structure.
	struct GDPoint {
		double x, y;
	};

	// Euclidean (L2) distance between two points.
	inline double getEuclideanDistance(const GDPoint* p1, const GDPoint* p2) {
		return std::sqrt(std::pow(p1->x - p2->x, 2) + std::pow(p1->y - p2->y, 2));
	}

	// Compute distance from a point to a line segment [fromPt, toPt].
	double getPointToLineDistance(const GDPoint* pt,
		const GDPoint* fromPt,
		const GDPoint* toPt,
		double unitGridRes,
		bool no_intersection_requirement)
	{
		double lineLength = getEuclideanDistance(fromPt, toPt);
		if (lineLength < 1e-6)
			return getEuclideanDistance(pt, fromPt);
		double U = ((pt->x - toPt->x) * (fromPt->x - toPt->x) +
			(pt->y - toPt->y) * (fromPt->y - toPt->y))
			/ (lineLength * lineLength);
		if (!no_intersection_requirement && (U < 0.0 || U > 1.0))
			return unitGridRes; // fallback value
		GDPoint intersection;
		intersection.x = toPt->x + U * (fromPt->x - toPt->x);
		intersection.y = toPt->y + U * (fromPt->y - toPt->y);
		double d1 = getEuclideanDistance(pt, &intersection);
		double d0 = getEuclideanDistance(pt, fromPt);
		double d2 = getEuclideanDistance(pt, toPt);
		return no_intersection_requirement ? std::min({ d1, d0, d2 }) : d1;
	}

	// ---------------------------------------------------------
	// Data Structures for MMLinks and GPS Points.
	// ---------------------------------------------------------
	// A link (named MMLink) is defined by its two endpoints.
	struct MMLink {
		GDPoint fromPt;
		GDPoint toPt;
		int external_from_node_id;
		int external_to_node_id;


		double link_distance; // precomputed total length
		bool bInsideFlag;     // true if the link is within a grid cell with GPS points
		double likelihood_distance; // computed matching cost
		int likely_trace_no;        // matching GPS trace id
		int hit_count;              // count of “hits” from intersection tests
		MMLink(const GDPoint& f, const GDPoint& t, const int& from_node_id, const int& to_node_id)
			: fromPt(f), toPt(t), external_from_node_id(from_node_id), external_to_node_id(to_node_id),
			link_distance(getEuclideanDistance(&f, &t)),
			bInsideFlag(false),
			likelihood_distance(1e9),
			likely_trace_no(-1),
			hit_count(0) {}
	};

	// A GPS trace point.
	struct GPSPoint {
		GDPoint pt;
		int trace_no;
		int global_time;   // in seconds
		GPSPoint() : trace_no(0), global_time(0) {}
	};

	// ---------------------------------------------------------
	// Grid Cell for Matching.
	// ---------------------------------------------------------
	struct GridCell {
		double x, y;  // cell origin (e.g., lower-left)
		std::vector<int> linkIndices; // indices of MMLinks overlapping this cell
		std::vector<int> OriginZoneNodeIndices; // indices of MMLinks overlapping this cell
		std::vector<int> DestinationZoneNodeIndices; // indices of MMLinks overlapping this cell
		std::vector<GPSPoint> gpsPoints; // GPS points falling in this cell
		bool dwell_flag = false;
		float start_min = 1e9; // earliest GPS time (minutes) in cell
		float end_min = 0;     // latest GPS time (minutes)
		float dwell_time = 0;  // computed dwell time (minutes)
	};

	// ---------------------------------------------------------
	// MatchingGrid: builds a grid for MMLinks, marks visited cells from GPS traces,
	// and processes each cell to compute a matching cost.
	// ---------------------------------------------------------
	class MatchingGrid {


	public:

		int m_FirstThruNode;

		MatchingGrid(double left, double right, double bottom, double top, double gridRes = DEFAULT_GRID_RESOLUTION, int first_through_node_no = 0)
			: left_(left), right_(right), bottom_(bottom), top_(top), gridRes_(gridRes), m_FirstThruNode(first_through_node_no)
		{
			cols_ = static_cast<int>(std::ceil((right_ - left_) / gridRes_)) + 1;
			rows_ = static_cast<int>(std::ceil((top_ - bottom_) / gridRes_)) + 1;
			grid_.resize(rows_, std::vector<GridCell>(cols_));
			for (int r = 0; r < rows_; ++r)
				for (int c = 0; c < cols_; ++c) {
					grid_[r][c].x = left_ + c * gridRes_;
					grid_[r][c].y = bottom_ + r * gridRes_;
				}
		}

		// Add this method to retrieve a grid cell by (x, y)
		const GridCell* getCell(double x, double y) const {
			int c = static_cast<int>((x - left_) / gridRes_);
			int r = static_cast<int>((y - bottom_) / gridRes_);
			if (r >= 0 && r < rows_ && c >= 0 && c < cols_)
				return &grid_[r][c];
			return nullptr;
		}

		// Insert a MMLink into every cell overlapping its bounding box.
		void insertLink(int linkIndex, const GDPoint& fromPt, const GDPoint& toPt, int internal_from_node_id, int internal_to_node_id) {
			double minX = std::min(fromPt.x, toPt.x);
			double maxX = std::max(fromPt.x, toPt.x);
			double minY = std::min(fromPt.y, toPt.y);
			double maxY = std::max(fromPt.y, toPt.y);
			int colStart = getColIndex(minX);
			int colEnd = getColIndex(maxX);
			int rowStart = getRowIndex(minY);
			int rowEnd = getRowIndex(maxY);
			for (int r = rowStart; r <= rowEnd; ++r)
				for (int c = colStart; c <= colEnd; ++c)
					if (isValidCell(r, c))
					{
						grid_[r][c].linkIndices.push_back(linkIndex);

						if(internal_from_node_id< m_FirstThruNode)  // connector, internal_from_node_id is the TAZ centroid 
							grid_[r][c].OriginZoneNodeIndices.push_back(internal_from_node_id);
						if (internal_to_node_id < m_FirstThruNode)  // connector, internal_to_node_id is the TAZ centroid 
							grid_[r][c].DestinationZoneNodeIndices.push_back(internal_to_node_id);
					}
		}

		// Insert a GPS point and update the cell's time bounds.
		void insertGPSPoint(const GPSPoint& gps, const GDPoint& pt) {
			int r = getRowIndex(pt.y);
			int c = getColIndex(pt.x);
			if (isValidCell(r, c)) {
				grid_[r][c].gpsPoints.push_back(gps);
				float timeMin = gps.global_time / 60.0f;
				grid_[r][c].start_min = std::min(grid_[r][c].start_min, timeMin);
				grid_[r][c].end_min = std::max(grid_[r][c].end_min, timeMin);
			}
		}

		// Insert a GPS point and update the cell's time bounds.
		void cleanGPSPoint() {
			for (int r = 0; r < rows_; ++r) {
				for (int c = 0; c < cols_; ++c) {
					grid_[r][c].gpsPoints.clear();
				}
			}
		}
		// Process visited cells:
		// 1) Mark dwell cells (if time span > 120 minutes).
		// 2) For each link in a cell, loop over GPS points and update a matching cost.
		// 3) For links with no matching GPS point, assign the closest one.
		void processCells(const std::vector<MMLink>& links,
			std::vector<double>& linkGeneralCost,
			std::vector<int>& linkMatchingTraceNo,
			double nonHitDistanceRatio = 10.0)
		{
			// First pass: mark dwell cells.
			for (int r = 0; r < rows_; ++r) {
				for (int c = 0; c < cols_; ++c) {
					GridCell& cell = grid_[r][c];
					if (!cell.gpsPoints.empty() && (cell.end_min - cell.start_min > 120)) {
						cell.dwell_flag = true;
						cell.dwell_time = cell.end_min - cell.start_min;
					}
				}
			}
			// Second pass: update link cost from GPS points.
			for (int r = 0; r < rows_; ++r) {
				for (int c = 0; c < cols_; ++c) {
					GridCell& cell = grid_[r][c];
					if (!cell.gpsPoints.empty() && !cell.linkIndices.empty()) {
						// Initialize each link’s cost to the grid resolution.
						for (int linkIdx : cell.linkIndices)
							linkGeneralCost[linkIdx] = gridRes_;
						// For each GPS point…
						for (const auto& gps : cell.gpsPoints) {
							for (int linkIdx : cell.linkIndices) {
								const MMLink& L = links[linkIdx];
								double d = getPointToLineDistance(&gps.pt, &L.fromPt, &L.toPt, 1.0, true);
								double cost = d; // Here, you might combine with other measures.
								if (cost < linkGeneralCost[linkIdx]) {
									linkGeneralCost[linkIdx] = cost;
									linkMatchingTraceNo[linkIdx] = gps.trace_no;
								}
							}
						}
						// Stage 2: for links with no matching GPS point, find the closest one.
						for (int linkIdx : cell.linkIndices) {
							if (linkMatchingTraceNo[linkIdx] < 0) {
								double min_d = 1e9;
								int bestTrace = -1;
								for (const auto& gps : cell.gpsPoints) {
									double d = getPointToLineDistance(&gps.pt, &links[linkIdx].fromPt, &links[linkIdx].toPt, 1.0, false);
									if (d < min_d) { min_d = d; bestTrace = gps.trace_no; }
								}
								linkMatchingTraceNo[linkIdx] = bestTrace;
							}
						}
					}
				}
			}
		}

		void printSummary() const {
			std::cout << "MatchingGrid: " << rows_ << " rows x " << cols_ << " cols, resolution "
				<< gridRes_ << "\n";
			std::cout << "Bounds: x[" << left_ << ", " << right_ << "], y[" << bottom_ << ", " << top_ << "]\n";
		}

	private:
		int getColIndex(double x) const { return static_cast<int>((x - left_) / gridRes_); }
		int getRowIndex(double y) const { return static_cast<int>((y - bottom_) / gridRes_); }
		bool isValidCell(int r, int c) const { return r >= 0 && r < rows_ && c >= 0 && c < cols_; }
		double left_, right_, bottom_, top_, gridRes_;
		int rows_, cols_;
		std::vector<std::vector<GridCell>> grid_;
	};

	// Define a structure to hold trace record information.
	struct TraceRecord {
		int trace_no = 0;
		std::string trace_id = "";
		std::string allowed_link_type_code = "";
		std::string blocked_link_type_code = "";
		double y_coord = 0.0;
		double x_coord = 0.0;
		int road_order = 0;
		double reference_speed = 0.0;
		std::string agent_id = "";
		int internal_road_order = 0;
		int o_node_id = -1;
		int d_node_id = -1;
	};

	// Here we change the parameter to a nested map:
	// Outer key: agent_id (std::string)
	// Inner key: trace_no (int)
	// Value: vector of GPSPoints that belong to that trace.
	int readGPSTraceFile(std::unordered_map<std::string, std::vector<GPSPoint>>& gpsTraces)

	{
		CDTACSVParser parser_trace;
		std::vector<TraceRecord> traceRecords;

		// Open the CSV file (the second parameter 'true' indicates that the file contains headers)
		if (parser_trace.OpenCSVFile("trace.csv", true))
		{
			// Read each record until end-of-file.
			while (parser_trace.ReadRecord())
			{
				TraceRecord rec;
				// Read the fields by name.
				parser_trace.GetValueByFieldName("trace_no", rec.trace_no);
				parser_trace.GetValueByFieldName("trace_id", rec.trace_id);
				parser_trace.GetValueByFieldName("allowed_link_type_code", rec.allowed_link_type_code);
				parser_trace.GetValueByFieldName("blocked_link_type_code", rec.blocked_link_type_code);
				parser_trace.GetValueByFieldName("y_coord", rec.y_coord);
				parser_trace.GetValueByFieldName("x_coord", rec.x_coord);
				parser_trace.GetValueByFieldName("road_order", rec.road_order);
				parser_trace.GetValueByFieldName("reference_speed", rec.reference_speed);
				parser_trace.GetValueByFieldName("agent_id", rec.agent_id);
				parser_trace.GetValueByFieldName("internal_road_order", rec.internal_road_order);
				parser_trace.GetValueByFieldName("o_node_id", rec.o_node_id);
				parser_trace.GetValueByFieldName("d_node_id", rec.d_node_id);

				traceRecords.push_back(rec);
			}
			parser_trace.CloseCSVFile();
		}
		else
		{
			std::cerr << "Error: Cannot open trace.csv" << std::endl;
			std::exit(1);
		}

		// Group the GPS traces by agent_id and trace_no
		for (const auto& rec : traceRecords)
		{
			GPSPoint gps;
			gps.pt.x = rec.x_coord;
			gps.pt.y = rec.y_coord;
			gps.trace_no = rec.trace_no;
			// If you have a global_time field in your CSV, you can set gps.global_time here.

			// Group by agent_id then trace_no:
			gpsTraces[rec.agent_id].push_back(gps);
		}

		return 0;
	}


	// Helper: Compute Euclidean distance between two points.
	inline double computeEuclideanDistance(const GDPoint& a, const GDPoint& b) {
		return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
	}

	// Given a set of candidate zone node IDs (from a grid cell) and a point,
	// return the candidate whose coordinate is closest.
	// zoneNodeCoords: a mapping from zone node id to its coordinate.
	int determineZoneNode(const std::vector<int>& candidateNodeIds, const GDPoint& avgPt,
		const std::map<int, GDPoint>& zoneNodeCoords)
	{
		double minDist = std::numeric_limits<double>::max();
		int bestNode = -1;
		for (int nodeId : candidateNodeIds) {
			auto it = zoneNodeCoords.find(nodeId);
			if (it != zoneNodeCoords.end()) {
				double d = computeEuclideanDistance(avgPt, it->second);
				if (d < minDist) {
					minDist = d;
					bestNode = nodeId;
				}
			}
		}
		return bestNode;
	}

	// This function uses the first few and last few GPS points of a trace
	// to determine the origin and destination zone node IDs.
	// Parameters:
	// - trace: vector of GPSPoint representing a single trace.
	// - grid: the MatchingGrid that contains the grid cells.
	// - zoneNodeCoords: mapping from zone node id to its GDPoint (x,y) coordinate.
	// - numPointsToAverage: number of points from the start/end to average.
	// - originNodeId & destinationNodeId: outputs.
	void determineOriginDestinationNodes(const std::vector<GPSPoint>& trace,
		const MatchingGrid& grid,
		const std::map<int, GDPoint>& zoneNodeCoords,
		int numPointsToAverage,
		int& originNodeId,
		int& destinationNodeId)
	{
		if (trace.empty()) {
			originNodeId = -1;
			destinationNodeId = -1;
			return;
		}

		// Average first numPointsToAverage GPS points for the origin.
		GDPoint avgOrigin = { 0, 0 };
		int countOrigin = std::min(numPointsToAverage, static_cast<int>(trace.size()));
		for (int i = 0; i < countOrigin; ++i) {
			avgOrigin.x += trace[i].pt.x;
			avgOrigin.y += trace[i].pt.y;
		}
		avgOrigin.x /= countOrigin;
		avgOrigin.y /= countOrigin;

		// Average last numPointsToAverage GPS points for the destination.
		GDPoint avgDest = { 0, 0 };
		int countDest = std::min(numPointsToAverage, static_cast<int>(trace.size()));
		for (int i = static_cast<int>(trace.size()) - countDest; i < static_cast<int>(trace.size()); ++i) {
			avgDest.x += trace[i].pt.x;
			avgDest.y += trace[i].pt.y;
		}
		avgDest.x /= countDest;
		avgDest.y /= countDest;

		// Retrieve the grid cell for the averaged origin and destination.
		const GridCell* originCell = grid.getCell(avgOrigin.x, avgOrigin.y);
		const GridCell* destCell = grid.getCell(avgDest.x, avgDest.y);

		// Determine the best matching zone node from the candidates in each cell.
		if (originCell && !originCell->OriginZoneNodeIndices.empty())
			originNodeId = determineZoneNode(originCell->OriginZoneNodeIndices, avgOrigin, zoneNodeCoords);
		else
			originNodeId = -1;

		if (destCell && !destCell->DestinationZoneNodeIndices.empty())
			destinationNodeId = determineZoneNode(destCell->DestinationZoneNodeIndices, avgDest, zoneNodeCoords);
		else
			destinationNodeId = -1;
	}


} // end namespace MapMatching



// -------------------------
// -------------------------
// Example main() demonstrating advanced grid-based matching using MMLink
// -------------------------
int mapmatchingAPI() {


	fopen_s(&summary_log_file, "summary_log_file.txt", "w");
	double* MainVolume, * SubVolume, * SDVolume, Lambda;
	int*** MDMinPathPredLink;

	read_settings_file();
	read_mode_type_file();
	fopen_s(&logfile, "TAP_log.csv", "w");  // Open the log file for writing.
	no_nodes = get_number_of_nodes_from_node_file(no_zones, FirstThruNode);
	number_of_links = get_number_of_links_from_link_file();

	printf("# of nodes= %d, # of zones = %d, First Through Node (Seq No) = %d, number of links = %d\n", no_nodes, no_zones,
		FirstThruNode, number_of_links);

	fprintf(summary_log_file, "no_nodes= %d, no_zones = %d, FirstThruNode (seq No) = %d, number_of_links = %d\n", no_nodes, no_zones,
		FirstThruNode, number_of_links);


	double system_wide_travel_time = 0;
	double system_least_travel_time = 0;

	//fprintf(logfile,
	//    "iteration_no,link_id,internal_from_node_id,internal_to_node_id,volume,capacity,voc,"
	//    "fftt,travel_time,delay\n");



	//Init(number_of_modes, no_zones);

	InitLinks();

	MDMinPathPredLink = (int***)Alloc_3D(number_of_modes, no_zones, no_nodes, sizeof(int)); 

	double** CostTo;

	CostTo = (double**)Alloc_2D(no_zones, no_nodes, sizeof(double));
	using namespace MapMatching;

	// Create some MMLinks.
	std::vector<MMLink> mmLinks;

	for (int k = 1; k <= number_of_links; k++)
	{
		double from_x, from_y, to_x, to_y; 
		from_x = g_node_vector[Link[k].internal_from_node_id].x;
		from_y = g_node_vector[Link[k].internal_to_node_id].y;
		to_x = g_node_vector[Link[k].internal_from_node_id].x;
		to_y = g_node_vector[Link[k].internal_to_node_id].y;

		mmLinks.push_back(MMLink({ from_x, from_y }, { to_x, to_y}, Link[k].internal_from_node_id, Link[k].internal_to_node_id));
	}



	// Compute overall boundaries from all MMLink endpoints.
	double left = 1e9, right = -1e9, bottom = 1e9, top = -1e9;
	for (const auto& L : mmLinks) {
		left = std::min({ left, L.fromPt.x, L.toPt.x });
		right = std::max({ right, L.fromPt.x, L.toPt.x });
		bottom = std::min({ bottom, L.fromPt.y, L.toPt.y });
		top = std::max({ top, L.fromPt.y, L.toPt.y });
	}

	// Build the matching grid.
	MatchingGrid mGrid(left, right, bottom, top, DEFAULT_GRID_RESOLUTION, FirstThruNode);
	mGrid.printSummary();
	for (size_t i = 0; i < mmLinks.size(); i++) {
		mGrid.insertLink(static_cast<int>(i), mmLinks[i].fromPt, mmLinks[i].toPt, mmLinks[i].external_from_node_id, mmLinks[i].external_to_node_id);
	}
	std::unordered_map<std::string, std::vector<GPSPoint>> gpsTraces; 
	std::unordered_map<std::string, int > gpsTraces_originNodeId;
	std::unordered_map<std::string, int > gpsTraces_destinationNodeId;
	std::unordered_map<std::string, int > gpsTraces_RouteId;
	readGPSTraceFile(gpsTraces);

	std::map<int, GDPoint> zoneNodeCoords;

	for (int i = 1; i <= no_zones; i++)  // zones-first dense centroid nodes
	{
		GDPoint pt;
		pt.x = g_node_vector[i].x;
		pt.y = g_node_vector[i].y;
		zoneNodeCoords[g_node_vector[i].node_id] = pt;
	}

	std::vector<std::vector<std::vector<int>>> ODRouteIndices(
		no_zones + 1,
		std::vector<std::vector<int>>(no_zones + 1)
	);
	// ------------------------
	// Insert GPS points into the grid.
	// ------------------------
	// Loop over each agent  for that agent.
	int 	max_route_size = 1;

	for (const auto& agentPair : gpsTraces) {

		mGrid.cleanGPSPoint(); 
		const std::string& agentId = agentPair.first;

		std::vector<GPSPoint> gpsPoints;
		gpsPoints = agentPair.second;
		// Optionally, you might want to print or log the agent and trace information:
		// std::cout << "Processing agent " << agentId << ", trace " << traceNo << std::endl;
		for (const auto& gps : gpsPoints) {
			mGrid.insertGPSPoint(gps, gps.pt);
		}

		////#pragma omp parallel for
//		int Orig = agentPair.second.;  // get origin zone id
		int numPointsToAverage = 3;
		int  originNodeId = 0;
		int destinationNodeId = 0;

		determineOriginDestinationNodes(gpsPoints,
			mGrid,
			zoneNodeCoords,
			numPointsToAverage,
			originNodeId,
			destinationNodeId);
		gpsTraces_originNodeId[agentId] = originNodeId;
		gpsTraces_destinationNodeId[agentId] = destinationNodeId;
		int originZoneSeq = NodeSeqFromExternal(originNodeId);
		int destinationZoneSeq = NodeSeqFromExternal(destinationNodeId);
		if (originZoneSeq < 1 || originZoneSeq > no_zones ||
			destinationZoneSeq < 1 || destinationZoneSeq > no_zones)
			continue;
		gpsTraces_RouteId[agentId] =
			(int)ODRouteIndices[originZoneSeq][destinationZoneSeq].size();
		ODRouteIndices[originZoneSeq][destinationZoneSeq].push_back(1);

		if (ODRouteIndices[originZoneSeq][destinationZoneSeq].size() > (size_t)max_route_size)
			max_route_size = (int)ODRouteIndices[originZoneSeq][destinationZoneSeq].size();

	}

	InitializeLinkIndices(number_of_modes, no_zones, max_route_size);;


	for (const auto& agentPair : gpsTraces) {

		mGrid.cleanGPSPoint();
		const std::string& agentId = agentPair.first;

		std::vector<GPSPoint> gpsPoints;
		gpsPoints = agentPair.second;
		// Optionally, you might want to print or log the agent and trace information:
		// std::cout << "Processing agent " << agentId << ", trace " << traceNo << std::endl;
		for (const auto& gps : gpsPoints) {
			mGrid.insertGPSPoint(gps, gps.pt);
		}
		// ------------------------
		// Prepare arrays for each MMLink’s cost and matching trace number.
		// ------------------------
		std::vector<double> linkGeneralCost(mmLinks.size(), DEFAULT_GRID_RESOLUTION);
		std::vector<int>    linkMatchingTraceNo(mmLinks.size(), -1);

		// Process the grid cells.
		mGrid.processCells(mmLinks, linkGeneralCost, linkMatchingTraceNo);

		// Print the matching cost for each MMLink.
		for (size_t i = 0; i < mmLinks.size(); ++i) {
			std::printf("MMLink %zu: cost = %.3f, matched trace = %d\n",
				i, linkGeneralCost[i], linkMatchingTraceNo[i]);
		}


		if (gpsTraces_originNodeId[agentId] >= 1 && gpsTraces_destinationNodeId[agentId] >= 1)
		{
			// external GPS node ids -> internal seq (renumbering boundary)
			int Orig = NodeSeqFromExternal(gpsTraces_originNodeId[agentId]);
			int Dest = NodeSeqFromExternal(gpsTraces_destinationNodeId[agentId]);
			if (Orig < 1 || Orig > no_zones || Dest < 1 || Dest > no_zones)
				continue;
			std::vector<int> currentLinkSequence; // Temporary vector to store link indices


			int m = 1;
			Minpath(m, Orig, MDMinPathPredLink[m][Orig], CostTo[Orig]);
			int CurrentNode = Dest;   // already internal seq

			// MinPathPredLink is coded as internal node id 
			// 
			//double total_travel_time = 0;
			//double total_length = 0;
			//double total_FFTT = 0;
			// back trace 
			while (CurrentNode != Orig)
			{
				int k = MDMinPathPredLink[m][Orig][CurrentNode];

				if (k <= 0 || k > number_of_links || k == INVALID)
				{
					printf("A problem in All_or_Nothing_Assign() Invalid pred for node seq no %d Orig zone= %d \n\n", CurrentNode, Orig);
					break;
				}
				CurrentNode = Link[k].internal_from_node_id;

				if (CurrentNode <= 0 || CurrentNode > no_nodes)
				{
					printf("A problem in All_or_Nothing_Assign() Invalid node seq no %d Orig zone = %d \n\n", CurrentNode, Orig);
					break;
				}

				currentLinkSequence.push_back(k); // Store the link index
			}

			int route_id = gpsTraces_RouteId[agentId];
			AddLinkSequence(m, Orig, Dest, route_id, currentLinkSequence);

		}

	}
				std::vector<double> m_theta; 


				OutputRouteDetails("route_assignment.csv", m_theta);
	
	Free_3D((void***)MDMinPathPredLink, number_of_modes, no_zones, no_nodes); 
	Free_2D((void**)CostTo, no_zones, no_nodes);
	return 0;
}


extern "C" PATH_ENGINE_API int DTA_AssignmentAPIWithStatus()
{
	return AssignmentAPI();
}

int AutoCalibrationAPI()
{
	const bool previous_force = g_force_auto_calibration;
	g_force_auto_calibration = true;
	const int status = AssignmentAPI();
	g_force_auto_calibration = previous_force;
	return status;
}

extern "C" PATH_ENGINE_API int DTA_AutoCalibrationAPIWithStatus()
{
	return AutoCalibrationAPI();
}

extern "C" PATH_ENGINE_API int DTA_AutoCalibrationOneShotAPIWithStatus()
{
	g_auto_calibration_one_shot = true;
	const int status = AutoCalibrationAPI();
	// All CSV/JSON streams are closed by this point. Flush C stdio diagnostics
	// and terminate without destructing hundreds of millions of nested paths.
	fflush(NULL);
	std::_Exit(status);
}

extern "C" PATH_ENGINE_API int DTA_SimulationAPIWithStatus()
{
	return SimulationAPI();
}

extern "C" PATH_ENGINE_API void DTA_AssignmentAPI()
{
	(void)DTA_AssignmentAPIWithStatus();
}

extern "C" PATH_ENGINE_API void DTA_AutoCalibrationAPI()
{
	(void)DTA_AutoCalibrationAPIWithStatus();
}

extern "C" PATH_ENGINE_API void DTA_SimulationAPI()
{
	(void)DTA_SimulationAPIWithStatus();
}

// Optional standalone executable entry point. Built only when BUILD_EXE is
// defined (see -DBUILD_EXE). Lets us run the kernel as a normal .exe in the
// working directory (reads node/link/settings/mode_type/demand from cwd and
// writes link_performance.csv etc.) without deploying a shared library.
#ifdef BUILD_EXE
int main(int argc, char** argv) {
	if (argc >= 2 && std::string(argv[1]) == "--auto-calibration-one-shot")
		return DTA_AutoCalibrationOneShotAPIWithStatus();
	if (argc >= 2 && std::string(argv[1]) == "--auto-calibration")
		return AutoCalibrationAPI();
	return AssignmentAPI();
}
#endif
