// pybind11 in-process binding for the TAPLite kernel.
//
// Compiles TAPLite.cpp WITHOUT BUILD_EXE (so its main() is excluded) and exposes the kernel's
// AssignmentAPI() as taplite_calibration._native.run_in_dir(path). The kernel reads/writes CSVs in the
// current working directory, so we chdir into `path` first.
//
// Build: see CMakeLists.txt in this folder (needs pybind11). The resulting _native module is
// packaged directly inside taplite-calibration; no separate TAPLite executable is used.

#include <pybind11/pybind11.h>
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef _WIN32
#include <direct.h>
#define portable_chdir _chdir
#else
#include <unistd.h>
#define portable_chdir chdir
#endif

int AssignmentAPI();    // defined in TAPLite.cpp (compiled here without BUILD_EXE)
int AutoCalibrationAPI();
int ProcessorCountValidationStatus(int requested_processors);

namespace py = pybind11;

static int run_in_dir(const std::string& path) {
    if (!path.empty())
        portable_chdir(path.c_str());
    int rc;
    {
        py::gil_scoped_release release;   // the assignment is a long C++ run; free the GIL
        rc = AssignmentAPI();
    }
    return rc;
}

static int run_auto_calibration_in_dir(const std::string& path) {
    if (!path.empty())
        portable_chdir(path.c_str());
    int rc;
    {
        py::gil_scoped_release release;
        rc = AutoCalibrationAPI();
    }
    return rc;
}

static py::dict openmp_status(int requested_threads) {
    if (requested_threads < 0)
        throw py::value_error("requested_threads must be >= 0");

    bool compiled = false;
    int openmp_version = 0;
    int max_threads = 1;
    int num_procs = 1;
    bool dynamic = false;
    int probe_team_size = 1;

#ifdef _OPENMP
    compiled = true;
    openmp_version = _OPENMP;
    max_threads = omp_get_max_threads();
    num_procs = omp_get_num_procs();
    dynamic = omp_get_dynamic() != 0;

    if (requested_threads > 0) {
#pragma omp parallel num_threads(requested_threads)
        {
#pragma omp single
            probe_team_size = omp_get_num_threads();
        }
    } else {
#pragma omp parallel
        {
#pragma omp single
            probe_team_size = omp_get_num_threads();
        }
    }
#endif

    py::dict status;
    status["compiled"] = compiled;
    status["openmp_version"] = openmp_version;
    status["max_threads"] = max_threads;
    status["num_procs"] = num_procs;
    status["dynamic"] = dynamic;
    status["requested_threads"] = requested_threads;
    status["probe_team_size"] = probe_team_size;
    return status;
}

PYBIND11_MODULE(_native, m) {
    m.doc() = "Bundled TAPLite assignment and refined auto-calibration kernel.";
    m.def("run_in_dir", &run_in_dir, py::arg("path") = "",
          "Run a static assignment, reading CSV inputs from `path` (or the current working\n"
          "directory) and writing link_performance.csv there. Returns the kernel exit code.\n"
          "NOTE: the kernel keeps global state, so run ONE assignment per process — for many\n"
          "runs use an isolated Python worker (taplite_calibration.assign does this by default).");
    m.def("run_auto_calibration_in_dir", &run_auto_calibration_in_dir,
          py::arg("path") = "",
          "Run equilibrium-coupled QVDF auto calibration in one initialized kernel "
          "session. Inner equilibria exchange route and link state in memory and only "
          "the final assignment and calibration outputs are written.");
    m.def("openmp_status", &openmp_status, py::arg("requested_threads") = 0,
          "Return native OpenMP build and runtime diagnostics without running an assignment.");
    m.def("_processor_count_validation_status", &ProcessorCountValidationStatus,
          py::arg("requested_processors"),
          "Return the kernel settings-validation status for a processor count.");
}
