# Bundled native-kernel example

This four-node, two-route network makes a real call to the bundled TAPLite C++
auto-calibration API. It is separate from `examples/quickstart`, whose smoke
backend tests the complete Python orchestration without performing assignment.

After installing the package from a platform wheel or source checkout:

```text
python examples/native-kernel/run_example.py --output native-example-output
```

The output folder receives the final `link_performance.csv`, calibration
history, link audit, summary, volume-constraint audit, and other kernel logs.
No TAPLite executable or adjacent repository is used.

