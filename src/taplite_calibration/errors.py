"""Package-specific exceptions with actionable user messages."""


class CalibrationError(RuntimeError):
    """Base error raised for a failed calibration workflow."""


class ConfigurationError(CalibrationError, ValueError):
    """The project configuration is incomplete or inconsistent."""


class CapabilityError(CalibrationError):
    """A required optional runtime capability is unavailable."""
