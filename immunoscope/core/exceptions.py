"""
Custom exceptions for ImmunoScope pipeline.

Defines a hierarchy of exceptions for better error handling and debugging.
"""


class ImmunoScopeError(Exception):
    """Base exception for all ImmunoScope errors."""
    pass


class InputValidationError(ImmunoScopeError):
    """
    Input validation error (user error).

    Raised when user-provided input is invalid (e.g., missing files,
    invalid parameters).

    Example:
        >>> raise InputValidationError(
        ...     param_name="topology",
        ...     value="/path/to/missing.tpr",
        ...     reason="file does not exist"
        ... )
    """

    def __init__(self, param_name: str, value: any, reason: str):
        message = f"Invalid {param_name}: {reason} (got: {value})"
        super().__init__(message)
        self.param_name = param_name
        self.value = value
        self.reason = reason


class ProcessingError(ImmunoScopeError):
    """
    Processing error (system error).

    Raised when processing fails due to system issues (e.g., external
    command failure, insufficient resources).

    Example:
        >>> raise ProcessingError(
        ...     step="PBC correction",
        ...     reason="gmx trjconv failed with exit code 1",
        ...     suggestion="check GROMACS installation"
        ... )
    """

    def __init__(self, step: str, reason: str, suggestion: str = None):
        message = f"Failed at {step}: {reason}"
        if suggestion:
            message += f"\nSuggestion: {suggestion}"
        super().__init__(message)
        self.step = step
        self.reason = reason
        self.suggestion = suggestion


class PipelineError(ImmunoScopeError):
    """
    Pipeline execution error.

    Raised when pipeline flow control encounters issues (e.g., missing
    intermediate results, node failure).

    Example:
        >>> raise PipelineError(
        ...     node_name="RMSDNode",
        ...     reason="trajectory_processed not found in context",
        ...     context_state={'system_id': '1ao7', 'has_errors': False}
        ... )
    """

    def __init__(self, node_name: str, reason: str, context_state: dict = None):
        message = f"Pipeline error at {node_name}: {reason}"
        if context_state:
            message += f"\nContext state: {context_state}"
        super().__init__(message)
        self.node_name = node_name
        self.reason = reason
        self.context_state = context_state


class ConfigurationError(ImmunoScopeError):
    """
    Configuration error.

    Raised when configuration is invalid (e.g., invalid YAML, missing
    required fields).

    Example:
        >>> raise ConfigurationError(
        ...     config_file="pipeline.yaml",
        ...     reason="missing required field 'nodes'"
        ... )
    """

    def __init__(self, config_file: str, reason: str):
        message = f"Configuration error in {config_file}: {reason}"
        super().__init__(message)
        self.config_file = config_file
        self.reason = reason


class TaskDiscoveryError(ImmunoScopeError):
    """
    Task discovery error.

    Raised when task discovery encounters a fatal error (e.g., root directory
    does not exist, not accessible).

    Note: Finding invalid/ambiguous tasks is NOT an error - those should be
    recorded in the discovery report.

    Example:
        >>> raise TaskDiscoveryError(
        ...     root="/path/to/nonexistent",
        ...     reason="directory does not exist"
        ... )
    """

    def __init__(self, root: str, reason: str):
        message = f"Task discovery failed for {root}: {reason}"
        super().__init__(message)
        self.root = root
        self.reason = reason


class TaskValidationError(ImmunoScopeError):
    """
    Task validation error.

    Raised when task validation encounters unexpected internal issues.
    This is for internal validation logic errors, not for tasks that
    fail validation (those go into the report).

    Example:
        >>> raise TaskValidationError(
        ...     task_id="system_001",
        ...     reason="unexpected validation state"
        ... )
    """

    def __init__(self, task_id: str, reason: str):
        message = f"Task validation error for {task_id}: {reason}"
        super().__init__(message)
        self.task_id = task_id
        self.reason = reason


class ManifestWriteError(ImmunoScopeError):
    """
    Manifest write error.

    Raised when manifest cannot be written to disk (e.g., directory not
    writable, disk full, serialization errors).

    Example:
        >>> raise ManifestWriteError(
        ...     output_path="/readonly/manifest.jsonl",
        ...     reason="permission denied"
        ... )
    """

    def __init__(self, output_path: str, reason: str):
        message = f"Failed to write manifest to {output_path}: {reason}"
        super().__init__(message)
        self.output_path = output_path
        self.reason = reason


class TaskFileParseError(ImmunoScopeError):
    """
    Task file parse error.

    Raised when task.yaml cannot be parsed (e.g., invalid YAML syntax,
    malformed structure).

    Example:
        >>> raise TaskFileParseError(
        ...     task_file="/path/to/task.yaml",
        ...     reason="invalid YAML syntax"
        ... )
    """

    def __init__(self, task_file: str, reason: str):
        message = f"Failed to parse task file {task_file}: {reason}"
        super().__init__(message)
        self.task_file = task_file
        self.reason = reason
