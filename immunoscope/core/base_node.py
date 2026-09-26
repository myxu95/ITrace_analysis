"""
Pipeline Node - Abstract base class for pipeline nodes.

All pipeline nodes (Layer 2) must inherit from this base class.
"""

from abc import ABC, abstractmethod
from typing import Optional
import logging

from .context import PipelineContext
from .exceptions import PipelineError

logger = logging.getLogger(__name__)


class PipelineNode(ABC):
    """
    Abstract base class for pipeline nodes.

    Pipeline nodes are lightweight wrappers around Core Modules (Layer 1).
    They are responsible for:
    - Validating inputs from context
    - Calling the underlying module
    - Recording outputs to context
    - Handling errors gracefully

    Subclasses must implement the execute() method.

    Example:
        >>> class RMSDNode(PipelineNode):
        ...     def execute(self, context: PipelineContext) -> PipelineContext:
        ...         # Validate inputs
        ...         self.validate_inputs(context)
        ...
        ...         # Call underlying module
        ...         calculator = RMSDCalculator()
        ...         result = calculator.calculate(...)
        ...
        ...         # Record outputs
        ...         context.add_result('rmsd', result.to_dict())
        ...
        ...         return context
    """

    def __init__(self, name: Optional[str] = None):
        """
        Initialize pipeline node.

        Args:
            name: Node name (defaults to class name)
        """
        self.name = name or self.__class__.__name__
        self.logger = logging.getLogger(f"immunoscope.pipeline.{self.name}")

    @abstractmethod
    def execute(self, context: PipelineContext) -> PipelineContext:
        """
        Execute node logic.

        This method must be implemented by subclasses.

        Args:
            context: Pipeline context

        Returns:
            Updated pipeline context

        Raises:
            PipelineError: If node execution fails
            Any other exception from underlying modules

        Example:
            >>> def execute(self, context):
            ...     self.logger.info(f"Executing {self.name}")
            ...     # ... do work ...
            ...     return context
        """
        pass

    def validate_inputs(self, context: PipelineContext):
        """
        Validate inputs from context.

        Subclasses can override this method to add custom validation.
        Default implementation does nothing.

        Args:
            context: Pipeline context

        Raises:
            PipelineError: If validation fails

        Example:
            >>> def validate_inputs(self, context):
            ...     if not context.trajectory_processed:
            ...         raise PipelineError(
            ...             node_name=self.name,
            ...             reason="trajectory_processed not found"
            ...         )
        """
        pass

    def __repr__(self) -> str:
        """String representation for debugging."""
        return f"{self.__class__.__name__}(name='{self.name}')"

    def __str__(self) -> str:
        """Human-readable string representation."""
        return self.name
