"""
Pipeline Base Class - Layer 3 of the architecture.

Pipelines orchestrate the execution of nodes in a specific order.
"""

from typing import List
import logging

from ..core.base_node import PipelineNode
from ..core.context import PipelineContext

logger = logging.getLogger(__name__)


class Pipeline:
    """
    Base class for all pipelines.

    A pipeline is a sequence of nodes that are executed in order.
    It handles:
    - Sequential execution of nodes
    - Error propagation
    - Flow control (should_stop flag)

    Example:
        >>> class MyPipeline(Pipeline):
        ...     def __init__(self):
        ...         self.nodes = [
        ...             PreprocessNode(),
        ...             RMSDNode(),
        ...         ]
        ...
        >>> pipeline = MyPipeline()
        >>> result = pipeline.execute(context)
    """

    def __init__(self, nodes: List[PipelineNode] = None):
        """
        Initialize pipeline.

        Args:
            nodes: List of pipeline nodes
        """
        self.nodes = nodes or []
        self.logger = logging.getLogger(f"immunoscope.pipeline.{self.__class__.__name__}")

    def execute(self, context: PipelineContext) -> PipelineContext:
        """
        Execute pipeline.

        Executes all nodes in sequence. If a node sets context.should_stop,
        or if an error occurs, execution stops.

        Args:
            context: Pipeline context

        Returns:
            Updated context after execution

        Example:
            >>> context = PipelineContext(
            ...     system_id="1ao7",
            ...     topology="data/1ao7/md.tpr",
            ...     trajectory_raw="data/1ao7/md.xtc"
            ... )
            >>> result_context = pipeline.execute(context)
        """
        self.logger.info(f"Starting pipeline for system {context.system_id}")

        for i, node in enumerate(self.nodes):
            # Check if should stop
            if context.should_stop:
                self.logger.warning(
                    f"Pipeline stopped at node {i}/{len(self.nodes)} ({node.name})"
                )
                break

            # Execute node
            try:
                self.logger.info(f"Executing node {i+1}/{len(self.nodes)}: {node.name}")
                context = node.execute(context)

            except Exception as e:
                error_msg = f"Node {node.name} failed with exception: {str(e)}"
                self.logger.exception(error_msg)
                context.add_error(error_msg)
                context.should_stop = True
                break

        # Log completion
        if context.has_errors():
            self.logger.error(
                f"Pipeline completed with {len(context.errors)} errors for {context.system_id}"
            )
        else:
            self.logger.info(
                f"Pipeline completed successfully for {context.system_id}"
            )

        return context

    def add_node(self, node: PipelineNode):
        """
        Add a node to the pipeline.

        Args:
            node: Pipeline node to add
        """
        self.nodes.append(node)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(nodes={len(self.nodes)})"
