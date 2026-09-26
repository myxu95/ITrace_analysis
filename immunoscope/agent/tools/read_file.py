"""Read file tool for ImmunoScope Agent."""

from __future__ import annotations
from pathlib import Path
from typing import ClassVar

from pydantic import BaseModel, Field

from immunoscope.agent.tool import Tool, ToolContext, ToolResult
from immunoscope.agent.tools import register_tool


class ReadFileInput(BaseModel):
    """Input schema for read_file tool."""

    file_path: str = Field(
        ...,
        description="Path to the file to read",
    )
    max_lines: int = Field(
        default=100,
        ge=1,
        le=1000,
        description="Maximum number of lines to read (default: 100, max: 1000)",
    )
    start_line: int = Field(
        default=0,
        ge=0,
        description="Line number to start reading from (0-indexed)",
    )


@register_tool
class ReadFileTool(Tool):
    """Read contents of a text file."""

    name: ClassVar[str] = "read_file"
    description: ClassVar[str] = (
        "Read the contents of a text file. Use this to inspect analysis results, "
        "configuration files, or data files. Supports reading specific line ranges. "
        "Good for CSV, JSON, log files, and other text formats. Returns file content "
        "with line numbers."
    )
    Input: ClassVar[type[BaseModel]] = ReadFileInput

    is_read_only: ClassVar[bool] = True
    is_concurrency_safe: ClassVar[bool] = True
    needs_permission: ClassVar[bool] = False
    timeout_seconds: ClassVar[float] = 10.0

    async def call(self, args: ReadFileInput, ctx: ToolContext) -> ToolResult:
        """Read file contents."""
        try:
            file_path = Path(args.file_path).resolve()

            if not file_path.exists():
                return ToolResult(
                    content=f"File does not exist: {args.file_path}",
                    is_error=True,
                )

            if not file_path.is_file():
                return ToolResult(
                    content=f"Path is not a file: {args.file_path}",
                    is_error=True,
                )

            # Check file size
            size_mb = file_path.stat().st_size / (1024 * 1024)
            if size_mb > 10:
                return ToolResult(
                    content=f"File too large to read: {size_mb:.1f} MB (max 10 MB)",
                    is_error=True,
                )

            # Read file
            try:
                with open(file_path, "r", encoding="utf-8") as f:
                    lines = f.readlines()
            except UnicodeDecodeError:
                return ToolResult(
                    content=f"File is not a text file (binary content detected)",
                    is_error=True,
                )

            total_lines = len(lines)
            start = args.start_line
            end = min(start + args.max_lines, total_lines)

            if start >= total_lines:
                return ToolResult(
                    content=f"Start line {start} exceeds file length ({total_lines} lines)",
                    is_error=True,
                )

            # Format with line numbers
            content_lines = []
            for i in range(start, end):
                content_lines.append(f"{i+1:6d}\t{lines[i].rstrip()}")

            content = "\n".join(content_lines)

            return ToolResult(content={
                "file_path": str(file_path),
                "total_lines": total_lines,
                "start_line": start,
                "end_line": end,
                "lines_read": end - start,
                "truncated": end < total_lines,
                "content": content,
            })

        except Exception as e:
            return ToolResult(
                content=f"Failed to read file: {str(e)}",
                is_error=True,
            )
