"""List files tool for ImmunoScope Agent."""

from __future__ import annotations
import os
from pathlib import Path
from typing import ClassVar

from pydantic import BaseModel, Field

from immunoscope.agent.tool import Tool, ToolContext, ToolResult
from immunoscope.agent.tools import register_tool


class ListFilesInput(BaseModel):
    """Input schema for list_files tool."""

    directory: str = Field(
        default=".",
        description="Directory path to list files from (default: current directory)",
    )
    pattern: str = Field(
        default="*",
        description="Glob pattern to filter files (e.g., '*.xtc', '*.pdb')",
    )
    recursive: bool = Field(
        default=False,
        description="Whether to list files recursively",
    )


@register_tool
class ListFilesTool(Tool):
    """List files in a directory with optional filtering."""

    name: ClassVar[str] = "list_files"
    description: ClassVar[str] = (
        "List files in a directory. Use this to explore available trajectory files, "
        "topology files, or analysis outputs. Supports glob patterns (*.xtc, *.pdb) "
        "and recursive listing. Returns file paths, sizes, and modification times."
    )
    Input: ClassVar[type[BaseModel]] = ListFilesInput

    is_read_only: ClassVar[bool] = True
    is_concurrency_safe: ClassVar[bool] = True
    needs_permission: ClassVar[bool] = False
    timeout_seconds: ClassVar[float] = 10.0

    async def call(self, args: ListFilesInput, ctx: ToolContext) -> ToolResult:
        """List files in directory."""
        try:
            directory = Path(args.directory).resolve()

            if not directory.exists():
                return ToolResult(
                    content=f"Directory does not exist: {args.directory}",
                    is_error=True,
                )

            if not directory.is_dir():
                return ToolResult(
                    content=f"Path is not a directory: {args.directory}",
                    is_error=True,
                )

            # List files
            if args.recursive:
                pattern = f"**/{args.pattern}"
            else:
                pattern = args.pattern

            files = []
            for path in directory.glob(pattern):
                if path.is_file():
                    stat = path.stat()
                    files.append({
                        "path": str(path.relative_to(directory)),
                        "absolute_path": str(path),
                        "size_bytes": stat.st_size,
                        "size_mb": round(stat.st_size / (1024 * 1024), 2),
                        "modified": stat.st_mtime,
                    })

            # Sort by name
            files.sort(key=lambda x: x["path"])

            return ToolResult(content={
                "directory": str(directory),
                "pattern": args.pattern,
                "recursive": args.recursive,
                "count": len(files),
                "files": files[:100],  # Limit to 100 files
                "truncated": len(files) > 100,
            })

        except Exception as e:
            return ToolResult(
                content=f"Failed to list files: {str(e)}",
                is_error=True,
            )
