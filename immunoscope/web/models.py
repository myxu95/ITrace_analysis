"""Pydantic models for ImmunoScope web jobs."""

from __future__ import annotations

import enum
from datetime import datetime

from pydantic import BaseModel, Field


class JobStatus(str, enum.Enum):
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class AnalysisProfile(str, enum.Enum):
    STANDARD = "standard"
    SAMPLING_COMPARE = "sampling_compare"
    CUSTOM = "custom"


class InputMode(str, enum.Enum):
    RAW = "raw"
    PREPARED = "prepared"


class AnalysisConfig(BaseModel):
    profile: AnalysisProfile = AnalysisProfile.STANDARD
    modules: list[str] = Field(default_factory=lambda: ["quality", "identity", "contact", "cluster", "report"])
    notes: str = ""
    input_mode: InputMode = InputMode.RAW
    auto_start: bool = True
    preprocess_method: str = "2step"
    stride: int = 5
    landscape_reducer: str = "pca"
    tica_lag: int = 10
    continue_on_error: bool = True


class JobInputPaths(BaseModel):
    structure: str | None = None
    topology: str | None = None
    trajectory: str | None = None
    prepared_structure: str | None = None
    processed_trajectory: str | None = None
    case_a: str | None = None
    case_b: str | None = None
    label_a: str | None = None
    label_b: str | None = None
    comparison_mode: str | None = None
    comparison_scope: str | None = None
    alignment_selection: str | None = None
    residue_mapping: str | None = None


class JobSummary(BaseModel):
    id: str
    name: str
    status: JobStatus
    config: AnalysisConfig
    input_files: list[str] = Field(default_factory=list)
    created_at: datetime
    started_at: datetime | None = None
    finished_at: datetime | None = None
    progress: int = 0
    error: str | None = None
    job_dir: str | None = None
    run_dir: str | None = None
    report_html: str | None = None
    result_index: str | None = None


class JobDetail(JobSummary):
    inputs: JobInputPaths = Field(default_factory=JobInputPaths)
    log: str = ""
    output_files: list[str] = Field(default_factory=list)
    run_summary: dict = Field(default_factory=dict)
    result_index_payload: dict = Field(default_factory=dict)


class AssistantQuery(BaseModel):
    question: str
    skill: str | None = None
    intent: str | None = None
    case_b: str | None = None
    label_a: str | None = None
    label_b: str | None = None


class AssistantSkillInvoke(BaseModel):
    skill: str
    question: str = ""
    intent: str | None = None
    parameters: dict = Field(default_factory=dict)


class CompareJobRequest(BaseModel):
    case_a: str
    case_b: str
    label_a: str = "Condition A"
    label_b: str = "Condition B"
    job_name: str = ""
    comparison_mode: str = "sampling"
    comparison_scope: str = "same-system"
    alignment_selection: str = "phla_core_ca"
    residue_mapping: str = "auto"
    comparison_context: str = ""
    auto_start: bool = True
