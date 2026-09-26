"""SQLite-backed job store for the ImmunoScope web app."""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import aiosqlite

from immunoscope.web.config import settings
from immunoscope.web.models import JobDetail, JobSummary


async def init_db() -> None:
    settings.DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    async with aiosqlite.connect(settings.DB_PATH) as db:
        await db.execute(
            """
            CREATE TABLE IF NOT EXISTS jobs (
                id TEXT PRIMARY KEY,
                payload TEXT NOT NULL,
                created_at TEXT NOT NULL
            )
            """
        )
        await db.commit()


class JobStore:
    async def save(self, job: JobSummary | JobDetail) -> None:
        payload = job.model_dump_json(exclude={"log", "output_files", "run_summary", "result_index_payload"})
        async with aiosqlite.connect(settings.DB_PATH) as db:
            await db.execute(
                "INSERT OR REPLACE INTO jobs (id, payload, created_at) VALUES (?, ?, ?)",
                (job.id, payload, job.created_at.isoformat()),
            )
            await db.commit()

    async def get(self, job_id: str) -> JobDetail | None:
        async with aiosqlite.connect(settings.DB_PATH) as db:
            cursor = await db.execute("SELECT payload FROM jobs WHERE id = ?", (job_id,))
            row = await cursor.fetchone()
        if not row:
            return None
        data = json.loads(row[0])
        job_dir = settings.JOBS_DIR / job_id
        output_dir = job_dir / "output"
        data["output_files"] = sorted(path.name for path in output_dir.iterdir() if path.is_file()) if output_dir.exists() else []
        log_path = job_dir / "job.log"
        data["log"] = log_path.read_text(encoding="utf-8", errors="ignore") if log_path.exists() else ""
        run_summary_path = job_dir / "run" / "run_summary.json"
        if run_summary_path.exists():
            try:
                data["run_summary"] = json.loads(run_summary_path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                data["run_summary"] = {}
        result_index_path = job_dir / "job_result_index.json"
        if result_index_path.exists():
            try:
                data["result_index_payload"] = json.loads(result_index_path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                data["result_index_payload"] = {}
        return JobDetail.model_validate(data)

    async def list_all(self) -> list[JobSummary]:
        async with aiosqlite.connect(settings.DB_PATH) as db:
            cursor = await db.execute("SELECT payload FROM jobs ORDER BY created_at DESC")
            rows = await cursor.fetchall()
        return [JobSummary.model_validate(json.loads(row[0])) for row in rows]

    async def update(self, job_id: str, **changes) -> JobDetail | None:
        job = await self.get(job_id)
        if not job:
            return None
        data = job.model_dump()
        data.update(changes)
        updated = JobDetail.model_validate(data)
        await self.save(updated)
        return updated

    async def delete(self, job_id: str) -> None:
        async with aiosqlite.connect(settings.DB_PATH) as db:
            await db.execute("DELETE FROM jobs WHERE id = ?", (job_id,))
            await db.commit()


job_store = JobStore()
