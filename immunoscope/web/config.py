"""Configuration for the ImmunoScope web server."""

from __future__ import annotations

from pathlib import Path

from pydantic_settings import BaseSettings


class WebSettings(BaseSettings):
    HOST: str = "127.0.0.1"
    PORT: int = 8890
    DEBUG: bool = True

    BASE_DIR: Path = Path(__file__).resolve().parent
    DATA_DIR: Path = Path("output/web").resolve()
    UPLOAD_DIR: Path = DATA_DIR / "uploads"
    JOBS_DIR: Path = DATA_DIR / "jobs"
    DB_PATH: Path = DATA_DIR / "immunoscope_web.db"

    MAX_UPLOAD_SIZE_MB: int = 2048

    class Config:
        env_prefix = "IMMUNOSCOPE_WEB_"


settings = WebSettings()
settings.UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
settings.JOBS_DIR.mkdir(parents=True, exist_ok=True)
