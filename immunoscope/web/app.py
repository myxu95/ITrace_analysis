"""FastAPI application for the ImmunoScope web interface."""

from __future__ import annotations

from fastapi import FastAPI
from fastapi.responses import FileResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware

from immunoscope.web.config import settings
from immunoscope.web.routers import jobs, agent_chat, upload, reports, trajectories, design, compare_views, knowledge, settings as settings_router
from immunoscope.web.services.job_store import init_db
from immunoscope.web.services.agent_store import get_agent_store
from immunoscope.agent.audit import init_audit_db
from immunoscope.agent.config import AGENT_AUDIT_DB_PATH
import os

app = FastAPI(title="ImmunoScope Web", version="0.1.0")

# CORS middleware for WebSocket support
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Configure appropriately for production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

STATIC_DIR = settings.BASE_DIR / "static"

app.include_router(jobs.router, prefix="/api")
app.include_router(agent_chat.router, prefix="/api")
app.include_router(upload.router, prefix="/api")
app.include_router(reports.router)
app.include_router(trajectories.router)
app.include_router(design.router)
app.include_router(compare_views.router)
app.include_router(knowledge.router)
app.include_router(settings_router.router, prefix="/api")
app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")

# Mount output directory for serving generated reports
OUTPUT_DIR = settings.BASE_DIR.parent / "output"
if OUTPUT_DIR.exists():
    app.mount("/reports", StaticFiles(directory=OUTPUT_DIR), name="reports")


@app.get("/")
async def index():
    return FileResponse(STATIC_DIR / "index.html")


@app.get("/agent")
async def agent_page():
    """Redirect to the unified web app Agent route."""
    return RedirectResponse(url="/#/agent")


@app.get("/agent/advanced")
async def agent_advanced_page():
    """Legacy advanced Agent route kept as a redirect."""
    return RedirectResponse(url="/#/agent")


@app.get("/comparison")
async def comparison_page():
    """Trajectory comparison interface."""
    comparison_html = STATIC_DIR / "comparison.html"
    if comparison_html.exists():
        return FileResponse(comparison_html)
    return FileResponse(STATIC_DIR / "index.html")


@app.get("/reports")
async def reports_page():
    """Reports browser interface."""
    reports_html = STATIC_DIR / "reports.html"
    if reports_html.exists():
        return FileResponse(reports_html)
    return FileResponse(STATIC_DIR / "index.html")


@app.on_event("startup")
async def startup() -> None:
    await init_db()
    # Initialize agent audit database
    os.makedirs(os.path.dirname(AGENT_AUDIT_DB_PATH), exist_ok=True)
    await init_audit_db(AGENT_AUDIT_DB_PATH)
    # Start agent session store cleanup
    store = get_agent_store()
    await store.start()


@app.on_event("shutdown")
async def shutdown() -> None:
    # Stop agent session store
    store = get_agent_store()
    await store.stop()
