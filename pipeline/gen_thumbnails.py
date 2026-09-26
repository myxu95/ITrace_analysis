"""Render a cartoon structure thumbnail per trajectory with headless Mol*.

For each trajectory this drives a headless Chrome instance (via the Chrome
DevTools Protocol) to load ``frontend/_thumb.html?id=<id>``, which renders the
first-frame structure as a clean cartoon on a white background, then captures
the canvas to ``<web_data>/<id>/thumb.png``. The card view on the list page
shows these thumbnails.

Requirements: ``google-chrome`` on PATH and the ``websockets`` package. The dev
server must be running (the render page fetches ``/data/<id>/topology.pdb`` from
it) -- pass its URL with ``--base-url``.

Usage:
    # with `uvicorn backend.app:app --port 8011` already running:
    python -m pipeline.gen_thumbnails
    python -m pipeline.gen_thumbnails --ids 1ao7_run2 2ak4_run2
"""
from __future__ import annotations

import argparse
import asyncio
import base64
import json
import shutil
import subprocess
import sys
import time
import urllib.request

import websockets

from . import config

CHROME_FLAGS = [
    "--headless=new", "--no-sandbox", "--use-gl=angle", "--use-angle=swiftshader",
    "--enable-webgl", "--ignore-gpu-blocklist", "--window-size=700,500",
    "--remote-allow-origins=*",
]


def _devtools_ws(port: int) -> str:
    for _ in range(60):
        try:
            data = json.load(urllib.request.urlopen(f"http://127.0.0.1:{port}/json"))
            pages = [t for t in data if t["type"] == "page"]
            if pages:
                return pages[0]["webSocketDebuggerUrl"]
        except Exception:
            pass
        time.sleep(0.25)
    raise RuntimeError("Chrome DevTools endpoint did not come up")


async def _render_all(ids: list[str], base_url: str, port: int) -> tuple[int, int]:
    ws = await websockets.connect(_devtools_ws(port), max_size=None)
    counter = 0
    ok = fail = 0

    async def cmd(method, params=None):
        nonlocal counter
        counter += 1
        mid = counter
        await ws.send(json.dumps({"id": mid, "method": method, "params": params or {}}))
        while True:
            msg = json.loads(await ws.recv())
            if msg.get("id") == mid:
                return msg

    async def ev(expr):
        r = await cmd("Runtime.evaluate", {"expression": expr, "returnByValue": True})
        return r.get("result", {}).get("result", {}).get("value")

    await cmd("Page.enable")
    await cmd("Runtime.enable")
    for i, tid in enumerate(ids):
        try:
            await cmd("Page.navigate", {"url": f"{base_url}/_thumb.html?id={tid}"})
            ready = None
            for _ in range(60):
                ready = await ev("window.__thumbReady")
                if ready:
                    break
                await asyncio.sleep(0.3)
            if ready is not True:
                fail += 1
                print(f"[{i + 1}/{len(ids)}] skip {tid} (ready={ready})")
                continue
            rect = await ev("(function(){var r=document.getElementById('app')"
                            ".getBoundingClientRect();return {x:r.x,y:r.y,w:r.width,h:r.height};})()")
            shot = await cmd("Page.captureScreenshot", {"format": "png", "clip": {
                "x": rect["x"], "y": rect["y"], "width": rect["w"], "height": rect["h"], "scale": 2}})
            data = shot.get("result", {}).get("data")
            if not data:
                fail += 1
                print(f"[{i + 1}/{len(ids)}] no screenshot {tid}")
                continue
            out_dir = config.WEB_DATA / tid
            out_dir.mkdir(parents=True, exist_ok=True)
            (out_dir / "thumb.png").write_bytes(base64.b64decode(data))
            ok += 1
            if (i + 1) % 20 == 0 or i == 0:
                print(f"[{i + 1}/{len(ids)}] ok {tid}")
        except Exception as exc:  # noqa: BLE001 - keep the batch going
            fail += 1
            print(f"[{i + 1}/{len(ids)}] error {tid}: {exc}")
    await ws.close()
    return ok, fail


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ids", nargs="*", help="trajectory ids (default: all in web_data)")
    ap.add_argument("--base-url", default="http://127.0.0.1:8011",
                    help="URL of the running dev server")
    ap.add_argument("--port", type=int, default=9222, help="Chrome remote-debugging port")
    args = ap.parse_args(argv)

    chrome = shutil.which("google-chrome") or shutil.which("chromium")
    if not chrome:
        print("google-chrome / chromium not found on PATH", file=sys.stderr)
        return 1

    ids = args.ids or [d.name for d in sorted(config.WEB_DATA.glob("*")) if d.is_dir()]

    proc = subprocess.Popen(
        [chrome, *CHROME_FLAGS, f"--remote-debugging-port={args.port}", "about:blank"],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    try:
        ok, fail = asyncio.run(_render_all(ids, args.base_url.rstrip("/"), args.port))
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()
    print(f"\nThumbnails: {ok} ok, {fail} failed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
