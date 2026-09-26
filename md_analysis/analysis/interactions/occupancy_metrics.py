"""
Interaction occupancy and continuity metric utilities.
"""

from __future__ import annotations

import json
from math import gcd
from typing import Iterable


def frame_set_to_segments(frames_seen: Iterable[int]) -> list[list[int]]:
    frames = sorted({int(frame) for frame in frames_seen})
    if not frames:
        return []
    step = _infer_frame_step(frames)
    segments: list[list[int]] = []
    start = frames[0]
    prev = frames[0]
    for frame in frames[1:]:
        if frame == prev + step:
            prev = frame
            continue
        segments.append([start, prev])
        start = frame
        prev = frame
    segments.append([start, prev])
    return segments


def _infer_frame_step(frames: list[int]) -> int:
    if len(frames) < 2:
        return 1
    diffs = [curr - prev for prev, curr in zip(frames, frames[1:]) if curr > prev]
    if not diffs:
        return 1
    step = diffs[0]
    for diff in diffs[1:]:
        step = gcd(step, diff)
    return max(step, 1)


def summarize_frames(frames_seen: Iterable[int]) -> dict[str, int | None]:
    frames = sorted({int(frame) for frame in frames_seen})
    if not frames:
        return {
            "n_segments": 0,
            "max_consecutive_frames": 0,
            "first_frame": None,
            "last_frame": None,
        }
    step = _infer_frame_step(frames)
    n_segments = 1
    current_run = 1
    max_run = 1
    for prev, curr in zip(frames, frames[1:]):
        if curr - prev == step:
            current_run += 1
        else:
            n_segments += 1
            max_run = max(max_run, current_run)
            current_run = 1
    max_run = max(max_run, current_run)
    return {
        "n_segments": n_segments,
        "max_consecutive_frames": max_run,
        "first_frame": int(frames[0]),
        "last_frame": int(frames[-1]),
    }


def summarize_segments(segments: list[list[int]]) -> dict[str, int | None]:
    if not segments:
        return {
            "n_segments": 0,
            "max_consecutive_frames": 0,
            "first_frame": None,
            "last_frame": None,
        }
    return {
        "n_segments": len(segments),
        "max_consecutive_frames": 0,
        "first_frame": int(segments[0][0]),
        "last_frame": int(segments[-1][1]),
    }


def encode_segments(segments: list[list[int]]) -> str:
    return json.dumps(segments, separators=(",", ":"))
