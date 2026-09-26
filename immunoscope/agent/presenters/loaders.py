"""Cached data loaders for JSON and CSV files."""

from pathlib import Path
from typing import Optional, Dict, Any, List
import json
import pandas as pd


class DataCache:
    """Simple cache for loaded data files."""

    def __init__(self):
        self._json_cache: Dict[Path, Dict] = {}
        self._csv_cache: Dict[Path, pd.DataFrame] = {}

    def get_json(self, path: Optional[Path]) -> Dict[str, Any]:
        """
        Load JSON file with caching.

        Args:
            path: Path to JSON file

        Returns:
            Parsed JSON dict, or empty dict if file missing/invalid
        """
        if not path or not path.exists():
            return {}

        if path in self._json_cache:
            return self._json_cache[path]

        try:
            with open(path) as f:
                data = json.load(f)
            self._json_cache[path] = data
            return data
        except Exception:
            # Return empty dict on error (mirrors QueryAnswerBuilder behavior)
            return {}

    def get_csv(
        self,
        path: Optional[Path],
        usecols: Optional[List[str]] = None
    ) -> pd.DataFrame:
        """
        Load CSV file with caching.

        Args:
            path: Path to CSV file
            usecols: Optional list of columns to load

        Returns:
            DataFrame, or empty DataFrame if file missing/invalid
        """
        if not path or not path.exists():
            return pd.DataFrame()

        # Cache key includes usecols for different column selections
        cache_key = (path, tuple(usecols) if usecols else None)

        if cache_key in self._csv_cache:
            return self._csv_cache[cache_key]

        try:
            df = pd.read_csv(path, usecols=usecols)
            self._csv_cache[cache_key] = df
            return df
        except Exception:
            # Return empty DataFrame on error
            return pd.DataFrame()

    def clear(self):
        """Clear all cached data."""
        self._json_cache.clear()
        self._csv_cache.clear()


# Global cache instance (one per process)
_cache = DataCache()


def read_json(path: Optional[Path]) -> Dict[str, Any]:
    """
    Read JSON file with caching.

    Args:
        path: Path to JSON file

    Returns:
        Parsed JSON dict, or empty dict if missing/invalid
    """
    return _cache.get_json(path)


def read_csv(
    path: Optional[Path],
    usecols: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    Read CSV file with caching.

    Args:
        path: Path to CSV file
        usecols: Optional list of columns to load

    Returns:
        DataFrame, or empty DataFrame if missing/invalid
    """
    return _cache.get_csv(path, usecols)


def clear_cache():
    """Clear the data cache (useful for testing)."""
    _cache.clear()
