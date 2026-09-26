"""Statistical comparison tools for two-sample analysis.

This module provides statistical tests and effect size calculations
for comparing distributions between two conditions (e.g., Standard MD vs REST2).

Author: ITrace Development Team
Date: 2026-05-03
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy import stats


@dataclass
class StatisticalComparisonResult:
    """Result of statistical comparison between two distributions."""

    mean_a: float
    mean_b: float
    std_a: float
    std_b: float
    median_a: float
    median_b: float
    ci_95_a: tuple[float, float]
    ci_95_b: tuple[float, float]
    p_value: float
    effect_size: float
    significant: bool
    test_method: str
    n_samples_a: int
    n_samples_b: int


class StatisticalComparator:
    """Statistical comparison tools for two-sample analysis."""

    def __init__(self, alpha: float = 0.05, n_bootstrap: int = 1000):
        """
        Initialize statistical comparator.

        Args:
            alpha: Significance level (default 0.05)
            n_bootstrap: Number of bootstrap iterations for CI estimation
        """
        self.alpha = alpha
        self.n_bootstrap = n_bootstrap

    def compare_distributions(
        self,
        data_a: np.ndarray,
        data_b: np.ndarray,
        test_method: Literal["wilcoxon", "ttest", "auto"] = "auto"
    ) -> StatisticalComparisonResult:
        """
        Compare two distributions with statistical tests.

        Args:
            data_a: Data from condition A
            data_b: Data from condition B
            test_method: Statistical test to use
                - "wilcoxon": Wilcoxon rank-sum test (non-parametric)
                - "ttest": Independent t-test (parametric)
                - "auto": Automatically choose based on normality

        Returns:
            StatisticalComparisonResult with all statistics

        Raises:
            ValueError: If data arrays are empty or invalid
        """
        # Validate inputs
        data_a = self._validate_data(data_a, "data_a")
        data_b = self._validate_data(data_b, "data_b")

        # Compute descriptive statistics
        mean_a = float(np.mean(data_a))
        mean_b = float(np.mean(data_b))
        std_a = float(np.std(data_a, ddof=1))
        std_b = float(np.std(data_b, ddof=1))
        median_a = float(np.median(data_a))
        median_b = float(np.median(data_b))

        # Compute confidence intervals
        ci_95_a = self._bootstrap_ci(data_a)
        ci_95_b = self._bootstrap_ci(data_b)

        # Choose test method
        if test_method == "auto":
            test_method = self._choose_test_method(data_a, data_b)

        # Perform statistical test
        if test_method == "wilcoxon":
            p_value = self._wilcoxon_test(data_a, data_b)
        elif test_method == "ttest":
            p_value = self._ttest(data_a, data_b)
        else:
            raise ValueError(f"Unknown test method: {test_method}")

        # Compute effect size (Cohen's d)
        effect_size = self._cohens_d(data_a, data_b)

        # Determine significance
        significant = p_value < self.alpha

        return StatisticalComparisonResult(
            mean_a=mean_a,
            mean_b=mean_b,
            std_a=std_a,
            std_b=std_b,
            median_a=median_a,
            median_b=median_b,
            ci_95_a=ci_95_a,
            ci_95_b=ci_95_b,
            p_value=p_value,
            effect_size=effect_size,
            significant=significant,
            test_method=test_method,
            n_samples_a=len(data_a),
            n_samples_b=len(data_b)
        )

    def _validate_data(self, data: np.ndarray, name: str) -> np.ndarray:
        """Validate and clean input data."""
        if data is None:
            raise ValueError(f"{name} is None")

        data = np.asarray(data, dtype=float)

        if data.size == 0:
            raise ValueError(f"{name} is empty")

        # Remove NaN and Inf
        valid_mask = np.isfinite(data)
        if not np.any(valid_mask):
            raise ValueError(f"{name} contains only NaN or Inf values")

        data = data[valid_mask]

        if len(data) < 2:
            raise ValueError(f"{name} has fewer than 2 valid samples")

        return data

    def _bootstrap_ci(
        self,
        data: np.ndarray,
        confidence: float = 0.95
    ) -> tuple[float, float]:
        """
        Compute bootstrap confidence interval for the mean.

        Args:
            data: Input data
            confidence: Confidence level (default 0.95)

        Returns:
            (lower_bound, upper_bound)
        """
        n = len(data)
        bootstrap_means = np.empty(self.n_bootstrap)

        rng = np.random.default_rng(seed=42)

        for i in range(self.n_bootstrap):
            sample = rng.choice(data, size=n, replace=True)
            bootstrap_means[i] = np.mean(sample)

        alpha = 1 - confidence
        lower = np.percentile(bootstrap_means, alpha / 2 * 100)
        upper = np.percentile(bootstrap_means, (1 - alpha / 2) * 100)

        return (float(lower), float(upper))

    def _choose_test_method(
        self,
        data_a: np.ndarray,
        data_b: np.ndarray
    ) -> Literal["wilcoxon", "ttest"]:
        """
        Automatically choose test method based on normality.

        Uses Shapiro-Wilk test for normality.
        If both distributions are normal, use t-test.
        Otherwise, use Wilcoxon rank-sum test.
        """
        # For small samples (n < 50), always use non-parametric
        if len(data_a) < 50 or len(data_b) < 50:
            return "wilcoxon"

        # Test normality
        try:
            _, p_a = stats.shapiro(data_a)
            _, p_b = stats.shapiro(data_b)

            # If both are normal (p > 0.05), use t-test
            if p_a > 0.05 and p_b > 0.05:
                return "ttest"
        except Exception:
            # If normality test fails, use non-parametric
            pass

        return "wilcoxon"

    def _wilcoxon_test(self, data_a: np.ndarray, data_b: np.ndarray) -> float:
        """
        Perform Wilcoxon rank-sum test (Mann-Whitney U test).

        Non-parametric test for comparing two independent samples.
        """
        try:
            statistic, p_value = stats.mannwhitneyu(
                data_a,
                data_b,
                alternative='two-sided'
            )
            return float(p_value)
        except Exception as e:
            # If test fails, return NaN
            return np.nan

    def _ttest(self, data_a: np.ndarray, data_b: np.ndarray) -> float:
        """
        Perform independent t-test.

        Parametric test for comparing two independent samples.
        Uses Welch's t-test (unequal variances).
        """
        try:
            statistic, p_value = stats.ttest_ind(
                data_a,
                data_b,
                equal_var=False  # Welch's t-test
            )
            return float(p_value)
        except Exception as e:
            # If test fails, return NaN
            return np.nan

    def _cohens_d(self, data_a: np.ndarray, data_b: np.ndarray) -> float:
        """
        Compute Cohen's d effect size.

        Cohen's d = (mean_b - mean_a) / pooled_std

        Interpretation:
            |d| < 0.2: negligible
            0.2 <= |d| < 0.5: small
            0.5 <= |d| < 0.8: medium
            |d| >= 0.8: large
        """
        mean_a = np.mean(data_a)
        mean_b = np.mean(data_b)
        std_a = np.std(data_a, ddof=1)
        std_b = np.std(data_b, ddof=1)

        n_a = len(data_a)
        n_b = len(data_b)

        # Pooled standard deviation
        pooled_std = np.sqrt(
            ((n_a - 1) * std_a**2 + (n_b - 1) * std_b**2) / (n_a + n_b - 2)
        )

        if pooled_std == 0:
            return 0.0

        cohens_d = (mean_b - mean_a) / pooled_std

        return float(cohens_d)

    def interpret_effect_size(self, effect_size: float) -> str:
        """
        Interpret Cohen's d effect size.

        Args:
            effect_size: Cohen's d value

        Returns:
            Interpretation string
        """
        abs_d = abs(effect_size)

        if abs_d < 0.2:
            return "negligible"
        elif abs_d < 0.5:
            return "small"
        elif abs_d < 0.8:
            return "medium"
        else:
            return "large"

    def format_p_value(self, p_value: float) -> str:
        """
        Format p-value for display.

        Args:
            p_value: P-value

        Returns:
            Formatted string (e.g., "p < 0.001", "p = 0.023")
        """
        if np.isnan(p_value):
            return "p = N/A"
        elif p_value < 0.001:
            return "p < 0.001"
        elif p_value < 0.01:
            return f"p < 0.01"
        elif p_value < 0.05:
            return f"p < 0.05"
        else:
            return f"p = {p_value:.3f}"

    def significance_stars(self, p_value: float) -> str:
        """
        Convert p-value to significance stars.

        Args:
            p_value: P-value

        Returns:
            Stars string ("***", "**", "*", "ns")
        """
        if np.isnan(p_value):
            return ""
        elif p_value < 0.001:
            return "***"
        elif p_value < 0.01:
            return "**"
        elif p_value < 0.05:
            return "*"
        else:
            return "ns"
