"""Residue-label canonicalization for the multi-agent pipeline (2026-06-02).

Problem this solves
-------------------
The four round-1 readers each see the SAME physical residue under a
*different label dialect*, because the presenter views they read render
residues differently:

  - the `hotspots` view shows ``ASP92`` with a separate ``Region`` column
    (``CDR3α``) — the chain letter is NOT shown, only the greek region —
    so the LLM tends to emit ``α-ASP92`` / ``α-D92``.
  - the `exposure` view shows ``ASP92 (D)`` / ``GLU98 (E)`` — chain letter
    IS shown, region is NOT — so the LLM emits ``ASP92 (D)``.

The orchestrator's `_aggregate_candidates` keyed candidates on the raw
display string, so ``α-ASP92`` and ``ASP92 (D)`` (the same residue) ended
up as two separate rows — defeating the whole point of deterministic
aggregation. This module collapses every dialect onto one canonical key.

Canonical key
-------------
A residue's identity is ``(chain, resid)``. ``resname`` is treated as
*decoration*: the same residue may be cited as ``D92`` (no resname) by one
agent and ``ASP92`` by another, so it must NOT participate in the key.
``resname`` is retained for display and (when present on both) used as a
consistency cross-check.

The chain is essential — ``resname+resid`` is NOT globally unique
(e.g. ``ASP30`` exists in two chains of a typical TCR-pHLA complex), so we
never drop it.

Greek ↔ chain-letter aliasing
-----------------------------
The hotspots dialect encodes the chain as a greek symbol (``α`` = TCR
alpha chain, ``β`` = TCR beta chain) while the exposure dialect uses the
raw PDB chain letter (``D`` / ``E``). The mapping between them is
*case-specific* (chain letters vary per PDB), so it is loaded from the
case's RRCS summary CSV (``tcr_chain`` ↔ ``chain_id`` columns) rather than
hardcoded. When that CSV is unavailable the parser still works for every
letter-based dialect; only the greek→letter collapse is skipped and greek
labels are keyed under a normalized greek token instead.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional

log = logging.getLogger("immunoscope.agent.multi_agent.residue_id")


# Greek symbols and their English chain words, both pointing at the same
# logical TCR chain. Used both to recognise greek tokens in a label and to
# look the case-specific chain letter up in the alias map.
_GREEK_TO_WORD: Dict[str, str] = {
    "α": "alpha",
    "β": "beta",
    "a": "alpha",   # only when used as an explicit greek stand-in, e.g. "a-D92"
    "b": "beta",
}
_WORD_TO_GREEK: Dict[str, str] = {"alpha": "α", "beta": "β"}

# 3-letter -> 1-letter is intentionally NOT done here: we keep resname as
# given (upper-cased) and never key on it, so AA-code normalization is moot.

# Common rrcs summary locations under a case dir, most-specific first.
_RRCS_SUMMARY_CANDIDATES = (
    "analysis/rrcs/analysis/interactions/rrcs/annotated_rrcs_pair_summary.csv",
    "analysis/rrcs/analysis/interactions/rrcs/rrcs_pair_summary.csv",
)


@dataclass(frozen=True)
class CanonicalResidue:
    """A residue reduced to its case-stable identity.

    ``chain`` is the resolved PDB chain letter when known, otherwise a
    normalized greek token (``"α"`` / ``"β"``) so greek-only labels still
    collapse together within a run. ``key`` is what the aggregator dedupes
    on; ``display`` is for human-facing tables.
    """

    chain: str
    resid: int
    resname: str = ""

    @property
    def key(self) -> tuple[str, int]:
        """Identity key — chain + resid only (resname is decoration)."""
        return (self.chain, self.resid)

    def display(self) -> str:
        """Stable, dialect-free display label, e.g. ``ASP92 (D)`` or ``92 (D)``."""
        head = f"{self.resname}{self.resid}" if self.resname else str(self.resid)
        return f"{head} ({self.chain})"


class ChainAliasMap:
    """Case-specific greek ↔ chain-letter resolver.

    Built from the case's RRCS summary CSV. Maps the words ``alpha`` /
    ``beta`` (and the greek symbols ``α`` / ``β``) to the PDB chain letter
    used for that TCR chain in this case. Best-effort: an empty map simply
    means greek tokens are not collapsed to letters.
    """

    def __init__(self, word_to_letter: Optional[Dict[str, str]] = None) -> None:
        self._word_to_letter: Dict[str, str] = dict(word_to_letter or {})

    @property
    def is_empty(self) -> bool:
        return not self._word_to_letter

    def resolve(self, token: str) -> Optional[str]:
        """Resolve a chain token (``α``/``β``/``alpha``/``beta``/letter) to a
        PDB chain letter, or ``None`` if it cannot be resolved.

        A token that already looks like a plain chain letter (single
        alpha char that is not a recognised greek stand-in) is returned
        verbatim — letter dialects need no map.
        """
        if not token:
            return None
        t = token.strip()
        low = t.lower()
        word = _GREEK_TO_WORD.get(t) or (low if low in ("alpha", "beta") else None)
        if word is not None:
            letter = self._word_to_letter.get(word)
            if letter:
                return letter
            # No case map → fall back to a normalized greek token so greek
            # labels at least collapse together among themselves.
            return _WORD_TO_GREEK.get(word, word)
        # Not greek: treat as a raw chain identifier (letter / digit chain).
        return t

    @classmethod
    def from_case_dir(cls, case_dir: str | Path) -> "ChainAliasMap":
        """Load the alpha/beta → chain-letter map from a case's RRCS CSV.

        Never raises — on any failure it returns an empty map and logs at
        INFO so the pipeline degrades to greek-only collapsing.
        """
        try:
            base = Path(case_dir)
        except Exception:  # noqa: BLE001
            return cls()
        for rel in _RRCS_SUMMARY_CANDIDATES:
            csv_path = base / rel
            if csv_path.is_file():
                mapping = cls._read_map(csv_path)
                if mapping:
                    return cls(mapping)
        log.info("ChainAliasMap: no usable rrcs summary under %s", case_dir)
        return cls()

    @staticmethod
    def _read_map(csv_path: Path) -> Dict[str, str]:
        try:
            import pandas as pd

            df = pd.read_csv(csv_path)
        except Exception as exc:  # noqa: BLE001
            log.info("ChainAliasMap: failed to read %s: %s", csv_path, exc)
            return {}

        mapping: Dict[str, str] = {}
        # The CSV carries (tcr_chain, chain_id) on each side of a pair.
        for chain_col, id_col in (
            ("tcr_chain_1", "chain_id_1"),
            ("tcr_chain_2", "chain_id_2"),
            ("tcr_chain", "chain_id"),
        ):
            if chain_col not in df.columns or id_col not in df.columns:
                continue
            sub = df[[chain_col, id_col]].dropna()
            for word, letter in zip(sub[chain_col], sub[id_col]):
                word = str(word).strip().lower()
                letter = str(letter).strip()
                if word in ("alpha", "beta") and letter:
                    # First mapping wins; conflicting maps are logged once.
                    if word in mapping and mapping[word] != letter:
                        log.warning(
                            "ChainAliasMap: %s maps to both %s and %s; keeping %s",
                            word, mapping[word], letter, mapping[word],
                        )
                    else:
                        mapping.setdefault(word, letter)
        return mapping


# 3-letter (and a handful of non-standard) residue names, used to tell a
# resname prefix (``ASP-D92``) apart from a chain-letter prefix (``D-ASP92``).
_STD_RESNAMES = frozenset(
    """ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR
    TRP TYR VAL SEC PYL MSE HID HIE HIP HSD HSE HSP CYX ASH GLH LYN""".split()
)

# A residue body: optional resname (1-4 letters), required resid.
_RESNAME_RESID = re.compile(r"^\s*(?P<resname>[A-Za-z]{1,4})?\s*(?P<resid>\d+)\s*$")
# Trailing "(X)" chain suffix, e.g. "ASP92 (D)".
_PAREN_CHAIN = re.compile(r"\(([^)]+)\)\s*$")
# Leading "<prefix>-<rest>" or "<prefix>:<rest>", e.g. "α-ASP92", "D-ASP92".
_PREFIX_SPLIT = re.compile(r"^\s*([A-Za-zαβ]+)\s*[-:]\s*(.+)$")


def _is_greek_token(tok: str) -> bool:
    return tok in _GREEK_TO_WORD or tok.lower() in ("alpha", "beta")


# Antibody/TCR region keywords that can appear as a trailing "(CDR3α)" /
# "(FRβ)" annotation on a residue label. These are decoration, never a chain —
# the RecommendationAgent emits labels like "α-GLY93 (CDR3α)" where the chain is
# the leading greek prefix and the paren is the loop name. Without this guard
# the paren would be misread as the chain and the residue silently dropped.
_REGION_PREFIXES = ("CDR", "FR", "FW")


def _is_region_token(tok: str) -> bool:
    return tok.strip().upper().startswith(_REGION_PREFIXES)


def parse_residue_label(
    label: str, alias_map: Optional[ChainAliasMap] = None
) -> Optional[CanonicalResidue]:
    """Parse any reader label dialect into a `CanonicalResidue`.

    Handles, among others::

        "α-ASP92"        -> chain α→(case letter), resid 92, resname ASP
        "ASP92 (D)"      -> chain D, resid 92, resname ASP
        "GLU98 (E)"      -> chain E, resid 98, resname GLU
        "D-ASP92"        -> chain D, resid 92, resname ASP
        "β-W97"          -> chain β→(case letter), resid 97, resname W
        "ASP-D92"        -> chain D, resid 92, resname ASP  (ResidueKey.label form)
        "D92"            -> chain D, resid 92
        "ASP92"          -> chain "", resid 92, resname ASP  (chain unknown)

    Returns ``None`` when no residue number can be recovered. A trailing
    ``[note]`` and surrounding backticks/whitespace are tolerated.
    """
    if not isinstance(label, str):
        return None
    alias = alias_map or ChainAliasMap()
    text = label.strip().strip("`").strip()
    # Drop a trailing bracketed annotation like "[rim_tune baseline row]".
    text = re.sub(r"\s*\[[^\]]*\]\s*$", "", text).strip()
    if not text:
        return None

    chain_token: Optional[str] = None
    resname: str = ""

    # 1) Trailing "(X)" suffix: a chain token ("(D)") is the chain source, but a
    #    region annotation ("(CDR3α)", "(FRβ)") is decoration — drop it and let
    #    the leading greek/letter prefix supply the chain instead.
    m = _PAREN_CHAIN.search(text)
    if m:
        paren_tok = m.group(1).strip()
        text = _PAREN_CHAIN.sub("", text).strip()
        if not _is_region_token(paren_tok):
            chain_token = paren_tok

    # 2) Leading "<prefix>-<rest>": prefix is greek/word, a chain letter, or a
    #    resname (the ResidueKey.label form "ASP-D92"). Disambiguate by type.
    m = _PREFIX_SPLIT.match(text)
    if m:
        prefix, rest = m.group(1).strip(), m.group(2).strip()
        if _is_greek_token(prefix):
            if chain_token is None:
                chain_token = prefix
            text = rest
        elif prefix.upper() in _STD_RESNAMES:
            # prefix is the residue name; the chain (if any) lives in `rest`.
            resname = prefix.upper()
            text = rest
        elif chain_token is None:
            # prefix is a plain chain letter ("D-ASP92").
            chain_token = prefix
            text = rest
        else:
            # We already have a chain from the suffix; treat prefix as resname
            # if plausible, else ignore it.
            text = rest

    return _resolve_body(text, resname, chain_token, alias)


def _resolve_body(
    text: str,
    resname: str,
    chain_token: Optional[str],
    alias: ChainAliasMap,
) -> Optional[CanonicalResidue]:
    """Resolve a residue body (``RESNAME###`` / ``###`` / ``D92``) into a key.

    ``chain_token`` carries any chain already split off a prefix/suffix.
    ``resname`` carries any residue name already split off a prefix.
    """
    m = _RESNAME_RESID.match(text)
    if m:
        body_name = (m.group("resname") or "").strip()
        resid_str = m.group("resid")
    else:
        digits = re.search(r"(\d+)\s*$", text)  # last resort: trailing integer
        if not digits:
            return None
        resid_str = digits.group(1)
        body_name = text[: digits.start()].strip()

    try:
        resid = int(resid_str)
    except ValueError:
        return None

    if body_name:
        up = body_name.upper()
        # Reject region keywords masquerading as residues ("CDR3", "FR2").
        if chain_token is None and resname == "" and (
            up.startswith("CDR") or up.startswith("FR") or up.startswith("FW")
        ):
            return None
        short_non_resname = up not in _STD_RESNAMES and len(body_name) <= 2
        if chain_token is None and short_non_resname:
            # "D92" dialect: short, non-resname leading letters = the chain.
            # This fires even when a resname was already split off a prefix
            # (the "ASP-D92" ResidueKey.label form: resname=ASP, body=D92).
            chain_token = body_name
        elif resname == "":
            # Keep the body as resname (covers 1-letter AA codes like "W" in
            # "β-W97"). When chain already came from a greek/letter prefix and
            # the body happens to be a redundant chain letter, this stores it
            # as a (harmless) resname — the identity key never uses resname.
            resname = up
        # else: body_name is a redundant resname when we already have one.

    chain = ""
    if chain_token is not None:
        chain = alias.resolve(chain_token) or chain_token.strip()

    return CanonicalResidue(chain=chain, resid=resid, resname=resname)


def canonical_key(
    label: str, alias_map: Optional[ChainAliasMap] = None
) -> Optional[tuple[str, int]]:
    """Convenience: parse ``label`` and return its identity key, or ``None``."""
    res = parse_residue_label(label, alias_map)
    return res.key if res is not None else None


__all__ = [
    "CanonicalResidue",
    "ChainAliasMap",
    "canonical_key",
    "parse_residue_label",
]
