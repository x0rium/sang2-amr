"""Structural verification of AMR candidates via ESMFold + Foldseek.

The FINAL filter in the pipeline:
  SANG2 (composition) → candidates → ESMFold (3D) → Foldseek (structural search)
  → AMR fold match? → confirmed / rejected

Replaces aa_resonance.py's crude AA-pair profile with actual 3D structural
comparison against all known protein folds.

Requires internet access (ESMFold API + Foldseek API).
"""

from __future__ import annotations

import requests
import time
import json
import numpy as np
from pathlib import Path

# AMR-related keywords for Foldseek hit classification
AMR_KEYWORDS = [
    "lactam", "carbapenem", "penicillin", "cephalosporin",
    "resistance", "efflux", "multidrug",
    "aminoglycoside", "acetyltransferase", "phosphotransferase",
    "tetracycline", "chloramphenicol", "vancomycin", "colistin",
    "transpeptidase", "metallo-beta", "oxacillinase",
    "dihydrofolate reductase", "sulfonamide",
]

NON_AMR_KEYWORDS = [
    "ntpase", "atpase", "kinase", "repressor", "recombinase",
    "recbcd", "exonuclease", "helicase", "topoisomerase",
    "ribosom", "polymerase", "transposase", "integrase",
    "chaperone", "protease", "lipase", "dehydrogenase",
]


def predict_structure(sequence: str, timeout: int = 120) -> str | None:
    """Predict 3D structure with ESMFold API. Returns PDB text or None."""
    if len(sequence) > 400:
        sequence = sequence[:400]
    try:
        resp = requests.post(
            "https://api.esmatlas.com/foldSequence/v1/pdb/",
            data=sequence,
            headers={"Content-Type": "text/plain"},
            timeout=timeout,
        )
        if resp.status_code == 200:
            return resp.text
    except requests.RequestException:
        pass
    return None


def plddt_from_pdb(pdb_text: str) -> float:
    """Extract mean pLDDT from ESMFold PDB (0-1 scale)."""
    scores = []
    for line in pdb_text.split("\n"):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            scores.append(float(line[60:66]))
    return float(np.mean(scores)) if scores else 0.0


def search_foldseek(pdb_text: str, databases: list[str] | None = None,
                    timeout: int = 300) -> list[dict]:
    """Search Foldseek for structural matches. Returns list of hit dicts."""
    if databases is None:
        databases = ["pdb100", "afdb50"]

    try:
        files = {"q": ("query.pdb", pdb_text, "application/octet-stream")}
        data = {"mode": "3diaa", "database[]": databases}
        resp = requests.post(
            "https://search.foldseek.com/api/ticket",
            files=files, data=data, timeout=60,
        )
        if resp.status_code != 200:
            return []

        ticket_id = resp.json().get("id", "")

        # Poll for completion
        deadline = time.time() + timeout
        while time.time() < deadline:
            time.sleep(5)
            status = requests.get(
                f"https://search.foldseek.com/api/ticket/{ticket_id}",
                timeout=30,
            ).json()
            if status.get("status") == "COMPLETE":
                break
            if status.get("status") == "ERROR":
                return []
        else:
            return []

        # Collect hits from all databases
        all_hits = []
        for db_idx in range(len(databases)):
            try:
                result = requests.get(
                    f"https://search.foldseek.com/api/result/{ticket_id}/{db_idx}",
                    timeout=30,
                ).json()
                for db_result in result.get("results", []):
                    for hit_group in db_result.get("alignments", []):
                        hits = hit_group if isinstance(hit_group, list) else [hit_group]
                        for hit in hits:
                            if isinstance(hit, dict):
                                all_hits.append({
                                    "target": hit.get("target", ""),
                                    "description": hit.get("tDescription", ""),
                                    "evalue": hit.get("eval", 999),
                                    "seqid": hit.get("seqId", 0),
                                })
            except Exception:
                continue

        return sorted(all_hits, key=lambda h: h["evalue"])

    except requests.RequestException:
        return []


def classify_hits(hits: list[dict], top_n: int = 10) -> dict:
    """Classify Foldseek hits as AMR-related or non-AMR.

    Returns:
      amr_score: fraction of top hits matching AMR keywords (0-1)
      nonamr_score: fraction matching non-AMR keywords (0-1)
      differential: amr_score - nonamr_score (-1 to +1)
      top_amr_hit: best AMR-related hit description (or None)
      top_hit: best overall hit description
    """
    top = hits[:top_n]
    if not top:
        return {"amr_score": 0, "nonamr_score": 0, "differential": 0,
                "top_amr_hit": None, "top_hit": None}

    amr_count = 0
    nonamr_count = 0
    top_amr_hit = None

    for h in top:
        desc = h["description"].lower()
        is_amr = any(k in desc for k in AMR_KEYWORDS)
        is_nonamr = any(k in desc for k in NON_AMR_KEYWORDS)

        if is_amr:
            amr_count += 1
            if top_amr_hit is None:
                top_amr_hit = h["description"]
        if is_nonamr:
            nonamr_count += 1

    n = len(top)
    amr_score = amr_count / n
    nonamr_score = nonamr_count / n

    return {
        "amr_score": amr_score,
        "nonamr_score": nonamr_score,
        "differential": amr_score - nonamr_score,
        "top_amr_hit": top_amr_hit,
        "top_hit": top[0]["description"] if top else None,
    }


def verify_candidate(aa_sequence: str) -> dict:
    """Full structural verification pipeline for one candidate.

    1. ESMFold → 3D structure
    2. pLDDT → is it a real protein?
    3. Foldseek → structural search
    4. Classify → AMR or non-AMR?

    Returns dict with all results.
    """
    result = {
        "length": len(aa_sequence),
        "structure_predicted": False,
        "plddt": 0.0,
        "well_folded": False,
        "foldseek_hits": 0,
        "amr_score": 0.0,
        "nonamr_score": 0.0,
        "differential": 0.0,
        "verdict": "UNKNOWN",
        "top_hit": None,
        "top_amr_hit": None,
    }

    # Step 1: Predict structure
    pdb_text = predict_structure(aa_sequence)
    if pdb_text is None:
        result["verdict"] = "STRUCTURE_FAILED"
        return result
    result["structure_predicted"] = True

    # Step 2: Check quality
    plddt = plddt_from_pdb(pdb_text)
    result["plddt"] = plddt
    result["well_folded"] = plddt > 0.7

    if plddt < 0.5:
        result["verdict"] = "DISORDERED"
        return result

    # Step 3: Structural search
    hits = search_foldseek(pdb_text)
    result["foldseek_hits"] = len(hits)

    if not hits:
        result["verdict"] = "NO_STRUCTURAL_MATCHES"
        return result

    # Step 4: Classify
    classification = classify_hits(hits, top_n=10)
    result.update(classification)

    # Verdict
    diff = classification["differential"]
    if diff > 0.1:
        result["verdict"] = "AMR_LIKELY"
    elif diff < -0.1:
        result["verdict"] = "NON_AMR"
    elif classification["top_amr_hit"]:
        result["verdict"] = "AMR_POSSIBLE"
    else:
        result["verdict"] = "UNKNOWN_FUNCTION"

    return result
