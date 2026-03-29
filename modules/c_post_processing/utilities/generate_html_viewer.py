#!/usr/bin/env python3
"""
generate_html_viewer.py — Auto-generate an interactive HTML viewer for RNA-seq
alignment result heatmaps stored under 3_POST_PROC/.

Usage:
    python generate_html_viewer.py <post_proc_dir> [--output <html_path>]

Path convention parsed (relative to post_proc_dir):
    {method}/Figure_Outputs/{analysis}/{reference}/{gene_group}/
        {processing_level}/{count_type}/{label_type}/{norm_scheme}/
        {row_orientation}/{sort_order}/{filename}.png
"""

import argparse
import datetime
import hashlib
import json
import os
import shutil
import sys
from pathlib import Path

# ── Path segment indices (after stripping post_proc_dir) ──────────────────────
#  0  method           e.g. M1_HISAT2_RefGuided
#  1  "Figure_Outputs" (literal, skipped)
#  2  analysis         e.g. I_Basic_Heatmap
#  3  reference        e.g. Eggplant_V4.1
#  4  gene_group       e.g. Selected_GRF_GIF_…
#  5  processing_level e.g. gene_level
#  6  count_type       e.g. tpm
#  7  label_type       e.g. Shortened_Name
#  8  norm_scheme      e.g. zscore
#  9  row_orientation  e.g. Genes_as_Rows
# 10  sort_order       e.g. Original_Order
# 11  filename         e.g. *.png


def scan_figures(post_proc_dir: Path) -> list[dict]:
    """Walk post_proc_dir and return list of parsed image records."""
    records = []
    # Walk only method subdirectories, skipping _viewer_cache/ which can contain
    # thousands of hardlinked PNGs that would all be discarded by the filter below.
    # On WSL2, each stat() for files in _viewer_cache costs 5-20ms.
    for subdir in sorted(post_proc_dir.iterdir()):
        if not subdir.is_dir() or subdir.name.startswith("_"):
            continue
        for png in sorted(subdir.rglob("*.png")):
            rel = png.relative_to(post_proc_dir)
            parts = rel.parts
            if len(parts) < 12 or parts[1] != "Figure_Outputs":
                continue
            records.append(
                {
                    "path": rel.as_posix(),
                    "method": parts[0],
                    "analysis": parts[2],
                    "reference": parts[3],
                    "gene_group": parts[4],
                    "processing_level": parts[5],
                    "count_type": parts[6],
                    "label_type": parts[7],
                    "norm_scheme": parts[8],
                    "row_orientation": parts[9],
                    "sort_order": parts[10],
                }
            )
    return records


def build_manifest(records: list[dict]) -> dict:
    # Single O(R) pass extracts all 9 dimension sets simultaneously
    # (was 9× O(R) with separate unique_sorted calls per key)
    dim_keys = ("method", "analysis", "reference", "gene_group",
                "processing_level", "count_type", "label_type",
                "norm_scheme", "row_orientation", "sort_order")
    dims: dict[str, set[str]] = {k: set() for k in dim_keys}
    for r in records:
        for k in dim_keys:
            dims[k].add(r[k])
    return {
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "total_images": len(records),
        "dimensions": {
            "methods": sorted(dims["method"]),
            "analyses": sorted(dims["analysis"]),
            "references": sorted(dims["reference"]),
            "gene_groups": sorted(dims["gene_group"]),
            "processing_levels": sorted(dims["processing_level"]),
            "count_types": sorted(dims["count_type"]),
            "label_types": sorted(dims["label_type"]),
            "norm_schemes": sorted(dims["norm_scheme"]),
            "row_orientations": sorted(dims["row_orientation"]),
            "sort_orders": sorted(dims["sort_order"]),
        },
        "images": records,
    }


def _win_long(p: Path) -> str:
    """Return a string path with \\\\?\\ prefix on Windows for long path support."""
    s = str(p)
    if sys.platform == "win32" and not s.startswith("\\\\?\\"):
        return "\\\\?\\" + s
    return s


def create_viewer_cache(post_proc_dir: Path, records: list[dict]) -> None:
    """Create a flat cache of short-named hardlinks for file:// compatibility.

    Windows has a 260-character MAX_PATH limit.  The pipeline's deep directory
    structure routinely produces paths of 370-400+ characters, which prevents
    browsers from loading images via ``file://`` URLs.

    This function creates ``_viewer_cache/`` inside *post_proc_dir* containing
    one hardlink (falling back to copy) per image, named by a 12-hex-digit MD5
    of the original relative path.  Manifest records are updated in-place so
    ``img.path`` points to the short path and ``img.original_path`` retains
    the original for display / download filename purposes.

    Incremental: existing cache entries are reused if the cached file already
    exists and the source file has not been modified since.  Stale entries
    (present in cache but not referenced by any current record) are removed.
    This avoids the O(R) full-rebuild cost on re-generation runs where only
    a few images have changed — significant on WSL2 where filesystem ops are
    expensive (~20-60 ms per hardlink/copy vs ~1 ms on native Linux).
    """
    cache_dir = post_proc_dir / "_viewer_cache"
    cache_dir.mkdir(exist_ok=True)

    # O(R) pass: compute short names and link/copy only new or stale entries
    used_names: set[str] = set()
    created = 0
    reused = 0
    for rec in records:
        original_rel = rec["path"]  # e.g. "M1_.../file.png"
        src = post_proc_dir / original_rel

        # Deterministic short name from path hash
        h = hashlib.md5(original_rel.encode("utf-8")).hexdigest()[:12]
        short_name = f"{h}.png"

        # Handle (very unlikely) hash collisions
        if short_name in used_names:
            for i in range(1, 1000):
                short_name = f"{h}_{i}.png"
                if short_name not in used_names:
                    break
            else:
                raise RuntimeError(f"Could not resolve hash collision for {original_rel}")
        used_names.add(short_name)

        dst = cache_dir / short_name

        # Incremental: skip if cached file exists and is not older than source.
        # Single stat() replaces exists() + stat() — 1 syscall instead of 2 per cached file.
        # O(C) where C = cached files; saves ~C stat syscalls on WSL2 cross-fs mounts.
        try:
            _dst_mtime = dst.stat().st_mtime
        except OSError:
            _dst_mtime = None
        if _dst_mtime is not None:
            try:
                if _dst_mtime >= src.stat().st_mtime:
                    rec["original_path"] = original_rel
                    rec["path"] = f"_viewer_cache/{short_name}"
                    reused += 1
                    continue
            except OSError:
                pass  # src stat failed — fall through to recreate
            # Stale: source is newer, remove old cached entry
            try:
                dst.unlink()
            except OSError:
                pass

        # Use \\?\ prefix on Windows to bypass 260-char MAX_PATH
        src_s = _win_long(src)
        dst_s = _win_long(dst)
        try:
            os.link(src_s, dst_s)
        except OSError:
            # Hardlink may fail (cross-device, FAT32, permissions) — fall back
            try:
                shutil.copy2(src_s, dst_s)
            except OSError as e:
                print(f"WARNING: Could not cache {original_rel}: {e}", file=sys.stderr)
                continue

        rec["original_path"] = original_rel
        rec["path"] = f"_viewer_cache/{short_name}"
        created += 1

    # Remove stale cache entries not referenced by any current record.
    # O(E) where E = existing cache files.
    removed = 0
    try:
        for cached in cache_dir.iterdir():
            if cached.name not in used_names:
                try:
                    cached.unlink()
                    removed += 1
                except OSError:
                    pass
    except OSError:
        pass

    print(f"  Cache: {reused} reused, {created} created, {removed} stale removed",
          file=sys.stderr)


# ── HTML template ─────────────────────────────────────────────────────────────

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>RNA-seq Alignment Results Viewer</title>
<style>
:root {
  --bg: #0d1117; --surface: #161b22; --surface2: #21262d;
  --border: #30363d; --accent: #238636; --accent2: #1f6feb;
  --text: #e6edf3; --muted: #8b949e; --warn: #d29922;
  --tag-bg: #1c2128; --hover: #2d333b; --radius: 6px;
  --method-colors: #79c0ff,#7ee787,#ffa657,#f78166,#d2a8ff;
}
* { box-sizing: border-box; margin: 0; padding: 0; }
body { background: var(--bg); color: var(--text); font-family: -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; font-size: 14px; min-height: 100vh; }

/* ── Header ── */
header { background: var(--surface); border-bottom: 1px solid var(--border); padding: 14px 20px; display: flex; align-items: center; gap: 16px; position: sticky; top: 0; z-index: 100; }
header h1 { font-size: 16px; font-weight: 600; color: var(--text); white-space: nowrap; }
.badge { background: var(--accent); color: #fff; font-size: 11px; padding: 2px 8px; border-radius: 20px; font-weight: 600; }
.spacer { flex: 1; }
#stats { font-size: 12px; color: var(--muted); }
#refresh-btn { background: var(--surface2); border: 1px solid var(--border); color: var(--text); padding: 5px 12px; border-radius: var(--radius); cursor: pointer; font-size: 12px; }
#refresh-btn:hover { background: var(--hover); }

/* ── Layout ── */
.layout { display: flex; height: calc(100vh - 53px); }

/* ── Sidebar ── */
aside { width: 260px; min-width: 220px; background: var(--surface); border-right: 1px solid var(--border); overflow-y: auto; flex-shrink: 0; }
.sidebar-section { border-bottom: 1px solid var(--border); padding: 14px; }
.sidebar-section h3 { font-size: 11px; font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em; color: var(--muted); margin-bottom: 10px; }
.filter-group { margin-bottom: 10px; }
.filter-group label { display: block; font-size: 12px; color: var(--muted); margin-bottom: 4px; }
.filter-group select { width: 100%; background: var(--surface2); border: 1px solid var(--border); color: var(--text); padding: 5px 8px; border-radius: var(--radius); font-size: 13px; appearance: none; cursor: pointer; }
.filter-group select:focus { outline: none; border-color: var(--accent2); }
.filter-group select option { background: var(--surface2); }
.chip-group { display: flex; flex-wrap: wrap; gap: 5px; }
.chip { background: var(--tag-bg); border: 1px solid var(--border); color: var(--muted); padding: 3px 10px; border-radius: 20px; font-size: 12px; cursor: pointer; transition: all .15s; user-select: none; }
.chip.active { background: var(--accent2); border-color: var(--accent2); color: #fff; }
.chip:hover:not(.active) { background: var(--hover); color: var(--text); }
.reset-btn { width: 100%; background: var(--surface2); border: 1px solid var(--border); color: var(--muted); padding: 7px; border-radius: var(--radius); cursor: pointer; font-size: 12px; margin-top: 4px; transition: all .15s; }
.reset-btn:hover { background: var(--hover); color: var(--text); }

/* ── Main content ── */
main { flex: 1; overflow: auto; padding: 16px; }

/* ── Grid ── */
.grid-wrap { overflow: auto; }
table.grid { border-collapse: collapse; }
table.grid th, table.grid td { border: 1px solid var(--border); }
table.grid th { background: var(--surface); padding: 8px 12px; font-size: 12px; font-weight: 600; white-space: nowrap; color: var(--muted); text-align: center; }
table.grid td { padding: 6px; vertical-align: top; background: var(--surface); overflow: hidden; }
table.grid td.row-header { background: var(--surface2); padding: 8px 12px; font-size: 12px; font-weight: 600; white-space: nowrap; color: var(--text); min-width: 160px; }
.method-label { display: flex; align-items: center; gap: 6px; }
.method-dot { width: 10px; height: 10px; border-radius: 50%; flex-shrink: 0; }

/* ── Cell ── */
.cell { position: relative; border-radius: var(--radius); overflow: hidden; cursor: pointer; transition: transform .15s; background: var(--surface2); }
.cell:hover { transform: scale(1.02); z-index: 2; }
.cell img { width: 100%; display: block; border-radius: var(--radius); }
.cell-empty { height: 80px; display: flex; align-items: center; justify-content: center; color: var(--muted); font-size: 11px; border: 1px dashed var(--border); border-radius: var(--radius); }
.cell-count { position: absolute; top: 4px; right: 4px; background: rgba(0,0,0,.6); color: var(--muted); font-size: 10px; padding: 1px 5px; border-radius: 10px; }
.cell-nav { position: absolute; bottom: 4px; left: 50%; transform: translateX(-50%); display: flex; gap: 3px; opacity: 0; transition: opacity .15s; }
.cell:hover .cell-nav { opacity: 1; }
.cell-nav button { background: rgba(0,0,0,.7); border: none; color: #fff; width: 18px; height: 18px; border-radius: 50%; cursor: pointer; font-size: 10px; display: flex; align-items: center; justify-content: center; }
.cell-nav button:hover { background: var(--accent2); }

/* ── No results ── */
.no-results { text-align: center; padding: 60px 20px; color: var(--muted); }
.no-results h2 { font-size: 18px; margin-bottom: 8px; }

/* ── Modal ── */
#modal { display: none; position: fixed; inset: 0; z-index: 1000; background: rgba(0,0,0,.85); backdrop-filter: blur(4px); align-items: center; justify-content: center; }
#modal.open { display: flex; }
.modal-box { background: var(--surface); border: 1px solid var(--border); border-radius: 10px; max-width: 90vw; max-height: 92vh; display: flex; flex-direction: column; overflow: hidden; }
.modal-header { padding: 12px 16px; border-bottom: 1px solid var(--border); display: flex; align-items: flex-start; gap: 12px; }
.modal-meta { flex: 1; }
.modal-title { font-size: 14px; font-weight: 600; margin-bottom: 6px; }
.modal-tags { display: flex; flex-wrap: wrap; gap: 4px; }
.modal-tag { background: var(--tag-bg); border: 1px solid var(--border); padding: 2px 8px; border-radius: 4px; font-size: 11px; color: var(--muted); }
.modal-actions { display: flex; gap: 8px; flex-shrink: 0; }
.modal-btn { background: var(--surface2); border: 1px solid var(--border); color: var(--text); padding: 6px 14px; border-radius: var(--radius); cursor: pointer; font-size: 13px; text-decoration: none; display: flex; align-items: center; gap: 5px; }
.modal-btn:hover { background: var(--hover); }
.modal-btn.primary { background: var(--accent2); border-color: var(--accent2); color: #fff; }
.modal-btn.primary:hover { background: #1a5bc6; }
.modal-close { background: none; border: none; color: var(--muted); cursor: pointer; font-size: 20px; line-height: 1; padding: 2px 4px; flex-shrink: 0; }
.modal-close:hover { color: var(--text); }
.modal-img-wrap { overflow: auto; flex: 1; padding: 12px; display: flex; align-items: center; justify-content: center; }
.modal-img-wrap img { max-width: 100%; max-height: 80vh; object-fit: contain; border-radius: var(--radius); }
.modal-nav { padding: 10px 16px; border-top: 1px solid var(--border); display: flex; align-items: center; justify-content: space-between; gap: 8px; }
.modal-nav-btn { background: var(--surface2); border: 1px solid var(--border); color: var(--text); padding: 5px 16px; border-radius: var(--radius); cursor: pointer; font-size: 13px; }
.modal-nav-btn:hover:not(:disabled) { background: var(--hover); }
.modal-nav-btn:disabled { opacity: 0.4; cursor: default; }
.modal-pos { font-size: 12px; color: var(--muted); }

/* ── Tabs ── */
.tabs { display: flex; gap: 2px; margin-bottom: 14px; border-bottom: 1px solid var(--border); }
.tab { padding: 7px 16px; cursor: pointer; font-size: 13px; color: var(--muted); border-bottom: 2px solid transparent; margin-bottom: -1px; transition: all .15s; }
.tab.active { color: var(--text); border-bottom-color: var(--accent2); }
.tab:hover:not(.active) { color: var(--text); }

/* ── Column resize ── */
.col-resize-handle { position: absolute; right: -2px; top: 0; bottom: 0; width: 5px; cursor: col-resize; z-index: 10; }
.col-resize-handle:hover, .col-resize-handle.active { background: rgba(31,111,235,0.6); }
body.col-resizing { cursor: col-resize !important; user-select: none; }
body.col-resizing * { cursor: col-resize !important; }

/* ── Zoom controls ── */
.zoom-controls { display: flex; align-items: center; gap: 4px; }
.zoom-btn { background: var(--surface2); border: 1px solid var(--border); color: var(--text); width: 26px; height: 26px; border-radius: var(--radius); cursor: pointer; font-size: 14px; display: flex; align-items: center; justify-content: center; }
.zoom-btn:hover { background: var(--hover); }
#zoom-level { font-size: 11px; color: var(--muted); min-width: 36px; text-align: center; }
</style>
</head>
<body>

<header>
  <h1>RNA-seq Alignment Results</h1>
  <span class="badge" id="img-count">0 images</span>
  <div class="spacer"></div>
  <span id="stats"></span>
  <div class="zoom-controls">
    <button class="zoom-btn" onclick="setZoom(state._zoom - 10)" title="Zoom out">&#8722;</button>
    <span id="zoom-level">100%</span>
    <button class="zoom-btn" onclick="setZoom(state._zoom + 10)" title="Zoom in">&#43;</button>
    <button class="zoom-btn" onclick="setZoom(100)" title="Reset zoom" style="font-size:11px;">1:1</button>
  </div>
  <button id="refresh-btn" onclick="location.reload()">&#8635; Refresh</button>
</header>

<div class="layout">
<aside>
  <div class="sidebar-section">
    <h3>Analysis Type</h3>
    <div class="chip-group" id="filter-analysis"></div>
  </div>
  <div class="sidebar-section">
    <h3>Gene Group</h3>
    <div class="chip-group" id="filter-gene-group"></div>
  </div>
  <div class="sidebar-section">
    <h3>Processing Level</h3>
    <div class="chip-group" id="filter-proc-level"></div>
  </div>
  <div class="sidebar-section">
    <h3>Count Type</h3>
    <div class="chip-group" id="filter-count-type"></div>
  </div>
  <div class="sidebar-section">
    <h3>Label Type</h3>
    <div class="chip-group" id="filter-label-type"></div>
  </div>
  <div class="sidebar-section">
    <h3>Normalization</h3>
    <div class="chip-group" id="filter-norm"></div>
  </div>
  <div class="sidebar-section">
    <h3>Row Orientation</h3>
    <div class="chip-group" id="filter-row-orient"></div>
  </div>
  <div class="sidebar-section">
    <h3>Sort Order</h3>
    <div class="chip-group" id="filter-sort"></div>
  </div>
  <div class="sidebar-section">
    <h3>Display</h3>
    <div class="chip-group">
      <span class="chip" id="toggle-hide-empty" onclick="toggleHideEmpty()">Hide empty columns</span>
    </div>
  </div>
  <div class="sidebar-section">
    <button class="reset-btn" onclick="resetFilters()">&#10006; Reset All Filters</button>
  </div>
</aside>

<main>
  <div class="grid-wrap" id="grid-wrap"></div>
</main>
</div>

<!-- Modal -->
<div id="modal">
  <div class="modal-box">
    <div class="modal-header">
      <div class="modal-meta">
        <div class="modal-title" id="modal-title"></div>
        <div class="modal-tags" id="modal-tags"></div>
      </div>
      <div class="modal-actions">
        <a id="modal-download" class="modal-btn primary" download>&#8659; Download</a>
        <button class="modal-close" onclick="closeModal()">&#215;</button>
      </div>
    </div>
    <div class="modal-img-wrap">
      <img id="modal-img" src="" alt="">
    </div>
    <div class="modal-nav">
      <button class="modal-nav-btn" id="modal-prev" onclick="navModal(-1)">&#8592; Prev</button>
      <span class="modal-pos" id="modal-pos"></span>
      <button class="modal-nav-btn" id="modal-next" onclick="navModal(1)">Next &#8594;</button>
    </div>
  </div>
</div>

<script>
// ── Data injected by generator ─────────────────────────────────────────────
const MANIFEST = __MANIFEST_JSON__;

// ── Method display config ──────────────────────────────────────────────────
const METHOD_META = {
  "M1_HISAT2_RefGuided": { short: "M1 · HISAT2 Ref-Guided", color: "#79c0ff" },
  "M2_HISAT2_DeNovo":    { short: "M2 · HISAT2 De Novo",    color: "#7ee787" },
  "M3_STAR_Align":       { short: "M3 · STAR Align",        color: "#ffa657" },
  "M4_Salmon_Saf":       { short: "M4 · Salmon SAF",        color: "#f78166" },
  "M5_RSEM_Bowtie2":     { short: "M5 · RSEM Bowtie2",      color: "#d2a8ff" },
};
const METHOD_ORDER = ["M1_HISAT2_RefGuided","M2_HISAT2_DeNovo","M3_STAR_Align","M4_Salmon_Saf","M5_RSEM_Bowtie2"];

// ── Accession groups — each becomes one grid column ────────────────────────
// References that share an accession are merged into one cell per method.
// Order here defines left-to-right column order.
const ACCESSION_GROUPS = [
  {
    id:    "V4.1",
    label: "V4.1",
    refs:  ["Eggplant_V4.1", "Eggplant_V4.1_transcripts.function"],
  },
  {
    id:    "GPE001970",
    label: "GPE001970",
    refs:  ["GPE001970_genome", "GPE001970_transcripts"],
  },
];

// Build ref → accession id lookup for fast access
const REF_TO_ACCESSION = {};
ACCESSION_GROUPS.forEach(ag => ag.refs.forEach(r => { REF_TO_ACCESSION[r] = ag.id; }));

// ── State ──────────────────────────────────────────────────────────────────
// Every filter is a Set. Empty Set = no restriction (show all).
// Multiple selections within a filter are OR'd; across filters are AND'd.
const state = {
  analyses_active:          new Set(),
  gene_groups_active:       new Set(),
  processing_levels_active: new Set(),
  count_types_active:       new Set(),
  label_types_active:       new Set(),
  norm_schemes_active:      new Set(),
  row_orientations_active:  new Set(),
  sort_orders_active:       new Set(),
  hide_empty_cols:          false,
  _colWidths:               {},
  _zoom:                    100,
};

// Cell image navigation state (per cell: method×reference)
const cellState = {};  // key → { images:[], idx:0 }

// Modal state
let modalImages = [];
let modalIdx = 0;

// ── Boot ───────────────────────────────────────────────────────────────────
document.addEventListener("DOMContentLoaded", () => {
  const d = MANIFEST.dimensions;

  // All filters start as null (no restriction) so every method is visible.
  // User can narrow down using the sidebar toggles.

  // Info banner when opened via file:// (cached short paths should work fine)
  if (location.protocol === "file:") {
    const banner = document.createElement("div");
    banner.style.cssText = "background:#1c3a6e;color:#79c0ff;padding:7px 16px;font-size:12px;text-align:center;position:sticky;top:53px;z-index:99;";
    banner.textContent = "Opened via file:// — for best experience you can also serve with: bash modules/other_tools/run_post_processing_html_viewer.sh --serve";
    document.body.insertBefore(banner, document.querySelector(".layout"));
  }

  buildChips("filter-analysis",   d.analyses,          "analyses_active");
  buildChips("filter-gene-group", d.gene_groups,       "gene_groups_active",       prettyGeneGroup);
  buildChips("filter-proc-level", d.processing_levels, "processing_levels_active");
  buildChips("filter-count-type", d.count_types,       "count_types_active");
  buildChips("filter-label-type", d.label_types,       "label_types_active");
  buildChips("filter-norm",       d.norm_schemes,      "norm_schemes_active");
  buildChips("filter-row-orient", d.row_orientations,  "row_orientations_active");
  buildChips("filter-sort",       d.sort_orders,       "sort_orders_active");

  document.getElementById("img-count").textContent =
    MANIFEST.total_images + " images";

  // Ctrl+scroll to zoom the grid
  document.getElementById("grid-wrap").addEventListener("wheel", e => {
    if (e.ctrlKey) {
      e.preventDefault();
      setZoom(state._zoom + (e.deltaY < 0 ? 10 : -10));
    }
  }, { passive: false });

  render();
});

// ── Filter helpers ─────────────────────────────────────────────────────────
// All chips are multi-select: each toggles independently within its Set.
// Empty Set = no restriction; non-empty = show only matching values.
// labelFn: optional function to format the chip label (defaults to prettyLabel).
function buildChips(containerId, values, stateKey, labelFn) {
  const fmt = labelFn || prettyLabel;
  const el = document.getElementById(containerId);
  values.forEach(v => {
    const c = document.createElement("span");
    c.className = "chip";
    c.textContent = fmt(v);
    c.dataset.value = v;
    c.onclick = () => {
      if (state[stateKey].has(v)) {
        state[stateKey].delete(v);
        c.classList.remove("active");
      } else {
        state[stateKey].add(v);
        c.classList.add("active");
      }
      render();
    };
    el.appendChild(c);
  });
}

function resetFilters() {
  state.analyses_active          = new Set();
  state.gene_groups_active       = new Set();
  state.processing_levels_active = new Set();
  state.count_types_active       = new Set();
  state.label_types_active       = new Set();
  state.norm_schemes_active      = new Set();
  state.row_orientations_active  = new Set();
  state.sort_orders_active       = new Set();
  state.hide_empty_cols          = false;
  ["filter-analysis","filter-gene-group","filter-proc-level","filter-count-type",
   "filter-label-type","filter-norm","filter-row-orient","filter-sort"].forEach(id => {
    document.getElementById(id).querySelectorAll(".chip")
      .forEach(c => c.classList.remove("active"));
  });
  document.getElementById("toggle-hide-empty").classList.remove("active");
  render();
}

function toggleHideEmpty() {
  state.hide_empty_cols = !state.hide_empty_cols;
  document.getElementById("toggle-hide-empty").classList.toggle("active", state.hide_empty_cols);
  render();
}

// ── Filter images ──────────────────────────────────────────────────────────
function filteredImages() {
  return MANIFEST.images.filter(img => {
    if (state.analyses_active.size > 0          && !state.analyses_active.has(img.analysis))                return false;
    if (state.gene_groups_active.size > 0       && !state.gene_groups_active.has(img.gene_group))           return false;
    if (state.processing_levels_active.size > 0 && !state.processing_levels_active.has(img.processing_level)) return false;
    if (state.count_types_active.size > 0       && !state.count_types_active.has(img.count_type))           return false;
    if (state.label_types_active.size > 0      && !state.label_types_active.has(img.label_type))           return false;
    if (state.norm_schemes_active.size > 0      && !state.norm_schemes_active.has(img.norm_scheme))         return false;
    if (state.row_orientations_active.size > 0  && !state.row_orientations_active.has(img.row_orientation)) return false;
    if (state.sort_orders_active.size > 0       && !state.sort_orders_active.has(img.sort_order))           return false;
    return true;
  });
}

// ── Render grid ────────────────────────────────────────────────────────────
function render() {
  const images = filteredImages();
  const wrap = document.getElementById("grid-wrap");
  document.getElementById("stats").textContent =
    images.length + " / " + MANIFEST.total_images + " shown";

  if (images.length === 0) {
    wrap.innerHTML = `<div class="no-results"><h2>No images match the current filters</h2><p>Try adjusting or resetting the filters.</p></div>`;
    return;
  }

  // Build method → accession_id → gene_group → images lookup
  const lookup = {};
  images.forEach(img => {
    const acc = REF_TO_ACCESSION[img.reference] || img.reference;
    if (!lookup[img.method]) lookup[img.method] = {};
    if (!lookup[img.method][acc]) lookup[img.method][acc] = {};
    const gg = img.gene_group;
    if (!lookup[img.method][acc][gg]) lookup[img.method][acc][gg] = [];
    lookup[img.method][acc][gg].push(img);
  });

  // Determine which accession columns have any data in the current filter
  const activeAccessions = images.length > 0
    ? ACCESSION_GROUPS.filter(ag =>
        MANIFEST.dimensions.methods.some(m => lookup[m] && lookup[m][ag.id]))
    : ACCESSION_GROUPS;
  const cols = activeAccessions.length > 0 ? activeAccessions : ACCESSION_GROUPS;

  // All gene groups present in filtered images (preserving encounter order, deduped)
  // O(n) Set-based dedup instead of O(n²) .includes() scan
  const geneGroupSet = [...new Set(images.map(img => img.gene_group))];
  if (geneGroupSet.length === 0) MANIFEST.dimensions.gene_groups.forEach(g => geneGroupSet.push(g));

  // Always render all known methods as rows (cells are "No data" when empty)
  // O(M) Set lookup instead of O(M×D) .includes() scan
  const _methodSet = new Set(MANIFEST.dimensions.methods);
  const methods = METHOD_ORDER.filter(m => _methodSet.has(m));

  // Per-accession gene group visibility: when hide_empty_cols is on, drop sub-columns
  // where every method has no data for that (accession × gene_group) pair.
  const ggPerAccession = {};
  cols.forEach(ag => {
    if (state.hide_empty_cols) {
      ggPerAccession[ag.id] = geneGroupSet.filter(gg =>
        methods.some(m => (((lookup[m] || {})[ag.id] || {})[gg] || []).length > 0)
      );
    } else {
      ggPerAccession[ag.id] = [...geneGroupSet];
    }
  });

  // Total leaf columns = sum of visible gene groups per accession
  const totalLeafCols = cols.reduce((sum, ag) => sum + ggPerAccession[ag.id].length, 0);

  const table = document.createElement("table");
  table.className = "grid";
  table.style.tableLayout = "fixed";

  // Build colgroup for column width control (enables drag-to-resize)
  const _colgroup = document.createElement("colgroup");
  const _colM = document.createElement("col");
  _colM.style.width = state._colWidths["__method__"] || "160px";
  _colM.dataset.colKey = "__method__";
  _colgroup.appendChild(_colM);
  cols.forEach(ag => {
    ggPerAccession[ag.id].forEach(gg => {
      const _col = document.createElement("col");
      const _ck = ag.id + "||" + gg;
      _col.style.width = state._colWidths[_ck] || "200px";
      _col.dataset.colKey = _ck;
      _colgroup.appendChild(_col);
    });
  });
  table.appendChild(_colgroup);

  // ── Header: 3 rows — Category / Accession / Gene Group ──
  const thead = table.createTHead();

  // Row 1: "Accessions" category label spanning all leaf columns
  const catRow = thead.insertRow();
  const thCorner = document.createElement("th");
  thCorner.rowSpan = 3;
  thCorner.style.cssText = "vertical-align:bottom;min-width:150px;";
  thCorner.textContent = "Method";
  catRow.appendChild(thCorner);
  const thCat = document.createElement("th");
  thCat.colSpan = totalLeafCols;
  thCat.style.cssText = "text-align:center;color:var(--accent2);letter-spacing:.06em;font-size:11px;border-bottom:none;padding-bottom:2px;text-transform:uppercase;";
  thCat.textContent = "Accessions";
  catRow.appendChild(thCat);

  // Row 2: one accession header spanning its visible gene_group sub-columns
  const accRow = thead.insertRow();
  cols.forEach(ag => {
    const visibleCount = ggPerAccession[ag.id].length;
    if (visibleCount === 0) return;  // skip accession with no visible columns
    const th = document.createElement("th");
    th.colSpan = visibleCount;
    th.style.cssText = "text-align:center;font-size:13px;font-weight:700;border-bottom:1px solid var(--border);";
    th.textContent = ag.label;
    accRow.appendChild(th);
  });

  // Row 3: gene group sub-column headers (per accession, only visible groups)
  const ggRow = thead.insertRow();
  cols.forEach(ag => {
    ggPerAccession[ag.id].forEach(gg => {
      const th = document.createElement("th");
      th.style.cssText = "text-align:center;min-width:180px;font-size:11px;color:var(--muted);font-weight:500;";
      th.textContent = prettyGeneGroup(gg);
      ggRow.appendChild(th);
    });
  });

  // Data rows
  const tbody = table.createTBody();
  methods.forEach(method => {
    const row = tbody.insertRow();
    // Row header
    const rh = document.createElement("td");
    rh.className = "row-header";
    const meta = METHOD_META[method] || { short: method, color: "#8b949e" };
    rh.innerHTML = `<div class="method-label"><span class="method-dot" style="background:${meta.color}"></span>${meta.short}</div>`;
    row.appendChild(rh);

    cols.forEach(ag => {
      ggPerAccession[ag.id].forEach(gg => {
        const td = row.insertCell();
        const imgs = ((lookup[method] || {})[ag.id] || {})[gg] || [];
        if (imgs.length === 0) {
          td.innerHTML = `<div class="cell-empty">No data</div>`;
        } else {
          const cellKey = method + "||" + ag.id + "||" + gg;
          if (!cellState[cellKey]) cellState[cellKey] = { images: imgs, idx: 0 };
          else { cellState[cellKey].images = imgs; }
          const cs = cellState[cellKey];
          cs.idx = Math.min(cs.idx, imgs.length - 1);
          td.innerHTML = buildCellHTML(cellKey, cs);
        }
      });
    });
  });

  wrap.innerHTML = "";
  wrap.appendChild(table);
  initColumnResize();
  if (state._zoom !== 100) table.style.zoom = state._zoom / 100;
}

// ── Column resize ─────────────────────────────────────────────────────────
function initColumnResize() {
  const table = document.querySelector("table.grid");
  if (!table || !table.tHead) return;
  const colEls = table.querySelectorAll("colgroup col");

  // Corner th (Method label) — first cell of first header row
  addResizeHandle(table.tHead.rows[0].cells[0], colEls[0]);

  // Gene group header row — last header row; cells map to col indices 1+
  const ggRow = table.tHead.rows[table.tHead.rows.length - 1];
  for (let i = 0; i < ggRow.cells.length; i++) {
    addResizeHandle(ggRow.cells[i], colEls[i + 1]);
  }
}

function addResizeHandle(th, col) {
  if (!th || !col) return;
  th.style.position = "relative";
  th.style.overflow = "hidden";
  const handle = document.createElement("div");
  handle.className = "col-resize-handle";
  th.appendChild(handle);

  handle.addEventListener("mousedown", e => {
    e.preventDefault();
    e.stopPropagation();
    const startX = e.pageX;
    const startW = th.offsetWidth;
    handle.classList.add("active");
    document.body.classList.add("col-resizing");

    const onMove = ev => {
      const w = Math.max(60, startW + ev.pageX - startX);
      col.style.width = w + "px";
      if (col.dataset.colKey) state._colWidths[col.dataset.colKey] = w + "px";
    };
    const onUp = () => {
      document.removeEventListener("mousemove", onMove);
      document.removeEventListener("mouseup", onUp);
      handle.classList.remove("active");
      document.body.classList.remove("col-resizing");
    };
    document.addEventListener("mousemove", onMove);
    document.addEventListener("mouseup", onUp);
  });
}

// ── Zoom ──────────────────────────────────────────────────────────────────
function setZoom(level) {
  state._zoom = Math.max(25, Math.min(300, level));
  const table = document.querySelector("table.grid");
  if (table) table.style.zoom = state._zoom / 100;
  document.getElementById("zoom-level").textContent = state._zoom + "%";
}

function buildCellHTML(cellKey, cs) {
  const img = cs.images[cs.idx];
  const count = cs.images.length;
  const esc = cellKey.replace(/'/g,"\\'");
  const navHTML = count > 1
    ? `<div class="cell-nav">
        <button title="Previous" onclick="event.stopPropagation();cycleCell('${esc}',-1)">&#8592;</button>
        <button title="Next"     onclick="event.stopPropagation();cycleCell('${esc}',1)">&#8594;</button>
       </div>`
    : "";
  return `<div class="cell" data-cell-key="${cellKey.replace(/"/g,'&quot;')}" onclick="openModal('${esc}')">
    <img src="${img.path}" alt="" loading="lazy">
    ${count > 1 ? `<div class="cell-count">${cs.idx+1}/${count}</div>` : ""}
    ${navHTML}
  </div>`;
}

function cycleCell(cellKey, dir) {
  const cs = cellState[cellKey];
  cs.idx = (cs.idx + dir + cs.images.length) % cs.images.length;
  // Find the cell div by data attribute and re-render its parent TD
  const el = document.querySelector(`.cell[data-cell-key="${CSS.escape ? CSS.escape(cellKey) : cellKey}"]`);
  if (el) el.parentElement.innerHTML = buildCellHTML(cellKey, cs);
}

// ── Modal ──────────────────────────────────────────────────────────────────
function openModal(cellKey) {
  const cs = cellState[cellKey];
  modalImages = cs.images;
  modalIdx = cs.idx;
  showModalAt(modalIdx);
  document.getElementById("modal").classList.add("open");
  document.addEventListener("keydown", onModalKey);
}

function showModalAt(idx) {
  const img = modalImages[idx];
  document.getElementById("modal-img").src = img.path;
  document.getElementById("modal-download").href = img.path;
  document.getElementById("modal-download").download = (img.original_path || img.path).split("/").pop();
  document.getElementById("modal-title").textContent = (img.original_path || img.path).split("/").pop().replace(/_/g," ").replace(/\.png$/i,"");
  const tags = document.getElementById("modal-tags");
  tags.innerHTML = ["method","analysis","reference","gene_group","processing_level","count_type","label_type","norm_scheme","row_orientation","sort_order"]
    .map(k => `<span class="modal-tag">${escHTML(prettyLabel(k))}: ${escHTML(prettyLabel(img[k]))}</span>`)
    .join("");
  document.getElementById("modal-pos").textContent = (idx+1) + " / " + modalImages.length;
  document.getElementById("modal-prev").disabled = idx === 0;
  document.getElementById("modal-next").disabled = idx === modalImages.length - 1;
}

function navModal(dir) {
  modalIdx = Math.max(0, Math.min(modalImages.length-1, modalIdx+dir));
  showModalAt(modalIdx);
}

function closeModal() {
  document.getElementById("modal").classList.remove("open");
  document.removeEventListener("keydown", onModalKey);
}

function onModalKey(e) {
  if (e.key === "ArrowRight") navModal(1);
  else if (e.key === "ArrowLeft") navModal(-1);
  else if (e.key === "Escape") closeModal();
}

document.getElementById("modal").addEventListener("click", e => {
  if (e.target === document.getElementById("modal")) closeModal();
});

// ── Util ───────────────────────────────────────────────────────────────────
function escHTML(s) {
  return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;').replace(/'/g,'&#39;');
}
function prettyLabel(s) {
  if (!s) return "";
  return s.replace(/_/g," ").replace(/([a-z])([A-Z])/g,"$1 $2");
}
// Shorten gene group label for the column header (strip dataset suffix)
function prettyGeneGroup(gg) {
  if (!gg) return "";
  // Strip trailing "_in_PRJNA…" dataset suffix to keep headers compact
  return gg.replace(/_in_PRJ[A-Z0-9]+(_\w+)?$/, "").replace(/_/g," ");
}
</script>
</body>
</html>
"""


def generate_html(manifest: dict) -> str:
    json_str = json.dumps(manifest, separators=(",", ":"))
    # Escape </ sequences to prevent </script> breakout in inline JSON
    json_str = json_str.replace("</", r"<\/")
    return HTML_TEMPLATE.replace("__MANIFEST_JSON__", json_str)


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Generate interactive HTML viewer for RNA-seq heatmaps.")
    parser.add_argument("post_proc_dir", help="Path to 3_POST_PROC directory")
    parser.add_argument("--output", default=None, help="Output HTML path (default: <post_proc_dir>/alignment_results_viewer.html)")
    parser.add_argument("--manifest-only", action="store_true", help="Only write viewer_manifest.json, skip HTML")
    args = parser.parse_args()

    post_proc_dir = Path(args.post_proc_dir).resolve()
    if not post_proc_dir.is_dir():
        print(f"ERROR: Directory not found: {post_proc_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Scanning {post_proc_dir} …", file=sys.stderr)
    records = scan_figures(post_proc_dir)
    if not records:
        print("WARNING: No PNG files found under Figure_Outputs/", file=sys.stderr)

    # Build short-path cache so file:// works on Windows (MAX_PATH = 260)
    if records:
        print(f"Creating viewer cache ({len(records)} images) …", file=sys.stderr)
        create_viewer_cache(post_proc_dir, records)

    manifest = build_manifest(records)
    print(f"Found {len(records)} images across {len(manifest['dimensions']['methods'])} methods, "
          f"{len(manifest['dimensions']['references'])} references.", file=sys.stderr)

    # Write manifest JSON (used by viewer for caching)
    manifest_path = post_proc_dir / "viewer_manifest.json"
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    print(f"Manifest written: {manifest_path}", file=sys.stderr)

    if args.manifest_only:
        return

    # Write HTML
    output_path = Path(args.output).resolve() if args.output else post_proc_dir / "alignment_results_viewer.html"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    html = generate_html(manifest)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"HTML viewer written: {output_path}", file=sys.stderr)
    print(str(output_path))


if __name__ == "__main__":
    main()
