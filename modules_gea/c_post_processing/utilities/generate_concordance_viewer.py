#!/usr/bin/env python3
"""
generate_concordance_viewer.py — Auto-generate an interactive HTML viewer for
concordance analysis results stored under II_RESULTS/4_CONCORDANCE_ANALYSIS/.

Usage:
    python generate_concordance_viewer.py <concordance_dir> [--output <html_path>]

Scans for:
    {mode_dir}/figures/*.png          — Analysis figures
    {mode_dir}/tables/*.csv           — Data tables
    {mode_dir}/concordance_report.md  — Markdown reports
"""

import argparse
import datetime
import json
import sys
from pathlib import Path

# ── Figure classification ────────────────────────────────────────────────────

FIGURE_PREFIXES = [
    # (stem_prefix, analysis_type_label, has_gene_group_suffix)
    ("method_concordance_heatmap_spearman", "Concordance Heatmap",     False),
    ("equivalent_gene_heatmap_",           "Equivalent Gene Heatmap",  True),
    ("equivalent_gene_scatter_",           "Equivalent Gene Scatter",  True),
    ("gene_concordance_heatmap_",          "Gene Concordance Heatmap", True),
]

MODE_LABELS = {
    "cross_methods_vs_methods_concordance":                     "Cross-Method",
    "cross_genomes_vs_genomes_concordance":                     "Cross-Genome",
    "cross_equivalent_gene_between_genomes_concordance":        "Cross-Equivalent Gene",
    "cross_within_gene_group_vs_within_gene_group_concordance": "Cross Within Gene-Group",
}


def classify_figure(filename: str) -> tuple[str, str]:
    """Return (analysis_type_label, gene_group|'—')."""
    stem = Path(filename).stem
    for prefix, label, has_gg in FIGURE_PREFIXES:
        if has_gg and stem.startswith(prefix):
            return label, stem[len(prefix):]
        if not has_gg and stem == prefix:
            return label, "\u2014"
    return "Other", "\u2014"


def mode_label(dirname: str) -> str:
    return MODE_LABELS.get(dirname, dirname.replace("_", " ").title())


# ── Scanning ─────────────────────────────────────────────────────────────────

def scan_concordance(base_dir: Path) -> tuple[list, list, list]:
    """Return (figures, tables, reports) from concordance output tree."""
    figures, tables, reports = [], [], []

    for mode_dir in sorted(base_dir.iterdir()):
        if not mode_dir.is_dir() or mode_dir.name in ("logs", "_viewer_cache"):
            continue

        mode = mode_label(mode_dir.name)
        mode_id = mode_dir.name

        # Figures
        fig_dir = mode_dir / "figures"
        if fig_dir.is_dir():
            for png in sorted(fig_dir.glob("*.png")):
                atype, gg = classify_figure(png.name)
                figures.append({
                    "path": png.relative_to(base_dir).as_posix(),
                    "mode": mode,
                    "mode_id": mode_id,
                    "analysis_type": atype,
                    "gene_group": gg,
                    "filename": png.name,
                })

        # Tables (inline CSV content — files are small)
        tbl_dir = mode_dir / "tables"
        if tbl_dir.is_dir():
            for csv_f in sorted(tbl_dir.glob("*.csv")):
                try:
                    content = csv_f.read_text(encoding="utf-8", errors="replace")
                except Exception as e:
                    print(f"  [WARN] Could not read {csv_f}: {e}")
                    content = ""
                tables.append({
                    "path": csv_f.relative_to(base_dir).as_posix(),
                    "mode": mode,
                    "mode_id": mode_id,
                    "filename": csv_f.name,
                    "content": content,
                })

        # Markdown report
        rpt = mode_dir / "concordance_report.md"
        if rpt.is_file():
            try:
                content = rpt.read_text(encoding="utf-8", errors="replace")
            except Exception as e:
                print(f"  [WARN] Could not read {rpt}: {e}")
                content = ""
            reports.append({
                "mode": mode,
                "mode_id": mode_id,
                "content": content,
            })

    return figures, tables, reports


def build_manifest(figures, tables, reports) -> dict:
    # Single O(F) pass extracts all 3 dimension sets simultaneously
    # (was 3x O(F) with separate unique_sorted calls per key)
    dim_keys = ("mode", "analysis_type", "gene_group")
    dims: dict[str, set[str]] = {k: set() for k in dim_keys}
    for f in figures:
        for k in dim_keys:
            dims[k].add(f[k])
    return {
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "total_figures": len(figures),
        "total_tables": len(tables),
        "total_reports": len(reports),
        "dimensions": {
            "modes":          sorted(dims["mode"]),
            "analysis_types": sorted(dims["analysis_type"]),
            "gene_groups":    sorted(dims["gene_group"]),
        },
        "figures": figures,
        "tables": tables,
        "reports": reports,
    }


# ── HTML template ────────────────────────────────────────────────────────────

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Concordance Analysis Viewer</title>
<style>
:root {
  --bg: #0d1117; --surface: #161b22; --surface2: #21262d;
  --border: #30363d; --accent: #238636; --accent2: #1f6feb;
  --text: #e6edf3; --muted: #8b949e; --warn: #d29922;
  --tag-bg: #1c2128; --hover: #2d333b; --radius: 6px;
}
* { box-sizing: border-box; margin: 0; padding: 0; }
body { background: var(--bg); color: var(--text); font-family: -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; font-size: 14px; min-height: 100vh; }

/* ── Header ── */
header { background: var(--surface); border-bottom: 1px solid var(--border); padding: 14px 20px; display: flex; align-items: center; gap: 16px; position: sticky; top: 0; z-index: 100; }
header h1 { font-size: 16px; font-weight: 600; color: var(--text); white-space: nowrap; }
.badge { background: var(--accent); color: #fff; font-size: 11px; padding: 2px 8px; border-radius: 20px; font-weight: 600; }
.badge-blue { background: var(--accent2); }
.spacer { flex: 1; }
#stats { font-size: 12px; color: var(--muted); }
.zoom-controls { display: flex; align-items: center; gap: 4px; }
.zoom-btn { background: var(--surface2); border: 1px solid var(--border); color: var(--text); width: 26px; height: 26px; border-radius: var(--radius); cursor: pointer; font-size: 14px; display: flex; align-items: center; justify-content: center; }
.zoom-btn:hover { background: var(--hover); }
#zoom-level { font-size: 11px; color: var(--muted); min-width: 36px; text-align: center; }
#refresh-btn { background: var(--surface2); border: 1px solid var(--border); color: var(--text); padding: 5px 12px; border-radius: var(--radius); cursor: pointer; font-size: 12px; }
#refresh-btn:hover { background: var(--hover); }

/* ── Layout ── */
.layout { display: flex; height: calc(100vh - 53px); }

/* ── Sidebar ── */
aside { width: 260px; min-width: 220px; background: var(--surface); border-right: 1px solid var(--border); overflow-y: auto; flex-shrink: 0; }
.sidebar-section { border-bottom: 1px solid var(--border); padding: 14px; }
.sidebar-section h3 { font-size: 11px; font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em; color: var(--muted); margin-bottom: 10px; }

/* ── Sidebar tabs ── */
.sidebar-tabs { display: flex; flex-direction: column; gap: 4px; }
.sidebar-tab { background: var(--surface2); border: 1px solid var(--border); color: var(--muted); padding: 8px 12px; border-radius: var(--radius); cursor: pointer; font-size: 13px; text-align: left; transition: all .15s; }
.sidebar-tab:hover { background: var(--hover); color: var(--text); }
.sidebar-tab.active { background: var(--accent2); border-color: var(--accent2); color: #fff; }

/* ── Chips ── */
.chip-group { display: flex; flex-wrap: wrap; gap: 5px; }
.chip { background: var(--tag-bg); border: 1px solid var(--border); color: var(--muted); padding: 3px 10px; border-radius: 20px; font-size: 12px; cursor: pointer; transition: all .15s; user-select: none; }
.chip.active { background: var(--accent2); border-color: var(--accent2); color: #fff; }
.chip:hover:not(.active) { background: var(--hover); color: var(--text); }
.reset-btn { width: 100%; background: var(--surface2); border: 1px solid var(--border); color: var(--muted); padding: 7px; border-radius: var(--radius); cursor: pointer; font-size: 12px; margin-top: 4px; transition: all .15s; }
.reset-btn:hover { background: var(--hover); color: var(--text); }

/* ── Nav list (reports / tables sidebar) ── */
.nav-list { display: flex; flex-direction: column; gap: 3px; }
.nav-item { background: var(--surface2); border: 1px solid var(--border); color: var(--muted); padding: 7px 10px; border-radius: var(--radius); cursor: pointer; font-size: 12px; transition: all .15s; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
.nav-item:hover { background: var(--hover); color: var(--text); }
.nav-item.active { background: var(--accent2); border-color: var(--accent2); color: #fff; }

/* ── Main content ── */
main { flex: 1; overflow: auto; padding: 20px; }

/* ── Figure gallery ── */
.mode-section { margin-bottom: 28px; }
.mode-title { font-size: 15px; font-weight: 600; color: var(--accent2); margin-bottom: 12px; padding-bottom: 6px; border-bottom: 1px solid var(--border); }
.fig-gallery { display: grid; grid-template-columns: repeat(auto-fill, minmax(280px, 1fr)); gap: 14px; }
.fig-card { background: var(--surface); border: 1px solid var(--border); border-radius: var(--radius); overflow: hidden; cursor: pointer; transition: transform .15s, box-shadow .15s; }
.fig-card:hover { transform: translateY(-2px); box-shadow: 0 4px 12px rgba(0,0,0,.4); }
.fig-card img { width: 100%; display: block; background: var(--surface2); }
.fig-card-meta { padding: 8px 10px; }
.fig-card-type { font-size: 12px; font-weight: 600; color: var(--text); margin-bottom: 2px; }
.fig-card-gg { font-size: 11px; color: var(--muted); }
.no-results { text-align: center; padding: 60px 20px; color: var(--muted); }
.no-results h2 { font-size: 18px; margin-bottom: 8px; }

/* ── Report viewer ── */
.report-wrap { max-width: 900px; }
.report-wrap h1 { font-size: 22px; font-weight: 700; margin: 24px 0 12px; color: var(--text); border-bottom: 1px solid var(--border); padding-bottom: 6px; }
.report-wrap h2 { font-size: 18px; font-weight: 600; margin: 20px 0 10px; color: var(--text); }
.report-wrap h3 { font-size: 15px; font-weight: 600; margin: 16px 0 8px; color: var(--text); }
.report-wrap h4 { font-size: 13px; font-weight: 600; margin: 12px 0 6px; color: var(--muted); }
.report-wrap p { margin: 6px 0; line-height: 1.6; color: var(--text); }
.report-wrap strong { color: var(--text); }
.report-wrap em { font-style: italic; color: var(--muted); }
.report-wrap code { background: var(--surface2); padding: 1px 5px; border-radius: 3px; font-size: 13px; }
.report-wrap hr { border: none; border-top: 1px solid var(--border); margin: 16px 0; }
.report-wrap ul { margin: 6px 0 6px 20px; line-height: 1.6; }
.report-wrap li { margin: 2px 0; color: var(--text); }
.report-wrap table { border-collapse: collapse; margin: 10px 0; width: 100%; }
.report-wrap th { background: var(--surface2); padding: 6px 10px; border: 1px solid var(--border); font-size: 12px; font-weight: 600; color: var(--muted); text-align: left; white-space: nowrap; }
.report-wrap td { padding: 5px 10px; border: 1px solid var(--border); font-size: 13px; }
.report-wrap .md-img-wrap { margin: 12px 0; text-align: center; }
.report-wrap .md-img { max-width: 100%; border-radius: var(--radius); cursor: pointer; transition: transform .15s; }
.report-wrap .md-img:hover { transform: scale(1.01); }

/* ── Data table viewer ── */
.tbl-section { margin-bottom: 24px; }
.tbl-header { background: var(--surface); border: 1px solid var(--border); border-radius: var(--radius) var(--radius) 0 0; padding: 10px 14px; cursor: pointer; display: flex; align-items: center; justify-content: space-between; transition: background .15s; }
.tbl-header:hover { background: var(--hover); }
.tbl-header h4 { font-size: 13px; font-weight: 600; color: var(--text); }
.tbl-header .tbl-toggle { color: var(--muted); font-size: 11px; }
.tbl-body { border: 1px solid var(--border); border-top: none; border-radius: 0 0 var(--radius) var(--radius); overflow-x: auto; max-height: 500px; overflow-y: auto; display: none; }
.tbl-body.open { display: block; }
.data-table { border-collapse: collapse; width: 100%; }
.data-table th { background: var(--surface2); padding: 6px 10px; border: 1px solid var(--border); font-size: 11px; font-weight: 600; color: var(--muted); text-align: left; white-space: nowrap; cursor: pointer; position: sticky; top: 0; z-index: 1; user-select: none; }
.data-table th:hover { background: var(--hover); color: var(--text); }
.data-table th .sort-arrow { margin-left: 4px; font-size: 9px; }
.data-table td { padding: 4px 10px; border: 1px solid var(--border); font-size: 12px; white-space: nowrap; }
.data-table tr:nth-child(even) td { background: var(--surface); }
.data-table tr:nth-child(odd) td { background: var(--bg); }
.data-table td.num-cell { text-align: right; font-variant-numeric: tabular-nums; }
.corr-cell { font-weight: 500; }

/* ── Modal ── */
#modal { display: none; position: fixed; inset: 0; z-index: 1000; background: rgba(0,0,0,.85); backdrop-filter: blur(4px); align-items: center; justify-content: center; }
#modal.open { display: flex; }
.modal-box { background: var(--surface); border: 1px solid var(--border); border-radius: 10px; max-width: 92vw; max-height: 92vh; display: flex; flex-direction: column; overflow: hidden; }
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
</style>
</head>
<body>

<header>
  <h1>Concordance Analysis</h1>
  <span class="badge" id="fig-count">0 figures</span>
  <span class="badge badge-blue" id="report-count">0 reports</span>
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
    <h3>View</h3>
    <div class="sidebar-tabs" id="sidebar-tabs"></div>
  </div>

  <!-- Shared filters (always visible) -->
  <div class="sidebar-section">
    <h3>Concordance Mode</h3>
    <div class="chip-group" id="filter-mode"></div>
  </div>
  <div class="sidebar-section">
    <h3>Gene Group</h3>
    <div class="chip-group" id="filter-gg"></div>
  </div>

  <!-- Figures-only filters -->
  <div id="fig-filters">
    <div class="sidebar-section">
      <h3>Analysis Type</h3>
      <div class="chip-group" id="filter-atype"></div>
    </div>
  </div>

  <div class="sidebar-section">
    <button class="reset-btn" onclick="resetFilters()">&#10006; Reset All Filters</button>
  </div>
</aside>

<main id="main-content">
  <div id="view-figures"></div>
  <div id="view-reports" style="display:none"></div>
  <div id="view-tables" style="display:none"></div>
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
// ── Data ──────────────────────────────────────────────────────────────────
const MANIFEST = __MANIFEST_JSON__;

// ── State ─────────────────────────────────────────────────────────────────
const state = {
  views: new Set(["figures"]),   // multi-select: any combo of figures/reports/tables
  modes_active:  new Set(),      // single-select: only one concordance mode at a time
  atypes_active: new Set(),
  gg_active:     new Set(),
  _zoom: 100,
};

// Modal
let modalImages = [];
let modalIdx = 0;

// ── Boot ──────────────────────────────────────────────────────────────────
document.addEventListener("DOMContentLoaded", () => {
  const d = MANIFEST.dimensions;

  // View tabs
  const tabs = [
    { id: "figures", label: "Figures (" + MANIFEST.total_figures + ")" },
    { id: "reports", label: "Reports (" + MANIFEST.total_reports + ")" },
    { id: "tables",  label: "Data Tables (" + MANIFEST.total_tables + ")" },
  ];
  const tabContainer = document.getElementById("sidebar-tabs");
  tabs.forEach(t => {
    const btn = document.createElement("button");
    btn.className = "sidebar-tab" + (t.id === "figures" ? " active" : "");
    btn.textContent = t.label;
    btn.onclick = () => toggleView(t.id);
    btn.dataset.view = t.id;
    tabContainer.appendChild(btn);
  });

  // Figure filter chips
  buildChipsSingle("filter-mode",  d.modes,    "modes_active");
  buildChips("filter-atype", d.analysis_types,  "atypes_active");
  buildChips("filter-gg",    d.gene_groups,     "gg_active",    prettyGeneGroup);

  // Badges
  document.getElementById("fig-count").textContent = MANIFEST.total_figures + " figures";
  document.getElementById("report-count").textContent = MANIFEST.total_reports + " reports";

  // Ctrl+scroll zoom
  document.getElementById("main-content").addEventListener("wheel", e => {
    if (e.ctrlKey) { e.preventDefault(); setZoom(state._zoom + (e.deltaY < 0 ? 10 : -10)); }
  }, { passive: false });

  renderAll();
});

// ── Helpers ───────────────────────────────────────────────────────────────
function prettyLabel(s) { return s ? s.replace(/_/g, " ") : ""; }
function prettyGeneGroup(s) {
  if (!s) return "";
  return s.replace(/_in_PRJ[A-Z0-9]+(_\w+)?$/, "").replace(/_/g, " ");
}

function renderAll() { renderFigures(); renderReports(); renderTables(); }

// Multi-select chips (OR within filter, AND across filters)
function buildChips(containerId, values, stateKey, labelFn) {
  const fmt = labelFn || prettyLabel;
  const el = document.getElementById(containerId);
  values.forEach(v => {
    const c = document.createElement("span");
    c.className = "chip";
    c.textContent = fmt(v);
    c.dataset.value = v;
    c.onclick = () => {
      if (state[stateKey].has(v)) { state[stateKey].delete(v); c.classList.remove("active"); }
      else { state[stateKey].add(v); c.classList.add("active"); }
      renderAll();
    };
    el.appendChild(c);
  });
}

// Single-select chips (only one active at a time, click again to deselect)
function buildChipsSingle(containerId, values, stateKey, labelFn) {
  const fmt = labelFn || prettyLabel;
  const el = document.getElementById(containerId);
  values.forEach(v => {
    const c = document.createElement("span");
    c.className = "chip";
    c.textContent = fmt(v);
    c.dataset.value = v;
    c.onclick = () => {
      if (state[stateKey].has(v)) {
        state[stateKey].clear(); c.classList.remove("active");
      } else {
        state[stateKey].clear();
        el.querySelectorAll(".chip").forEach(ch => ch.classList.remove("active"));
        state[stateKey].add(v); c.classList.add("active");
      }
      updateDependentChips();
      renderAll();
    };
    el.appendChild(c);
  });
}

// Build lookup: mode → Set of analysis types that exist for that mode
const ATYPES_BY_MODE = {};
(function() {
  MANIFEST.figures.forEach(f => {
    if (!ATYPES_BY_MODE[f.mode]) ATYPES_BY_MODE[f.mode] = new Set();
    ATYPES_BY_MODE[f.mode].add(f.analysis_type);
  });
})();

// Show/hide analysis type chips based on selected concordance mode
function updateDependentChips() {
  const el = document.getElementById("filter-atype");
  const selectedMode = state.modes_active.size === 1 ? [...state.modes_active][0] : null;
  const allowed = selectedMode ? (ATYPES_BY_MODE[selectedMode] || new Set()) : null;

  el.querySelectorAll(".chip").forEach(c => {
    const v = c.dataset.value;
    if (allowed && !allowed.has(v)) {
      c.style.display = "none";
      // Clear if was active but now hidden
      if (state.atypes_active.has(v)) { state.atypes_active.delete(v); c.classList.remove("active"); }
    } else {
      c.style.display = "";
    }
  });
}

function resetFilters() {
  state.modes_active  = new Set();
  state.atypes_active = new Set();
  state.gg_active     = new Set();
  ["filter-mode","filter-atype","filter-gg"].forEach(id => {
    document.getElementById(id).querySelectorAll(".chip").forEach(c => c.classList.remove("active"));
  });
  updateDependentChips();
  renderAll();
}

// ── View toggling (multi-select) ──────────────────────────────────────────
function toggleView(view) {
  if (state.views.has(view)) {
    if (state.views.size > 1) state.views.delete(view);  // keep at least one active
  } else {
    state.views.add(view);
  }
  document.querySelectorAll("#sidebar-tabs .sidebar-tab").forEach(b => {
    b.classList.toggle("active", state.views.has(b.dataset.view));
  });
  const has = id => state.views.has(id);
  document.getElementById("fig-filters").style.display = has("figures") ? "" : "none";
  document.getElementById("view-figures").style.display = has("figures") ? "" : "none";
  document.getElementById("view-reports").style.display = has("reports") ? "" : "none";
  document.getElementById("view-tables").style.display  = has("tables")  ? "" : "none";
  if (has("reports")) renderReports();
  if (has("tables"))  renderTables();
}

// ── Figures ───────────────────────────────────────────────────────────────
function filteredFigures() {
  return MANIFEST.figures.filter(f => {
    if (state.modes_active.size  > 0 && !state.modes_active.has(f.mode))           return false;
    if (state.atypes_active.size > 0 && !state.atypes_active.has(f.analysis_type)) return false;
    if (state.gg_active.size     > 0 && !state.gg_active.has(f.gene_group))        return false;
    return true;
  });
}

function renderFigures() {
  const figs = filteredFigures();
  const wrap = document.getElementById("view-figures");
  document.getElementById("stats").textContent = figs.length + " / " + MANIFEST.total_figures + " shown";

  if (figs.length === 0) {
    wrap.innerHTML = '<div class="no-results"><h2>No figures match the current filters</h2><p>Try adjusting or resetting the filters.</p></div>';
    return;
  }

  // Group by mode
  const byMode = {};
  figs.forEach(f => { if (!byMode[f.mode]) byMode[f.mode] = []; byMode[f.mode].push(f); });

  // O(n) Map build — avoids O(n²) .indexOf() scan inside render loop
  const figIndexMap = new Map(figs.map((f, i) => [f, i]));

  // O(F) array accumulation + single .join('') — avoids O(F²) worst-case string concat
  const _parts = [];
  for (const mode of Object.keys(byMode).sort()) {
    _parts.push('<div class="mode-section">');
    _parts.push('<div class="mode-title">' + mode + '</div>');
    _parts.push('<div class="fig-gallery">');
    byMode[mode].forEach((f, i) => {
      const globalIdx = figIndexMap.get(f);
      _parts.push('<div class="fig-card" onclick="openModal(' + globalIdx + ')">');
      _parts.push('<img src="' + f.path + '" alt="" loading="lazy">');
      _parts.push('<div class="fig-card-meta">');
      _parts.push('<div class="fig-card-type">' + f.analysis_type + '</div>');
      if (f.gene_group !== "\u2014") _parts.push('<div class="fig-card-gg">' + prettyGeneGroup(f.gene_group) + '</div>');
      _parts.push('</div></div>');
    });
    _parts.push('</div></div>');
  }
  const html = _parts.join('');

  wrap.innerHTML = html;
  if (state._zoom !== 100) {
    wrap.querySelectorAll(".fig-gallery").forEach(g => { g.style.zoom = state._zoom / 100; });
  }
}

// ── Reports ──────────────────────────────────────────────────────────────
function renderReports() {
  const wrap = document.getElementById("view-reports");
  if (MANIFEST.reports.length === 0) {
    wrap.innerHTML = '<div class="no-results"><h2>No reports available</h2></div>';
    return;
  }
  const reports = state.modes_active.size > 0
    ? MANIFEST.reports.filter(r => state.modes_active.has(r.mode))
    : MANIFEST.reports;
  if (reports.length === 0) {
    wrap.innerHTML = '<div class="no-results"><h2>No reports match the selected mode</h2></div>';
    return;
  }
  // O(R) array accumulation + single .join('')
  const _rParts = [];
  reports.forEach(r => {
    _rParts.push('<div class="mode-section"><div class="mode-title">' + r.mode + '</div>');
    _rParts.push('<div class="report-wrap">' + mdToHtml(r.content, r.mode_id) + '</div></div>');
  });
  wrap.innerHTML = _rParts.join('');
}

// ── Minimal Markdown → HTML ──────────────────────────────────────────────
function mdToHtml(md, modeId) {
  const lines = md.split("\n");
  // O(L) array accumulation + single .join('') — avoids O(L²) worst-case string concat
  const _p = [];
  let inTable = false, inList = false, tableRows = [];

  for (let i = 0; i < lines.length; i++) {
    const line = lines[i];
    const trimmed = line.trim();

    // Headers
    const hm = trimmed.match(/^(#{1,4})\s+(.+)/);
    if (hm) {
      if (inList) { _p.push("</ul>\n"); inList = false; }
      const lvl = hm[1].length;
      _p.push("<h" + lvl + ">" + mdInline(hm[2]) + "</h" + lvl + ">\n");
      continue;
    }

    // HR
    if (/^---+$/.test(trimmed)) { _p.push("<hr>\n"); continue; }

    // Table row
    if (trimmed.startsWith("|") && trimmed.endsWith("|")) {
      if (/^\|[\s:|-]+\|$/.test(trimmed)) continue; // separator
      if (!inTable) { inTable = true; tableRows = []; }
      tableRows.push(trimmed);
      const next = (i + 1 < lines.length) ? lines[i + 1].trim() : "";
      if (!next.startsWith("|")) {
        _p.push(mdTable(tableRows));
        inTable = false;
        tableRows = [];
      }
      continue;
    }

    // List
    if (/^[-*]\s/.test(trimmed)) {
      if (!inList) { _p.push("<ul>\n"); inList = true; }
      _p.push("<li>" + mdInline(trimmed.replace(/^[-*]\s+/, "")) + "</li>\n");
      const next = (i + 1 < lines.length) ? lines[i + 1].trim() : "";
      if (!/^[-*]\s/.test(next)) { _p.push("</ul>\n"); inList = false; }
      continue;
    }

    // Image
    const imgM = trimmed.match(/^!\[([^\]]*)\]\(([^)]+)\)/);
    if (imgM) {
      // Rewrite path: figures/X.png → modeId/figures/X.png
      let src = imgM[2];
      if (src.startsWith("figures/")) src = modeId + "/" + src;
      _p.push('<div class="md-img-wrap"><img src="' + src + '" alt="' + imgM[1] + '" class="md-img" onclick="openModalSingle(this.src)"></div>\n');
      continue;
    }

    // Empty
    if (trimmed === "") { continue; }

    // Paragraph
    if (inList) { _p.push("</ul>\n"); inList = false; }
    _p.push("<p>" + mdInline(trimmed) + "</p>\n");
  }
  if (inList) _p.push("</ul>\n");
  return _p.join('');
}

function mdInline(t) {
  return t
    .replace(/\*\*([^*]+)\*\*/g, "<strong>$1</strong>")
    .replace(/\*([^*]+)\*/g, "<em>$1</em>")
    .replace(/`([^`]+)`/g, "<code>$1</code>");
}

function mdTable(rows) {
  if (rows.length === 0) return "";
  const _p = ["<table>"];
  rows.forEach((row, ri) => {
    const cells = row.split("|").slice(1, -1).map(c => c.trim());
    const tag = ri === 0 ? "th" : "td";
    _p.push("<tr>" + cells.map(c => "<" + tag + ">" + mdInline(c) + "</" + tag + ">").join("") + "</tr>");
  });
  _p.push("</table>\n");
  return _p.join("");
}

// ── Data Tables ──────────────────────────────────────────────────────────
function renderTables() {
  const wrap = document.getElementById("view-tables");
  if (MANIFEST.tables.length === 0) {
    wrap.innerHTML = '<div class="no-results"><h2>No data tables available</h2></div>';
    return;
  }

  const tables = state.modes_active.size > 0
    ? MANIFEST.tables.filter(t => state.modes_active.has(t.mode))
    : MANIFEST.tables;
  if (tables.length === 0) {
    wrap.innerHTML = '<div class="no-results"><h2>No tables match the selected mode</h2></div>';
    return;
  }

  // Group by mode
  const byMode = {};
  tables.forEach(t => { if (!byMode[t.mode]) byMode[t.mode] = []; byMode[t.mode].push(t); });

  // O(T) array accumulation + single .join('')
  const _tParts = [];
  let tblIdx = 0;
  for (const mode of Object.keys(byMode).sort()) {
    _tParts.push('<div class="mode-section"><div class="mode-title">' + mode + '</div>');
    byMode[mode].forEach(t => {
      const id = "tbl-body-" + (tblIdx++);
      _tParts.push('<div class="tbl-section">');
      _tParts.push('<div class="tbl-header" onclick="toggleTable(\'' + id + '\',this)">');
      _tParts.push('<h4>' + prettyLabel(t.filename.replace(".csv", "")) + '</h4>');
      _tParts.push('<span class="tbl-toggle">&#9654;</span></div>');
      _tParts.push('<div class="tbl-body" id="' + id + '">' + csvToHtmlTable(t.content, t.filename) + '</div>');
      _tParts.push('</div>');
    });
    _tParts.push('</div>');
  }
  wrap.innerHTML = _tParts.join('');
}

function toggleTable(id, headerEl) {
  const body = document.getElementById(id);
  const isOpen = body.classList.toggle("open");
  headerEl.querySelector(".tbl-toggle").textContent = isOpen ? "\u25BC" : "\u25B6";
}

function csvToHtmlTable(csv, filename) {
  const lines = csv.trim().split("\n").filter(l => l.trim());
  if (lines.length === 0) return "<p>Empty table</p>";

  const isCorr = /correlation|spearman/i.test(filename);
  const rows = lines.map(l => parseCSVRow(l));

  // O(R×C) array accumulation + single .join('') — avoids O(R²×C²) worst-case string concat
  const _cp = ['<table class="data-table">'];
  rows.forEach((cells, ri) => {
    _cp.push("<tr>");
    cells.forEach((c, ci) => {
      if (ri === 0) {
        _cp.push('<th onclick="sortDataTable(this,' + ci + ')">' + c + ' <span class="sort-arrow"></span></th>');
      } else {
        const num = parseFloat(c);
        const isNum = !isNaN(num) && c.trim() !== "";
        let style = "";
        if (isCorr && isNum && num >= 0 && num <= 1 && ci > 0) {
          const g = Math.round(60 + num * 140);
          style = ' style="background:rgba(' + (255 - g) + ',' + g + ',100,0.15)" class="corr-cell num-cell"';
        } else if (isNum) {
          style = ' class="num-cell"';
        }
        _cp.push("<td" + style + ">" + c + "</td>");
      }
    });
    _cp.push("</tr>");
  });
  _cp.push("</table>");
  return _cp.join('');
}

function parseCSVRow(line) {
  // Handle quoted fields with commas
  const cells = [];
  let current = "", inQuotes = false;
  for (let i = 0; i < line.length; i++) {
    const ch = line[i];
    if (ch === '"') { inQuotes = !inQuotes; }
    else if (ch === ',' && !inQuotes) { cells.push(current.trim()); current = ""; }
    else { current += ch; }
  }
  cells.push(current.trim());
  return cells;
}

function sortDataTable(th, colIdx) {
  const table = th.closest("table");
  const headerRow = table.rows[0];
  const rows = Array.from(table.rows).slice(1);

  // Determine sort direction
  const currentDir = th.dataset.sortDir || "none";
  const newDir = currentDir === "asc" ? "desc" : "asc";

  // Clear all sort arrows
  headerRow.querySelectorAll(".sort-arrow").forEach(s => s.textContent = "");
  th.querySelector(".sort-arrow").textContent = newDir === "asc" ? " \u25B2" : " \u25BC";
  th.dataset.sortDir = newDir;

  // Pre-parse sort keys once O(n) — avoids O(n log n) parseFloat calls inside comparator
  const keys = rows.map(r => {
    const text = r.cells[colIdx] ? r.cells[colIdx].textContent.trim() : "";
    const num = parseFloat(text);
    return { row: r, text, num, isNum: !isNaN(num) && text !== "" };
  });

  keys.sort((a, b) => {
    if (a.isNum && b.isNum) return newDir === "asc" ? a.num - b.num : b.num - a.num;
    return newDir === "asc" ? a.text.localeCompare(b.text) : b.text.localeCompare(a.text);
  });

  // Batch DOM append via DocumentFragment — single reflow instead of N reflows
  const frag = document.createDocumentFragment();
  keys.forEach(k => frag.appendChild(k.row));
  table.appendChild(frag);
}

// ── Modal ─────────────────────────────────────────────────────────────────
function openModal(idx) {
  modalImages = filteredFigures();
  modalIdx = idx;
  showModalAt(modalIdx);
  document.getElementById("modal").classList.add("open");
  document.addEventListener("keydown", onModalKey);
}

function openModalSingle(src) {
  // For report images — single image modal
  modalImages = [{ path: src, mode: "", analysis_type: "", gene_group: "", filename: src.split("/").pop() }];
  modalIdx = 0;
  showModalAt(0);
  document.getElementById("modal").classList.add("open");
  document.addEventListener("keydown", onModalKey);
}

function showModalAt(idx) {
  const f = modalImages[idx];
  document.getElementById("modal-img").src = f.path;
  document.getElementById("modal-download").href = f.path;
  document.getElementById("modal-download").download = f.filename || f.path.split("/").pop();
  document.getElementById("modal-title").textContent = (f.filename || f.path.split("/").pop()).replace(/_/g, " ").replace(".png", "");
  document.getElementById("modal-tags").innerHTML =
    ["mode","analysis_type","gene_group"].filter(k => f[k] && f[k] !== "\u2014")
    .map(k => '<span class="modal-tag">' + prettyLabel(k) + ': ' + prettyLabel(f[k]) + '</span>').join("");
  document.getElementById("modal-pos").textContent = (idx + 1) + " / " + modalImages.length;
  document.getElementById("modal-prev").disabled = idx === 0;
  document.getElementById("modal-next").disabled = idx === modalImages.length - 1;
}

function navModal(dir) {
  modalIdx = Math.max(0, Math.min(modalImages.length - 1, modalIdx + dir));
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

// ── Zoom ──────────────────────────────────────────────────────────────────
function setZoom(level) {
  state._zoom = Math.max(25, Math.min(300, level));
  document.getElementById("zoom-level").textContent = state._zoom + "%";
  // Apply to active view content
  const z = state._zoom / 100;
  document.querySelectorAll(".fig-gallery").forEach(g => { g.style.zoom = z; });
  document.querySelectorAll(".report-wrap").forEach(g => { g.style.zoom = z; });
  document.querySelectorAll(".data-table").forEach(g => { g.style.zoom = z; });
}
</script>
</body>
</html>
"""


# ── HTML generation ──────────────────────────────────────────────────────────

def generate_html(manifest: dict) -> str:
    json_str = json.dumps(manifest, separators=(",", ":"))
    # Escape </ sequences to prevent </script> breakout in inline JSON
    json_str = json_str.replace("</", r"<\/")
    return HTML_TEMPLATE.replace("__MANIFEST_JSON__", json_str)


# ── Entry point ──────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Generate interactive HTML viewer for concordance analysis results.")
    parser.add_argument("concordance_dir",
                        help="Path to II_RESULTS/4_CONCORDANCE_ANALYSIS directory")
    parser.add_argument("--output", default=None,
                        help="Output HTML path (default: <dir>/concordance_viewer.html)")
    parser.add_argument("--manifest-only", action="store_true",
                        help="Only write concordance_manifest.json, skip HTML")
    args = parser.parse_args()

    base_dir = Path(args.concordance_dir).resolve()
    if not base_dir.is_dir():
        print(f"ERROR: Directory not found: {base_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Scanning {base_dir} \u2026", file=sys.stderr)
    figures, tables, reports = scan_concordance(base_dir)

    if not figures and not tables and not reports:
        print("WARNING: No concordance outputs found.", file=sys.stderr)

    manifest = build_manifest(figures, tables, reports)
    print(f"Found {len(figures)} figures, {len(tables)} tables, "
          f"{len(reports)} reports across {len(manifest['dimensions']['modes'])} modes.",
          file=sys.stderr)

    # Write manifest JSON
    manifest_path = base_dir / "concordance_manifest.json"
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    print(f"Manifest written: {manifest_path}", file=sys.stderr)

    if args.manifest_only:
        return

    # Write HTML
    output_path = Path(args.output).resolve() if args.output else base_dir / "concordance_viewer.html"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    html = generate_html(manifest)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"HTML viewer written: {output_path}", file=sys.stderr)
    print(str(output_path))


if __name__ == "__main__":
    main()
