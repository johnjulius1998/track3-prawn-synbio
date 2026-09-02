#!/usr/bin/env python3
"""
monitor_dashboard.py — Live Progress Dashboard for Track 3 Pipeline
=====================================================================
Serves a self-refreshing HTML dashboard on http://localhost:8080
showing batch processing progress, sample stats, and system resources.

USAGE: python src/utils/monitor_dashboard.py [--port 8080]
      Then open http://localhost:8080 in a browser.
"""

import argparse
import json
import os
import time
from pathlib import Path
from http.server import HTTPServer, BaseHTTPRequestHandler

PROJECT_ROOT = Path(__file__).resolve().parents[3]


def get_pipeline_status():
    """Read current pipeline state and return as dict."""
    status = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "samples": {},
        "batch_log": [],
        "system": {},
        "microbial": {},
        "edges": {},
    }

    # ---- Batch log ----
    log_path = PROJECT_ROOT / "results/logs/batch_process.log"
    if log_path.exists():
        with open(log_path) as f:
            lines = f.readlines()
            status["batch_log"] = [l.strip() for l in lines[-30:]]  # Last 30 lines
            # Parse current sample
            for line in reversed(lines):
                if "Processing" in line and "===" in line:
                    status["current_action"] = line.strip()
                    break
                elif "Downloading" in line:
                    status["current_action"] = line.strip()
                    break
    else:
        status["batch_log"] = ["Batch not yet started."]
        status["current_action"] = "Idle"

    # ---- Quantified samples ----
    quant_dir = PROJECT_ROOT / "data/interim/host_counts"
    total_expected = 20
    completed = []
    if quant_dir.exists():
        for d in sorted(quant_dir.iterdir()):
            if d.is_dir() and (d / "quant.sf").exists():
                srr = d.name
                # Get sample metadata
                meta_path = PROJECT_ROOT / "data/raw/sra/PRJNA875278/metadata.tsv"
                meta_info = {}
                if meta_path.exists():
                    with open(meta_path) as f:
                        header = f.readline().strip().split("\t")
                        for line in f:
                            parts = line.strip().split("\t")
                            if parts[0] == srr:
                                meta_info = dict(zip(header, parts))
                                break
                # Get gene count from quant.sf
                try:
                    import pandas as pd
                    qf = pd.read_csv(d / "quant.sf", sep="\t", nrows=5)
                    n_genes = len(pd.read_csv(d / "quant.sf", sep="\t", usecols=["TPM"]))
                    n_expressed = int((pd.read_csv(d / "quant.sf", sep="\t", usecols=["TPM"])["TPM"] > 0).sum())
                except Exception:
                    n_genes = "?"
                    n_expressed = "?"

                completed.append({
                    "srr": srr,
                    "tissue": meta_info.get("tissue", "?"),
                    "sex": meta_info.get("sex", "?"),
                    "weight_gain": meta_info.get("weight_gain", "?"),
                    "genes": n_genes,
                    "expressed": n_expressed,
                })

    status["samples"] = {
        "completed": len(completed),
        "total": total_expected,
        "percent": round(100 * len(completed) / total_expected, 1),
        "list": completed,
    }

    # ---- System resources ----
    try:
        import subprocess
        mem = subprocess.run(["free", "-h"], capture_output=True, text=True).stdout
        disk = subprocess.run(["df", "-h", str(PROJECT_ROOT)], capture_output=True, text=True).stdout
        status["system"] = {
            "memory": mem.strip().split("\n")[1] if mem else "?",
            "disk": disk.strip().split("\n")[-1] if disk else "?",
        }
    except Exception:
        status["system"] = {"memory": "?", "disk": "?"}

    # ---- Microbial shortlist stats ----
    shortlist_path = PROJECT_ROOT / "results/shortlist/microbial_taxa.csv"
    if shortlist_path.exists():
        import pandas as pd
        df = pd.read_csv(shortlist_path)
        status["microbial"] = {
            "total_taxa": len(df),
            "jumper_associated": int(df["direction_of_effect"].str.contains("Jumper", na=False).sum()),
            "laggard_associated": int(df["direction_of_effect"].str.contains("Laggard", na=False).sum()),
            "contaminants": int((df.get("contaminant_risk") == "HIGH").sum()),
        }

    # ---- Network edges ----
    edges_path = PROJECT_ROOT / "results/shortlist/network_edges.csv"
    if edges_path.exists():
        import pandas as pd
        df = pd.read_csv(edges_path)
        status["edges"] = {
            "total": len(df),
            "function_overlap": int((df["edge_basis"] == "predicted_function_overlap").sum()),
            "phenotype_concordance": int((df["edge_basis"] == "phenotype_concordance").sum()),
        }

    return status


HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Track 3 Pipeline Monitor</title>
<meta http-equiv="refresh" content="30">
<style>
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{ font-family: 'Segoe UI', system-ui, sans-serif; background: #0d1117; color: #c9d1d9; padding: 20px; }}
h1 {{ color: #58a6ff; margin-bottom: 10px; }}
h2 {{ color: #8b949e; font-size: 14px; margin-bottom: 15px; font-weight: normal; }}
.grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 15px; margin-bottom: 20px; }}
.card {{ background: #161b22; border: 1px solid #30363d; border-radius: 8px; padding: 15px; }}
.card h3 {{ color: #58a6ff; margin-bottom: 10px; font-size: 14px; text-transform: uppercase; letter-spacing: 0.5px; }}
.stat {{ display: flex; justify-content: space-between; padding: 4px 0; border-bottom: 1px solid #21262d; font-size: 13px; }}
.stat .val {{ color: #f0f6fc; font-weight: bold; }}
.bar-bg {{ background: #21262d; border-radius: 4px; height: 20px; margin: 8px 0; overflow: hidden; }}
.bar-fg {{ background: linear-gradient(90deg, #238636, #3fb950); height: 100%; border-radius: 4px; transition: width 0.5s; }}
.log {{ background: #0d1117; border: 1px solid #30363d; border-radius: 4px; padding: 8px; font-family: 'Cascadia Code', 'Fira Code', monospace; font-size: 12px; max-height: 300px; overflow-y: auto; white-space: pre-wrap; color: #8b949e; }}
.sample-row {{ display: flex; justify-content: space-between; padding: 2px 0; font-size: 12px; border-bottom: 1px solid #21262d; }}
.sample-row .srr {{ color: #58a6ff; }}
.male {{ color: #58a6ff; }} .female {{ color: #f778ba; }}
.warn {{ color: #d2991d; }} .good {{ color: #3fb950; }} .bad {{ color: #f85149; }}
.footer {{ text-align: center; color: #484f58; font-size: 11px; margin-top: 20px; }}
@keyframes pulse {{ 0%,100% {{opacity:1}} 50% {{opacity:0.5}} }}
.live {{ animation: pulse 2s infinite; color: #3fb950; }}
</style>
</head>
<body>
<h1>🦐 Track 3 Host-Microbe Pipeline</h1>
<h2>Last update: {timestamp} <span class="live">● LIVE</span></h2>

<div class="grid">
<div class="card">
<h3>📊 Host RNA-seq Progress</h3>
<div class="stat"><span>Samples quantified</span><span class="val">{samples_completed} / {samples_total}</span></div>
<div class="bar-bg"><div class="bar-fg" style="width:{samples_percent}%"></div></div>
<div class="stat"><span>Current action</span><span class="val">{current_action}</span></div>
</div>

<div class="card">
<h3>💻 System</h3>
<div class="stat"><span>Memory</span><span class="val">{memory}</span></div>
<div class="stat"><span>Disk</span><span class="val">{disk}</span></div>
</div>

<div class="card">
<h3>🦠 Microbial Layer</h3>
<div class="stat"><span>Taxa shortlisted</span><span class="val">{microbial_taxa}</span></div>
<div class="stat"><span>Jumper-associated</span><span class="val">{microbial_jumper}</span></div>
<div class="stat"><span>Laggard-associated</span><span class="val">{microbial_laggard}</span></div>
<div class="stat"><span>Contaminants flagged</span><span class="val">{microbial_contam}</span></div>
</div>

<div class="card">
<h3>🔗 Integration Edges</h3>
<div class="stat"><span>Total edges</span><span class="val">{edges_total}</span></div>
<div class="stat"><span>Function overlap</span><span class="val">{edges_func}</span></div>
<div class="stat"><span>Phenotype concordance</span><span class="val">{edges_pheno}</span></div>
</div>
</div>

<h3 style="margin-bottom:8px;">📋 Quantified Samples</h3>
<div class="card" style="margin-bottom:15px;">
<div class="sample-row" style="font-weight:bold;color:#8b949e;">
<span>SRR ID</span><span>Tissue</span><span>Sex</span><span>Weight Gain</span><span>Genes</span>
</div>
{sample_rows}
</div>

<h3 style="margin-bottom:8px;">📜 Batch Log (last 30 lines)</h3>
<div class="log">{batch_log}</div>

<div class="footer">Track 3 Host-Microbe Integration &bull; Auto-refresh every 30s &bull; {timestamp}</div>
</body>
</html>"""


class DashboardHandler(BaseHTTPRequestHandler):
    def do_GET(self):
        if self.path == "/api/status":
            # JSON API endpoint
            status = get_pipeline_status()
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.end_headers()
            self.wfile.write(json.dumps(status, indent=2).encode())
        else:
            # HTML dashboard
            status = get_pipeline_status()
            s = status
            html = HTML_TEMPLATE.format(
                timestamp=s["timestamp"],
                samples_completed=s["samples"]["completed"],
                samples_total=s["samples"]["total"],
                samples_percent=s["samples"]["percent"],
                current_action=s.get("current_action", "Idle")[:80],
                memory=s["system"].get("memory", "?"),
                disk=s["system"].get("disk", "?"),
                microbial_taxa=s["microbial"].get("total_taxa", "?"),
                microbial_jumper=s["microbial"].get("jumper_associated", "?"),
                microbial_laggard=s["microbial"].get("laggard_associated", "?"),
                microbial_contam=s["microbial"].get("contaminants", "?"),
                edges_total=s["edges"].get("total", "?"),
                edges_func=s["edges"].get("function_overlap", "?"),
                edges_pheno=s["edges"].get("phenotype_concordance", "?"),
                sample_rows="\n".join(
                    f'<div class="sample-row"><span class="srr">{r["srr"]}</span>'
                    f'<span>{r["tissue"][:16]}</span>'
                    f'<span class="{r["sex"]}">{r["sex"]}</span>'
                    f'<span>{r["weight_gain"]}</span>'
                    f'<span>{r["expressed"]}</span></div>'
                    for r in s["samples"]["list"]
                ) if s["samples"]["list"] else '<div class="sample-row"><span colspan="5">No samples quantified yet</span></div>',
                batch_log="\n".join(s["batch_log"]) if s["batch_log"] else "No log data",
            )
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.end_headers()
            self.wfile.write(html.encode())

    def log_message(self, format, *args):
        pass  # Suppress console log spam


def main():
    parser = argparse.ArgumentParser(description="Track 3 Pipeline Monitor Dashboard")
    parser.add_argument("--port", type=int, default=8080, help="Port to serve on (default: 8080)")
    args = parser.parse_args()

    print(f"\n{'='*60}")
    print(f"  🦐 Track 3 Pipeline Monitor")
    print(f"  Open: http://localhost:{args.port}")
    print(f"  API:  http://localhost:{args.port}/api/status")
    print(f"  Press Ctrl+C to stop")
    print(f"{'='*60}\n")

    server = HTTPServer(("0.0.0.0", args.port), DashboardHandler)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down...")
        server.shutdown()


if __name__ == "__main__":
    main()
