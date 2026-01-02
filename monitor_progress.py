#!/usr/bin/env python3
"""Monitor search progress by watching runs directory."""
import os
import time
import json
from pathlib import Path
from datetime import datetime

def find_latest_run():
    """Find most recent run directory."""
    runs_dir = Path("runs")
    if not runs_dir.exists():
        return None
    run_dirs = [d for d in runs_dir.iterdir() if d.is_dir()]
    if not run_dirs:
        return None
    return max(run_dirs, key=lambda x: x.stat().st_mtime)

def check_progress():
    """Check current progress."""
    run_dir = find_latest_run()
    if not run_dir:
        return {"status": "no_runs", "message": "No runs started yet"}
    
    checkpoint = run_dir / "intersection_checkpoint.sobj"
    report = run_dir / "hit_report_intersection.json"
    
    result = {
        "run_dir": str(run_dir),
        "has_checkpoint": checkpoint.exists(),
        "has_report": report.exists(),
    }
    
    if report.exists():
        try:
            with open(report, 'r') as f:
                data = json.load(f)
            result["final_dim"] = data.get("final_dim")
            result["invariant_mode"] = data.get("invariant_mode")
            result["boundaries"] = len(data.get("boundaries", []))
            result["success"] = data.get("final_dim") == 1
            result["message"] = f"Final dimension: {result['final_dim']}"
        except Exception as e:
            result["error"] = str(e)
    elif checkpoint.exists():
        result["message"] = "Search in progress (checkpoint exists)"
    else:
        result["message"] = "Search starting..."
    
    return result

if __name__ == '__main__':
    progress = check_progress()
    print(json.dumps(progress, indent=2))













