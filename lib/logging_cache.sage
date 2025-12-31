# Lightweight logging and caching utilities for Sage-based pipelines
from sageall import *
import os, json, hashlib, time, pathlib

def ensure_dir(path):
    if not path:
        return
    os.makedirs(path, exist_ok=True)

def _now_timestamp():
    return time.strftime("%Y-%m-%d_%H%M%S")

def new_run_dir(base_dir="runs", seed=None, mode=None):
    run_id = _now_timestamp()
    if seed is not None:
        run_id += f"_seed{int(seed)}"
    if mode is not None:
        run_id += f"_mode{str(mode)}"
    run_dir = os.path.join(os.getcwd(), base_dir, run_id)
    ensure_dir(run_dir)
    return run_dir

def sha256_of_bytes(b):
    return hashlib.sha256(b).hexdigest()

def sha256_of_jsonable(obj):
    blob = json.dumps(obj, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return sha256_of_bytes(blob)

def write_json(path, obj):
    ensure_dir(os.path.dirname(path) or ".")
    with open(path, "w") as f:
        json.dump(obj, f, indent=2)

def read_json(path, default=None):
    try:
        with open(path, "r") as f:
            return json.load(f)
    except Exception:
        return default

def save_sobj(obj, path):
    ensure_dir(os.path.dirname(path) or ".")
    save(obj, path)

def load_sobj(path):
    return load(path)

def snapshot_config(run_dir, config_dict):
    write_json(os.path.join(run_dir, "config.json"), config_dict)

def log_line(message, run_dir=None):
    print(message)
    if run_dir:
        ensure_dir(run_dir)
        with open(os.path.join(run_dir, "search_live.log"), "a") as f:
            f.write(message.rstrip() + "\n")








