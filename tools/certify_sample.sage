import json
import datetime
from sage.all import *

def json_converter(obj):
    if hasattr(obj, 'list'): # Vector
        return obj.list()
    if isinstance(obj, (Integer, int)):
        return int(obj)
    if isinstance(obj, (Rational, float)):
        return float(obj) # Warning: precision loss for Rationals
    if hasattr(obj, 'str'):
        return str(obj)
    return str(obj)

def serialize_sample(data):
    """
    Serialize sample data to a JSON-compatible dict.
    Converts Sage types to Python types.
    For exact rational numbers, we might want to store them as strings "p/q".
    """
    out = {}
    for k, v in data.items():
        if isinstance(v, list):
            out[k] = [serialize_value(x) for x in v]
        elif isinstance(v, dict):
            out[k] = {str(sk): serialize_value(sv) for sk, sv in v.items()}
        else:
            out[k] = serialize_value(v)
    return out

def serialize_value(v):
    if hasattr(v, 'list'): # Vector
        return [str(x) for x in v]
    if v in QQ:
        return str(v)
    if isinstance(v, (Integer, int)):
        return int(v)
    return str(v)

def save_certificate(data, filename_prefix="cert"):
    """
    Save certificate to results/certificates/
    """
    import os
    os.makedirs("results/certificates", exist_ok=True)
    
    timestamp = datetime.datetime.now().isoformat().replace(":", "-")
    filename = f"results/certificates/{filename_prefix}_{timestamp}.json"
    
    serialized = serialize_sample(data)
    
    with open(filename, 'w') as f:
        json.dump(serialized, f, indent=2)
    
    return filename






