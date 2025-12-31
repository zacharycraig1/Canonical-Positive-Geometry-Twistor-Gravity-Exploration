#!/usr/bin/env sage
"""Read PDF using Sage's capabilities"""
import subprocess
import sys

try:
    # Try using pdftotext if available
    result = subprocess.run(['pdftotext', 'factorization_dim8_constraints_report.pdf', '-'], 
                          capture_output=True, text=True, timeout=10)
    if result.returncode == 0:
        print(result.stdout)
    else:
        print("pdftotext not available or failed", file=sys.stderr)
        # Try alternative: read as binary and look for text patterns
        with open('factorization_dim8_constraints_report.pdf', 'rb') as f:
            content = f.read()
            # Look for text streams (basic extraction)
            import re
            # PDF text is often in streams between stream and endstream
            # This is a very basic extraction
            text_matches = re.findall(rb'stream\s+(.*?)\s+endstream', content, re.DOTALL)
            for match in text_matches[:10]:  # First 10 streams
                try:
                    # Try to decode as text
                    decoded = match.decode('latin-1', errors='ignore')
                    # Filter for readable text
                    readable = ''.join(c if 32 <= ord(c) < 127 or c in '\n\t' else ' ' for c in decoded)
                    if len(readable.strip()) > 20:  # Only print substantial text
                        print(readable[:500])  # First 500 chars
                except:
                    pass
except Exception as e:
    print(f"Error: {e}", file=sys.stderr)
    sys.exit(1)



