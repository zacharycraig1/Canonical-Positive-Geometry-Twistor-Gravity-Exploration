#!/usr/bin/env sage
"""Extract readable text from PDF"""
import re

# Read PDF as binary
with open('factorization_dim8_constraints_report.pdf', 'rb') as f:
    content = f.read()

# Look for text patterns in PDF
# PDF text streams are often between stream and endstream markers
# Try to find readable ASCII text
text_patterns = []

# Method 1: Look for text between parentheses (PDF text objects)
paren_text = re.findall(rb'\(([^)]{10,200})\)', content)
for match in paren_text[:50]:  # First 50 matches
    try:
        decoded = match.decode('ascii', errors='ignore')
        if len(decoded.strip()) > 5 and any(c.isalpha() for c in decoded):
            text_patterns.append(decoded)
    except:
        pass

# Method 2: Look for readable ASCII sequences
ascii_pattern = re.findall(rb'[ -~]{20,100}', content)
for match in ascii_pattern[:30]:
    try:
        text = match.decode('ascii', errors='ignore')
        if any(c.isalpha() for c in text) and 'factorization' in text.lower() or 'dim' in text.lower() or 'constraint' in text.lower():
            text_patterns.append(text)
    except:
        pass

# Print found text
print("="*70)
print("EXTRACTED TEXT FROM PDF:")
print("="*70)
for i, text in enumerate(set(text_patterns), 1):
    if len(text.strip()) > 10:
        print(f"\n{i}. {text[:200]}")








