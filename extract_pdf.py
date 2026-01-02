#!/usr/bin/env python3
"""Extract text from factorization PDF"""
import sys

try:
    import PyPDF2
    with open('factorization_dim8_constraints_report.pdf', 'rb') as pdf_file:
        reader = PyPDF2.PdfReader(pdf_file)
        text = ''
        for page in reader.pages:
            text += page.extract_text() + '\n'
        print(text)
except ImportError:
    try:
        import pdfplumber
        with pdfplumber.open('factorization_dim8_constraints_report.pdf') as pdf:
            text = ''
            for page in pdf.pages:
                text += page.extract_text() + '\n'
            print(text)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
except Exception as e:
    print(f"Error: {e}", file=sys.stderr)
    sys.exit(1)












