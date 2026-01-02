import json
import os

def check_results():
    try:
        with open('channel_classification_results.json', 'r') as f:
            data = json.load(f)
            
        print(f"Loaded {len(data)} results.")
        
        classifications = {}
        
        for entry in data:
            ch = entry['channel']
            cls = entry['classification']
            
            if ch not in classifications:
                classifications[ch] = []
            classifications[ch].append(cls)
            
        print("\nSummary by Channel:")
        for ch, clss in classifications.items():
            finite_count = sum(1 for c in clss if "Finite" in c)
            pole_count = sum(1 for c in clss if "Pole" in c)
            zero_count = sum(1 for c in clss if "Zero" in c)
            unknown_count = sum(1 for c in clss if "Unknown" in c)
            
            status = "MIXED"
            if finite_count == len(clss): status = "FINITE"
            elif pole_count == len(clss): status = "POLE"
            
            print(f"  {ch}: {status} (Finite: {finite_count}, Pole: {pole_count}, Zero: {zero_count}, Unknown: {unknown_count})")
            
    except FileNotFoundError:
        print("Results file not found.")

if __name__ == "__main__":
    check_results()







