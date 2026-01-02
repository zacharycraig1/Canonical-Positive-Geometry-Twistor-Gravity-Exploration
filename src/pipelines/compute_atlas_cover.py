import json
import sys
import os

def solve_set_cover():
    input_path = "RESULTS/atlas_sweep_all_pairs_n6.json"
    if not os.path.exists(input_path):
        print(f"Error: {input_path} not found.")
        return

    print(f"Loading data from {input_path}...")
    with open(input_path, 'r') as f:
        data = json.load(f)

    # Organize data: for each root set, what channels does it cover?
    # And for each channel, is it covered?
    
    # Map: tuple(roots) -> set(channels)
    roots_to_channels = {}
    all_channels = set()
    
    for entry in data:
        roots = tuple(sorted(entry['roots']))
        channel = entry['channel']
        match = entry['match']
        
        all_channels.add(channel)
        
        if roots not in roots_to_channels:
            roots_to_channels[roots] = set()
            
        if match:
            roots_to_channels[roots].add(channel)
            
    print(f"Total Channels: {len(all_channels)}")
    print(f"Total Root Sets: {len(roots_to_channels)}")
    
    # Check if a solution exists
    union_all = set()
    for s in roots_to_channels.values():
        union_all.update(s)
        
    missing = all_channels - union_all
    if missing:
        print(f"CRITICAL WARNING: The following channels are NOT covered by ANY chart: {missing}")
        # Identify why? 
        # For now, continue to cover what we can
    else:
        print("All channels are coverable.")

    # Greedy Set Cover
    covered_channels = set()
    cover_solution = []
    
    uncovered = all_channels.copy()
    
    # Heuristic: Pick chart that covers most uncovered channels
    while uncovered:
        best_root = None
        best_count = -1
        best_new_cover = set()
        
        candidates = []
        
        for roots, channels in roots_to_channels.items():
            new_cover = channels.intersection(uncovered)
            count = len(new_cover)
            candidates.append((count, roots, new_cover))
            
            if count > best_count:
                best_count = count
                best_root = roots
                best_new_cover = new_cover
        
        if best_count <= 0:
            print("Cannot cover remaining channels (should not happen if check passed).")
            break
            
        print(f"Selected Chart {best_root} covering {best_count} new channels: {sorted(list(best_new_cover))}")
        
        cover_solution.append({
            "roots": list(best_root),
            "covers": list(best_new_cover) # This is just what it 'newly' covers in greedy step
        })
        
        covered_channels.update(best_new_cover)
        uncovered -= best_new_cover
        
    print(f"Found cover with {len(cover_solution)} charts.")
    
    # Build complete coverage map (each channel covered by which chart in solution?)
    # A channel might be covered by multiple in the solution.
    # We assign it to the first one in the list (or keep all?)
    # For the sum, we might want 'partitions of unity' or just 'sum over all'? 
    # The hypothesis is Sum_{R} c_R Omega_R.
    # If we sum over just the cover, we hope it works. 
    
    # Output
    out_path = "RESULTS/atlas_cover_n6.json"
    with open(out_path, 'w') as f:
        json.dump(cover_solution, f, indent=2)
    print(f"Saved cover to {out_path}")

if __name__ == "__main__":
    solve_set_cover()



