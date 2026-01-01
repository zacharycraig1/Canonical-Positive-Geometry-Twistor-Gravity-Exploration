import subprocess
import sys
import os

def run_sage_script(script_path):
    """Runs a sage script and returns True if successful."""
    print(f"Running {script_path}...")
    try:
        # Use 'sage' command if available, otherwise try to run with python (assuming inside sage shell)
        # On Windows Cursor environment, we often run scripts with python if sage is in path or python IS sage
        # Let's try executing with sys.executable first if it's sage, or 'sage' command.
        
        # Checking if we are in a sage environment or if we need to call sage
        cmd = ['sage', script_path]
        
        # If 'sage' is not in path, and we are running inside sage's python, we might just call python
        # But 'sage' command handles imports better usually.
        # Fallback to python if sage command fails? 
        # In this specific environment, user likely runs via `sage -python` or similar. 
        # But let's assume `sage` is available as a command or we use the current python executable 
        # if it is the sage python.
        
        # Adjust for Windows/specific env if needed. 
        # For now, let's try running it as a subprocess with the same interpreter 
        # IF the file ends in .py, otherwise use 'sage'.
        # The files are .sage, so they need pre-parsing or running with sage.
        
        # However, the user provided files like `src/scripts/physics_pullback_n4.sage`
        # BUT there are also `.sage.py` versions generated. 
        # We should prefer running the .sage file with `sage`.
        
        if sys.platform == 'win32':
             # On Windows, 'sage' might be a batch file or shell script, or not in path directly.
             # We'll try using the local 'sage.ps1' if it exists.
             sage_ps1 = os.path.join(os.getcwd(), 'sage.ps1')
             if os.path.exists(sage_ps1):
                 # Use powershell to run sage.ps1
                 # We pass the arguments to sage.ps1 which passes them to sage in docker
                 cmd = ['powershell', '-ExecutionPolicy', 'Bypass', '-File', sage_ps1, script_path]
             else:
                 cmd = ['sage', script_path]
        
        result = subprocess.run(cmd, check=True, text=True, capture_output=False)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running {script_path}: {e}")
        return False
    except FileNotFoundError:
        print("Command 'sage' not found. Trying 'python' on the .sage.py equivalent if it exists...")
        py_path = script_path + ".py"
        if os.path.exists(py_path):
             try:
                 subprocess.run([sys.executable, py_path], check=True)
                 return True
             except subprocess.CalledProcessError as e:
                 print(f"Error running {py_path}: {e}")
                 return False
        else:
            print(f"Could not find 'sage' command or {py_path}")
            return False

def main():
    scripts = [
        "src/scripts/physics_pullback_n4.sage",
        "src/scripts/physics_pullback_n5.sage",
        "src/scripts/physics_pullback_n6.sage"
    ]
    
    failures = []
    
    for script in scripts:
        if not run_sage_script(script):
            failures.append(script)
            
    if failures:
        print("\nCORE TESTS FAILED:")
        for f in failures:
            print(f"- {f}")
        sys.exit(1)
    else:
        print("\nALL CORE TESTS PASSED. Certificate Valid.")
        sys.exit(0)

if __name__ == "__main__":
    main()

