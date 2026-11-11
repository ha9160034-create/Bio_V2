import os
import sys
import time
import subprocess
import traceback
import urllib.request
from pathlib import Path

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")

IMAGE_NAME = "biobio-app:latest"
CONTAINER_NAME = "biobio_desktop_instance"
PORT = 8501
URL = f"http://localhost:{PORT}"

def main():
    print("==================================================")
    print("   BIOBIO - Molecular Docking & MD Suite")
    print("==================================================")

    chk = subprocess.run(["docker", "info"], capture_output=True, text=True, encoding="utf-8", errors="replace")
    if chk.returncode != 0:
        print("\n[!] Please open Docker Desktop first and wait until it is running.")
        input("\nPress Enter to close...")
        return

    subprocess.run(["docker", "rm", "-f", CONTAINER_NAME], capture_output=True)

    base_dir = Path(sys.executable).parent if getattr(sys, "frozen", False) else Path(__file__).parent
    docking_dir = (base_dir / "docking_runs").resolve()
    md_dir = (base_dir / "md_runs").resolve()
    docking_dir.mkdir(exist_ok=True)
    md_dir.mkdir(exist_ok=True)

    print("\n[*] Starting BIOBIO Engine (GROMACS + Vina + ACPYPE)...")
    run_cmd = [
        "docker", "run", "-d", "--rm",
        "--name", CONTAINER_NAME,
        "-p", f"{PORT}:8501",
        "-v", f"{docking_dir}:/app/docking_runs",
        "-v", f"{md_dir}:/app/md_runs",
        IMAGE_NAME
    ]
    res = subprocess.run(run_cmd, capture_output=True, text=True, encoding="utf-8", errors="replace")
    if res.returncode != 0:
        print(f"\n[!] Failed to start container:\n{res.stderr}")
        input("\nPress Enter to close...")
        return

    print("[*] Waiting for BIOBIO interface to become ready...")
    ready = False
    for _ in range(45):
        try:
            with urllib.request.urlopen(URL, timeout=1) as r:
                if r.status == 200:
                    ready = True
                    break
        except Exception:
            time.sleep(1)

    if not ready:
        logs = subprocess.run(["docker", "logs", CONTAINER_NAME], capture_output=True, text=True, encoding="utf-8", errors="replace")
        print(f"\n[!] Container logs:\n{logs.stdout}\n{logs.stderr}")
        input("\nPress Enter to close...")
        return

    print("[+] BIOBIO is running! Opening Desktop Window...")
    edge_paths = [
        r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
        r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
    ]
    edge_exe = next((p for p in edge_paths if os.path.exists(p)), None)
    if edge_exe:
        subprocess.Popen([edge_exe, f"--app={URL}", "--window-size=1400,900"])
    else:
        import webbrowser
        webbrowser.open(URL)

    print("\n==================================================")
    print(" BIOBIO is now running in its own Desktop Window.")
    print(" Keep this window open while working.")
    print(" Press Enter here at any time to stop and exit.")
    print("==================================================")
    input()
    print("\n[*] Stopping BIOBIO Engine safely...")
    subprocess.run(["docker", "rm", "-f", CONTAINER_NAME], capture_output=True)

if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        input("\nPress Enter to exit...")
