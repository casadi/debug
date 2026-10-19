# Runs each open_issue_repros/<n>/repro.py in a subprocess and prints a summary.
import os, sys, subprocess, time
here = os.path.dirname(os.path.abspath(__file__))
only = sys.argv[1].split(",") if len(sys.argv) > 1 else None
rows = []
for n in sorted(d for d in os.listdir(here) if d.isdigit()):
    if only and n not in only: continue
    print(f"\n######## #{n}", flush=True)
    t = time.time()
    try:
        r = subprocess.run([sys.executable, os.path.join(here, n, "repro.py")], capture_output=True, text=True, timeout=1500, cwd=os.path.join(here, n))
        out = r.stdout + r.stderr
        verdict = "BUG PRESENT" if "BUG PRESENT" in r.stdout else "BUG ABSENT" if "BUG ABSENT" in r.stdout else f"NO VERDICT (rc={r.returncode})"
    except subprocess.TimeoutExpired as e:
        out = (e.stdout or b"").decode(errors="replace") if isinstance(e.stdout, bytes) else (e.stdout or "")
        verdict = "TIMEOUT"
    print(out[-6000:], flush=True)
    rows.append((n, verdict, time.time() - t))
print("\n======== SUMMARY")
for n, v, dt in rows: print(f"RESULT #{n}: {v} ({dt:.0f}s)")
