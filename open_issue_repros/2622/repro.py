# Issue #2622: Jupyter kernel dies on `import casadi` (macOS, Anaconda). Runs import in a fresh subprocess,
# and additionally through an ipykernel (jupyter_client) if available.
import subprocess, sys
code = "import faulthandler; faulthandler.enable(); import casadi; x=casadi.MX.sym('x'); print(casadi.Function('f',[x],[casadi.sin(x)])(1.0)); print('IMPORT_OK')"
p = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
print(p.stdout, p.stderr[-2000:])
ok = p.returncode == 0 and "IMPORT_OK" in p.stdout
try:
    from jupyter_client.manager import start_new_kernel
    km, kc = start_new_kernel(kernel_name="python3")
    msg = kc.execute_interactive("import casadi; print('KERNEL_OK', casadi.__version__)", timeout=120)
    print("kernel status:", msg["content"]["status"]); ok = ok and msg["content"]["status"] == "ok" and km.is_alive()
    kc.stop_channels(); km.shutdown_kernel(now=True)
except ImportError:
    print("jupyter_client not installed; tested plain subprocess import only")
except Exception as e:
    print("kernel error:", e); ok = False
print("BUG ABSENT" if ok else "BUG PRESENT")
