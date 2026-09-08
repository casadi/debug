"""Windows Opti lifetime probe with real Ctrl+C during native MUMPS work."""
import argparse
import collections
import ctypes
import gc
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent


def system_memory():
    from ctypes import wintypes as w
    class PerformanceInfo(ctypes.Structure):
        _fields_ = [('cb', w.DWORD)] + [(name, ctypes.c_size_t) for name in (
            'CommitTotal', 'CommitLimit', 'CommitPeak', 'PhysicalTotal',
            'PhysicalAvailable', 'SystemCache', 'KernelTotal', 'KernelPaged',
            'KernelNonpaged', 'PageSize')] + [(name, w.DWORD) for name in (
            'HandleCount', 'ProcessCount', 'ThreadCount')]
    info = PerformanceInfo()
    info.cb = ctypes.sizeof(info)
    api = ctypes.WinDLL('psapi', use_last_error=True).GetPerformanceInfo
    api.argtypes = [ctypes.POINTER(PerformanceInfo), w.DWORD]
    api.restype = w.BOOL
    if not api(ctypes.byref(info), info.cb):
        raise ctypes.WinError(ctypes.get_last_error())
    return {name: getattr(info, name) * info.PageSize for name in (
        'CommitTotal', 'CommitLimit', 'CommitPeak', 'PhysicalAvailable')}


def memory():
    import psutil
    process = psutil.Process()
    info = process.memory_info()
    return dict(rss=info.rss, private=info.private, pagefile=info.pagefile,
                handles=process.num_handles(), threads=process.num_threads())


def send_ctrl_c(pid):
    """Join only the worker's console; exclude the CI runner from the broadcast."""
    from ctypes import wintypes as w
    kernel = ctypes.WinDLL('kernel32', use_last_error=True)
    kernel.AttachConsole.argtypes = [w.DWORD]
    kernel.AttachConsole.restype = w.BOOL
    kernel.SetConsoleCtrlHandler.argtypes = [ctypes.c_void_p, w.BOOL]
    kernel.SetConsoleCtrlHandler.restype = w.BOOL
    kernel.GenerateConsoleCtrlEvent.argtypes = [w.DWORD, w.DWORD]
    kernel.GenerateConsoleCtrlEvent.restype = w.BOOL
    kernel.FreeConsole()
    if not kernel.AttachConsole(pid):
        raise ctypes.WinError(ctypes.get_last_error())
    if not kernel.SetConsoleCtrlHandler(None, True):
        raise ctypes.WinError(ctypes.get_last_error())
    if not kernel.GenerateConsoleCtrlEvent(0, 0):
        raise ctypes.WinError(ctypes.get_last_error())
    time.sleep(0.01)
    kernel.FreeConsole()


def worker(args):
    import casadi as ca
    log = open(args.output.with_suffix('.worker.log'), 'w', buffering=1)
    os.dup2(log.fileno(), 1)
    os.dup2(log.fileno(), 2)
    marker = ctypes.CDLL(str(HERE / 'marker.dll')) if args.mode == 'native' else None
    if marker:
        marker.arm_interrupt.argtypes = []
        marker.arm_interrupt.restype = None
    output = args.output.open('w', buffering=1)
    def record(**row):
        output.write(json.dumps(row) + '\n')
    record(version=ca.__version__, revision=ca.CasadiMeta.git_revision(),
           module=ca.__file__, python=sys.version, platform=platform.platform(),
           pid=os.getpid(), mode=args.mode, factor=args.factor, count=args.count)

    def cycle():
        o = ca.Opti()
        x = o.variable(100)
        o.minimize(ca.sumsqr(1-x[:-1]) + 100*ca.sumsqr(x[1:]-x[:-1]**2))
        o.subject_to(o.bounded(-3, x, 3))
        o.set_initial(x, 0)
        o.solver('ipopt', {'print_time': False}, {
            'print_level': 0, 'sb': 'yes', 'linear_solver': 'mumps',
            'max_iter': 1 if args.mode == 'failure' else 300})
        if marker: marker.arm_interrupt()
        caught = None
        try:
            solution = o.solve()
            if args.mode == 'success':
                assert solution.value(o.f) < 1e-8
        except (RuntimeError, KeyboardInterrupt) as ex:
            caught = type(ex).__name__
            if args.mode == 'success': raise
        status = o.stats()['return_status']
        allowed = {'success': {'Solve_Succeeded'},
                   'failure': {'Maximum_Iterations_Exceeded'},
                   'native': {'NonIpopt_Exception_Thrown', 'User_Requested_Stop'}}
        assert status in allowed[args.mode], (status, caught)
        if args.mode != 'success': assert caught is not None
        return status

    counts = collections.Counter()
    start = time.monotonic()
    record(iteration=0, elapsed=0, **memory())
    for i in range(1, args.count+1):
        counts[cycle()] += 1
        if i == 1 or i % 25 == 0 or i == args.count:
            gc.collect()
            sample = memory()
            record(iteration=i, elapsed=time.monotonic()-start, **sample)
            if sample['private'] > 2*1024**3:
                raise RuntimeError('Private committed memory exceeded 2 GiB safety limit')
    record(status_counts=dict(counts), completed=args.count)


def run(args):
    args.output.mkdir(parents=True, exist_ok=True)
    manifest = []
    for mode, factor in [('success', 0), ('failure', 0), ('native', 1), ('native', 3)]:
        name = mode + (f'-factor{factor}' if factor else '')
        output = (args.output / f'{name}.jsonl').resolve()
        command = [sys.executable, str(HERE / 'probe.py'), 'worker',
                   '--mode', mode, '--factor', str(factor), '--count', str(args.count),
                   '--output', str(output)]
        env = os.environ.copy()
        env.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        if mode == 'native':
            env.update(PROBE_PYTHON=sys.executable, PROBE_SCRIPT=str(HERE / 'probe.py'),
                       PROBE_FACTOR=str(factor), PROBE_COUNT=str(args.count),
                       PROBE_INJECTIONS=str(output.with_suffix('.injections.json')))
            command = [args.gdb, '-q', '-batch', '-x', str(HERE / 'native.gdb'),
                       '--args'] + command
        before = system_memory()
        print(f'Starting {name}: {args.count} cycles', flush=True)
        with output.with_suffix('.driver.log').open('w') as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                                    env=env, timeout=args.timeout)
        after = system_memory()
        rows = [json.loads(line) for line in output.read_text().splitlines()] if output.exists() else []
        complete = rows and rows[-1].get('completed') == args.count
        injection = None
        if mode == 'native' and output.with_suffix('.injections.json').exists():
            injection = json.loads(output.with_suffix('.injections.json').read_text())
            complete = complete and injection['injections'] == args.count and not injection['errors']
        elif mode == 'native':
            complete = False
        entry = dict(case=name, exit_code=result.returncode, complete=bool(complete),
                     system_before=before, system_after=after, injections=injection)
        manifest.append(entry)
        (args.output / 'manifest.json').write_text(json.dumps(manifest, indent=2))
        print(json.dumps(entry), flush=True)
    summarize(args.output)
    assert all(r['exit_code'] == 0 and r['complete'] for r in manifest), manifest


def summarize(directory):
    lines = ['| Case | Cycles | RSS at 250 / final (MiB) | Private at 250 / final (MiB) | Late private slope (KiB/cycle) |',
             '|---|---:|---:|---:|---:|']
    for path in sorted(directory.glob('*.jsonl')):
        rows = [json.loads(line) for line in path.read_text().splitlines()]
        samples = [r for r in rows if 'iteration' in r]
        if not samples: continue
        final = samples[-1]
        baseline = next((r for r in samples if r['iteration'] == 250), samples[0])
        late = [r for r in samples if r['iteration'] >= max(250, final['iteration']//2)]
        slope = 0.0
        if len(late) >= 2:
            mx = sum(r['iteration'] for r in late)/len(late)
            my = sum(r['private'] for r in late)/len(late)
            slope = sum((r['iteration']-mx)*(r['private']-my) for r in late)/sum((r['iteration']-mx)**2 for r in late)/1024
        lines.append(f"| {path.stem} | {final['iteration']} | {baseline['rss']/2**20:.3f} / {final['rss']/2**20:.3f} | {baseline['private']/2**20:.3f} / {final['private']/2**20:.3f} | {slope:.4f} |")
    report = '\n'.join(lines)+'\n\nMemory growth is diagnostic, not an automatic leak verdict. System commit includes other runner processes.\n'
    (directory / 'summary.md').write_text(report)
    print(report, flush=True)
    if os.environ.get('GITHUB_STEP_SUMMARY'):
        with open(os.environ['GITHUB_STEP_SUMMARY'], 'a') as f: f.write(report)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['run', 'worker', 'send'])
    parser.add_argument('--mode', choices=['success', 'failure', 'native'], default='success')
    parser.add_argument('--factor', type=int, default=1)
    parser.add_argument('--count', type=int, default=1000)
    parser.add_argument('--output', type=Path, default=Path('results'))
    parser.add_argument('--gdb', default='gdb')
    parser.add_argument('--timeout', type=int, default=600)
    parser.add_argument('--pid', type=int)
    args = parser.parse_args()
    if args.action == 'send': send_ctrl_c(args.pid)
    elif args.action == 'worker': worker(args)
    else: run(args)


if __name__ == '__main__': main()
