# Issue 3605: Windows Opti memory experiment

Runs published CasADi 3.6.5 and 3.7.2 wheels on Windows Server 2022,
Python 3.11, using Ipopt/MUMPS with one BLAS/OpenMP thread. Each case creates,
solves, and releases 1,000 independent 100-variable bounded Rosenbrock Opti
problems in one Python process. The four cases are successful solves,
iteration-limit failures, and Ctrl+C during the first or third numeric
factorization. A three-cycle pilot validates each path before the long run.

For native interruptions, GDB arms a breakpoint in MUMPS's LDLT routine
after the requested factorization driver entry. At that breakpoint a helper
joins the worker's separate console and calls GenerateConsoleCtrlEvent with
CTRL_C_EVENT. GDB passes the resulting console exception to the application.
The helper ignores its own Ctrl+C. The CI shell does not share this console.
The first native stack is logged; every injection is counted. Every solve's
return status and each worker's completion count are checked.

JSONL samples record current working set, private committed bytes, the process
pagefile counter, handles, and threads after scope exit and garbage collection.
System commit is recorded before and after each process exits. Workflow logs
also record allocated pagefile size and pagefile usage before and after the
experiment. System counters include other runner activity. They are not a
substitute for per-process memory trends. No heap trimming is requested.

Artifacts contain raw measurements, native stacks, worker diagnostics,
injection counts, environment versions, and a Markdown summary with the late
private-memory slope. A plateau is evidence against sustained growth for this
workload, not proof that all workloads and interruption locations are leak-free.
Windows Server 2022 is not necessarily the reporter's Windows version.

This isolated branch adapts the already registered workflow
`.github/workflows/linux-valgrind-corruptor.yml`. Dispatch with:

```sh
gh workflow run linux-valgrind-corruptor.yml --repo casadi/debug --ref issue-3605-windows-memory
```

Push experiment commits with `[skip ci]` to avoid starting unrelated push
workflows. The debug repository's main branch retains its original workflow.
