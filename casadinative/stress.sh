#!/bin/bash
# Loop casadi's Julia suite against CasADiNative; report crashes.
# MODE: full | noconic | onlyconic | scribble | gmalloc
set -u
N=$1; MODE=$2
mkdir -p stub && echo 'import CasADiNative' > stub/CasADiNative.jl
rm -rf t && cp -R cn/test/casadi/test/julia t
case $MODE in
  noconic) rm t/conic.jl ;;
  onlyconic) find t -name '*.jl' ! -name runtests.jl ! -name conic.jl ! -name '_*' -delete ;;
  scribble) export MallocScribble=1 MallocPreScribble=1 MallocGuardEdges=1 MallocCheckHeapStart=1000 MallocCheckHeapEach=1000 ;;
  gmalloc) export DYLD_INSERT_LIBRARIES=/usr/lib/libgmalloc.dylib ;;
esac
ls t
crashes=0
for i in $(seq $N); do
  CASADI_JL=$PWD/stub julia --startup-file=no --project=env -e 'using CasADiNative, Test; m=Module(:T); Core.eval(m, :(include(p)=Base.include($m,p))); Base.include(m, joinpath("t","runtests.jl"))' > run$i.log 2>&1
  rc=$?
  echo "run $i rc=$rc"
  if [ $rc -ne 0 ]; then
    crashes=$((crashes+1))
    (grep -v nope_nope run$i.log | grep -n -B3 -A40 -E 'signal|GuardMalloc|malloc:|Malloc' | head -120 || true)
    (grep -E 'Error During|did not pass' run$i.log || true)
  fi
done
echo "MODE=$MODE crashes=$crashes/$N"
