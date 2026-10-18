#!/bin/bash
# Loop casadi's Julia suite against CasADiNative; report crashes.
set -u
N=$1; MODE=$2
mkdir -p stub && echo 'import CasADiNative' > stub/CasADiNative.jl
pre=""; [ "$MODE" = nofinal ] && pre="GC.enable_finalizers(false);"
crashes=0
for i in $(seq $N); do
  CASADI_JL=$PWD/stub julia --startup-file=no --project=env -e "$pre using CasADiNative, Test; m=Module(:T); Core.eval(m, :(include(p)=Base.include(\$m,p))); Base.include(m, joinpath(\"cn\",\"test\",\"casadi\",\"test\",\"julia\",\"runtests.jl\"))" > run$i.log 2>&1
  rc=$?
  echo "run $i rc=$rc"
  if [ $rc -ne 0 ]; then
    crashes=$((crashes+1))
    (grep -n -A60 'signal' run$i.log | grep -v nope_nope | head -80 || true)
    (grep -E 'Error During|did not pass' run$i.log || true)
  fi
done
echo "MODE=$MODE crashes=$crashes/$N"
