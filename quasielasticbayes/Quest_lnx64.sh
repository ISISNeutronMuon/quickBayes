#!/bin/sh

f2py3 -c --fcompiler=gnu95 --compiler=unix \
 -m Quest_lnx64 \
  Quest_main.f90 only: quest : Quest_subs.f90 BlrRes.f90 Bayes.f90 Four.f90 Util.f90 Simopt.f90
