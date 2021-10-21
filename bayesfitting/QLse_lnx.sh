#!/bin/sh

f2py3 -c --fcompiler=gnu95 --compiler=unix \
 -m QLse_lnx64 \
  QLse_main.f90 only: qlstexp : QLse_subs.f90 BlrRes.f90 Bayes.f90 Four.f90 Util.f90 Simopt.f90
