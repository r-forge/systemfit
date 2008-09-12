#!/bin/sh
echo "Sweave(\"systemfit.Rnw\")" | LC_ALL="C" R --no-save --no-restore

echo "Stangle(\"systemfit.Rnw\", annotate=FALSE, split=TRUE,
   prefix=FALSE)" | LC_ALL="C" R --no-save --no-restore

echo "Sweave(\"systemfit_reliability.Snw\")" | LC_ALL="C" R --no-save --no-restore
echo "Stangle(\"systemfit_reliability.Snw\", annotate=FALSE, split=TRUE,
   prefix=FALSE)" | LC_ALL="C" R --no-save --no-restore

echo "Sweave(\"systemfit_timings.Snw\")" | R --no-save --no-restore
echo "Stangle(\"systemfit_timings.Snw\", annotate=FALSE, split=TRUE,
   prefix=FALSE)" | LC_ALL="C" R --no-save --no-restore

echo "Sweave(\"systemfit_sem.Snw\")" | R --no-save --no-restore
echo "Stangle(\"systemfit_sem.Snw\", annotate=FALSE, split=TRUE,
   prefix=FALSE)" | LC_ALL="C" R --no-save --no-restore
