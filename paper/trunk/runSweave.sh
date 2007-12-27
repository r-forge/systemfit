#!/bin/sh
echo "Sweave(\"systemfit_code.Snw\")" | R --no-save --no-restore

echo "Sweave(\"systemfit_usage.Snw\")" | R --no-save --no-restore
echo "Stangle(\"systemfit_usage.Snw\", annotate=FALSE, split=TRUE,
   prefix=FALSE)" | R --no-save --no-restore

echo "Sweave(\"systemfit_reliability.Snw\")" | R --no-save --no-restore
echo "Stangle(\"systemfit_reliability.Snw\", annotate=FALSE, split=TRUE,
   prefix=FALSE)" | R --no-save --no-restore

echo "Sweave(\"systemfit_timings.Snw\")" | R --no-save --no-restore
echo "Stangle(\"systemfit_timings.Snw\", annotate=FALSE, split=TRUE,
   prefix=FALSE)" | R --no-save --no-restore

echo "Sweave(\"systemfit_sem.Snw\")" | R --no-save --no-restore
echo "Stangle(\"systemfit_sem.Snw\", annotate=FALSE, split=TRUE,
   prefix=FALSE)" | R --no-save --no-restore
