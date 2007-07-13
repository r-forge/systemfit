#!/bin/sh
echo "Sweave(\"systemfit_code.Snw\")" | R --no-save --no-restore

echo "Sweave(\"systemfit_usage.Snw\")" | R --no-save --no-restore

echo "Sweave(\"systemfit_reliability.Snw\")" | R --no-save --no-restore

echo "Sweave(\"systemfit_sem.Snw\")" | R --no-save --no-restore
