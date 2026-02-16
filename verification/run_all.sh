#!/bin/bash
set -e
echo "FIRST PROOF BENCHMARK — COMPLETE VERIFICATION"
echo "==============================================="
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"
for dir in p1 p2 p3_p4 p5 p6 p7 p9 p10; do
    echo ""
    echo "--- $dir ---"
    cd "$SCRIPT_DIR/$dir"
    for f in *.py; do
        [ -f "$f" ] || continue
        echo "  Running $f..."
        python3 "$f" 2>&1 | tail -5
    done
done
echo ""
echo "ALL VERIFICATIONS COMPLETE"
