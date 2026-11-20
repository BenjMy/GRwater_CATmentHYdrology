#!/usr/bin/env bash
# Deprecated: the batch runner is now scenario_runner.py (same flags, same run order, same
# skip-if-done check). This wrapper only keeps old command lines working.
exec python3 "$(dirname "${BASH_SOURCE[0]}")/scenario_runner.py" "$@"
