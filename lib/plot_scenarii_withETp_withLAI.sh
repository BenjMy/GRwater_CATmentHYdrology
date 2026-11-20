#!/usr/bin/env bash
# Kept for backward compatibility: the post-processor is now post_process.py (same options).
exec python "$(dirname "$0")/post_process.py" "$@"
