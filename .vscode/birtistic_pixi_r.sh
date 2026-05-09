#!/bin/bash
# R wrapper script for birtistic-pixi environment
cd "$(dirname "$0")/.." || exit 1

# Prefer radian for syntax highlighting and better REPL UX; fallback to base R.
if pixi run which radian >/dev/null 2>&1; then
	exec pixi run radian "$@"
else
	exec pixi run R "$@"
fi
