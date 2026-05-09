#!/bin/bash
# Python wrapper script for birtistic-pixi environment
cd "$(dirname "$0")/.." || exit 1
exec pixi run python "$@"
