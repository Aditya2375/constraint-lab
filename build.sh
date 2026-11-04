#!/bin/sh
set -eu
em++ engine.cpp -std=c++20 -O3 -Wall -Wextra -Werror --no-entry -s WASM=1 -s MODULARIZE=1 -s EXPORT_ES6=1 -s EXPORT_NAME=createEngine -s ENVIRONMENT=web,node -s ALLOW_MEMORY_GROWTH=0 -s INITIAL_MEMORY=16777216 -s EXPORTED_FUNCTIONS='["_reset","_advance","_body_count","_joint_count","_body_value","_joint_value","_metric","_configure","_spawn","_grab","_drag","_release"]' -o docs/engine.js
chmod 644 docs/engine.js docs/engine.wasm
