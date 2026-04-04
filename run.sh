#!/bin/bash

export JULIA_DEPOT_PATH="/tmp/.julia:/.julia"
export JULIA_PKG_DEVDIR="/tmp"
exec /usr/local/julia/bin/julia \
    --sysimage=/usr/local/julia/bin/julia_base.so \
    --project=/app \
    /app/swot.jl "$@"
