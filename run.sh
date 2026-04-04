#!/bin/bash

export JULIA_DEPOT_PATH="/.julia"
exec /usr/local/julia/bin/julia \
    --sysimage=/usr/local/julia/bin/julia_base.so \
    --project=/app \
    /app/swot.jl "$@"
