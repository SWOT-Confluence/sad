# SAD algorithm image.
# Derives from the pre-built base image that has all Julia dependencies
# compiled into a sysimage at /usr/local/julia/bin/julia_base.so.
# Only rebuild this when algorithm source changes
# Should be fast since no compilation happens at build time.

FROM ghcr.io/kandread/sad-julia-base:1.0

COPY Sad.jl/Project.toml  /app/Project.toml
COPY Sad.jl/Manifest.toml /app/Manifest.toml
COPY Sad.jl/src/Sad.jl        /app/Sad.jl
COPY Sad.jl/src/preprocess.jl /app/preprocess.jl
COPY Sad.jl/src/priors.jl     /app/priors.jl
COPY Sad.jl/src/gvf.jl        /app/gvf.jl
COPY Sad.jl/src/rejection.jl  /app/rejection.jl
COPY Sad.jl/src/inference.jl  /app/inference.jl
COPY swot.jl       /app/swot.jl
COPY --chmod=755 run.sh /app/run.sh
COPY ./sos_read    /app/sos_read/

LABEL version="1.0" \
    description="Containerized SAD algorithm." \
    "confluence.contact"="ntebaldi@umass.edu" \
    "algorithm.contact"="kandread@umass.edu"

ENTRYPOINT ["/app/run.sh"]
