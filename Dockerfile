# Stage 0 - Create from julia image and install OS packages
FROM julia:latest as stage0
RUN apt update && apt -y install bzip2 build-essential

# Stage 1 - Install SAD dependencies
FROM stage0 as stage1
RUN mkdir -p /usr/local/bin/julia_pkgs
ENV JULIA_LOAD_PATH="/usr/local/bin/julia_pkgs:$JULIA_LOAD_PATH"
ENV JULIA_DEPOT_PATH="/usr/local/bin/julia_pkgs:$JULIA_DEPOT_PATH"
COPY deps.jl /app/deps.jl
RUN julia /app/deps.jl \
	&& find /usr/local/bin/julia_pkgs -type d -exec chmod 755 {} \; \
	&& find /usr/local/bin/julia_pkgs -type f -exec chmod 644 {} \;

# Stage 2 - Copy SWOT script
FROM stage1 as stage2
COPY swot.jl /app/swot.jl

# Stage 3 - Execute algorithm
FROM stage2 as stage3
LABEL version="1.0" \
	description="Containerized SAD algorithm." \
	"confluence.contact"="ntebaldi@umass.edu" \
	"algorithm.contact"="kandread@umass.edu"
ENTRYPOINT ["/usr/local/julia/bin/julia", "/app/swot.jl"]