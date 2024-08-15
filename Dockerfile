# Stage 0 - Create from julia image and install OS packages
FROM julia:latest as stage0
RUN apt update && apt -y install bzip2 build-essential

# STAGE 1 - Python and python packages for S3 functionality
FROM stage0 as stage1
RUN apt update && apt -y install python3 python3-dev python3-pip python3-venv python3-boto3

# Stage 2 - Install SAD dependencies
FROM stage1 as stage2
RUN git submodule update --init --recursive
RUN mkdir -p /usr/local/bin/julia_pkgs
ENV JULIA_LOAD_PATH="/usr/local/bin/julia_pkgs:$JULIA_LOAD_PATH"
ENV JULIA_DEPOT_PATH="/usr/local/bin/julia_pkgs:$JULIA_DEPOT_PATH"
ENV PYTHON="/usr/bin/python3"
COPY deps.jl /app/deps.jl
RUN julia /app/deps.jl \
	&& find /usr/local/bin/julia_pkgs -type d -exec chmod 755 {} \; \
	&& find /usr/local/bin/julia_pkgs -type f -exec chmod 644 {} \;

# Stage 3 - Copy SWOT script
FROM stage2 as stage3
COPY swot.jl /app/swot.jl
COPY ./sos_read /app/sos_read/

# Stage 4 - Execute algorithm
FROM stage3 as stage4
LABEL version="1.0" \
	description="Containerized SAD algorithm." \
	"confluence.contact"="ntebaldi@umass.edu" \
	"algorithm.contact"="kandread@umass.edu"
ENTRYPOINT ["/usr/local/julia/bin/julia", "/app/swot.jl"]