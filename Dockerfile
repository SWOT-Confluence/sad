# SAD algorithm image.
# Derives from the pre-built base image that has all dependencies installed.
#
# Stage 0 - Create from Python3.12 image
FROM python:3.12-slim as stage0

# Stage 1 - Debian dependencies
FROM stage0 as stage1
RUN apt update \
    && DEBIAN_FRONTEND=noninteractive apt install -y curl zip python3-dev build-essential libhdf5-serial-dev netcdf-bin libnetcdf-dev

# Stage 2 - Create virtual environment and install dependencies
FROM stage1 as stage2
COPY requirements.txt /app/requirements.txt
RUN /usr/local/bin/python3 -m venv /app/env
RUN /app/env/bin/pip install -r /app/requirements.txt

# Stage 3 - Copy SAD code
FROM stage2 as stage3
COPY src/preprocess.py  /app/preprocess.py
COPY src/priors.py      /app/priors.py
COPY src/gvf.py         /app/gvf.py
COPY src/rejection.py   /app/rejection.py
COPY src/infer.py       /app/infer.py
COPY src/utils.py       /app/utils.py
COPY ./sos_read /app/sos_read/

# Stage 4 - Execute algorithm
FROM stage3 as stage4
COPY swot.py        /app/swot.py
LABEL version="1.0" \
    description="Containerized SAD algorithm." \
    "confluence.contact"="ntebaldi@umass.edu" \
    "algorithm.contact"="kandread@umass.edu"
ENV JAX_ENABLE_X64="1"
ENTRYPOINT ["/app/env/bin/python3", "/app/swot.py"]
