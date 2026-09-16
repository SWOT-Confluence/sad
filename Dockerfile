# Stage 0 - Create from Python3.12 image
FROM python:3.12-slim as stage0

# Stage 1 - Debian dependencies
FROM stage0 as stage1
RUN apt update \
	&& DEBIAN_FRONTEND=noninteractive apt install -y curl git zip python3-dev build-essential libhdf5-serial-dev netcdf-bin libnetcdf-dev

# Stage 2 - Create virtual environment and install dependencies
FROM stage1 as stage2
COPY requirements.txt /app/requirements.txt
RUN /usr/local/bin/python3 -m venv /app/env
RUN /app/env/bin/pip install -r /app/requirements.txt

# Stage 3 - Copy SoS reader
FROM stage2 as stage3
COPY ./sos_read /app/sos_read/

# Stage 4 - Execute algorithm
FROM stage3 as stage4
COPY run_sad.py /app/run_sad.py
LABEL version="2.0" \
	description="Containerized SAD algorithm (SADnm)." \
	"confluence.contact"="ntebaldi@umass.edu" \
	"algorithm.contact"="kandread@umass.edu"
ENTRYPOINT ["/app/env/bin/python3", "/app/run_sad.py"]
