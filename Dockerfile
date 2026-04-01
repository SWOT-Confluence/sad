# Stage 0 - Base Python image
FROM python:3.12-slim AS stage0
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libxml2 \
    libhdf5-dev \
    libnetcdf-dev \
    && rm -rf /var/lib/apt/lists/*

# Stage 1 - Python dependencies
FROM stage0 AS stage1
COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt

# Stage 2 - Copy algorithm source
FROM stage1 AS stage2
COPY src/preprocess.py  /app/preprocess.py
COPY src/priors.py      /app/priors.py
COPY src/gvf.py         /app/gvf.py
COPY src/rejection.py   /app/rejection.py
COPY src/infer.py       /app/infer.py
COPY src/utils.py       /app/utils.py
COPY swot.py        /app/swot.py
COPY ./sos_read     /app/sos_read/

# Stage 3 - Final image
FROM stage2 AS stage3
LABEL version="1.0" \
    description="Containerized SAD algorithm." \
    "confluence.contact"="ntebaldi@umass.edu" \
    "algorithm.contact"="kandread@umass.edu"

ENV PYTHONPATH="/app:${PYTHONPATH}"
# Enable JAX float64 globally via environment variable
ENV JAX_ENABLE_X64="1"

WORKDIR /app
ENTRYPOINT ["python", "/app/swot.py"]
