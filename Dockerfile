# SAD algorithm image.
# Derives from the pre-built base image that has all dependencies installed.
#

FROM kandread/sad-base:1.0

COPY src/preprocess.py  /app/preprocess.py
COPY src/priors.py      /app/priors.py
COPY src/gvf.py         /app/gvf.py
COPY src/rejection.py   /app/rejection.py
COPY src/infer.py       /app/infer.py
COPY src/utils.py       /app/utils.py
COPY src/swot.py        /app/swot.py
COPY ./sos_read     /app/sos_read/

LABEL version="1.0" \
    description="Containerized SAD algorithm." \
    "confluence.contact"="ntebaldi@umass.edu" \
    "algorithm.contact"="kandread@umass.edu"

ENV PYTHONPATH="/app:${PYTHONPATH}"
ENV JAX_ENABLE_X64="1"

WORKDIR /app
ENTRYPOINT ["python3", "/app/swot.py"]
