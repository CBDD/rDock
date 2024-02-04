FROM python:3.12 AS base

RUN DEBIAN_FRONTEND=noninteractive apt update -y \
    && apt install -y --no-install-recommends \
        swig \
        openbabel \
        libopenbabel-dev \
    && rm -rf /var/lib/apt/lists/*

# this is a bit of a mess but it is a quick way to be able to
# pip install openbable, we will improve it later
RUN ln -s /usr/include/openbabel3 /usr/local/include/openbabel3
RUN ln -s /lib/x86_64-linux-gnu/libopenbabel.so /usr/local/lib/libopenbabel.so


COPY requirements.txt requirements.txt
RUN python -m pip install --upgrade pip
RUN python -m pip install -r requirements.txt

FROM base AS development

COPY requirements-dev.txt requirements-dev.txt
RUN python -m pip install -r requirements-dev.txt

FROM base AS production

COPY . /rdock-utils
WORKDIR /rdock-utils
RUN python -m pip install .
