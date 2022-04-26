FROM debian:bullseye-20220418-slim AS build
RUN \
  apt-get update -y \
  && apt-get install -y --no-install-recommends \
  gcc=4:10.2.1-1 \
  g++=4:10.2.1-1 \
  make=4.3-4.1 \
  openbabel=3.1.1+dfsg-6 \
  libcppunit-dev=1.15.1-2 \
  libpopt-dev=1.18-2 \
  && rm -rf /var/lib/apt/lists/*
ENV RBT_FILE rdock
ENV RBT_ROOT /$RBT_FILE
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:$RBT_ROOT/lib
ENV PATH $PATH:$RBT_ROOT/bin
WORKDIR /rdock
COPY build /rdock/build
COPY import /rdock/import
COPY include /rdock/include
COPY lib /rdock/lib
COPY src /rdock/src
COPY bin /rdock/bin
COPY data /rdock/data
WORKDIR /rdock/build
RUN make -j 4
# hadolint ignore=DL3059
RUN useradd -m -s /bin/bash -G 0 rdock
USER rdock
WORKDIR /home/rdock