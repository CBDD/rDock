FROM debian:bullseye-slim AS build

RUN apt update -y && apt install -y gcc g++ make openbabel libpopt-dev libcppunit-dev && rm -rf /var/lib/apt/lists/*

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

RUN useradd -m -s /bin/bash -G 0 rdock
USER rdock
WORKDIR /home/rdock
