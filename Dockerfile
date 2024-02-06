FROM ubuntu:22.04

RUN apt-get update && apt-get install -y --no-install-recommends make g++ libpopt-dev libpopt0
RUN apt-get install -y --no-install-recommends git clang-tidy clang-format ssh
