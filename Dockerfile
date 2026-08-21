FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential cmake git patch python3 python3-dev python3-pip python3-venv \
    openjdk-17-jdk maven zstd zlib1g-dev llvm-17 clang-17 libc++-17-dev \
    libc++abi-17-dev libtbb-dev ca-certificates curl && \
    rm -rf /var/lib/apt/lists/*
WORKDIR /workspace
COPY requirements.txt requirements-lsf.txt ./
RUN python3 -m venv .venv && .venv/bin/pip install --upgrade pip && \
    .venv/bin/pip install -r requirements.txt -r requirements-lsf.txt
COPY . .
RUN ./setup.sh --no-system --skip-python
ENV AUTOCSF_DOCKER_DESCRIPTION="Docker Desktop linux/amd64"
ENTRYPOINT ["/bin/bash"]
