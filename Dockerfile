FROM ubuntu:22.04

# Install system dependencies including all scientific libraries
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-dev \
    bowtie2 \
    freebayes \
    samtools \
    build-essential \
    libhdf5-dev \
    llvm \
    pkg-config \
    wget \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Install meteor2 - ONE command that just works!
RUN pip3 install --no-cache-dir meteor

CMD ["/bin/bash"]
