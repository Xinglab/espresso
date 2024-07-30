FROM debian:bullseye

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       ca-certificates \
       curl \
       git \
       make \
    && rm -rf /var/lib/apt/lists/* \
    # Install conda to /conda
    && mkdir /conda \
    && cd /conda \
    && curl -L 'https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh' -O \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /conda/install \
    && /conda/install/bin/conda init \
    # The mamba installer solves the package dependencies for ESPRESSO
    # using much less RAM
    && /conda/install/bin/conda install conda-libmamba-solver \
    && /conda/install/bin/conda config --set solver libmamba \
    && git clone https://github.com/Xinglab/espresso.git /espresso \
    && cd /espresso \
    # && git checkout {commit} \
    && cd /espresso/snakemake \
    && ./install

# Make conda installed programs available on PATH
ENV PATH /espresso/snakemake/conda_env/bin:${PATH}

# Set defaults for running the image.
# The ENTRYPOINT AND CMD are empty to be compatible with
# CWL and WDL implementations that cannot override those values
WORKDIR /espresso
ENTRYPOINT []
CMD []
