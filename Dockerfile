FROM mambaorg/micromamba:2.0.3-ubuntu22.04 AS app

# build and run as root users since micromamba image has 'mambauser' set as the $USER
USER root
# set workdir to default for building; set to /data at the end
WORKDIR /

ARG CLOCI_VER="0.3.1"

# Label the image with metadata that might be important to the user
LABEL base.image="mambaorg/micromamba:2.0.3-ubuntu22.04"
LABEL dockerfile.version="1"
LABEL software="cloci"
LABEL software.version="${CLOCI_VER}"
LABEL description="Co-occurrence Locus and Orthologous Cluster Identifier"
LABEL website="https://github.com/xonq/cloci/"
LABEL license="https://github.com/xonq/cloci/blob/master/LICENSE"
LABEL maintainer="Zachary Konkel"
LABEL maintainer.email="konkelzach@protonmail.com"

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install --no-install-recommends -y \
  wget \
  ca-certificates \
  procps && \
  apt-get autoclean && rm -rf /var/lib/apt/lists/*

COPY env.yaml /tmp/env.yaml

RUN micromamba install \
  --yes \
  --name base \
  -c conda-forge -c bioconda \
  --file /tmp/env.yaml \
  && micromamba clean --all -y \
  && rm /tmp/env.yaml

RUN wget https://github.com/xonq/cloci/archive/refs/tags/${CLOCI_VER}.tar.gz \
  && tar -xzvf ${CLOCI_VER}.tar.gz

ENV PYTHONPATH="/cloci-${CLOCI_VER}/cloci/:${PYTHONPATH}" 
RUN mv cloci-${CLOCI_VER}/cloci/main.py cloci-${CLOCI_VER}/cloci/cloci

ENV PATH="/opt/conda/bin/:/cloci-${CLOCI_VER}/cloci/:${PATH}" \
  LC_ALL=C.UTF-8

# final working directory in "app" layer is /data for passing data in/out of container
WORKDIR /cloci-${CLOCI_VER}/
