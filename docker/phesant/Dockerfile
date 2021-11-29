FROM rocker/r-ver:3.3.1

COPY install_packages.R install_packages.R

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    git \ 
    build-essential \
    libcurl4-gnutls-dev \
    libxml2-dev \
    libssl-dev \
    zlib1g-dev && \
    git clone https://github.com/MRCIEU/PHESANT.git && \
    Rscript install_packages.R

