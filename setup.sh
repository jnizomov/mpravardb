#!/bin/bash

# add CRAN to apt sources
apt-key adv --keyserver keyserver.ubuntu.com \
    --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
printf '\ndeb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/\n' \
    | tee -a /etc/apt/sources.list
add-apt-repository -y ppa:c2d4u.team/c2d4u4.0+

# install system requirements & R
export DEBIAN_FRONTEND=noninteractive
apt-get -y update
apt-get -y upgrade
apt-get -yq install \
    software-properties-common \
    libopenblas-dev \
    libsodium-dev \
    r-base \
    r-base-dev \
    r-cran-shiny \
    r-cran-rmarkdown \
    r-cran-plotly \
    r-cran-ggplot2 \
    r-cran-jsonlite \
    r-cran-forecast