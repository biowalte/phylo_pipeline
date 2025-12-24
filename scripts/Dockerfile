FROM rocker/tidyverse:devel

RUN apt-get update && apt-get install -y \
    libglpk-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('ape', 'phytools', 'ggtree'), repos='https://cloud.r-project.org')"
RUN R -e "install.packages('BiocManager'); BiocManager::install('ggtree', update=FALSE, ask=FALSE)"
WORKDIR /data
COPY plot_phylogeny.R /usr/local/bin/plot_phylogeny.R
RUN chmod +x /usr/local/bin/plot_phylogeny.R
ENTRYPOINT ["plot_phylogeny.R"]
