FROM r-base:4.4.1

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    git \
    libcurl4-openssl-dev \
    libssl-dev \
  && rm -rf /var/lib/apt/lists/*

RUN R -e 'install.packages(c("tibble", "dplyr", "BiocManager"))' \
  && R -e 'BiocManager::install(c("IRanges", "Rsamtools", "S4Vectors", "GenomicRanges"))'

RUN cd /tmp \
  && git clone https://gitlab.com/jdotlim/assertify.git \
  && R CMD build assertify \
  && R CMD INSTALL assertify*.tar.gz \
  && rm -rf assertify*

COPY src/gelpers/ gelpers

RUN R CMD build gelpers \
  && R CMD INSTALL gelpers*.tar.gz \
  && rm -rf gelpers*

RUN cd opt \
  && mkdir -p gcnv/scripts

COPY src/scripts/*.R gcnv/scripts/
