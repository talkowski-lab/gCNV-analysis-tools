FROM r-base:4.4.1

ARG HTSLIB_VERSION="1.21"
ARG HTSLIB_URI="https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2"

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    git \
    izlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libdeflate-dev \
    libcurl4-openssl-dev \
    libssl-dev \
  && rm -rf /var/lib/apt/lists/*

RUN curl -L -o htslib.tar.bz2 "${HTSLIB_URI}" \
  && tar -jxf htslib.tar.bz2 \
  && cd "htslib-${HTSLIB_VERSION}" \
  && ./configure --prefix=/usr/local \
       --enable-libcurl \
       --enable-gcs \
       --enable-s3 \
       --with-libdeflate \
  && make \
  && make install \
  && cd .. \
  && rm -rf "htslib-${HTSLIB_VERSION}" htslib.tar.bz2

RUN R -e 'install.packages(c("data.table", "R.utils", "BiocManager"))' \
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

RUN mkdir -p /opt/gcnv/scripts

COPY src/scripts/*.R /opt/gcnv/scripts/
