FROM rstudio/plumber

#system libraries of general use
RUN apt-get update && apt-get install -y \
    libssl-dev \
    build-essential \
    gfortran \
    libjpeg-dev \
    xorg \
    libx11-dev \
    libglu1-mesa-dev \
    libfreetype6-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libmariadbclient-dev \
    libbz2-dev \
    liblzma-dev

# Install R packages
RUN R -e "install.packages('Morpho')"
RUN R -e "install.packages('plumber')"
RUN R -e "install.packages('dplyr')"
RUN R -e "install.packages('dbplyr')"
RUN R -e "install.packages('devtools')"
RUN R -e "install.packages('jsonlite')"
RUN R -e "devtools::install_github('hlorenzo/ddsPLS')"
RUN R -e "devtools::install_github('J0vid/Jovid')"
RUN R -e "install.packages('BiocManager')"
RUN R -e 'BiocManager::install(c("GenomicFeatures", "AnnotationDbi", "org.Mm.eg.db"))'
RUN R -e "install.packages('future')"

# select port
EXPOSE 3838

WORKDIR /srv/shiny-server/MGP_API/

# Run R file
CMD ["./plumber.R"]
