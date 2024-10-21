# Use the official R image as a base
FROM r-base:4.3.0

# Set the working directory inside the container
WORKDIR /usr/src/app

# Install necessary system dependencies
RUN apt-get update && apt-get install -y \
    procps \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install required R packages
RUN R -e "install.packages(c('ggplot2', 'dplyr', 'tidyr', 'cowplot', 'scales'), repos='http://cran.rstudio.com/')"

# Copy the R script to the working directory
COPY bedgraph-visualizer.R .



# Example usage:
# docker build -t my-r-script .
# docker run --rm -v $(pwd)/output:/usr/src/app/output my-r-script BAF.bedgraph.gz LRR.bedgraph.gz regions.bed output
