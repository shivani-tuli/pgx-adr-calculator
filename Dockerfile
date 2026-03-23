# ============================================================
# PGx-ADR Risk Calculator — Docker Image
# Includes: R 4.4 + Shiny + Java 17 + PharmCAT
# ============================================================
FROM rocker/shiny:4.4.0

LABEL maintainer="shivani-tuli"
LABEL description="PGx-ADR Risk Calculator with PharmCAT star-allele calling"

# Install system dependencies + Java 17
RUN apt-get update && apt-get install -y --no-install-recommends \
    openjdk-17-jre-headless \
    python3 \
    python3-pip \
    python3-venv \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set JAVA_HOME
ENV JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64
ENV PATH="${JAVA_HOME}/bin:${PATH}"

# Install R packages
RUN R -e 'install.packages(c( \
    "shiny", "shinydashboard", "DT", "plotly", "ggplot2", \
    "dplyr", "tidyr", "jsonlite", "httr", "scales" \
), repos="https://cran.r-project.org")'

# Create app directory
RUN mkdir -p /srv/shiny-server/pgx-adr-calculator
WORKDIR /srv/shiny-server/pgx-adr-calculator

# Copy application files
COPY app.R .
COPY R/ R/
COPY data/ data/
COPY scripts/ scripts/
COPY pharmcat/ pharmcat/

# Download PharmCAT if not already present
RUN if [ ! -f pharmcat/pharmcat.jar ]; then \
    echo "Downloading PharmCAT..." && \
    mkdir -p pharmcat && \
    curl -L -o pharmcat/pharmcat.jar \
      "https://github.com/PharmGKB/PharmCAT/releases/latest/download/pharmcat.jar"; \
    fi

# Set up PharmCAT Python environment
RUN if [ -f pharmcat/requirements.txt ]; then \
    python3 -m venv pharmcat/.venv && \
    pharmcat/.venv/bin/pip install -r pharmcat/requirements.txt; \
    fi

# Verify Java works
RUN java -version

# Shiny server config: single app mode
RUN echo 'run_as shiny;\n\
server {\n\
  listen 3838;\n\
  location / {\n\
    site_dir /srv/shiny-server/pgx-adr-calculator;\n\
    log_dir /var/log/shiny-server;\n\
    directory_index on;\n\
  }\n\
}' > /etc/shiny-server/shiny-server.conf

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
