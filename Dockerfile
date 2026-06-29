FROM python:3.11-slim-bookworm

# Install standard dependencies, curl/wget, and compilation utilities
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ncbi-blast+ \
    prodigal \
    git \
    wget \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Manually pull, extract, and deploy the authentic pre-compiled NCBI AMRFinderPlus Linux binary
RUN wget https://github.com/ncbi/amr/releases/download/amrfinder_v3.12.8/amrfinder_binaries_v3.12.8.tar.gz \
    && tar -xzf amrfinder_binaries_v3.12.8.tar.gz \
    && mv amrfinder /usr/local/bin/amrfinder \
    && mv amrfinder_update /usr/local/bin/amrfinder_update \
    && rm amrfinder_binaries_v3.12.8.tar.gz

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

# Expose the mandatory Hugging Face Space port
EXPOSE 7860

# Use explicit path parameters and absolute execution definitions to prevent shell blocking
ENTRYPOINT ["sh", "-c", "mkdir -p /app/amr_db && amrfinder_update -d /app/amr_db || true; streamlit run FishAMR_Link_v15_Supreme.py --server.port=7860 --server.address=0.0.0.0"]

