FROM python:3.11-slim-bookworm

# Clean update and package installation with your complete bio-tool suite
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ncbi-blast+ \
    ncbi-amrfinderplus \
    prodigal \
    git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

# Expose the mandatory Hugging Face Space port
EXPOSE 7860

ENTRYPOINT ["streamlit", "run", "FishAMR_Link_v15_Supreme.py", "--server.port=7860", "--server.address=0.0.0.0"]
