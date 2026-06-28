# Use an official lightweight Python engine baseline image
FROM python:3.11-slim-bookworm

# Establish root administrative privileges to install native biological binary systems
USER root

# Install system utilities along with the required compiled biological programs
RUN apt-get update && apt-get install -y --no-not-get --no-install-recommends \
    build-essential \
    ncbi-blast+ \
    prodigal \
    git \
    && rm -rf /var/lib/apt/lists/*

# Set the operational directory inside the virtual system container
WORKDIR /app

# Create a secure, unprivileged system user account for Hugging Face compatibility
RUN useradd -m -u 1000 user
RUN chown -R user:user /app

# Switch context execution permanently to the unprivileged account space
USER user
ENV HOME=/home/user \
    PATH=/home/user/.local/bin:$PATH

# Copy over the Python library asset manifest record
COPY --chown=user:user requirements.txt .

# Compile and install Python tools inside user-space boundaries
RUN pip install --no-cache-dir --user -r requirements.txt

# Migrate your local source directory code completely into the application workspace
COPY --chown=user:user . .

# Expose Streamlit's default network communication gateway port
EXPOSE 8501

# Run the containerized interactive web application engine on startup
ENTRYPOINT ["streamlit", "run", "FishAMR_Link_v15_Supreme.py", "--server.port=8501", "--server.address=0.0.0.0"]
