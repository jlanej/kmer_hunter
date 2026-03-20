# ── Build stage: install Python deps ────────────────────────────────────────
FROM python:3.11-slim AS builder

RUN pip install --no-cache-dir --upgrade pip \
 && pip install --no-cache-dir \
      plotly==5.22.0 \
      pandas==2.2.2 \
      numpy==1.26.4

# ── Runtime image ────────────────────────────────────────────────────────────
FROM python:3.11-slim

# Install only runtime system libraries (bwa for fast exact-match search; curl for optional reference download)
RUN apt-get update \
 && apt-get install -y --no-install-recommends bwa ca-certificates curl samtools \
 && rm -rf /var/lib/apt/lists/*

# Copy installed Python packages from builder
COPY --from=builder /usr/local/lib/python3.11/site-packages \
                    /usr/local/lib/python3.11/site-packages

# Copy the kmer_hunter script
COPY kmer_hunter.py /usr/local/bin/kmer_hunter.py
RUN chmod +x /usr/local/bin/kmer_hunter.py

# Default cache location (override with --cache-dir or mount a volume)
ENV KMER_HUNTER_CACHE=/cache
RUN mkdir -p /cache /data

WORKDIR /data

ENTRYPOINT ["python", "/usr/local/bin/kmer_hunter.py"]
CMD ["--help"]
