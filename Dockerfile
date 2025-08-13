# ---- Minimal, fast base ----
    FROM python:3.11-slim

    # Quality-of-life tools (tiny) — no heavy deps
    RUN apt-get update && apt-get install -y --no-install-recommends \
        bash coreutils procps \
     && rm -rf /var/lib/apt/lists/*
    
    ENV PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1
    
    # Install Bvalcalc from PyPI (use the TestPyPI variant below if needed)
    RUN pip install --no-cache-dir \
    --index-url https://test.pypi.org/simple \
    --extra-index-url https://pypi.org/simple \
    Bvalcalc  
    
    # Workdir for your mounted files
    WORKDIR /work
    
    # Default entrypoint runs the CLI
    ENTRYPOINT ["Bvalcalc"]
    