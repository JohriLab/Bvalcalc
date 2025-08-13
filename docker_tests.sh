docker run --rm -it \
  --memory=8g --memory-swap=8g \
  -v "$PWD":/work -w /work \
  python:3.11-slim \
  bash -lc "pip install --no-cache-dir --index-url https://test.pypi.org/simple --extra-index-url https://pypi.org/simple Bvalcalc && \
            Bvalcalc --gene --pop_params DrosophilaParams_phastcons.py  "
