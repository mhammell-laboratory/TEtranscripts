# syntax=docker/dockerfile:1

ARG PYTHON_VERSION=3.12

FROM python:${PYTHON_VERSION}-slim-bookworm AS builder

ENV PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /build
COPY . .

# Build application and dependency wheels once, leaving compilers and source
# files out of the runtime image.
RUN python -m pip wheel --wheel-dir /wheels .


FROM python:${PYTHON_VERSION}-slim-bookworm AS runtime

ARG VCS_REF=unknown

LABEL org.opencontainers.image.title="TEtranscripts" \
      org.opencontainers.image.description="Transposable-element quantification and Python-native differential analysis" \
      org.opencontainers.image.source="https://github.com/molikd/TEtranscripts" \
      org.opencontainers.image.revision="${VCS_REF}" \
      org.opencontainers.image.licenses="GPL-3.0-only"

ENV PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_NO_CACHE_DIR=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

COPY --from=builder /wheels /wheels
RUN python -m pip install /wheels/* \
    && rm -rf /wheels \
    && groupadd --gid 10001 tetranscripts \
    && useradd --uid 10001 --gid 10001 --create-home --shell /usr/sbin/nologin tetranscripts \
    && install -d --owner tetranscripts --group tetranscripts /data

WORKDIR /data
USER tetranscripts

ENTRYPOINT ["TEtranscripts"]
CMD ["--help"]


# Optional compatibility image for reproducing the historical R analysis or
# running legacy/native comparisons.  The default final image below continues
# to inherit directly from ``runtime`` and therefore remains R-free.
FROM runtime AS legacy

USER root
RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        r-base \
        r-bioc-deseq \
        r-bioc-deseq2 \
    && rm -rf /var/lib/apt/lists/*

LABEL org.opencontainers.image.description="TEtranscripts compatibility image with R, DESeq, and DESeq2" \
      org.opencontainers.image.variant="legacy-deseq"

USER tetranscripts


# Keep the no-target Docker build as the small Python-native production image.
FROM runtime AS final
