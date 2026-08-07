FROM python:3.12-slim

WORKDIR /home/genomics
COPY . /home/genomics

# Install the checked-out source. Differential analysis and normalization are
# implemented in Python and no longer require R or Bioconductor packages.
RUN pip install --no-cache-dir .
