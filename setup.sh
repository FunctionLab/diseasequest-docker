#!/bin/bash

# Build Docker image
docker build -t dq .

# Create output directory
mkdir -p outputs
chmod 777 -R outputs

# Download diseaseQUEST data
mkdir -p data
echo "Downloading test data file..."
wget --progress=bar:force -O data/dq_test_data.tar.gz http://wisp.princeton.edu/static/data/dq_test_data.tar.gz
tar xzvf data/dq_test_data.tar.gz -C data/ --strip-components=1
