#!/bin/bash

# Exit on any error
set -e

# Check if python3 is installed
if ! command -v python3 &> /dev/null; then
    echo "Python3 could not be found. Please install Python3."
    exit 1
fi

if command -v python3 &>/dev/null; then
    PYTHON_CMD=python3
elif command -v python &>/dev/null; then
    PYTHON_CMD=python

# Create a virtual environment
echo "Creating virtual environment..."
$PYTHON_CMD -m venv venv


# On Windows
#.\venv\Scripts\activate

# Install required packages
echo "Installing required packages..."
pip install --upgrade pip  # Upgrade pip to the latest version
pip install numpy matplotlib scipy pandas healpy astropy casatools casatasks ipykernel jupyter ipython tqdm


# Add the virtual environment as a new Jupyter kernel
$PYTHON_CMD -m ipykernel install --user --name=venv --display-name "Hirax Python env"


# Verification of the installation
echo "Verifying installation..."
pip list

echo "Installation complete!"



