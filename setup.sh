#!/bin/bash

# GPCR B1 Family Pipeline - Automated Setup Script for Linux/macOS
# This script sets up the complete environment and installs all dependencies

set -e  # Exit on any error

echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║  🚀 GPCR B1 Family Pipeline - Automatic Setup             ║"
echo "║  Version 1.0                                              ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

# Check if Python is installed
echo "📋 Checking for Python..."
if ! command -v python3 &> /dev/null; then
    echo "❌ Python 3 is not installed!"
    echo "   Please install Python 3.7+ from https://www.python.org/downloads/"
    exit 1
fi

PYTHON_VERSION=$(python3 --version 2>&1)
echo "✅ Found: $PYTHON_VERSION"
echo ""

# Create virtual environment
echo "🔧 Creating virtual environment..."
if [ -d "venv" ]; then
    echo "   ⚠️  venv folder already exists. Skipping creation..."
else
    python3 -m venv venv
    echo "✅ Virtual environment created"
fi
echo ""

# Activate virtual environment
echo "🔓 Activating virtual environment..."
source venv/bin/activate
echo "✅ Virtual environment activated"
echo ""

# Upgrade pip
echo "📦 Upgrading pip..."
pip install --upgrade pip --quiet
echo "✅ pip upgraded"
echo ""

# Install dependencies
echo "📚 Installing dependencies from requirements.txt..."
if [ -f "requirements.txt" ]; then
    pip install -r requirements.txt
    echo "✅ All dependencies installed successfully"
else
    echo "❌ requirements.txt not found!"
    exit 1
fi
echo ""

# Verify installation
echo "🔍 Verifying installation..."
echo ""

# Check BioPython
python -c "import Bio; print(f'✅ BioPython {Bio.__version__} installed')" || {
    echo "❌ BioPython verification failed"
    exit 1
}

# Check input file
if [ -f "Galaxy14_PHMMER.txt" ]; then
    echo "✅ Input file Galaxy14_PHMMER.txt found"
else
    echo "⚠️  Input file Galaxy14_PHMMER.txt not found"
    echo "   You'll need to provide this file before running the pipeline"
fi

# Check script
if [ -f "master_pipeline_gpcr.py" ]; then
    echo "✅ Pipeline script master_pipeline_gpcr.py found"
else
    echo "❌ Pipeline script master_pipeline_gpcr.py not found"
    exit 1
fi

echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║  ✅ Setup Complete!                                        ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "📝 Next steps:"
echo ""
echo "1. Activate the virtual environment (next time you open terminal):"
echo "   $ source venv/bin/activate"
echo ""
echo "2. Edit the script and set your email (REQUIRED):"
echo "   Open master_pipeline_gpcr.py and change:"
echo "   Entrez.email = 'your_email@example.com'"
echo ""
echo "3. Run the pipeline:"
echo "   $ python master_pipeline_gpcr.py"
echo ""
echo "4. Find results in FINAL_OUTPUT/ folder"
echo ""
echo "📚 For more information, see:"
echo "   - README.md - Full documentation"
echo "   - SETUP_GUIDE.md - Detailed setup instructions"
echo ""
echo "Good luck! 🧬"
echo ""
