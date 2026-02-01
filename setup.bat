@echo off
REM GPCR B1 Family Pipeline - Automated Setup Script for Windows
REM This script sets up the complete environment and installs all dependencies

setlocal enabledelayedexpansion

cls
echo.
echo ╔════════════════════════════════════════════════════════════╗
echo ║  🚀 GPCR B1 Family Pipeline - Automatic Setup             ║
echo ║  Version 1.0                                              ║
echo ╚════════════════════════════════════════════════════════════╝
echo.

REM Check if Python is installed
echo 📋 Checking for Python...
python --version >nul 2>&1
if errorlevel 1 (
    echo ❌ Python is not installed!
    echo.
    echo    Please install Python 3.7+ from https://www.python.org/downloads/
    echo    IMPORTANT: Check "Add Python to PATH" during installation
    echo.
    pause
    exit /b 1
)

for /f "tokens=2" %%i in ('python --version 2^>^&1') do set PYTHON_VERSION=%%i
echo ✅ Found: Python %PYTHON_VERSION%
echo.

REM Create virtual environment
echo 🔧 Creating virtual environment...
if exist "venv" (
    echo    ⚠️  venv folder already exists. Skipping creation...
) else (
    python -m venv venv
    if errorlevel 1 (
        echo ❌ Failed to create virtual environment
        pause
        exit /b 1
    )
    echo ✅ Virtual environment created
)
echo.

REM Activate virtual environment
echo 🔓 Activating virtual environment...
call venv\Scripts\activate.bat
if errorlevel 1 (
    echo ❌ Failed to activate virtual environment
    pause
    exit /b 1
)
echo ✅ Virtual environment activated
echo.

REM Upgrade pip
echo 📦 Upgrading pip...
python -m pip install --upgrade pip --quiet
echo ✅ pip upgraded
echo.

REM Install dependencies
echo 📚 Installing dependencies from requirements.txt...
if exist "requirements.txt" (
    pip install -r requirements.txt
    if errorlevel 1 (
        echo ❌ Failed to install dependencies
        pause
        exit /b 1
    )
    echo ✅ All dependencies installed successfully
) else (
    echo ❌ requirements.txt not found!
    pause
    exit /b 1
)
echo.

REM Verify installation
echo 🔍 Verifying installation...
echo.

REM Check BioPython
python -c "import Bio; print(f'✅ BioPython {Bio.__version__} installed')" 2>nul
if errorlevel 1 (
    echo ❌ BioPython verification failed
    pause
    exit /b 1
)

REM Check input file
if exist "Galaxy14_PHMMER.txt" (
    echo ✅ Input file Galaxy14_PHMMER.txt found
) else (
    echo ⚠️  Input file Galaxy14_PHMMER.txt not found
    echo    You'll need to provide this file before running the pipeline
)

REM Check script
if exist "master_pipeline_gpcr.py" (
    echo ✅ Pipeline script master_pipeline_gpcr.py found
) else (
    echo ❌ Pipeline script master_pipeline_gpcr.py not found
    pause
    exit /b 1
)

echo.
echo ╔════════════════════════════════════════════════════════════╗
echo ║  ✅ Setup Complete!                                        ║
echo ╚════════════════════════════════════════════════════════════╝
echo.
echo 📝 Next steps:
echo.
echo 1. Activate the virtual environment (next time you open command prompt):
echo    ^> venv\Scripts\activate
echo.
echo 2. Edit the script and set your email (REQUIRED):
echo    Open master_pipeline_gpcr.py and change:
echo    Entrez.email = 'your_email@example.com'
echo.
echo 3. Run the pipeline:
echo    ^> python master_pipeline_gpcr.py
echo.
echo 4. Find results in FINAL_OUTPUT\ folder
echo.
echo 📚 For more information, see:
echo    - README.md - Full documentation
echo    - SETUP_GUIDE.md - Detailed setup instructions
echo.
echo Good luck! 🧬
echo.
pause
