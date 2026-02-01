# 🚀 Complete Setup Guide - GPCR B1 Family Pipeline

**This guide will walk you through EVERYTHING step-by-step to get the pipeline running on your machine.**

## Table of Contents
1. [System Requirements](#system-requirements)
2. [Installation by Operating System](#installation-by-operating-system)
3. [Verification Checklist](#verification-checklist)
4. [Common Issues & Solutions](#common-issues--solutions)

---

## System Requirements

Before you start, make sure you have:
- ✅ **Python 3.7 or higher** installed
- ✅ **Internet connection** (for NCBI database access)
- ✅ **Git** (optional, but recommended)
- ✅ **A text editor** (VS Code, Notepad++, etc.)

### Check Your Python Version

Open a terminal/command prompt and type:

```bash
python --version
```

or

```bash
python3 --version
```

You should see something like:
```
Python 3.12.3
```

❌ **If you don't have Python**: Download from https://www.python.org/downloads/
- Make sure to check **"Add Python to PATH"** during installation!

---

## Installation by Operating System

### 🐧 **LINUX / macOS**

#### Step 1: Clone or Download the Repository

**Option A: Using Git (Recommended)**
```bash
git clone https://github.com/JB5000/gpcr-b1family-pipeline.git
cd gpcr-b1family-pipeline
```

**Option B: Manual Download**
1. Go to https://github.com/JB5000/gpcr-b1family-pipeline
2. Click green "Code" button → "Download ZIP"
3. Extract the ZIP file
4. Open terminal and navigate to the folder:
```bash
cd ~/Downloads/gpcr-b1family-pipeline
```

#### Step 2: Create Virtual Environment

A **virtual environment** is a isolated Python workspace for this project. It keeps everything clean and organized.

```bash
python3 -m venv venv
```

This creates a folder called `venv` with all the isolated Python files.

#### Step 3: Activate Virtual Environment

```bash
source venv/bin/activate
```

✅ **Success indicator**: Your terminal prompt should now show `(venv)` at the beginning:
```
(venv) user@computer:~/gpcr-b1family-pipeline$
```

#### Step 4: Install Dependencies

```bash
pip install -r requirements.txt
```

This command reads the `requirements.txt` file and installs **BioPython** automatically.

**Expected output:**
```
Collecting biopython>=1.81
  Downloading biopython-1.84-cp312-cp312-linux_x86_64.whl ...
Installing collected packages: biopython
Successfully installed biopython-1.84
```

✅ **Done!** Your environment is ready.

#### Step 5: Deactivate Virtual Environment (When done)

```bash
deactivate
```

Your terminal prompt will return to normal (no `(venv)` prefix).

---

### 🪟 **WINDOWS**

#### Step 1: Clone or Download the Repository

**Option A: Using Git**
```bash
git clone https://github.com/JB5000/gpcr-b1family-pipeline.git
cd gpcr-b1family-pipeline
```

**Option B: Manual Download**
1. Go to https://github.com/JB5000/gpcr-b1family-pipeline
2. Click green "Code" button → "Download ZIP"
3. Extract the ZIP file to your preferred location
4. Open **Command Prompt** (cmd.exe) or **PowerShell** and navigate:
```bash
cd C:\Users\YourUsername\Downloads\gpcr-b1family-pipeline
```

#### Step 2: Create Virtual Environment

```bash
python -m venv venv
```

This creates a folder called `venv`.

#### Step 3: Activate Virtual Environment

**Command Prompt (cmd.exe):**
```bash
venv\Scripts\activate
```

**PowerShell:**
```bash
venv\Scripts\Activate.ps1
```

✅ **Success indicator**: Your prompt should show `(venv)`:
```
(venv) C:\Users\YourUsername\gpcr-b1family-pipeline>
```

#### Step 4: Install Dependencies

```bash
pip install -r requirements.txt
```

**Expected output:**
```
Collecting biopython>=1.81
  Downloading biopython-1.84-cp312-win_amd64.whl ...
Installing collected packages: biopython
Successfully installed biopython-1.84
```

✅ **Done!** Your environment is ready.

#### Step 5: Deactivate Virtual Environment (When done)

```bash
deactivate
```

Your prompt returns to normal.

---

### 🍎 **macOS (Apple Silicon M1/M2/M3)**

⚠️ **Special note**: Apple Silicon Macs may have specific requirements.

#### Step 1-2: Same as Linux/macOS above

```bash
git clone https://github.com/JB5000/gpcr-b1family-pipeline.git
cd gpcr-b1family-pipeline
python3 -m venv venv
```

#### Step 3: Activate & Install

```bash
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

If you encounter issues with BioPython, try:
```bash
pip install --pre biopython
```

---

## Verification Checklist

### ✅ Test Your Installation

After completing the steps above, **verify everything works**:

#### 1. Check Python Version
```bash
python --version
```
Expected: `Python 3.7+`

#### 2. Check Virtual Environment is Active
Your terminal should show `(venv)` prefix.

#### 3. Check BioPython is Installed
```bash
python -c "import Bio; print(Bio.__version__)"
```
Expected output: Something like `1.84`

#### 4. Check Input File Exists
```bash
ls Galaxy14_PHMMER.txt
```
or on Windows:
```bash
dir Galaxy14_PHMMER.txt
```
Expected: File listed (no error)

#### 5. Run a Test
```bash
python master_pipeline_gpcr.py
```

You should see:
```
Total unique accessions: 240

[10%] 24/240 | elapsed: 0:00:50 | ETA: 0:07:33
...
PIPELINE COMPLETE ✅
Results written to: FINAL_OUTPUT/
```

**If all these work**, you're ready to use the pipeline! 🎉

---

## Common Issues & Solutions

### ❌ "python: command not found"

**Problem**: Python is not installed or not in PATH.

**Solution**:
1. Download Python: https://www.python.org/downloads/
2. **IMPORTANT**: During installation, check ✅ **"Add Python to PATH"**
3. Restart your terminal/command prompt
4. Try `python --version` again

---

### ❌ "ModuleNotFoundError: No module named 'Bio'"

**Problem**: BioPython is not installed.

**Solution**:
1. Make sure virtual environment is **activated** (you should see `(venv)` in prompt)
2. Run:
```bash
pip install -r requirements.txt
```

If that doesn't work:
```bash
pip install biopython
```

---

### ❌ "No such file or directory: venv/bin/activate"

**Problem**: Virtual environment wasn't created properly.

**Solution**:
1. Delete the `venv` folder (if it exists)
2. Recreate it:
```bash
python3 -m venv venv
```

3. Reactivate:
```bash
source venv/bin/activate  # Linux/macOS
# or
venv\Scripts\activate     # Windows
```

---

### ❌ "No such file or directory: Galaxy14_PHMMER.txt"

**Problem**: Input file is missing or you're in the wrong folder.

**Solution**:
1. Make sure you're in the correct directory:
```bash
ls                    # Linux/macOS - lists files
dir                   # Windows - lists files
```

2. You should see `Galaxy14_PHMMER.txt`, `master_pipeline_gpcr.py`, `README.md`, etc.

3. If you don't see the input file, download it from GitHub or copy it to this directory.

---

### ❌ "NCBI connection timeout"

**Problem**: Network issues or NCBI is temporarily down.

**Solution**:
1. Check your internet connection
2. Try again in a few minutes
3. In `master_pipeline_gpcr.py`, increase the `SLEEP` value to be more respectful of NCBI:
```python
SLEEP = 0.5  # Instead of 0.34
```

---

### ❌ "ModuleNotFoundError: biopython (on Windows with PowerShell)"

**Problem**: PowerShell execution policy may prevent activation.

**Solution**: Use Command Prompt instead:
```bash
cmd.exe
cd gpcr-b1family-pipeline
venv\Scripts\activate
pip install -r requirements.txt
```

---

## Quick Start Scripts

### 🐧 Linux/macOS Automatic Setup

Create a file called `setup.sh`:

```bash
#!/bin/bash
echo "🚀 Setting up GPCR B1 Family Pipeline..."
echo ""
echo "✓ Creating virtual environment..."
python3 -m venv venv

echo "✓ Activating virtual environment..."
source venv/bin/activate

echo "✓ Installing dependencies..."
pip install --upgrade pip
pip install -r requirements.txt

echo ""
echo "✅ Setup complete!"
echo ""
echo "To activate the environment in the future, run:"
echo "  source venv/bin/activate"
echo ""
echo "To run the pipeline:"
echo "  python master_pipeline_gpcr.py"
```

Then run:
```bash
chmod +x setup.sh
./setup.sh
```

---

### 🪟 Windows Automatic Setup

Create a file called `setup.bat`:

```batch
@echo off
echo.
echo 🚀 Setting up GPCR B1 Family Pipeline...
echo.
echo ✓ Creating virtual environment...
python -m venv venv

echo ✓ Activating virtual environment...
call venv\Scripts\activate

echo ✓ Installing dependencies...
pip install --upgrade pip
pip install -r requirements.txt

echo.
echo ✅ Setup complete!
echo.
echo To activate the environment in the future, run:
echo   venv\Scripts\activate
echo.
echo To run the pipeline:
echo   python master_pipeline_gpcr.py
echo.
pause
```

Then double-click `setup.bat` or run:
```bash
setup.bat
```

---

## Understanding Virtual Environments

### 💡 What is a Virtual Environment?

Think of it like a **separate Python installation** just for this project:
- **Without venv**: All Python packages go to one system-wide location → potential conflicts
- **With venv**: Each project has its own isolated packages → clean, organized, safe

### 🎯 Why Use Virtual Environments?

✅ **Project Isolation**: Different projects can use different package versions  
✅ **Clean System**: Doesn't mess with your system Python  
✅ **Reproducibility**: Same setup works on any computer  
✅ **Collaboration**: Easy to share requirements.txt with teammates  

### 📋 What's in requirements.txt?

The `requirements.txt` file lists all packages needed:
```
biopython>=1.81
```

This means: "Install BioPython version 1.81 or higher"

When someone runs `pip install -r requirements.txt`, all these packages are automatically installed.

---

## Troubleshooting Checklist

Before asking for help, check:

- [ ] Python 3.7+ is installed (`python --version`)
- [ ] You're in the correct folder (can you see Galaxy14_PHMMER.txt?)
- [ ] Virtual environment is created (`venv` folder exists)
- [ ] Virtual environment is activated (`(venv)` shows in prompt)
- [ ] Dependencies are installed (ran `pip install -r requirements.txt`)
- [ ] Internet is connected (needed for NCBI queries)
- [ ] You configured Entrez.email in the script

---

## Getting Help

If you run into issues:

1. **Check the README.md** - It has general information
2. **Check this guide** - It has common issues
3. **Check Python/BioPython docs**:
   - https://www.python.org/
   - https://biopython.org/
4. **Open a GitHub Issue**: https://github.com/JB5000/gpcr-b1family-pipeline/issues

---

## Next Steps

Once your setup is complete:

1. **Read the README.md** for full documentation
2. **Edit master_pipeline_gpcr.py** to set your email:
   ```python
   Entrez.email = "your_email@example.com"  # CHANGE THIS!
   ```
3. **Run the pipeline**:
   ```bash
   python master_pipeline_gpcr.py
   ```
4. **Check the output** in `FINAL_OUTPUT/` folder

---

**Happy bioinformatics! 🧬**

Last Updated: February 2026
