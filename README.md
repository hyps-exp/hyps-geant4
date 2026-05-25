<!--
# hyps-geant4

A Geant4-based Monte Carlo simulation and ROOT-based offline analyzer for the HYPS experiment.
> **Note:**
> This codebase originates from a previous experiment (E40) and contains legacy code.
> Major update and refactoring are expected in the near future. *(As of May 2026)*

## Overview
This repository contains:
- **simulation** - A Geant4-based Monte Calro simulation that generates output in `.dat` format.
- **analyzer** - A ROOT-based offline analyzer that reads simulation output (`.dat`) and performs physics analysis.

## Requirements
This software is designed to run on the KEKCC, where Geant4 and ROOT are available as pre-installed moduels.
It has not been tested on other environments.

### simulation
| Dependency | Version |
|---|---|
| Geant4 | 11.2.2 |
| CMake | 3.16 (minimum required) |
| GCC | 8 or later (C++17 required) |

### analyzer
| Dependency | Version |
|---|---|
| ROOT | 6.32.04 (tested) |
| Geant4 | 11.2.2 (for bundled CLHEP) |
| Make | 4.3 (tested)|
| GCC | 11.5.0 (tested) |

---

## Installation
### PATH settings
Please add following line in your `.bashrc` on KEKCC.
```bash
. /sw/packages/geant4/11.2.2/bin/geant4.sh
```


## Repository Structure
```
|--- simulation/
     |--- CMakeLists.txt
     |--- build.sh
     |--- clean.sh
     |--- main.cc
     |--- include/
     |--- src/
     |--- init_vis.mac # Macro for visualization initialization
     |--- vis.mac      # Macro for visualization settings (called by init_vis.mac)
     |--- run.mac      # Macro for batch mode execution
     |--- script/      # Scripts for job submission
     |--- tools/       # Utilities for geometry calculation (not used in HYPS)

|--- analyzer/
     |---
```
-->

# hyps-geant4

A Geant4-based Monte Carlo simulator and ROOT-based offline analyzer for the HYPS experiment.

> **Note:**
> This codebase originates from a previous experiment (E40) and contains legacy code.
> Major updates and refactoring are expected in the near future. *(As of May 2026)*

---

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Getting Started](#getting-started)
  - [Prerequisites](#prerequisites)
  - [Build](#build)
  - [Run](#run)
- [Usage](#usage)
- [Repository Structure](#repository-structure)

---

## Overview

This repository contains:

- **`simulation/`** — A Geant4-based Monte Carlo simulator that generates output in `.dat` format.
- **`analyzer/`** — A ROOT-based offline analyzer that reads simulation output (`.dat`) and performs physics analysis.

---

## Requirements

### simulation

| Dependency | Version |
|---|---|
| Geant4 | 11.2.2 |
| CMake | 3.16 (minimum required) |
| GCC | 8 or later (C++17 required) |

### analyzer

| Dependency | Version |
|---|---|
| ROOT | 6.32.04 (tested) |
| Geant4 | 11.2.2 (for bundled CLHEP) |
| Make | 4.3 (tested) |
| GCC | 11.5.0 (tested) |

---

## Getting Started

This software is designed to run on the KEKCC, where Geant4 and ROOT are available as pre-installed modules.
It has not been tested on other environments.

### Prerequisites

Add the following line to your `~/.bashrc`:

```bash
. /sw/packages/geant4/11.2.2/bin/geant4.sh
```

This single line is sufficient to set up the environment. No other configuration is needed.

#### ⚠️ Warning: Conflicting Old Configuration

If you previously set up the HYPS Geant4 environment manually, your `~/.bashrc` may contain the following lines:

```bash
. /sw/packages/geant4/11.2.2/bin/geant4.sh
. /sw/packages/geant4/11.2.2/share/Geant4/geant4make/geant4make.sh
export CLHEPSYS=/home/had/miwaq/cern/clhep/install
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CLHEPSYS/lib
export PATH=$CLHEPSYS/bin:$PATH
# export CPLUS_INCLUDE_PATH=/usr/include/qt5:/usr/include/qt5/QtCore:$CPLUS_INCLUDE_PATH
```

**These entries must be removed or commented out before using the current setup.**

Leaving them in place can cause the following issues:

- **CLHEP version conflict** — The system Geant4 already bundles CLHEP internally. A user-local `CLHEPSYS` pointing to a separate installation will take precedence and may cause link errors or undefined behavior at runtime.
- **`LD_LIBRARY_PATH` conflict** — Appending a user-local CLHEP lib path can override shared libraries loaded by Geant4, leading to subtle crashes or mismatched symbol errors that are difficult to diagnose.
- **`CPLUS_INCLUDE_PATH` conflict** — This was intended for Qt5-based visualization. Since `geant4.sh` already handles Qt5 configuration, this setting is redundant and may cause the compiler to pick up incorrect Qt5 headers.

To fix, open `~/.bashrc` and remove or comment out (with `#`) all lines from the old block.

Note that simply running `source ~/.bashrc` is **not sufficient** — environment variables already exported in the current shell session will remain active. You must fully reset your environment:

- **Standard SSH session** — Disconnect and reconnect.
- **Inside a tmux session** — Kill the session and start a new one:

  ```bash
  tmux kill-session
  ```

  Then reconnect via SSH and create a new tmux session.

### Build

```bash
git clone https://github.com/hyps-exp/hyps-geant4.git
cd hyps-geant4/simulation
./build.sh
```

The executable will be placed in `simulation/bin/`. To rebuild after modifying source files, run `./build.sh` again.

If you want to clean the build artifacts (`bin/` and `.build/`):

```bash
./clean.sh
```

### Run

**Interactive (visualization) mode:**

```bash
./G4Hyps
```

**Batch mode:**

```bash
./G4Hyps run.mac
```

---

## Usage

For details on macro files, geometry configuration, and the analyzer workflow, see the Wiki.

---

## Repository Structure

```
hyps-geant4/
|--- simulation/
|    |--- CMakeLists.txt
|    |--- build.sh
|    |--- clean.sh
|    |--- main.cc
|    |--- include/
|    |--- src/
|    |--- init_vis.mac    # Macro for visualization initialization
|    |--- vis.mac         # Macro for visualization settings (called by init_vis.mac)
|    |--- run.mac         # Macro for batch mode execution
|    |--- script/         # Scripts for job submission
|    |--- tools/          # Utilities for geometry calculation (not used in HYPS)
|--- analyzer/
|    |--- script/         # Scripts for job submission (Condor, bsub)
|    |--- src/            # Source files
|--- param-sim/
     |--- CFTEFF/
     |--- DCGEO/
     |--- E40/
     |--- conf/
     |--- fieldmap/
     |--- misc/
     |--- util/
```