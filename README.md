# HyperC: Scalable Distributed Workflow for Long Reads Self-Correction

HyperC is a high-performance, distributed workflow designed to accelerate long-read self-correction pipelines. It integrates with existing tools like [CONSENT](https://github.com/paolap/mg-consent) to efficiently process third-generation sequencing (TGS) data at scale using a hybrid MPI + OpenMP parallelization strategy.

Developed to address the limitations of slow correction processes in genomic studies, HyperC enables fast and scalable analysis suitable for population-scale datasets.

---

## Features

- **Distributed parallelization** using MPI across compute nodes.
- **Intra-node multithreading** using OpenMP.
- **Optimized I/O**: compressed FASTQ loading (via Zstd) and balanced PAF partitioning.
- **Modular design**: easily plug in any MSA-based correction tool.
- **No need for parallel programming knowledge** to run.

---

## Use Case

Using real datasets (like NA12878 from the Nanopore WGS Consortium), it significantly reduced runtime, scaling efficiently across multiple nodes.

---

## Repository Structure

- `main.cpp` – Entry point for the HyperC workflow.
- `CMakeLists.txt` – Build configuration for the project.
- `compile.sh`, `compile_and_run.sh` – Scripts for compilation and execution.
- `consent.h` – Interface for integrating CONSENT's correction module (an example of correction module).
- `robin_hood.h` – Efficient hash map implementation used for data structures (needed for consent.h code).
- `utils.h` – Utility functions for data handling and decompression.
- `.gitmodules` – Contains references to submodules:
  - `spoa`: MSA library based on partial order alignment.
  - `Complete-Striped-Smith-Waterman-Library`: Optimized Smith-Waterman alignment.

---

## Installation

### Prerequisites

- GCC ≥ 10 with OpenMP support
- MPI (e.g. OpenMPI)
- CMake ≥ 3.10
- Zstandard (`libzstd`)

### Clone the Repository

Make sure to clone the repository **with submodules** to include required dependencies:

```bash
git clone --recurse-submodules https://github.com/riccardoc95/hipesc.git
```

### Build

```bash
./compile.sh
```

Or, for building and running:

```bash
./compile_and_run.sh
```

---

## Running the Pipeline

To run the pipeline, you need:

* A **FASTQ** file with long reads
* A corresponding **PAF** file (generated via [Minimap2](https://github.com/lh3/minimap2))

### Example (SLURM):

```bash
srun -n 4 ./hyperc /path/to/input.fq /path/to/input.paf
```

Each MPI rank will process a portion of the input, and correction results will be written to separate output files.

---

## Plug in Your Own Correction Module

To integrate another MSA-based correction module:

1. Implement a C/C++ function that accepts a target read and its overlapping reads.
2. Modify `consent.h` to wrap the new module.
3. Recompile with `compile.sh`.

HyperC will handle all job distribution and parallelization transparently.

---

## License

[MIT License](LICENSE)

---

## 🔗 Reference

If you use HyperC in your research, please cite:

> Di Rocco, L., Ferraro Petrillo, U., Ceccaroni, R., & Brutti, P. (202*). *A Scalable Distributed Workflow for Accelerating Long Reads Self-Correction*.

