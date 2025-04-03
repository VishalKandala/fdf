# **FDF Swarm**

---

## **Table of Contents**

- [Overview](#overview)
- [Directory Structure](#directory-structure)
- [Dependencies](#dependencies)
- [Installation and Building](#installation-and-building)
  - [1. Setting Up PETSc](#1-setting-up-petsc)
  - [2. Cloning the Repository](#2-cloning-the-repository)
  - [3. Building the Project](#3-building-the-project)
- [Usage](#usage)
  - [Running Executables](#running-executables)
  - [Available Executables](#available-executables)
- [Makefile Targets](#makefile-targets)
  - [Default Target](#default-target)
  - [Building Specific Executables](#building-specific-executables)
  - [Custom Targets](#custom-targets)
- [Scripts](#scripts)
- [Generating TAGS File](#generating-tags-file)
- [Cleaning Up](#cleaning-up)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)
- [Contact Information](#contact-information)

---

## **Overview**

This project is a computational framework that utilizes the PETSc library for solving complex scientific simulations. The primary focus of the project is particle swarm management, grid interpolation, and walking search algorithms.

---

## **Directory Structure**

```
.
├── include/            # Header files for modularized functionality
│   ├── grid.h
│   ├── ParticleSwarm.h
│   ├── walkingsearch.h
│   ├── logging.h
│   ├── variables.h
├── src/                # Source files implementing project logic
│   ├── grid.c
│   ├── ParticleSwarm.c
│   ├── walkingsearch.c
│   ├── logging.c
├── scripts/            # Preprocessing Python scripts
├── bin/                # Compiled executables
├── obj/                # Object files generated during compilation
├── Makefile            # Build automation file
├── README.md           # Project documentation
```

---

## **Dependencies**

- PETSc (Portable, Extensible Toolkit for Scientific Computation)
- GCC (GNU Compiler Collection)
- MPI (Message Passing Interface)

---

## **Installation and Building**

### 1. Setting Up PETSc

Ensure PETSc is installed and configured:
```bash
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-linux-c-opt
```

### 2. Cloning the Repository

Clone this repository into your local directory:
```bash
git clone https://github.com/your-repo/project-name.git
cd project-name
```

### 3. Building the Project

Use the provided `Makefile` to compile the project:
```bash
make
```
This will generate executables in the `bin/` directory.

---

## **Usage**

### Running Executables

Run the primary executable `swarm_interp`:
```bash
bin/swarm_interp
```

### Available Executables

- **`swarm_interp`**: Main simulation executable.
- Additional executables can be found in the `bin/` directory.

---

## **Makefile Targets**

### Default Target
The default target builds all executables:
```bash
make
```

### Building Specific Executables
To build a specific target, run:
```bash
make swarm_interp
```

### Custom Targets

- **`make clean`**: Removes object files.
- **`make clean_all`**: Cleans everything, including executables and object files.
- **`make tags`**: Generates a TAGS file for code navigation.

---

## **Scripts**

Preprocessing scripts are located in the `scripts/` directory. Run them as:
```bash
python3 scripts/your_script.py
```

---

## **Generating TAGS File**

Use the `make tags` command to generate a TAGS file for code navigation:
```bash
make tags
```

---

## **Cleaning Up**

To clean object files and executables, use:
```bash
make clean
```

For a complete cleanup:
```bash
make clean_all
```

---

## **Troubleshooting**

- Ensure PETSc is correctly configured and environment variables `PETSC_DIR` and `PETSC_ARCH` are set.
- Use `make clean_all` if compilation issues persist and rebuild.

---

## **Contributing**

Contributions are welcome! Please follow the standard Git workflow:
1. Fork the repository.
2. Create a feature branch.
3. Submit a pull request.

---

## **License**

This project is licensed under the MIT License. See the LICENSE file for details.

---

## **Contact Information**

For issues or inquiries, contact:
- Email: vishalkandala@tamu.edu
- GitHub: https://github.com/VishalKandala
*/


1.  **Storage IS Node-Indexed:** The memory allocation for `ucat` uses the `fda` DMDA, which is based on the **node** indices `(i,j,k)` ranging from `0..IM`, `0..JM`, `0..KM`. So, `ucat[k][j][i]` accesses the data slot associated with **node (i,j,k)**. This part is unambiguous from the allocation (`VecDuplicate` from `Csi` based on `fda`).

2.  **Value *Stored* Represents Cell Center (for Interior Indices):**
    *   Look at `Contra2Cart`. It calculates the Cartesian velocity `(det0/det, etc.)` by averaging fluxes and metrics surrounding **cell `C(i,j,k)`**.
    *   It then **stores** this calculated value, which physically represents the velocity at the center of cell `C(i,j,k)`, into the array slot `ucat[k][j][i]`.
    *   **The Convention:** The code adopts the convention that the *value conceptually belonging to cell C(i,j,k)* is stored in the memory slot associated with *node index (i,j,k)*, but **only when `(i,j,k)` corresponds to an interior cell** (roughly indices `1` to `IM-1`, etc., matching the loop bounds `lxs..lxe` in `Contra2Cart`). The node `(i,j,k)` acts as the "anchor" or "marker" storage location for the value associated with cell `C(i,j,k)`.

3.  **Why this Convention?**
    *   It keeps the `ucat` data structure aligned with `coor` and `ucont` (all using `fda`).
    *   When stencils in `Convection` or `Viscous` access neighbors like `ucat[k][j][i+1]` while operating relative to cell `C(i,j,k)`, they correctly retrieve the value *representing* cell `C(i+1,j,k)` because that value was stored at node index `(i+1,j,k)`.

4.  **Boundary Indices Store Something Else:**
    *   As established, `Contra2Cart` does *not* calculate values for boundary indices `i=0, IM`, `j=0, JM`, `k=0, KM`.
    *   `FormBCS` explicitly populates `ucat[k][j][0]`, `ucat[k][j][IM]`, etc., using extrapolation formulas. These values stored at the boundary node indices **do not represent** the physical velocity at the center of the adjacent boundary cells (`C(0,j,k)` or `C(IM-1,j,k)`). They represent the **extrapolated ghost values** needed by the computational stencils operating in the first/last interior cells.

**Resolving the Contradiction:**

*   **ucat IS stored using node indexing `0..IM`.**
*   **For INTERIOR node indices `i` (from `1` to `IM-1`, etc.):** The value `ucat[k][j][i]` *physically represents* the state variable (Cartesian velocity) associated with the **center of cell `C(i,j,k)`**.
*   **For BOUNDARY node indices `i` (`0` and `IM`, etc.):** The value `ucat[k][j][i]` stores an **extrapolated "ghost" value** determined by the boundary condition routines (`FormBCS`), not the physical state at the center of cell `C(0,j,k)` or `C(IM-1,j,k)`.

**Revised Table:**

| Conceptual Cell                 | Cell Index `(i,j,k)` Range        | Type      | `ucat` Storage Index `[k][j][i]` Used | Physical Meaning of Value Stored in `ucat[k][j][i]`             | `P`, `nvert` Storage Index `[k][j][i]` |
| :------------------------------ | :-------------------------------- | :-------- | :------------------------------------ | :---------------------------------------------------------- | :--------------------------------------- |
| **Interior Cells**              | `1 <= i <= IM-1` <br> `1 <= j <= JM-1` <br> `1 <= k <= KM-1` | Interior  | `[k][j][i]`                         | Velocity representing center of **Cell `C(i,j,k)`**           | `[k][j][i]`                              |
| **Boundary Cells (e.g., i=0)**  | `i=0` <br> `1 <= j <= JM-1` <br> `1 <= k <= KM-1` | Boundary /<br>Ghost Role | `[k][j][0]` *(Node Index)*          | **Extrapolated "ghost" value** based on BC at `i=0` face     | `[k][j][0]`                              |
| **Boundary Cells (e.g., i=IM-1)**| `i=IM-1` <br> `1 <= j <= JM-1` <br> `1 <= k <= KM-1` | Boundary /<br>Ghost Role | `[k][j][IM-1]` *(Node Index)*     | Velocity representing center of **Cell `C(IM-1,j,k)`**      | `[k][j][IM-1]`                           |
| *Boundary Nodes (e.g., i=IM)*   | *(N/A - Not a Cell)*              | Node      | `[k][j][IM]` *(Node Index)*         | **Extrapolated "ghost" value** based on BC at `i=IM` face | *(N/A)*                                  |
| **Edge Cells (e.g., i=0, j=0)** | `i=0, j=0` <br> `1 <= k <= KM-1`   | Boundary /<br>Ghost Role | `[k][0][0]` *(Node Index)*          | **Extrapolated "ghost" value** based on BCs at `i=0, j=0` | `[k][0][0]`                              |
| **Corner Cells (e.g., i=0, j=0, k=0)** | `i=0, j=0, k=0`                | Boundary /<br>Ghost Role | `[0][0][0]` *(Node Index)*          | **Extrapolated "ghost" value** based on BCs at `i=0, j=0, k=0`| `[0][0][0]`                              |

This refined table clarifies that the *meaning* of the value stored in `ucat[k][j][i]` depends on whether `(i,j,k)` is an interior or boundary node index.
