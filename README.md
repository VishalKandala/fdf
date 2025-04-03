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

`ucont`, the face-centered contravariant velocity/flux, for the `IM=JM=KM=5` case.

**1. Total Size of `ucont`:**

*   Basis: Node-based `fda` DMDA (`M=6, N=6, P=6`).
*   Total node locations: `M * N * P = 6 * 6 * 6 = 216`.
*   Degrees of Freedom (DOF): 3 (for `.x`, `.y`, `.z` components stored at each node index).
*   **Total allocated `Cmpnts` structs:** 216.
*   **Total allocated scalar values:** `216 * 3 = 648`.
*   Accessible array indices: `ucont[k][j][i]` where `i,j,k` range `0..5`.

**2. Range Storing Meaningful "Interior" Face Flux Values:**

*   These are the fluxes across faces located strictly *between* physical cells (not boundary faces).
*   **U-flux (`.x`) across Interior i-Faces `Fi(i,j,k)`:**
    *   Face index `i` range: `1` to `IM-1 = 4`.
    *   Cell indices `j, k` range: `0` to `JM-1 = 4`, `0` to `KM-1 = 4`.
    *   Storage: `ucont[k][j][i].x` where `i=1..4`, `j=0..4`, `k=0..4`.
    *   Number: `(IM-1) * JM * KM = 4 * 5 * 5 = 100` locations.
*   **V-flux (`.y`) across Interior j-Faces `Fj(i,j,k)`:**
    *   Face index `j` range: `1` to `JM-1 = 4`.
    *   Cell indices `i, k` range: `0` to `IM-1 = 4`, `0` to `KM-1 = 4`.
    *   Storage: `ucont[k][j][i].y` where `i=0..4`, `j=1..4`, `k=0..4`.
    *   Number: `IM * (JM-1) * KM = 5 * 4 * 5 = 100` locations.
*   **W-flux (`.z`) across Interior k-Faces `Fk(i,j,k)`:**
    *   Face index `k` range: `1` to `KM-1 = 4`.
    *   Cell indices `i, j` range: `0` to `IM-1 = 4`, `0` to `JM-1 = 4`.
    *   Storage: `ucont[k][j][i].z` where `i=0..4`, `j=0..4`, `k=1..4`.
    *   Number: `IM * JM * (KM-1) = 5 * 5 * 4 = 100` locations.
*   **Total Meaningful Interior Flux Scalars:** `100 + 100 + 100 = 300`.

**3. Range Storing Meaningful "Boundary" Face Flux Values:**

*   These are the fluxes across the faces that coincide with the physical boundaries of the domain. Their values are set by `FormBCS`.
*   **U-flux (`.x`) across Boundary i-Faces `Fi(0,j,k)` and `Fi(IM,j,k)`:**
    *   Face indices `i=0` and `i=IM=5`.
    *   Cell indices `j, k` range: `0` to `JM-1 = 4`, `0` to `KM-1 = 4`.
    *   Storage: `ucont[k][j][0].x` and `ucont[k][j][5].x` where `j=0..4`, `k=0..4`.
    *   Number: `2 * JM * KM = 2 * 5 * 5 = 50` locations.
*   **V-flux (`.y`) across Boundary j-Faces `Fj(i,0,k)` and `Fj(i,JM,k)`:**
    *   Face indices `j=0` and `j=JM=5`.
    *   Cell indices `i, k` range: `0` to `IM-1 = 4`, `0` to `KM-1 = 4`.
    *   Storage: `ucont[k][0][i].y` and `ucont[k][5][i].y` where `i=0..4`, `k=0..4`.
    *   Number: `2 * IM * KM = 2 * 5 * 5 = 50` locations.
*   **W-flux (`.z`) across Boundary k-Faces `Fk(i,j,0)` and `Fk(i,j,KM)`:**
    *   Face indices `k=0` and `k=KM=5`.
    *   Cell indices `i, j` range: `0` to `IM-1 = 4`, `0` to `JM-1 = 4`.
    *   Storage: `ucont[0][j][i].z` and `ucont[5][j][i].z` where `i=0..4`, `j=0..4`.
    *   Number: `2 * IM * JM = 2 * 5 * 5 = 50` locations.
*   **Total Meaningful Boundary Flux Scalars:** `50 + 50 + 50 = 150`.

**4. Unused Storage Locations:**

*   **Total Allocated Scalars:** 648.
*   **Total Meaningful Flux Scalars (Interior + Boundary):** `300 + 150 = 450`.
*   **Number of Unused Scalar Values:** `648 - 450 = 198`.
*   **Explanation:** At each of the 216 node storage locations `ucont[k][j][i]`, only one of the three components (`.x`, `.y`, or `.z`) is used to store the physically meaningful flux across the corresponding face (`Fi(i)`, `Fj(j)`, or `Fk(k)`). The other two components at that specific `[k][j][i]` index are unused by this convention. (`216 locations * 2 unused components/location = 432` unused scalar slots? No, that's wrong).
*   **Correct Calculation:** Total slots = 216. Each holds 3 scalars (648 total). Meaningful fluxes = 450 scalars. Unused = 648 - 450 = 198 scalar slots. This means out of the 216 `Cmpnts` structs, some components are unused. For example, `ucont[k][j][i].y` and `ucont[k][j][i].z` have no meaning relative to the i-face `Fi(i)` whose flux is stored in `ucont[k][j][i].x`.

**Summary Table for `ucont`**

| Feature                         | Description                                      | Storage Index `[k][j][i]` Range & Component Used                | Number of Meaningful Scalar Locations |
| :------------------------------ | :----------------------------------------------- | :-------------------------------------------------------------- | :------------------------------------ |
| **Total Allocation Size**       | Full node-indexed grid (`fda`)                   | `i,j,k = 0..5` (for `.x`, `.y`, `.z`)                           | 648 (216 `Cmpnts`)                  |
| **Interior i-Face Flux Storage**| U-flux across interior face `Fi(i,j,k)`          | `i=1..4`, `j=0..4`, `k=0..4` --> **`.x`** component used         | 100                                   |
| **Interior j-Face Flux Storage**| V-flux across interior face `Fj(i,j,k)`          | `i=0..4`, `j=1..4`, `k=0..4` --> **`.y`** component used         | 100                                   |
| **Interior k-Face Flux Storage**| W-flux across interior face `Fk(i,j,k)`          | `i=0..4`, `j=0..4`, `k=1..4` --> **`.z`** component used         | 100                                   |
| **Boundary i-Face Flux Storage**| U-flux across boundary face `Fi(0/IM,j,k)`       | `i=0,5`, `j=0..4`, `k=0..4` --> **`.x`** component used        | 50                                    |
| **Boundary j-Face Flux Storage**| V-flux across boundary face `Fj(i,0/JM,k)`       | `i=0..4`, `j=0,5`, `k=0..4` --> **`.y`** component used        | 50                                    |
| **Boundary k-Face Flux Storage**| W-flux across boundary face `Fk(i,j,0/KM)`       | `i=0..4`, `j=0..4`, `k=0,5` --> **`.z`** component used        | 50                                    |
| **Total Meaningful Flux Values**| Sum of Interior and Boundary Flux Scalars        | N/A                                                             | 450                                   |
| **Unused Scalar Slots**         | Allocated scalar slots not holding meaningful flux | Components `.y`,`.z` at i-face indices; `.x`,`.z` at j-face; etc. | 198                                   |

