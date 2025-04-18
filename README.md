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

**Grid Structure** 
**Assumptions:**

*   `IM`, `JM`, `KM`: Number of *allocated cells*.
*   `ucat[k][j][i]`: Cell-centered value for cell with origin node `(i,j,k)`.
*   `coor[k][j][i]`: Coordinate `(x,y,z)` of **node** `(i,j,k)`. Node indices: `i=0..IM`, `j=0..JM`, `k=0..KM`.
*   `ucont[k][j][i].x`: x-flux component associated with the **constant-i face** passing through node `(i,j,k)`. This face separates cell `(i-1, j, k)` from cell `(i, j, k)`. Face index `i` runs `0..IM`.
*   `ucont[k][j][i].y`: y-flux component associated with the **constant-j face** passing through node `(i,j,k)`. This face separates cell `(i, j-1, k)` from cell `(i, j, k)`. Face index `j` runs `0..JM`.
*   `ucont[k][j][i].z`: z-flux component associated with the **constant-k face** passing through node `(i,j,k)`. This face separates cell `(i, j, k-1)` from cell `(i, j, k)`. Face index `k` runs `0..KM`.
*   Coordinate origin `(0,0,0)` is typically at node `(0,0,0)`.

**Explanation of Face Indices:**

*   **`Ucont.x` (i-faces):**
    *   `i=0`: Boundary face before the first ghost cell `(0,j,k)`.
    *   `i=1`: Face between ghost cell `(0,j,k)` and first physical cell `(1,j,k)`. Needed by cell `(1,j,k)`.
    *   `i=2..IM-2`: Internal faces *within* the physical domain. Needed by adjacent physical cells.
    *   `i=IM-1`: Face between last physical cell `(IM-2,j,k)` and ghost cell `(IM-1,j,k)`. Needed by cell `(IM-2,j,k)`.
    *   `i=IM`: Boundary face after the last ghost cell `(IM-1,j,k)`.
*   **`Ucont.y` (j-faces):** Analogous structure with index `j` running from `0` to `JM`. Physical domain requires faces `j=1..JM-1`.
*   **`Ucont.z` (k-faces):** Analogous structure with index `k` running from `0` to `KM`. Physical domain requires faces `k=1..KM-1`.
*   `IM`, `JM`, `KM`: Number of *allocated cells* in the i, j, k directions, respectively.
*   `ucat[k][j][i]`: Stores the cell-centered value for the cell whose origin node (lower-left-back corner) is `(i,j,k)`.
*   `coor[k][j][i]`: Stores the physical coordinate `(x,y,z)` of the grid **node** `(i,j,k)`.
*   The physical coordinate system origin `(0,0,0)` is typically located at node `(0,0,0)`, i.e., `coor[0][0][0] == (0,0,0)`.

**Grid Layout Table:**

| Feature           | Description                          | Cell Indices (for `ucat`)                     | Node Indices (for `coor`, defining cells) | Face Indices: `Ucont.x` (i-faces)                 | Face Indices: `Ucont.y` (j-faces)                 | Face Indices: `Ucont.z` (k-faces)                 | Physical Coordinate Range (Approx.)               |
| :---------------- | :----------------------------------- | :-------------------------------------------- | :---------------------------------------- | :------------------------------------------------ | :------------------------------------------------ | :------------------------------------------------ | :------------------------------------------------ |
| **Allocated Grid**  | Entire grid in memory              | `i=0..IM-1` `j=0..JM-1` `k=0..KM-1`           | `i=0..IM` `j=0..JM` `k=0..KM`             | `i=0..IM` (`j=0..JM-1`, `k=0..KM-1`)              | `j=0..JM` (`i=0..IM-1`, `k=0..KM-1`)              | `k=0..KM` (`i=0..IM-1`, `j=0..JM-1`)              | `coor[0][0][0]` to `coor[KM][JM][IM]`             |
| **Physical Domain** | Region for solving physics         | `i=1..IM-2` `j=1..JM-2` `k=1..KM-2`           | Nodes `(1,1,1)` to `(IM-1,JM-1,KM-1)`     | `i=1..IM-1` (`j=1..JM-2`, `k=1..KM-2`)            | `j=1..JM-1` (`i=1..IM-2`, `k=1..KM-2`)            | `k=1..KM-1` (`i=1..IM-2`, `j=1..JM-2`)            | `coor[1][1][1]` to `coor[KM-1][JM-1][IM-1]`       |
| **Ghost Layer**   | Cells/Faces for boundary data      | `i=0` or `IM-1` OR `j=0` or `JM-1` OR `k=0` or `KM-1` | Nodes `i=0,IM` OR `j=0,JM` OR `k=0,KM`      | Faces `i=0, IM` (Set by BCs/MPI)                | Faces `j=0, JM` (Set by BCs/MPI)                | Faces `k=0, KM` (Set by BCs/MPI)                | Space outside physical node range             |

**Explanation Notes:**

*   **Cell Indices:** Refer to the cell whose origin node (minimum i,j,k corner) has these indices. Used for `ucat`, `nvert`.
*   **Node Indices:** Refer to the actual grid point locations. Used for `coor`. These nodes define the corners and extent of the cells.
*   **Face Indices:**
    *   `Ucont.x` lives on constant-`i` faces. Index `i` indicates the face passing through nodes with that `i` coordinate. The range in parentheses shows the relevant `j,k` node indices spanned by that face within the specified domain (allocated/physical).
    *   `Ucont.y` lives on constant-`j` faces. Index `j` indicates the face passing through nodes with that `j` coordinate.
    *   `Ucont.z` lives on constant-`k` faces. Index `k` indicates the face passing through nodes with that `k` coordinate.
*   **Physical Coordinates:** Assumes `coor[0][0][0]` is the physical origin `(0,0,0)`. The physical domain then starts roughly one grid step away (`coor[1][1][1]`).


**Key Takeaways:**

1.  **Index Shift:** The **physical** cell indices are shifted inwards by one compared to the allocated indices (`1..IM-2` vs `0..IM-1`).
2.  **Consistent Storage:** A cell-centered value for cell `(i,j,k)` (whether physical or ghost) is always stored at `ucat[k][j][i]`. A node coordinate for node `(i,j,k)` is always stored at `coor[k][j][i]`.
3.  **Coordinate Mapping:** The node `(0,0,0)` (corner of the first ghost cell `(0,0,0)`) usually corresponds to the physical origin `(0,0,0)`. The node `(1,1,1)` (corner of the first physical cell `(1,1,1)`) is then located one grid step away from the origin.
4.  **Boundary Conditions:** Functions like `FormBCS` operate on the **ghost layer** (e.g., setting `ucont` at indices `i=0`, `i=IM-1`, etc.) using information extrapolated or copied from the adjacent **physical cells** (e.g., using `ucat` at indices `i=1`, `i=IM-2`).
5.  **DMDA Role:** PETSc's `DMGlobalToLocalBegin/End` automatically handles filling the ghost cell values in `ucat`, `ucont`, etc., with data from neighboring processes when running in parallel. Boundary condition functions handle filling the ghost values at the true domain boundaries.
Okay, let's add the face-centered variables (`Ucont`, and implicitly metrics like `csi`, `eta`, `zet`) to the table, keeping the Arakawa C-grid structure in mind.

**Key Points:**

*   The `FormBCS` function explicitly sets `ucont` values on the **Ghost Layer faces** (`i=0`, `i=IM`, `j=0`, `j=JM`, `k=0`, `k=KM` - although the code often uses indices like `i=mx-2` to access the *face* data structure index corresponding to the *last* physical cell boundary).
*   Computations within the physical domain cells (`i=1..IM-2`, etc.) will typically read `ucont` values from the faces bounding them (i.e., face indices `i=1..IM-1`, `j=1..JM-1`, `k=1..KM-1`).




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
