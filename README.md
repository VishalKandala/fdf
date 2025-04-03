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

Okay, let's break down the `ucat` array for `IM=JM=KM=5`.

**1. Total Size of `ucat`:**

*   As established, `ucat` uses the node-based `fda` DMDA.
*   Number of nodes: `M=IM+1=6`, `N=JM+1=6`, `P=KM+1=6`.
*   Total node locations: `M * N * P = 6 * 6 * 6 = 216`.
*   Degrees of Freedom (DOF): 3 (for `ucat.x`, `ucat.y`, `ucat.z`).
*   **Total allocated `Cmpnts` structs:** 216.
*   **Total allocated scalar values:** `216 * 3 = 648`.
*   The accessible array indices are `ucat[k][j][i]` where `i` ranges `0..5`, `j` ranges `0..5`, `k` ranges `0..5`.

**2. Range Storing Corresponding "Interior" Cell Center `ucat` Values:**

*   **Conceptual Cells:** The interior physical cells are indexed `i=1..IM-1`, `j=1..JM-1`, `k=1..KM-1`. For `IM=5`, this is `i=1..4`, `j=1..4`, `k=1..4`.
*   **Storage Convention:** The code stores the `ucat` value calculated by `Contra2Cart` (which represents the state at the center of cell `C(i,j,k)`) in the storage slot corresponding to **node index `(i,j,k)`**.
*   **Index Range:**
    *   `i` from `1` to `4` (`IM-1`)
    *   `j` from `1` to `4` (`JM-1`)
    *   `k` from `1` to `4` (`KM-1`)
*   **Number of Locations:** `(IM-1 - 1 + 1) * (JM-1 - 1 + 1) * (KM-1 - 1 + 1) = 4 * 4 * 4 = 64` locations.
*   **Meaning:** `ucat[k][j][i]` within this range holds the Cartesian velocity computed by `Contra2Cart`, representing the physical state at the center of cell `C(i,j,k)`.

**3. Range Storing Boundary Extrapolated "Ghost" `ucat` Values:**

*   **Storage Indices:** These values are stored at the node indices that lie on the boundaries of the `0..IM`, `0..JM`, `0..KM` node index space.
*   **Index Range:** `ucat[k][j][i]` where:
    *   `i = 0` OR `i = 5` (`IM`)
    *   OR `j = 0` OR `j = 5` (`JM`)
    *   OR `k = 0` OR `k = 5` (`KM`)
*   **Number of Locations:** This is the total number of node locations minus the number of interior node locations used for cell values: `216 - 64 = 152` locations.
*   **Meaning:** The values `ucat[k][j][i]` within this range **do not** represent the physical state at the center of the adjacent boundary cell (e.g., `ucat[0][0][0]` doesn't hold the value for the center of `C(0,0,0)`). Instead, they store the **extrapolated values calculated by `FormBCS`** (e.g., `ucat[0][0][0] = 2*ubcs[0][0][0] - ucat[1][0][0]`). These act as the necessary "ghost" values required by computational stencils operating in the first/last layer of interior cells (e.g., cells `C(1,j,k)` need the value stored at `ucat[k][j][0]`).

**4. Unused Storage Locations?**

*   Based on this convention, **all 216 node locations** `ucat[k][j][i]` (where `i,j,k` range `0..5`) are used to store meaningful data for the numerical scheme:
    *   Indices `1..4` store interior cell values.
    *   Indices `0` and `5` store extrapolated boundary/ghost values.
*   Therefore, there are **0 unused `Cmpnts` struct storage locations**. All allocated slots hold either a value representing an interior cell or an extrapolated value representing the boundary condition's effect for the interior stencils. (It's possible some *components* within the `Cmpnts` struct at boundary locations might be less critical depending on the BC, but the extrapolation in `FormBCS` calculates all three).

**Summary Table:**

| Feature                         | Description                                                                                                | Index Range (`i`, `j`, `k`) | Number of Locations (`Cmpnts` structs) |
| :------------------------------ | :--------------------------------------------------------------------------------------------------------- | :-------------------------- | :------------------------------------- |
| **Total Allocation Size**       | Full node-indexed grid (`fda`)                                                                             | `0..5`, `0..5`, `0..5`      | 216                                    |
| **Interior Cell Value Storage** | Slot `[k][j][i]` stores value representing center of cell `C(i,j,k)`                                       | `1..4`, `1..4`, `1..4`      | 64                                     |
| **Boundary/Ghost Value Storage**| Slot `[k][j][i]` stores extrapolated value from `FormBCS` based on physical BC at the corresponding boundary | `i/j/k` is `0` or `5`       | 152                                    |
| **Unused Locations**            | Locations not used for either purpose                                                                      | None                        | 0                                      |

`ucont`, the face-centered contravariant velocity/flux, for the `IM=JM=KM=5` case.

Storage versus the conceptual meaning for `ucont` (the face-centered contravariant velocity/flux).

**1. Storage Allocation:**

*   **Code Evidence:** Like `ucat`, `ucont` is created in `MG_Initial` by duplicating a vector (`user->Csi`) that is based on the `fda` DMDA (`DMCreateGlobalVector(user[bi].fda, &(user[bi].Csi))`).
*   **Structure:** This means `ucont` is allocated using the **node-based indexing** scheme (`0..IM`, `0..JM`, `0..KM`) corresponding to the `M x N x P` grid layout. Each storage location `(i,j,k)` holds `DOF=3` scalar values (accessed as `ucont[k][j][i].x`, `.y`, `.z`).
*   **Total Slots:** `(IM+1) * (JM+1) * (KM+1)` locations, each holding 3 scalars.

**2. Conceptual Meaning:**

*   **Physical Quantity:** `ucont` represents the contravariant velocity components (U, V, W), which are equivalent to the **volume flux normal to the cell faces** in computational space (when multiplied by the Jacobian, which is often implicitly included in the way metrics are defined or used).
*   **Location:** These fluxes physically occur **at the center of the cell faces**.
    *   The U-component (related to `ucont.x`) exists on faces of constant computational coordinate `i` (the i-faces, `Fi`).
    *   The V-component (related to `ucont.y`) exists on faces of constant computational coordinate `j` (the j-faces, `Fj`).
    *   The W-component (related to `ucont.z`) exists on faces of constant computational coordinate `k` (the k-faces, `Fk`).

**3. The Mapping Convention (Storage Index vs. Physical Face):**

The code uses a specific convention to store these face-centered flux components within the node-indexed `ucont` array:

*   **`ucont[k][j][i].x` stores the U-flux** (normal to the i-face) associated with the **face `Fi(i, j, k)`**. This is the face located *at* computational coordinate `i`, separating cell `C(i-1,j,k)` from `C(i,j,k)`. This mapping holds for face indices `i` from `0` to `IM`.
*   **`ucont[k][j][i].y` stores the V-flux** (normal to the j-face) associated with the **face `Fj(i, j, k)`**. This is the face located *at* computational coordinate `j`, separating cell `C(i,j-1,k)` from `C(i,j,k)`. This mapping holds for face indices `j` from `0` to `JM`.
*   **`ucont[k][j][i].z` stores the W-flux** (normal to the k-face) associated with the **face `Fk(i, j, k)`**. This is the face located *at* computational coordinate `k`, separating cell `C(i,j,k-1)` from `C(i,j,k)`. This mapping holds for face indices `k` from `0` to `KM`.

**Code Evidence for the Mapping:**

*   **`Contra2Cart` Averaging:**
    *   `q[0] = 0.5 * (ucont[k][j][i-1].x + ucont[k][j][i].x)` averages the U-flux from face `Fi(i-1)` (stored at index `i-1`) and face `Fi(i)` (stored at index `i`) to approximate the U-flux at the center of cell `C(i)`.
    *   Similar logic applies to `q[1]` using `.y` components at indices `j-1` and `j`, and `q[2]` using `.z` components at indices `k-1` and `k`.
*   **`Convection` Flux Calculation:**
    *   `ucon = ucont[k][j][i].x * 0.5` uses the value stored at index `i` (representing face `Fi(i)`) when calculating the advective transport across that specific i-face.
    *   The loops for `fp2` (j-flux) use `ucont[k][j][i].y` with `j` as the loop index, and `fp3` (k-flux) uses `ucont[k][j][i].z` with `k` as the loop index, confirming the component-direction association.
*   **`FormBCS` Boundary Setting:**
    *   Sets `ucont[k][j][0].x` for BC at face `Fi(0)`.
    *   Sets `ucont[k][0][i].y` for BC at face `Fj(0)`.
    *   Sets `ucont[0][j][i].z` for BC at face `Fk(0)`.
    *   (And similarly for upper boundary faces `Fi(IM)`, `Fj(JM)`, `Fk(KM)`).

**4. Boundary vs. Interior for `ucont`:**

*   This mapping convention applies identically to both **boundary faces** and **interior faces**.
*   `ucont[k][j][0].x` is the flux across the physical boundary face `Fi(0)`. Its *value* is determined by `FormBCS`.
*   `ucont[k][j][1].x` is the flux across the first interior face `Fi(1)`. Its *value* is determined by the flow solver (momentum update/projection).
*   `ucont[k][j][IM].x` is the flux across the other physical boundary face `Fi(IM)`. Its *value* is determined by `FormBCS`.

**5. Unused Components:**

*   Crucially, at a specific storage location `ucont[k][j][i]`, only *one* of the `.x`, `.y`, or `.z` components has a direct physical meaning according to this convention.
    *   At `ucont[k][j][i]`, only `.x` represents the flux (across face `Fi(i)`). The values stored in `.y` and `.z` at this location are irrelevant to face `Fi(i)`.
    *   Similarly, `ucont[k][j][i].y` is meaningful for face `Fj(j)`, while `.x` and `.z` stored at `[k][j][i]` are not directly related to that face.
*   These "unused" components at each storage location might contain zero, leftover values, or potentially temporary data depending on the solver implementation, but they don't represent the primary contravariant fluxes associated with that index according to the convention used in `Contra2Cart`, `Convection`, and `FormBCS`.

**Markdown Table Summary for `ucont`**

| Physical Face                   | Face Index `(i,j,k)` Range         | Type      | `ucont` Storage Index `[k][j][i]` Used | Component Used | Physical Meaning of Stored Value                                     |
| :------------------------------ | :--------------------------------- | :-------- | :------------------------------------ | :------------- | :------------------------------------------------------------------- |
| **Interior i-Face `Fi(i,j,k)`** | `1 <= i <= IM-1` <br> `0 <= j <= JM-1` <br> `0 <= k <= KM-1` | Interior  | `[k][j][i]`                         | `.x`           | U-flux across face `Fi(i,j,k)` (Calculated by solver)              |
| **Boundary i-Face `Fi(0,j,k)`** | `i=0` <br> `0 <= j <= JM-1` <br> `0 <= k <= KM-1` | Boundary  | `[k][j][0]`                         | `.x`           | U-flux across face `Fi(0,j,k)` (Set by `FormBCS`)                   |
| **Boundary i-Face `Fi(IM,j,k)`**| `i=IM` <br> `0 <= j <= JM-1` <br> `0 <= k <= KM-1` | Boundary  | `[k][j][IM]`                        | `.x`           | U-flux across face `Fi(IM,j,k)` (Set by `FormBCS`)                  |
| **Interior j-Face `Fj(i,j,k)`** | `0 <= i <= IM-1` <br> `1 <= j <= JM-1` <br> `0 <= k <= KM-1` | Interior  | `[k][j][i]`                         | `.y`           | V-flux across face `Fj(i,j,k)` (Calculated by solver)              |
| **Boundary j-Face `Fj(i,0,k)`** | `0 <= i <= IM-1` <br> `j=0` <br> `0 <= k <= KM-1` | Boundary  | `[k][0][i]`                         | `.y`           | V-flux across face `Fj(i,0,k)` (Set by `FormBCS`)                   |
| **Boundary j-Face `Fj(i,JM,k)`**| `0 <= i <= IM-1` <br> `j=JM` <br> `0 <= k <= KM-1` | Boundary  | `[k][JM][i]`                        | `.y`           | V-flux across face `Fj(i,JM,k)` (Set by `FormBCS`)                  |
| **Interior k-Face `Fk(i,j,k)`** | `0 <= i <= IM-1` <br> `0 <= j <= JM-1` <br> `1 <= k <= KM-1` | Interior  | `[k][j][i]`                         | `.z`           | W-flux across face `Fk(i,j,k)` (Calculated by solver)              |
| **Boundary k-Face `Fk(i,j,0)`** | `0 <= i <= IM-1` <br> `0 <= j <= JM-1` <br> `k=0` | Boundary  | `[0][j][i]`                         | `.z`           | W-flux across face `Fk(i,j,0)` (Set by `FormBCS`)                   |
| **Boundary k-Face `Fk(i,j,KM)`**| `0 <= i <= IM-1` <br> `0 <= j <= JM-1` <br> `k=KM` | Boundary  | `[KM][j][i]`                        | `.z`           | W-flux across face `Fk(i,j,KM)` (Set by `FormBCS`)                  |
| *Unused Components*             | *Any index*                        | N/A       | `[k][j][i]`                         | `.y`, `.z`     | *No direct meaning for flux across face `Fi(i,j,k)`*               |
|                                 |                                    |           | `[k][j][i]`                         | `.x`, `.z`     | *No direct meaning for flux across face `Fj(i,j,k)`*               |
|                                 |                                    |           | `[k][j][i]`                         | `.x`, `.y`     | *No direct meaning for flux across face `Fk(i,j,k)`*               |

This table clarifies that while `ucont[k][j][i]` stores 3 values, only one component (`.x`, `.y`, or `.z`) at that index `(i,j,k)` corresponds to the physical flux across the respective face (`Fi(i)`, `Fj(j)`, or `Fk(k)`).

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

![image](https://github.com/user-attachments/assets/46890cac-7c6e-4326-8fdb-0aaf7a5dcfc1)
![image](https://github.com/user-attachments/assets/f202afcf-d434-4f5b-8565-b215042199ad)
![image](https://github.com/user-attachments/assets/1dfad9bb-01ff-48d9-9c0e-7c5f8d348dc8)
![image](https://github.com/user-attachments/assets/ff4ea445-7831-468f-9e82-f636190340c1)
![image](https://github.com/user-attachments/assets/a2255c03-6a7c-4a5a-9db6-ae098f59d318)
![image](https://github.com/user-attachments/assets/72ed5af8-e0d3-4d9a-a019-0952f6b44a2c)

