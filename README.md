Of course! Here is the provided information formatted as a clean and professional Markdown `README.md` file.

---

# 🌀 PICurv: A Parallel Particle-in-Cell Solver for Curvilinear LES

> **PICurv** is a parallel particle-in-cell (PIC) solver designed for low-Mach-number turbulent flows in complex geometries. Built using PETSc’s `DMDA` and `DMSwarm` infrastructure, it features:
> * A **fractional-step incompressible Navier–Stokes solver** with **dynamic Smagorinsky LES**
> * Support for **body-fitted curvilinear grids** and **immersed boundaries (CURVIB)**
> * Fully **MPI-parallel particle tracking** with two-way Eulerian–Lagrangian coupling
> * Trilinear interpolation (grid → particle) and conservative face-based projection (particle → grid)
>
> The solver has been validated on canonical geometries such as curved pipes and demonstrates scalable performance for **millions of particles** on both structured Cartesian and curvilinear meshes.

---

## 🔧 Key Features

*   **Parallelized 3D particle-in-cell solver** with PETSc
*   **Low-Mach number incompressible Navier–Stokes** with LES modeling
*   **Two-way coupled Lagrangian particle transport**
*   **Trilinear interpolation and deposition kernels**
*   **Geometric flexibility** via immersed boundary methods (CURVIB)
*   **Modular architecture** with plug-and-play routines for scalar transport, stochastic mixing, and particle diagnostics

---

## 🚀 Getting Started

### Dependencies

*   PETSc 3.20.3 or newer (built with MPI and `DMSWARM` support)

### Building the Project

To build the solver and postprocessor:

```bash
make inttest       # Builds main solver executable `inttest`
make postprocess   # Builds postprocessor executable `postprocess`
```

Ensure PETSc paths are configured correctly in your Makefile or environment.

### Directory Structure

```text
.
├── src/           # Source files (*.c)
├── include/       # Header files (*.h)
├── scripts/       # Python utilities
├── params/        # Input files for a run (grid, config, BCs)
├── run/           # Example run directory (created manually)
│   ├── control.dat
│   ├── params/
│   │   ├── config.dat
│   │   ├── bcs.dat
│   │   └── grid.dat (optional if using generated grid)
│   ├── inttest -> ../inttest (symlink)
│   └── postprocess -> ../postprocess (symlink)
```

---

## ▶️ Running the Code

1.  **Create a new `run/` directory**.
2.  **Place a `control.dat`** file inside. See below for a sample.
3.  **Create symbolic links** to the built executables:
    ```bash
    ln -s ../inttest inttest
    ln -s ../postprocess postprocess
    ```
4.  Ensure a `params/` subdirectory exists within `run/` with:
    *   `config.dat`
    *   `bcs.dat`
    *   (optional) `grid.dat` if using curvilinear meshes
5.  Run the solver:
    ```bash
    mpirun -np 4 ./inttest -f control.dat
    ```

---

## 📄 Sample `control.dat` Configuration

```plaintext
###########################################################################################################################
### I/O & Management  ###
-tio 1
-logfreq 10
-func_config_file "params/config.dat"
-grid_file "params/grid_bent_tube.dat"
-bcs_file "params/bcs.dat"
-Setup_Only 0

### PARTICLES
-numParticles 100
-pinit 0

### SIMULATION
-dt 0.1
-totalsteps 10
-rstart 0
-ti 0.0
-finit 2
-ucont_x 0.0
-ucont_y 0.0
-ucont_z 1.0

### GRID
-grid 0
-xMin 0.00
-xMax 1.00
-yMin 0.00
-yMax 1.00
-zMin 0.00
-zMax 1.00
-nblk 1
-im 20
-jm 20
-km 280
-r_x 1.0
-r_y 1.0
-r_z 1.0

### PARALLELIZATION
-dm_processors_x 2
-dm_processors_y 2
#-dm_processors_z 2
###########################################################################################################################
```

---

## 🗖️ Grid File Format (`grid.dat`)

To use a pre-generated body-fitted or curvilinear mesh, create a file `grid.dat` in the `params/` directory with the following structure:

```text
FDFGRID                # format tag
1                      # number of blocks
20 20 280              # number of cells in i, j, k for the block
x0 y0 z0               # node 0 coordinates
x1 y1 z1
x2 y2 z2
...
```

The coordinates specify node positions along the primary directions, and the grid is assumed to be structured and logically rectangular.

---

## 📂 Function-Level Logging (Sample `config.dat`)

Use `params/config.dat` to selectively enable logging of specific functions when running with `export LOG_LEVEL=DEBUG`.

Example:

```text
# ============================================================================
#               Function Logging Allow-List for Phase 1
# ============================================================================
main
AdvanceSimulation_TEST
PerformInitialSetup_TEST

# Boundary condition routines (uncomment to debug):
#Create_InletConstantVelocity
#Apply_InletConstantVelocity
#Initialize_InletConstantVelocity
```

This allows clean and controlled tracing of solver behavior for debugging or testing.

---

## 🧭 Grid Indexing and Layout (Arakawa C-grid)

PICurv follows a node-centered Arakawa C-grid indexing strategy with ghost cells and directional velocity storage. The figures below illustrate this storage scheme:

### Boundary Cell C(0,0) Layout

![Boundary Cell Layout](bd_cell_storage.png)

### Interior Cell C(1,1) Layout

![Interior Cell Layout](interior_cell_storage.png)

### Schematic: Physical and Ghost Cells

![Schematic Grid Layout](schematic%20examle.png)

--- ## Physical Locations of Key Vector and Scalar Quantities This table outlines the physical locations of various fields used in the CFD solver, relative to computational cell `(ic,jc,kc)` and its bounding nodes and faces. | Variable Name (Array Access Example)         | Data Type     | Physical Quantity                                   | Primary Physical Location                                     | Associated DMDA | Notes                                                                                                                               | | :------------------------------------------- | :------------ | :-------------------------------------------------- | :------------------------------------------------------------ | :-------------- | :---------------------------------------------------------------------------------------------------------------------------------- | | **Cell-Centered Quantities**                 |               |                                                     | **Geometric Center of Cell `(ic,jc,kc)`**                     |                 | Stored with array indices `[kc][jc][ic]`                                                                                              | | `P[kc][jc][ic]`                              | `PetscScalar` | Pressure                                            | Center of Cell `(ic,jc,kc)`                                   | `da`            |                                                                                                                                     | | `ucat[kc][jc][ic]`                           | `Cmpnts`      | Cartesian Velocity (u,v,w)                          | Center of Cell `(ic,jc,kc)`                                   | `fda`           | Derived from face-based `ucont`.                                                                                                      | | `Aj[kc][jc][ic]`                             | `PetscScalar` | Jacobian Inverse (1/J)                              | Center of Cell `(ic,jc,kc)`                                   | `da`            | Primary Jacobian calculation.                                                                                                       | | `nvert[kc][jc][ic]`                          | `PetscScalar` | Cell Status / Blanking Flag                         | Center of Cell `(ic,jc,kc)`                                   | `da`            | Used for IBM, overset grids.                                                                                                        | | `K_Omega[kc][jc][ic]` (if cell-centered)    | `Cmpnts`      | Turbulence k & ω                                    | Center of Cell `(ic,jc,kc)`                                   | `fda2`          | For RANS models.                                                                                                                    | | `Cent[kc][jc][ic]`                           | `Cmpnts`      | Physical Coords (x,y,z) of Cell Center              | Center of Cell `(ic,jc,kc)`                                   | `fda`           | Calculated by averaging node coordinates.                                                                                             | | `GridSpace[kc][jc][ic]`                      | `Cmpnts`      | Characteristic Grid Spacings (dx,dy,dz-like)        | Center of Cell `(ic,jc,kc)`                                   | `fda`           | Represents distances between face centers in computational directions.                                                                | | **Node-Based Quantities**                    |               |                                                     | **Geometric Grid Node**                                       |                 |                                                                                                                                     | | `coor[kn][jn][in]`                           | `Cmpnts`      | Physical Coords (x,y,z) of Node                     | Grid Node `(in,jn,kn)`                                        | `fda`           | `coor[kc][jc][ic]` is the Bottom-Left-Front node of cell `(ic,jc,kc)`.                                                              | | **Face-Based Quantities (Staggered)**        |               |                                                     | **On the Surface of a Cell Face**                             |                 |                                                                                                                                     | | `ucont[k_any][j_any][iface].x`               | `PetscScalar` | U-Contravariant Velocity (ξ-component)              | On I-Face `iface`                                             | `fda`           | Normal velocity component to I-Face. `csi` is collocated.                                                                           | | `csi[k_any][j_any][iface]`                   | `Cmpnts`      | Metric Vector ∇ξ (ξ_x, ξ_y, ξ_z)                    | On I-Face `iface`                                             | `fda` (`lCsi`)  | Related to area vector of I-Face.                                                                                                   | | `ucont[k_any][jface][i_any].y`               | `PetscScalar` | V-Contravariant Velocity (η-component)              | On J-Face `jface`                                             | `fda`           | Normal velocity component to J-Face. `eta` is collocated.                                                                           | | `eta[k_any][jface][i_any]`                   | `Cmpnts`      | Metric Vector ∇η (η_x, η_y, η_z)                    | On J-Face `jface`                                             | `fda` (`lEta`)  | Related to area vector of J-Face.                                                                                                   | | `ucont[kface][j_any][i_any].z`               | `PetscScalar` | W-Contravariant Velocity (ζ-component)              | On K-Face `kface`                                             | `fda`           | Normal velocity component to K-Face. `zet` is collocated.                                                                           | | `zet[kface][j_any][i_any]`                   | `Cmpnts`      | Metric Vector ∇ζ (ζ_x, ζ_y, ζ_z)                    | On K-Face `kface`                                             | `fda` (`lZet`)  | Related to area vector of K-Face.                                                                                                   | | **Face-Center Based Quantities (Interpolated/Calculated)** |               |                                       | **Geometric Center of a Cell Face**                         |                 |                                                                                                                                     | | `Centx[k_any][j_any][iface]`                 | `Cmpnts`      | Physical Coords (x,y,z) of I-Face Center            | Geometric Center of I-Face `iface`                            | `fda`           |                                                                                                                                     | | `IAj[k_any][j_any][iface]`                   | `PetscScalar` | Jacobian Inverse (1/J) at I-Face Center             | Geometric Center of I-Face `iface`                            | `da`            |                                                                                                                                     | | `ICsi[k_any][j_any][iface]` (and IEta, IZet) | `Cmpnts`      | Metrics (∇ξ, ∇η, ∇ζ) at I-Face Center               | Geometric Center of I-Face `iface`                            | `fda`           |                                                                                                                                     | | `Centy[k_any][jface][i_any]`                 | `Cmpnts`      | Physical Coords (x,y,z) of J-Face Center            | Geometric Center of J-Face `jface`                            | `fda`           |                                                                                                                                     | | `JAj[k_any][jface][i_any]`                   | `PetscScalar` | Jacobian Inverse (1/J) at J-Face Center             | Geometric Center of J-Face `jface`                            | `da`            |                                                                                                                                     | | `JCsi[k_any][jface][i_any]` (and JEta, JZet) | `Cmpnts`      | Metrics (∇ξ, ∇η, ∇ζ) at J-Face Center               | Geometric Center of J-Face `jface`                            | `fda`           |                                                                                                                                     | | `Centz[kface][j_any][i_any]`                 | `Cmpnts`      | Physical Coords (x,y,z) of K-Face Center            | Geometric Center of K-Face `kface`                            | `fda`           |                                                                                                                                     | | `KAj[kface][j_any][i_any]`                   | `PetscScalar` | Jacobian Inverse (1/J) at K-Face Center             | Geometric Center of K-Face `kface`                            | `da`            |                                                                                                                                     | | `KCsi[kface][j_any][i_any]` (and KEta, KZet) | `Cmpnts`      | Metrics (∇ξ, ∇η, ∇ζ) at K-Face Center               | Geometric Center of K-Face `kface`                            | `fda`           |                                                                                                                                     | **Notes on `ucont`, `csi`, `eta`, `zet` Indexing Relative to Cell `(ic,jc,kc)`:** *   **Left I-Face (Face `ic`):** *   U-velocity: `ucont[kc][jc][ic].x` *   Metrics: `csi[kc][jc][ic]` *   **Right I-Face (Face `ic+1`):** *   U-velocity: `ucont[kc][jc][ic+1].x` *   Metrics: `csi[kc][jc][ic+1]` *   **Bottom J-Face (Face `jc`):** *   V-velocity: `ucont[kc][jc][ic].y` (Assuming the `jc` index corresponds to the J-Face `jc` when accessing the `.y` component for cell `ic,kc`) *   Metrics: `eta[kc][jc][ic]` *   **Top J-Face (Face `jc+1`):** *   V-velocity: `ucont[kc][jc+1][ic].y` *   Metrics: `eta[kc][jc+1][ic]` *   **Front K-Face (Face `kc`):** *   W-velocity: `ucont[kc][jc][ic].z` (Assuming the `kc` index corresponds to the K-Face `kc` when accessing the `.z` component for cell `ic,jc`) *   Metrics: `zet[kc][jc][ic]` *   **Back K-Face (Face `kc+1`):** *   W-velocity: `ucont[kc+1][jc][ic].z` *   Metrics: `zet[kc+1][jc][ic]` The "k_any", "j_any", "i_any" in the table for face-based quantities indicate that the other two indices typically align with a cell or a path along the face, while the defining index (`iface`, `jface`, `kface`) specifies the particular face sheet. The exact mapping for the non-defining indices depends on how the loops iterate over the faces but usually aligns with the `k,j,i` of one of the cells bounded by that face. **Grid Structure** **Assumptions:** *   `IM`, `JM`, `KM`: Number of *allocated cells*. *   `ucat[k][j][i]`: Cell-centered value for cell with origin node `(i,j,k)`. *   `coor[k][j][i]`: Coordinate `(x,y,z)` of **node** `(i,j,k)`. Node indices: `i=0..IM`, `j=0..JM`, `k=0..KM`. *   `ucont[k][j][i].x`: x-flux component associated with the **constant-i face** passing through node `(i,j,k)`. This face separates cell `(i-1, j, k)` from cell `(i, j, k)`. Face index `i` runs `0..IM`. *   `ucont[k][j][i].y`: y-flux component associated with the **constant-j face** passing through node `(i,j,k)`. This face separates cell `(i, j-1, k)` from cell `(i, j, k)`. Face index `j` runs `0..JM`. *   `ucont[k][j][i].z`: z-flux component associated with the **constant-k face** passing through node `(i,j,k)`. This face separates cell `(i, j, k-1)` from cell `(i, j, k)`. Face index `k` runs `0..KM`. *   Coordinate origin `(0,0,0)` is typically at node `(0,0,0)`. **Explanation of Face Indices:** *   **`Ucont.x` (i-faces):** *   `i=0`: Boundary face before the first ghost cell `(0,j,k)`. *   `i=1`: Face between ghost cell `(0,j,k)` and first physical cell `(1,j,k)`. Needed by cell `(1,j,k)`. *   `i=2..IM-2`: Internal faces *within* the physical domain. Needed by adjacent physical cells. *   `i=IM-1`: Face between last physical cell `(IM-2,j,k)` and ghost cell `(IM-1,j,k)`. Needed by cell `(IM-2,j,k)`. *   `i=IM`: Boundary face after the last ghost cell `(IM-1,j,k)`. *   **`Ucont.y` (j-faces):** Analogous structure with index `j` running from `0` to `JM`. Physical domain requires faces `j=1..JM-1`. *   **`Ucont.z` (k-faces):** Analogous structure with index `k` running from `0` to `KM`. Physical domain requires faces `k=1..KM-1`. *   `IM`, `JM`, `KM`: Number of *allocated cells* in the i, j, k directions, respectively. *   `ucat[k][j][i]`: Stores the cell-centered value for the cell whose origin node (lower-left-back corner) is `(i,j,k)`. *   `coor[k][j][i]`: Stores the physical coordinate `(x,y,z)` of the grid **node** `(i,j,k)`. *   The physical coordinate system origin `(0,0,0)` is typically located at node `(0,0,0)`, i.e., `coor[0][0][0] == (0,0,0)`. **Grid Layout Table:** | Feature           | Description                          | Cell Indices (for `ucat`)                     | Node Indices (for `coor`, defining cells) | Face Indices: `Ucont.x` (i-faces)                 | Face Indices: `Ucont.y` (j-faces)                 | Face Indices: `Ucont.z` (k-faces)                 | Physical Coordinate Range (Approx.)               | | :---------------- | :----------------------------------- | :-------------------------------------------- | :---------------------------------------- | :------------------------------------------------ | :------------------------------------------------ | :------------------------------------------------ | :------------------------------------------------ | | **Allocated Grid**  | Entire grid in memory              | `i=0..IM-1` `j=0..JM-1` `k=0..KM-1`           | `i=0..IM` `j=0..JM` `k=0..KM`             | `i=0..IM` (`j=0..JM-1`, `k=0..KM-1`)              | `j=0..JM` (`i=0..IM-1`, `k=0..KM-1`)              | `k=0..KM` (`i=0..IM-1`, `j=0..JM-1`)              | `coor[0][0][0]` to `coor[KM][JM][IM]`             | | **Physical Domain** | Region for solving physics         | `i=1..IM-2` `j=1..JM-2` `k=1..KM-2`           | Nodes `(1,1,1)` to `(IM-1,JM-1,KM-1)`     | `i=1..IM-1` (`j=1..JM-2`, `k=1..KM-2`)            | `j=1..JM-1` (`i=1..IM-2`, `k=1..KM-2`)            | `k=1..KM-1` (`i=1..IM-2`, `j=1..JM-2`)            | `coor[1][1][1]` to `coor[KM-1][JM-1][IM-1]`       | | **Ghost Layer**   | Cells/Faces for boundary data      | `i=0` or `IM-1` OR `j=0` or `JM-1` OR `k=0` or `KM-1` | Nodes `i=0,IM` OR `j=0,JM` OR `k=0,KM`      | Faces `i=0, IM` (Set by BCs/MPI)                | Faces `j=0, JM` (Set by BCs/MPI)                | Faces `k=0, KM` (Set by BCs/MPI)                | Space outside physical node range             | **Explanation Notes:** *   **Cell Indices:** Refer to the cell whose origin node (minimum i,j,k corner) has these indices. Used for `ucat`, `nvert`. *   **Node Indices:** Refer to the actual grid point locations. Used for `coor`. These nodes define the corners and extent of the cells. *   **Face Indices:** *   `Ucont.x` lives on constant-`i` faces. Index `i` indicates the face passing through nodes with that `i` coordinate. The range in parentheses shows the relevant `j,k` node indices spanned by that face within the specified domain (allocated/physical). *   `Ucont.y` lives on constant-`j` faces. Index `j` indicates the face passing through nodes with that `j` coordinate. *   `Ucont.z` lives on constant-`k` faces. Index `k` indicates the face passing through nodes with that `k` coordinate. *   **Physical Coordinates:** Assumes `coor[0][0][0]` is the physical origin `(0,0,0)`. The physical domain then starts roughly one grid step away (`coor[1][1][1]`). **Key Takeaways:** 1.  **Index Shift:** The **physical** cell indices are shifted inwards by one compared to the allocated indices (`1..IM-2` vs `0..IM-1`). 2.  **Consistent Storage:** A cell-centered value for cell `(i,j,k)` (whether physical or ghost) is always stored at `ucat[k][j][i]`. A node coordinate for node `(i,j,k)` is always stored at `coor[k][j][i]`. 3.  **Coordinate Mapping:** The node `(0,0,0)` (corner of the first ghost cell `(0,0,0)`) usually corresponds to the physical origin `(0,0,0)`. The node `(1,1,1)` (corner of the first physical cell `(1,1,1)`) is then located one grid step away from the origin. 4.  **Boundary Conditions:** Functions like `FormBCS` operate on the **ghost layer** (e.g., setting `ucont` at indices `i=0`, `i=IM-1`, etc.) using information extrapolated or copied from the adjacent **physical cells** (e.g., using `ucat` at indices `i=1`, `i=IM-2`). 5.  **DMDA Role:** PETSc's `DMGlobalToLocalBegin/End` automatically handles filling the ghost cell values in `ucat`, `ucont`, etc., with data from neighboring processes when running in parallel. Boundary condition functions handle filling the ghost values at the true domain boundaries. Okay, let's add the face-centered variables (`Ucont`, and implicitly metrics like `csi`, `eta`, `zet`) to the table, keeping the Arakawa C-grid structure in mind. **Key Points:** *   The `FormBCS` function explicitly sets `ucont` values on the **Ghost Layer faces** (`i=0`, `i=IM`, `j=0`, `j=JM`, `k=0`, `k=KM` - although the code often uses indices like `i=mx-2` to access the *face* data structure index corresponding to the *last* physical cell boundary). *   Computations within the physical domain cells (`i=1..IM-2`, etc.) will typically read `ucont` values from the faces bounding them (i.e., face indices `i=1..IM-1`, `j=1..JM-1`, `k=1..KM-1`).

---

## 📜 Citation / Acknowledgements

This work is supported by the **National Science Foundation (NSF Award No. 2309630)** and computational resources from the **High-Performance Research Computing (HPRC) center at Texas A\&M University**.

---

Let me know if you’d like to add performance benchmarks, case studies, or a “Contributing” section next.
