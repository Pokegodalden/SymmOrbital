# SymmOrbital — Simplified Usage Guide

Welcome to **SymmOrbital**! This quick guide walks you through everything you can do in the application in plain, easy-to-follow steps—no prior computational chemistry setup required.

---

## 1. Pick or Search for a Molecule

At the top of the page, you have two fast ways to load a molecule:

1. **Search the Live PubChem Bar**:
   - Type any chemical name (e.g., `acetone`, `caffeine`, `sulfur tetrafluoride`), chemical formula (`H2O`, `XeF6`, `C6H6`), SMILES string, or PubChem CID (`962`) into the large search bar and press **Search PubChem** (or hit `Enter`).
   - SymmOrbital automatically fetches 3D coordinates from NCBI PubChem—or reconstructs accurate 3D VSEPR geometry if only a flat 2D record exists.
2. **Browse Benchmark Molecules**:
   - Click the **Benchmark Molecules (16 NIST PES Standards)** button right below the search bar to expand curated molecules grouped by category:
     - **Diatomics & Radicals** ($\mathrm{CO}, \mathrm{NO}$)
     - **Main-Group & VSEPR** ($\mathrm{H_2O}, \mathrm{NH_3}, \mathrm{CH_4}, \mathrm{CO_2}, \mathrm{BF_3}, \mathrm{SF_6}$)
     - **Hypervalent & Lone-Pair VSEPR** ($\mathrm{XeF_6}, \mathrm{XeF_4}, \mathrm{SF_4}, \mathrm{ClF_3}, \mathrm{PF_5}$)
     - **Organic $\pi$-Systems** ($\mathrm{C_2H_4}, \mathrm{CH_2O}, \mathrm{C_6H_6}$)

---

## 2. Sketch or Modify Molecules in 2D (Left Pane)

Use the **2D Chemical Sketcher** on the left to build custom molecules or tweak an existing structure:

- **Open Drawing Tools & Elements**:
  - Click **Tools & Elements** in the top-right corner of the 2D Sketcher pane to reveal bond tools, ring stamps, and element buttons.
- **Add Bonds & Stereochemistry**:
  - Choose **— Single**, **═ Double**, **≡ Triple**, **▲ Wedge** (pointing toward you, $+z$), or **┄ Dash** (pointing away from you, $-z$).
  - Click an empty spot on the grid to place a new bond, click an existing atom to sprout a new bonded atom, or click an existing bond to change its bond order or convert it into a **Wedge** or **Dash** bond.
- **Change Chemical Elements**:
  - Select an element (**C, H, O, N, B, F, S, P, Cl, Xe**) from the element palette, or hover over any atom on the canvas and press its letter key on your keyboard.
- **Quick Keyboard Shortcuts** (while hovering over the 2D canvas):
  - `1`, `2`, `3`: Switch to Single, Double, or Triple bond tool
  - `D` or `4`: Switch to **Dashed Bond** (`┄ Dash`) tool
  - `W`: Switch to **Wedge Bond** (`▲ Wedge`) tool
  - `V`: Move / Select tool (drag atoms to reposition them)
  - `Del` or `Backspace`: Delete the hovered atom or bond
  - `Ctrl+Z` / `Ctrl+Shift+Z`: Undo / Redo

> **Automatic 3D Geometry Reconstruction**: Even though you draw on a flat 2D canvas, SymmOrbital automatically calculates central-atom lone pairs ($AX_nE_m$) and reconstructs true 3D VSEPR geometry (tetrahedral, trigonal pyramidal, seesaw, T-shaped, octahedral, etc.) in the 3D viewer!

---

## 3. Explore Molecular Orbitals & Energy Levels (Middle Pane)

The **Molecular Orbital Energy Correlation Diagram** in the middle shows how individual Atomic Orbitals (left and right columns) combine to form Molecular Orbitals (center column):

- **Inspect Any Orbital**:
  - **Click** any horizontal energy bar in the diagram (or any row in the table below it) to lock that orbital into the 3D viewer.
  - **Hover** over an energy bar to preview its wavefunction in 3D and highlight the contributing atomic orbitals.
- **Key Frontier Orbitals**:
  - **HOMO** (Highest Occupied Molecular Orbital) / **SOMO** (Singly Occupied MO for radicals): Highlighted in blue.
  - **LUMO** (Lowest Unoccupied Molecular Orbital): Highlighted in amber.
  - **$\Delta E$ Gap**: Shown in the top badge—larger gaps generally indicate greater electronic stability.
- **Interactive Controls** (`Options` button in the MO header):
  - **Charge ($q$)**: Add or remove electrons (`-2` to `+2`) to see how ionization or reduction changes the HOMO, LUMO, and bond order.
  - **Ground vs. Excited State ($h\nu$)**: Simulate promoting an electron from the HOMO to the LUMO.
  - **SALC Lines**: Toggle the dashed correlation lines connecting atomic orbitals to molecular orbitals.

---

## 4. Interact with the 3D Orbital & Symmetry Stage (Right Pane)

The **3D Orbital & Symmetry Stage** renders the molecule's 3D structure alongside its quantum-mechanical wavefunction ($\psi_i(\mathbf{r})$):

- **Mouse & Touch Controls**:
  - **Left-Click + Drag**: Rotate the 3D molecule and orbital in any direction.
  - **Scroll Wheel / Pinch**: Zoom in and out.
  - **Right-Click + Drag**: Pan the camera.
- **Understand the Colors**:
  - **Blue Lobes ($\psi > 0$)**: Positive wavefunction phase.
  - **Red/Amber Lobes ($\psi < 0$)**: Negative wavefunction phase.
  - Where blue and red lobes meet without overlapping, there is a **nodal plane** (zero electron density).
- **View Options** (`View Options` button in the 3D header):
  - **Nodal Plane**: Toggle the translucent nodal surface slicing through the orbital.
  - **Symmetry Axes**: Toggle 3D rotation axes ($C_n$) and mirror planes ($\sigma$) directly inside the 3D scene.
  - **Lobe Opacity**: Adjust the transparency slider to see the underlying atomic skeleton clearly.
  - **Animate $C_n$**: Spin the molecule around its principal symmetry axis.

---

## 5. Bottom Reference Dock (Character Tables & Experimental Verification)

At the bottom of the screen, click any tab to expand the reference drawer (click the active tab again to collapse it):

1. **Point Group & Character Table**:
   - Displays the molecule's Schoenflies point group ($C_{2v}, C_{3v}, C_{4v}, D_{3h}, D_{4h}, O_h$, etc.) and Mulliken character table.
   - **Click any symmetry row** (e.g., $A_1, B_1, E, T_{1u}$) to highlight all molecular orbitals belonging to that symmetry species.
   - Click a symmetry operation button (like **Rotate $C_n$** or **Reflect $\sigma$**) to watch that operation animate live on the 3D molecule.
2. **NIST CCCBDB & PES Audit**:
   - Compares SymmOrbital's calculated orbital energies against real laboratory photoelectron spectroscopy (PES) ionization energies from **NIST CCCBDB**.
3. **PubChem API Inspector**:
   - Shows the live IUPAC name, SMILES string, PubChem CID, molecular weight, and whether 3D coordinates came from a PubChem MMFF94s 3D conformer or VSEPR 3D reconstruction.

---

## 6. Workspace Customization (Top Bar)

- **Dark / Light Mode**: Click the **Moon / Sun icon** in the top navigation bar to switch between Light and Dark themes.
- **Pro 3-Pane vs. Classic Layout**: Click **Classic Layout** in the top-right corner if you prefer a vertically stacked view with larger individual panels, or switch back to **Pro 3-Pane** to see 2D sketching, MO energy levels, and 3D orbitals side by side.
