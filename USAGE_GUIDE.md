# SymmOrbital — Simplified Usage Guide

Welcome to **SymmOrbital**! This quick-start guide walks you through every workspace panel, interactive control, chemical drawing tool, and 3D orbital visualization feature in the application.

---

## 1. Searching or Selecting a Molecule (Top Bar & Benchmark Library)

At the top of the screen, you can load any molecule or coordination complex in two ways:

1. **Live PubChem Search**:
   - Type any chemical name (e.g., `hexacyanoferrate`, `chromium hexacarbonyl`, `hydrogen fluoride`, `caffeine`, `benzene`), formula (e.g., `HF`, `HCl`, `Fe(CO)5`, `SiH4`), SMILES string, or PubChem CID into the top search box and press **Search PubChem**.
   - SymmOrbital fetches the structure from the **PubChem PUG REST API**, automatically infers any dative metal–ligand coordination bonds (such as the 6 $\mathrm{Fe\text{--}C}$ bonds in **hexacyanoferrate $[\mathrm{Fe(CN)_6}]^{4-}$** or $\mathrm{Cr\text{--}C}$ bonds in **$\mathrm{Cr(CO)_6}$**), strips spectator counterions, reconstructs full 3D VSEPR / coordination geometry, detects its **3D Schoenflies Point Group**, and builds its Molecular Orbital (MO) correlation diagram.
2. **Benchmark Molecules Drawer (6 Category Tabs)**:
   - Click **Benchmark Molecules** in the top bar to open the curated **NIST PES, CCCBDB, Coordination & Biomolecular Standards** drawer:
     - **d & f Coordination Complexes ($T_d / D_{4h} / O_h$)**: $\mathrm{TiCl_4}$ ($3d^0$ tetrahedral), $[\mathrm{PtCl_4}]^{2-}$ ($5d^8$ square planar), $[\mathrm{CrCl_6}]^{3-}$ ($3d^3$ high-spin octahedral quartet with Hund's rule $2t_{2g}^3$ SOMOs), and $[\mathrm{CeCl_6}]^{3-}$ ($4f^1$ lanthanide octahedral complex with an 8-lobe $4f_{xyz}$ $1a_{2u}$ SOMO).
     - **Hypervalent & Lone-Pair VSEPR**: $\mathrm{XeF_6}$, $\mathrm{XeF_4}$, $\mathrm{XeF_2}$, $\mathrm{XeOF_4}$, $\mathrm{SF_4}$, $\mathrm{ClF_3}$, $\mathrm{PF_5}$, $\mathrm{SF_6}$.
     - **Diatomics & Radicals ($C_{\infty v} / D_{\infty h}$)**: $\mathrm{CO}$, $\mathrm{NO}$ (11e⁻ radical), $\mathrm{O_2}$ (triplet diradical), $\mathrm{N_2}$ ($s\text{-}p$ mixing), $\mathrm{CO_2}$.
     - **Main-Group VSEPR**: $\mathrm{H_2O}$, $\mathrm{NH_3}$, $\mathrm{CH_4}$, $\mathrm{BF_3}$, $\mathrm{PF_5}$, $\mathrm{SF_6}$, $\mathrm{XeF_2}$, $\mathrm{XeOF_4}$.
     - **Conjugated $\pi$-Systems**: $\mathrm{C_2H_4}$ (Ethylene), $\mathrm{CH_2O}$ (Formaldehyde), $\mathrm{C_6H_6}$ (Benzene), $\mathrm{C_{40}H_{56}}$ ($\beta$-Carotene).
     - **Large Molecules**: **$\beta$-Carotene ($\mathrm{C_{40}H_{56}}$, 96 atoms)** and **Adenosine Triphosphate / ATP ($\mathrm{C_{10}H_{16}N_5O_{13}P_3}$, 47 atoms)**.

---

## 2. Drawing Custom Molecules & Using the Full 2D Chemical Sketcher (Left Pane)

Use the **2D Sketcher** on the left to draw or modify any chemical structure, or click **`Full 2D Sketcher ↗`** in the lower-right corner of the **Tools & Elements** box to open the full **PubChem-Style Draw Structure Modal**:

- **Two-Click Atom Bonding & Carbon Placement (`1`, `2`, `3` Bond Modes)**:
  - Select **Single (`1`)**, **Double (`2`)**, or **Triple (`3`)** bond mode.
  - **First Click on an Atom**: Selects that atom as the bond start point (highlighted with a blue selection ring).
  - **Second Click on Another Atom**: Immediately adds (or updates to) a single, double, or triple bond between the two selected atoms.
  - **Second Click on Empty Space**: If your second click is not on an existing atom, a standalone **Carbon ($\mathrm{C}$) atom** is placed at that canvas position without adding a bond.
  - **Click-and-Drag**: You can also click and drag from any existing atom to sprout a new bonded atom at the release point, or click an existing bond to cycle or update its bond order.
- **All 86 Elements ($\mathrm{H}$ through $\mathrm{Rn}$)**:
  - Pick common elements from the quick palette, click **Periodic Table**, or open **`Full 2D Sketcher ↗`** to access all **86 elements**—including $3d/4d/5d$ transition metals ($\mathrm{Sc\text{--}Zn}$, $\mathrm{Y\text{--}Cd}$, $\mathrm{Lu\text{--}Hg}$) and $4f$ lanthanides ($\mathrm{La\text{--}Lu}$).
- **Full 2D Sketcher Modal Features (`Full 2D Sketcher ↗`)**:
  - **Top Live String Bar**: Type or paste any **SMILES**, **Molecular Formula** (e.g., `SiH4`, `HBr`, `NdF3`), or **MDL Molfile** and click **Apply** to load it onto the canvas.
  - **Transform & Query Tools**: **New**, **Undo (`Udo`)**, **Clean (`Cln`)**, **Simplify (`Sfy`)**, **Delete (`Del`)**, **Query (`Qry`)**, **Move (`✥`)**, **Rotate (`↻` by $30^\circ$)**, and **Flip Vertical / Horizontal**.
  - **Rings, Chains & Functional Groups**: Stamp 3- to 8-membered rings, Benzene, Naphthalene, 2- to 6-carbon chains, branched `t-Bu` / neopentyl groups, or `-CHO`, `-CO2H`, `-NO2`, and `-SO3H` groups with one click.
  - **Charges & Export**: Increment/decrement atomic formal charges (`⊕`, `⊖`, `⊙`), adjust hydrogen display (**Keep AsIs**, **Add Hydrogens**, **Remove Hydrogens**), and export your structure to **MDL V2000 Molfile**, **SMILES**, or **Formula** before clicking **Accept Structure**.

---

## 3. Exploring the MO Energy Correlation Diagram (Center Pane)

The **MO Energy Correlation Diagram** shows how Central Metal / Primary Atomic Orbitals ($s, p, d, 4f$) and Peripheral Ligand **SALCs** combine to form **Molecular Orbitals ($\psi_i$)**:

- **Selecting & Stepping Through Orbitals**:
  - Click any **MO energy rung** in the center column (or use the **`↑` / `↓` buttons** and keyboard **Up / Down arrow keys**) to inspect its 3D orbital wavefunction, energy ($\text{eV}$), literature-standard Mulliken symmetry label ($2\sigma^+, 1\pi, 2t_{2g}, 1a_{2u}, \dots$), and LCAO formula.
  - Degenerate $E$ doublets and $T$ triplets ($t_{2g}, t_{1u}, t_{2u}$) are grouped side-by-side on the same energy rung and filled spin-parallel according to **Hund's rule**.
- **State & Options Menu** (top-right of the MO pane):
  - **Molecular Charge ($q$)**: Add or remove electrons to observe how **Bond Order**, **HOMO/LUMO**, and **Magnetism** (Diamagnetic vs. Paramagnetic) update in real time.
  - **Electronic State**: Toggle between **Ground State** and **Excited State ($h\nu$)** to promote an electron from the HOMO to the LUMO.
  - **Full MO Level Table (`All Levels →`)**: Click **`All Levels →`** at the bottom-right of the LCAO formula bar to open a complete table of all molecular orbitals, energies, symmetries, and electron occupancies.

---

## 4. 3D Orbital Stage & Continuous AO $\to$ MO Wavefunction Morphing (Right Pane)

The **3D Orbital & Symmetry Stage** renders the selected Molecular Orbital alongside the molecule's 3D ball-and-stick framework:

- **AO $\to$ MO Continuous 3D Wavefunction Morphing (`Morph MO` & `0%–100%` Slider)**:
  - Click **`Morph MO`** in the top-right header of the 3D viewport to animate a smooth transition from **Separate Atomic Orbitals (`0%`)** to the **Merged Molecular Orbital Isosurface (`100%`)**, or drag the **`0%–100%` slider** manually:
    - **`0%` (Separate AOs)**: Displays the unmixed constituent $s, p, d, 4f$ atomic orbitals centered on their respective nuclei, weighted by their LCAO coefficients $c_\mu$.
    - **`1%–99%` (AO $\to$ MO Overlap)**: Watch constructive in-phase lobes fuse across covalent and coordinate bonds into delocalized $\sigma, \pi, \delta$ envelopes while destructive out-of-phase lobes carve out nodal planes.
    - **`100%` (Merged MO Isosurface)**: Displays the smooth, Taubin-faired, watertight 3D Molecular Orbital isosurface ($\psi(\mathbf{r}) = \sum_\mu c_\mu \phi_\mu(\mathbf{r})$).
- **Camera & 3D Display Controls**:
  - **Left-click + drag** to rotate in 3D, **Right-click + drag** to pan, and **Scroll** to zoom.
  - Open **3D Options** (or press **`H`** to toggle the bottom 3D HUD pill) to control **Nodal Planes**, **Symmetry Axes ($C_n, \sigma_v, \sigma_h, \sigma_d$)**, **Spin $C_n$**, **Reset Camera**, **Lobe Opacity (`15%–95%`)**, and **Lobe Scale (`0.5×–1.6×`)**.

---

## 5. Analytical Inspector Dock (Bottom Tabs)

Switch between the tabs along the bottom bar to inspect group-theoretical and experimental reference data:

1. **Symmetry & Character Table + $\Gamma_{\text{red}}$ Decomposition Calculator**:
   - View the complete group theory character table for the active molecule's point group, including linear/rotational, quadratic ($d$-orbital), and cubic ($4f$-orbital) basis annotations.
   - Use the **Reducible Representation ($\Gamma_{\text{red}}$) Decomposition Calculator** to decompose any preset basis (**$\Gamma_\sigma$ Bonds**, **$\Gamma_{\text{val}}$ Valence AOs**, **$\Gamma_d$ (5 $d$ Orbitals)**, **$\Gamma_f$ (7 $4f$ Orbitals)**, **$\Gamma_{3N}$ Total**, **$\Gamma_{\text{vib}}$ Vibrational**, **$\Gamma_\pi$**) or custom character vector $\chi(R)$ via the Great Orthogonality Theorem.
2. **NIST WebBook PES & CCCBDB Audit**:
   - Compare computed bond lengths ($r_e$), bond angles ($\theta_e$), first vertical ionization energies ($\text{IE}_1 = -\varepsilon_{\text{HOMO}}$), and orbital eigenvalues against experimental Photoelectron Spectroscopy (PES) literature benchmarks.
3. **PubChem PUG REST Inspector**:
   - Inspect live REST API URLs, conformer type (`3D-SDF` vs. `2D-SDF → 3D VSEPR Reconstructed`), molecular descriptors, and authoritative NCBI PubChem record links.
