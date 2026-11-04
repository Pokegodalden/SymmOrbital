# SymmOrbital — Simplified Usage Guide

Welcome to **SymmOrbital**! This quick-start guide walks you through every workspace panel, interactive control, and 3D visualization tool in the application.

---

## 1. Searching or Selecting a Molecule (Top Bar & Benchmark Library)

At the top of the screen, you can load any molecule in two ways:

1. **Live PubChem Search**:
   - Type any chemical name (e.g., `caffeine`, `ferrocene`, `aspirin`), formula, or PubChem CID into the top search box and press **Search**.
   - SymmOrbital fetches the 3D molecular coordinates directly from the **PubChem PUG REST API**, detects its 3D **Schoenflies Point Group**, and builds its Molecular Orbital (MO) diagram automatically.
2. **22 Benchmark Molecules (6 Category Tabs)**:
   - Click **Benchmark Molecules** in the top bar to open the curated **22 NIST PES, CCCBDB & Biomolecular Standards** drawer.
   - Switch between the 6 category tabs to explore classic and advanced electronic structures:
     - **Linear ($D_{\infty h} / C_{\infty v}$)**: $\mathrm{H_2}$, $\mathrm{N_2}$ ($s\text{-}p$ mixing), $\mathrm{O_2}$ (triplet diradical), $\mathrm{CO}$, $\mathrm{CO_2}$, $\mathrm{C_2H_2}$.
     - **Bent & Trigonal ($C_{2v} / D_{3h} / C_{3v} / D_{2h}$)**: $\mathrm{H_2O}$, $\mathrm{H_2CO}$, $\mathrm{B_2H_6}$ (3c-2e bridging bonds), $\mathrm{NH_3}$, $\mathrm{BF_3}$.
     - **3D Polyhedral ($T_d / O_h$)**: $\mathrm{CH_4}$, $\mathrm{P_4}$ (strained cluster), $\mathrm{SF_6}$.
     - **Hypervalent ($D_{\infty h} / D_{3h} / C_{4v}$)**: $\mathrm{XeF_2}$, $\mathrm{PF_5}$, $\mathrm{XeOF_4}$.
     - **Conjugated $\pi$-Systems ($D_{2h} / C_{2h} / D_{6h}$)**: $\mathrm{C_2H_4}$ (Ethylene, full 12-MO manifold), $\mathrm{C_4H_6}$ (1,3-Butadiene), $\mathrm{C_6H_6}$ (Benzene).
     - **Large molecules ($D_{2h} / C_s$)**: **$\beta$-Carotene ($\mathrm{C_{40}H_{56}}$, 96 atoms)** showcasing delocalized 11-double-bond polyene $\pi \to \pi^*$ frontier orbitals, and **Adenosine Triphosphate / ATP ($\mathrm{C_{10}H_{16}N_5O_{13}P_3}$, 47 atoms)** showcasing adenine $\pi$ vs. triphosphate chain orbitals.

---

## 2. Drawing Custom Molecules (Left Pane: 2D Sketcher)

Use the **2D Sketcher** on the left to draw or modify any chemical structure:

- **Drawing Tools**:
  - **Select / Move (`V`)**: Click and drag atoms to rearrange geometry.
  - **Atom (`A`)**: Choose an element from the quick-pick bar ($\mathrm{C, N, O, F, P, S, Cl, H, Xe}$) or click **Periodic Table** for all elements, then click the canvas to place or replace atoms.
  - **Bonds (`1`, `2`, `3`)**: Draw Single, Double, or Triple covalent bonds between atoms.
  - **Stereo Wedges (`W` / `D` or `4`)**: Draw **Solid Up Wedges (`W`)** or **Hashed Down Wedges (`D` / `4`)** to define 3D tetrahedral stereocenters, trigonal pyramidal lone pairs, or non-planar geometries.
  - **Rings**: Stamp Cyclopropane, Cyclobutane, Cyclopentane, Cyclohexane, or Benzene rings with one click.
  - **Add Hydrogens (`+H`)**: Automatically satisfies open valence shells on all heavy atoms.
- **Automatic 2D $\to$ 3D VSEPR Reconstruction**: Every edit on the 2D canvas immediately reconstructs a 3D VSEPR geometry, updates the point group badge, and recalculates the MO diagram.

---

## 3. Exploring the MO Energy Correlation Diagram (Center Pane)

The **MO Energy Correlation Diagram** displays how Left/Central Atomic Orbitals and Right/Ligand Symmetry-Adapted Linear Combinations (**SALCs**) combine to form **Molecular Orbitals ($\psi_i$)** in the center column:

- **Selecting & Stepping Through Orbitals**:
  - Click any **MO energy rung** in the center column to inspect its 3D orbital wavefunction, energy ($\text{eV}$), Mulliken symmetry label ($a_{1g}, e_{1g}, t_2, \pi_u, \dots$), and LCAO formula at the bottom of the diagram.
  - Use the **`↑` and `↓` buttons** in the MO Diagram header (or your keyboard's **Up / Down arrow keys**) to step smoothly up and down the energy ladder from core $\sigma$ orbitals through the **HOMO–LUMO frontier** to virtual $\pi^*$ and $\sigma^*$ orbitals.
- **State & Options Menu** (top-right of the MO pane):
  - **Molecular Charge ($q$)**: Add or remove electrons ($-2$ to $+2$) to observe how **Bond Order**, **HOMO/LUMO**, and **Magnetism** (Diamagnetic vs. Paramagnetic) change in real time.
  - **Electronic State**: Switch between **Ground State** and **Excited State ($h\nu$)** to promote an electron from the HOMO to the LUMO.
  - **SALC Correlation Lines**: Toggle the dashed symmetry-matching connector lines between AOs/SALCs and MOs.
  - **Full MO Level Table (`All Levels →`)**: Toggle inside **State & Options** or click **`All Levels →`** at the bottom-right of the LCAO formula bar to open a complete table of all molecular orbitals, energies, symmetries, and electron occupancies.

---

## 4. 3D Orbital Stage & Continuous AO $\to$ MO Wavefunction Morphing (Right Pane)

The **3D Orbital & Symmetry Stage** renders the selected Molecular Orbital alongside the molecule's 3D ball-and-stick framework:

- **AO $\to$ MO Continuous 3D Wavefunction Morphing (`Morph MO` & `0%–100%` Slider)**:
  - Click **`Morph MO`** in the top-right header of the 3D viewport to play a smooth, continuous 3D quantum animation from **Separate Atomic Orbitals (`0%`)** to the **Merged Molecular Orbital Isosurface (`100%`)**, or drag the **`0%–100%` slider** manually to inspect any intermediate stage:
    - **`0%` (Separate AOs)**: Displays the unmixed constituent $s, p, d$ atomic orbitals centered on their respective nuclei, weighted by their LCAO coefficients $c_\mu$.
    - **`1%–99%` (AO $\to$ MO Overlap)**: Watch constructive in-phase lobes (same color/sign) reach across bonds and fuse together into delocalized $\sigma$ or $\pi$ envelopes, while destructive out-of-phase lobes (opposite color/sign) repel and carve out nodal planes.
    - **`100%` (Merged MO Isosurface)**: Displays the fully formed, watertight 3D Molecular Orbital isosurface ($\psi(\mathbf{r}) = \sum_\mu c_\mu \phi_\mu(\mathbf{r})$).
- **Camera Controls**:
  - **Left-click + drag** to rotate in 3D, **Right-click + drag** to pan, and **Scroll** to zoom in or out.
- **3D Options Menu** (top-right of the 3D viewport):
  - **Nodal Plane**: Toggle translucent purple nodal surfaces where the wavefunction $\psi = 0$ changes phase sign.
  - **Symmetry Axes**: Overlay principal rotation axes ($C_n$) and mirror planes ($\sigma_v, \sigma_h, \sigma_d$).
  - **Spin $C_n$**: Continuously rotate the molecule around its principal symmetry axis to verify whether the orbital is symmetric or antisymmetric.
  - **Reset Camera**: Snap the 3D camera back to the default orientation and zoom.
  - **Lobe Opacity (`15%–95%`) & Lobe Scale (`0.5×–1.6×`)**: Adjust the translucency and spatial size of the positive (`+`, red) and negative (`−`, blue) orbital lobes. All lobes remain 100% closed and unclipped even at maximum `1.6×` scale.
- **Show / Hide 3D Controls (`H` Key)**:
  - Press **`H`** on your keyboard (or click the floating pill at the bottom of the 3D stage) to expand or collapse the quick 3D camera and display controls bar.

---

## 5. Analytical Inspector Dock (Bottom Collapsible Drawer)

Click any of the three tabs along the bottom bar to expand or collapse the analytical reference drawer:

1. **Point Group & Character Table + $\Gamma_{\text{red}}$ Decomposition Calculator**:
   - View the complete group theory character table for the active molecule's point group, with the selected MO's irreducible representation highlighted.
   - **Reducible Representation ($\Gamma_{\text{red}}$) Decomposition Calculator**: Scroll to the calculator inside this tab to decompose any reducible representation into its lower-case Mulliken irreducible direct sum ($\Gamma_{\text{red}} = \bigoplus_i a_i \Gamma_i$) via the Great Orthogonality Theorem. Click a preset basis (**$\Gamma_\sigma$ Bonds**, **$\Gamma_{\text{val}}$ Valence AOs**, **$\Gamma_{3N}$ Total**, **$\Gamma_{\text{vib}}$ Vibrational**, **$\Gamma_\pi$ Out-of-Plane**) or type custom characters $\chi(R)$ for each symmetry class.
2. **NIST CCCBDB & PES Audit**:
   - Compare computed bond lengths ($r_e$), bond angles ($\theta_e$), first vertical ionization energies ($\text{IE}_1 = -\varepsilon_{\text{HOMO}}$), and orbital eigenvalues against experimental Photoelectron Spectroscopy (PES) literature values.
3. **PubChem API Inspector**:
   - Inspect live REST API URLs, conformer IDs, response status, and raw JSON payloads returned by PubChem.
