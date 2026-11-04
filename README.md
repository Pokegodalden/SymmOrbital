<div align="center">
  <h1>SymmOrbital</h1>
  <p><strong>Molecular Orbital & Symmetry Analyzer — Live PubChem API, Interactive 2D/3D Chemical Sketcher, Continuous AO &rarr; MO 3D Wavefunction Morphing & NIST Benchmark Suite</strong></p>
  <p><a href="./USAGE_GUIDE.md"><strong>&#128214; Simplified Usage Guide</strong></a></p>
</div>

---

## Overview

**SymmOrbital** is a computational and educational chemistry platform designed for visualizing **Molecular Orbital (MO) Energy Correlation Diagrams**, **Linear Combination of Atomic Orbitals (LCAO-MO)**, **Continuous 3D Atomic-to-Molecular Orbital Wavefunction Morphing**, **Symmetry-Adapted Linear Combinations (SALCs)**, and **3D Schoenflies Point Group Symmetry Elements**.

The application combines a **live PubChem REST API structural resolver**, a **2D/3D chemical structure editor** with stereochemical wedge/dash bonds and VSEPR 3D coordinate reconstruction, a **distance-scaled Extended Hückel / LCAO-MO engine** (supporting $s, p, d$ valence basis sets, hypervalent 3c-4e bonding, conjugated $\pi$-systems, and large macromolecular biomarkers), a **Marching Tetrahedra 3D scalar-field isosurface engine** for real-time AO $\to$ MO topological morphing, and **curated NIST Photoelectron Spectroscopy (PES) & CCCBDB benchmark datasets** for 22 reference molecules.

---

## Key Features

### 1. Live PubChem PUG REST API Integration
- **Universal Compound Search**: Query over 115 million chemical structures by common name, IUPAC name, chemical formula, or PubChem CID.
- **3D & 2D Coordinate Resolution**: Automatically fetches MMFF94-minimized 3D conformers (`record_type=3d`) from PubChem, with seamless fallback to 2D structural records and automatic VSEPR 3D coordinate generation when 3D conformers are unavailable.
- **API Payload Inspector**: Inspect live PubChem request URLs, response latency, conformer metadata, and raw JSON payloads directly inside the bottom **PubChem API Inspector** tab.

### 2. Interactive 2D Chemical Structure Editor & 2D $\to$ 3D VSEPR Reconstructor
- **Full Drawing Suite**: Select/Move, Atom Placement ($\mathrm{C, N, O, F, P, S, Cl, Br, I, H, B, Si, Xe, Fe, Pt}$ + Periodic Table modal), Single/Double/Triple covalent bonds, **Up Wedge (`W`)** and **Down Hashed Wedge (`D` / `4`)** stereochemical bonds, Eraser, and Ring Templates (Cyclopropane, Cyclobutane, Cyclopentane, Cyclohexane, Benzene).
- **Stereocenter-Aware 3D VSEPR Builder**: Converts 2D skeletal drawings and wedge/dash stereocenters into 3D Cartesian coordinates ($\text{\AA}$) on the fly, automatically updating the Schoenflies point group, character table, and MO correlation diagram.

### 3. Distance-Scaled Extended Hückel & Symmetry-Adapted LCAO-MO Engine
- **Curated NIST PES Benchmarks**: Exact experimental/CCCBDB vertical ionization and orbital energy levels, Mulliken irreducible representations ($\Gamma_i$), and symmetry-adapted LCAO wavefunctions for **22 benchmark molecules** spanning linear, bent, trigonal, tetrahedral, hypervalent, octahedral, conjugated $\pi$-systems, and large biomolecules.
- **Large Molecule & Conjugated Polyene Support**:
  - **Frontier & Delocalized $\pi$-Subspace Prioritization**: Automatically retains all conjugated $p_\pi$ and heteroatom lone-pair orbitals alongside top-weighted $\sigma/\sigma^*$ framework orbitals so large systems like **$\beta$-Carotene ($\mathrm{C_{40}H_{56}}$, 96 atoms)** and **Adenosine Triphosphate ($\mathrm{C_{10}H_{16}N_5O_{13}P_3}$, 47 atoms)** preserve their physically meaningful HOMO–LUMO frontier gap and delocalized polyene/phosphate bands.
  - **Occupied/Virtual Window Slicing**: Centers the interactive correlation diagram around the frontier HOMO–LUMO manifold (`10 occupied + 10 virtual` levels) so energy rungs remain clear and legible regardless of molecular size.
- **Electronic State & Charge Controls**:
  - Adjust molecular **Charge ($q \in [-2, +2]$)** to model radical cations, anions, and isoelectronic series.
  - Toggle **Ground State vs. Excited State ($h\nu$, HOMO $\to$ LUMO promotion)** with automatic spin-multiplicity ($2S+1$), magnetic state (Diamagnetic vs. Paramagnetic), and **Bond Order** recalculation.
  - Step cleanly through molecular orbitals using the header **`↑` / `↓` orbital stepping buttons** (or keyboard arrow keys) and toggle the **Full MO Level Table (`All Levels →`)** for complete eigenvalue and occupancy inspection.

### 4. Continuous 3D AO $\to$ MO Wavefunction Morphing & Symmetry Stage
- **Real-Time Scalar-Field Wavefunction Morphing (`Morph MO` & `0%–100%` Slider)**:
  - Play/Pause or scrub the interactive **`0%–100%` morph slider** in the 3D viewport header to watch **Separate Atomic Orbitals (`0%`)** continuously overlap, interfere, and fuse into the **Merged Molecular Orbital Isosurface (`100%`)**.
  - Powered by a **6-tetrahedra Kuhn Marching Tetrahedra** isosurface extractor over a **watertight boundary-windowed 3D scalar field**:
    $$F(\mathbf{r}, t) = (1 - t)\,F_{\text{AO}}(\mathbf{r}) + t\,F_{\text{MO}}(\mathbf{r})$$
    where $F_{\text{AO}}(\mathbf{r}) = \operatorname{sign}(\psi_{\mu^*})\sqrt{\sum_\mu (c_\mu \phi_\mu(\mathbf{r}))^2}$ keeps individual constituent atomic orbitals spatially distinct at $t = 0$, and $F_{\text{MO}}(\mathbf{r}) = \sum_\mu c_\mu \phi_\mu(\mathbf{r})$ evaluates the true LCAO-MO interference sum at $t = 1$.
  - **Constructive & Destructive Quantum Interference**: Bonding combinations with matching signs ($c_\mu c_\nu > 0$) lower their saddle barrier and topologically bridge across covalent bonds, whereas antibonding combinations ($c_\mu c_\nu < 0$) carve out zero-crossing nodal planes between nuclei.
  - **Tangent-Plane Laplacian Surface Fairing & Analytical $\nabla F$ Normals**: Projects local Laplacian smoothing onto the isosurface tangent plane $(\mathbf{I} - \mathbf{n}\mathbf{n}^\top)\Delta\mathbf{p}$ and evaluates central-difference analytical gradients $\nabla F(\mathbf{r})$ for smooth Phong shading and 100% closed, unclipped 2-manifold lobes across all scale settings (`0.5×–1.6×`).
- **3D Options & Symmetry Overlays**:
  - Toggle **Nodal Planes**, **Symmetry Axes ($C_n, \sigma_v, \sigma_h, \sigma_d$)**, **Continuous $C_n$ Rotation (`Spin Cn`)**, **Reset Camera**, **Lobe Opacity (`15%–95%`)**, and **Lobe Scale (`0.5×–1.6×`)**.
  - Press **`H`** (or click **Show/Hide 3D Controls**) to toggle the floating bottom 3D HUD pill.

### 5. Group Theory Character Tables, $\Gamma_{\text{red}}$ Decomposition & NIST Audit
- **Point Group & Character Table Inspector**: Full Schoenflies character tables ($C_1, C_s, C_i, C_2, C_{2v}, C_{3v}, C_{4v}, C_{2h}, D_2, D_{2h}, D_{3h}, D_{4h}, D_{6h}, D_{2d}, D_{3d}, C_{\infty v}, D_{\infty h}, T_d, O_h$) with conjugacy classes, Mulliken symbols, linear/rotational ($x, y, z, R_x, R_y, R_z$) and quadratic ($x^2, z^2, xy, xz, yz$) basis functions, and dynamically highlighted active MO irreps.
- **Interactive Reducible Representation ($\Gamma_{\text{red}}$) Decomposition Calculator**:
  - Decompose any reducible representation into its direct sum of irreducible representations ($\Gamma_{\text{red}} = \bigoplus_i a_i \Gamma_i$) using the **Great Orthogonality Theorem**:
    $$a_i = \frac{1}{h} \sum_R N_R \, \chi_{\text{red}}(R) \, \chi_i^*(R)$$
  - Select from built-in molecular basis presets ($\Gamma_\sigma$, $\Gamma_{\text{val}}$, $\Gamma_{3N}$, $\Gamma_{\text{trans+rot}}$, $\Gamma_{\text{vib}}$, $\Gamma_\pi$) or enter custom characters $\chi(R)$ per conjugacy class with automatic lower-case Mulliken direct-sum formatting and step-by-step inner-product work.
- **NIST CCCBDB & PES Audit Tab**: Compares computed geometry parameters ($r_e, \theta_e$), first vertical ionization potentials ($\text{IE}_1 = -\varepsilon_{\text{HOMO}}$), and Koopmans' theorem orbital energies against experimental Photoelectron Spectroscopy (PES) literature references (Turner, Kimura, Potts, Siegbahn, Bieri).

---

## Curated Benchmark Molecules (22 NIST PES, CCCBDB & Biomolecular Standards)

| Category | Molecule | Formula | Point Group | Highlighted Electronic & Symmetry Features |
| :--- | :--- | :--- | :---: | :--- |
| **Linear** | Dihydrogen | $\mathrm{H_2}$ | $D_{\infty h}$ | Canonical 2-center $1\sigma_g / 1\sigma_u^*$ covalent bond |
| **Linear** | Dinitrogen | $\mathrm{N_2}$ | $D_{\infty h}$ | $s\text{-}p$ mixing inversion ($3\sigma_g$ HOMO above $1\pi_u$) |
| **Linear** | Dioxygen | $\mathrm{O_2}$ | $D_{\infty h}$ | Paramagnetic triplet diradical ($^3\Sigma_g^-$, half-filled $1\pi_g^*$) |
| **Linear** | Carbon Monoxide | $\mathrm{CO}$ | $C_{\infty v}$ | Heteronuclear polarization; carbon-centered $5\sigma$ HOMO lone pair |
| **Linear** | Carbon Dioxide | $\mathrm{CO_2}$ | $D_{\infty h}$ | 3-center 4-electron $\pi$-system ($1\pi_u, 1\pi_g$ non-bonding HOMO) |
| **Linear** | Acetylene | $\mathrm{C_2H_2}$ | $D_{\infty h}$ | Degenerate cylindrical $1\pi_u$ triple-bond manifold |
| **Bent / Planar** | Water | $\mathrm{H_2O}$ | $C_{2v}$ | Non-equivalent oxygen lone pairs ($1b_1$ pure $p_x$ vs. $3a_1$ $sp^2$-hybrid) |
| **Bent / Planar** | Formaldehyde | $\mathrm{H_2CO}$ | $C_{2v}$ | Classic $n \to \pi^*$ ($2b_2 \to 2b_1$) carbonyl chromophore |
| **Bent / Planar** | Diborane | $\mathrm{B_2H_6}$ | $D_{2h}$ | Electron-deficient 3-center 2-electron ($b_{2u}, a_g$) bridging $\mathrm{B\text{--}H\text{--}B}$ bonds |
| **Pyramidal / Planar** | Ammonia | $\mathrm{NH_3}$ | $C_{3v}$ | $3a_1$ pyramidal lone pair and degenerate $1e$ $\mathrm{N\text{--}H}$ bonding pair |
| **Pyramidal / Planar** | Boron Trifluoride | $\mathrm{BF_3}$ | $D_{3h}$ | Lewis acidity and $1a_2''$ out-of-plane $\pi$-backbonding |
| **3D Polyhedral** | Methane | $\mathrm{CH_4}$ | $T_d$ | Triply degenerate $1t_2$ photoemission bands ($1a_1 + 1t_2$, no $d$-orbital involvement) |
| **3D Polyhedral** | White Phosphorus | $\mathrm{P_4}$ | $T_d$ | Strained tetrahedral cage with bent "banana" $\sigma$-cluster orbitals ($2a_1, 2t_2, 1e$) |
| **3D Polyhedral** | Sulfur Hexafluoride | $\mathrm{SF_6}$ | $O_h$ | Hypervalent 12-electron octahedral bonding ($1e_g$ non-bonding ligand SALCs) |
| **Hypervalent** | Xenon Difluoride | $\mathrm{XeF_2}$ | $D_{\infty h}$ | Classic Pimentel-Rundle 3-center 4-electron ($3\sigma_u$ non-bonding) linear hypervalency |
| **Hypervalent** | Phosphorus Pentafluoride | $\mathrm{PF_5}$ | $D_{3h}$ | Trigonal bipyramidal axial vs. equatorial bonding and Berry pseudorotation |
| **Hypervalent** | Xenon Oxytetrafluoride | $\mathrm{XeOF_4}$ | $C_{4v}$ | Square pyramidal $d^2sp^3$ / 3c-4e basal $\mathrm{XeF_4}$ vs. axial $\mathrm{Xe{=}O}$ $\pi$-bond |
| **Conjugated $\pi$** | Ethylene | $\mathrm{C_2H_4}$ | $D_{2h}$ | Full 12-valence-MO $D_{2h}$ manifold and frontier $1b_{3u}(\pi) \to 1b_{2g}(\pi^*)$ gap |
| **Conjugated $\pi$** | 1,3-Butadiene | $\mathrm{C_4H_6}$ | $C_{2h}$ | 4-center conjugated $\pi$-manifold ($1a_u, 1b_g, 2a_u^*, 2b_g^*$) |
| **Conjugated $\pi$** | Benzene | $\mathrm{C_6H_6}$ | $D_{6h}$ | Hückel $6\pi$ aromaticity ($1a_{2u}$, degenerate $1e_{1g}$ HOMO & $1e_{2u}$ LUMO) |
| **Large Molecules** | $\beta$-Carotene | $\mathrm{C_{40}H_{56}}$ | $D_{2h}$ | 96-atom conjugated 11-double-bond polyene ($22\pi$ electron bathochromic $\pi \to \pi^*$ chromophore) |
| **Large Molecules** | Adenosine Triphosphate | $\mathrm{C_{10}H_{16}N_5O_{13}P_3}$ | $C_s$ | 47-atom bioenergetic coenzyme with adenine $\pi$-ring HOMO and triphosphate anhydride LUMO |

---

## Tech Stack & Verification

- **Frontend**: React 19, TypeScript, Tailwind CSS 4, Lucide Icons, KaTeX mathematical typesetting
- **3D Rendering**: Three.js (`WebGLRenderer`, custom Marching Tetrahedra scalar-field isosurface generator, tangent-plane Laplacian surface fairing, analytical $\nabla F$ vertex normals)
- **Chemical Data**: Live PubChem PUG REST API + Curated NIST CCCBDB / PES Benchmark Library
- **Automated Test Suite**: Run `npx tsx --test src/chemistry/moEngine.test.ts` to verify all 22 benchmark eigenvalues, Mulliken symmetries, watertight 3D Marching Tetrahedra isosurfaces (`0` open boundary edges), and PubChem fallback calculations.
