<div align="center">
  <h1>SymmOrbital</h1>
  <p><strong>Molecular Orbital & Symmetry Analyzer — Live PubChem API, Full 86-Element 2D/3D Chemical Sketcher, s–p–d–f Extended Hückel Engine, Continuous AO &rarr; MO 3D Wavefunction Morphing & NIST Benchmark Suite</strong></p>
  <p><a href="./USAGE_GUIDE.md"><strong>&#128214; Simplified Usage Guide</strong></a></p>
</div>

---

## Overview

**SymmOrbital** is a computational and educational chemistry platform designed for visualizing **Molecular Orbital (MO) Energy Correlation Diagrams**, **Linear Combination of Atomic Orbitals (LCAO-MO)**, **Continuous 3D Atomic-to-Molecular Orbital Wavefunction Morphing**, **Symmetry-Adapted Linear Combinations (SALCs)**, and **3D Schoenflies Point Group Symmetry Elements**.

The application combines a **live PubChem REST API structural resolver**, a **dual 2D/3D chemical structure editor** (including a full **PubChem-style Draw Structure modal** supporting all **86 elements from $\mathrm{H}$ to $\mathrm{Rn}$**), an **automatic coordination-complex bond & 3D VSEPR geometry reconstructor**, a **distance-scaled $s\text{--}p\text{--}d\text{--}f$ Extended Hückel / LCAO-MO engine** (supporting main-group, hypervalent 3c-4e, transition-metal $3d/4d/5d$, and lanthanide $4f$ coordination complexes with Hund's rule spin-parallel filling), a **Taubin-faired Marching Tetrahedra 3D scalar-field isosurface engine** for real-time AO $\to$ MO topological morphing, and **curated NIST Photoelectron Spectroscopy (PES) & CCCBDB benchmark datasets** for 22 reference molecules and complexes.

---

## Key Features

### 1. Live PubChem PUG REST API Integration & Coordination Complex Resolver
- **Universal Compound Search**: Query over 115 million chemical structures by common name, IUPAC name, chemical formula, SMILES, or PubChem CID (e.g., `hexacyanoferrate`, `chromium hexacarbonyl`, `nitric oxide`, `hydrogen fluoride`, `benzene`).
- **3D & 2D Coordinate Resolution with Principal-Axis Alignment**: Automatically fetches MMFF94-minimized 3D conformers (`record_type=3d`) from PubChem, aligns tilted diatomic/linear conformers (such as $\mathrm{HF}$ and $\mathrm{HCl}$) along the canonical molecular axis, and seamlessly falls back to 2D structural records with automatic 3D VSEPR/coordination coordinate reconstruction when 3D conformers are unavailable.
- **Automatic Dative Metal–Ligand Bond Inference (`resolveCoordinationBonds`)**: Detects disconnected coordination complexes in PubChem 2D-SDFs, MDL Molfiles, and multi-component SMILES (such as **Hexacyanoferrate(II) $[\mathrm{Fe(CN)_6}]^{4-}$**, **Ferricyanide $[\mathrm{Fe(CN)_6}]^{3-}$**, **$\mathrm{Cr(CO)_6}$**, **$\mathrm{Fe(CO)_5}$**, **$\mathrm{Ni(CO)_4}$**, **$[\mathrm{Ni(CN)_4}]^{2-}$**, **$[\mathrm{Fe(CN)_5(NO)}]^{2-}$**, and **$[\mathrm{PtCl_6}]^{2-}$**), bonds each ligand via its $\sigma$-donor atom ($\mathrm{C}$ for $-\mathrm{CN}$/$-\mathrm{CO}$/$-\mathrm{CS}$; $\mathrm{N}$ for $-\mathrm{NO}$/$-\mathrm{NH_3}$), strips spectator alkali/alkaline-earth counterions ($\mathrm{K^+, Na^+, Li^+}$), and reconstructs full 3D Octahedral ($O_h$), Trigonal Bipyramidal ($D_{3h}$), Square Planar ($D_{4h}$), Square Pyramidal ($C_{4v}$), or Tetrahedral ($T_d$) geometries.
- **API Payload Inspector**: Inspect live PubChem request URLs, response latency, conformer metadata, and raw JSON payloads directly inside the bottom **PubChem PUG REST Inspector** tab.

### 2. Interactive 2D Chemical Sketcher & Full PubChem-Style Draw Structure Modal
- **Quick 2D Workbench Sketcher + `Full 2D Sketcher ↗` Modal**:
  - Use the left-pane **2D Sketcher** for rapid editing or click **`Full 2D Sketcher ↗`** in the lower-right of the **Tools & Elements** panel to launch the full-featured **PubChem Draw Structure** modal.
  - **All 86 Elements ($\mathrm{H}$ through $\mathrm{Rn}$)**: Full Periodic Table coverage including alkali/alkaline-earth metals, main-group $p$-block elements, $3d/4d/5d$ transition metals ($\mathrm{Sc\text{--}Zn}$, $\mathrm{Y\text{--}Cd}$, $\mathrm{Lu\text{--}Hg}$), and $4f$ lanthanides ($\mathrm{La\text{--}Lu}$) via interactive dropdown selectors.
  - **Two-Click Atom-to-Atom Bonding, Sprouting & Replacement**: In **Single (`1`), Double (`2`), or Triple (`3`) bond mode**, click a first atom to select it (highlighted with a blue ring):
    - **Second click on another atom**: Adds or updates a single, double, or triple bond between the two atoms.
    - **Second click on empty canvas space**: Places a new atom of the **currently selected element** at the exact spot of the second click and connects it to the first atom with the selected bond.
    - **Second click on the same atom**: Replaces that atom in-place with the **currently selected element**.
    - **Click on empty canvas space (when no atom is selected)**: Places a standalone atom of the **currently selected element** (defaulting to $\mathrm{C}$) without forming a bond.
  - **Complete Structural Drawing & Import/Export Toolkit**: Includes **Up Wedge (`W`)** and **Down Hashed Wedge (`D` / `4`)** stereochemical bonds, Wavy/Crossed/Aromatic query bonds, 3- to 8-membered rings + Benzene & Naphthalene, aliphatic & branched carbon chains, functional group stamps (`-CHO`, `-CO2H`, `-NO2`, `-SO3H`), formal charges (`+`, `-`, `•`), transform tools (Clean, Simplify, Move, Rotate, Flip H/V), and live **SMILES / Formula / MDL V2000 Molfile** import and export.

### 3. Full $s\text{--}p\text{--}d\text{--}f$ Extended Hückel & Symmetry-Adapted LCAO-MO Engine
- **s, p, d, and 4f Valence Basis Sets**:
  - **5-Orbital Real $d$ Basis** ($d_{z^2}, d_{x^2-y^2}, d_{xy}, d_{xz}, d_{yz}$) for transition metals and hypervalent $p$-block centers with exact 3D Slater-Koster directional overlap transformations ($S_{dd\sigma}, S_{dd\pi}, S_{dd\delta}, S_{sd\sigma}, S_{pd\sigma}, S_{pd\pi}$).
  - **7-Orbital Real Cubic $4f$ Basis** ($f_{z^3}, f_{xz^2}, f_{yz^2}, f_{xyz}, f_{z(x^2-y^2)}, f_{x(x^2-3y^2)}, f_{y(3x^2-y^2)}$) for lanthanides ($\mathrm{La\text{--}Lu}$) alongside $5d/6s/6p$ shells, capturing cubic crystal-field splitting ($\Gamma_f = a_{2u} \oplus t_{1u} \oplus t_{2u}$ in $O_h$).
  - **Hund's Rule Spin-Parallel Filling**: Automatically distributes electrons spin-parallel ($1, 1, 1, \dots$) across degenerate $E$ doublets, $T$ triplets, and $d/f$ subshells before spin-pairing (e.g., high-spin $3d^3$ $t_{2g}^3$ in $[\mathrm{CrCl_6}]^{3-}$ and triplet $1\pi_g^*$ in $\mathrm{O_2}$).
  - **Literature-Standard Mulliken Shell Numbering**: Accounts for inner core shells when numbering valence $\sigma^+$ and $\pi$ irreps (producing exact experimental labels such as $2\sigma^+, 3\sigma^+, 1\pi, 4\sigma^{+*}$ for $\mathrm{HF}$, $\mathrm{HCl}$, $\mathrm{HBr}$, and $\mathrm{HI}$).
- **Large Molecule & Conjugated Polyene Support**:
  - Preserves delocalized $\pi/\pi^*$ frontier bands and avoids slicing 2-fold or 3-fold degenerate multiplets across large systems like **$\beta$-Carotene ($\mathrm{C_{40}H_{56}}$, 96 atoms)** and **Adenosine Triphosphate ($\mathrm{C_{10}H_{16}N_5O_{13}P_3}$, 47 atoms)**.
- **Electronic State & Charge Controls**:
  - Adjust molecular **Charge ($q \in [-4, +4]$)** and toggle **Ground State vs. Excited State ($h\nu$, HOMO $\to$ LUMO promotion)** with automatic spin-multiplicity ($2S+1$), magnetic state (Diamagnetic vs. Paramagnetic), and **Bond Order** recalculation.

### 4. Continuous 3D AO $\to$ MO Wavefunction Morphing & Symmetry Stage
- **Real-Time Scalar-Field Wavefunction Morphing (`Morph MO` & `0%–100%` Slider)**:
  - Play/Pause or scrub the interactive **`0%–100%` morph slider** in the 3D viewport header to watch **Separate Atomic Orbitals (`0%`)** continuously overlap, interfere, and fuse into the **Merged Molecular Orbital Isosurface (`100%`)**:
    $$F(\mathbf{r}, t) = (1 - t)\,F_{\text{AO}}(\mathbf{r}) + t\,F_{\text{MO}}(\mathbf{r})$$
- **High-Resolution Anisotropic Grid & Volume-Preserving Taubin $\lambda/\mu$ Surface Fairing**:
  - Uses an adaptive `48³–54³` anisotropic 3D grid tailored to the molecular bounding box (`margin = 2.8 Å`), dynamic peak-normalized isovalues (`32%` of $|\psi|_{\max}$), and two-stage **Taubin $\lambda/\mu$ fairing** ($\lambda = +0.45, \mu = -0.47$) with analytical central-difference $\nabla F(\mathbf{r})$ normals so large polyenes, 6-lobe/8-lobe $4f$ orbitals, and octahedral coordination complexes render with smooth, hole-free, watertight 2-manifold lobes.
- **3D Options & Symmetry Overlays**:
  - Toggle **Nodal Planes**, **Symmetry Axes ($C_n, \sigma_v, \sigma_h, \sigma_d$)**, **Continuous $C_n$ Rotation (`Spin Cn`)**, **Reset Camera**, **Lobe Opacity (`15%–95%`)**, and **Lobe Scale (`0.5×–1.6×`)**.

### 5. Group Theory Character Tables, $\Gamma_{\text{red}}$ Decomposition & NIST Audit
- **Point Group & Character Table Inspector**: Full Schoenflies character tables ($C_1, C_s, C_i, C_2, C_{2v}, C_{3v}, C_{4v}, C_{2h}, D_2, D_{2h}, D_{3h}, D_{4h}, D_{6h}, D_{2d}, D_{3d}, C_{\infty v}, D_{\infty h}, T_d, O_h$) with linear/rotational, quadratic ($d$-orbital), and cubic ($f$-orbital) basis functions.
- **Interactive Reducible Representation ($\Gamma_{\text{red}}$) Decomposition Calculator**:
  - Decompose any reducible representation into its direct sum of irreducible representations ($\Gamma_{\text{red}} = \bigoplus_i a_i \Gamma_i$) using the **Great Orthogonality Theorem**, including built-in presets for $\Gamma_\sigma$, $\Gamma_{\text{val}}$, $\Gamma_{3N}$, $\Gamma_{\text{vib}}$, $\Gamma_\pi$, **$\Gamma_d$ (5 $d$ Orbitals)**, and **$\Gamma_f$ (7 $4f$ Orbitals)**.
- **NIST CCCBDB & PES Audit Tab**: Compares computed geometry parameters ($r_e, \theta_e$), first vertical ionization potentials ($\text{IE}_1 = -\varepsilon_{\text{HOMO}}$), and Koopmans' theorem orbital energies against experimental Photoelectron Spectroscopy (PES) literature references.

---

## Curated Benchmark Molecules (22 NIST PES, CCCBDB, Coordination & Biomolecular Standards)

| Category | Molecule / Complex | Formula | Point Group | Highlighted Electronic & Symmetry Features |
| :--- | :--- | :---: | :---: | :--- |
| **d & f Complexes** | Titanium Tetrachloride | $\mathrm{TiCl_4}$ | $T_d$ | $3d^0$ tetrahedral complex; $1t_1$ Cl lone-pair HOMO and $2e / 3t_2$ $\mathrm{Ti}\;3d$ LUMO splitting |
| **d & f Complexes** | Tetrachloroplatinate(II) | $[\mathrm{PtCl_4}]^{2-}$ | $D_{4h}$ | $5d^8$ square-planar complex; occupied $b_{2g}, e_g, a_{1g}$ $5d$ manifold and $2b_{1g}^*$ ($5d_{x^2-y^2}$) LUMO |
| **d & f Complexes** | Hexachlorochromate(III) | $[\mathrm{CrCl_6}]^{3-}$ | $O_h$ | $3d^3$ octahedral high-spin quartet ($S=3/2$); three singly occupied $2t_{2g}$ SOMOs and $3e_g^*$ LUMO |
| **d & f Complexes** | Hexachlorocerate(III) | $[\mathrm{CeCl_6}]^{3-}$ | $O_h$ | $4f^1$ lanthanide octahedral complex; singly occupied $1a_{2u}$ ($4f_{xyz}$) SOMO and $2t_{1u}^*, 1t_{2u}$ $4f$ LUMOs |
| **Diatomics & Radicals** | Carbon Monoxide | $\mathrm{CO}$ | $C_{\infty v}$ | Heteronuclear polarization; carbon-centered $5\sigma$ HOMO lone pair and $2\pi^*$ LUMO |
| **Diatomics & Radicals** | Nitric Oxide | $\mathrm{NO}$ | $C_{\infty v}$ | 11-valence-electron open-shell radical ($\mathrm{BO} = 2.5$); singly occupied $2\pi^*$ SOMO |
| **Diatomics & Radicals** | Dioxygen | $\mathrm{O_2}$ | $D_{\infty h}$ | Paramagnetic triplet diradical ($^3\Sigma_g^-$); Hund's rule half-filled $1\pi_g^*(x), 1\pi_g^*(y)$ SOMOs |
| **Diatomics & Radicals** | Dinitrogen | $\mathrm{N_2}$ | $D_{\infty h}$ | $s\text{-}p$ mixing inversion ($3\sigma_g$ HOMO above $1\pi_u$ $\pi$-bonding pair) |
| **Diatomics & Radicals** | Carbon Dioxide | $\mathrm{CO_2}$ | $D_{\infty h}$ | 3-center 4-electron $\pi$-system ($1\pi_u$ bonding and $1\pi_g$ non-bonding oxygen HOMO) |
| **Main-Group VSEPR** | Water | $\mathrm{H_2O}$ | $C_{2v}$ | Non-equivalent oxygen lone pairs ($1b_1$ pure $p_x$ HOMO vs. $3a_1$ $sp^2$-hybrid) |
| **Main-Group VSEPR** | Ammonia | $\mathrm{NH_3}$ | $C_{3v}$ | $3a_1$ pyramidal lone pair HOMO and degenerate $1e$ $\mathrm{N\text{--}H}$ bonding pair |
| **Main-Group VSEPR** | Methane | $\mathrm{CH_4}$ | $T_d$ | Triply degenerate $1t_2$ photoemission bands ($2a_1 + 1t_2$ tetrahedral $\sigma$ manifold) |
| **Main-Group VSEPR** | Boron Trifluoride | $\mathrm{BF_3}$ | $D_{3h}$ | Trigonal planar Lewis acid with $1a_2''$ out-of-plane $\pi$ donor and $2a_2''^*$ $2p_z$ LUMO |
| **Hypervalent** | Xenon Hexafluoride | $\mathrm{XeF_6}$ | $C_{3v}$ | Fluxional $7$-pair ($AX_6E$) distorted octahedron with stereochemically active $5s^2$ lone pair |
| **Hypervalent** | Xenon Tetrafluoride | $\mathrm{XeF_4}$ | $D_{4h}$ | Square-planar $AX_4E_2$ hypervalent noble-gas fluoride with axial $5p_z$ lone pair |
| **Hypervalent** | Xenon Difluoride | $\mathrm{XeF_2}$ | $D_{\infty h}$ | Classic Pimentel-Rundle 3-center 4-electron ($3\sigma_u$ non-bonding) linear hypervalency |
| **Hypervalent** | Xenon Oxytetrafluoride | $\mathrm{XeOF_4}$ | $C_{4v}$ | Square pyramidal $AX_5E$ basal $\mathrm{XeF_4}$ 3c-4e bonds vs. axial $\mathrm{Xe{=}O}$ $\pi$-bond |
| **Hypervalent** | Sulfur Tetrafluoride | $\mathrm{SF_4}$ | $C_{2v}$ | Seesaw ($AX_4E$) geometry with equatorial lone pair and axial 3c-4e $\mathrm{F\text{--}S\text{--}F}$ bond |
| **Hypervalent** | Chlorine Trifluoride | $\mathrm{ClF_3}$ | $C_{2v}$ | T-shaped ($AX_3E_2$) interhalogen with two equatorial lone pairs |
| **Hypervalent** | Phosphorus Pentafluoride | $\mathrm{PF_5}$ | $D_{3h}$ | Trigonal bipyramidal axial vs. equatorial bonding and Berry pseudorotation |
| **Hypervalent** | Sulfur Hexafluoride | $\mathrm{SF_6}$ | $O_h$ | Hypervalent 12-electron octahedral bonding ($1t_{1g}$ HOMO, $2e_g^*$ / $3t_{1u}^*$ LUMO) |
| **Conjugated $\pi$** | Ethylene | $\mathrm{C_2H_4}$ | $D_{2h}$ | Full 12-valence-MO $D_{2h}$ manifold and frontier $1b_{3u}(\pi) \to 1b_{2g}(\pi^*)$ gap |
| **Conjugated $\pi$** | Formaldehyde | $\mathrm{CH_2O}$ | $C_{2v}$ | Classic $n \to \pi^*$ ($2b_2 \to 2b_1^*$) carbonyl chromophore |
| **Conjugated $\pi$** | Benzene | $\mathrm{C_6H_6}$ | $D_{6h}$ | Hückel $6\pi$ aromaticity ($1a_{2u}$, degenerate $1e_{1g}$ HOMO & $1e_{2u}^*$ LUMO) |
| **Large Molecules** | $\beta$-Carotene | $\mathrm{C_{40}H_{56}}$ | $D_{2h}$ | 96-atom conjugated 11-double-bond polyene ($22\pi$ electron bathochromic $\pi \to \pi^*$ chromophore) |
| **Large Molecules** | Adenosine Triphosphate | $\mathrm{C_{10}H_{16}N_5O_{13}P_3}$ | $C_s$ | 47-atom bioenergetic coenzyme with adenine $\pi$-ring HOMO and triphosphate anhydride LUMO |

---

## Tech Stack & Verification

- **Frontend**: React 19, TypeScript, Tailwind CSS 4, Lucide Icons, KaTeX mathematical typesetting
- **3D Rendering**: Three.js (`WebGLRenderer`, Marching Tetrahedra scalar-field isosurface generator, volume-preserving Taubin $\lambda/\mu$ surface fairing, analytical $\nabla F$ vertex normals)
- **Chemical Data**: Live PubChem PUG REST API + Curated NIST CCCBDB / PES Benchmark Library
- **Automated Test Suite**: Run `npx tsx --test src/chemistry/moEngine.test.ts` to verify all benchmark eigenvalues, Mulliken symmetries, $s\text{--}p\text{--}d\text{--}f$ Slater-Koster overlaps, coordination bond inference (`[Fe(CN)6]4-`, `Cr(CO)6`, `Fe(CO)5`), and watertight 3D Marching Tetrahedra isosurfaces.
