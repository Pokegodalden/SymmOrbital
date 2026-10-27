[Simplified Usage Guide](./USAGE_GUIDE.md)

# SymmOrbital

> **Interactive quantum chemistry workbench for 2D chemical sketching, VSEPR-aware 3D geometry reconstruction, real-time LCAO-MO energy diagrams, 3D WebGL orbital isosurfaces, and Schoenflies point-group symmetry analysis—benchmarked against NIST CCCBDB & PubChem.**

![TypeScript](https://img.shields.io/badge/TypeScript-5.8-3178C6?logo=typescript&logoColor=white)
![React](https://img.shields.io/badge/React-19-61DAFB?logo=react&logoColor=0F172A)
![Three.js](https://img.shields.io/badge/Three.js-WebGL-000000?logo=threedotjs&logoColor=white)
![KaTeX](https://img.shields.io/badge/KaTeX-LaTeX_Math-008080?logo=latex&logoColor=white)
![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)

**SymmOrbital** bridges the gap between static textbook Molecular Orbital (MO) diagrams and research-grade quantum chemistry visualizers. Sketch any molecule in 2D—or query **PubChem** by name, formula, SMILES, or CID—and immediately inspect its 3D $\psi_i(\mathbf{r})$ wavefunction phase lobes, Symmetry-Adapted Linear Combinations (SALCs), Mulliken character tables, and experimental photoelectron spectroscopy benchmarks.

---

## Key Features

### 1. Synchronized 3-Pane Scientific Workbench (with Dark Mode & Reversible Layout)
- **Pane 1 — Interactive 2D Chemical Sketcher**: Draw single (`—`), double (`═`), triple (`≡`), stereo wedge (`▲`), and **dashed (`┄`)** bonds, stamp rings ($\mathrm{C_6}$, $\mathrm{C_4}$, $\mathrm{C_3}$), or use rapid keyboard shortcuts (`1`/`2`/`3` for bond order, `D`/`4` for dashed bond, `W` for wedge bond, `C`/`N`/`O`/`H`/`F`/`S`/`P`/`Cl`/`Xe` for elements, `Del` to erase).
- **Pane 2 — Interactive LCAO-MO Correlation Diagram**: Full Fragment AO $\leftrightarrow$ Molecular Orbital $\leftrightarrow$ SALC energy ladder with HOMO/SOMO–LUMO gap ($\Delta E$) telemetry, ionic charge state ($q \in [-2, +2]$), ground vs. excited state ($h\nu$) electron promotion, and complete KaTeX-rendered LCAO expansions.
- **Pane 3 — 3D WebGL Orbital Isosurface & Symmetry Stage**: Real-time Three.js visualization of positive ($\psi > 0$) and negative ($\psi < 0$) wavefunction phase lobes, nodal surfaces, principal rotation axes ($C_n$), mirror planes ($\sigma_v, \sigma_h, \sigma_d$), and animated point-group symmetry operations.
- **Dark Mode & Layout Toggle**: Switch seamlessly between Light and Dark themes via the header sun/moon toggle, or toggle between the equal-width **Pro 3-Pane** workbench and the **Classic v1** stacked view.

### 2. VSEPR 3D Coordinate Reconstructor & Lone-Pair Point-Group Engine
- **Automatic 2D $\to$ 3D VSEPR Coordinate Reconstruction (`reconstruct3DGeometryFromVSEPR`)**: Converts flat 2D sketches and PubChem `2D-SDF` fallbacks (such as hypervalent fluorides where PubChem lacks MMFF94s 3D conformers) into accurate 3D stereochemical coordinates across steric numbers $2\text{--}7$.
- **Lone-Pair-Aware $AX_nE_m$ Symmetry Solver (`getCentralAtomLonePairs`)**: Computes central-atom non-bonding lone pairs $E = \lfloor\max(0, V - \sum b_i)/2\rfloor$ to distinguish true 3D molecular point groups from flat 2D projections:
  - $AX_6E_1$ ($\mathrm{XeF_6}$): Fluxional monocapped octahedron ($C_{3v}$) instead of $O_h$ or planar $D_{6h}$
  - $AX_4E_2$ ($\mathrm{XeF_4}$): Square planar ($D_{4h}$)
  - $AX_4E_1$ ($\mathrm{SF_4}$): Seesaw ($C_{2v}$)
  - $AX_3E_2$ ($\mathrm{ClF_3}$): T-shaped ($C_{2v}$)
  - $AX_5E_1$ ($\mathrm{BrF_5}$ / $\mathrm{IF_5}$): Square pyramidal ($C_{4v}$)
  - $AX_5E_0$ ($\mathrm{PF_5}$): Trigonal bipyramidal ($D_{3h}$)
- **Interactive Character Tables**: Supports $C_s, C_{2v}, C_{3v}, C_{4v}, D_{2h}, D_{3h}, D_{4h}, D_{6h}, D_{\infty h}, C_{\infty v}, T_d, O_h$. Hovering or clicking any Mulliken irreducible representation ($A_1, B_1, E, T_2, \Pi, \Sigma^+$, etc.) highlights all matching MOs across the correlation diagram.

### 3. Live PubChem Integration & NIST CCCBDB Audit Suite
- **Unified PubChem Omnibar**: Resolve live 3D coordinates, IUPAC names, SMILES, and molecular metadata directly from the NCBI PubChem PUG-REST API, complete with 3D geometry provenance indicators (`PubChem MMFF94s 3D-SDF Conformer` vs. `2D-SDF → 3D VSEPR Reconstructed`).
- **Experimental Database Audit**: Built-in verification suite comparing computed Extended Hückel / Koopmans' theorem orbital energies ($\mathrm{IP}_i \approx -\varepsilon_i$) against experimental **NIST CCCBDB** and **He(I) Photoelectron Spectroscopy (PES)** vertical ionization potentials across all 16 benchmark molecules, complete with Mean Absolute Error (MAE) and Pearson $r^2$ correlation metrics.

### 4. Chemical Edge-Case & Open-Shell Radical Handling
- **Subvalent & Multiple-Bond Valences**: Accurately preserves unsaturated and dative bonding topologies in $\mathrm{CO}$ ($C\equiv O$), $\mathrm{NO}$, $\mathrm{NO_2}$, $\mathrm{O_3}$, and $\mathrm{SO_2}$ without erroneous implicit-hydrogen padding.
- **Hund's Rule & Open-Shell Radicals**: Supports degenerate $\pi / e / t_2$ manifold half-filling and odd-electron SOMO radicals (such as Nitric Oxide, $\mathrm{NO}$, with a $2\pi^*$ SOMO and bond order $2.5$).
- **s–p Mixing Correction**: Captures $2\sigma_g / 1\pi_u$ and $5\sigma / 1\pi$ orbital ordering inversions in heteronuclear and homonuclear diatomics.

---

## Curated Benchmark Molecules (16 NIST PES Standards)

| Molecule | Formula | Point Group | VSEPR / Topology | Frontier Highlight |
| :--- | :--- | :--- | :--- | :--- |
| **Water** | $\mathrm{H_2O}$ | $C_{2v}$ | $AX_2E_2$ (Bent) | Non-bonding oxygen $1b_1$ ($2p_x$) lone pair HOMO |
| **Ammonia** | $\mathrm{NH_3}$ | $C_{3v}$ | $AX_3E_1$ (Trigonal Pyramidal) | $3a_1$ pyramidal lone pair & degenerate $1e$ ($\mathrm{N\text{--}H}$) bonding |
| **Methane** | $\mathrm{CH_4}$ | $T_d$ | $AX_4E_0$ (Tetrahedral) | Triply degenerate $1t_2$ HOMO — PES proof against naive $sp^3$ equivalence |
| **Benzene** | $\mathrm{C_6H_6}$ | $D_{6h}$ | Planar $\pi$-Ring | Frost circle $\pi$-manifold ($1a_{2u}, 1e_{1g}, 1e_{2u}, 1b_{2g}$) |
| **Ethylene** | $\mathrm{C_2H_4}$ | $D_{2h}$ | Planar $\pi$-Alkene | Classic $\pi(1b_{3u}) \to \pi^*(1b_{2g})$ frontier transition |
| **Formaldehyde** | $\mathrm{CH_2O}$ | $C_{2v}$ | Planar Carbonyl | Oxygen $2b_2$ ($n_O$) HOMO and $2b_1$ ($\pi^*_{\mathrm{CO}}$) LUMO |
| **Carbon Dioxide** | $\mathrm{CO_2}$ | $D_{\infty h}$ | $AX_2E_0$ (Linear) | Doubly degenerate $1\pi_g$ non-bonding terminal oxygen HOMO |
| **Carbon Monoxide** | $\mathrm{CO}$ | $C_{\infty v}$ | Heteronuclear Diatomic | Carbon-centered $5\sigma$ HOMO via $s\text{--}p$ mixing |
| **Nitric Oxide** | $\mathrm{NO}$ | $C_{\infty v}$ | Open-Shell Radical | $11e^-$ radical with singly occupied $2\pi^*$ SOMO ($\mathrm{BO}=2.5$) |
| **Boron Trifluoride** | $\mathrm{BF_3}$ | $D_{3h}$ | $AX_3E_0$ (Trigonal Planar) | Vacant boron $2a_2''$ ($2p_z$) Lewis-acid acceptor LUMO |
| **Sulfur Hexafluoride** | $\mathrm{SF_6}$ | $O_h$ | $AX_6E_0$ (Octahedral) | Octahedral ligand-field $1t_{1g}$ HOMO & $2e_g^*$ antibonding LUMO |
| **Xenon Hexafluoride** | $\mathrm{XeF_6}$ | $C_{3v}$ | $AX_6E_1$ (Distorted Octahedral) | Stereochemically active $5s^2$ lone pair in $4a_1^*$ HOMO & $4e^*$ LUMO |
| **Xenon Tetrafluoride** | $\mathrm{XeF_4}$ | $D_{4h}$ | $AX_4E_2$ (Square Planar) | Axial $5p_z$ lone-pair $2a_{1g}$ HOMO & $3e_u^*$ 3c–4e LUMO |
| **Sulfur Tetrafluoride** | $\mathrm{SF_4}$ | $C_{2v}$ | $AX_4E_1$ (Seesaw) | Equatorial $3a_1$ lone-pair HOMO & axial $2b_1^*$ 3c–4e LUMO |
| **Chlorine Trifluoride** | $\mathrm{ClF_3}$ | $C_{2v}$ | $AX_3E_2$ (T-Shaped) | Equatorial lone pair $3a_1$ HOMO & hypervalent axial $2b_1^*$ LUMO |
| **Phosphorus Pentafluoride** | $\mathrm{PF_5}$ | $D_{3h}$ | $AX_5E_0$ (Trigonal Bipyramidal) | Non-bonding fluorine $1e''$ HOMO & $3a_1'^*$ antibonding LUMO |

---

## Quick Start

### Prerequisites
- **Node.js** `>= 18.x`
- **npm** `>= 9.x`

### Installation & Local Development

```bash
# Install dependencies
npm install

# Start the development server on http://localhost:3000
npm run dev
```

### Running the Quantum Chemistry & Edge-Case Test Suite

SymmOrbital includes an automated test suite verifying point-group classification, VSEPR 3D coordinate reconstruction, implicit-hydrogen valence rules, open-shell Hund's rule occupancy, and experimental HOMO/SOMO symmetries across all 16 benchmark molecules:

```bash
npx tsx --test src/chemistry/moEngine.test.ts
```

### Production Build & GitHub Pages Deployment

```bash
npm run build
```

This repository includes an automated GitHub Actions workflow (`.github/workflows/deploy-pages.yml`) and relative asset base configuration (`vite.config.ts`) for zero-config **GitHub Pages** hosting. Push to `main` and set **Settings → Pages → Source** to **GitHub Actions**.

---

## Architecture & Directory Structure

```text
src/
├── chemistry/
│   ├── moEngine.ts            # Extended Hückel LCAO-MO solver, 3D VSEPR reconstructor & point-group classifier
│   ├── moEngine.test.ts       # Unit tests for valences, VSEPR special cases, radicals & Hund's rule
│   ├── moleculesData.ts       # 16 canonical benchmark molecules with calibrated AO/MO/SALC sets
│   ├── symmetryData.ts        # Schoenflies point groups (including C4v/D4h/Oh), symmetry elements & character tables
│   ├── databaseBenchmarks.ts  # Experimental NIST CCCBDB & He(I) PES vertical ionization potentials
│   ├── pubchemService.ts      # Live NCBI PubChem PUG-REST 3D/2D structure & property resolver
│   └── latexUtils.ts          # Chemical formula, Schoenflies, Mulliken irrep & LCAO LaTeX formatter
├── components/
│   ├── ChemicalSketcher2D.tsx # Interactive 2D vector sketcher with single/double/triple/wedge/dash & ring tools
│   ├── OrbitalVisualizer3D.tsx# Three.js WebGL 3D molecular & wavefunction phase isosurface stage
│   ├── MODiagramPanel.tsx     # Interactive SVG AO-MO-SALC correlation ladder & level table
│   ├── SymmetryGroupPanel.tsx # Interactive Mulliken character table & symmetry operation trigger
│   ├── DatabaseVerificationPanel.tsx # NIST CCCBDB vs. Extended Hückel Koopmans' audit table
│   └── Latex.tsx              # Fast, accessible KaTeX math & chemical renderer
├── backup_v1/
│   └── AppClassicV1.tsx       # Reversible V1 Classic stacked layout snapshot
└── App.tsx                    # V2 Pro 3-Pane Synchronized Workbench with Light/Dark theme support
```

---

## License

Licensed under the **Apache License 2.0**.
