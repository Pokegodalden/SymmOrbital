# SymmOrbital

> **Interactive quantum chemistry workbench for 2D chemical sketching, real-time LCAO-MO energy diagrams, 3D WebGL orbital isosurfaces, and Schoenflies point-group symmetry analysis—benchmarked against NIST CCCBDB & PubChem.**

![TypeScript](https://img.shields.io/badge/TypeScript-5.8-3178C6?logo=typescript&logoColor=white)
![React](https://img.shields.io/badge/React-19-61DAFB?logo=react&logoColor=0F172A)
![Three.js](https://img.shields.io/badge/Three.js-WebGL-000000?logo=threedotjs&logoColor=white)
![KaTeX](https://img.shields.io/badge/KaTeX-LaTeX_Math-008080?logo=latex&logoColor=white)
![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)

**SymmOrbital** bridges the gap between static textbook Molecular Orbital (MO) diagrams and research-grade quantum chemistry visualizers. Sketch any molecule in 2D—or query **PubChem** by name, formula, SMILES, or CID—and immediately inspect its 3D $\psi_i(\mathbf{r})$ wavefunction phase lobes, Symmetry-Adapted Linear Combinations (SALCs), Mulliken character tables, and experimental photoelectron spectroscopy benchmarks.

---

## Key Features

### 1. Synchronized 3-Pane Scientific Workbench
- **Pane 1 — Interactive 2D Chemical Sketcher**: Draw single, double, triple, wedge, and dashed stereochemical bonds, place rings ($\mathrm{C_6}$, $\mathrm{C_5}$, $\mathrm{C_4}$, $\mathrm{C_3}$), or use rapid keyboard shortcuts (`1`/`2`/`3` for bond order, `C`/`N`/`O`/`H`/`F`/`S` for elements, `Del` to erase).
- **Pane 2 — 3D WebGL Orbital Isosurface & Symmetry Stage**: Real-time Three.js visualization of positive ($\psi > 0$) and negative ($\psi < 0$) wavefunction phase lobes, nodal surfaces, principal rotation axes ($C_n$), mirror planes ($\sigma_v, \sigma_h, \sigma_d$), and animated point-group symmetry operations.
- **Pane 3 — Interactive LCAO-MO Correlation Diagram**: Full Fragment AO $\leftrightarrow$ Molecular Orbital $\leftrightarrow$ SALC energy ladder with HOMO/SOMO–LUMO gap ($\Delta E$) telemetry, ground vs. excited state ($h\nu$) electron promotion, and complete KaTeX-rendered LCAO expansions.

### 2. Automated Point-Group & Mulliken Character Table Engine
- Detects Schoenflies point groups ($C_s, C_{2v}, C_{3v}, D_{2h}, D_{3h}, D_{6h}, D_{\infty h}, C_{\infty v}, T_d, O_h$) directly from 3D coordinates.
- Interactive character tables where hovering or clicking any Mulliken irreducible representation ($A_1, B_1, E, T_2, \Pi, \Sigma^+$, etc.) highlights all matching MOs across the correlation diagram.

### 3. Live PubChem Integration & NIST CCCBDB Audit Suite
- **Unified PubChem Omnibar**: Resolve live 3D coordinates, IUPAC names, SMILES, and molecular metadata directly from the NCBI PubChem PUG-REST API.
- **Experimental Database Audit**: Built-in verification suite comparing computed Extended Hückel / Koopmans' theorem orbital energies ($\mathrm{IP}_i \approx -\varepsilon_i$) against experimental **NIST CCCBDB** and ** Turner He(I) Photoelectron Spectroscopy (PES)** vertical ionization potentials, complete with Mean Absolute Error (MAE) and Pearson $r^2$ correlation metrics.

### 4. Chemical Edge-Case & Open-Shell Radical Handling
- **Subvalent & Multiple-Bond Valences**: Accurately preserves unsaturated and dative bonding topologies in $\mathrm{CO}$ ($C\equiv O$), $\mathrm{NO}$, $\mathrm{NO_2}$, $\mathrm{O_3}$, and $\mathrm{SO_2}$ without erroneous implicit-hydrogen padding.
- **Hund's Rule & Open-Shell Radicals**: Supports degenerate $\pi / e / t_2$ manifold half-filling and odd-electron SOMO radicals (such as Nitric Oxide, $\mathrm{NO}$, with a $2\pi^*$ SOMO and bond order $2.5$).
- **s–p Mixing Correction**: Captures $2\sigma_g / 1\pi_u$ and $5\sigma / 1\pi$ orbital ordering inversions in heteronuclear and homonuclear diatomics.

---

## Curated Benchmark Molecules

| Molecule | Formula | Point Group | Frontier Highlight |
| :--- | :--- | :--- | :--- |
| **Water** | $\mathrm{H_2O}$ | $C_{2v}$ | Non-bonding oxygen $1b_1$ ($2p_x$) lone pair HOMO |
| **Ammonia** | $\mathrm{NH_3}$ | $C_{3v}$ | $3a_1$ pyramidal lone pair & degenerate $1e$ ($\mathrm{N\text{--}H}$) bonding |
| **Methane** | $\mathrm{CH_4}$ | $T_d$ | Triply degenerate $1t_2$ HOMO — PES proof against naive $sp^3$ equivalence |
| **Benzene** | $\mathrm{C_6H_6}$ | $D_{6h}$ | Frost circle $\pi$-manifold ($1a_{2u}, 1e_{1g}, 1e_{2u}, 1b_{2g}$) |
| **Ethylene** | $\mathrm{C_2H_4}$ | $D_{2h}$ | Classic $\pi(1b_{3u}) \to \pi^*(1b_{2g})$ frontier transition |
| **Formaldehyde** | $\mathrm{CH_2O}$ | $C_{2v}$ | Oxygen $2b_2$ ($n_O$) HOMO and $2b_1$ ($\pi^*_{\mathrm{CO}}$) LUMO |
| **Carbon Dioxide** | $\mathrm{CO_2}$ | $D_{\infty h}$ | Doubly degenerate $1\pi_g$ non-bonding terminal oxygen HOMO |
| **Carbon Monoxide** | $\mathrm{CO}$ | $C_{\infty v}$ | Carbon-centered $5\sigma$ HOMO via $s\text{--}p$ mixing |
| **Nitric Oxide** | $\mathrm{NO}$ | $C_{\infty v}$ | Open-shell $11e^-$ radical with singly occupied $2\pi^*$ SOMO ($\mathrm{BO}=2.5$) |
| **Boron Trifluoride** | $\mathrm{BF_3}$ | $D_{3h}$ | Vacant boron $2a_2''$ ($2p_z$) Lewis-acid acceptor LUMO |
| **Sulfur Hexafluoride** | $\mathrm{SF_6}$ | $O_h$ | Octahedral ligand-field $1t_{1g}$ HOMO & $2e_g^*$ antibonding LUMO |

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

SymmOrbital includes an automated test suite verifying point-group classification, implicit-hydrogen valence rules, open-shell Hund's rule occupancy, and experimental HOMO/SOMO symmetries:

```bash
npx tsx --test src/chemistry/moEngine.test.ts
```

### Production Build

```bash
npm run build
```

---

## Architecture & Directory Structure

```text
src/
├── chemistry/
│   ├── moEngine.ts            # Extended Hückel LCAO-MO solver, 3D geometry & point-group classifier
│   ├── moEngine.test.ts       # Unit tests for valences, diatomics, radicals & Hund's rule
│   ├── moleculesData.ts       # 11 canonical benchmark molecules with calibrated AO/MO/SALC sets
│   ├── symmetryData.ts        # Schoenflies point groups, symmetry elements & character tables
│   ├── databaseBenchmarks.ts  # Experimental NIST CCCBDB & He(I) PES vertical ionization potentials
│   ├── pubchemService.ts      # Live NCBI PubChem PUG-REST 3D/2D structure & property resolver
│   └── latexUtils.ts          # Chemical formula, Schoenflies, Mulliken irrep & LCAO LaTeX formatter
├── components/
│   ├── ChemicalSketcher2D.tsx # Interactive 2D vector sketcher with SMILES & ring/bond tools
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
