# Integration of $d$ and $f$ Valence Orbitals into the Quantum MO & Symmetry Engine

Extend **SymmOrbital** from an $s\text{--}p$ valence solver to a full **$s\text{--}p\text{--}d\text{--}f$ Linear Combination of Atomic Orbitals (LCAO-MO)** quantum engine, supporting the 5 real cubic $d$ orbitals ($d_{z^2}, d_{x^2-y^2}, d_{xy}, d_{xz}, d_{yz}$) and 7 real cubic $f$ orbitals ($f_{z^3}, f_{xz^2}, f_{yz^2}, f_{xyz}, f_{z(x^2-y^2)}, f_{x(x^2-3y^2)}, f_{y(3x^2-y^2)}$) across Extended Hückel diagonalization, Slater–Koster two-center overlap integrals, point-group character tables, and 3D Marching Tetrahedra AO-to-MO morphing.

---

## User Review & Critical Decisions

> [!IMPORTANT]
> The following architectural decisions incorporate your confirmed choices from Phase 1 and define the exact physical scope of the $d$- and $f$-orbital expansion:

- **Confirmed Decision 1 — Active Basis Set Scope**:
  - **Transition Metals ($3d, 4d, 5d$)**: Include the full 9-orbital $(n-1)d + ns + np$ valence manifold ($5\times d + 1\times s + 3\times p$).
  - **Hypervalent $p$-Block Elements ($\text{Period } \ge 3$, e.g., $\text{P, S, Cl, As, Se, Br, Te, I, Xe}$)**: Automatically activate outer $nd$ polarization/acceptor orbitals when the atom's coordination number exceeds 3 or total bond order exceeds 4 (such as $\text{SF}_6$, $\text{PF}_5$, $\text{XeF}_4$, $\text{SO}_4^{2-}$, $\text{POCl}_3$), complementing the 3-center 4-electron ($3\text{c--}4\text{e}$) bonding picture with $d_\pi\text{--}p_\pi$ backbonding and $e_g / t_{2g}$ symmetry participation.
  - **Lanthanides ($\text{La--Lu}$, $Z = 57\text{--}71$)**: Include the full 16-orbital $4f + 5d + 6s + 6p$ valence basis ($7\times 4f + 5\times 5d + 1\times 6s + 3\times 6p$).
- **Confirmed Decision 2 — Real Cubic $4f$ Quantum Solver & Hund's Rule**:
  - Model all 7 real cubic $4f$ spherical harmonics explicitly in the Extended Hückel Hamiltonian with analytical two-center Slater–Koster $\sigma, \pi, \delta, \phi$ coupling and crystal-field splitting, populated via Aufbau + Hund's Rule of Maximum Multiplicity across degenerate $f$-shells ($a_{2u}, t_{1u}, t_{2u}$ in $O_h$).
- **Confirmed Decision 3 — Benchmark Library Expansion**:
  - Add a dedicated **"Transition Metal & Lanthanide Complexes"** filter tab to the Benchmark Molecules bar featuring canonical $d$- and $f$-block coordination complexes ($\text{TiCl}_4$, $[\text{PtCl}_4]^{2-}$, $\text{Cr(CO)}_6$, $\text{FeCl}_4^-$, and $\text{CeCl}_6^{2-}$ / $\text{EuCl}_3$).

---

## 1. Overview & Core Concept

- **What It Does**:
  1. Expands the quantum chemistry engine to construct and diagonalize $s, p, d, f$ Extended Hückel Hamiltonians for any transition metal, lanthanide, or hypervalent main-group structure drawn in the 2D Sketcher, imported from PubChem, or selected from the Benchmark bar.
  2. Computes directional two-center Slater–Koster overlap integrals ($ss\sigma, sp\sigma, pp\sigma, pp\pi, sd\sigma, pd\sigma, pd\pi, dd\sigma, dd\pi, dd\delta, sf\sigma, pf\sigma, pf\pi, df\sigma, df\pi, df\delta, ff\sigma$) and atomic Valence State Ionization Energies ($\text{VSIE}_d, \text{VSIE}_f$) for all $d$- and $f$-block elements.
  3. Evaluates 3D spatial wavefunctions and AO-to-MO Marching Tetrahedra morphing isosurfaces for all 5 real $d$ orbitals (cloverleafs and $d_{z^2}$ collar) and all 7 real cubic $f$ orbitals ($f_{z^3}$ double-collar, 8-lobed $f_{xyz}$ / $f_{z(x^2-y^2)}$, and 6-lobed planar $f_{x(x^2-3y^2)}$ / $f_{y(3x^2-y^2)}$).
  4. Upgrades the MO Energy Diagram, Fragment Correlation panels, and Point-Group Character Tables to display ligand-field splittings ($t_{2g}$ vs. $e_g$ in $O_h$; $e$ vs. $t_2$ in $T_d$; $a_{2u}, t_{1u}, t_{2u}$ $4f$-manifolds) and $\delta$ / $\phi$ bonding classifications.
- **Target Audience / Persona**: Chemistry students, educators, and researchers exploring ligand-field theory, $d$-orbital splitting in octahedral/tetrahedral/square-planar complexes, $\pi$-backbonding, and lanthanide $4f$ orbital topologies.
- **Key Value**: Bridges main-group MO theory and coordination chemistry in a single real-time interactive workbench without sacrificing 60 FPS 3D rendering.

---

## 2. User Experience & Visual Design

- **Key User Flows**:
  1. **Benchmark Exploration**: The user clicks the new **"d & f Coordination Complexes"** tab in the Benchmark Molecules bar and selects a complex such as Titanium Tetrachloride ($\text{TiCl}_4$, $T_d$), Hexacarbonylchromium(0) ($\text{Cr(CO)}_6$, $O_h$), Tetrachloroplatinate(II) ($[\text{PtCl}_4]^{2-}$, $D_{4h}$), or Hexachlorocerate(IV) ($[\text{CeCl}_6]^{2-}$, $O_h$).
  2. **Ligand-Field MO Diagram Inspection**: The central MO Energy Diagram displays the metal $(n-1)d$, $ns$, $np$, and (for lanthanides) $4f$ fragment levels on the left, the ligand $\sigma/\pi$ SALCs on the right, and the resulting ligand-field split molecular orbitals ($t_{2g}, e_g^*, t_{1u}, a_{2u}$) in the center with spin arrows obeying Hund's rule.
  3. **3D $d$/$f$ Isosurface & Morphing Playback**: Selecting any $d$- or $f$-derived MO renders its multi-lobed 3D isosurface (4-lobe cloverleafs, toroidal $d_{z^2}$ / $f_{z^3}$ collars, or 6-/8-lobe cubic $f$ manifolds). Clicking **Morph MO** smoothly animates the separate metal $d/f$ and ligand $s/p$ atomic orbitals fusing into the delocalized coordination MO.
  4. **Custom 2D Drawing & PubChem Import**: Selecting any transition metal ($\text{Sc--Zn, Y--Cd, Hf--Hg}$) or lanthanide ($\text{La--Lu}$) in the Full 2D Sketcher and accepting the structure automatically invokes the $d/f$-enabled VSEPR 3D geometry reconstructor and Extended Hückel solver.
- **Visual Identity & Theme**:
  - Preserves the scientific workbench aesthetic (`Plus Jakarta Sans` UI headers, `JetBrains Mono` / `KaTeX` mathematical formulas, 60-30-10 neutral slate canvas with semantic blue/red or teal/amber phase lobes).
  - Adds clean LaTeX rendering for all $d$ and $f$ basis labels ($3d_{z^2}, 3d_{x^2-y^2}, 4f_{z^3}, 4f_{xyz}, 4f_{z(x^2-y^2)}, 4f_{x(x^2-3y^2)}$) and $\delta / \delta^*$ interaction badges.
- **Interactive Feedback & Motion**:
  - Adaptive Marching Tetrahedra grid resolution automatically increases angular sampling around high-order $d$ and $f$ nodal cones so 6-lobe and 8-lobe $f$ orbitals remain crisp, symmetric, and 100% watertight during 60 FPS morph scrubbing.

---

## 3. Key Product Decisions & Trade-Offs

- **Decision 1 — Real Cubic $f$-Orbital Basis Set vs. General Spherical Set**:
  - *Chosen Approach*: Use the **7 real cubic $f$ orbitals** ($f_{z^3}, f_{xz^2}, f_{yz^2}, f_{xyz}, f_{z(x^2-y^2)}, f_{x(x^2-3y^2)}, f_{y(3x^2-y^2)}$).
  - *Why*: In octahedral ($O_h$), tetrahedral ($T_d$), and cubic ligand fields, the cubic $f$ set partitions cleanly into irreducible representations ($a_{2u} \rightarrow f_{xyz}$; $t_{1u} \rightarrow f_{x^3}, f_{y^3}, f_{z^3}$; $t_{2u} \rightarrow f_{z(x^2-y^2)}, f_{x(y^2-z^2)}, f_{y(z^2-x^2)}$), making character-table highlighting and symmetry assignment exact.
  - *Alternatives Considered*: Complex spherical harmonics $Y_3^m(\theta, \phi)$, which cannot be rendered as real signed $+\psi / -\psi$ lobes without taking real/imaginary linear combinations.
- **Decision 2 — Directional Cosine Tensor Evaluation for $d$ and $f$ Slater–Koster Overlaps**:
  - *Chosen Approach*: Evaluate two-center overlaps between arbitrary $s, p, d, f$ basis functions by projecting both orbitals onto a local diatomic frame $(\hat{e}_\sigma, \hat{e}_{\pi 1}, \hat{e}_{\pi 2})$ aligned with the internuclear vector $\vec{R}_{AB}$ via numerical Gauss–Legendre / analytical spherical-harmonic rotation.
  - *Why*: Writing out 16×16 analytical polynomial tables in global $(x,y,z)$ coordinates is error-prone under arbitrary 3D rotations; local-frame projection guarantees exact rotational invariance ($SO(3)$ symmetry) so degenerate $t_{2g}$, $e_g$, $t_{1u}$, and $t_{2u}$ levels never split artificially when a complex is rotated.
- **Decision 3 — Adaptive Radial Exponent & Grid Resolution for Multi-Nodal $d/f$ Surfaces**:
  - *Chosen Approach*: Scale the angular polynomial prefactors by $(r^l e^{-\zeta_l r})$ with $l=2$ for $d$ and $l=3$ for $f$, and bump the Marching Tetrahedra grid density slightly when the selected MO has dominant $f$-orbital weight ($\sum c_f^2 > 0.25$).
  - *Why*: An $f$ orbital has 3 angular nodal surfaces (up to 8 lobes around a single nucleus). Slightly finer voxel sampling around the metal center prevents adjacent opposite-phase lobes from merging or pinching while keeping keyframe generation fast.

---

## 4. Technical Architecture & Data Strategy *(Technical Reference)*

### Architecture & Data Flow Diagram

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                     SymmOrbital s-p-d-f Quantum Pipeline                     │
└──────────────────────────────────────────────────────────────────────────────┘
                                       │
       ┌───────────────────────────────┼───────────────────────────────┐
       ▼                               ▼                               ▼
┌────────────────────────────┐ ┌────────────────────────────┐ ┌────────────────────────────┐
│  1. Element & Basis Store  │ │ 2. 3D VSEPR & Coordination │ │ 3. Point-Group & Character │
│  - vsieS, vsieP, vsieD,    │ │    Geometry Reconstructor  │ │    Tables (C2v..Oh, Dinfh) │
│    vsieF (eV) for Z=1..86  │ │  - Octahedral (Oh, ML6)    │ │  - Cubic d & f basis funcs │
│  - Active orbital selector │ │  - Square Planar (D4h,ML4) │ │  - Delta (δ) & Phi (φ)     │
│    (TM: d+s+p, Ln: f+d+s+p │ │  - Tetrahedral (Td, ML4)   │ │    irreducible reps        │
│     Hypervalent p: s+p+d)  │ │  - Trigonal Bipyramid(D3h) │ │  - Γ_d and Γ_f presets     │
└──────────────┬─────────────┘ └──────────────┬─────────────┘ └──────────────┬─────────────┘
               │                              │                              │
               └──────────────────────────────┼──────────────────────────────┘
                                              ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│           4. Extended Hückel s-p-d-f Solver & Local-Frame Overlaps           │
│  - Constructs N×N Hamiltonian H_ii = VSIE(n,l), H_ij = ½ K(H_ii+H_jj) S_ij   │
│  - SO(3)-invariant local internuclear frame (ê_σ, ê_π1, ê_π2) overlap S_ij   │
│  - Classifies σ, σ*, π, π*, δ, δ*, n (d/f non-bonding) & Mulliken irreps     │
│  - Applies Aufbau + Hund's Rule across degenerate d (t2g/eg) & f shells      │
└──────────────────────────────────────┬───────────────────────────────────────┘
                                       │
                  ┌────────────────────┴────────────────────┐
                  ▼                                         ▼
┌─────────────────────────────────────────┐ ┌─────────────────────────────────────────┐
│     5. MO Diagram & Inspector UI        │ │   6. 3D Wavefunction & Isosurface Stage │
│  - Metal (n-1)d / 4f vs Ligand SALCs    │ │  - Real cubic l=2 (5 d) & l=3 (7 f)     │
│  - LaTeX formatting for d & f subshells │ │    angular harmonics in ψ(x,y,z)        │
│  - Benchmark Tab: "d & f Complexes"     │ │  - Watertight Marching Tetrahedra morph │
│    (TiCl4, [PtCl4]2-, Cr(CO)6, [CeCl6]) │ │  - Analytic lobe builder for 'lobes'    │
└─────────────────────────────────────────┘ └─────────────────────────────────────────┘
```

### Data Model & Basis Specification

1. **Expanded `OrbitalBasisType` Union**:
   - **$s$ & $p$ (4 orbitals)**: `'s' | 'px' | 'py' | 'pz'`
   - **$d$ (5 real cubic orbitals, $l = 2$)**:
     - `'dz2'` ($\propto 2z^2 - x^2 - y^2$ or principal symmetry axis)
     - `'dx2y2'` ($\propto x^2 - y^2$)
     - `'dxy'` ($\propto xy$), `'dxz'` ($\propto xz$), `'dyz'` ($\propto yz$)
   - **$f$ (7 real cubic orbitals, $l = 3$)**:
     - `'fz3'` ($\propto z(2z^2 - 3x^2 - 3y^2)$)
     - `'fxz2'` ($\propto x(4z^2 - x^2 - y^2)$)
     - `'fyz2'` ($\propto y(4z^2 - x^2 - y^2)$)
     - `'fxyz'` ($\propto xyz$)
     - `'fz_x2y2'` ($\propto z(x^2 - y^2)$)
     - `'fx_x2_3y2'` ($\propto x(x^2 - 3y^2)$)
     - `'fy_3x2_y2'` ($\propto y(3x^2 - y^2)$)
2. **Element Parameterization (`ELEMENT_PROPERTIES`)**:
   - Add `vsieD?: number` and `vsieF?: number` (in eV) and `block: 's' | 'p' | 'd' | 'f'` across all 86 elements:
     - $3d$ series ($\text{Sc}\ -8.5\text{ eV}$ to $\text{Zn}\ -14.8\text{ eV}$), $4d$ series ($\text{Y}\ -8.2\text{ eV}$ to $\text{Cd}\ -13.9\text{ eV}$), $5d$ series ($\text{Hf}\ -8.8\text{ eV}$ to $\text{Hg}\ -14.2\text{ eV}$).
     - $4f$ lanthanide series ($\text{La}\ -11.2\text{ eV}$ to $\text{Lu}\ -16.8\text{ eV}$, with $5d$ at $-7.5\text{ to }-8.8\text{ eV}$ and $6s$ at $-6.0\text{ to }-6.8\text{ eV}$).
     - Hypervalent Period 3–6 $p$-block $nd$ polarization energies ($\text{P, S, Cl, As, Se, Br, Te, I, Xe}$ at $-4.5\text{ to }-6.5\text{ eV}$).
3. **Expanded `MolecularOrbital['interactionType']`**:
   - Add `'δ'` and `'δ*'` alongside `'σ' | 'π' | 'σ*' | 'π*' | 'n (Lone Pair)'` so metal–metal and metal–ligand $\delta$-symmetry orbitals are explicitly badged.

### Interactive Component & State Mapping

- **Coordination VSEPR & $d^n$ Geometry Rules**:
  - Update `reconstruct3DGeometryFromVSEPR` and `detectPointGroup` so transition-metal and lanthanide complexes adopt their experimental coordination geometries:
    - $\text{ML}_6$ ($d^0\text{--}d^{10}, f^n$, e.g., $\text{Cr(CO)}_6, \text{TiF}_6^{2-}, \text{CeCl}_6^{2-}$) $\rightarrow$ Octahedral ($O_h$).
    - $\text{ML}_4$ with $d^8$ metal ($\text{Pt, Pd, Ni, Au}$, e.g., $[\text{PtCl}_4]^{2-}$) $\rightarrow$ Square Planar ($D_{4h}$); other $\text{ML}_4$ ($\text{TiCl}_4, \text{FeCl}_4^-$) $\rightarrow$ Tetrahedral ($T_d$).
- **3D Wavefunction, Marching Tetrahedra & Analytic Lobe Rendering**:
  - Update `evaluateOrbitalWavefunction` and `moIsosurface.ts` (`evalSingleBasis`) with normalized real cubic $l=2$ and $l=3$ harmonic polynomials so both `Merged MO (Iso)` morphing and `Point Cloud` modes render all $5\times d$ and $7\times f$ orbitals accurately.
  - Update `OrbitalVisualizer3D.tsx` (`mode === 'lobes'`) to construct multi-lobed analytic geometries for all 5 $d$ orbitals (4-lobe cloverleafs in $xy, xz, yz, x^2-y^2$ plus $d_{z^2}$ dumbbell+torus) and all 7 $f$ orbitals (axial dumbbell+double-torus for $f_{z^3}$, 6-lobe planar stars for $f_{x(x^2-3y^2)} / f_{y(3x^2-y^2)}$, and 8-lobe octants for $f_{xyz} / f_{z(x^2-y^2)}$).
- **LaTeX & Character Table Integration**:
  - Update `latexUtils.ts` so $d$ and $f$ orbital labels format cleanly into mathematical subscripts/superscripts (e.g., `3dx2y2` $\rightarrow 3d_{x^2-y^2}$, `4fz3` $\rightarrow 4f_{z^3}$, `4fx_x2_3y2` $\rightarrow 4f_{x(x^2-3y^2)}$).
  - Enrich `symmetryData.ts` character tables ($O_h, T_d, D_{4h}, D_{3h}, C_{4v}, C_{3v}, C_{2v}$) with cubic $f$-orbital basis assignments and $\Gamma_d / \Gamma_f$ reducible representation presets.
