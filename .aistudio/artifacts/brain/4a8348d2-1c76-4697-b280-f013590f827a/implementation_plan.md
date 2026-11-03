# 3D Atomic-to-Molecular Orbital (AO $\rightarrow$ MO) Scalar-Field Morphing Animation

This feature introduces a real-time 3D **Scalar-Field Marching Cubes Morphing Engine** into the 3D Orbital & Symmetry Stage, enabling users to animate or scrub the continuous quantum-mechanical interference transition from isolated **Atomic Orbitals ($t = 0\%$)** to the fully fused **Molecular Orbital Isosurface ($t = 100\%$)** at 60 FPS.

---

## User Review & Critical Decisions

> [!IMPORTANT]
> The following architectural and UX decisions incorporate your confirmed selections from Phase 1:

- **Confirmed Decision 1 — Scalar-Field Marching Cubes Morphing**: The 3D stage will evaluate a continuous parameterized 3D scalar field $F(\mathbf{r}, t)$ that interpolates between non-interfering atom-centered basis orbitals at $t = 0$ and the full constructive/destructive LCAO superposition $\psi_{\text{MO}}(\mathbf{r}) = \sum_i c_i \phi_i(\mathbf{r} - \mathbf{R}_i)$ at $t = 1$, extracting watertight $+\psi$ and $-\psi$ triangle meshes with smooth surface normals.
- **Confirmed Decision 2 — Play/Pause Button + Interactive 0–100% Morph Scrubber**: The 3D stage will provide an accessible **Play / Pause** toggle paired with a live **0% (Separate AOs) $\leftrightarrow$ 100% (Merged MO)** range slider so users can either play the smooth fusion animation or manually scrub through the orbital necking and nodal pinch-off stages.
- **Confirmed Decision 3 — Pre-Cached Keyframes with Adaptive Grid Resolution**: To guarantee zero frame drops on both compact molecules ($\mathrm{H_2O}, \mathrm{C_2H_4}, \mathrm{C_6H_6}$) and large macromolecules ($\beta\text{-Carotene}, \mathrm{ATP}$), the engine will pre-compute and cache keyframe isosurface geometries using an adaptive voxel grid resolution scaled to the molecular footprint.

---

## 1. Overview & Core Concept

- **What It Does**: Bridges the pedagogical gap between individual atomic basis functions ($s, p_x, p_y, p_z, d_{z^2}$) and delocalized molecular orbitals ($\sigma, \sigma^*, \pi, \pi^*$). At $t = 0\%$, lobes sit strictly centered on their parent nuclei without cross-atom overlap. As $t$ advances toward $100\%$, in-phase lobes ($\text{sign}(c_i \phi_i) = \text{sign}(c_j \phi_j)$) swell toward the inter-nuclear region and fuse into a continuous bonding envelope, while out-of-phase lobes ($\text{sign}(c_i \phi_i) \neq \text{sign}(c_j \phi_j)$) repel and carve out a sharp nodal plane ($\psi = 0$).
- **Target Audience / Persona**: Chemistry students, educators, and researchers inspecting how LCAO coefficients and point-group symmetry govern covalent bond formation, $\pi$-conjugation, and antibonding nodes.
- **Key Value**: Lets users directly observe **constructive interference** (neck fusion in bonding MOs) and **destructive interference** (nodal cleavage in antibonding MOs) rather than viewing only static, overlapping atom-centered teardrops.

---

## 2. User Experience & Visual Design

- **Key User Flows**:
  1. **Quick Morph Playback (Header & Bottom Dock)**: The 3D stage header and collapsible bottom control dock include a **Play / Pause Morph** button and a **Morph ($0\% \rightarrow 100\%$)** slider with tabular percentage readout (`tabular-nums`).
  2. **Automatic Mode Activation**: Clicking **Play** or dragging the **Morph** slider automatically switches the 3D render style into the **Merged MO (Morph)** surface mode (or users can select `Merged MO` directly in the render mode selector alongside `LCAO Lobes`, `ψ Grid Cloud`, and `Mesh`).
  3. **Interactive Scrubbing**: Dragging the slider between `0%` (`Separate AOs`) and `100%` (`Merged MO`) immediately updates the 3D isosurface mesh in the viewport while preserving camera orientation, lighting, symmetry axes/planes, and nodal plane overlays.
  4. **Seamless Loop / Completion**: Clicking **Play** smoothly animates $t$ from $0\% \rightarrow 100\%$ over $2.4\text{ s}$ with an ease-in-out curve (holding briefly at $100\%$ so the user can inspect the merged MO). If clicked again while already at $100\%$, it replays from $0\% \rightarrow 100\%$.
- **Visual Identity & Theme**:
  - *Aesthetic Direction*: High-precision scientific instrument viewport matching the existing studio PBR aesthetic.
  - *Color Palette & Mood*: Retains the dual-phase $+\psi$ / $-\psi$ physical material system (`#2563EB` Royal Blue / `#DC2626` Crimson Rose in standard mode; `#0D9488` Teal / `#D97706` Amber in colorblind-safe mode) with `MeshPhysicalMaterial` clearcoat and studio `RoomEnvironment` reflections.
  - *Typography & Hierarchy*: Monospace tabular numerals (`font-mono tabular-nums`) for morph percentage (`0%`–`100%`) and lobe scale readouts; clean single-line button labels (`whitespace-nowrap shrink-0`).
  - *Viewport Ergonomics*: Compact inline Play/Pause + Scrubber control in the 3D stage header bar for immediate one-click access without opening the drawer, mirrored inside the collapsible bottom control dock so the central 3D canvas remains 100% unobstructed.

---

## 3. Key Product Decisions & Trade-Offs

- **Decision 1: Interference-Parameterized Scalar Field $F(\mathbf{r}, t)$**
  - *Chosen Approach*: For any point $\mathbf{r}$ and morph parameter $t \in [0, 1]$, let $\psi_i(\mathbf{r}) = c_i\,\phi_i(\mathbf{r} - \mathbf{R}_i)$ be each atomic orbital contribution. Define the dominant single-atom field $\psi_{\text{dom}}(\mathbf{r}) = \arg\max_{\psi_i} |\psi_i(\mathbf{r})|$ (which has zero cross-atom constructive fusion and slightly contracted radial decay so lobes remain distinct at $t = 0$) and the full LCAO superposition $\psi_{\text{MO}}(\mathbf{r}) = \sum_i \psi_i(\mathbf{r})$. The morphing field blends:
    $$F(\mathbf{r}, t) = (1 - t)\,\psi_{\text{dom}}(\mathbf{r}) + t\,\psi_{\text{MO}}(\mathbf{r})$$
  - *Why*: At $t = 0$, each voxel is governed strictly by its nearest dominant atomic orbital, producing clean, separate atom-centered $s/p/d$ lobes. As $t \rightarrow 1$, the cross-terms $\sum_{j \neq i} \psi_j(\mathbf{r})$ turn on continuously, causing in-phase lobes to bridge across the bond axis and out-of-phase lobes to pinch off along nodal surfaces.
  - *Alternatives Considered*: Translating atoms physically apart at $t = 0$ was rejected because moving nuclei alters the molecular geometry and misaligns point-group symmetry elements.
- **Decision 2: Adaptive Pre-Cached Keyframe Ring Buffer**
  - *Chosen Approach*: Pre-extract and cache $K$ discrete Marching Cubes `BufferGeometry` pairs ($+\psi$ and $-\psi$ meshes across $t \in [0, 1]$) on demand for the active `(moleculeId, moId)` pair, using an adaptive grid resolution ($28^3$ for small/medium molecules, $22^3$ for large molecules like $\beta\text{-Carotene}$ and $\mathrm{ATP}$) with smooth vertex normal interpolation.
  - *Why*: Eliminates per-frame CPU voxel evaluation during playback and scrubbing, guaranteeing 60 FPS on hardware GPUs and smooth frame-paced rendering on software WebGL without memory leaks.

---

## 4. Technical Architecture & Data Strategy

### Architecture & Component Diagram

```
┌─────────────────────────────────────────────────────────────────────────┐
│                        OrbitalVisualizer3D.tsx                          │
│                                                                         │
│  ┌───────────────────────────────────────────────────────────────────┐  │
│  │ Top Stage Header & Collapsible Bottom Dock                        │  │
│  │  • Render Mode: [LCAO Lobes] [Merged MO] [ψ Cloud] [Mesh]         │  │
│  │  • Morph Controls: [▶ Play / ⏸ Pause]  [0% ════●════ 100%]        │  │
│  └───────────────────────────────┬───────────────────────────────────┘  │
│                                  │ morphProgress (0.0 → 1.0)            │
│                                  ▼                                      │
│  ┌───────────────────────────────────────────────────────────────────┐  │
│  │ AO-to-MO Scalar-Field Polygonizer & Keyframe Cache                │  │
│  │  1. Adaptive Bounding Grid (28³ standard / 22³ macromolecules)    │  │
│  │  2. Dual-Field Evaluator:                                         │  │
│  │     • ψ_dom(r): Isolated dominant AO term at t = 0                │  │
│  │     • ψ_MO(r):  Full LCAO superposition Σ c_i φ_i(r) at t = 1     │  │
│  │  3. Marching Cubes Isosurface Extractor (+ψ_iso & -ψ_iso meshes)  │  │
│  │  4. Keyframe Geometry Cache (disposed cleanly on MO/Mol change)   │  │
│  └───────────────────────────────┬───────────────────────────────────┘  │
│                                  │ THREE.BufferGeometry (+ψ, -ψ)        │
│                                  ▼                                      │
│  ┌───────────────────────────────────────────────────────────────────┐  │
│  │ Three.js Studio PBR Scene (contentGroup)                          │  │
│  │  • Atoms & Bonds (CPK + Metallic Collars)                         │  │
│  │  • Morphing Orbital Surface (+ψ Blue/Teal, -ψ Rose/Amber)         │  │
│  │  • Nodal Planes & 3D Point-Group Symmetry Elements                │  │
│  └───────────────────────────────────────────────────────────────────┘  │
└─────────────────────────────────────────────────────────────────────────┘
```

### Interactive Component & State Mapping

1. **Render Mode & Morph State (`OrbitalVisualizer3D.tsx`)**:
   - Adds `'merged'` (`Merged MO`) to `renderMode` (`'lobes' | 'merged' | 'cloud' | 'wireframe'`).
   - Adds `morphProgress` state (`0` to `1`, defaulting to `1.0` when selecting `Merged MO`, or animating `0.0 -> 1.0` when `isMorphPlaying` is triggered).
   - Adds `isMorphPlaying` boolean state driven by a `requestAnimationFrame` (and `afterDraw` software-WebGL safe) animation loop that advances `morphProgress` smoothly and stops cleanly at `1.0` or when paused.
2. **Scalar-Field Marching Cubes Polygonizer (`src/chemistry/moIsosurface.ts`)**:
   - Pure TypeScript module exporting `buildMorphIsosurfaceKeyframe(mo, atoms, morphT, options)` and a lightweight LRU keyframe cache for the active orbital.
   - Evaluates active LCAO basis terms on a 3D bounding grid around the active orbital's participating atoms (Automatically tightens bounding box around active LCAO centers so localized orbitals in ATP and delocalized polyene chains in $\beta\text{-Carotene}$ both achieve optimal voxel density).
   - Extracts smooth indexed triangle meshes (`positions`, `normals`) for $+\psi_{\text{iso}}$ and $-\psi_{\text{iso}}$ using standard 256-case Marching Cubes edge/triangle lookup tables and linear edge interpolation.
3. **Memory & Lifecycle Hygiene**:
   - Cached keyframe `BufferGeometry` instances are reused across slider scrubbing and play loops for the active MO and disposed whenever the molecule or selected MO changes or the component unmounts.
