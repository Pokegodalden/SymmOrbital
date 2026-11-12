import test from 'node:test';
import assert from 'node:assert/strict';
import { CANONICAL_MOLECULES } from './moleculesData.ts';
import {
  POINT_GROUPS,
  decomposeReducibleRepresentation,
} from './symmetryData.ts';
import {
  analyzeMolecularGraph,
  computeMolecularFormula,
  detectPointGroup,
  evaluateOrbitalWavefunction,
  jacobiEigenvalues,
  parseChemicalInputToGraph,
  reconstruct3DGeometryFromVSEPR,
} from './moEngine.ts';
import {
  formulaToLatex,
  irrepToLatex,
  lcaoExpansionToLatex,
  moLabelToLatex,
  schoenfliesToLatex,
} from './latexUtils.ts';
import {
  NIST_CCCBDB_BENCHMARKS,
  evaluateEngineAgainstBenchmarks,
} from './databaseBenchmarks.ts';
import {
  buildMorphIsosurfaceKeyframe,
  evaluateMorphScalarFieldAtPoint,
} from './moIsosurface.ts';

test('LaTeX formatting utilities convert chemical formulas, Schoenflies groups, irreps, and MO labels accurately', () => {
  assert.equal(formulaToLatex('H2O'), '\\mathrm{H_{2}O}');
  assert.equal(formulaToLatex('C6H6'), '\\mathrm{C_{6}H_{6}}');
  assert.equal(schoenfliesToLatex('C2v'), 'C_{2v}');
  assert.equal(schoenfliesToLatex('D∞h'), 'D_{\\infty h}');
  assert.equal(schoenfliesToLatex('Td'), 'T_d');
  assert.equal(irrepToLatex('B1'), 'B_{1}');
  assert.equal(irrepToLatex('E1g'), 'E_{1g}');
  assert.equal(irrepToLatex("A2''"), "A_{2}^{''}");
  assert.equal(irrepToLatex('Σg+'), '\\Sigma_g^+');
  assert.equal(irrepToLatex('b1'), 'b_{1}');
  assert.equal(irrepToLatex('e1g'), 'e_{1g}');
  assert.equal(irrepToLatex("a2''"), "a_{2}^{''}");
  assert.equal(irrepToLatex('σg+'), '\\sigma_g^+');
  assert.equal(moLabelToLatex('1b1'), '1b_{1}');
  assert.equal(moLabelToLatex('4a1*'), '4a_{1}^{*}');
  assert.equal(moLabelToLatex('1e1g (π2)'), '1e_{1g}\\;(\\pi_{2})');
  assert.ok(moLabelToLatex("21a'' (n_O, γ-PO3)").includes('\\gamma '));

  const water = CANONICAL_MOLECULES.find((m) => m.formula === 'H2O')!;
  const homo = water.molecularOrbitals.find((m) => m.isHOMO)!;
  const lcaoTex = lcaoExpansionToLatex(homo.label, homo.lcao);
  assert.ok(lcaoTex.includes('\\psi\\!\\left(1b_{1}\\right)'), `Unexpected LCAO LaTeX: ${lcaoTex}`);
});

test('NIST Chemistry WebBook PES & NIST CCCBDB benchmark audit achieves 100% point-group parity and < 0.05 eV calibrated MAE across all 22 benchmarks', () => {
  const water = CANONICAL_MOLECULES[0];
  const audit = evaluateEngineAgainstBenchmarks(water);

  assert.equal(audit.rows.length, 22);
  assert.equal(audit.pointGroupAccuracyPercent, 100);
  assert.ok(
    audit.meanAbsoluteErrorCalibratedEv < 0.05,
    `Expected calibrated MAE < 0.05 eV vs NIST WebBook PES, got ${audit.meanAbsoluteErrorCalibratedEv}`
  );
  assert.ok(
    audit.meanGeometryErrorAngstrom < 0.025,
    `Expected geometry MAE < 0.025 Å vs NIST CCCBDB, got ${audit.meanGeometryErrorAngstrom}`
  );

  // Verify Benzene (C6H6) HOMO matches NIST WebBook PES (9.24 eV)
  const benzeneRow = audit.rows.find((r) => r.formula === 'C6H6')!;
  assert.equal(-benzeneRow.calibratedHomoEnergyEv, 9.24);
  assert.equal(benzeneRow.nistWebBookIpEv, 9.24);
  assert.equal(NIST_CCCBDB_BENCHMARKS.C6H6.experimentalPointGroup, 'D6h');

  // Verify Nitric Oxide (NO) SOMO matches NIST WebBook PES (9.26 eV) and C∞v symmetry
  const noRow = audit.rows.find((r) => r.formula === 'NO')!;
  assert.equal(-noRow.calibratedHomoEnergyEv, 9.26);
  assert.equal(noRow.nistWebBookIpEv, 9.26);
  assert.equal(NIST_CCCBDB_BENCHMARKS.NO.experimentalPointGroup, 'C∞v');

  // Verify Xenon Hexafluoride (XeF6) matches C3v monocapped octahedron and NIST PES IP1 = 12.19 eV
  const xef6Row = audit.rows.find((r) => r.formula === 'XeF6')!;
  assert.equal(xef6Row.detectedPointGroup, 'C3v');
  assert.equal(-xef6Row.calibratedHomoEnergyEv, 12.19);

  // Verify Xenon Difluoride (XeF2) matches D∞h 3c–4e benchmark (r = 1.977 Å, BO = 1.0, IP1 = 12.35 eV)
  const xef2Row = audit.rows.find((r) => r.formula === 'XeF2')!;
  assert.equal(xef2Row.detectedPointGroup, 'D∞h');
  assert.equal(-xef2Row.calibratedHomoEnergyEv, 12.35);
  assert.equal(xef2Row.computedBondLength, 1.977);
  const xef2Mol = CANONICAL_MOLECULES.find((m) => m.formula === 'XeF2')!;
  assert.equal(xef2Mol.bondOrder, 1.0);

  // Verify Xenon Oxytetrafluoride (XeOF4) matches C4v heteroleptic AX4YE1 benchmark (r(Xe=O) = 1.703 Å, IP1 = 12.90 eV)
  const xeof4Row = audit.rows.find((r) => r.formula === 'XeOF4')!;
  assert.equal(xeof4Row.detectedPointGroup, 'C4v');
  assert.equal(-xeof4Row.calibratedHomoEnergyEv, 12.90);
  assert.equal(xeof4Row.computedBondLength, 1.703);

  // Verify Dioxygen (O2) triplet diradical (1πg* SOMO pair at -12.07 eV, BO = 2.0)
  const o2Mol = analyzeMolecularGraph(
    CANONICAL_MOLECULES.find((m) => m.formula === 'O2')!.atoms,
    CANONICAL_MOLECULES.find((m) => m.formula === 'O2')!.bonds
  );
  const o2Somos = o2Mol.molecularOrbitals.filter((m) => m.degenerateGroup === 'o2_1pig');
  assert.equal(o2Somos.length, 2);
  assert.equal(o2Somos[0].occupancy, 1);
  assert.equal(o2Somos[1].occupancy, 1);
  assert.equal(o2Mol.bondOrder, 2.0);

  // Verify Large Molecules: β-Carotene (C40H56, D2h, IP1 = 6.50 eV, 22-center delocalized π HOMO)
  // and Adenosine Triphosphate (ATP, C10H16N5O13P3, Cs, IP1 = 8.40 eV)
  const betaCaroteneRow = audit.rows.find((r) => r.formula === 'C40H56')!;
  assert.equal(betaCaroteneRow.detectedPointGroup, 'D2h');
  assert.equal(-betaCaroteneRow.calibratedHomoEnergyEv, 6.50);
  const betaCaroteneMol = CANONICAL_MOLECULES.find((m) => m.formula === 'C40H56')!;
  const bcHomo = betaCaroteneMol.molecularOrbitals.find((m) => m.isHOMO)!;
  assert.equal(bcHomo.lcao.length, 22);

  const atpRow = audit.rows.find((r) => r.formula === 'C10H16N5O13P3')!;
  assert.equal(atpRow.detectedPointGroup, 'Cs');
  assert.equal(-atpRow.calibratedHomoEnergyEv, 8.40);

  // Verify Ethylene (C2H4, D2h) has its complete 12-MO valence manifold (6 occupied = 12e-, 6 virtual)
  // matching Γ_val = 3Ag + B2g + 2B3g + 3B1u + 2B2u + B3u and all 6 experimental NIST PES bands
  const c2h4Mol = CANONICAL_MOLECULES.find((m) => m.formula === 'C2H4')!;
  assert.equal(c2h4Mol.molecularOrbitals.length, 12);
  const c2h4Occupied = c2h4Mol.molecularOrbitals.filter((m) => m.occupancy === 2);
  const c2h4Virtual = c2h4Mol.molecularOrbitals.filter((m) => m.occupancy === 0);
  assert.equal(c2h4Occupied.length, 6);
  assert.equal(c2h4Virtual.length, 6);
  const c2h4TotalElec = c2h4Mol.molecularOrbitals.reduce((s, m) => s + m.occupancy, 0);
  assert.equal(c2h4TotalElec, 12);
  const d2hDecomp = decomposeReducibleRepresentation(
    POINT_GROUPS.D2h,
    [12, 0, 0, 0, 0, 0, 4, 8]
  );
  assert.equal(d2hDecomp.formula, '3ag + b2g + 2b3g + 3b1u + 2b2u + b3u');
  assert.equal(NIST_CCCBDB_BENCHMARKS.C2H4.pesBands.length, 6);
});

test('Point group detection, VSEPR 3D reconstruction, and edge-case implicit-H rules preserve all 22 canonical molecules including C40H56 and C10H16N5O13P3', () => {
  const expected: Record<string, string> = {
    H2O: 'C2v',
    NH3: 'C3v',
    C6H6: 'D6h',
    CH4: 'Td',
    C2H4: 'D2h',
    CH2O: 'C2v',
    CO2: 'Dinfh',
    BF3: 'D3h',
    SF6: 'Oh',
    CO: 'Cinfv',
    NO: 'Cinfv',
    XeF6: 'C3v',
    XeF4: 'D4h',
    SF4: 'C2v',
    ClF3: 'C2v',
    PF5: 'D3h',
    XeF2: 'Dinfh',
    O2: 'Dinfh',
    N2: 'Dinfh',
    XeOF4: 'C4v',
    C40H56: 'D2h',
    C10H16N5O13P3: 'Cs',
    TiCl4: 'Td',
    PtCl4: 'D4h',
    CrCl6: 'Oh',
    CeCl6: 'Oh',
    C2H4Cl3Pt: 'C2v',
  };

  for (const mol of CANONICAL_MOLECULES) {
    const detected = detectPointGroup(mol.atoms, mol.bonds);
    assert.equal(
      detected,
      expected[mol.formula],
      `Expected ${mol.formula} to be classified as ${expected[mol.formula]}, got ${detected}`
    );

    // Ensure analyzeMolecularGraph and computeMolecularFormula(..., true) never mutate CO -> CHO or NO -> HNO
    const analyzed = analyzeMolecularGraph(mol.atoms, mol.bonds);
    assert.equal(
      analyzed.formula,
      mol.formula,
      `Expected analyzeMolecularGraph(${mol.formula}) to have formula ${mol.formula}, got ${analyzed.formula}`
    );
    const sketchFormula = computeMolecularFormula(mol.atoms, mol.bonds, true);
    assert.equal(
      sketchFormula.formula,
      mol.formula,
      `Expected computeMolecularFormula(${mol.formula}, true) to be ${mol.formula}, got ${sketchFormula.formula}`
    );
  }

  // Test Nitric Oxide (NO) drawn with a single, double, or triple bond in sketcher:
  // all must resolve to NO (11 valence e-, BO = 2.5, 2pi* SOMO occupancy = 1), never HNO or H2NO
  for (const order of [1, 2, 3] as const) {
    const noDrawn = analyzeMolecularGraph(
      [
        { id: 'N1', element: 'N', x2d: -0.58, y2d: 0, x3d: -0.58, y3d: 0, z3d: 0 },
        { id: 'O1', element: 'O', x2d: 0.57, y2d: 0, x3d: 0.57, y3d: 0, z3d: 0 },
      ],
      [{ id: 'b1', source: 'N1', target: 'O1', order }]
    );
    assert.equal(noDrawn.formula, 'NO');
    assert.equal(noDrawn.valenceElectrons, 11);
    assert.equal(noDrawn.bondOrder, 2.5);
    const somo = noDrawn.molecularOrbitals.find((m) => m.isHOMO)!;
    assert.equal(somo.occupancy, 1);
  }

  // Test Hund's rule on reduced NO- (12 valence e-, isoelectronic with O2):
  // both degenerate 2pi*(x) and 2pi*(y) orbitals must be singly occupied (1, 1), not (2, 0)
  const noAnion = analyzeMolecularGraph(
    [
      { id: 'N1', element: 'N', x2d: -0.58, y2d: 0, x3d: -0.58, y3d: 0, z3d: 0 },
      { id: 'O1', element: 'O', x2d: 0.57, y2d: 0, x3d: 0.57, y3d: 0, z3d: 0 },
    ],
    [{ id: 'b1', source: 'N1', target: 'O1', order: 2 }],
    { charge: -1 }
  );
  const piStarShell = noAnion.molecularOrbitals.filter((m) => m.degenerateGroup === 'no_2pi');
  assert.equal(piStarShell.length, 2);
  assert.equal(piStarShell[0].occupancy, 1);
  assert.equal(piStarShell[1].occupancy, 1);

  // Bent homonuclear triatomic Ozone (O3, angle 116.8 deg, open chain) must be C2v and stay O3 (never H2O3)
  const ozoneAtoms = [
    { id: 'O1', element: 'O' as const, x2d: 0, y2d: 0.4, x3d: 0, y3d: 0.4, z3d: 0 },
    { id: 'O2', element: 'O' as const, x2d: -1.08, y2d: -0.25, x3d: -1.08, y3d: -0.25, z3d: 0 },
    { id: 'O3', element: 'O' as const, x2d: 1.08, y2d: -0.25, x3d: 1.08, y3d: -0.25, z3d: 0 },
  ];
  const ozoneBonds = [
    { id: 'b1', source: 'O1', target: 'O2', order: 2 as const },
    { id: 'b2', source: 'O1', target: 'O3', order: 1 as const },
  ];
  assert.equal(detectPointGroup(ozoneAtoms, ozoneBonds), 'C2v');
  const ozoneAnalyzed = analyzeMolecularGraph(ozoneAtoms, ozoneBonds);
  assert.equal(ozoneAnalyzed.formula, 'O3');
  assert.equal(ozoneAnalyzed.valenceElectrons, 18);

  // Nitrogen Dioxide (NO2) radical (17 valence e-) must remain NO2 (never HNO2) AND be classified as bent C2v (never Dinfh)
  const no2Atoms = [
    { id: 'N1', element: 'N' as const, x2d: 0, y2d: 0, x3d: 0, y3d: 0, z3d: 0 },
    { id: 'O1', element: 'O' as const, x2d: -1.2, y2d: 0, x3d: -1.2, y3d: 0, z3d: 0 },
    { id: 'O2', element: 'O' as const, x2d: 1.2, y2d: 0, x3d: 1.2, y3d: 0, z3d: 0 },
  ];
  const no2Bonds = [
    { id: 'b1', source: 'N1', target: 'O1', order: 2 as const },
    { id: 'b2', source: 'N1', target: 'O2', order: 1 as const },
  ];
  assert.equal(detectPointGroup(no2Atoms, no2Bonds), 'C2v');
  const no2Analyzed = analyzeMolecularGraph(no2Atoms, no2Bonds);
  assert.equal(no2Analyzed.formula, 'NO2');
  assert.equal(no2Analyzed.pointGroupId, 'C2v');
  assert.equal(no2Analyzed.valenceElectrons, 17);

  // Xenon Oxytetrafluoride (XeOF4, AX5E1 -> C4v Square Pyramidal) must place unique Oxygen atom on axial +y axis even when drawn last in 2D
  const xeof4Atoms = [
    { id: 'Xe1', element: 'Xe' as const, x2d: 0, y2d: 0, x3d: 0, y3d: 0, z3d: 0 },
    { id: 'F1', element: 'F' as const, x2d: -1.3, y2d: -0.5, x3d: -1.3, y3d: 0.5, z3d: 0 },
    { id: 'F2', element: 'F' as const, x2d: 1.3, y2d: -0.5, x3d: 1.3, y3d: 0.5, z3d: 0 },
    { id: 'F3', element: 'F' as const, x2d: -1.3, y2d: 0.5, x3d: -1.3, y3d: -0.5, z3d: 0 },
    { id: 'F4', element: 'F' as const, x2d: 1.3, y2d: 0.5, x3d: 1.3, y3d: -0.5, z3d: 0 },
    { id: 'O1', element: 'O' as const, x2d: 0, y2d: -1.4, x3d: 0, y3d: 1.4, z3d: 0 },
  ];
  const xeof4Bonds = [
    { id: 'b1', source: 'Xe1', target: 'F1', order: 1 as const },
    { id: 'b2', source: 'Xe1', target: 'F2', order: 1 as const },
    { id: 'b3', source: 'Xe1', target: 'F3', order: 1 as const },
    { id: 'b4', source: 'Xe1', target: 'F4', order: 1 as const },
    { id: 'b5', source: 'Xe1', target: 'O1', order: 2 as const },
  ];
  const reconstructedXeOF4 = reconstruct3DGeometryFromVSEPR(xeof4Atoms, xeof4Bonds);
  const oReconstructed = reconstructedXeOF4.find((a) => a.id === 'O1')!;
  assert.equal(oReconstructed.x3d, 0);
  assert.equal(oReconstructed.z3d, 0);
  assert.ok(oReconstructed.y3d > 1.5, `Expected VSEPR-reconstructed axial O1 on +y axis, got y3d=${oReconstructed.y3d}`);

  const xeof4Analyzed = analyzeMolecularGraph(xeof4Atoms, xeof4Bonds);
  assert.equal(xeof4Analyzed.pointGroupId, 'C4v');
  const oAtom3d = xeof4Analyzed.atoms.find((a) => a.id === 'O1')!;
  assert.equal(oAtom3d.x3d, 0);
  assert.equal(oAtom3d.z3d, 0);
  assert.ok(oAtom3d.y3d > 1.5, `Expected axial O1 on +y axis, got y3d=${oAtom3d.y3d}`);

  // Sulfur Dioxide (SO2, 18 valence e-) must remain SO2 (never H2SO2)
  const so2Atoms = [
    { id: 'S1', element: 'S' as const, x2d: 0, y2d: 0.4, x3d: 0, y3d: 0.4, z3d: 0 },
    { id: 'O1', element: 'O' as const, x2d: -1.2, y2d: -0.25, x3d: -1.2, y3d: -0.25, z3d: 0 },
    { id: 'O2', element: 'O' as const, x2d: 1.2, y2d: -0.25, x3d: 1.2, y3d: -0.25, z3d: 0 },
  ];
  const so2Bonds = [
    { id: 'b1', source: 'S1', target: 'O1', order: 2 as const },
    { id: 'b2', source: 'S1', target: 'O2', order: 2 as const },
  ];
  const so2Analyzed = analyzeMolecularGraph(so2Atoms, so2Bonds);
  assert.equal(so2Analyzed.formula, 'SO2');
  assert.equal(so2Analyzed.valenceElectrons, 18);
});

test('Great Orthogonality Theorem reducible representation decomposition works for C2v and D6h', () => {
  const c2v = POINT_GROUPS.C2v;
  // 6 valence basis orbitals of H2O: characters [6, 0, 2, 4] -> 3a1 + b1 + 2b2
  const decomp = decomposeReducibleRepresentation(c2v, [6, 0, 2, 4]);
  assert.equal(decomp.formula, '3a1 + b1 + 2b2');

  const d6h = POINT_GROUPS.D6h;
  const benzenePiPreset = d6h.defaultReduciblePresets[0];
  const piDecomp = decomposeReducibleRepresentation(d6h, benzenePiPreset.characters);
  assert.ok(piDecomp.formula.includes('e1g') && piDecomp.formula.includes('a2u'));
});

test('Slater-Koster 4-orbital Extended Hückel solver and Jacobi diagonalizer conserve trace and valence occupancy', () => {
  const H = [
    [-11.4, -2.4, 0, -2.4],
    [-2.4, -11.4, -2.4, 0],
    [0, -2.4, -11.4, -2.4],
    [-2.4, 0, -2.4, -11.4],
  ];
  const { eigenvalues } = jacobiEigenvalues(H);
  const sumEig = eigenvalues.reduce((a, b) => a + b, 0);
  const trace = H[0][0] + H[1][1] + H[2][2] + H[3][3];
  assert.ok(Math.abs(sumEig - trace) < 1e-4);

  // Custom 1,3-Butadiene graph
  const butadieneAtoms = [
    { id: 'C1', element: 'C' as const, x2d: -1.8, y2d: 0.3, x3d: -1.8, y3d: -0.3, z3d: 0 },
    { id: 'C2', element: 'C' as const, x2d: -0.6, y2d: -0.3, x3d: -0.6, y3d: 0.3, z3d: 0 },
    { id: 'C3', element: 'C' as const, x2d: 0.6, y2d: 0.3, x3d: 0.6, y3d: -0.3, z3d: 0 },
    { id: 'C4', element: 'C' as const, x2d: 1.8, y2d: -0.3, x3d: 1.8, y3d: 0.3, z3d: 0 },
  ];
  const butadieneBonds = [
    { id: 'b1', source: 'C1', target: 'C2', order: 2 as const },
    { id: 'b2', source: 'C2', target: 'C3', order: 1 as const },
    { id: 'b3', source: 'C3', target: 'C4', order: 2 as const },
  ];
  const mol = analyzeMolecularGraph(butadieneAtoms, butadieneBonds, {
    customName: '1,3-Butadiene',
  });
  assert.equal(mol.formula, 'C4H6');
  const homo = mol.molecularOrbitals.find((m) => m.isHOMO);
  const lumo = mol.molecularOrbitals.find((m) => m.isLUMO);
  assert.ok(homo && lumo && lumo.energy > homo.energy);
});

test('3D LCAO wavefunction evaluator reflects nodal symmetry of H2O 1b1 HOMO across z=0', () => {
  const water = CANONICAL_MOLECULES.find((m) => m.formula === 'H2O')!;
  const homo1b1 = water.molecularOrbitals.find((m) => m.id === 'mo_1b1')!;
  const psiPosZ = evaluateOrbitalWavefunction(homo1b1, water.atoms, 0.0, 0.35, 0.8);
  const psiNegZ = evaluateOrbitalWavefunction(homo1b1, water.atoms, 0.0, 0.35, -0.8);
  const psiNode = evaluateOrbitalWavefunction(homo1b1, water.atoms, 0.0, 0.35, 0.0);

  assert.ok(psiPosZ > 0.05);
  assert.ok(psiNegZ < -0.05);
  assert.ok(Math.abs(psiNode) < 1e-6);
});

test('Scalar-field Marching Cubes AO-to-MO morphing separates AOs at t=0 and fuses constructive bonding lobes at t=1', () => {
  const ethylene = CANONICAL_MOLECULES.find((m) => m.formula === 'C2H4')!;
  const piHomo = ethylene.molecularOrbitals.find((m) => m.isHOMO)!;

  // At the C=C π-cloud midpoint above the molecular plane (x = 0, y = 0, z = 0.55),
  // t = 0 suppresses cross-atom overlap whereas t = 1 constructively sums both C1(2pz) and C2(2pz) lobes
  const midSeparate = evaluateMorphScalarFieldAtPoint(
    piHomo,
    ethylene.atoms,
    0,
    0,
    0.55,
    0.0
  );
  const midMerged = evaluateMorphScalarFieldAtPoint(
    piHomo,
    ethylene.atoms,
    0,
    0,
    0.55,
    1.0
  );
  assert.ok(
    midMerged > midSeparate * 2.5,
    `Expected merged MO midpoint amplitude (${midMerged}) > 2.5x separate AO midpoint (${midSeparate})`
  );

  // Extract 0% (Separate AOs) and 100% (Merged MO) isosurface meshes
  const meshSep = buildMorphIsosurfaceKeyframe(piHomo, ethylene.atoms, 0.0, 1.0);
  const meshMerged = buildMorphIsosurfaceKeyframe(piHomo, ethylene.atoms, 1.0, 1.0);

  assert.ok(meshSep.positive.positions.length > 0 && meshSep.negative.positions.length > 0);
  assert.ok(meshMerged.positive.positions.length > 0 && meshMerged.negative.positions.length > 0);
  assert.equal(meshMerged.positive.positions.length, meshMerged.positive.normals.length);

  // Check unit length of first normal vector
  const nx = meshMerged.positive.normals[0];
  const ny = meshMerged.positive.normals[1];
  const nz = meshMerged.positive.normals[2];
  assert.ok(Math.abs(Math.hypot(nx, ny, nz) - 1.0) < 1e-3);

  // Verify every extracted positive-phase vertex lies strictly on crossed edges near F(p, 1.0) ≈ 0.20
  const posVerts = meshMerged.positive.positions;
  for (let i = 0; i < posVerts.length; i += 3) {
    const val = evaluateMorphScalarFieldAtPoint(
      piHomo,
      ethylene.atoms,
      posVerts[i],
      posVerts[i + 1],
      posVerts[i + 2],
      1.0
    );
    assert.ok(
      Math.abs(val - 0.20) < 0.06,
      `Expected vertex ${i / 3} to lie on isosurface F≈0.20, got ${val}`
    );
  }

  // Verify 100% watertightness (0 open boundary edges -> zero side clipping) across morph states and high lobeScale
  const countOpenBoundaryEdges = (indices: Uint32Array): number => {
    const edgeCounts = new Map<string, number>();
    for (let i = 0; i < indices.length; i += 3) {
      const tri = [indices[i], indices[i + 1], indices[i + 2]];
      for (let e = 0; e < 3; e++) {
        const a = tri[e];
        const b = tri[(e + 1) % 3];
        const key = a < b ? `${a}_${b}` : `${b}_${a}`;
        edgeCounts.set(key, (edgeCounts.get(key) || 0) + 1);
      }
    }
    let openEdges = 0;
    for (const count of edgeCounts.values()) {
      if (count === 1) openEdges++;
    }
    return openEdges;
  };

  assert.equal(
    countOpenBoundaryEdges(meshSep.positive.indices),
    0,
    'Separate AO positive mesh must be 100% watertight with 0 clipped boundary edges'
  );
  assert.equal(
    countOpenBoundaryEdges(meshMerged.positive.indices),
    0,
    'Merged MO positive mesh must be 100% watertight with 0 clipped boundary edges'
  );
  assert.equal(
    countOpenBoundaryEdges(meshMerged.negative.indices),
    0,
    'Merged MO negative mesh must be 100% watertight with 0 clipped boundary edges'
  );

  // Also verify at maximum user Lobe Scale (1.6x) that no side clipping occurs
  const meshMaxScale = buildMorphIsosurfaceKeyframe(piHomo, ethylene.atoms, 1.0, 1.6);
  assert.equal(
    countOpenBoundaryEdges(meshMaxScale.positive.indices),
    0,
    'Max-scale (1.6x) Merged MO mesh must have 0 clipped boundary edges'
  );
});

test('HF and HCl (including tilted 3D PubChem SDF geometries) produce literature-standard Mulliken labels (2σ+, 3σ+, 1π, 4σ+*) and accurate 3D bonding/antibonding topologies', () => {
  const hf = analyzeMolecularGraph(
    [
      { id: 'F1', element: 'F', x2d: -0.5, y2d: 0, x3d: 0, y3d: 0, z3d: 0 },
      { id: 'H1', element: 'H', x2d: 0.5, y2d: 0, x3d: 0.4358, y3d: -0.1477, z3d: -0.8188 },
    ],
    [{ id: 'b1', source: 'F1', target: 'H1', order: 1 }],
    { customName: 'Hydrogen Fluoride' }
  );

  const hcl = analyzeMolecularGraph(
    [
      { id: 'Cl1', element: 'Cl', x2d: -0.5, y2d: 0, x3d: 0, y3d: 0, z3d: 0 },
      { id: 'H1', element: 'H', x2d: 0.5, y2d: 0, x3d: 0.6085, y3d: -0.2063, z3d: -1.1433 },
    ],
    [{ id: 'b1', source: 'Cl1', target: 'H1', order: 1 }],
    { customName: 'Hydrogen Chloride' }
  );

  for (const mol of [hf, hcl]) {
    assert.equal(mol.pointGroupId, 'Cinfv');
    // Verify canonical x-axis alignment (y3d = 0, z3d = 0)
    for (const a of mol.atoms) {
      assert.ok(Math.abs(a.y3d) < 1e-6 && Math.abs(a.z3d) < 1e-6);
    }

    const labels = mol.molecularOrbitals.map((m) => m.label);
    assert.deepEqual(labels, ['2σ+', '3σ+', '1π(x)', '1π(y)', '4σ+*']);

    const mo3Sigma = mol.molecularOrbitals[1];
    const mo4SigmaStar = mol.molecularOrbitals[4];
    assert.equal(mo3Sigma.bondingType, 'bonding');
    assert.equal(mo4SigmaStar.bondingType, 'antibonding');
    assert.ok(mo4SigmaStar.isLUMO);

    // Verify 4σ+* has both positive and negative 3D phase lobes separated by an internuclear nodal plane
    const mesh4s = buildMorphIsosurfaceKeyframe(mo4SigmaStar, mol.atoms, 1.0, 1.0);
    assert.ok(
      mesh4s.positive.positions.length > 300 && mesh4s.negative.positions.length > 300,
      `Expected both positive and negative lobes in ${mol.formula} 4σ+*, got +${mesh4s.positive.positions.length / 3} / -${mesh4s.negative.positions.length / 3}`
    );
  }
});

test('Full 86-element Periodic Table (H to Rn), principal shell assignment, and Molfile/SMILES/Formula sketcher round-tripping work accurately', async () => {
  const { ELEMENT_PROPERTIES } = await import('./moleculesData.ts');
  const {
    getPrincipalShell,
    generateMolfileFromGraph,
    parseChemicalInputToGraph,
  } = await import('./moEngine.ts');

  const allSymbols = Object.keys(ELEMENT_PROPERTIES);
  assert.equal(allSymbols.length, 86, 'Should support all 86 elements from H (Z=1) to Rn (Z=86)');

  assert.equal(getPrincipalShell('H'), '1');
  assert.equal(getPrincipalShell('Li'), '2');
  assert.equal(getPrincipalShell('Si'), '3');
  assert.equal(getPrincipalShell('Ge'), '4');
  assert.equal(getPrincipalShell('Br'), '4');
  assert.equal(getPrincipalShell('Sn'), '5');
  assert.equal(getPrincipalShell('Pb'), '6');
  assert.equal(getPrincipalShell('Lu'), '6');

  // Test Formula parsing and VSEPR + MO analysis for Silane (SiH4, Td)
  const sih4Parsed = parseChemicalInputToGraph('SiH4', 'Formula')!;
  assert.ok(sih4Parsed && sih4Parsed.atoms.length === 5);
  const sih4 = analyzeMolecularGraph(sih4Parsed.atoms, sih4Parsed.bonds);
  assert.equal(sih4.formula, 'SiH4');
  assert.equal(sih4.pointGroupId, 'Td');
  assert.equal(sih4.valenceElectrons, 8);
  assert.ok(sih4.leftFragmentTitle.includes('3s, 3p'));

  // Test Formula parsing and MO analysis for Hydrogen Bromide (HBr, Cinfv, 4s/4p valence)
  const hbrParsed = parseChemicalInputToGraph('HBr', 'Formula')!;
  const hbr = analyzeMolecularGraph(hbrParsed.atoms, hbrParsed.bonds);
  assert.equal(hbr.formula, 'HBr');
  assert.equal(hbr.pointGroupId, 'Cinfv');
  assert.equal(hbr.valenceElectrons, 8);
  assert.deepEqual(
    hbr.molecularOrbitals.map((m) => m.label),
    ['2σ+', '3σ+', '1π(x)', '1π(y)', '4σ+*']
  );

  // Test MDL V2000 Molfile export & re-import round-trip with formal charges
  const molfile = generateMolfileFromGraph(
    [
      { id: 'N1', element: 'N', x2d: 0, y2d: 0, x3d: 0, y3d: 0, z3d: 0, formalCharge: 1 },
      { id: 'O1', element: 'O', x2d: 1.2, y2d: 0, x3d: 1.2, y3d: 0, z3d: 0, formalCharge: -1 },
    ],
    [{ id: 'b1', source: 'N1', target: 'O1', order: 1, stereo: 'wedge' }],
    'Charged N-O Test'
  );
  assert.ok(molfile.includes('V2000') && molfile.includes('M  END'));
  const reimported = parseChemicalInputToGraph(molfile, 'MDL Molfile')!;
  assert.equal(reimported.atoms.length, 2);
  assert.equal(reimported.atoms[0].formalCharge, 1);
  assert.equal(reimported.atoms[1].formalCharge, -1);
  assert.equal(reimported.bonds[0].stereo, 'wedge');
});

test('Full 5-orbital d and 7-orbital 4f Slater-Koster basis, Hunds rule t2g3 filling, 4f isosurfaces, and Oh Γ_d / Γ_f decomposition work accurately', async () => {
  const { runExtendedHuckelSolver, parseChemicalInputToGraph } = await import('./moEngine.ts');

  // 1. Verify Great Orthogonality Theorem decomposition for Γ_d and Γ_f in Oh
  const oh = POINT_GROUPS.Oh;
  const gammaD = decomposeReducibleRepresentation(oh, [5, -1, 1, -1, 1, 5, -1, -1, 1, 1]);
  assert.equal(gammaD.formula, 'eg + t2g');
  const gammaF = decomposeReducibleRepresentation(oh, [7, 1, -1, -1, -1, -7, 1, -1, 1, 1]);
  assert.equal(gammaF.formula, 'a2u + t1u + t2u');

  // 2. Verify canonical d and f coordination complexes (TiCl4, PtCl4, CrCl6, CeCl6)
  const ticl4 = CANONICAL_MOLECULES.find((m) => m.formula === 'TiCl4')!;
  const ptcl4 = CANONICAL_MOLECULES.find((m) => m.formula === 'PtCl4')!;
  const crcl6 = CANONICAL_MOLECULES.find((m) => m.formula === 'CrCl6')!;
  const cecl6 = CANONICAL_MOLECULES.find((m) => m.formula === 'CeCl6')!;

  assert.equal(ticl4.pointGroupId, 'Td');
  assert.equal(ptcl4.pointGroupId, 'D4h');
  assert.equal(crcl6.pointGroupId, 'Oh');
  assert.equal(cecl6.pointGroupId, 'Oh');

  // CrCl63- (d3, S = 3/2) must have three singly occupied t2g SOMOs (1, 1, 1)
  const crT2g = crcl6.molecularOrbitals.filter((m) => m.degenerateGroup === 'crcl6_2t2g');
  assert.equal(crT2g.length, 3);
  assert.deepEqual(
    crT2g.map((m) => m.occupancy),
    [1, 1, 1]
  );

  // CeCl63- (4f1) must have singly occupied 1a2u (4fxyz) SOMO and watertight 8-lobe 3D isosurface
  const ceHomo = cecl6.molecularOrbitals.find((m) => m.isHOMO)!;
  assert.equal(ceHomo.occupancy, 1);
  assert.equal(ceHomo.lcao[0].orbital, 'fxyz');

  const fxyzMesh = buildMorphIsosurfaceKeyframe(ceHomo, cecl6.atoms, 1.0, 1.0);
  assert.ok(
    fxyzMesh.positive.positions.length > 600 && fxyzMesh.negative.positions.length > 600,
    'CeCl6 4fxyz isosurface must generate rich positive and negative octant lobes'
  );

  // 3. Verify custom EHT solver includes d orbitals for transition metals (FeCl4) and f+d orbitals for lanthanides (NdF3)
  const fecl4Graph = parseChemicalInputToGraph('FeCl4', 'Formula')!;
  const fecl4Eht = runExtendedHuckelSolver(fecl4Graph.atoms, fecl4Graph.bonds);
  const hasFe3d = fecl4Eht.basis.some((b) =>
    ['dz2', 'dx2y2', 'dxy', 'dxz', 'dyz'].includes(b.orbital)
  );
  assert.ok(hasFe3d, 'Expected Extended Hückel solver to include 3d basis functions for FeCl4');

  const ndf3Graph = parseChemicalInputToGraph('NdF3', 'Formula')!;
  const ndf3Anal = analyzeMolecularGraph(ndf3Graph.atoms, ndf3Graph.bonds);
  const hasNd4f = ndf3Anal.molecularOrbitals.some((mo) =>
    mo.lcao.some((c) =>
      ['fz3', 'fxz2', 'fyz2', 'fxyz', 'fz_x2y2', 'fx_x2_3y2', 'fy_3x2_y2'].includes(c.orbital)
    )
  );
  assert.ok(hasNd4f, 'Expected Extended Hückel solver to include 4f basis functions for NdF3');
  assert.ok(ndf3Anal.leftFragmentTitle.includes('4f'), 'Expected lanthanide leftFragmentTitle to list 4f shell');
});

test('Disconnected coordination complexes (hexacyanoferrate [Fe(CN)6]4-, Cr(CO)6, Fe(CO)5, [Ni(CN)4]2-, Ni(CO)4) infer metal-ligand σ bonds and reconstruct 3D Oh / D3h / D4h / Td symmetry', async () => {
  const { parseChemicalInputToGraph, resolveCoordinationBonds } = await import('./moEngine.ts');

  // 1. Hexacyanoferrate ([Fe(CN)6]4-, PubChem CID 9552077):
  // PubChem 2D-SDF stores Fe1 + 6 disconnected C≡N ligands at z=0 with NO Fe-C bonds.
  const rawFeCN6Atoms = [
    { id: 'Fe1', element: 'Fe' as const, x2d: 5.13, y2d: -1.16, x3d: 5.13, y3d: 1.16, z3d: 0 },
    { id: 'N1', element: 'N' as const, x2d: 2.87, y2d: -0.5, x3d: 2.87, y3d: 0.5, z3d: 0 },
    { id: 'N2', element: 'N' as const, x2d: 2.87, y2d: -3.5, x3d: 2.87, y3d: 3.5, z3d: 0 },
    { id: 'N3', element: 'N' as const, x2d: 2.87, y2d: 2.5, x3d: 2.87, y3d: -2.5, z3d: 0 },
    { id: 'N4', element: 'N' as const, x2d: 5.95, y2d: -0.5, x3d: 5.95, y3d: 0.5, z3d: 0 },
    { id: 'N5', element: 'N' as const, x2d: 5.95, y2d: -3.5, x3d: 5.95, y3d: 3.5, z3d: 0 },
    { id: 'N6', element: 'N' as const, x2d: 5.95, y2d: 2.5, x3d: 5.95, y3d: -2.5, z3d: 0 },
    { id: 'C1', element: 'C' as const, x2d: 2.0, y2d: 0.0, x3d: 2.0, y3d: 0.0, z3d: 0 },
    { id: 'C2', element: 'C' as const, x2d: 2.0, y2d: -3.0, x3d: 2.0, y3d: 3.0, z3d: 0 },
    { id: 'C3', element: 'C' as const, x2d: 2.0, y2d: 3.0, x3d: 2.0, y3d: -3.0, z3d: 0 },
    { id: 'C4', element: 'C' as const, x2d: 5.08, y2d: 0.0, x3d: 5.08, y3d: 0.0, z3d: 0 },
    { id: 'C5', element: 'C' as const, x2d: 5.08, y2d: -3.0, x3d: 5.08, y3d: 3.0, z3d: 0 },
    { id: 'C6', element: 'C' as const, x2d: 5.08, y2d: 3.0, x3d: 5.08, y3d: -3.0, z3d: 0 },
  ];
  const rawFeCN6Bonds = [
    { id: 'sdf_b_1', source: 'N1', target: 'C1', order: 3 as const },
    { id: 'sdf_b_2', source: 'N2', target: 'C2', order: 3 as const },
    { id: 'sdf_b_3', source: 'N3', target: 'C3', order: 3 as const },
    { id: 'sdf_b_4', source: 'N4', target: 'C4', order: 3 as const },
    { id: 'sdf_b_5', source: 'N5', target: 'C5', order: 3 as const },
    { id: 'sdf_b_6', source: 'N6', target: 'C6', order: 3 as const },
  ];

  const resolvedFeCN6 = resolveCoordinationBonds(rawFeCN6Atoms, rawFeCN6Bonds, true);
  assert.equal(resolvedFeCN6.bonds.length, 12, 'Expected 6 C≡N bonds + 6 inferred Fe-C coordination bonds');
  const feBonds = resolvedFeCN6.bonds.filter((b) => b.source === 'Fe1' || b.target === 'Fe1');
  assert.equal(feBonds.length, 6, 'Expected all 6 cyanide ligands to bond to Fe1 via Carbon');
  for (const b of feBonds) {
    const partnerId = b.source === 'Fe1' ? b.target : b.source;
    assert.ok(partnerId.startsWith('C'), `Expected Fe1 to bond to C donor atom, got ${partnerId}`);
  }

  const fecn6Analyzed = analyzeMolecularGraph(rawFeCN6Atoms, rawFeCN6Bonds, {
    customName: 'Hexacyanoferrate(II)',
    charge: -4,
  });
  assert.equal(fecn6Analyzed.pointGroupId, 'Oh', 'Hexacyanoferrate must have octahedral (Oh) symmetry');
  assert.equal(fecn6Analyzed.formula, 'C6FeN6', 'Coordinated cyanide carbons must not gain implicit hydrogens');
  const maxAbsZ = fecn6Analyzed.atoms.reduce((m, a) => Math.max(m, Math.abs(a.z3d)), 0);
  assert.ok(maxAbsZ > 1.8, `Expected 3D octahedral axial ligands along ±z (max |z3d| = ${maxAbsZ})`);

  // Verify collinear Fe-C-N軸 on all 6 octahedral branches
  const feCenter = fecn6Analyzed.atoms.find((a) => a.id === 'Fe1')!;
  assert.ok(Math.hypot(feCenter.x3d, feCenter.y3d, feCenter.z3d) < 1e-4, 'Fe1 must sit at 3D origin');
  for (let k = 1; k <= 6; k++) {
    const cAtom = fecn6Analyzed.atoms.find((a) => a.id === `C${k}`)!;
    const nAtom = fecn6Analyzed.atoms.find((a) => a.id === `N${k}`)!;
    const rC = Math.hypot(cAtom.x3d, cAtom.y3d, cAtom.z3d);
    const rN = Math.hypot(nAtom.x3d, nAtom.y3d, nAtom.z3d);
    assert.ok(rC > 1.7 && rC < 2.2, `Expected Fe-C bond distance ~1.93 Å, got ${rC}`);
    assert.ok(rN > rC + 1.0, `Expected terminal N atom radially outward from C (rN=${rN}, rC=${rC})`);
  }

  // 2. Disconnected PubChem SMILES for Cr(CO)6, Fe(CO)5, [Ni(CN)4]2-, and Ni(CO)4
  const crCO6 = parseChemicalInputToGraph('[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+].[Cr]', 'SMILES')!;
  const crCO6Entry = analyzeMolecularGraph(crCO6.atoms, crCO6.bonds);
  assert.equal(crCO6Entry.pointGroupId, 'Oh', 'Cr(CO)6 must be Octahedral (Oh)');
  assert.equal(crCO6Entry.bonds.filter((b) => b.source === 'Cr1' || b.target === 'Cr1').length, 6);

  const feCO5 = parseChemicalInputToGraph('[C]=O.[C]=O.[C]=O.[C]=O.[C]=O.[Fe]', 'SMILES')!;
  const feCO5Entry = analyzeMolecularGraph(feCO5.atoms, feCO5.bonds);
  assert.equal(feCO5Entry.pointGroupId, 'D3h', 'Fe(CO)5 must be Trigonal Bipyramidal (D3h)');
  assert.equal(feCO5Entry.bonds.filter((b) => b.source === 'Fe1' || b.target === 'Fe1').length, 5);

  const niCN4 = parseChemicalInputToGraph('[C-]#N.[C-]#N.[C-]#N.[C-]#N.[Ni+2]', 'SMILES')!;
  const niCN4Entry = analyzeMolecularGraph(niCN4.atoms, niCN4.bonds, { charge: -2 });
  assert.equal(niCN4Entry.pointGroupId, 'D4h', '[Ni(CN)4]2- (d8) must be Square Planar (D4h)');

  const niCO4 = parseChemicalInputToGraph('[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+].[Ni]', 'SMILES')!;
  const niCO4Entry = analyzeMolecularGraph(niCO4.atoms, niCO4.bonds);
  assert.equal(niCO4Entry.pointGroupId, 'Td', 'Ni(CO)4 (d10) must be Tetrahedral (Td)');
});

test('Multi-center polycyclic aromatics (Coronene D6h, Triphenylene D3h, Naphthalene D2h, Phenanthrene C2v), heterocycles (Borazine D3h, Pyridine C2v), 3D cages (Cubane Oh, Adamantane Td), and polyynes (Diacetylene Dinfh, Cyanoacetylene Cinfv) classify accurately', () => {
  // 1. Coronene (C24H12): 6 inner hub carbons at r = 1.40 Å, 6 radial spoke carbons at r = 2.80 Å, 12 outer rim carbons at r = 3.70 Å
  const coroneneAtoms = [];
  const coroneneBonds = [];
  for (let k = 0; k < 6; k++) {
    const phi = (k * Math.PI) / 3;
    coroneneAtoms.push({
      id: `C_in_${k}`,
      element: 'C' as const,
      x2d: Number((1.4 * Math.cos(phi)).toFixed(3)),
      y2d: Number((1.4 * Math.sin(phi)).toFixed(3)),
      x3d: Number((1.4 * Math.cos(phi)).toFixed(3)),
      y3d: Number((1.4 * Math.sin(phi)).toFixed(3)),
      z3d: 0,
    });
    coroneneAtoms.push({
      id: `C_sp_${k}`,
      element: 'C' as const,
      x2d: Number((2.8 * Math.cos(phi)).toFixed(3)),
      y2d: Number((2.8 * Math.sin(phi)).toFixed(3)),
      x3d: Number((2.8 * Math.cos(phi)).toFixed(3)),
      y3d: Number((2.8 * Math.sin(phi)).toFixed(3)),
      z3d: 0,
    });
    for (const delta of [-Math.PI / 9, Math.PI / 9]) {
      const tag = delta < 0 ? 'a' : 'b';
      const ang = phi + delta;
      coroneneAtoms.push({
        id: `C_out_${k}_${tag}`,
        element: 'C' as const,
        x2d: Number((3.7 * Math.cos(ang)).toFixed(3)),
        y2d: Number((3.7 * Math.sin(ang)).toFixed(3)),
        x3d: Number((3.7 * Math.cos(ang)).toFixed(3)),
        y3d: Number((3.7 * Math.sin(ang)).toFixed(3)),
        z3d: 0,
      });
    }
  }
  for (let k = 0; k < 6; k++) {
    const next = (k + 1) % 6;
    coroneneBonds.push(
      { id: `b_in_${k}`, source: `C_in_${k}`, target: `C_in_${next}`, order: (k % 2 === 0 ? 2 : 1) as 1 | 2 },
      { id: `b_rad_${k}`, source: `C_in_${k}`, target: `C_sp_${k}`, order: 1 as const },
      { id: `b_sp_b_${k}`, source: `C_sp_${k}`, target: `C_out_${k}_b`, order: (k % 2 === 0 ? 1 : 2) as 1 | 2 },
      { id: `b_rim_${k}`, source: `C_out_${k}_b`, target: `C_out_${next}_a`, order: (k % 2 === 0 ? 2 : 1) as 1 | 2 },
      { id: `b_sp_a_${next}`, source: `C_out_${next}_a`, target: `C_sp_${next}`, order: (next % 2 === 0 ? 2 : 1) as 1 | 2 }
    );
  }
  assert.equal(detectPointGroup(coroneneAtoms, coroneneBonds), 'D6h', 'Coronene (C24H12) must be D6h');
  const coroneneEntry = analyzeMolecularGraph(coroneneAtoms, coroneneBonds, { customName: 'Coronene' });
  assert.equal(coroneneEntry.pointGroupId, 'D6h');
  assert.equal(coroneneEntry.formula, 'C24H12');
  const coroneneHomo = coroneneEntry.molecularOrbitals.find((m) => m.isHOMO)!;
  assert.equal(coroneneHomo.irrepId, 'E2u', 'Coronene (24π e-) HOMO must have doubly degenerate e2u symmetry');

  // 2. Borazine (B3N3H6, inorganic benzene): alternating B and N hexagon -> D3h (not D6h or Cs)
  const borazineAtoms = [];
  const borazineBonds = [];
  for (let k = 0; k < 6; k++) {
    const phi = (k * Math.PI) / 3;
    const el = (k % 2 === 0 ? 'B' : 'N') as 'B' | 'N';
    borazineAtoms.push({
      id: `${el}${k + 1}`,
      element: el,
      x2d: Number((1.44 * Math.cos(phi)).toFixed(3)),
      y2d: Number((1.44 * Math.sin(phi)).toFixed(3)),
      x3d: Number((1.44 * Math.cos(phi)).toFixed(3)),
      y3d: Number((1.44 * Math.sin(phi)).toFixed(3)),
      z3d: 0,
    });
  }
  for (let k = 0; k < 6; k++) {
    const a1 = borazineAtoms[k].id;
    const a2 = borazineAtoms[(k + 1) % 6].id;
    borazineBonds.push({ id: `bb_${k}`, source: a1, target: a2, order: 1 as const });
  }
  assert.equal(detectPointGroup(borazineAtoms, borazineBonds), 'D3h', 'Borazine (B3N3H6) must be D3h');

  // 3. Naphthalene (C10H8): two fused benzene rings -> D2h
  const naphAtoms = [
    { id: 'C1', element: 'C' as const, x2d: 0, y2d: 0.71, x3d: 0, y3d: 0.71, z3d: 0 },
    { id: 'C2', element: 'C' as const, x2d: 0, y2d: -0.71, x3d: 0, y3d: -0.71, z3d: 0 },
    { id: 'C3', element: 'C' as const, x2d: 1.23, y2d: 1.4, x3d: 1.23, y3d: 1.4, z3d: 0 },
    { id: 'C4', element: 'C' as const, x2d: 2.45, y2d: 0.71, x3d: 2.45, y3d: 0.71, z3d: 0 },
    { id: 'C5', element: 'C' as const, x2d: 2.45, y2d: -0.71, x3d: 2.45, y3d: -0.71, z3d: 0 },
    { id: 'C6', element: 'C' as const, x2d: 1.23, y2d: -1.4, x3d: 1.23, y3d: -1.4, z3d: 0 },
    { id: 'C7', element: 'C' as const, x2d: -1.23, y2d: 1.4, x3d: -1.23, y3d: 1.4, z3d: 0 },
    { id: 'C8', element: 'C' as const, x2d: -2.45, y2d: 0.71, x3d: -2.45, y3d: 0.71, z3d: 0 },
    { id: 'C9', element: 'C' as const, x2d: -2.45, y2d: -0.71, x3d: -2.45, y3d: -0.71, z3d: 0 },
    { id: 'C10', element: 'C' as const, x2d: -1.23, y2d: -1.4, x3d: -1.23, y3d: -1.4, z3d: 0 },
  ];
  const naphBonds = [
    { id: 'nb1', source: 'C1', target: 'C2', order: 2 as const },
    { id: 'nb2', source: 'C1', target: 'C3', order: 1 as const },
    { id: 'nb3', source: 'C3', target: 'C4', order: 2 as const },
    { id: 'nb4', source: 'C4', target: 'C5', order: 1 as const },
    { id: 'nb5', source: 'C5', target: 'C6', order: 2 as const },
    { id: 'nb6', source: 'C6', target: 'C2', order: 1 as const },
    { id: 'nb7', source: 'C1', target: 'C7', order: 1 as const },
    { id: 'nb8', source: 'C7', target: 'C8', order: 2 as const },
    { id: 'nb9', source: 'C8', target: 'C9', order: 1 as const },
    { id: 'nb10', source: 'C9', target: 'C10', order: 2 as const },
    { id: 'nb11', source: 'C10', target: 'C2', order: 1 as const },
  ];
  assert.equal(detectPointGroup(naphAtoms, naphBonds), 'D2h', 'Naphthalene (C10H8) must be D2h');

  // 4. Pyridine (C5H5N): 1 N + 5 C hexagon -> C2v
  const pyrAtoms = borazineAtoms.map((a, idx) => ({
    ...a,
    id: idx === 0 ? 'N1' : `C${idx}`,
    element: (idx === 0 ? 'N' : 'C') as 'N' | 'C',
  }));
  const pyrBonds = borazineBonds.map((b, idx) => ({
    ...b,
    source: pyrAtoms[idx].id,
    target: pyrAtoms[(idx + 1) % 6].id,
  }));
  assert.equal(detectPointGroup(pyrAtoms, pyrBonds), 'C2v', 'Pyridine (C5H5N) must be C2v');

  // 5. Cubane (C8H8): 8 vertices of a cube at (±0.78, ±0.78, ±0.78) -> Oh
  const cubaneAtoms = [];
  const cubaneBonds = [];
  let cIdx = 1;
  for (const sx of [-0.78, 0.78]) {
    for (const sy of [-0.78, 0.78]) {
      for (const sz of [-0.78, 0.78]) {
        cubaneAtoms.push({
          id: `C${cIdx++}`,
          element: 'C' as const,
          x2d: sx,
          y2d: sy,
          x3d: sx,
          y3d: sy,
          z3d: sz,
        });
      }
    }
  }
  for (let i = 0; i < cubaneAtoms.length; i++) {
    for (let j = i + 1; j < cubaneAtoms.length; j++) {
      const a = cubaneAtoms[i];
      const b = cubaneAtoms[j];
      if (Math.abs(Math.hypot(a.x3d - b.x3d, a.y3d - b.y3d, a.z3d - b.z3d) - 1.56) < 0.05) {
        cubaneBonds.push({ id: `cb_${i}_${j}`, source: a.id, target: b.id, order: 1 as const });
      }
    }
  }
  assert.equal(detectPointGroup(cubaneAtoms, cubaneBonds), 'Oh', 'Cubane (C8H8) must be Oh');

  // 6. White Phosphorus (P4) tetrahedral cage -> Td
  const p4Atoms = [
    { id: 'P1', element: 'P' as const, x2d: 0.78, y2d: 0.78, x3d: 0.78, y3d: 0.78, z3d: 0.78 },
    { id: 'P2', element: 'P' as const, x2d: -0.78, y2d: -0.78, x3d: -0.78, y3d: -0.78, z3d: 0.78 },
    { id: 'P3', element: 'P' as const, x2d: -0.78, y2d: 0.78, x3d: -0.78, y3d: 0.78, z3d: -0.78 },
    { id: 'P4', element: 'P' as const, x2d: 0.78, y2d: -0.78, x3d: 0.78, y3d: -0.78, z3d: -0.78 },
  ];
  const p4Bonds = [
    { id: 'pb1', source: 'P1', target: 'P2', order: 1 as const },
    { id: 'pb2', source: 'P1', target: 'P3', order: 1 as const },
    { id: 'pb3', source: 'P1', target: 'P4', order: 1 as const },
    { id: 'pb4', source: 'P2', target: 'P3', order: 1 as const },
    { id: 'pb5', source: 'P2', target: 'P4', order: 1 as const },
    { id: 'pb6', source: 'P3', target: 'P4', order: 1 as const },
  ];
  assert.equal(detectPointGroup(p4Atoms, p4Bonds), 'Td', 'P4 tetrahedral cage must be Td');

  // 7. Linear polyatomics: Diacetylene (C4H2 -> Dinfh) and Cyanoacetylene (HC3N -> Cinfv)
  const c4h2Atoms = [
    { id: 'H1', element: 'H' as const, x2d: -3.0, y2d: 0, x3d: -3.0, y3d: 0, z3d: 0 },
    { id: 'C1', element: 'C' as const, x2d: -1.9, y2d: 0, x3d: -1.9, y3d: 0, z3d: 0 },
    { id: 'C2', element: 'C' as const, x2d: -0.7, y2d: 0, x3d: -0.7, y3d: 0, z3d: 0 },
    { id: 'C3', element: 'C' as const, x2d: 0.7, y2d: 0, x3d: 0.7, y3d: 0, z3d: 0 },
    { id: 'C4', element: 'C' as const, x2d: 1.9, y2d: 0, x3d: 1.9, y3d: 0, z3d: 0 },
    { id: 'H2', element: 'H' as const, x2d: 3.0, y2d: 0, x3d: 3.0, y3d: 0, z3d: 0 },
  ];
  const c4h2Bonds = [
    { id: 'db1', source: 'H1', target: 'C1', order: 1 as const },
    { id: 'db2', source: 'C1', target: 'C2', order: 3 as const },
    { id: 'db3', source: 'C2', target: 'C3', order: 1 as const },
    { id: 'db4', source: 'C3', target: 'C4', order: 3 as const },
    { id: 'db5', source: 'C4', target: 'H2', order: 1 as const },
  ];
  assert.equal(detectPointGroup(c4h2Atoms, c4h2Bonds), 'Dinfh', 'Diacetylene (C4H2) must be Dinfh');

  const hc3nAtoms = c4h2Atoms.slice(0, 5).map((a, idx) => (idx === 4 ? { ...a, id: 'N1', element: 'N' as const } : a));
  const hc3nBonds = [
    { id: 'hb1', source: 'H1', target: 'C1', order: 1 as const },
    { id: 'hb2', source: 'C1', target: 'C2', order: 1 as const },
    { id: 'hb3', source: 'C2', target: 'C3', order: 1 as const },
    { id: 'hb4', source: 'C3', target: 'N1', order: 3 as const },
  ];
  hc3nBonds[1].order = 3;
  assert.equal(detectPointGroup(hc3nAtoms, hc3nBonds), 'Cinfv', 'Cyanoacetylene (HC3N) must be Cinfv');
});

test('Long-chain alkanes (dodecane, hexadecane, eicosane), fatty acids (stearic acid), branched terpenoids (pristane), and cycloalkanes (cyclohexane) reconstruct true 3D tetrahedral sp3 carbon geometry (~109.5°)', () => {
  const computeMeanTetrahedralAngles = (entry: ReturnType<typeof analyzeMolecularGraph>) => {
    const sp3Carbons = entry.atoms.filter((a) => {
      if (a.element !== 'C') return false;
      const conn = entry.bonds.filter((b) => b.source === a.id || b.target === a.id);
      return conn.length === 4 && conn.every((b) => b.order === 1);
    });
    const hchAngles: number[] = [];
    const cccAngles: number[] = [];
    for (const c of sp3Carbons) {
      const nbs = entry.bonds
        .filter((b) => b.source === c.id || b.target === c.id)
        .map((b) => entry.atoms.find((x) => x.id === (b.source === c.id ? b.target : b.source))!);
      const hNbs = nbs.filter((x) => x.element === 'H');
      const cNbs = nbs.filter((x) => x.element === 'C');
      if (hNbs.length >= 2) {
        const u1 = [hNbs[0].x3d - c.x3d, hNbs[0].y3d - c.y3d, hNbs[0].z3d - c.z3d];
        const u2 = [hNbs[1].x3d - c.x3d, hNbs[1].y3d - c.y3d, hNbs[1].z3d - c.z3d];
        const l1 = Math.hypot(u1[0], u1[1], u1[2]) || 1;
        const l2 = Math.hypot(u2[0], u2[1], u2[2]) || 1;
        const dot = (u1[0] * u2[0] + u1[1] * u2[1] + u1[2] * u2[2]) / (l1 * l2);
        hchAngles.push((Math.acos(Math.max(-1, Math.min(1, dot))) * 180) / Math.PI);
      }
      if (cNbs.length >= 2) {
        const u1 = [cNbs[0].x3d - c.x3d, cNbs[0].y3d - c.y3d, cNbs[0].z3d - c.z3d];
        const u2 = [cNbs[1].x3d - c.x3d, cNbs[1].y3d - c.y3d, cNbs[1].z3d - c.z3d];
        const l1 = Math.hypot(u1[0], u1[1], u1[2]) || 1;
        const l2 = Math.hypot(u2[0], u2[1], u2[2]) || 1;
        const dot = (u1[0] * u2[0] + u1[1] * u2[1] + u1[2] * u2[2]) / (l1 * l2);
        cccAngles.push((Math.acos(Math.max(-1, Math.min(1, dot))) * 180) / Math.PI);
      }
    }
    const meanHCH = hchAngles.reduce((a, b) => a + b, 0) / Math.max(1, hchAngles.length);
    const meanCCC = cccAngles.reduce((a, b) => a + b, 0) / Math.max(1, cccAngles.length);
    const maxCarbonZ = Math.max(...sp3Carbons.map((a) => Math.abs(a.z3d)));
    const maxAllZ = Math.max(...entry.atoms.map((a) => Math.abs(a.z3d)));
    return { sp3Count: sp3Carbons.length, meanHCH, meanCCC, maxCarbonZ, maxAllZ };
  };

  // 1. Dodecane (C12H26, 38 atoms) & Eicosane (C20H42, 62 atoms)
  const dodecane = analyzeMolecularGraph(
    parseChemicalInputToGraph('CCCCCCCCCCCC', 'SMILES')!.atoms,
    parseChemicalInputToGraph('CCCCCCCCCCCC', 'SMILES')!.bonds
  );
  assert.equal(dodecane.formula, 'C12H26');
  assert.equal(dodecane.atoms.length, 38);
  const dodStats = computeMeanTetrahedralAngles(dodecane);
  assert.equal(dodStats.sp3Count, 12);
  assert.ok(dodStats.maxCarbonZ > 0.24, `Expected 3D sp3 carbon puckering > 0.24 Å, got ${dodStats.maxCarbonZ}`);
  assert.ok(dodStats.maxAllZ > 1.2, `Expected 3D tetrahedral C-H depth > 1.2 Å, got ${dodStats.maxAllZ}`);
  assert.ok(
    Math.abs(dodStats.meanHCH - 109.47) < 1.0 && Math.abs(dodStats.meanCCC - 109.47) < 2.0,
    `Expected tetrahedral ~109.5° angles in Dodecane, got HCH=${dodStats.meanHCH.toFixed(2)}°, CCC=${dodStats.meanCCC.toFixed(2)}°`
  );

  const eicosaneGraph = parseChemicalInputToGraph('CCCCCCCCCCCCCCCCCCCC', 'SMILES')!;
  const eicosane = analyzeMolecularGraph(eicosaneGraph.atoms, eicosaneGraph.bonds);
  assert.equal(eicosane.formula, 'C20H42');
  assert.equal(eicosane.atoms.length, 62);
  const eicStats = computeMeanTetrahedralAngles(eicosane);
  assert.equal(eicStats.sp3Count, 20);
  assert.ok(Math.abs(eicStats.meanHCH - 109.47) < 1.0 && Math.abs(eicStats.meanCCC - 109.47) < 2.0);

  // 2. Stearic Acid (C18H36O2, 56 atoms): 17 sp3 chain carbons + 1 sp2 carboxylic acid carbon
  const stearicGraph = parseChemicalInputToGraph('CCCCCCCCCCCCCCCCCC(=O)O', 'SMILES')!;
  const stearic = analyzeMolecularGraph(stearicGraph.atoms, stearicGraph.bonds);
  assert.equal(stearic.formula, 'C18H36O2');
  assert.equal(stearic.atoms.length, 56);
  const stStats = computeMeanTetrahedralAngles(stearic);
  assert.equal(stStats.sp3Count, 17, 'Expected 17 sp3 alkyl carbons and 1 sp2 carbonyl carbon in Stearic Acid');
  assert.ok(stStats.maxAllZ > 1.2 && Math.abs(stStats.meanHCH - 109.47) < 1.0);

  // 3. Cyclohexane (C1CCCCC1, C6H12, 18 atoms): 3D chair conformation with ~109.5° CCC and HCH angles
  const chxGraph = parseChemicalInputToGraph('C1CCCCC1', 'SMILES')!;
  const chx = analyzeMolecularGraph(chxGraph.atoms, chxGraph.bonds);
  assert.equal(chx.formula, 'C6H12');
  assert.equal(chx.atoms.length, 18);
  const chxStats = computeMeanTetrahedralAngles(chx);
  assert.equal(chxStats.sp3Count, 6);
  assert.ok(
    Math.abs(chxStats.meanCCC - 109.47) < 1.5 && Math.abs(chxStats.meanHCH - 109.47) < 1.0,
    `Expected 3D chair cyclohexane ~109.5° angles, got CCC=${chxStats.meanCCC.toFixed(2)}°, HCH=${chxStats.meanHCH.toFixed(2)}°`
  );

  // 4. Branched isoprenoid chain Pristane (C19H40, 59 atoms): methyl branches project out of the plane (maxCarbonZ > 1.0 Å)
  const pristaneGraph = parseChemicalInputToGraph('CC(C)CCCC(C)CCCC(C)CCCC(C)C', 'SMILES')!;
  const pristane = analyzeMolecularGraph(pristaneGraph.atoms, pristaneGraph.bonds);
  assert.equal(pristane.formula, 'C19H40');
  const prStats = computeMeanTetrahedralAngles(pristane);
  assert.ok(prStats.maxCarbonZ > 1.0, `Expected out-of-plane tetrahedral methyl branches > 1.0 Å, got ${prStats.maxCarbonZ}`);
  assert.ok(Math.abs(prStats.meanHCH - 109.47) < 1.0);
});

test("Complexes with complex ligands (Wilkinson's catalyst, Vaska's complex, Pd(PPh3)4, Cisplatin, [Co(en)3]3+, [Ru(bpy)3]2+, Fe(acac)3, [Fe(C2O4)3]3-, Ni(dmgH)2, Ferrocene, Nitroprusside dihydrate, Zeise's salt) connect via true donor atoms and assign accurate 3D geometry & point groups", () => {
  const pph3 = 'C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3';

  // 1. Wilkinson's Catalyst: Chlorotris(triphenylphosphine)rhodium(I), [RhCl(PPh3)3]
  const wilkinsonGraph = parseChemicalInputToGraph(
    `${pph3}.${pph3}.${pph3}.[Cl-].[Rh+]`,
    'SMILES'
  )!;
  const wilkinson = analyzeMolecularGraph(wilkinsonGraph.atoms, wilkinsonGraph.bonds, {
    customName: "Wilkinson's Catalyst",
  });
  const rhAtom = wilkinson.atoms.find((a) => a.element === 'Rh')!;
  assert.ok(rhAtom, "Wilkinson's catalyst must contain Rh center");
  const rhPartners = wilkinson.bonds
    .filter((b) => b.source === rhAtom.id || b.target === rhAtom.id)
    .map((b) => wilkinson.atoms.find((a) => a.id === (b.source === rhAtom.id ? b.target : b.source))!.element)
    .sort();
  assert.deepEqual(
    rhPartners,
    ['Cl', 'P', 'P', 'P'],
    "Rh in Wilkinson's catalyst must bond exclusively to Cl + 3 P donors (not phenyl carbons)"
  );
  assert.equal(
    wilkinson.pointGroupId,
    'C2v',
    "Wilkinson's catalyst (d8 square-planar ML3X) must have C2v symmetry"
  );
  const wilkinsonMaxZ = Math.max(...wilkinson.atoms.map((a) => Math.abs(a.z3d)));
  assert.ok(
    wilkinsonMaxZ > 1.5,
    `Expected 3D tetrahedral propeller phenyl rings on -PPh3 (max |z3d| = ${wilkinsonMaxZ})`
  );
  const hasRh4d = wilkinson.molecularOrbitals.some((mo) =>
    mo.lcao.some(
      (c) =>
        c.atomId === rhAtom.id &&
        ['dz2', 'dx2y2', 'dxy', 'dxz', 'dyz'].includes(c.orbital)
    )
  );
  assert.ok(hasRh4d, "Wilkinson's catalyst MO basis must retain Rh 4d orbitals");

  // 2. Vaska's Complex: trans-[IrCl(CO)(PPh3)2]
  const vaskaGraph = parseChemicalInputToGraph(
    `[C-]#[O+].${pph3}.${pph3}.[Cl-].[Ir+]`,
    'SMILES'
  )!;
  const vaska = analyzeMolecularGraph(vaskaGraph.atoms, vaskaGraph.bonds);
  const irAtom = vaska.atoms.find((a) => a.element === 'Ir')!;
  const irPartners = vaska.bonds
    .filter((b) => b.source === irAtom.id || b.target === irAtom.id)
    .map((b) => vaska.atoms.find((a) => a.id === (b.source === irAtom.id ? b.target : b.source))!.element)
    .sort();
  assert.deepEqual(irPartners, ['C', 'Cl', 'P', 'P'], "Vaska's complex Ir must bond to C(O), Cl, and 2 P");
  assert.equal(vaska.pointGroupId, 'C2v', "Vaska's complex (trans-ML2XY) must have C2v symmetry");

  // 3. Tetrakis(triphenylphosphine)palladium(0): Pd(PPh3)4 (d10 tetrahedral -> Td)
  const pdPPh3_4Graph = parseChemicalInputToGraph(
    `${pph3}.${pph3}.${pph3}.${pph3}.[Pd]`,
    'SMILES'
  )!;
  const pdPPh3_4 = analyzeMolecularGraph(pdPPh3_4Graph.atoms, pdPPh3_4Graph.bonds);
  const pdAtom = pdPPh3_4.atoms.find((a) => a.element === 'Pd')!;
  const pdPartners = pdPPh3_4.bonds
    .filter((b) => b.source === pdAtom.id || b.target === pdAtom.id)
    .map((b) => pdPPh3_4.atoms.find((a) => a.id === (b.source === pdAtom.id ? b.target : b.source))!.element);
  assert.deepEqual(pdPartners, ['P', 'P', 'P', 'P'], 'Pd(PPh3)4 must connect all 4 PPh3 ligands via P');
  assert.equal(pdPPh3_4.pointGroupId, 'Td', 'Pd(PPh3)4 (d10 M0L4) must have Td symmetry');

  // 4. Cisplatin: cis-[PtCl2(NH3)2]
  const cisplatinGraph = parseChemicalInputToGraph('N.N.Cl[Pt]Cl', 'SMILES')!;
  const cisplatin = analyzeMolecularGraph(cisplatinGraph.atoms, cisplatinGraph.bonds);
  const ptAtom = cisplatin.atoms.find((a) => a.element === 'Pt')!;
  const ptPartners = cisplatin.bonds
    .filter((b) => b.source === ptAtom.id || b.target === ptAtom.id)
    .map((b) => cisplatin.atoms.find((a) => a.id === (b.source === ptAtom.id ? b.target : b.source))!.element)
    .sort();
  assert.deepEqual(ptPartners, ['Cl', 'Cl', 'N', 'N']);
  assert.equal(cisplatin.pointGroupId, 'C2v', 'Cisplatin must have C2v symmetry');

  // 5. Tris-bidentate octahedral propeller chelates: [Co(en)3]3+, [Ru(bpy)3]2+, Fe(acac)3, [Fe(C2O4)3]3- -> D3
  const coEn3Graph = parseChemicalInputToGraph('C(CN)N.C(CN)N.C(CN)N.[Co+3]', 'SMILES')!;
  const coEn3 = analyzeMolecularGraph(coEn3Graph.atoms, coEn3Graph.bonds, { charge: 3 });
  assert.equal(coEn3.pointGroupId, 'D3', '[Co(en)3]3+ must have D3 tris-chelate propeller symmetry');
  const coAtom = coEn3.atoms.find((a) => a.element === 'Co')!;
  assert.equal(
    coEn3.bonds.filter((b) => b.source === coAtom.id || b.target === coAtom.id).length,
    6,
    '[Co(en)3]3+ must form 6 Co-N chelate bonds'
  );

  const bpy = 'C1=CC=NC(=C1)C2=CC=CC=N2';
  const ruBpy3Graph = parseChemicalInputToGraph(`${bpy}.${bpy}.${bpy}.[Ru+2]`, 'SMILES')!;
  const ruBpy3 = analyzeMolecularGraph(ruBpy3Graph.atoms, ruBpy3Graph.bonds, { charge: 2 });
  assert.equal(ruBpy3.pointGroupId, 'D3', '[Ru(bpy)3]2+ must have D3 tris-chelate propeller symmetry');
  const ruAtom = ruBpy3.atoms.find((a) => a.element === 'Ru')!;
  assert.equal(
    ruBpy3.bonds.filter((b) => b.source === ruAtom.id || b.target === ruAtom.id).length,
    6,
    '[Ru(bpy)3]2+ must form 6 Ru-N chelate bonds'
  );

  const feAcac3Graph = parseChemicalInputToGraph(
    'CC(=CC(=O)C)[O-].CC(=CC(=O)C)[O-].CC(=CC(=O)C)[O-].[Fe+3]',
    'SMILES'
  )!;
  const feAcac3 = analyzeMolecularGraph(feAcac3Graph.atoms, feAcac3Graph.bonds);
  assert.equal(feAcac3.pointGroupId, 'D3', 'Fe(acac)3 must have D3 symmetry');

  const feOx3Graph = parseChemicalInputToGraph(
    'C(=O)(C(=O)[O-])[O-].C(=O)(C(=O)[O-])[O-].C(=O)(C(=O)[O-])[O-].[Fe+3]',
    'SMILES'
  )!;
  const feOx3 = analyzeMolecularGraph(feOx3Graph.atoms, feOx3Graph.bonds, { charge: -3 });
  assert.equal(feOx3.pointGroupId, 'D3', '[Fe(C2O4)3]3- must have D3 symmetry');

  // 6. Square-planar chelates: Ni(dmgH)2 -> D2h, Oxaliplatin -> C2v
  const niDmg2Graph = parseChemicalInputToGraph(
    'CC(=NO)C(=NO)C.CC(=NO)C(=NO)C.[Ni]',
    'SMILES'
  )!;
  const niDmg2 = analyzeMolecularGraph(niDmg2Graph.atoms, niDmg2Graph.bonds);
  assert.equal(niDmg2.pointGroupId, 'D2h', 'Ni(dmgH)2 must have D2h symmetry');

  const oxaliplatinGraph = parseChemicalInputToGraph(
    'C1CCC(C(C1)N)N.C(=O)(C(=O)O)O.[Pt]',
    'SMILES'
  )!;
  const oxaliplatin = analyzeMolecularGraph(oxaliplatinGraph.atoms, oxaliplatinGraph.bonds);
  assert.equal(oxaliplatin.pointGroupId, 'C2v', 'Oxaliplatin must have C2v symmetry');

  // 7. Metallocene sandwich complexes: Ferrocene Fe(C5H5)2 -> D5d
  const ferroceneGraph = parseChemicalInputToGraph(
    '[CH-]1C=CC=C1.[CH-]1C=CC=C1.[Fe+2]',
    'SMILES'
  )!;
  const ferrocene = analyzeMolecularGraph(ferroceneGraph.atoms, ferroceneGraph.bonds);
  assert.equal(ferrocene.pointGroupId, 'D5d', 'Ferrocene Fe(C5H5)2 must have D5d staggered sandwich symmetry');
  const feSandwich = ferrocene.atoms.find((a) => a.element === 'Fe')!;
  assert.equal(
    ferrocene.bonds.filter((b) => b.source === feSandwich.id || b.target === feSandwich.id).length,
    10,
    'Ferrocene must have 10 η5 Fe-C coordination bonds'
  );

  // 8. Hydrated & η2-alkene complexes: Sodium nitroprusside dihydrate (-> C4v) & Zeise's salt (-> C2v)
  const nitroprussideGraph = parseChemicalInputToGraph(
    '[C-]#N.[C-]#N.[C-]#N.[C-]#N.[C-]#N.N#[O+].O.O.[Fe-2]',
    'SMILES'
  )!;
  const nitroprusside = analyzeMolecularGraph(nitroprussideGraph.atoms, nitroprussideGraph.bonds, {
    charge: -2,
  });
  assert.equal(
    nitroprusside.pointGroupId,
    'C4v',
    'Nitroprusside dihydrate must strip uncoordinated lattice .2H2O and yield C4v [Fe(CN)5(NO)]2-'
  );

  const zeiseGraph = parseChemicalInputToGraph('C=C.O.Cl[Pt-](Cl)Cl', 'SMILES')!;
  const zeise = analyzeMolecularGraph(zeiseGraph.atoms, zeiseGraph.bonds, { charge: -1 });
  assert.equal(
    zeise.pointGroupId,
    'C2v',
    "Zeise's salt [PtCl3(C2H4)]- must strip lattice .H2O, coordinate ethylene, and yield C2v symmetry"
  );
});

test("Dewar-Chatt-Duncanson (DCD) model for transition-metal η²-alkene complexes (Zeise's salt [PtCl3(η²-C2H4)]⁻, [PdCl3(η²-C2H4)]⁻, and custom sketches)", () => {
  // 1. Disconnected Zeise's salt monohydrate SMILES: C=C.O.Cl[Pt-](Cl)Cl.[K+]
  const zeiseParsed = parseChemicalInputToGraph('C=C.O.Cl[Pt-](Cl)Cl.[K+]', 'SMILES')!;
  const zeise = analyzeMolecularGraph(zeiseParsed.atoms, zeiseParsed.bonds);
  assert.equal(zeise.pointGroupId, 'C2v', "Zeise's salt must have C2v symmetry");

  const pt = zeise.atoms.find((a) => a.element === 'Pt')!;
  const carbons = zeise.atoms.filter((a) => a.element === 'C');
  const hydrogens = zeise.atoms.filter((a) => a.element === 'H');
  assert.equal(carbons.length, 2, 'Ethylene ligand must have 2 carbons');
  assert.equal(hydrogens.length, 4, 'Ethylene ligand must have 4 hydrogens');

  // Verify η² dihapto coordination: Pt is bonded to BOTH C1 and C2 of ethylene (plus 3 Cl = 5 bonds total)
  const ptBonds = zeise.bonds.filter((b) => b.source === pt.id || b.target === pt.id);
  assert.equal(
    ptBonds.length,
    5,
    'Pt in Zeise’s salt must have 5 coordination bonds (3 Pt-Cl + 2 η² Pt-C bonds)'
  );
  for (const c of carbons) {
    const hasPtCBond = ptBonds.some((b) => b.source === c.id || b.target === c.id);
    assert.ok(hasPtCBond, `Alkene carbon ${c.id} must be η²-coordinated to Pt`);
  }

  // Verify 3D Dewar-Chatt-Duncanson geometry:
  // (a) C=C bond axis is perpendicular to the PtCl3 xy plane (|Δz(C1,C2)| ≈ 1.376 Å, elongated from 1.339 Å)
  const [c1, c2] = carbons;
  const dzCC = Math.abs(c2.z3d - c1.z3d);
  const rCC = Math.hypot(c2.x3d - c1.x3d, c2.y3d - c1.y3d, c2.z3d - c1.z3d);
  assert.ok(
    dzCC > 1.30 && Math.abs(rCC - 1.376) < 0.02,
    `Expected perpendicular elongated C=C bond (~1.376 Å along z), got r(C=C)=${rCC.toFixed(3)} Å, |Δz|=${dzCC.toFixed(3)} Å`
  );
  // (b) All 4 C-H bonds bend back away from Pt (metallacyclopropane pyramidalization)
  const meanCx = 0.5 * (c1.x3d + c2.x3d);
  const ptToSign = Math.sign(meanCx - pt.x3d) || 1;
  for (const h of hydrogens) {
    assert.ok(
      (h.x3d - meanCx) * ptToSign > 0.2,
      `Hydrogen ${h.id} must bend back away from Pt (x_H=${h.x3d}, x_C=${meanCx})`
    );
  }

  // Verify Dewar-Chatt-Duncanson Fragment Orbitals
  const fragIds = new Set(zeise.fragmentOrbitals.map((f) => f.id));
  assert.ok(fragIds.has('frag_dcd_m_dpi'), 'Must include Metal d_π Back-Donor fragment orbital');
  assert.ok(fragIds.has('frag_dcd_m_dsigma'), 'Must include Metal d_σ Acceptor fragment orbital');
  assert.ok(fragIds.has('frag_dcd_alkene_pi'), 'Must include Alkene π(C=C) σ-Donor fragment orbital');
  assert.ok(
    fragIds.has('frag_dcd_alkene_pistar'),
    'Must include Alkene π*(C=C) π-Acceptor fragment orbital'
  );

  // Verify Dewar-Chatt-Duncanson Molecular Orbitals (μ/σ-Donation + d_π -> π* Backbonding)
  const sigmaDonMo = zeise.molecularOrbitals.find((m) => m.label.includes('DCD μ/σ-Don'))!;
  assert.ok(sigmaDonMo, 'Must include DCD μ/σ-Donation bonding MO');
  assert.equal(sigmaDonMo.bondingType, 'bonding');
  assert.equal(sigmaDonMo.interactionType, 'σ');
  assert.ok(sigmaDonMo.correlatedFragmentIds.includes('frag_dcd_alkene_pi'));

  const piBackMo = zeise.molecularOrbitals.find((m) => m.label.includes('DCD π-Back'))!;
  assert.ok(piBackMo, 'Must include DCD π-Backbonding MO');
  assert.equal(piBackMo.bondingType, 'bonding');
  assert.equal(piBackMo.interactionType, 'π');
  assert.ok(piBackMo.correlatedFragmentIds.includes('frag_dcd_alkene_pistar'));
  // In the π-backbonding MO, C1 and C2 2px coefficients have opposite signs (π* antibonding character on C=C)
  const c1Term = piBackMo.lcao.find((t) => t.atomId === c1.id && t.orbital === 'px')!;
  const c2Term = piBackMo.lcao.find((t) => t.atomId === c2.id && t.orbital === 'px')!;
  assert.ok(
    c1Term && c2Term && c1Term.coeff * c2Term.coeff < 0,
    'DCD π-backbonding MO must have opposite-sign 2p coefficients on C1 and C2 (alkene π* character)'
  );

  // 2. Dynamic non-canonical DCD complex: Palladium analogue [PdCl3(η²-C2H4)]⁻ and 2D sketcher single-bond completion
  const pdSingleBondAtoms = [
    { id: 'Pd1', element: 'Pd' as const, x2d: 0, y2d: 0, x3d: 0, y3d: 0, z3d: 0, formalCharge: -1 },
    { id: 'Cl1', element: 'Cl' as const, x2d: -1.5, y2d: 0, x3d: -1.5, y3d: 0, z3d: 0 },
    { id: 'Cl2', element: 'Cl' as const, x2d: 0, y2d: -1.5, x3d: 0, y3d: 1.5, z3d: 0 },
    { id: 'Cl3', element: 'Cl' as const, x2d: 0, y2d: 1.5, x3d: 0, y3d: -1.5, z3d: 0 },
    { id: 'C1', element: 'C' as const, x2d: 1.5, y2d: 0, x3d: 1.5, y3d: 0, z3d: 0 },
    { id: 'C2', element: 'C' as const, x2d: 2.2, y2d: 1.0, x3d: 2.2, y3d: -1.0, z3d: 0 },
  ];
  // User only draws a single Pd1-C1 bond to C1=C2 in the 2D sketcher
  const pdSingleBonds = [
    { id: 'b1', source: 'Pd1', target: 'Cl1', order: 1 as const },
    { id: 'b2', source: 'Pd1', target: 'Cl2', order: 1 as const },
    { id: 'b3', source: 'Pd1', target: 'Cl3', order: 1 as const },
    { id: 'b4', source: 'Pd1', target: 'C1', order: 1 as const },
    { id: 'b5', source: 'C1', target: 'C2', order: 2 as const },
  ];
  const pdZeise = analyzeMolecularGraph(pdSingleBondAtoms, pdSingleBonds);
  assert.equal(
    pdZeise.pointGroupId,
    'C2v',
    '[PdCl3(η²-C2H4)]⁻ drawn with a single Pd-C bond must complete η² hapticity and yield C2v symmetry'
  );
  assert.ok(
    pdZeise.bonds.some(
      (b) => (b.source === 'Pd1' && b.target === 'C2') || (b.source === 'C2' && b.target === 'Pd1')
    ),
    'resolveCoordinationBonds must automatically complete the second η² Pd1-C2 bond'
  );
  assert.ok(
    pdZeise.molecularOrbitals.some((m) => m.label.includes('DCD μ/σ-Don')),
    '[PdCl3(η²-C2H4)]⁻ must generate DCD μ/σ-Donation MO'
  );
  assert.ok(
    pdZeise.molecularOrbitals.some((m) => m.label.includes('DCD π-Back')),
    '[PdCl3(η²-C2H4)]⁻ must generate DCD π-Backbonding MO'
  );
});

test('Complexes with larger alkene, chelating diene (Ni(COD)2, Pt(COD)Cl2, Pt(nbd)Cl2), bis-arene (Cr(C6H6)2), piano-stool ((C6H6)Cr(CO)3), and multi-alkene (Ni(C2H4)3) ligands resolve all η²/η⁴/η⁶ bonds, 3D geometries, and point groups', () => {
  // 1. Bis(1,5-cyclooctadiene)nickel(0) Ni(COD)2 -> D2d tetrahedral diolefin chelate with 8 Ni-C bonds
  const niCod2Parsed = parseChemicalInputToGraph('C1CC=CCCC=C1.C1CC=CCCC=C1.[Ni]', 'SMILES')!;
  const niCod2 = analyzeMolecularGraph(niCod2Parsed.atoms, niCod2Parsed.bonds);
  assert.equal(
    niCod2.pointGroupId,
    'D2d',
    'Bis(1,5-cyclooctadiene)nickel(0) Ni(COD)2 must have D2d symmetry'
  );
  const niAtom = niCod2.atoms.find((a) => a.element === 'Ni')!;
  const niBonds = niCod2.bonds.filter((b) => b.source === niAtom.id || b.target === niAtom.id);
  assert.equal(
    niBonds.length,
    8,
    'Ni(COD)2 must have 8 Ni-C coordination bonds (both C=C double bonds on each of the two COD rings)'
  );
  assert.ok(
    niCod2.fragmentOrbitals.some((f) => f.id === 'frag_dcd_m_dpi'),
    'Ni(COD)2 must include Dewar-Chatt-Duncanson fragment orbitals'
  );
  assert.ok(
    niCod2.molecularOrbitals.some((m) => m.label.includes('DCD π-Back')),
    'Ni(COD)2 must include Dewar-Chatt-Duncanson π-Backbonding MO'
  );

  // Verify 3D geometry of Ni(COD)2 has no collapsed atoms and physical C-C bond lengths
  const check3DGeometry = (mol: ReturnType<typeof analyzeMolecularGraph>, label: string) => {
    const heavy = mol.atoms.filter((a) => a.element !== 'H');
    let minHeavyDist = Infinity;
    for (let i = 0; i < heavy.length; i++) {
      for (let j = i + 1; j < heavy.length; j++) {
        const d = Math.hypot(
          heavy[i].x3d - heavy[j].x3d,
          heavy[i].y3d - heavy[j].y3d,
          (heavy[i].z3d || 0) - (heavy[j].z3d || 0)
        );
        if (d < minHeavyDist) minHeavyDist = d;
      }
    }
    let maxCCDist = 0;
    for (const b of mol.bonds) {
      const a1 = mol.atoms.find((a) => a.id === b.source)!;
      const a2 = mol.atoms.find((a) => a.id === b.target)!;
      if (a1.element === 'C' && a2.element === 'C') {
        const d = Math.hypot(
          a1.x3d - a2.x3d,
          a1.y3d - a2.y3d,
          (a1.z3d || 0) - (a2.z3d || 0)
        );
        if (d > maxCCDist) maxCCDist = d;
      }
    }
    assert.ok(
      minHeavyDist >= 1.12,
      `${label}: expected minHeavyDist >= 1.12 Å, got ${minHeavyDist.toFixed(3)} Å`
    );
    assert.ok(
      maxCCDist <= 1.60,
      `${label}: expected max C-C bond length <= 1.60 Å, got ${maxCCDist.toFixed(3)} Å`
    );
  };

  check3DGeometry(niCod2, 'Ni(COD)2');

  // 2. Square-planar chelating diene complexes: Pt(COD)Cl2, Pd(COD)Cl2, Pt(nbd)Cl2 -> C2v with 6 M-L bonds
  for (const [name, smi] of [
    ['Pt(COD)Cl2', 'C1CC=CCCC=C1.Cl[Pt]Cl'],
    ['Pd(COD)Cl2', 'C1CC=CCCC=C1.Cl[Pd]Cl'],
    ['Pt(nbd)Cl2', 'C1C2C=CC1C=C2.Cl[Pt]Cl'],
  ] as const) {
    const parsed = parseChemicalInputToGraph(smi, 'SMILES')!;
    const mol = analyzeMolecularGraph(parsed.atoms, parsed.bonds);
    assert.equal(mol.pointGroupId, 'C2v', `${name} must have C2v symmetry`);
    const m = mol.atoms.find((a) => a.element === 'Pt' || a.element === 'Pd')!;
    const mBonds = mol.bonds.filter((b) => b.source === m.id || b.target === m.id);
    assert.equal(mBonds.length, 6, `${name} must have 6 metal-ligand bonds (4 M-C + 2 M-Cl)`);
    check3DGeometry(mol, name);
  }

  // 3. Bis(benzene)chromium Cr(η⁶-C6H6)2 -> D6h with 12 Cr-C bonds
  const crBz2Parsed = parseChemicalInputToGraph(
    '[CH-]1[CH-][CH-][CH-][CH-][CH-]1.[CH-]1[CH-][CH-][CH-][CH-][CH-]1.[Cr]',
    'SMILES'
  )!;
  const crBz2 = analyzeMolecularGraph(crBz2Parsed.atoms, crBz2Parsed.bonds);
  assert.equal(crBz2.pointGroupId, 'D6h', 'Bis(benzene)chromium Cr(C6H6)2 must have D6h symmetry');
  const crAtom = crBz2.atoms.find((a) => a.element === 'Cr')!;
  const crBonds = crBz2.bonds.filter((b) => b.source === crAtom.id || b.target === crAtom.id);
  assert.equal(crBonds.length, 12, 'Cr(C6H6)2 must have 12 Cr-C coordination bonds');
  check3DGeometry(crBz2, 'Cr(C6H6)2');

  // 4. Piano-stool complexes: (C6H6)Cr(CO)3 -> C3v, (C4H4)Fe(CO)3 -> Cs, (C4H6)Fe(CO)3 -> Cs
  const bzCrCo3Parsed = parseChemicalInputToGraph(
    'c1ccccc1.[C-]#[O+].[C-]#[O+].[C-]#[O+].[Cr]',
    'SMILES'
  )!;
  const bzCrCo3 = analyzeMolecularGraph(bzCrCo3Parsed.atoms, bzCrCo3Parsed.bonds);
  assert.equal(bzCrCo3.pointGroupId, 'C3v', '(C6H6)Cr(CO)3 must have C3v symmetry');
  check3DGeometry(bzCrCo3, '(C6H6)Cr(CO)3');

  const cbdFeCo3Parsed = parseChemicalInputToGraph(
    'C1=CC=C1.[C-]#[O+].[C-]#[O+].[C-]#[O+].[Fe]',
    'SMILES'
  )!;
  const cbdFeCo3 = analyzeMolecularGraph(cbdFeCo3Parsed.atoms, cbdFeCo3Parsed.bonds);
  assert.equal(cbdFeCo3.pointGroupId, 'Cs', '(C4H4)Fe(CO)3 must have Cs symmetry');
  check3DGeometry(cbdFeCo3, '(C4H4)Fe(CO)3');

  const bdFeCo3Parsed = parseChemicalInputToGraph(
    'C=CC=C.[C-]#[O+].[C-]#[O+].[C-]#[O+].[Fe]',
    'SMILES'
  )!;
  const bdFeCo3 = analyzeMolecularGraph(bdFeCo3Parsed.atoms, bdFeCo3Parsed.bonds);
  assert.equal(bdFeCo3.pointGroupId, 'Cs', '(C4H6)Fe(CO)3 must have Cs symmetry');
  check3DGeometry(bdFeCo3, '(C4H6)Fe(CO)3');

  // 5. Tris(ethylene)nickel(0) Ni(C2H4)3 -> D3h with 6 Ni-C bonds
  const niEth3Parsed = parseChemicalInputToGraph('C=C.C=C.C=C.[Ni]', 'SMILES')!;
  const niEth3 = analyzeMolecularGraph(niEth3Parsed.atoms, niEth3Parsed.bonds);
  assert.equal(niEth3.pointGroupId, 'D3h', 'Ni(C2H4)3 must have D3h symmetry');
  check3DGeometry(niEth3, 'Ni(C2H4)3');

  // 6. Larger monodentate η²-alkene complexes (propene, 1-hexene, styrene, cyclohexene, norbornene) -> C2v
  for (const [name, smi] of [
    ['[PtCl3(propene)]-', 'CC=C.Cl[Pt-](Cl)Cl'],
    ['[PtCl3(1-hexene)]-', 'CCCCC=C.Cl[Pt-](Cl)Cl'],
    ['[PtCl3(styrene)]-', 'C=Cc1ccccc1.Cl[Pt-](Cl)Cl'],
    ['[PtCl3(cyclohexene)]-', 'C1CCC=CC1.Cl[Pt-](Cl)Cl'],
    ['[PtCl3(norbornene)]-', 'C1CC2CC1C=C2.Cl[Pt-](Cl)Cl'],
  ] as const) {
    const parsed = parseChemicalInputToGraph(smi, 'SMILES')!;
    const mol = analyzeMolecularGraph(parsed.atoms, parsed.bonds, { charge: -1 });
    assert.equal(mol.pointGroupId, 'C2v', `${name} must have C2v coordination symmetry`);
    const pt = mol.atoms.find((a) => a.element === 'Pt')!;
    const ptBonds = mol.bonds.filter((b) => b.source === pt.id || b.target === pt.id);
    assert.equal(ptBonds.length, 5, `${name} must have 5 Pt-L bonds (3 Pt-Cl + 2 η² Pt-C)`);
    check3DGeometry(mol, name);
  }
});


