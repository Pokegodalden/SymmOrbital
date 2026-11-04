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
    { customName: 'Hydrogen Fluoride', iupacName: 'fluorane', smiles: 'F' }
  );

  const hcl = analyzeMolecularGraph(
    [
      { id: 'Cl1', element: 'Cl', x2d: -0.5, y2d: 0, x3d: 0, y3d: 0, z3d: 0 },
      { id: 'H1', element: 'H', x2d: 0.5, y2d: 0, x3d: 0.6085, y3d: -0.2063, z3d: -1.1433 },
    ],
    [{ id: 'b1', source: 'Cl1', target: 'H1', order: 1 }],
    { customName: 'Hydrogen Chloride', iupacName: 'chlorane', smiles: 'Cl' }
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

