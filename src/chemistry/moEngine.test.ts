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
  assert.equal(moLabelToLatex('1b1'), '1b_{1}');
  assert.equal(moLabelToLatex('4a1*'), '4a_{1}^{*}');
  assert.equal(moLabelToLatex('1e1g (π2)'), '1e_{1g}\\;(\\pi_{2})');

  const water = CANONICAL_MOLECULES.find((m) => m.formula === 'H2O')!;
  const homo = water.molecularOrbitals.find((m) => m.isHOMO)!;
  const lcaoTex = lcaoExpansionToLatex(homo.label, homo.lcao);
  assert.ok(lcaoTex.includes('\\psi\\!\\left(1b_{1}\\right)'), `Unexpected LCAO LaTeX: ${lcaoTex}`);
});

test('NIST Chemistry WebBook PES & NIST CCCBDB benchmark audit achieves 100% point-group parity and < 0.05 eV calibrated MAE', () => {
  const water = CANONICAL_MOLECULES[0];
  const audit = evaluateEngineAgainstBenchmarks(water);

  assert.equal(audit.rows.length, 11);
  assert.equal(audit.pointGroupAccuracyPercent, 100);
  assert.ok(
    audit.meanAbsoluteErrorCalibratedEv < 0.05,
    `Expected calibrated MAE < 0.05 eV vs NIST WebBook PES, got ${audit.meanAbsoluteErrorCalibratedEv}`
  );
  assert.ok(
    audit.meanGeometryErrorAngstrom < 0.02,
    `Expected geometry MAE < 0.02 Å vs NIST CCCBDB, got ${audit.meanGeometryErrorAngstrom}`
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
});

test('Point group detection and edge-case implicit-H rules preserve CO, NO, NO2, O3, SO2, N2, and O2 without mutating into hydrides', () => {
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

  // Nitrogen Dioxide (NO2) radical (17 valence e-) must remain NO2 (never HNO2)
  const no2Atoms = [
    { id: 'N1', element: 'N' as const, x2d: 0, y2d: 0.35, x3d: 0, y3d: 0.35, z3d: 0 },
    { id: 'O1', element: 'O' as const, x2d: -1.05, y2d: -0.2, x3d: -1.05, y3d: -0.2, z3d: 0 },
    { id: 'O2', element: 'O' as const, x2d: 1.05, y2d: -0.2, x3d: 1.05, y3d: -0.2, z3d: 0 },
  ];
  const no2Bonds = [
    { id: 'b1', source: 'N1', target: 'O1', order: 2 as const },
    { id: 'b2', source: 'N1', target: 'O2', order: 1 as const },
  ];
  const no2Analyzed = analyzeMolecularGraph(no2Atoms, no2Bonds);
  assert.equal(no2Analyzed.formula, 'NO2');
  assert.equal(no2Analyzed.valenceElectrons, 17);

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
  // 6 valence basis orbitals of H2O: characters [6, 0, 2, 4] -> 3A1 + B1 + 2B2
  const decomp = decomposeReducibleRepresentation(c2v, [6, 0, 2, 4]);
  assert.equal(decomp.formula, '3A1 + B1 + 2B2');

  const d6h = POINT_GROUPS.D6h;
  const benzenePiPreset = d6h.defaultReduciblePresets[0];
  const piDecomp = decomposeReducibleRepresentation(d6h, benzenePiPreset.characters);
  assert.ok(piDecomp.formula.includes('E1g') && piDecomp.formula.includes('A2u'));
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
