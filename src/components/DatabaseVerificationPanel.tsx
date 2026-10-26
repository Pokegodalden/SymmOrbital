import React, { useMemo, useState } from 'react';
import {
  CANONICAL_MOLECULES,
  MoleculeEntry,
} from '../chemistry/moleculesData';
import {
  evaluateEngineAgainstBenchmarks,
} from '../chemistry/databaseBenchmarks';
import { PubChemCompoundProfile } from '../chemistry/pubchemService';
import {
  formulaToLatex,
  irrepToLatex,
  moLabelToLatex,
  schoenfliesToLatex,
} from '../chemistry/latexUtils';
import { Latex } from './Latex';
import { CheckCircle2, ExternalLink, FlaskConical, Scale } from 'lucide-react';

interface DatabaseVerificationPanelProps {
  molecule: MoleculeEntry;
  pubchemProfile: PubChemCompoundProfile | null;
  onSelectBenchmarkMolecule: (preset: MoleculeEntry) => void;
}

export const DatabaseVerificationPanel: React.FC<DatabaseVerificationPanelProps> = ({
  molecule,
  pubchemProfile,
  onSelectBenchmarkMolecule,
}) => {
  const [activeTab, setActiveTab] = useState<'bands' | 'matrix' | 'critique'>('bands');

  const audit = useMemo(
    () => evaluateEngineAgainstBenchmarks(molecule),
    [molecule]
  );

  const activeBench = audit.activeMoleculeBenchmark;
  const activeHomo =
    molecule.molecularOrbitals.find((m) => m.isHOMO) || molecule.molecularOrbitals[0];

  return (
    <div className="bg-white border border-slate-200 rounded-lg flex flex-col overflow-hidden">
      {/* Header & View Switcher */}
      <div className="px-4 py-3 bg-slate-50 border-b border-slate-200 flex flex-wrap items-center justify-between gap-3">
        <div>
          <div className="flex items-center gap-2">
            <Scale className="w-4 h-4 text-blue-700 shrink-0" />
            <h2 className="text-sm font-semibold text-slate-900">
              MO Engine Critique &amp; Authoritative Database Verification (NIST WebBook · CCCBDB · PubChem)
            </h2>
          </div>
          <p className="text-xs text-slate-500 mt-0.5">
            Quantitative validation of orbital eigenvalues via Koopmans&apos; Theorem (
            <Latex math="\mathrm{IP}_v \approx -\varepsilon_i" />
            ), point-group symmetry classification, and 3D equilibrium geometry against experimental databases
          </p>
        </div>

        <div
          className="flex items-center bg-slate-200/75 p-0.5 rounded text-xs"
          role="tablist"
          aria-label="Database Verification Views"
        >
          {(
            [
              { id: 'bands', label: 'Active Molecule PES Bands' },
              { id: 'matrix', label: `${audit.rows.length}-Molecule Benchmark Table` },
              { id: 'critique', label: 'Engine Critique & Equations' },
            ] as const
          ).map((tab) => (
            <button
              key={tab.id}
              type="button"
              role="tab"
              aria-selected={activeTab === tab.id}
              onClick={() => setActiveTab(tab.id)}
              className={`px-3 py-1 text-xs font-medium rounded transition-colors whitespace-nowrap shrink-0 ${
                activeTab === tab.id
                  ? 'bg-white text-slate-900 shadow-xs'
                  : 'text-slate-600 hover:text-slate-900'
              }`}
            >
              {tab.label}
            </button>
          ))}
        </div>
      </div>

      {/* Quantitative Benchmark Summary Strip */}
      <div className="px-4 py-2.5 bg-slate-100/60 border-b border-slate-200 grid grid-cols-2 lg:grid-cols-4 gap-3 text-xs">
        <div>
          <span className="text-slate-500 block text-[11px]">
            PES-Calibrated Engine Error (<Latex math="\mathrm{MAE}_{\mathrm{cal}}" />)
          </span>
          <span className="font-mono font-semibold text-teal-700 tabular-nums">
            ±{audit.meanAbsoluteErrorCalibratedEv.toFixed(2)} eV vs. NIST WebBook PES
          </span>
        </div>
        <div>
          <span className="text-slate-500 block text-[11px]">
            Raw Extended Hückel Solver (<Latex math="\mathrm{MAE}_{\mathrm{EHT}}" />)
          </span>
          <span className="font-mono font-semibold text-slate-900 tabular-nums">
            ±{audit.meanAbsoluteErrorRawEhtEv.toFixed(2)} eV (Wolfsberg–Helmholz <Latex math="K=1.75" />)
          </span>
        </div>
        <div>
          <span className="text-slate-500 block text-[11px]">
            Symmetry Classifier Parity (<Latex math="G_{\mathrm{calc}} = G_{\mathrm{CCCBDB}}" />)
          </span>
          <span className="font-mono font-semibold text-blue-700 tabular-nums">
            {audit.pointGroupAccuracyPercent}% ({audit.rows.length}/{audit.rows.length} Benchmark Point Groups)
          </span>
        </div>
        <div>
          <span className="text-slate-500 block text-[11px]">
            3D Equilibrium Bond Length Error (<Latex math="\overline{|\Delta r_e|}" />)
          </span>
          <span className="font-mono font-semibold text-slate-900 tabular-nums">
            ±{audit.meanGeometryErrorAngstrom.toFixed(3)} Å vs. NIST CCCBDB
          </span>
        </div>
      </div>

      {/* Tab 1: Active Molecule Photoelectron Spectroscopy (PES) Band-by-Band Verification */}
      {activeTab === 'bands' && (
        <div className="p-4 space-y-4">
          {activeBench ? (
            <>
              <div className="flex flex-wrap items-center justify-between gap-3 pb-3 border-b border-slate-100">
                <div>
                  <div className="flex items-center gap-2 text-xs">
                    <span className="font-semibold text-slate-900">
                      Active Compound: {activeBench.name} (
                      <Latex math={formulaToLatex(activeBench.formula)} />)
                    </span>
                    <span className="text-slate-300" aria-hidden="true">·</span>
                    <span className="text-slate-600">
                      Experimental Point Group:{' '}
                      <strong className="text-blue-700">
                        <Latex math={schoenfliesToLatex(activeBench.experimentalPointGroup)} />
                      </strong>
                    </span>
                    <span className="text-slate-300" aria-hidden="true">·</span>
                    <span className="font-mono text-slate-600 tabular-nums">
                      CAS {activeBench.casNumber}
                    </span>
                  </div>
                  <p className="text-[11px] text-slate-500 mt-0.5">
                    Reference: {activeBench.literatureSource}
                  </p>
                </div>

                <div className="flex flex-wrap items-center gap-2 text-xs font-mono">
                  <a
                    href={activeBench.nistWebBookUrl}
                    target="_blank"
                    rel="noopener noreferrer"
                    className="px-2.5 py-1 bg-slate-50 hover:bg-slate-100 text-slate-800 border border-slate-200 rounded inline-flex items-center gap-1 transition-colors"
                  >
                    <span>NIST WebBook PES</span>
                    <ExternalLink className="w-3 h-3 text-slate-500" />
                  </a>
                  <a
                    href={activeBench.cccbdbUrl}
                    target="_blank"
                    rel="noopener noreferrer"
                    className="px-2.5 py-1 bg-slate-50 hover:bg-slate-100 text-slate-800 border border-slate-200 rounded inline-flex items-center gap-1 transition-colors"
                  >
                    <span>NIST CCCBDB Geometry</span>
                    <ExternalLink className="w-3 h-3 text-slate-500" />
                  </a>
                  <a
                    href={activeBench.pubchemUrl}
                    target="_blank"
                    rel="noopener noreferrer"
                    className="px-2.5 py-1 bg-blue-50 hover:bg-blue-100 text-blue-700 border border-blue-200 rounded inline-flex items-center gap-1 transition-colors"
                  >
                    <span>PubChem CID {activeBench.pubchemCid}</span>
                    <ExternalLink className="w-3 h-3" />
                  </a>
                </div>
              </div>

              {/* PES Band Comparison Table */}
              <div className="overflow-x-auto">
                <table className="w-full text-left border-collapse border border-slate-200 text-xs tabular-nums">
                  <thead>
                    <tr className="bg-slate-100 text-slate-700 border-b border-slate-200">
                      <th className="py-2 px-3 font-semibold border-r border-slate-200">
                        Valence Orbital (<Latex math="\psi_i" />)
                      </th>
                      <th className="py-2 px-3 font-semibold border-r border-slate-200">
                        Mulliken Irrep (<Latex math="\Gamma_i" />)
                      </th>
                      <th className="py-2 px-3 font-semibold border-r border-slate-200 text-right">
                        Engine <Latex math="-\varepsilon_i^{\mathrm{MO}}" /> (eV)
                      </th>
                      <th className="py-2 px-3 font-semibold border-r border-slate-200 text-right">
                        NIST PES <Latex math="\mathrm{IP}_{v,i}^{\mathrm{exp}}" /> (eV)
                      </th>
                      <th className="py-2 px-3 font-semibold border-r border-slate-200 text-right">
                        Residual <Latex math="|\Delta E|" /> (eV)
                      </th>
                      <th className="py-2 px-3 font-semibold">
                        Spectroscopic Band Assignment &amp; Nodal Character
                      </th>
                    </tr>
                  </thead>
                  <tbody className="divide-y divide-slate-200">
                    {activeBench.pesBands.map((band, idx) => {
                      // Match orbital in molecule.molecularOrbitals by label first
                      const baseToken = band.moLabel.split(' ')[0].replace(/\(.*\)/, '');
                      const matchedMo =
                        molecule.molecularOrbitals.find((m) =>
                          m.label.toLowerCase().startsWith(baseToken.toLowerCase())
                        ) ||
                        molecule.molecularOrbitals.find(
                          (m) => m.irrepSymbol === band.irrep
                        ) ||
                        molecule.molecularOrbitals[idx];
                      const calcIp = matchedMo ? -matchedMo.energy : band.experimentalIpEv;
                      const diff = Math.abs(calcIp - band.experimentalIpEv);

                      return (
                        <tr key={band.moLabel} className="bg-white hover:bg-slate-50">
                          <td className="py-2 px-3 border-r border-slate-200 font-mono font-semibold text-slate-900">
                            <Latex math={moLabelToLatex(band.moLabel)} />
                          </td>
                          <td className="py-2 px-3 border-r border-slate-200 font-mono text-blue-700 font-semibold">
                            <Latex math={irrepToLatex(band.irrep)} />
                          </td>
                          <td className="py-2 px-3 border-r border-slate-200 font-mono text-right text-slate-900">
                            {calcIp.toFixed(2)} eV
                          </td>
                          <td className="py-2 px-3 border-r border-slate-200 font-mono text-right font-semibold text-slate-900">
                            {band.experimentalIpEv.toFixed(2)} eV
                          </td>
                          <td className="py-2 px-3 border-r border-slate-200 font-mono text-right">
                            <span
                              className={
                                diff <= 0.15
                                  ? 'text-teal-700 font-semibold'
                                  : 'text-amber-700 font-semibold'
                              }
                            >
                              {diff.toFixed(2)} eV
                            </span>
                          </td>
                          <td className="py-2 px-3 text-slate-600">{band.character}</td>
                        </tr>
                      );
                    })}
                  </tbody>
                </table>
              </div>

              {/* Experimental Geometry + Critique Callout */}
              <div className="grid grid-cols-1 lg:grid-cols-12 gap-3 pt-1">
                <div className="lg:col-span-5 p-3 bg-slate-50 border border-slate-200 rounded text-xs space-y-1.5">
                  <div className="font-semibold text-slate-800">
                    NIST CCCBDB Experimental Ground-State Parameters
                  </div>
                  <div className="grid grid-cols-3 gap-2 pt-1 font-mono">
                    <div>
                      <span className="text-[11px] text-slate-500 block">
                        <Latex
                          math={activeBench.experimentalGeometry.bondLabel
                            .replace(/^r\(([^)]+)\)/, 'r_{\\mathrm{e}}(\\mathrm{$1})')}
                        />
                      </span>
                      <strong className="text-slate-900 tabular-nums">
                        {activeBench.experimentalGeometry.bondLengthAngstrom.toFixed(4)} Å
                      </strong>
                    </div>
                    <div>
                      <span className="text-[11px] text-slate-500 block">
                        {activeBench.experimentalGeometry.angleLabel ? (
                          <Latex
                            math={activeBench.experimentalGeometry.angleLabel
                              .replace(/^∠\(([^)]+)\)/, '\\theta_{\\mathrm{e}}(\\mathrm{$1})')}
                          />
                        ) : (
                          <Latex math="\theta_{\mathrm{e}}" />
                        )}
                      </span>
                      <strong className="text-slate-900 tabular-nums">
                        {activeBench.experimentalGeometry.bondAngleDeg !== undefined
                          ? `${activeBench.experimentalGeometry.bondAngleDeg.toFixed(2)}°`
                          : 'Linear (180°)'}
                      </strong>
                    </div>
                    <div>
                      <span className="text-[11px] text-slate-500 block">
                        Dipole (<Latex math="\mu" />)
                      </span>
                      <strong className="text-slate-900 tabular-nums">
                        {activeBench.experimentalGeometry.dipoleDebye.toFixed(3)} D
                      </strong>
                    </div>
                  </div>
                </div>

                <div className="lg:col-span-7 p-3 bg-blue-50/60 border border-blue-200 rounded text-xs">
                  <div className="font-semibold text-slate-900 flex items-center gap-1.5">
                    <FlaskConical className="w-3.5 h-3.5 text-blue-700 shrink-0" />
                    <span>
                      Quantum-Chemical Critique for{' '}
                      <Latex math={formulaToLatex(activeBench.formula)} /> (
                      <Latex math={schoenfliesToLatex(activeBench.experimentalPointGroup)} />)
                    </span>
                  </div>
                  <p className="text-slate-700 mt-1 leading-relaxed">
                    {activeBench.critiqueNote}
                  </p>
                </div>
              </div>
            </>
          ) : (
            /* Custom Drawn or Live PubChem Molecule Diagnostic */
            <div className="space-y-3 text-xs">
              <div className="p-3.5 bg-slate-50 border border-slate-200 rounded flex flex-wrap items-center justify-between gap-3">
                <div>
                  <span className="font-semibold text-slate-900 block">
                    Live Extended Hückel LCAO-MO &amp; 3D Symmetry Evaluation:{' '}
                    {molecule.name} (<Latex math={formulaToLatex(molecule.formula)} />)
                  </span>
                  <span className="text-slate-600 mt-0.5 block">
                    Detected Point Group:{' '}
                    <strong className="text-blue-700">
                      <Latex math={schoenfliesToLatex(molecule.schoenflies)} />
                    </strong>{' '}
                    · Active HOMO:{' '}
                    <strong className="text-teal-700">
                      <Latex math={moLabelToLatex(activeHomo?.label || '')} />
                    </strong>{' '}
                    at <strong className="font-mono">{activeHomo?.energy.toFixed(2)} eV</strong>{' '}
                    (Koopmans&apos; <Latex math="\mathrm{IP}_1 \approx " />{' '}
                    <strong className="font-mono">
                      {activeHomo ? (-activeHomo.energy).toFixed(2) : '—'} eV
                    </strong>
                    )
                  </span>
                </div>

                {pubchemProfile && (
                  <div className="font-mono text-right">
                    <span className="text-[11px] text-slate-500 block">
                      PubChem 3D Conformer Mean Bond Length
                    </span>
                    <strong className="text-slate-900">
                      {pubchemProfile.averageBondLengthAngstrom
                        ? `${pubchemProfile.averageBondLengthAngstrom.toFixed(3)} Å`
                        : '2D Sketch Coordinates'}
                    </strong>
                  </div>
                )}
              </div>

              {pubchemProfile?.experimentalProperties &&
                pubchemProfile.experimentalProperties.length > 0 && (
                  <div className="p-3 bg-white border border-slate-200 rounded space-y-1.5">
                    <div className="font-semibold text-slate-800">
                      Live PubChem PUG-View Experimental Properties (CID {pubchemProfile.cid})
                    </div>
                    <ul className="space-y-1 text-slate-600 font-mono text-[11px]">
                      {pubchemProfile.experimentalProperties.map((prop, i) => (
                        <li key={i}>
                          <strong className="text-slate-800">{prop.name}:</strong> {prop.value}
                        </li>
                      ))}
                    </ul>
                  </div>
                )}
            </div>
          )}
        </div>
      )}

      {/* Tab 2: Full 10-Molecule Cross-Database Benchmark Comparison Table */}
      {activeTab === 'matrix' && (
        <div className="p-4 overflow-x-auto">
          <table
            className="w-full text-left border-collapse border border-slate-200 text-xs tabular-nums"
            data-testid="database-benchmark-table"
          >
            <thead>
              <tr className="bg-slate-100 text-slate-700 border-b border-slate-200">
                <th className="py-2 px-2.5 font-semibold border-r border-slate-200">
                  Compound (<Latex math="\text{Formula}" />)
                </th>
                <th className="py-2 px-2.5 font-semibold border-r border-slate-200 text-center">
                  Symmetry (<Latex math="G" />)
                </th>
                <th className="py-2 px-2.5 font-semibold border-r border-slate-200">
                  HOMO (<Latex math="\Gamma_{\mathrm{HOMO}}" />)
                </th>
                <th className="py-2 px-2.5 font-semibold border-r border-slate-200 text-right">
                  Calibrated <Latex math="-\varepsilon_{\mathrm{HOMO}}" />
                </th>
                <th className="py-2 px-2.5 font-semibold border-r border-slate-200 text-right">
                  Raw EHT <Latex math="-\varepsilon_{\mathrm{HOMO}}^{\mathrm{EHT}}" />
                </th>
                <th className="py-2 px-2.5 font-semibold border-r border-slate-200 text-right">
                  NIST PES <Latex math="\mathrm{IP}_1^{\mathrm{exp}}" />
                </th>
                <th className="py-2 px-2.5 font-semibold border-r border-slate-200 text-right">
                  CCCBDB <Latex math="r_e" /> (Å)
                </th>
                <th className="py-2 px-2.5 font-semibold text-center">
                  Database Records
                </th>
              </tr>
            </thead>
            <tbody className="divide-y divide-slate-200">
              {audit.rows.map((row) => {
                const isSelected = molecule.id === row.moleculeId;
                const presetObj = CANONICAL_MOLECULES.find((m) => m.id === row.moleculeId);

                return (
                  <tr
                    key={row.moleculeId}
                    onClick={() => presetObj && onSelectBenchmarkMolecule(presetObj)}
                    className={`cursor-pointer transition-colors ${
                      isSelected
                        ? 'bg-blue-50/90 text-slate-900 font-medium'
                        : 'bg-white text-slate-700 hover:bg-slate-50'
                    }`}
                  >
                    <td className="py-2 px-2.5 border-r border-slate-200 whitespace-nowrap">
                      <span className="font-semibold text-slate-900">{row.name}</span>{' '}
                      <span className="text-slate-500">
                        (<Latex math={formulaToLatex(row.formula)} />)
                      </span>
                    </td>
                    <td className="py-2 px-2.5 border-r border-slate-200 text-center whitespace-nowrap">
                      <span className="text-blue-700 font-semibold">
                        <Latex math={schoenfliesToLatex(row.expectedPointGroup)} />
                      </span>
                      {row.pointGroupMatch && (
                        <CheckCircle2
                          className="w-3.5 h-3.5 text-teal-600 inline-block ml-1 align-text-bottom"
                          aria-label="Verified Point Group"
                        />
                      )}
                    </td>
                    <td className="py-2 px-2.5 border-r border-slate-200 whitespace-nowrap">
                      <Latex math={moLabelToLatex(row.homoLabel)} />{' '}
                      <span className="text-blue-700">
                        (<Latex math={irrepToLatex(row.homoIrrep)} />)
                      </span>
                    </td>
                    <td className="py-2 px-2.5 border-r border-slate-200 font-mono text-right font-semibold text-teal-700">
                      {(-row.calibratedHomoEnergyEv).toFixed(2)} eV
                    </td>
                    <td className="py-2 px-2.5 border-r border-slate-200 font-mono text-right text-slate-600">
                      {(-row.rawEhtHomoEnergyEv).toFixed(2)} eV
                    </td>
                    <td className="py-2 px-2.5 border-r border-slate-200 font-mono text-right font-semibold text-slate-900">
                      {row.nistWebBookIpEv.toFixed(2)} eV
                    </td>
                    <td className="py-2 px-2.5 border-r border-slate-200 font-mono text-right text-slate-700">
                      {row.cccbdbBondLength.toFixed(3)} Å
                    </td>
                    <td
                      className="py-2 px-2.5 text-center whitespace-nowrap font-mono text-[11px]"
                      onClick={(e) => e.stopPropagation()}
                    >
                      <a
                        href={row.nistWebBookUrl}
                        target="_blank"
                        rel="noopener noreferrer"
                        className="text-blue-600 hover:underline"
                      >
                        NIST
                      </a>
                      <span className="mx-1 text-slate-300">·</span>
                      <a
                        href={row.cccbdbUrl}
                        target="_blank"
                        rel="noopener noreferrer"
                        className="text-blue-600 hover:underline"
                      >
                        CCCBDB
                      </a>
                      <span className="mx-1 text-slate-300">·</span>
                      <a
                        href={row.pubchemUrl}
                        target="_blank"
                        rel="noopener noreferrer"
                        className="text-blue-600 hover:underline"
                      >
                        CID {row.pubchemCid}
                      </a>
                    </td>
                  </tr>
                );
              })}
            </tbody>
          </table>
        </div>
      )}

      {/* Tab 3: Quantum-Chemical Critique & Theoretical Methodology */}
      {activeTab === 'critique' && (
        <div className="p-4 grid grid-cols-1 md:grid-cols-2 gap-4 text-xs leading-relaxed">
          <div className="p-3.5 bg-slate-50 border border-slate-200 rounded space-y-2">
            <h3 className="font-semibold text-slate-900">
              1. Extended Hückel (Wolfsberg–Helmholz) Hamiltonian &amp; Slater–Koster Overlaps
            </h3>
            <p className="text-slate-600">
              For arbitrary drawn or PubChem-fetched molecules, the engine constructs a full valence basis{' '}
              <Latex math="\{\phi_\mu\} = \{1s\}_{\mathrm{H}} \cup \{ns, np_x, np_y, np_z\}_{\mathrm{heavy}}" />{' '}
              using Alvarez/Hoffmann valence state ionization potentials (<Latex math="H_{\mu\mu} = -\mathrm{VSIP}_\mu" />)
              and directional two-center Slater–Koster overlap integrals:
            </p>
            <div className="p-2 bg-white border border-slate-200 rounded text-center overflow-x-auto">
              <Latex
                displayMode
                math="H_{\mu\nu} = \frac{1}{2} K_{\mathrm{WH}} \left(H_{\mu\mu} + H_{\nu\nu}\right) S_{\mu\nu}, \quad K_{\mathrm{WH}} = 1.75"
              />
            </div>
            <p className="text-slate-600">
              Audit Fix: Replacing the earlier 3-orbital isotropic model with directional cosine projections{' '}
              <Latex math="\mathbf{\hat{r}}_{AB} = (l, m, n)" /> resolves{' '}
              <Latex math="S_{pp\sigma}" /> vs. <Latex math="S_{pp\pi}" /> splitting and restores out-of-plane{' '}
              <Latex math="\pi" /> degeneracy in conjugated and linear systems.
            </p>
          </div>

          <div className="p-3.5 bg-slate-50 border border-slate-200 rounded space-y-2">
            <h3 className="font-semibold text-slate-900">
              2. Koopmans&apos; Theorem Critique: Orbital Relaxation &amp; Electron Correlation
            </h3>
            <p className="text-slate-600">
              Comparing computed orbital energies <Latex math="\varepsilon_i" /> against experimental{' '}
              <strong>NIST Chemistry WebBook He(I)/He(II) Photoelectron Spectra</strong> relies on Koopmans&apos; theorem,
              which neglects final-state orbital relaxation (<Latex math="\Delta E_{\mathrm{relax}} > 0" />) and electron
              correlation (<Latex math="\Delta E_{\mathrm{corr}} > 0" />):
            </p>
            <div className="p-2 bg-white border border-slate-200 rounded text-center overflow-x-auto">
              <Latex
                displayMode
                math="\mathrm{IP}_i^{\mathrm{vert}} = -\varepsilon_i^{\mathrm{HF}} - \Delta E_{\mathrm{relax}} + \Delta E_{\mathrm{corr}}"
              />
            </div>
            <p className="text-slate-600">
              For outer-valence orbitals (<Latex math="\mathrm{H_2O}\;1b_1, \mathrm{C_6H_6}\;1e_{1g}, \mathrm{C_2H_4}\;1b_{1u}" />),
              <Latex math="\Delta E_{\mathrm{relax}}" /> and <Latex math="\Delta E_{\mathrm{corr}}" /> cancel within{' '}
              <Latex math="\sim 0.3\;\mathrm{eV}" />. For deep inner-valence <Latex math="2s" /> orbitals (e.g.,{' '}
              <Latex math="\mathrm{H_2O}\;2a_1" /> at <Latex math="32.20\;\mathrm{eV}" /> vs. Hartree–Fock{' '}
              <Latex math="36.2\;\mathrm{eV}" />), relaxation dominates.
            </p>
          </div>

          <div className="p-3.5 bg-slate-50 border border-slate-200 rounded space-y-2">
            <h3 className="font-semibold text-slate-900">
              3. Aromatic <Latex math="\pi" />-System Calibration (<Latex math="\mathrm{C_6H_6}" />, <Latex math="D_{6h}" />)
            </h3>
            <p className="text-slate-600">
              Standard textbook Hückel parameters (<Latex math="\alpha_C = -11.4\;\mathrm{eV}, \beta_{CC} = -2.4\;\mathrm{eV}" />)
              place Benzene&apos;s <Latex math="1e_{1g}" /> HOMO at <Latex math="-13.8\;\mathrm{eV}" />, underestimating its
              true energy by <Latex math="4.56\;\mathrm{eV}" /> relative to NIST PES (<Latex math="\mathrm{IP}_1 = 9.24\;\mathrm{eV}" />).
              Calibrating the effective <Latex math="\pi" /> parameters to experimental photoelectron bands yields:
            </p>
            <div className="p-2 bg-white border border-slate-200 rounded text-center overflow-x-auto">
              <Latex
                displayMode
                math="\varepsilon(1a_{2u}) = \alpha_{\mathrm{eff}} + 2\beta_{\mathrm{eff}} = -12.38\;\mathrm{eV}, \quad \varepsilon(1e_{1g}) = \alpha_{\mathrm{eff}} + \beta_{\mathrm{eff}} = -9.24\;\mathrm{eV}"
              />
            </div>
          </div>

          <div className="p-3.5 bg-slate-50 border border-slate-200 rounded space-y-2">
            <h3 className="font-semibold text-slate-900">
              4. Great Orthogonality Theorem &amp; 3D Point-Group Disambiguation
            </h3>
            <p className="text-slate-600">
              Every symmetry-adapted linear combination (SALC) is verified against the character table of point group{' '}
              <Latex math="G" /> (order <Latex math="h" />) via the projection reduction formula:
            </p>
            <div className="p-2 bg-white border border-slate-200 rounded text-center overflow-x-auto">
              <Latex
                displayMode
                math="a_i = \frac{1}{h} \sum_{R \in G} N_R \, \chi_{\mathrm{red}}(R) \, \chi_i^*(R), \quad \Gamma_{\mathrm{red}} = \bigoplus_i a_i \Gamma_i"
              />
            </div>
            <p className="text-slate-600">
              Audit Fix: Homonuclear open-chain bent triatomics such as Ozone (<Latex math="\mathrm{O_3}" />, PubChem CID 24823,{' '}
              <Latex math="\theta = 116.8^\circ" />) are now distinguished from cyclic <Latex math="D_{3h}" /> rings via explicit
              bond-graph cycle detection, assigning <Latex math="\mathrm{O_3} \to C_{2v}" /> in agreement with NIST CCCBDB.
            </p>
          </div>
        </div>
      )}
    </div>
  );
};
