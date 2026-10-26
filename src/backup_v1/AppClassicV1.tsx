/**
 * @license
 * SPDX-License-Identifier: Apache-2.0
 *
 * BACKUP SNAPSHOT (V1 Classic Stacked Layout)
 * Saved prior to the V2 Pro 3-Column Synchronized Scientific Workbench upgrade.
 * You can switch to this exact Classic V1 view at any time using the "Classic Layout (v1)"
 * button in the top navigation bar, or restore it as the default by exporting AppClassicV1
 * from /src/App.tsx.
 */

import React, { useMemo, useState } from 'react';
import {
  AtomNode,
  BondEdge,
  CANONICAL_MOLECULES,
  MolecularOrbital,
  MoleculeEntry,
} from '../chemistry/moleculesData';
import { POINT_GROUPS } from '../chemistry/symmetryData';
import { analyzeMolecularGraph } from '../chemistry/moEngine';
import {
  PubChemCompoundProfile,
  lookupPubChemCompound,
} from '../chemistry/pubchemService';
import {
  formulaToLatex,
  irrepToLatex,
  moLabelToLatex,
  schoenfliesToLatex,
} from '../chemistry/latexUtils';
import { Latex } from '../components/Latex';
import { ChemicalSketcher2D } from '../components/ChemicalSketcher2D';
import { MODiagramPanel } from '../components/MODiagramPanel';
import { OrbitalVisualizer3D } from '../components/OrbitalVisualizer3D';
import { SymmetryGroupPanel } from '../components/SymmetryGroupPanel';
import { DatabaseVerificationPanel } from '../components/DatabaseVerificationPanel';
import { ExternalLink, Layers, Moon, Search, Sun } from 'lucide-react';

interface AppClassicV1Props {
  onSwitchToProWorkbench?: () => void;
  darkMode?: boolean;
  onToggleDarkMode?: () => void;
}

export function AppClassicV1({
  onSwitchToProWorkbench,
  darkMode = false,
  onToggleDarkMode,
}: AppClassicV1Props) {
  const defaultWater = CANONICAL_MOLECULES[0];

  const [atoms, setAtoms] = useState<AtomNode[]>(defaultWater.atoms);
  const [bonds, setBonds] = useState<BondEdge[]>(defaultWater.bonds);
  const [customTitle, setCustomTitle] = useState<string | undefined>(undefined);
  const [customCid, setCustomCid] = useState<number | undefined>(undefined);
  const [charge, setCharge] = useState<number>(0);
  const [excitedState, setExcitedState] = useState<boolean>(false);
  const [pointGroupOverride, setPointGroupOverride] = useState<string | null>(null);
  const [selectedMoId, setSelectedMoId] = useState<string>('mo_1b1');
  const [highlightedSymmetryElementId, setHighlightedSymmetryElementId] = useState<string | null>(null);

  // Live PubChem API Lookup state
  const [pubchemQuery, setPubchemQuery] = useState<string>('');
  const [pubchemProfile, setPubchemProfile] = useState<PubChemCompoundProfile | null>(null);
  const [pubchemLoading, setPubchemLoading] = useState<boolean>(false);
  const [pubchemError, setPubchemError] = useState<string | null>(null);

  // Compute live MoleculeEntry (matches canonical benchmark or runs Extended Hückel LCAO-MO + Symmetry solver)
  const molecule: MoleculeEntry = useMemo(() => {
    return analyzeMolecularGraph(atoms, bonds, {
      customName: customTitle,
      customCid,
      charge,
      excitedState,
    });
  }, [atoms, bonds, customTitle, customCid, charge, excitedState]);

  const activePointGroup = useMemo(() => {
    const pgKey = pointGroupOverride || molecule.pointGroupId;
    return POINT_GROUPS[pgKey] || POINT_GROUPS.C2v;
  }, [molecule.pointGroupId, pointGroupOverride]);

  const selectedMo: MolecularOrbital = useMemo(() => {
    const found = molecule.molecularOrbitals.find((m) => m.id === selectedMoId);
    if (found) return found;
    const homo = molecule.molecularOrbitals.find((m) => m.isHOMO);
    return homo || molecule.molecularOrbitals[0];
  }, [molecule, selectedMoId]);

  const handleSelectPresetMolecule = (preset: MoleculeEntry) => {
    setAtoms(preset.atoms);
    setBonds(preset.bonds);
    setCustomTitle(undefined);
    setCustomCid(preset.pubchemCid);
    setCharge(0);
    setExcitedState(false);
    setPointGroupOverride(null);
    setHighlightedSymmetryElementId(null);
    setPubchemError(null);
    const homo = preset.molecularOrbitals.find((m) => m.isHOMO) || preset.molecularOrbitals[0];
    if (homo) setSelectedMoId(homo.id);
  };

  const handleGraphChange = (
    nextAtoms: AtomNode[],
    nextBonds: BondEdge[],
    nextTitle?: string
  ) => {
    setAtoms(nextAtoms);
    setBonds(nextBonds);
    setCustomTitle(nextTitle);
    setPointGroupOverride(null);
  };

  const handlePubChemLookup = async (queryText: string) => {
    const q = queryText.trim();
    if (!q) return;
    setPubchemLoading(true);
    setPubchemError(null);

    try {
      const lower = q.toLowerCase();
      const canonical = CANONICAL_MOLECULES.find(
        (m) =>
          m.name.toLowerCase() === lower ||
          m.iupacName.toLowerCase() === lower ||
          m.id.toLowerCase() === lower ||
          String(m.pubchemCid) === lower
      );

      const profile = await lookupPubChemCompound(q);
      setPubchemProfile(profile);

      if (canonical) {
        handleSelectPresetMolecule(canonical);
      } else if (profile.atoms.length > 0) {
        setAtoms(profile.atoms);
        setBonds(profile.bonds);
        setCustomTitle(profile.title || profile.iupacName);
        setCustomCid(profile.cid);
        setCharge(profile.charge || 0);
        setExcitedState(false);
        setPointGroupOverride(null);
      }
    } catch (err: any) {
      setPubchemError(err?.message || 'Unable to resolve compound from PubChem.');
    } finally {
      setPubchemLoading(false);
    }
  };

  return (
    <div className="min-h-screen flex flex-col bg-[#F8FAFC] text-slate-900">
      {/* Strict 3-Zone Top Navigation Bar Contract */}
      <header className="flex items-center justify-between gap-8 px-6 py-3.5 bg-white border-b border-slate-200 sticky top-0 z-30">
        <a
          href="#top"
          className="text-lg font-bold tracking-tight text-slate-900 whitespace-nowrap shrink-0"
        >
          SymmOrbital
        </a>

        <nav
          className="hidden md:flex items-center gap-6 text-sm font-medium text-slate-600"
          aria-label="Primary Workspace Navigation"
        >
          <a
            href="#sketcher-3d"
            className="hover:text-slate-900 hover:underline underline-offset-4 transition-colors whitespace-nowrap shrink-0"
          >
            Structure &amp; 3D
          </a>
          <a
            href="#mo-diagram"
            className="hover:text-slate-900 hover:underline underline-offset-4 transition-colors whitespace-nowrap shrink-0"
          >
            MO Diagram
          </a>
          <a
            href="#symmetry-irreps"
            className="hover:text-slate-900 hover:underline underline-offset-4 transition-colors whitespace-nowrap shrink-0"
          >
            Symmetry &amp; Irreps
          </a>
          <a
            href="#database-audit"
            className="hover:text-slate-900 hover:underline underline-offset-4 transition-colors whitespace-nowrap shrink-0"
          >
            Database Audit
          </a>
        </nav>

        <div className="flex items-center gap-2.5 shrink-0">
          {onToggleDarkMode && (
            <button
              type="button"
              onClick={onToggleDarkMode}
              title={darkMode ? 'Switch to Light Mode' : 'Switch to Dark Mode'}
              aria-label={darkMode ? 'Switch to Light Mode' : 'Switch to Dark Mode'}
              className="p-2 text-xs font-medium text-slate-700 bg-slate-100 border border-slate-300 rounded-lg hover:bg-slate-200 transition-colors flex items-center justify-center shrink-0"
            >
              {darkMode ? (
                <Sun className="w-4 h-4 text-amber-400" />
              ) : (
                <Moon className="w-4 h-4 text-slate-600" />
              )}
            </button>
          )}
          {onSwitchToProWorkbench && (
            <button
              type="button"
              onClick={onSwitchToProWorkbench}
              className="px-3 py-2 text-xs font-medium text-slate-700 bg-slate-100 border border-slate-300 rounded-lg hover:bg-slate-200 transition-colors flex items-center gap-1.5 whitespace-nowrap shrink-0"
              title="Switch to the Pro 3-Pane Synchronized Workbench Layout"
            >
              <Layers className="w-3.5 h-3.5 text-blue-600" />
              <span>Switch to Pro Workbench (v2)</span>
            </button>
          )}
          <button
            type="button"
            onClick={() => handleSelectPresetMolecule(defaultWater)}
            className="px-4 py-2 text-xs font-medium text-white bg-blue-600 rounded-lg hover:bg-blue-700 transition-colors whitespace-nowrap shrink-0"
          >
            Reset to Water (<Latex math="\mathrm{H_2O}" />)
          </button>
        </div>
      </header>

      {/* Main Scientific Workbench Container */}
      <main id="top" className="flex-1 max-w-[1440px] w-full mx-auto px-4 sm:px-6 py-5 space-y-5">
        {/* Quick-Try Canonical Molecule & Symmetry Group Bar */}
        <section
          aria-label="Canonical Molecules and Active Telemetry"
          className="bg-white border border-slate-200 rounded-lg p-3.5 flex flex-col gap-3"
        >
          <div className="flex flex-wrap items-center justify-between gap-3">
            <div>
              <h1 className="text-base font-bold text-slate-900">
                Interactive Chemical Sketcher, Molecular Orbital Diagram &amp; Point-Group Symmetry Workbench
              </h1>
              <p className="text-xs text-slate-500 mt-0.5">
                Draw any 2D chemical structure or select a benchmark molecule below to compute its LCAO-MO energy levels{' '}
                <Latex math="\varepsilon_i" />, 3D orbital wave functions{' '}
                <Latex math="\psi_i(\mathbf{r})" />, and irreducible representations{' '}
                <Latex math="\Gamma_i" />
              </p>
            </div>

            {/* Active Molecule Summary Metadata */}
            <div className="flex flex-wrap items-center gap-2 text-xs font-mono text-slate-700">
              <span>
                Compound:{' '}
                <strong className="text-slate-900">
                  {molecule.name} (<Latex math={formulaToLatex(molecule.formula)} />)
                </strong>
              </span>
              <span className="text-slate-300" aria-hidden="true">·</span>
              <span>
                Point Group:{' '}
                <strong className="text-blue-700">
                  <Latex math={schoenfliesToLatex(activePointGroup.schoenflies)} />
                </strong>
              </span>
              <span className="text-slate-300" aria-hidden="true">·</span>
              <span>
                Active MO:{' '}
                <strong className="text-teal-700">
                  {selectedMo ? (
                    <>
                      <Latex math={moLabelToLatex(selectedMo.label)} /> (
                      <Latex math={irrepToLatex(selectedMo.irrepSymbol)} />)
                    </>
                  ) : (
                    '—'
                  )}
                </strong>
              </span>
              {molecule.pubchemCid > 0 && (
                <>
                  <span className="text-slate-300" aria-hidden="true">·</span>
                  <a
                    href={`https://pubchem.ncbi.nlm.nih.gov/compound/${molecule.pubchemCid}`}
                    target="_blank"
                    rel="noopener noreferrer"
                    className="text-blue-600 hover:underline inline-flex items-center gap-1"
                  >
                    <span>CID {molecule.pubchemCid}</span>
                    <ExternalLink className="w-3 h-3" />
                  </a>
                </>
              )}
            </div>
          </div>

          {/* Preset Molecule Selector Buttons */}
          <div
            className="flex flex-wrap items-center gap-1.5 pt-2 border-t border-slate-100"
            role="group"
            aria-label="Preset Benchmark Molecules"
          >
            <span className="text-xs font-medium text-slate-500 mr-1 whitespace-nowrap">
              Benchmark Molecules:
            </span>
            {CANONICAL_MOLECULES.map((preset) => {
              const isCurrent = molecule.formula === preset.formula && !customTitle;
              return (
                <button
                  key={preset.id}
                  type="button"
                  onClick={() => handleSelectPresetMolecule(preset)}
                  className={`px-3 py-1.5 text-xs font-medium rounded-md border transition-colors whitespace-nowrap shrink-0 ${
                    isCurrent
                      ? 'bg-slate-900 text-white border-slate-900 shadow-xs'
                      : 'bg-slate-50 text-slate-700 border-slate-200 hover:bg-slate-100 hover:text-slate-900'
                  }`}
                >
                  <Latex math={formulaToLatex(preset.formula)} /> ({preset.name}) ·{' '}
                  <Latex math={schoenfliesToLatex(preset.schoenflies)} />
                </button>
              );
            })}
          </div>
        </section>

        {/* Top Dual-Pane Stage: 2D Chemical Structure Sketcher (Left) + 3D Orbital & Symmetry Visualizer (Right) */}
        <section
          id="sketcher-3d"
          aria-label="2D Structure Sketcher and 3D Orbital Visualizer"
          className="grid grid-cols-1 lg:grid-cols-12 gap-5 items-stretch"
        >
          <div className="lg:col-span-6 flex flex-col">
            <ChemicalSketcher2D
              atoms={atoms}
              bonds={bonds}
              onChangeGraph={handleGraphChange}
              onTriggerPubChemSearch={handlePubChemLookup}
              isSearchingPubChem={pubchemLoading}
            />
          </div>

          <div className="lg:col-span-6 flex flex-col">
            <OrbitalVisualizer3D
              molecule={molecule}
              selectedMo={selectedMo}
              pointGroup={activePointGroup}
              highlightedSymmetryElementId={highlightedSymmetryElementId}
            />
          </div>
        </section>

        {/* Middle Stage: Interactive MO Energy Level & Correlation Diagram */}
        <section id="mo-diagram" aria-label="Interactive Molecular Orbital Diagram">
          <MODiagramPanel
            molecule={molecule}
            selectedMoId={selectedMo?.id || ''}
            onSelectMo={(mo) => setSelectedMoId(mo.id)}
            charge={charge}
            onChangeCharge={setCharge}
            excitedState={excitedState}
            onToggleExcitedState={() => setExcitedState((prev) => !prev)}
          />
        </section>

        {/* Point-Group Character Table & Irreducible Representation Decomposer */}
        <section id="symmetry-irreps" aria-label="Symmetry Point Group and Irreducible Representations">
          <SymmetryGroupPanel
            molecule={molecule}
            pointGroup={activePointGroup}
            selectedMo={selectedMo}
            highlightedSymmetryElementId={highlightedSymmetryElementId}
            onHighlightSymmetryElement={setHighlightedSymmetryElementId}
            onSelectPointGroupOverride={(pgId) => setPointGroupOverride(pgId)}
          />
        </section>

        {/* MO Engine Critique & Authoritative Database Verification */}
        <section
          id="database-audit"
          aria-label="MO Engine Critique and Database Benchmark Verification"
        >
          <DatabaseVerificationPanel
            molecule={molecule}
            pubchemProfile={pubchemProfile}
            onSelectBenchmarkMolecule={handleSelectPresetMolecule}
          />
        </section>

        {/* PubChem PUG REST Live Structure & Property Inspector */}
        <section
          id="pubchem-inspector"
          aria-label="PubChem Database Structure Lookup"
          className="bg-white border border-slate-200 rounded-lg p-4 space-y-3"
        >
          <div className="flex flex-wrap items-center justify-between gap-3">
            <div>
              <h2 className="text-sm font-semibold text-slate-900">
                PubChem PUG REST &amp; PUG-View Compound &amp; 3D Conformer Lookup
              </h2>
              <p className="text-xs text-slate-500">
                Fetch live 3D SDF atomic coordinates, IUPAC identifiers, and experimental properties from NCBI PubChem to populate the MO &amp; Symmetry Engine
              </p>
            </div>

            <form
              onSubmit={(e) => {
                e.preventDefault();
                handlePubChemLookup(pubchemQuery);
              }}
              className="flex items-center gap-2"
            >
              <label htmlFor="pubchem-main-search" className="sr-only">
                Search PubChem by Name, CID, or SMILES
              </label>
              <input
                id="pubchem-main-search"
                type="text"
                value={pubchemQuery}
                onChange={(e) => setPubchemQuery(e.target.value)}
                placeholder="Search PubChem (e.g. ozone, pyridine, 962)..."
                className="w-64 px-3 py-1.5 text-xs font-mono bg-slate-50 border border-slate-200 rounded-md focus:outline-none focus:border-blue-600 text-slate-900"
              />
              <button
                type="submit"
                disabled={pubchemLoading || !pubchemQuery.trim()}
                className="px-3.5 py-1.5 text-xs font-medium text-white bg-slate-900 rounded-md hover:bg-slate-800 disabled:opacity-50 transition-colors flex items-center gap-1.5 whitespace-nowrap shrink-0"
              >
                <Search className="w-3.5 h-3.5" />
                <span>{pubchemLoading ? 'Querying PubChem...' : 'Query PubChem API'}</span>
              </button>
            </form>
          </div>

          {pubchemError && (
            <div className="p-3 bg-rose-50 border border-rose-200 rounded text-xs text-rose-800 font-mono">
              ▲ PubChem Query Error: {pubchemError}
            </div>
          )}

          {pubchemProfile ? (
            <div className="grid grid-cols-1 md:grid-cols-4 gap-4 pt-2 border-t border-slate-100 text-xs">
              <div>
                <span className="text-slate-500 block text-[11px]">Compound / IUPAC Name</span>
                <strong className="text-slate-900 font-medium block mt-0.5">
                  {pubchemProfile.title} ({pubchemProfile.iupacName})
                </strong>
                <span className="text-[11px] text-slate-500 font-mono mt-1 block">
                  Provenance: Experimentally determined / PubChem PUG REST
                </span>
              </div>
              <div>
                <span className="text-slate-500 block text-[11px]">Formula &amp; Molecular Weight</span>
                <span className="font-mono font-semibold text-slate-900 tabular-nums block mt-0.5">
                  <Latex math={formulaToLatex(pubchemProfile.molecularFormula)} /> ·{' '}
                  {pubchemProfile.molecularWeight.toFixed(2)} g/mol
                </span>
                <span className="font-mono text-[11px] text-slate-500 block mt-0.5 truncate">
                  InChIKey: {pubchemProfile.inchikey}
                </span>
              </div>
              <div>
                <span className="text-slate-500 block text-[11px]">Canonical SMILES &amp; Polarity</span>
                <span className="font-mono font-semibold text-blue-700 block mt-0.5 truncate">
                  {pubchemProfile.canonicalSmiles}
                </span>
                <span className="font-mono text-[11px] text-slate-600 block mt-0.5 tabular-nums">
                  TPSA: {pubchemProfile.tpsa ?? 0} <Latex math="\text{\AA}^2" /> · XLogP:{' '}
                  {pubchemProfile.xlogp ?? 'N/A'}
                </span>
              </div>
              <div className="flex flex-col justify-between items-start md:items-end">
                <span className="text-slate-500 text-[11px]">Authoritative Record</span>
                <a
                  href={pubchemProfile.pubchemRecordUrl}
                  target="_blank"
                  rel="noopener noreferrer"
                  className="px-3 py-1.5 bg-blue-50 text-blue-700 border border-blue-200 rounded font-mono text-xs font-medium hover:bg-blue-100 transition-colors inline-flex items-center gap-1.5 whitespace-nowrap"
                >
                  <span>View PubChem CID {pubchemProfile.cid}</span>
                  <ExternalLink className="w-3.5 h-3.5" />
                </a>
              </div>
            </div>
          ) : (
            <div className="pt-2 border-t border-slate-100 flex flex-wrap items-center justify-between gap-2 text-xs text-slate-500">
              <span>
                Active Provenance:{' '}
                <strong className="text-slate-700">
                  {customTitle
                    ? 'Predicted by Extended Hückel LCAO-MO & Jacobi Eigenvalue Solver'
                    : 'Curated NIST WebBook Photoelectron Spectroscopy & CCCBDB Group Theory Benchmark'}
                </strong>
              </span>
              {molecule.pubchemCid > 0 && (
                <a
                  href={`https://pubchem.ncbi.nlm.nih.gov/compound/${molecule.pubchemCid}`}
                  target="_blank"
                  rel="noopener noreferrer"
                  className="font-mono text-blue-600 hover:underline inline-flex items-center gap-1"
                >
                  <span>PubChem Reference Record CID {molecule.pubchemCid}</span>
                  <ExternalLink className="w-3 h-3" />
                </a>
              )}
            </div>
          )}
        </section>
      </main>

      {/* Quiet Editorial Footer */}
      <footer className="px-6 py-4 bg-white border-t border-slate-200 text-xs text-slate-500 flex flex-wrap items-center justify-between gap-4">
        <span>
          SymmOrbital — Molecular Orbital Theory, LCAO Wavefunction &amp; Point-Group Symmetry Workbench (Classic v1 Layout)
        </span>
        <span>
          Chemical Data Provenance: NIST Chemistry WebBook PES · NIST CCCBDB · NCBI PubChem PUG REST API
        </span>
      </footer>
    </div>
  );
}
