/**
 * @license
 * SPDX-License-Identifier: Apache-2.0
 */

import React, { useMemo, useState } from 'react';
import { AtomNode, BondEdge, MolecularOrbital } from '../chemistry/moleculesData';
import {
  LowdinEHTResult,
  computeCOOP,
  computePDOS,
} from '../chemistry/quantumEngineV3';
import { Latex } from './Latex';
import { Activity, BarChart3, GitCommitHorizontal } from 'lucide-react';

interface PopulationBondingPanelProps {
  eht: LowdinEHTResult;
  atoms: AtomNode[];
  bonds: BondEdge[];
  molecularOrbitals: MolecularOrbital[];
  selectedMoId: string;
  onSelectMoIndex?: (moIndex: number) => void;
  darkMode?: boolean;
}

const ELEMENT_SERIES_COLORS = [
  '#2563EB', // blue
  '#0D9488', // teal
  '#D97706', // amber
  '#7C3AED', // violet
  '#DC2626', // red
  '#0891B2', // cyan
];

const SHELL_COLORS: Record<string, string> = {
  s: '#2563EB',
  p: '#0D9488',
  d: '#D97706',
  f: '#7C3AED',
};

export const PopulationBondingPanel: React.FC<PopulationBondingPanelProps> = ({
  eht,
  atoms,
  bonds,
  molecularOrbitals,
  selectedMoId,
  onSelectMoIndex,
  darkMode = false,
}) => {
  const [pdosMode, setPdosMode] = useState<'element' | 'shell'>('element');
  const [sigmaEv, setSigmaEv] = useState<number>(0.65);
  const [selectedPair, setSelectedPair] = useState<[number, number] | null>(() => {
    if (eht.mayerBonds.length > 0) {
      return [eht.mayerBonds[0].atomAIndex, eht.mayerBonds[0].atomBIndex];
    }
    return null;
  });
  const [hoverPdosIdx, setHoverPdosIdx] = useState<number | null>(null);
  const [hoverCoopIdx, setHoverCoopIdx] = useState<number | null>(null);

  // Keep selectedPair valid when molecule changes
  const activePair: [number, number] | null = useMemo(() => {
    if (!selectedPair) return null;
    if (selectedPair[0] < atoms.length && selectedPair[1] < atoms.length) {
      return selectedPair;
    }
    if (eht.mayerBonds.length > 0) {
      return [eht.mayerBonds[0].atomAIndex, eht.mayerBonds[0].atomBIndex];
    }
    return null;
  }, [selectedPair, atoms.length, eht.mayerBonds]);

  const activePairLabel = useMemo(() => {
    if (!activePair) return 'All Bonds';
    const a = atoms[activePair[0]];
    const b = atoms[activePair[1]];
    if (!a || !b) return 'All Bonds';
    return `${a.element}${activePair[0] + 1}–${b.element}${activePair[1] + 1}`;
  }, [activePair, atoms]);

  const pdos = useMemo(
    () => computePDOS(eht, { sigmaEv, numPoints: 140 }),
    [eht, sigmaEv]
  );

  const coop = useMemo(
    () =>
      computeCOOP(
        eht,
        bonds,
        activePair ? { atomAIndex: activePair[0], atomBIndex: activePair[1] } : null,
        { sigmaEv, numPoints: 140 }
      ),
    [eht, bonds, activePair, sigmaEv]
  );

  const isOpenShell = useMemo(
    () => eht.mullikenAtoms.some((a) => Math.abs(a.spinPopulation) > 1e-3),
    [eht.mullikenAtoms]
  );

  const hasDF = useMemo(
    () =>
      eht.mullikenAtoms.some(
        (a) => a.shellPopulations.d > 1e-3 || a.shellPopulations.f > 1e-3
      ),
    [eht.mullikenAtoms]
  );

  const selectedMoEnergy = useMemo(() => {
    const found = molecularOrbitals.find((m) => m.id === selectedMoId);
    return found ? found.energy : null;
  }, [molecularOrbitals, selectedMoId]);

  const selectNearestMoByEnergy = (energyEv: number) => {
    if (!onSelectMoIndex || eht.energies.length === 0) return;
    let bestIdx = 0;
    let bestDiff = Infinity;
    for (let i = 0; i < eht.energies.length; i++) {
      const diff = Math.abs(eht.energies[i] - energyEv);
      if (diff < bestDiff) {
        bestDiff = diff;
        bestIdx = i;
      }
    }
    onSelectMoIndex(bestIdx);
  };

  // Shared SVG geometry
  const W = 520;
  const H = 220;
  const padL = 44;
  const padR = 16;
  const padT = 16;
  const padB = 34;
  const plotW = W - padL - padR;
  const plotH = H - padT - padB;

  const axisStroke = darkMode ? '#475569' : '#CBD5E1';
  const gridStroke = darkMode ? '#1E293B' : '#F1F5F9';
  const textFill = darkMode ? '#94A3B8' : '#475569';
  const strongFill = darkMode ? '#E2E8F0' : '#0F172A';

  // PDOS scales
  const eMin = pdos.eMin;
  const eMax = pdos.eMax;
  const eSpan = Math.max(1e-3, eMax - eMin);
  const pdosMax = Math.max(0.2, ...pdos.points.map((pt) => pt.total));
  const xForE = (e: number) => padL + ((e - eMin) / eSpan) * plotW;
  const yForPdos = (v: number) => padT + plotH - (v / (pdosMax * 1.12)) * plotH;

  const pdosSeries = useMemo(() => {
    if (pdosMode === 'element') {
      return pdos.elements.map((el, idx) => ({
        key: el,
        label: el,
        color: ELEMENT_SERIES_COLORS[idx % ELEMENT_SERIES_COLORS.length],
        values: pdos.points.map((pt) => pt.byElement[el] || 0),
      }));
    }
    return (['s', 'p', 'd', 'f'] as const)
      .filter((sh) => pdos.points.some((pt) => pt.byShell[sh] > 1e-3))
      .map((sh) => ({
        key: sh,
        label: `${sh}-shell`,
        color: SHELL_COLORS[sh],
        values: pdos.points.map((pt) => pt.byShell[sh]),
      }));
  }, [pdos, pdosMode]);

  const buildPolyline = (xs: number[], ys: number[]) =>
    xs.map((x, i) => `${x.toFixed(1)},${ys[i].toFixed(1)}`).join(' ');

  // COOP scales
  const coopMaxAbs = Math.max(0.05, ...coop.points.map((pt) => Math.abs(pt.value)));
  const yForCoop = (v: number) =>
    padT + plotH * 0.5 - (v / (coopMaxAbs * 1.15)) * (plotH * 0.46);
  const yCoopZero = yForCoop(0);

  const coopPosPath = useMemo(() => {
    if (coop.points.length === 0) return '';
    const pts: string[] = [];
    const x0 = xForE(coop.points[0].energy);
    pts.push(`M ${x0.toFixed(1)} ${yCoopZero.toFixed(1)}`);
    for (let i = 0; i < coop.points.length; i++) {
      const x = xForE(coop.points[i].energy);
      const y = yForCoop(Math.max(0, coop.points[i].value));
      pts.push(`L ${x.toFixed(1)} ${y.toFixed(1)}`);
    }
    const xN = xForE(coop.points[coop.points.length - 1].energy);
    pts.push(`L ${xN.toFixed(1)} ${yCoopZero.toFixed(1)} Z`);
    return pts.join(' ');
  }, [coop, eMin, eSpan, coopMaxAbs]);

  const coopNegPath = useMemo(() => {
    if (coop.points.length === 0) return '';
    const pts: string[] = [];
    const x0 = xForE(coop.points[0].energy);
    pts.push(`M ${x0.toFixed(1)} ${yCoopZero.toFixed(1)}`);
    for (let i = 0; i < coop.points.length; i++) {
      const x = xForE(coop.points[i].energy);
      const y = yForCoop(Math.min(0, coop.points[i].value));
      pts.push(`L ${x.toFixed(1)} ${y.toFixed(1)}`);
    }
    const xN = xForE(coop.points[coop.points.length - 1].energy);
    pts.push(`L ${xN.toFixed(1)} ${yCoopZero.toFixed(1)} Z`);
    return pts.join(' ');
  }, [coop, eMin, eSpan, coopMaxAbs]);

  const energyTicks = useMemo(() => {
    const step = eSpan > 40 ? 10 : eSpan > 20 ? 5 : 2;
    const start = Math.ceil(eMin / step) * step;
    const ticks: number[] = [];
    for (let t = start; t <= eMax; t += step) {
      ticks.push(Math.round(t * 10) / 10);
    }
    return ticks;
  }, [eMin, eMax, eSpan]);

  const fermiX = xForE(eht.fermiEnergy);
  const selectedMoX =
    selectedMoEnergy !== null && selectedMoEnergy >= eMin && selectedMoEnergy <= eMax
      ? xForE(selectedMoEnergy)
      : null;

  const pdosXs = pdos.points.map((pt) => xForE(pt.energy));
  const coopXs = coop.points.map((pt) => xForE(pt.energy));

  return (
    <div
      data-testid="population-bonding-panel"
      className="bg-white border border-slate-200 rounded-lg p-4 space-y-4"
    >
      {/* Header & Summary Telemetry */}
      <div className="flex flex-wrap items-center justify-between gap-3 pb-3 border-b border-slate-200">
        <div>
          <h2 className="text-sm font-semibold text-slate-900 flex items-center gap-2">
            <BarChart3 className="w-4 h-4 text-teal-600 shrink-0" />
            <span>Mulliken Populations, Mayer Bond Orders, PDOS &amp; COOP</span>
          </h2>
          <p className="text-xs text-slate-500 mt-0.5">
            Exact Löwdin symmetric orthogonalization (
            <Latex math="\mathbf{S}^{-1/2}\mathbf{H}\mathbf{S}^{-1/2}\mathbf{C}' = \mathbf{C}'\mathbf{E}" />
            ) with real Clementi–Raimondi STO overlap integrals
          </p>
        </div>

        <div className="flex flex-wrap items-center gap-2 text-xs font-mono text-slate-600 tabular-nums">
          <span>
            Basis: <strong className="text-slate-900">{eht.basis.length} STOs</strong>
          </span>
          <span className="text-slate-300" aria-hidden="true">·</span>
          <span>
            <Latex math={`E_F = ${eht.fermiEnergy.toFixed(2)}\\;\\mathrm{eV}`} />
          </span>
          <span className="text-slate-300" aria-hidden="true">·</span>
          <span>
            <Latex math={`|\\vec{\\mu}| = ${eht.dipoleMagnitudeDebye.toFixed(2)}\\;\\mathrm{D}`} />
          </span>
          <span className="text-slate-300" aria-hidden="true">·</span>
          <span>
            <Latex
              math={`\\sum q_A = ${eht.mullikenAtoms
                .reduce((s, a) => s + a.netCharge, 0)
                .toFixed(2)}`}
            />
          </span>
        </div>
      </div>

      {/* Top Row: Mulliken Atomic Table + Mayer Bond Order Table */}
      <div className="grid grid-cols-1 lg:grid-cols-12 gap-4">
        {/* Mulliken Atomic Charges */}
        <div className="lg:col-span-6 border border-slate-200 rounded-lg overflow-hidden flex flex-col">
          <div className="px-3 py-2 bg-slate-50 border-b border-slate-200 flex items-center justify-between gap-2">
            <span className="text-xs font-semibold text-slate-800">
              Mulliken Atomic Charges &amp; Shell Breakdown
            </span>
            <span className="text-[11px] font-mono text-slate-500 tabular-nums">
              <Latex math="q_A = Z_A^{\mathrm{val}} - \sum_{\mu\in A}(\mathbf{PS})_{\mu\mu}" />
            </span>
          </div>

          <div className="overflow-x-auto max-h-60 overflow-y-auto">
            <table className="w-full text-xs font-mono tabular-nums">
              <thead className="bg-slate-50 text-slate-500 border-b border-slate-200 sticky top-0">
                <tr>
                  <th className="text-left px-3 py-1.5 font-medium">Atom</th>
                  <th className="text-right px-2 py-1.5 font-medium">Net q</th>
                  <th className="text-left px-2 py-1.5 font-medium w-20">Polarity</th>
                  <th className="text-right px-2 py-1.5 font-medium">Gross N</th>
                  <th className="text-right px-2 py-1.5 font-medium">s</th>
                  <th className="text-right px-2 py-1.5 font-medium">p</th>
                  {hasDF && <th className="text-right px-2 py-1.5 font-medium">d/f</th>}
                  {isOpenShell && <th className="text-right px-2 py-1.5 font-medium">Spin</th>}
                </tr>
              </thead>
              <tbody className="divide-y divide-slate-100">
                {eht.mullikenAtoms.map((ma) => {
                  const mag = Math.min(1, Math.abs(ma.netCharge) / 1.0);
                  const isNeg = ma.netCharge < -0.005;
                  const isPos = ma.netCharge > 0.005;
                  return (
                    <tr key={ma.atomId} className="hover:bg-slate-50">
                      <td className="px-3 py-1.5 font-semibold text-slate-900">
                        {ma.element}
                        <span className="text-slate-400 font-normal">{ma.atomIndex + 1}</span>
                      </td>
                      <td
                        className={`px-2 py-1.5 text-right font-semibold ${
                          isNeg
                            ? 'text-rose-600'
                            : isPos
                              ? 'text-blue-600'
                              : 'text-slate-600'
                        }`}
                      >
                        {ma.netCharge >= 0 ? '+' : ''}
                        {ma.netCharge.toFixed(3)}
                      </td>
                      <td className="px-2 py-1.5">
                        <div className="h-1.5 w-16 bg-slate-100 rounded-xs overflow-hidden flex items-center">
                          <div className="w-1/2 flex justify-end">
                            {isNeg && (
                              <div
                                className="h-1.5 bg-rose-500"
                                style={{ width: `${Math.round(mag * 100)}%` }}
                              />
                            )}
                          </div>
                          <div className="w-px h-2 bg-slate-300 shrink-0" />
                          <div className="w-1/2 flex justify-start">
                            {isPos && (
                              <div
                                className="h-1.5 bg-blue-500"
                                style={{ width: `${Math.round(mag * 100)}%` }}
                              />
                            )}
                          </div>
                        </div>
                      </td>
                      <td className="px-2 py-1.5 text-right text-slate-700">
                        {ma.grossPopulation.toFixed(2)}
                      </td>
                      <td className="px-2 py-1.5 text-right text-slate-600">
                        {ma.shellPopulations.s.toFixed(2)}
                      </td>
                      <td className="px-2 py-1.5 text-right text-slate-600">
                        {ma.shellPopulations.p.toFixed(2)}
                      </td>
                      {hasDF && (
                        <td className="px-2 py-1.5 text-right text-slate-600">
                          {(ma.shellPopulations.d + ma.shellPopulations.f).toFixed(2)}
                        </td>
                      )}
                      {isOpenShell && (
                        <td className="px-2 py-1.5 text-right text-teal-700">
                          {ma.spinPopulation >= 0 ? '+' : ''}
                          {ma.spinPopulation.toFixed(2)}
                        </td>
                      )}
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </div>

          <div className="px-3 py-2 bg-slate-50 border-t border-slate-200 flex flex-wrap items-center justify-between gap-2 text-[11px] font-mono text-slate-600 tabular-nums">
            <span>
              Dipole Vector:{' '}
              <strong className="text-slate-900">
                ({eht.dipoleVectorDebye[0].toFixed(2)}, {eht.dipoleVectorDebye[1].toFixed(2)},{' '}
                {eht.dipoleVectorDebye[2].toFixed(2)}) D
              </strong>
            </span>
            <span>
              <span className="text-rose-600 font-semibold">δ−</span> nucleophilic ·{' '}
              <span className="text-blue-600 font-semibold">δ+</span> electrophilic
            </span>
          </div>
        </div>

        {/* Mayer Bond Orders */}
        <div className="lg:col-span-6 border border-slate-200 rounded-lg overflow-hidden flex flex-col">
          <div className="px-3 py-2 bg-slate-50 border-b border-slate-200 flex items-center justify-between gap-2">
            <span className="text-xs font-semibold text-slate-800 flex items-center gap-1.5">
              <GitCommitHorizontal className="w-3.5 h-3.5 text-blue-600" />
              <span>Mayer Bond Orders (Click bond to inspect COOP)</span>
            </span>
            <button
              type="button"
              onClick={() => setSelectedPair(null)}
              className={`px-2 py-0.5 text-[11px] font-mono rounded border transition-colors whitespace-nowrap shrink-0 ${
                activePair === null
                  ? 'bg-slate-900 text-white border-slate-900'
                  : 'bg-white text-slate-600 border-slate-300 hover:bg-slate-100'
              }`}
            >
              All Bonds
            </button>
          </div>

          <div className="overflow-x-auto max-h-60 overflow-y-auto">
            {eht.mayerBonds.length > 0 ? (
              <table className="w-full text-xs font-mono tabular-nums">
                <thead className="bg-slate-50 text-slate-500 border-b border-slate-200 sticky top-0">
                  <tr>
                    <th className="text-left px-3 py-1.5 font-medium">Bond Pair</th>
                    <th className="text-right px-2 py-1.5 font-medium">R (Å)</th>
                    <th className="text-right px-2 py-1.5 font-medium">Lewis</th>
                    <th className="text-right px-2 py-1.5 font-medium">Mayer BO</th>
                    <th className="text-left px-3 py-1.5 font-medium">Character</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-100">
                  {eht.mayerBonds.slice(0, 24).map((mb) => {
                    const isSelected =
                      activePair !== null &&
                      ((activePair[0] === mb.atomAIndex && activePair[1] === mb.atomBIndex) ||
                        (activePair[0] === mb.atomBIndex && activePair[1] === mb.atomAIndex));
                    const charLabel =
                      mb.mayerOrder >= 2.4
                        ? 'Triple / strong π'
                        : mb.mayerOrder >= 1.5
                          ? 'Double / conjugated'
                          : mb.mayerOrder >= 0.75
                            ? 'Single σ covalent'
                            : mb.mayerOrder >= 0.35
                              ? 'Polar / coordinate'
                              : 'Weak / non-bonded';

                    return (
                      <tr
                        key={`${mb.atomAId}-${mb.atomBId}`}
                        onClick={() => setSelectedPair([mb.atomAIndex, mb.atomBIndex])}
                        className={`cursor-pointer transition-colors ${
                          isSelected ? 'bg-blue-50 text-slate-900' : 'hover:bg-slate-50'
                        }`}
                      >
                        <td className="px-3 py-1.5 font-semibold text-slate-900">
                          {mb.elementA}
                          <span className="text-slate-400 font-normal">{mb.atomAIndex + 1}</span>
                          {'–'}
                          {mb.elementB}
                          <span className="text-slate-400 font-normal">{mb.atomBIndex + 1}</span>
                        </td>
                        <td className="px-2 py-1.5 text-right text-slate-600">
                          {mb.distanceAngstrom.toFixed(2)}
                        </td>
                        <td className="px-2 py-1.5 text-right text-slate-500">
                          {mb.lewisOrder > 0 ? mb.lewisOrder.toFixed(1) : '—'}
                        </td>
                        <td className="px-2 py-1.5 text-right font-semibold text-teal-700">
                          {mb.mayerOrder.toFixed(3)}
                        </td>
                        <td className="px-3 py-1.5 text-slate-600 font-sans text-[11px]">
                          {charLabel}
                        </td>
                      </tr>
                    );
                  })}
                </tbody>
              </table>
            ) : (
              <div className="p-4 text-xs text-slate-500">
                No significant covalent bond orders detected above threshold.
              </div>
            )}
          </div>

          <div className="px-3 py-2 bg-slate-50 border-t border-slate-200 text-[11px] font-mono text-slate-500">
            <Latex math="B_{AB} = \sum_{\mu\in A}\sum_{\nu\in B}(\mathbf{PS})_{\mu\nu}(\mathbf{PS})_{\nu\mu}" />
          </div>
        </div>
      </div>

      {/* Bottom Row: Side-by-Side PDOS and COOP SVG Curves */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
        {/* PDOS Chart */}
        <div className="border border-slate-200 rounded-lg p-3 space-y-2">
          <div className="flex flex-wrap items-center justify-between gap-2">
            <div>
              <h3 className="text-xs font-semibold text-slate-900 flex items-center gap-1.5">
                <Activity className="w-3.5 h-3.5 text-blue-600" />
                <span>Projected Density of States (PDOS)</span>
              </h3>
              <p className="text-[11px] text-slate-500">
                Gaussian-broadened orbital states · Click any energy to select nearest MO
              </p>
            </div>

            <div className="flex items-center gap-1.5">
              <div className="flex items-center bg-slate-100 p-0.5 rounded text-[11px]">
                <button
                  type="button"
                  onClick={() => setPdosMode('element')}
                  className={`px-2 py-0.5 rounded font-medium transition-colors whitespace-nowrap shrink-0 ${
                    pdosMode === 'element'
                      ? 'bg-white text-slate-900 shadow-2xs'
                      : 'text-slate-600 hover:text-slate-900'
                  }`}
                >
                  Element
                </button>
                <button
                  type="button"
                  onClick={() => setPdosMode('shell')}
                  className={`px-2 py-0.5 rounded font-medium transition-colors whitespace-nowrap shrink-0 ${
                    pdosMode === 'shell'
                      ? 'bg-white text-slate-900 shadow-2xs'
                      : 'text-slate-600 hover:text-slate-900'
                  }`}
                >
                  s/p/d/f Shell
                </button>
              </div>

              <select
                aria-label="Gaussian broadening sigma"
                value={sigmaEv}
                onChange={(e) => setSigmaEv(Number(e.target.value))}
                className="text-[11px] font-mono bg-slate-50 border border-slate-200 rounded px-1.5 py-0.5 text-slate-700"
              >
                <option value={0.4}>σ=0.4 eV</option>
                <option value={0.65}>σ=0.65 eV</option>
                <option value={1.1}>σ=1.1 eV</option>
              </select>
            </div>
          </div>

          <div className="flex flex-wrap items-center gap-3 text-[11px] font-mono text-slate-600">
            <span className="inline-flex items-center gap-1">
              <span className="w-3 h-0.5 border-t border-dashed border-slate-500 inline-block" />
              Total DOS
            </span>
            {pdosSeries.map((s) => (
              <span key={s.key} className="inline-flex items-center gap-1">
                <span
                  className="w-2.5 h-2.5 rounded-xs inline-block"
                  style={{ backgroundColor: s.color }}
                />
                {s.label}
              </span>
            ))}
            {hoverPdosIdx !== null && pdos.points[hoverPdosIdx] && (
              <span className="ml-auto text-slate-900 font-semibold tabular-nums">
                E = {pdos.points[hoverPdosIdx].energy.toFixed(2)} eV · DOS ={' '}
                {pdos.points[hoverPdosIdx].total.toFixed(2)}
              </span>
            )}
          </div>

          <svg
            viewBox={`0 0 ${W} ${H}`}
            className="w-full h-52 select-none cursor-crosshair"
            onMouseMove={(e) => {
              const rect = e.currentTarget.getBoundingClientRect();
              const sx = ((e.clientX - rect.left) / rect.width) * W;
              const frac = Math.max(0, Math.min(1, (sx - padL) / plotW));
              setHoverPdosIdx(Math.round(frac * (pdos.points.length - 1)));
            }}
            onMouseLeave={() => setHoverPdosIdx(null)}
            onClick={() => {
              if (hoverPdosIdx !== null && pdos.points[hoverPdosIdx]) {
                selectNearestMoByEnergy(pdos.points[hoverPdosIdx].energy);
              }
            }}
          >
            {[0.25, 0.5, 0.75, 1.0].map((f) => {
              const y = padT + plotH * (1 - f);
              return (
                <line
                  key={f}
                  x1={padL}
                  y1={y}
                  x2={W - padR}
                  y2={y}
                  stroke={gridStroke}
                  strokeWidth={1}
                />
              );
            })}

            <line x1={padL} y1={padT} x2={padL} y2={padT + plotH} stroke={axisStroke} strokeWidth={1} />
            <line
              x1={padL}
              y1={padT + plotH}
              x2={W - padR}
              y2={padT + plotH}
              stroke={axisStroke}
              strokeWidth={1}
            />

            {energyTicks.map((t) => {
              const x = xForE(t);
              return (
                <g key={t}>
                  <line
                    x1={x}
                    y1={padT + plotH}
                    x2={x}
                    y2={padT + plotH + 4}
                    stroke={axisStroke}
                    strokeWidth={1}
                  />
                  <text
                    x={x}
                    y={H - 12}
                    textAnchor="middle"
                    fontSize={10}
                    fill={textFill}
                    className="font-mono"
                  >
                    {t}
                  </text>
                </g>
              );
            })}
            <text
              x={padL + plotW * 0.5}
              y={H - 1}
              textAnchor="middle"
              fontSize={10}
              fill={textFill}
              className="font-mono"
            >
              Energy (eV)
            </text>

            {pdosSeries.map((s) => {
              const ys = s.values.map(yForPdos);
              const areaD = `M ${pdosXs[0]?.toFixed(1) ?? padL} ${(padT + plotH).toFixed(1)} ${pdosXs
                .map((x, i) => `L ${x.toFixed(1)} ${ys[i].toFixed(1)}`)
                .join(' ')} L ${pdosXs[pdosXs.length - 1]?.toFixed(1) ?? padL} ${(padT + plotH).toFixed(1)} Z`;
              return (
                <g key={s.key}>
                  <path d={areaD} fill={s.color} fillOpacity={0.14} />
                  <polyline
                    fill="none"
                    stroke={s.color}
                    strokeWidth={1.75}
                    points={buildPolyline(pdosXs, ys)}
                  />
                </g>
              );
            })}

            <polyline
              fill="none"
              stroke={strongFill}
              strokeOpacity={0.55}
              strokeWidth={1.25}
              strokeDasharray="4 3"
              points={buildPolyline(
                pdosXs,
                pdos.points.map((pt) => yForPdos(pt.total))
              )}
            />

            {selectedMoX !== null && (
              <line
                x1={selectedMoX}
                y1={padT}
                x2={selectedMoX}
                y2={padT + plotH}
                stroke="#0D9488"
                strokeWidth={1.25}
                strokeOpacity={0.75}
              />
            )}

            <line
              x1={fermiX}
              y1={padT}
              x2={fermiX}
              y2={padT + plotH}
              stroke="#DC2626"
              strokeWidth={1.25}
              strokeDasharray="3 3"
            />
            <text
              x={Math.min(W - padR - 22, fermiX + 4)}
              y={padT + 10}
              fontSize={9.5}
              fill="#DC2626"
              className="font-mono font-semibold"
            >
              E_F
            </text>

            {hoverPdosIdx !== null && pdos.points[hoverPdosIdx] && (
              <line
                x1={xForE(pdos.points[hoverPdosIdx].energy)}
                y1={padT}
                x2={xForE(pdos.points[hoverPdosIdx].energy)}
                y2={padT + plotH}
                stroke={strongFill}
                strokeOpacity={0.45}
                strokeWidth={1}
              />
            )}
          </svg>
        </div>

        {/* COOP Chart */}
        <div className="border border-slate-200 rounded-lg p-3 space-y-2">
          <div className="flex flex-wrap items-center justify-between gap-2">
            <div>
              <h3 className="text-xs font-semibold text-slate-900">
                COOP Bonding Curve ·{' '}
                <span className="font-mono text-blue-700">{activePairLabel}</span>
              </h3>
              <p className="text-[11px] text-slate-500">
                <span className="text-teal-700 font-medium">+ Bonding</span> vs.{' '}
                <span className="text-amber-700 font-medium">− Antibonding</span> overlap population
              </p>
            </div>

            <div className="text-xs font-mono tabular-nums text-right">
              <span className="text-slate-500">ICOOP: </span>
              <strong className={coop.icoopAtFermi >= 0 ? 'text-teal-700' : 'text-amber-700'}>
                {coop.icoopAtFermi >= 0 ? '+' : ''}
                {coop.icoopAtFermi.toFixed(3)}
              </strong>
            </div>
          </div>

          <div className="flex flex-wrap items-center justify-between gap-2 text-[11px] font-mono text-slate-600">
            <div className="flex items-center gap-3">
              <span className="inline-flex items-center gap-1">
                <span className="w-2.5 h-2.5 rounded-xs bg-teal-600 inline-block" />
                Bonding (+COOP)
              </span>
              <span className="inline-flex items-center gap-1">
                <span className="w-2.5 h-2.5 rounded-xs bg-amber-600 inline-block" />
                Antibonding (−COOP)
              </span>
            </div>
            {hoverCoopIdx !== null && coop.points[hoverCoopIdx] && (
              <span className="text-slate-900 font-semibold tabular-nums">
                E = {coop.points[hoverCoopIdx].energy.toFixed(2)} eV · COOP ={' '}
                {coop.points[hoverCoopIdx].value >= 0 ? '+' : ''}
                {coop.points[hoverCoopIdx].value.toFixed(3)}
              </span>
            )}
          </div>

          <svg
            viewBox={`0 0 ${W} ${H}`}
            className="w-full h-52 select-none cursor-crosshair"
            onMouseMove={(e) => {
              const rect = e.currentTarget.getBoundingClientRect();
              const sx = ((e.clientX - rect.left) / rect.width) * W;
              const frac = Math.max(0, Math.min(1, (sx - padL) / plotW));
              setHoverCoopIdx(Math.round(frac * (coop.points.length - 1)));
            }}
            onMouseLeave={() => setHoverCoopIdx(null)}
            onClick={() => {
              if (hoverCoopIdx !== null && coop.points[hoverCoopIdx]) {
                selectNearestMoByEnergy(coop.points[hoverCoopIdx].energy);
              }
            }}
          >
            <line
              x1={padL}
              y1={yCoopZero}
              x2={W - padR}
              y2={yCoopZero}
              stroke={axisStroke}
              strokeWidth={1.25}
            />
            <line x1={padL} y1={padT} x2={padL} y2={padT + plotH} stroke={axisStroke} strokeWidth={1} />
            <line
              x1={padL}
              y1={padT + plotH}
              x2={W - padR}
              y2={padT + plotH}
              stroke={axisStroke}
              strokeWidth={1}
            />

            <text x={padL - 6} y={padT + 10} textAnchor="end" fontSize={9.5} fill="#0D9488" className="font-mono">
              +B
            </text>
            <text x={padL - 6} y={yCoopZero + 3} textAnchor="end" fontSize={9.5} fill={textFill} className="font-mono">
              0
            </text>
            <text x={padL - 6} y={padT + plotH - 4} textAnchor="end" fontSize={9.5} fill="#D97706" className="font-mono">
              −AB
            </text>

            <path d={coopPosPath} fill="#0D9488" fillOpacity={0.28} />
            <path d={coopNegPath} fill="#D97706" fillOpacity={0.28} />
            <polyline
              fill="none"
              stroke={strongFill}
              strokeWidth={1.5}
              strokeOpacity={0.85}
              points={buildPolyline(
                coopXs,
                coop.points.map((pt) => yForCoop(pt.value))
              )}
            />

            {energyTicks.map((t) => {
              const x = xForE(t);
              return (
                <g key={t}>
                  <line
                    x1={x}
                    y1={padT + plotH}
                    x2={x}
                    y2={padT + plotH + 4}
                    stroke={axisStroke}
                    strokeWidth={1}
                  />
                  <text
                    x={x}
                    y={H - 12}
                    textAnchor="middle"
                    fontSize={10}
                    fill={textFill}
                    className="font-mono"
                  >
                    {t}
                  </text>
                </g>
              );
            })}
            <text
              x={padL + plotW * 0.5}
              y={H - 1}
              textAnchor="middle"
              fontSize={10}
              fill={textFill}
              className="font-mono"
            >
              Energy (eV)
            </text>

            {selectedMoX !== null && (
              <line
                x1={selectedMoX}
                y1={padT}
                x2={selectedMoX}
                y2={padT + plotH}
                stroke="#2563EB"
                strokeWidth={1.25}
                strokeOpacity={0.75}
              />
            )}

            <line
              x1={fermiX}
              y1={padT}
              x2={fermiX}
              y2={padT + plotH}
              stroke="#DC2626"
              strokeWidth={1.25}
              strokeDasharray="3 3"
            />
            <text
              x={Math.min(W - padR - 22, fermiX + 4)}
              y={padT + 10}
              fontSize={9.5}
              fill="#DC2626"
              className="font-mono font-semibold"
            >
              E_F
            </text>

            {hoverCoopIdx !== null && coop.points[hoverCoopIdx] && (
              <line
                x1={xForE(coop.points[hoverCoopIdx].energy)}
                y1={padT}
                x2={xForE(coop.points[hoverCoopIdx].energy)}
                y2={padT + plotH}
                stroke={strongFill}
                strokeOpacity={0.45}
                strokeWidth={1}
              />
            )}
          </svg>
        </div>
      </div>
    </div>
  );
};
