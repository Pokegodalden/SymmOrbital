/**
 * @license
 * SPDX-License-Identifier: Apache-2.0
 */

import React, { useEffect, useMemo, useState } from 'react';
import { AtomNode, BondEdge } from '../chemistry/moleculesData';
import {
  WalshCoordinateType,
  WalshFrame,
  runWalshDiagramScan,
} from '../chemistry/walshEngine';
import { formulaToLatex, schoenfliesToLatex } from '../chemistry/latexUtils';
import { ChemText, Latex } from './Latex';
import { Pause, Play, RotateCcw, TrendingUp } from 'lucide-react';

export interface WalshLiveOverride {
  atoms: AtomNode[];
  frame: WalshFrame;
  coordLabel: string;
  valueFormatted: string;
}

interface WalshDiagramPanelProps {
  atoms: AtomNode[];
  bonds: BondEdge[];
  charge: number;
  moleculeName: string;
  moleculeFormula: string;
  selectedMoIndex: number;
  onWalshGeometryOverride: (override: WalshLiveOverride | null) => void;
  onSelectMoIndex?: (moIndex: number) => void;
  darkMode?: boolean;
}

const CURVE_PALETTE = [
  '#64748B',
  '#2563EB',
  '#0D9488',
  '#0284C7',
  '#7C3AED',
  '#D97706',
  '#DC2626',
  '#059669',
  '#4F46E5',
  '#EA580C',
];

export const WalshDiagramPanel: React.FC<WalshDiagramPanelProps> = ({
  atoms,
  bonds,
  charge,
  moleculeName,
  moleculeFormula,
  selectedMoIndex,
  onWalshGeometryOverride,
  onSelectMoIndex,
  darkMode = false,
}) => {
  const [preferredType, setPreferredType] = useState<WalshCoordinateType | undefined>(undefined);
  const [frameIndex, setFrameIndex] = useState<number>(0);
  const [isPlaying, setIsPlaying] = useState<boolean>(false);
  const [liveSync3D, setLiveSync3D] = useState<boolean>(false);

  const scan = useMemo(
    () => runWalshDiagramScan(atoms, bonds, charge, preferredType, 19),
    [atoms, bonds, charge, preferredType]
  );

  // Locate equilibrium frame index whenever scan changes
  const eqFrameIdx = useMemo(() => {
    let best = 0;
    let bestDiff = Infinity;
    for (let i = 0; i < scan.frames.length; i++) {
      const d = Math.abs(scan.frames[i].value - scan.coordinate.equilibrium);
      if (d < bestDiff) {
        bestDiff = d;
        best = i;
      }
    }
    return best;
  }, [scan]);

  // Reset to equilibrium when underlying molecule changes
  useEffect(() => {
    setFrameIndex(eqFrameIdx);
    setIsPlaying(false);
    setLiveSync3D(false);
    onWalshGeometryOverride(null);
  }, [atoms, bonds, charge]);

  // Reset frame index when coordinate type changes
  useEffect(() => {
    setFrameIndex(eqFrameIdx);
  }, [eqFrameIdx]);

  // Push live geometry override into 3D viewport when user scrubs or plays
  const activeFrame = scan.frames[Math.min(frameIndex, scan.frames.length - 1)] || scan.frames[0];

  useEffect(() => {
    if (!activeFrame || !liveSync3D) {
      onWalshGeometryOverride(null);
      return;
    }
    onWalshGeometryOverride({
      atoms: activeFrame.atoms,
      frame: activeFrame,
      coordLabel: scan.coordinate.label,
      valueFormatted: `${activeFrame.value.toFixed(scan.coordinate.unit === '×' ? 2 : 0)}${scan.coordinate.unit}`,
    });
  }, [activeFrame, liveSync3D, scan.coordinate.label, scan.coordinate.unit]);

  // Cleanup override on unmount
  useEffect(() => {
    return () => onWalshGeometryOverride(null);
  }, []);

  // Playback loop
  useEffect(() => {
    if (!isPlaying) return;
    setLiveSync3D(true);
    let dir = 1;
    const id = window.setInterval(() => {
      setFrameIndex((prev) => {
        const next = prev + dir;
        if (next >= scan.frames.length - 1) {
          dir = -1;
          return scan.frames.length - 1;
        }
        if (next <= 0) {
          dir = 1;
          return 0;
        }
        return next;
      });
    }, 180);
    return () => window.clearInterval(id);
  }, [isPlaying, scan.frames.length]);

  const handleScrub = (idx: number) => {
    setIsPlaying(false);
    setLiveSync3D(true);
    setFrameIndex(Math.max(0, Math.min(scan.frames.length - 1, idx)));
  };

  const handleResetEquilibrium = () => {
    setIsPlaying(false);
    setLiveSync3D(false);
    setFrameIndex(eqFrameIdx);
    onWalshGeometryOverride(null);
  };

  // SVG dimensions for Correlation Diagram (full-width scientific aspect)
  const W = 1120;
  const H = 250;
  const padL = 56;
  const padR = 92;
  const padT = 18;
  const padB = 34;
  const plotW = W - padL - padR;
  const plotH = H - padT - padB;

  const { eMin, eMax } = useMemo(() => {
    let lo = Infinity;
    let hi = -Infinity;
    for (const f of scan.frames) {
      for (const moIdx of scan.trackedMoIndices) {
        const e = f.energies[moIdx];
        if (e !== undefined) {
          if (e < lo) lo = e;
          if (e > hi) hi = e;
        }
      }
    }
    if (!Number.isFinite(lo)) return { eMin: -25, eMax: 5 };
    const pad = Math.max(1.5, (hi - lo) * 0.12);
    return { eMin: lo - pad, eMax: hi + pad };
  }, [scan]);

  const xForFrame = (idx: number) =>
    padL + (idx / Math.max(1, scan.frames.length - 1)) * plotW;
  const yForEnergy = (e: number) =>
    padT + plotH - ((e - eMin) / Math.max(1e-3, eMax - eMin)) * plotH;

  // Total energy subplot dimensions
  const HTot = 110;
  const padTTot = 14;
  const padBTot = 26;
  const plotHTot = HTot - padTTot - padBTot;

  const { totMin, totMax, minTotFrameIdx } = useMemo(() => {
    let lo = Infinity;
    let hi = -Infinity;
    let minIdx = 0;
    scan.frames.forEach((f, idx) => {
      if (f.totalElectronicEnergy < lo) {
        lo = f.totalElectronicEnergy;
        minIdx = idx;
      }
      if (f.totalElectronicEnergy > hi) hi = f.totalElectronicEnergy;
    });
    const span = Math.max(0.25, hi - lo);
    return {
      totMin: lo - span * 0.15,
      totMax: hi + span * 0.15,
      minTotFrameIdx: minIdx,
    };
  }, [scan.frames]);

  const yForTot = (e: number) =>
    padTTot + plotHTot - ((e - totMin) / Math.max(1e-3, totMax - totMin)) * plotHTot;

  const axisStroke = darkMode ? '#475569' : '#CBD5E1';
  const gridStroke = darkMode ? '#1E293B' : '#F1F5F9';
  const textFill = darkMode ? '#94A3B8' : '#475569';

  // De-collide right-edge MO labels when orbitals become degenerate (e.g. 3a1 + 1b1 -> π_u at 180°)
  const labelYByMoIdx = useMemo(() => {
    const items = scan.trackedMoIndices.map((moIdx) => {
      const lastE = scan.frames[scan.frames.length - 1]?.energies[moIdx] ?? 0;
      return { moIdx, y: yForEnergy(lastE) };
    });
    items.sort((a, b) => a.y - b.y);
    const minGap = 12;
    for (let i = 1; i < items.length; i++) {
      if (items[i].y - items[i - 1].y < minGap) {
        items[i].y = items[i - 1].y + minGap;
      }
    }
    const map: Record<number, number> = {};
    for (const it of items) map[it.moIdx] = it.y;
    return map;
  }, [scan, eMin, eMax]);

  const homoE = activeFrame?.energies[activeFrame.homoIndex];
  const lumoE = activeFrame?.energies[activeFrame.lumoIndex];
  const gapEv =
    homoE !== undefined && lumoE !== undefined && lumoE > homoE
      ? lumoE - homoE
      : 0;

  return (
    <div
      data-testid="walsh-diagram-panel"
      className="bg-white border border-slate-200 rounded-lg p-4 space-y-3.5"
    >
      {/* Header Row */}
      <div className="flex flex-wrap items-center justify-between gap-3 pb-3 border-b border-slate-200">
        <div>
          <h2 className="text-sm font-semibold text-slate-900 flex items-center gap-2">
            <TrendingUp className="w-4 h-4 text-blue-600 shrink-0" />
            <span>
              Interactive Walsh Diagram · <ChemText text={moleculeName} /> (
              <Latex math={formulaToLatex(moleculeFormula)} />)
            </span>
          </h2>
          <p className="text-xs text-slate-500 mt-0.5">{scan.coordinate.description}</p>
        </div>

        {/* Perturbation Coordinate Selector */}
        <div className="flex flex-wrap items-center gap-1.5">
          {scan.availableCoordinates.map((coord) => {
            const active = coord.type === scan.coordinate.type;
            return (
              <button
                key={coord.type}
                type="button"
                onClick={() => {
                  setPreferredType(coord.type);
                  setIsPlaying(false);
                }}
                className={`px-2.5 py-1 text-xs font-medium rounded border transition-colors whitespace-nowrap shrink-0 ${
                  active
                    ? 'bg-slate-900 text-white border-slate-900'
                    : 'bg-slate-50 text-slate-700 border-slate-200 hover:bg-slate-100'
                }`}
              >
                {coord.label}
              </button>
            );
          })}
        </div>
      </div>

      {/* Interactive Scrubber & Live Geometry Telemetry */}
      <div className="bg-slate-50 border border-slate-200 rounded-lg p-3 flex flex-col lg:flex-row lg:items-center justify-between gap-3">
        <div className="flex items-center gap-2 flex-1 min-w-0">
          <button
            type="button"
            onClick={() => setIsPlaying((p) => !p)}
            data-testid="walsh-play-button"
            className="px-2.5 py-1.5 text-xs font-medium text-white bg-blue-600 hover:bg-blue-700 rounded flex items-center gap-1.5 whitespace-nowrap shrink-0 transition-colors"
          >
            {isPlaying ? <Pause className="w-3.5 h-3.5" /> : <Play className="w-3.5 h-3.5" />}
            <span>{isPlaying ? 'Pause' : 'Animate'}</span>
          </button>

          <button
            type="button"
            onClick={handleResetEquilibrium}
            data-testid="walsh-reset-button"
            className="px-2.5 py-1.5 text-xs font-medium text-slate-700 bg-white border border-slate-300 hover:bg-slate-100 rounded flex items-center gap-1.5 whitespace-nowrap shrink-0 transition-colors"
          >
            <RotateCcw className="w-3.5 h-3.5" />
            <span>Equilibrium</span>
          </button>

          <div className="flex items-center gap-2 flex-1 min-w-[160px] pl-1">
            <span className="text-[11px] font-mono text-slate-500 whitespace-nowrap shrink-0">
              {scan.coordinate.min}
              {scan.coordinate.unit}
            </span>
            <input
              type="range"
              aria-label="Walsh geometry perturbation coordinate"
              data-testid="walsh-scrubber"
              min={0}
              max={Math.max(0, scan.frames.length - 1)}
              step={1}
              value={frameIndex}
              onChange={(e) => handleScrub(Number(e.target.value))}
              className="flex-1 accent-blue-600 cursor-pointer"
            />
            <span className="text-[11px] font-mono text-slate-500 whitespace-nowrap shrink-0">
              {scan.coordinate.max}
              {scan.coordinate.unit}
            </span>
          </div>
        </div>

        {/* Live Telemetry Strip (Zero-Pill) */}
        {activeFrame && (
          <div className="flex flex-wrap items-center gap-2 text-xs font-mono text-slate-700 tabular-nums shrink-0">
            <span>
              Coord:{' '}
              <strong className="text-blue-700" data-testid="walsh-coord-readout">
                {activeFrame.value.toFixed(scan.coordinate.unit === '×' ? 2 : 0)}
                {scan.coordinate.unit}
              </strong>
            </span>
            <span className="text-slate-300" aria-hidden="true">·</span>
            <span>
              Group:{' '}
              <strong className="text-slate-900" data-testid="walsh-point-group">
                <Latex math={schoenfliesToLatex(activeFrame.schoenflies)} />
              </strong>
            </span>
            <span className="text-slate-300" aria-hidden="true">·</span>
            <span>
              Gap: <strong className="text-teal-700">{gapEv.toFixed(2)} eV</strong>
            </span>
            <span className="text-slate-300" aria-hidden="true">·</span>
            <span>
              <Latex math={`|\\vec{\\mu}| = ${activeFrame.dipoleDebye.toFixed(2)}\\;\\mathrm{D}`} />
            </span>
            {liveSync3D && (
              <>
                <span className="text-slate-300" aria-hidden="true">·</span>
                <span className="text-blue-700 font-semibold">3D Stage Live</span>
              </>
            )}
          </div>
        )}
      </div>

      {/* Main Walsh Orbital Correlation SVG */}
      <div className="border border-slate-200 rounded-lg p-3 space-y-2">
        <div className="flex flex-wrap items-center justify-between gap-2 text-[11px] font-mono text-slate-500">
          <span>
            Frontier MO Energies <Latex math="\varepsilon_i(q)" /> vs. {scan.coordinate.label} (Click curve to select MO · Click plot to jump geometry)
          </span>
          <span className="tabular-nums">
            <span className="text-teal-700 font-semibold">● HOMO</span> ·{' '}
            <span className="text-amber-700 font-semibold">● LUMO</span> ·{' '}
            <span className="text-blue-600 font-semibold">┃ Active Geometry</span>
          </span>
        </div>

        <svg
          viewBox={`0 0 ${W} ${H}`}
          className="w-full h-64 select-none cursor-pointer"
          onClick={(e) => {
            const rect = e.currentTarget.getBoundingClientRect();
            const sx = ((e.clientX - rect.left) / rect.width) * W;
            const frac = Math.max(0, Math.min(1, (sx - padL) / plotW));
            handleScrub(Math.round(frac * (scan.frames.length - 1)));
          }}
        >
          {/* Horizontal grid lines */}
          {[0, 0.25, 0.5, 0.75, 1].map((f) => {
            const eVal = eMin + f * (eMax - eMin);
            const y = yForEnergy(eVal);
            return (
              <g key={f}>
                <line x1={padL} y1={y} x2={W - padR} y2={y} stroke={gridStroke} strokeWidth={1} />
                <text
                  x={padL - 8}
                  y={y + 3.5}
                  textAnchor="end"
                  fontSize={10}
                  fill={textFill}
                  className="font-mono"
                >
                  {eVal.toFixed(1)}
                </text>
              </g>
            );
          })}

          {/* Axes */}
          <line x1={padL} y1={padT} x2={padL} y2={padT + plotH} stroke={axisStroke} strokeWidth={1} />
          <line
            x1={padL}
            y1={padT + plotH}
            x2={W - padR}
            y2={padT + plotH}
            stroke={axisStroke}
            strokeWidth={1}
          />

          {/* Equilibrium vertical guide */}
          <line
            x1={xForFrame(eqFrameIdx)}
            y1={padT}
            x2={xForFrame(eqFrameIdx)}
            y2={padT + plotH}
            stroke={axisStroke}
            strokeWidth={1}
            strokeDasharray="3 3"
          />

          {/* MO Correlation Curves */}
          {scan.trackedMoIndices.map((moIdx, ci) => {
            const isHomo = activeFrame && moIdx === activeFrame.homoIndex;
            const isLumo = activeFrame && moIdx === activeFrame.lumoIndex;
            const isSelected = moIdx === selectedMoIndex;
            const color = isHomo
              ? '#0D9488'
              : isLumo
                ? '#D97706'
                : CURVE_PALETTE[ci % CURVE_PALETTE.length];
            const pts = scan.frames
              .map((f, fi) => `${xForFrame(fi).toFixed(1)},${yForEnergy(f.energies[moIdx] ?? 0).toFixed(1)}`)
              .join(' ');
            const labelY = labelYByMoIdx[moIdx] ?? yForEnergy(
              scan.frames[scan.frames.length - 1]?.energies[moIdx] ?? 0
            );
            const activeY = yForEnergy(activeFrame?.energies[moIdx] ?? 0);

            return (
              <g
                key={moIdx}
                onClick={(ev) => {
                  ev.stopPropagation();
                  onSelectMoIndex?.(moIdx);
                }}
              >
                <polyline
                  fill="none"
                  stroke={color}
                  strokeWidth={isSelected ? 3 : isHomo || isLumo ? 2.25 : 1.5}
                  strokeOpacity={isSelected || isHomo || isLumo ? 1 : 0.72}
                  points={pts}
                />
                <circle
                  cx={xForFrame(frameIndex)}
                  cy={activeY}
                  r={isSelected ? 4.5 : isHomo || isLumo ? 3.5 : 2.5}
                  fill={color}
                />
                <text
                  x={W - padR + 8}
                  y={labelY + 3}
                  fontSize={10}
                  fill={color}
                  className="font-mono font-semibold"
                >
                  {isHomo ? `MO${moIdx + 1} HOMO` : isLumo ? `MO${moIdx + 1} LUMO` : `MO${moIdx + 1}`}
                </text>
              </g>
            );
          })}

          {/* Active Geometry Cursor Line */}
          <line
            x1={xForFrame(frameIndex)}
            y1={padT}
            x2={xForFrame(frameIndex)}
            y2={padT + plotH}
            stroke="#2563EB"
            strokeWidth={1.75}
            strokeDasharray="4 2"
          />

          {/* X-axis coordinate labels */}
          {[0, 0.25, 0.5, 0.75, 1].map((f) => {
            const idx = Math.round(f * (scan.frames.length - 1));
            const fr = scan.frames[idx];
            if (!fr) return null;
            const x = xForFrame(idx);
            return (
              <g key={f}>
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
                  {fr.value.toFixed(scan.coordinate.unit === '×' ? 2 : 0)}
                  {scan.coordinate.unit} ({fr.schoenflies})
                </text>
              </g>
            );
          })}

          <text
            x={14}
            y={padT + plotH * 0.5}
            textAnchor="middle"
            fontSize={10}
            fill={textFill}
            transform={`rotate(-90 14 ${padT + plotH * 0.5})`}
            className="font-mono"
          >
            Energy (eV)
          </text>
        </svg>

        {/* Total Electronic Energy Subplot */}
        <div className="pt-2 border-t border-slate-100">
          <div className="flex items-center justify-between text-[11px] font-mono text-slate-500 mb-1">
            <span>
              Total Valence Electronic Energy <Latex math="E_{\mathrm{elec}} = \sum_i n_i \varepsilon_i" />
            </span>
            <span className="tabular-nums text-slate-700">
              Minimum at{' '}
              <strong>
                {scan.frames[minTotFrameIdx]?.value.toFixed(
                  scan.coordinate.unit === '×' ? 2 : 0
                )}
                {scan.coordinate.unit}
              </strong>{' '}
              ({scan.frames[minTotFrameIdx]?.totalElectronicEnergy.toFixed(2)} eV)
            </span>
          </div>

          <svg
            viewBox={`0 0 ${W} ${HTot}`}
            className="w-full h-24 select-none cursor-pointer"
            onClick={(e) => {
              const rect = e.currentTarget.getBoundingClientRect();
              const sx = ((e.clientX - rect.left) / rect.width) * W;
              const frac = Math.max(0, Math.min(1, (sx - padL) / plotW));
              handleScrub(Math.round(frac * (scan.frames.length - 1)));
            }}
          >
            <line
              x1={padL}
              y1={padTTot}
              x2={padL}
              y2={padTTot + plotHTot}
              stroke={axisStroke}
              strokeWidth={1}
            />
            <line
              x1={padL}
              y1={padTTot + plotHTot}
              x2={W - padR}
              y2={padTTot + plotHTot}
              stroke={axisStroke}
              strokeWidth={1}
            />

            <text
              x={padL - 8}
              y={padTTot + 8}
              textAnchor="end"
              fontSize={9.5}
              fill={textFill}
              className="font-mono"
            >
              {totMax.toFixed(1)}
            </text>
            <text
              x={padL - 8}
              y={padTTot + plotHTot}
              textAnchor="end"
              fontSize={9.5}
              fill={textFill}
              className="font-mono"
            >
              {totMin.toFixed(1)}
            </text>

            <polyline
              fill="none"
              stroke="#2563EB"
              strokeWidth={2}
              points={scan.frames
                .map(
                  (f, fi) =>
                    `${xForFrame(fi).toFixed(1)},${yForTot(f.totalElectronicEnergy).toFixed(1)}`
                )
                .join(' ')}
            />

            {/* E_min marker */}
            <circle
              cx={xForFrame(minTotFrameIdx)}
              cy={yForTot(scan.frames[minTotFrameIdx]?.totalElectronicEnergy ?? totMin)}
              r={4}
              fill="#0D9488"
            />

            {/* Active cursor */}
            <line
              x1={xForFrame(frameIndex)}
              y1={padTTot}
              x2={xForFrame(frameIndex)}
              y2={padTTot + plotHTot}
              stroke="#2563EB"
              strokeWidth={1.5}
              strokeDasharray="4 2"
            />
            {activeFrame && (
              <circle
                cx={xForFrame(frameIndex)}
                cy={yForTot(activeFrame.totalElectronicEnergy)}
                r={4}
                fill="#2563EB"
              />
            )}
          </svg>
        </div>
      </div>
    </div>
  );
};
