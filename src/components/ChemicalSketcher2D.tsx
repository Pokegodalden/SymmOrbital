import React, { useEffect, useRef, useState } from 'react';
import {
  AtomNode,
  BondEdge,
  ElementSymbol,
} from '../chemistry/moleculesData';
import {
  computeMolecularFormula,
  generateSmilesFromGraph,
} from '../chemistry/moEngine';
import { formulaToLatex } from '../chemistry/latexUtils';
import { Latex } from './Latex';
import {
  ChevronDown,
  ChevronUp,
  Eraser,
  Hexagon,
  MousePointer,
  PanelLeftClose,
  Redo2,
  RotateCcw,
  SlidersHorizontal,
  Square,
  Undo2,
} from 'lucide-react';

export type SketchTool =
  | 'select'
  | 'bond_1'
  | 'bond_2'
  | 'bond_3'
  | 'wedge'
  | 'dash'
  | 'erase'
  | 'ring_6'
  | 'ring_4'
  | 'ring_3';

interface GraphSnapshot {
  atoms: AtomNode[];
  bonds: BondEdge[];
  title?: string;
}

interface ChemicalSketcher2DProps {
  atoms: AtomNode[];
  bonds: BondEdge[];
  onChangeGraph: (atoms: AtomNode[], bonds: BondEdge[], customTitle?: string) => void;
  onTriggerPubChemSearch?: (query: string) => void;
  isSearchingPubChem?: boolean;
  onCollapseSketcher?: () => void;
  hideBottomSearch?: boolean;
}

const ELEMENTS: ElementSymbol[] = ['C', 'H', 'O', 'N', 'B', 'F', 'S', 'P', 'Cl', 'Xe'];

export const ChemicalSketcher2D: React.FC<ChemicalSketcher2DProps> = ({
  atoms,
  bonds,
  onChangeGraph,
  onTriggerPubChemSearch,
  isSearchingPubChem = false,
  onCollapseSketcher,
  hideBottomSearch = false,
}) => {
  const [activeTool, setActiveTool] = useState<SketchTool>('bond_1');
  const [activeElement, setActiveElement] = useState<ElementSymbol>('C');
  const [showHeaderOptions, setShowHeaderOptions] = useState<boolean>(false);
  const [draggingAtomId, setDraggingAtomId] = useState<string | null>(null);
  const [hoveredAtomId, setHoveredAtomId] = useState<string | null>(null);
  const [hoveredBondId, setHoveredBondId] = useState<string | null>(null);
  const [canvasFocused, setCanvasFocused] = useState<boolean>(false);
  const [smilesInput, setSmilesInput] = useState<string>('');
  const [history, setHistory] = useState<GraphSnapshot[]>([]);
  const [redoStack, setRedoStack] = useState<GraphSnapshot[]>([]);
  const svgRef = useRef<SVGSVGElement | null>(null);

  // Coordinate conversion between SVG pixels (460 x 340) and chemical Angstrom space (-3.5..3.5)
  const SVG_WIDTH = 460;
  const SVG_HEIGHT = 340;
  const SCALE = 58; // pixels per Angstrom

  const toScreenX = (x2d: number) => SVG_WIDTH / 2 + x2d * SCALE;
  const toScreenY = (y2d: number) => SVG_HEIGHT / 2 + y2d * SCALE;
  const toChemX = (px: number) => Number(((px - SVG_WIDTH / 2) / SCALE).toFixed(2));
  const toChemY = (py: number) => Number(((py - SVG_HEIGHT / 2) / SCALE).toFixed(2));

  const commitGraphChange = (
    nextAtoms: AtomNode[],
    nextBonds: BondEdge[],
    customTitle?: string
  ) => {
    setHistory((prev) => [...prev.slice(-24), { atoms, bonds }]);
    setRedoStack([]);
    onChangeGraph(nextAtoms, nextBonds, customTitle);
  };

  const handleUndo = () => {
    if (history.length === 0) return;
    const prev = history[history.length - 1];
    setRedoStack((r) => [...r, { atoms, bonds }]);
    setHistory((h) => h.slice(0, -1));
    onChangeGraph(prev.atoms, prev.bonds, prev.title);
  };

  const handleRedo = () => {
    if (redoStack.length === 0) return;
    const next = redoStack[redoStack.length - 1];
    setHistory((h) => [...h, { atoms, bonds }]);
    setRedoStack((r) => r.slice(0, -1));
    onChangeGraph(next.atoms, next.bonds, next.title);
  };

  // Keyboard shortcuts when hovering over the 2D sketcher canvas
  useEffect(() => {
    if (!canvasFocused) return;
    const onKeyDown = (e: KeyboardEvent) => {
      if (
        document.activeElement?.tagName === 'INPUT' ||
        document.activeElement?.tagName === 'TEXTAREA' ||
        document.activeElement?.tagName === 'SELECT'
      ) {
        return;
      }

      if ((e.metaKey || e.ctrlKey) && e.key.toLowerCase() === 'z') {
        e.preventDefault();
        if (e.shiftKey) handleRedo();
        else handleUndo();
        return;
      }

      const upper = e.key.toUpperCase();
      if (['C', 'H', 'O', 'N', 'B', 'F', 'S', 'P'].includes(upper)) {
        const el = upper as ElementSymbol;
        setActiveElement(el);
        if (hoveredAtomId) {
          const nextAtoms = atoms.map((a) =>
            a.id === hoveredAtomId ? { ...a, element: el } : a
          );
          commitGraphChange(nextAtoms, bonds);
        }
        return;
      }

      if (e.key === '1') setActiveTool('bond_1');
      else if (e.key === '2') setActiveTool('bond_2');
      else if (e.key === '3') setActiveTool('bond_3');
      else if (e.key === '4' || e.key.toLowerCase() === 'd') setActiveTool('dash');
      else if (e.key.toLowerCase() === 'w') setActiveTool('wedge');
      else if (e.key.toLowerCase() === 'v') setActiveTool('select');
      else if (e.key === 'Backspace' || e.key === 'Delete') {
        if (hoveredAtomId) {
          const nextAtoms = atoms.filter((a) => a.id !== hoveredAtomId);
          const nextBonds = bonds.filter(
            (b) => b.source !== hoveredAtomId && b.target !== hoveredAtomId
          );
          commitGraphChange(nextAtoms, nextBonds);
        } else if (hoveredBondId) {
          commitGraphChange(
            atoms,
            bonds.filter((b) => b.id !== hoveredBondId)
          );
        }
      }
    };

    window.addEventListener('keydown', onKeyDown);
    return () => window.removeEventListener('keydown', onKeyDown);
  }, [canvasFocused, hoveredAtomId, hoveredBondId, atoms, bonds, history, redoStack]);

  const getNextAtomId = (el: ElementSymbol, currentAtoms: AtomNode[]) => {
    const existingNums = currentAtoms
      .filter((a) => a.element === el)
      .map((a) => parseInt(a.id.replace(/\D+/g, ''), 10) || 0);
    const nextNum = (existingNums.length > 0 ? Math.max(...existingNums) : 0) + 1;
    return `${el}${nextNum}`;
  };

  const getSvgCoords = (e: React.MouseEvent<SVGSVGElement>): { x: number; y: number } => {
    if (!svgRef.current) return { x: 0, y: 0 };
    const rect = svgRef.current.getBoundingClientRect();
    const scaleX = SVG_WIDTH / (rect.width || SVG_WIDTH);
    const scaleY = SVG_HEIGHT / (rect.height || SVG_HEIGHT);
    const px = (e.clientX - rect.left) * scaleX;
    const py = (e.clientY - rect.top) * scaleY;
    return { x: toChemX(px), y: toChemY(py) };
  };

  const handleCanvasClick = (e: React.MouseEvent<SVGSVGElement>) => {
    if (e.target !== svgRef.current && (e.target as Element).tagName !== 'rect') return;
    const { x, y } = getSvgCoords(e);

    if (activeTool === 'ring_6' || activeTool === 'ring_4' || activeTool === 'ring_3') {
      const nSides = activeTool === 'ring_6' ? 6 : activeTool === 'ring_4' ? 4 : 3;
      const radius = nSides === 6 ? 1.35 : nSides === 4 ? 1.1 : 0.95;
      const newAtoms = [...atoms];
      const newBonds = [...bonds];
      const ringIds: string[] = [];

      for (let k = 0; k < nSides; k++) {
        const angle = -Math.PI / 2 + (k * 2 * Math.PI) / nSides;
        const ax = Number((x + Math.cos(angle) * radius).toFixed(2));
        const ay = Number((y + Math.sin(angle) * radius).toFixed(2));
        const id = getNextAtomId(activeElement === 'H' ? 'C' : activeElement, newAtoms);
        ringIds.push(id);
        newAtoms.push({
          id,
          element: activeElement === 'H' ? 'C' : activeElement,
          x2d: ax,
          y2d: ay,
          x3d: ax,
          y3d: -ay,
          z3d: 0,
        });
      }

      for (let k = 0; k < nSides; k++) {
        const nextK = (k + 1) % nSides;
        const order: 1 | 2 = (nSides === 6 || nSides === 4) && k % 2 === 0 ? 2 : 1;
        newBonds.push({
          id: `b_ring_${Date.now()}_${k}`,
          source: ringIds[k],
          target: ringIds[nextK],
          order,
        });
      }

      commitGraphChange(
        newAtoms,
        newBonds,
        nSides === 6 ? 'Benzene Derivative' : nSides === 4 ? 'Cyclobutadiene' : 'Cyclopropane'
      );
      return;
    }

    if (
      activeTool === 'bond_1' ||
      activeTool === 'bond_2' ||
      activeTool === 'bond_3' ||
      activeTool === 'wedge' ||
      activeTool === 'dash'
    ) {
      const id1 = getNextAtomId(activeElement, atoms);
      if (atoms.length === 0) {
        const a1: AtomNode = {
          id: id1,
          element: activeElement,
          x2d: Number((x - 0.6).toFixed(2)),
          y2d: y,
          x3d: Number((x - 0.6).toFixed(2)),
          y3d: -y,
          z3d: 0,
        };
        const id2 = getNextAtomId(activeElement, [...atoms, a1]);
        const a2: AtomNode = {
          id: id2,
          element: activeElement,
          x2d: Number((x + 0.6).toFixed(2)),
          y2d: y,
          x3d: Number((x + 0.6).toFixed(2)),
          y3d: -y,
          z3d: 0,
        };
        const order: 1 | 2 | 3 = activeTool === 'bond_2' ? 2 : activeTool === 'bond_3' ? 3 : 1;
        const stereo = activeTool === 'wedge' ? 'wedge' : activeTool === 'dash' ? 'dash' : 'none';
        commitGraphChange([a1, a2], [
          { id: `b_${Date.now()}`, source: a1.id, target: a2.id, order, stereo },
        ]);
      } else {
        let nearest: AtomNode | null = null;
        let minDist = Infinity;
        for (const a of atoms) {
          const d = Math.hypot(a.x2d - x, a.y2d - y);
          if (d < minDist) {
            minDist = d;
            nearest = a;
          }
        }
        const newAtom: AtomNode = {
          id: id1,
          element: activeElement,
          x2d: x,
          y2d: y,
          x3d: x,
          y3d: -y,
          z3d: activeTool === 'wedge' ? 0.6 : activeTool === 'dash' ? -0.6 : 0,
        };
        if (nearest && minDist < 2.4) {
          const order: 1 | 2 | 3 = activeTool === 'bond_2' ? 2 : activeTool === 'bond_3' ? 3 : 1;
          const stereo = activeTool === 'wedge' ? 'wedge' : activeTool === 'dash' ? 'dash' : 'none';
          commitGraphChange(
            [...atoms, newAtom],
            [
              ...bonds,
              {
                id: `b_${Date.now()}`,
                source: nearest.id,
                target: newAtom.id,
                order,
                stereo,
              },
            ]
          );
        } else {
          commitGraphChange([...atoms, newAtom], bonds);
        }
      }
    }
  };

  const handleAtomClick = (e: React.MouseEvent, atom: AtomNode) => {
    e.stopPropagation();

    if (activeTool === 'erase') {
      const nextAtoms = atoms.filter((a) => a.id !== atom.id);
      const nextBonds = bonds.filter((b) => b.source !== atom.id && b.target !== atom.id);
      commitGraphChange(nextAtoms, nextBonds);
      return;
    }

    if (activeTool === 'select') {
      if (atom.element !== activeElement) {
        const nextAtoms = atoms.map((a) =>
          a.id === atom.id ? { ...a, element: activeElement } : a
        );
        commitGraphChange(nextAtoms, bonds);
      }
      return;
    }

    if (
      activeTool === 'bond_1' ||
      activeTool === 'bond_2' ||
      activeTool === 'bond_3' ||
      activeTool === 'wedge' ||
      activeTool === 'dash'
    ) {
      const connected = bonds.filter((b) => b.source === atom.id || b.target === atom.id);
      let angle = -Math.PI / 6;
      if (connected.length === 1) {
        const otherId = connected[0].source === atom.id ? connected[0].target : connected[0].source;
        const other = atoms.find((a) => a.id === otherId);
        if (other) {
          const baseAngle = Math.atan2(atom.y2d - other.y2d, atom.x2d - other.x2d);
          angle = baseAngle + Math.PI / 3;
        }
      } else if (connected.length >= 2) {
        angle = -Math.PI / 2 + connected.length * 1.1;
      }

      const bondLen = activeElement === 'H' ? 0.95 : 1.25;
      const nx = Number((atom.x2d + Math.cos(angle) * bondLen).toFixed(2));
      const ny = Number((atom.y2d + Math.sin(angle) * bondLen).toFixed(2));
      const newId = getNextAtomId(activeElement, atoms);
      const newAtom: AtomNode = {
        id: newId,
        element: activeElement,
        x2d: nx,
        y2d: ny,
        x3d: nx,
        y3d: -ny,
        z3d: activeTool === 'wedge' ? 0.65 : activeTool === 'dash' ? -0.65 : 0,
      };
      const order: 1 | 2 | 3 = activeTool === 'bond_2' ? 2 : activeTool === 'bond_3' ? 3 : 1;
      const stereo = activeTool === 'wedge' ? 'wedge' : activeTool === 'dash' ? 'dash' : 'none';

      commitGraphChange(
        [...atoms, newAtom],
        [
          ...bonds,
          {
            id: `b_${Date.now()}`,
            source: atom.id,
            target: newAtom.id,
            order,
            stereo,
          },
        ]
      );
    }
  };

  const handleBondClick = (e: React.MouseEvent, bond: BondEdge) => {
    e.stopPropagation();
    if (activeTool === 'erase') {
      commitGraphChange(
        atoms,
        bonds.filter((b) => b.id !== bond.id)
      );
      return;
    }

    const updatedBonds = bonds.map((b) => {
      if (b.id !== bond.id) return b;
      if (activeTool === 'bond_2') return { ...b, order: 2 as const, stereo: 'none' as const };
      if (activeTool === 'bond_3') return { ...b, order: 3 as const, stereo: 'none' as const };
      if (activeTool === 'wedge') return { ...b, order: 1 as const, stereo: 'wedge' as const };
      if (activeTool === 'dash') return { ...b, order: 1 as const, stereo: 'dash' as const };
      const nextOrder: 1 | 2 | 3 = b.order === 1 ? 2 : b.order === 2 ? 3 : 1;
      return { ...b, order: nextOrder, stereo: 'none' as const };
    });

    const updatedAtoms =
      activeTool === 'wedge' || activeTool === 'dash'
        ? atoms.map((a) =>
            a.id === bond.target
              ? { ...a, z3d: activeTool === 'wedge' ? 0.65 : -0.65 }
              : a
          )
        : atoms;

    commitGraphChange(updatedAtoms, updatedBonds);
  };

  const handleMouseMove = (e: React.MouseEvent<SVGSVGElement>) => {
    if (!draggingAtomId) return;
    const { x, y } = getSvgCoords(e);
    const updated = atoms.map((a) =>
      a.id === draggingAtomId
        ? { ...a, x2d: x, y2d: y, x3d: x, y3d: -y }
        : a
    );
    onChangeGraph(updated, bonds);
  };

  const handleCleanGeometry = () => {
    if (atoms.length === 0) return;
    const cx = atoms.reduce((s, a) => s + a.x2d, 0) / atoms.length;
    const cy = atoms.reduce((s, a) => s + a.y2d, 0) / atoms.length;
    const centered = atoms.map((a) => ({
      ...a,
      x2d: Number((a.x2d - cx).toFixed(2)),
      y2d: Number((a.y2d - cy).toFixed(2)),
      x3d: Number((a.x2d - cx).toFixed(2)),
      y3d: Number((-(a.y2d - cy)).toFixed(2)),
    }));
    commitGraphChange(centered, bonds);
  };

  const applyTemplate = (
    templateId: 'butadiene' | 'cyclobutadiene' | 'ozone' | 'hcn' | 'so2'
  ) => {
    if (templateId === 'butadiene') {
      commitGraphChange(
        [
          { id: 'C1', element: 'C', x2d: -1.8, y2d: 0.35, x3d: -1.8, y3d: -0.35, z3d: 0 },
          { id: 'C2', element: 'C', x2d: -0.6, y2d: -0.35, x3d: -0.6, y3d: 0.35, z3d: 0 },
          { id: 'C3', element: 'C', x2d: 0.6, y2d: 0.35, x3d: 0.6, y3d: -0.35, z3d: 0 },
          { id: 'C4', element: 'C', x2d: 1.8, y2d: -0.35, x3d: 1.8, y3d: 0.35, z3d: 0 },
        ],
        [
          { id: 'b1', source: 'C1', target: 'C2', order: 2 },
          { id: 'b2', source: 'C2', target: 'C3', order: 1 },
          { id: 'b3', source: 'C3', target: 'C4', order: 2 },
        ],
        '1,3-Butadiene'
      );
    } else if (templateId === 'cyclobutadiene') {
      commitGraphChange(
        [
          { id: 'C1', element: 'C', x2d: -0.8, y2d: -0.8, x3d: -0.8, y3d: 0.8, z3d: 0 },
          { id: 'C2', element: 'C', x2d: 0.8, y2d: -0.8, x3d: 0.8, y3d: 0.8, z3d: 0 },
          { id: 'C3', element: 'C', x2d: 0.8, y2d: 0.8, x3d: 0.8, y3d: -0.8, z3d: 0 },
          { id: 'C4', element: 'C', x2d: -0.8, y2d: 0.8, x3d: -0.8, y3d: -0.8, z3d: 0 },
        ],
        [
          { id: 'b1', source: 'C1', target: 'C2', order: 2 },
          { id: 'b2', source: 'C2', target: 'C3', order: 1 },
          { id: 'b3', source: 'C3', target: 'C4', order: 2 },
          { id: 'b4', source: 'C4', target: 'C1', order: 1 },
        ],
        'Cyclobutadiene (D4h)'
      );
    } else if (templateId === 'ozone') {
      commitGraphChange(
        [
          { id: 'O1', element: 'O', x2d: 0, y2d: -0.42, x3d: 0, y3d: 0.42, z3d: 0 },
          { id: 'O2', element: 'O', x2d: -1.08, y2d: 0.25, x3d: -1.08, y3d: -0.25, z3d: 0 },
          { id: 'O3', element: 'O', x2d: 1.08, y2d: 0.25, x3d: 1.08, y3d: -0.25, z3d: 0 },
        ],
        [
          { id: 'b1', source: 'O1', target: 'O2', order: 2 },
          { id: 'b2', source: 'O1', target: 'O3', order: 1 },
        ],
        'Ozone (O3, C2v)'
      );
    } else if (templateId === 'so2') {
      commitGraphChange(
        [
          { id: 'S1', element: 'S', x2d: 0, y2d: -0.45, x3d: 0, y3d: 0.45, z3d: 0 },
          { id: 'O1', element: 'O', x2d: -1.24, y2d: 0.26, x3d: -1.24, y3d: -0.26, z3d: 0 },
          { id: 'O2', element: 'O', x2d: 1.24, y2d: 0.26, x3d: 1.24, y3d: -0.26, z3d: 0 },
        ],
        [
          { id: 'b1', source: 'S1', target: 'O1', order: 2 },
          { id: 'b2', source: 'S1', target: 'O2', order: 2 },
        ],
        'Sulfur Dioxide (SO2, C2v)'
      );
    } else if (templateId === 'hcn') {
      commitGraphChange(
        [
          { id: 'H1', element: 'H', x2d: -1.3, y2d: 0, x3d: -1.3, y3d: 0, z3d: 0 },
          { id: 'C1', element: 'C', x2d: 0, y2d: 0, x3d: 0, y3d: 0, z3d: 0 },
          { id: 'N1', element: 'N', x2d: 1.2, y2d: 0, x3d: 1.2, y3d: 0, z3d: 0 },
        ],
        [
          { id: 'b1', source: 'H1', target: 'C1', order: 1 },
          { id: 'b2', source: 'C1', target: 'N1', order: 3 },
        ],
        'Hydrogen Cyanide (C∞v)'
      );
    }
  };

  const formulaInfo = computeMolecularFormula(atoms, bonds, true);
  const currentSmiles = generateSmilesFromGraph(atoms, bonds);

  return (
    <div className="bg-white border border-slate-200 rounded-lg flex flex-col overflow-hidden h-full">
      {/* Header Bar with Collapsible Tools & Options Toggle (Collapsed by Default) */}
      <div className="px-3.5 py-2.5 bg-slate-50 border-b border-slate-200 flex items-center justify-between gap-2 shrink-0">
        <div className="flex items-center gap-2 min-w-0">
          <h2 className="text-sm font-semibold text-slate-900 truncate">
            2D Chemical Sketcher
          </h2>
          <span className="text-slate-300" aria-hidden="true">·</span>
          <span className="font-mono text-xs font-semibold text-blue-700 whitespace-nowrap">
            <Latex math={formulaToLatex(formulaInfo.formula)} />
          </span>
        </div>

        <div className="flex items-center gap-1.5 shrink-0">
          <button
            type="button"
            onClick={() => setShowHeaderOptions((v) => !v)}
            aria-expanded={showHeaderOptions}
            className={`px-2.5 py-1 text-xs font-medium rounded border transition-colors flex items-center gap-1 whitespace-nowrap ${
              showHeaderOptions
                ? 'bg-slate-900 text-white border-slate-900'
                : 'bg-white text-slate-700 border-slate-200 hover:bg-slate-100'
            }`}
          >
            <SlidersHorizontal className="w-3 h-3" />
            <span>Tools &amp; Elements</span>
            {showHeaderOptions ? (
              <ChevronUp className="w-3.5 h-3.5" />
            ) : (
              <ChevronDown className="w-3.5 h-3.5" />
            )}
          </button>

          {onCollapseSketcher && (
            <button
              type="button"
              onClick={onCollapseSketcher}
              title="Minimize 2D Sketcher Pane"
              aria-label="Minimize 2D Sketcher Pane"
              className="p-1.5 text-slate-600 bg-white border border-slate-200 rounded hover:bg-slate-100 transition-colors"
            >
              <PanelLeftClose className="w-3.5 h-3.5" />
            </button>
          )}
        </div>
      </div>

      {/* Collapsible Drawing Tools, Element Palette & Templates (Collapsed by Default) */}
      {showHeaderOptions && (
        <div className="border-b border-slate-200 bg-slate-50/90 divide-y divide-slate-200/80 shrink-0">
          {/* Row 1: Bond/Ring Tools + Undo/Redo/Clear */}
          <div
            className="px-3 py-1.5 flex flex-wrap items-center justify-between gap-1.5"
            role="toolbar"
            aria-label="2D Structure Drawing Tools"
          >
            <div className="flex flex-wrap items-center gap-1">
              <button
                type="button"
                onClick={() => setActiveTool('select')}
                title="Select / Move Atom or Relabel Element (Shortcut: V)"
                className={`px-2 py-1 text-[11px] font-medium rounded transition-colors flex items-center gap-1 whitespace-nowrap shrink-0 ${
                  activeTool === 'select'
                    ? 'bg-blue-600 text-white'
                    : 'bg-white text-slate-700 border border-slate-200 hover:bg-slate-100'
                }`}
              >
                <MousePointer className="w-3 h-3" />
                <span>Move</span>
              </button>

              <button
                type="button"
                onClick={() => setActiveTool('bond_1')}
                title="Single Bond (Shortcut: 1)"
                className={`px-2 py-1 text-[11px] font-mono font-medium rounded transition-colors whitespace-nowrap shrink-0 ${
                  activeTool === 'bond_1'
                    ? 'bg-blue-600 text-white'
                    : 'bg-white text-slate-700 border border-slate-200 hover:bg-slate-100'
                }`}
              >
                — Single
              </button>

              <button
                type="button"
                onClick={() => setActiveTool('bond_2')}
                title="Double Bond (Shortcut: 2)"
                className={`px-2 py-1 text-[11px] font-mono font-medium rounded transition-colors whitespace-nowrap shrink-0 ${
                  activeTool === 'bond_2'
                    ? 'bg-blue-600 text-white'
                    : 'bg-white text-slate-700 border border-slate-200 hover:bg-slate-100'
                }`}
              >
                ═ Double
              </button>

              <button
                type="button"
                onClick={() => setActiveTool('bond_3')}
                title="Triple Bond (Shortcut: 3)"
                className={`px-2 py-1 text-[11px] font-mono font-medium rounded transition-colors whitespace-nowrap shrink-0 ${
                  activeTool === 'bond_3'
                    ? 'bg-blue-600 text-white'
                    : 'bg-white text-slate-700 border border-slate-200 hover:bg-slate-100'
                }`}
              >
                ≡ Triple
              </button>

              <button
                type="button"
                onClick={() => setActiveTool('wedge')}
                title="Stereo Wedge Bond (+z, Shortcut: W)"
                className={`px-2 py-1 text-[11px] font-medium rounded transition-colors whitespace-nowrap shrink-0 ${
                  activeTool === 'wedge'
                    ? 'bg-blue-600 text-white'
                    : 'bg-white text-slate-700 border border-slate-200 hover:bg-slate-100'
                }`}
              >
                ▲ Wedge
              </button>

              <button
                type="button"
                onClick={() => setActiveTool('dash')}
                title="Dashed Bond / Stereo Hash (-z, Shortcut: D or 4)"
                className={`px-2 py-1 text-[11px] font-mono font-medium rounded transition-colors whitespace-nowrap shrink-0 ${
                  activeTool === 'dash'
                    ? 'bg-blue-600 text-white'
                    : 'bg-white text-slate-700 border border-slate-200 hover:bg-slate-100'
                }`}
              >
                ┄ Dash
              </button>

              <button
                type="button"
                onClick={() => setActiveTool('ring_6')}
                title="Stamp Benzene Ring (C6)"
                className={`px-2 py-1 text-[11px] font-medium rounded transition-colors flex items-center gap-1 whitespace-nowrap shrink-0 ${
                  activeTool === 'ring_6'
                    ? 'bg-blue-600 text-white'
                    : 'bg-white text-slate-700 border border-slate-200 hover:bg-slate-100'
                }`}
              >
                <Hexagon className="w-3 h-3" />
                <span>Ring-6</span>
              </button>

              <button
                type="button"
                onClick={() => setActiveTool('ring_4')}
                title="Stamp 4-Membered Ring (C4)"
                className={`px-2 py-1 text-[11px] font-medium rounded transition-colors flex items-center gap-1 whitespace-nowrap shrink-0 ${
                  activeTool === 'ring_4'
                    ? 'bg-blue-600 text-white'
                    : 'bg-white text-slate-700 border border-slate-200 hover:bg-slate-100'
                }`}
              >
                <Square className="w-3 h-3" />
                <span>Ring-4</span>
              </button>

              <button
                type="button"
                onClick={() => setActiveTool('erase')}
                title="Erase Atom or Bond (Shortcut: Del)"
                className={`px-2 py-1 text-[11px] font-medium rounded transition-colors flex items-center gap-1 whitespace-nowrap shrink-0 ${
                  activeTool === 'erase'
                    ? 'bg-rose-600 text-white'
                    : 'bg-white text-slate-700 border border-slate-200 hover:bg-rose-50 hover:text-rose-700'
                }`}
              >
                <Eraser className="w-3 h-3" />
                <span>Erase</span>
              </button>
            </div>

            <div className="flex items-center gap-1">
              <button
                type="button"
                onClick={handleUndo}
                disabled={history.length === 0}
                title="Undo (Ctrl+Z)"
                aria-label="Undo last structure edit"
                className="p-1 text-slate-600 bg-white border border-slate-200 rounded hover:bg-slate-100 disabled:opacity-40 transition-colors"
              >
                <Undo2 className="w-3.5 h-3.5" />
              </button>
              <button
                type="button"
                onClick={handleRedo}
                disabled={redoStack.length === 0}
                title="Redo (Ctrl+Shift+Z)"
                aria-label="Redo structure edit"
                className="p-1 text-slate-600 bg-white border border-slate-200 rounded hover:bg-slate-100 disabled:opacity-40 transition-colors"
              >
                <Redo2 className="w-3.5 h-3.5" />
              </button>
              <button
                type="button"
                onClick={handleCleanGeometry}
                title="Center 2D Structure Coordinates"
                className="px-2 py-1 text-[11px] font-medium text-slate-700 bg-white border border-slate-200 rounded hover:bg-slate-100 transition-colors whitespace-nowrap"
              >
                Center
              </button>
              <button
                type="button"
                onClick={() => commitGraphChange([], [], 'Empty Canvas')}
                title="Clear Sketcher Canvas"
                className="px-2 py-1 text-[11px] font-medium text-rose-700 bg-white border border-slate-200 rounded hover:bg-rose-50 transition-colors flex items-center gap-1 whitespace-nowrap"
              >
                <RotateCcw className="w-3 h-3" />
                <span>Clear</span>
              </button>
            </div>
          </div>

          {/* Row 2: Periodic Element Selector + Quick Templates */}
          <div className="px-3 py-1.5 bg-slate-100/70 flex flex-wrap items-center justify-between gap-1.5">
            <div className="flex items-center gap-1 flex-wrap" role="radiogroup" aria-label="Active Chemical Element">
              {ELEMENTS.map((el) => {
                const isSelected = activeElement === el;
                return (
                  <button
                    key={el}
                    type="button"
                    onClick={() => setActiveElement(el)}
                    title={`Select ${el} (or press '${el}' while hovering an atom)`}
                    className={`w-6 h-6 text-[11px] font-mono font-semibold rounded transition-colors flex items-center justify-center shrink-0 ${
                      isSelected
                        ? 'bg-slate-900 text-white shadow-xs'
                        : 'bg-white text-slate-700 border border-slate-200 hover:border-slate-400'
                    }`}
                    aria-label={`Element ${el}`}
                    aria-pressed={isSelected}
                  >
                    {el}
                  </button>
                );
              })}
            </div>

            <div className="flex items-center gap-1 text-[11px] flex-wrap">
              <button
                type="button"
                onClick={() => applyTemplate('ozone')}
                className="px-1.5 py-0.5 text-[11px] font-medium text-slate-700 bg-white border border-slate-200 rounded hover:bg-slate-50 whitespace-nowrap shrink-0"
              >
                <Latex math="\mathrm{O_3}" />
              </button>
              <button
                type="button"
                onClick={() => applyTemplate('so2')}
                className="px-1.5 py-0.5 text-[11px] font-medium text-slate-700 bg-white border border-slate-200 rounded hover:bg-slate-50 whitespace-nowrap shrink-0"
              >
                <Latex math="\mathrm{SO_2}" />
              </button>
              <button
                type="button"
                onClick={() => applyTemplate('butadiene')}
                className="px-1.5 py-0.5 text-[11px] font-medium text-slate-700 bg-white border border-slate-200 rounded hover:bg-slate-50 whitespace-nowrap shrink-0"
              >
                Butadiene
              </button>
              <button
                type="button"
                onClick={() => applyTemplate('hcn')}
                className="px-1.5 py-0.5 text-[11px] font-medium text-slate-700 bg-white border border-slate-200 rounded hover:bg-slate-50 whitespace-nowrap shrink-0"
              >
                HCN
              </button>
            </div>
          </div>
        </div>
      )}

      {/* Main Interactive 2D SVG Chemical Drawing Canvas */}
      <div
        className="relative bg-white select-none flex-1 min-h-0 flex flex-col justify-between"
        onMouseEnter={() => setCanvasFocused(true)}
        onMouseLeave={() => {
          setCanvasFocused(false);
          setDraggingAtomId(null);
        }}
      >
        <svg
          ref={svgRef}
          viewBox={`0 0 ${SVG_WIDTH} ${SVG_HEIGHT}`}
          className="w-full flex-1 min-h-[280px] cursor-crosshair block"
          onClick={handleCanvasClick}
          onMouseMove={handleMouseMove}
          onMouseUp={() => setDraggingAtomId(null)}
          data-testid="chemical-sketcher-svg"
        >
          <defs>
            <pattern id="chem-grid" width="29" height="29" patternUnits="userSpaceOnUse">
              <path
                d="M 29 0 L 0 0 0 29"
                fill="none"
                stroke="#F1F5F9"
                strokeWidth="1"
              />
            </pattern>
          </defs>

          <rect width={SVG_WIDTH} height={SVG_HEIGHT} fill="url(#chem-grid)" />

          {/* Subtle Symmetry Crosshair Axes */}
          <line
            x1={SVG_WIDTH / 2}
            y1={0}
            x2={SVG_WIDTH / 2}
            y2={SVG_HEIGHT}
            stroke="#E2E8F0"
            strokeDasharray="4 4"
            strokeWidth="1"
          />
          <line
            x1={0}
            y1={SVG_HEIGHT / 2}
            x2={SVG_WIDTH}
            y2={SVG_HEIGHT / 2}
            stroke="#E2E8F0"
            strokeDasharray="4 4"
            strokeWidth="1"
          />

          {/* Render Bonds */}
          {bonds.map((bond) => {
            const s = atoms.find((a) => a.id === bond.source);
            const t = atoms.find((a) => a.id === bond.target);
            if (!s || !t) return null;

            const x1 = toScreenX(s.x2d);
            const y1 = toScreenY(s.y2d);
            const x2 = toScreenX(t.x2d);
            const y2 = toScreenY(t.y2d);

            const dx = x2 - x1;
            const dy = y2 - y1;
            const len = Math.hypot(dx, dy) || 1;
            const nx = -dy / len;
            const ny = dx / len;
            const isHovered = hoveredBondId === bond.id;
            const strokeColor = isHovered
              ? activeTool === 'erase'
                ? '#DC2626'
                : '#2563EB'
              : '#1E293B';

            return (
              <g
                key={bond.id}
                className="cursor-pointer"
                onMouseEnter={() => setHoveredBondId(bond.id)}
                onMouseLeave={() => setHoveredBondId(null)}
                onClick={(e) => handleBondClick(e, bond)}
              >
                <line
                  x1={x1}
                  y1={y1}
                  x2={x2}
                  y2={y2}
                  stroke="transparent"
                  strokeWidth="14"
                />

                {bond.order === 1 && (!bond.stereo || bond.stereo === 'none') && (
                  <line
                    x1={x1}
                    y1={y1}
                    x2={x2}
                    y2={y2}
                    stroke={strokeColor}
                    strokeWidth={isHovered ? 3.2 : 2.4}
                    strokeLinecap="round"
                  />
                )}

                {bond.order === 1 && bond.stereo === 'wedge' && (
                  <polygon
                    points={`${x1},${y1} ${x2 + nx * 5},${y2 + ny * 5} ${x2 - nx * 5},${y2 - ny * 5}`}
                    fill={strokeColor}
                  />
                )}

                {bond.order === 1 && bond.stereo === 'dash' && (
                  <line
                    x1={x1}
                    y1={y1}
                    x2={x2}
                    y2={y2}
                    stroke={strokeColor}
                    strokeWidth={3}
                    strokeDasharray="3 4"
                  />
                )}

                {bond.order === 2 && (
                  <>
                    <line
                      x1={x1 + nx * 3.2}
                      y1={y1 + ny * 3.2}
                      x2={x2 + nx * 3.2}
                      y2={y2 + ny * 3.2}
                      stroke={strokeColor}
                      strokeWidth={2.2}
                      strokeLinecap="round"
                    />
                    <line
                      x1={x1 - nx * 3.2}
                      y1={y1 - ny * 3.2}
                      x2={x2 - nx * 3.2}
                      y2={y2 - ny * 3.2}
                      stroke={strokeColor}
                      strokeWidth={2.2}
                      strokeLinecap="round"
                    />
                  </>
                )}

                {bond.order === 3 && (
                  <>
                    <line
                      x1={x1}
                      y1={y1}
                      x2={x2}
                      y2={y2}
                      stroke={strokeColor}
                      strokeWidth={2.0}
                      strokeLinecap="round"
                    />
                    <line
                      x1={x1 + nx * 4.5}
                      y1={y1 + ny * 4.5}
                      x2={x2 + nx * 4.5}
                      y2={y2 + ny * 4.5}
                      stroke={strokeColor}
                      strokeWidth={2.0}
                      strokeLinecap="round"
                    />
                    <line
                      x1={x1 - nx * 4.5}
                      y1={y1 - ny * 4.5}
                      x2={x2 - nx * 4.5}
                      y2={y2 - ny * 4.5}
                      stroke={strokeColor}
                      strokeWidth={2.0}
                      strokeLinecap="round"
                    />
                  </>
                )}
              </g>
            );
          })}

          {/* Render Atoms */}
          {atoms.map((atom) => {
            const cx = toScreenX(atom.x2d);
            const cy = toScreenY(atom.y2d);
            const isHovered = hoveredAtomId === atom.id;
            const elemColor =
              atom.element === 'O'
                ? '#DC2626'
                : atom.element === 'N'
                  ? '#2563EB'
                  : atom.element === 'F' || atom.element === 'Cl'
                    ? '#059669'
                    : atom.element === 'S'
                      ? '#CA8A04'
                      : atom.element === 'B'
                        ? '#DB2777'
                        : atom.element === 'Xe'
                          ? '#0891B2'
                          : atom.element === 'H'
                            ? '#475569'
                            : '#0F172A';

            return (
              <g
                key={atom.id}
                transform={`translate(${cx}, ${cy})`}
                className="cursor-pointer"
                onMouseEnter={() => setHoveredAtomId(atom.id)}
                onMouseLeave={() => setHoveredAtomId(null)}
                onMouseDown={(e) => {
                  if (activeTool === 'select') {
                    e.stopPropagation();
                    setDraggingAtomId(atom.id);
                  }
                }}
                onClick={(e) => handleAtomClick(e, atom)}
              >
                <circle
                  r={atom.element === 'H' ? 12 : 15}
                  fill="#FFFFFF"
                  stroke={
                    isHovered
                      ? activeTool === 'erase'
                        ? '#DC2626'
                        : '#2563EB'
                      : '#CBD5E1'
                  }
                  strokeWidth={isHovered ? 2.2 : 1.4}
                />
                <text
                  textAnchor="middle"
                  dominantBaseline="central"
                  fill={elemColor}
                  className="font-mono text-xs font-bold select-none"
                >
                  {atom.element}
                </text>
                <text
                  y={atom.element === 'H' ? 19 : 23}
                  textAnchor="middle"
                  fill="#94A3B8"
                  className="font-mono text-[9px] select-none"
                >
                  {atom.id}
                </text>
              </g>
            );
          })}
        </svg>

        {atoms.length === 0 && (
          <div className="absolute inset-0 pointer-events-none flex flex-col items-center justify-center text-center p-4">
            <p className="text-xs font-medium text-slate-500">
              Click anywhere on the grid to place a bond or ring, or use the search bar above.
            </p>
          </div>
        )}

        {/* Bottom Canvas Overlay Info Strip */}
        <div className="px-3 py-2 bg-slate-50 border-t border-slate-200 flex flex-wrap items-center justify-between gap-2 text-[11px] shrink-0">
          <div className="flex items-center gap-1.5 text-slate-600 font-mono">
            <span>
              Active: <strong className="text-slate-900">{activeElement}</strong>
            </span>
            <span aria-hidden="true">·</span>
            <span>
              Atoms: <strong className="text-slate-900 tabular-nums">{atoms.length}</strong>
            </span>
            <span aria-hidden="true">·</span>
            <span>
              Bonds: <strong className="text-slate-900 tabular-nums">{bonds.length}</strong>
            </span>
            {currentSmiles && (
              <>
                <span aria-hidden="true">·</span>
                <span className="truncate max-w-[120px]" title={currentSmiles}>
                  <strong className="text-blue-700">{currentSmiles}</strong>
                </span>
              </>
            )}
          </div>
          <span className="text-[10px] text-slate-500">
            1/2/3 Bond · D/W Dash/Wedge · C/N/O/H Element · Del Erase
          </span>
        </div>
      </div>

      {/* Optional Quick Chemical Name / SMILES / PubChem Fetcher Bar (used in Classic v1 mode) */}
      {onTriggerPubChemSearch && !hideBottomSearch && (
        <form
          onSubmit={(e) => {
            e.preventDefault();
            if (smilesInput.trim()) {
              onTriggerPubChemSearch(smilesInput.trim());
            }
          }}
          className="px-2.5 py-2 bg-white border-t border-slate-200 flex items-center gap-1.5 shrink-0"
        >
          <label
            htmlFor="pubchem-quick-input"
            className="text-[11px] font-medium text-slate-600 whitespace-nowrap shrink-0"
          >
            PubChem:
          </label>
          <input
            id="pubchem-quick-input"
            type="text"
            value={smilesInput}
            onChange={(e) => setSmilesInput(e.target.value)}
            placeholder="Name, formula, or CID..."
            className="flex-1 min-w-0 px-2 py-1 text-xs font-mono bg-slate-50 border border-slate-200 rounded focus:outline-none focus:border-blue-600 text-slate-900"
          />
          <button
            type="submit"
            disabled={isSearchingPubChem || !smilesInput.trim()}
            className="px-2.5 py-1 text-xs font-medium text-white bg-blue-600 rounded hover:bg-blue-700 disabled:opacity-50 transition-colors whitespace-nowrap shrink-0"
          >
            {isSearchingPubChem ? '...' : 'Load'}
          </button>
        </form>
      )}
    </div>
  );
};
