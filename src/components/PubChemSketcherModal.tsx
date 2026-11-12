import React, { useEffect, useRef, useState } from 'react';
import {
  AtomNode,
  BondEdge,
  ELEMENT_PROPERTIES,
  ElementSymbol,
} from '../chemistry/moleculesData';
import {
  computeMolecularFormula,
  expandImplicitHydrogens,
  generateMolfileFromGraph,
  generateSmilesFromGraph,
  parseChemicalInputToGraph,
} from '../chemistry/moEngine';

export type PubChemTool =
  | 'move'
  | 'erase'
  | 'query'
  | 'bond_single'
  | 'bond_double'
  | 'bond_triple'
  | 'bond_wedge'
  | 'bond_dash'
  | 'bond_wavy'
  | 'bond_crossed'
  | 'bond_any'
  | 'bond_sa'
  | 'bond_da'
  | 'bond_sd'
  | 'ring_3'
  | 'ring_4'
  | 'ring_5'
  | 'ring_6'
  | 'ring_7'
  | 'ring_8'
  | 'ring_benzene'
  | 'ring_bicyclic'
  | 'charge_plus'
  | 'charge_minus'
  | 'charge_neutral'
  | 'chain_ethyl'
  | 'chain_propyl'
  | 'chain_ipr'
  | 'chain_butyl'
  | 'chain_tbu'
  | 'chain_neo'
  | 'chain_cross'
  | 'fg_cho'
  | 'fg_co2h'
  | 'fg_no2'
  | 'fg_so3h'
  | 'element';

interface PubChemSketcherModalProps {
  isOpen: boolean;
  initialAtoms: AtomNode[];
  initialBonds: BondEdge[];
  onAccept: (atoms: AtomNode[], bonds: BondEdge[], customTitle?: string) => void;
  onCancel: () => void;
}

const QUERY_ELEMENTS: ElementSymbol[] = ['C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I', 'Si', 'B', 'Se'];
const ROW4_TM: ElementSymbol[] = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'];
const ROW5_TM: ElementSymbol[] = ['Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd'];
const ROW6_LN_TM: ElementSymbol[] = [
  'Lu', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
  'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
];

export function getReadableAtomColor(el: ElementSymbol): string {
  if (el === 'C') return '#0F172A';
  if (el === 'H') return '#475569';
  if (el === 'O') return '#DC2626';
  if (el === 'N') return '#2563EB';
  if (el === 'S') return '#B45309';
  if (el === 'P') return '#D97706';
  if (el === 'F' || el === 'Cl') return '#059669';
  if (el === 'Br') return '#991B1B';
  if (el === 'I') return '#6B21A8';
  const raw = ELEMENT_PROPERTIES[el]?.cpkColor || '#1E293B';
  // Ensure light CPK colors remain readable on a white 2D sketcher background
  if (raw.toUpperCase() === '#FFFFFF' || raw.toUpperCase() === '#D9FFFF' || raw.toUpperCase() === '#FFFF30') {
    return '#334155';
  }
  return raw;
}

export const PubChemSketcherModal: React.FC<PubChemSketcherModalProps> = ({
  isOpen,
  initialAtoms,
  initialBonds,
  onAccept,
  onCancel,
}) => {
  const [draftAtoms, setDraftAtoms] = useState<AtomNode[]>([]);
  const [draftBonds, setDraftBonds] = useState<BondEdge[]>([]);
  const [history, setHistory] = useState<{ atoms: AtomNode[]; bonds: BondEdge[] }[]>([]);
  const [activeTool, setActiveTool] = useState<PubChemTool>('bond_single');
  const [activeElement, setActiveElement] = useState<ElementSymbol>('C');

  // Periodic table dropdown states matching PubChem Sketcher
  const [queryEl, setQueryEl] = useState<ElementSymbol>('C');
  const [tmRow4El, setTmRow4El] = useState<ElementSymbol>('Sc');
  const [tmRow5El, setTmRow5El] = useState<ElementSymbol>('Y');
  const [tmRow6El, setTmRow6El] = useState<ElementSymbol>('Lu');

  // Top bar & bottom export/hydrogen states
  const [networkMode, setNetworkMode] = useState<'Broadband' | 'Standard'>('Broadband');
  const [topFormat, setTopFormat] = useState<'SMILES' | 'Formula' | 'MDL Molfile' | 'InChI'>('SMILES');
  const [topInputText, setTopInputText] = useState<string>('');
  const [exportFormat, setExportFormat] = useState<'MDL Molfile' | 'SMILES' | 'Formula'>('MDL Molfile');
  const [hydrogenMode, setHydrogenMode] = useState<'Keep AsIs' | 'Add Implicit H' | 'Remove Explicit H'>('Keep AsIs');
  const [showExportPanel, setShowExportPanel] = useState<boolean>(false);
  const [exportText, setExportText] = useState<string>('');
  const [showHelpBanner, setShowHelpBanner] = useState<boolean>(false);
  const [importedFileName, setImportedFileName] = useState<string>('No file chosen');

  // Canvas interaction state
  const [draggingAtomId, setDraggingAtomId] = useState<string | null>(null);
  const [selectedAtomId, setSelectedAtomId] = useState<string | null>(null);
  const [cursorChemPos, setCursorChemPos] = useState<{ x: number; y: number } | null>(null);
  const [draggingCanvas, setDraggingCanvas] = useState<{ startX: number; startY: number } | null>(null);
  const [hoveredAtomId, setHoveredAtomId] = useState<string | null>(null);
  const [hoveredBondId, setHoveredBondId] = useState<string | null>(null);
  const svgRef = useRef<SVGSVGElement | null>(null);
  const fileInputRef = useRef<HTMLInputElement | null>(null);

  useEffect(() => {
    if (isOpen) {
      const clonedAtoms = initialAtoms.map((a) => ({ ...a }));
      const clonedBonds = initialBonds.map((b) => ({ ...b }));
      setDraftAtoms(clonedAtoms);
      setDraftBonds(clonedBonds);
      setHistory([]);
      setSelectedAtomId(null);
      setShowExportPanel(false);
      setShowHelpBanner(false);
    }
  }, [isOpen, initialAtoms, initialBonds]);

  useEffect(() => {
    setSelectedAtomId(null);
  }, [activeTool]);

  useEffect(() => {
    if (selectedAtomId && !draftAtoms.some((a) => a.id === selectedAtomId)) {
      setSelectedAtomId(null);
    }
  }, [draftAtoms, selectedAtomId]);

  // Keep top text box synchronized with current graph & format selection
  useEffect(() => {
    if (!isOpen) return;
    if (topFormat === 'SMILES') {
      setTopInputText(generateSmilesFromGraph(draftAtoms, draftBonds));
    } else if (topFormat === 'Formula') {
      setTopInputText(computeMolecularFormula(draftAtoms, draftBonds, true).plainFormula);
    } else if (topFormat === 'InChI') {
      const f = computeMolecularFormula(draftAtoms, draftBonds, true).plainFormula;
      setTopInputText(f && f !== 'Empty' ? `InChI=1S/${f}` : '');
    } else if (topFormat === 'MDL Molfile') {
      setTopInputText(generateSmilesFromGraph(draftAtoms, draftBonds));
    }
  }, [draftAtoms, draftBonds, topFormat, isOpen]);

  if (!isOpen) return null;

  const SVG_W = 540;
  const SVG_H = 420;
  const minX = draftAtoms.length > 0 ? Math.min(...draftAtoms.map((a) => a.x2d)) : -2.5;
  const maxX = draftAtoms.length > 0 ? Math.max(...draftAtoms.map((a) => a.x2d)) : 2.5;
  const minY = draftAtoms.length > 0 ? Math.min(...draftAtoms.map((a) => a.y2d)) : -2.0;
  const maxY = draftAtoms.length > 0 ? Math.max(...draftAtoms.map((a) => a.y2d)) : 2.0;
  const spanX = Math.max(5.5, maxX - minX + 2.6);
  const spanY = Math.max(4.5, maxY - minY + 2.4);
  const cx2d = draftAtoms.length > 0 ? (minX + maxX) / 2 : 0;
  const cy2d = draftAtoms.length > 0 ? (minY + maxY) / 2 : 0;
  const SCALE = Math.max(16, Math.min(64, Math.min((SVG_W - 64) / spanX, (SVG_H - 64) / spanY)));

  const toScreenX = (x2d: number) => SVG_W / 2 + (x2d - cx2d) * SCALE;
  const toScreenY = (y2d: number) => SVG_H / 2 + (y2d - cy2d) * SCALE;
  const toChemX = (px: number) => Number((cx2d + (px - SVG_W / 2) / SCALE).toFixed(2));
  const toChemY = (py: number) => Number((cy2d + (py - SVG_H / 2) / SCALE).toFixed(2));

  const pushState = (nextAtoms: AtomNode[], nextBonds: BondEdge[]) => {
    setHistory((prev) => [...prev.slice(-30), { atoms: draftAtoms, bonds: draftBonds }]);
    setDraftAtoms(nextAtoms);
    setDraftBonds(nextBonds);
  };

  const handleUndo = () => {
    if (history.length === 0) return;
    const prev = history[history.length - 1];
    setHistory((h) => h.slice(0, -1));
    setDraftAtoms(prev.atoms);
    setDraftBonds(prev.bonds);
  };

  const getNextAtomId = (el: ElementSymbol, currentAtoms: AtomNode[]) => {
    const nums = currentAtoms
      .filter((a) => a.element === el)
      .map((a) => parseInt(a.id.replace(/\D+/g, ''), 10) || 0);
    const nextNum = (nums.length > 0 ? Math.max(...nums) : 0) + 1;
    return `${el}${nextNum}`;
  };

  const selectElementTool = (el: ElementSymbol) => {
    setActiveElement(el);
    // If not currently in a bond placement mode, switch to element tool;
    // if already in bond mode, keep bond mode active so clicking empty space places `el`
    const isBondTool =
      activeTool === 'bond_single' ||
      activeTool === 'bond_double' ||
      activeTool === 'bond_triple' ||
      activeTool === 'bond_wedge' ||
      activeTool === 'bond_dash' ||
      activeTool === 'bond_wavy' ||
      activeTool === 'bond_crossed' ||
      activeTool === 'bond_any' ||
      activeTool === 'bond_sa' ||
      activeTool === 'bond_da' ||
      activeTool === 'bond_sd';
    if (!isBondTool) {
      setActiveTool('element');
    }
  };

  const handleCleanStructure = () => {
    if (draftAtoms.length === 0) return;
    const mx = draftAtoms.reduce((s, a) => s + a.x2d, 0) / draftAtoms.length;
    const my = draftAtoms.reduce((s, a) => s + a.y2d, 0) / draftAtoms.length;
    const cleaned = draftAtoms.map((a) => {
      const nx = Number((a.x2d - mx).toFixed(2));
      const ny = Number((a.y2d - my).toFixed(2));
      return { ...a, x2d: nx, y2d: ny, x3d: nx, y3d: -ny };
    });
    pushState(cleaned, draftBonds);
  };

  const handleRotateStructure = () => {
    if (draftAtoms.length === 0) return;
    const mx = draftAtoms.reduce((s, a) => s + a.x2d, 0) / draftAtoms.length;
    const my = draftAtoms.reduce((s, a) => s + a.y2d, 0) / draftAtoms.length;
    const rad = Math.PI / 6; // 30 degrees
    const c = Math.cos(rad);
    const s = Math.sin(rad);
    const rotated = draftAtoms.map((a) => {
      const dx = a.x2d - mx;
      const dy = a.y2d - my;
      const nx = Number((mx + dx * c - dy * s).toFixed(2));
      const ny = Number((my + dx * s + dy * c).toFixed(2));
      return { ...a, x2d: nx, y2d: ny, x3d: nx, y3d: -ny };
    });
    pushState(rotated, draftBonds);
  };

  const handleFlipVertical = () => {
    if (draftAtoms.length === 0) return;
    const my = draftAtoms.reduce((s, a) => s + a.y2d, 0) / draftAtoms.length;
    const flipped = draftAtoms.map((a) => {
      const ny = Number((2 * my - a.y2d).toFixed(2));
      return { ...a, y2d: ny, y3d: -ny };
    });
    pushState(flipped, draftBonds);
  };

  const handleFlipHorizontal = () => {
    if (draftAtoms.length === 0) return;
    const mx = draftAtoms.reduce((s, a) => s + a.x2d, 0) / draftAtoms.length;
    const flipped = draftAtoms.map((a) => {
      const nx = Number((2 * mx - a.x2d).toFixed(2));
      return { ...a, x2d: nx, x3d: nx };
    });
    pushState(flipped, draftBonds);
  };

  const handleHydrogenAction = (mode: 'Keep AsIs' | 'Add Implicit H' | 'Remove Explicit H') => {
    setHydrogenMode(mode);
    if (mode === 'Add Implicit H') {
      const expanded = expandImplicitHydrogens(draftAtoms, draftBonds, false);
      pushState(expanded.atoms, expanded.bonds);
    } else if (mode === 'Remove Explicit H') {
      const heavyAtoms = draftAtoms.filter((a) => a.element !== 'H');
      if (heavyAtoms.length > 0) {
        const heavySet = new Set(heavyAtoms.map((a) => a.id));
        const nextBonds = draftBonds.filter((b) => heavySet.has(b.source) && heavySet.has(b.target));
        pushState(heavyAtoms, nextBonds);
      }
    }
  };

  const handleExportClick = () => {
    if (exportFormat === 'MDL Molfile') {
      setExportText(generateMolfileFromGraph(draftAtoms, draftBonds));
    } else if (exportFormat === 'SMILES') {
      setExportText(generateSmilesFromGraph(draftAtoms, draftBonds));
    } else {
      setExportText(computeMolecularFormula(draftAtoms, draftBonds, true).plainFormula);
    }
    setShowExportPanel((v) => !v);
  };

  const handleFileImport = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;
    setImportedFileName(file.name);
    const reader = new FileReader();
    reader.onload = () => {
      const content = String(reader.result || '');
      const parsed = parseChemicalInputToGraph(
        content,
        file.name.endsWith('.mol') || file.name.endsWith('.sdf') || content.includes('V2000')
          ? 'MDL Molfile'
          : 'SMILES'
      );
      if (parsed && parsed.atoms.length > 0) {
        pushState(parsed.atoms, parsed.bonds);
      }
    };
    reader.readAsText(file);
  };

  const handleTopInputSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    const parsed = parseChemicalInputToGraph(topInputText, topFormat);
    if (parsed && parsed.atoms.length > 0) {
      pushState(parsed.atoms, parsed.bonds);
    }
  };

  // Helper to stamp a ring, spiro/bicyclic, chain, or functional group at (x, y) or attached to anchorAtom
  const stampStructureAt = (
    tool: PubChemTool,
    x: number,
    y: number,
    anchorAtom?: AtomNode
  ): boolean => {
    const newAtoms = [...draftAtoms];
    const newBonds = [...draftBonds];
    const addAtom = (el: ElementSymbol, ax: number, ay: number, formalCharge?: number): AtomNode => {
      const id = getNextAtomId(el, newAtoms);
      const node: AtomNode = {
        id,
        element: el,
        x2d: Number(ax.toFixed(2)),
        y2d: Number(ay.toFixed(2)),
        x3d: Number(ax.toFixed(2)),
        y3d: Number((-ay).toFixed(2)),
        z3d: 0,
        ...(formalCharge ? { formalCharge } : {}),
      };
      newAtoms.push(node);
      return node;
    };
    const addBond = (
      sId: string,
      tId: string,
      order: 1 | 2 | 3 = 1,
      stereo: BondEdge['stereo'] = 'none'
    ) => {
      newBonds.push({
        id: `b_pc_${Date.now()}_${newBonds.length}`,
        source: sId,
        target: tId,
        order,
        stereo,
      });
    };

    // Outward direction from anchorAtom if present
    let dirAngle = -Math.PI / 6;
    if (anchorAtom) {
      const connected = draftBonds.filter(
        (b) => b.source === anchorAtom.id || b.target === anchorAtom.id
      );
      if (connected.length > 0) {
        let vx = 0;
        let vy = 0;
        for (const b of connected) {
          const otherId = b.source === anchorAtom.id ? b.target : b.source;
          const other = draftAtoms.find((a) => a.id === otherId);
          if (other) {
            vx += anchorAtom.x2d - other.x2d;
            vy += anchorAtom.y2d - other.y2d;
          }
        }
        if (Math.hypot(vx, vy) > 0.05) {
          dirAngle = Math.atan2(vy, vx);
        } else {
          dirAngle = Math.PI / 3;
        }
      }
    }

    // 1. Rings (3..8 and Benzene)
    if (
      tool === 'ring_3' ||
      tool === 'ring_4' ||
      tool === 'ring_5' ||
      tool === 'ring_6' ||
      tool === 'ring_7' ||
      tool === 'ring_8' ||
      tool === 'ring_benzene'
    ) {
      const nSides =
        tool === 'ring_3'
          ? 3
          : tool === 'ring_4'
            ? 4
            : tool === 'ring_5'
              ? 5
              : tool === 'ring_7'
                ? 7
                : tool === 'ring_8'
                  ? 8
                  : 6;
      const isAromatic = tool === 'ring_benzene';
      const radius = nSides === 3 ? 0.85 : nSides === 4 ? 1.0 : nSides >= 7 ? 1.45 : 1.25;
      const ringEl = activeElement === 'H' ? 'C' : activeElement;
      const cx = anchorAtom ? anchorAtom.x2d + Math.cos(dirAngle) * (radius + 1.1) : x;
      const cy = anchorAtom ? anchorAtom.y2d + Math.sin(dirAngle) * (radius + 1.1) : y;
      const ringNodes: AtomNode[] = [];

      for (let k = 0; k < nSides; k++) {
        const ang = dirAngle + Math.PI + (k * 2 * Math.PI) / nSides;
        ringNodes.push(
          addAtom(ringEl, cx + Math.cos(ang) * radius, cy + Math.sin(ang) * radius)
        );
      }
      for (let k = 0; k < nSides; k++) {
        const nextK = (k + 1) % nSides;
        const order: 1 | 2 = isAromatic && k % 2 === 0 ? 2 : 1;
        addBond(ringNodes[k].id, ringNodes[nextK].id, order);
      }
      if (anchorAtom) {
        addBond(anchorAtom.id, ringNodes[0].id, 1);
      }
      pushState(newAtoms, newBonds);
      return true;
    }

    // 2. Bicyclic / Spiro ring system (Naphthalene-like fused bicyclic)
    if (tool === 'ring_bicyclic') {
      const ringEl = activeElement === 'H' ? 'C' : activeElement;
      const r = 1.15;
      const cx1 = (anchorAtom ? anchorAtom.x2d + 1.8 : x) - r * 0.866;
      const cy1 = anchorAtom ? anchorAtom.y2d : y;
      const cx2 = cx1 + r * 1.732;
      const ring1: AtomNode[] = [];
      for (let k = 0; k < 6; k++) {
        const ang = Math.PI / 6 + (k * Math.PI) / 3;
        ring1.push(addAtom(ringEl, cx1 + Math.cos(ang) * r, cy1 + Math.sin(ang) * r));
      }
      for (let k = 0; k < 6; k++) {
        addBond(ring1[k].id, ring1[(k + 1) % 6].id, 1);
      }
      // Second fused 6-ring sharing ring1[0] and ring1[5]
      const r2a = addAtom(ringEl, cx2 + Math.cos(Math.PI / 6) * r, cy1 + Math.sin(Math.PI / 6) * r);
      const r2b = addAtom(ringEl, cx2 + Math.cos(-Math.PI / 6) * r, cy1 + Math.sin(-Math.PI / 6) * r);
      const r2c = addAtom(ringEl, cx2 + Math.cos(-Math.PI / 2) * r, cy1 + Math.sin(-Math.PI / 2) * r);
      const r2d = addAtom(ringEl, cx2 + Math.cos(Math.PI / 2) * r, cy1 + Math.sin(Math.PI / 2) * r);
      addBond(ring1[0].id, r2d.id, 1);
      addBond(r2d.id, r2a.id, 1);
      addBond(r2a.id, r2b.id, 1);
      addBond(r2b.id, r2c.id, 1);
      addBond(r2c.id, ring1[5].id, 1);
      if (anchorAtom) addBond(anchorAtom.id, ring1[2].id, 1);
      pushState(newAtoms, newBonds);
      return true;
    }

    // 3. Chains & Branched Chains
    if (
      tool === 'chain_ethyl' ||
      tool === 'chain_propyl' ||
      tool === 'chain_butyl'
    ) {
      const nCarbons = tool === 'chain_ethyl' ? 2 : tool === 'chain_propyl' ? 3 : 4;
      let prev = anchorAtom || addAtom('C', x - 1.2, y);
      for (let k = 0; k < nCarbons; k++) {
        const stepAng = dirAngle + (k % 2 === 0 ? -Math.PI / 6 : Math.PI / 6);
        const nextNode = addAtom(
          'C',
          prev.x2d + Math.cos(stepAng) * 1.25,
          prev.y2d + Math.sin(stepAng) * 1.25
        );
        addBond(prev.id, nextNode.id, 1);
        prev = nextNode;
      }
      pushState(newAtoms, newBonds);
      return true;
    }

    if (tool === 'chain_ipr' || tool === 'chain_tbu' || tool === 'chain_neo' || tool === 'chain_cross') {
      const root =
        anchorAtom
          ? addAtom('C', anchorAtom.x2d + Math.cos(dirAngle) * 1.25, anchorAtom.y2d + Math.sin(dirAngle) * 1.25)
          : addAtom('C', x, y);
      if (anchorAtom) addBond(anchorAtom.id, root.id, 1);

      const nBranches = tool === 'chain_ipr' ? 2 : tool === 'chain_cross' && !anchorAtom ? 4 : 3;
      for (let b = 0; b < nBranches; b++) {
        const ang =
          nBranches === 4
            ? (b * Math.PI) / 2
            : dirAngle + ((b - (nBranches - 1) / 2) * Math.PI) / 2.4;
        const leaf = addAtom('C', root.x2d + Math.cos(ang) * 1.22, root.y2d + Math.sin(ang) * 1.22);
        addBond(root.id, leaf.id, 1);
      }
      pushState(newAtoms, newBonds);
      return true;
    }

    // 4. Functional Groups: CHO, CO2H, NO2, SO3H
    if (tool === 'fg_cho') {
      const cNode = anchorAtom
        ? addAtom('C', anchorAtom.x2d + Math.cos(dirAngle) * 1.25, anchorAtom.y2d + Math.sin(dirAngle) * 1.25)
        : addAtom('C', x, y);
      if (anchorAtom) addBond(anchorAtom.id, cNode.id, 1);
      const oNode = addAtom(
        'O',
        cNode.x2d + Math.cos(dirAngle - Math.PI / 3) * 1.2,
        cNode.y2d + Math.sin(dirAngle - Math.PI / 3) * 1.2
      );
      const hNode = addAtom(
        'H',
        cNode.x2d + Math.cos(dirAngle + Math.PI / 3) * 1.0,
        cNode.y2d + Math.sin(dirAngle + Math.PI / 3) * 1.0
      );
      addBond(cNode.id, oNode.id, 2);
      addBond(cNode.id, hNode.id, 1);
      pushState(newAtoms, newBonds);
      return true;
    }

    if (tool === 'fg_co2h') {
      const cNode = anchorAtom
        ? addAtom('C', anchorAtom.x2d + Math.cos(dirAngle) * 1.25, anchorAtom.y2d + Math.sin(dirAngle) * 1.25)
        : addAtom('C', x, y);
      if (anchorAtom) addBond(anchorAtom.id, cNode.id, 1);
      const oDouble = addAtom(
        'O',
        cNode.x2d + Math.cos(dirAngle - Math.PI / 3) * 1.2,
        cNode.y2d + Math.sin(dirAngle - Math.PI / 3) * 1.2
      );
      const oSingle = addAtom(
        'O',
        cNode.x2d + Math.cos(dirAngle + Math.PI / 3) * 1.25,
        cNode.y2d + Math.sin(dirAngle + Math.PI / 3) * 1.25
      );
      const hNode = addAtom(
        'H',
        oSingle.x2d + Math.cos(dirAngle) * 0.95,
        oSingle.y2d + Math.sin(dirAngle) * 0.95
      );
      addBond(cNode.id, oDouble.id, 2);
      addBond(cNode.id, oSingle.id, 1);
      addBond(oSingle.id, hNode.id, 1);
      pushState(newAtoms, newBonds);
      return true;
    }

    if (tool === 'fg_no2') {
      const nNode = anchorAtom
        ? addAtom('N', anchorAtom.x2d + Math.cos(dirAngle) * 1.25, anchorAtom.y2d + Math.sin(dirAngle) * 1.25, 1)
        : addAtom('N', x, y, 1);
      if (anchorAtom) addBond(anchorAtom.id, nNode.id, 1);
      const o1 = addAtom(
        'O',
        nNode.x2d + Math.cos(dirAngle - Math.PI / 3) * 1.2,
        nNode.y2d + Math.sin(dirAngle - Math.PI / 3) * 1.2
      );
      const o2 = addAtom(
        'O',
        nNode.x2d + Math.cos(dirAngle + Math.PI / 3) * 1.2,
        nNode.y2d + Math.sin(dirAngle + Math.PI / 3) * 1.2,
        -1
      );
      addBond(nNode.id, o1.id, 2);
      addBond(nNode.id, o2.id, 1);
      pushState(newAtoms, newBonds);
      return true;
    }

    if (tool === 'fg_so3h') {
      const sNode = anchorAtom
        ? addAtom('S', anchorAtom.x2d + Math.cos(dirAngle) * 1.35, anchorAtom.y2d + Math.sin(dirAngle) * 1.35)
        : addAtom('S', x, y);
      if (anchorAtom) addBond(anchorAtom.id, sNode.id, 1);
      const o1 = addAtom(
        'O',
        sNode.x2d + Math.cos(dirAngle - Math.PI / 2) * 1.22,
        sNode.y2d + Math.sin(dirAngle - Math.PI / 2) * 1.22
      );
      const o2 = addAtom(
        'O',
        sNode.x2d + Math.cos(dirAngle + Math.PI / 2) * 1.22,
        sNode.y2d + Math.sin(dirAngle + Math.PI / 2) * 1.22
      );
      const o3 = addAtom(
        'O',
        sNode.x2d + Math.cos(dirAngle) * 1.28,
        sNode.y2d + Math.sin(dirAngle) * 1.28
      );
      const hNode = addAtom(
        'H',
        o3.x2d + Math.cos(dirAngle) * 0.95,
        o3.y2d + Math.sin(dirAngle) * 0.95
      );
      addBond(sNode.id, o1.id, 2);
      addBond(sNode.id, o2.id, 2);
      addBond(sNode.id, o3.id, 1);
      addBond(o3.id, hNode.id, 1);
      pushState(newAtoms, newBonds);
      return true;
    }

    return false;
  };

  const getBondSpecFromTool = (
    tool: PubChemTool
  ): { order: 1 | 2 | 3; stereo: BondEdge['stereo'] } | null => {
    switch (tool) {
      case 'bond_single':
      case 'bond_sa':
        return { order: 1, stereo: 'none' };
      case 'bond_double':
      case 'bond_da':
      case 'bond_sd':
        return { order: 2, stereo: 'none' };
      case 'bond_triple':
        return { order: 3, stereo: 'none' };
      case 'bond_wedge':
        return { order: 1, stereo: 'wedge' };
      case 'bond_dash':
        return { order: 1, stereo: 'dash' };
      case 'bond_wavy':
      case 'bond_any':
        return { order: 1, stereo: 'wavy' };
      case 'bond_crossed':
        return { order: 2, stereo: 'crossed' };
      default:
        return null;
    }
  };

  const getSvgCoords = (e: React.MouseEvent<SVGSVGElement>): { x: number; y: number } => {
    if (!svgRef.current) return { x: 0, y: 0 };
    const rect = svgRef.current.getBoundingClientRect();
    const px = ((e.clientX - rect.left) / (rect.width || SVG_W)) * SVG_W;
    const py = ((e.clientY - rect.top) / (rect.height || SVG_H)) * SVG_H;
    return { x: toChemX(px), y: toChemY(py) };
  };

  const handleCanvasClick = (e: React.MouseEvent<SVGSVGElement>) => {
    if (e.target !== svgRef.current && (e.target as Element).tagName !== 'rect') return;
    const { x, y } = getSvgCoords(e);

    if (stampStructureAt(activeTool, x, y)) {
      return;
    }

    if (activeTool === 'element') {
      const id = getNextAtomId(activeElement, draftAtoms);
      const newAtom: AtomNode = {
        id,
        element: activeElement,
        x2d: x,
        y2d: y,
        x3d: x,
        y3d: -y,
        z3d: 0,
      };
      setSelectedAtomId(null);
      pushState([...draftAtoms, newAtom], draftBonds);
      return;
    }

    const bondSpec = getBondSpecFromTool(activeTool);
    if (bondSpec) {
      const id1 = getNextAtomId(activeElement, draftAtoms);
      const newAtom: AtomNode = {
        id: id1,
        element: activeElement,
        x2d: x,
        y2d: y,
        x3d: x,
        y3d: -y,
        z3d: bondSpec.stereo === 'wedge' ? 0.65 : bondSpec.stereo === 'dash' ? -0.65 : 0,
      };
      const sourceAtom = selectedAtomId
        ? draftAtoms.find((a) => a.id === selectedAtomId)
        : undefined;

      if (sourceAtom) {
        // First click on atom, second click on empty canvas space:
        // Add an atom of the selected element at (x, y) and connect with first atom using the selected bond
        const newBond: BondEdge = {
          id: `b_pc_${Date.now()}`,
          source: sourceAtom.id,
          target: id1,
          order: bondSpec.order,
          stereo: bondSpec.stereo,
        };
        setSelectedAtomId(null);
        pushState([...draftAtoms, newAtom], [...draftBonds, newBond]);
      } else {
        // No first atom selected: place standalone atom of the selected element at (x, y)
        setSelectedAtomId(null);
        pushState([...draftAtoms, newAtom], draftBonds);
      }
    }
  };

  const handleAtomClick = (e: React.MouseEvent, atom: AtomNode) => {
    e.stopPropagation();

    if (activeTool === 'erase') {
      if (selectedAtomId === atom.id) setSelectedAtomId(null);
      pushState(
        draftAtoms.filter((a) => a.id !== atom.id),
        draftBonds.filter((b) => b.source !== atom.id && b.target !== atom.id)
      );
      return;
    }

    if (activeTool === 'charge_plus') {
      pushState(
        draftAtoms.map((a) =>
          a.id === atom.id ? { ...a, formalCharge: Math.min(3, (a.formalCharge ?? 0) + 1) } : a
        ),
        draftBonds
      );
      return;
    }

    if (activeTool === 'charge_minus') {
      pushState(
        draftAtoms.map((a) =>
          a.id === atom.id ? { ...a, formalCharge: Math.max(-3, (a.formalCharge ?? 0) - 1) } : a
        ),
        draftBonds
      );
      return;
    }

    if (activeTool === 'charge_neutral') {
      pushState(
        draftAtoms.map((a) => (a.id === atom.id ? { ...a, formalCharge: 0 } : a)),
        draftBonds
      );
      return;
    }

    if (activeTool === 'element' || activeTool === 'query') {
      const targetEl = activeTool === 'query' ? queryEl : activeElement;
      const newId = getNextAtomId(targetEl, draftAtoms.filter((a) => a.id !== atom.id));
      const updatedAtoms = draftAtoms.map((a) =>
        a.id === atom.id ? { ...a, id: newId, element: targetEl } : a
      );
      const updatedBonds = draftBonds.map((b) => ({
        ...b,
        source: b.source === atom.id ? newId : b.source,
        target: b.target === atom.id ? newId : b.target,
      }));
      pushState(updatedAtoms, updatedBonds);
      return;
    }

    if (stampStructureAt(activeTool, atom.x2d, atom.y2d, atom)) {
      return;
    }

    const bondSpec = getBondSpecFromTool(activeTool);
    if (bondSpec) {
      // First click: select the first atom
      if (!selectedAtomId) {
        setSelectedAtomId(atom.id);
        return;
      }

      const sourceAtom = draftAtoms.find((a) => a.id === selectedAtomId);
      if (!sourceAtom) {
        setSelectedAtomId(atom.id);
        return;
      }

      // Second click is another atom: add or update bond between the 2 selected atoms
      if (atom.id !== sourceAtom.id) {
        const existingBond = draftBonds.find(
          (b) =>
            (b.source === sourceAtom.id && b.target === atom.id) ||
            (b.source === atom.id && b.target === sourceAtom.id)
        );

        if (existingBond) {
          const updatedBonds = draftBonds.map((b) =>
            b.id === existingBond.id
              ? {
                  ...b,
                  source: sourceAtom.id,
                  target: atom.id,
                  order: bondSpec.order,
                  stereo: bondSpec.stereo,
                }
              : b
          );
          setSelectedAtomId(null);
          pushState(draftAtoms, updatedBonds);
        } else {
          const newBond: BondEdge = {
            id: `b_pc_${Date.now()}`,
            source: sourceAtom.id,
            target: atom.id,
            order: bondSpec.order,
            stereo: bondSpec.stereo,
          };
          setSelectedAtomId(null);
          pushState(draftAtoms, [...draftBonds, newBond]);
        }
        return;
      }

      // Second click is the same atom: replace the current atom with the selected element
      const newId =
        atom.element === activeElement
          ? atom.id
          : getNextAtomId(
              activeElement,
              draftAtoms.filter((a) => a.id !== atom.id)
            );
      const updatedAtoms = draftAtoms.map((a) =>
        a.id === atom.id ? { ...a, id: newId, element: activeElement } : a
      );
      const updatedBonds = draftBonds.map((b) => ({
        ...b,
        source: b.source === atom.id ? newId : b.source,
        target: b.target === atom.id ? newId : b.target,
      }));
      setSelectedAtomId(null);
      pushState(updatedAtoms, updatedBonds);
    }
  };

  const handleBondClick = (e: React.MouseEvent, bond: BondEdge) => {
    e.stopPropagation();
    setSelectedAtomId(null);
    if (activeTool === 'erase') {
      pushState(
        draftAtoms,
        draftBonds.filter((b) => b.id !== bond.id)
      );
      return;
    }
    const bondSpec = getBondSpecFromTool(activeTool);
    const updatedBonds = draftBonds.map((b) => {
      if (b.id !== bond.id) return b;
      if (bondSpec) {
        return { ...b, order: bondSpec.order, stereo: bondSpec.stereo };
      }
      const nextOrder: 1 | 2 | 3 = b.order === 1 ? 2 : b.order === 2 ? 3 : 1;
      return { ...b, order: nextOrder, stereo: 'none' as const };
    });
    pushState(draftAtoms, updatedBonds);
  };

  const handleMouseMove = (e: React.MouseEvent<SVGSVGElement>) => {
    const coords = getSvgCoords(e);
    setCursorChemPos(coords);
    if (draggingAtomId) {
      const { x, y } = coords;
      setDraftAtoms((prev) =>
        prev.map((a) =>
          a.id === draggingAtomId ? { ...a, x2d: x, y2d: y, x3d: x, y3d: -y } : a
        )
      );
    } else if (draggingCanvas) {
      const { x, y } = coords;
      const dx = Number((x - draggingCanvas.startX).toFixed(2));
      const dy = Number((y - draggingCanvas.startY).toFixed(2));
      if (Math.abs(dx) > 0.05 || Math.abs(dy) > 0.05) {
        setDraftAtoms((prev) =>
          prev.map((a) => ({
            ...a,
            x2d: Number((a.x2d + dx).toFixed(2)),
            y2d: Number((a.y2d + dy).toFixed(2)),
            x3d: Number((a.x2d + dx).toFixed(2)),
            y3d: Number((-(a.y2d + dy)).toFixed(2)),
          }))
        );
        setDraggingCanvas({ startX: x, startY: y });
      }
    }
  };

  const toolBtnClass = (tool: PubChemTool) =>
    `h-7 flex items-center justify-center text-[11px] font-medium border border-slate-400 select-none transition-colors ${
      activeTool === tool
        ? 'bg-amber-400 text-slate-950 font-bold shadow-inner'
        : 'bg-gradient-to-b from-slate-100 to-slate-200 text-slate-900 hover:from-slate-200 hover:to-slate-300'
    }`;

  const elemBtnClass = (el: ElementSymbol) =>
    `h-7 flex items-center justify-center text-xs font-semibold border border-slate-400 select-none transition-colors ${
      activeElement === el
        ? 'bg-amber-400 text-slate-950 font-bold shadow-inner'
        : 'bg-gradient-to-b from-[#e8f0f8] to-[#d0e0f0] text-slate-900 hover:from-[#d8e8f8] hover:to-[#bfd4ec]'
    }`;

  return (
    <div
      className="fixed inset-0 z-50 bg-slate-900/60 backdrop-blur-[1px] flex items-center justify-center p-3 overflow-y-auto"
      role="dialog"
      aria-modal="true"
      aria-label="Draw Structure - Full 2D Chemical Sketcher"
    >
      <div className="bg-[#d6d6d6] border border-slate-500 shadow-2xl w-full max-w-[980px] rounded-xs overflow-hidden text-slate-900">
        {/* Top Modal Header Bar */}
        <div className="px-4 py-2.5 flex items-center justify-between bg-[#d6d6d6]">
          <h2 className="text-base font-normal tracking-wide text-slate-800 uppercase">
            DRAW STRUCTURE
          </h2>
          <button
            type="button"
            onClick={onCancel}
            aria-label="Close Draw Structure Modal"
            className="text-slate-700 hover:text-slate-950 text-xl leading-none px-1.5 py-0.5 font-semibold"
          >
            ✕
          </button>
        </div>

        {/* Inner Sketcher Frame */}
        <div className="px-4 pb-3">
          <div className="border border-slate-500 bg-[#e4ecf5]">
            {/* Top Format & Live String Bar */}
            <form
              onSubmit={handleTopInputSubmit}
              className="grid grid-cols-[120px_110px_1fr] border-b border-slate-500 bg-[#dce6f2]"
            >
              <div className="p-0.5 border-r border-slate-400">
                <select
                  aria-label="Sketcher Network Mode"
                  value={networkMode}
                  onChange={(e) => setNetworkMode(e.target.value as 'Broadband' | 'Standard')}
                  className="w-full h-6 bg-white border border-slate-400 text-xs px-1 text-slate-900"
                >
                  <option value="Broadband">Broadband</option>
                  <option value="Standard">Standard</option>
                </select>
              </div>
              <div className="p-0.5 border-r border-slate-400">
                <select
                  aria-label="Chemical String Format"
                  value={topFormat}
                  onChange={(e) =>
                    setTopFormat(
                      e.target.value as 'SMILES' | 'Formula' | 'MDL Molfile' | 'InChI'
                    )
                  }
                  className="w-full h-6 bg-white border border-slate-400 text-xs px-1 text-slate-900"
                >
                  <option value="SMILES">SMILES</option>
                  <option value="Formula">Formula</option>
                  <option value="MDL Molfile">MDL Molfile</option>
                  <option value="InChI">InChI</option>
                </select>
              </div>
              <div className="p-0.5 flex items-center gap-1 bg-white">
                <input
                  type="text"
                  aria-label="Live SMILES or Formula Input"
                  value={topInputText}
                  onChange={(e) => setTopInputText(e.target.value)}
                  placeholder="Type SMILES or Formula and press Enter..."
                  className="w-full h-6 px-2 text-xs font-mono text-slate-900 bg-white focus:outline-none"
                />
                <button
                  type="submit"
                  title="Parse SMILES/Formula into Canvas"
                  className="px-2 h-5 text-[10px] font-semibold bg-slate-200 hover:bg-slate-300 border border-slate-400 text-slate-800 mr-0.5 shrink-0"
                >
                  Apply
                </button>
              </div>
            </form>

            {/* Main 2-Column Body: Left Tool/Periodic Grid + Right Interactive White Canvas */}
            <div className="grid grid-cols-1 md:grid-cols-[332px_1fr]">
              {/* Left Tool & Periodic Table Panel */}
              <div className="border-r border-slate-500 bg-[#dce6f2] select-none">
                {/* Row 1: Action & Transform Tools (New, Udo, Cln, Sfy, Del, Qry, Move, Rotate, FlipV, FlipH) */}
                <div className="grid grid-cols-11">
                  <button
                    type="button"
                    onClick={() => pushState([], [])}
                    title="New / Clear Canvas"
                    className="col-span-2 h-7 text-[11px] font-medium border border-slate-400 bg-gradient-to-b from-slate-100 to-slate-200 hover:from-slate-200 hover:to-slate-300"
                  >
                    New
                  </button>
                  <button
                    type="button"
                    onClick={handleUndo}
                    disabled={history.length === 0}
                    title="Undo Last Edit"
                    className="h-7 text-[11px] font-medium border border-slate-400 bg-gradient-to-b from-slate-100 to-slate-200 hover:from-slate-200 hover:to-slate-300 disabled:opacity-50"
                  >
                    Udo
                  </button>
                  <button
                    type="button"
                    onClick={handleCleanStructure}
                    title="Clean & Center 2D Coordinates"
                    className="h-7 text-[11px] font-medium border border-slate-400 bg-gradient-to-b from-slate-100 to-slate-200 hover:from-slate-200 hover:to-slate-300"
                  >
                    Cln
                  </button>
                  <button
                    type="button"
                    onClick={handleCleanStructure}
                    title="Simplify & Standardize Structure"
                    className="h-7 text-[11px] font-medium border border-slate-400 bg-gradient-to-b from-slate-100 to-slate-200 hover:from-slate-200 hover:to-slate-300"
                  >
                    Sfy
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('erase')}
                    title="Delete Atom or Bond"
                    className={toolBtnClass('erase')}
                  >
                    Del
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('query')}
                    title="Query Atom Tool"
                    className={toolBtnClass('query')}
                  >
                    Qry
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('move')}
                    title="Translate / Move Atom or Structure"
                    className={toolBtnClass('move')}
                  >
                    ✥
                  </button>
                  <button
                    type="button"
                    onClick={handleRotateStructure}
                    title="Rotate Structure by 30°"
                    className="h-7 text-xs font-bold border border-slate-400 bg-gradient-to-b from-slate-100 to-slate-200 hover:from-slate-200 hover:to-slate-300"
                  >
                    ↻
                  </button>
                  <button
                    type="button"
                    onClick={handleFlipVertical}
                    title="Flip Structure Vertically"
                    className="h-7 text-xs font-bold border border-slate-400 bg-gradient-to-b from-slate-100 to-slate-200 hover:from-slate-200 hover:to-slate-300"
                  >
                    ⇕
                  </button>
                  <button
                    type="button"
                    onClick={handleFlipHorizontal}
                    title="Flip Structure Horizontally"
                    className="h-7 text-xs font-bold border border-slate-400 bg-gradient-to-b from-slate-100 to-slate-200 hover:from-slate-200 hover:to-slate-300"
                  >
                    ⇔
                  </button>
                </div>

                {/* Row 2: 11 Bond Tools */}
                <div className="grid grid-cols-11">
                  <button
                    type="button"
                    onClick={() => setActiveTool('bond_single')}
                    title="Single Bond"
                    className={toolBtnClass('bond_single')}
                  >
                    <span className="font-bold text-sm">—</span>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('bond_double')}
                    title="Double Bond"
                    className={toolBtnClass('bond_double')}
                  >
                    <span className="font-bold text-sm">═</span>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('bond_triple')}
                    title="Triple Bond"
                    className={toolBtnClass('bond_triple')}
                  >
                    <span className="font-bold text-sm">≡</span>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('bond_wedge')}
                    title="Up Stereo Wedge Bond (+z)"
                    className={toolBtnClass('bond_wedge')}
                  >
                    <svg viewBox="0 0 20 12" className="w-4 h-3 fill-current">
                      <polygon points="2,6 18,2 18,10" />
                    </svg>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('bond_dash')}
                    title="Down Stereo Hashed Bond (-z)"
                    className={toolBtnClass('bond_dash')}
                  >
                    <svg viewBox="0 0 20 12" className="w-4 h-3 stroke-current" strokeWidth="1.8">
                      <line x1="3" y1="5" x2="3" y2="7" />
                      <line x1="7" y1="4" x2="7" y2="8" />
                      <line x1="11" y1="3" x2="11" y2="9" />
                      <line x1="15" y1="2" x2="15" y2="10" />
                    </svg>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('bond_wavy')}
                    title="Wavy / Unspecified Stereo Bond"
                    className={toolBtnClass('bond_wavy')}
                  >
                    <span className="text-[11px] font-bold">∿</span>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('bond_crossed')}
                    title="Cis/Trans Crossed Double Bond"
                    className={toolBtnClass('bond_crossed')}
                  >
                    <span className="text-[11px] font-bold">✕</span>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('bond_any')}
                    title="Any Bond Query"
                    className={toolBtnClass('bond_any')}
                  >
                    <span className="text-[11px] font-bold">⋀</span>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('bond_sa')}
                    title="Single or Aromatic Bond"
                    className={toolBtnClass('bond_sa')}
                  >
                    <span className="text-[10px]">S/A</span>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('bond_da')}
                    title="Double or Aromatic Bond"
                    className={toolBtnClass('bond_da')}
                  >
                    <span className="text-[10px]">D/A</span>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('bond_sd')}
                    title="Single or Double Bond"
                    className={toolBtnClass('bond_sd')}
                  >
                    <span className="text-[10px]">S/D</span>
                  </button>
                </div>

                {/* Row 3: 8 Ring & Spiro Templates + 3 Charge Tools */}
                <div className="grid grid-cols-11">
                  <button
                    type="button"
                    onClick={() => setActiveTool('ring_3')}
                    title="Cyclopropane (3-Ring)"
                    className={toolBtnClass('ring_3')}
                  >
                    <svg viewBox="0 0 20 20" className="w-3.5 h-3.5 stroke-current fill-none" strokeWidth="1.6">
                      <polygon points="10,3 17,16 3,16" />
                    </svg>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('ring_4')}
                    title="Cyclobutane (4-Ring)"
                    className={toolBtnClass('ring_4')}
                  >
                    <svg viewBox="0 0 20 20" className="w-3.5 h-3.5 stroke-current fill-none" strokeWidth="1.6">
                      <rect x="4" y="4" width="12" height="12" />
                    </svg>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('ring_5')}
                    title="Cyclopentane (5-Ring)"
                    className={toolBtnClass('ring_5')}
                  >
                    <svg viewBox="0 0 20 20" className="w-3.5 h-3.5 stroke-current fill-none" strokeWidth="1.6">
                      <polygon points="10,2.5 17.5,8 14.5,16.8 5.5,16.8 2.5,8" />
                    </svg>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('ring_6')}
                    title="Cyclohexane (6-Ring)"
                    className={toolBtnClass('ring_6')}
                  >
                    <svg viewBox="0 0 20 20" className="w-3.5 h-3.5 stroke-current fill-none" strokeWidth="1.6">
                      <polygon points="10,2.2 16.8,6.1 16.8,13.9 10,17.8 3.2,13.9 3.2,6.1" />
                    </svg>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('ring_7')}
                    title="Cycloheptane (7-Ring)"
                    className={toolBtnClass('ring_7')}
                  >
                    <svg viewBox="0 0 20 20" className="w-3.5 h-3.5 stroke-current fill-none" strokeWidth="1.5">
                      <polygon points="10,2 16.2,5 17.8,11.8 13.5,17.2 6.5,17.2 2.2,11.8 3.8,5" />
                    </svg>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('ring_8')}
                    title="Cyclooctane (8-Ring)"
                    className={toolBtnClass('ring_8')}
                  >
                    <svg viewBox="0 0 20 20" className="w-3.5 h-3.5 stroke-current fill-none" strokeWidth="1.5">
                      <polygon points="7,2.2 13,2.2 17.8,7 17.8,13 13,17.8 7,17.8 2.2,13 2.2,7" />
                    </svg>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('ring_benzene')}
                    title="Benzene (Aromatic 6-Ring)"
                    className={toolBtnClass('ring_benzene')}
                  >
                    <svg viewBox="0 0 20 20" className="w-3.5 h-3.5 stroke-current fill-none" strokeWidth="1.5">
                      <polygon points="10,2.2 16.8,6.1 16.8,13.9 10,17.8 3.2,13.9 3.2,6.1" />
                      <circle cx="10" cy="10" r="3.8" />
                    </svg>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('ring_bicyclic')}
                    title="Fused Bicyclic / Spiro Ring System"
                    className={toolBtnClass('ring_bicyclic')}
                  >
                    <svg viewBox="0 0 20 20" className="w-3.5 h-3.5 stroke-current fill-none" strokeWidth="1.4">
                      <rect x="2.5" y="2.5" width="15" height="15" />
                      <line x1="7.5" y1="2.5" x2="7.5" y2="17.5" />
                      <line x1="12.5" y1="2.5" x2="12.5" y2="17.5" />
                      <line x1="2.5" y1="7.5" x2="17.5" y2="7.5" />
                      <line x1="2.5" y1="12.5" x2="17.5" y2="12.5" />
                    </svg>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('charge_plus')}
                    title="Increment Formal Positive Charge (+1)"
                    className={toolBtnClass('charge_plus')}
                  >
                    ⊕
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('charge_minus')}
                    title="Decrement Formal Negative Charge (-1)"
                    className={toolBtnClass('charge_minus')}
                  >
                    ⊖
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('charge_neutral')}
                    title="Neutralize Formal Charge / Radical (0)"
                    className={toolBtnClass('charge_neutral')}
                  >
                    ⊙
                  </button>
                </div>

                {/* Row 4: 7 Alkyl / Branched Chains + 4 Functional Groups (CHO, CO2H, NO2, SO3H) */}
                <div className="grid grid-cols-11">
                  <button
                    type="button"
                    onClick={() => setActiveTool('chain_ethyl')}
                    title="Ethyl Chain (C2)"
                    className={toolBtnClass('chain_ethyl')}
                  >
                    <svg viewBox="0 0 20 14" className="w-4 h-3 stroke-current fill-none" strokeWidth="1.6">
                      <polyline points="2,11 10,4 18,11" />
                    </svg>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('chain_propyl')}
                    title="Propyl Chain (C3)"
                    className={toolBtnClass('chain_propyl')}
                  >
                    <svg viewBox="0 0 20 14" className="w-4 h-3 stroke-current fill-none" strokeWidth="1.6">
                      <polyline points="2,10 7.5,4 12.5,10 18,4" />
                    </svg>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('chain_ipr')}
                    title="Isopropyl Branch (i-Pr)"
                    className={toolBtnClass('chain_ipr')}
                  >
                    <svg viewBox="0 0 20 16" className="w-4 h-3.5 stroke-current fill-none" strokeWidth="1.6">
                      <line x1="10" y1="2" x2="10" y2="9" />
                      <line x1="10" y1="9" x2="3" y2="14" />
                      <line x1="10" y1="9" x2="17" y2="14" />
                    </svg>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('chain_butyl')}
                    title="Butyl Zigzag Chain (C4)"
                    className={toolBtnClass('chain_butyl')}
                  >
                    <svg viewBox="0 0 22 14" className="w-4 h-3 stroke-current fill-none" strokeWidth="1.6">
                      <polyline points="1.5,10 6,4 11,10 16,4 20.5,10" />
                    </svg>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('chain_tbu')}
                    title="tert-Butyl Group (t-Bu)"
                    className={toolBtnClass('chain_tbu')}
                  >
                    <svg viewBox="0 0 20 16" className="w-4 h-3.5 stroke-current fill-none" strokeWidth="1.6">
                      <polyline points="2,9 8,5 14,9" />
                      <line x1="8" y1="5" x2="13" y2="2" />
                      <line x1="14" y1="9" x2="14" y2="15" />
                    </svg>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('chain_neo')}
                    title="Neopentyl / Branched Group"
                    className={toolBtnClass('chain_neo')}
                  >
                    <svg viewBox="0 0 20 16" className="w-4 h-3.5 stroke-current fill-none" strokeWidth="1.6">
                      <polyline points="2,13 8,9 14,13" />
                      <line x1="8" y1="9" x2="8" y2="2" />
                      <line x1="8" y1="9" x2="15" y2="5" />
                    </svg>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('chain_cross')}
                    title="Tetrahedral Quaternary Carbon C(C)4"
                    className={toolBtnClass('chain_cross')}
                  >
                    <svg viewBox="0 0 20 20" className="w-3.5 h-3.5 stroke-current fill-none" strokeWidth="1.6">
                      <line x1="10" y1="2" x2="10" y2="18" />
                      <line x1="2" y1="10" x2="18" y2="10" />
                    </svg>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('fg_cho')}
                    title="Aldehyde Group (-CHO)"
                    className={toolBtnClass('fg_cho')}
                  >
                    <span className="text-[9px] font-bold">CHO</span>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('fg_co2h')}
                    title="Carboxylic Acid Group (-CO2H)"
                    className={toolBtnClass('fg_co2h')}
                  >
                    <span className="text-[9px] font-bold">CO₂H</span>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('fg_no2')}
                    title="Nitro Group (-NO2)"
                    className={toolBtnClass('fg_no2')}
                  >
                    <span className="text-[9px] font-bold">NO₂</span>
                  </button>
                  <button
                    type="button"
                    onClick={() => setActiveTool('fg_so3h')}
                    title="Sulfonic Acid Group (-SO3H)"
                    className={toolBtnClass('fg_so3h')}
                  >
                    <span className="text-[9px] font-bold">SO₃H</span>
                  </button>
                </div>

                {/* Periodic Table Rows (Periods 1 to 6) */}
                {/* Period 1: H, [empty], ?, [? dropdown], ... He */}
                <div className="grid grid-cols-11">
                  <button
                    type="button"
                    onClick={() => selectElementTool('H')}
                    className={elemBtnClass('H')}
                  >
                    H
                  </button>
                  <div className="h-7 border border-slate-400 bg-[#e8eff7]" />
                  <button
                    type="button"
                    onClick={() => selectElementTool(queryEl)}
                    title={`Select ${queryEl}`}
                    className={elemBtnClass(queryEl)}
                  >
                    ?
                  </button>
                  <div className="col-span-2 h-7 border border-slate-400 bg-[#e8eff7] flex items-center px-0.5">
                    <select
                      aria-label="Wildcard Element Selector"
                      value={queryEl}
                      onChange={(e) => {
                        const el = e.target.value as ElementSymbol;
                        setQueryEl(el);
                        selectElementTool(el);
                      }}
                      className="w-full h-5 bg-white border border-slate-400 text-[10px] text-slate-900"
                    >
                      {QUERY_ELEMENTS.map((el) => (
                        <option key={el} value={el}>
                          ? ({el})
                        </option>
                      ))}
                    </select>
                  </div>
                  <div className="h-7 border border-slate-400 bg-[#e8eff7]" />
                  <div className="h-7 border border-slate-400 bg-[#e8eff7]" />
                  <div className="h-7 border border-slate-400 bg-[#e8eff7]" />
                  <div className="h-7 border border-slate-400 bg-[#e8eff7]" />
                  <div className="h-7 border border-slate-400 bg-[#e8eff7]" />
                  <button
                    type="button"
                    onClick={() => selectElementTool('He')}
                    className={elemBtnClass('He')}
                  >
                    He
                  </button>
                </div>

                {/* Period 2: Li, Be, [3-col gap], B, C, N, O, F, Ne */}
                <div className="grid grid-cols-11">
                  <button
                    type="button"
                    onClick={() => selectElementTool('Li')}
                    className={elemBtnClass('Li')}
                  >
                    Li
                  </button>
                  <button
                    type="button"
                    onClick={() => selectElementTool('Be')}
                    className={elemBtnClass('Be')}
                  >
                    Be
                  </button>
                  <div className="col-span-3 h-7 border border-slate-400 bg-[#e8eff7]" />
                  {(['B', 'C', 'N', 'O', 'F', 'Ne'] as ElementSymbol[]).map((el) => (
                    <button
                      key={el}
                      type="button"
                      onClick={() => selectElementTool(el)}
                      className={elemBtnClass(el)}
                    >
                      {el}
                    </button>
                  ))}
                </div>

                {/* Period 3: Na, Mg, [3-col gap], Al, Si, P, S, Cl, Ar */}
                <div className="grid grid-cols-11">
                  <button
                    type="button"
                    onClick={() => selectElementTool('Na')}
                    className={elemBtnClass('Na')}
                  >
                    Na
                  </button>
                  <button
                    type="button"
                    onClick={() => selectElementTool('Mg')}
                    className={elemBtnClass('Mg')}
                  >
                    Mg
                  </button>
                  <div className="col-span-3 h-7 border border-slate-400 bg-[#e8eff7]" />
                  {(['Al', 'Si', 'P', 'S', 'Cl', 'Ar'] as ElementSymbol[]).map((el) => (
                    <button
                      key={el}
                      type="button"
                      onClick={() => selectElementTool(el)}
                      className={elemBtnClass(el)}
                    >
                      {el}
                    </button>
                  ))}
                </div>

                {/* Period 4: K, Ca, Sc, [Sc..Zn dropdown], Ga, Ge, As, Se, Br, Kr */}
                <div className="grid grid-cols-11">
                  <button
                    type="button"
                    onClick={() => selectElementTool('K')}
                    className={elemBtnClass('K')}
                  >
                    K
                  </button>
                  <button
                    type="button"
                    onClick={() => selectElementTool('Ca')}
                    className={elemBtnClass('Ca')}
                  >
                    Ca
                  </button>
                  <button
                    type="button"
                    onClick={() => selectElementTool(tmRow4El)}
                    className={elemBtnClass(tmRow4El)}
                  >
                    {tmRow4El}
                  </button>
                  <div className="col-span-2 h-7 border border-slate-400 bg-[#e8eff7] flex items-center px-0.5">
                    <select
                      aria-label="3d Transition Metal Selector"
                      value={tmRow4El}
                      onChange={(e) => {
                        const el = e.target.value as ElementSymbol;
                        setTmRow4El(el);
                        selectElementTool(el);
                      }}
                      className="w-full h-5 bg-white border border-slate-400 text-[10px] text-slate-900"
                    >
                      {ROW4_TM.map((el) => (
                        <option key={el} value={el}>
                          {el}
                        </option>
                      ))}
                    </select>
                  </div>
                  {(['Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'] as ElementSymbol[]).map((el) => (
                    <button
                      key={el}
                      type="button"
                      onClick={() => selectElementTool(el)}
                      className={elemBtnClass(el)}
                    >
                      {el}
                    </button>
                  ))}
                </div>

                {/* Period 5: Rb, Sr, Y, [Y..Cd dropdown], In, Sn, Sb, Te, I, Xe */}
                <div className="grid grid-cols-11">
                  <button
                    type="button"
                    onClick={() => selectElementTool('Rb')}
                    className={elemBtnClass('Rb')}
                  >
                    Rb
                  </button>
                  <button
                    type="button"
                    onClick={() => selectElementTool('Sr')}
                    className={elemBtnClass('Sr')}
                  >
                    Sr
                  </button>
                  <button
                    type="button"
                    onClick={() => selectElementTool(tmRow5El)}
                    className={elemBtnClass(tmRow5El)}
                  >
                    {tmRow5El}
                  </button>
                  <div className="col-span-2 h-7 border border-slate-400 bg-[#e8eff7] flex items-center px-0.5">
                    <select
                      aria-label="4d Transition Metal Selector"
                      value={tmRow5El}
                      onChange={(e) => {
                        const el = e.target.value as ElementSymbol;
                        setTmRow5El(el);
                        selectElementTool(el);
                      }}
                      className="w-full h-5 bg-white border border-slate-400 text-[10px] text-slate-900"
                    >
                      {ROW5_TM.map((el) => (
                        <option key={el} value={el}>
                          {el}
                        </option>
                      ))}
                    </select>
                  </div>
                  {(['In', 'Sn', 'Sb', 'Te', 'I', 'Xe'] as ElementSymbol[]).map((el) => (
                    <button
                      key={el}
                      type="button"
                      onClick={() => selectElementTool(el)}
                      className={elemBtnClass(el)}
                    >
                      {el}
                    </button>
                  ))}
                </div>

                {/* Period 6: Cs, Ba, Lu, [La..Lu/Hf..Hg dropdown], Tl, Pb, Bi, Po, At, Rn */}
                <div className="grid grid-cols-11">
                  <button
                    type="button"
                    onClick={() => selectElementTool('Cs')}
                    className={elemBtnClass('Cs')}
                  >
                    Cs
                  </button>
                  <button
                    type="button"
                    onClick={() => selectElementTool('Ba')}
                    className={elemBtnClass('Ba')}
                  >
                    Ba
                  </button>
                  <button
                    type="button"
                    onClick={() => selectElementTool(tmRow6El)}
                    className={elemBtnClass(tmRow6El)}
                  >
                    {tmRow6El}
                  </button>
                  <div className="col-span-2 h-7 border border-slate-400 bg-[#e8eff7] flex items-center px-0.5">
                    <select
                      aria-label="Lanthanide and 5d Metal Selector"
                      value={tmRow6El}
                      onChange={(e) => {
                        const el = e.target.value as ElementSymbol;
                        setTmRow6El(el);
                        selectElementTool(el);
                      }}
                      className="w-full h-5 bg-white border border-slate-400 text-[10px] text-slate-900"
                    >
                      {ROW6_LN_TM.map((el) => (
                        <option key={el} value={el}>
                          {el}
                        </option>
                      ))}
                    </select>
                  </div>
                  {(['Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn'] as ElementSymbol[]).map((el) => (
                    <button
                      key={el}
                      type="button"
                      onClick={() => selectElementTool(el)}
                      className={elemBtnClass(el)}
                    >
                      {el}
                    </button>
                  ))}
                </div>

                {/* Row 11: Export + Format Dropdown + Done */}
                <div className="grid grid-cols-[62px_1fr_62px] border-t border-slate-400">
                  <button
                    type="button"
                    onClick={handleExportClick}
                    className="h-7 text-[11px] font-medium border border-slate-400 bg-gradient-to-b from-[#e8f0f8] to-[#d0e0f0] hover:from-[#d8e8f8] hover:to-[#bfd4ec]"
                  >
                    Export
                  </button>
                  <div className="h-7 border border-slate-400 bg-[#e8eff7] flex items-center px-1">
                    <select
                      aria-label="Export Structure Format"
                      value={exportFormat}
                      onChange={(e) =>
                        setExportFormat(
                          e.target.value as 'MDL Molfile' | 'SMILES' | 'Formula'
                        )
                      }
                      className="w-full h-5 bg-white border border-slate-400 text-xs text-slate-900"
                    >
                      <option value="MDL Molfile">MDL Molfile</option>
                      <option value="SMILES">SMILES</option>
                      <option value="Formula">Formula</option>
                    </select>
                  </div>
                  <button
                    type="button"
                    onClick={() => onAccept(draftAtoms, draftBonds)}
                    className="h-7 text-[11px] font-medium border border-slate-400 bg-gradient-to-b from-[#e8f0f8] to-[#d0e0f0] hover:from-[#d8e8f8] hover:to-[#bfd4ec]"
                  >
                    Done
                  </button>
                </div>

                {/* Row 12: Hydrogen + Mode Dropdown + Help */}
                <div className="grid grid-cols-[62px_1fr_90px]">
                  <button
                    type="button"
                    onClick={() => handleHydrogenAction(hydrogenMode)}
                    className="h-7 text-[11px] font-medium border border-slate-400 bg-gradient-to-b from-[#e8f0f8] to-[#d0e0f0] hover:from-[#d8e8f8] hover:to-[#bfd4ec]"
                  >
                    Hydrogen
                  </button>
                  <div className="h-7 border border-slate-400 bg-[#e8eff7] flex items-center px-1">
                    <select
                      aria-label="Hydrogen Handling Mode"
                      value={hydrogenMode}
                      onChange={(e) =>
                        handleHydrogenAction(
                          e.target.value as 'Keep AsIs' | 'Add Implicit H' | 'Remove Explicit H'
                        )
                      }
                      className="w-full h-5 bg-white border border-slate-400 text-xs text-slate-900"
                    >
                      <option value="Keep AsIs">Keep AsIs</option>
                      <option value="Add Implicit H">Add Implicit H</option>
                      <option value="Remove Explicit H">Remove Explicit H</option>
                    </select>
                  </div>
                  <button
                    type="button"
                    onClick={() => setShowHelpBanner((v) => !v)}
                    className="h-7 text-xs font-bold underline border border-slate-400 bg-[#e8eff7] hover:bg-slate-200 text-slate-900"
                  >
                    Help
                  </button>
                </div>

                {/* Row 13: Import + Choose File */}
                <div className="grid grid-cols-[62px_1fr]">
                  <button
                    type="button"
                    onClick={() => fileInputRef.current?.click()}
                    className="h-7 text-[11px] font-medium border border-slate-400 bg-gradient-to-b from-[#e8f0f8] to-[#d0e0f0] hover:from-[#d8e8f8] hover:to-[#bfd4ec]"
                  >
                    Import
                  </button>
                  <div className="h-7 border border-slate-400 bg-[#e8eff7] flex items-center gap-1.5 px-1.5 overflow-hidden">
                    <input
                      ref={fileInputRef}
                      type="file"
                      accept=".mol,.sdf,.smi,.txt"
                      onChange={handleFileImport}
                      className="hidden"
                    />
                    <button
                      type="button"
                      onClick={() => fileInputRef.current?.click()}
                      className="px-2 py-0.5 text-[11px] bg-slate-100 hover:bg-slate-200 border border-slate-400 text-slate-900 shrink-0"
                    >
                      Choose File
                    </button>
                    <span className="text-[11px] text-slate-700 truncate">
                      {importedFileName}
                    </span>
                  </div>
                </div>
              </div>

              {/* Right Interactive 2D White Canvas */}
              <div className="bg-white relative flex flex-col justify-between h-[364px]">
                <svg
                  ref={svgRef}
                  viewBox={`0 0 ${SVG_W} ${SVG_H}`}
                  className="w-full h-full flex-1 cursor-crosshair select-none block"
                  onClick={handleCanvasClick}
                  onMouseDown={(e) => {
                    if (
                      activeTool === 'move' &&
                      (e.target === svgRef.current || (e.target as Element).tagName === 'rect')
                    ) {
                      const { x, y } = getSvgCoords(e);
                      setDraggingCanvas({ startX: x, startY: y });
                    }
                  }}
                  onMouseMove={handleMouseMove}
                  onMouseUp={() => {
                    setDraggingAtomId(null);
                    setDraggingCanvas(null);
                  }}
                  data-testid="pubchem-modal-svg"
                >
                  <rect width={SVG_W} height={SVG_H} fill="#FFFFFF" />

                  {/* Render Bonds */}
                  {draftBonds.map((bond) => {
                    const s = draftAtoms.find((a) => a.id === bond.source);
                    const t = draftAtoms.find((a) => a.id === bond.target);
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

                    // Wavy path generator
                    const buildWavyPath = () => {
                      const steps = 8;
                      let d = `M ${x1} ${y1}`;
                      for (let k = 1; k <= steps; k++) {
                        const f = k / steps;
                        const mx = x1 + dx * f + (k % 2 === 0 ? 1 : -1) * nx * 3.5;
                        const my = y1 + dy * f + (k % 2 === 0 ? 1 : -1) * ny * 3.5;
                        d += ` L ${mx.toFixed(1)} ${my.toFixed(1)}`;
                      }
                      return d;
                    };

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
                            strokeWidth={isHovered ? 3.2 : 2.2}
                            strokeLinecap="round"
                          />
                        )}

                        {bond.order === 1 && bond.stereo === 'wedge' && (
                          <polygon
                            points={`${x1},${y1} ${x2 + nx * 5.5},${y2 + ny * 5.5} ${x2 - nx * 5.5},${y2 - ny * 5.5}`}
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
                            strokeWidth={3.2}
                            strokeDasharray="3 4"
                          />
                        )}

                        {bond.order === 1 && bond.stereo === 'wavy' && (
                          <path
                            d={buildWavyPath()}
                            fill="none"
                            stroke={strokeColor}
                            strokeWidth={2.0}
                          />
                        )}

                        {bond.order === 2 && bond.stereo === 'crossed' && (
                          <>
                            <line
                              x1={x1 + nx * 3.5}
                              y1={y1 + ny * 3.5}
                              x2={x2 - nx * 3.5}
                              y2={y2 - ny * 3.5}
                              stroke={strokeColor}
                              strokeWidth={2.0}
                            />
                            <line
                              x1={x1 - nx * 3.5}
                              y1={y1 - ny * 3.5}
                              x2={x2 + nx * 3.5}
                              y2={y2 + ny * 3.5}
                              stroke={strokeColor}
                              strokeWidth={2.0}
                            />
                          </>
                        )}

                        {bond.order === 2 && bond.stereo !== 'crossed' && (
                          <>
                            <line
                              x1={x1 + nx * 3.2}
                              y1={y1 + ny * 3.2}
                              x2={x2 + nx * 3.2}
                              y2={y2 + ny * 3.2}
                              stroke={strokeColor}
                              strokeWidth={2.1}
                              strokeLinecap="round"
                            />
                            <line
                              x1={x1 - nx * 3.2}
                              y1={y1 - ny * 3.2}
                              x2={x2 - nx * 3.2}
                              y2={y2 - ny * 3.2}
                              stroke={strokeColor}
                              strokeWidth={2.1}
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

                  {/* Interactive Preview Line when an Atom is Selected in Bond Mode */}
                  {selectedAtomId &&
                    getBondSpecFromTool(activeTool) !== null &&
                    (() => {
                      const srcAtom = draftAtoms.find((a) => a.id === selectedAtomId);
                      if (!srcAtom) return null;
                      const targetAtom =
                        hoveredAtomId && hoveredAtomId !== selectedAtomId
                          ? draftAtoms.find((a) => a.id === hoveredAtomId)
                          : null;
                      const tx = targetAtom
                        ? toScreenX(targetAtom.x2d)
                        : cursorChemPos
                          ? toScreenX(cursorChemPos.x)
                          : null;
                      const ty = targetAtom
                        ? toScreenY(targetAtom.y2d)
                        : cursorChemPos
                          ? toScreenY(cursorChemPos.y)
                          : null;
                      if (tx === null || ty === null) return null;
                      const sx = toScreenX(srcAtom.x2d);
                      const sy = toScreenY(srcAtom.y2d);
                      return (
                        <line
                          x1={sx}
                          y1={sy}
                          x2={tx}
                          y2={ty}
                          stroke={targetAtom ? '#2563EB' : '#94A3B8'}
                          strokeWidth={targetAtom ? 2.4 : 1.6}
                          strokeDasharray="4 3"
                          className="pointer-events-none"
                        />
                      );
                    })()}

                  {/* Render Atoms */}
                  {draftAtoms.map((atom) => {
                    const cx = toScreenX(atom.x2d);
                    const cy = toScreenY(atom.y2d);
                    const isHovered = hoveredAtomId === atom.id;
                    const isSelected = selectedAtomId === atom.id;
                    const color = getReadableAtomColor(atom.element);
                    const chg = atom.formalCharge ?? 0;
                    const chgLabel =
                      chg > 0 ? (chg === 1 ? '+' : `${chg}+`) : chg < 0 ? (chg === -1 ? '−' : `${Math.abs(chg)}−`) : '';

                    return (
                      <g
                        key={atom.id}
                        transform={`translate(${cx}, ${cy})`}
                        className="cursor-pointer"
                        onMouseEnter={() => setHoveredAtomId(atom.id)}
                        onMouseLeave={() => setHoveredAtomId(null)}
                        onMouseDown={(e) => {
                          if (activeTool === 'move') {
                            e.stopPropagation();
                            setDraggingAtomId(atom.id);
                          }
                        }}
                        onClick={(e) => handleAtomClick(e, atom)}
                      >
                        {isSelected && (
                          <circle
                            r={18}
                            fill="none"
                            stroke="#2563EB"
                            strokeWidth="1.8"
                            strokeDasharray="3 2"
                          />
                        )}
                        <circle
                          r={14}
                          fill={isSelected ? '#EFF6FF' : '#FFFFFF'}
                          stroke={
                            isSelected
                              ? '#2563EB'
                              : isHovered
                                ? activeTool === 'erase'
                                  ? '#DC2626'
                                  : '#2563EB'
                                : '#CBD5E1'
                          }
                          strokeWidth={isSelected || isHovered ? 2.2 : 1.2}
                        />
                        <text
                          textAnchor="middle"
                          dominantBaseline="central"
                          fill={color}
                          className="font-mono text-xs font-bold select-none"
                        >
                          {atom.element}
                        </text>
                        {chgLabel && (
                          <text
                            x={10}
                            y={-8}
                            textAnchor="middle"
                            fill={chg > 0 ? '#DC2626' : '#2563EB'}
                            className="font-mono text-[10px] font-bold select-none"
                          >
                            {chgLabel}
                          </text>
                        )}
                      </g>
                    );
                  })}
                </svg>

                {/* Optional Export Textarea Drawer */}
                {showExportPanel && (
                  <div className="p-2.5 bg-slate-50 border-t border-slate-300 text-xs">
                    <div className="flex items-center justify-between mb-1">
                      <span className="font-semibold text-slate-700">
                        Exported {exportFormat}:
                      </span>
                      <button
                        type="button"
                        onClick={() => setShowExportPanel(false)}
                        className="text-slate-500 hover:text-slate-800 text-[11px]"
                      >
                        Close
                      </button>
                    </div>
                    <textarea
                      readOnly
                      value={exportText}
                      rows={4}
                      className="w-full p-1.5 font-mono text-[11px] bg-white border border-slate-300 rounded-xs text-slate-800"
                    />
                  </div>
                )}

                {/* Optional Quick Help Strip */}
                {showHelpBanner && (
                  <div className="px-3 py-2 bg-amber-50 border-t border-amber-200 text-[11px] text-slate-700 flex items-center justify-between">
                    <span>
                      In Single/Double/Triple bond mode, click an atom first then click a 2nd atom to bond them, click empty canvas space to sprout a bonded {activeElement} at that spot, or click the same atom twice to replace it with {activeElement}. Clicking empty space with no atom selected places a standalone {activeElement}.
                    </span>
                    <button
                      type="button"
                      onClick={() => setShowHelpBanner(false)}
                      className="ml-2 font-semibold text-slate-900 underline shrink-0"
                    >
                      Dismiss
                    </button>
                  </div>
                )}
              </div>
            </div>
          </div>
        </div>

        {/* Bottom Modal Footer: Full-Width Accept Structure & Cancel Buttons */}
        <div className="px-4 pb-4 pt-1 grid grid-cols-2 gap-4">
          <button
            type="button"
            onClick={() => onAccept(draftAtoms, draftBonds)}
            className="w-full py-2.5 bg-[#0071bc] hover:bg-[#005a96] text-white text-sm font-medium rounded-xs shadow-xs transition-colors cursor-pointer"
          >
            Accept Structure
          </button>
          <button
            type="button"
            onClick={onCancel}
            className="w-full py-2.5 bg-[#1b2b38] hover:bg-[#111c24] text-white text-sm font-medium rounded-xs shadow-xs transition-colors cursor-pointer"
          >
            Cancel
          </button>
        </div>
      </div>
    </div>
  );
};
