/**
 * @license
 * SPDX-License-Identifier: Apache-2.0
 */

import React, { useEffect, useMemo, useRef, useState } from 'react';
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { RoomEnvironment } from 'three/addons/environments/RoomEnvironment.js';
import {
  AtomNode,
  ELEMENT_PROPERTIES,
  MolecularOrbital,
  MoleculeEntry,
} from '../chemistry/moleculesData';
import {
  detectHapticRingCoordinations,
  evaluateOrbitalWavefunction,
} from '../chemistry/moEngine';
import { buildMorphIsosurfaceKeyframe } from '../chemistry/moIsosurface';
import { PointGroupData, SymmetryElement3D } from '../chemistry/symmetryData';
import { LowdinEHTResult } from '../chemistry/quantumEngineV3';
import { GaussianCubeData } from '../chemistry/fileIO';
import {
  ColoredIsosurfaceMeshData,
  Slice2DResult,
  SliceOrientation,
  VolumetricFieldMode,
  buildVolumetricSurfaceV3,
  compute2DSlicePlane,
} from '../chemistry/volumetricV3';
import {
  formulaToLatex,
  irrepToLatex,
  moLabelToLatex,
  schoenfliesToLatex,
} from '../chemistry/latexUtils';
import { Latex } from './Latex';
import {
  ChevronDown,
  ChevronUp,
  Grid2x2,
  Pause,
  Play,
  Sun,
} from 'lucide-react';

interface OrbitalVisualizer3DV3Props {
  molecule: MoleculeEntry;
  displayAtoms?: AtomNode[];
  selectedMo: MolecularOrbital;
  pointGroup: PointGroupData;
  alignedSymmetryElements?: SymmetryElement3D[];
  eht: LowdinEHTResult;
  fieldMode: VolumetricFieldMode;
  onChangeFieldMode: (mode: VolumetricFieldMode) => void;
  importedCube?: GaussianCubeData | null;
  highlightedSymmetryElementId?: string | null;
  walshBanner?: string | null;
  onInspectOpenShellO2?: () => void;
  onIonizeRadical?: () => void;
}

const MORPH_KEYFRAME_STEPS = 40;

function createChargeSprite(text: string, isNeg: boolean, isPos: boolean): THREE.Sprite {
  const canvas = document.createElement('canvas');
  canvas.width = 128;
  canvas.height = 48;
  const ctx = canvas.getContext('2d');
  if (ctx) {
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.fillStyle = isNeg
      ? 'rgba(190, 18, 60, 0.90)'
      : isPos
        ? 'rgba(29, 78, 216, 0.90)'
        : 'rgba(51, 65, 85, 0.88)';
    ctx.beginPath();
    ctx.roundRect(8, 6, 112, 36, 8);
    ctx.fill();
    ctx.strokeStyle = 'rgba(255,255,255,0.45)';
    ctx.lineWidth = 2;
    ctx.stroke();
    ctx.fillStyle = '#FFFFFF';
    ctx.font = 'bold 22px "IBM Plex Mono", monospace';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.fillText(text, 64, 24);
  }
  const tex = new THREE.CanvasTexture(canvas);
  tex.colorSpace = THREE.SRGBColorSpace;
  tex.needsUpdate = true;
  const mat = new THREE.SpriteMaterial({ map: tex, transparent: true, depthTest: false });
  const sprite = new THREE.Sprite(mat);
  sprite.scale.set(0.88, 0.33, 1);
  sprite.renderOrder = 10;
  return sprite;
}

export const OrbitalVisualizer3DV3: React.FC<OrbitalVisualizer3DV3Props> = ({
  molecule,
  displayAtoms,
  selectedMo,
  pointGroup,
  alignedSymmetryElements,
  eht,
  fieldMode,
  onChangeFieldMode,
  importedCube = null,
  highlightedSymmetryElementId = null,
  walshBanner = null,
  onInspectOpenShellO2,
  onIonizeRadical,
}) => {
  const atoms = displayAtoms || molecule.atoms;

  const canvasMountRef = useRef<HTMLDivElement | null>(null);
  const sliceInsetCanvasRef = useRef<HTMLCanvasElement | null>(null);
  const requestRenderRef = useRef<(() => void) | null>(null);
  const resetCameraRef = useRef<((preset?: 'iso' | 'top' | 'front' | 'side') => void) | null>(null);
  const rebuildSceneRef = useRef<(() => void) | null>(null);
  const updateMorphFastRef = useRef<((tProgress: number) => void) | null>(null);

  const [lobeScale, setLobeScale] = useState<number>(1.0);
  const [orbitalOpacity, setOrbitalOpacity] = useState<number>(0.72);
  const [renderMode, setRenderMode] = useState<'merged' | 'wireframe'>('merged');
  const [morphProgress, setMorphProgress] = useState<number>(1.0);
  const [isMorphPlaying, setIsMorphPlaying] = useState<boolean>(false);
  const [showSymmetryElements, setShowSymmetryElements] = useState<boolean>(true);
  const [showNodalPlanes, setShowNodalPlanes] = useState<boolean>(false);
  const [showChargeLabels, setShowChargeLabels] = useState<boolean>(false);
  const [showDipoleArrow, setShowDipoleArrow] = useState<boolean>(true);
  const [darkStage, setDarkStage] = useState<boolean>(true);
  const [colorblindPalette, setColorblindPalette] = useState<boolean>(false);
  const [dockCollapsed, setDockCollapsed] = useState<boolean>(true);

  // 2D Contour & Phase Slice Plane State
  const [showSlicePlane, setShowSlicePlane] = useState<boolean>(false);
  const [sliceOrientation, setSliceOrientation] = useState<SliceOrientation>('xz');
  const [sliceOffset, setSliceOffset] = useState<number>(0);
  const [sliceIsolines, setSliceIsolines] = useState<boolean>(true);
  const [sliceProbe, setSliceProbe] = useState<{ u: number; v: number; val: number } | null>(null);
  const [spinEmptyState, setSpinEmptyState] = useState<boolean>(false);
  const [mepRange, setMepRange] = useState<[number, number] | null>(null);

  const morphStartLinearRef = useRef<number>(0);

  // Compute 2D slice plane when enabled
  const sliceData: Slice2DResult | null = useMemo(() => {
    if (!showSlicePlane) return null;
    return compute2DSlicePlane({
      fieldMode,
      atoms,
      selectedMo,
      eht,
      importedCube,
      orientation: sliceOrientation,
      offset: sliceOffset,
      showIsolines: sliceIsolines,
      resolution: 72,
    });
  }, [
    showSlicePlane,
    fieldMode,
    atoms,
    selectedMo,
    eht,
    importedCube,
    sliceOrientation,
    sliceOffset,
    sliceIsolines,
  ]);

  // Draw 2D slice into inset preview canvas
  useEffect(() => {
    const cvs = sliceInsetCanvasRef.current;
    if (!cvs || !sliceData) return;
    const res = sliceData.width;
    cvs.width = res;
    cvs.height = sliceData.height;
    const ctx = cvs.getContext('2d');
    if (!ctx) return;
    const imgData = ctx.createImageData(res, sliceData.height);
    // Flip vertically so +v is at top of 2D inset
    for (let y = 0; y < sliceData.height; y++) {
      const srcRow = (sliceData.height - 1 - y) * res * 4;
      const dstRow = y * res * 4;
      for (let x = 0; x < res * 4; x++) {
        imgData.data[dstRow + x] = sliceData.rgba[srcRow + x];
      }
    }
    ctx.putImageData(imgData, 0, 0);
  }, [sliceData]);

  const stateRef = useRef({
    molecule,
    atoms,
    selectedMo,
    pointGroup,
    alignedSymmetryElements,
    eht,
    fieldMode,
    importedCube,
    highlightedSymmetryElementId,
    lobeScale,
    orbitalOpacity,
    renderMode,
    morphProgress,
    showSymmetryElements,
    showNodalPlanes,
    showChargeLabels,
    showDipoleArrow,
    darkStage,
    colorblindPalette,
    showSlicePlane,
    sliceData,
  });

  stateRef.current = {
    molecule,
    atoms,
    selectedMo,
    pointGroup,
    alignedSymmetryElements,
    eht,
    fieldMode,
    importedCube,
    highlightedSymmetryElementId,
    lobeScale,
    orbitalOpacity,
    renderMode,
    morphProgress,
    showSymmetryElements,
    showNodalPlanes,
    showChargeLabels,
    showDipoleArrow,
    darkStage,
    colorblindPalette,
    showSlicePlane,
    sliceData,
  };

  useEffect(() => {
    setIsMorphPlaying(false);
  }, [molecule.id, fieldMode]);

  // Smooth AO -> MO morph playback (only in MO field mode)
  useEffect(() => {
    if (!isMorphPlaying || fieldMode !== 'mo') return;

    let rafId = 0;
    const startWallTime = performance.now();
    const startU = morphStartLinearRef.current;
    const totalDurationMs = 2200;
    const remainingDurationMs = Math.max(200, (1 - startU) * totalDurationMs);

    const stepAnimation = (now: number) => {
      const elapsed = Math.max(0, now - startWallTime);
      const frac = Math.min(1, elapsed / remainingDurationMs);
      const u = startU + frac * (1 - startU);
      const easedT = 0.5 * (1 - Math.cos(Math.PI * u));

      stateRef.current.morphProgress = easedT;
      updateMorphFastRef.current?.(easedT);
      setMorphProgress(easedT);

      if (frac >= 1.0) {
        setIsMorphPlaying(false);
        return;
      }
      rafId = requestAnimationFrame(stepAnimation);
    };

    rafId = requestAnimationFrame(stepAnimation);
    return () => cancelAnimationFrame(rafId);
  }, [isMorphPlaying, fieldMode]);

  const handleTogglePlayMorph = () => {
    if (fieldMode !== 'mo') onChangeFieldMode('mo');
    if (!isMorphPlaying) {
      const startVal = morphProgress >= 0.96 ? 0 : morphProgress;
      morphStartLinearRef.current =
        startVal <= 0
          ? 0
          : startVal >= 1
            ? 0
            : Math.acos(Math.max(-1, Math.min(1, 1 - 2 * startVal))) / Math.PI;
      if (startVal === 0) {
        stateRef.current.morphProgress = 0;
        updateMorphFastRef.current?.(0);
        setMorphProgress(0);
      }
      setIsMorphPlaying(true);
    } else {
      setIsMorphPlaying(false);
    }
  };

  const handleScrubMorph = (nextValue: number) => {
    setIsMorphPlaying(false);
    stateRef.current.morphProgress = nextValue;
    updateMorphFastRef.current?.(nextValue);
    setMorphProgress(nextValue);
  };

  // Persistent WebGL Renderer Lifecycle
  useEffect(() => {
    const container = canvasMountRef.current;
    if (!container) return;

    const isSoftwareGl = (() => {
      const gl = document.createElement('canvas').getContext('webgl2');
      if (!gl) return true;
      const ext = gl.getExtension('WEBGL_debug_renderer_info');
      const param = ext ? ext.UNMASKED_RENDERER_WEBGL : gl.RENDERER;
      const gpu = String(gl.getParameter(param));
      gl.getExtension('WEBGL_lose_context')?.loseContext();
      return /SwiftShader|llvmpipe|softpipe|Software|Basic Render/i.test(gpu);
    })();

    const width = container.clientWidth || 560;
    const height = container.clientHeight || 360;

    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0x0b0f19);
    scene.fog = new THREE.FogExp2(0x0b0f19, 0.032);

    const camera = new THREE.PerspectiveCamera(42, width / height, 0.1, 250);
    camera.position.set(0, 2.2, 6.2);

    const renderer = new THREE.WebGLRenderer({
      antialias: !isSoftwareGl,
      powerPreference: 'high-performance',
      preserveDrawingBuffer: true,
    });
    renderer.setSize(width, height);
    renderer.outputColorSpace = THREE.SRGBColorSpace;
    renderer.toneMapping = THREE.ACESFilmicToneMapping;
    renderer.toneMappingExposure = 1.12;
    renderer.setPixelRatio(isSoftwareGl ? 1 : Math.min(window.devicePixelRatio, 2));
    renderer.shadowMap.enabled = !isSoftwareGl;
    renderer.shadowMap.type = THREE.PCFSoftShadowMap;

    container.innerHTML = '';
    container.appendChild(renderer.domElement);

    const pmremGenerator = new THREE.PMREMGenerator(renderer);
    scene.environment = pmremGenerator.fromScene(new RoomEnvironment(), 0.04).texture;
    pmremGenerator.dispose();

    const hemiLight = new THREE.HemisphereLight(0xdbeafe, 0x1e293b, 0.85);
    scene.add(hemiLight);

    const keyLight = new THREE.DirectionalLight(0xfffbeb, 1.45);
    keyLight.position.set(5, 7, 6);
    if (!isSoftwareGl) {
      keyLight.castShadow = true;
      keyLight.shadow.mapSize.set(1024, 1024);
      keyLight.shadow.bias = -0.0005;
    }
    scene.add(keyLight);

    const rimLight = new THREE.DirectionalLight(0x38bdf8, 0.65);
    rimLight.position.set(-5, -3, -5);
    scene.add(rimLight);

    const groundGroup = new THREE.Group();
    const shadowGeo = new THREE.CircleGeometry(3.4, 32);
    const shadowMat = new THREE.MeshBasicMaterial({
      color: 0x020617,
      transparent: true,
      opacity: 0.35,
      depthWrite: false,
    });
    const shadowDisc = new THREE.Mesh(shadowGeo, shadowMat);
    shadowDisc.rotation.x = -Math.PI / 2;
    shadowDisc.position.y = -2.25;
    groundGroup.add(shadowDisc);

    const ringGeo = new THREE.RingGeometry(2.1, 2.14, 48);
    const ringMat = new THREE.MeshBasicMaterial({
      color: 0x334155,
      side: THREE.DoubleSide,
      transparent: true,
      opacity: 0.45,
      depthWrite: false,
    });
    const floorRing = new THREE.Mesh(ringGeo, ringMat);
    floorRing.rotation.x = -Math.PI / 2;
    floorRing.position.y = -2.23;
    groundGroup.add(floorRing);
    scene.add(groundGroup);

    const controls = new OrbitControls(camera, renderer.domElement);
    controls.enableDamping = !isSoftwareGl;
    controls.dampingFactor = 0.06;
    controls.minDistance = 1.6;
    controls.maxDistance = 42.0;

    const contentGroup = new THREE.Group();
    scene.add(contentGroup);

    const morphKeyframeCache = new Map<
      string,
      { posGeo: THREE.BufferGeometry; negGeo: THREE.BufferGeometry }
    >();
    let activePosMesh: THREE.Mesh | null = null;
    let activeNegMesh: THREE.Mesh | null = null;
    let activeMorphSignature = '';
    let lastFramedMolId = '';

    const clearMorphCache = () => {
      for (const entry of morphKeyframeCache.values()) {
        entry.posGeo.dispose();
        entry.negGeo.dispose();
      }
      morphKeyframeCache.clear();
    };

    const buildBufferGeoFromMeshData = (
      data: ColoredIsosurfaceMeshData
    ): THREE.BufferGeometry => {
      const geo = new THREE.BufferGeometry();
      if (data.positions.length > 0 && data.indices.length > 0) {
        geo.setAttribute('position', new THREE.BufferAttribute(data.positions, 3));
        geo.setAttribute('normal', new THREE.BufferAttribute(data.normals, 3));
        if (data.colors && data.colors.length === data.positions.length) {
          geo.setAttribute('color', new THREE.BufferAttribute(data.colors, 3));
        }
        geo.setIndex(new THREE.BufferAttribute(data.indices, 1));
      }
      return geo;
    };

    const getMorphKeyframe = (
      mo: MolecularOrbital,
      atomList: AtomNode[],
      scaleMult: number,
      tProgress: number
    ) => {
      const clampedT = Math.max(0, Math.min(1, tProgress));
      const stepIdx = Math.round(clampedT * MORPH_KEYFRAME_STEPS);
      const atomSig = atomList
        .slice(0, 4)
        .map((a) => `${a.x3d.toFixed(2)},${a.y3d.toFixed(2)},${a.z3d.toFixed(2)}`)
        .join('|');
      const sig = `${mo.id}|${atomList.length}|${atomSig}|${scaleMult.toFixed(2)}`;
      if (sig !== activeMorphSignature) {
        clearMorphCache();
        activeMorphSignature = sig;
      }
      const key = `${sig}|${stepIdx}`;
      let entry = morphKeyframeCache.get(key);
      if (!entry) {
        const res = buildMorphIsosurfaceKeyframe(
          mo,
          atomList,
          stepIdx / MORPH_KEYFRAME_STEPS,
          scaleMult
        );
        entry = {
          posGeo: buildBufferGeoFromMeshData(res.positive),
          negGeo: buildBufferGeoFromMeshData(res.negative),
        };
        morphKeyframeCache.set(key, entry);
      }
      return entry;
    };

    let renderScheduled = false;
    let isDisposed = false;
    const maxAnisotropy = renderer.capabilities.getMaxAnisotropy();

    const doRender = () => {
      renderScheduled = false;
      if (isDisposed) return;
      controls.update();
      renderer.render(scene, camera);
    };

    const requestRender = () => {
      if (isDisposed || renderScheduled) return;
      renderScheduled = true;
      requestAnimationFrame(doRender);
    };
    requestRenderRef.current = requestRender;

    const updateMorphFast = (tProgress: number) => {
      if (!activePosMesh || !activeNegMesh) return;
      const { selectedMo: mo, atoms: atomList, lobeScale: scaleMult, fieldMode: fm } = stateRef.current;
      if (!mo || fm !== 'mo') return;
      const kf = getMorphKeyframe(mo, atomList, scaleMult, tProgress);
      activePosMesh.geometry = kf.posGeo;
      activeNegMesh.geometry = kf.negGeo;
      requestRender();
    };
    updateMorphFastRef.current = updateMorphFast;

    const disposeContentGroup = () => {
      const cachedGeos = new Set<THREE.BufferGeometry>();
      for (const e of morphKeyframeCache.values()) {
        cachedGeos.add(e.posGeo);
        cachedGeos.add(e.negGeo);
      }
      while (contentGroup.children.length > 0) {
        const child = contentGroup.children[0];
        contentGroup.remove(child);
        child.traverse((obj: any) => {
          if (obj.geometry && !cachedGeos.has(obj.geometry)) {
            obj.geometry.dispose();
          }
          if (obj.material) {
            const mats = Array.isArray(obj.material) ? obj.material : [obj.material];
            for (const m of mats) {
              if (m.map) m.map.dispose();
              m.dispose();
            }
          }
        });
      }
      activePosMesh = null;
      activeNegMesh = null;
    };

    const buildMolecularAndOrbitalScene = () => {
      disposeContentGroup();

      const {
        molecule: mol,
        atoms: atomList,
        selectedMo: mo,
        pointGroup: pg,
        alignedSymmetryElements: alignedSym,
        eht: ehtData,
        fieldMode: fm,
        importedCube: cubeData,
        highlightedSymmetryElementId: hiSymId,
        lobeScale: scaleMult,
        orbitalOpacity: opacityVal,
        renderMode: mode,
        morphProgress: tProgress,
        showSymmetryElements: showSym,
        showNodalPlanes: showNodes,
        showChargeLabels: showCharges,
        showDipoleArrow: showDipole,
        darkStage: isDark,
        colorblindPalette: cbMode,
        showSlicePlane: sliceOn,
        sliceData: sliceRes,
      } = stateRef.current;

      // Compute molecular center & bounding radius
      let cx = 0,
        cy = 0,
        cz = 0;
      for (const a of atomList) {
        cx += a.x3d;
        cy += a.y3d;
        cz += a.z3d;
      }
      const nA = Math.max(1, atomList.length);
      cx /= nA;
      cy /= nA;
      cz /= nA;

      let maxR = 1.6;
      for (const a of atomList) {
        const r = Math.hypot(a.x3d - cx, a.y3d - cy, a.z3d - cz);
        if (r > maxR) maxR = r;
      }
      const rMol = Math.max(1.6, maxR + 0.6);

      // Update stage theme
      const bgColorHex = isDark ? 0x0b0f19 : 0xf1f5f9;
      scene.background = new THREE.Color(bgColorHex);
      scene.fog = new THREE.FogExp2(bgColorHex, Math.max(0.008, 0.035 / Math.max(1, rMol / 2.6)));
      shadowMat.opacity = isDark ? 0.35 : 0.08;
      ringMat.color.setHex(isDark ? 0x334155 : 0x94a3b8);
      groundGroup.position.set(cx, cy - rMol - 0.85, cz);
      const ringScale = Math.max(1, rMol / 2.1);
      floorRing.scale.setScalar(ringScale);
      shadowDisc.scale.setScalar(ringScale);

      const hapticRings = detectHapticRingCoordinations(atomList, mol.bonds);
      const isZAlignedSandwich =
        hapticRings.length > 0 &&
        hapticRings.every((hg) => {
          const metal = atomList.find((a) => a.id === hg.metalId);
          const ringCs = hg.ringCarbonIds
            .map((id) => atomList.find((a) => a.id === id))
            .filter((a): a is AtomNode => Boolean(a));
          if (!metal || ringCs.length < 4) return false;
          const rcx = ringCs.reduce((s, a) => s + a.x3d, 0) / ringCs.length;
          const rcy = ringCs.reduce((s, a) => s + a.y3d, 0) / ringCs.length;
          const rcz = ringCs.reduce((s, a) => s + a.z3d, 0) / ringCs.length;
          const dz = Math.abs(rcz - metal.z3d);
          const dAll = Math.hypot(rcx - metal.x3d, rcy - metal.y3d, dz) || 1;
          return dz / dAll > 0.85;
        });

      // Re-frame camera cleanly when switching molecule (include formula + atom count so consecutive PubChem loads re-frame)
      const frameKey = `${mol.id}|${mol.formula}|${atomList.length}`;
      if (lastFramedMolId !== frameKey) {
        lastFramedMolId = frameKey;
        const dist = Math.max(6.2, rMol * 2.95);
        if (isZAlignedSandwich) {
          camera.position.set(cx + dist * 0.82, cy + dist * 0.35, cz + dist * 0.45);
        } else {
          camera.position.set(cx + dist * 0.15, cy + dist * 0.34, cz + dist * 0.94);
        }
        controls.target.set(cx, cy, cz);
        controls.update();
      }

      const vUp = new THREE.Vector3(0, 1, 0);

      // 1. Render Atoms
      const atomMap = new Map<string, AtomNode>();
      for (const a of atomList) {
        atomMap.set(a.id, a);
        const props = ELEMENT_PROPERTIES[a.element] || ELEMENT_PROPERTIES.C;
        const rad = Math.max(0.22, Math.min(0.56, props.covalentRadius * 0.42));
        const sphereGeo = new THREE.SphereGeometry(rad, 24, 18);
        const sphereMat = new THREE.MeshStandardMaterial({
          color: new THREE.Color(props.cpkColor),
          roughness: 0.28,
          metalness: props.atomicNumber >= 21 && props.atomicNumber <= 30 ? 0.45 : 0.12,
        });
        const mesh = new THREE.Mesh(sphereGeo, sphereMat);
        mesh.position.set(a.x3d, a.y3d, a.z3d);
        mesh.castShadow = !isSoftwareGl;
        mesh.receiveShadow = !isSoftwareGl;
        contentGroup.add(mesh);
      }

      // 2. Render Bonds (with textbook haptic centroid axes + inner aromatic tori for η⁴/η⁵/η⁶ rings)
      const hapticMetalBondIds = new Set<string>();
      const hapticPerimeterBondIds = new Set<string>();
      for (const hg of hapticRings) {
        hg.metalCarbonBondIds.forEach((id) => hapticMetalBondIds.add(id));
        hg.ringPerimeterBondIds.forEach((id) => hapticPerimeterBondIds.add(id));
      }

      const bondMat = new THREE.MeshStandardMaterial({
        color: isDark ? 0x94a3b8 : 0x64748b,
        roughness: 0.42,
        metalness: 0.18,
      });

      for (const b of mol.bonds) {
        if (hapticMetalBondIds.has(b.id)) continue;
        const a1 = atomMap.get(b.source);
        const a2 = atomMap.get(b.target);
        if (!a1 || !a2) continue;
        const p1 = new THREE.Vector3(a1.x3d, a1.y3d, a1.z3d);
        const p2 = new THREE.Vector3(a2.x3d, a2.y3d, a2.z3d);
        const diff = p2.clone().sub(p1);
        const len = diff.length();
        if (len < 0.05) continue;
        const mid = p1.clone().add(p2).multiplyScalar(0.5);
        const dir = diff.clone().normalize();
        const perp = new THREE.Vector3()
          .crossVectors(dir, Math.abs(dir.y) < 0.9 ? vUp : new THREE.Vector3(1, 0, 0))
          .normalize();

        const effectiveOrder = hapticPerimeterBondIds.has(b.id) ? 1 : b.order;
        const offsets =
          effectiveOrder >= 2.5
            ? [-0.11, 0, 0.11]
            : effectiveOrder >= 1.5
              ? [-0.075, 0.075]
              : [0];
        const cylRad = effectiveOrder >= 1.5 ? 0.055 : 0.08;
        const cylGeo = new THREE.CylinderGeometry(cylRad, cylRad, len, 12);

        for (const off of offsets) {
          const stick = new THREE.Mesh(cylGeo, bondMat);
          stick.position.copy(mid).addScaledVector(perp, off);
          stick.quaternion.setFromUnitVectors(vUp, dir);
          contentGroup.add(stick);
        }
      }

      // Render crisp segmented metal-to-centroid haptic coordination axes & inscribed aromatic ring tori
      if (hapticRings.length > 0) {
        const hapticMat = new THREE.MeshStandardMaterial({
          color: isDark ? 0xcbd5e1 : 0x475569,
          roughness: 0.35,
          metalness: 0.32,
        });
        for (const hg of hapticRings) {
          const metal = atomMap.get(hg.metalId);
          const rPts = hg.ringCarbonIds
            .map((id) => atomMap.get(id))
            .filter((a): a is AtomNode => Boolean(a))
            .map((a) => new THREE.Vector3(a.x3d, a.y3d, a.z3d));
          if (!metal || rPts.length < 4) continue;

          const pMetal = new THREE.Vector3(metal.x3d, metal.y3d, metal.z3d);
          const pCent = new THREE.Vector3();
          for (const p of rPts) pCent.add(p);
          pCent.divideScalar(rPts.length);

          // Segmented dashed coordination axis from metal to ring centroid
          const axisVec = pCent.clone().sub(pMetal);
          const axisLen = axisVec.length();
          if (axisLen > 0.25) {
            const axisDir = axisVec.clone().normalize();
            const segCount = 5;
            const segLen = (axisLen / segCount) * 0.54;
            const dashGeo = new THREE.CylinderGeometry(0.044, 0.044, segLen, 10);
            for (let s = 0; s < segCount; s++) {
              const frac = (s + 0.5) / segCount;
              const segMesh = new THREE.Mesh(dashGeo, hapticMat);
              segMesh.position.copy(pMetal).addScaledVector(axisVec, frac);
              segMesh.quaternion.setFromUnitVectors(vUp, axisDir);
              contentGroup.add(segMesh);
            }
          }

          // Inscribed aromatic ring torus lying in the plane of the ring carbons
          const avgR =
            rPts.reduce((sum, p) => sum + p.distanceTo(pCent), 0) / rPts.length;
          if (avgR > 0.3) {
            const e1 = rPts[0].clone().sub(pCent);
            const e2 = rPts[1].clone().sub(pCent);
            const ringNorm = new THREE.Vector3().crossVectors(e1, e2);
            if (ringNorm.lengthSq() > 1e-6) {
              ringNorm.normalize();
              const torusGeo = new THREE.TorusGeometry(avgR * 0.62, 0.032, 10, 36);
              const torusMesh = new THREE.Mesh(torusGeo, hapticMat);
              torusMesh.position.copy(pCent);
              torusMesh.quaternion.setFromUnitVectors(new THREE.Vector3(0, 0, 1), ringNorm);
              contentGroup.add(torusMesh);
            }
          }
        }
      }

      // 3. Render Volumetric Field (MO / Density / MEP / Spin / Cube)
      const posHex = cbMode ? 0x0284c7 : 0x2563eb;
      const negHex = cbMode ? 0xea580c : 0xe11d48;
      const densityHex = 0x0d9488;
      const spinPosHex = 0x10b981;
      const spinNegHex = 0x8b5cf6;

      if (fm === 'mo' && mo) {
        setSpinEmptyState(false);
        setMepRange(null);

        const posMat = new THREE.MeshPhysicalMaterial({
          color: posHex,
          roughness: 0.24,
          metalness: 0.02,
          clearcoat: 0.25,
          transparent: true,
          opacity: opacityVal,
          wireframe: mode === 'wireframe',
          side: THREE.DoubleSide,
          depthWrite: false,
        });
        const negMat = new THREE.MeshPhysicalMaterial({
          color: negHex,
          roughness: 0.24,
          metalness: 0.02,
          clearcoat: 0.25,
          transparent: true,
          opacity: opacityVal,
          wireframe: mode === 'wireframe',
          side: THREE.DoubleSide,
          depthWrite: false,
        });

        const kf = getMorphKeyframe(mo, atomList, scaleMult, tProgress);
        const posMesh = new THREE.Mesh(kf.posGeo, posMat);
        posMesh.renderOrder = 3;
        const negMesh = new THREE.Mesh(kf.negGeo, negMat);
        negMesh.renderOrder = 3;
        activePosMesh = posMesh;
        activeNegMesh = negMesh;
        contentGroup.add(posMesh, negMesh);

        // Optional Nodal Planes in MO mode
        if (showNodes && mo.nodalPlanes > 0) {
          const planeSpan = Math.max(4.6, (rMol + 1.0) * 1.55);
          const nodePlaneGeo = new THREE.PlaneGeometry(planeSpan, Math.min(planeSpan, 6.4));
          const nodePlaneMat = new THREE.MeshBasicMaterial({
            color: isDark ? 0x94a3b8 : 0x475569,
            transparent: true,
            opacity: 0.14,
            side: THREE.DoubleSide,
            depthWrite: false,
          });
          const nodeBorderGeo = new THREE.EdgesGeometry(nodePlaneGeo);
          const nodeBorderMat = new THREE.LineBasicMaterial({
            color: isDark ? 0x38bdf8 : 0x0284c7,
            transparent: true,
            opacity: 0.45,
          });

          const addNodePlane = (normal: [number, number, number]) => {
            const pGroup = new THREE.Group();
            pGroup.position.set(cx, cy, cz);
            const mesh = new THREE.Mesh(nodePlaneGeo, nodePlaneMat);
            mesh.renderOrder = 1;
            const edges = new THREE.LineSegments(nodeBorderGeo, nodeBorderMat);
            edges.renderOrder = 1;
            pGroup.add(mesh, edges);
            pGroup.quaternion.setFromUnitVectors(
              new THREE.Vector3(0, 0, 1),
              new THREE.Vector3(...normal).normalize()
            );
            contentGroup.add(pGroup);
          };

          if (mo.interactionType === 'π' || mo.interactionType === 'π*') {
            addNodePlane([0, 0, 1]);
            if (mo.nodalPlanes >= 2) addNodePlane([1, 0, 0]);
          } else {
            const psiAlongX = Math.abs(evaluateOrbitalWavefunction(mo, atomList, 0.75, 0, 0));
            const psiAlongY = Math.abs(evaluateOrbitalWavefunction(mo, atomList, 0, 0.75, 0));
            const isYZeroNode = mo.nodalPlanes === 1 && psiAlongX > 0.02 && psiAlongY < 1e-4;
            addNodePlane(isYZeroNode ? [0, 1, 0] : [1, 0, 0]);
            if (mo.nodalPlanes >= 2) addNodePlane([0, 1, 0]);
          }
        }
      } else if (fm === 'density' || fm === 'mep' || fm === 'spin' || fm === 'cube') {
        // Density / MEP / Spin / Cube Volumetric Mode
        const surf = buildVolumetricSurfaceV3({
          fieldMode: fm,
          atoms: atomList,
          eht: ehtData,
          importedCube: cubeData,
          lobeScale: scaleMult,
        });

        setSpinEmptyState(fm === 'spin' && !surf.hasData);
        setMepRange(
          fm === 'mep' && surf.hasData && surf.mepRange
            ? [surf.mepRange.min, surf.mepRange.max]
            : null
        );

        if (surf.positive.positions.length > 0 && surf.positive.indices.length > 0) {
          const geo = buildBufferGeoFromMeshData(surf.positive);
          const isMep = fm === 'mep' && !!surf.positive.colors;
          const mat = isMep
            ? new THREE.MeshStandardMaterial({
                // Calibrated mid-albedo multiplier keeps saturated linear vertex colors inside
                // ACES Filmic's peak-chroma range instead of clipping toward pastel peach/lavender.
                color: 0xa6a6a6,
                vertexColors: true,
                roughness: 0.62,
                metalness: 0.0,
                envMapIntensity: 0.08,
                transparent: true,
                opacity: Math.min(0.88, opacityVal + 0.12),
                wireframe: mode === 'wireframe',
                side: mode === 'wireframe' ? THREE.DoubleSide : THREE.FrontSide,
                depthWrite: false,
              })
            : new THREE.MeshPhysicalMaterial({
                color:
                  fm === 'density'
                    ? densityHex
                    : fm === 'spin'
                      ? spinPosHex
                      : posHex,
                roughness: 0.24,
                metalness: 0.02,
                clearcoat: 0.2,
                transparent: true,
                opacity: opacityVal,
                wireframe: mode === 'wireframe',
                side: THREE.DoubleSide,
                depthWrite: false,
              });
          const mesh = new THREE.Mesh(geo, mat);
          mesh.renderOrder = 3;
          contentGroup.add(mesh);
        }

        if (surf.negative.positions.length > 0 && surf.negative.indices.length > 0) {
          const geo = buildBufferGeoFromMeshData(surf.negative);
          const mat = new THREE.MeshPhysicalMaterial({
            color: fm === 'spin' ? spinNegHex : negHex,
            roughness: 0.24,
            metalness: 0.02,
            clearcoat: 0.2,
            transparent: true,
            opacity: opacityVal,
            wireframe: mode === 'wireframe',
            side: THREE.DoubleSide,
            depthWrite: false,
          });
          const mesh = new THREE.Mesh(geo, mat);
          mesh.renderOrder = 3;
          contentGroup.add(mesh);
        }
      }

      // 4. Interactive 2D Contour & Phase Slice Plane in 3D Scene
      if (sliceOn && sliceRes) {
        const res = sliceRes.width;
        const offCanvas = document.createElement('canvas');
        offCanvas.width = res;
        offCanvas.height = sliceRes.height;
        const ctx = offCanvas.getContext('2d');
        if (ctx) {
          const img = ctx.createImageData(res, sliceRes.height);
          img.data.set(sliceRes.rgba);
          ctx.putImageData(img, 0, 0);
        }
        const tex = new THREE.CanvasTexture(offCanvas);
        tex.colorSpace = THREE.SRGBColorSpace;
        tex.anisotropy = maxAnisotropy;
        tex.needsUpdate = true;

        const sliceGeo = new THREE.PlaneGeometry(sliceRes.spanU, sliceRes.spanV);
        const sliceMat = new THREE.MeshBasicMaterial({
          map: tex,
          transparent: true,
          opacity: 0.88,
          side: THREE.DoubleSide,
          polygonOffset: true,
          polygonOffsetFactor: -1,
          polygonOffsetUnits: -1,
        });
        const sliceMesh = new THREE.Mesh(sliceGeo, sliceMat);
        sliceMesh.renderOrder = 4;

        // Orient plane to match u/v/normal definition from compute2DSlicePlane
        if (sliceRes.orientation === 'xy') {
          sliceMesh.position.set(sliceRes.center[0], sliceRes.center[1], sliceRes.center[2] + sliceRes.offset);
        } else if (sliceRes.orientation === 'xz') {
          sliceMesh.rotation.x = -Math.PI / 2;
          sliceMesh.position.set(sliceRes.center[0], sliceRes.center[1] + sliceRes.offset, sliceRes.center[2]);
        } else {
          sliceMesh.rotation.y = Math.PI / 2;
          sliceMesh.position.set(sliceRes.center[0] + sliceRes.offset, sliceRes.center[1], sliceRes.center[2]);
        }

        const borderGeo = new THREE.EdgesGeometry(sliceGeo);
        const borderMat = new THREE.LineBasicMaterial({
          color: isDark ? 0x38bdf8 : 0x0284c7,
          transparent: true,
          opacity: 0.85,
        });
        const border = new THREE.LineSegments(borderGeo, borderMat);
        sliceMesh.add(border);
        contentGroup.add(sliceMesh);
      }

      // 5. Optional 3D Mulliken Partial Charge Sprites
      if (showCharges && ehtData.mullikenAtoms.length === atomList.length) {
        for (let i = 0; i < atomList.length; i++) {
          const a = atomList[i];
          const ma = ehtData.mullikenAtoms[i];
          const label = `${ma.netCharge >= 0 ? '+' : ''}${ma.netCharge.toFixed(2)}`;
          const sp = createChargeSprite(label, ma.netCharge < -0.02, ma.netCharge > 0.02);
          sp.position.set(a.x3d, a.y3d + 0.56, a.z3d);
          contentGroup.add(sp);
        }
      }

      // 6. Optional 3D Molecular Dipole Moment Arrow (μ)
      if (showDipole && ehtData.dipoleMagnitudeDebye > 0.08) {
        const [dx, dy, dz] = ehtData.dipoleVectorDebye;
        const vec = new THREE.Vector3(dx, dy, dz);
        const dir = vec.clone().normalize();
        const arrowLen = Math.min(rMol * 1.15, Math.max(1.1, ehtData.dipoleMagnitudeDebye * 0.65));
        const arrow = new THREE.ArrowHelper(
          dir,
          new THREE.Vector3(cx, cy, cz),
          arrowLen,
          0xf59e0b,
          0.34,
          0.2
        );
        contentGroup.add(arrow);
      }

      // 7. 3D Point-Group Symmetry Elements (using molecule-aligned axes when available)
      const elementsToDraw =
        alignedSym && alignedSym.length > 0 ? alignedSym : pg.elements3D;
      if (showSym && elementsToDraw.length > 0) {
        const symGroup = new THREE.Group();
        symGroup.position.set(cx, cy, cz);
        const symAxisSpan = Math.max(5.0, (rMol + 1.0) * 1.9);
        const symPlaneSpan = Math.max(4.6, (rMol + 1.0) * 1.5);

        for (const el of elementsToDraw) {
          const isHighlighted = hiSymId === el.id;
          const baseColor = new THREE.Color(el.color);

          if (el.type === 'axis' && el.axis) {
            const dir = new THREE.Vector3(...el.axis).normalize();
            const rodGeo = new THREE.CylinderGeometry(
              isHighlighted ? 0.045 : 0.024,
              isHighlighted ? 0.045 : 0.024,
              symAxisSpan,
              12
            );
            const rodMat = new THREE.MeshStandardMaterial({
              color: baseColor,
              emissive: baseColor,
              emissiveIntensity: isHighlighted ? 0.6 : 0.25,
              roughness: 0.3,
            });
            const rod = new THREE.Mesh(rodGeo, rodMat);
            rod.quaternion.setFromUnitVectors(vUp, dir);
            symGroup.add(rod);

            const tipGeo = new THREE.ConeGeometry(0.11, 0.24, Math.max(3, el.fold || 4));
            const tip = new THREE.Mesh(tipGeo, rodMat);
            tip.position.copy(dir.clone().multiplyScalar(symAxisSpan / 2));
            tip.quaternion.setFromUnitVectors(vUp, dir);
            symGroup.add(tip);
          } else if (el.type === 'plane' && el.axis) {
            const normal = new THREE.Vector3(...el.axis).normalize();
            const planeGeo = new THREE.PlaneGeometry(symPlaneSpan, Math.min(symPlaneSpan, 6.4));
            const planeMat = new THREE.MeshBasicMaterial({
              color: baseColor,
              transparent: true,
              opacity: isHighlighted ? 0.26 : 0.08,
              side: THREE.DoubleSide,
              depthWrite: false,
            });
            const planeMesh = new THREE.Mesh(planeGeo, planeMat);
            planeMesh.renderOrder = 1;
            planeMesh.quaternion.setFromUnitVectors(new THREE.Vector3(0, 0, 1), normal);
            symGroup.add(planeMesh);
          } else if (el.type === 'center') {
            const invGeo = new THREE.OctahedronGeometry(isHighlighted ? 0.22 : 0.15, 0);
            const invMat = new THREE.MeshStandardMaterial({
              color: 0xf59e0b,
              emissive: 0xd97706,
              emissiveIntensity: 0.5,
              roughness: 0.2,
              metalness: 0.8,
            });
            symGroup.add(new THREE.Mesh(invGeo, invMat));
          }
        }
        contentGroup.add(symGroup);
      }

      controls.update();
      renderer.render(scene, camera);
    };

    rebuildSceneRef.current = buildMolecularAndOrbitalScene;
    buildMolecularAndOrbitalScene();

    if (!isSoftwareGl) {
      renderer.compile(scene, camera);
    }

    const onControlChange = () => requestRender();
    controls.addEventListener('change', onControlChange);

    resetCameraRef.current = (preset = 'iso') => {
      const atomList = stateRef.current.atoms;
      let cx = 0,
        cy = 0,
        cz = 0;
      for (const a of atomList) {
        cx += a.x3d;
        cy += a.y3d;
        cz += a.z3d;
      }
      const nA = Math.max(1, atomList.length);
      cx /= nA;
      cy /= nA;
      cz /= nA;
      let maxR = 1.6;
      for (const a of atomList) {
        const r = Math.hypot(a.x3d - cx, a.y3d - cy, a.z3d - cz);
        if (r > maxR) maxR = r;
      }
      const dist = Math.max(6.2, (maxR + 0.6) * 2.95);
      controls.target.set(cx, cy, cz);
      if (preset === 'top') camera.position.set(cx, cy + dist, cz + 0.01);
      else if (preset === 'front') camera.position.set(cx, cy, cz + dist);
      else if (preset === 'side') camera.position.set(cx + dist, cy, cz);
      else camera.position.set(cx + dist * 0.15, cy + dist * 0.34, cz + dist * 0.94);
      controls.update();
      requestRender();
    };

    const ro = new ResizeObserver(() => {
      if (!container || isDisposed) return;
      const w = container.clientWidth;
      const h = container.clientHeight;
      if (w > 20 && h > 20) {
        camera.aspect = w / h;
        camera.updateProjectionMatrix();
        renderer.setSize(w, h);
        requestRender();
      }
    });
    ro.observe(container);

    return () => {
      isDisposed = true;
      ro.disconnect();
      controls.removeEventListener('change', onControlChange);
      controls.dispose();
      disposeContentGroup();
      clearMorphCache();
      shadowGeo.dispose();
      shadowMat.dispose();
      ringGeo.dispose();
      ringMat.dispose();
      renderer.dispose();
      if (renderer.domElement.parentNode === container) {
        container.removeChild(renderer.domElement);
      }
    };
  }, []);

  // Rebuild 3D meshes whenever chemistry or visual controls change
  useEffect(() => {
    rebuildSceneRef.current?.();
  }, [
    molecule,
    atoms,
    selectedMo,
    pointGroup,
    alignedSymmetryElements,
    eht,
    fieldMode,
    importedCube,
    highlightedSymmetryElementId,
    lobeScale,
    orbitalOpacity,
    renderMode,
    showSymmetryElements,
    showNodalPlanes,
    showChargeLabels,
    showDipoleArrow,
    darkStage,
    colorblindPalette,
    showSlicePlane,
    sliceData,
  ]);

  const fieldTabs: { id: VolumetricFieldMode; label: string }[] = [
    { id: 'mo', label: 'MO ψ' },
    { id: 'density', label: 'Density ρ' },
    { id: 'mep', label: 'MEP' },
    { id: 'spin', label: 'Spin' },
    ...(importedCube ? [{ id: 'cube' as const, label: 'Cube' }] : []),
  ];

  return (
    <div
      data-testid="orbital-visualizer-v3"
      className="bg-white border border-slate-200 rounded-lg flex flex-col h-full overflow-hidden"
    >
      {/* Compact Header with Volumetric Field Mode Switcher */}
      <div className="px-3 py-2 border-b border-slate-200 flex flex-wrap items-center justify-between gap-2 bg-white">
        <div className="flex items-center gap-1.5 min-w-0">
          <span className="text-xs font-semibold text-slate-900 whitespace-nowrap">
            3D Field
          </span>
          <span className="text-slate-300" aria-hidden="true">·</span>
          <span className="text-[11px] font-mono text-slate-600 truncate">
            <Latex math={formulaToLatex(molecule.formula)} />{' '}
            (<Latex math={schoenfliesToLatex(pointGroup.schoenflies)} />)
          </span>
        </div>

        {/* Volumetric Field Selector Tabs */}
        <div
          className="flex items-center gap-0.5 bg-slate-100 p-0.5 rounded"
          role="tablist"
          aria-label="3D Volumetric Field Mode"
        >
          {fieldTabs.map((tab) => {
            const active = fieldMode === tab.id;
            return (
              <button
                key={tab.id}
                type="button"
                role="tab"
                aria-selected={active}
                data-testid={`field-mode-${tab.id}`}
                onClick={() => onChangeFieldMode(tab.id)}
                className={`px-2 py-0.5 text-[11px] font-mono font-medium rounded transition-colors whitespace-nowrap shrink-0 ${
                  active
                    ? 'bg-white text-slate-900 shadow-2xs'
                    : 'text-slate-600 hover:text-slate-900'
                }`}
              >
                {tab.label}
              </button>
            );
          })}
        </div>
      </div>

      {/* Main 3D Viewport Stage */}
      <div className="relative flex-1 min-h-0 bg-[#0b0f19]">
        <div ref={canvasMountRef} className="w-full h-full" />

        {/* Top-Left Compact Telemetry Overlay (Solid high-contrast, no backdrop-blur) */}
        <div className="absolute top-2.5 left-2.5 pointer-events-none flex flex-col gap-1 max-w-[65%]">
          <div
            className={`px-2.5 py-1 rounded text-[11px] font-mono border ${
              darkStage
                ? 'bg-slate-950/90 text-slate-100 border-slate-700'
                : 'bg-white/95 text-slate-900 border-slate-300'
            }`}
          >
            {fieldMode === 'mo' ? (
              <span>
                <Latex math={moLabelToLatex(selectedMo.label)} /> (
                <Latex math={irrepToLatex(selectedMo.irrepSymbol)} />) ·{' '}
                {selectedMo.energy.toFixed(2)} eV
              </span>
            ) : fieldMode === 'density' ? (
              <span>
                Total Density <Latex math="\rho(\mathbf{r}) = \sum n_i|\psi_i|^2" />
              </span>
            ) : fieldMode === 'mep' ? (
              <span>
                MEP <Latex math="V(\mathbf{r})" /> on <Latex math="\rho" /> isosurface
              </span>
            ) : fieldMode === 'spin' ? (
              <span>
                Spin Density <Latex math="\rho_\alpha(\mathbf{r}) - \rho_\beta(\mathbf{r})" />
              </span>
            ) : (
              <span>Imported Gaussian .cube ({importedCube?.title || 'Volume'})</span>
            )}
          </div>

          {walshBanner && (
            <div className="px-2.5 py-1 rounded text-[11px] font-mono bg-blue-700 text-white border border-blue-500">
              Walsh Live: {walshBanner}
            </div>
          )}
        </div>

        {/* Top-Right Quick Viewport Controls: Slice Toggle + Stage Theme + Camera */}
        <div className="absolute top-2.5 right-2.5 flex items-center gap-1">
          <button
            type="button"
            data-testid="toggle-2d-slice"
            onClick={() => setShowSlicePlane((s) => !s)}
            title="Toggle 2D Contour & Phase Slice Plane"
            className={`px-2 py-1 text-[11px] font-mono font-medium rounded border flex items-center gap-1 transition-colors whitespace-nowrap shrink-0 ${
              showSlicePlane
                ? 'bg-blue-600 text-white border-blue-500'
                : darkStage
                  ? 'bg-slate-900/90 text-slate-200 border-slate-700 hover:bg-slate-800'
                  : 'bg-white/95 text-slate-800 border-slate-300 hover:bg-slate-100'
            }`}
          >
            <Grid2x2 className="w-3 h-3" />
            <span>2D Slice</span>
          </button>

          <button
            type="button"
            onClick={() => setDarkStage((d) => !d)}
            title="Toggle Dark / Light 3D Stage"
            className={`p-1 rounded border transition-colors ${
              darkStage
                ? 'bg-slate-900/90 text-amber-300 border-slate-700 hover:bg-slate-800'
                : 'bg-white/95 text-slate-700 border-slate-300 hover:bg-slate-100'
            }`}
          >
            <Sun className="w-3.5 h-3.5" />
          </button>

          <button
            type="button"
            onClick={() => resetCameraRef.current?.('iso')}
            className={`px-2 py-1 text-[11px] font-mono rounded border transition-colors ${
              darkStage
                ? 'bg-slate-900/90 text-slate-200 border-slate-700 hover:bg-slate-800'
                : 'bg-white/95 text-slate-800 border-slate-300 hover:bg-slate-100'
            }`}
          >
            Reset
          </button>
        </div>

        {/* Spin Density Closed-Shell Empty State */}
        {spinEmptyState && (
          <div
            data-testid="spin-empty-banner"
            className="absolute inset-x-3 bottom-3 bg-slate-950/95 border border-slate-700 rounded-lg p-3 text-xs text-slate-200 space-y-2"
          >
            <div className="font-semibold text-white">
              Closed-Shell Singlet: Zero Net Spin Density (<Latex math="\rho_\alpha - \rho_\beta \equiv 0" />)
            </div>
            <p className="text-[11px] text-slate-400">
              All occupied molecular orbitals are doubly paired. Load an open-shell triplet/radical or remove one electron to inspect unpaired spin density.
            </p>
            <div className="flex flex-wrap items-center gap-2">
              {onInspectOpenShellO2 && (
                <button
                  type="button"
                  onClick={onInspectOpenShellO2}
                  className="px-2.5 py-1 text-[11px] font-medium bg-teal-600 hover:bg-teal-700 text-white rounded transition-colors"
                >
                  Inspect O₂ (³Σg⁻ Triplet)
                </button>
              )}
              {onIonizeRadical && (
                <button
                  type="button"
                  onClick={onIonizeRadical}
                  className="px-2.5 py-1 text-[11px] font-medium bg-slate-800 hover:bg-slate-700 text-slate-100 border border-slate-600 rounded transition-colors"
                >
                  Ionize (+1 Radical Cation)
                </button>
              )}
            </div>
          </div>
        )}

        {/* MEP Calibrated Electrostatic Potential Legend */}
        {fieldMode === 'mep' && mepRange && (
          <div
            data-testid="mep-legend"
            className="absolute bottom-2 left-2 bg-slate-950/92 border border-slate-700 rounded px-2 py-1 text-[9.5px] font-mono text-slate-200 space-y-0.5"
          >
            <div className="flex items-center justify-between gap-2">
              <span className="text-rose-400">δ−</span>
              <span className="text-slate-400">MEP</span>
              <span className="text-blue-400">δ+</span>
            </div>
            <div
              className="h-1.5 w-28 rounded-xs"
              style={{
                background:
                  'linear-gradient(90deg, #D6080F 0%, #F77300 25%, #F2CC0D 41%, #0FAD47 50%, #00B8C7 59%, #0070F0 75%, #0F24DB 100%)',
              }}
            />
            <div className="flex items-center justify-between text-slate-400 tabular-nums">
              <span>{mepRange[0].toFixed(2)}</span>
              <span>+{mepRange[1].toFixed(2)}</span>
            </div>
          </div>
        )}

        {/* 2D Contour & Phase Slice Inset with Live Value Probing */}
        {showSlicePlane && sliceData && (
          <div
            data-testid="slice-2d-inset"
            className="absolute bottom-2 right-2 bg-slate-950/95 border border-slate-700 rounded-md p-1.5 text-[9.5px] font-mono text-slate-200 w-[132px] space-y-1"
          >
            <div className="flex items-center justify-between">
              <span className="font-semibold text-sky-400 uppercase">
                {sliceOrientation} ({sliceOffset >= 0 ? '+' : ''}
                {sliceOffset.toFixed(1)}Å)
              </span>
              <span className="text-slate-400">
                ±{(sliceData.spanU * 0.5).toFixed(1)}
              </span>
            </div>
            <canvas
              ref={sliceInsetCanvasRef}
              className="w-[118px] h-[118px] rounded-xs border border-slate-700 cursor-crosshair block"
              onMouseMove={(e) => {
                const rect = e.currentTarget.getBoundingClientRect();
                const fx = Math.max(0, Math.min(0.999, (e.clientX - rect.left) / rect.width));
                const fy = Math.max(0, Math.min(0.999, 1 - (e.clientY - rect.top) / rect.height));
                const ix = Math.floor(fx * sliceData.width);
                const iy = Math.floor(fy * sliceData.height);
                const val = sliceData.values[iy * sliceData.width + ix] ?? 0;
                const u = (fx - 0.5) * sliceData.spanU;
                const v = (fy - 0.5) * sliceData.spanV;
                setSliceProbe({ u, v, val });
              }}
              onMouseLeave={() => setSliceProbe(null)}
            />
            <div className="text-[9.5px] text-slate-300 tabular-nums truncate">
              {sliceProbe
                ? `(${sliceProbe.u.toFixed(1)},${sliceProbe.v.toFixed(1)}) ${sliceProbe.val >= 0 ? '+' : ''}${sliceProbe.val.toFixed(3)}`
                : `Peak |f|=${sliceData.absMax.toFixed(3)}`}
            </div>
          </div>
        )}
      </div>

      {/* Inline 2D Slice Plane Controls Strip (visible only when 2D Slice is active) */}
      {showSlicePlane && (
        <div className="px-3 py-1.5 bg-slate-50 border-t border-slate-200 flex flex-wrap items-center justify-between gap-2 text-[11px] font-mono">
          <div className="flex items-center gap-1">
            <span className="text-slate-500">Plane:</span>
            {(['xy', 'xz', 'yz'] as SliceOrientation[]).map((ori) => (
              <button
                key={ori}
                type="button"
                onClick={() => setSliceOrientation(ori)}
                className={`px-1.5 py-0.5 rounded uppercase border transition-colors ${
                  sliceOrientation === ori
                    ? 'bg-blue-600 text-white border-blue-600'
                    : 'bg-white text-slate-700 border-slate-300 hover:bg-slate-100'
                }`}
              >
                {ori}
              </button>
            ))}
          </div>

          <div className="flex items-center gap-1.5 flex-1 min-w-[120px] max-w-[200px]">
            <span className="text-slate-500 shrink-0">Offset:</span>
            <input
              type="range"
              aria-label="2D slice plane offset"
              min={-2.5}
              max={2.5}
              step={0.1}
              value={sliceOffset}
              onChange={(e) => setSliceOffset(Number(e.target.value))}
              className="w-full accent-blue-600"
            />
          </div>

          <label className="inline-flex items-center gap-1 text-slate-700 cursor-pointer">
            <input
              type="checkbox"
              checked={sliceIsolines}
              onChange={(e) => setSliceIsolines(e.target.checked)}
              className="accent-blue-600"
            />
            <span>Isolines</span>
          </label>
        </div>
      )}

      {/* Compact Bottom Controls Bar (AO->MO Morph + Collapsible Display Options) */}
      <div className="px-3 py-2 bg-white border-t border-slate-200 space-y-2">
        <div className="flex flex-wrap items-center justify-between gap-2">
          {fieldMode === 'mo' ? (
            <div className="flex items-center gap-1.5 flex-1 min-w-[160px]">
              <button
                type="button"
                onClick={handleTogglePlayMorph}
                className="px-2 py-1 text-[11px] font-medium bg-blue-600 hover:bg-blue-700 text-white rounded flex items-center gap-1 whitespace-nowrap shrink-0 transition-colors"
              >
                {isMorphPlaying ? <Pause className="w-3 h-3" /> : <Play className="w-3 h-3" />}
                <span>AO→MO</span>
              </button>
              <input
                type="range"
                aria-label="AO to MO fusion morph"
                min={0}
                max={1}
                step={0.025}
                value={morphProgress}
                onChange={(e) => handleScrubMorph(Number(e.target.value))}
                className="flex-1 accent-blue-600"
              />
              <span className="text-[11px] font-mono text-slate-600 tabular-nums w-8 text-right">
                {Math.round(morphProgress * 100)}%
              </span>
            </div>
          ) : (
            <div className="flex items-center gap-2 flex-1 min-w-[160px] text-[11px] font-mono text-slate-600">
              <span>Iso Scale:</span>
              <input
                type="range"
                aria-label="Volumetric isosurface scale"
                min={0.5}
                max={1.6}
                step={0.05}
                value={lobeScale}
                onChange={(e) => setLobeScale(Number(e.target.value))}
                className="flex-1 accent-blue-600 max-w-[140px]"
              />
              <span className="tabular-nums">{lobeScale.toFixed(2)}×</span>
            </div>
          )}

          <button
            type="button"
            onClick={() => setDockCollapsed((c) => !c)}
            className="px-2 py-1 text-[11px] font-medium text-slate-700 bg-slate-100 border border-slate-200 rounded hover:bg-slate-200 flex items-center gap-1 whitespace-nowrap shrink-0 transition-colors"
          >
            <span>Display</span>
            {dockCollapsed ? <ChevronUp className="w-3 h-3" /> : <ChevronDown className="w-3 h-3" />}
          </button>
        </div>

        {!dockCollapsed && (
          <div className="pt-2 border-t border-slate-100 grid grid-cols-2 gap-x-3 gap-y-2 text-[11px] text-slate-700">
            <label className="flex items-center justify-between gap-2">
              <span className="font-mono text-slate-500">Opacity</span>
              <input
                type="range"
                min={0.25}
                max={0.95}
                step={0.05}
                value={orbitalOpacity}
                onChange={(e) => setOrbitalOpacity(Number(e.target.value))}
                className="w-24 accent-blue-600"
              />
            </label>

            <label className="flex items-center justify-between gap-2">
              <span className="font-mono text-slate-500">Mesh</span>
              <select
                value={renderMode}
                onChange={(e) => setRenderMode(e.target.value as 'merged' | 'wireframe')}
                className="text-[11px] font-mono bg-slate-50 border border-slate-200 rounded px-1.5 py-0.5"
              >
                <option value="merged">Smooth Glass</option>
                <option value="wireframe">Wireframe</option>
              </select>
            </label>

            <label className="inline-flex items-center gap-1.5 cursor-pointer">
              <input
                type="checkbox"
                checked={showSymmetryElements}
                onChange={(e) => setShowSymmetryElements(e.target.checked)}
                className="accent-blue-600"
              />
              <span>Symmetry Axes/Planes</span>
            </label>

            <label className="inline-flex items-center gap-1.5 cursor-pointer">
              <input
                type="checkbox"
                checked={showDipoleArrow}
                onChange={(e) => setShowDipoleArrow(e.target.checked)}
                className="accent-blue-600"
              />
              <span>Dipole Vector (μ)</span>
            </label>

            <label className="inline-flex items-center gap-1.5 cursor-pointer">
              <input
                type="checkbox"
                checked={showChargeLabels}
                onChange={(e) => setShowChargeLabels(e.target.checked)}
                className="accent-blue-600"
              />
              <span>Mulliken Charges (q)</span>
            </label>

            <label className="inline-flex items-center gap-1.5 cursor-pointer">
              <input
                type="checkbox"
                checked={showNodalPlanes}
                onChange={(e) => setShowNodalPlanes(e.target.checked)}
                className="accent-blue-600"
              />
              <span>Nodal Planes (ψ=0)</span>
            </label>

            <label className="inline-flex items-center gap-1.5 cursor-pointer">
              <input
                type="checkbox"
                checked={colorblindPalette}
                onChange={(e) => setColorblindPalette(e.target.checked)}
                className="accent-blue-600"
              />
              <span>Colorblind Phase</span>
            </label>

            <div className="flex items-center gap-1 justify-end">
              {(['top', 'front', 'side'] as const).map((p) => (
                <button
                  key={p}
                  type="button"
                  onClick={() => resetCameraRef.current?.(p)}
                  className="px-1.5 py-0.5 text-[10px] font-mono bg-slate-100 hover:bg-slate-200 rounded uppercase"
                >
                  {p}
                </button>
              ))}
            </div>
          </div>
        )}
      </div>
    </div>
  );
};
