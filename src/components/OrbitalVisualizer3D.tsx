import React, { useEffect, useRef, useState } from 'react';
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { RoomEnvironment } from 'three/addons/environments/RoomEnvironment.js';
import {
  ELEMENT_PROPERTIES,
  MolecularOrbital,
  MoleculeEntry,
} from '../chemistry/moleculesData';
import { evaluateOrbitalWavefunction } from '../chemistry/moEngine';
import { buildMorphIsosurfaceKeyframe } from '../chemistry/moIsosurface';
import { PointGroupData } from '../chemistry/symmetryData';
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
  Pause,
  Play,
  Sun,
} from 'lucide-react';

interface OrbitalVisualizer3DProps {
  molecule: MoleculeEntry;
  selectedMo: MolecularOrbital;
  pointGroup: PointGroupData;
  highlightedSymmetryElementId?: string | null;
}

const MORPH_KEYFRAME_STEPS = 40; // 41 high-density keyframes (2.5% increments: 0/40 .. 40/40)

export const OrbitalVisualizer3D: React.FC<OrbitalVisualizer3DProps> = ({
  molecule,
  selectedMo,
  pointGroup,
  highlightedSymmetryElementId = null,
}) => {
  const canvasMountRef = useRef<HTMLDivElement | null>(null);
  const requestRenderRef = useRef<(() => void) | null>(null);
  const resetCameraRef = useRef<((preset?: 'iso' | 'top' | 'front' | 'side') => void) | null>(null);
  const rebuildSceneRef = useRef<(() => void) | null>(null);
  const updateMorphFastRef = useRef<((tProgress: number) => void) | null>(null);

  const [lobeScale, setLobeScale] = useState<number>(1.0);
  const [orbitalOpacity, setOrbitalOpacity] = useState<number>(0.68);
  const [renderMode, setRenderMode] = useState<'lobes' | 'merged' | 'cloud' | 'wireframe'>('merged');
  const [morphProgress, setMorphProgress] = useState<number>(1.0);
  const [isMorphPlaying, setIsMorphPlaying] = useState<boolean>(false);
  const [showSymmetryElements, setShowSymmetryElements] = useState<boolean>(true);
  const [showNodalPlanes, setShowNodalPlanes] = useState<boolean>(true);
  const [darkStage, setDarkStage] = useState<boolean>(true);
  const [colorblindPalette, setColorblindPalette] = useState<boolean>(false);
  const [dockCollapsed, setDockCollapsed] = useState<boolean>(true);

  const morphStartLinearRef = useRef<number>(0);

  // Store mutable props in ref so WebGL lifecycle never tears down on prop change
  const stateRef = useRef({
    molecule,
    selectedMo,
    pointGroup,
    highlightedSymmetryElementId,
    lobeScale,
    orbitalOpacity,
    renderMode,
    morphProgress,
    isMorphPlaying,
    showSymmetryElements,
    showNodalPlanes,
    darkStage,
    colorblindPalette,
  });

  stateRef.current = {
    molecule,
    selectedMo,
    pointGroup,
    highlightedSymmetryElementId,
    lobeScale,
    orbitalOpacity,
    renderMode,
    morphProgress,
    isMorphPlaying,
    showSymmetryElements,
    showNodalPlanes,
    darkStage,
    colorblindPalette,
  };

  // Stop active playback when user switches molecule
  useEffect(() => {
    setIsMorphPlaying(false);
  }, [molecule.id]);

  // Smooth 60 FPS cosine ease-in-out playback loop for AO -> MO morphing animation
  useEffect(() => {
    if (!isMorphPlaying) return;

    let rafId = 0;
    const startWallTime = performance.now();
    const startU = morphStartLinearRef.current;
    const totalDurationMs = 2400; // 2.4s silky fusion transition
    const remainingDurationMs = Math.max(200, (1 - startU) * totalDurationMs);

    const stepAnimation = (now: number) => {
      const elapsed = Math.max(0, now - startWallTime);
      const frac = Math.min(1, elapsed / remainingDurationMs);
      const u = startU + frac * (1 - startU);
      // Cosine ease-in-out curve for natural acceleration and gentle settling into the merged MO
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

    return () => {
      cancelAnimationFrame(rafId);
    };
  }, [isMorphPlaying]);

  const handleTogglePlayMorph = () => {
    if (renderMode !== 'merged' && renderMode !== 'wireframe') {
      setRenderMode('merged');
    }
    if (!isMorphPlaying) {
      // If already near 100% (Merged MO), restart from 0% (Separate AOs)
      const startVal = morphProgress >= 0.96 ? 0 : morphProgress;
      // Invert cosine ease-in-out so resuming mid-scrub starts seamlessly without a jump
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
    if (renderMode !== 'merged' && renderMode !== 'wireframe') {
      setRenderMode('merged');
    }
    stateRef.current.morphProgress = nextValue;
    updateMorphFastRef.current?.(nextValue);
    setMorphProgress(nextValue);
  };

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

    const width = container.clientWidth || 640;
    const height = container.clientHeight || 420;

    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0x0b0f19);
    scene.fog = new THREE.FogExp2(0x0b0f19, 0.035);

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

    // Studio Environment Reflection (RoomEnvironment)
    const pmremGenerator = new THREE.PMREMGenerator(renderer);
    scene.environment = pmremGenerator.fromScene(new RoomEnvironment(), 0.04).texture;
    pmremGenerator.dispose();

    // 3-Point Warm/Cool Studio Lighting
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

    // Contact shadow & Subtle Reference Coordinate Floor
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

    // OrbitControls
    const controls = new OrbitControls(camera, renderer.domElement);
    controls.enableDamping = !isSoftwareGl;
    controls.dampingFactor = 0.08;
    controls.minDistance = 2.2;
    controls.maxDistance = 16;

    // Root dynamic content group
    const contentGroup = new THREE.Group();
    scene.add(contentGroup);

    // Shared geometries for performance
    const sphereGeo = new THREE.SphereGeometry(1, 28, 22);
    const torusCollarGeo = new THREE.TorusGeometry(1.03, 0.07, 12, 28);
    const cylinderGeo = new THREE.CylinderGeometry(1, 1, 1, 18);

    // Deformed teardrop lobe geometry for discrete LCAO lobes mode
    const lobeGeo = new THREE.SphereGeometry(0.65, 26, 20);
    const posAttr = lobeGeo.attributes.position;
    const vTemp = new THREE.Vector3();
    for (let i = 0; i < posAttr.count; i++) {
      vTemp.fromBufferAttribute(posAttr, i);
      const ny = (vTemp.y + 0.65) / 1.3;
      const pinch = Math.pow(ny, 0.55) * 1.15;
      vTemp.x *= pinch;
      vTemp.z *= pinch;
      vTemp.y = ny * 1.35;
      posAttr.setXYZ(i, vTemp.x, vTemp.y, vTemp.z);
    }
    lobeGeo.computeVertexNormals();

    // Pre-cached Marching Cubes / Tetrahedra keyframe geometries for active orbital
    let activeCacheScopeKey = '';
    let prewarmTimerId = 0;
    let activePosIsoMesh: THREE.Mesh | null = null;
    let activeNegIsoMesh: THREE.Mesh | null = null;

    const morphKeyframeCache = new Map<
      number,
      { posGeo: THREE.BufferGeometry; negGeo: THREE.BufferGeometry }
    >();

    const clearMorphKeyframeCache = () => {
      if (prewarmTimerId) {
        clearTimeout(prewarmTimerId);
        prewarmTimerId = 0;
      }
      for (const entry of morphKeyframeCache.values()) {
        entry.posGeo.dispose();
        entry.negGeo.dispose();
      }
      morphKeyframeCache.clear();
    };

    const buildKeyframeStep = (
      mol: MoleculeEntry,
      mo: MolecularOrbital,
      scaleMult: number,
      stepIdx: number
    ) => {
      const existing = morphKeyframeCache.get(stepIdx);
      if (existing) return existing;

      const quantizedT = stepIdx / MORPH_KEYFRAME_STEPS;
      const isoData = buildMorphIsosurfaceKeyframe(mo, mol.atoms, quantizedT, scaleMult);

      const posGeo = new THREE.BufferGeometry();
      if (isoData.positive.positions.length > 0) {
        posGeo.setAttribute(
          'position',
          new THREE.BufferAttribute(isoData.positive.positions, 3)
        );
        posGeo.setAttribute(
          'normal',
          new THREE.BufferAttribute(isoData.positive.normals, 3)
        );
        posGeo.setIndex(new THREE.BufferAttribute(isoData.positive.indices, 1));
        posGeo.computeBoundingSphere();
      }

      const negGeo = new THREE.BufferGeometry();
      if (isoData.negative.positions.length > 0) {
        negGeo.setAttribute(
          'position',
          new THREE.BufferAttribute(isoData.negative.positions, 3)
        );
        negGeo.setAttribute(
          'normal',
          new THREE.BufferAttribute(isoData.negative.normals, 3)
        );
        negGeo.setIndex(new THREE.BufferAttribute(isoData.negative.indices, 1));
        negGeo.computeBoundingSphere();
      }

      const pair = { posGeo, negGeo };
      morphKeyframeCache.set(stepIdx, pair);
      return pair;
    };

    const scheduleKeyframePrewarm = (
      mol: MoleculeEntry,
      mo: MolecularOrbital,
      scaleMult: number
    ) => {
      if (prewarmTimerId) {
        clearTimeout(prewarmTimerId);
        prewarmTimerId = 0;
      }
      let nextStep = 0;
      const warmNextSlice = () => {
        const start = performance.now();
        while (nextStep <= MORPH_KEYFRAME_STEPS && performance.now() - start < 8) {
          if (!morphKeyframeCache.has(nextStep)) {
            buildKeyframeStep(mol, mo, scaleMult, nextStep);
          }
          nextStep++;
        }
        if (nextStep <= MORPH_KEYFRAME_STEPS) {
          prewarmTimerId = window.setTimeout(warmNextSlice, 12);
        } else {
          prewarmTimerId = 0;
        }
      };
      prewarmTimerId = window.setTimeout(warmNextSlice, 24);
    };

    const getOrCreateMorphKeyframeGeometries = (
      mol: MoleculeEntry,
      mo: MolecularOrbital,
      scaleMult: number,
      tProgress: number
    ) => {
      const atomSig = mol.atoms
        .slice(0, 8)
        .map((a) => `${a.id}:${a.x3d.toFixed(2)},${a.y3d.toFixed(2)},${a.z3d.toFixed(2)}`)
        .join('|');
      const scopeKey = `${mol.id}_${mol.atoms.length}_${atomSig}__${mo.id}_${mo.label}__${scaleMult.toFixed(2)}`;
      if (scopeKey !== activeCacheScopeKey) {
        clearMorphKeyframeCache();
        activeCacheScopeKey = scopeKey;
        scheduleKeyframePrewarm(mol, mo, scaleMult);
      }

      const stepIdx = Math.max(
        0,
        Math.min(MORPH_KEYFRAME_STEPS, Math.round(tProgress * MORPH_KEYFRAME_STEPS))
      );
      return buildKeyframeStep(mol, mo, scaleMult, stepIdx);
    };

    let frameId = 0;
    let renderBusy = false;
    let renderDirty = false;

    function afterDraw(cb: () => void) {
      const gl = renderer.getContext() as WebGL2RenderingContext;
      if (!gl || !gl.fenceSync) {
        frameId = window.setTimeout(cb, 80);
        return;
      }
      const sync = gl.fenceSync(gl.SYNC_GPU_COMMANDS_COMPLETE, 0);
      gl.flush();
      const poll = () => {
        if (!sync) {
          cb();
          return;
        }
        const state = gl.getSyncParameter(sync, gl.SYNC_STATUS);
        if (state === gl.UNSIGNALED) {
          frameId = window.setTimeout(poll, 5);
        } else {
          gl.deleteSync(sync);
          frameId = window.setTimeout(cb, 80);
        }
      };
      frameId = window.setTimeout(poll, 5);
    }

    function requestRender() {
      if (renderBusy) {
        renderDirty = true;
        return;
      }
      renderBusy = true;
      frameId = requestAnimationFrame(() => {
        renderDirty = false;
        controls.update();
        renderer.render(scene, camera);
        const done = () => {
          renderBusy = false;
          if (renderDirty) requestRender();
        };
        if (isSoftwareGl) afterDraw(done);
        else done();
      });
    }

    requestRenderRef.current = requestRender;

    // Fast O(1) morph geometry swap without rebuilding atoms, bonds, or symmetry elements
    updateMorphFastRef.current = (tProgress: number) => {
      const {
        molecule: mol,
        selectedMo: mo,
        lobeScale: scaleMult,
        renderMode: mode,
      } = stateRef.current;

      if (
        !mo ||
        (mode !== 'merged' && mode !== 'wireframe') ||
        !activePosIsoMesh ||
        !activeNegIsoMesh
      ) {
        return;
      }

      const { posGeo, negGeo } = getOrCreateMorphKeyframeGeometries(
        mol,
        mo,
        scaleMult,
        tProgress
      );

      activePosIsoMesh.geometry = posGeo;
      activePosIsoMesh.visible =
        !!posGeo.attributes.position && posGeo.attributes.position.count > 0;

      activeNegIsoMesh.geometry = negGeo;
      activeNegIsoMesh.visible =
        !!negGeo.attributes.position && negGeo.attributes.position.count > 0;

      controls.update();
      renderer.render(scene, camera);
    };

    let lastFramedMoleculeId = '';
    let currentCamDist = 6.2;

    resetCameraRef.current = (preset = 'iso') => {
      const d = currentCamDist;
      if (preset === 'top') {
        camera.position.set(0, d * 1.05, 0.01);
      } else if (preset === 'front') {
        camera.position.set(0, d * 0.08, d * 1.05);
      } else if (preset === 'side') {
        camera.position.set(d * 1.05, d * 0.08, 0.01);
      } else {
        camera.position.set(d * 0.34, d * 0.32, d * 0.94);
      }
      controls.target.set(0, 0, 0);
      controls.update();
      requestRender();
    };

    const buildMolecularAndOrbitalScene = () => {
      activePosIsoMesh = null;
      activeNegIsoMesh = null;
      while (contentGroup.children.length > 0) {
        const child = contentGroup.children[0];
        contentGroup.remove(child);
      }

      const {
        molecule: mol,
        selectedMo: mo,
        pointGroup: pg,
        highlightedSymmetryElementId: hiSymId,
        lobeScale: scaleMult,
        orbitalOpacity: opacityVal,
        renderMode: mode,
        morphProgress: tMorph,
        showSymmetryElements: showSym,
        showNodalPlanes: showNodes,
        darkStage: isDark,
        colorblindPalette: isCb,
      } = stateRef.current;

      const rMol = Math.max(
        2.2,
        ...mol.atoms.map((a) => Math.hypot(a.x3d, a.y3d, a.z3d))
      );
      const minY3d =
        mol.atoms.length > 0 ? Math.min(...mol.atoms.map((a) => a.y3d)) : -1.2;

      const rOrbital = rMol + 1.25 * Math.max(1.0, scaleMult);
      currentCamDist = Math.max(7.4, rOrbital * 2.05);
      controls.maxDistance = Math.max(20, currentCamDist * 2.6);

      const floorScale = Math.max(1.15, rOrbital / 2.8);
      shadowDisc.scale.setScalar(floorScale);
      floorRing.scale.setScalar(floorScale);
      const floorY = Math.min(-2.65, minY3d - 1.65);
      shadowDisc.position.y = floorY - 0.02;
      floorRing.position.y = floorY;

      const bgHex = isDark ? 0x0b0f19 : 0xf1f5f9;
      scene.background = new THREE.Color(bgHex);
      scene.fog = new THREE.FogExp2(bgHex, Math.min(0.032, 0.18 / currentCamDist));

      const molKey = `${mol.id}_${mol.atoms.length}`;
      if (molKey !== lastFramedMoleculeId) {
        lastFramedMoleculeId = molKey;
        camera.position.set(
          currentCamDist * 0.28,
          currentCamDist * 0.30,
          currentCamDist * 0.95
        );
        controls.target.set(0, 0, 0);
        controls.update();
      }

      // 1. Build 3D Composite Atoms
      const atomGroup = new THREE.Group();
      for (const atom of mol.atoms) {
        const props = ELEMENT_PROPERTIES[atom.element] || ELEMENT_PROPERTIES.C;
        const radius = Math.max(0.22, Math.min(0.46, props.covalentRadius * 0.45));
        const singleAtomGroup = new THREE.Group();
        singleAtomGroup.position.set(atom.x3d, atom.y3d, atom.z3d);

        const matMain = new THREE.MeshStandardMaterial({
          color: new THREE.Color(props.cpkColor),
          roughness: 0.25,
          metalness: 0.12,
        });
        const sphereMesh = new THREE.Mesh(sphereGeo, matMain);
        sphereMesh.scale.setScalar(radius);
        singleAtomGroup.add(sphereMesh);

        const matCollar = new THREE.MeshStandardMaterial({
          color: 0x64748b,
          roughness: 0.55,
          metalness: 0.75,
        });
        const collarMesh = new THREE.Mesh(torusCollarGeo, matCollar);
        collarMesh.scale.setScalar(radius);
        singleAtomGroup.add(collarMesh);

        atomGroup.add(singleAtomGroup);
      }
      contentGroup.add(atomGroup);

      // 2. Build 3D Composite Chemical Bonds
      const bondGroup = new THREE.Group();
      const bondMat = new THREE.MeshStandardMaterial({
        color: isDark ? 0x94a3b8 : 0x475569,
        roughness: 0.4,
        metalness: 0.45,
      });

      const vStart = new THREE.Vector3();
      const vEnd = new THREE.Vector3();
      const vDir = new THREE.Vector3();
      const vUp = new THREE.Vector3(0, 1, 0);

      for (const bond of mol.bonds) {
        const a1 = mol.atoms.find((a) => a.id === bond.source);
        const a2 = mol.atoms.find((a) => a.id === bond.target);
        if (!a1 || !a2) continue;

        vStart.set(a1.x3d, a1.y3d, a1.z3d);
        vEnd.set(a2.x3d, a2.y3d, a2.z3d);
        vDir.subVectors(vEnd, vStart);
        const dist = vDir.length();
        if (dist < 0.05) continue;

        const mid = new THREE.Vector3().addVectors(vStart, vEnd).multiplyScalar(0.5);
        const quat = new THREE.Quaternion().setFromUnitVectors(
          vUp,
          vDir.clone().normalize()
        );

        const offsets =
          bond.order === 1
            ? [0]
            : bond.order === 2
              ? [-0.09, 0.09]
              : [-0.12, 0, 0.12];
        const cylRadius = bond.order === 1 ? 0.075 : 0.048;

        for (const off of offsets) {
          const bMesh = new THREE.Mesh(cylinderGeo, bondMat);
          bMesh.scale.set(cylRadius, dist, cylRadius);
          bMesh.position.copy(mid);
          bMesh.position.z += off;
          bMesh.quaternion.copy(quat);
          bondGroup.add(bMesh);
        }
      }
      contentGroup.add(bondGroup);

      // 3. Build 3D Molecular Orbital Representation
      if (mo) {
        const orbitalGroup = new THREE.Group();
        const posHex = isCb ? 0x0d9488 : 0x2563eb;
        const posEmissive = isCb ? 0x0f766e : 0x1d4ed8;
        const negHex = isCb ? 0xd97706 : 0xdc2626;
        const negEmissive = isCb ? 0xb45309 : 0xb91c1c;

        const posMat = new THREE.MeshPhysicalMaterial({
          color: posHex,
          emissive: posEmissive,
          emissiveIntensity: 0.22,
          roughness: 0.2,
          metalness: 0.05,
          clearcoat: isSoftwareGl ? 0 : 0.6,
          transparent: true,
          opacity: opacityVal,
          wireframe: mode === 'wireframe',
          depthWrite: false,
          side: THREE.DoubleSide,
        });

        const negMat = new THREE.MeshPhysicalMaterial({
          color: negHex,
          emissive: negEmissive,
          emissiveIntensity: 0.22,
          roughness: 0.2,
          metalness: 0.05,
          clearcoat: isSoftwareGl ? 0 : 0.6,
          transparent: true,
          opacity: opacityVal,
          wireframe: mode === 'wireframe',
          depthWrite: false,
          side: THREE.DoubleSide,
        });

        if (mode === 'merged' || mode === 'wireframe') {
          // Scalar-Field Marching Cubes / Tetrahedra Isosurface (t = 0 Separate AOs -> t = 1 Merged MO)
          const { posGeo, negGeo } = getOrCreateMorphKeyframeGeometries(
            mol,
            mo,
            scaleMult,
            tMorph
          );

          const posIsoMesh = new THREE.Mesh(posGeo, posMat);
          posIsoMesh.renderOrder = 2;
          posIsoMesh.visible =
            !!posGeo.attributes.position && posGeo.attributes.position.count > 0;
          orbitalGroup.add(posIsoMesh);
          activePosIsoMesh = posIsoMesh;

          const negIsoMesh = new THREE.Mesh(negGeo, negMat);
          negIsoMesh.renderOrder = 2;
          negIsoMesh.visible =
            !!negGeo.attributes.position && negGeo.attributes.position.count > 0;
          orbitalGroup.add(negIsoMesh);
          activeNegIsoMesh = negIsoMesh;
        } else if (mode === 'lobes') {
          const isDelocalized = mo.lcao.length > 10;
          const delocBoost = isDelocalized ? 1.28 : 1.0;

          for (const term of mo.lcao) {
            const atom = mol.atoms.find((a) => a.id === term.atomId);
            if (!atom || Math.abs(term.coeff) < 0.04) continue;

            const weight =
              Math.pow(Math.abs(term.coeff), isDelocalized ? 0.58 : 0.72) *
              scaleMult *
              delocBoost;
            const isPositiveCoeff = term.coeff >= 0;
            const primaryMat = isPositiveCoeff ? posMat : negMat;
            const oppositeMat = isPositiveCoeff ? negMat : posMat;

            const atomOrigin = new THREE.Vector3(atom.x3d, atom.y3d, atom.z3d);

            if (term.orbital === 's') {
              const sMesh = new THREE.Mesh(sphereGeo, primaryMat);
              sMesh.position.copy(atomOrigin);
              sMesh.scale.setScalar(0.95 * weight);
              sMesh.renderOrder = 2;
              orbitalGroup.add(sMesh);
            } else if (
              term.orbital === 'px' ||
              term.orbital === 'py' ||
              term.orbital === 'pz'
            ) {
              const axisVec =
                term.orbital === 'px'
                  ? new THREE.Vector3(1, 0, 0)
                  : term.orbital === 'py'
                    ? new THREE.Vector3(0, 1, 0)
                    : new THREE.Vector3(0, 0, 1);

              const posLobe = new THREE.Mesh(lobeGeo, primaryMat);
              posLobe.position.copy(atomOrigin);
              posLobe.quaternion.setFromUnitVectors(vUp, axisVec);
              posLobe.scale.setScalar(0.95 * weight);
              posLobe.renderOrder = 2;
              orbitalGroup.add(posLobe);

              const negLobe = new THREE.Mesh(lobeGeo, oppositeMat);
              negLobe.position.copy(atomOrigin);
              negLobe.quaternion.setFromUnitVectors(vUp, axisVec.clone().negate());
              negLobe.scale.setScalar(0.95 * weight);
              negLobe.renderOrder = 2;
              orbitalGroup.add(negLobe);
            } else {
              const addDirectedLobe = (
                dx: number,
                dy: number,
                dz: number,
                mat: THREE.Material,
                scaleFactor = 0.86
              ) => {
                const lobe = new THREE.Mesh(lobeGeo, mat);
                lobe.position.copy(atomOrigin);
                const dir = new THREE.Vector3(dx, dy, dz).normalize();
                lobe.quaternion.setFromUnitVectors(vUp, dir);
                lobe.scale.setScalar(scaleFactor * weight);
                lobe.renderOrder = 2;
                orbitalGroup.add(lobe);
              };

              if (term.orbital === 'dz2') {
                addDirectedLobe(0, 0, 1, primaryMat, 1.02);
                addDirectedLobe(0, 0, -1, primaryMat, 1.02);
                const torusBelly = new THREE.Mesh(
                  new THREE.TorusGeometry(0.55 * weight, 0.19 * weight, 14, 28),
                  oppositeMat
                );
                torusBelly.position.copy(atomOrigin);
                torusBelly.renderOrder = 2;
                orbitalGroup.add(torusBelly);
              } else if (term.orbital === 'dx2y2') {
                addDirectedLobe(1, 0, 0, primaryMat, 0.9);
                addDirectedLobe(-1, 0, 0, primaryMat, 0.9);
                addDirectedLobe(0, 1, 0, oppositeMat, 0.9);
                addDirectedLobe(0, -1, 0, oppositeMat, 0.9);
              } else if (term.orbital === 'dxy') {
                addDirectedLobe(1, 1, 0, primaryMat, 0.9);
                addDirectedLobe(-1, -1, 0, primaryMat, 0.9);
                addDirectedLobe(1, -1, 0, oppositeMat, 0.9);
                addDirectedLobe(-1, 1, 0, oppositeMat, 0.9);
              } else if (term.orbital === 'dxz') {
                addDirectedLobe(1, 0, 1, primaryMat, 0.9);
                addDirectedLobe(-1, 0, -1, primaryMat, 0.9);
                addDirectedLobe(1, 0, -1, oppositeMat, 0.9);
                addDirectedLobe(-1, 0, 1, oppositeMat, 0.9);
              } else if (term.orbital === 'dyz') {
                addDirectedLobe(0, 1, 1, primaryMat, 0.9);
                addDirectedLobe(0, -1, -1, primaryMat, 0.9);
                addDirectedLobe(0, 1, -1, oppositeMat, 0.9);
                addDirectedLobe(0, -1, 1, oppositeMat, 0.9);
              } else if (term.orbital === 'fz3') {
                addDirectedLobe(0, 0, 1, primaryMat, 1.04);
                addDirectedLobe(0, 0, -1, oppositeMat, 1.04);
                const collarGeo = new THREE.TorusGeometry(0.52 * weight, 0.15 * weight, 12, 26);
                const topCollar = new THREE.Mesh(collarGeo, oppositeMat);
                topCollar.position.copy(atomOrigin).add(new THREE.Vector3(0, 0, 0.34 * weight));
                topCollar.renderOrder = 2;
                orbitalGroup.add(topCollar);
                const botCollar = new THREE.Mesh(collarGeo, primaryMat);
                botCollar.position.copy(atomOrigin).add(new THREE.Vector3(0, 0, -0.34 * weight));
                botCollar.renderOrder = 2;
                orbitalGroup.add(botCollar);
              } else if (term.orbital === 'fxyz') {
                for (const sx of [-1, 1]) {
                  for (const sy of [-1, 1]) {
                    for (const sz of [-1, 1]) {
                      addDirectedLobe(
                        sx,
                        sy,
                        sz,
                        sx * sy * sz > 0 ? primaryMat : oppositeMat,
                        0.82
                      );
                    }
                  }
                }
              } else if (term.orbital === 'fz_x2y2') {
                for (const sz of [-1, 1]) {
                  addDirectedLobe(1, 0, sz, sz > 0 ? primaryMat : oppositeMat, 0.82);
                  addDirectedLobe(-1, 0, sz, sz > 0 ? primaryMat : oppositeMat, 0.82);
                  addDirectedLobe(0, 1, sz, sz > 0 ? oppositeMat : primaryMat, 0.82);
                  addDirectedLobe(0, -1, sz, sz > 0 ? oppositeMat : primaryMat, 0.82);
                }
              } else if (term.orbital === 'fx_x2_3y2' || term.orbital === 'fy_3x2_y2') {
                const offset = term.orbital === 'fy_3x2_y2' ? Math.PI / 6 : 0;
                for (let k = 0; k < 6; k++) {
                  const phi = offset + (k * Math.PI) / 3;
                  addDirectedLobe(
                    Math.cos(phi),
                    Math.sin(phi),
                    0,
                    k % 2 === 0 ? primaryMat : oppositeMat,
                    0.84
                  );
                }
              } else if (term.orbital === 'fxz2' || term.orbital === 'fyz2') {
                const isX = term.orbital === 'fxz2';
                for (let k = 0; k < 6; k++) {
                  const theta = (k * Math.PI) / 3;
                  const c = Math.cos(theta);
                  const s = Math.sin(theta);
                  addDirectedLobe(
                    isX ? c : 0,
                    isX ? 0 : c,
                    s,
                    k % 2 === 0 ? primaryMat : oppositeMat,
                    0.84
                  );
                }
              }
            }
          }
        } else {
          // Volumetric 3D Wavefunction Point-Cloud Isosurface Sampling
          const posPoints: number[] = [];
          const negPoints: number[] = [];
          const bound = Math.max(2.6, Math.min(15.0, rMol + 0.8));
          const step = Math.max(0.28, bound / 9.5);
          const threshold = (rMol > 5 ? 0.04 : 0.08) / Math.max(0.4, scaleMult);

          for (let x = -bound; x <= bound; x += step) {
            for (let y = -bound; y <= bound; y += step) {
              for (let z = -bound; z <= bound; z += step) {
                const psi = evaluateOrbitalWavefunction(mo, mol.atoms, x, y, z);
                if (psi > threshold) {
                  posPoints.push(x, y, z);
                } else if (psi < -threshold) {
                  negPoints.push(x, y, z);
                }
              }
            }
          }

          const addCloud = (pts: number[], colorHex: number) => {
            if (pts.length === 0) return;
            const count = Math.min(480, pts.length / 3);
            const inst = new THREE.InstancedMesh(
              new THREE.SphereGeometry(0.12 * scaleMult * Math.max(1, rMol / 4), 8, 8),
              new THREE.MeshStandardMaterial({
                color: colorHex,
                roughness: 0.3,
                transparent: true,
                opacity: opacityVal,
                depthWrite: false,
              }),
              count
            );
            const dummy = new THREE.Object3D();
            for (let i = 0; i < count; i++) {
              dummy.position.set(pts[i * 3], pts[i * 3 + 1], pts[i * 3 + 2]);
              dummy.updateMatrix();
              inst.setMatrixAt(i, dummy.matrix);
            }
            inst.computeBoundingSphere();
            orbitalGroup.add(inst);
          };

          addCloud(posPoints, posHex);
          addCloud(negPoints, negHex);
        }

        // 4. Nodal Plane Visualization (ψ = 0 surfaces)
        if (showNodes && mo.nodalPlanes > 0) {
          const planeSpan = Math.max(5.4, (rMol + 1.2) * 1.7);
          const nodePlaneGeo = new THREE.PlaneGeometry(planeSpan, Math.min(planeSpan, 7.2));
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
            opacity: 0.5,
          });

          const addNodePlane = (normal: [number, number, number]) => {
            const pGroup = new THREE.Group();
            const mesh = new THREE.Mesh(nodePlaneGeo, nodePlaneMat);
            mesh.renderOrder = 1;
            const edges = new THREE.LineSegments(nodeBorderGeo, nodeBorderMat);
            edges.renderOrder = 1;
            pGroup.add(mesh);
            pGroup.add(edges);
            pGroup.quaternion.setFromUnitVectors(
              new THREE.Vector3(0, 0, 1),
              new THREE.Vector3(...normal).normalize()
            );
            orbitalGroup.add(pGroup);
          };

          if (mo.interactionType === 'π' || mo.interactionType === 'π*') {
            addNodePlane([0, 0, 1]);
            if (mo.nodalPlanes >= 2) {
              addNodePlane([1, 0, 0]);
            }
          } else {
            const psiAlongX = Math.abs(evaluateOrbitalWavefunction(mo, mol.atoms, 0.75, 0, 0));
            const psiAlongY = Math.abs(evaluateOrbitalWavefunction(mo, mol.atoms, 0, 0.75, 0));
            const isYZeroNode = mo.nodalPlanes === 1 && psiAlongX > 0.02 && psiAlongY < 1e-4;
            addNodePlane(isYZeroNode ? [0, 1, 0] : [1, 0, 0]);
            if (mo.nodalPlanes >= 2) {
              addNodePlane([0, 1, 0]);
            }
          }
        }

        contentGroup.add(orbitalGroup);
      }

      // 5. Render 3D Point-Group Symmetry Elements (Cn Axes, σ Mirror Planes, Inversion Center i)
      if (showSym && pg) {
        const symGroup = new THREE.Group();
        const symAxisSpan = Math.max(5.8, (rMol + 1.2) * 2.1);
        const symPlaneSpan = Math.max(5.2, (rMol + 1.2) * 1.65);
        for (const el of pg.elements3D) {
          const isHighlighted = hiSymId === el.id;
          const baseColor = new THREE.Color(el.color);

          if (el.type === 'axis' && el.axis) {
            const dir = new THREE.Vector3(...el.axis).normalize();
            const axisLen = symAxisSpan;
            const rodGeo = new THREE.CylinderGeometry(
              isHighlighted ? 0.045 : 0.025,
              isHighlighted ? 0.045 : 0.025,
              axisLen,
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

            const tipGeo = new THREE.ConeGeometry(
              0.12,
              0.26,
              Math.max(3, el.fold || 4)
            );
            const tip = new THREE.Mesh(tipGeo, rodMat);
            tip.position.copy(dir.clone().multiplyScalar(axisLen / 2));
            tip.quaternion.setFromUnitVectors(vUp, dir);
            symGroup.add(tip);
          } else if (el.type === 'plane' && el.axis) {
            const normal = new THREE.Vector3(...el.axis).normalize();
            const planeGeo = new THREE.PlaneGeometry(symPlaneSpan, Math.min(symPlaneSpan, 7.0));
            const planeMat = new THREE.MeshBasicMaterial({
              color: baseColor,
              transparent: true,
              opacity: isHighlighted ? 0.28 : 0.09,
              side: THREE.DoubleSide,
              depthWrite: false,
            });
            const planeMesh = new THREE.Mesh(planeGeo, planeMat);
            planeMesh.renderOrder = 1;
            planeMesh.quaternion.setFromUnitVectors(
              new THREE.Vector3(0, 0, 1),
              normal
            );
            symGroup.add(planeMesh);
          } else if (el.type === 'center') {
            const invGeo = new THREE.OctahedronGeometry(isHighlighted ? 0.24 : 0.16, 0);
            const invMat = new THREE.MeshStandardMaterial({
              color: 0xf59e0b,
              emissive: 0xd97706,
              emissiveIntensity: 0.5,
              roughness: 0.2,
              metalness: 0.8,
            });
            const invMesh = new THREE.Mesh(invGeo, invMat);
            symGroup.add(invMesh);
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

    let lastWidth = width;
    let lastHeight = height;
    const onResize = () => {
      if (!container) return;
      const w = container.clientWidth || 640;
      const h = container.clientHeight || 420;
      if (Math.abs(w - lastWidth) < 2 && Math.abs(h - lastHeight) < 2) return;
      lastWidth = w;
      lastHeight = h;
      camera.aspect = w / h;
      camera.updateProjectionMatrix();
      renderer.setSize(w, h);
      controls.update();
      renderer.render(scene, camera);
    };
    window.addEventListener('resize', onResize);

    const onKeyDown = (e: KeyboardEvent) => {
      if (
        document.activeElement?.tagName === 'INPUT' ||
        document.activeElement?.tagName === 'TEXTAREA'
      ) {
        return;
      }
      if (e.key === 'h' || e.key === 'H') {
        setDockCollapsed((prev) => !prev);
      }
    };
    window.addEventListener('keydown', onKeyDown);

    return () => {
      window.removeEventListener('resize', onResize);
      window.removeEventListener('keydown', onKeyDown);
      controls.removeEventListener('change', onControlChange);
      controls.dispose();
      cancelAnimationFrame(frameId);
      clearTimeout(frameId);
      clearMorphKeyframeCache();
      sphereGeo.dispose();
      torusCollarGeo.dispose();
      cylinderGeo.dispose();
      lobeGeo.dispose();
      renderer.dispose();
    };
  }, []);

  // Trigger full scene rebuild only when structural or shading props change (morphProgress updates orbital geometry in-place via updateMorphFastRef)
  useEffect(() => {
    if (rebuildSceneRef.current) {
      rebuildSceneRef.current();
    }
  }, [
    molecule,
    selectedMo,
    pointGroup,
    highlightedSymmetryElementId,
    lobeScale,
    orbitalOpacity,
    renderMode,
    showSymmetryElements,
    showNodalPlanes,
    darkStage,
    colorblindPalette,
  ]);

  const morphPercent = Math.round(morphProgress * 100);

  return (
    <div className="bg-white border border-slate-200 rounded-lg flex flex-col overflow-hidden relative h-full">
      {/* Top Viewport Header Bar */}
      <div className="px-3.5 py-2 bg-slate-50 border-b border-slate-200 flex flex-wrap items-center justify-between gap-2 shrink-0">
        <div className="flex items-center gap-2 min-w-0">
          <h2 className="text-sm font-semibold text-slate-900 truncate">
            3D Orbital &amp; Symmetry Stage
          </h2>
          <span className="text-slate-300" aria-hidden="true">·</span>
          <span className="font-mono text-xs font-semibold text-blue-700 whitespace-nowrap">
            {selectedMo ? (
              <>
                <Latex math={moLabelToLatex(selectedMo.label)} /> (
                <Latex math={irrepToLatex(selectedMo.irrepSymbol)} />)
              </>
            ) : (
              ''
            )}
          </span>
        </div>

        {/* AO -> MO Morph Playback Controls & 3D Options Toggle */}
        <div className="flex items-center gap-2 shrink-0">
          <div
            className="flex items-center gap-1.5 bg-white border border-slate-200 rounded px-2 py-0.5"
            title="Animate or scrub continuous scalar-field Marching Cubes fusion from Separate AOs (0%) to Merged MO (100%)"
          >
            <button
              type="button"
              data-testid="morph-play-button"
              onClick={handleTogglePlayMorph}
              className={`px-2 py-0.5 text-[11px] font-medium rounded transition-colors flex items-center gap-1 whitespace-nowrap shrink-0 ${
                isMorphPlaying
                  ? 'bg-blue-600 text-white'
                  : 'bg-slate-900 text-white hover:bg-slate-800'
              }`}
            >
              {isMorphPlaying ? (
                <>
                  <Pause className="w-3 h-3" />
                  <span>Pause</span>
                </>
              ) : (
                <>
                  <Play className="w-3 h-3" />
                  <span>Morph MO</span>
                </>
              )}
            </button>

            <input
              id="header-morph-slider"
              data-testid="morph-slider"
              type="range"
              aria-label="AO to Merged MO Morph Percentage"
              min={0}
              max={1}
              step={0.01}
              value={morphProgress}
              onChange={(e) => handleScrubMorph(parseFloat(e.target.value))}
              className="w-20 accent-blue-600 cursor-pointer"
            />

            <span
              data-testid="morph-percent-label"
              className="font-mono text-[11px] font-semibold text-slate-700 tabular-nums w-8 text-right whitespace-nowrap"
            >
              {morphPercent}%
            </span>
          </div>

          <button
            type="button"
            onClick={() => setDockCollapsed((c) => !c)}
            aria-expanded={!dockCollapsed}
            className={`px-2.5 py-1 text-xs font-medium rounded border transition-colors flex items-center gap-1 whitespace-nowrap shrink-0 ${
              !dockCollapsed
                ? 'bg-slate-900 text-white border-slate-900'
                : 'bg-white text-slate-700 border-slate-200 hover:bg-slate-100'
            }`}
          >
            <span>3D Options</span>
            {!dockCollapsed ? (
              <ChevronUp className="w-3.5 h-3.5" />
            ) : (
              <ChevronDown className="w-3.5 h-3.5" />
            )}
          </button>
        </div>
      </div>

      {/* Main 3D WebGL Viewport Container */}
      <div className="relative w-full flex-1 min-h-[320px] select-none overflow-hidden">
        {/* Dedicated empty ref mount for Three.js WebGL canvas */}
        <div
          ref={canvasMountRef}
          className="w-full h-full block"
          data-testid="orbital-3d-canvas"
          data-orbital={selectedMo?.label || 'none'}
          data-point-group={pointGroup?.id || 'C2v'}
          data-render-mode={renderMode}
          data-morph-percent={morphPercent}
        />

        {/* Top-Left Minimal Spec Overlay (100% unobstructed center) */}
        <div
          className={`absolute top-3 left-3 pointer-events-none px-3 py-1.5 rounded border text-xs font-mono ${
            darkStage
              ? 'bg-slate-900/90 border-slate-700 text-slate-100'
              : 'bg-white/95 border-slate-200 text-slate-800'
          }`}
        >
          <div>
            <strong>
              <Latex math={formulaToLatex(molecule.formula)} />
            </strong>{' '}
            · Point Group{' '}
            <strong className="text-cyan-400">
              <Latex math={schoenfliesToLatex(molecule.schoenflies)} />
            </strong>
          </div>
          <div className="text-[11px] opacity-80">
            MO: <Latex math={moLabelToLatex(selectedMo?.label || '')} /> ·{' '}
            {selectedMo?.bondingType} · {selectedMo?.nodalPlanes} Nodal Plane
            {selectedMo?.nodalPlanes === 1 ? '' : 's'}
          </div>
          {(renderMode === 'merged' || renderMode === 'wireframe') && (
            <div className="text-[10.5px] text-sky-300 mt-0.5 tabular-nums">
              Surface:{' '}
              {morphPercent === 0
                ? 'Separate AOs (0%)'
                : morphPercent === 100
                  ? 'Merged MO Isosurface (100%)'
                  : `AO → MO Fusion (${morphPercent}%)`}
            </div>
          )}
        </div>

        {/* Collapsible Bottom Control Dock (Zero Viewport Occlusion) */}
        <div
          className="absolute left-1/2 -translate-x-1/2 bottom-3 z-10 w-[94%] max-w-[680px] flex flex-col items-center pointer-events-auto"
        >
          <button
            type="button"
            onClick={() => setDockCollapsed((c) => !c)}
            className="px-3 py-1 text-[11px] font-medium bg-slate-900 text-slate-200 border border-slate-700 rounded-t-md hover:bg-slate-800 transition-colors flex items-center gap-1 whitespace-nowrap"
          >
            {dockCollapsed ? (
              <>
                <ChevronUp className="w-3.5 h-3.5" />
                <span>Show 3D Controls (H)</span>
              </>
            ) : (
              <>
                <ChevronDown className="w-3.5 h-3.5" />
                <span>Hide 3D Controls (H)</span>
              </>
            )}
          </button>

          {!dockCollapsed && (
            <div className="w-full px-3.5 py-2.5 bg-slate-900/95 border border-slate-700 rounded-lg text-slate-100 shadow-lg flex flex-col gap-2.5 text-xs">
              {/* Row 1: AO -> MO Morph Scrubber + Lobe Scale Slider */}
              <div className="flex flex-wrap items-center justify-between gap-3 pb-2 border-b border-slate-800">
                <div className="flex items-center gap-2">
                  <button
                    type="button"
                    onClick={handleTogglePlayMorph}
                    className={`px-2.5 py-1 text-[11px] font-medium rounded transition-colors flex items-center gap-1 whitespace-nowrap shrink-0 ${
                      isMorphPlaying
                        ? 'bg-blue-600 text-white'
                        : 'bg-slate-800 text-slate-100 border border-slate-700 hover:bg-slate-700'
                    }`}
                  >
                    {isMorphPlaying ? (
                      <>
                        <Pause className="w-3 h-3" />
                        <span>Pause Morph</span>
                      </>
                    ) : (
                      <>
                        <Play className="w-3 h-3" />
                        <span>Play AO → MO</span>
                      </>
                    )}
                  </button>
                  <label
                    htmlFor="dock-morph-slider"
                    className="text-[11px] text-slate-300 font-mono tabular-nums whitespace-nowrap"
                  >
                    Separate AOs (0%) → Merged MO ({morphPercent}%):
                  </label>
                  <input
                    id="dock-morph-slider"
                    type="range"
                    min={0}
                    max={1}
                    step={0.01}
                    value={morphProgress}
                    onChange={(e) => handleScrubMorph(parseFloat(e.target.value))}
                    className="w-28 accent-blue-500 cursor-pointer"
                  />
                </div>

                {/* Lobe Scale Slider */}
                <div className="flex items-center gap-2">
                  <label htmlFor="lobe-scale-slider" className="text-[11px] text-slate-300 font-mono tabular-nums whitespace-nowrap">
                    Scale: {lobeScale.toFixed(1)}×
                  </label>
                  <input
                    id="lobe-scale-slider"
                    type="range"
                    min={0.5}
                    max={1.6}
                    step={0.1}
                    value={lobeScale}
                    onChange={(e) => setLobeScale(parseFloat(e.target.value))}
                    className="w-20 accent-blue-500 cursor-pointer"
                  />
                </div>
              </div>

              {/* Row 2: Render Mode Selector, Symmetry/Nodal Toggles & Camera Presets */}
              <div className="flex flex-wrap items-center justify-between gap-2">
                <div className="flex items-center gap-1 bg-slate-800 p-0.5 rounded" role="group" aria-label="Orbital Render Style">
                  {(
                    [
                      { id: 'merged', label: 'Merged MO' },
                      { id: 'lobes', label: 'LCAO Lobes' },
                      { id: 'cloud', label: 'ψ Grid Cloud' },
                      { id: 'wireframe', label: 'Mesh' },
                    ] as const
                  ).map((m) => (
                    <button
                      key={m.id}
                      type="button"
                      onClick={() => setRenderMode(m.id)}
                      className={`px-2 py-1 text-[11px] font-medium rounded transition-colors whitespace-nowrap shrink-0 ${
                        renderMode === m.id
                          ? 'bg-blue-600 text-white'
                          : 'text-slate-300 hover:text-white'
                      }`}
                    >
                      {m.label}
                    </button>
                  ))}
                </div>

                <div className="flex items-center gap-1.5">
                  <button
                    type="button"
                    onClick={() => setShowSymmetryElements((s) => !s)}
                    className={`px-2 py-1 text-[11px] font-medium rounded border transition-colors whitespace-nowrap shrink-0 ${
                      showSymmetryElements
                        ? 'bg-teal-600/90 text-white border-teal-500'
                        : 'bg-slate-800 text-slate-300 border-slate-700 hover:text-white'
                    }`}
                  >
                    Symmetry Axes/Planes
                  </button>

                  <button
                    type="button"
                    onClick={() => setShowNodalPlanes((n) => !n)}
                    className={`px-2 py-1 text-[11px] font-medium rounded border transition-colors whitespace-nowrap shrink-0 ${
                      showNodalPlanes
                        ? 'bg-blue-600/90 text-white border-blue-500'
                        : 'bg-slate-800 text-slate-300 border-slate-700 hover:text-white'
                    }`}
                  >
                    Nodal Planes
                  </button>
                </div>

                <div className="flex items-center gap-1">
                  <button
                    type="button"
                    onClick={() => resetCameraRef.current?.('iso')}
                    title="Reset 3D Isometric Camera"
                    className="px-2 py-1 text-[11px] font-medium bg-slate-800 text-slate-200 border border-slate-700 rounded hover:bg-slate-700 transition-colors whitespace-nowrap shrink-0"
                  >
                    Reset View
                  </button>
                  <button
                    type="button"
                    onClick={() => resetCameraRef.current?.('top')}
                    title="Top-Down Principal Axis View"
                    className="px-2 py-1 text-[11px] font-medium bg-slate-800 text-slate-200 border border-slate-700 rounded hover:bg-slate-700 transition-colors whitespace-nowrap shrink-0"
                  >
                    Top (z)
                  </button>
                  <button
                    type="button"
                    onClick={() => setColorblindPalette((p) => !p)}
                    title="Toggle Colorblind-Safe Phase Palette (Teal/Amber vs Blue/Rose)"
                    className="px-1.5 py-1 text-[11px] font-medium bg-slate-800 text-slate-200 border border-slate-700 rounded hover:bg-slate-700 transition-colors whitespace-nowrap shrink-0"
                  >
                    {colorblindPalette ? 'Teal/Amber' : 'Blue/Rose'}
                  </button>
                  <button
                    type="button"
                    onClick={() => setDarkStage((d) => !d)}
                    title="Toggle Light / Dark 3D Stage"
                    aria-label="Toggle Light or Dark 3D Stage"
                    className="p-1.5 bg-slate-800 text-slate-200 border border-slate-700 rounded hover:bg-slate-700 transition-colors shrink-0"
                  >
                    <Sun className="w-3.5 h-3.5" />
                  </button>
                </div>
              </div>
            </div>
          )}
        </div>
      </div>
    </div>
  );
};

