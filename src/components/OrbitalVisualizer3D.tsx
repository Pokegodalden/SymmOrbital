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
  Sun,
} from 'lucide-react';

interface OrbitalVisualizer3DProps {
  molecule: MoleculeEntry;
  selectedMo: MolecularOrbital;
  pointGroup: PointGroupData;
  highlightedSymmetryElementId?: string | null;
}

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

  const [lobeScale, setLobeScale] = useState<number>(1.0);
  const [orbitalOpacity, setOrbitalOpacity] = useState<number>(0.68);
  const [renderMode, setRenderMode] = useState<'lobes' | 'cloud' | 'wireframe'>('lobes');
  const [showSymmetryElements, setShowSymmetryElements] = useState<boolean>(true);
  const [showNodalPlanes, setShowNodalPlanes] = useState<boolean>(true);
  const [darkStage, setDarkStage] = useState<boolean>(true);
  const [colorblindPalette, setColorblindPalette] = useState<boolean>(false);
  const [dockCollapsed, setDockCollapsed] = useState<boolean>(true);

  // Store mutable props in ref so WebGL lifecycle never tears down on prop change
  const stateRef = useRef({
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
  });

  stateRef.current = {
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

    const camera = new THREE.PerspectiveCamera(42, width / height, 0.1, 100);
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

    // Deformed teardrop lobe geometry for p/d orbital lobes
    const lobeGeo = new THREE.SphereGeometry(0.65, 26, 20);
    const posAttr = lobeGeo.attributes.position;
    const vTemp = new THREE.Vector3();
    for (let i = 0; i < posAttr.count; i++) {
      vTemp.fromBufferAttribute(posAttr, i);
      // Shift center along +Y so lobe emanates from nucleus at y=0 to y=1.3
      const ny = (vTemp.y + 0.65) / 1.3; // 0 at base, 1 at tip
      const pinch = Math.pow(ny, 0.55) * 1.15;
      vTemp.x *= pinch;
      vTemp.z *= pinch;
      vTemp.y = ny * 1.35;
      posAttr.setXYZ(i, vTemp.x, vTemp.y, vTemp.z);
    }
    lobeGeo.computeVertexNormals();

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

    resetCameraRef.current = (preset = 'iso') => {
      if (preset === 'top') {
        camera.position.set(0, 6.5, 0.01);
      } else if (preset === 'front') {
        camera.position.set(0, 0.1, 6.5);
      } else if (preset === 'side') {
        camera.position.set(6.5, 0.1, 0.01);
      } else {
        camera.position.set(2.2, 2.0, 5.8);
      }
      controls.target.set(0, 0, 0);
      controls.update();
      requestRender();
    };

    const buildMolecularAndOrbitalScene = () => {
      // Clear previous children & dispose dynamic materials
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
        showSymmetryElements: showSym,
        showNodalPlanes: showNodes,
        darkStage: isDark,
        colorblindPalette: isCb,
      } = stateRef.current;

      const bgHex = isDark ? 0x0b0f19 : 0xf1f5f9;
      scene.background = new THREE.Color(bgHex);
      scene.fog = new THREE.FogExp2(bgHex, 0.032);

      // 1. Build 3D Composite Atoms
      const atomGroup = new THREE.Group();
      for (const atom of mol.atoms) {
        const props = ELEMENT_PROPERTIES[atom.element] || ELEMENT_PROPERTIES.C;
        const radius = Math.max(0.22, Math.min(0.46, props.covalentRadius * 0.45));
        const singleAtomGroup = new THREE.Group();
        singleAtomGroup.position.set(atom.x3d, atom.y3d, atom.z3d);

        // Primary CPK sphere
        const matMain = new THREE.MeshStandardMaterial({
          color: new THREE.Color(props.cpkColor),
          roughness: 0.25,
          metalness: 0.12,
        });
        const sphereMesh = new THREE.Mesh(sphereGeo, matMain);
        sphereMesh.scale.setScalar(radius);
        singleAtomGroup.add(sphereMesh);

        // Secondary metallic equatorial collar ring for visual craft
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

      // 3. Build 3D Molecular Orbital Lobes
      if (mo) {
        const orbitalGroup = new THREE.Group();
        const posHex = isCb ? 0x0d9488 : 0x2563eb;
        const posEmissive = isCb ? 0x0f766e : 0x1d4ed8;
        const negHex = isCb ? 0xd97706 : 0xdc2626;
        const negEmissive = isCb ? 0xb45309 : 0xb91c1c;

        const posMat = new THREE.MeshPhysicalMaterial({
          color: posHex, // Positive wavefunction phase (+ψ)
          emissive: posEmissive,
          emissiveIntensity: 0.18,
          roughness: 0.18,
          metalness: 0.05,
          clearcoat: isSoftwareGl ? 0 : 0.6,
          transparent: true,
          opacity: opacityVal,
          wireframe: mode === 'wireframe',
          depthWrite: false,
          side: THREE.DoubleSide,
        });

        const negMat = new THREE.MeshPhysicalMaterial({
          color: negHex, // Negative wavefunction phase (−ψ)
          emissive: negEmissive,
          emissiveIntensity: 0.18,
          roughness: 0.18,
          metalness: 0.05,
          clearcoat: isSoftwareGl ? 0 : 0.6,
          transparent: true,
          opacity: opacityVal,
          wireframe: mode === 'wireframe',
          depthWrite: false,
          side: THREE.DoubleSide,
        });

        if (mode === 'lobes' || mode === 'wireframe') {
          for (const term of mo.lcao) {
            const atom = mol.atoms.find((a) => a.id === term.atomId);
            if (!atom || Math.abs(term.coeff) < 0.08) continue;

            const weight = Math.pow(Math.abs(term.coeff), 0.72) * scaleMult;
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

              // Positive lobe along +axisVec
              const posLobe = new THREE.Mesh(lobeGeo, primaryMat);
              posLobe.position.copy(atomOrigin);
              posLobe.quaternion.setFromUnitVectors(vUp, axisVec);
              posLobe.scale.setScalar(0.95 * weight);
              posLobe.renderOrder = 2;
              orbitalGroup.add(posLobe);

              // Negative lobe along -axisVec
              const negLobe = new THREE.Mesh(lobeGeo, oppositeMat);
              negLobe.position.copy(atomOrigin);
              negLobe.quaternion.setFromUnitVectors(vUp, axisVec.clone().negate());
              negLobe.scale.setScalar(0.95 * weight);
              negLobe.renderOrder = 2;
              orbitalGroup.add(negLobe);
            } else if (term.orbital === 'dz2') {
              // Axial lobes along ±Y and equatorial torus in XZ
              const posLobe = new THREE.Mesh(lobeGeo, primaryMat);
              posLobe.position.copy(atomOrigin);
              posLobe.scale.setScalar(1.05 * weight);
              orbitalGroup.add(posLobe);

              const posLobe2 = new THREE.Mesh(lobeGeo, primaryMat);
              posLobe2.position.copy(atomOrigin);
              posLobe2.rotation.x = Math.PI;
              posLobe2.scale.setScalar(1.05 * weight);
              orbitalGroup.add(posLobe2);

              const torusBelly = new THREE.Mesh(
                new THREE.TorusGeometry(0.55 * weight, 0.2 * weight, 14, 28),
                oppositeMat
              );
              torusBelly.position.copy(atomOrigin);
              torusBelly.rotation.x = Math.PI / 2;
              orbitalGroup.add(torusBelly);
            }
          }
        } else {
          // Volumetric 3D Wavefunction Point-Cloud Isosurface Sampling
          const posPoints: number[] = [];
          const negPoints: number[] = [];
          const step = 0.28;
          const bound = 2.6;
          const threshold = 0.08 / Math.max(0.4, scaleMult);

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
            const count = Math.min(400, pts.length / 3);
            const inst = new THREE.InstancedMesh(
              new THREE.SphereGeometry(0.12 * scaleMult, 8, 8),
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
          const nodePlaneGeo = new THREE.PlaneGeometry(4.2, 4.2);
          const nodePlaneMat = new THREE.MeshBasicMaterial({
            color: isDark ? 0x94a3b8 : 0x475569,
            transparent: true,
            opacity: 0.16,
            side: THREE.DoubleSide,
            depthWrite: false,
          });
          const nodeBorderGeo = new THREE.EdgesGeometry(nodePlaneGeo);
          const nodeBorderMat = new THREE.LineBasicMaterial({
            color: isDark ? 0x38bdf8 : 0x0284c7,
            transparent: true,
            opacity: 0.55,
          });

          const addNodePlane = (normal: [number, number, number]) => {
            const pGroup = new THREE.Group();
            const mesh = new THREE.Mesh(nodePlaneGeo, nodePlaneMat);
            const edges = new THREE.LineSegments(nodeBorderGeo, nodeBorderMat);
            pGroup.add(mesh);
            pGroup.add(edges);
            pGroup.quaternion.setFromUnitVectors(
              new THREE.Vector3(0, 0, 1),
              new THREE.Vector3(...normal).normalize()
            );
            orbitalGroup.add(pGroup);
          };

          if (mo.interactionType === 'π' || mo.interactionType === 'π*') {
            // Horizontal π nodal plane (z = 0)
            addNodePlane([0, 0, 1]);
            if (mo.nodalPlanes >= 2) {
              addNodePlane([1, 0, 0]);
            }
          } else {
            addNodePlane([1, 0, 0]);
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
        for (const el of pg.elements3D) {
          const isHighlighted = hiSymId === el.id;
          const baseColor = new THREE.Color(el.color);

          if (el.type === 'axis' && el.axis) {
            const dir = new THREE.Vector3(...el.axis).normalize();
            const axisLen = 4.8;
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

            // End-cap n-fold polygon marker
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
            const planeGeo = new THREE.PlaneGeometry(3.8, 3.8);
            const planeMat = new THREE.MeshBasicMaterial({
              color: baseColor,
              transparent: true,
              opacity: isHighlighted ? 0.28 : 0.1,
              side: THREE.DoubleSide,
              depthWrite: false,
            });
            const planeMesh = new THREE.Mesh(planeGeo, planeMat);
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

      requestRender();
    };

    rebuildSceneRef.current = buildMolecularAndOrbitalScene;
    buildMolecularAndOrbitalScene();

    if (!isSoftwareGl) {
      renderer.compile(scene, camera);
    }

    const onControlChange = () => requestRender();
    controls.addEventListener('change', onControlChange);

    const onResize = () => {
      if (!container) return;
      const w = container.clientWidth || 640;
      const h = container.clientHeight || 420;
      camera.aspect = w / h;
      camera.updateProjectionMatrix();
      renderer.setSize(w, h);
      requestRender();
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
      sphereGeo.dispose();
      torusCollarGeo.dispose();
      cylinderGeo.dispose();
      lobeGeo.dispose();
      renderer.dispose();
    };
  }, []);

  // Trigger scene rebuild whenever props or orbital controls change
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

  return (
    <div className="bg-white border border-slate-200 rounded-lg flex flex-col overflow-hidden relative h-full">
      {/* Top Viewport Header Bar */}
      <div className="px-3.5 py-2.5 bg-slate-50 border-b border-slate-200 flex items-center justify-between gap-2 shrink-0">
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

        {/* Collapsible 3D Options Button in Header (Collapsed by Default) */}
        <div className="flex items-center gap-2 shrink-0">
          <div className="hidden sm:flex items-center gap-1.5 text-xs font-mono">
            <span className={colorblindPalette ? 'text-teal-700 font-semibold' : 'text-blue-700 font-semibold'}>
              <Latex math="+\psi" />
            </span>
            <span className="text-slate-300" aria-hidden="true">/</span>
            <span className={colorblindPalette ? 'text-amber-700 font-semibold' : 'text-rose-700 font-semibold'}>
              <Latex math="-\psi" />
            </span>
          </div>

          <button
            type="button"
            onClick={() => setDockCollapsed((c) => !c)}
            aria-expanded={!dockCollapsed}
            className={`px-2.5 py-1 text-xs font-medium rounded border transition-colors flex items-center gap-1 whitespace-nowrap ${
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
            <div className="w-full px-3.5 py-2.5 bg-slate-900/95 border border-slate-700 rounded-lg text-slate-100 shadow-lg flex flex-wrap items-center justify-between gap-3 text-xs">
              {/* Lobe Scale Slider */}
              <div className="flex items-center gap-2">
                <label htmlFor="lobe-scale-slider" className="text-[11px] text-slate-300 font-mono whitespace-nowrap">
                  Lobe Scale: {lobeScale.toFixed(1)}×
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

              {/* Render Mode Selector */}
              <div className="flex items-center gap-1 bg-slate-800 p-0.5 rounded" role="group" aria-label="Orbital Render Style">
                {(
                  [
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

              {/* Symmetry & Nodal Toggles */}
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

              {/* Camera Presets & Stage Theme */}
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
          )}
        </div>
      </div>
    </div>
  );
};
