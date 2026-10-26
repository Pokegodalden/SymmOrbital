# SymmOrbital V1 Classic Layout Backup

This folder preserves the complete **V1 Classic Stacked Layout** (`AppClassicV1.tsx`) prior to the **V2 Professional 3-Pane Synchronized Workbench** redesign.

## How to Switch or Revert

1. **Instant In-App Reversibility (Zero Code Changes)**:
   - Click the **"Classic Layout (v1)"** / **"Pro Workbench (v2)"** toggle button in the top navigation bar at any time. Your preference is saved in `localStorage` under `symmorbital_layout_mode`.

2. **Source-Level Reversibility**:
   - In `/src/App.tsx`, change the default `layoutMode` initial state from `'pro'` to `'classic'`, or re-export `AppClassicV1` directly as default:
     ```tsx
     export { AppClassicV1 as default } from './backup_v1/AppClassicV1';
     ```
