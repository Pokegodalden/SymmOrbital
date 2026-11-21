/**
 * @license
 * SPDX-License-Identifier: Apache-2.0
 */

import React, { useEffect, useState } from 'react';
import { AppWorkbenchV3 } from './components/AppWorkbenchV3';
import { AppWorkbenchV2 } from './backup_v2/AppWorkbenchV2';
import { AppClassicV1 } from './backup_v1/AppClassicV1';
import { AppVersion } from './components/VersionSwitcherDropdown';

const VERSION_STORAGE_KEY = 'symmorbital_version';

export default function App() {
  // Active application release ('v3' default Computational Workbench, 'v2' Pro 3-Pane, 'v1' Classic Stacked)
  const [appVersion, setAppVersion] = useState<AppVersion>(() => {
    try {
      const saved = window.localStorage.getItem(VERSION_STORAGE_KEY);
      if (saved === 'v3' || saved === 'v2' || saved === 'v1') {
        return saved;
      }
      // One-time graceful default to v3.0 on upgrade
      return 'v3';
    } catch {
      return 'v3';
    }
  });

  useEffect(() => {
    try {
      window.localStorage.setItem(VERSION_STORAGE_KEY, appVersion);
    } catch {
      // ignore storage errors in restricted iframes
    }
  }, [appVersion]);

  // Global Dark / Light Mode state synchronized across v3.0, v2.0, and v1.0
  const [darkMode, setDarkMode] = useState<boolean>(() => {
    try {
      return window.localStorage.getItem('symmorbital_theme') === 'dark';
    } catch {
      return false;
    }
  });

  useEffect(() => {
    document.documentElement.classList.toggle('dark', darkMode);
    try {
      window.localStorage.setItem('symmorbital_theme', darkMode ? 'dark' : 'light');
    } catch {
      // ignore storage errors in restricted iframes
    }
  }, [darkMode]);

  const toggleDarkMode = () => setDarkMode((prev) => !prev);

  if (appVersion === 'v1') {
    return (
      <AppClassicV1
        currentVersion={appVersion}
        onSelectVersion={setAppVersion}
        onSwitchToProWorkbench={() => setAppVersion('v2')}
        darkMode={darkMode}
        onToggleDarkMode={toggleDarkMode}
      />
    );
  }

  if (appVersion === 'v2') {
    return (
      <AppWorkbenchV2
        currentVersion={appVersion}
        onSelectVersion={setAppVersion}
        darkMode={darkMode}
        onToggleDarkMode={toggleDarkMode}
      />
    );
  }

  return (
    <AppWorkbenchV3
      currentVersion={appVersion}
      onSelectVersion={setAppVersion}
      darkMode={darkMode}
      onToggleDarkMode={toggleDarkMode}
    />
  );
}
