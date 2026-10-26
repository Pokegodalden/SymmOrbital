import {createRoot} from 'react-dom/client';
import App from './App.tsx';
import './index.css';

// Ensure third-party browser extensions that overwrite window.fetch do not throw
// when window.fetch is defined with only a getter in the host environment.
try {
  const desc = Object.getOwnPropertyDescriptor(window, 'fetch');
  if (desc && desc.get && !desc.set && desc.configurable !== false) {
    let currentFetch = window.fetch.bind(window);
    Object.defineProperty(window, 'fetch', {
      configurable: true,
      enumerable: desc.enumerable ?? true,
      get() {
        return currentFetch;
      },
      set(fn) {
        if (typeof fn === 'function') {
          currentFetch = fn;
        }
      },
    });
  }
} catch {
  // Ignore if descriptor is non-configurable
}

createRoot(document.getElementById('root')!).render(<App />);
