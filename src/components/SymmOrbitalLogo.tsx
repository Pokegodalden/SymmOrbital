import React from 'react';

interface SymmOrbitalLogoProps {
  className?: string;
  size?: number;
}

/**
 * SymmOrbital Brand Logo (Option 4: "The LCAO-MO Correlation Diamond")
 * Depicts left/right Atomic Orbital (AO) energy rungs splitting via dashed SALC lines
 * into a lower bonding HOMO rung (blue) and upper antibonding LUMO rung (amber),
 * framing a central dual-phase p-orbital wavefunction.
 */
export const SymmOrbitalLogo: React.FC<SymmOrbitalLogoProps> = ({
  className = 'w-8 h-8',
  size = 30,
}) => {
  return (
    <svg
      xmlns="http://www.w3.org/2000/svg"
      viewBox="0 0 120 120"
      width={size}
      height={size}
      style={{ width: size, height: size, minWidth: size, minHeight: size, display: 'inline-block' }}
      fill="none"
      className={className}
      role="img"
      aria-label="SymmOrbital Logo"
      data-testid="symmorbital-brand-logo"
    >
      <defs>
        <linearGradient id="symmLogoBg" x1="6" y1="6" x2="114" y2="114" gradientUnits="userSpaceOnUse">
          <stop offset="0%" stopColor="#1E293B" />
          <stop offset="100%" stopColor="#0F172A" />
        </linearGradient>
      </defs>
      {/* High-Contrast Background Badge (visible in both Light and Dark mode headers) */}
      <rect
        x="5"
        y="5"
        width="110"
        height="110"
        rx="24"
        fill="url(#symmLogoBg)"
        stroke="#3B82F6"
        strokeWidth="4"
      />
      {/* Left AO Energy Rung */}
      <line
        x1="15"
        y1="60"
        x2="34"
        y2="60"
        stroke="#CBD5E1"
        strokeWidth="5.5"
        strokeLinecap="round"
      />
      {/* Right AO Energy Rung */}
      <line
        x1="86"
        y1="60"
        x2="105"
        y2="60"
        stroke="#CBD5E1"
        strokeWidth="5.5"
        strokeLinecap="round"
      />
      {/* Top Antibonding LUMO (psi*) Rung */}
      <line
        x1="44"
        y1="24"
        x2="76"
        y2="24"
        stroke="#FBBF24"
        strokeWidth="6"
        strokeLinecap="round"
      />
      {/* Bottom Bonding HOMO (psi) Rung */}
      <line
        x1="44"
        y1="96"
        x2="76"
        y2="96"
        stroke="#38BDF8"
        strokeWidth="6"
        strokeLinecap="round"
      />
      {/* SALC Correlation Dashed Lines Forming a Diamond */}
      <line
        x1="34"
        y1="60"
        x2="44"
        y2="24"
        stroke="#94A3B8"
        strokeWidth="2.8"
        strokeDasharray="4 4"
      />
      <line
        x1="86"
        y1="60"
        x2="76"
        y2="24"
        stroke="#94A3B8"
        strokeWidth="2.8"
        strokeDasharray="4 4"
      />
      <line
        x1="34"
        y1="60"
        x2="44"
        y2="96"
        stroke="#94A3B8"
        strokeWidth="2.8"
        strokeDasharray="4 4"
      />
      <line
        x1="86"
        y1="60"
        x2="76"
        y2="96"
        stroke="#94A3B8"
        strokeWidth="2.8"
        strokeDasharray="4 4"
      />
      {/* Central Dual-Phase p-Orbital Inside the Correlation Diamond */}
      <path
        d="M60 60 C44 48, 44 32, 60 32 C76 32, 76 48, 60 60 Z"
        fill="#3B82F6"
        stroke="#93C5FD"
        strokeWidth="1.5"
      />
      <path
        d="M60 60 C44 72, 44 88, 60 88 C76 88, 76 72, 60 60 Z"
        fill="#F43F5E"
        stroke="#FDA4AF"
        strokeWidth="1.5"
      />
      <circle cx="60" cy="60" r="4.5" fill="#FFFFFF" />
    </svg>
  );
};
