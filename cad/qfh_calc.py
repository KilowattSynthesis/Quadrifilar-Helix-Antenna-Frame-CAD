"""Quadrifilar Helix/Helicoidal (QFH) Antenna geometry calculator.

Based on https://jcoppens.com/ant/qfh/calc.en.php

Translated from qfhcalc.js by John Coppens (jcoppens@usa.net)
Original: Copyright (C) 2000 John Coppens, GNU GPL v2+

Usage:

```python
results = calculate_qfh(
    freq=137.5,
    wdiam=7, wrad=15,
    ratio=0.44, turns=0.5,
    nrwavel=1
)
print_results(results)
```
"""

import math
from dataclasses import dataclass

# ---------------------------------------------------------------------------
# Correction factor tables
# ---------------------------------------------------------------------------

# Velocity factor correction based on conductor diameter (lambda correction).
_DELTA_L_TABLE = [
    1.045,
    1.053,
    1.060,
    1.064,
    1.068,
    1.070,
    1.070,
    1.071,
    1.071,
    1.070,
    1.070,
    1.070,
    1.070,
    1.069,
    1.069,
    1.068,
    1.067,
]

# End-effect correction based on conductor diameter.
_DELTA_F_TABLE = [
    1.013,
    1.014,
    1.015,
    1.016,
    1.017,
    1.018,
    1.020,
    1.022,
    1.025,
    1.027,
    1.030,
    1.033,
    1.036,
    1.041,
    1.044,
    1.049,
    1.054,
]


def _interpolate(table: list[float], diam: float) -> float:
    """Linearly interpolate a value from a lookup table.

    For example, look up using conductor diameter.
    """
    idx = int(diam)
    idx = max(0, min(idx, len(table) - 2))  # clamp to valid range
    return table[idx] + (table[idx + 1] - table[idx]) * (diam - idx)


def delta_l(diam: float) -> float:
    """Wavelength correction factor for conductor diameter."""
    return _interpolate(_DELTA_L_TABLE, diam)


def delta_f(diam: float) -> float:
    """End-effect correction factor for conductor diameter."""
    return _interpolate(_DELTA_F_TABLE, diam)


# ---------------------------------------------------------------------------
# Core geometry
# ---------------------------------------------------------------------------


def wavelength(freq_mhz: float) -> float:
    """Free-space wavelength in mm for a given frequency in MHz."""
    return 300_000.0 / freq_mhz


def compensated_wavelength(wavel: float, wdiam: float) -> float:
    """Wavelength corrected for conductor diameter effect.

    Conductor diameter is clamped to 15 mm (table limit).
    """
    wd_eff = min(wdiam, 15.0)
    return wavel * delta_l(wd_eff)


def bending_correction(wrad: float) -> float:
    """Length correction per bend to account for the rounded corner radius.

    Each loop has 4 bends; multiply by 4 for total correction.
    """
    return 2 * wrad - (math.pi * wrad / 2)


def optimal_diameter(wavelc: float) -> float:
    """Optimal conductor diameter for a given compensated wavelength (mm)."""
    return 0.0088 * wavelc


def _loop_geometry(
    total_compensated: float,
    ratio: float,
    turns: float,
    wrad: float,
    wdiam: float,
) -> dict:
    """Compute geometric dimensions for one loop (large or small).

    Parameters
    ----------
    total_compensated : float
        Total compensated wire length for this loop (mm).
    ratio : float
        Width-to-height (diameter/height) ratio.
    turns : float
        Number of turns (twist), e.g. 0.5 for 180°.
    wrad : float
        Bending radius (mm).
    wdiam : float
        Conductor outer diameter (mm).

    Returns
    -------
    dict with keys: rad, vert, height, idiam,
                    rad_comp, vert_comp

    """
    divisor = 1 + math.sqrt(1 / ratio**2 + (turns * math.pi) ** 2)
    rad = 0.5 * total_compensated / divisor  # horizontal half-separation
    vert = (total_compensated - 2 * rad) / 2  # vertical tube length
    height = rad / ratio  # antenna height

    return {
        "rad": rad,
        "vert": vert,
        "height": height,
        "idiam": rad - wdiam,  # internal diameter of cylinder
        "rad_comp": rad - 2 * wrad,  # horizontal separator (no bends)
        "vert_comp": vert - 2 * wrad,  # vertical separator (no bends)
    }


# ---------------------------------------------------------------------------
# Main calculation
# ---------------------------------------------------------------------------


@dataclass
class QfhInputSpec:
    """Input to the QFH calculator.

    Parameters
    ----------
    freq     : Design frequency (MHz). Default 137.5 MHz (NOAA weather sats).
    wdiam    : Conductor outer diameter (mm).
    wrad     : Bending radius — center of bend to conductor center (mm).
    ratio    : Diameter-to-height ratio (0.44 typical; 0.3-0.4 for better
               horizon coverage).
    turns    : Helix twist in fractions of a full turn (0.5 = 180°).
    nrwavel  : Loop circumference expressed in wavelengths (1, 1.5, or 2).

    """

    freq: float  # MHz
    wdiam: float  # mm  conductor outer diameter
    wrad: float  # mm  bending radius

    # Extremely-default settings:
    ratio: float = 0.44  # width/height ratio
    turns: float = 0.5  # number of turns (0.25 / 0.5 / 0.75 / 1.0)
    nrwavel: float = 1.0  # loop circumference in wavelengths (normally 1)

    def to_pretty_str(self, prefix: str = "") -> str:
        """Return a human-readable strrepresentation of the input parameters.

        Spans multiple lines for better readability.
        """
        return "\n".join(
            [
                prefix + f"Frequency: {self.freq} MHz",
                prefix + f"Conductor diameter: {self.wdiam} mm",
                prefix + f"Bending radius: {self.wrad} mm",
                prefix + f"Diameter/height ratio: {self.ratio}",
                prefix + f"Turns: {self.turns}",
                prefix + f"Loop length: {self.nrwavel} wavelengths",
            ]
        )


@dataclass
class LoopResult:
    """Result for one of two loops in a QFH antenna."""

    total: float  # Total wire length before compensation (mm)
    total_comp: float  # Total compensated wire length (mm)
    vert: float  # Vertical tube length (mm)
    vert_comp: float  # Compensated vertical tube (mm)
    height: float  # Antenna height H (mm)
    idiam: float  # Internal cylinder diameter Di (mm)
    rad: float  # Horizontal separator D (mm)
    rad_comp: float  # Compensated horizontal separator Dc (mm)


@dataclass
class QfhResult:
    """Overall result for the geometry of a QFH antenna."""

    # Wavelength data
    wavelength: float  # Free-space wavelength (mm)
    wavelength_comp: float  # Compensated wavelength (mm)
    bending_correction: float  # Per-bend correction (mm)
    optimal_diam: float  # Optimal conductor diameter (mm)

    # Loop dimensions
    large_loop: LoopResult
    small_loop: LoopResult


def _make_loop(
    total: float, total_c: float, geo: dict[str, float]
) -> LoopResult:
    return LoopResult(
        total=total,
        total_comp=total_c,
        vert=geo["vert"],
        vert_comp=geo["vert_comp"],
        height=geo["height"],
        idiam=geo["idiam"],
        rad=geo["rad"],
        rad_comp=geo["rad_comp"],
    )


def calculate_qfh(qfh_input_spec: QfhInputSpec) -> QfhResult:
    """Calculate all QFH antenna dimensions.

    Returns
    -------
    QfhResult dataclass with all intermediate and final dimensions.

    """
    wavel = wavelength(qfh_input_spec.freq)
    wavelc = compensated_wavelength(wavel, qfh_input_spec.wdiam)
    bcorr = bending_correction(qfh_input_spec.wrad)
    optd = optimal_diameter(wavelc)

    # Large loop: 1.026 x compensated wavelength
    total1 = wavelc * qfh_input_spec.nrwavel * 1.026
    total1c = total1 + 4 * bcorr
    geo1 = _loop_geometry(
        total1c,
        qfh_input_spec.ratio,
        qfh_input_spec.turns,
        qfh_input_spec.wrad,
        qfh_input_spec.wdiam,
    )

    # Small loop: 0.975 x compensated wavelength
    total2 = wavelc * qfh_input_spec.nrwavel * 0.975
    total2c = total2 + 4 * bcorr
    geo2 = _loop_geometry(
        total2c,
        qfh_input_spec.ratio,
        qfh_input_spec.turns,
        qfh_input_spec.wrad,
        qfh_input_spec.wdiam,
    )

    return QfhResult(
        wavelength=wavel,
        wavelength_comp=wavelc,
        bending_correction=bcorr,
        optimal_diam=optd,
        large_loop=_make_loop(total1, total1c, geo1),
        small_loop=_make_loop(total2, total2c, geo2),
    )


# ---------------------------------------------------------------------------
# Pretty-printer
# ---------------------------------------------------------------------------


def _fmt(value: float, decimals: int = 1) -> str:
    return f"{value:.{decimals}f}"


def print_results(input_spec: QfhInputSpec, r: QfhResult) -> None:
    """Print QFH calculation results in a readable format."""
    sep = "-" * 48

    print(sep)
    print("  QFH Antenna Calculator")
    print(sep)
    print("INPUT PARAMETERS:")
    print(input_spec.to_pretty_str(prefix="  "))

    print(sep)
    print("RESULTS:")
    print(f"  Wavelength             : {_fmt(r.wavelength)} mm")
    print(f"  Compensated wavelength : {_fmt(r.wavelength_comp)} mm")
    print(f"  Bending correction     : {_fmt(r.bending_correction)} mm")
    print(f"  Optimal conductor diam : {_fmt(r.optimal_diam)} mm")

    for label, loop in (("LARGE", r.large_loop), ("SMALL", r.small_loop)):
        print()
        print(f"  {label} LOOP")
        print(f"  {'Total length':<34}: {_fmt(loop.total)} mm")
        print(
            f"  {'Total compensated length':<34}: {_fmt(loop.total_comp)} mm"
        )
        print(f"  {'Vertical tube':<34}: {_fmt(loop.vert)} mm")
        print(
            f"  {'Compensated vertical tube':<34}: {_fmt(loop.vert_comp)} mm"
        )
        print(f"  {'Antenna height (H)':<34}: {_fmt(loop.height)} mm")
        print(f"  {'Internal diameter (Di)':<34}: {_fmt(loop.idiam)} mm")
        print(f"  {'Horizontal separator (D)':<34}: {_fmt(loop.rad)} mm")
        print(
            f"  {'Compensated horiz. sep. (Dc)':<34}: {_fmt(loop.rad_comp)} mm"
        )

    print(sep)


# ---------------------------------------------------------------------------
# Example
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # Default: 137.5 MHz NOAA satellite receive antenna.
    input_spec = QfhInputSpec(
        freq=137.5,
        wdiam=7.0,
        wrad=15.0,
        # Extremely-default settings:
        ratio=0.44,
        turns=0.5,
        nrwavel=1.0,
    )
    results = calculate_qfh(input_spec)
    print_results(input_spec, r=results)
