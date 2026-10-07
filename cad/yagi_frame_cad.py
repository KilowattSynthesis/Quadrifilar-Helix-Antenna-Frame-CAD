"""436 MHz handheld/mast 2-element beam (reflector + folded-dipole driver).

Frame for 10 mm copper tape. Layout, in the assembled pose:

* +Z is the firing direction (boom points up, mast at the bottom).
* Elements run along X and lie in the XZ plane, in FRONT of the boom
  (0 <= y <= bar_depth).
* The boom sits BEHIND the elements (-boom_thickness <= y <= 0) so it never
  interrupts the tape.
* Reflector tape: on the +Z (forward-facing) face of the reflector bar.
* Driven tape: wrapped around the OUTER perimeter of the folded-dipole frame
  (top face of top bar, both ends, bottom face of bottom bar). The bottom run
  is split at x = 0 by the feed gap.

Each element is printed as two IDENTICAL halves that half-lap over the boom
and are clamped to it with two M4 bolts (use nylon, especially at the feed).
Print two of each half.

Starting dimensions only: tape on plastic, flat conductors and your print
material all shift resonance. Model it (e.g. 4nec2) and/or verify with a VNA.
The reflector has trim margin and slotted mounting holes so you can tune
length and spacing after printing.
"""

from dataclasses import dataclass
from pathlib import Path

import build123d as bd
from build123d_ease import show
from loguru import logger

SPEED_OF_LIGHT_MM_PER_S = 299_792_458_000.0
MIN_REFLECTOR_TO_SOCKET_CLEARANCE = 10.0
MIN_SOCKET_WALL = 2.0


@dataclass
class Spec:
    """Specification for the 436 MHz 2-element beam frame."""

    frequency_mhz: float = 436.0

    # --- Electrical dimensions (copper tape) ---
    # Reflector tape length, tip to tip (~0.51 λ).
    reflector_length: float = 352.0
    # Folded dipole outer length; the tape wraps the frame ends (~0.47 λ).
    driven_length: float = 326.0
    # Folded dipole outer height = spacing between its two conductors.
    driven_height: float = 40.0
    # Reflector centreline to folded-dipole centreline (~0.125 λ).
    # Closer spacing -> lower feed impedance.
    # ~0.12 λ is a starting guess for 50 ohm.
    element_spacing: float = 86.0
    # Gap in the bottom conductor where the coax is soldered.
    feed_gap: float = 4.0

    # --- Tape ---
    tape_width: float = 10.0
    tape_channel_clearance: float = 0.6
    tape_channel_depth: float = 0.4
    # Extra bar length beyond the nominal reflector length, per end, for
    # trimming.
    reflector_trim_margin: float = 8.0
    # Small cross-groove marking the nominal reflector tape end.
    trim_mark_width: float = 0.8

    # --- Element bars ---
    bar_depth: float = 14.0  # Y; the tape-carrying faces are this wide.
    bar_height: float = 12.0  # Z (bar thickness in the element plane).
    frame_corner_radius: float = (
        2.0  # Outer corner radius so the tape wraps smoothly.
    )

    # --- Centre joint (half-lap over the boom) ---
    lap_length: float = 30.0
    bolt_hole_diameter: float = 4.4  # M4 clearance.
    bolt_x_offset: float = 7.0
    reflector_slot_travel: float = (
        16.0  # +/- 8 mm of reflector spacing adjustment.
    )

    # --- Boom ---
    boom_width: float = 24.0  # X
    boom_thickness: float = 20.0  # Y
    # Top of mast socket to reflector centreline.
    reflector_standoff: float = 45.0
    boom_top_margin: float = 6.0
    coax_hole_diameter: float = 6.5
    zip_slot_width: float = 3.5
    zip_slot_height: float = 6.0

    # --- Mast socket ---
    mast_pipe_od: float | None = 42.2  # 1-1/4" PVC pipe — OD = 1.660 in
    mast_clearance: float = 0.6
    socket_wall: float = 4.0
    socket_depth: float = 50.0
    socket_cap: float = 6.0
    mast_bolt_diameter: float = 5.5  # M5 cross-bolt.
    gusset_height: float = 25.0

    def __post_init__(self) -> None:
        """Post initialization checks."""
        assert self.frequency_mhz > 0
        assert (
            self.tape_width + self.tape_channel_clearance < self.bar_depth
        ), "Tape channel must fit on the bar face"
        assert self.driven_height > 2 * self.bar_height + 4, (
            "Folded dipole frame needs an inner opening"
        )
        assert self.driven_length < self.reflector_length, (
            "Reflector should be longer than the driven element"
        )
        assert self.lap_length > self.boom_width, (
            "Lap joint should span the boom"
        )
        assert (
            self.bolt_x_offset + self.bolt_hole_diameter / 2
            < min(self.lap_length, self.boom_width) / 2 - 1.5
        ), "Bolt holes too close to the edge of the lap/boom"
        assert (
            self.bolt_x_offset
            > self.bolt_hole_diameter / 2 + self.feed_gap / 2
        ), "Bolt holes collide at the centre"
        assert self.feed_gap < self.lap_length
        reflector_top = self.reflector_standoff + self.bar_height / 2
        driven_bottom = (
            self.reflector_standoff
            + self.element_spacing
            - self.driven_height / 2
        )
        assert driven_bottom - reflector_top > self.coax_hole_diameter + 12, (
            "Element spacing too small for the coax hole"
        )
        assert (
            self.reflector_standoff
            - self.bar_height / 2
            - self.reflector_slot_travel / 2
            > MIN_REFLECTOR_TO_SOCKET_CLEARANCE
        ), "Reflector too close to the mast socket"
        if self.mast_pipe_od is not None:
            assert self.mast_pipe_od > 0
            assert self.socket_wall >= MIN_SOCKET_WALL

    @property
    def wavelength(self) -> float:
        """Free-space wavelength in mm."""
        return SPEED_OF_LIGHT_MM_PER_S / (self.frequency_mhz * 1e6)

    @property
    def reflector_z(self) -> float:
        """Z of the reflector centreline."""
        return self.reflector_standoff

    @property
    def driven_z(self) -> float:
        """Z of the folded-dipole centreline."""
        return self.reflector_standoff + self.element_spacing

    @property
    def boom_top_z(self) -> float:
        """Z of the top of the boom."""
        return self.driven_z + self.driven_height / 2 + self.boom_top_margin

    @property
    def channel_width(self) -> float:
        """Width of the tape-locating channel."""
        return self.tape_width + self.tape_channel_clearance


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _y_prism(
    sketch: bd.Sketch, y0: float, depth: float
) -> bd.Part | bd.Compound:
    """Extrude an XY-sketch (X = X, sketch Y = world Z) along +Y from y0."""
    return bd.Part(None) + bd.Pos(0, y0, 0) * (
        bd.Rot(X=-90) * bd.extrude(sketch, amount=depth)
    )


def _y_hole(
    spec: Spec, x: float, z: float, diameter: float
) -> bd.Part | bd.Compound:
    """Through-hole along Y covering the boom and element bars."""
    length = spec.boom_thickness + spec.bar_depth + 4
    return bd.Part(None) + bd.Pos(
        x, (spec.bar_depth - spec.boom_thickness) / 2, z
    ) * (bd.Rot(X=90) * bd.Cylinder(radius=diameter / 2, height=length))


def _right_half_region(spec: Spec) -> bd.Part | bd.Compound:
    """Region kept by the right-hand half of a half-lapped element.

    Full section for x >= lap/2; front half (y >= depth/2) across the lap.
    Rotating this 180 deg about Z through y = depth/2 gives the left half, so
    both halves are the same printed part.
    """
    big = 2000.0
    outer = bd.Part(None) + bd.Pos(
        spec.lap_length / 2 + big / 2, spec.bar_depth / 2, 0
    ) * bd.Box(big, spec.bar_depth + 2, big)
    lap = bd.Pos(0, spec.bar_depth * 3 / 4 + 0.5, 0) * bd.Box(
        spec.lap_length, spec.bar_depth / 2 + 1, big
    )
    return outer + lap


def _split(
    spec: Spec, element: bd.Part | bd.Compound
) -> tuple[bd.Part | bd.Compound, bd.Part | bd.Compound]:
    """Split an element into (right, left) half-lap halves."""
    region = _right_half_region(spec)
    return element & region, element - region


def _for_printing(part: bd.Part | bd.Compound) -> bd.Part | bd.Compound:
    """Lay an element half on its front face (+Y) and drop it onto the bed."""
    p = bd.Rot(X=-90) * part  # +Y face -> -Z (bed)
    bb = p.bounding_box()
    return bd.Part(None) + (
        bd.Pos(-bb.center().X, -bb.center().Y, -bb.min.Z) * p
    )


# ---------------------------------------------------------------------------
# Elements
# ---------------------------------------------------------------------------


def _element_bolt_holes(spec: Spec, z: float) -> bd.Part | bd.Compound:
    holes = bd.Part(None)
    for x in (-spec.bolt_x_offset, spec.bolt_x_offset):
        holes += _y_hole(spec, x, z, spec.bolt_hole_diameter)
    return holes


def reflector_full(spec: Spec) -> bd.Part | bd.Compound:
    """Reflector bar in its assembled position (before splitting)."""
    z = spec.reflector_z
    bar_len = spec.reflector_length + 2 * spec.reflector_trim_margin

    p = bd.Part(None)
    p += bd.Pos(0, spec.bar_depth / 2, z) * bd.Box(
        bar_len, spec.bar_depth, spec.bar_height
    )

    # Tape channel along the forward (+Z) face, full bar length.
    top = z + spec.bar_height / 2
    p -= bd.Pos(0, spec.bar_depth / 2, top) * bd.Box(
        bar_len + 2, spec.channel_width, 2 * spec.tape_channel_depth
    )

    # Trim marks at the nominal tape ends.
    for x in (-spec.reflector_length / 2, spec.reflector_length / 2):
        p -= bd.Pos(x, spec.bar_depth / 2, top) * bd.Box(
            spec.trim_mark_width,
            spec.bar_depth + 2,
            2 * (spec.tape_channel_depth + 0.4),
        )

    p -= _element_bolt_holes(spec, z)
    return p


def driven_full(spec: Spec) -> bd.Part | bd.Compound:
    """Folded-dipole frame in its assembled position (before splitting)."""
    z = spec.driven_z
    outer_sk = bd.RectangleRounded(
        spec.driven_length, spec.driven_height, spec.frame_corner_radius
    )
    inner_sk = bd.Rectangle(
        spec.driven_length - 2 * spec.bar_height,
        spec.driven_height - 2 * spec.bar_height,
    )

    frame_sk = bd.Sketch((outer_sk - inner_sk).faces())
    p = bd.Part(None) + bd.Pos(0, 0, z) * _y_prism(frame_sk, 0, spec.bar_depth)

    # Tape channel: a band around the outer perimeter.
    inset_sk = bd.offset(outer_sk, amount=-spec.tape_channel_depth)
    assert isinstance(inset_sk, bd.Sketch)  # Type checking.
    band = bd.Part(None) + bd.Pos(0, 0, z) * (
        _y_prism(outer_sk, -1, spec.bar_depth + 2)
        - _y_prism(inset_sk, -2, spec.bar_depth + 4)
    )
    band &= bd.Pos(0, spec.bar_depth / 2, z) * bd.Box(
        spec.driven_length + 10, spec.channel_width, spec.driven_height + 10
    )
    # Leave a raised divider at the feed gap so the two tape ends can't touch.
    band -= bd.Pos(0, spec.bar_depth / 2, z - spec.driven_height / 2) * bd.Box(
        spec.feed_gap, spec.bar_depth + 2, 4 * spec.tape_channel_depth
    )
    p -= band

    top_bar_z = z + spec.driven_height / 2 - spec.bar_height / 2
    bottom_bar_z = z - spec.driven_height / 2 + spec.bar_height / 2
    p -= _element_bolt_holes(spec, top_bar_z)
    p -= _element_bolt_holes(spec, bottom_bar_z)
    return p


# ---------------------------------------------------------------------------
# Boom + mast socket
# ---------------------------------------------------------------------------


def boom(spec: Spec) -> bd.Part | bd.Compound:
    """Boom with integrated mast socket, in its assembled position."""
    yc = -spec.boom_thickness / 2
    p = bd.Part(None)
    p += bd.Pos(0, yc, spec.boom_top_z / 2) * bd.Box(
        spec.boom_width, spec.boom_thickness, spec.boom_top_z
    )

    if spec.mast_pipe_od is not None:
        bore_r = (spec.mast_pipe_od + spec.mast_clearance) / 2
        outer_r = bore_r + spec.socket_wall
        socket_h = spec.socket_depth + spec.socket_cap
        p += bd.Pos(0, yc, -socket_h / 2) * bd.Cylinder(outer_r, socket_h)
        # Gusset blending the socket into the boom.
        gusset = bd.Pos(0, yc, spec.gusset_height / 2) * bd.Cone(
            outer_r, spec.boom_width / 2, spec.gusset_height
        )
        # Keep the gusset behind the element plane so it can't touch the
        # reflector.
        gusset &= bd.Pos(0, -500, 0) * bd.Box(1000, 1000, 1000)
        p += gusset

        # Bore, with a 45 deg conical roof so it prints without bridging.
        p -= bd.Pos(
            0, yc, -socket_h + spec.socket_depth / 2 - 0.01
        ) * bd.Cylinder(bore_r, spec.socket_depth)
        cone_h = min(bore_r, spec.socket_cap + spec.socket_depth / 2)
        p -= bd.Pos(
            0, yc, -socket_h + spec.socket_depth + cone_h / 2 - 0.01
        ) * bd.Cone(bore_r, 0.01, cone_h)

        # M5 cross-bolt through socket and pipe.
        p -= bd.Pos(0, yc, -socket_h + spec.socket_depth / 2) * (
            bd.Rot(Y=90)
            * bd.Cylinder(spec.mast_bolt_diameter / 2, 2 * outer_r + 4)
        )

    # Reflector mounting: vertical slots for spacing adjustment.
    slot_sk = bd.SlotCenterToCenter(
        spec.reflector_slot_travel, spec.bolt_hole_diameter, rotation=90
    )
    for x in (-spec.bolt_x_offset, spec.bolt_x_offset):
        p -= bd.Pos(x, 0, spec.reflector_z) * _y_prism(
            slot_sk, -spec.boom_thickness - 1, spec.boom_thickness + 2
        )

    # Driven-element bolt holes.
    top_bar_z = spec.driven_z + spec.driven_height / 2 - spec.bar_height / 2
    bottom_bar_z = spec.driven_z - spec.driven_height / 2 + spec.bar_height / 2
    p -= _element_bolt_holes(spec, top_bar_z)
    p -= _element_bolt_holes(spec, bottom_bar_z)

    # Coax pass-through from the back of the boom to just below the feed gap.
    coax_z = (
        spec.driven_z
        - spec.driven_height / 2
        - spec.coax_hole_diameter / 2
        - 3
    )
    p -= _y_hole(spec, 0, coax_z, spec.coax_hole_diameter)

    # Zip-tie slots (through X) to hold the coax against the back of the boom.
    reflector_top = spec.reflector_z + spec.bar_height / 2
    zip_zs = [
        spec.reflector_z
        - spec.bar_height / 2
        - spec.reflector_slot_travel / 2
        - 8,
        (reflector_top + spec.reflector_slot_travel / 2 + coax_z) / 2,
    ]
    for zz in zip_zs:
        p -= bd.Pos(
            0, -spec.boom_thickness + 2 + spec.zip_slot_width / 2, zz
        ) * bd.Box(
            spec.boom_width + 2, spec.zip_slot_width, spec.zip_slot_height
        )

    return p


# ---------------------------------------------------------------------------
# Top level
# ---------------------------------------------------------------------------


def assembly(spec: Spec) -> bd.Compound:
    """Everything in its assembled position.

    Purpose: for checking fit, not printing.
    """
    refl_r, refl_l = _split(spec, reflector_full(spec))
    drv_r, drv_l = _split(spec, driven_full(spec))
    return bd.Compound(children=[boom(spec), refl_r, refl_l, drv_r, drv_l])


def reflector_half(spec: Spec) -> bd.Part | bd.Compound:
    """One reflector half, laid flat for printing. Print 2."""
    return _for_printing(_split(spec, reflector_full(spec))[0])


def driven_half(spec: Spec) -> bd.Part | bd.Compound:
    """One folded-dipole half, laid flat for printing. Print 2."""
    return _for_printing(_split(spec, driven_full(spec))[0])


def log_electrical_summary(spec: Spec) -> None:
    """Log the electrical dimensions in wavelengths."""
    lam = spec.wavelength
    logger.info(f"Wavelength at {spec.frequency_mhz} MHz: {lam:.1f} mm")
    logger.info(
        f"Reflector:     {spec.reflector_length:.0f} mm = "
        f"{spec.reflector_length / lam:.3f} λ"
    )
    logger.info(
        f"Driven (FD):   {spec.driven_length:.0f} mm = "
        f"{spec.driven_length / lam:.3f} λ"
    )
    logger.info(
        f"FD conductor spacing: {spec.driven_height:.0f} mm = "
        f"{spec.driven_height / lam:.3f} λ"
    )
    logger.info(
        f"Element spacing: {spec.element_spacing:.0f} mm = "
        f"{spec.element_spacing / lam:.3f} λ"
    )


if __name__ == "__main__":
    spec = Spec()
    log_electrical_summary(spec)

    parts = {
        "assembly": show(assembly(spec)),
        "boom": show(boom(spec)),
        "reflector_half_print_x2": show(reflector_half(spec)),
        "driven_half_print_x2": show(driven_half(spec)),
    }

    logger.info("Showing CAD model(s)")

    (
        export_folder := Path(__file__).parent.parent
        / "build"
        / Path(__file__).stem
    ).mkdir(exist_ok=True, parents=True)
    for name, part in parts.items():
        assert isinstance(part, bd.Part | bd.Solid | bd.Compound), (
            f"{name} is not an expected type ({type(part)})"
        )
        if not part.is_manifold:
            logger.warning(f'Part "{name}" is not manifold')

        bd.export_stl(part, str(export_folder / f"{name}.stl"))
        bd.export_step(part, str(export_folder / f"{name}.step"))
