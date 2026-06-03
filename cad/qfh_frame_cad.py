import math
from dataclasses import dataclass
from pathlib import Path

import build123d as bd
from build123d_ease import show
from loguru import logger

from cad.qfh_calc import QfhInputSpec, QfhResult, calculate_qfh


@dataclass
class PartSpec:
    """Mechanical specification for the QFH antenna support structure.

    All arguments are in millimeters (mm) except where noted. The RF geometry
    is read from the ``qfh`` argument, which is the output of the QFH
    calculator step.

    Parameters
    ----------
    qfh:
        Calculated RF geometry from :func:`calculate_qfh`.  All helix radii
        and heights are read from this object.
    wire_diameter:
        Outer diameter of the conductor that will be threaded through the
        channels (mm).  Drives channel bore and C-slot dimensions.
    cylinder_height:
        Total height of the support cylinder (mm).  Should be slightly
        taller than the tallest helix so both helices are fully enclosed.
    extrusion_width:
        Wall thickness of the elliptic support tube (mm).
        Typically set to one or two printer extrusion widths.
    pedestal_height:
        Thickness of the base safety-grid plate (mm).
    mesh_size:
        Overall footprint of the base grid (mm x mm).
    mesh_bars:
        Number of bars in each axis of the base grid.
    mesh_bar_width:
        Width of each bar in the base grid (mm).
    channel_outer_scale:
        ``r_outer = wire_diameter * channel_outer_scale``.  Controls the
        C-clip body size and wall thickness.
    channel_bore_clearance:
        ``r_bore = wire_diameter/2 + channel_bore_clearance``.  Clearance
        so the wire slides in without binding.
    channel_opening_ratio:
        Slot opening = ``wire_diameter * channel_opening_ratio``.  Must be
        < 1.0 so the opening is narrower than the wire and retains it.
        Smaller values grip more firmly but require more insertion force.

    """

    qfh: QfhResult

    # Mechanical / print parameters
    wire_diameter: float = 1.0  # conductor outer diameter
    cylinder_height: float = 120.0  # support cylinder height
    extrusion_width: float = 0.6  # tube wall thickness

    # Base safety-grid parameters
    pedestal_height: float = 1.0
    mesh_size: float = 100.0  # grid footprint (square)
    mesh_bars: int = 5  # bars per axis
    mesh_bar_width: float = 1.2

    # Wire-channel C-clip snap-fit geometry
    channel_outer_scale: float = 1.2
    channel_bore_clearance: float = 0.1
    channel_opening_ratio: float = 0.7

    def __post_init__(self) -> None:
        """Validate spec parameters."""
        required = max(self.qfh.large_loop.height, self.qfh.small_loop.height)
        assert self.cylinder_height >= required, (
            f"cylinder_height ({self.cylinder_height:.1f} mm) must be >= "
            f"max helix height ({required:.1f} mm)"
        )
        assert 0.0 < self.channel_opening_ratio < 1.0, (
            "channel_opening_ratio must be in (0, 1)"
        )

    # ------------------------------------------------------------------
    # Derived quantities
    # ------------------------------------------------------------------

    @property
    def cylh2(self) -> float:
        """Half the cylinder height (the vertical mid-point)."""
        return self.cylinder_height / 2.0

    @property
    def d1(self) -> float:
        """Large-helix full diameter (mm).

        In the QFH calculator, ``LoopResult.rad`` is the *horizontal
        separator D* (the full diameter of the helix coil, not its radius).
        """
        return self.qfh.large_loop.rad

    @property
    def hh1(self) -> float:
        """Large-helix height (mm)."""
        return self.qfh.large_loop.height

    @property
    def d2(self) -> float:
        """Small-helix full diameter (mm)."""
        return self.qfh.small_loop.rad

    @property
    def hh2(self) -> float:
        """Small-helix height (mm)."""
        return self.qfh.small_loop.height

    # Wire-hole heights (centre of cylinder +/- half helix height)
    @property
    def hwire11(self) -> float:
        """Lower wire-hole height for the large helix (mm)."""
        return self.cylh2 - self.hh1 / 2.0

    @property
    def hwire12(self) -> float:
        """Lower wire-hole height for the small helix (mm)."""
        return self.cylh2 - self.hh2 / 2.0

    @property
    def hwire21(self) -> float:
        """Upper wire-hole height for the large helix (mm)."""
        return self.cylh2 + self.hh1 / 2.0

    @property
    def hwire22(self) -> float:
        """Upper wire-hole height for the small helix (mm)."""
        return self.cylh2 + self.hh2 / 2.0

    # Elliptic support-tube geometry
    @property
    def tube_r_outer(self) -> float:
        """X semi-axis of the tube outer wall (mm)."""
        return self.d1 / 2.0 - self.wire_diameter / 2.0

    @property
    def tube_r_inner(self) -> float:
        """X semi-axis of the tube inner wall (mm)."""
        return self.tube_r_outer - self.extrusion_width

    @property
    def tube_scale_y(self) -> float:
        """Y/X axis ratio of the elliptic tube (D2/D1)."""
        return self.d2 / self.d1

    # Channel helix placement: centroid on tube outer wall
    @property
    def channel_helix_radius(self) -> float:
        """Radius of the circular helix path for each wire channel (mm).

        Placing the channel centroid exactly at ``tube_r_outer`` ensures the
        channel body always intersects the tube wall regardless of azimuthal
        angle, which is required for ``fuse(glue=True)`` to work correctly.
        """
        return self.tube_r_outer

    # Channel cross-section derived quantities
    @property
    def channel_r_outer(self) -> float:
        """Outer radius of the C-clip body (mm)."""
        return self.wire_diameter * self.channel_outer_scale

    @property
    def channel_r_bore(self) -> float:
        """Inner bore radius of the C-clip (mm)."""
        return self.wire_diameter / 2.0 + self.channel_bore_clearance

    @property
    def channel_opening(self) -> float:
        """Width of the C-clip slot opening (mm)."""
        return self.wire_diameter * self.channel_opening_ratio

    @property
    def channel_slot_x(self) -> float:
        """X coordinate where the slot begins, derived for exact opening width.

        Neat.
        """
        return math.sqrt(
            self.channel_r_bore**2 - (self.channel_opening / 2.0) ** 2
        )

    @property
    def channel_slot_rect_cx(self) -> float:
        """X centre of the slot-cutting rectangle (mm)."""
        return (
            self.channel_slot_x
            + (self.channel_r_outer - self.channel_slot_x) / 2.0
        )

    @property
    def channel_slot_rect_w(self) -> float:
        """Width of the slot-cutting rectangle (mm)."""
        return (self.channel_r_outer - self.channel_slot_x) * 2.0

    @property
    def channel_slot_rect_h(self) -> float:
        """Height of the slot-cutting rectangle (mm)."""
        return self.channel_r_bore * 2.0 + 0.1

    @property
    def mesh_space(self) -> float:
        """Clear gap between grid bars (mm)."""
        return (
            self.mesh_size - self.mesh_bar_width * self.mesh_bars
        ) / self.mesh_bars


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------


def _channel_face_at_zero(spec: PartSpec) -> bd.Face:
    """Build the C-clip cross-section face at phi=0 (centroid on +X axis).

    The face is a circular annulus with a slot cut from the outward (+X) face.
    The slot is narrower than the wire diameter, so the wire snaps in and is
    retained.  The face is centred at ``(channel_helix_radius, 0)`` in the XY
    plane and can be rotated to the desired start azimuth before extruding.

    """
    r = spec.channel_helix_radius
    with bd.BuildSketch(bd.Plane.XY) as sk, bd.Locations((r, 0)):
        bd.Circle(spec.channel_r_outer)
        bd.Circle(spec.channel_r_bore, mode=bd.Mode.SUBTRACT)
        with bd.Locations((spec.channel_slot_rect_cx, 0)):
            bd.Rectangle(
                spec.channel_slot_rect_w,
                spec.channel_slot_rect_h,
                mode=bd.Mode.SUBTRACT,
            )
    return sk.face()


def _make_channel(
    spec: PartSpec,
    phi_start_deg: float,
    height: float,
    z_offset: float,
) -> bd.Solid:
    """Return one helical wire-channel via ``extrude_linear_with_rotation``.

    The C-clip face is rotated to *phi_start_deg* then extruded upward by
    *height* with a -180° twist (one half-turn, clockwise = left-hand winding,
    matching the negative twist in the original OpenSCAD extrusion).

    Arguments:
    ---------
    spec:           Part specification.
    phi_start_deg:  Start azimuth on the tube surface (degrees).
                    0°/180° for the large-helix pair; 90°/270° for the small.
    height:         Axial span: ``spec.hh1`` or ``spec.hh2`` (mm).
    z_offset:      Bottom of the channel: ``spec.hwire11`` or ``spec.hwire12``.

    Notes:
    -----
    ``extrude_linear_with_rotation`` replaces the previous spline-path
    multisection sweep.  It is simpler, more reliable, and keeps the slot
    orientation correct throughout by construction: the face rotates uniformly
    with the extrusion twist, so the outward-facing slot stays outward at every
    height.

    The channel centroid radius equals ``tube_r_outer`` so the body always
    intersects the tube wall, enabling ``fuse(glue=True)`` to join them.

    """
    face = _channel_face_at_zero(spec)
    if phi_start_deg != 0:
        face = face.rotate(bd.Axis.Z, phi_start_deg)

    solid = bd.Solid.extrude_linear_with_rotation(
        section=face,
        center=(0, 0),
        normal=(0, 0, height),
        angle=-180.0,  # CW half-turn = left-hand winding
    )
    return solid.translate((0, 0, z_offset))


# ---------------------------------------------------------------------------
# Public model function
# ---------------------------------------------------------------------------


def qfh_antenna(spec: PartSpec) -> bd.Compound:
    """Create the QFH antenna support structure.

    The model consists of:

    * an elliptic support tube (height = ``cylinder_height``),
    * four helical wire channels (two large-loop, two small-loop),
    * horizontal through-holes where the wire enters and exits,
    * open slots at the top of the upper holes,
    * a safety-grid base plate.

    """
    logger.debug(
        (
            "QFH dimensions: "
            "D1={:.1f} mm, D2={:.1f} mm, H1={:.1f} mm, H2={:.1f} mm"
        ),
        spec.d1,
        spec.d2,
        spec.hh1,
        spec.hh2,
    )

    # ------------------------------------------------------------------
    # Elliptic support tube (straight extrusion, no twist needed).
    # ------------------------------------------------------------------
    logger.debug("Building elliptic support tube …")
    with bd.BuildSketch(bd.Plane.XY) as sk_tube:
        bd.Ellipse(spec.tube_r_outer, spec.tube_r_outer * spec.tube_scale_y)
        bd.Ellipse(
            spec.tube_r_inner,
            spec.tube_r_inner * spec.tube_scale_y,
            mode=bd.Mode.SUBTRACT,
        )
    with bd.BuildPart() as tube_ctx:
        bd.add(sk_tube.sketch)
        bd.extrude(amount=spec.cylinder_height)
    result = tube_ctx.part
    assert result

    # ------------------------------------------------------------------
    # Wire channels: extrude_linear_with_rotation produces a solid that
    # twists uniformly so the C-slot stays radially outward throughout.
    # fuse(glue=True) joins touching-but-not-overlapping surfaces cleanly.
    # ------------------------------------------------------------------
    logger.debug("Building and fusing helical wire channels …")

    channel_specs = [
        (0, spec.hh1, spec.hwire11),
        (180, spec.hh1, spec.hwire11),
        (90, spec.hh2, spec.hwire12),
        (270, spec.hh2, spec.hwire12),
    ]
    for phi, height, z_off in channel_specs:
        channel = _make_channel(spec, phi, height, z_off)
        result = result.fuse(channel, glue=True)
        assert isinstance(result, bd.Solid | bd.Part | bd.Compound)

    # ------------------------------------------------------------------
    # Base safety grid.
    # ------------------------------------------------------------------
    logger.debug("Building base safety grid …")
    grid_offset = spec.mesh_size / 2.0
    with bd.BuildPart() as grid_ctx:
        for i in range(1, spec.mesh_bars):
            xc = (
                i * (spec.mesh_bar_width + spec.mesh_space)
                - grid_offset
                + spec.mesh_bar_width / 2.0
            )
            with bd.Locations((xc, 0, spec.pedestal_height / 2.0)):
                bd.Box(
                    spec.mesh_bar_width, spec.mesh_size, spec.pedestal_height
                )
        for i in range(1, spec.mesh_bars):
            yc = (
                i * (spec.mesh_bar_width + spec.mesh_space)
                - grid_offset
                + spec.mesh_bar_width / 2.0
            )
            with bd.Locations((0, yc, spec.pedestal_height / 2.0)):
                bd.Box(
                    spec.mesh_size, spec.mesh_bar_width, spec.pedestal_height
                )

    grid_ctx_part = grid_ctx.part
    assert grid_ctx_part
    result = result.fuse(grid_ctx_part, glue=True)

    # ------------------------------------------------------------------
    # Wire entry/exit holes and upper open slots.
    # ------------------------------------------------------------------
    logger.debug("Subtracting wire holes and upper slots …")
    hole_length = 3.0 * spec.hh1
    _ca = (bd.Align.CENTER, bd.Align.CENTER, bd.Align.CENTER)

    for z, rot in [
        (spec.hwire11, (0, 90, 0)),
        (spec.hwire12, (90, 0, 0)),
        (spec.hwire21, (0, 90, 0)),
        (spec.hwire22, (90, 0, 0)),
    ]:
        assert isinstance(result, bd.Solid | bd.Part | bd.Compound)
        result = result - bd.Cylinder(
            spec.wire_diameter / 2.0,
            hole_length,
            rotation=rot,
            align=_ca,
        ).translate((0, 0, z))

    result = result - bd.Box(
        spec.cylinder_height,
        spec.wire_diameter,
        spec.cylinder_height,
        align=_ca,
    ).translate((0, 0, spec.hwire21 + spec.cylh2))
    result = result - bd.Box(
        spec.wire_diameter,
        spec.cylinder_height,
        spec.cylinder_height,
        align=_ca,
    ).translate((0, 0, spec.hwire22 + spec.cylh2))

    assert result
    return result


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    """Generate the QFH antenna support structure and export."""
    input_spec = QfhInputSpec(
        frequency_hz=913.0e6,  # MHz
        wire_diameter=1.5,  # conductor outer diameter (mm)
        wire_bending_radius=1.5,  # bending radius (mm)
        ratio=0.44,  # width / height ratio
        turns=0.5,  # half-turn helix
        num_wavelengths=1.0,  # one wavelength per loop
    )
    qfh_result = calculate_qfh(input_spec)

    parts = {
        "QFH_Antenna_913_MHz": show(qfh_antenna(PartSpec(qfh=qfh_result))),
        # "QFH_Antenna_436_MHz": show(
        #     qfh_antenna(
        #         PartSpec(
        #             qfh=calculate_qfh(
        #                 QfhInputSpec(
        #                     frequency_hz=436.0e6,
        #                     wire_diameter=1.5,
        #                     wire_bending_radius=1.5,
        #                 )
        #             )
        #         )
        #     )
        # ),
    }

    logger.info("Showing CAD model(s)")

    (export_folder := Path(__file__).parent / "build").mkdir(exist_ok=True)

    for name, part in parts.items():
        assert isinstance(part, bd.Part | bd.Solid | bd.Compound), (
            f"{name} is not an expected type ({type(part)})"
        )
        if not part.is_manifold:
            logger.warning('Part "{}" is not manifold', name)

        bd.export_stl(part, str(export_folder / f"{name}.stl"))
        bd.export_step(part, str(export_folder / f"{name}.step"))
        logger.info('Exported "{}" to {}', name, export_folder)


if __name__ == "__main__":
    main()
