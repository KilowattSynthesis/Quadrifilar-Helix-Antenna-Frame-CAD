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
    helix_spline_pts:
        Number of interpolation points used when building each elliptic helix
        spline.  Higher values give a smoother groove at the cost of longer
        build time.  64 is a good default; raise to 120+ for final export.

    """

    qfh: QfhResult

    # Mechanical / print parameters
    wire_diameter: float = 1.0  # conductor outer diameter
    cylinder_height: float = 120.0  # support cylinder height
    extrusion_width: float = 0.6  # tube wall thickness

    # Base safety-grid parameters.
    pedestal_height: float = 1.0
    mesh_size: float = 100.0  # grid footprint (square)
    mesh_bars: int = 5  # bars per axis
    mesh_bar_width: float = 1.2  # bar width

    helix_spline_pts: int = 64  # spline path resolution
    helix_n_sections: int = 8  # number of cross-section planes along the helix

    # Note: helix_n_sections > 8 may trigger a BRepOffsetAPI_MakePipeShell
    # failure in the underlying OCCT kernel for this profile shape; 8 gives
    # adequate smoothness.

    def __post_init__(self) -> None:
        """Validate the cylinder is tall enough to contain both helices."""
        required = max(self.qfh.large_loop.height, self.qfh.small_loop.height)
        assert self.cylinder_height >= required, (
            f"cylinder_height ({self.cylinder_height:.1f} mm) must be >= "
            f"max helix height ({required:.1f} mm)"
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
        """Small-helix full diameter (mm).

        See :attr:`d1` - ``LoopResult.rad`` is the full diameter.
        """
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

    # Elliptic helix-path semi-axes
    @property
    def helix_a(self) -> float:
        """X semi-axis of the wire-channel centreline ellipse (mm).

        Exactly ``tube_r_outer + wire_diameter/2``: the channel centre sits
        one wire radius proud of the tube outer wall along X.
        """
        return self.tube_r_outer + self.wire_diameter / 2.0

    @property
    def helix_b(self) -> float:
        """Y semi-axis of the wire-channel centreline ellipse (mm).

        Exactly ``tube_r_outer * tube_scale_y + wire_diameter/2``: the same
        constant normal offset from the tube outer wall along Y.
        """
        return self.tube_r_outer * self.tube_scale_y + self.wire_diameter / 2.0

    @property
    def mesh_space(self) -> float:
        """Clear gap between grid bars (mm)."""
        return (
            self.mesh_size - self.mesh_bar_width * self.mesh_bars
        ) / self.mesh_bars


# ---------------------------------------------------------------------------
# Internal geometry helpers
# ---------------------------------------------------------------------------


def _wirechannel_sketch(plane: bd.Plane, wire_diameter: float) -> bd.Sketch:
    """C-shaped groove cross-section, normal to *plane*.

    The profile matches the original OpenSCAD ``wirechannel()`` module:

    * outer disc  - radius ``wire_diameter x 0.8``
    * bore        - radius ``wire_diameter / 2``   (subtracted)
    * open slot   - rectangle on the outward face  (subtracted)

    The *plane* x-axis must point radially outward so the open slot faces
    away from the antenna axis.  See :func:`_radial_plane_at_edge_start`.

    """
    with bd.BuildSketch(plane) as cs:
        bd.Circle(wire_diameter * 0.8)
        bd.Circle(wire_diameter / 2.0, mode=bd.Mode.SUBTRACT)
        with bd.Locations((wire_diameter * 1.5, 0)):
            bd.Rectangle(
                wire_diameter * 3.0, wire_diameter, mode=bd.Mode.SUBTRACT
            )
    return cs.sketch


def _elliptic_helix_points(
    a: float,
    b: float,
    height: float,
    phi_start_deg: float,
    n_pts: int,
) -> list[tuple[float, float, float]]:
    """Compute points along a left-hand (CW) elliptic half-turn helix.

    The helix traces the ellipse  ``x = a·cos(φ)``, ``y = b·sin(φ)``
    starting at parametric angle *phi_start_deg* and winding clockwise
    (decreasing φ) for exactly half a turn (π radians), rising linearly
    from ``z = 0`` to ``z = height``.

    Using *phi_start_deg* (0°, 90°, 180°, 270°) to place all four channel
    centrelines on the *same* ellipse ensures each stays at a constant
    ``wire_diameter / 2`` offset from the tube outer wall at every azimuth.

    """
    phi0 = math.radians(phi_start_deg)
    return [
        (
            a * math.cos(phi0 - math.pi * i / n_pts),
            b * math.sin(phi0 - math.pi * i / n_pts),
            height * i / n_pts,
        )
        for i in range(n_pts + 1)
    ]


def _helix_wire_solid(
    spec: PartSpec,
    phi_start_deg: float,
    height: float,
    z_offset: float,
) -> bd.Solid | bd.Part:
    """Return one helical wire-channel solid.

    Uses a **multisection sweep**: the C-groove cross-section is placed at
    ``spec.helix_n_sections + 1`` positions along the elliptic helix path,
    each independently oriented so its x-axis points radially outward.
    The sections are then swept together along the spline path.

    Parameters
    ----------
    spec:
        Part specification (provides ``helix_a``, ``helix_b``,
        ``wire_diameter``, ``helix_spline_pts``, ``helix_n_sections``).
    phi_start_deg:
        Start angle on the elliptic helix path (degrees).
        0° / 180° for the large-helix pair; 90° / 270° for the small pair.
    height:
        Axial span of the helix: ``spec.hh1`` or ``spec.hh2`` (mm).
    z_offset:
        Vertical translation: ``spec.hwire11`` or ``spec.hwire12`` (mm).

    Notes
    -----
    **Why multisection sweep is required**
        A single-section ``sweep()`` uses parallel transport (TRANSFORMED
        mode) to propagate the initial profile orientation along the path.
        For an elliptic helix the Frenet frame precesses, so the C-slot —
        which starts radially outward — rotates by up to ~150° by the time
        it reaches the far end, ending up buried between the channel wall
        and the tube body.  Placing independently-oriented sections at
        regular intervals and using ``multisection=True`` forces the slot
        to remain radially outward at every point.

    **Elliptic helix path**
        The support tube is elliptic (semi-axes ``tube_r_outer`` x Y).
        A circular path at radius ``d1/2`` drifts up to ~1.6 mm from the
        tube wall at the minor-axis positions.  The elliptic path with
        semi-axes ``helix_a`` and ``helix_b`` gives an exact
        ``wire_diameter/2`` offset at every azimuth.

    **Left-hand (CW) winding**
        Matches the negative twist of the original OpenSCAD extrusion.

    """
    pts = _elliptic_helix_points(
        spec.helix_a,
        spec.helix_b,
        height,
        phi_start_deg,
        spec.helix_spline_pts,
    )

    with bd.BuildLine() as bl:
        bd.Spline(*pts)
    path_wire = bl.wire()
    path_edge = bl.edges()[0]

    # Build one cross-section at each equally-spaced station, each with its
    # x-axis pointing radially outward at that point on the elliptic path.
    sections: list[bd.Sketch] = []
    for i in range(spec.helix_n_sections + 1):
        t = i / spec.helix_n_sections
        phi = math.radians(phi_start_deg) - math.pi * t
        pos = bd.Vector(
            spec.helix_a * math.cos(phi),
            spec.helix_b * math.sin(phi),
            height * t,
        )
        tan = path_edge.tangent_at(t)
        radial = bd.Vector(pos.X, pos.Y, 0).normalized()
        radial_perp = (radial - tan * radial.dot(tan)).normalized()
        plane = bd.Plane(origin=pos, x_dir=radial_perp, z_dir=tan)
        sections.append(_wirechannel_sketch(plane, spec.wire_diameter))

    with bd.BuildPart() as wp:
        bd.sweep(sections=sections, path=path_wire, multisection=True)

    solid = wp.part
    assert solid
    return solid.translate((0, 0, z_offset))


# ---------------------------------------------------------------------------
# Public model function
# ---------------------------------------------------------------------------


def qfh_antenna(spec: PartSpec) -> bd.Solid | bd.Part:
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
    logger.debug(
        "Elliptic helix path: a={:.3f} mm, b={:.3f} mm",
        spec.helix_a,
        spec.helix_b,
    )

    # ------------------------------------------------------------------
    # Pre-build helical wire-channel solids.
    # They are built outside BuildPart because bd.add() fuses external
    # pre-built solids incorrectly when they don't overlap the existing
    # context body (the channels are flush against the tube wall, not
    # embedded in it).  They are merged via .fuse() after the base body
    # is complete.
    # ------------------------------------------------------------------
    logger.debug("Building helical wire channels …")

    wire_channels = [
        _helix_wire_solid(
            spec, phi_start_deg=0, height=spec.hh1, z_offset=spec.hwire11
        ),
        _helix_wire_solid(
            spec, phi_start_deg=180, height=spec.hh1, z_offset=spec.hwire11
        ),
        _helix_wire_solid(
            spec, phi_start_deg=90, height=spec.hh2, z_offset=spec.hwire12
        ),
        _helix_wire_solid(
            spec, phi_start_deg=270, height=spec.hh2, z_offset=spec.hwire12
        ),
    ]

    # ------------------------------------------------------------------
    # Build the base body: support tube + safety grid + all subtractions.
    # Everything that can be done inside one BuildPart context is done
    # here; channels are fused in afterwards.
    # ------------------------------------------------------------------
    hole_length = 3.0 * spec.hh1
    _ca = (bd.Align.CENTER, bd.Align.CENTER, bd.Align.CENTER)
    grid_offset = spec.mesh_size / 2.0

    logger.debug("Building elliptic support tube and base grid …")

    with bd.BuildPart() as bp:
        # ---- elliptic support tube -----------------------------------
        with bd.BuildSketch(bd.Plane.XY):
            bd.Ellipse(
                spec.tube_r_outer, spec.tube_r_outer * spec.tube_scale_y
            )
            bd.Ellipse(
                spec.tube_r_inner,
                spec.tube_r_inner * spec.tube_scale_y,
                mode=bd.Mode.SUBTRACT,
            )
        bd.extrude(amount=spec.cylinder_height)

        # ---- base safety grid ----------------------------------------
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

        # ---- wire entry / exit holes ---------------------------------
        logger.debug("Subtracting wire holes and upper slots …")
        with bd.Locations((0, 0, spec.hwire11)):
            bd.Cylinder(
                spec.wire_diameter / 2.0,
                hole_length,
                rotation=(0, 90, 0),
                align=_ca,
                mode=bd.Mode.SUBTRACT,
            )
        with bd.Locations((0, 0, spec.hwire12)):
            bd.Cylinder(
                spec.wire_diameter / 2.0,
                hole_length,
                rotation=(90, 0, 0),
                align=_ca,
                mode=bd.Mode.SUBTRACT,
            )
        with bd.Locations((0, 0, spec.hwire21)):
            bd.Cylinder(
                spec.wire_diameter / 2.0,
                hole_length,
                rotation=(0, 90, 0),
                align=_ca,
                mode=bd.Mode.SUBTRACT,
            )
        with bd.Locations((0, 0, spec.hwire22)):
            bd.Cylinder(
                spec.wire_diameter / 2.0,
                hole_length,
                rotation=(90, 0, 0),
                align=_ca,
                mode=bd.Mode.SUBTRACT,
            )

        # ---- upper open slots ----------------------------------------
        with bd.Locations((0, 0, spec.hwire21 + spec.cylh2)):
            bd.Box(
                spec.cylinder_height,
                spec.wire_diameter,
                spec.cylinder_height,
                mode=bd.Mode.SUBTRACT,
            )
        with bd.Locations((0, 0, spec.hwire22 + spec.cylh2)):
            bd.Box(
                spec.wire_diameter,
                spec.cylinder_height,
                spec.cylinder_height,
                mode=bd.Mode.SUBTRACT,
            )

    # ------------------------------------------------------------------
    # Fuse the pre-built wire channels into the base body.
    # Using .fuse() rather than bd.add() inside BuildPart because the
    # channels sit flush against (not embedded in) the tube wall, so
    # bd.add() would incorrectly report zero intersection and discard
    # the larger solid.
    # ------------------------------------------------------------------
    logger.debug("Fusing wire channels into base body …")
    result = bp.part
    assert isinstance(result, bd.Part | bd.Solid)
    for channel in wire_channels:
        result = result.fuse(channel)
        assert isinstance(result, bd.Part | bd.Solid)

    assert result
    return result


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    """Generate the QFH antenna support structure and export."""
    # --- RF calculator ---------------------------------------------------
    input_spec = QfhInputSpec(
        freq=913.0,  # MHz
        wdiam=1.5,  # conductor outer diameter (mm)
        wrad=1.5,  # bending radius (mm)
        ratio=0.44,  # width / height ratio
        turns=0.5,  # half-turn helix
        nrwavel=1.0,  # one wavelength per loop
    )
    qfh_result = calculate_qfh(input_spec)

    # --- build & display -------------------------------------------------
    parts = {
        "QFH_Antenna_913_MHz": show(qfh_antenna(PartSpec(qfh=qfh_result))),
        # "QFH_Antenna_436_MHz": show(
        #     qfh_antenna(
        #         PartSpec(
        #             qfh=calculate_qfh(
        #                 QfhInputSpec(
        #                     freq=436.0,
        #                     wdiam=1.5,
        #                     wrad=1.5,
        #                 )
        #             )
        #         )
        #     )
        # ),
    }

    logger.info("Showing CAD model(s)")

    # --- export ----------------------------------------------------------
    (export_folder := Path(__file__).parent.parent / "build").mkdir(
        exist_ok=True
    )

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
