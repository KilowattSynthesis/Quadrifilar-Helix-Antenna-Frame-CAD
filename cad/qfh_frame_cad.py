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

    def __post_init__(self) -> None:
        """Validate the cylinder is tall enough to contain both helices."""
        ll = self.qfh.large_loop
        sl = self.qfh.small_loop
        required = max(ll.height, sl.height)
        assert self.cylinder_height >= required, (
            f"cylinder_height ({self.cylinder_height:.1f} mm) must be >= "
            f"max helix height ({required:.1f} mm)"
        )

    # ------------------------------------------------------------------
    # Derived quantities (computed from qfh + mechanical params)
    # ------------------------------------------------------------------

    @property
    def cylh2(self) -> float:
        """Half the cylinder height (the vertical mid-point)."""
        return self.cylinder_height / 2.0

    # Large-helix dimensions
    @property
    def d1(self) -> float:
        """Large-helix diameter (mm).

        In the QFH calculator, ``LoopResult.rad`` is the *horizontal
        separator D* (the full diameter of the helix coil, not its radius).
        """
        return self.qfh.large_loop.rad

    @property
    def hh1(self) -> float:
        """Large-helix height (mm)."""
        return self.qfh.large_loop.height

    # Small-helix dimensions
    @property
    def d2(self) -> float:
        """Small-helix diameter (mm).

        See :attr:`d1` - ``LoopResult.rad`` is the full diameter.
        """
        return self.qfh.small_loop.rad

    @property
    def hh2(self) -> float:
        """Small-helix height (mm)."""
        return self.qfh.small_loop.height

    # Wire-hole heights (centre of the cylinder ± half the helix height)
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

    # Helix sweep parameters
    @property
    def pitch1(self) -> float:
        """Pitch (mm/turn) for the large helix - 0.5 turns over hh1."""
        return 2.0 * self.hh1

    @property
    def pitch2(self) -> float:
        """Pitch (mm/turn) for the small helix - 0.5 turns over hh2."""
        return 2.0 * self.hh2

    @property
    def radius1(self) -> float:
        """Centreline radius of the large helix (mm)."""
        return self.d1 / 2.0

    @property
    def radius2(self) -> float:
        """Centreline radius of the small helix (mm)."""
        return self.d2 / 2.0

    # Elliptic support-tube radii
    @property
    def tube_r_outer(self) -> float:
        """Outer semi-axis of the support tube (mm)."""
        return self.d1 / 2.0 - self.wire_diameter / 2.0

    @property
    def tube_r_inner(self) -> float:
        """Inner semi-axis of the support tube (mm)."""
        return self.tube_r_outer - self.extrusion_width

    @property
    def tube_scale_y(self) -> float:
        """Y-axis scale factor to make the tube elliptic (D2/D1)."""
        return self.d2 / self.d1

    # Base grid derived quantities
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

    """
    with bd.BuildSketch(plane) as cs:
        bd.Circle(wire_diameter * 0.8)
        bd.Circle(wire_diameter / 2.0, mode=bd.Mode.SUBTRACT)
        with bd.Locations((wire_diameter * 1.5, 0)):
            bd.Rectangle(
                wire_diameter * 3.0, wire_diameter, mode=bd.Mode.SUBTRACT
            )
    return cs.sketch


def _helix_wire_solid(  # noqa: D417, PLR0913
    radius: float,
    pitch: float,
    height: float,
    wire_diameter: float,
    rot_z: float = 0.0,
    z_offset: float = 0.0,
) -> bd.Part:
    """Return one helical wire-channel solid.

    Builds a C-section profile swept along a right-hand 0.5-turn helix,
    then rotates it *rot_z* degrees around Z and translates it to *z_offset*.

    Parameters
    ----------
    radius:       Helix centreline radius (mm).
    pitch:        Helix pitch in mm/full-turn (``2 x height`` for 0.5 turns).
    height:       Axial span of the helix (mm).
    wire_diameter: Conductor outer diameter (mm) - sizes the C-groove.
    rot_z:        Azimuthal rotation before placement (degrees).
                  Use 0/180 for the large-helix pair, 90/270 for the small.
    z_offset:     Vertical translation to apply after rotation (mm).

    """
    # --- helix path ---
    with bd.BuildLine() as bl:
        bd.Helix(pitch=pitch, height=height, radius=radius)
    edge = bl.edges()[0]

    # --- cross-section at the helix start, perpendicular to the tangent ---
    start_plane = bd.Plane(
        origin=edge.location_at(0).position,
        z_dir=edge.tangent_at(0),
    )
    cs = _wirechannel_sketch(start_plane, wire_diameter)

    # --- sweep ---
    with bd.BuildPart() as wp:
        with bd.BuildLine():
            bd.Helix(pitch=pitch, height=height, radius=radius)
        bd.add(cs)
        bd.sweep()

    solid = wp.part
    assert solid
    if rot_z:
        solid = solid.rotate(bd.Axis.Z, rot_z)
    return solid.translate((0, 0, z_offset))


# ---------------------------------------------------------------------------
# Public model function
# ---------------------------------------------------------------------------


def qfh_antenna(spec: PartSpec) -> bd.Part:
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
    # Pre-build helical wire-channel solids (each requires its own helix
    # path, so they are created outside the main BuildPart context).
    # ------------------------------------------------------------------
    logger.debug("Building helical wire channels …")

    wire_channels = [
        _helix_wire_solid(
            spec.radius1,
            spec.pitch1,
            spec.hh1,
            spec.wire_diameter,
            rot_z=0,
            z_offset=spec.hwire11,
        ),
        _helix_wire_solid(
            spec.radius1,
            spec.pitch1,
            spec.hh1,
            spec.wire_diameter,
            rot_z=180,
            z_offset=spec.hwire11,
        ),
        _helix_wire_solid(
            spec.radius2,
            spec.pitch2,
            spec.hh2,
            spec.wire_diameter,
            rot_z=90,
            z_offset=spec.hwire12,
        ),
        _helix_wire_solid(
            spec.radius2,
            spec.pitch2,
            spec.hh2,
            spec.wire_diameter,
            rot_z=270,
            z_offset=spec.hwire12,
        ),
    ]

    # ------------------------------------------------------------------
    # Assemble everything inside a single BuildPart so the return type
    # stays bd.Part.
    # ------------------------------------------------------------------
    hole_length = 3.0 * spec.hh1  # long enough to pierce any cross-section

    # Alias for center alignment.
    _ca = (bd.Align.CENTER, bd.Align.CENTER, bd.Align.CENTER)

    with bd.BuildPart() as bp:
        # ---- elliptic support tube -----------------------------------
        logger.debug("Building elliptic support tube …")
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

        # ---- helical wire channels -----------------------------------
        for channel in wire_channels:
            bd.add(channel)

        # ---- base safety grid ----------------------------------------
        logger.debug("Building base safety grid …")
        grid_offset = spec.mesh_size / 2.0  # grid origin is at (-50, -50)
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

        # ---- wire entry / exit holes (lower pair) --------------------
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

        # ---- wire exit holes (upper pair) ----------------------------
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

        # ---- upper open slots (wire exits clear above cylinder top) --
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

    p = bp.part
    assert p
    return p


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    """Generate the QFH antenna support structure.

    Render, and export it to STL and STEP files.
    """
    # --- RF calculator ---------------------------------------------------
    input_spec = QfhInputSpec(
        freq=913.0,  # MHz
        wdiam=1.5,  # conductor outer diameter (mm)
        wrad=1.5,  # bending radius (mm)
        # Extremely-default settings:
        ratio=0.44,  # width / height ratio
        turns=0.5,  # half-turn helix
        nrwavel=1.0,  # one wavelength per loop
    )
    qfh_result = calculate_qfh(input_spec)

    # --- build & display -------------------------------------------------
    parts = {
        "QFH_Antenna_913_MHz": show(qfh_antenna(PartSpec(qfh=qfh_result))),
        "QFH_Antenna_436_MHz": show(
            qfh_antenna(
                PartSpec(
                    qfh=calculate_qfh(
                        QfhInputSpec(
                            freq=436.0,
                            wdiam=1.5,
                            wrad=1.5,
                        )
                    )
                )
            )
        ),
    }

    logger.info("Showing CAD model(s)")

    # --- export ----------------------------------------------------------
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
