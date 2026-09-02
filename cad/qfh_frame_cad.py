"""Build123d model of a 3D-printable QFH antenna frame for foil-tape radiators.

Design summary
--------------
Each of the two bifilar loops is a single **twisted blade**: a flat rectangular
slab, ``rad`` long and ``tape_land_width`` thick, extruded from z=0 to the loop
height while rotating ``turns`` revolutions.  The two blades cross at 90
deg and fuse into one printable frame.

The conductor is **self-adhesive foil tape stuck to the outside of the
frame** -- no wire channels, no threading.  For one loop, the tape path is::

    PCB pad -> underside of the bottom bar -> around the bottom corner
            -> up the outer end face (helical, 180 deg) -> over the top face
            -> down the opposite end face -> underside of the other bottom bar
            -> PCB pad

The blade's end faces sit at exactly +/-``rad``/2, so the tape centreline lands
on the design radius of the RF geometry (``LoopResult.rad`` is the
centre-to-centre diameter D).

Redundancy: **zip-tie holes** are drilled through the blade's wide faces near
each end and near the top/bottom.  A tie threads through a hole, wraps
around the nearby edge *over the tape*, and cinches -- so the tape is held
mechanically as well as by its adhesive.

At the bottom the frame carries a hub that does two jobs:

* a **mast sleeve** that slips over a 1.5 in (38.1 mm) OD PVC pipe, with radial
  screw holes for fastening into the pipe, and
* a **balun PCB mount** -- three M3 bosses on a 24 mm bolt circle matching the
  QFHBAL01 board (https://github.com/ODZ-UJF-AV-CR/QFHBAL01).  The PCB hangs
  under the hub inside the sleeve bore, on a flat underside that the four
  tape ends run straight across to reach the board's pads.  The board also
  sets how far the mast pipe can be pushed in.

Made with Claude Opus.
"""

import math
from dataclasses import dataclass
from pathlib import Path

import build123d as bd
from build123d_ease import show
from loguru import logger

from cad.qfh_calc import QfhInputSpec, QfhResult, calculate_qfh

MM_PER_INCH = 25.4


@dataclass
class PartSpec:
    """Mechanical specification for the QFH antenna support structure.

    All lengths are in millimeters (mm).

    The RF geometry is read from the ``qfh`` argument, which is the output of
    the QFH geometry calculation step.  Everything else here is mechanical:
    how wide the tape land is, where the zip-tie holes go, and how the frame
    mounts to the mast and to the balun PCB.
    """

    qfh: QfhResult

    # --- Blade (tape land) -------------------------------------------------
    # Thickness of the twisted blade == width of the flat face the foil tape
    # is stuck to.  Use tape a little narrower than this.
    tape_land_width: float = 8.0

    # --- Zip-tie holes -----------------------------------------------------
    tie_hole_diameter: float = 3.2  # Fits a standard 2.5 x 1.0 mm zip tie.
    # Distance from the blade's outer end face to the hole centre.  The tie
    # wraps around the end, so its loop is 2*inset + tape_land_width long.
    tie_hole_inset: float = 5.0
    # Ties along each helical end face.
    tie_holes_per_helix: int = 6
    tie_hole_z_margin: float = 15.0  # Keep-out at the very top/bottom.
    # Ties that wrap the horizontal top/bottom faces, where the tape turns
    # the corner and runs inboard.  Given as fractions along the bar's bare
    # span (from the hub, or from the axis at the top, out to the end).
    tie_hole_bar_radius_fractions: tuple[float, ...] = (0.35, 0.75)
    tie_hole_bar_z: float = 8.0  # Height above bottom / below top.

    # --- Mast sleeve -------------------------------------------------------
    # Set to None for no mast sleeve (e.g. antennas too small to straddle a
    # pipe); the hub then reduces to a plate carrying just the PCB bosses.
    mast_pipe_od: float | None = 1.5 * MM_PER_INCH  # 38.1 mm: PVC pipe OD.
    mast_bore_clearance: float = 0.4  # Diametral slip fit.
    mast_sleeve_wall: float = 4.0
    mast_sleeve_length: float = 40.0  # Hangs below the frame.
    mast_screw_hole_diameter: float = 3.4  # M3 clearance, into the pipe.
    mast_screw_count: int = 3
    mast_screw_z_from_sleeve_end: float = 12.0

    # Hub plate: fills the 45 deg gaps between the bottom bars so the sleeve
    # and the PCB bosses have something to hang from.
    hub_plate_thickness: float = 6.0
    # Minimum bare bottom bar left outboard of the hub, so the hub never
    # shrouds the blade's end face (the helical tape land).
    hub_edge_margin: float = 4.0

    # --- Balun PCB mount (QFHBAL01) ---------------------------------------
    # Board has 3x 3.2 mm holes on a 24 mm bolt circle, at 90 deg spacing with
    # the fourth quadrant taken by the RF connector.  The board is fitted
    # rotated 45 deg to the frame's bars so its X-shaped pads line up with the
    # four tape ends; that puts the bosses between the bars.
    pcb_screw_circle_diameter: float = 24.0
    pcb_screw_angles_deg: tuple[float, ...] = (45.0, 135.0, 225.0)
    pcb_boss_diameter: float = 7.0
    pcb_standoff_height: float = 6.0  # Gap from frame underside to PCB.
    pcb_screw_hole_diameter: float = 2.7  # M3 thread-forming into plastic.
    pcb_screw_hole_depth: float = 8.0

    def __post_init__(self) -> None:
        """Validate spec parameters."""
        min_half_len = (
            min(self.qfh.small_loop.rad, self.qfh.large_loop.rad) / 2.0
        )
        # Leave the blade ends (the tape lands) clear of the hub.
        if self.hub_radius > min_half_len - self.hub_edge_margin:
            msg = (
                f"Hub (r={self.hub_radius:.1f}) crowds the tape lands at "
                f"r={min_half_len:.1f}. Use a smaller mast_pipe_od, or set "
                f"mast_pipe_od=None for no mast sleeve."
            )
            raise ValueError(msg)

        if self.mast_pipe_od is not None and (
            self.pcb_boss_outer_radius > self.mast_bore_radius
        ):
            msg = (
                f"PCB bosses (to r={self.pcb_boss_outer_radius:.1f}) do not "
                f"fit inside the mast bore "
                f"(r={self.mast_bore_radius:.1f}); the board must hang "
                f"inside the sleeve."
            )
            raise ValueError(msg)

        if self.pcb_screw_hole_depth > (
            self.pcb_standoff_height + self.hub_plate_thickness
        ):
            msg = "PCB screw holes would break through the top of the hub."
            raise ValueError(msg)

        if self.tie_hole_inset < self.tie_hole_diameter:
            msg = "Zip-tie holes are too close to the blade's end face."
            raise ValueError(msg)

    @property
    def mast_bore_radius(self) -> float:
        """Radius of the hole the PVC pipe slides into."""
        assert self.mast_pipe_od is not None
        return (self.mast_pipe_od + self.mast_bore_clearance) / 2.0

    @property
    def pcb_boss_outer_radius(self) -> float:
        """Radius reached by the outside of the PCB mounting bosses."""
        return (self.pcb_screw_circle_diameter + self.pcb_boss_diameter) / 2.0

    @property
    def hub_radius(self) -> float:
        """Outer radius of the mast sleeve / hub plate."""
        if self.mast_pipe_od is None:
            return self.pcb_boss_outer_radius + 1.5
        return self.mast_bore_radius + self.mast_sleeve_wall


def _polar(radius: float, angle_deg: float) -> tuple[float, float]:
    """Cartesian (x, y) for a polar coordinate, in mm and degrees."""
    return (
        radius * math.cos(math.radians(angle_deg)),
        radius * math.sin(math.radians(angle_deg)),
    )


def _tie_hole(
    *,
    radius: float,
    angle_deg: float,
    z: float,
    diameter: float,
    length: float,
) -> bd.Part:
    """One zip-tie hole: axis horizontal, across the blade's wide faces.

    The blade at height ``z`` lies along ``angle_deg``; the hole is drilled
    perpendicular to it, so a tie threaded through can wrap around the nearby
    edge and over the tape.
    """
    cyl = bd.Cylinder(
        radius=diameter / 2.0, height=length, rotation=(90, 0, 0)
    )
    return cyl.rotate(bd.Axis.Z, angle_deg).translate(
        (*_polar(radius, angle_deg), z)
    )


def _draw_twisted_blade(
    *,
    loop_diameter: float,
    loop_height: float,
    spec: PartSpec,
) -> bd.Part | bd.Compound:
    """One loop's twisted blade, with its zip-tie holes.

    The blade's two end faces are the helical tape lands; they sit at
    +/-``loop_diameter``/2, i.e. exactly on the RF design radius.
    """
    thickness = spec.tape_land_width
    turns = spec.qfh.input_spec.turns
    half_len = loop_diameter / 2.0

    section = bd.Rectangle(loop_diameter, thickness).face()
    assert section is not None

    p = bd.Part(None)
    p += bd.Solid.extrude_linear_with_rotation(
        section=section,
        center=(0, 0),
        normal=(0, 0, loop_height),  # Distance.
        angle=(360 * turns),
    )

    # Twist is CCW with height: the blade at height z lies along this angle.
    def blade_angle_deg(z: float) -> float:
        return 360.0 * turns * z / loop_height

    hole_len = thickness * 3.0
    holes = bd.Part(None)

    # Ties along the two helical end faces (the main tape runs).
    n = spec.tie_holes_per_helix
    z_lo, z_hi = spec.tie_hole_z_margin, loop_height - spec.tie_hole_z_margin
    for i in range(n):
        z = z_lo + (z_hi - z_lo) * (i / (n - 1) if n > 1 else 0.5)
        ang = blade_angle_deg(z)
        for end_sign in (0.0, 180.0):
            holes += _tie_hole(
                radius=half_len - spec.tie_hole_inset,
                angle_deg=ang + end_sign,
                z=z,
                diameter=spec.tie_hole_diameter,
                length=hole_len,
            )

    # Ties that wrap the flat bottom / top faces, where the tape turns
    # inboard.  At the bottom the hub covers the inner part of the bar, so a
    # tie can only wrap outboard of it; at the top the bar is bare.
    r_out = half_len - spec.tie_hole_inset
    bar_spans = {
        spec.tie_hole_bar_z: (spec.hub_radius + spec.tie_hole_inset, r_out),
        loop_height - spec.tie_hole_bar_z: (0.0, r_out),
    }
    for z, (r_lo, r_hi) in bar_spans.items():
        if r_hi - r_lo < spec.tie_hole_diameter * 3:
            continue  # No bare bar to wrap a tie around.
        ang = blade_angle_deg(z)
        for frac in spec.tie_hole_bar_radius_fractions:
            for end_sign in (0.0, 180.0):
                holes += _tie_hole(
                    radius=r_lo + (r_hi - r_lo) * frac,
                    angle_deg=ang + end_sign,
                    z=z,
                    diameter=spec.tie_hole_diameter,
                    length=hole_len,
                )

    p -= holes
    return p


def _draw_hub(spec: PartSpec) -> bd.Part | bd.Compound:
    """Mast sleeve + hub plate + balun-PCB bosses, all at/below the bottom.

    The hub plate (z = 0 .. hub_plate_thickness) ties the two bottom bars
    together and closes off the 45 deg gaps.  The sleeve, if any, hangs below
    z = 0 and slides over the mast pipe.  The PCB bosses hang below z = 0 too,
    *inside* the sleeve bore -- so they have to be added after the bore is
    cut, and the PCB itself becomes the mast pipe's insertion depth stop.
    """
    outer_r = spec.hub_radius
    has_sleeve = spec.mast_pipe_od is not None

    p = bd.Part(None)

    # Hub plate, sitting on top of z=0 so its underside is flush with the
    # bottom bars (the tape runs along that flat underside).  It also caps
    # the top of the bore, which is what the PCB bosses hang from.
    p += bd.Pos(Z=spec.hub_plate_thickness / 2.0) * bd.Cylinder(
        radius=outer_r, height=spec.hub_plate_thickness
    )

    if has_sleeve:
        bore_r = spec.mast_bore_radius
        r_mid = (bore_r + outer_r) / 2.0

        p += bd.Pos(Z=-spec.mast_sleeve_length / 2.0) * bd.Cylinder(
            radius=outer_r, height=spec.mast_sleeve_length
        )

        # Bore for the mast pipe (also the space the PCB hangs in).
        p -= bd.Pos(Z=-spec.mast_sleeve_length / 2.0) * bd.Cylinder(
            radius=bore_r, height=spec.mast_sleeve_length
        )

        # Radial screw holes, for screwing the sleeve to the mast pipe.
        z_screw = -(
            spec.mast_sleeve_length - spec.mast_screw_z_from_sleeve_end
        )
        for i in range(spec.mast_screw_count):
            # Offset off the bar axes, clear of the tape runs above.
            ang = 360.0 * i / spec.mast_screw_count + 45.0
            p -= (
                bd.Cylinder(
                    radius=spec.mast_screw_hole_diameter / 2.0,
                    height=spec.mast_sleeve_wall * 3.0,
                    rotation=(0, 90, 0),  # Axis along X, radial once rotated.
                )
                .rotate(bd.Axis.Z, ang)
                .translate((*_polar(r_mid, ang), z_screw))
            )

    # PCB standoff bosses, hanging below the hub underside, inside the bore.
    for ang in spec.pcb_screw_angles_deg:
        p += bd.Pos(
            *_polar(spec.pcb_screw_circle_diameter / 2.0, ang),
            Z=-spec.pcb_standoff_height / 2.0,
        ) * bd.Cylinder(
            radius=spec.pcb_boss_diameter / 2.0,
            height=spec.pcb_standoff_height,
        )

    # Blind, thread-forming screw holes up into the bosses for the PCB.
    for ang in spec.pcb_screw_angles_deg:
        p -= bd.Pos(
            *_polar(spec.pcb_screw_circle_diameter / 2.0, ang),
            Z=-spec.pcb_standoff_height + spec.pcb_screw_hole_depth / 2.0,
        ) * bd.Cylinder(
            radius=spec.pcb_screw_hole_diameter / 2.0,
            height=spec.pcb_screw_hole_depth,
        )

    return p


def qfh_antenna_frame(spec: PartSpec) -> bd.Part | bd.Compound:
    """Create the QFH antenna support structure."""
    p = bd.Part(None)

    # Large loop blade: bottom bar along +/-X.
    p += _draw_twisted_blade(
        loop_diameter=spec.qfh.large_loop.rad,
        loop_height=spec.qfh.large_loop.height,
        spec=spec,
    )

    # Small loop blade: a quarter turn around, bottom bar along +/-Y.
    p += _draw_twisted_blade(
        loop_diameter=spec.qfh.small_loop.rad,
        loop_height=spec.qfh.small_loop.height,
        spec=spec,
    ).rotate(bd.Axis.Z, 90)

    # The hub lives entirely at and below z=0, while the blades run upward
    # from z=0, so fusing it in last cannot backfill any of its holes.
    p += _draw_hub(spec)

    return p


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    """Generate the QFH antenna support structure and export."""
    parts = {
        # 436 MHz: ~100 mm across, so it straddles a 1.5 in PVC mast.
        "QFH_Antenna_436_MHz": show(
            qfh_antenna_frame(
                PartSpec(
                    qfh=calculate_qfh(
                        QfhInputSpec(
                            frequency_hz=436.0e6,
                            wire_diameter=1.5,
                            wire_bending_radius=3.0,
                        )
                    )
                )
            )
        ),
        # 913 MHz: only ~46 mm across -- narrower than a 1.5 in pipe -- so
        # there is no mast sleeve; it just carries the balun PCB.
        "QFH_Antenna_913_MHz": (
            qfh_antenna_frame(
                PartSpec(
                    qfh=calculate_qfh(
                        QfhInputSpec(
                            frequency_hz=913.0e6,
                            wire_diameter=1.5,  # Conductor outer dia (mm).
                            wire_bending_radius=3.0,  # Bending radius (mm).
                            ratio=0.44,  # Width / height ratio.
                            turns=0.5,  # Half-turn helix.
                            num_wavelengths=1.0,  # One wavelength per loop.
                        )
                    ),
                    mast_pipe_od=None,
                )
            )
        ),
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
            logger.warning('Part "{}" is not manifold', name)

        bd.export_stl(part, str(export_folder / f"{name}.stl"))
        bd.export_step(part, str(export_folder / f"{name}.step"))
        logger.info('Exported "{}" to {}', name, export_folder)


if __name__ == "__main__":
    main()
