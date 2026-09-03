"""Build123d model of a 3D-printable QFH antenna frame for foil-tape radiators.

Design summary
--------------
Each of the two bifilar loops is a single **twisted blade**: a slab ``rad``
long, extruded from z=0 to the loop height while rotating ``turns``
revolutions.  Its section is thin (``blade_core_thickness``) for most of its
length and flares out over the last ``tape_land_flare_length`` at each end to
the full ``tape_land_width``, so the material goes where the tape needs it.
The flare only widens the blade radially, so the top and bottom of every
blade is a full-width land instead, giving the tape something to sit on where
it turns inboard.  The bottom land is the thicker of the two: it runs up to
the outside top of the PCB enclosure, so the base of the antenna is one solid
full-width slab flush with the hub.  The two blades cross at 90
deg and fuse into one printable frame.

The conductor is **self-adhesive foil tape stuck to the outside of the
frame** -- no wire channels, no threading.  For one loop, the tape path is::

    PCB pad -> wire through the feed-through -> underside of the bottom bar
            -> around the bottom corner
            -> up the outer end face (helical, 180 deg) -> over the top face
            -> down the opposite end face -> underside of the other bottom bar
            -> feed-through -> PCB pad

The blade's end faces sit at exactly +/-``rad``/2, so the tape centreline lands
on the design radius of the RF geometry (``LoopResult.rad`` is the
centre-to-centre diameter D).

Redundancy: **zip-tie holes** are drilled through the blade's wide faces near
each end and near the top/bottom.  A tie threads through a hole, wraps
around the nearby edge *over the tape*, and cinches -- so the tape is held
mechanically as well as by its adhesive.

The two loops are different heights, so at the top the taller blade stands
right where the shorter loop's tape has to cross the axis.  A **crossover
window** through the taller blade, sitting on the shorter blade's top face,
lets that tape run straight through.  It is shortened automatically if the
loop heights are too close to leave a bridge above it.

At the bottom the frame carries a hub that does two jobs:

* a **mast sleeve** that slips over a 1.5 in (38.1 mm) OD PVC pipe, with radial
  screw holes for fastening into the pipe, and
* a **balun PCB mount** -- four M3 bosses on a 24 mm bolt circle matching the
  QFHBAL01 board (https://github.com/ODZ-UJF-AV-CR/QFHBAL01), which has three
  holes plus an RF connector in the fourth quadrant, so one boss goes unused
  and the board can be fitted in any rotation.  The PCB hangs
  under the hub inside the sleeve bore, which the sleeve wall closes off
  from the tape runs outside it.  Each of the four bars gets one 1.5 mm
  2.5 mm **wire feed-through**, tucked just under the surface the tape runs
  along: the tape ends outside and a short copper wire runs through to the
  board's pad, which is far easier to weatherproof than an open slot.  The
  board also sets how far the pipe can be pushed in.

Taller frames do not fit a print bed, so ``qfh_antenna_frame_sections``
cuts the frame into as few equal horizontal **sections** as will fit under
``max_print_height`` (default 200 mm, counting the balls that stand proud of
each cut).  Locating pins are 3 mm **balls** sitting on the cut plane, half
proud of one section and half dished into the other: a pair on the axis,
inside a boss that also ties the two blades together, and more out on both
blades' arms.  By default the balls belong to the section *above* each cut,
which leaves the section below with a dish that opens upward; the sections
carrying balls want flipping on the print bed.  A ball reaches only its own
radius clear of the plane, which is what lets pins go out on the arms at all
-- a pin standing any real height there would twist straight out the side of
the blade as the blade turns.
Zip-tie anchor holes flank each cut further out, where they have the lever
arm to resist bending.  Set ``max_print_height=None`` for one piece.

Made with Claude Opus.
"""

import itertools
import math
from dataclasses import dataclass
from pathlib import Path

import build123d as bd
from build123d_ease import show
from loguru import logger

from cad.qfh_calc import QfhInputSpec, QfhResult, calculate_qfh

MM_PER_INCH = 25.4

# Azimuths of the four bottom bars: the large loop's blade lies along +/-X at
# z=0 and the small loop's, a quarter turn around, along +/-Y.  The four tape
# runs head inboard along these.
BAR_ANGLES_DEG = (0.0, 90.0, 180.0, 270.0)

# Which build ``main`` puts in the viewer, and how far apart its sections are
# stood.  Everything is still exported; this only picks what is displayed.
SHOWN_PART_PREFIX = "QFH_Antenna_436_MHz"
SHOWN_PART_GAP = 20.0


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
    # Width of the flat outer end face the foil tape is stuck to.  Use tape a
    # little narrower than this.
    tape_land_width: float = 10.0
    # Thickness of the blade away from the ends.  The section is this thin
    # for most of its length and flares out to the full tape land only where
    # the tape actually needs it, which is most of the material saved.
    blade_core_thickness: float = 5.0
    # Length, at each end, over which the core flares out to the tape land.
    tape_land_flare_length: float = 10.0
    # The flare only widens the blade radially, so the flat top and bottom
    # faces would be left at core width, and tape turning inboard across them
    # would overhang.  These pads take the top and bottom of each blade out to
    # the full tape land width, across the whole blade: a
    # tape_land_width x rad x tape_pad_thickness slab at each end.  They
    # occupy the blade's own last few mm, so the loop height is unchanged.
    #
    # The bottom one is thicker: it runs all the way up to the outside top of
    # the PCB enclosure (see ``bottom_land_thickness``), so the whole base of
    # the antenna is one full-width slab flush with the hub.
    tape_pad_thickness: float = 3.0

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

    # --- Top tape crossover gap -------------------------------------------
    # The two loops are different heights, so at the top the taller blade
    # stands right where the shorter loop's tape has to cross the axis.  Cut
    # a window through the taller blade, sitting on the shorter blade's top
    # face, for that tape to pass through.
    top_tape_gap_width: float = 20.0
    top_tape_gap_height: float = 10.0
    # Material left above the window, at the very top of the taller blade.
    # The window is shortened if it would leave less than this.
    top_tape_gap_min_bridge: float = 3.0

    # --- Splitting into printable sections --------------------------------
    # Tall antennas do not fit on a print bed, so the frame is cut into as
    # few equal horizontal sections as will fit.  Set to None to keep it in
    # one piece.  The pins that stand proud of each cut are counted against
    # this, so a section plus its pins really does fit.
    max_print_height: float | None = 200.0

    # Each joint is a boss on the axis, where both blades cross, carrying two
    # pins.  Below the cut it is a cone that grows to full diameter, so the
    # lower section has no overhang under it; above the cut it is a plain
    # cylinder, whose top face points up and prints fine.
    joint_boss_diameter: float = 18.0
    joint_boss_height_below: float = 12.0
    joint_boss_height_above: float = 8.0

    # The pins are balls sitting on the cut plane: the lower section grows a
    # protruding dome, the upper one a matching dished socket.  A ball reaches
    # only its own radius above the plane, which is what lets pins go out on
    # the arms at all -- a pin standing any real height there would twist
    # straight out the side of the blade, since the blade turns as it rises.
    # They also self-centre, so the joint pulls itself into line as it closes.
    joint_pin_diameter: float = 3.0
    joint_pin_clearance: float = 0.3  # Diametral, on the socket.
    # The pair flanking the axis, inside the boss.  Two is enough to fix the
    # orientation: each blade is symmetric under a half turn about Z, so the
    # only other orientation they allow is the identical one.
    joint_pin_offset: float = 5.0  # From the axis.
    # And more out on both blades' arms, as fractions of the arm's length.
    joint_pin_radius_fractions: tuple[float, ...] = (0.45, 0.8)
    # Which side of a cut carries the balls.  With them on the section above
    # (the default), the section below gets a dish that opens upward -- about
    # the easiest feature there is to print -- while the balls themselves hang
    # under the upper section's mating face, so that section wants flipping on
    # the bed.  Set True to put the balls on the section below instead: they
    # then point up, and the dish becomes a downward-opening cavity whose
    # ceiling has to bridge.
    joint_balls_on_lower_section: bool = False
    # A whisker taken off the outer pole of every ball and socket.  A
    # sphere's mesh collapses to a point at its poles, and the outer one
    # sits on the exposed dome, where it leaves a zero-area triangle and an
    # STL that fails a watertight check.  Cutting the pole away removes it,
    # and gives the printed tip a real flat instead of a point.
    joint_pin_tip_flat: float = 0.15

    # Zip-tie anchors: a hole through each blade above and below the cut, out
    # on the arms.  A tie threaded through both wraps the blade between them,
    # and cinching it pulls the joint shut.  Out on the arms it has the lever
    # arm to resist the bending the joint actually sees.
    joint_tie_hole_z_offset: float = 8.0
    joint_tie_radius_fractions: tuple[float, ...] = (0.65,)

    # --- Feed-throughs into the PCB housing -------------------------------
    # The mast sleeve's wall closes the PCB housing off from the tape runs
    # outside it.  Rather than open it up with a window, each bar gets one
    # small radial hole: the tape ends outside, and a short copper wire
    # passes through to the board's pad.  A 1.5 mm hole is far easier to seal
    # than an open slot, which is the point.
    pcb_wire_hole_diameter: float = 2.5
    # Material left between the frame underside -- the surface the tape runs
    # along -- and the near edge of the hole.  The hole's depth is derived
    # from this, so it stays put against that surface if the diameter changes.
    pcb_wire_hole_edge_margin: float = 1.0

    # --- Mast sleeve -------------------------------------------------------
    # Set to None for no mast sleeve (e.g. antennas too small to straddle a
    # pipe); the hub then reduces to a plate carrying just the PCB bosses.
    mast_pipe_od: float | None = 1.5 * MM_PER_INCH  # 38.1 mm: PVC pipe OD.
    mast_bore_clearance: float = 0.4  # Diametral slip fit.
    mast_sleeve_wall: float = 8.0  # Doubled: this is what grips the mast.
    mast_sleeve_length: float = 40.0  # Hangs below the frame.
    mast_screw_hole_diameter: float = 3.4  # M3 clearance, into the pipe.
    mast_screw_count: int = 6
    # One ring of holes per entry, each measured up from the sleeve's open
    # (bottom) end.  Two rows -- one near each end of the sleeve -- grip the
    # pipe without relying on a single ring to resist tilting.
    mast_screw_rows_z_from_sleeve_end: tuple[float, ...] = (12.0, 28.0)

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
    #
    # A fourth boss is modelled at 315 deg for symmetry: the board has no hole
    # there, so it goes unused and simply backs up the connector quadrant.  It
    # also means the board can be fitted in any of the four rotations.
    pcb_screw_circle_diameter: float = 24.0
    pcb_screw_angles_deg: tuple[float, ...] = (45.0, 135.0, 225.0, 315.0)
    pcb_boss_diameter: float = 7.0
    pcb_standoff_height: float = 6.0  # Gap from frame underside to PCB.
    pcb_screw_hole_diameter: float = 2.7  # M3 thread-forming into plastic.
    pcb_screw_hole_depth: float = 8.0

    def __post_init__(self) -> None:
        """Validate spec parameters."""
        self._validate_blade()
        self._validate_hub()
        self._validate_pcb_mount()
        self._validate_tape_path()
        self._validate_sections()

    @property
    def min_half_length(self) -> float:
        """Half-length of the shorter of the two blades."""
        return min(self.qfh.small_loop.rad, self.qfh.large_loop.rad) / 2.0

    def _validate_blade(self) -> None:
        if self.blade_core_thickness > self.tape_land_width:
            msg = (
                f"Blade core ({self.blade_core_thickness:.1f} mm) is thicker "
                f"than the tape land ({self.tape_land_width:.1f} mm) it is "
                f"supposed to flare out to."
            )
            raise ValueError(msg)

        blade_length = self.min_half_length * 2.0
        if 2.0 * self.tape_land_flare_length >= blade_length:
            msg = (
                f"The two {self.tape_land_flare_length:.1f} mm flares meet in "
                f"the middle of the shortest blade "
                f"({blade_length:.1f} mm long)."
            )
            raise ValueError(msg)

        if self.tie_hole_inset < self.tie_hole_diameter:
            msg = "Zip-tie holes are too close to the blade's end face."
            raise ValueError(msg)

        if self.bottom_land_thickness >= self.tie_hole_bar_z:
            msg = (
                f"Tape pads ({self.bottom_land_thickness:.1f} mm) reach "
                f"the bar zip-tie holes at z={self.tie_hole_bar_z:.1f} mm, "
                f"so a tie could not wrap around them."
            )
            raise ValueError(msg)

    def _validate_hub(self) -> None:
        # Leave the blade ends (the tape lands) clear of the hub.
        if self.hub_radius > self.min_half_length - self.hub_edge_margin:
            msg = (
                f"Hub (r={self.hub_radius:.1f}) crowds the tape lands at "
                f"r={self.min_half_length:.1f}. Use a smaller mast_pipe_od, "
                f"or set mast_pipe_od=None for no mast sleeve."
            )
            raise ValueError(msg)

        if self.top_tape_gap_headroom <= self.top_tape_gap_bridge:
            msg = (
                f"The two loops' tops are only "
                f"{self.top_tape_gap_headroom:.1f} mm apart -- no room for a "
                f"tape crossover gap that still leaves a "
                f"{self.top_tape_gap_bridge:.1f} mm bridge above it."
            )
            raise ValueError(msg)

    def _validate_pcb_mount(self) -> None:
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

    def _validate_sections(self) -> None:
        if self.max_print_height is None:
            return  # Kept in one piece.

        joint_height = (
            self.joint_boss_height_below + self.joint_boss_height_above
        )
        if self.max_print_height <= joint_height * 2:
            msg = (
                f"max_print_height ({self.max_print_height:.1f} mm) is not "
                f"usefully bigger than a joint ({joint_height:.1f} mm)."
            )
            raise ValueError(msg)

        # An arm ball has to fit inside the blade's thin core.
        if self.joint_pin_diameter > self.blade_core_thickness - 1.0:
            msg = (
                f"Joint balls ({self.joint_pin_diameter:.1f} mm) leave under "
                f"0.5 mm of wall in the "
                f"{self.blade_core_thickness:.1f} mm blade core."
            )
            raise ValueError(msg)

        # A cut has to land on plain twisting blade: clear of the hub below,
        # and clear of the crossover window near the top.
        window_lo = min(
            self.qfh.large_loop.height, self.qfh.small_loop.height
        )
        window_hi = window_lo + self.top_tape_gap_height_used
        for cut_z in self.section_cut_heights:
            boss_lo = cut_z - self.joint_boss_height_below
            boss_hi = cut_z + self.joint_boss_height_above
            if boss_lo <= self.hub_plate_thickness:
                msg = (
                    f"A section cut at z={cut_z:.1f} mm puts its joint into "
                    f"the hub. Adjust max_print_height."
                )
                raise ValueError(msg)
            if boss_lo < window_hi and boss_hi > window_lo:
                msg = (
                    f"A section cut at z={cut_z:.1f} mm puts its joint "
                    f"across the tape crossover window "
                    f"(z={window_lo:.1f}..{window_hi:.1f} mm). Adjust "
                    f"max_print_height."
                )
                raise ValueError(msg)

    def _validate_tape_path(self) -> None:
        if self.mast_pipe_od is None:
            return  # No sleeve wall between the tape and the board.

        if self.pcb_wire_hole_edge_margin <= 0.0:
            msg = (
                "Feed-through would break through the frame underside that "
                "the tape runs along."
            )
            raise ValueError(msg)

        hole_bottom = self.pcb_wire_hole_z + self.pcb_wire_hole_diameter / 2.0
        if hole_bottom >= self.pcb_standoff_height:
            msg = (
                f"Feed-through reaches z=-{hole_bottom:.1f} mm, at or below "
                f"the board itself (z=-{self.pcb_standoff_height:.1f} mm); "
                f"it must come through above the board."
            )
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
    def frame_z_min(self) -> float:
        """Lowest z of the whole frame: the mast sleeve's open end."""
        if self.mast_pipe_od is None:
            return 0.0
        return -self.mast_sleeve_length

    @property
    def frame_z_max(self) -> float:
        """Highest z of the whole frame: the taller blade's top."""
        return max(self.qfh.large_loop.height, self.qfh.small_loop.height)

    @property
    def section_count(self) -> int:
        """How many printable sections the frame is cut into."""
        if self.max_print_height is None:
            return 1
        # The balls stand proud of each cut by their radius, so a section's
        # real print height is its slice plus that.  Budget for it rather
        # than discovering it at the printer.
        budget = self.max_print_height - self.joint_pin_diameter / 2.0
        total = self.frame_z_max - self.frame_z_min
        return max(1, math.ceil(total / budget))

    @property
    def section_cut_heights(self) -> tuple[float, ...]:
        """The z heights the frame is cut at, bottom-most first."""
        n = self.section_count
        z_lo, z_hi = self.frame_z_min, self.frame_z_max
        return tuple(
            z_lo + (z_hi - z_lo) * k / n for k in range(1, n)
        )

    @property
    def pcb_wire_hole_z(self) -> float:
        """Depth of the feed-through's axis below the frame underside."""
        return (
            self.pcb_wire_hole_edge_margin + self.pcb_wire_hole_diameter / 2.0
        )

    @property
    def sleeve_mid_radius(self) -> float:
        """Mid-wall radius of the mast sleeve."""
        return (self.mast_bore_radius + self.hub_radius) / 2.0

    @property
    def top_tape_gap_headroom(self) -> float:
        """Height of the taller blade standing above the shorter one."""
        return abs(self.qfh.large_loop.height - self.qfh.small_loop.height)

    @property
    def bottom_land_thickness(self) -> float:
        """Thickness of the full-width land along the bottom of each blade.

        Taken up to the outside top of the PCB enclosure -- the hub plate's
        upper face -- so the base of the antenna is one solid full-width slab
        flush with the hub, rather than a thin land on a thin core.
        """
        return max(self.tape_pad_thickness, self.hub_plate_thickness)

    @property
    def top_tape_gap_bridge(self) -> float:
        """Material kept above the crossover window.

        Never less than the tape pad, so the window's roof is the pad's
        underside and the bridge is the full-width pad rather than thin core.
        """
        return max(self.top_tape_gap_min_bridge, self.tape_pad_thickness)

    @property
    def top_tape_gap_height_used(self) -> float:
        """Gap height actually cut, shortened to keep a top bridge."""
        return min(
            self.top_tape_gap_height,
            self.top_tape_gap_headroom - self.top_tape_gap_bridge,
        )

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


def _blade_section(*, loop_diameter: float, spec: PartSpec) -> bd.Face:
    r"""Build the 2D profile that gets twisted up to make one blade.

    A thin core that flares out at both ends, so material is spent only where
    it earns its place: the full ``tape_land_width`` appears at the outer end
    faces, where the foil tape sticks, and the long middle stays down at
    ``blade_core_thickness``::

        +--_                              _--+   <-- tape land, full width
        |    \____________________________/    |
        |     ____________________________     |   <-- thin core
        +--_/                            \_--+

    Note that the flat top and bottom faces, where the tape turns inboard
    toward the PCB, are only as wide as the core -- tape running along them
    overhangs a little.
    """
    half_len = loop_diameter / 2.0
    t_end = spec.tape_land_width / 2.0
    t_core = spec.blade_core_thickness / 2.0
    x_core = half_len - spec.tape_land_flare_length

    section = bd.Polygon(
        (-half_len, -t_end),
        (-x_core, -t_core),
        (x_core, -t_core),
        (half_len, -t_end),
        (half_len, t_end),
        (x_core, t_core),
        (-x_core, t_core),
        (-half_len, t_end),
        align=None,
    ).face()
    assert section is not None
    return section


def _twisted_segment(
    *,
    section: bd.Face,
    loop_height: float,
    turns: float,
    z_start: float,
    z_end: float,
) -> bd.Solid:
    """Extrude one section over a z range, following the blade's twist.

    The segment starts at whatever angle the blade has reached at ``z_start``
    and twists at the blade's own rate, so segments stack up into a single
    continuous blade.
    """
    height = z_end - z_start
    segment = bd.Solid.extrude_linear_with_rotation(
        section=section,
        center=(0, 0),
        normal=(0, 0, height),
        angle=(360.0 * turns * height / loop_height),
    )
    return segment.rotate(
        bd.Axis.Z, 360.0 * turns * z_start / loop_height
    ).translate((0, 0, z_start))


def _draw_twisted_blade(
    *,
    loop_diameter: float,
    loop_height: float,
    spec: PartSpec,
    crossover_z: float | None = None,
) -> bd.Part | bd.Compound:
    """One loop's twisted blade, with its zip-tie holes.

    The blade's two end faces are the helical tape lands; they sit at
    +/-``loop_diameter``/2, i.e. exactly on the RF design radius.
    """
    turns = spec.qfh.input_spec.turns
    half_len = loop_diameter / 2.0
    pad_t = spec.tape_pad_thickness
    bottom_t = spec.bottom_land_thickness

    # The blade is a stack of three twisted segments rather than one solid
    # with pads laid over it: the flared core in the middle, and a full-width
    # tape land at each end.  Stacking them keeps every joint a plain
    # face-to-face union -- overlapping a pad onto the core instead makes
    # their end faces touch tangentially, which is what turns the boolean
    # degenerate.
    pad_section = bd.Rectangle(loop_diameter, spec.tape_land_width).face()
    assert pad_section is not None
    core_section = _blade_section(loop_diameter=loop_diameter, spec=spec)
    core_top = loop_height - pad_t

    def segment(section: bd.Face, z_start: float, z_end: float) -> bd.Solid:
        return _twisted_segment(
            section=section,
            loop_height=loop_height,
            turns=turns,
            z_start=z_start,
            z_end=z_end,
        )

    core = bd.Part(None) + segment(core_section, bottom_t, core_top)

    # Window for the other loop's tape.  Cut it from the bare core, letting
    # the cutter overshoot the core's top face rather than stopping flush
    # with it, so the top land can then close the window off as its roof
    # without the two ever sharing a coincident cut face.
    if crossover_z is not None:
        z_top = crossover_z + spec.top_tape_gap_height_used
        core -= _top_tape_gap(
            spec=spec,
            blade_height=loop_height,
            z_bottom=crossover_z,
            z_top=(
                core_top + pad_t if z_top >= core_top - 1e-6 else z_top
            ),
        )

    p = bd.Part(None)
    p += segment(pad_section, 0.0, bottom_t)
    p += core
    p += segment(pad_section, core_top, loop_height)

    # Twist is CCW with height: the blade at height z lies along this angle.
    def blade_angle_deg(z: float) -> float:
        return 360.0 * turns * z / loop_height

    hole_len = spec.tape_land_width * 3.0
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


def _top_tape_gap(
    *, spec: PartSpec, blade_height: float, z_bottom: float, z_top: float
) -> bd.Part:
    """Window through the taller blade for the shorter loop's top tape.

    The shorter loop's tape crosses the axis along the top face of its own
    blade, at ``z_bottom``; the taller blade stands right in the way.  The
    window sits directly on that face so the tape runs straight through it.
    """
    z_mid = (z_bottom + z_top) / 2.0

    # Follow the taller blade's twist at the window's mid-height.
    angle_deg = 360.0 * spec.qfh.input_spec.turns * z_mid / blade_height

    return (
        bd.Box(
            length=spec.top_tape_gap_width,
            width=spec.tape_land_width * 3.0,  # Punch clean through.
            height=z_top - z_bottom,
        )
        .rotate(bd.Axis.Z, angle_deg)
        .translate((0.0, 0.0, z_mid))
    )


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
        r_mid = spec.sleeve_mid_radius

        p += bd.Pos(Z=-spec.mast_sleeve_length / 2.0) * bd.Cylinder(
            radius=outer_r, height=spec.mast_sleeve_length
        )

        # Bore for the mast pipe (also the space the PCB hangs in).
        p -= bd.Pos(Z=-spec.mast_sleeve_length / 2.0) * bd.Cylinder(
            radius=bore_r, height=spec.mast_sleeve_length
        )

        # Radial screw holes, for screwing the sleeve to the mast pipe.
        # Each row is staggered half a step from the last so the holes don't
        # all line up along the same four bar lines.
        step = 360.0 / spec.mast_screw_count
        for row, z_from_end in enumerate(
            spec.mast_screw_rows_z_from_sleeve_end
        ):
            z_screw = -(spec.mast_sleeve_length - z_from_end)
            row_offset = 45.0 + (step / 2.0) * (row % 2)
            for i in range(spec.mast_screw_count):
                # Offset off the bar axes, clear of the tape runs above.
                ang = step * i + row_offset
                p -= (
                    bd.Cylinder(
                        radius=spec.mast_screw_hole_diameter / 2.0,
                        height=spec.mast_sleeve_wall * 3.0,
                        # Axis along X, radial once rotated.
                        rotation=(0, 90, 0),
                    )
                    .rotate(bd.Axis.Z, ang)
                    .translate((*_polar(r_mid, ang), z_screw))
                )

        # One small feed-through per bar: the tape stops outside the
        # sleeve and a short copper wire runs through to the board's pad.
        for ang in BAR_ANGLES_DEG:
            p -= (
                bd.Cylinder(
                    radius=spec.pcb_wire_hole_diameter / 2.0,
                    height=spec.mast_sleeve_wall * 3.0,
                    rotation=(0, 90, 0),  # Axis along X, radial once rotated.
                )
                .rotate(bd.Axis.Z, ang)
                .translate((*_polar(r_mid, ang), -spec.pcb_wire_hole_z))
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
    large_h = spec.qfh.large_loop.height
    small_h = spec.qfh.small_loop.height

    # The taller blade blocks the shorter loop's tape where it crosses the
    # axis along its top face, so it gets a window at the shorter one's top.
    large_is_taller = large_h >= small_h

    # Large loop blade: bottom bar along +/-X.
    large_blade = _draw_twisted_blade(
        loop_diameter=spec.qfh.large_loop.rad,
        loop_height=large_h,
        spec=spec,
        crossover_z=small_h if large_is_taller else None,
    )

    # Small loop blade: a quarter turn around, bottom bar along +/-Y.
    small_blade = _draw_twisted_blade(
        loop_diameter=spec.qfh.small_loop.rad,
        loop_height=small_h,
        spec=spec,
        crossover_z=None if large_is_taller else large_h,
    ).rotate(bd.Axis.Z, 90)

    used = spec.top_tape_gap_height_used
    if used < spec.top_tape_gap_height:
        logger.warning(
            "Top tape gap shortened to {:.1f} mm (asked {:.1f}); the loops' "
            "tops are only {:.1f} mm apart and a {:.1f} mm bridge is kept.",
            used,
            spec.top_tape_gap_height,
            spec.top_tape_gap_headroom,
            spec.top_tape_gap_bridge,
        )

    p = bd.Part(None)
    p += large_blade
    p += small_blade

    # The hub lives entirely at and below z=0, while the blades run upward
    # from z=0, so fusing it in last cannot backfill any of its holes.
    p += _draw_hub(spec)

    return p


# ---------------------------------------------------------------------------
# Splitting into printable sections
# ---------------------------------------------------------------------------


def _joint_boss(*, spec: PartSpec, cut_z: float) -> bd.Part | bd.Compound:
    """Boss straddling one cut, on the axis where both blades cross.

    Below the cut it is a cone growing to full diameter, so the lower section
    has no overhang beneath it; above the cut it is a plain cylinder, whose
    top face points upward and so prints fine on the upper section.
    """
    r_big = spec.joint_boss_diameter / 2.0
    r_small = spec.blade_core_thickness / 2.0

    p = bd.Part(None)
    p += bd.Pos(Z=cut_z - spec.joint_boss_height_below / 2.0) * bd.Cone(
        bottom_radius=r_small,
        top_radius=r_big,
        height=spec.joint_boss_height_below,
    )
    p += bd.Pos(Z=cut_z + spec.joint_boss_height_above / 2.0) * bd.Cylinder(
        radius=r_big, height=spec.joint_boss_height_above
    )
    return p


def _joint_pin_positions(
    *, spec: PartSpec, cut_z: float
) -> list[tuple[float, float]]:
    """Where the balls sit on one cut plane.

    A pair flanking the axis inside the boss, plus more out on both blades'
    arms, each at that blade's own angle where the cut crosses it.
    """
    turns = spec.qfh.input_spec.turns
    positions: list[tuple[float, float]] = []

    # The central pair, laid along the large blade.
    axis_angle_deg = 360.0 * turns * cut_z / spec.qfh.large_loop.height
    positions += [
        _polar(sign * spec.joint_pin_offset, axis_angle_deg)
        for sign in (1.0, -1.0)
    ]

    # Out on the arms, on both blades.
    for loop, blade_rot in (
        (spec.qfh.large_loop, 0.0),
        (spec.qfh.small_loop, 90.0),
    ):
        if cut_z > loop.height:
            continue  # Past the end of this blade.
        angle_deg = 360.0 * turns * cut_z / loop.height + blade_rot
        positions += [
            _polar(loop.rad / 2.0 * frac, angle_deg + end_sign)
            for frac in spec.joint_pin_radius_fractions
            for end_sign in (0.0, 180.0)
        ]

    return positions


def _joint_pins(
    *, spec: PartSpec, cut_z: float, socket: bool
) -> bd.Part | bd.Compound:
    """Build the balls at one cut, or the sockets they seat into.

    Each ball is centred on the cut plane, so half of it stands proud of one
    section and the matching half is dished out of the other.  Ball and socket
    always use the same half, the one ``joint_balls_on_lower_section`` puts
    the protruding dome in.
    """
    diameter = spec.joint_pin_diameter + (
        spec.joint_pin_clearance if socket else 0.0
    )
    radius = diameter / 2.0
    sign = 1.0 if spec.joint_balls_on_lower_section else -1.0

    # Trim the outer pole away (see joint_pin_tip_flat).  The inner one is
    # buried in the section the feature belongs to, so it needs no such
    # treatment.  Taking the same slice off ball and socket alike leaves the
    # clearance between them untouched.
    cap = bd.Pos(
        Z=sign * (2.0 * radius - spec.joint_pin_tip_flat)
    ) * bd.Box(4.0 * radius, 4.0 * radius, 2.0 * radius)

    p = bd.Part(None)
    for x, y in _joint_pin_positions(spec=spec, cut_z=cut_z):
        p += bd.Pos(x, y, cut_z) * (bd.Sphere(radius=radius) - cap)
    return p


def _joint_tie_holes(
    *, spec: PartSpec, cut_z: float
) -> bd.Part | bd.Compound:
    """Zip-tie anchor holes flanking one cut, out on both blades' arms.

    A tie threaded through the hole below the cut and the one above wraps the
    blade between them, so cinching it pulls the joint shut.
    """
    turns = spec.qfh.input_spec.turns
    holes = bd.Part(None)

    for loop, blade_rot in (
        (spec.qfh.large_loop, 0.0),
        (spec.qfh.small_loop, 90.0),
    ):
        half_len = loop.rad / 2.0
        for dz in (
            -spec.joint_tie_hole_z_offset,
            spec.joint_tie_hole_z_offset,
        ):
            z = cut_z + dz
            if not 0.0 <= z <= loop.height:
                continue  # Past the end of this blade.
            angle_deg = 360.0 * turns * z / loop.height + blade_rot
            for frac in spec.joint_tie_radius_fractions:
                for end_sign in (0.0, 180.0):
                    holes += _tie_hole(
                        radius=half_len * frac,
                        angle_deg=angle_deg + end_sign,
                        z=z,
                        diameter=spec.tie_hole_diameter,
                        length=spec.tape_land_width * 3.0,
                    )
    return holes


def qfh_antenna_frame_sections(spec: PartSpec) -> list[bd.Part | bd.Compound]:
    """Build the frame, split into printable sections.

    Returns the sections in assembled position, bottom-most first.  A frame
    that already fits comes back as a single-item list.

    Joint features are applied to the whole frame before it is cut, so each
    cut divides them between the two sections it makes.  Only the pins go on
    afterwards, since they stand proud of the cut plane and would otherwise
    be sliced straight off.
    """
    cut_zs = spec.section_cut_heights
    frame = qfh_antenna_frame(spec)
    if not cut_zs:
        return [frame]

    for cut_z in cut_zs:
        frame += _joint_boss(spec=spec, cut_z=cut_z)
        frame -= _joint_tie_holes(spec=spec, cut_z=cut_z)

    # Slice into sections.
    bounds = [spec.frame_z_min, *cut_zs, spec.frame_z_max]
    span = 4.0 * spec.qfh.large_loop.rad  # Comfortably wider than the frame.
    sections: list[bd.Part | bd.Compound] = []
    for z_lo, z_hi in itertools.pairwise(bounds):
        slab = bd.Pos(Z=(z_lo + z_hi) / 2.0) * bd.Box(
            length=span, width=span, height=z_hi - z_lo
        )
        sections.append(frame & slab)

    # Balls straddle the cut plane, so unlike the boss and the tie holes they
    # cannot be applied before the split: each half belongs to one section
    # only.
    ball_below = spec.joint_balls_on_lower_section
    for i, cut_z in enumerate(cut_zs):
        male, female = (i, i + 1) if ball_below else (i + 1, i)
        sections[male] += _joint_pins(spec=spec, cut_z=cut_z, socket=False)
        sections[female] -= _joint_pins(spec=spec, cut_z=cut_z, socket=True)

    logger.info(
        "Split into {} sections at z = {}; tallest is {:.1f} mm "
        "(max print height {:.1f} mm)",
        len(sections),
        ", ".join(f"{z:.1f}" for z in cut_zs),
        max(s.bounding_box().size.Z for s in sections),
        spec.max_print_height,
    )
    if not ball_below:
        logger.info(
            "Balls are on the section above each cut, so they hang below "
            "that section's mating face: print sections 2..{} flipped, or "
            "the whole mating face is a {:.2f} mm overhang. Section 1 has "
            "only an upward-opening dish and prints as modelled.",
            len(sections),
            spec.joint_pin_diameter / 2.0 - spec.joint_pin_tip_flat,
        )
    return sections


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    """Generate the QFH antenna support structure and export."""
    specs = {
        # 436 MHz: ~100 mm across, so it straddles a 1.5 in PVC mast.  At
        # 267 mm tall it does not fit a print bed, so it comes out in two
        # sections that pin and zip-tie together.
        "QFH_Antenna_436_MHz": PartSpec(
            qfh=calculate_qfh(
                QfhInputSpec(
                    frequency_hz=436.0e6,
                    wire_diameter=1.5,
                    wire_bending_radius=3.0,
                )
            )
        ),
        # 913 MHz: only ~46 mm across -- narrower than a 1.5 in pipe -- so
        # there is no mast sleeve; it just carries the balun PCB.  At 109 mm
        # tall it prints in one piece.
        "QFH_Antenna_913_MHz": PartSpec(
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
        ),
    }

    parts: dict[str, bd.Part | bd.Compound | bd.Solid] = {}
    for name, spec in specs.items():
        sections = qfh_antenna_frame_sections(spec)
        if len(sections) == 1:
            parts[name] = sections[0]
            continue
        for i, section in enumerate(sections, start=1):
            parts[f"{name}_section_{i}_of_{len(sections)}"] = section

    # Show just the 436 MHz build, its sections stood side by side on the
    # bed.  They are laid out rather than left stacked so both are visible at
    # once instead of forming one 267 mm tower, and shown in a single call --
    # show() replaces what is displayed, so calling it per part would leave
    # only the last one on screen.
    shown = {
        name: part
        for name, part in parts.items()
        if name.startswith(SHOWN_PART_PREFIX)
    }
    pitch = (
        max(part.bounding_box().size.X for part in shown.values())
        + SHOWN_PART_GAP
    )
    logger.info("Showing {}", ", ".join(shown))
    show(
        *(
            part.translate(
                (i * pitch, 0.0, -part.bounding_box().min.Z)
            )
            for i, part in enumerate(shown.values())
        )
    )

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

        # Drop each section onto the print bed rather than leaving it at its
        # assembled height.
        printable = part.translate((0, 0, -part.bounding_box().min.Z))

        bd.export_stl(printable, str(export_folder / f"{name}.stl"))
        bd.export_step(printable, str(export_folder / f"{name}.step"))
        logger.info('Exported "{}" to {}', name, export_folder)


if __name__ == "__main__":
    main()
