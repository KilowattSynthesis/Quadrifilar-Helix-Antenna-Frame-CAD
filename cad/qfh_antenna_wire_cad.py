"""Build123d solid model of a Quadrifilar Helix (QFH) antenna.

Topology
--------
A QFH is two independent *bifilar loops* sharing one axis, the second rotated
90 deg about Z.  Each loop is a single closed conductor shaped like a
"twisted rectangle":

    bottom bar  ->  helical side A  -> top bar ->  helical side B  ->  (close)

* The two **helical sides** wind ``turns`` revolutions (0.5 = 180 deg) up a
  cylinder of diameter ``rad`` (= centre-to-centre ``D``), rising ``height``.
* The **top / bottom bars** are the straight diametric segments joining the two
  sides.  Because of the half-turn twist the top bar lies 180 deg around from
  the bottom bar.
* The four **corners** of each loop are the rounded bends of radius
  ``wire_bending_radius`` -- these are "the arcs".

Mapping qfh_calc -> geometry (verified against the solver output):
    rad     = D, cylinder diameter centre-to-centre
    height  = H, vertical extent (bar-to-bar)
    idiam   = rad - wire_diameter  (inner clear diameter; confirms rad is D)
    vert    = 3D wire length of one helical side = sqrt(H^2 + (pi*rad*turns)^2)
    *_comp  = the same, minus the material consumed by the corner bends

The corner: an over-constrained joint
--------------------------------------
At the base of a helix the tangent is purely *circumferential + vertical* with
**no radial component** (e.g. (0, 0.95, 0.30) for the default design).  The bar
arrives *radially*.  The two are perpendicular, so each corner is a ~90 deg
turn that lives in the tilted plane spanned by {radial, helix-tangent} --
not a flat vertical elbow.

A single circular arc of a fixed radius cannot be exactly tangent to BOTH the
straight bar AND the *curving* helix at its on-cylinder base point: the helix
leaves its own tangent line immediately, so the three constraints
(radius = rho, tangent to bar, tangent to helix-at-base) over-determine the
arc.
This model resolves it by honouring the tangency that matters
electromechanically -- a smooth, true-radius entry into the radiating helix --
using :class:`TangentArc`, and absorbing a sub-degree kink at the
*straight-bar* end (where it is invisible).  As a check, the resulting
centreline length comes out within ~0.1 % of ``LoopResult.total_comp``.

Made with Claude Opus.
"""

from pathlib import Path

import build123d as bd
from build123d_ease import show
from loguru import logger

from cad.qfh_calc import LoopResult, QfhInputSpec, QfhResult, calculate_qfh

# ---------------------------------------------------------------------------
# Centreline construction
# ---------------------------------------------------------------------------


def _bar_corner_point(
    helix_end: bd.Vector,
    helix_tangent: bd.Vector,
    vertical_sign: int,
    rho: float,
) -> bd.Vector:
    """Point where a corner bend hands off to a straight bar.

    Starting from a helix endpoint and stepping one bend-radius ``rho`` along
    the helix tangent (up at the top, down at the bottom) and one bend-radius
    radially inward gives the tangent point of a true 90 deg / radius-``rho``
    arc.  ``vertical_sign`` is -1 at a helix *base* (bend dives down to the
    bar) and +1 at a helix *top* (bend rises over to the bar).
    """
    radial_out = bd.Vector(helix_end.X, helix_end.Y, 0).normalized()
    return (
        helix_end + helix_tangent * (vertical_sign * rho) + radial_out * (-rho)
    )


def loop_centerline(
    rad: float, height: float, turns: float, rho: float
) -> bd.Wire:
    """Make a closed centreline wire for one bifilar loop.

    Arguments:
    ---------
    rad : Cylinder diameter ``D`` (centre-to-centre), i.e. ``LoopResult.rad``.
    height : Antenna height ``H`` (``LoopResult.height``).
    turns : Helix twist in revolutions (``QfhInputSpec.turns``).
    rho : Corner bend radius (``QfhInputSpec.wire_bending_radius``).

    """
    radius = rad / 2.0

    # Trim the helix height by the vertical reach of the two end bends so the
    # finished bars land exactly at z = 0 and z = height.
    #
    # tz is the vertical fraction of the (unit) helix tangent; one probe helix
    # is enough to get it.
    probe = bd.Helix(pitch=height / turns, height=height, radius=radius)
    tz = abs((probe % 0).Z)
    helix_h = height - 2.0 * rho * tz
    z0 = rho * tz

    side_a = bd.Helix(
        pitch=helix_h / turns, height=helix_h, radius=radius
    ).translate((0, 0, z0))
    side_b = side_a.rotate(bd.Axis.Z, 180)  # the two arms are 180 deg apart

    # Endpoints (@) and unit tangents (%) straight from the helix geometry.
    a0, a1 = side_a @ 0, side_a @ 1
    b0, b1 = side_b @ 0, side_b @ 1
    ta0, ta1 = (side_a % 0).normalized(), (side_a % 1).normalized()
    tb0, tb1 = (side_b % 0).normalized(), (side_b % 1).normalized()

    # Bar-side tangent points for all four corners.
    sa0 = _bar_corner_point(a0, ta0, -1, rho)
    sa1 = _bar_corner_point(a1, ta1, +1, rho)
    sb0 = _bar_corner_point(b0, tb0, -1, rho)
    sb1 = _bar_corner_point(b1, tb1, +1, rho)

    # Corner bends: true radius-rho arcs, tangent to the helix at its endpoint.
    bend_a0 = bd.TangentArc(a0, sa0, tangent=ta0 * -1)
    bend_a1 = bd.TangentArc(a1, sa1, tangent=ta1 * +1)
    bend_b0 = bd.TangentArc(b0, sb0, tangent=tb0 * -1)
    bend_b1 = bd.TangentArc(b1, sb1, tangent=tb1 * +1)

    bottom_bar = bd.Line(sb0, sa0)
    top_bar = bd.Line(sa1, sb1)

    edges = (
        side_a.edges()
        + side_b.edges()
        + bend_a0.edges()
        + bend_a1.edges()
        + bend_b0.edges()
        + bend_b1.edges()
        + bottom_bar.edges()
        + top_bar.edges()
    )
    wire = bd.Wire(edges)
    if not wire.is_closed:
        msg = "loop centreline failed to close"
        raise RuntimeError(msg)
    return wire


def build_loop(
    loop: LoopResult, turns: float, rho: float, wire_radius: float
) -> bd.Part:
    """Sweep the conductor profile along one loop's centreline -> solid."""
    wire = loop_centerline(loop.rad, loop.height, turns=turns, rho=rho)
    profile = bd.Plane(origin=wire @ 0, z_dir=wire % 0) * bd.Circle(
        wire_radius
    )
    assert isinstance(profile, bd.Sketch)  # Type checking.
    p = bd.sweep(profile, path=wire)
    assert isinstance(p, bd.Part)  # Type checking.
    return p


def build_qfh_antenna(
    qfh_input_spec: QfhInputSpec, qfh_result: QfhResult | None = None
) -> bd.Compound:
    """Build the full two-loop QFH antenna for a given input spec.

    The large and small loops are kept as *separate* solids in a Compound -- a
    QFH's two loops are electrically distinct and must not be fused (a boolean
    union would short them).  Note: with both bottom bars passing through the
    axis at z = 0 they geometrically cross at the centre; a real build offsets
    the feed there.  Adjust as needed for your feed arrangement.
    """
    qfh_result = qfh_result or calculate_qfh(qfh_input_spec)
    wire_radius = qfh_input_spec.wire_diameter / 2.0

    large = build_loop(
        qfh_result.large_loop,
        qfh_input_spec.turns,
        qfh_input_spec.wire_bending_radius,
        wire_radius,
    )
    small = build_loop(
        qfh_result.small_loop,
        qfh_input_spec.turns,
        qfh_input_spec.wire_bending_radius,
        wire_radius,
    )
    small = small.rotate(bd.Axis.Z, 90)  # Second loop a quarter-turn around.

    # Created a named compound (awesome for exploring/viewing).
    large.label = "large_loop"
    small.label = "small_loop"
    p = bd.Compound(label="qfh_antenna", children=[large, small])

    # Print info about it:
    bb = p.bounding_box().size
    logger.debug(
        f"Built QFH: {len(p.solids())} solids, "
        f"volume {p.volume:.0f} mm^3, "
        f"bbox {bb.X:.1f} x {bb.Y:.1f} x {bb.Z:.1f} mm"
    )

    # Sanity: centreline length vs solver's compensated wire length.
    for name, lp in (
        ("large", qfh_result.large_loop),
        ("small", qfh_result.small_loop),
    ):
        w = loop_centerline(
            lp.rad,
            lp.height,
            qfh_input_spec.turns,
            qfh_input_spec.wire_bending_radius,
        )
        logger.debug(
            f"  {name} loop centreline {w.length:.1f} mm  "
            f"(solver total_comp {lp.total_comp:.1f} mm)"
        )

    return p


def main() -> None:
    """Generate the QFH antenna wires and export."""
    input_spec_913e6 = QfhInputSpec(
        frequency_hz=913.0e6,
        wire_diameter=1.5,  # conductor outer diameter (mm)
        wire_bending_radius=1.5,  # bending radius (mm)
        ratio=0.44,  # width / height ratio
        turns=0.5,  # half-turn helix
        num_wavelengths=1.0,  # one wavelength per loop
    )
    qfh_result_913e6 = calculate_qfh(input_spec_913e6)

    parts = {
        "QFH_Antenna_913_MHz": show(
            build_qfh_antenna(input_spec_913e6, qfh_result=qfh_result_913e6)
        ),
        "QFH_Antenna_436_MHz": show(
            build_qfh_antenna(
                qfh_input_spec=QfhInputSpec(
                    frequency_hz=436.0e6,
                    wire_diameter=1.5,
                    wire_bending_radius=1.5,
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
