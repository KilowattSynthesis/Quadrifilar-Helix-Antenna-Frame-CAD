from dataclasses import dataclass
from pathlib import Path

import build123d as bd
from build123d_ease import show
from loguru import logger

from cad.qfh_antenna_wire_cad import build_qfh_antenna
from cad.qfh_calc import QfhInputSpec, QfhResult, calculate_qfh


@dataclass
class PartSpec:
    """Mechanical specification for the QFH antenna support structure.

    All lengths are in millimeters (mm).

    The RF geometry is read from the ``qfh`` argument, which is the output of
    the QFH geometry calculation step.

    """

    qfh: QfhResult

    # Mechanical / print parameters:
    wire_diameter: float = 2.0  # Conductor outer diameter.

    channel_thickness_on_sides: float = 2.0

    def __post_init__(self) -> None:
        """Validate spec parameters."""


def _draw_loop_wire_channel(
    *,
    loop_diameter: float,
    loop_height: float,
    wire_diameter: float,
    channel_thickness_on_sides: float,
    turns: float,
) -> bd.Part | bd.Compound:

    # Create smaller loop (in X axis).
    small_loop_sk = bd.Sketch()
    # Line connecting the circles.
    small_loop_sk += bd.Rectangle(
        loop_diameter + wire_diameter,
        wire_diameter + channel_thickness_on_sides * 2,
    )
    for x in (-1, 1):
        # Round wire.
        small_loop_sk -= bd.Pos(X=x * (loop_diameter / 2)) * bd.Circle(
            radius=wire_diameter / 2
        )

        # Square channel extension to surface.
        small_loop_sk -= bd.Pos(
            X=x * (loop_diameter / 2 + wire_diameter / 2)
        ) * bd.Rectangle(wire_diameter, wire_diameter)

    small_loop_sk_face = small_loop_sk.face()
    assert small_loop_sk_face is not None

    p = bd.Part(None)
    p += bd.Solid.extrude_linear_with_rotation(
        section=small_loop_sk_face,
        center=(0, 0),
        normal=(0, 0, loop_height),  # Distance.
        angle=(360 * turns),
    )

    return p


def qfh_antenna_frame(spec: PartSpec) -> bd.Compound:
    """Create the QFH antenna support structure."""
    p = bd.Part(None)

    # Small loop channels.
    p += _draw_loop_wire_channel(
        loop_diameter=spec.qfh.small_loop.rad,
        loop_height=spec.qfh.small_loop.height,
        wire_diameter=spec.wire_diameter,
        channel_thickness_on_sides=spec.channel_thickness_on_sides,
        turns=spec.qfh.input_spec.turns,
    ).rotate(bd.Axis.Z, 90)
    # .translate((0, 0, spec.qfh.large_loop.height-spec.qfh.small_loop.height))

    # Large loop channels.
    p += _draw_loop_wire_channel(
        loop_diameter=spec.qfh.large_loop.rad,
        loop_height=spec.qfh.large_loop.height,
        wire_diameter=spec.wire_diameter,
        channel_thickness_on_sides=spec.channel_thickness_on_sides,
        turns=spec.qfh.input_spec.turns,
    )

    alignment_spin_deg_small = spec.qfh.small_loop.frame_alignment_spin_deg(
        turns=spec.qfh.input_spec.turns,
        wire_bending_radius=spec.qfh.input_spec.wire_bending_radius,
    )
    alignment_spin_deg_large = spec.qfh.large_loop.frame_alignment_spin_deg(
        turns=spec.qfh.input_spec.turns,
        wire_bending_radius=spec.qfh.input_spec.wire_bending_radius,
    )
    # Surprisingly, the angles are slightly different for the center connector.
    alignment_spin_deg_avg = (
        alignment_spin_deg_small + alignment_spin_deg_large
    ) / 2

    antenna_wire = build_qfh_antenna(
        qfh_input_spec=spec.qfh.input_spec,
        qfh_result=spec.qfh,
    ).rotate(axis=bd.Axis.Z, angle=alignment_spin_deg_avg)

    p -= antenna_wire

    return p


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    """Generate the QFH antenna support structure and export."""
    input_spec = QfhInputSpec(
        frequency_hz=913.0e6,  # MHz
        wire_diameter=1.5,  # conductor outer diameter (mm)
        wire_bending_radius=3.0,  # bending radius (mm)
        ratio=0.44,  # width / height ratio
        turns=0.5,  # half-turn helix
        num_wavelengths=1.0,  # one wavelength per loop
    )
    qfh_result = calculate_qfh(input_spec)

    parts = {
        "QFH_Antenna_913_MHz": show(
            qfh_antenna_frame(PartSpec(qfh=qfh_result))
        ),
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
