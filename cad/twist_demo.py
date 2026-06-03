import build123d as bd
from build123d_ease import show

if __name__ == "__main__":
    hex_shape = bd.RegularPolygon(radius=1, side_count=6)
    hex_shape_face = hex_shape.face()
    assert hex_shape_face

    twist_extrude = bd.Solid.extrude_linear_with_rotation(
        section=hex_shape_face,
        center=(0, 0),
        normal=(0, 0, 5),  # Distance to extrude.
        angle=(360 / 5),
    )

    show(twist_extrude)
