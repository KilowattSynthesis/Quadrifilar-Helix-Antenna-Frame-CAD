# Quadrifilar-Helix-Antenna-Frame-CAD

A 3d-printable QFH antenna frame, made with Build123d. Customizable for any frequency.

## Goals

* Customizable, reliable frame for building a QFH antenna at home.
* Focus on ease-of-sourcing with tools like 3D printers and JLCPCB PCB fabriaction.
* Fully self-contained. All calculations are completed in this repo. No reliance on external geometry calculators.

## Variants/Components

### Main

* `qfh_tape_frame_cad.py`: Main CAD model using foil tape on the surface of this frame. Works well.

### Work-in-Progress and/or Ancillary

* `qfh_antenna_wire_cad.py`: Render of the ideal copper wire/pipe conductor itself, without a frame.
* `qfh_frame_cad.py`: Work-in-progress variant using copper wire/pipe as the conductor.
* `qfh_calc.py`: Calculations for antenna geometry. No CAD result.

## Rough Build Guide

1. Print antenna.
2. Assemble antenna with copper foil.
3. Figure out which orientation the [QFHBAL01 PCB](https://github.com/ODZ-UJF-AV-CR/QFHBAL01) should go in, and clip off the standoff that physically interferes with the PCB's SMA connector.
    * Follow the table in its README for the Standard/Anti-standard config.
    * Recommendation: Populate transformer in galvanic isolation mode to protect circuitry against static buildup.
4. Install the [QFHBAL01 PCB](https://github.com/ODZ-UJF-AV-CR/QFHBAL01) in antenna. Solder wires from foil to PCB.
5. Screw/glue antenna onto a 1.75" OD pipe.

## See Also

* Geometry Calculator: https://jcoppens.com/ant/qfh/calc.en.php
* https://github.com/cernohorsky/QFH-Antenna-868MHz/blob/master/QFH-Antenna-868MHz.scad
* Balun PCB: https://github.com/ODZ-UJF-AV-CR/QFHBAL01

### Related, but less helpful

* https://github.com/cernohorsky/QFH-Antenna-868MHz
* https://www.thingiverse.com/thing:634205
* https://usradioguy.com/wp-content/uploads/2020/05/20200307-How-To-Build-A-QFH.pdf
* https://network.satnogs.org/stations/4704/

## AI Note

This project is made heavily with AI, especially for the complex geometries around the helical elements.
