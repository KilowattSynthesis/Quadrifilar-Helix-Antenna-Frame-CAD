# Quadrifilar-Helix-Antenna-Frame-CAD

A 3d-printable QFH antenna frame, made with Build123d. Customizable for any frequency.

## Goals

* Customizable, reliable frame for building a QFH antenna at home.
* Focus on ease-of-sourcing with tools like 3D printers and JLCPCB PCB fabriaction.
* Fully self-contained. All calculations are completed in this repo. No reliance on external geometry calculators.

## Design

The frame is two **twisted blades** -- one per bifilar loop -- crossing at 90
deg and printed as a single part.

* **The conductor is foil tape**, stuck to the flat outer end faces of the
  blades. There are no wire channels to thread. Each blade end face sits
  exactly on the RF design radius, so the tape centreline lands where the
  solver says the conductor should be.
* **Zip-tie holes** run through each blade near its ends and near the top and
  bottom bars. A tie goes through a hole, wraps around the nearby edge *over
  the tape*, and cinches -- mechanical backup for the adhesive.
* **Top crossover window** through the taller blade, sitting on the shorter
  blade's top face, so the shorter loop's tape can cross the axis. It is
  shortened automatically if the two loop heights are too close to leave a
  `top_tape_gap_min_bridge` (3 mm) bridge above it -- with the stock 0.44
  ratio there is only 11.2 mm of headroom at 436 MHz and 5.3 mm at 913 MHz,
  so the 20 mm wide window comes out 8.2 mm and 2.3 mm tall respectively.
  Foil tape is thin, so a short window is still enough to pass it.
* **Mast sleeve** at the bottom, bored for a 1.5 in (38.1 mm) OD PVC pipe,
  with radial screw holes for fastening into the pipe.
* **Balun PCB mount**: three M3 bosses on a 24 mm bolt circle, matching the
  [QFHBAL01](https://github.com/ODZ-UJF-AV-CR/QFHBAL01) board. The board hangs
  under the hub inside the sleeve bore, and doubles as the mast pipe's
  insertion depth stop.

Tape path, per loop:

```text
PCB pad -> underside of the bottom bar -> around the bottom corner
        -> up the outer end face (helical, 180 deg) -> over the top face
        -> down the opposite end face -> underside of the other bottom bar
        -> PCB pad
```

### Assembly

1. Print the frame upright (as modelled). No supports needed.
2. Run foil tape along all four tape paths above. Keep it centred on the
   narrow end faces.
3. Pass the shorter loop's top tape through the crossover window in the
   taller blade.
4. Thread a zip tie through each hole, over the tape, and cinch.
5. Fit the balun PCB **rotated 45 deg to the frame's bars**, so its X-shaped
   pads line up with the four tape ends. Three M3 screws from below thread
   into the bosses. Solder the tape ends to the pads.
6. Slide the frame onto the PVC mast until the pipe meets the PCB, then drive
   the three mast screws.

Everything is parametric -- see `PartSpec` in `cad/qfh_frame_cad.py` for tape
land width, tie hole spacing, mast diameter, standoff height, and so on. Set
`mast_pipe_od=None` for antennas too small to straddle a pipe (the bundled
913 MHz build is only ~46 mm across, so it has no sleeve).

## See Also

* Geometry Calculator: https://jcoppens.com/ant/qfh/calc.en.php
* https://github.com/cernohorsky/QFH-Antenna-868MHz/blob/master/QFH-Antenna-868MHz.scad
* Balun PCB: https://github.com/ODZ-UJF-AV-CR/QFHBAL01

### Related, but less helpful

* https://github.com/cernohorsky/QFH-Antenna-868MHz
* https://www.thingiverse.com/thing:634205
* https://usradioguy.com/wp-content/uploads/2020/05/20200307-How-To-Build-A-QFH.pdf
* https://network.satnogs.org/stations/4704/
