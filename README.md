# Quadrifilar-Helix-Antenna-Frame-CAD

A 3d-printable QFH antenna frame, made with Build123d. Customizable for any frequency.

## Goals

* Customizable, reliable frame for building a QFH antenna at home.
* Focus on ease-of-sourcing with tools like 3D printers and JLCPCB PCB fabriaction.
* Fully self-contained. All calculations are completed in this repo. No reliance on external geometry calculators.

## Design

The frame is two **twisted blades** -- one per bifilar loop -- crossing at 90
deg and printed as a single part.

Each blade's section is a thin 5 mm core that flares out over the last 10 mm
at each end to a full 10 mm, so the tape gets a wide land to stick to without
paying for that thickness along the whole blade -- about 45% less section area
than a plain 10 mm slab. The flare only widens the blade radially, so the top
and bottom of every blade is a full-width land across its whole length,
giving the tape a 10 mm surface where it turns inboard. The bottom land is
the thicker of the two -- it runs up to the outside top of the PCB enclosure
(8.1 mm: the 6 mm hub plate, raised by the 2.1 mm `hub_rise` below it), so
the base of the antenna is one solid full-width slab flush with the hub,
which is where the frame carries the most load.

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
  with radial screw holes for fastening into the pipe. The wall is 8 mm --
  this is the joint the whole antenna hangs off, so it is deliberately heavy.
* **Balun PCB mount**: four M3 bosses on a 24 mm bolt circle, matching the
  [QFHBAL01](https://github.com/ODZ-UJF-AV-CR/QFHBAL01) board. That board has
  three holes and an RF connector in the fourth quadrant, so one boss goes
  unused -- it just backs up the connector side, and lets the board be fitted
  in any of the four rotations. The board hangs under the hub inside the
  sleeve bore, and doubles as the mast pipe's insertion depth stop.
* **Wire feed-throughs**: the sleeve wall closes the PCB housing off from
  the tape runs outside it, so each of the four bars gets a single hole,
  bored to the feed wire plus 0.7 mm (2.2 mm for the stock 1.5 mm wire). The
  tape ends outside and a short copper wire passes through to the board's
  pad -- much easier to weatherproof than an open slot.
* **Straight feed wires**, which is what sets the hub's height. The tape is
  stuck to the *outside* of the frame, so its exposed face -- the one the
  wire solders to -- sits `tape_thickness` (0.9 mm) below the bar's
  underside. A `feed_wire_diameter` (1.5 mm) wire lying on that face has its
  axis 1.65 mm below the underside, and it keeps that height the whole way:
  the feed-through is bored concentric with it, and the board's top face is
  set one wire radius lower again, at 2.4 mm below the underside, so the
  wire arrives lying flat on the board. Nothing to bend, nothing to hold
  while the iron is on it -- a straight wire soldered flat at both ends.
* **Raised hub**: the board's height is fixed by that wire run, so
  `pcb_standoff_height` (4.5 mm) lifts the *hub* rather than lowering the
  board. The whole hub -- plate, sleeve, bore and bosses -- rises
  `hub_rise` = 2.1 mm up into the frame. That leaves the standoff free to be
  whatever the tallest part on top of the board needs, and takes the mast
  mounting up inside the antenna instead of leaving it all hanging beneath.
  The bore rises with it and takes a 2.1 mm bite out of the bottom of the
  blades where they cross it -- under the hub, well inboard of any tape.

* **Printable sections**: anything taller than `max_print_height` (200 mm by
  default) is cut into as few equal horizontal sections as will fit, counting
  the balls that stand proud of each cut. Locating pins are 3 mm **balls**
  sitting on the cut plane -- half proud of one section, half dished into the
  other -- placed both on the axis and out on both blades' arms.
  Balls self-centre as the joint closes, and reaching only their own radius
  above the plane is what lets pins go out on the arms at all: a pin standing
  any real height there would twist straight out the side of the blade. The
  axial pair sits in a boss that also ties the two blades together, a cone
  below the cut and a cylinder above it so both sections print without an
  overhang under the joint. Zip-tie anchors flank each cut further out. Set `max_print_height=None` to keep the frame in one
  piece. The bundled 436 MHz build is 267 mm tall and comes out in two
  sections; the 913 MHz one is 109 mm and stays whole.

Tape path, per loop:

```text
PCB pad -> wire through the feed-through -> underside of the bottom bar
        -> around the bottom corner
        -> up the outer end face (helical, 180 deg) -> over the top face
        -> down the opposite end face -> underside of the other bottom bar
        -> feed-through -> PCB pad
```

### Assembly

1. Print each section. The balls sit on the section above each cut, so the
   bottom-most section prints as modelled -- it has only an upward-opening
   dish -- while every section above it should be **flipped on the bed**, so
   its balls point up. Printed the other way up its mating face is a 1.35 mm
   overhang resting on ten ball tips. Set
   `joint_balls_on_lower_section=True` to move the balls to the other side
   instead, which lets every section print as modelled at the cost of a
   bridged dish ceiling. No supports needed either way.
2. Stack the sections: the balls drop into their sockets and pull the joint
   into line, and there is only one way they go together -- a blade is
   symmetric under a half turn about Z, so the other orientation the pins
   allow is the identical one. Cinch a zip tie through each pair of anchor
   holes flanking the joint; the tie wraps the blade between them and pulls
   the joint shut.
3. Run foil tape along all four tape paths above, across the joints. Keep it
   centred on the narrow end faces.
4. Pass the shorter loop's top tape through the crossover window in the
   taller blade.
5. Thread a zip tie through each hole, over the tape, and cinch.
6. End each tape run outside the sleeve and solder a short 1.5 mm copper
   wire to it, passed through that bar's feed-through into the PCB housing.
   The wire is straight: it lies flat on the tape, runs level through the
   hole, and comes out level with the board's top face.
7. Fit the balun PCB **rotated 45 deg to the frame's bars**, so its X-shaped
   pads line up with the four tape ends. Three M3 screws from below thread
   into the bosses; the fourth boss has no matching hole and stays empty.
   With the board pulled tight against its standoffs, each of the four wires
   is already lying flat on the board -- solder them to the pads as they
   sit.
8. Slide the frame onto the PVC mast until the pipe meets the PCB, then drive
   the three mast screws.

Everything is parametric -- see `PartSpec` in `cad/qfh_tape_frame_cad.py` for tape
land width, tie hole spacing, mast diameter, standoff height, max print
height, joint pin sizing, and so on. Set
`mast_pipe_od=None` for antennas too small to straddle a pipe (the bundled
913 MHz build is only ~46 mm across, so it has no sleeve).

### Tuning to a measured build

The calculator assumes bare conductor in free space. A real build runs foil
tape over printed plastic, which slows the wave and drags resonance *below*
the design frequency, so the first print of a new design lands low. Sweep
S11, then set `empirical_tuning_factor = f_measured / f_design` in the
`QfhInputSpec`: it scales every conductor length uniformly, so resonance
moves up by exactly `1 / factor` and everything else about the geometry is
unchanged.

The bundled 436 MHz build measured **408 MHz** on its first (uncorrected)
print -- 6.4% low -- so it now carries `empirical_tuning_factor = 408/436`
(0.9358). That takes the large loop from 226.8 mm tall / 99.8 mm across to
212.2 mm / 93.4 mm.

## See Also

* Geometry Calculator: https://jcoppens.com/ant/qfh/calc.en.php
* https://github.com/cernohorsky/QFH-Antenna-868MHz/blob/master/QFH-Antenna-868MHz.scad
* Balun PCB: https://github.com/ODZ-UJF-AV-CR/QFHBAL01

### Related, but less helpful

* https://github.com/cernohorsky/QFH-Antenna-868MHz
* https://www.thingiverse.com/thing:634205
* https://usradioguy.com/wp-content/uploads/2020/05/20200307-How-To-Build-A-QFH.pdf
* https://network.satnogs.org/stations/4704/
