# wrap_pbc.tcl -- SURF 2026
# Load a LAMMPS trajectory, wrap PBCs, and set up a clean representation
# of ions at a charged planar surface.
#
# Usage in VMD:  vmd -e wrap_pbc.tcl -args system.data prod.dcd

if { $argc < 2 } {
    puts "usage: vmd -e wrap_pbc.tcl -args <data> <dcd>"
    exit 1
}

set datafile [lindex $argv 0]
set trajfile [lindex $argv 1]

# Load topology + trajectory
mol new     $datafile type lammpsdata waitfor all
mol addfile $trajfile type dcd        waitfor all

# Wrap all atoms into the primary cell, centered on the box
package require pbctools
pbc wrap -all -center com -centersel "all"
pbc box  -center com

# Clean default representation
mol delrep 0 top

# Water: transparent lines
mol representation Lines 1.0
mol selection      "type 1 2"
mol color          Name
mol addrep top

# Ions: VDW
mol representation VDW 0.8 20
mol selection      "type 3 or type 4 or type 5"
mol color          Type
mol addrep top

# Electrode / wall (assumes type 6 — adjust for your system)
mol representation Surf 1.4
mol selection      "type 6"
mol color          ColorID 6
mol addrep top

display projection Orthographic
axes location off
color Display Background white
