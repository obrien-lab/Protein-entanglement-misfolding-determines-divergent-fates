# Entanglement type: gain
package require topotools
display rendermode GLSL
axes location off

color Display {Background} white

mol new /storage/home/yuj179/mygroup/protein_Ubq/new/Tq/NU_E/Q9H6R3/setup/Q9H6R3.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol delrep 0 top
mol representation NewCartoon 0.300000 10.000000 4.100000 0
mol color ColorID 6
mol selection {resid 1 to 686}
mol material AOChalky
mol addrep top
mol representation NewCartoon 0.300000 10.000000 4.100000 0
mol color ColorID 4
mol selection {resid 181 to 342}
mol material Opaque
mol addrep top
mol representation NewCartoon 0.300000 10.000000 4.100000 0
mol color ColorID 12
mol selection {resid 157 to 176}
mol material Opaque
mol addrep top
mol representation VDW 1.000000 12.000000
mol color Name
mol selection {not resid 1 to 686 and not water}
mol material Opaque
mol addrep top

mol new ./state_8.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol delrep 0 top
mol representation NewCartoon 0.300000 10.000000 4.100000 0
mol color ColorID 6
mol selection {all}
mol material AOChalky
mol addrep top
mol representation NewCartoon 0.350000 10.000000 4.100000 0
mol color ColorID 1
mol selection {resid 181 to 342}
mol material Opaque
mol addrep top
mol representation NewCartoon 0.350000 10.000000 4.100000 0
mol color ColorID 0
mol selection {resid 157 to 176}
mol material Opaque
mol addrep top
mol representation VDW 1.000000 12.000000
mol color ColorID 3
mol selection {resid 181 342 and name CA}
mol material Opaque
mol addrep top

set sel [atomselect top "resid 181 342 and name CA"]
set idx [$sel get index]
topo addbond [lindex $idx 0] [lindex $idx 1]
mol representation Bonds 0.300000 12.000000
mol color ColorID 3
mol selection {resid 181 342 and name CA}
mol material Opaque
mol addrep top

set sel1 [atomselect 0 "resid 1 to 686 and not (resid 181 to 342 157 to 176) and name CA"]
set sel2 [atomselect 1 "resid 1 to 686 and not (resid 181 to 342 157 to 176) and name CA"]
set trans_mat [measure fit $sel1 $sel2]
set move_sel [atomselect 0 "all"]
$move_sel move $trans_mat
