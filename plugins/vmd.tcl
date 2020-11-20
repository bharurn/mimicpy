
##### MiMiCPy VMD settings script
##
proc prepqm {top {sele atomselect0} {molid 0} {inp None} {mdp None} {ndx index.ndx} {out cpmd.inp}} {
	set a [molinfo $molid get a]
	set b [molinfo $molid get b]
	set c [molinfo $molid get c]
	set alpha [molinfo $molid get alpha]
	set beta [molinfo $molid get beta]
	set gamma [molinfo $molid get gamma]
	
	# get all params from selection
	set name [$sele get name]
	set type [$sele get type]
	set index [$sele get index]
	set mass [$sele get mass]
	set element [$sele get element]
	set resname [$sele get resname]
	set resid [$sele get resid]
	set x [$sele get x]
	set y [$sele get y]
	set z [$sele get z]
	
	# execute mimicpy main_vmd script & disp logger output
	puts $[exec mimicpy_vmd $top $inp $mdp $ndx $out $name $type $index\
	 $mass $element $resname $resid $x $y $z $a $b $c $alpha\
	  $beta $gamma]
}
##
##################################