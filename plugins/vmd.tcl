
##### MiMiCPy VMD settings script
##
proc prepqm {mpt {sele atomselect0} {molid 0} {cpmd cpmd.inp} {ndx index.ndx}} {
	if { [file exists $mpt] == 0} {
		puts "$mpt_file not found!"
		return
	}
	
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
	puts $[exec mimicpy_vmd $mpt $cpmd $ndx $name $type $index\
	 $mass $element $resname $resid $x $y $z $a $b $:c $alpha\
	  $beta $gamma]
}
##
##################################