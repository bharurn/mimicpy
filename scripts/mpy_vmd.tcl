namespace eval globalvars {
  # namespace to hold global vars
  variable mpt 
  variable gro
  variable a
  variable b
  variable c
  variable alpha
  variable beta
  variable gamma
}

proc mpyload {mpt_file gro_file} {
	# convenience function to load both mpt and gro
	set ::globalvars::mpt $mpt_file
	if { [file exists $mpt_file] == 0} {
		puts "$mpt_file not found!"
		return
	}
	
	set ::globalvars::gro $gro_file
	
	set molid [mol load gro $gro_file]
	
	set ::globalvars::a [molinfo $molid get a]
	set ::globalvars::b [molinfo $molid get b]
	set ::globalvars::c [molinfo $molid get c]
	set ::globalvars::alpha [molinfo $molid get alpha]
	set ::globalvars::beta [molinfo $molid get beta]
	set ::globalvars::gamma [molinfo $molid get gamma]
	
	return $molid
}

proc prepqm {sele {cpmd cpmd.inp} {ndx index.ndx}} {
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
	
	# execute python script & disp output
	puts $[exec python mpy_vmd.py $::globalvars::mpt $::globalvars::gro $cpmd $ndx $name $type $index\
	 $mass $element $resname $resid $x $y $z $::globalvars::a $::globalvars::b $::globalvars::c $::globalvars::alpha\
	  $::globalvars::beta $::globalvars::gamma]
}