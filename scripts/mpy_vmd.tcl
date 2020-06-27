namespace eval loadedFiles {
  # namespace to hold global varsmpt and gro file names
  variable mpt 
  variable gro
}

proc mpyload {mpt_file gro_file} {
	# convenience function to load both mpt and gro
	set ::loadedFiles::mpt $mpt_file
	
	set ::loadedFiles::gro $gro_file
	
	set molid [mol load gro $gro_file]
	puts $molid
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
	
	# execute vmd.py
	set out [exec python test.py $::loadedFiles::mpt $::loadedFiles::gro $cpmd $ndx $name $type $index $mass $element $resname $resid $x $y $z]
	puts $out # display messages from python
}