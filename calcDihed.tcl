set fl [open "fdihed.dat" w]
set nframe [molinfo "top" get numframes]
for {set iframe 0} {$iframe < $nframe} {incr iframe} {
  set framedihed [measure dihed {3 0 5 7} frame $iframe]
  puts $fl "$iframe $framedihed"
}
close $fl
}

set f2 [open "sdihed.dat" w]
set nframe [molinfo "top" get numframes]
for {set iframe 0} {$iframe < $nframe} {incr iframe} {
  measure dihed {7 5 9 12} frame $iframe
  puts $f2 "$iframe $framedihed"
}
close $f2
}
