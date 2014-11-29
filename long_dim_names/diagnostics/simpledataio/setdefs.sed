:mline
/\\$/{
 N
 s,\\\n,,
 b mline
}
t clear
:clear
s/^[   ]*#[  ]*define[   ][  ]*\([^  (][^  (]*([^)]*)\)[   ]*\(.*\)/PP_DEFINE\1=\2/g
t quote
s/^[   ]*#[  ]*define[   ][  ]*\([^  ][^   ]*\)[   ]*\(.*\)/PP_DEFINE\1=\2/g
t quote
b any
:quote
s/[  `~#$^&*(){}\\|;'\''"<>?]/\\&/g
s/\[/\\&/g
s/\]/\\&/g
s/\$/$$/g
H
:any
${
  g
  s/^\n//
  s/\n/ /g
  p
}
