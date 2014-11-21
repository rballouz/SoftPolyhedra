#!/bin/tcsh
#
# makes movie from already existing pngs that are ss.*.r.png, argument is frame rate
#
if !($#argv) set argv = 30
echo $1
set pngs = `\ls -1d *.png`
rm -f ffmpeg*.png >& /dev/null
@ i = 1
foreach png ($pngs)
    ln -s $png ffmpeg`printf %04d $i`.png
    @ i++
end
echo $#pngs

ffmpeg -sameq -r $1 -b 9600 -y -i ffmpeg%04d.png movie.mp4 << EOF

EOF
rm -f ffmpeg*.png >& /dev/null
exit 0

#############
#############

syntax:
echo syntax: $0