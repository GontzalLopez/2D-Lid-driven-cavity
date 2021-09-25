reset session

set terminal png 1 size 500,500 font "Verdana,10" 
set title "u velocity [m/s]"
set xlabel "x (meter)"
set ylabel "y (meter)"
set palette rgbformulae 22,13,-31
set xrange [0:1]
set yrange [0:1]

#set pm3d map # no interpolation
#set pm3d map interpolate 1,3 # no interpolation in y, 3 interpolation in x
	#set pm3d map interpolate 0,0 # automatic interpolation
	#splot "dataCavity.dat" u 1:2:3

plot "dataCavity.dat" u 1:2:3 with image
#set output 'contour_u.png'

set terminal png 2 size 500,500 font "Verdana,10" 

set title "Velocity components at the mid-cross section y=0.5"
set xlabel "u,v [m/s]"
set ylabel "x [m]"
set xrange [0:1]
set yrange [-0.25:0.25]
set grid
plot "dataCavity_prof.dat" u 1:2 with points pt 7 title "u componet", \
     "dataCavity_prof.dat" u 1:3 with points pt 7 title "v component"
#set output 'profiles.png'

set terminal wxt 3 size 500,500 font "Verdana,10" 

set title "u velocity at the mid-cross section x=0.5"
set xlabel "u [m/s]"
set ylabel "y [m]"
set xrange [0:1]
set yrange [-0.3:1]
set grid
set key left
plot "dataCavity_uExp.dat" u 1:2 with points pt 7 title "Numerical", \
     "exp_data.dat" u 1:2 with points pt 7 title "Experimental"
#set output 'numVSexp.png'







