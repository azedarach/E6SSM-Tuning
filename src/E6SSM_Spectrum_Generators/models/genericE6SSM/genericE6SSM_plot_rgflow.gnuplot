set terminal x11

set title "genericE6SSM renormalization group flow"
set xlabel "renormalization scale / GeV"
set logscale x

if (!exists("filename")) filename='genericE6SSM_rgflow.dat'

plot for [i=2:175+1] filename using 1:(column(i)) title columnhead(i)
