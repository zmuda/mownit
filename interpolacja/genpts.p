plot 	'out.tmp' using 1:2 title 'polynomial' with lines linewidth 3 linecolor rgb "red",'out.tmp' using 1:3 title 'lagrange' with lines linewidth 2 linecolor rgb "black", 'out.tmp' using 1:4 title 'newton' with lines linecolor rgb "blue", 'outPts.tmp' using 1:2 title 'input' linecolor rgb "green", 'out.tmp' using 1:5 title 'cspline' with lines linecolor rgb "yellow", 'out.tmp' using 1:6 title 'akima' with lines linecolor rgb "pink"
pause -1
set terminal svg
set output "output.svg"
replot
