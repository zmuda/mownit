plot 	'out.tmp' using 1:2 title 'polynomial' with lines,	       		'out.tmp' using 1:3 title 'lagrange' with lines,	      		'out.tmp' using 1:4 title 'newton' with lines,	  		'outPts.tmp' using 1:2 title 'input',			   		'out.tmp' using 1:5 title 'cspline' with lines,		   		'out.tmp' using 1:6 title 'akima' with lines
pause -1
