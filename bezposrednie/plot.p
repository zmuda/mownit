LU(x) = -1063.52 + 17.151*x + -0.0803813*x**2 + 0.000460895*x**3
Cholesky(x)= -1353.46 + 19.5649*x + -0.0649049*x**2 + 0.000303124*x**3
plot 'tmp.out' using 1:2:3 title 'LU' linecolor rgb 'blue' with yerrorbars,'tmp.out' using 1:4:5 title 'Cholesky' linecolor rgb 'green' with yerrorbars, Cholesky(x) with lines linecolor rgb 'green',LU(x) with lines linecolor rgb 'blue'

pause -1
set terminal svg
set output"output.svg"
replot