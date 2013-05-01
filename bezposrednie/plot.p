LU(x) = -2260.7 + 33.3602*x + -0.140903*x**2 + 0.000533983*x**3
Cholesky(x)= -1095.02 + 16.6564*x + -0.0595828*x**2 + 0.000316429*x**3
plot 'tmp.out' using 1:2:3 title 'LU' linecolor rgb 'blue' with yerrorbars,'tmp.out' using 1:4:5 title 'Cholesky' linecolor rgb 'green' with yerrorbars, Cholesky(x) with lines linecolor rgb 'green',LU(x) with lines linecolor rgb 'blue'

pause -1
set terminal svg
set output"output.svg"
replot