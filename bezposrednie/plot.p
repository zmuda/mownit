LU(x) = -1721.85 + 24.7211*x + -0.0981444*x**2 + 0.000483919*x**3
Cholesky(x)= -990.221 + 14.9992*x + -0.0497547*x**2 + 0.000297139*x**3
plot 'tmp.out' using 1:2:3 title 'LU' linecolor rgb 'blue' with yerrorbars,'tmp.out' using 1:4:5 title 'Cholesky' linecolor rgb 'green' with yerrorbars, Cholesky(x) with lines linecolor rgb 'green',LU(x) with lines linecolor rgb 'blue'

pause -1
set terminal svg
set output"output.svg"
replot