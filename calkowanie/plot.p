set terminal svg
set output "komp.svg"
plot [0: 1] 3*(x*x-2*x)

set output "osc.svg"
plot [0: 3.14] x * sin(30 * x) * cos(x)
set output "osob.svg"
plot [0: 1] log(x) / sqrt(x)
set output "adap.svg"
plot [0.01: 100] sin(x*2)/x




