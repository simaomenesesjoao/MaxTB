
using PyPlot

x = LinRange(0,1,100)
y = x.^2

p = plot(x,y)
display(p)
readline()
