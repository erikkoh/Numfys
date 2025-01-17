using Plots

time = range(0,100,1000)
x = sin.(10*time)
plot(time,x)