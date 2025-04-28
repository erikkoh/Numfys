

x = 10_000_000
y = 1.0
x1 = Float16(x)
y1 = Float16(y)
x2 = Float32(x)
y2 = Float32(y)
x3 = Float64(x)
y3 = Float64(y)

z1 = x1 + y1
z2 = x2 + y2
z3 = x3 + y3
true_value = BigFloat(x) + BigFloat(y)
println("Float16: ", z1, " Float32: ", z2, " Float64: ", z3, " BigFloat: ", true_value)

a = 0.1 +0.2
b = 0.3 
println("0.1 + 0.2 == 0.3: ", a == b)
println(" abs(a-b) < 1e-10: ", abs(a-b) < 1e-10)


y = 1.0 + 10^6*10^(-16)
x = 1.0 
for i in 1:10^6
    global x = x + 10^(-6)
end

println("x: ", x, " y: ", y, " x - y = ", x - y)	


facotrial_list = [i for i in 1:171 ]
fac_171 = Float64(1.0)
for i in facotrial_list
    global fac_171 *= i
end
println("Facotrial of 171: ", fac_171)
fac_170 = Float64(1.0)
for i in facotrial_list[1:end-1]
    global fac_170 *= i
end
println("Facotrial of 170: ", fac_170)
println("Directly 171! - 170!: ", fac_171 - fac_170)
println("Using 171*(170!- 1/171*170!): ", 171*(fac_170 - fac_170/171))
println("Stirling: ", exp(171*log(171) - 171) - exp(170*log(170) - 170))
println("171! using Stirling: ", exp(171*log(171) - 171 + log(2*pi*171) / 2))


a = Float32(1.0000000000010)
println("deltaa", Float32(10^-12))
b = Float32(1.0)
println("a: ", Float32(a), " b: ", Float32(b), " a - b = ", sqrt(Float32(a)) - sqrt(Float32(b)), " Using BigFloat: ", sqrt(BigFloat(a)) - sqrt(BigFloat(b)))
println((sqrt(Float32(10^12 + 1)) - sqrt(Float32(b*10^12)))*10^-6)
println((Float32(10^-12))/(sqrt(Float32(a))+sqrt(Float32(b))))
