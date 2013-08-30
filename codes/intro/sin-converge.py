import math

X = [0.5, 0.25, 0.125, 0.0625]

for x in X:
    print x, math.sin(x) - (x - x**3/6.0)

