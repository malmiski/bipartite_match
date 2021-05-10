from random import random
from math import cos, sin, sqrt
import sys
def gen_circle_rand(x,y,r):
    rad = r*sqrt(random())
    theta = random()*2*3.14159
    return (x + rad*cos(theta), y + rad*sin(theta))
def generate(center_a_x, center_a_y, radius_1, center_b_x, center_b_y, radius_2, n, name):
    a = []
    b = []
    for i in range(int(n)):
        a.append(gen_circle_rand(center_a_x, center_a_y, radius_1))
        b.append(gen_circle_rand(center_b_x, center_b_y, radius_2))
    f = open(name, "w")
    f.write(str(a) + '\n' + str(b))
    f.close()

if(len(sys.argv) < 8):
    print("Too few args")
for i in range(1, 8):
    sys.argv[i] = float(sys.argv[i])
generate(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8])
print("Done")
