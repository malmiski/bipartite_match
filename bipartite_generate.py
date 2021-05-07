from random import random
f = open("bipartite_set_4.txt", "a")
a = []
b = []
for i in range(2000):
  a.append((random()*100, random()*400))
  b.append((random()*300, random()*100))
f.write(str(a) + '\n' + str(b))
f.close()
