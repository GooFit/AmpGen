import os
import random

f = open("AuxBlinding.opt","w")

random.seed(2019)

blindingoffset = round(random.uniform(1,10)*10,4)

for i in range(15): f.write("\n")
f.write("#          Empty lines just prevent accidental unblinding. \n")
for i in range(100): f.write("\n")

f.write("BLINDINGOFFSET 5 ")
f.write(str(blindingoffset)+ " 0. \n")
f.close()
