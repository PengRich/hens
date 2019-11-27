import os
i = 1
while True:
    os.system("mpirun -np 2 debug/para_sade_4sp1")
    i += 1
    if i==3:
        break
