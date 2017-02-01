import os

cmd = 'mpirun -n 300 ./driver_nomr --infile final-final/by-draw/output{:04d}.json --sim nomr'

for i in range(0, 500):

    print('Running %d' % i)
    os.system(cmd.format(i))

