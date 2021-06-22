import os
import random

def add_carbon(filename):
    random.seed(os.urandom(32))
    num_points = 0
    with open(filename, 'r') as source_file:
        line = source_file.readline().strip()
        num_points = int(line.split()[0])
    with open(filename, 'a') as dest_file:
        dest_file.write('\n')
        for i in range(1, num_points):
            factors = []
            for j in range(i):
                factors.append(random.uniform(0.5, 2.0))
            dest_file.write(f'{" ".join(format(x, ".5f") for x in factors)}\n')

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        add_carbon(sys.argv[1].strip())
    else:
        print('file path required.')