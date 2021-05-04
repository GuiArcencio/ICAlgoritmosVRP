def transformFile(filename):
    num_points = 0
    capacity = 0
    points = []
    best_known = 0
    with open(filename, 'r') as source_file:
        line = source_file.readline().strip()
        num_points = int(line.split()[0]) + 1
        best_known = float(line.split()[1])

        capacity = int(source_file.readline().strip())
        line = source_file.readline().strip()
        points.append((0, float(line.split()[0]), float(line.split()[1])))

        for k in range(num_points - 1):
            line = source_file.readline().strip()
            points.append((int(line.split()[3]), float(line.split()[1]), float(line.split()[2])))

    with open(filename, 'w') as dest_file:
        dest_file.write(f'{num_points} {num_points} {capacity}\n')
        for k in range(num_points):
            dest_file.write(f'{points[k][0]} {points[k][1]:.1f} {points[k][2]:.1f}\n')
        dest_file.write(f'\nBest known: {best_known:.2f}\n')

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        transformFile(sys.argv[1].strip())
    else:
        print('file path required.')