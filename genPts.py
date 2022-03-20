import random
with open("Groundtruth.vtk", 'w') as f:
    f.write('# vtk DataFile Version 3.0\n')
    f.write('vtk output\n')
    f.write('ASCII\n')
    f.write('DATASET UNSTRUCTURED_GRID\n')
    f.write('POINTS 1000000 double\n')
    for i in range(500000):
        x = round(random.uniform(0.0, 10.0),5)
        y = round(random.uniform(0.0, 10.0),5)
        z = round(random.uniform(0.0, 10.0),5)

        # write start and end
        f.write(f'{x} {y} {z}\n')
        f.write(f'{x} {y} {z}\n')

