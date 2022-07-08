from kneed import KneeLocator
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style

radius = []
index = []
i = 0
with open("D:\\Projects\\Spadas\\spadas_4-15\\index\\dss\\index\\beijing\\radius.txt", "r") as f:
    for line in f.readlines():
        line = line.strip('\n')
        radius.append(float(line))
        index.append(i)
        i += 1

if __name__ == '__main__':
    # print(index)
    # print(radius)
    my_kneed = KneeLocator(index, radius, S=1.0, curve='convex', direction='increasing')
    # my_kneed.plot_knee_normalized()
    plt.show()
    print(my_kneed.elbow)
    print(my_kneed.elbow_y)
    style.use('seaborn-whitegrid')
    plt.plot(index, radius)
    plt.annotate(text='Knee Point', xy=(my_kneed.knee, my_kneed.knee_y))
    plt.show()
