import matplotlib.pyplot as plt
import pandas as pd

def show_graph(filename):
    data = pd.read_csv(filename, header=None)
    data2 = pd.read_csv('utils/1.csv', header=None)

    plt.style.use('ggplot')
    plt.plot(data[0], data[1], 'go')
    plt.plot(data2[0], data2[1], 'ro')
    plt.xlabel('distance traveled')
    plt.ylabel('carbon emission')

    plt.show()

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        show_graph(sys.argv[1].strip())
    else:
        print('file path required.')