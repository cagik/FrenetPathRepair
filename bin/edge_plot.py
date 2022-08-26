from turtle import color
import matplotlib.pyplot as plt;

def readData(str):
    data = []
    with open(str, "r") as f:
        for line in f.readlines():
            line = line.strip('\n')
            data.append(float(line))
    data.pop(0)
    return data

def main():
    edgel = readData("./data/ref_edge_l.txt")
    edger = readData("./data/ref_edge_r.txt")
    edgedl = readData("./data/update_edge_l.txt")
    edgedr = readData("./data/update_edge_r.txt")
    res_s = readData("./data/refine_s.txt")
    plt.plot(edgel, color = 'r')
    plt.plot(edger, color = 'r')
    plt.plot(edgedl, color = 'b')
    plt.plot(edgedr, color = 'b')
    plt.plot(res_s, color = 'g')

    plt.show()
main()
