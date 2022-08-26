from turtle import color
import matplotlib.pyplot as plt;


def PathReader(dir1, color1):
    f1 = open(dir1, "r")
    line1 = f1.readline()
    x, y = [], []
    while line1:
      linedate1 = line1.split(",")
      if line1 == "\n":
        break
      x.append(float(linedate1[0]))
      y.append(float(linedate1[1]))
      line1 = f1.readline()
    x.pop(0)
    y.pop(0)
    plt.plot(x, y, color1)


def main():
    PathReader("./data/rightBound.txt", "g")
    PathReader("./data/leftBound.txt", "b")
    PathReader("./data/refPath.txt", "r")
    PathReader("./data/refinePath.txt", "purple")
    plt.show()
main()
