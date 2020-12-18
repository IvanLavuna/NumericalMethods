import matplotlib.pyplot as plt
import json

X_init = []
Y_init = []
X = []
Y = []
data = None


if __name__ == '__main__':

    with open("./../../data/json/output_points.json", "r") as file:
        data = json.load(file)
    if data:
        for item in data:
            X.append(item["x"])
            Y.append(item["y"])

    with open("./../../data/json/input_points.json", "r") as file:
        data = json.load(file)
    if data:
        for item in data:
            X_init.append(item["x"])
            Y_init.append(item["y"])

    for i in range(len(X_init)):
        plt.scatter(X_init[i], Y_init[i], s=100)

    plt.plot(X, Y)
    plt.xlabel('x - axis')
    plt.ylabel('y - axis')
    plt.title('Lagrange interpolate function')

    plt.show()
