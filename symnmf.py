import sys
import numpy as np
import symnmf as s


def read_data_points(data_points):
    data_list = np.loadtxt(data_points,delimiter=',')
    return data_list.tolist()
    
    
def initial_H(W, k, n):
    m = np.mean(W)
    upper = 2 * np.sqrt(m/k)
    H_init = np.random.uniform(0, upper, size=(n, k))
    return H_init


def output_matrix(matrix):
    for r in matrix:
        print(','.join(f"{value:.4f}" for value in r))


def symnmf1(dataPoints,k):
    W = s.norm(dataPoints, D)
    H_init = initial_H(W, k)
    H_final = s.symnmf(H_init, W)
    return H_final

def parseArgs(args):
    if len(args) == 4:
        try:
            k = int(args[1])
        except:
            print("Invalid number of clusters!")
            exit()
        goal = args[2]
        if (goal != "symnmf" and goal != "sym" and goal != "ddg" and goal != "norm"):
            print("Invalid goal name!")
            exit()

        filePath = args[3]
        return k, goal, filePath



if __name__ == '__main__':
    k, goal, filePath = parseArgs(sys.argv)
    dataPoints = read_data_points(filePath)   
    print(dataPoints)
    n = len(dataPoints)
    try:
        if (goal == "symnmf"):
            np.random.seed(0)
            W = s.norm(dataPoints)
            H_init = initial_H(W, k,n)   
            H_init = H_init.tolist()    
            for line in H_init:
                print(line)
            print("---------AND NOW-------")
            H_final = s.symnmf(H_init, W)
            output_matrix(H_final)

        elif (goal == "sym"):
            A = s.sym(dataPoints)
            output_matrix(A)

        elif (goal == "ddg"):
            D = s.ddg(dataPoints)
            print(D)
        elif (goal == "norm"):
            W = s.norm(dataPoints)    
            output_matrix(W)

    except:
        print("last one Error")
        print("An Error Has Occurred")
        exit()
