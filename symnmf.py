import sys
import numpy as np
import symnmf

def read_data_points(data_points): 
    return np.loadtxt(data_points)

def output_matrix(matrix):
    np.matrix(matrix)
    
def initial_H(W,k):
    m = np.mean(W)
    upper = 2 * np.sqrt(m/k)
    H_init = np.random.uniform(0,upper,size=(W.shape[0],k))
    return H_init
   
def output_matrix(matrix):
    for r in matrix:
        print(','.join(f"{value:.4f}" for value in r))
        
            
def parseArgs(args):
    if len(args) == 3: 
        try:
            k = int(args[0]) 
        except:
            print("Invalid number of clusters!")    
            exit()
            
        goal = args[1]
        if(goal != "symnmf" and goal != "sym" and goal != "ddg" and goal != "norm"):
            print("Invalid goal name!")    
            exit()
            
        filePath = args[2]
                
        
   
    return k, goal,filePath  
    
if __name__ == '__main__':
    k, goal,filePath  = parseArgs(sys.argv)
    dataPoints = np.loadtxt(filePath)
    
    try:
        if(goal == "symnmf"):
           #beginning of the code ?
            np.random.seed(0)
            A = symnmf.sym(dataPoints)
            D = symnmf.ddg(A,n)
            W = symnmf.norm(D,A,n)
            H_init = initial_H(W,k)
            H_final = symnmf.calcOptimalDecompMatrix(H_init,W,n,k)
            output_matrix(H_final)
            
        elif (goal == "sym"):
            A = symnmf.sym(dataPoints)
            output_matrix(A)
            
        elif (goal == "ddg"):
            A = symnmf.sym(dataPoints)
            D = symnmf.ddg(dataPoints)
            output_matrix(D)
        elif (goal == "norm"):
            A = symnmf.sym(dataPoints)
            D = symnmf.ddg(dataPoints)
            W = symnmf.calcNormalizedSymilarityMatrix(dataPoints,D)
            output_matrix(W)
        
    except:
        print("An Error Has Occurred")
        exit()
