import numpy as np
from sklearn.metrics import silhouette_score
import kmeans
import symnmf as s
import sys

def read_data_points(data_points):
    return np.loadtxt(data_points,delimiter=',')


def initial_H(W, k, n):
    m = np.mean(W)
    upper = 2 * np.sqrt(m/k)
    H_init = np.random.uniform(0, upper, size=(n, k))
    return H_init

def Symnmf_analysis(data ,k):
    n, d = data.shape
    W = s.norm(data)
    H_init = initial_H(W, k,n)   
    H_init = H_init.tolist()
    H_final = s.symnmf(H_init, W)
    clusters = np.array(H_final)
    return clusters.argmax(axis=1)

def Kmeans_analysis(data ,k):
    n,d = data.shape
    clusters = kmeans.kMeans(data,k,n,d)
    return clusters

def main():
    data_path = sys.argv[2]
    k = int(sys.argv[1])
    data_points = read_data_points(data_path)
    #finding clusters
    symnmf_clusters = Symnmf_analysis(data_points,k)
    symnmf_Kmeans = Kmeans_analysis(data_points,k)
    
    #finding score
    symnmf_score = silhouette_score(data_points,symnmf_clusters)
    Kmeans_score = silhouette_score(data_points,symnmf_Kmeans)
    
    print(f"nmf: {symnmf_score:.4f}")
    print(f"kmeans: {Kmeans_score:.4f}")
    
if __name__ == "__main__":
    main()
    