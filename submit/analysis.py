import numpy as np
from sklearn.metrics import silhouette_score
import kmeans
import symnmf
import sys

def read_data_points(data_points):
    data_list = np.loadtxt(data_points,delimiter=',')
    return data_list.tolist()



def Symnmf_analysis(data ,k):
    clusters = np.array(symnmf.symnmf(data,k))
    return clusters.argmax(axis=1)


def Kmeans_analysis(data ,k):
    n = len(data)
    d = len(data[0])
    clusters = kmeans.kMeans(data,k,n,d)
    return clusters

def main():
    data_path = sys.argv[2]
    k = int(sys.argv[1])
    data_points = read_data_points(data_path)
    #finding clusters
    symnmf_clusters = Symnmf_analysis(data_points,k)
    kmeans_clusters = Kmeans_analysis(data_points,k)

    #finding score
    symnmf_score = silhouette_score(data_points,symnmf_clusters)
    Kmeans_score = silhouette_score(data_points,kmeans_clusters)
    
    print(f"nmf: {symnmf_score:.4f}")
    print(f"kmeans: {Kmeans_score:.4f}")
    
if __name__ == "__main__":
    main()
    