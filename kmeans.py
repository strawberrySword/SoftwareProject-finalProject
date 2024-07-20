import sys

def kMeans(dataPoints, k, n, d, iter=200 ):
    # initialize centroids to first k data points
    epsilon = 0.001
    centroids = []
    for i in range(k):
        centroids.append({'center': dataPoints[i].copy(), 'currentCenter': [0]*d, 'size': 0})
        
    i = 0
    maxDelta = epsilon + 1
    while(i < iter and maxDelta > epsilon):
        for j in range(n):
            closestCluster = findClosestCluster(dataPoints[j], centroids)
            centroids[closestCluster]['size'] += 1
            centroids[closestCluster]['currentCenter'] = sumList(centroids[closestCluster]['currentCenter'], dataPoints[j])
        
        maxDelta = 0 
        for u in centroids:
            u['currentCenter'] = [x / (u['size']) for x in u['currentCenter']]
            delta = calcEclideanDistance(u['currentCenter'], u['center'])
            if(delta > maxDelta):
                maxDelta = delta
            u['center'] = u['currentCenter'].copy()
            u['currentCenter'] = [0]*d
            u['size'] = 0
        i += 1
        
    for u in centroids:
        formatted = [ '%.4f' % elem for elem in u['center'] ]
        print(','.join(formatted))

def sumList(a,b):
    return [x + y for x,y in zip(a,b)]

def findClosestCluster(datapoint, centroids):
    minDistance = calcEclideanDistance(centroids[0]['center'], datapoint)
    closestCluster = 0
    for index, u in enumerate(centroids):
        distance = calcEclideanDistance(u['center'], datapoint)
        if(distance < minDistance):
            minDistance = distance
            closestCluster = index
    return closestCluster

def calcEclideanDistance(u, v):
    squareSum = 0
    for i in range(len(u)):
        squareSum += (u[i] - v[i])**2
    return squareSum**0.5

def parseArgs(args):
    try:
        n = int(args[2])
        if(n < 2):
            exit()
    except:
        print("Invalid number of points!")
        exit()
    try:
        k = int(args[1])
        if(k < 2 or k > n):
            exit()    
    except:
        print("Invalid number of clusters!")    
        exit()
    try:
        d = int(args[3])
        if(d<=0):
            exit()    
    except:
        print("Invalid dimension of point!")
        exit()
    if len(args) == 5:
        iter = 200
        filePath = args[4]
    if len(args) == 6:
        try:
            iter = int(args[4])
            if(iter < 2 or iter >= 1000):
                exit()
        except:
            print("Invalid maximum iteration!")
            exit()
        filePath = args[5]
    return k, n, d, iter, filePath

def parseDataPoints(filePath):
    f = open(filePath, 'r')
    raw = f.read()
    lines = raw.split('\n')
    dataPoints = [l.split(',') for l in lines if len(l) > 0 ]
    dataPoints = [[float(x) for x in dp] for dp in dataPoints]

    f.close()
    return dataPoints

if __name__ == '__main__':
    k, n, d, iter, filePath = parseArgs(sys.argv)
    list = parseDataPoints(filePath)

    try:
        kMeans(list, k, n, d, iter)
    except:
        print("An Error Has Occurred")
        exit()