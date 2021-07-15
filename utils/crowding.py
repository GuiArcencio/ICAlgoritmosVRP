import pandas as pd

def generateCrowding(filename):
    solutions = pd.read_csv(filename, header=None)
    S = []
    for i, sol in solutions.iterrows():
        S.append((sol[0], sol[1], i))
    
    distance = [0.0]*solutions.shape[0]

    for k in range(2): # número de funções objetivo
        S.sort(key = lambda x: x[k], reverse=True)
        
        for j in range(1, len(S) - 1):
            delta_j = (S[j+1][k] - S[j-1][k]) / (S[len(S)-1][k] - S[0][k])
            distance[S[j][2]] += delta_j

        distance[S[0][2]] = float("inf")
        distance[S[len(S)-1][2]] = float("inf")

    print(distance)
        

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        generateCrowding(sys.argv[1].strip())
    else:
        print('file path required.')