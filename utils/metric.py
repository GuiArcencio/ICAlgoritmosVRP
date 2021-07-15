import pandas as pd

def generateCrowding(filename):
    solutions = pd.read_csv(filename, header=None)
    S = []
    for _, sol in solutions.iterrows():
        S.append((sol[0], sol[1]))

    ideal_point = (min(S, key = lambda x: x[0])[0] , min(S, key = lambda x: x[1])[1])
    nadir_point = (max(S, key = lambda x: x[0])[0] , max(S, key = lambda x: x[1])[1])

    DM = 0.0
    for h in range(2):
        S.sort(key=lambda x: x[h])
        distances = []
        for i in range(len(S) - 1):
            distances.append(S[i+1][h] - S[i][h])
        mean_distance = sum(distances) / len(distances)
        sd_distance = 0.0
        for d in distances:
            sd_distance += (d - mean_distance)**2
        sd_distance /= len(distances) - 1
        sd_distance = sd_distance ** 0.5

        obj_range = nadir_point[h] - ideal_point[h]

        DM += (sd_distance / mean_distance) / (obj_range / abs(ideal_point[h] - nadir_point[h]))

    DM /= len(S)
    print(DM)
        

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        generateCrowding(sys.argv[1].strip())
    else:
        print('file path required.')