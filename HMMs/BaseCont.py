import math
import matplotlib.pyplot as plt

baseIDX = {"A": 0, "C": 1, "G": 2, "T": 3}


def main():
    spaciiFA = "MSpacii.fa"
    pathogenFA = "pathogen.fa"
    spaciiFA_T = "MSpacii_training.fa"
    pathogenFA_T = "pathogen_training.fa"

    spaciiID2seq = getSeq(spaciiFA)
    pathogenID2seq = getSeq(pathogenFA)

    spaciiTrainModel = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
    pathTrainModel = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]

    spaciiTrainModel = trainModel(spaciiTrainModel, spaciiFA_T)
    pathTrainModel = trainModel(pathTrainModel, pathogenFA_T)

    markovScoresSpacii = []
    markovScoresPath = []

    for ID in spaciiID2seq.keys():
        markovScoresSpacii.append(getLogLike(spaciiTrainModel, pathTrainModel, spaciiID2seq[ID]))

    for ID in pathogenID2seq.keys():
        markovScoresPath.append(getLogLike(spaciiTrainModel, pathTrainModel, pathogenID2seq[ID]))

    #### ----- output -----
    plt.hist([markovScoresPath, markovScoresSpacii], bins=20, label=['pathogen', 'spacii'], rwidth=1, density=True)
    plt.xlabel('Log-likelihood Ratio Score')
    plt.ylabel('Density')
    plt.legend()
    plt.show()
    scoresOutputText(markovScoresSpacii, markovScoresPath)
    #### ---- output -----


def scoresOutputText(markovScoresSpacii, markovScoresPath):
    f = open("results.tab", "w")
    f.write("SpaciiScores\tpathogenScores\n")
    for i in range(len(markovScoresSpacii)):
        f.write(str(markovScoresSpacii[i]) + "\t" + str(markovScoresPath[i]) + "\n")
    f.close()


def getLogLike(model1, model2, seq):
    Pmod1 = 1
    Pmod2 = 1

    for i in range(1, len(seq)):
        prev_base = baseIDX[seq[i - 1]]
        curr_base = baseIDX[seq[i]]
        Pmod1 *= model1[prev_base][curr_base]
        Pmod2 *= model2[prev_base][curr_base]

    if Pmod1 > 0 and Pmod2 > 0:
        score = math.log(Pmod1) - math.log(Pmod2)
    else:
        score = float('-inf')

    return score


def trainModel(model, data):
    seq_data = getSeq(data)
    dinuc_counts = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]

    for seq in seq_data.values():
        for i in range(1, len(seq)):
            prev_base = baseIDX[seq[i - 1]]
            curr_base = baseIDX[seq[i]]
            dinuc_counts[prev_base][curr_base] += 1

    for i in range(4):
        total_count = sum(dinuc_counts[i])
        if total_count > 0:
            for j in range(4):
                model[i][j] = dinuc_counts[i][j] / total_count

    print("Trained Model:")
    for row in model:
        print(row)
    return model


def getSeq(filename):
    f = open(filename)
    id2seq = {}
    currkey = ''

    for line in f:
        if line.find(">") == 0:
            currkey = line.rstrip()[1:]
            id2seq[currkey] = ''
        else:
            id2seq[currkey] = id2seq[currkey] + line.rstrip()

    return id2seq



main()