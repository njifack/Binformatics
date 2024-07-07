import math

modelCodons = ['TTT', 'TTC', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA',
               'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTT', 'GTC',
               'GTA', 'GTG', 'TCT', 'TCC', 'TCA', 'TCG', 'AGT',
               'AGC', 'CCT', 'CCC', 'CCA', 'CCG', 'ACT', 'ACC',
               'ACA', 'ACG', 'GCT', 'GCC', 'GCA', 'GCG', 'TAT',
               'TAC', 'CAT', 'CAC', 'CAA', 'CAG', 'AAT', 'AAC',
               'AAA', 'AAG', 'GAT', 'GAC', 'GAA', 'GAG', 'TGT',
               'TGC', 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG',
               'GGT', 'GGC', 'GGA', 'GGG', 'TGG', 'TAA', 'TAG',
               'TGA']

def scoreModels():
    codingMatrix = getProbs("./codingModel.tab")
    noncodingMatrix = getProbs("./noncodingModel.tab")

    ancestorSeq = getSeq("./Ancestor.fa")
    spaciiSeq = getSeq("./Spacii.fa")

    with open("Evanessa.fa", "w") as output_file:
        allIDs = list(ancestorSeq.keys())

        for ID in allIDs:
            cScore = 0
            nScore = 0

            ancestorSeg = ancestorSeq[ID]
            spaciiSeg = spaciiSeq[ID]

            for i in range(len(ancestorSeg) - 2):  # Adjusted range to ensure we don't go out of bounds
                ancestorCodon = ancestorSeg[i:i + 3]
                spaciiCodon = spaciiSeg[i:i + 3]

                if ancestorCodon in modelCodons and spaciiCodon in modelCodons:
                    codonIndex = modelCodons.index(ancestorCodon)
                    if codingMatrix[codonIndex][modelCodons.index(spaciiCodon)] != 0:
                        cScore += math.log(codingMatrix[codonIndex][modelCodons.index(spaciiCodon)])

                    if noncodingMatrix[codonIndex][modelCodons.index(spaciiCodon)] != 0:
                        nScore += math.log(noncodingMatrix[codonIndex][modelCodons.index(spaciiCodon)])

            if cScore > nScore:
                output_file.write(f">{ID} is coding ({cScore} vs {nScore})\n")
            else:
                output_file.write(f">{ID} is NOT coding ({cScore} vs {nScore})\n")


def getProbs(f1):
    with open(f1) as f:
        pMatrix = []
        for line in f:
            tmp = line.rstrip().split("\t")
            tmp = [float(i) for i in tmp]
            pMatrix.append(tmp)
    return pMatrix


def getSeq(filename):
    id2seq = {}
    with open(filename) as f:
        currkey = ""
        for line in f:
            if line.startswith(">"):
                currkey = line[1:].split("|")[0].strip()
                id2seq[currkey] = ""
            else:
                id2seq[currkey] += line.strip()
    return id2seq


scoreModels()
