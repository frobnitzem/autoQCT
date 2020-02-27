# Tools for running combinatorial search through a group
# of atoms in order to identify the included molecules.

import itertools
from functools import reduce
import os
import string
import molecule as mo


def clusterHydroxide(pdFrame, clusterSize, clusterLength):
    oxygenList = pdFrame.index[pdFrame['type'] == 'O'].tolist()
    hydrogenList = pdFrame.index[pdFrame['type'] == 'H'].tolist()

    hydrogenCombs = list(itertools.combinations(hydrogenList, 2))

    waterHydrogenIndex = []
    waterOxygenIndex = []
    waterOxygen = []
    waterHydrogen = []

    waterHydrogenPairIndex = []
    boolCount = []
    isWaterBool = []

    for oxygen in oxygenList:

        for n in range(len(hydrogenCombs)):
            atomH1 = pdFrame.loc[hydrogenCombs[n][0]]
            atomH2 = pdFrame.loc[hydrogenCombs[n][1]]
            atomO = pdFrame.loc[oxygen]

            if mo.isWater(atomH1, atomO, atomH2):
                isWaterBool.append(True)
                waterHydrogenIndex.append(hydrogenCombs[n][0])
                waterHydrogenIndex.append(hydrogenCombs[n][1])
                waterHydrogenPairIndex.append([hydrogenCombs[n][0], hydrogenCombs[n][1]])
                waterOxygenIndex.append(oxygen)

    if sum(isWaterBool) == (clusterSize - 2) / 3:
        boolCount.append(True)
    else:
        boolCount.append(False)

    hydroxideOxygenIndex = list(set(oxygenList) - set(waterOxygenIndex))
    hydroxideHydrogenIndex = list(set(hydrogenList) - set(waterHydrogenIndex))

    hydroxideHydrogen = pdFrame.loc[hydroxideHydrogenIndex[0]]
    hydroxideOxygen = pdFrame.loc[hydroxideOxygenIndex[0]]

    if mo.isHydroxide(hydroxideHydrogen, hydroxideOxygen):
        boolCount.append(True)
    else:
        boolCount.append(False)

    if mo.isClustered(hydroxideOxygen, waterHydrogen, clusterLength):
        boolCount.append(True)
    else:
        boolCount.append(False)

    for pair, single in zip(waterHydrogenPairIndex, waterOxygenIndex):
        waterOxygen.append(pdFrame.loc[single])
        waterHydrogen.append([pdFrame.loc[pair[0]], pdFrame.loc[pair[1]]])

    return [reduce((lambda x, y: x * y), boolCount), waterHydrogen, waterOxygen, hydroxideOxygen, hydroxideHydrogen]


def printClusterHydroxide(waterHydrogen, waterOxygen, hydroxideOxygen, hydroxideHydrogen, cluster, itr):
    if not os.path.exists("full"):
        os.mkdir("full")

    fullFileName = str("full/" + str(itr) + ".com")
    file = open(fullFileName, "w")
    file.write("t = " + str(cluster) + "\n")
    file.write('FULL CLUSTER' + "\n")

    for m in range(len(waterOxygen)):
        atomH1 = waterHydrogen[m][0]
        atomH2 = waterHydrogen[m][1]
        atomO = waterOxygen[m]
        file.write(str(atomO.type + "\t" + str(atomO.x) + "\t" + str(atomO.y) + "\t" + str(atomO.z) + "\n"))
        file.write(str(atomH1.type + "\t" + str(atomH1.x) + "\t" + str(atomH1.y) + "\t" + str(atomH1.z) + "\n"))
        file.write(str(atomH2.type + "\t" + str(atomH2.x) + "\t" + str(atomH2.y) + "\t" + str(atomH2.z) + "\n"))

    file.write(str(hydroxideOxygen.type + "\t" + str(hydroxideOxygen.x) + "\t" + str(hydroxideOxygen.y) + "\t" + str(
        hydroxideOxygen.z) + "\n"))
    file.write(str(
        hydroxideHydrogen.type + "\t" + str(hydroxideHydrogen.x) + "\t" + str(hydroxideHydrogen.y) + "\t" + str(
            hydroxideHydrogen.z) + "\n"))
    file.write("")
    file.close()

    clusterConfigurations = list(string.ascii_uppercase[0:len(waterOxygen)])

    ####### HYDROXIDE #######
    if not os.path.exists("OH"):
        os.mkdir("OH")

    OHFileName = str("OH/" + str(itr) + ".com")
    file = open(OHFileName, "w")
    file.write("t = " + str(cluster) + "\n")
    file.write("OH CLUSTER" + "\n")
    file.write(str(hydroxideOxygen.type + "\t" + str(hydroxideOxygen.x) + "\t" + str(hydroxideOxygen.y) + "\t" + str(
        hydroxideOxygen.z) + "\n"))
    file.write(str(
        hydroxideHydrogen.type + "\t" + str(hydroxideHydrogen.x) + "\t" + str(hydroxideHydrogen.y) + "\t" + str(
            hydroxideHydrogen.z) + "\n"))
    file.write("")
    file.close()

    ###### 1 Cluster ######
    for m in range(len(waterOxygen)):
        if not os.path.exists(clusterConfigurations[m]):
            os.mkdir(clusterConfigurations[m])

        FileName1 = clusterConfigurations[m] + "/" + str(itr) + ".com"
        file = open(FileName1, "w")
        file.write("t = " + str(cluster) + "\n")
        file.write(clusterConfigurations[m] + " Cluster" + "\n")
        atomH1 = waterHydrogen[m][0]
        atomH2 = waterHydrogen[m][1]
        atomO = waterOxygen[m]
        file.write(str(atomO.type + "\t" + str(atomO.x) + "\t" + str(atomO.y) + "\t" + str(atomO.z) + "\n"))
        file.write(str(atomH1.type + "\t" + str(atomH1.x) + "\t" + str(atomH1.y) + "\t" + str(atomH1.z) + "\n"))
        file.write(str(atomH2.type + "\t" + str(atomH2.x) + "\t" + str(atomH2.y) + "\t" + str(atomH2.z) + "\n"))
        file.write(
            str(hydroxideOxygen.type + "\t" + str(hydroxideOxygen.x) + "\t" + str(hydroxideOxygen.y) + "\t" + str(
                hydroxideOxygen.z) + "\n"))
        file.write(str(
            hydroxideHydrogen.type + "\t" + str(hydroxideHydrogen.x) + "\t" + str(hydroxideHydrogen.y) + "\t" + str(
                hydroxideHydrogen.z) + "\n"))
        file.write("")
        file.close()

        if (len(waterHydrogen) > 2):
            Nminus1 = [x for i, x in enumerate(clusterConfigurations) if i != m]
            waterHydrogen_N1 = [x for i, x in enumerate(waterHydrogen) if i != m]
            waterOxygen_N1 = [x for i, x in enumerate(waterOxygen) if i != m]

            ####### N - 1
            Nminus1 = ''.join(Nminus1)
            if not os.path.exists(Nminus1):
                os.mkdir(Nminus1)

            FileNameN1 = Nminus1 + "/" + str(itr) + ".com"
            file = open(FileNameN1, "w")
            file.write("t = " + str(cluster) + "\n")
            file.write(Nminus1 + " Cluster" + "\n")

            for l in range(len(waterHydrogen_N1)):
                atomH1_N1 = waterHydrogen_N1[l][0]
                atomH2_N1 = waterHydrogen_N1[l][1]
                atomO_N1 = waterOxygen_N1[l]
                file.write(
                    str(atomO_N1.type + "\t" + str(atomO_N1.x) + "\t" + str(atomO_N1.y) + "\t" + str(
                        atomO_N1.z) + "\n"))
                file.write(str(
                    atomH1_N1.type + "\t" + str(atomH1_N1.x) + "\t" + str(atomH1_N1.y) + "\t" + str(
                        atomH1_N1.z) + "\n"))
                file.write(str(
                    atomH2_N1.type + "\t" + str(atomH2_N1.x) + "\t" + str(atomH2_N1.y) + "\t" + str(
                        atomH2_N1.z) + "\n"))

            file.write(
                str(hydroxideOxygen.type + "\t" + str(hydroxideOxygen.x) + "\t" + str(hydroxideOxygen.y) + "\t" + str(
                    hydroxideOxygen.z) + "\n"))
            file.write(str(
                hydroxideHydrogen.type + "\t" + str(hydroxideHydrogen.x) + "\t" + str(hydroxideHydrogen.y) + "\t" + str(
                    hydroxideHydrogen.z) + "\n"))
            file.write("")
            file.close()

    return None


def clusterIon(pdFrame, clusterSize, clusterLength, ion):
    oxygenList = pdFrame.index[pdFrame['type'] == 'O'].tolist()
    hydrogenList = pdFrame.index[pdFrame['type'] == 'H'].tolist()

    ionList = pdFrame.index[pdFrame['type'] == ion].tolist()

    hydrogenCombs = list(itertools.combinations(hydrogenList, 2))

    waterHydrogenIndex = []
    waterOxygenIndex = []
    waterOxygen = []
    waterHydrogen = []

    waterHydrogenPairIndex = []
    boolCount = []
    isWaterBool = []

    for oxygen in oxygenList:

        for n in range(len(hydrogenCombs)):
            atomH1 = pdFrame.loc[hydrogenCombs[n][0]]
            atomH2 = pdFrame.loc[hydrogenCombs[n][1]]
            atomO = pdFrame.loc[oxygen]

            if mo.isWater(atomH1, atomO, atomH2):
                isWaterBool.append(True)
                waterHydrogenIndex.append(hydrogenCombs[n][0])
                waterHydrogenIndex.append(hydrogenCombs[n][1])
                waterHydrogenPairIndex.append([hydrogenCombs[n][0], hydrogenCombs[n][1]])
                waterOxygenIndex.append(oxygen)

    if sum(isWaterBool) == (clusterSize - 1) / 3:
        boolCount.append(True)
    else:
        boolCount.append(False)

    for pair, single in zip(waterHydrogenPairIndex, waterOxygenIndex):
        waterOxygen.append(pdFrame.loc[single])
        waterHydrogen.append([pdFrame.loc[pair[0]], pdFrame.loc[pair[1]]])

    ionX = pdFrame.loc[ionList[0]]

    if mo.isClustered(ionX, waterHydrogen, clusterLength):
        boolCount.append(True)
    else:
        boolCount.append(False)

    return [reduce((lambda x, y: x * y), boolCount), waterHydrogen, waterOxygen, ionX]


def printClusterIon(waterHydrogen, waterOxygen, ionFrame, cluster, itr):
    if not os.path.exists("full"):
        os.mkdir("full")

    fullFileName = str("full/" + str(itr) + ".com")
    file = open(fullFileName, "w")
    file.write("t = " + str(cluster) + "\n")
    file.write('FULL CLUSTER' + "\n")

    for m in range(len(waterOxygen)):
        atomH1 = waterHydrogen[m][0]
        atomH2 = waterHydrogen[m][1]
        atomO = waterOxygen[m]
        file.write(str(atomO.type + "\t" + str(atomO.x) + "\t" + str(atomO.y) + "\t" + str(atomO.z) + "\n"))
        file.write(str(atomH1.type + "\t" + str(atomH1.x) + "\t" + str(atomH1.y) + "\t" + str(atomH1.z) + "\n"))
        file.write(str(atomH2.type + "\t" + str(atomH2.x) + "\t" + str(atomH2.y) + "\t" + str(atomH2.z) + "\n"))

    file.write(str(ionFrame.type + "\t" + str(ionFrame.x) + "\t" + str(ionFrame.y) + "\t" + str(ionFrame.z) + "\n"))
    file.write("")
    file.close()

    clusterConfigurations = list(string.ascii_uppercase[0:len(waterOxygen)])

    # ####### ION #######
    # SINGLE ION SHOULDN"T BE PRINTED
    # if not os.path.exists("OH"):
    #     os.mkdir("OH")
    #
    # OHFileName = str("OH/" + str(itr) + ".com")
    # file = open(OHFileName, "w")
    # file.write("t = " + str(cluster) + "\n")
    # file.write("I CLUSTER" + "\n")
    # file.write(str(ionFrame.type + "\t" + str(ionFrame.x) + "\t" + str(ionFrame.y) + "\t" + str(
    #     ionFrame.z) + "\n"))
    # file.write("")
    # file.close()

    ###### 1 Cluster ######
    for m in range(len(waterOxygen)):
        if not os.path.exists(clusterConfigurations[m]):
            os.mkdir(clusterConfigurations[m])

        FileName1 = clusterConfigurations[m] + "/" + str(itr) + ".com"
        file = open(FileName1, "w")
        file.write("t = " + str(cluster) + "\n")
        file.write(clusterConfigurations[m] + " Cluster" + "\n")
        atomH1 = waterHydrogen[m][0]
        atomH2 = waterHydrogen[m][1]
        atomO = waterOxygen[m]
        file.write(str(atomO.type + "\t" + str(atomO.x) + "\t" + str(atomO.y) + "\t" + str(atomO.z) + "\n"))
        file.write(str(atomH1.type + "\t" + str(atomH1.x) + "\t" + str(atomH1.y) + "\t" + str(atomH1.z) + "\n"))
        file.write(str(atomH2.type + "\t" + str(atomH2.x) + "\t" + str(atomH2.y) + "\t" + str(atomH2.z) + "\n"))
        file.write(str(ionFrame.type + "\t" + str(ionFrame.x) + "\t" + str(ionFrame.y) + "\t" + str(
            ionFrame.z) + "\n"))
        file.write("")
        file.close()

        if (len(waterHydrogen) > 2):
            Nminus1 = [x for i, x in enumerate(clusterConfigurations) if i != m]
            waterHydrogen_N1 = [x for i, x in enumerate(waterHydrogen) if i != m]
            waterOxygen_N1 = [x for i, x in enumerate(waterOxygen) if i != m]

            ####### N - 1
            Nminus1 = ''.join(Nminus1)
            if not os.path.exists(Nminus1):
                os.mkdir(Nminus1)

            FileNameN1 = Nminus1 + "/" + str(itr) + ".com"
            file = open(FileNameN1, "w")
            file.write("t = " + str(cluster) + "\n")
            file.write(Nminus1 + " Cluster" + "\n")

            for l in range(len(waterHydrogen_N1)):
                atomH1_N1 = waterHydrogen_N1[l][0]
                atomH2_N1 = waterHydrogen_N1[l][1]
                atomO_N1 = waterOxygen_N1[l]
                file.write(
                    str(atomO_N1.type + "\t" + str(atomO_N1.x) + "\t" + str(atomO_N1.y) + "\t" + str(
                        atomO_N1.z) + "\n"))
                file.write(str(
                    atomH1_N1.type + "\t" + str(atomH1_N1.x) + "\t" + str(atomH1_N1.y) + "\t" + str(
                        atomH1_N1.z) + "\n"))
                file.write(str(
                    atomH2_N1.type + "\t" + str(atomH2_N1.x) + "\t" + str(atomH2_N1.y) + "\t" + str(
                        atomH2_N1.z) + "\n"))

            file.write(
                str(ionFrame.type + "\t" + str(ionFrame.x) + "\t" + str(ionFrame.y) + "\t" + str(
                    ionFrame.z) + "\n"))
            file.write("")
            file.close()

    return None


def clusterProton(pdFrame, clusterSize, clusterLength):
    oxygenList = pdFrame.index[pdFrame['type'] == 'O'].tolist()
    hydrogenList = pdFrame.index[pdFrame['type'] == 'H'].tolist()

    hydrogenCombs = list(itertools.combinations(hydrogenList, 2))

    waterHydrogenIndex = []
    waterOxygenIndex = []
    waterOxygen = []
    waterHydrogen = []

    waterHydrogenPairIndex = []
    boolCount = []
    isWaterBool = []

    for oxygen in oxygenList:

        for n in range(len(hydrogenCombs)):
            atomH1 = pdFrame.loc[hydrogenCombs[n][0]]
            atomH2 = pdFrame.loc[hydrogenCombs[n][1]]
            atomO = pdFrame.loc[oxygen]

            if mo.isWater(atomH1, atomO, atomH2):
                isWaterBool.append(True)
                waterHydrogenIndex.append(hydrogenCombs[n][0])
                waterHydrogenIndex.append(hydrogenCombs[n][1])
                waterHydrogenPairIndex.append([hydrogenCombs[n][0], hydrogenCombs[n][1]])
                waterOxygenIndex.append(oxygen)

    if sum(isWaterBool) == (clusterSize - 1) / 3:
        boolCount.append(True)
    else:
        boolCount.append(False)

    protonHydrogenIndex = list(set(hydrogenList) - set(waterHydrogenIndex))

    protonHydrogen = pdFrame.loc[protonHydrogenIndex[0]]

    if mo.isProton(protonHydrogen):
        boolCount.append(True)
    else:
        boolCount.append(False)

    if mo.isClustered(protonHydrogen, waterHydrogen, clusterLength):
        boolCount.append(True)
    else:
        boolCount.append(False)

    for pair, single in zip(waterHydrogenPairIndex, waterOxygenIndex):
        waterOxygen.append(pdFrame.loc[single])
        waterHydrogen.append([pdFrame.loc[pair[0]], pdFrame.loc[pair[1]]])

    return [reduce((lambda x, y: x * y), boolCount), waterHydrogen, waterOxygen, protonHydrogen]


def printClusterProton(waterHydrogen, waterOxygen, protonHydrogen, cluster, itr):
    if not os.path.exists("full"):
        os.mkdir("full")

    fullFileName = str("full/" + str(itr) + ".com")
    file = open(fullFileName, "w")
    file.write("t = " + str(cluster) + "\n")
    file.write('FULL CLUSTER' + "\n")

    for m in range(len(waterOxygen)):
        atomH1 = waterHydrogen[m][0]
        atomH2 = waterHydrogen[m][1]
        atomO = waterOxygen[m]
        file.write(str(atomO.type + "\t" + str(atomO.x) + "\t" + str(atomO.y) + "\t" + str(atomO.z) + "\n"))
        file.write(str(atomH1.type + "\t" + str(atomH1.x) + "\t" + str(atomH1.y) + "\t" + str(atomH1.z) + "\n"))
        file.write(str(atomH2.type + "\t" + str(atomH2.x) + "\t" + str(atomH2.y) + "\t" + str(atomH2.z) + "\n"))

    file.write(str(protonHydrogen.type + "\t" + str(protonHydrogen.x) + "\t" + str(protonHydrogen.y) + "\t" + str(
        protonHydrogen.z) + "\n"))
    file.write("")
    file.close()

    clusterConfigurations = list(string.ascii_uppercase[0:len(waterOxygen)])

    # ####### HYDROXIDE #######
    # SINGLE ION SHOULDN"T BE PRINTED
    # if not os.path.exists("OH"):
    #     os.mkdir("OH")
    #
    # OHFileName = str("OH/" + str(itr) + ".com")
    # file = open(OHFileName, "w")
    # file.write("t = " + str(cluster) + "\n")
    # file.write("OH CLUSTER" + "\n")
    # file.write(str(protonHydrogen.type + "\t" + str(protonHydrogen.x) + "\t" + str(protonHydrogen.y) + "\t" + str(
    #     protonHydrogen.z) + "\n"))
    # file.write("")
    # file.close()

    ###### 1 Cluster ######
    for m in range(len(waterOxygen)):
        if not os.path.exists(clusterConfigurations[m]):
            os.mkdir(clusterConfigurations[m])

        FileName1 = clusterConfigurations[m] + "/" + str(itr) + ".com"
        file = open(FileName1, "w")
        file.write("t = " + str(cluster) + "\n")
        file.write(clusterConfigurations[m] + " Cluster" + "\n")
        atomH1 = waterHydrogen[m][0]
        atomH2 = waterHydrogen[m][1]
        atomO = waterOxygen[m]
        file.write(str(atomO.type + "\t" + str(atomO.x) + "\t" + str(atomO.y) + "\t" + str(atomO.z) + "\n"))
        file.write(str(atomH1.type + "\t" + str(atomH1.x) + "\t" + str(atomH1.y) + "\t" + str(atomH1.z) + "\n"))
        file.write(str(atomH2.type + "\t" + str(atomH2.x) + "\t" + str(atomH2.y) + "\t" + str(atomH2.z) + "\n"))
        file.write(str(protonHydrogen.type + "\t" + str(protonHydrogen.x) + "\t" + str(protonHydrogen.y) + "\t" + str(
            protonHydrogen.z) + "\n"))
        file.write("")
        file.close()

        if (len(waterHydrogen) > 2):
            Nminus1 = [x for i, x in enumerate(clusterConfigurations) if i != m]
            waterHydrogen_N1 = [x for i, x in enumerate(waterHydrogen) if i != m]
            waterOxygen_N1 = [x for i, x in enumerate(waterOxygen) if i != m]

            ####### N - 1
            Nminus1 = ''.join(Nminus1)
            if not os.path.exists(Nminus1):
                os.mkdir(Nminus1)

            FileNameN1 = Nminus1 + "/" + str(itr) + ".com"
            file = open(FileNameN1, "w")
            file.write("t = " + str(cluster) + "\n")
            file.write(Nminus1 + " Cluster" + "\n")

            for l in range(len(waterHydrogen_N1)):
                atomH1_N1 = waterHydrogen_N1[l][0]
                atomH2_N1 = waterHydrogen_N1[l][1]
                atomO_N1 = waterOxygen_N1[l]
                file.write(
                    str(atomO_N1.type + "\t" + str(atomO_N1.x) + "\t" + str(atomO_N1.y) + "\t" + str(
                        atomO_N1.z) + "\n"))
                file.write(str(
                    atomH1_N1.type + "\t" + str(atomH1_N1.x) + "\t" + str(atomH1_N1.y) + "\t" + str(
                        atomH1_N1.z) + "\n"))
                file.write(str(
                    atomH2_N1.type + "\t" + str(atomH2_N1.x) + "\t" + str(atomH2_N1.y) + "\t" + str(
                        atomH2_N1.z) + "\n"))

            file.write(
                str(protonHydrogen.type + "\t" + str(protonHydrogen.x) + "\t" + str(protonHydrogen.y) + "\t" + str(
                    protonHydrogen.z) + "\n"))
            file.write("")
            file.close()

    return None
