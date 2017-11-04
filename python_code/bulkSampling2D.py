import numpy as np
import pickle
import sys
import pandas as pd
import time
import copy
import math
import os
from collections import Counter


def getRealTree(tum_id, bulk_id, phylogeny_read, sample_IDs, phylo = True):

    phylogeny = phylogeny_read
    sampledIDs = sample_IDs

    # choose to build a tree from all the cells (phylogeny) or from the sampled cells (sampledIDs) only
    if phylo:
        # quick loop over phylogeny to get tip nodes only and store as a first column in tree
        tree = []
        tip_labels = []
        tip_types = []
        for i in range(len(phylogeny.id)):
            if phylogeny.d1[i] == 0 and phylogeny.d2[i] == 0:
                tree.append([phylogeny.id[i]])
                #tree.append(['n'+str(phylogeny.id[i])] if phylogeny.type[i] == 1 else ['s'+str(phylogeny.id[i])])
                tip_types.append('n' if phylogeny.type[i] == 1 else 's')
                tip_labels.append(phylogeny.id[i])
    else:
        # quick loop over sampledIDs to get tip nodes only and store as a first column in tree
        tree = []
        tip_labels = []
        tip_types = []
        for i in range(len(sampledIDs)):
            if int(phylogeny[phylogeny.id == sampledIDs[i]].d1) == 0 and int(phylogeny[phylogeny.id == sampledIDs[i]].d2) == 0:
                tree.append([sampledIDs[i]])
                tip_types.append('n' if int(phylogeny[phylogeny.id == sampledIDs[i]].type) == 1 else 's')
                tip_labels.append(sampledIDs[i])

    # second loop over phylogeny (from bottom to up) to group internal nodes by ancestors
    start = time.time()
    for i in range(len(phylogeny.id)-1,-1,-1):
        print(i)
        if phylogeny.d1[i] > 0 or phylogeny.d2[i] > 0:
            for j in range(len(tree)):
                curBranch = pd.Series(tree[j])
                if sum(curBranch.isin([phylogeny.d1[i], phylogeny.d2[i]])) > 0:
                    tree[j].append(phylogeny.id[i])
    end = time.time()
    print('tree building finishied in: ', end - start)

    # get number of branching points per leaf from tree nested list
    generations = [len(x) for x in tree]
    maxGen = max(generations)

    # convert tree array into data.frame (dimension = tips vs generations)
    treeMatrix = np.zeros((len(tree), maxGen))
    for i in range(len(tree)):
        curLength = len(tree[i])
        if curLength== maxGen:
            # insert directly into the data.frame the reversed branch
            treeMatrix[i] = list(reversed(tree[i]))
        else:
            #if not enough nodes, frist replicate the node before the last node by maxGen - branchLength and then insert
            replicate = maxGen-curLength+1
            tree[i] = list(reversed(tree[i]))
            curInsert = [tree[i][0:curLength-2]]
            curInsert.append(replicate*[tree[i][curLength-2]])
            curInsert.append([tree[i][curLength-1]])
            curInsert = [x for y in curInsert for x in y]
            treeMatrix[i] = curInsert


    colNames = ["c" + str(i) for i in range(treeMatrix.shape[1])]

    treeMatrix = pd.DataFrame(treeMatrix, columns = colNames)

    treeMatrix.to_csv('treeMatrix{}_bulk{}.txt'.format(tum_id, bulk_id), sep = ' ', index = False)

    tip_labels = pd.DataFrame(np.vstack((tip_labels, tip_types)))
    tip_labels.to_csv('tip_labels{}_bulk{}.txt'.format(tum_id, bulk_id), sep = ' ', index = False)


#sampledIDs = allNonzeroIDs; phylogeny = phylogeny_read; cellMutations = cellMutations_read

def getSamples(sampledIDs, phylogeny, cellMutations):

    leavesList = []
    leavesMutList = []
    cellTypesList = []
    tipNodeTypeList = []
    for i in range(len(sampledIDs)):
        fatherID = int(sampledIDs[i])
        currType = 'n'+str(i) if int(phylogeny[phylogeny.id == fatherID].type) == 1 else 's'+str(i)
        tipNodeTypeList.append(currType)
        if int(phylogeny[phylogeny.id == fatherID].dead) == 0: #and phylogeny.d1[fatherID-1] == 0 and phylogeny.d2[fatherID-1] == 0:
            newLeaves = np.array([])
            newLeaves = np.append(newLeaves, fatherID)

            newMuts = np.array([])
            nextMut = 0 if (cellMutations[fatherID-1] is None or len(cellMutations[fatherID-1]) < 2) else list(cellMutations[fatherID-1][1])
            newMuts = np.append(newMuts, nextMut)

            #newTypes = np.array([])
            nextType = ['n'] if int(phylogeny[phylogeny.id == fatherID].type) == 1 else ['s']
            #newTypes = np.append(newTypes, nextType*len(nextMut))

            while (fatherID-1 != 0):

                newFatherID = int(phylogeny[phylogeny.id == fatherID].f)

                newLeaves = np.append(newLeaves, newFatherID)

                nextMut = 0 if (cellMutations[newFatherID-1] is None or len(cellMutations[newFatherID-1]) < 2) else list(cellMutations[newFatherID-1][1])
                newMuts = np.append(newMuts, nextMut)

                #nextType = ['n'] if int(phylogeny[phylogeny.id == fatherID].type) == 1 else ['s']
                #newTypes = np.append(newTypes, nextType*len(nextMut))

                fatherID = newFatherID

            newLeaves = newLeaves.astype(int).tolist()
            newMuts = newMuts.astype(int).tolist()
            #cellTypesList.append(newTypes)
            cellTypesList.append(len(newMuts)*nextType)
            leavesList.append(newLeaves)
            leavesMutList.append(newMuts)
        else:
            print(fatherID)

    return {'sampledLeaves':leavesList, 'sampledMutations':leavesMutList, 'sampledIDs':sampledIDs, 'sampledTypes':cellTypesList, 'tipNodeCellTypes':tipNodeTypeList}


def getUniqMutList(sampledMutList):

    #uniqMutList = [sampledMutList[y][0] for y in range(len(sampledMutList))]
    uniqMutList = np.array([item for sublist in sampledMutList for item in sublist])
    uniqMutList = np.unique(uniqMutList)

    return uniqMutList


def getMutations(uniqMutations, sampMutations, sampTypes, outFileName):

    result = np.zeros((uniqMutations.size, len(sampMutations)+1))
    result[:,0] = uniqMutations
    for i in range(uniqMutations.size):
        for j in range(len(sampMutations)):
            result[i][j+1] = 1 if uniqMutations[i] in set(sampMutations[j]) else 0

    #sampleNames = ["s" + str(i) for i in range(len(sampMutations))]
    sampleNames = sampTypes
    colNames = np.concatenate((["muts"], sampleNames))

    result = pd.DataFrame(result, columns = colNames)

    result = result.astype(int)

    result.to_csv(outFileName, sep = ' ', index = False)

    return result


def addNoise(mutData, seqDepth, minFreq):

    result = np.array([])

    for index, row in mutData.iterrows():

        samp_vaf = float(sum(np.random.binomial(1, p = .5*row['freq'], size = seqDepth)))

        if samp_vaf > minFreq:
            result = np.append(result, [row['count'],samp_vaf/seqDepth, row['mut_id'], row['clone']])

    result = np.reshape(result, (int(len(result)/4), 4))
    result = pd.DataFrame(result, columns = ['count', 'freq', 'mut_id', 'clone'])

    return result

def processBulkSampling(idsToBeSampled_in, phylogeny_read, cellMutations_read):

    # get bulk samples in a dictionary
    samples_bulk = getSamples(idsToBeSampled_in, phylogeny_read, cellMutations_read)

    # get sampled mutations list
    bulk_mutlist = samples_bulk['sampledMutations']

    # get sampled cell types
    bulk_typelist = samples_bulk['sampledTypes']

    # aggregate all mutations into one array list
    bulk_mutations = [item for sublist in bulk_mutlist for item in sublist]
    bulk_types = [item for sublist in bulk_typelist for item in sublist]

    bulk_mut_type = pd.DataFrame({'mut_id': bulk_mutations, 'clone': bulk_types})

    bulk_mut_type = bulk_mut_type.groupby(['mut_id', 'clone']).size().reset_index(name='count')

    bulk_mut_type['freq'] = round(bulk_mut_type['count']/len(idsToBeSampled_in),4)

    return(bulk_mut_type)


def sampSingleCells(samp_id, scIdsToBeSampled, addGenotyping, outFileName, cellIDs_read, phylogeny_read, cellMutations_read):

    singleCellIDs = scIdsToBeSampled.astype(int)

    singleCellSamples = getSamples(singleCellIDs, phylogeny_read, cellMutations_read)

    singleCellMutlist = singleCellSamples['sampledMutations']
    singleCellTypes = singleCellSamples['tipNodeCellTypes']

    # get mutations binary matrix for SC samples
    singleCellUniqMuts = getUniqMutList(singleCellMutlist)
    getMutations(singleCellUniqMuts, singleCellMutlist, singleCellTypes, 'mutations_SC'+outFileName+'_{}.txt'.format(samp_id))

    # export same SC cells for plateau analysis
    plateau_SC_maxlen = len(max(singleCellMutlist, key = len))

    plateau_SC_df = np.zeros(shape = (plateau_SC_maxlen, len(singleCellMutlist)))

    for i in range(len(singleCellMutlist)):
        for j in range(len(singleCellMutlist[i])):
            plateau_SC_df[j][i] = singleCellMutlist[i][j]

    plateau_SC_df = pd.DataFrame(plateau_SC_df, columns = singleCellTypes)
    plateau_SC_df.to_csv('plateau_SC'+outFileName+'_{}.txt'.format(samp_id), sep = ' ', index = False)

    if addGenotyping:

        genoCellIDs = np.random.choice(cellIDs_read[np.nonzero(cellIDs_read)], size = 100, replace = False)
        genoCellIDs = genoCellIDs.astype(int)

        genoSamples = getSamples(genoCellIDs, phylogeny_read, cellMutations_read)

        genoMutlist = genoSamples['sampledMutations']
        genoTypes = genoSamples['tipNodeCellTypes']

        #genoUniqMuts = getUniqMutList(genoMutlist)

        getMutations(singleCellUniqMuts, genoMutlist, genoTypes, 'mutations_SC_genotyped100'+outFileName+'_{}.txt'.format(samp_id))


def sampBulkCrossSectors(centre, radius, cellIDs_read, phylogeny_read, cellMutations_read, seqDepth, minFreq, samp_id, outFileName):


    m = cellIDs_read.shape[0]

    sampledCellCount = []
    for c in range(centre.shape[0]):

        # difining each region (building sides = twice the radius) and sampling uniformly nsample number of points from each
        centxEnd = m - 1 if centre[c][0] + radius > m - 1 else int(centre[c][0] + radius)
        centyEnd = m - 1 if centre[c][1] + radius > m - 1 else int(centre[c][1] + radius)

        centxStart = 0 if centre[c][0] - radius < 0 else int(centre[c][0] - radius)
        centyStart = 0 if centre[c][1] - radius < 0 else int(centre[c][1] - radius)

        #print(centxStart, centxEnd, centyStart, centyEnd)

        # collecting all the cell id-s from the above defined region
        idsToBeSampled = np.array([])

        for i in np.arange(centxStart, centxEnd):
            for j in np.arange(centyStart, centyEnd):
                idsToBeSampled = np.append(idsToBeSampled, cellIDs_read[i][j])

        # delete all non-positive(i.e. empty) id-s
        idsToBeSampled = idsToBeSampled[np.nonzero(idsToBeSampled)]
        idsToBeSampled = idsToBeSampled.astype(int)

        sampledCellCount.append(len(idsToBeSampled))

        samp_bulk_muts_freq = processBulkSampling(idsToBeSampled,  phylogeny_read, cellMutations_read)

        samp_bulk_muts_freq.to_csv(outFileName+'_{}_{}.txt'.format(samp_id, c), sep = ' ', index = False)

        samp_bulk_muts_freq_TEMP = copy.deepcopy(samp_bulk_muts_freq)
        samp_bulk_muts_freq_noise = addNoise(samp_bulk_muts_freq_TEMP, seqDepth, minFreq)

        samp_bulk_muts_freq_noise.to_csv(outFileName+'_noise_{}_{}.txt'.format(samp_id, c), sep = ' ', index = False)

    sampledCellCount = pd.DataFrame(sampledCellCount, columns = ['count'])
    sampledCellCount.to_csv('CellCount_'+outFileName+'_{}.txt'.format(samp_id), sep = ' ', index = False)


def sampBulkInnerOuter(cellIDs_read, phylogeny_read, cellMutations_read, seqDepth, minFreq, samp_id):

    sampledCellCount = []
    m = cellIDs_read.shape[0]

    allNonzeroIDs = cellIDs_read[np.nonzero(cellIDs_read)]

    innerIDsToBeSampled = np.array([])

    for i in np.arange(int(m/4), int(3*m/4)):
        for j in np.arange(int(m/4), int(3*m/4)):
            innerIDsToBeSampled = np.append(innerIDsToBeSampled, cellIDs_read[i][j])

    innerIDsToBeSampled = innerIDsToBeSampled[np.nonzero(innerIDsToBeSampled)]
    innerIDsToBeSampled = innerIDsToBeSampled.astype(int)

    outerIDsToBeSampled = [x for x in allNonzeroIDs if x not in innerIDsToBeSampled]

    sampledCellCount.append(len(innerIDsToBeSampled))
    sampledCellCount.append(len(outerIDsToBeSampled))

    inner_bulk_muts_freq = processBulkSampling(innerIDsToBeSampled, phylogeny_read, cellMutations_read)
    outer_bulk_muts_freq = processBulkSampling(outerIDsToBeSampled, phylogeny_read, cellMutations_read)

    inner_bulk_muts_freq.to_csv('innerBulk_{}.txt'.format(samp_id), sep = ' ', index = False)
    outer_bulk_muts_freq.to_csv('outerBulk_{}.txt'.format(samp_id), sep = ' ', index = False)

    inner_bulk_muts_freq_TEMP = copy.deepcopy(inner_bulk_muts_freq)
    inner_bulk_muts_freq_noise = addNoise(inner_bulk_muts_freq_TEMP, seqDepth, minFreq)
    inner_bulk_muts_freq_noise.to_csv('innerBulk_noise_{}.txt'.format(samp_id), sep = ' ', index = False)

    outer_bulk_muts_freq_TEMP = copy.deepcopy(outer_bulk_muts_freq)
    outer_bulk_muts_freq_noise = addNoise(outer_bulk_muts_freq_TEMP, seqDepth, minFreq)
    outer_bulk_muts_freq_noise.to_csv('outerBulk_noise_{}.txt'.format(samp_id), sep = ' ', index = False)

    sampledCellCount = pd.DataFrame(sampledCellCount, columns = ['count'])
    sampledCellCount.to_csv('CellCount_InnerOuter_{}.txt'.format(samp_id), sep = ' ', index = False)


def sampSectors(samp_id, seqDepth, minFreq, cellIDs_read, phylogeny_read, cellMutations_read):

    m = cellIDs_read.shape[0]

    singleCellIDs = np.array([])
    sampledCellCount = []

    grid = np.zeros((m,m))
    grid = grid.astype(int)

    radius = m
    startAngle = np.arange(0,360,45)
    percentange = 25

    for a, stAngle in enumerate(startAngle):
        print(stAngle)

        idsToBeSampled = np.array([])

        endAngle = 360/percentange + stAngle
        for i in np.arange(0, m):

            x = i - int(m/2)
            for j in np.arange(0, m):

                y = j - int(m/2)
                angle = math.degrees(np.arctan2(y, x))+180
                polarRadius = math.sqrt(x**2 + y**2)

                if angle >= stAngle and angle <= endAngle and polarRadius < radius:
                    idsToBeSampled = np.append(idsToBeSampled, cellIDs_read[i][j])
                    grid[i][j] = 3

        idsToBeSampled = idsToBeSampled[np.nonzero(idsToBeSampled)]
        idsToBeSampled = idsToBeSampled.astype(int)

        sampledCellCount.append(len(idsToBeSampled))

        samp_bulk_muts_freq = processBulkSampling(idsToBeSampled,  phylogeny_read, cellMutations_read)

        samp_bulk_muts_freq.to_csv('sectorsBulk2D_{}_{}.txt'.format(samp_id, a), sep = ' ', index = False)

        samp_bulk_muts_freq_TEMP = copy.deepcopy(samp_bulk_muts_freq)
        samp_bulk_muts_freq_noise = addNoise(samp_bulk_muts_freq_TEMP, seqDepth, minFreq)

        samp_bulk_muts_freq_noise.to_csv('sectorsBulk2D_noise'+str(seqDepth)+'_{}_{}.txt'.format(samp_id, a), sep = ' ', index = False)

        singleCellIDs = np.append(singleCellIDs, np.random.choice(idsToBeSampled, size = 1, replace = False))

    sampledCellCount = pd.DataFrame(sampledCellCount, columns = ['count'])
    sampledCellCount.to_csv('CellCount_SectorsBulk2D_{}.txt'.format(samp_id), sep = ' ', index = False)

    grid = pd.DataFrame(grid)
    grid.to_csv('sampled_grid2D_sectors_{}.txt'.format(samp_id), sep = ' ', index = False)

    sampSingleCells(samp_id, singleCellIDs, True, '_2Dsectors', cellIDs_read, phylogeny_read, cellMutations_read)



def sampStripes(samp_id, seqDepth, minFreq, cellIDs_read, phylogeny_read, cellMutations_read):

    m = cellIDs_read.shape[0]

    #xBase = [m/3, m/2, 2*m/3]
    #width = [2,1,2]

    xBase = [m/6, m/3, m/2, 2*m/3, 5*m/6]
    width = [m/20, m/23, m/27, m/23, m/20]

    #xStart = [int(i-m/20) for i in xBase]
    #xEnd = [int(i+m/20) for i in xBase]

    #singleCellIDs = np.array([])
    sampledCellCount = []

    grid = np.zeros((m,m))
    grid = grid.astype(int)

    for c, val in enumerate(xBase):

        idsToBeSampled = np.array([])

        #for i in np.arange(xStart[c], xEnd[c]):
        #for i in np.arange(int(val - width[c]/2), int(val + width[c]/2)):
        for i in np.arange(int(val - width[c]), int(val + width[c])):
            for j in np.arange(0,m):
                if cellIDs_read[i][j] > 0:
                #if phylogeny_read.dead[cellIDs_read[i][j]-1] == 0:

                    idsToBeSampled = np.append(idsToBeSampled, cellIDs_read[i][j])
                    grid[i][j] = 3

        idsToBeSampled = idsToBeSampled[np.nonzero(idsToBeSampled)]
        idsToBeSampled = idsToBeSampled.astype(int)

        sampledCellCount.append(len(idsToBeSampled))

        #getRealTree(tum_id = samp_id, bulk_id = c, phylogeny_read = phylogeny_read, sample_IDs = idsToBeSampled, phylo = False)
        sampSingleCells(samp_id, idsToBeSampled, False, '_2D_bulk'+str(c), cellIDs_read, phylogeny_read, cellMutations_read)

        samp_bulk_muts_freq = processBulkSampling(idsToBeSampled,  phylogeny_read, cellMutations_read)

        samp_bulk_muts_freq.to_csv('stripesBulk2D_{}_{}.txt'.format(samp_id, c), sep = ' ', index = False)

        samp_bulk_muts_freq_TEMP = copy.deepcopy(samp_bulk_muts_freq)
        samp_bulk_muts_freq_noise = addNoise(samp_bulk_muts_freq_TEMP, seqDepth, minFreq)

        samp_bulk_muts_freq_noise.to_csv('stripesBulk2D_noise'+str(seqDepth)+'_{}_{}.txt'.format(samp_id, c), sep = ' ', index = False)

        #singleCellIDs = np.append(singleCellIDs, np.random.choice(idsToBeSampled, size = 2, replace = False))

    sampledCellCount = pd.DataFrame(sampledCellCount, columns = ['count'])
    sampledCellCount.to_csv('CellCount_StripesBulk2D_{}.txt'.format(samp_id), sep = ' ', index = False)

    grid = pd.DataFrame(grid)
    grid.to_csv('sampled_grid2D_stripes_{}.txt'.format(samp_id), sep = ' ', index = False)

    #sampSingleCells(samp_id, singleCellIDs, True, '_2Dstripes', cellIDs_read, phylogeny_read, cellMutations_read)


def sampSquares(samp_id, seqDepth, minFreq, radius, cellIDs_read, phylogeny_read, cellMutations_read):

    # getting grid size
    m = cellIDs_read.shape[0]

    # MRS bulk and single cell sampling:
    x = y = np.array([m/4, m/2, 3*m/4])
    x = y = y.astype(int)

    centre = np.array([])

    for i in range(3):
        for j in range(3):
            centre = np.append(centre, [x[i], y[j]])


    # centre array contains the central coordinates for the 9 squares to be sampled from (both as bulk and SCs)
    centre = centre.reshape(9,2).astype(int)

    #singleCellIDs = np.array([])
    sampledCellCount = []

    grid = np.zeros((m,m))
    grid = grid.astype(int)

    for c, val in enumerate(centre):

        print(c, val)

        idsToBeSampled = np.array([])

        for i in np.arange(val[0]-radius, val[0]+radius):
            for j in np.arange(val[1]-radius, val[1]+radius):
                if cellIDs_read[i][j] > 0:
                    idsToBeSampled = np.append(idsToBeSampled, cellIDs_read[i][j])
                    grid[i][j] = 3

        idsToBeSampled = idsToBeSampled[np.nonzero(idsToBeSampled)]
        idsToBeSampled = idsToBeSampled.astype(int)

        sampledCellCount.append(len(idsToBeSampled))

        #getRealTree(tum_id = samp_id, bulk_id = c, phylogeny_read = phylogeny_read, sample_IDs = idsToBeSampled, phylo = False)

        treeCellIDs = np.random.choice(idsToBeSampled, size = 200, replace = False)

        sampSingleCells(samp_id, treeCellIDs, False, '_2D_bulk'+str(c), cellIDs_read, phylogeny_read, cellMutations_read)

        #bulk_muts_freq = processBulkSampling(idsToBeSampled, phylogeny_read, cellMutations_read)

        #bulk_muts_freq.to_csv('squaresBulk2D_radius{}_{}_{}.txt'.format(radius, samp_id, c), sep = ' ', index = False)

        #bulk_muts_freq_TEMP = copy.deepcopy(bulk_muts_freq)

        #bulk_muts_freq_noise = addNoise(bulk_muts_freq_TEMP, seqDepth, minFreq)

        #bulk_muts_freq_noise.to_csv('squaresBulk2D_noise'+str(seqDepth)+'_radius{}_{}_{}.txt'.format(radius, samp_id, c), sep = ' ', index = False)

        #singleCellIDs = np.append(singleCellIDs, np.random.choice(idsToBeSampled, size = 1, replace = False))

    sampledCellCount = pd.DataFrame(sampledCellCount, columns = ['count'])
    sampledCellCount.to_csv('CellCount_SquaresBulk2D_radius{}_{}.txt'.format(radius, samp_id), sep = ' ', index = False)

    grid = pd.DataFrame(grid)
    grid.to_csv('sampled_grid2D_squares_radius{}_{}.txt'.format(radius, samp_id), sep = ' ', index = False)

    #sampSingleCells(samp_id, singleCellIDs, True, '_2Dsquares', cellIDs_read, phylogeny_read, cellMutations_read)


def sampWholeBulk(samp_id, seqDepth, minFreq, cellIDsFile, phyloFile, cellMutFile):

    # reading data: cellIDs, phylogeny and cellMutations
    with open(cellIDsFile, 'rb') as cfile:
        cellIDs_read = pickle.load(cfile, encoding='latin1')

    with open(phyloFile, 'rb') as pfile:
        phylogeny_read = pickle.load(pfile, encoding='latin1')

    with open(cellMutFile, 'rb') as mfile:
        cellMutations_read = pickle.load(mfile, encoding='latin1')

    allNonzeroIDs = cellIDs_read[np.nonzero(cellIDs_read)]

    treeCellIDs = np.random.choice(allNonzeroIDs, size = 10, replace = False)

    sampSingleCells(samp_id, treeCellIDs, False, '_2D_whole', cellIDs_read, phylogeny_read, cellMutations_read)

    bulk_muts_freq = processBulkSampling(allNonzeroIDs, phylogeny_read, cellMutations_read)

    bulk_muts_freq.to_csv('wholeBulk_{}.txt'.format(samp_id), sep = ' ', index = False)

    bulk_muts_freq_TEMP = copy.deepcopy(bulk_muts_freq)

    bulk_muts_freq_noise = addNoise(bulk_muts_freq_TEMP, seqDepth, minFreq)

    bulk_muts_freq_noise.to_csv('wholeBulk_noise{}x_{}.txt'.format(seqDepth, samp_id), sep = ' ', index = False)


def runSampling(samp_id, seqDepth, minFreq, radius_squares, cellIDsFile, phyloFile, cellMutFile):
    # reading data: cellIDs, phylogeny and cellMutations
    with open(cellIDsFile, 'rb') as cfile:
        cellIDs_read = pickle.load(cfile, encoding='latin1')

    with open(phyloFile, 'rb') as pfile:
        phylogeny_read = pickle.load(pfile, encoding='latin1')

    with open(cellMutFile, 'rb') as mfile:
        cellMutations_read = pickle.load(mfile, encoding='latin1')

    allNonzeroIDs = cellIDs_read[np.nonzero(cellIDs_read)]; print(len(allNonzeroIDs))

    sampSingleCells(samp_id, np.random.choice(allNonzeroIDs, size = 20, replace = False), True, '_2D', cellIDs_read, phylogeny_read, cellMutations_read);

    #sampSingleCells(samp_id, allNonzeroIDs, False, '_2D_whole', cellIDs_read, phylogeny_read, cellMutations_read)


    # sample 9 squares
    #sampSquares(samp_id, seqDepth, minFreq, radius_squares, cellIDs_read, phylogeny_read, cellMutations_read)

    # sample 8 sectors
    #sampSectors(samp_id, seqDepth, minFreq, cellIDs_read, phylogeny_read, cellMutations_read)

    # sample 5 stripes
    #sampStripes(samp_id, seqDepth, minFreq, cellIDs_read, phylogeny_read, cellMutations_read)

    # getting grid size
    #m = cellIDs_read.shape[0]

    # direction: north -> east -> south -> west
    #crossCentres = np.vstack(([int(m/2)-1, int(3*m/4)-1], [int(3*m/4)-1, int(m/2)-1], [int(m/2)-1, int(m/4)-1], [[int(m/4)-1, int(m/2)-1]]))
    #sectorCentres = np.vstack(([int(3*m/4)-1, int(3*m/4)-1],[int(3*m/4)-1, int(m/4)-1], [int(m/4)-1, int(m/4)-1], [int(m/4)-1, int(3*m/4)-1]))

    # sample four bulks along cross
    #sampBulkCrossSectors(crossCentres, m/8, cellIDs_read, phylogeny_read, cellMutations_read, seqDepth, minFreq, samp_id, 'crossSquaresBulk')

    # sample four sectors
    #sampBulkCrossSectors(sectorCentres, m/4, cellIDs_read, phylogeny_read, cellMutations_read, seqDepth, minFreq, samp_id, 'entireSectorsBulk')

    # sample inner and outer (surrounding) bulks
    #sampBulkInnerOuter(cellIDs_read, phylogeny_read, cellMutations_read, seqDepth, minFreq, samp_id)

#sampWholeBulk(samp_id = sys.argv[1], seqDepth = int(sys.argv[2]), minFreq = float(sys.argv[3]), cellIDsFile = sys.argv[4], phyloFile = sys.argv[5], cellMutFile = sys.argv[6])
#runSampling(samp_id = sys.argv[1], seqDepth = int(sys.argv[2]), minFreq = float(sys.argv[3]), cellIDsFile = sys.argv[4], phyloFile = sys.argv[5], cellMutFile = sys.argv[6], radius_squares = int(sys.argv[7]))

#exec mutations_SC_genotyped100_2D_1.txt.nex; set increase=auto; hsearch start=stepwise addseq=random nreps=50 swap=TBR; outgroup normal; filter binary=yes; savetree file=mutations_SC_genotyped100_2D_1.txt.nex.tre brlen=yes;


samp_id = sid = 1; seqDepth = 200; minFreq = 3; bulkSize = 100;
cellIDsFile = 'cellIDs'+str(sid)+'_pickle.txt';
phyloFile = 'phylogeny'+str(sid)+'_pickle.txt';
cellMutFile = 'cellMutations'+str(sid)+'_pickle.txt';

# running the functions
start = time.time()
os.chdir('/Users/kchkhaidze/Documents/Projects/simulations/2D')
#ids=[x+10 for x in range(7)]
ids = [1]
for i in range(len(ids)):
    sid = ids[i]
    print(sid)
    sampWholeBulk(samp_id = sid, seqDepth = 200, minFreq = 3, cellIDsFile = 'cellIDs'+str(sid)+'_pickle.txt', phyloFile = 'phylogeny'+str(sid)+'_pickle.txt', cellMutFile = 'cellMutations'+str(sid)+'_pickle.txt')
    #runSampling(samp_id = sid, seqDepth = 100, minFreq = .01, radius_squares = 20, cellIDsFile = 'cellIDs'+str(sid)+'_pickle.txt', phyloFile = 'phylogeny'+str(sid)+'_pickle.txt', cellMutFile = 'cellMutations'+str(sid)+'_pickle.txt')
end = time.time()
print(end - start)
