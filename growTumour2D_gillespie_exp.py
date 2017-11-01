import numpy as np
import pandas as pd
import math
import pickle
import sys


def simulateTumorGrowth2D(sd, grid_size, mu, wildBirthRate, mutantBirthRate, alpha, beta, cellAggression, cloneStartTime):

	## defining variables
	# generation count
	generation = 1

	# 2D grid to store cells' unique IDs
	cellIDs = np.zeros((grid_size, grid_size))
	cellIDs = cellIDs.astype(int)

	# 2D grid of 0/1 indicating the state of each cell - dead/alive
	cellDeadLive = np.zeros((grid_size, grid_size))
	cellDeadLive = cellDeadLive.astype(int)

	# 2D grid to store cell types - 1/2 wild/mutant
	liveCellTypes = np.zeros((grid_size, grid_size))
	liveCellTypes = liveCellTypes.astype(int)

	# 2D matrix to store phylogeny of cells (columns: id, father, daughter1, daughter2, dead/alive, cellType)
	phylogeny = np.zeros((1, 6))
	phylogeny = phylogeny.astype(int)
	phylogeny = pd.DataFrame(phylogeny, columns = ['id', 'f', 'd1','d2','dead', 'type'])

	# list to store lists of mutations for each cell
	cellMutations = []

	# 2D matrix to store generation id and total mutation per generation
	globalMutations = np.zeros((1,2))
	globalMutations = globalMutations.astype(int)

	### initializing variables ###
	centre = int(np.floor((grid_size+1)/2))

	# register first wild type cell
	firstCellID = 1
	cellIDs[centre-1][centre-1] = firstCellID
	cellDeadLive[centre-1][centre-1] = 1  
	liveCellTypes[centre-1][centre-1] = 1 

	# her corresponding entries in the phylogeny 
	phylogeny.id[0] = 1

    # start with zero number of mutations
	totalMutCount = 0

	# generate random number of mutations for the first cell
	mutNum = np.random.poisson(mu, 1)[0]
	cellMutations.append([firstCellID])
	cellMutations[firstCellID-1].append([totalMutCount + k for k in range(mutNum)])
	totalMutCount += mutNum
	
	totalIdCount = 1

	# initialize the condition for the while loop to commence and be running unless the condition changes
	finish = 0

    # store total live cells per generation to plot growth curve
	totalLiveCells = [1]
    
    # store time steps passed per generation
	time = [0]
	timeIncrement = 0

	converted = False
    
    #random.seed(sd)
	while(finish != 1):

		# increase generation count (generation is one cell division event i.e. one reaction)
		generation += 1
		# count total number of cells in wild type
		wildCellCount = liveCellTypes[liveCellTypes == 1 ].size
		
		# get next time interval for wild reaction
		wildReaction = np.log(1/np.random.uniform())/(wildCellCount*wildBirthRate)

		
		if time[len(time)-1] >= cloneStartTime and not converted:

			print('converting: blue ----------> red')

			firstMutantCellID = np.random.choice(cellIDs[np.nonzero(cellIDs)], size = 1, replace = False)

			print(firstMutantCellID)

			liveCellTypes[np.where(cellIDs == firstMutantCellID[0])] = 2

			phylogeny.loc[phylogeny.id == firstMutantCellID[0], 'type'] = 2

			converted = True


		if converted:

			# count total number of cells in mutant type
			mutantCellCount = liveCellTypes[liveCellTypes == 2].size

			# get next time interval for mutant reaction
			mutantReaction = np.log(1/np.random.uniform())/(mutantCellCount*mutantBirthRate)

			# the reaction with the smallest time interval is chosen for the next step
			if min(wildReaction, mutantReaction) == wildReaction:
				cellType = 1
				timeIncrement = wildReaction
			else:
				cellType = 2
				timeIncrement = mutantReaction
		else:
			
			cellType = 1
			timeIncrement = wildReaction

		

        # increment time
		time.append(time[len(time)-1] + timeIncrement)

		# pick up uniformly which particular cell to divide next within a chosen reaction 
		curTypeCellCoords = np.where(liveCellTypes == cellType)
		curTypeCellCount = curTypeCellCoords[0].size

		if curTypeCellCount > 0:

			whichCoord = np.random.choice(curTypeCellCount)

			# set next cell's i.e. next father's coordinates
			nextFatherX = curTypeCellCoords[0][whichCoord]
			nextFatherY = curTypeCellCoords[1][whichCoord]

		    # and get he's ID from cellIDs grid
			fatherID = cellIDs[nextFatherX][nextFatherY]

			# check if both daughter cells die or at least one survives
			if np.random.uniform(0,1,1)[0] > alpha:
				# at least one daughter cell survived
				# store father cell's coordinates for the first daughter
				xtemp = np.where(cellIDs == fatherID)[0]
				ytemp = np.where(cellIDs == fatherID)[1]

				xtemp = xtemp[0]
				ytemp = ytemp[0]
				 
				# store neighbourhood coordinates of the father cell
				neighbourCoords = sum([[xtemp, ytemp, 0,\
										xtemp + 1, ytemp, 0,\
										xtemp - 1, ytemp, 0,\

										xtemp, ytemp - 1, 0,\
										xtemp + 1, ytemp - 1, 0,\
										xtemp - 1, ytemp - 1, 0,\
										
										xtemp, ytemp + 1, 0,\
										xtemp + 1, ytemp + 1, 0,\
										xtemp - 1, ytemp + 1, 0]], [])

				neighbourCoords = np.array(neighbourCoords)
				neighbourCoords.shape = (9, 3)
				freeCoords = np.zeros((1,3))
				neighbourFullSpace = 0
				spaceCorrection = 0

		        # count how many cells in the father cell's neighbourhood are occupied or are on the boundary of the initial 3D grid
				for j in range(neighbourCoords.shape[0]):
					if (neighbourCoords[j][0] >= grid_size or neighbourCoords[j][1] >= grid_size)\
					 or (neighbourCoords[j][0] < 0 or neighbourCoords[j][1] < 0):
						spaceCorrection += 1
						next
					else:
						if cellDeadLive[neighbourCoords[j][0]][neighbourCoords[j][1]] == 0:
							freeCoords = np.vstack([freeCoords, [neighbourCoords[j][0], neighbourCoords[j][1], neighbourCoords[j][2]]])
						else:
							neighbourFullSpace += 1
		        
		            # remove the empty first row (needed for vstack)
				freeCoords = np.delete(freeCoords, 0, 0)

				# if there is a space in the father cell neighbourhood: 
				if neighbourFullSpace < neighbourCoords.shape[0] - spaceCorrection:
					# check if any of the spaces from squares and diagonals are empty (ie not in the freeCoords array)
					pickedSpace = np.random.permutation(len(freeCoords))[0]
					xtest = int(freeCoords[pickedSpace][0])
					ytest = int(freeCoords[pickedSpace][1])

                      # create IDs for the daughter cells
					currD1Id = totalIdCount + 1
					currD2Id = totalIdCount + 2

					# store daughters' info in the phylogeny matrix
					phylogeny.loc[phylogeny.id == fatherID, 'd1'] = currD1Id
					phylogeny.loc[phylogeny.id == fatherID, 'd2'] = currD2Id
					phylogeny.loc[phylogeny.id == fatherID, 'type'] = cellType
					phylogeny = phylogeny.append(pd.DataFrame([[currD1Id, fatherID, 0, 0, 0, cellType], [currD2Id, fatherID, 0, 0, 0, cellType]], columns = ['id', 'f', 'd1','d2','dead', 'type']))
					
					# store both daughters' IDs, their states and types in the corresponding 3D grids
					# first daughter is registered over the father's space
					cellIDs[xtemp][ytemp] = currD1Id
					cellDeadLive[xtemp][ytemp] = 1
					liveCellTypes[xtemp][ytemp] = 1 if cellType == 1 else 2

					# second daughter is registered on the randomly chosen free space in a father cell neighbourhood
					cellIDs[xtest][ytest] = currD2Id
					cellDeadLive[xtest][ytest]= 1
					liveCellTypes[xtest][ytest] = 1 if cellType == 1 else 2

					# and increase the total ID count by two
					totalIdCount += 2

					# generate and store random number of mutations for the first daughter
					mutNum = np.random.poisson(mu, 1)[0]
					cellMutations.append([currD1Id])
					cellMutations[currD1Id-1].append([totalMutCount + k for k in range(mutNum)])
					totalMutCount += mutNum

					# check if the second daughter cell dies with probability = beta(local death)
					deathTest = np.random.uniform(0,1,1)[0]
					if deathTest < (1-beta):
						# dies: erase the second daughter from all 2D grids and phylogeny matrix
						cellIDs[xtest][ytest] = 0
						cellDeadLive[xtest][ytest] = 0
						liveCellTypes[xtest][ytest] = 0
						phylogeny.loc[phylogeny.id == currD2Id, 'dead'] = -1
						xtest = xtemp
						ytest = ytemp
						cellMutations.append([currD2Id])
						cellMutations[currD2Id-1] = None
					else:
						# survives: generate and store random number of mutations for the second daughter
						mutNum = np.random.poisson(mu, 1)[0]
						cellMutations.append([currD2Id])
						cellMutations[currD2Id-1].append([totalMutCount + k for k in range(mutNum)])
						totalMutCount += mutNum
				else:
					# if there is no free space in the father cell's neighbourhood, with the probability cellAggression shift all the cells 
					# from father cell till the first  empty space going in a random direction by one unit towards the grid bourders 
					aggression = np.random.uniform(0,1,1)[0]
					if aggression < cellAggression:
						x = 2*np.random.uniform(0,1,1)[0] - 1
						y = 2*np.random.uniform(0,1,1)[0] - 1
						r = math.sqrt(x**2 + y**2)
						dx = x/r
						dy = y/r
						shift = 1
						if round(xtemp + shift*dx) >= grid_size or round(ytemp + shift*dy) >= grid_size\
                                 or round(xtemp + shift*dx) < 0 or round(ytemp + shift*dy) < 0:
							shift = None
						else:
							while cellDeadLive[int(round(xtemp + shift*dx))][int(round(ytemp + shift*dy))]!= 0:
								shift += 1
								if round(xtemp + shift*dx) >= grid_size or round(ytemp + shift*dy) >= grid_size\
								or round(xtemp + shift*dx) < 0 or round(ytemp + shift*dy) < 0:
									shift = None
									break
						if shift is not None:
							for s in range(shift, 1, -1):
								if round(xtemp + s*dx) < grid_size and round(ytemp + s*dy) < grid_size\
								and round(xtemp + s*dx) >= 0 and round(ytemp + s*dy) >= 0 :
									if round(xtemp + (s-1)*dx) < grid_size and round(ytemp + (s-1)*dy) < grid_size\
									and round(xtemp + (s-1)*dx) >= 0 and round(ytemp + (s-1)*dy) >= 0:
										cellIDs[int(round(xtemp + s*dx))][int(round(ytemp + s*dy))] = cellIDs[int(round(xtemp + (s-1)*dx))][int(round(ytemp + (s-1)*dy))]
										cellDeadLive[int(round(xtemp + s*dx))][int(round(ytemp + s*dy))] = cellDeadLive[int(round(xtemp + (s-1)*dx))][int(round(ytemp + (s-1)*dy))]
										liveCellTypes[int(round(xtemp + s*dx))][int(round(ytemp + s*dy))] = liveCellTypes[int(round(xtemp + (s-1)*dx))][int(round(ytemp + (s-1)*dy))]
							
							# if shifting finishes successfully, take over the space of the 'first' cell (meaning closest to the father cell) shifted 
							xtest = int(round(xtemp + shift*dx))
							ytest = int(round(ytemp + shift*dy))

							# create IDs for the daughter cells
							currD1Id = totalIdCount + 1
							currD2Id = totalIdCount + 2

							# store daughters' info in the phylogeny matrix
							phylogeny.loc[phylogeny.id == fatherID, 'd1'] = currD1Id
							phylogeny.loc[phylogeny.id == fatherID, 'd2'] = currD2Id
							phylogeny.loc[phylogeny.id == fatherID, 'type'] = cellType
							phylogeny = phylogeny.append(pd.DataFrame([[currD1Id, fatherID, 0, 0, 0, cellType], [currD2Id, fatherID, 0, 0, 0, cellType]], columns = ['id', 'f', 'd1','d2','dead', 'type']))
							
							# store both daughters' IDs and their states in the corresponding 3D grids
					            # first daughter is registered over the father's space
							cellIDs[xtemp][ytemp] = currD1Id
							cellDeadLive[xtemp][ytemp] = 1
							liveCellTypes[xtemp][ytemp] = 1 if cellType == 1 else 2

							# second daughter is registered on the 'first' shifted cell's place
							cellIDs[int(round(xtemp + dx))][int(round(ytemp + dy))] = currD2Id
							cellDeadLive[int(round(xtemp + dx))][int(round(ytemp + dy))] = 1
							liveCellTypes[int(round(xtemp + dx))][int(round(ytemp + dy))] = 1 if cellType == 1 else 2
							
							# increase total ID count by two
							totalIdCount += 2

							# generate and store random number of mutations for the first daughter
							mutNum = np.random.poisson(mu, 1)[0]
							cellMutations.append([currD1Id])
							cellMutations[currD1Id-1].append([totalMutCount + k for k in range(mutNum)])
							totalMutCount += mutNum

							# check if the second daughter cell dies with the probability = beta (local death)
							deathTest = np.random.uniform(0,1,1)[0]
							if deathTest < (1-beta):
								# dies: erase the second daughter from all 2D grids and phylogeny matrix
								cellIDs[int(round(xtemp + dx))][int(round(ytemp + dy))] = 0
								cellDeadLive[int(round(xtemp + dx))][int(round(ytemp + dy))] = 0
								liveCellTypes[int(round(xtemp + dx))][int(round(ytemp + dy))] = 0
								phylogeny.loc[phylogeny.id == currD2Id, 'dead'] = -1
								cellMutations.append([currD2Id])
								cellMutations[currD2Id-1] = None
							else:
								# survives: generate and store random number of mutations for the second daughter
								mutNum = np.random.poisson(mu, 1)[0]
								cellMutations.append([currD2Id])
								cellMutations[currD2Id-1].append([totalMutCount + k for k in range(mutNum)])
								totalMutCount += mutNum
			
		# count current live cells
		curLiveCells = np.size(cellDeadLive[np.nonzero(cellDeadLive)])

		totalLiveCells.append(curLiveCells)

		# check if the currently uccupied space is on the grid boundary and if so, stop the growth
		if xtest == grid_size or ytest == grid_size or xtest == 1 or ytest == 1:
			finish = 1 # finishing the growth
		

	# store current generation count and the total number of mutations within it
	globalMutations = np.vstack((globalMutations, [generation, totalMutCount]))

	globalMutations = np.delete(globalMutations, 0, 0)

	globalMutations = pd.DataFrame(globalMutations, columns = ['gen', 'mutCount']) 

	totalLiveCells = pd.DataFrame(totalLiveCells, columns = ['cell_count'])

	timeSteps = pd.DataFrame(time)

	timeSteps.to_csv('timeSteps2D_{}.txt'.format(sd), sep = ' ', index = False)

	cellTypes = pd.DataFrame(liveCellTypes)

	with open('cellMutations{}.txt'.format(sd), 'w') as cellMutFile:
		for i in cellMutations:
			cellMutFile.write("%s\n" % i)

	with open('cellMutations{}_pickle.txt'.format(sd), 'wb') as mfile:
		pickle.dump(cellMutations, mfile)

	with open('phylogeny{}_pickle.txt'.format(sd), 'wb') as pfile:
		pickle.dump(phylogeny, pfile)

	with open('cellIDs{}_pickle.txt'.format(sd), 'wb') as cfile:
		pickle.dump(cellIDs, cfile)

	with open('cellDeadLive{}_pickle.txt'.format(sd), 'wb') as cfile:
		pickle.dump(cellDeadLive, cfile)

	with open('cellTypes{}_pickle.txt'.format(sd), 'wb') as tfile:
		pickle.dump(cellTypes, tfile)

	phylogeny.to_csv('phylogeny{}.txt'.format(sd), sep = ' ', index = False)

	globalMutations.to_csv('globalMutations{}.txt'.format(sd), sep = ' ', index = False)

	totalLiveCells.to_csv('totalLiveCells2D_{}.txt'.format(sd), sep = ' ', index = False)

	cellTypes.to_csv('cellTypes{}.txt'.format(sd), sep = ' ', index = False)


#simulateTumorGrowth2D(sd = int(sys.argv[1]), grid_size = int(sys.argv[2]), mu = int(sys.argv[3]), wildBirthRate = int(sys.argv[4]), mutantBirthRate = int(sys.argv[5]), alpha = float(sys.argv[6]), beta = float(sys.argv[7]), cellAggression = float(sys.argv[8]), cloneStartTime = float(sys.argv[9]))

#sd = 1; grid_size = 10; mu = 5;
#wildBirthRate = 1; mutantBirthRate = 2; alpha = 0; beta = 1;
#cellAggression = 1; cloneStartTime = 2;



simulateTumorGrowth2D(sd = 1, grid_size = 50, mu = 10, wildBirthRate = 1, mutantBirthRate = 2, alpha = 0, beta = 1, cellAggression = 1, cloneStartTime = 5)

