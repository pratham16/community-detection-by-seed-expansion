import networkx as nx
import sys
import argparse
from collections import defaultdict
import time
import codecs
from multiprocessing import Pool
import pprgrow_min_cond


INF = float('inf')

def pprgrow(args):
	seed,G,stopping,nruns,alpha,maxexpand,fast = args
	expandseq = [2,3,4,5,10,15]
	expands = list()
	curmod = 1
	while len(expands) < nruns:
		temp = [curmod*i for i in expandseq]
		for i in temp:
			expands.append(i)
		curmod *= 10

	expands = expands[:nruns]
	maxdeg = max(G.degree(G.nodes()).values())
	bestcond = INF
	bestset = list()
	bestexpand = 0.0
	bestceil = 0.0
	if fast==True:
		expands = [1000]
	print maxdeg
	for ei in range(len(expands)):
		if fast==True:
			curexpand = expands[ei]
		else:
			curexpand = expands[ei]*len(seed)+maxdeg
		assert len(seed)>0.0
		if curexpand > maxexpand:
			continue
		if stopping=='cond':
			start_time = time.time()
			curset, cond = pprgrow_min_cond.pprgrow(G,seed,alpha,curexpand)
			end_time = time.time()
			print end_time - start_time
			if cond < bestcond:
				bestcond = cond
				bestset = curset
				bestexpand = curexpand

	return curset


def growclusters(G,seeds,expansion,stopping,nworkers,nruns,alpha,maxexpand,fast):
	if maxexpand == INF:
		maxexpand = G.number_of_edges()

	n = G.number_of_nodes()
	ns = len(seeds)
	communities = list()

	if nworkers==1:
		for i in range(ns):
			seed = seeds[i]
			if expansion=='ppr':
				curset = pprgrow((seed,G,stopping,nruns,alpha,maxexpand,fast))
			else:
				print 'Method not implemented yet'

			#H[curset,i] = 1.0
			communities.append(curset)
			print 'Seed',i,'Done'

	else:
		print 'Initiating parallel seed expansion'
		slen = len(seeds)
		args = zip(seeds,[G]*slen,[stopping]*slen,[nruns]*slen,[alpha]*slen,\
			[maxexpand]*slen,[fast]*slen)
		p = Pool(nworkers)
		if expansion=='ppr':
			communities = p.map(pprgrow, args)

	return communities


# TO DO: Write function to remove duplicate communities
def remove_duplicates(G, communities,delta):
	# Create node2com dictionary
	node2com = defaultdict(list)
	com_id = 0
	for comm in communities:
		for node in comm:
			node2com[node].append(com_id)
		com_id += 1

	deleted = dict()
	i = 0
	for i in range(len(communities)):
		comm = communities[i]
		if deleted.get(i,0) == 0:
			nbrnodes = nx.node_boundary(G, comm)
			for nbr in nbrnodes:
				nbrcomids = node2com[nbr]
				for nbrcomid in nbrcomids:
					if i!=nbrcomid and deleted.get(i,0)==0 and deleted.get(nbrcomid,0)==0:
						nbrcom = communities[nbrcomid]
						distance = 1.0 - (len(set(comm) & set(nbrcom))*1.0 / (min(len(comm),len(nbrcom))*1.0))

						if distance <= delta:
							# Near duplicate communities found.
							# Discard current community
							# Followed the idea of Lee et al. in GCE
							deleted[i] = 1
							for node in comm:
								node2com[node].remove(i)
	for i in range(len(communities)):
		if deleted.get(i,0)==1:
			communities[i] = []

	communities = filter(lambda c: c!=[], communities) # Discard empty communities
	return communities

def neighbor_inflation(G,seeds):
	# Seed = union(seedNode, egonet(seedNode))
	for i in range(len(seeds)):
		seed = seeds[i]
		egonet = list()
		for s in seed:
			egonet.append(s)
			[egonet.append(k) for k in G.neighbors(s)]
		seeds[i] = list(set(egonet))
		# print sorted([int(k) for k in seeds[i]])

	return seeds


def __main():
	# graph_file = sys.argv[1]
	# seed_file = sys.argv[2]

	parser = argparse.ArgumentParser()
	parser.add_argument('graph_file',type=str,help='Input Graph File Path')
	parser.add_argument('seed_file',type=str,help='Input Seeds File Path')
	parser.add_argument('--ninf',help='Neighbourhood Inflation parameter',action='store_false')
	parser.add_argument('--expansion',type=str,help='Seed expansion: PPR or VPPR',default='ppr')
	parser.add_argument('--stopping',type=str,help='Stopping criteria',default='cond')
	parser.add_argument('--nworkers',type=int,help='Number of Workers',default=1)
	parser.add_argument('--nruns',type=int,help='Maximum number of runs',default=13)
	parser.add_argument('--alpha',type=float,help='alpha value for Personalized PageRank expansion',default=0.99)
	parser.add_argument('--maxexpand',type=float,help='Maximum expansion allowed for approximate PPR',default=INF)
	parser.add_argument('--delta',type=float,help='Minimum distance parameter for near duplicate communities',default=0.2)
	args = parser.parse_args()

	G = nx.read_edgelist(args.graph_file)
	print "Graph Loaded"
	with codecs.open(args.seed_file,'r',encoding='utf-8') as f:
		seeds = f.readlines()
		seeds = [x.strip() for x in seeds]
		seeds = [x.split() for x in seeds]
	print "Seeds Loaded\n"

	# print seeds
	if args.ninf==True:
		seeds = neighbor_inflation(G,seeds)
	print "Initiating Seed Expansion------"
	communities = growclusters(G,seeds,args.expansion,args.stopping,args.nworkers,args.nruns,args.alpha,args.maxexpand,False)
	print "Seed Expansion Finished.\n"
	print "Initiating removal of near duplicate communities."
	communities = remove_duplicates(G,communities,args.delta)
	print "Duplicate communities removed\n"

	
	print "Writing communities to output file:"
	for c in communities:
		print c

if __name__ == "__main__":
	__main()
