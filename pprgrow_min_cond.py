import Queue
from operator import itemgetter
import networkx as nx

__all__ = ['compute_local_pagerank', 'cluster_from_sweep', 'pprgrow']

INF = float('inf')
def compute_local_pagerank(G, r, p, alpha, epsilon, max_push_count, q):
	for node, res in r.iteritems():
		if res > epsilon*G.degree(node):
			q.put(node)

	push_count = 0
	while q.empty()==False and push_count < max_push_count:
		push_count += 1
		u = q.get()
		du = G.degree(u)*1.0
		moving_prob = r[u] - 0.5*epsilon*du
		r[u] = 0.5*epsilon*du
		p[u] = p.get(u,0.0) + (1-alpha)*moving_prob

		neighbor_update = alpha*moving_prob/du

		for nbr in G.neighbors(u):
			nbrdeg = G.degree(nbr)
			rxold = r.get(nbr,0.0)
			rxnew = rxold + neighbor_update
			r[nbr] = rxnew
		if rxnew > epsilon*nbrdeg and rxold <= epsilon*nbrdeg:
			q.put(nbr)

	return push_count

def cluster_from_sweep(G, p, community):
	p_sorted = sorted(p.iteritems(),key=itemgetter(1),reverse=True)

	conductance = list()
	volume = list()
	cutsize = list()

	i=0
	rank = dict()
	for node, ppr in p_sorted:
		rank[node] = i
		i += 1

	total_degree,curcutsize,curvolume,i = 2.0*G.number_of_edges(),0,0,0
	for node, ppr in p_sorted:
		deg = G.degree(node)*1.0
		change = deg*1.0
		for nbr in G.neighbors(node):
			if rank.get(nbr,-1) != -1 and rank.get(nbr,-1) < rank[node]:
				change -= 2

		curcutsize += change
		curvolume += deg

		cutsize.append(curcutsize)
		volume.append(curvolume)
		if curvolume==0 or total_degree-curvolume==0:
			conductance.append(1.0)
		else:
			conductance.append(curcutsize/min(curvolume,total_degree-curvolume))


	# print conductance
	lastind = len(conductance)
	mincond = INF
	mincondind = 0
	for i in range(lastind):
		if conductance[i] <= mincond:
			mincond = conductance[i]
			mincondind = i
	i=0
	for i in range(min(mincondind+1,len(p_sorted))):
		community.append(p_sorted[i][0])

	return mincond


def pprgrow(G, seed, alpha, targetvol):
	p = dict()
	r = dict()
	q = Queue.Queue()
	community = list()

	assert targetvol > 0.0
	assert alpha > 0.0 , alpha < 1.0

	for s in seed:
		r[s] = 1.0 / (len(seed)*1.0)

	pr_eps = 1.0 / max(10.0*targetvol,100.0)
	maxsteps = 1.0 / (pr_eps*(1.0-alpha))
	maxsteps = min(maxsteps, 0.5*(2.0**32-1.0))

	nsteps = compute_local_pagerank(G, r, p, alpha, pr_eps, int(maxsteps), q)
	if nsteps == 0:
		p = r

	# Scale the probabilities by their degree
	for node, pr in p.iteritems():
		pr *= 1.0/max(G.degree(node), 1.0)

	mincond = 2.0**32 - 1.0
	mincond = cluster_from_sweep(G, p, community)
	return community,mincond


#def __main():
	

#if __name__ == "__main__":
#	__main()
