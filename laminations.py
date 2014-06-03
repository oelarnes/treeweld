import random
import sys

sys.setrecursionlimit(2001)

def make_lam_rich(s, j):
	done = False
	i = j
	lam = []
	lam_deg = {}
	node_lists = []
	if j == 0:
		deg = 0
		node_list = []
	else:
		deg = 1
		node_list = [i]
	while not done:
		if i not in s:
			done = True
		else:
			temp = make_lam_rich(s, i+1)
			lam = lam + [(i+1, temp[2]+1)] + temp[0]
			lam_deg = dict(lam_deg.items() + temp[1].items())
			i = temp[2] + 1
			node_list.append(i)
			deg += 1
			node_lists.extend(temp[3])
	for node in node_list:
		lam_deg[node] = deg
	node_lists.append(node_list)
	return lam, lam_deg, i, node_lists
	
def make_lam(s):
	lam = make_lam_rich(s,0)
	return lam[0], lam[1].values(), lam[3]
	
def convert_to_set(lam):
	s = []
	for p in lam:
		if p[0] < p[1]:
			s.append(p[0]-1)
		else:
			s.append(p[1]-1)
	return s
	
def fix(L):
    n = 2*len(L)
    if n==0:
        return L
    p0 = first_pair(L)
    temp = [p for p in L if p[0] not in p0] 
    rem = [j+1 for j in range(n) if j+1 not in p0]
    inv_rem = {rem[j]:j+1 for j in range(n-2)}
    temp = fix([(inv_rem[p[0]], inv_rem[p[1]]) for p in temp])
    temp = [(rem[p[0]-1],rem[p[1]-1]) for p in temp]
    return [p0] + temp

def first_pair(L):
    n = 2*len(L)
    for p in L:
        if (p[1]-p[0])%n == 1:
            return (p[1],p[0])
        if (p[0]-p[1])%n == 1:
            return (p[0],p[1])

def expand(L, sub_div):
	temp = fix(L)
	new = []
	for pair in temp:
		new.extend([((pair[0]-1)*sub_div + j, pair[1]*sub_div - j + 1)\
			for j in range(1,sub_div+1)])
	return new
	
def random_lam(n, sub_div=1):
	b = random.sample(range(2*n), n)
	# find the root of the bridge walk
	val = 0
	min_val = 0
	min_ind = 0
	for i in range(2 * n):
		if i in b:
			val += 1
		else:
			val -= 1
			if val < min_val:
				min_ind = i + 1
				min_val = val
	b = [(x - min_ind) % (2 * n) for x in b]
 	lam = make_lam(b)
 	return lam
	
def all_lams(n):
	subs = all_subsets(n)
	lams = []
	for s in subs:
		lams.append(make_lam(s))
	return lams
	
def all_subsets(n):
	if n == 0:
		return [[]]
	subs = []
	for j in range(n):
		subs1 = [[el + 1 for el in s] for s in all_subsets(j)]
		subs2 = [[el + 2 * (j + 1) for el in s] for s in all_subsets(n - j - 1)]
		for s1 in subs1:
			subs += [[0] + s1 + s2 for s2 in subs2]
	return subs
	
inp = raw_input("Enter a number of edges to generate a random tree. Hit enter to generate a tree from a file edges.dat. Max 2000 vertices: ")

if inp != "":
	n = int(inp)
	k = raw_input("Enter the number of trees or type all: ")
	if k == "all":
		lams = all_lams(n)
		k = len(lams)
	else:
		k = int(k)
		lams = []
		for j in range(k):
			lams.append(random_lam(n))
else:
	lam = []
	f = open("edges.dat", "r")
	lam = f.read().split('\n')[:-1]
	lam = [(int(p.split(' ')[0]), int(p.split(' ')[1])) for p in lam]
	print lam
	print convert_to_set(lam)
	lams = [make_lam(convert_to_set(lam))]
	n = len(lam)
				
sub_div = int(raw_input("Enter the number of welding subdivisions: "))
f = open("edges_sub.dat", "w")
f.write(str(n) + ' ' + str(sub_div) + '\n')
for lam in lams:
	degs = lam[1]
	for j in range(n):
		f.write(str(degs[2 * j]) + ' ' + str(degs[2 * j + 1]) + '\n')
	vert = lam[2]
	for j in range(n + 1):
		f.write(str(vert[j][0]) + ' ' + str(vert[j][-1]) + '\n')
	lam = expand(lam[0], sub_div)
	for p in lam:
		f.write(str(p[0]) + ' ' + str(p[1]) + '\n')
		

	