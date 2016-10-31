#neighbor join algorithm lab

def readFromFile(filename):
    """
    parse an n-by-n matrix of distances and a header list.
    Input: name of file to parse
    Output: distance matrix D; headers list
    """
    D = []
    headers = []
    with open(filename, 'r') as fin:
        for line in fin:
            row = line.strip().split()
            if len(headers) == 0:
                headers = row
            else:
                row = [int(a) for a in row[1:]]
                D.append(row)
    return D,headers


def make_dict(D, headers):
    d = {}
    for i in range(len(headers)):
        sub_d = {}
        for j in range(len(headers)):
            sub_d[headers[j]] = D[i][j]
        d[headers[i]] = sub_d
    return d


def min_of_layered_dict(d):
    """
    turns the distance matrix into a layered dictionary
    """
    m = float('inf')
    mk1 = None
    mk2 = None
    for k in d:
        sub_d = d[k]
        mk_sub_d = min(sub_d,key=sub_d.get)
        if sub_d[mk_sub_d] < m:
            m = sub_d[mk_sub_d]
            mk1 = k
            mk2 = mk_sub_d

    return mk1, mk2

    
def neighbor_join_v2(D):
    """
    Neighbor join using dictionary representation instead of matrix
    """
    num_otus = len(D)
    interior_points = 0
    T = Tree()
    N = [] #curated list of nodes that are "in" D so that I don't have to prune and update the dictionary when adding/removing rows from the Distance matrix
    for k in D:
        N.append(k)

    while interior_points < num_otus - 2:
        n = len(D)
        Q = {}
        for i in N:
            sub_q = {}
            for j in N:
                sub_q[j] = get_q_v2(i,j,D)
            Q[i] = sub_q
                
        merge1,merge2 = min_of_layered_dict(Q)

        bi = get_new_dist_v2(merge1, merge2, D)
        bj = D[merge1][merge2] - bi

        new_header = "(" + str(merge1) + ", " + str(merge2) + ")"

        
        T.add(merge1,new_header,bi,interior_points)
        T.add(merge2,new_header,bj,interior_points)

        N.remove(merge1)
        N.remove(merge2)

        D[new_header] = {}

        
        for i in N:
            new_d = .5 * (D[merge1][i] + D[merge2][i] - D[merge1][merge2])
            D[new_header][i] = new_d
            D[i][new_header] = new_d
        D[new_header][new_header] = 0

        N.append(new_header)

        interior_points += 1

    return T








def neighbor_join(D):
    num_otus = len(D)
    interior_points = 0
    #initialize a tree here
    T = Tree()
    
    while interior_points < num_otus - 2:
        n = len(D)
        sub_q = [0 for i in range(n)]
        Q = [sub_q[:] for i in range(n)]
        print("Length of D is: ",len(D))

        for i in range(len(D)):
            print(D[i])

        for i in range(len(D)):
            for j in range(len(D[i])):
                Q[i][j] = get_q_val(i,j,D)

        mv = float('inf')
        mi = 0
        mj = 0

        for i in range(len(Q)):
            for j in range(len(Q[i])):
                if Q[i][j] < mv:
                    mi = i
                    mj = j
                    mv = Q[i][j]

        print("Length of Q is:",len(Q))
        for i in range(len(Q)):
            print(Q[i])
            
        bi = get_new_dist(mi, mj, D)
        bj = D[mi][mj] - bi

        T.add(mi, ("(" + str(mi) + ", " + str(mj) + ")"), bi, interior_points)
        T.add(mj, ("(" + str(mi) + ", " + str(mj) + ")"), bj, interior_points)

        D = make_merged_D(mi,mj,D)
        interior_points += 1

    return T
        


def make_merged_D(i,j,D):
    i_row = D[i]
    j_row = D[j]
    dij = D[i][j]
    u_row = []

    new_D = []
    for k in range(len(D)):
        if k != i and k != j:
            dik = D[i][k]
            djk = D[j][k]
            duk = (.5)*(dik + djk - dij)
            row_copy = []
            for m in range(len(D[k])):
                if m != i and m != j:
                    row_copy.append(D[k][m])
            u_row.append(duk)
            row_copy.append(duk)
            new_D.append(row_copy)
    u_row.append(0)
    new_D.append(u_row)

    return new_D
        
        

def get_new_dist(i,j,D):
    n = len(D)
    i_sum = 0
    j_sum = 0

    for dist in D[i]:
        i_sum += dist

    for dist in D[j]:
        j_sum += dist

    b = (.5) * D[i][j] + (1 / float(2 * (n - 2))) * (i_sum - j_sum)

    return b
            
        
def get_q_val(i,j, D):
    n = len(D)
    i_sum = 0
    j_sum = 0
    for dist in D[i]:
        i_sum += dist

    for dist in D[j]:
        j_sum += dist

    if i != j:
        return ((n - 2) * D[i][j]) - i_sum - j_sum

    else:
        return float('inf')


class Tree:
    def __init__(self):
        self.edge_ls = []

    def add(self, s, t, weight,iteration):
        self.edge_ls.append(Edge(s,t,weight,iteration))

    def __str__(self):
        s = ""
        for e in self.edge_ls:
            s += str(e) + "\n"
        return s


        
        
class Edge(tuple):
    def __new__(cls, s, t, weight, iteration):
        return tuple.__new__(cls, (frozenset([s,t]), weight, iteration))

    def nodes(self):
        return self[0]

    def weight(self):
        return self[1]

    def iteration(self):
        return self[2]

    def __str__(self):
        return str(self.nodes()) + " weight:" + str(self.weight()) + " iteration:" + str(self.iteration()) 
    
    def __setattr__(self, *ignored):
        raise NotImplementedError

    def __delattr__(self, *ignored):
        raise NotImplementedError

    
        


#######################
#v2 helper functions
#######################

def get_q_v2(i,j, D):
    n = len(D)
    i_sum = 0
    j_sum = 0
    for dist in D[i].values():
        i_sum += dist

    for dist in D[j].values():
        j_sum += dist

    if i != j:
        return ((n - 2) * D[i][j]) - i_sum - j_sum

    else:
        return float('inf')



def get_new_dist_v2(i,j,D):
    n = len(D)
    i_sum = 0
    j_sum = 0

    for dist in D[i].values():
        i_sum += dist

    for dist in D[j].values():
        j_sum += dist

    b = (.5) * D[i][j] + (1 / float(2 * (n - 2))) * (i_sum - j_sum)

    return b





def main():
    D,header = readFromFile("class-example.txt")
    #T = neighbor_join(D)
    #print(T)
    big_dict = make_dict(D,header)

    T = neighbor_join_v2(big_dict)
    print(T)
    

if __name__ == "__main__":
    main()
