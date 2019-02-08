import numpy as np
import scipy.sparse
import scipy.sparse.linalg as SPLA
from functools import reduce
import matplotlib.pyplot as plt

pauli_Z = np.array([[1, 0],[0, -1]])
pauli_X = np.array([[0, 1],[1, 0]])
I = np.eye(2)

def random_clause(num_variables, k=3):
    '''Generate a random clause for k-SAT. Returns a list of three
    integers
    '''
    clause = np.random.choice(num_variables, k, replace=False) + 1
    clause = [c * ((-1)**np.random.randint(2)) for c in clause]
    return clause

def random_SAT_instance(num_variables, num_clauses, k=3):
    '''
    Generate a random instance of K-SAT
    '''
    return [random_clause(num_variables, k) for i in range(num_clauses)]

def make_SAT_Hamiltonian(sat_instance, num_bits, pauli=pauli_Z):
    '''Creates a Hamiltonian from a SAT instance. pauli specifies which
    basis to use
    returns a sparse matrix
    '''    
#    num_bits = np.max(np.abs(sat_instance))
    H = scipy.sparse.csr_matrix((2**num_bits, 2**num_bits))
    for clause in sat_instance:
        term_list = [I] * num_bits
        for i in clause:
            term_list[abs(i) - 1] = 0.5 * (I + pauli * np.sign(i))
        H += reduce(lambda A, B: scipy.sparse.kron(A, B, format="csr"), term_list)
    return H
        
    

def find_energy_gap(H, num_vals=6):
    '''finds the energy gap of a given Hamiltonian, i.e. the difference
    between the ground state and the first excited state.  Eigenvalue
    finder can't avoid degenerate eigenvalues (I think), so it will
    have to look for all lowest eigs until

    '''
    w, v = SPLA.eigsh(H, k=num_vals, which="SA")
    w_sorted = sorted(w)
    #print(w_sorted)
    E_0 = w_sorted[0]
    pos = np.where(np.array(w_sorted) > (E_0 + 1e-9))
    if len(pos[0]):
        return (w_sorted[pos[0][0]] - E_0)
    else:
        raise ValueError("GS too degenerate or dE < tol!")

def dense_energy_gap(H):
    '''Finds the energy gap for a small Hamiltonian written as a dense matrix'''
    w, v = np.linalg.eigh(H)
    w_sorted = np.sort(w)
    E_0 = w_sorted[0]
    pos = np.where(w_sorted > (E_0 + 1e-9))
    if len(pos[0]):
        return (w_sorted[pos[0][0]] - E_0)
    else:
        # H is fully degenerate otherwise
        return 0

def make_ZX_factory(num_variables, num_clauses_1, num_clauses_2, k_1=3, k_2=3):
    '''Creates a function that takes no arguments and spits out
    Hamiltonians with 
    '''
    def factory():
        instance_1 = random_SAT_instance(num_variables, num_clauses_1, k=k_1)
        instance_2 = random_SAT_instance(num_variables, num_clauses_2, k=k_2)
        if num_clauses_2 > 0 and num_clauses_1 > 0:
            H = make_SAT_Hamiltonian(instance_1, num_variables, pauli_Z) + make_SAT_Hamiltonian(instance_2, num_variables, pauli_X)
        elif num_clauses_2 > 0:
            H = make_SAT_Hamiltonian(instance_2, num_variables, pauli_X)
        elif num_clauses_1 > 0:
            H = make_SAT_Hamiltonian(instance_1, num_variables, pauli_Z)
        else:
            raise ValueError("Empty Hamiltonian")
        return H
    return factory

def symmetric_ZX_factory(num_variables, num_clauses_1, num_clauses_2, k=3):
    '''Creates a Hamiltonian factory that produces ``symmetric'' SAT
    Hamiltonians: if the number of Z and X clauses are equal, they
    coincide, otherwise the first min(num_clauses_1, num_clauses_2)
    coincide.
    '''
    def factory():
        big_instance = random_SAT_instance(num_variables, np.max([num_clauses_1, num_clauses_2]), k=k)
        instance_1 = big_instance[:num_clauses_1]
        instance_2 = big_instance[:num_clauses_2]
        if num_clauses_2 > 0 and num_clauses_1 > 0:
            H = make_SAT_Hamiltonian(instance_1, num_variables, pauli_Z) + make_SAT_Hamiltonian(instance_2, num_variables, pauli_X)
        elif num_clauses_2 > 0:
            H = make_SAT_Hamiltonian(instance_2, num_variables, pauli_X)
        elif num_clauses_1 > 0:
            H = make_SAT_Hamiltonian(instance_1, num_variables, pauli_Z)
        else:
            raise ValueError("Empty Hamiltonian")
        return H
    return factory
    
def average_energy_gap(ham_factory, num_instances=100, num_eigs=10):
    '''
    Calculate the average energy gap for the ZX clause Hamiltonians
    '''
    dE_data = []
    failures = 0
    for i in range(num_instances):        
        H = ham_factory()
#        print(H)
        try:
            dE = find_energy_gap(H, num_eigs)
            #dE = 10
            dE_data.append(dE)
        except ValueError:
            failures += 1

    print('{} failures'.format(failures))

    if dE_data:
        return np.mean(dE_data)
    else:
        raise ValueError("No convergencies")

#        big_instance = random_SAT_instance(num_variables, np.max(num_clauses_1, num_clauses_2), k = k_1)
#        instance_1 = big_instance[:num_clauses_1]
#        instance_2 = big_instance[:num_clauses_2]
        
##        instance_1 = random_SAT_instance(num_variables, num_clauses_1, k=k_1)
##        instance_2 = random_SAT_instance(num_variables, num_clauses_2, k=k_2)
#        if num_clauses_2 > 0 and num_clauses_1 > 0:
#            H = make_SAT_Hamiltonian(instance_1, num_variables, pauli_Z) + make_SAT_Hamiltonian(instance_2, num_variables, pauli_X)
#        elif num_clauses_2 > 0:
#            H = make_SAT_Hamiltonian(instance_2, num_variables, pauli_X)
#        elif num_clauses_1 > 0:
#            H = make_SAT_Hamiltonian(instance_1, num_variables, pauli_Z)
#        else:
#            raise ValueError("Empty Hamiltonian")
#        #print(np.shape(H))
#        #print(H.nnz)
#        try:
#            N = np.shape(H)[0]
#            if N <= num_eigs:
#               # print(N)
#                dE = dense_energy_gap(H.todense())
#            else:
#                dE = find_energy_gap(H, num_eigs)                
#            #n_eigs = np.min(np.shape(H)[0] - 1, num_eigs) 
#            dE_data.append(dE)
#
#            
#        except ValueError:
#            failures += 1
#
#    print('Failed to find gap {} times'.format(failures))
#    if dE_data:
#        return np.mean(dE_data)
#    else:
#        raise ValueError("No convergencies")
#
    
if (__name__=='__main__'):
    pass
    #num_variables = 6
    #num_instances = 25
    #nz = 10
    #nx = 10
    #zs = [4 * i for i in range(nz)]
    #xs = [4 * i for i in range(nx)]
    #data = np.zeros((nz, nx))
    #for i, z in enumerate(zs):
    #    for j, x in enumerate(xs):
    #        print(z, x, end=' ')
    #        try:
    #            data[i, j] = average_energy_gap(num_variables, z, 3, x, 3,
    #                                            num_instances=num_instances,
    #                                            num_eigs=50)
    #        except ValueError as e:
    #            print(e)
    #            data[i, j] = np.nan
    #print(data)
    f = symmetric_ZX_factory(6, 50, 10)
    dE = average_energy_gap(f, num_eigs=20)
    print(dE)
    #print(f())
    #quit()
