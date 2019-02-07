import numpy as np
import scipy.sparse
import scipy.sparse.linalg as SPLA
from functools import reduce

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

def make_SAT_Hamiltonian(sat_instance, pauli=pauli_Z):
    '''Creates a Hamiltonian from a SAT instance. pauli specifies which
    basis to use
    returns a sparse matrix
    '''
    num_bits = np.max(np.abs(sat_instance))
    H = scipy.sparse.csr_matrix((2**num_bits, 2**num_bits))
    for clause in sat_instance:
        term_list = [I] * num_bits
        for i in clause:
            term_list[abs(i) - 1] = 0.5 * (I + pauli * np.sign(i))
        H += reduce(lambda A, B: scipy.sparse.kron(A, B, format="csr"), term_list)
    return H
        
    

def find_energy_gap(H, num_vals=6):
    '''finds the energy gap of a given Hamiltonian, i.e. the difference
    between the ground state and the first excited state.
    Eigenvalue finder can't avoid degenerate eigenvalues (I think), so it will have to look for all lowest eigs until 
    '''
    w, v = SPLA.eigsh(H, k=num_vals, which="SA")
    w_sorted = sorted(w)
    print(w_sorted)
    E_0 = w_sorted[0]
    pos = np.where(np.array(w_sorted) > (E_0 + 1e-9))
    if len(pos[0]):
        return (w_sorted[pos[0][0]] - E_0)
    else:
        raise ValueError("GS too degenerate or dE < tol!")

def average_energy_gap(num_variables_1, num_clauses_1, k_1,
                       num_variables_2, num_clauses_2, k_2, num_instances=50):
    '''
    Calculate the average energy gap for the ZX clause Hamiltonians
    '''
    dE_data = []
    for i in num_instances:
        instance_1 = random_SAT_instance(num_variables_1, k_1, num_clauses_1)
        instance_2 = random_SAT_instance(num_variables_2, k_2, num_clauses_2)
        H = make_SAT_Hamiltonian(instance_1, 'Z') + make_SAT_Hamiltonian(instance_2, 'X')
        dE = find_energy_gap(H)
        dE_data.append(dE)
    return np.mean(dE_data)

    
if (__name__=='__main__'):
    pass
