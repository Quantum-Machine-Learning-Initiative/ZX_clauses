import numpy as np

def random_clause(num_variables, k):
    '''Generate a random clause for k-SAT. Returns a list of three
    integers
    '''
    pass

def random_SAT_instance(num_variables, k, num_clauses):
    '''
    Generate a random instance of K-SAT
    '''

    # you can also change [] to () and make a generator expression but
    # that's for later
    return [random_clause(num_variables, k) for i in num_clauses]

def make_SAT_Hamiltonian(sat_instance, pauli):
    '''Creates a Hamiltonian from a SAT instance. pauli specifies which
    basis to use
    returns a sparse matrix
    '''
    pass

def find_energy_gap(H):
    '''finds the energy gap of a given Hamiltonian, i.e. the difference
    between the ground state and the first excited state
    '''
    pass

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
