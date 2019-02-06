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



if (__name__=='__main__'):
    pass
