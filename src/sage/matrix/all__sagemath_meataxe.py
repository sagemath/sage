def create_matrix(n, **kwargs):
    if 'base_ring' in kwargs:
        return matrix(n, base_ring=kwargs['base_ring'])
    else:
        return matrix(n, **kwargs)

# Example usage:
M = create_matrix(5, base_ring=QQ)
print(M)

# Alternatively,
M = create_matrix(5, ring=QQ)
print(M)
