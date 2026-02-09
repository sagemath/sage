from sage.all import *
import time

def profile_me():
    # get_memory_usage() is a built-in Sage function
    # It returns a float representing memory in MB
    try:
        start_mem = get_memory_usage()
        print(f"Starting Memory: {start_mem:.2f} MB")
    except NameError:
        print("Sage memory tools not found, skipping memory check.")
        start_mem = 0

    # The Workload: Large Polynomial factorization
    R, x = ZZ['x'].objgen()
    p1 = (x**500 + 3*x**250 + 1)**10
    p2 = (x**300 + 2*x**150 + 1)**12
    
    print("Performing heavy factorization...")
    start_time = time.time()
    
    # This triggers the Cython/C code
    product = p1 * p2
    factors = product.factor()
    
    end_time = time.time()
    
    if start_mem > 0:
        end_mem = get_memory_usage()
        print("-" * 30)
        print(f"Time Taken: {end_time - start_time:.4f} seconds")
        print(f"Memory Increment: {end_mem - start_mem:.2f} MB")
        print("-" * 30)

if __name__ == "__main__":
    profile_me()