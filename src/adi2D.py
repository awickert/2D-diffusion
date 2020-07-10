# Python prototype to test method
import numpy as np
from scipy.sparse import csc_matrix, spdiags
from scipy.sparse.linalg import spsolve
from matplotlib import pyplot as plt

tmax = 100. # seconds

D = 1E0
dx = 1.
dt = 1.
K = D * dx / dt**2

ncols = 51
nrows = 51
n = nrows*ncols
T = np.zeros((nrows, ncols))
T[ncols//2, nrows//2] = 1.

# A from Ax=b: Columns
A_upper_right_columns = -K * np.ones(nrows)
A_main_diagonal_columns = 2*K * np.ones(nrows) + 1
A_lower_left_columns = -K * np.ones(nrows)

diagonals = np.vstack(( A_lower_left_columns,
                        A_main_diagonal_columns,
                        A_upper_right_columns ))
offsets = np.array([-1, 0, 1])
A_LHSmatrix_columns = spdiags(diagonals, offsets, nrows, nrows, format='csr')

# A from Ax=b: Rows
A_upper_right_rows = -K * np.ones(ncols)
A_main_diagonal_rows = 2*K * np.ones(ncols) + 1
A_lower_left_rows = -K * np.ones(ncols)

diagonals = np.vstack(( A_lower_left_rows,
                        A_main_diagonal_rows,
                        A_upper_right_rows ))
offsets = np.array([-1, 0, 1])
A_LHSmatrix_rows = spdiags(diagonals, offsets, ncols, ncols, format='csr')


plt.figure()

t = 0
while t < tmax:

    # Loop over columns:
    for j in range(ncols):
        B_RHSmatrix = T[:,j]
        T[:,j] = spsolve(A_LHSmatrix_rows, B_RHSmatrix)
        
    # Loop over rows
    for i in range(nrows):
        B_RHSmatrix = T[i,:]
        T[i,:] = spsolve(A_LHSmatrix_rows, B_RHSmatrix)
    
    t += dt
    
    if t % 2 == 0:
        plt.cla()
        plt.imshow(T)
        plt.pause(0.01)

