import numpy as np

def Reduction(matrix):
    for i in range(matrix.shape[0]):
        min_row = np.min(matrix[i, :])
        matrix[i, :] -= min_row

    for j in range(matrix.shape[1]):
        min_col = np.min(matrix[:, j])
        matrix[:, j] -= min_col

    return matrix

def adjustMatrix(matrix, covered_rows, covered_cols):
    dimensions = matrix.shape
    n = dimensions[0]
    smallest = np.inf

    print("smallest uncovered number should be 0.31")

    for i in range(n):
        for j in range(n):
            if not covered_rows[i] and not covered_cols[j]:
                if matrix[i, j] < smallest:
                    smallest = matrix[i, j]
    print("We get: ", smallest)
    #smallest element is correct...Good!!!!
    
    for x in range(n):
        if not covered_rows[x]:
            matrix[x, :] -= smallest

    for y in range(n):
        if covered_cols[y]:
            matrix[:, y] += smallest            

    return matrix



def coverZeros(matrix):
    r,c  = matrix.shape
    covered_rows = [False]*r
    covered_cols = [False]*c

    # zero positions in an array
    zeros= [(i,j) for i in range(r) for j in range(c) if matrix[i][j] ==0]

    while len(zeros) > 0:
        # uncovered zeros in each row
        row_zeros = [sum(1 for x,y in zeros if x ==i and not covered_cols[y]) for i in range(r)]
        #uncovered zeros in each colum
        col_zeros = [sum(1 for x,y in zeros if y ==j and not covered_rows[x]) for j in range(c)]
        
        max_row = max(row_zeros)
        max_col = max(col_zeros)

        if(max_row >= max_col):
            row_idx = row_zeros.index(max_row)
            covered_rows[row_idx] = True
        else:
            col_idx = col_zeros.index(max_col)
            covered_cols[col_idx] =True
        zeros = [(x,y) for x, y in zeros if not (covered_rows[x] or covered_cols[y])]

    lines_used = sum(covered_rows) + sum(covered_cols)
    return covered_rows, covered_cols, lines_used

def finalAssignments(matrix):
    dup_matrix = matrix.copy()
    n = dup_matrix.shape[0]

    zeros_in_rows = np.sum(dup_matrix == 0, axis=1)

    assignment = np.full((n,), -1)

    while np.any(zeros_in_rows > 0):
        min_zeros = np.argmin(zeros_in_rows + (zeros_in_rows == 0)*n)#we get the index
        for j in range(n): #found it, now find the first 0
            if dup_matrix[min_zeros, j] == 0:
                
                assignment[min_zeros] = j
                
                dup_matrix[min_zeros, :] = np.inf 
                dup_matrix[:, j] = np.inf
                break
            # continue
        zeros_in_rows = np.sum(dup_matrix == 0, axis=1)

    return assignment
    
def hungarianMethod(matrix):
    matrix  = Reduction(matrix)
    print("after reduction")
    print(matrix)

    covered_rows, covered_cols, lines_used = coverZeros(matrix)
    # print("lines used: ", lines_used)
    # print("rows:", covered_rows)
    # print("cols:", covered_cols)
    #works till here...good!

    n= matrix.shape[0]

    while lines_used < n:
        matrix = adjustMatrix(matrix, covered_rows, covered_cols)
        print("not enough lines used... adjusted matrix")
        print(matrix)

        covered_rows, covered_cols , lines_used = coverZeros(matrix)
        print("lines used: ", lines_used)
        print("rows:", covered_rows)
        print("cols:", covered_cols)
    
    print("optimal matrix: ")
    print(matrix)

    return finalAssignments(matrix)




def constructMatrix(teammate_positions, formation_positions):
    num_players = len(teammate_positions)
    cost_matrix = np.zeros((num_players, num_players))

    for i in range(num_players):
        teammate_pos = np.array(teammate_positions[i])
        for j in range(num_players):
            formation_pos = np.array(formation_positions[j])
            cost_matrix[i, j] = np.linalg.norm(teammate_pos - formation_pos)
    cost_matrix = np.round(cost_matrix, 2)
    return cost_matrix

def role_assignment(teammate_positions, formation_positions): 
    
    cost_matrix = constructMatrix(teammate_positions, formation_positions)
    print("Cost matrix before reduction:\n", cost_matrix)

    assignment = hungarianMethod(cost_matrix)
    print("assignemnet:\n", assignment)

    point_preferences = {}
    for row, col in enumerate(assignment):
        x, y = formation_positions[col]
        point_preferences[row + 1] = np.array([x, y])

    print("Point preferences:", point_preferences)
    return point_preferences

def pass_receiver_selector(player_unum, teammate_positions, final_target):
    pass_receiver_unum = player_unum + 1

    if pass_receiver_unum != 12:
        target = teammate_positions[pass_receiver_unum - 1]
    else:
        target = final_target

    return target

# Sample data
teammate_positions = [
    [-14., 0],
    [-9.00671841, -5.0212432],
    [-9, 0],
    [-8.9, 5.0],
    [-5.05693005, -5.02106853],
    [-5.0216475, -0.03376451],
    [-4.9, 4.9],
    [-1.06363536, -5.95482965],
    [-1.03039866, -2.47499713],
    [-1.01752962, 2.52498091],
    [-1.02189553, 6.02497721]
]

formation_positions = [
    [-13, 0],
    [-10, -2],
    [-11, 3],
    [-8, 0],
    [-3, 0],
    [0, 1],
    [2, 0],
    [3, 3],
    [9, 1],
    [12, 0],
    [8, 0]
]

# Execute role assignment
assigned_roles = role_assignment(teammate_positions, formation_positions)
