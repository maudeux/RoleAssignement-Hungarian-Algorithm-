import numpy as np

def Reduction(matrix):
    # Row reduction
    for i in range(matrix.shape[0]):
        min_row = np.min(matrix[i, :])
        matrix[i, :] -= min_row

    # Column reduction
    for j in range(matrix.shape[1]):
        min_col = np.min(matrix[:, j])
        matrix[:, j] -= min_col

    return matrix

def adjustMatrix(matrix, covered_rows, covered_cols):
    n = matrix.shape[0]
    smallest = float('inf')

    # Find smallest uncovered value
    for i in range(n):
        for j in range(n):
            if not covered_rows[i] and not covered_cols[j]:
                if matrix[i, j] < smallest:
                    smallest = matrix[i, j]
    
    # Subtract from uncovered rows
    for i in range(n):
        if not covered_rows[i]:
            matrix[i, :] -= smallest
    
    # Add to covered columns
    for j in range(n):
        if covered_cols[j]:
            matrix[:, j] += smallest

    return matrix

def coverZeros(matrix):
    r, c = matrix.shape
    covered_rows = [False] * r
    covered_cols = [False] * c
    
    # Find all zero positions
    zeros = [(i, j) for i in range(r) for j in range(c) if abs(matrix[i][j]) < 1e-10]
    
    while zeros:
        # Count uncovered zeros in each row and column
        row_zeros = [sum(1 for x, y in zeros if x == i and not covered_cols[y]) for i in range(r)]
        col_zeros = [sum(1 for x, y in zeros if y == j and not covered_rows[x]) for j in range(c)]
        
        # Choose the maximum between rows and columns
        max_row = max(row_zeros)
        max_col = max(col_zeros)
        
        if max_row >= max_col:
            row_idx = row_zeros.index(max_row)
            covered_rows[row_idx] = True
        else:
            col_idx = col_zeros.index(max_col)
            covered_cols[col_idx] = True
            
        # Update remaining zeros
        zeros = [(x, y) for x, y in zeros if not (covered_rows[x] or covered_cols[y])]
    
    return covered_rows, covered_cols, sum(covered_rows) + sum(covered_cols)

def finalAssignments(matrix):
    n = matrix.shape[0]
    assignment = np.full(n, -1)
    
    # Create masks for assigned rows and columns
    assigned_rows = np.zeros(n, dtype=bool)
    assigned_cols = np.zeros(n, dtype=bool)
    
    # First pass: Assign unique zeros where possible
    for i in range(n):
        zero_positions = np.where(abs(matrix[i, :]) < 1e-10)[0]
        for j in zero_positions:
            if not assigned_cols[j]:
                assignment[i] = j
                assigned_rows[i] = True
                assigned_cols[j] = True
                break
    
    # Second pass: Handle remaining unassigned positions
    for i in range(n):
        if assignment[i] == -1:
            # Find the minimum value in unassigned columns
            available_cols = ~assigned_cols
            if np.any(available_cols):
                min_val_col = np.argmin(matrix[i, available_cols])
                actual_col = np.where(available_cols)[0][min_val_col]
                assignment[i] = actual_col
                assigned_cols[actual_col] = True
    
    return assignment

def hungarianMethod(matrix):
    # Make a copy to avoid modifying the original
    working_matrix = matrix.copy()
    
    # Step 1: Reduction
    working_matrix = Reduction(working_matrix)
    
    # Step 2: Cover zeros and adjust matrix until optimal
    covered_rows, covered_cols, lines_used = coverZeros(working_matrix)
    
    while lines_used < matrix.shape[0]:
        working_matrix = adjustMatrix(working_matrix, covered_rows, covered_cols)
        covered_rows, covered_cols, lines_used = coverZeros(working_matrix)
    
    # Step 3: Make final assignments
    return finalAssignments(working_matrix)

def constructMatrix(teammate_positions, formation_positions):
    num_players = len(teammate_positions)
    cost_matrix = np.zeros((num_players, num_players))
    
    for i in range(num_players):
        teammate_pos = np.array(teammate_positions[i])
        for j in range(num_players):
            formation_pos = np.array(formation_positions[j])
            cost_matrix[i, j] = np.linalg.norm(teammate_pos - formation_pos)
    
    return np.round(cost_matrix, 2)

def role_assignment(teammate_positions, formation_positions):
    cost_matrix = constructMatrix(teammate_positions, formation_positions)
    print("Cost matrix before reduction:\n", cost_matrix)
    
    assignment = hungarianMethod(cost_matrix)
    print("assignment:\n", assignment)
    
    point_preferences = {}
    for row, col in enumerate(assignment):
        if col != -1:  # Only create preference if there's a valid assignment
            x, y = formation_positions[int(col)]
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
