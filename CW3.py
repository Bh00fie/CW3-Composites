import numpy as np

#Example done in class using Rule of Mixture and HT Method
#Data given:
S3 = 0.2
S12 = 1
WFibre = 0.6 #Fibre weight Fraction
thickness = 0.001 #m For each layer
# thickness = 0.0002 #m For each layer (Question 2)
# thickness = 0.000125 #m For each Layer (Question 5)

# Laminate Configuration for Question 1
laminate_config = [0, 90, 90, 0]
# laminate_config = [-45, +45, 0, 90, 90, 0, 0, 90, 90, 0, +45, -45] # (Question 2)
# laminate_config = [90, 45,-45, 0] # (Question 5)
nLayers = len(laminate_config)

print(f'\nThe laminate configuration for this laminate is: {laminate_config}')
print(f'The thickness of the laminate is: {nLayers*thickness} m or {nLayers*thickness*1e3} mm, each laminate in the configuration has the same thickness, thus: {thickness} m')

# Calculating the bottom and top thickness of the laminate
position = [-nLayers * thickness / 2 + i * thickness for i in range(nLayers + 1)]

# Creating an empty list to store Qbar, A1, B1, D1 matrices for each angle
Qbar_list = []
A1_list = []
B1_list = []
D1_list = []

# Use a loop to keep track of the current layer
for i in range(nLayers):  
    angle = laminate_config[i] # Changing index/angle for every loop iteration
    theta = angle * (np.pi / 180)  # Calculate theta based on the angle converted 
    
    # Defining Carbon Fibre Properties
    pFibre = 1870 #km/m^3 - Density
    EFibre = 310 #GN/m^2 - Young's Modulus
    vFibre = 0.23 # - Poisson Ratio
    #Epoxy
    pMatrix = 1200 #km/m^3 - Density
    EMatrix = 2.4 #GN/m^2 - Young's Modulus
    vMatrix = 0.33 # - Poisson Ratio

    WMatrix = 1 - WFibre # Calculating the Weight Fraction value of the Matrix
    VFibre = np.round((WFibre / pFibre) / ((WFibre / pFibre) + (WMatrix / pMatrix)), 3)  # Calculating the Volume Fraction of the Fibre, rounded to 3 digits using np.round
    VMatrix = 1 - VFibre # Calculating the Volume Fraction of the Matrix

    # RULE OF MIXTURE
    #Young Modulus Calculation
    E1 = EMatrix*VMatrix + EFibre*VFibre  # Young's Modulus in the 1 direction
    E2_RoM = (EFibre*EMatrix)/((VFibre*EMatrix)+(VMatrix*EFibre)) # Young's Modulus in the 2 direction

    #Shear Modulus Calculation
    GFibre = EFibre/(2*(1+vFibre)) # Shear Modulus of the Fibre
    GMatrix = EMatrix/(2*(1+vMatrix)) # Shear Modulus of the Matrix
    G12_RoM = (GFibre*GMatrix)/((GFibre*VMatrix) + (GMatrix*VFibre)) # Shear Modulus of the Laminate
    
    #Major Poisson Ratio
    v12 = vFibre*VFibre + vMatrix*VMatrix
    #Minor Poisson Ratio
    v21_RoM = v12 *(E2_RoM/E1)  # Minor Poisson Ration
    
    # HALPIN TSAI 
    nE = ((EFibre/EMatrix)-1)/((EFibre/EMatrix) + S3) 
    E2 = EMatrix * ((1+(S3*nE*VFibre))/(1-(nE*VFibre))) # Young's Modulus in the 2 direction
     
    nG = ((GFibre/GMatrix) - 1)/((GFibre/GMatrix) + 1)
    G12 = GMatrix * ((1 + S12*nG*VFibre)/(1 - nG*VFibre)) # Shear Modulus of the Laminate

    v21 = v12 *(E2/E1)  # Minor Poisson Ration

    # # Question 2 - Material Properties
    # E1 = 54 #GPa
    # E2 = 18 #GPa
    # v12 = 0.28
    # v21 = v12 *(E2/E1)  # Minor Poisson Ration
    # G12 = 6 #GPa
    
    # Qmatrix
    Q11 = E1 / (1 - v12 * v21)
    Q22 = E2 / (1 - v12 * v21)
    Q12 = v12 * E2 / (1 - v12 * v21)
    Q66 = G12
    Q = np.array([[Q11, Q12, 0],
                [Q12, Q22, 0],
                [0, 0, Q66]])
    # Q = np.round(Q, 3)
    # print(f'\n Q in GN/m^2: \n{Q}')

    # Tmatrix
    T11 = np.cos(theta) ** 2
    T12 = np.sin(theta) ** 2
    T13 = 2 * np.sin(theta) * np.cos(theta)
    T21 = np.sin(theta) ** 2
    T22 = np.cos(theta) ** 2
    T23 = -2 * np.sin(theta) * np.cos(theta)
    T31 = -np.sin(theta) * np.cos(theta)
    T32 = np.sin(theta) * np.cos(theta)
    T33 = np.cos(theta) ** 2 - np.sin(theta) ** 2

    T = np.array([[T11, T12, T13],
                [T21, T22, T23],
                [T31, T32, T33]]) 

    # Inverse of T matrix
    T_inv = np.linalg.inv(T)

    # Transpose of the inverse matrix
    T_inv_transpose = np.transpose(T_inv)

    # Qbar Matrix calculate dot product with T_Inv, Q, T_inv_transpose
    Qbar = T_inv @ Q @ T_inv_transpose*1e9
    Qbar = np.round(Qbar, 3)

    # Converting the NumPy array to a regular Python list
    # Qbar_list.append(Qbar.tolist())
    Qbar_list.append(Qbar)

    # A Matrix
    A1 = Qbar * (position[i + 1] - position[i])
    A1_list.append(A1)

    # B Matrix    
    B1 = (1/2)*Qbar * ((position[i + 1]**2) - (position[i]**2))
    B1_list.append(B1)
    
    # D Matrix    
    D1 = (1/3)*Qbar * ((position[i + 1]**3) - (position[i]**3))
    D1_list.append(D1)

print('\nTHE MATERIAL PROPERTIES ARE:')   
print(f"Fibre Young's Modulus: {EFibre} GN/m^2")
print(f"Fibre Poisson's Ratio: {vFibre} Kg/m^3")
print(f"Fibre Volume Factor: {VFibre}")
print(f"Matrix Young's Modulus: {EMatrix} GN/m^2")
print(f"Matrix Poisson's Ratio: {vMatrix} Kg/m^3")
print(f"Matrix Volume Factor: {VMatrix}")
print(f"Fibre Shear Stress: {GFibre} GN/m^2")
print(f"Matrix Shear Stress: {GMatrix} GN/m^2 \n")
print('-----------------')
print('\n Using the RULE OF MIXTURE:')
print(f"E11 is: {E1} GN/m^2")
print(f"E22 is: {E2_RoM} GN/m^2")
print(f"v12 is: {v12}")
print(f"v21 is: {v21_RoM}")
print(f"G12 is: {G12_RoM} GN/m^2")
print('\n Using HALPIN-TSAI:')
print(f"E11 is: {E1} GN/m^2")
print(f"E22 is: {E2} GN/m^2")
print(f"v12 is: {v12}")
print(f"v21 is: {v21}")
print(f"G12 is: {G12} GN/m^2 \n")
print('----------------- \n')

# Set numpy print options to display numbers without scientific notation, NOT ROUNDING
np.set_printoptions(suppress=True, formatter={'float': lambda x: '0' if x == 0.0 else '{:0.3f}'.format(x)})

# Print the Qbar matrices in a human-readable format
for i, Qbar_matrix in enumerate(Qbar_list):
    print(f"Qbar for angle {laminate_config[i]} degrees in N/m^2:")
    for row in Qbar_matrix:
        print(row)
    print()
    
# Defining one single A, B, D matrix for the laminate (through summing A, B, D for each lamina)
A = np.sum(A1_list, axis=0)
B = np.sum(B1_list, axis=0)
D = np.sum(D1_list, axis=0)
# A = np.sum(A1_list, axis=0)/1e9 # To get values in GN/m
# B = np.sum(B1_list, axis=0)/1e9 # To get values in GN
# D = np.sum(D1_list, axis=0)/1e9 # To get values in GN*m

print(f'The A Matrix in N/m is: \n {A}')
print(f'\nThe B Matrix in N is: \n {B}')
print(f'\nThe D Matrix in Nm is: \n {D}')
print(' \n----------------- \n')

# Concatenating the A, B, and D matrices in the specified order
ABBD = np.block([[A, B], [B, D]])
print(f"ABBD Matrix: \n{ABBD} \n ") # Printing the ABBD matrix

# Inverse of ABBD Matrix
ABBD_inv = np.linalg.inv(ABBD) 
print(f"ABBD_inv Matrix: \n{ABBD_inv}") # Printing the ABBD_inv matrix
print('\n ----------------- \n')

# Define the different directions
directions = ["Ex", "Ey", "Gxy"]
print("Engineering Properties:")
for direction in directions:
    NM = [0] * 6     # Initialize the NM Matrix
    
    if direction == "Ex":
        NM[0] = 1
    elif direction == "Ey":
        NM[1] = 1
    elif direction == "Gxy":
        NM[2] = 1

    NM = np.array(NM).reshape(6, 1)

    result = np.dot(ABBD_inv, NM) # Calculate the strain in different directions -> [ϵx, ϵy, γxy, κx, κy, κxy] vector
    # print(result)

    ex, ey, yxy, kx, ky, kxy = result # Assigning each strain to its value in the result matrix

    if direction == "Ex":
        Sx = ((A[0, 0] * ex + A[0, 1] * ey + A[0, 2] * yxy + B[0, 0] * kx + B[0, 1] * ky + B[0, 2] * kxy)) / (nLayers * thickness)
        Ex = Sx / ex
        print(f"Ex: {Ex} N/m^2")
    elif direction == "Ey":
        Sy = ((A[1, 0] * ex + A[1, 1] * ey + A[1, 2] * yxy + B[1, 0] * kx + B[1, 1] * ky + B[1, 2] * kxy)) / (nLayers * thickness)
        Ey = Sy / ey
        print(f"Ey: {Ey} N/m^2")
    elif direction == "Gxy":
        Sxy = ((A[2, 0] * ex + A[2, 1] * ey + A[2, 2] * yxy + B[2, 0] * kx + B[2, 1] * ky + B[2, 2] * kxy)) / (nLayers * thickness)
        Gxy = Sxy / yxy
        print(f"Gxy: {Gxy} N/m^2 \n")