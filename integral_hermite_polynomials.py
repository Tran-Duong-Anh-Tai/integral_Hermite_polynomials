
import numpy as np
import math

# Output file
outfile_W = open('Wijkl.txt', 'w')
outfile_Y = open('Yijkl.txt', 'w')
outfile_U = open('Uijklmn.txt', 'w')

# Creating dictionary Wijkl, Yijkl, Uijklmn
Wijkl   = {} 
Yijkl   = {}
Uijklmn = {} 

# Input value Morb
Morb = 6

# Function to lookup Wijkl
def calc_W(Wijkl,index):
    if max(index) > 0 and min(index) >= 0 :
        index.sort(reverse=True)
        return Wijkl.get(''.join(map(str,index)),0)
    elif index[0] == index[1] == index[2] == index[3] == 0  :
        return 1 / math.sqrt(2 * math.pi)
    else:
        return 0

# Function to lookup Uijklmn
def calc_U(Uijklmn,index):
    if max(index) > 0 and min(index) >= 0 :
        index.sort(reverse=True)
        return Uijklmn.get(''.join(map(str,index)),0)
    elif index[0] == index[1] == index[2] == index[3] == index[4] == index[5] == 0 :
        return 1 / (math.pi*math.sqrt(3))
    else:
        return 0
        
# Printing W_{0,0,0,0}
outfile_W.write(f"W_{(0,0,0,0)} = {1 / math.sqrt(2 * math.pi)} \n")
print(f"W_{(0,0,0,0)} = {1 / math.sqrt(2 * math.pi)} \n")

# Calculating W_ijkl
for i in range(1, Morb + 1, 1):
    for j in range(0, i + 1, 1):
        for k in range(0, j + 1, 1):
            for l in range(0, k + 1, 1):
                if ( i + j + k + l ) % 2 == 0: 
                    # Creating index
                    index = str(i) + str(j) + str(k) + str(l)

                    # Calculating the value of the integral using a recursive formula
                    Wijkl[index] = 0.5 * (-math.sqrt((i - 1) / i) * calc_W(Wijkl, [i - 2, j, k, l]) +
                        math.sqrt(j / i) * calc_W(Wijkl, [i - 1, j - 1, k, l]) +
                        math.sqrt(k / i) * calc_W(Wijkl, [i - 1, j, k - 1, l]) +
                        math.sqrt(l / i) * calc_W(Wijkl, [i - 1, j, k, l - 1]) )

                    # Printing the calculated value
                    print(f"W_{(i,j,k,l)} = {Wijkl[index]} \n")

                    # Writing the calculated value
                    outfile_W.write(f"W_{(i,j,k,l)} = {Wijkl[index]}\n")


# Calculating Y_ijkl
for i in range(0, Morb + 1,1):
    for j in range(0, i + 1,1):
        for k in range(0, j + 1,1):
            for l in range(0, k + 1,1):
                if ( i + j + k + l ) % 2 == 0 and i!= j and k!=l : 
                    # Creating index
                    index = str(i) + str(j) + str(k) + str(l)

                    # Calculate the value of the integral using a recursive formula
                    Yijkl[index]= 2 * (math.sqrt(i*k) * calc_W(Wijkl, [i - 1, j, k-1, l]) -
                        math.sqrt(i*l) * calc_W(Wijkl, [i - 1, j, k, l - 1]) -
                        math.sqrt(j*k) * calc_W(Wijkl, [i, j - 1, k - 1, l]) +
                        math.sqrt(j*l) * calc_W(Wijkl, [i, j - 1, k, l - 1]))

                    # Printing the calculated value
                    print(f"Y_{(i,j,k,l)} = {Yijkl[index]} \n")

                    # Writing the calculated value
                    outfile_Y.write(f"Y_{(i,j,k,l)} = {Yijkl[index]}\n")

# Printing U_{0,0,0,0,0,0}
outfile_U.write(f"U_{0,0,0,0,0,0} = {1 / (math.pi*math.sqrt(3))} \n")
print(f"U_{0,0,0,0,0,0} = {1 / (math.pi*math.sqrt(3))} \n")

# Calculating U_ijklmn
for i in range(1, Morb + 1,1):
    for j in range(0, i + 1,1):
        for k in range(0, j + 1,1):
            for l in range(0, k + 1,1):
                for m in range(0, l + 1,1):
                    for n in range (0, m + 1,1):
                        if ( i + j + k + l + m + n) % 2 == 0: 
                            # Creating index
                            index = str(i) + str(j) + str(k) + str(l) + str(m) + str(n)

                            # Calculate the value of the integral using a recursive formula
                            Uijklmn[index] = (1/3) * (-2*math.sqrt((i - 1)/i) * calc_U(Uijklmn, [i - 2, j, k, l, m, n]) +
                                math.sqrt(j/i) * calc_U(Uijklmn, [i - 1, j - 1, k, l, m, n]) +
                                math.sqrt(k/i) * calc_U(Uijklmn, [i - 1, j, k - 1, l, m, n]) +
                                math.sqrt(l/i) * calc_U(Uijklmn, [i - 1, j, k, l - 1, m, n]) +
                                math.sqrt(m/i) * calc_U(Uijklmn, [i - 1, j, k, l, m - 1, n]) + 
                                math.sqrt(n/i) * calc_U(Uijklmn, [i - 1, j, k, l, m, n - 1])) 

                            # Printing the calculated value
                            print(f"U_{i,j,k,l,m,n} = {Uijklmn[index]} \n")

                            # Writing the calculated value
                            outfile_U.write(f"U_{i,j,k,l,m,n} = {Uijklmn[index]}\n")

# Close the output file
outfile_W.close()
outfile_Y.close()
outfile_U.close()