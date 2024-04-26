import numpy as np

def keccak_state_to_arr(state_array, w):
    keccak_res_array = np.arange(25*w)
    for i in range(len(keccak_res_array)):
        x = (i//w)%5
        y = (i//(w*5))
        z = i%w
        # reverse the state while filling an array
        keccak_res_array[len(keccak_res_array)-(i+1)] = state_array[x][y][z]
    return keccak_res_array

def kec_bin_to_hex(state_array):
    str_res = ''
    hex_num = ''
    for i in range(0, len(state_array), 4):
        hex_num = hex(int((bin(state_array[i])[2:] + bin(state_array[i+1])[2:] + \
                           bin(state_array[i+2])[2:] + bin(state_array[i+3])[2:]), 2))[2:]
        str_res += hex_num
    return str_res

# create keccak design
def theta(in_state, w):
    res_state = np.zeros((5, 5, w), dtype=np.uint16)
    col_right = 0
    col_left = 0
    for k in range(w):  # lane
        for i in range(5):  # row
            col_left = in_state[i-1][0][k] ^ in_state[i-1][1][k] ^ in_state[i-1][2][k] ^ \
                in_state[i-1][3][k] ^ in_state[i-1][4][k]
            col_right = in_state[(i+1)%5][0][k-1] ^ in_state[(i+1)%5][1][k-1] ^ in_state[(i+1)%5][2][k-1] \
                ^ in_state[(i+1)%5][3][k-1] ^ in_state[(i+1)%5][4][k-1]
            for j in range(5):  # column 
                res_state[i][j][k] = in_state[i][j][k] ^ col_left ^ col_right
    return res_state

def phi(in_state, w):
    res_state = np.zeros((5, 5, w), dtype=np.uint16)
    for i in range(5): # x-axis
        for j in range(5): # y-axis
            res_state[j][(2*i+3*j)%5] = in_state[i][j] # 4,4 <= 1,4
    return res_state

def rho(in_state, w):
    res_state = np.zeros((5, 5, w), dtype=np.uint16)
    #rot_mtrx = [
    #    [0, 1, 190, 28, 91],
    #    [36, 300, 6, 55, 276],
    #    [3, 10, 171, 153, 231],
    #    [105, 45, 15, 21, 136],
    #    [210, 66, 253, 120, 78]
    #]
    rot_mtrx = [
        [0, 36, 3, 105, 210],
        [1, 300, 10, 45, 66],
        [190, 6, 171, 15, 253],
        [28, 55, 153, 21, 120],
        [91, 276, 231, 136, 78]
    ]
    #print(rot_mtrx)
    for i in range(5): # x-axis
        for j in range(5): # y-axis
            res_state[i][j][rot_mtrx[i][j]%w:] = in_state[i][j][0:w-(rot_mtrx[i][j]%w)]
            res_state[i][j][0:rot_mtrx[i][j]%w:] = in_state[i][j][w-(rot_mtrx[i][j]%w):]
    return res_state

def chi(in_state, w):
    res_state = np.zeros((5, 5, w), dtype=np.uint16)
    for i in range(5): # x-axis
        for j in range(5):  # y-axis
            for k in range(w):  # z-axis
                res_state[i][j][k] = in_state[i][j][k] ^ (~in_state[(i+1)%5][j][k] & in_state[(i+2)%5][j][k])
    return res_state

def iota(in_state, w, cur_round):
    res_state = np.zeros((5, 5, w), dtype=np.uint16)
    res_state = in_state
    rnd_cnst = [
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0],
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0],
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,1,1],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1],
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,1,1],
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,1],
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,1],
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1],
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0],
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0],
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1],
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0]
    ]
    
    rnd_cnst_reversed = rnd_cnst[cur_round][:64-w-1:-1]
    res_state[0][0] = res_state[0][0] ^ rnd_cnst_reversed
    return res_state

def keccak_f(in_state, round_num, w):
    res_state = np.zeros((5, 5, w), dtype=np.uint16)
    res_state = in_state

    # Test intermediates
    interm_array = np.arange(25*w)
    for i in range(round_num):
        #print("Round " + str(i) + ":")
        #print("Theta:")
        res_state = theta(res_state, w)
        interm_array = keccak_state_to_arr(res_state, w)
        #print(kec_bin_to_hex(interm_array))
        #print("Rho:")
        res_state = rho(res_state, w)
        interm_array = keccak_state_to_arr(res_state, w)
        #print(kec_bin_to_hex(interm_array))
        #print("Phi:")
        res_state = phi(res_state, w)
        interm_array = keccak_state_to_arr(res_state, w)
        #print(kec_bin_to_hex(interm_array))
        #print("Chi:")
        res_state = chi(res_state, w)
        interm_array = keccak_state_to_arr(res_state, w)
        #print(kec_bin_to_hex(interm_array))
        #print("Iota:")
        res_state = iota(res_state, w, i)
        interm_array = keccak_state_to_arr(res_state, w)
        #print(kec_bin_to_hex(interm_array))
    return res_state

# Main
# Test vectors for 400 bits permutation:
# input:
# 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
# output:
# 37 E5 D6 D5 E7 DB F3 AA C7 9B 7D CA B2 86 EC FD 2C 69 5B 4E B1 67 AD 15 F7 A7 6F A6 FF 67 8A 3F 99 2F C2 E2 6B 65 31 5F A6 5B 29 CA 24 C2 5C B8 7C 09

# Hardware equivalent [399:0] MSB - LSB:
# in: b38e714c71ff58ee633e6e146d955d8159a15246e7760ec9dc7560520337bf8feff078705bd1eca0e89f14f50fa940ac09f5
# out: 097cb85cc224ca295ba65f31656be2c22f993f8a67ffa66fa7f715ad67b14e5b692cfdec86b2ca7d9bc7aaf3dbe7d5d6e537

#w = 8 #keccak 200
f = 400 # keccak-f[200 || 400 || 1600]
w = 16 #keccak 400
l = 4 # 2**l = w

#str_dec = "0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
# Encoding order: MSB - to - LSB (big endian?) HDL-like
str_dec = "b38e714c71ff58ee633e6e146d955d8159a15246e7760ec9dc7560520337bf8feff078705bd1eca0e89f14f50fa940ac09f5"

rounds = 12 + 2*l
keccak_array = np.arange(25*w)
keccak_res_array = np.arange(25*w)
keccak_state = np.zeros((5, 5, w), dtype=np.uint16)

#arrange the state
str_dec_bin_b = bin(int(str_dec,16)).zfill(f)
str_dec_bin = str_dec_bin_b.replace('b','0') # binary encoding

# Reversed HDL-like [0:399] LSB - MSB (reverse in bits)
for i in range(25*w):
    keccak_array[i] = str_dec_bin[len(str_dec_bin)-i-1]

#print(kec_bin_to_hex(keccak_array)) # reversed input

x = y = z = 0
for i in range(len(keccak_array)):
    x = (i//w)%5
    y = (i//(w*5))
    z = i%w
    keccak_state[x][y][z] = keccak_array[i] #s[w(5y+x)+z]

# keccak_state is according to the documentation LSB - 1st
keccak_result = keccak_f(keccak_state, rounds, w)

keccak_res_array = keccak_state_to_arr(keccak_result, w)

str_res = kec_bin_to_hex(keccak_res_array)
print(str_res)
