# Linear congruential generator (LCG) in RISC-V assembly, written for RARS.
#
# Same algorithm as python/lcg.py: x = (a*x + c) mod m, printed until the
# seed recurs. Parameters are the small example m=17, a=7, c=3, seed=6.
#
# Register map:
#   x6  (t1)  - m, the modulus
#   x7  (t2)  - a, the multiplier
#   x28 (t3)  - c, the increment
#   x29 (t4)  - seed (used at the end to detect the cycle closing)
#   x30 (t5)  - x, the current generator state
#   x31 (t6)  - address of the ' ' separator byte
#
# Syscalls (RARS convention): a0=1 + x11=int   -> print integer
#                              a0=4 + x11=addr -> print null-free string/byte
#                              a0=10           -> exit

.data
space: .byte ' '

.text
li   x7,  7          # a = 7
li   x28, 3           # c = 3
li   x29, 6           # seed = 6
mv   x30, x29         # x = seed
la   x31, space

li   a0, 1            # print the seed first
addi x11, x29, 0
ecall
li   x6, 17           # m = 17
li   a0, 4
addi x11, x31, 0
ecall

loop:
    mul  x30, x30, x7      # x = x * a
    add  x30, x30, x28     # x = x + c
    rem  x30, x30, x6      # x = x mod m
    beq  x30, x29, exit    # stop once the seed recurs (full period reached)
    li   a0, 1
    addi x11, x30, 0
    ecall
    li   a0, 4
    addi x11, x31, 0
    ecall
    j loop

exit:
    li a0, 10
    ecall
