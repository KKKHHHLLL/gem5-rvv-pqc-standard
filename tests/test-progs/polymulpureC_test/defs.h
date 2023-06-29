/******************************************************************************

 * BIKE -- Bit Flipping Key Encapsulation

 *

 * Copyright (c) 2017 Nir Drucker, Shay Gueron, Rafael Misoczki, Tobias Oder, Tim Gueneysu

 * (drucker.nir@gmail.com, shay.gueron@gmail.com, rafaelmisoczki@google.com, tobias.oder@rub.de, tim.gueneysu@rub.de)

 *

 * Permission to use this code for BIKE is granted.

 * Redistribution and use in source and binary forms, with or without

 * modification, are permitted provided that the following conditions are met:

 *

 * * Redistributions of source code must retain the above copyright notice,

 *   this list of conditions and the following disclaimer.

 *

 * * Redistributions in binary form must reproduce the above copyright

 *   notice, this list of conditions and the following disclaimer in the

 *   documentation and/or other materials provided with the distribution.

 *

 * * The names of the contributors may not be used to endorse or promote

 *   products derived from this software without specific prior written

 *   permission.

 *

 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ""AS IS"" AND ANY

 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE

 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR

 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS CORPORATION OR

 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,

 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,

 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR

 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF

 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING

 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS

 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 ******************************************************************************/

#include <stdint.h>

#ifndef __DEFS_H_INCLUDED__

#define __DEFS_H_INCLUDED__



////////////////////////////////////////////

//         BIKE main parameters

///////////////////////////////////////////



// UNCOMMENT TO SELECT THE NIST SECURITY LEVEL 1, 3 OR 5:

#define PARAM64 // NIST LEVEL 1

//#define PARAM96 // NIST LEVEL 3

//#define PARAM128 // NIST LEVEL 5



// UNCOMMENT TO ENABLE BANDWIDTH OPTIMISATION FOR BIKE-3:

//#define BANDWIDTH_OPTIMIZED



// BIKE shared-secret size:

#define ELL_BITS  256ULL

#define ELL_SIZE (ELL_BITS/8)



////////////////////////////////////////////

// Implicit Parameters (do NOT edit below)

///////////////////////////////////////////



// select the max between a and b:

#define MAX(a,b) ((a)>(b))?(a):(b)



// LEVEL-5 Security parameters:

#ifdef PARAM128

#define LEVEL 5

#define R_BITS 40973ULL

#define DV     137ULL

#define T1     264ULL

#define VAR_TH_FCT(x) (MAX(17.8785 + 0.00402312 * (x), 69))

// Parameters for BGF Decoder:

#define tau 3

#define NbIter 5

#define BLOCK_BITS 65536

// LEVEL-3 Security parameters:

#elif defined(PARAM96)

#define LEVEL 3

#define R_BITS 24659ULL

#define DV     103ULL

#define T1     199ULL

#define VAR_TH_FCT(x) (MAX(15.2588 + 0.005265 * (x), 52))

// Parameters for BGF Decoder:

#define tau 3

#define NbIter 5

#define BLOCK_BITS 32768

// LEVEL-1 security parameters:

#elif defined(PARAM64)

#define LEVEL 1

#define R_BITS 12323ULL

#define DV     71ULL

#define T1     134ULL

#define VAR_TH_FCT(x) (MAX(13.530 + 0.0069722 * (x), 36))

// Parameters for BGF Decoder:

#define tau 3

#define NbIter 5

#define BLOCK_BITS 16384

#endif



// Divide by the divider and round up to next integer:

#define DIVIDE_AND_CEIL(x, divider)  ((x/divider) + (x % divider == 0 ? 0 : 1ULL))



// Round the size to the nearest byte.

// SIZE suffix, is the number of bytes (uint8_t).

#define N_BITS   (R_BITS*2)

#define R_SIZE   DIVIDE_AND_CEIL(R_BITS, 8ULL)

#define R_SIZE_ DIVIDE_AND_CEIL(R_SIZE, 4ULL)

#define N_SIZE   DIVIDE_AND_CEIL(N_BITS, 8ULL)

#define R_DQWORDS DIVIDE_AND_CEIL(R_SIZE, 16ULL)



////////////////////////////////////////////

//             Debug

///////////////////////////////////////////



#ifndef VERBOSE

#define VERBOSE 0

#endif



#if (VERBOSE == 3)

#define MSG(...)     { printf(__VA_ARGS__); }

#define DMSG(...)    MSG(__VA_ARGS__)

#define EDMSG(...)   MSG(__VA_ARGS__)

#define SEDMSG(...)  MSG(__VA_ARGS__)

#elif (VERBOSE == 2)

#define MSG(...)     { printf(__VA_ARGS__); }

#define DMSG(...)    MSG(__VA_ARGS__)

#define EDMSG(...)   MSG(__VA_ARGS__)

#define SEDMSG(...)

#elif (VERBOSE == 1)

#define MSG(...)     { printf(__VA_ARGS__); }

#define DMSG(...)    MSG(__VA_ARGS__)

#define EDMSG(...)

#define SEDMSG(...)

#else

#define MSG(...)     { printf(__VA_ARGS__); }

#define DMSG(...)

#define EDMSG(...)

#define SEDMSG(...)

#endif



////////////////////////////////////////////

//              Printing

///////////////////////////////////////////



// Show timer results in cycles.

#define RDTSC



//#define PRINT_IN_BE

//#define NO_SPACE

//#define NO_NEWLINE



////////////////////////////////////////////

//              Testing

///////////////////////////////////////////

#define NUM_OF_CODE_TESTS       100ULL

#define NUM_OF_ENCRYPTION_TESTS 100ULL



////////////////////////////////////////////

//              MY BIKE CODE RELATED

///////////////////////////////////////////

//for RVV relization we now only consider using 32-bit chunks

#define SEW 32

#define SEW_LOG 5

#define T_TYPE uint32_t

#define U_TYPE uint64_t

#endif



#define B SEW  //chunk size

#define S SEW  //divstep step, S is less equal than B

#define R_SIZE_PRIME ((R_BITS/SEW) + (R_BITS % SEW == 0 ? 0 : 1ULL))//array length to store the binary polynomial coefficients

#define CHUNK_MASK ((1ULL << SEW) - 1)

#define FINAL_CHUNK_MASK ((1ULL << (R_BITS & (SEW - 1))) - 1)



#define IN

#define OUT



/////////////////////////////////////////////////

//              MY BIKE ITI POLY INVERSE RELATED

/////////////////////////////////////////////////

#define K_SQR_THR (64)



// k-squaring is computed by a permutation of bits of the input polynomial,

// as defined in [1](Observation 1). The required parameter for the permutation

// is l = (2^k)^-1 % R.

// Therefore, there are two sets of parameters for every exponentiation:

//   - exp0_k and exp1_k

//   - exp0_l and exp1_l



// Exponentiation 0 computes f^2^2^(i-1) for 0 < i < MAX_I.

// Exponentiation 1 computes f^2^((r-2) % 2^i) for 0 < i < MAX_I,

// only when the i-th bit of (r-2) is 1. Therefore, the value 0 in

// exp1_k[i] and exp1_l[i] means that exp1 is skipped in i-th iteration.



// To quickly generate all the required parameters in Sage:

//   r = DESIRED_R

//   max_i = floor(log(r-2, 2)) + 1

//   exp0_k = [2^i for i in range(max_i)]

//   exp0_l = [inverse_mod((2^k) % r, r) for k in exp0_k]

//   exp1_k = [(r-2)%(2^i) if ((r-2) & (1<<i)) else 0 for i in range(max_i)]

//   exp1_l = [inverse_mod((2^k) % r, r) if k != 0 else 0 for k in exp1_k]



#if(LEVEL == 1)

// The parameters below are hard-coded for R=12323

// MAX_I = floor(log(r-2)) + 1

#  define MAX_I (14)

#  define EXP0_K_VALS 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192

#  define EXP0_L_VALS 6162, 3081, 3851, 5632, 22, 484, 119, 1838, 1742, 3106, 10650, 1608, 10157, 8816

#  define EXP1_K_VALS 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 33, 4129

#  define EXP1_L_VALS 0, 0, 0, 0, 0, 6162, 0, 0, 0, 0, 0, 0, 242, 5717



#elif(LEVEL == 3)

// The parameters below are hard-coded for R=24659

// MAX_I = floor(log(r-2)) + 1

#  define MAX_I (15)

#  define EXP0_K_VALS 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384

#  define EXP0_L_VALS 12330, 6165, 7706, 3564, 2711, 1139, 15053, 1258, 4388, 20524, 9538, 6393, 10486, 1715, 6804

#  define EXP1_K_VALS 0, 0, 0, 0, 1, 0, 17, 0, 0, 0, 0, 0, 0, 81, 8273

#  define EXP1_L_VALS 0, 0, 0, 0, 12330, 0, 13685, 0, 0, 0, 0, 0, 0, 23678, 19056



#else

// The parameters below are hard-coded for R=40973

// MAX_I = floor(log(r-2)) + 1

#  define MAX_I (16)

#  define EXP0_K_VALS 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768

#  define EXP0_L_VALS 20487, 30730, 28169, 9443, 13001, 12376, 8302, 6618, 38760, 21582, 1660, 10409, 14669, 30338, 17745, 7520

#  define EXP1_K_VALS 0, 1, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 11, 0, 8203

#  define EXP1_L_VALS 0, 20487, 0, 15365, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6302, 0, 10058



#endif



#define BIT(len)       (1ULL << (len))

#define MASK(len)      (BIT(len) - 1)

#define BYTES_IN_QWORD 0x8



#define R_BLOCKS        DIVIDE_AND_CEIL(R_BITS, BLOCK_BITS)

#define R_QWORDS DIVIDE_AND_CEIL(R_BITS, (8 * BYTES_IN_QWORD))

#define R_PADDED        (R_BLOCKS * BLOCK_BITS)

#define R_PADDED_BYTES  (R_PADDED / 8)

#define R_PADDED_QWORDS (R_PADDED / 64)



#define LAST_R_QWORD_LEAD  (R_BITS & MASK(6))

#define LAST_R_QWORD_TRAIL (64 - LAST_R_QWORD_LEAD)

#define LAST_R_QWORD_MASK  MASK(LAST_R_QWORD_LEAD)



#define LAST_R_BYTE_LEAD  (R_BITS & MASK(3))

#define LAST_R_BYTE_TRAIL (8 - LAST_R_BYTE_LEAD)

#define LAST_R_BYTE_MASK  MASK(LAST_R_BYTE_LEAD)



#define REG_T uint64_t

#define REG_QWORDS (sizeof(REG_T) / sizeof(uint64_t)) // NOLINT



#  define LOAD(mem)       (mem)[0]

#  define STORE(mem, val) (mem)[0] = val



#  define SLLI_I64(a, imm) ((a) << (imm))

#  define SRLI_I64(a, imm) ((a) >> (imm))



#define LSB3(x) ((x)&7)



// The secure buffer size required for Karatsuba is computed by:

//    size(n) = 3*n/2 + size(n/2) = 3*sum_{i}{n/2^i} < 3n

#define SECURE_BUFFER_QWORDS (3 * R_PADDED_QWORDS)

//__TYPES_H_INCLUDED__
