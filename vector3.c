/****************************************************************************
 * examples/acceleration_calibration/vector.c
 *
 *   Copyright (C) 2017 Giorgio Groß. All rights reserved.
 *   Author: Giorgio Groß <giorgio.gross@robodev.eu>
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name NuttX nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

/****************************************************************************
 * Included Files
 ****************************************************************************/

#include "vector3.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/****************************************************************************
 * Pre-processor Definitions
 ****************************************************************************/

 #define DECIMAL_PLACES3(x) abs(((int)(((x)-((int)x))*1000)))

// determine where matrices memory is starting
#define m0(M) (float*)(&(M))

// determine where vector memory is starting
#define v0(v) (float*)(&(v))

// traverse a 3x3 matrix of type matrix3x3_t pointed to by base_p (start pointer) 
// and curr_p (pointer to current index) horizontally 
#define TRAVERSE_HORIZONTAL_3X3(curr_p, base_p) \
        (curr_p - base_p) == 8 ? \
        base_p : \
        (curr_p - base_p + 3) % (8 + ((curr_p - base_p + 3) == 8)) + base_p

/****************************************************************************
 * Public Data
 ****************************************************************************/

vector3_t X = {
  .x = 1,
  .y = 0,
  .z = 0
};

vector3_t Y = {
  .x = 0,
  .y = 1,
  .z = 0
};

vector3_t Z = {
  .x = 0,
  .y = 0,
  .z = 1
};

/****************************************************************************
 * Public Functions
 ****************************************************************************/

float enclosed_angle (vector3_t U, vector3_t V)
{
    return acosf(vector3_product(U, V) / (vector3_length(U) * vector3_length(V)));
}

float vector3_product (vector3_t U, vector3_t V)
{
    return U.x * V.x + U.y * V.y + U.z * V.z;
}

void vector3_sum (vector3_t* T, vector3_t U, vector3_t V)
{
    int c;
    for(c = 0; c < 3; c++)
    {
        *((float*)T + c) = *((float*)&U + c) + *((float*)&V + c);
    }
}

float vector3_length (vector3_t U)
{
    return sqrtf(vector3_product(U, U));
}

void vector3_normalize (vector3_t* U)
{
    float length = vector3_length(*U);
    multiply_vector3_scalar(U, *U, 1/length);
}

void multiply3x3_matrix_matrix (matrix3x3_t* res, matrix3x3_t A, matrix3x3_t B)
{
    float* a = m0(A);
    float* b = m0(B);
    float* r = m0(*res);

    int c;
    vector3_t row;
    for(c = 1; c <= 9; c++)
    {   
        // copy row of A into vector 'row'
        int i;
        for(i = 0; i < 3; i++) *(v0(row) + i) = *(a + 3 * i);
        *r = multiply_vector3_vector3(row, *((vector3_t*) b));

        a += ((c - 1) % 3) == 0;  // incs only if c==3 or c==6 / iterates through row vectors of A
        b = m0(B) + (c % 3) * 3;  // iterates through columns of B
        r = TRAVERSE_HORIZONTAL_3X3(r, m0(*res));
    }
}

void multiply3x3_matrix_vector (vector3_t* res, matrix3x3_t A, vector3_t V)
{
    float* a = m0(A);
    float* r = v0(*res);

    int c;
    vector3_t row;
    for(c = 0; c < 3; c++)
    {
        a += (c>0);
        r += (c>0);

        row.x = *a;
        row.y = *(a + 3);
        row.z = *(a + 6);

        *r = multiply_vector3_vector3(V, row);
    }
}

float multiply_vector3_vector3 (vector3_t U, vector3_t V)
{
    float* u = v0(U);
    float* v = v0(V);
    float res = 0;

    int c;
    for(c = 0; c < 3; c++)
    {
        u += (c>0);
        v += (c>0);

        res += (*u) * (*v);
    }
    return res;
}

void multiply_vector3_scalar (vector3_t* T, vector3_t U, float s)
{
    int c;
    for(c = 0; c < 3; c++)
    {
        *((float*)T + c) = *((float*)&U + c) * s;
    }
}

void orthonormalize3x3(matrix3x3_t* A)
{
    float scalarProduct;
    vector3_t summand2;
    vector3_t summand3;
    // eventually generalize for NxN...

    // Gram-Schmidt Algorithm.
    // Step 1
    vector3_normalize(&(A->v1));

    // Step 2
    scalarProduct = vector3_product(A->v1, A->v2);
    multiply_vector3_scalar(&summand2, A->v1, -scalarProduct);  // negative to subtract by sum in the next line
    vector3_sum(&A->v2, A->v2, summand2);
    vector3_normalize(&A->v2);

    // Step 3
    scalarProduct = vector3_product(A->v2, A->v3);
    multiply_vector3_scalar(&summand3, A->v2, -scalarProduct);  // negative to subtract by sum below
    
    scalarProduct = vector3_product(A->v1, A->v3);
    multiply_vector3_scalar(&summand2, A->v1, -scalarProduct);  // negative to subtract by sum below
    
    vector3_sum(&A->v3, A->v3, summand2);
    vector3_sum(&A->v3, A->v3, summand3);
    vector3_normalize(&A->v3);
}

void transpose3x3(matrix3x3_t* A)
{
    float* horPntr = m0(*A);
    float* vertPntr = m0(*A);
    float cachedVal = 0;
    
    // swap 1,1-1,1 ; 2,2-2,2 ; 3,3-3,3 ; 1,2-2,1 ; 2,3-3,2 ; 1,3-3,1
    int c;
    for(c = 0; c < 6; c++)
    {
        cachedVal = *horPntr;
        *horPntr = *vertPntr;
        *vertPntr = cachedVal;

        int i = 0;
        for(i = 0; i < 4; i++)
        {
            horPntr = TRAVERSE_HORIZONTAL_3X3(horPntr, m0(*A));
            vertPntr = (vertPntr - m0(*A) + 1) % 9 + m0(*A);
        }

    }
}

float det3x3(matrix3x3_t A)
{
    float sum = 0;
    int c;

    // adding part
    float* a = m0(A);
    for(c = 0; c < 3; c++)
    {
        sum += *a * *((a - m0(A) + 4) % 9 + m0(A)) * *((a - m0(A) + 8) % 9 + m0(A));
        a = TRAVERSE_HORIZONTAL_3X3(a, m0(A));
    }
    
    // subtracting part
    a = m0(A);
    for(c = 0; c < 3; c++)
    {
        sum -= *a * *((a - m0(A) + 7) % 9 + m0(A)) * *((a - m0(A) + 14) % 9 + m0(A));
        a = TRAVERSE_HORIZONTAL_3X3(a, m0(A));
    }

    return sum;
}

void print3x3_matrix(matrix3x3_t A)
{
    float* p = m0(A);
    int i;
    for(i = 0; i < 9; i++)
    {
        (i != 0) && ((i % 3) == 0) == 1? printf("\n") : printf("");
        // workaround for broken printf floating point support in nuttx
        printf("%c%2d.%03d ", *p < 0 ? '-' : ' ', (int)abs(*p), DECIMAL_PLACES3(*p));
        p = TRAVERSE_HORIZONTAL_3X3(p, m0(A));
    }
    printf("\n");
}