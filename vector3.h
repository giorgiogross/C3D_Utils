/****************************************************************************
 * examples/acceleration_calibration/vector.h
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

#ifndef INCLUDE_VECTOR3_H
#define INCLUDE_VECTOR3_H

/****************************************************************************
 * Public Types
 ****************************************************************************/

typedef struct __attribute__((packed))
{
  float x, y, z;
} vector3_t;

typedef struct __attribute__((packed))
{
    vector3_t v1, v2, v3;
} matrix3x3_t;

/****************************************************************************
 * Public Data
 ****************************************************************************/

extern vector3_t X, Y, Z;

/****************************************************************************
 * Public Function Prototypes
 ****************************************************************************/

/****************************************************************************
 * Name: enclosed_angle
 *
 * Description:
 *   Calculate the enclosed angle of two 3d vectors.
 *
 * Input Parameters:
 *   U - First vector
 *   V - Second vector
 *
 * Returned Value:
 *   Angle in radians
 *
 ****************************************************************************/

float enclosed_angle (vector3_t U, vector3_t V);

/****************************************************************************
 * Name: vector_product
 *
 * Description:
 *   Calculate the dot product of two 3d vectors.
 *
 * Input Parameters:
 *   U - First vector
 *   V - Second vector
 *
 * Returned Value:
 *   The dot product
 *
 ****************************************************************************/

float vector3_product (vector3_t U, vector3_t V);

/****************************************************************************
 * Name: vector_sum
 *
 * Description:
 *   Add two 3d vectors.
 *
 * Input Parameters:
 *   T - Pointer to target vector
 *   U - First vector
 *   V - Second vector
 *
 * Returned Value:
 *   void
 *
 ****************************************************************************/

void vector3_sum (vector3_t* T, vector3_t U, vector3_t V);

/****************************************************************************
 * Name: vector_length
 *
 * Description:
 *   Calculate the length of a 3d vector.
 *
 * Input Parameters:
 *   U - Vector
 *
 * Returned Value:
 *   Length of U
 *
 ****************************************************************************/

float vector3_length (vector3_t U);

/****************************************************************************
 * Name: vector_normalize
 *
 * Description:
 *   Normalize U to a length of 1.
 *
 * Input Parameters:
 *   U - Pointer to vector to be normalized
 *
 * Returned Value:
 *   void
 *
 ****************************************************************************/

void vector3_normalize (vector3_t* U);

/****************************************************************************
 * Name: multiply3x3_matrix_matrix
 *
 * Description:
 *   Multiply two 3x3 matrices. 
 *   Vectors of matrices are assumed to be their columns.
 *
 * Input Parameters:
 *   T - Pointer to target matrix
 *   A - First matrix
 *   B - Second matrix
 *
 * Returned Value:
 *   Angle in radians
 *
 ****************************************************************************/

void multiply3x3_matrix_matrix (matrix3x3_t* T, matrix3x3_t A, matrix3x3_t B);

/****************************************************************************
 * Name: multiply3x3_matrix_vector
 *
 * Description:
 *   Multiply a 3x3 matrix with a 3d vector. 
 *   Vectors of the matrix are assumed to be its columns.
 *
 * Input Parameters:
 *   T - Pointer to target vector
 *   A - 3x3 matrix
 *   V - 3d vector
 *
 * Returned Value:
 *   void
 *
 ****************************************************************************/

void multiply3x3_matrix_vector (vector3_t* T, matrix3x3_t A, vector3_t V);

/****************************************************************************
 * Name: multiply_vector_vector
 *
 * Description:
 *   Multiply two 3d vectors.
 *
 * Input Parameters:
 *   U - First vector
 *   V - Second vector
 *
 * Returned Value:
 *   The product
 *
 ****************************************************************************/

float multiply_vector3_vector3 (vector3_t U, vector3_t V);

/****************************************************************************
 * Name: multiply_vector_scalar
 *
 * Description:
 *   Multiply a 3d vector with a scalar
 *
 * Input Parameters:
 *   T - Pointer to target vector
 *   V - Vector to be multiplied
 *   s - Scalar
 *
 * Returned Value:
 *   void
 *
 ****************************************************************************/

void multiply_vector3_scalar (vector3_t *T, vector3_t U, float s);

/****************************************************************************
 * Name: orthonormalize3x3
 *
 * Description:
 *   Transform a 3x3 matrix into a orthonormalized matrix by using the 
 *   Gram-Schmidt algorithm
 *
 * Input Parameters:
 *   A - Pointer to matrix to be transformed
 *
 * Returned Value:
 *   void
 *
 ****************************************************************************/

void orthonormalize3x3(matrix3x3_t* A);

/****************************************************************************
 * Name: transpose3x3
 *
 * Description:
 *   Transpose a 3x3 matrix.
 *    Vectors of the matrix are assumed to be its columns.
 *
 * Input Parameters:
 *   A - Pointer to matrix to be transposed
 *
 * Returned Value:
 *   Angle in radians
 *
 ****************************************************************************/

void transpose3x3(matrix3x3_t* A);

/****************************************************************************
 * Name: det3x3
 *
 * Description:
 *   Calculate the determinant of a 3x3 matrix with Sarrus.
 *   Vectors of the matrix are assumed to be its columns.
 *
 * Input Parameters:
 *   A - The matrix to calculate the determinant of
 *
 * Returned Value:
 *   The determinant
 *
 ****************************************************************************/

float det3x3(matrix3x3_t A);

/****************************************************************************
 * Name: print3x3_matrix
 *
 * Description:
 *   Output a 3x3 matrix.
 *   Vectors of the matrix are assumed to be its columns.
 *
 * Input Parameters:
 *   A - Pointer to matrix to be printed
 *
 * Returned Value:
 *   void
 *
 ****************************************************************************/

void print3x3_matrix(matrix3x3_t A);

#endif /* INCLUDE_VECTOR3_H */