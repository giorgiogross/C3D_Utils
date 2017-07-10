/****************************************************************************
    MIT License

    Copyright (c) 2017 Giorgio Groß

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
    
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
} matrice3x3_t;

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
 * Name: multiply3x3_matrice_matrice
 *
 * Description:
 *   Multiply two 3x3 matrices. 
 *   Vectors of matrices are assumed to be their columns.
 *
 * Input Parameters:
 *   T - Pointer to target matrice
 *   A - First matrice
 *   B - Second matrice
 *
 * Returned Value:
 *   Angle in radians
 *
 ****************************************************************************/

void multiply3x3_matrice_matrice (matrice3x3_t* T, matrice3x3_t A, matrice3x3_t B);

/****************************************************************************
 * Name: multiply3x3_matrice_vector
 *
 * Description:
 *   Multiply a 3x3 matrice with a 3d vector. 
 *   Vectors of the matrice are assumed to be its columns.
 *
 * Input Parameters:
 *   T - Pointer to target vector
 *   A - 3x3 matrice
 *   V - 3d vector
 *
 * Returned Value:
 *   void
 *
 ****************************************************************************/

void multiply3x3_matrice_vector (vector3_t* T, matrice3x3_t A, vector3_t V);

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
 *   Transform a 3x3 matrice into a orthonormalized matrice by using the 
 *   Gram-Schmidt algorithm
 *
 * Input Parameters:
 *   A - Pointer to matrice to be transformed
 *
 * Returned Value:
 *   void
 *
 ****************************************************************************/

void orthonormalize3x3(matrice3x3_t* A);

/****************************************************************************
 * Name: transpose3x3
 *
 * Description:
 *   Transpose a 3x3 matrice.
 *    Vectors of the matrice are assumed to be its columns.
 *
 * Input Parameters:
 *   A - Pointer to matrice to be transposed
 *
 * Returned Value:
 *   Angle in radians
 *
 ****************************************************************************/

void transpose3x3(matrice3x3_t* A);

/****************************************************************************
 * Name: det3x3
 *
 * Description:
 *   Calculate the determinant of a 3x3 matrice with Sarrus.
 *   Vectors of the matrice are assumed to be its columns.
 *
 * Input Parameters:
 *   A - The matrice to calculate the determinant of
 *
 * Returned Value:
 *   The determinant
 *
 ****************************************************************************/

float det3x3(matrice3x3_t A);

/****************************************************************************
 * Name: print3x3_matrice
 *
 * Description:
 *   Output a 3x3 matrice.
 *   Vectors of the matrice are assumed to be its columns.
 *
 * Input Parameters:
 *   A - Pointer to matrice to be printed
 *
 * Returned Value:
 *   void
 *
 ****************************************************************************/

void print3x3_matrice(matrice3x3_t A);

#endif /* INCLUDE_VECTOR3_H */