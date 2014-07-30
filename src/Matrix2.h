#pragma once

/********************************************************************
    Author:         Paul Watkinson
    Website:        http://paulwatkinson.co.uk/

    Description:    A 2 dimensional matrix class
*********************************************************************/

#include "Vector2.h"

//////////////////////////////////////////////////////////////////////////
/// A matrix class for representing 2 dimensional matrices
//////////////////////////////////////////////////////////////////////////
template <typename T>
class Matrix2_t {
    public:
        union {
            struct { Vector2_t<T> column[2]; };
            struct { T linear[4]; };
        };

        /// Default Constructor
        Matrix2_t(void) { }
        
        /// Initialise this matrix with the given values
        Matrix2_t(const T v00, const T v01, const T v10, const T v11) {
            column[0] = Vector2_t<T>(v00, v01);
            column[1] = Vector2_t<T>(v10, v11);
        }

        /// Initialise this matrix with the given 2D vectors
        Matrix2_t(const Vector2_t<T> a, const Vector2_t<T> b) {
            column[0] = a; column[1] = b;
        }

        /// Initialise this matrix with the given matrix
        Matrix2_t(const Matrix2_t & other) {
            memcpy((T *)&linear, other.linear, sizeof(T) * 4);
        }

        /// Initialise this matrix with the given array
        Matrix2_t(T * values) {
            memcpy((T *)&linear, values, sizeof(T) * 4);
        }

        /// Destructor
        virtual ~Matrix2_t(void) { }

        /// Set this matrix's values to the given values
        Matrix2_t & Set(const T v00, const T v01, const T v10, const T v11) {
            column[0] = Vector2_t<T>(v00, v01);
            column[1] = Vector2_t<T>(v10, v11);

            return (*this);
        }

        /// Set this matrix's values to the given 2D vectors
        Matrix2_t & Set(const Vector2_t<T> a, const Vector2_t<T> b) {
            column[0] = a; column[1] = b;
            return (*this);
        }

        /// Set this matrix's values to the given array's values
        Matrix2_t & Set(T * values) {
            memcpy((T *)&linear, values, sizeof(T) * 4);
            return (*this);
        }

        /// Set this matrix's values to the given matrix's values
        Matrix2_t & Set(const Matrix2_t & other) {
            memcpy((T * )&linear, other.linear, sizeof(T) * 4);
            return (*this);
        }

        /// Set this matrix's values to the given matrix's values
        Matrix2_t & operator  = (const Matrix2_t & other) {
            memcpy((T * )&linear, other.linear, sizeof(T) * 4);
            return (*this);
        }

        /// Add all the values of the other matrix to this one and set this matrix to the result
        Matrix2_t & operator += (const Matrix2_t & other) {
            column[0] += other.column[0];
            column[1] += other.column[1];

            return (*this);
        }

        /// Negate all the values of the other matrix with this one and set this matrix to the result
        Matrix2_t & operator -= (const Matrix2_t & other) {
            column[0] -= other.column[0];
            column[1] -= other.column[1];

            return (*this);
        }

        /// Multiply this matrix with a given type and set this matrix to the result
        template <typename U>
        Matrix2_t & operator *= (const U & value) {
            return (*this) = ((*this) * value);
        }

        /// Divide this matrix with a given type and set this matrix to the result
        template <typename U>
        Matrix2_t & operator /= (const U & value) {
            return (*this) = ((*this) / value);
        }

        /// Multiply this matrix by the inverse of another and set this matrix to the result
        Matrix2_t & operator /= (const Matrix2_t & other) {
            return ((*this) *= other.Inverse());
        }

        /// Increment all values by one
        Matrix2_t & operator ++ (void) {
            column[0]++;
            column[1]++;

            return (*this);
        }

        /// Decrement all values by one
        Matrix2_t & operator -- (void) {
            column[0]--;
            column[1]--;

            return (*this);
        }

        /// Add this matrix to another
        Matrix2_t operator + (const Matrix2_t & other) const {
            return Matrix2_t<T>(
                column[0] + other[0],
                column[1] + other[1]
            );
        }

        /// Negate this matrix by another
        Matrix2_t operator - (const Matrix2_t & other) const {
            return Matrix2_t<T>(
                column[0] - other[0],
                column[1] - other[1]
            );
        }

        /// Invert the sign of all values of this matrix
        Matrix2_t operator - (void) const {
            return Matrix2_t<T>(
                -column[0],
                -column[1]
            );
        }

        /// Multiply this matrix with another
        Matrix2_t operator * (const Matrix2_t & other) const {
            return Matrix2_t<T>(
                (column[0][0] * other[0][0]) + (column[1][0] * other[0][1]),
                (column[0][1] * other[0][0]) + (column[1][1] * other[0][1]),
                (column[0][0] * other[1][0]) + (column[1][0] * other[1][1]),
                (column[0][1] * other[1][0]) + (column[1][1] * other[1][1])
            );
        }
        
        /// Multiply this matrix with a given value
        Matrix2_t operator * (const T & n) const {
            return Matrix2_t<T>(
                column[0] * n,
                column[1] * n
            );
        }

        /// Divide this matrix by a given value
        Matrix2_t operator / (const T & n) const {
            return Matrix2_t<T>(
                column[0] / n,
                column[1] / n
            );
        }

        /// Check if this matrix is equal to another
        bool operator == (const Matrix2_t & other) {
            for (unsigned int i = 0; i < 4; ++i) {
                if (linear[i] != other[i]) { return false; }
            }

            return true;
        }

        /// Check if this matrix is not equal to another
        bool operator != (const Matrix2_t & other) {
            return !((*this) == other);
        }

        /// Access a column within this matrix
        Vector2_t & operator [] (unsigned int n) {
            if (n > 1) { throw "Out of range."; }
            return column[n];
        }

        /// Access a column within this matrix
        const Vector2_t & operator [] (unsigned int n) const {
            if (n > 1) { throw "Out of range."; }
            return column[n];
        }

        /// Get the determinant of this matrix
        T Determinant(void) const {
            return ((column[0][0] * column[1][1]) - (column[1][0] * column[0][1]));
        }

        /// Get the inverse of this matrix
        Matrix2_t Inverse(void) const {
            T det = Determinant();
            
            if (det == 0) { return Matrix2_t<T>(); }

            return Matrix2_t<T>(
                 column[1][1],
                -column[0][1],
                -column[1][0],
                 column[0][0]
            ) * (static_cast<T>(1) / det);
        }

        /// Invert the matrix
        Matrix2_t & Invert(void) const {
            return (*this) = Inverse();
        }

        /// Get the transpose of this matrix
        Matrix2_t Transposed(void) const {
            return Matrix2_t<T>(
                linear[0], linear[2],
                linear[1], linear[3]
            );
        }

        /// Transpose this matrix
        Matrix2_t & Transpose(void) {
            return (*this) = Transposed();
        }

        /// Get this matrix scaled by a given vector
        Matrix2_t Scaled(const Vector2_t<T> & vec) {
            return (*this) * Matrix2_t<T>(
                vec.x, 0.0f,
                0.0f, vec.y
            );
        }

        /// Scale this matrix by a given vector
        Matrix2_t & Scale(const Vector2_t<T> & vec) {
            return (*this) = Scaled(vec);
        }

        /// Get a zero initialised matrix
        static Matrix2_t Matrix2_t::Zero(void) {
            return Matrix2_t<T>();
        }

        /// Get an identity matrix
        static Matrix2_t Matrix2_t::Identity(void) {
            return Matrix2_t<T>(
                1.0f, 0.0f,
                0.0f, 1.0f
            );
        }
};

typedef Matrix2_t<float>  Matrix2f;
typedef Matrix2_t<double> Matrix2d;
typedef Matrix2_t<int>    Matrix2i;