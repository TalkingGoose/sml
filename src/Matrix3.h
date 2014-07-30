#pragma once

/********************************************************************
    Author:         Paul Watkinson
    Website:        http://paulwatkinson.co.uk/

    Description:    A 3 dimensional matrix class
*********************************************************************/

#include "Vector3.h"
#include "Matrix2.h"
#include "Quaternion.h"

//////////////////////////////////////////////////////////////////////////
/// A matrix class for representing 3 dimensional matrices
//////////////////////////////////////////////////////////////////////////
template <typename T>
class Matrix3_t {
    public:
        union {
            struct { Vector3_t<T> column[3]; };
            struct { T linear[9]; };
        };

        /// Default Constructor
        Matrix3_t(void) { }

        /// Initialise this matrix with the given values
        Matrix3_t(const T v00, const T v01, const T v02, 
                  const T v10, const T v11, const T v12,
                  const T v20, const T v21, const T v22) {
            column[0] = Vector3_t<T>(v00, v01, v02);
            column[1] = Vector3_t<T>(v10, v11, v12);
            column[2] = Vector3_t<T>(v20, v21, v22);
        }

        /// Initialise this matrix with the given 3D vectors
        Matrix3_t(const Vector3_t<T> a, const Vector3_t<T> b, const Vector3_t<T> c) {
            column[0] = a; column[1] = b; column[2] = c;
        }

        /// Initialise this matrix with the given matrix
        Matrix3_t(const Matrix3_t & other) {
            memcpy((T *)&linear, other.linear, sizeof(T) * 9);
        }

        /// Initialise this matrix with the given array
        Matrix3_t(T * values) {
            memcpy((T *)&linear, values, sizeof(T) * 9);
        }

        /// Destructor
        virtual ~Matrix3_t(void) { }

        /// Set this matrix's values to the given values
        Matrix3_t & Set(const T v00, const T v01, const T v02, 
                        const T v10, const T v11, const T v12,
                        const T v20, const T v21, const T v22) {
            column[0] = Vector3_t<T>(v00, v01, v02);
            column[1] = Vector3_t<T>(v10, v11, v12);
            column[2] = Vector3_t<T>(v20, v21, v22);

            return (*this);
        }

        /// Set this matrix's values to the given 3D vectors
        Matrix3_t & Set(const Vector3_t<T> a, const Vector3_t<T> b, const Vector3_t<T> c) {
            column[0] = a; column[1] = b; column[2] = c;
            return (*this);
        }

        /// Set this matrix's values to the given array
        Matrix3_t & Set(T * values) {
            memcpy((T *)&linear, values, sizeof(T) * 9);
            return (*this);
        }

        /// Set this matrix's values to the given matrix's values
        Matrix3_t & Set(const Matrix3_t & other) {
            memcpy((T * )&linear, other.linear, sizeof(T) * 9);
            return (*this);
        }

        /// Set this matrix's values to the given matrix's values
        Matrix3_t & operator  = (const Matrix3_t & other) {
            memcpy((T * )&linear, other.linear, sizeof(T) * 9);
            return (*this);
        }

        /// Add all the values of the other matrix to this one and set this matrix to the result
        Matrix3_t & operator += (const Matrix3_t & other) {
            column[0] += other.column[0];
            column[1] += other.column[1];
            column[2] += other.column[2];

            return (*this);
        }

        /// Negate all the values of the other matrix with this one and set this matrix to the result
        Matrix3_t & operator -= (const Matrix3_t & other) {
            column[0] -= other.column[0];
            column[1] -= other.column[1];
            column[2] -= other.column[2];

            return (*this);
        }

        /// Multiply this matrix with a given type and set this matrix to the result
        template<typename U>
        Matrix3_t & operator *= (const U & value) {
            return (*this) = ((*this) * value);
        }

        /// Divide this matrix with a given type and set this matrix to the result
        template<typename U>
        Matrix3_t & operator /= (const U & value) {
            return (*this) = ((*this) / value);
        }

        /// Multiply this matrix by the inverse of another and set this matrix to the result
        Matrix3_t & operator /= (const Matrix3_t & other) {
            return ((*this) *= other.Inverse());
        }

        /// Increment all values by one
        Matrix3_t & operator ++ (void) {
            column[0]++;
            column[1]++;
            column[2]++;

            return (*this);
        }

        /// Decrement all values by one
        Matrix3_t & operator -- (void) {
            column[0]--;
            column[1]--;
            column[2]--;

            return (*this);
        }

        /// Add this matrix to another
        Matrix3_t operator + (const Matrix3_t & other) const {
            return Matrix3_t<T>(
                column[0] + other[0],
                column[1] + other[1],
                column[2] + other[2]
            );
        }

        /// Negate this matrix by another
        Matrix3_t operator - (const Matrix3_t & other) const {
            return Matrix3_t<T>(
                column[0] - other[0],
                column[1] - other[1],
                column[2] - other[2]
            );
        }

        /// Invert the sign of all values of this matrix
        Matrix3_t operator - (void) const {
            return Matrix3_t<T>(
                -column[0],
                -column[1],
                -column[2]
            );
        }

        /// Multiply this matrix with another
        Matrix3_t operator * (const Matrix3_t & other) const {
            return Matrix3_t<T>(
                (column[0][0] * other[0][0]) + (column[1][0] * other[0][1]) + (column[2][0] * other[0][2]),
                (column[0][1] * other[0][0]) + (column[1][1] * other[0][1]) + (column[2][1] * other[0][2]),
                (column[0][2] * other[0][0]) + (column[1][2] * other[0][1]) + (column[2][2] * other[0][2]),

                (column[0][0] * other[1][0]) + (column[1][0] * other[1][1]) + (column[2][0] * other[1][2]),
                (column[0][1] * other[1][0]) + (column[1][1] * other[1][1]) + (column[2][1] * other[1][2]),
                (column[0][2] * other[1][0]) + (column[1][2] * other[1][1]) + (column[2][2] * other[1][2]),

                (column[0][0] * other[2][0]) + (column[1][0] * other[2][1]) + (column[2][0] * other[2][2]),
                (column[0][1] * other[2][0]) + (column[1][1] * other[2][1]) + (column[2][1] * other[2][2]),
                (column[0][2] * other[2][0]) + (column[1][2] * other[2][1]) + (column[2][2] * other[2][2])
            );
        }

        /// Multiply this matrix with a given value
        Matrix3_t operator * (const T & n) const {
            return Matrix3_t<T>(
                column[0] * n,
                column[1] * n,
                column[2] * n
            );
        }

        /// Divide this matrix by a given value
        Matrix3_t operator / (const T & n) const {
            return Matrix3_t<T>(
                column[0] / n,
                column[1] / n,
                column[2] / n
            );
        }

        /// Check if this matrix is equal to another
        bool operator == (const Matrix3_t & other) {
            for (unsigned int i = 0; i < 9; ++i) {
                if (linear[i] != other[i]) { return false; }
            }

            return true;
        }

        /// Check if this matrix is not equal to another
        bool operator != (const Matrix3_t & other) {
            return !((*this) == other);
        }        

        /// Access a column within this matrix
        Vector3_t<T> & operator [] (unsigned int n) {
            if (n > 2) { throw "Out of range."; }
            return column[n];
        }

        /// Access a column within this matrix
        const Vector3_t<T> & operator [] (unsigned int n) const {
            if (n > 2) { throw "Out of range."; }
            return column[n];
        }

        /// Get the determinant of this matrix
        T Determinant(void) const {
            return ((linear[0] * ((linear[4] * linear[8]) - (linear[5] * linear[7]))) - 
                    (linear[1] * ((linear[8] * linear[3]) - (linear[5] * linear[6])))) + 
                    (linear[2] * ((linear[3] * linear[7]) - (linear[4] * linear[6])));
        }

        /// Get the inverse of this matrix
        Matrix3_t Inverse(void) const {
            T det = Determinant();

            if (det == 0) { return Matrix3_t<T>(); }

            return Matrix3_t<T>(
                 Matrix2_t<T>(linear[4], linear[7], linear[5], linear[8]).Determinant(),
                -Matrix2_t<T>(linear[1], linear[7], linear[2], linear[8]).Determinant(),
                 Matrix2_t<T>(linear[1], linear[4], linear[2], linear[5]).Determinant(),

                -Matrix2_t<T>(linear[3], linear[6], linear[5], linear[8]).Determinant(),
                 Matrix2_t<T>(linear[0], linear[6], linear[2], linear[8]).Determinant(),
                -Matrix2_t<T>(linear[0], linear[3], linear[2], linear[5]).Determinant(),

                 Matrix2_t<T>(linear[3], linear[6], linear[4], linear[7]).Determinant(),
                -Matrix2_t<T>(linear[0], linear[6], linear[1], linear[7]).Determinant(),
                 Matrix2_t<T>(linear[0], linear[3], linear[1], linear[4]).Determinant()
            ) * (static_cast<T>(1) / det);
        }

        /// Invert the matrix
        Matrix3_t & Invert(void) const {
            return (*this) = Inverse();
        }

        /// Get the transpose of this matrix
        Matrix3_t Transposed(void) const {
            return Matrix3_t<T>(
                linear[0], linear[3], linear[6],
                linear[1], linear[4], linear[7],
                linear[2], linear[5], linear[8]
            );
        }

        /// Transpose this matrix
        Matrix3_t & Transpose(void) {
            return (*this) = Transposed();
        }

        /// Get this matrix scaled by a given vector
        Matrix3_t Scaled(const Vector3_t<T> & vec) {
            return (*this) * Matrix3_t<T>(
                vec.x, 0.0f, 0.0f,
                0.0f, vec.y, 0.0f,
                0.0f, 0.0f, vec.z
            );
        }

        /// Scale this matrix by a given vector
        Matrix3_t & Scale(const Vector3_t<T> & vec) {
            return (*this) = Scaled(vec);
        }

        /// Get a zero initialised matrix
        static Matrix3_t Matrix3_t::Zero(void) {
            return Matrix3_t<T>();
        }

        /// Get an identity matrix
        static Matrix3_t Matrix3_t::Identity(void) {
            return Matrix3_t<T>(
                1.0f, 0.0f, 0.0f,
                0.0f, 1.0f, 0.0f,
                0.0f, 0.0f, 1.0f
            );
        }

        /// Multiply this matrix by a quaternion
        Matrix3_t operator * (const Quaternion & quat) const {
            T n = static_cast<T>(quat.SquaredMagnitude());
            T s = (n == 0 ? 0 : 2 / n);

            T xx = static_cast<T>(s * pow(quat.x, 2.0f)); T yy = static_cast<T>(s * pow(quat.y, 2.0f)); T zz = static_cast<T>(s * pow(quat.z, 2.0f));
            T xy = static_cast<T>(s * quat.x * quat.y);   T yz = static_cast<T>(s * quat.y * quat.z);   T xz = static_cast<T>(s * quat.x * quat.z);
            T wx = static_cast<T>(s * quat.w * quat.x);   T wy = static_cast<T>(s * quat.w * quat.y);   T wz = static_cast<T>(s * quat.w * quat.z);

            return (*this) * Matrix3_t<T>(
                1 - (yy + zz), xy + wz, xz - wy,
                xy - wz, 1 - (xx + zz), yz + wx,
                xz + wy, yz - wx, 1 - (xx + yy)
            );
        }
};

typedef Matrix3_t<float>  Matrix3f;
typedef Matrix3_t<double> Matrix3d;
typedef Matrix3_t<int>    Matrix3i;