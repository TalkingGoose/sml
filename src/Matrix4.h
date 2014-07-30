#pragma once

/********************************************************************
    Author:         Paul Watkinson
    Website:        http://paulwatkinson.co.uk/

    Description:    A 4 dimensional matrix class
*********************************************************************/

#include "Vector4.h"
#include "Matrix3.h"

//////////////////////////////////////////////////////////////////////////
/// A matrix class for representing 4 dimensional matrices
//////////////////////////////////////////////////////////////////////////
template <typename T>
class Matrix4_t {
    public:
        union {
            struct { Vector4_t<T> column[4]; };
            struct { T linear[16]; };
        };

        /// Default Constructor
        Matrix4_t(void) { }

        /// Initialise this matrix with the given values
        Matrix4_t(const T v00, const T v01, const T v02, const T v03,
                  const T v10, const T v11, const T v12, const T v13,
                  const T v20, const T v21, const T v22, const T v23,
                  const T v30, const T v31, const T v32, const T v33) {
            column[0] = Vector4_t<T>(v00, v01, v02, v03);
            column[1] = Vector4_t<T>(v10, v11, v12, v13);
            column[2] = Vector4_t<T>(v20, v21, v22, v23);
            column[3] = Vector4_t<T>(v30, v31, v32, v33);
        }

        /// Initialise this matrix with the given 4D vectors
        Matrix4_t(const Vector4_t<T> a, const Vector4_t<T> b, 
                  const Vector4_t<T> c, const Vector4_t<T> d) {
            column[0] = a; column[1] = b; column[2] = c; column[3] = d;
        }

        /// Initialise this matrix with the given matrix
        Matrix4_t(const Matrix4_t & other) {
            memcpy((T *)&linear, other.linear, sizeof(T) * 16);
        }

        /// Initialise this matrix with the given array
        Matrix4_t(T * values) {
            memcpy((T *)&linear, values, sizeof(T) * 16);
        }

        /// Destructor
        virtual ~Matrix4_t(void) { }

        /// Set this matrix's values to the given values
        Matrix4_t & Set(const T v00, const T v01, const T v02, const T v03,
                           const T v10, const T v11, const T v12, const T v13,
                           const T v20, const T v21, const T v22, const T v23,
                           const T v30, const T v31, const T v32, const T v33) {
            column[0] = Vector4_t<T>(v00, v01, v02, v03);
            column[1] = Vector4_t<T>(v10, v11, v12, v13);
            column[2] = Vector4_t<T>(v20, v21, v22, v23);
            column[3] = Vector4_t<T>(v30, v31, v32, v33);

            return (*this);
        }

        /// Set this matrix's values to the given 4D vectors
        Matrix4_t & Set(const Vector4_t<T> a, const Vector4_t<T> b, 
                        const Vector4_t<T> c, const Vector4_t<T> d) {
            column[0] = a; column[1] = b; column[2] = c; column[3] = d;
            return (*this);
        }

        /// Set this matrix's values to the given array
        Matrix4_t & Set(T * values) {
            memcpy((T *)&linear, values, sizeof(T) * 16);
            return (*this);
        }

        /// Set this matrix's values to the given matrix's values
        Matrix4_t & Set(const Matrix4_t & other) {
            memcpy((T * )&linear, other.linear, sizeof(T) * 16);
            return (*this);
        }

        /// Set this matrix's values to the given matrix's values
        Matrix4_t & operator  = (const Matrix4_t & other) {
            memcpy((T * )&linear, other.linear, sizeof(T) * 16);
            return (*this);
        }

        /// Add all the values of the other matrix to this one and set this matrix to the result
        Matrix4_t & operator += (const Matrix4_t & other) {
            column[0] += other.column[0];
            column[1] += other.column[1];
            column[2] += other.column[2];
            column[3] += other.column[3];

            return (*this);
        }

        /// Negate all the values of the other matrix with this one and set this matrix to the result
        Matrix4_t & operator -= (const Matrix4_t & other) {
            column[0] -= other.column[0];
            column[1] -= other.column[1];
            column[2] -= other.column[2];
            column[3] -= other.column[3];

            return (*this);
        }

        /// Multiply this matrix with a given type and set this matrix to the result
        template<typename U>
        Matrix4_t & operator *= (const U & value) {
            return (*this) = ((*this) * value);
        }

        /// Divide this matrix with a given type and set this matrix to the result
        template<typename U>
        Matrix4_t & operator /= (const U & value) {
            return (*this) = ((*this) / value);
        }

        /// Multiply this matrix by the inverse of another and set this matrix to the result
        Matrix4_t & operator /= (const Matrix4_t & other) {
            return ((*this) *= other.Inverse());
        }

        /// Increment all values by one
        Matrix4_t & operator ++ (void) {
            column[0]++;
            column[1]++;
            column[2]++;
            column[3]++;

            return (*this);
        }

        /// Decrement all values by one
        Matrix4_t & operator -- (void) {
            column[0]--;
            column[1]--;
            column[2]--;
            column[3]--;

            return (*this);
        }

        /// Add this matrix to another
        Matrix4_t operator + (const Matrix4_t & other) const {
            return Matrix4_t<T>(
                column[0] + other[0],
                column[1] + other[1],
                column[2] + other[2],
                column[3] + other[3]
            );
        }

        /// Negate this matrix by another
        Matrix4_t operator - (const Matrix4_t & other) const {
            return Matrix4_t<T>(
                column[0] - other[0],
                column[1] - other[1],
                column[2] - other[2],
                column[3] - other[3]
            );
        }

        /// Invert the sign of all values of this matrix
        Matrix4_t operator - (void) const {
            return Matrix4_t<T>(
                -column[0],
                -column[1],
                -column[2],
                -column[3]
            );
        }

        /// Multiply this matrix with another
        Matrix4_t operator * (const Matrix4_t & other) const {
            return Matrix4_t<T>(
                (column[0][0] * other[0][0]) + (column[1][0] * other[0][1]) + (column[2][0] * other[0][2]) + (column[3][0] * other[0][3]),
                (column[0][1] * other[0][0]) + (column[1][1] * other[0][1]) + (column[2][1] * other[0][2]) + (column[3][1] * other[0][3]),
                (column[0][2] * other[0][0]) + (column[1][2] * other[0][1]) + (column[2][2] * other[0][2]) + (column[3][2] * other[0][3]),
                (column[0][3] * other[0][0]) + (column[1][3] * other[0][1]) + (column[2][3] * other[0][2]) + (column[3][3] * other[0][3]),

                (column[0][0] * other[1][0]) + (column[1][0] * other[1][1]) + (column[2][0] * other[1][2]) + (column[3][0] * other[1][3]),
                (column[0][1] * other[1][0]) + (column[1][1] * other[1][1]) + (column[2][1] * other[1][2]) + (column[3][1] * other[1][3]),
                (column[0][2] * other[1][0]) + (column[1][2] * other[1][1]) + (column[2][2] * other[1][2]) + (column[3][2] * other[1][3]),
                (column[0][3] * other[1][0]) + (column[1][3] * other[1][1]) + (column[2][3] * other[1][2]) + (column[3][3] * other[1][3]),

                (column[0][0] * other[2][0]) + (column[1][0] * other[2][1]) + (column[2][0] * other[2][2]) + (column[3][0] * other[2][3]),
                (column[0][1] * other[2][0]) + (column[1][1] * other[2][1]) + (column[2][1] * other[2][2]) + (column[3][1] * other[2][3]),
                (column[0][2] * other[2][0]) + (column[1][2] * other[2][1]) + (column[2][2] * other[2][2]) + (column[3][2] * other[2][3]),
                (column[0][3] * other[2][0]) + (column[1][3] * other[2][1]) + (column[2][3] * other[2][2]) + (column[3][3] * other[2][3]),

                (column[0][0] * other[3][0]) + (column[1][0] * other[3][1]) + (column[2][0] * other[3][2]) + (column[3][0] * other[3][3]),
                (column[0][1] * other[3][0]) + (column[1][1] * other[3][1]) + (column[2][1] * other[3][2]) + (column[3][1] * other[3][3]),
                (column[0][2] * other[3][0]) + (column[1][2] * other[3][1]) + (column[2][2] * other[3][2]) + (column[3][2] * other[3][3]),
                (column[0][3] * other[3][0]) + (column[1][3] * other[3][1]) + (column[2][3] * other[3][2]) + (column[3][3] * other[3][3])
            );
        }


        /// Multiply this matrix with a given value
        Matrix4_t operator * (const T & n) const {
            return Matrix4_t<T>(
                column[0] * n,
                column[1] * n,
                column[2] * n,
                column[3] * n
            );
        }

        /// Divide this matrix by a given value
        Matrix4_t operator / (const T & n) const {
            return Matrix4_t<T>(
                column[0] / n,
                column[1] / n,
                column[2] / n,
                column[3] / n
            );
        }

        /// Check if this matrix is equal to another
        bool operator == (const Matrix4_t & other) {
            for (unsigned int i = 0; i < 16; ++i) {
                if (linear[i] != other[i]) { return false; }
            }

            return true;
        }

        /// Check if this matrix is not equal to another
        bool operator != (const Matrix4_t & other) {
            return !((*this) == other);
        }

        /// Access a column within this matrix
        Vector4_t<T> & operator [] (unsigned int n) {
            if (n > 3) { throw "Out of range."; }
            return column[n];
        }

        /// Access a column within this matrix
        const Vector4_t<T> & operator [] (unsigned int n) const {
            if (n > 3) { throw "Out of range."; }
            return column[n];
        }

        /// Get the determinant of this matrix
        // http://www.researchgate.net/publication/257068341_New_method_to_compute_the_determinant_of_a_4x4_matrix/file/9c9605244ab33c5a9f.pdf
        T Determinant(void) const {
            return (column[0][0] * column[1][1] * column[2][2] * column[3][3]) +
                   (column[0][0] * column[3][1] * column[1][2] * column[2][3]) +
                   (column[0][0] * column[2][1] * column[3][2] * column[1][3]) -
                   (column[0][0] * column[3][1] * column[2][2] * column[1][3]) -
                   (column[0][0] * column[1][1] * column[3][2] * column[2][3]) -
                   (column[0][0] * column[2][1] * column[1][2] * column[3][3]) -
                   
                   (column[1][0] * column[0][1] * column[2][2] * column[3][3]) -
                   (column[1][0] * column[2][1] * column[3][2] * column[0][3]) -
                   (column[1][0] * column[3][1] * column[0][2] * column[2][3]) +
                   (column[1][0] * column[3][1] * column[2][2] * column[0][3]) +
                   (column[1][0] * column[0][1] * column[3][2] * column[2][3]) +
                   (column[1][0] * column[2][1] * column[0][2] * column[3][3]) +
                   
                   (column[2][0] * column[0][1] * column[1][2] * column[3][3]) +
                   (column[2][0] * column[1][1] * column[3][2] * column[0][3]) +
                   (column[2][0] * column[3][1] * column[0][2] * column[1][3]) -
                   (column[2][0] * column[3][1] * column[1][2] * column[0][3]) -
                   (column[2][0] * column[0][1] * column[3][2] * column[1][3]) -
                   (column[2][0] * column[1][1] * column[0][2] * column[3][3]) -
                   
                   (column[3][0] * column[0][1] * column[1][2] * column[2][3]) -
                   (column[3][0] * column[1][1] * column[2][2] * column[0][3]) -
                   (column[3][0] * column[2][1] * column[0][2] * column[1][3]) +
                   (column[3][0] * column[2][1] * column[1][2] * column[0][3]) +
                   (column[3][0] * column[0][1] * column[2][2] * column[1][3]) +
                   (column[3][0] * column[1][1] * column[0][2] * column[2][3]);
        }

        /// Get the inverse of this matrix
        Matrix4_t Inverse(void) const {
            T det = Determinant();

            if (det == 0) { return Matrix4_t<T>(); }

            return Matrix4_t<T>(
                 Matrix3_t<T>(column[1][1], column[1][2], column[1][3], column[2][1], column[2][2], column[2][3], column[3][1], column[3][2], column[3][3]).Determinant(),
                -Matrix3_t<T>(column[0][1], column[0][2], column[0][3], column[2][1], column[2][2], column[2][3], column[3][1], column[3][2], column[3][3]).Determinant(),
                 Matrix3_t<T>(column[0][1], column[0][2], column[0][3], column[1][1], column[1][2], column[1][3], column[3][1], column[3][2], column[3][3]).Determinant(),
                -Matrix3_t<T>(column[0][1], column[0][2], column[0][3], column[1][1], column[1][2], column[1][3], column[2][1], column[2][2], column[2][3]).Determinant(),

                -Matrix3_t<T>(column[1][0], column[1][2], column[1][3], column[2][0], column[2][2], column[2][3], column[3][0], column[3][2], column[3][3]).Determinant(),
                 Matrix3_t<T>(column[0][0], column[0][2], column[0][3], column[2][0], column[2][2], column[2][3], column[3][0], column[3][2], column[3][3]).Determinant(),
                -Matrix3_t<T>(column[0][0], column[0][2], column[0][3], column[1][0], column[1][2], column[1][3], column[3][0], column[3][2], column[3][3]).Determinant(),
                 Matrix3_t<T>(column[0][0], column[0][2], column[0][3], column[1][0], column[1][2], column[1][3], column[2][0], column[2][2], column[2][3]).Determinant(),


                 Matrix3_t<T>(column[1][0], column[1][1], column[1][3], column[2][0], column[2][1], column[2][3], column[3][0], column[3][1], column[3][3]).Determinant(),
                -Matrix3_t<T>(column[0][0], column[0][1], column[0][3], column[2][0], column[2][1], column[2][3], column[3][0], column[3][1], column[3][3]).Determinant(),
                 Matrix3_t<T>(column[0][0], column[0][1], column[0][3], column[1][0], column[1][1], column[1][3], column[3][0], column[3][1], column[3][3]).Determinant(),
                -Matrix3_t<T>(column[0][0], column[0][1], column[0][3], column[1][0], column[1][1], column[1][3], column[2][0], column[2][1], column[2][3]).Determinant(),

                -Matrix3_t<T>(column[1][0], column[1][1], column[1][2], column[2][0], column[2][1], column[2][2], column[3][0], column[3][1], column[3][2]).Determinant(),
                 Matrix3_t<T>(column[0][0], column[0][1], column[0][2], column[2][0], column[2][1], column[2][2], column[3][0], column[3][1], column[3][2]).Determinant(),
                -Matrix3_t<T>(column[0][0], column[0][1], column[0][2], column[1][0], column[1][1], column[1][2], column[3][0], column[3][1], column[3][2]).Determinant(),
                 Matrix3_t<T>(column[0][0], column[0][1], column[0][2], column[1][0], column[1][1], column[1][2], column[2][0], column[2][1], column[2][2]).Determinant()
            ) * (static_cast<T>(1) / det);
        }

        /// Invert the matrix
        Matrix4_t & Invert(void) const {
            return (*this) = Inverse();
        }

        /// Get the transpose of this matrix
        Matrix4_t Transposed(void) const {
            return Matrix4_t<T>(
                linear[0],  linear[4],  linear[8],  linear[12],
                linear[1],  linear[5],  linear[9],  linear[13],
                linear[2],  linear[6],  linear[10], linear[14],
                linear[3],  linear[7],  linear[11], linear[15]
            );
        }

        /// Transpose this matrix
        Matrix4_t & Transpose(void) {
            return (*this) = Transposed();
        }
        
        // Get this matrix translated by the given vector
        Matrix4_t Translated(const Vector3_t<T> & vec) {
            return (*this) * Matrix4_t<T>(
                1.0f, 0.0f, 0.0f, 0.0f,
                0.0f, 1.0f, 0.0f, 0.0f,
                0.0f, 0.0f, 1.0f, 0.0f,
                vec.x, vec.y, vec.z, 1.0f
            );
        }

        /// Translate this matrix
        Matrix4_t & Translate(const Vector3_t<T> & vec) {
            return (*this) = Translated(vec);
        }

        /// Get this matrix scaled by a given vector
        Matrix4_t Scaled(const Vector3_t<T> & vec) {
            return (*this) * Matrix4_t<T>(
                vec.x, 0.0f, 0.0f, 0.0f,
                0.0f, vec.y, 0.0f, 0.0f,
                0.0f, 0.0f, vec.z, 0.0f,
                0.0f, 0.0f, 0.0f, 1.0f
            );
        }

        /// Scale this matrix by a given vector
        Matrix4_t & Scale(const Vector3_t<T> & vec) {
            return (*this) = Scaled(vec);
        }

        /// Get a zero initialised matrix
        static Matrix4_t Matrix4_t::Zero(void) {
            return Matrix4_t<T>();
        }

        /// Get an identity matrix
        static Matrix4_t Matrix4_t::Identity(void) {
            return Matrix4_t<T>(
                1.0f, 0.0f, 0.0f, 0.0f,
                0.0f, 1.0f, 0.0f, 0.0f,
                0.0f, 0.0f, 1.0f, 0.0f,
                0.0f, 0.0f, 0.0f, 1.0f
            );
        }

        /// Multiply this matrix by a quaternion
        Matrix4_t operator * (const Quaternion & quat) const {
            T n = static_cast<T>(quat.SquaredMagnitude());
            T s = (n == 0 ? 0 : 2 / n);

            T xx = static_cast<T>(s * pow(quat.x, 2.0f)); T yy = static_cast<T>(s * pow(quat.y, 2.0f)); T zz = static_cast<T>(s * pow(quat.z, 2.0f));
            T xy = static_cast<T>(s * quat.x * quat.y);   T yz = static_cast<T>(s * quat.y * quat.z);   T xz = static_cast<T>(s * quat.x * quat.z);
            T wx = static_cast<T>(s * quat.w * quat.x);   T wy = static_cast<T>(s * quat.w * quat.y);   T wz = static_cast<T>(s * quat.w * quat.z);

            return (*this) * Matrix4_t<T>(
                1 - (yy + zz), xy + wz, xz - wy, 0.0f,
                xy - wz, 1 - (xx + zz), yz + wx, 0.0f,
                xz + wy, yz - wx, 1 - (xx + yy), 0.0f,
                0.0f, 0.0f, 0.0f, 1.0f
                );
        }
};

typedef Matrix4_t<float>  Matrix4f;
typedef Matrix4_t<double> Matrix4d;
typedef Matrix4_t<int>    Matrix4i;