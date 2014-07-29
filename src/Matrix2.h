#pragma once

/********************************************************************
    Author:         Paul Watkinson
    Website:        http://paulwatkinson.co.uk/

    Description:    A 2 dimensional matrix class
*********************************************************************/

#include "Vector2.h"

template <typename T>
class Matrix2_t {
    public:
        union {
            struct { Vector2_t<T> column[2]; };
            struct { T linear[4]; };
        };

        Matrix2_t(void) { }
        
        Matrix2_t(const T v00, const T v01, const T v10, const T v11) {
            column[0] = Vector2_t<T>(v00, v01);
            column[1] = Vector2_t<T>(v10, v11);
        }

        Matrix2_t(const Vector2_t<T> a, const Vector2_t<T> b) {
            column[0] = a; column[1] = b;
        }

        Matrix2_t(T * values) {
            memcpy((T *)&linear, values, sizeof(T) * 4);
        }

        virtual ~Matrix2_t(void) { }

        Matrix2_t & Set(const T v00, const T v01, const T v10, const T v11) {
            column[0] = Vector2_t<T>(v00, v01);
            column[1] = Vector2_t<T>(v10, v11);

            return (*this);
        }

        Matrix2_t & Set(const Vector2_t<T> a, const Vector2_t<T> b) {
            column[0] = a; column[1] = b;
            return (*this);
        }

        Matrix2_t & Set(T * values) {
            memcpy((T *)&linear, values, sizeof(T) * 4);
            return (*this);
        }

        Matrix2_t & Set(const Matrix2_t & other) {
            memcpy((T * )&linear, other.linear, sizeof(T) * 4);
            return (*this);
        }

        Matrix2_t & operator  = (const Matrix2_t & other) {
            memcpy((T * )&linear, other.linear, sizeof(T) * 4);
            return (*this);
        }

        Matrix2_t & operator += (const Matrix2_t & other) {
            column[0] += other.column[0];
            column[1] += other.column[1];

            return (*this);
        }

        Matrix2_t & operator -= (const Matrix2_t & other) {
            column[0] -= other.column[0];
            column[1] -= other.column[1];

            return (*this);
        }

        template <typename U>
        Matrix2_t & operator *= (const U & value) {
            return (*this) = ((*this) * value);
        }

        template <typename U>
        Matrix2_t & operator /= (const U & value) {
            return (*this) = ((*this) / value);
        }

        Matrix2_t & operator /= (const Matrix2_t & other) {
            return ((*this) *= other.Inverse());
        }

        Matrix2_t & operator ++ (void) {
            column[0]++;
            column[1]++;

            return (*this);
        }

        Matrix2_t & operator -- (void) {
            column[0]--;
            column[1]--;

            return (*this);
        }

        Matrix2_t operator + (const Matrix2_t & other) const {
            return Matrix2_t<T>(
                column[0] + other[0],
                column[1] + other[1]
            );
        }

        Matrix2_t operator - (const Matrix2_t & other) const {
            return Matrix2_t<T>(
                column[0] - other[0],
                column[1] - other[1]
            );
        }

        Matrix2_t operator - (void) const {
            return Matrix2_t<T>(
                -column[0],
                -column[1]
            );
        }

        Matrix2_t operator * (const Matrix2_t & other) const {
            return Matrix2_t<T>(
                (column[0][0] * other[0][0]) + (column[1][0] * other[0][1]),
                (column[0][1] * other[0][0]) + (column[1][1] * other[0][1]),
                (column[0][0] * other[1][0]) + (column[1][0] * other[1][1]),
                (column[0][1] * other[1][0]) + (column[1][1] * other[1][1])
            );
        }
        
        Matrix2_t operator * (const T & n) const {
            return Matrix2_t<T>(
                column[0] * n,
                column[1] * n
            );
        }

        Matrix2_t operator / (const T & n) const {
            return Matrix2_t<T>(
                column[0] / n,
                column[1] / n
            );
        }

        bool operator == (const Matrix2_t & other) {
            for (unsigned int i = 0; i < 4; ++i) {
                if (linear[i] != other[i]) { return false; }
            }

            return true;
        }

        bool operator != (const Matrix2_t & other) {
            return !((*this) == other);
        }

        Vector2_t & operator [] (unsigned int n) {
            if (n > 1) { throw "Out of range."; }
            return column[n];
        }

        const Vector2_t & operator [] (unsigned int n) const {
            if (n > 1) { throw "Out of range."; }
            return column[n];
        }

        T Determinant(void) const {
            return ((column[0][0] * column[1][1]) - (column[1][0] * column[0][1]));
        }

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

        Matrix2_t & Invert(void) const {
            return (*this) = Inverse();
        }

        Matrix2_t Transposed(void) const {
            return Matrix2_t<T>(
                linear[0], linear[2],
                linear[1], linear[3]
            );
        }

        Matrix2_t & Transpose(void) {
            return (*this) = Transposed();
        }

        Matrix2_t Scaled(const Vector2_t<T> & vec) {
            return (*this) * Matrix2_t<T>(
                vec.x, 0.0f,
                0.0f, vec.y
            );
        }

        Matrix2_t & Scale(const Vector2_t<T> & vec) {
            return (*this) = Scaled(vec);
        }

        static Matrix2_t Matrix2_t::Zero(void) {
            return Matrix2_t<T>();
        }

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