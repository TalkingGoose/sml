#pragma once

/********************************************************************
    Author:         Paul Watkinson
    Website:        http://paulwatkinson.co.uk/

    Description:    A 4 dimensional vector class
*********************************************************************/

#include "Vector3.h"

/// A vector class for representing 4 dimensional vectors
template <typename T>
struct Vector4_t {
    public:
        T x, y, z, w;

        Vector4_t(void) {
            x = y = z = w = static_cast<T>(0);
        }
        
        /// <summary>Construct a 4 dimensional vector given the X, Y, Z and W values.</summary>
        /// <param name="a">The X value to be set.</param>
        /// <param name="b">The Y value to be set.</param>
        /// <param name="c">The Z value to be set.</param>
        /// <param name="d">The W value to be set.</param>
        Vector4_t(T a, T b, T c, T d) {
            x = a; y = b; z = c; w = d;
        }

        /// <summary>Construct a 4 dimensional vector given a 3 dimensional vector and a W value.</summary>
        /// <param name="other">The 3 dimensional vector.</param>
        /// <param name="d">The W value to be set.</param>
        Vector4_t(const Vector3_t<T> & other, T d) {
            this->x = other.x;
            this->y = other.y;
            this->z = other.z;
            this->w = d;
        }

        /// <summary>Set the X, Y, Z and W values of the 4 dimensional vector.</summary>
        /// <param name="a">The X value to be set.</param>
        /// <param name="b">The Y value to be set.</param>
        /// <param name="c">The Z value to be set.</param>
        /// <param name="d">The W value to be set.</param>
        Vector4_t & Set(T a, T b, T c, T d) {
            this->x = a;
            this->y = b;
            this->z = c;
            this->w = d;

            return (*this);
        }

        /// <summary>Set the  X, Y, Z and W values of the 4 dimensional vector using given a 3 dimensional 
        /// vector and a W value.</summary>
        /// <param name="other">The 3 dimensional vector of which the X, Y and Z value will be set.</param>
        /// <param name="d">The W value to be set.</param>
        Vector4_t & Set(const Vector3_t<T> other, T d = 0.0l) {
            this->x = other.x;
            this->y = other.y;
            this->z = other.z;
            this->w = d;

            return (*this);
        }

        /// <summary>Set the X, Y, Z and W values of the 4 dimensional vector to that of anothers.</summary>
        /// <param name="other">The 4 dimensional vector to be set the X, Y, Z and W values with.</param>
        Vector4_t & Set(const Vector4_t & other) {
            this->x = other.x;
            this->y = other.y;
            this->z = other.z;
            this->w = other.w;

            return (*this);
        }

        Vector4_t & operator = (const Vector3_t<T> & other) {
            this->x = other.x;
            this->y = other.y;
            this->z = other.z;
            this->w = 0.0f;

            return (*this);
        }

        Vector4_t & operator = (const Vector4_t & other) {
            this->x = other.x;
            this->y = other.y;
            this->z = other.z;
            this->w = other.w;

            return (*this);
        }
        
        Vector4_t & operator += (const Vector3_t<T> & other) {
            this->x += other.x;
            this->y += other.y;
            this->z += other.z;

            return (*this);
        }

        Vector4_t & operator += (const Vector4_t & other) {
            this->x += other.x;
            this->y += other.y;
            this->z += other.z;
            this->w += other.w;

            return (*this);
        }

        Vector4_t & operator -= (const Vector3_t<T> & other) {
            this->x -= other.x;
            this->y -= other.y;
            this->w -= other.z;

            return (*this);
        }

        Vector4_t & operator -= (const Vector4_t & other) {
            this->x -= other.x;
            this->y -= other.y;
            this->z -= other.z;
            this->w -= other.w;

            return (*this);
        }

        Vector4_t & operator *= (const Vector3_t<T> & other) {
            this->x *= other.x;
            this->y *= other.y;
            this->z *= other.z;

            return (*this);
        }

        Vector4_t & operator *= (const Vector4_t & other) {
            this->x *= other.x;
            this->y *= other.y;
            this->z *= other.z;
            this->w *= other.w;

            return (*this);
        }

        Vector4_t & operator *= (T n) {
            this->x *= n;
            this->y *= n;
            this->z *= n;
            this->w *= n;

            return (*this);
        }

        Vector4_t & operator /= (const Vector3_t<T> & other) {
            if (other.x != 0) {
                this->x /= other.x;
            }

            if (other.y != 0) {
                this->y /= other.y;
            }

            if (other.z != 0) {
                this->z /= other.z;
            }

            return (*this);
        }

        Vector4_t & operator /= (const Vector4_t & other) {
            if (other.x != 0) {
                x /= other.x;
            }

            if (other.y != 0) {
                y /= other.y;
            }

            if (other.z != 0) {
                z /= other.z;
            }

            if (other.w != 0) {
                w /= other.w;
            }

            return (*this);
        }

        Vector4_t & operator /= (T n) {
            x /= n;
            y /= n;
            z /= n;

            return (*this);
        }

        Vector4_t & operator ++ (void) {
            x++; y++; z++;
            return (*this);
        }

        Vector4_t & operator -- (void) {
            x--; y--; z--;
            return (*this);
        }

        Vector4_t operator + (const Vector3_t<T> & other) const {
            return Vector4_t(this->x + other.x, this->y + other.y, this->z + other.z, this->w);
        }

        Vector4_t operator + (const Vector4_t & other) const {
            return Vector4_t(this->x + other.x, this->y + other.y, this->z + other.z, this->w + other.w);
        }

        Vector4_t operator - (const Vector3_t<T> & other) const {
            return Vector4_t(this->x - other.x, this->y - other.y, this->z - other.z, this->w);
        }

        Vector4_t operator - (const Vector4_t & other) const {
            return Vector4_t(this->x - other.x, this->y - other.y, this->z - other.z, this->w - other.w);
        }

        Vector4_t operator - (void) const {
            return Vector4_t(-this->x, -this->y, -this->z, -this->w);
        }

        Vector4_t operator * (const Vector3_t<T> & other) const {
            return Vector4_t(this->x * other.x, this->y * other.y, this->z * other.z, this->w);
        }

        Vector4_t operator * (const Vector4_t & other) const {
            return Vector4_t(this->x * other.x, this->y * other.y, this->z * other.z, this->w * other.w);
        }

        Vector4_t operator * (T n) const {
            return Vector4_t(this->x * n, this->y * n, this->z * n, this->w * n);
        }

        Vector4_t operator / (const Vector3_t<T> & other) const {
            return Vector4_t(this->x / other.x, this->y / other.y, this->z / other.z, this->w);
        }

        Vector4_t operator / (const Vector4_t & other) const {
            return Vector4_t(this->x / other.x, this->y / other.y, this->z / other.z, this->w / other.w);
        }

        Vector4_t operator / (T n) const {
            return Vector4_t(this->x / n, this->y / n, this->z / n, this->w / n);
        }

        bool operator ! (void) const {
            return (this->x != 0 && this->y != 0 && this->z != 0 && this->w != 0);
        }

        bool operator == (const Vector4_t & other) const {
            return (this->x == other.x && this->y == other.y && this->z == other.z && this->w == other.w);
        }

        bool operator != (const Vector4_t & other) const {
            return (this->x != other.x || this->y != other.y || this->z != other.z || this->w != other.w);
        }

        T & operator [] (unsigned int n) {
            if (n > 3) { throw "Out of bounds."; }
            return *(&x + n);
        }

        const T & operator [] (const unsigned int n) const {
            if (n > 3) { throw "Out of bounds."; }
            return *(&x + n);
        }

        T SquaredMagnitude(void) const {
            return pow(this->x, 2.0f) + pow(this->y, 2.0f) + pow(this->z, 2.0f) + (this->w, 2.0f);
        }

        T Magnitude(void) const {
            return sqrt(this->SquaredMagnitude());
        }

        Vector4_t & Normalise(void) {
            if (!(*this)) {
                return (*this);
            }

            return (*this /= Magnitude());
        }

        Vector4_t Normalised(void) const {
            T m = Magnitude();
            return Vector4_t(this->x / m, this->y / m, this->z / m, this->w / m);
        }

        /// <summary>Returns the sum of the X, Y, Z and W components.</summary>
        T Sum(void) const {
            return x + y + z + w;
        }

        Vector4_t DotProduct(const Vector4_t & other) const {
            return ((*this) + other).Sum();
        }

        Vector4_t CrossProduct(const Vector4_t & other) const {
            return Vector4_t(
                (this->y * other.z) - (this->z * other.y),
                (this->z * other.x) - (this->x * other.z),
                (this->x * other.y) - (this->y * other.x),
                static_cast<T>(0.0f)
            );
        }
};

template <typename T>
inline Vector4_t<T> operator * (T n, const Vector4_t<T> & v) {
    return v * n;
}

typedef Vector4_t<float>  Vector4f;
typedef Vector4_t<double> Vector4d;
typedef Vector4_t<int>    Vector4i;