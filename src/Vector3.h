#pragma once

/********************************************************************
    Author:         Paul Watkinson
    Website:        http://paulwatkinson.co.uk/

    Description:    A 3 dimensional vector class
*********************************************************************/

#include <stdexcept>
#include "Vector2.h"

/// A vector class for representing 3 dimensional vectors
template<typename T>
struct Vector3_t {
    public:
        T x, y, z;

        Vector3_t(void) {
            x = y = z = static_cast<T>(0);
        }
        
        /// <summary>Construct a 3 dimensional vector given the X, Y and Z values.</summary>
        /// <param name="a">The X value to be set.</param>
        /// <param name="b">The Y value to be set.</param>
        /// <param name="c">The Z value to be set.</param>
        Vector3_t(T a, T b, T c) {
            x = a; y = b; z = c;
        }

        /// <summary>Construct a 3 dimensional vector given a 2 dimensional vector and a Z value.</summary>
        /// <param name="other">The 2 dimensional vector.</param>
        /// <param name="c">The Z value to be set.</param>
        Vector3_t(const Vector2_t<T> & other, T c) {
            this->x = other.x;
            this->y = other.y;
            this->z = c;
        }

        /// <summary>Set the X, Y and Z values of the 3 dimensional vector.</summary>
        /// <param name="a">The X value to be set.</param>
        /// <param name="b">The Y value to be set.</param>
        /// <param name="c">The Z value to be set.</param>
        Vector3_t & Set(T a, T b, T c) {
            this->x = a;
            this->y = b;
            this->z = c;

            return (*this);
        }

        /// <summary>Set the X, Y and Z values of the 3 dimensional vector using given a 2 dimensional 
        /// vector and a Z value.</summary>
        /// <param name="other">The 2 dimensional vector of which the X and Y value will be set.</param>
        /// <param name="c">The Z value to be set.</param>
        Vector3_t & Set(const Vector2_t<T> other, T c = 0.0f) {
            this->x = other.x;
            this->y = other.y;
            this->z = c;

            return (*this);
        }

        /// <summary>Set the X, Y and Z values of the 3 dimensional vector to that of anothers.</summary>
        /// <param name="other">The 3 dimensional vector to be set the X, Y and Z values with.</param>
        Vector3_t & Set(const Vector3_t & other) {
            this->x = other.x;
            this->y = other.y;
            this->z = other.z;

            return (*this);
        }

        Vector3_t & operator = (const Vector2_t<T> & other) {
            this->x = other.x;
            this->y = other.y;
            this->z = 0.0f;

            return (*this);
        }

        Vector3_t & operator = (const Vector3_t & other) {
            this->x = other.x;
            this->y = other.y;
            this->z = other.z;

            return (*this);
        }
        
        Vector3_t & operator += (const Vector2_t<T> & other) {
            this->x += other.x;
            this->y += other.y;

            return (*this);
        }

        Vector3_t & operator += (const Vector3_t & other) {
            this->x += other.x;
            this->y += other.y;
            this->z += other.z;

            return (*this);
        }

        Vector3_t & operator -= (const Vector2_t<T> & other) {
            this->x -= other.x;
            this->y -= other.y;

            return (*this);
        }

        Vector3_t & operator -= (const Vector3_t & other) {
            this->x -= other.x;
            this->y -= other.y;
            this->z -= other.z;

            return (*this);
        }

        Vector3_t & operator *= (const Vector2_t<T> & other) {
            this->x *= other.x;
            this->y *= other.y;

            return (*this);
        }

        Vector3_t & operator *= (const Vector3_t & other) {
            this->x *= other.x;
            this->y *= other.y;
            this->z *= other.z;

            return (*this);
        }

        Vector3_t & operator *= (T n) {
            this->x *= n;
            this->y *= n;
            this->z *= n;

            return (*this);
        }

        Vector3_t & operator /= (const Vector2_t<T> & other) {
            if (other.x != 0) {
                this->x /= other.x;
            }

            if (other.y != 0) {
                this->y /= other.y;
            }

            return (*this);
        }

        Vector3_t & operator /= (const Vector3_t & other) {
            if (other.x != 0) {
                x /= other.x;
            }

            if (other.y != 0) {
                y /= other.y;
            }

            if (other.z != 0) {
                z /= other.z;
            }

            return (*this);
        }

        Vector3_t & operator /= (T n) {
            x /= n;
            y /= n;
            z /= n;

            return (*this);
        }

        Vector3_t & operator ++ (void) {
            x++; y++; z++;
            return (*this);
        }

        Vector3_t & operator -- (void) {
            x--; y--; z--;
            return (*this);
        }

        Vector3_t operator + (const Vector2_t<T> & other) const {
            return Vector3_t(this->x + other.x, this->y + other.y, this->z);
        }

        Vector3_t operator + (const Vector3_t & other) const {
            return Vector3_t(this->x + other.x, this->y + other.y, this->z + other.z);
        }

        Vector3_t operator - (const Vector2_t<T> & other) const {
            return Vector3_t(this->x - other.x, this->y - other.y, this->z);
        }

        Vector3_t operator - (const Vector3_t & other) const {
            return Vector3_t(this->x - other.x, this->y - other.y, this->z - other.z);
        }

        Vector3_t operator - (void) const {
            return Vector3_t(-this->x, -this->y, -this->z);
        }

        Vector3_t operator * (const Vector2_t<T> & other) const {
            return Vector3_t(this->x * other.x, this->y * other.y, this->z);
        }

        Vector3_t operator * (const Vector3_t & other) const {
            return Vector3_t(this->x * other.x, this->y * other.y, this->z * other.z);
        }

        Vector3_t operator * (T n) const {
            return Vector3_t(this->x * n, this->y * n, this->z * n);
        }

        Vector3_t operator / (const Vector2_t<T> & other) const {
            return Vector3_t(this->x / other.x, this->y / other.y, this->z);
        }

        Vector3_t operator / (const Vector3_t & other) const {
            return Vector3_t(this->x / other.x, this->y / other.y, this->z / other.z);
        }

        Vector3_t operator / (T n) const {
            return Vector3_t(this->x / n, this->y / n, this->z / n);
        }

        bool operator ! (void) const {
            return (this->x != 0 && this->y != 0 && this->z != 0);
        }

        bool operator == (const Vector3_t & other) const {
            return (this->x == other.x && this->y == other.y && this->z == other.z);
        }

        bool operator != (const Vector3_t & other) const {
            return (this->x != other.x || this->y != other.y || this->z != other.z);
        }

        T & operator [] (unsigned int n) {
            if (n > 2) { throw "Out of bounds."; }
            return *(&x + n);
        }

        const T & operator [] (const unsigned int n) const {
            if (n > 2) { throw "Out of bounds."; }
            return *(&x + n);
        }

        T SquaredMagnitude(void) const {
            return pow(this->x, 2.0f) + pow(this->y, 2.0f) + pow(this->z, 2.0f);
        }

        T Magnitude(void) const {
            return sqrt(this->SquaredMagnitude());
        }

        Vector3_t & Normalise(void) {
            if (!(*this)) {
                return (*this);
            }

            return (*this /= Magnitude());
        }

        Vector3_t Normalised(void) const {
            T m = Magnitude();
            return Vector3_t(this->x / m, this->y / m, this->z / m);
        }

        /** All of the following rotation functions were calculated using the wikipedia page:
         ** http://en.wikipedia.org/wiki/Rotation_matrix */

        /// <summary>Rotates the vector by the given radians about the X axis.</summary>
        /// <param name="angle">The angle in radians.</param>
        /// <returns>The 3 dimensional vector rotated by a supplied radians.</returns>
        Vector3_t & RotateAboutX(T angle) {
            T s = sin(angle);
            T c = cos(angle);

            Set(x, (y * c) - (z * s), (y * s) + (z * c));

            return (*this);
        }

        /// <summary>Rotates the vector by the given radians about the Y axis.</summary>
        /// <param name="angle">The angle in radians.</param>
        /// <returns>The 3 dimensional vector rotated by a supplied radians.</returns>
        Vector3_t & RotateAboutY(T angle) {
            T s = sin(angle);
            T c = cos(angle);

            Set((x * c) + (z * s), y, (z * c) - (x * s));

            return (*this);
        }

        /// <summary>Rotates the vector by the given radians about the Z axis.</summary>
        /// <param name="angle">The angle in radians.</param>
        /// <returns>The 3 dimensional vector rotated by a supplied radians.</returns>
        Vector3_t & RotateAboutZ(T angle) {
            T s = sin(angle);
            T c = cos(angle);

            Set((x * c) - (y * s), (x * s) + (y * c), z);

            return (*this);
        }
        
        /// <summary>Rotates the vector by the given radians about the X, Y and Z axis.</summary>
        /// <param name="xAngle">The angle in radians to rotate around the X axis.</param>
        /// <param name="yAngle">The angle in radians to rotate around the Y axis.</param>
        /// <param name="zAngle">The angle in radians to rotate around the Z axis.</param>
        /// <returns>The 3 dimensional vector rotated by a supplied radians.</returns>
        Vector3_t & RotateByXYZ(T xAngle, T yAngle, T zAngle) {
            RotateAboutX(xAngle);
            RotateAboutY(yAngle);
            RotateAboutZ(zAngle);

            return (*this);
        }

        /// <summary>Rotates the vector by the given radians about the given axis.</summary>
        /// <param name="axis">The axis in which to rotate around.</param>
        /// <param name="angle">The angle in radians to rotate.</param>
        /// <returns>The 3 dimensional vector rotated by the given radians.</returns>
        Vector3_t & RotateAboutAxis(const Vector3_t & axis, T angle) {
            Vector3_t u = axis.Normalised();

            T s = sin(angle);
            T c = cos(angle);
            T mc = 1 - c;

            Set((x * (c + (pow(u.x, 2) * mc))) + (y * ((u.x * u.y * mc) - (u.z * s))) + (z * ((u.x * u.z * mc) + (u.y * s))),
                (x * ((u.y * u.x * mc) + (u.z * s))) + (y * (c + (pow(u.y, 2) * mc))) + (z * ((u.y * u.z * mc) - (u.x * s))), 
                (x * ((u.z * u.x * mc) - (u.y * s))) + (y * ((u.z * u.y * mc) + (u.x * s))) + (z * (c + (pow(u.z, 2) * mc))));

            return (*this);
        }

        /// <summary>Returns the sum of the X, Y and Z components.</summary>
        T Sum(void) const {
            return x + y + z;
        }

        /// <summary>Calculates the dot product when given another 3 dimensional vector.</summary>
        /// <param name="other">The 3 dimensional vector in which to calculate the dot product with.</param>
        /// <returns>The dot product between two, 3 dimensional vectors.</returns>
        T DotProduct(const Vector3_t & other) const {
            return ((*this) + other).Sum();
        }

        /** The cross product function was calculated using the wikipedia page: 
         ** http://en.wikipedia.org/wiki/Cross_product */
        /// <summary>Calculates the cross product when given another 3 dimensional vector.</summary>
        /// <param name="other">The 3 dimensional vector in which to calculate the cross product with.</param>
        /// <returns>The cross product between two, 3 dimensional vectors.</returns>
        Vector3_t CrossProduct(const Vector3_t & other) const {
            return Vector3_t(
                (this->y * other.z) - (this->z * other.y),
                (this->z * other.x) - (this->x * other.z),
                (this->x * other.y) - (this->y * other.x)
            );
        }

        static Vector3_t Lerp(const Vector3_t & vecA, const Vector3_t & vecB, T p) {
            return ((vecA * (1 - p)) + (vecB * p));
        }
};

template<typename T>
inline Vector3_t<T> operator * (T n, const Vector3_t<T> & v) {
    return v * n;
}

typedef Vector3_t<float>  Vector3f;
typedef Vector3_t<double> Vector3d;
typedef Vector3_t<int>    Vector3i;