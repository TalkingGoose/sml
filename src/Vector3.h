#pragma once

/********************************************************************
    Author:         Paul Watkinson
    Website:        http://paulwatkinson.co.uk/

    Description:    A 3 dimensional vector class
*********************************************************************/

#include "Vector2.h"

//////////////////////////////////////////////////////////////////////////
/// A vector class for representing 3 dimensional vectors
//////////////////////////////////////////////////////////////////////////
template<typename T>
struct Vector3_t {
    public:
        T x, y, z;

        /// Default constructor
        Vector3_t(void) {
            x = y = z = static_cast<T>(0);
        }
        
        /// Destructor
        ~Vector3_t(void) { }

        /// Initialise the values of the vector to those of another
        template <typename U>
        Vector3_t(Vector3_t<U> & other) {
            this->x = static_cast<T>(other.x);
            this->y = static_cast<T>(other.y);
            this->z = static_cast<T>(other.z);
        }

        /// Initialise the vector with the given X, Y and Z values
        Vector3_t(T a, T b, T c) {
            x = a; y = b; z = c;
        }

        /// Initialise the 3D vector with a 2D vector and a Z value
        Vector3_t(const Vector2_t<T> & other, T c) {
            this->x = other.x;
            this->y = other.y;
            this->z = c;
        }

        /// Set the vector's X, Y and Z values to those given
        Vector3_t & Set(T a, T b, T c) {
            this->x = a;
            this->y = b;
            this->z = c;

            return (*this);
        }

        /// Set the 3D vector's values to those of a 2D vector's and a Z value
        Vector3_t & Set(const Vector2_t<T> other, T c = 0.0f) {
            this->x = other.x;
            this->y = other.y;
            this->z = c;

            return (*this);
        }

        /// Set the values of the vector to those of another
        Vector3_t & Set(const Vector3_t & other) {
            this->x = other.x;
            this->y = other.y;
            this->z = other.z;

            return (*this);
        }

        /// Set the values of the 3D vector to those of a 2D vector
        Vector3_t & operator = (const Vector2_t<T> & other) {
            this->x = other.x;
            this->y = other.y;
            this->z = 0.0f;

            return (*this);
        }

        /// Set the values of the 3D vector to those of another
        Vector3_t & operator = (const Vector3_t & other) {
            this->x = other.x;
            this->y = other.y;
            this->z = other.z;

            return (*this);
        }
        
        /// Add the values of the 3D vector and a 2D vector's and set the 3D vector to the result
        Vector3_t & operator += (const Vector2_t<T> & other) {
            this->x += other.x;
            this->y += other.y;

            return (*this);
        }

        /// Add the values of the vector and another's and set the vector to the result
        Vector3_t & operator += (const Vector3_t & other) {
            this->x += other.x;
            this->y += other.y;
            this->z += other.z;

            return (*this);
        }

        /// Negate the values of the 3D vector and a 2D vector's and set the 3D vector to the result
        Vector3_t & operator -= (const Vector2_t<T> & other) {
            this->x -= other.x;
            this->y -= other.y;

            return (*this);
        }

        /// Negate the values of the vector and another's and set the vector to the result
        Vector3_t & operator -= (const Vector3_t & other) {
            this->x -= other.x;
            this->y -= other.y;
            this->z -= other.z;

            return (*this);
        }

        /// Multiply the values of the 3D vector and a 2D vector's and set the 3D vector to the result
        Vector3_t & operator *= (const Vector2_t<T> & other) {
            this->x *= other.x;
            this->y *= other.y;

            return (*this);
        }

        /// Multiply the values of the vector and another's and set the vector to the result
        Vector3_t & operator *= (const Vector3_t & other) {
            this->x *= other.x;
            this->y *= other.y;
            this->z *= other.z;

            return (*this);
        }

        /// Multiply the values of the vector by a value and set the vector to the result
        Vector3_t & operator *= (T n) {
            this->x *= n;
            this->y *= n;
            this->z *= n;

            return (*this);
        }

        /// Divide the values of the 3D vector by a 2D vector's and set the 3D vector to the result
        Vector3_t & operator /= (const Vector2_t<T> & other) {
            if (other.x != 0) {
                this->x /= other.x;
            }

            if (other.y != 0) {
                this->y /= other.y;
            }

            return (*this);
        }

        /// Divide the values of the vector by another's and set the vector to the result
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

        /// Divide the values of the vector by a value and set the vector to the result
        Vector3_t & operator /= (T n) {
            if (!n) { return (*this); }
            x /= n;
            y /= n;
            z /= n;

            return (*this);
        }

        /// Increment all values of the vector by one
        Vector3_t & operator ++ (void) {
            x++; y++; z++;
            return (*this);
        }

        /// Decrement all values of the vector by one
        Vector3_t & operator -- (void) {
            x--; y--; z--;
            return (*this);
        }

        /// Add the 3D vector to a 2D vector
        Vector3_t operator + (const Vector2_t<T> & other) const {
            return Vector3_t(this->x + other.x, this->y + other.y, this->z);
        }

        /// Add the vector to another
        Vector3_t operator + (const Vector3_t & other) const {
            return Vector3_t(this->x + other.x, this->y + other.y, this->z + other.z);
        }

        /// Negate the 3D vector by a 2D vector
        Vector3_t operator - (const Vector2_t<T> & other) const {
            return Vector3_t(this->x - other.x, this->y - other.y, this->z);
        }

        /// Negate the vector by another
        Vector3_t operator - (const Vector3_t & other) const {
            return Vector3_t(this->x - other.x, this->y - other.y, this->z - other.z);
        }

        /// Invert the sign of all values of the vector
        Vector3_t operator - (void) const {
            return Vector3_t(-this->x, -this->y, -this->z);
        }

        /// Multiply the 3D vector by a 2D vector
        Vector3_t operator * (const Vector2_t<T> & other) const {
            return Vector3_t(this->x * other.x, this->y * other.y, this->z);
        }

        /// Multiply the vector by another
        Vector3_t operator * (const Vector3_t & other) const {
            return Vector3_t(this->x * other.x, this->y * other.y, this->z * other.z);
        }

        /// Multiply the vector by a value
        Vector3_t operator * (T n) const {
            return Vector3_t(this->x * n, this->y * n, this->z * n);
        }

        /// Divide the 3D vector by a 2D vector
        Vector3_t operator / (const Vector2_t<T> & other) const {
            return Vector3_t(this->x / other.x, this->y / other.y, this->z);
        }

        /// Divide the vector by another
        Vector3_t operator / (const Vector3_t & other) const {
            return Vector3_t(this->x / other.x, this->y / other.y, this->z / other.z);
        }

        /// Divide the vector by a value
        Vector3_t operator / (T n) const {
            return Vector3_t(this->x / n, this->y / n, this->z / n);
        }

        /// Check if the vector has all non-zero values
        bool operator ! (void) const {
            return (this->x != 0 && this->y != 0 && this->z != 0);
        }

        /// Check if the vector is equal to another
        bool operator == (const Vector3_t & other) const {
            return (this->x == other.x && this->y == other.y && this->z == other.z);
        }

        /// Check if the vector is not equal to another
        bool operator != (const Vector3_t & other) const {
            return (this->x != other.x || this->y != other.y || this->z != other.z);
        }

        /// Access a value of the vector
        T & operator [] (unsigned int n) {
            if (n >= 3) { throw "Out of bounds."; }
            return *(&x + n);
        }

        /// Access a value of the vector
        const T & operator [] (const unsigned int n) const {
            if (n >= 3) { throw "Out of bounds."; }
            return *(&x + n);
        }

        /// Get the squared magnitude of the vector
        T SquaredMagnitude(void) const {
            return pow(this->x, 2.0f) + pow(this->y, 2.0f) + pow(this->z, 2.0f);
        }

        /// Get the magnitude of the vector
        T Magnitude(void) const {
            return sqrt(this->SquaredMagnitude());
        }

        /// Normalise the vector
        Vector3_t & Normalise(void) {
            if (!(*this)) {
                return (*this);
            }

            return (*this /= Magnitude());
        }

        /// Get the normalised vector
        Vector3_t Normalised(void) const {
            T m = Magnitude();
            return Vector3_t(this->x / m, this->y / m, this->z / m);
        }

        /** All of the following rotation functions were calculated using the wikipedia page:
         ** http://en.wikipedia.org/wiki/Rotation_matrix */

        /// Rotate the vector by the given radians about the X axis
        Vector3_t & RotateAboutX(T angle) {
            T s = sin(angle);
            T c = cos(angle);

            Set(x, (y * c) - (z * s), (y * s) + (z * c));

            return (*this);
        }

        /// Rotate the vector by the given radians about the Y axis
        Vector3_t & RotateAboutY(T angle) {
            T s = sin(angle);
            T c = cos(angle);

            Set((x * c) + (z * s), y, (z * c) - (x * s));

            return (*this);
        }

        /// Rotate the vector by the given radians about the Z axis
        Vector3_t & RotateAboutZ(T angle) {
            T s = sin(angle);
            T c = cos(angle);

            Set((x * c) - (y * s), (x * s) + (y * c), z);

            return (*this);
        }
        
        /// Rotate the vector by the given radians about the X, Y and Z axis
        Vector3_t & RotateByXYZ(T xAngle, T yAngle, T zAngle) {
            RotateAboutX(xAngle);
            RotateAboutY(yAngle);
            RotateAboutZ(zAngle);

            return (*this);
        }

        // Rotate the vector by the given radians about the given axis
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

        /// Return the sum of the X, Y and Z components
        T Sum(void) const {
            return x + y + z;
        }

        /// Calculate the dot product when given another 3D vector
        T DotProduct(const Vector3_t & other) const {
            return ((*this) + other).Sum();
        }

        /** The cross product function was calculated using the wikipedia page: 
         ** http://en.wikipedia.org/wiki/Cross_product */
        /// Calculate the cross product when given another 3D vector
        Vector3_t CrossProduct(const Vector3_t & other) const {
            return Vector3_t(
                (this->y * other.z) - (this->z * other.y),
                (this->z * other.x) - (this->x * other.z),
                (this->x * other.y) - (this->y * other.x)
            );
        }

        /// Linearly interpolate the values of the two vectors with the given percentage value
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