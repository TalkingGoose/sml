#pragma once

/********************************************************************
    Author:         Paul Watkinson
    Website:        http://paulwatkinson.co.uk/

    Description:    A 4 dimensional vector class
*********************************************************************/

#include "Vector3.h"

//////////////////////////////////////////////////////////////////////////
/// A vector class for representing 4 dimensional vectors
//////////////////////////////////////////////////////////////////////////
template <typename T>
struct Vector4_t {
    public:
        T x, y, z, w;

        /// Default constructor
        Vector4_t(void) {
            x = y = z = w = static_cast<T>(0);
        }

        /// Destructor
        ~Vector4_t(void) { }
        
        /// Initialise the values of the vector to those of another
        template <typename U>
        Vector4_t(Vector4_t<U> & other) {
            this->x = static_cast<T>(other.x);
            this->y = static_cast<T>(other.y);
            this->z = static_cast<T>(other.z);
            this->w = static_cast<T>(other.w);
        }

        /// Initialise the vector with the given X, Y, Z and W values
        Vector4_t(T a, T b, T c, T d) {
            x = a; y = b; z = c; w = d;
        }

        /// Initialise the 4D vector with a 3D vector and a W value
        Vector4_t(const Vector3_t<T> & other, T d) {
            this->x = other.x;
            this->y = other.y;
            this->z = other.z;
            this->w = d;
        }

        /// Set the vector's X, Y, Z and W values to those given
        Vector4_t & Set(T a, T b, T c, T d) {
            this->x = a;
            this->y = b;
            this->z = c;
            this->w = d;

            return (*this);
        }

        /// Set the 4D vector's values to those of a 3D vector's and a W value
        Vector4_t & Set(const Vector3_t<T> other, T d = 0.0l) {
            this->x = other.x;
            this->y = other.y;
            this->z = other.z;
            this->w = d;

            return (*this);
        }

        /// Set the values of the vector to those of another
        Vector4_t & Set(const Vector4_t & other) {
            this->x = other.x;
            this->y = other.y;
            this->z = other.z;
            this->w = other.w;

            return (*this);
        }

        /// Attempt to get the values of another unknown type and set the vectors values to them
        /// WARNING: Could cause undefined functionality
        template<typename U>
        Vector4_t & Set(const U & unknown) {
            memcpy(&x, &unknown.x, sizeof(T) * 3);

            return (*this);
        }

        /// Set the values of the 4D vector to those of a 3D vector
        Vector4_t & operator = (const Vector3_t<T> & other) {
            this->x = other.x;
            this->y = other.y;
            this->z = other.z;
            this->w = 0.0f;

            return (*this);
        }

        /// Set the values of the 4D vector to those of another
        Vector4_t & operator = (const Vector4_t & other) {
            this->x = other.x;
            this->y = other.y;
            this->z = other.z;
            this->w = other.w;

            return (*this);
        }

        //////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////// 
        
        /// Add the values of the 3D vector and a 2D vector's and set the 3D vector to the result
        Vector4_t & operator += (const Vector3_t<T> & other) {
            this->x += other.x;
            this->y += other.y;
            this->z += other.z;

            return (*this);
        }

        /// Add the values of the vector and another's and set the vector to the result
        Vector4_t & operator += (const Vector4_t & other) {
            this->x += other.x;
            this->y += other.y;
            this->z += other.z;
            this->w += other.w;

            return (*this);
        }

        /// Negate the values of the 3D vector and a 2D vector's and set the 3D vector to the result
        Vector4_t & operator -= (const Vector3_t<T> & other) {
            this->x -= other.x;
            this->y -= other.y;
            this->w -= other.z;

            return (*this);
        }

        /// Negate the values of the vector and another's and set the vector to the result
        Vector4_t & operator -= (const Vector4_t & other) {
            this->x -= other.x;
            this->y -= other.y;
            this->z -= other.z;
            this->w -= other.w;

            return (*this);
        }

        /// Multiply the values of the 3D vector and a 2D vector's and set the 3D vector to the result
        Vector4_t & operator *= (const Vector3_t<T> & other) {
            this->x *= other.x;
            this->y *= other.y;
            this->z *= other.z;

            return (*this);
        }

        /// Multiply the values of the vector and another's and set the vector to the result
        Vector4_t & operator *= (const Vector4_t & other) {
            this->x *= other.x;
            this->y *= other.y;
            this->z *= other.z;
            this->w *= other.w;

            return (*this);
        }

        /// Multiply the values of the vector by a value and set the vector to the result
        Vector4_t & operator *= (T n) {
            this->x *= n;
            this->y *= n;
            this->z *= n;
            this->w *= n;

            return (*this);
        }

        /// Divide the values of the 3D vector by a 2D vector's and set the 3D vector to the result
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

        /// Divide the values of the vector by another's and set the vector to the result
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

        /// Divide the values of the vector by a value and set the vector to the result
        Vector4_t & operator /= (T n) {
            if (!n) { return (*this); }
            x /= n;
            y /= n;
            z /= n;
            w /= n;

            return (*this);
        }

        /// Increment all values of the vector by one
        Vector4_t & operator ++ (void) {
            x++; y++; z++; w++;
            return (*this);
        }

        /// Decrement all values of the vector by one
        Vector4_t & operator -- (void) {
            x--; y--; z--; w--;
            return (*this);
        }

        /// Add the 3D vector to a 2D vector
        Vector4_t operator + (const Vector3_t<T> & other) const {
            return Vector4_t(this->x + other.x, this->y + other.y, this->z + other.z, this->w);
        }

        /// Add the vector to another
        Vector4_t operator + (const Vector4_t & other) const {
            return Vector4_t(this->x + other.x, this->y + other.y, this->z + other.z, this->w + other.w);
        }

        /// Negate the 3D vector by a 2D vector
        Vector4_t operator - (const Vector3_t<T> & other) const {
            return Vector4_t(this->x - other.x, this->y - other.y, this->z - other.z, this->w);
        }

        /// Negate the vector by another
        Vector4_t operator - (const Vector4_t & other) const {
            return Vector4_t(this->x - other.x, this->y - other.y, this->z - other.z, this->w - other.w);
        }

        /// Invert the sign of all values of the vector
        Vector4_t operator - (void) const {
            return Vector4_t(-this->x, -this->y, -this->z, -this->w);
        }

        /// Multiply the 3D vector by a 2D vector
        Vector4_t operator * (const Vector3_t<T> & other) const {
            return Vector4_t(this->x * other.x, this->y * other.y, this->z * other.z, this->w);
        }

        /// Multiply the vector by another
        Vector4_t operator * (const Vector4_t & other) const {
            return Vector4_t(this->x * other.x, this->y * other.y, this->z * other.z, this->w * other.w);
        }

        /// Multiply the vector by a value
        Vector4_t operator * (T n) const {
            return Vector4_t(this->x * n, this->y * n, this->z * n, this->w * n);
        }

        /// Divide the 3D vector by a 2D vector
        Vector4_t operator / (const Vector3_t<T> & other) const {
            return Vector4_t(this->x / other.x, this->y / other.y, this->z / other.z, this->w);
        }

        /// Divide the vector by another
        Vector4_t operator / (const Vector4_t & other) const {
            return Vector4_t(this->x / other.x, this->y / other.y, this->z / other.z, this->w / other.w);
        }

        /// Divide the vector by a value
        Vector4_t operator / (T n) const {
            return Vector4_t(this->x / n, this->y / n, this->z / n, this->w / n);
        }

        /// Check if the vector has all non-zero values
        bool operator ! (void) const {
            return (this->x != 0 && this->y != 0 && this->z != 0 && this->w != 0);
        }

        /// Check if the vector is equal to another
        bool operator == (const Vector4_t & other) const {
            return (this->x == other.x && this->y == other.y && this->z == other.z && this->w == other.w);
        }

        /// Check if the vector is not equal to another
        bool operator != (const Vector4_t & other) const {
            return (this->x != other.x || this->y != other.y || this->z != other.z || this->w != other.w);
        }

        /// Access a value of the vector
        T & operator [] (unsigned int n) {
            if (n >= 4) { throw "Out of bounds."; }
            return *(&x + n);
        }

        /// Access a value of the vector
        const T & operator [] (const unsigned int n) const {
            if (n >= 4) { throw "Out of bounds."; }
            return *(&x + n);
        }

        /// Get the squared magnitude of the vector
        T SquaredMagnitude(void) const {
            return pow(this->x, 2.0f) + pow(this->y, 2.0f) + pow(this->z, 2.0f) + (this->w, 2.0f);
        }

        /// Get the magnitude of the vector
        T Magnitude(void) const {
            return sqrt(this->SquaredMagnitude());
        }

        /// Normalise the vector
        Vector4_t & Normalise(void) {
            if (!(*this)) {
                return (*this);
            }

            return (*this /= Magnitude());
        }

        /// Get the normalised vector
        Vector4_t Normalised(void) const {
            T m = Magnitude();
            return Vector4_t(this->x / m, this->y / m, this->z / m, this->w / m);
        }

        /// Return the sum of the X, Y, Z and W components
        T Sum(void) const {
            return x + y + z + w;
        }

        /// Calculate the dot product when given another 3D vector
        Vector4_t DotProduct(const Vector4_t & other) const {
            return ((*this) + other).Sum();
        }

        /// Calculate the cross product when given another 3D vector
        Vector4_t CrossProduct(const Vector4_t & other) const {
            return Vector4_t(
                (this->y * other.z) - (this->z * other.y),
                (this->z * other.x) - (this->x * other.z),
                (this->x * other.y) - (this->y * other.x),
                static_cast<T>(0.0f)
            );
        }

        /// Linearly interpolate the values of the two vectors with the given percentage value
        static Vector4_t Lerp(const Vector4_t & vecA, const Vector4_t & vecB, T p) {
            return ((vecA * (1 - p)) + (vecB * p));
        }
};

template <typename T>
inline Vector4_t<T> operator * (T n, const Vector4_t<T> & v) {
    return v * n;
}

typedef Vector4_t<float>  Vector4f;
typedef Vector4_t<double> Vector4d;
typedef Vector4_t<int>    Vector4i;