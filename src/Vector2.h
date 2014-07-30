#pragma once

/********************************************************************
    Author:         Paul Watkinson
    Website:        http://paulwatkinson.co.uk/

    Description:    A 2 dimensional vector class
*********************************************************************/

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdexcept>

//////////////////////////////////////////////////////////////////////////
/// A vector class for representing 2 dimensional vectors
//////////////////////////////////////////////////////////////////////////
template <typename T>
struct Vector2_t {
    public:
        T x, y;
        
        /// Default constructor
        Vector2_t(void) {
            x = y = static_cast<T>(0);
        }

        /// Destructor
        ~Vector2_t(void) { }

        /// Initialise the values of the vector to those of another
        template <typename U>
        Vector2_t(Vector2_t<U> & other) {
            this->x = static_cast<T>(other.x);
            this->y = static_cast<T>(other.y);
        }

        /// Attempt to get the values of another unknown type and set the vector's values to them
        /// WARNING: Could cause undefined functionality
        template<typename U>
        Vector2_t(const U & unknown) {
            memcpy(&x, &unknown.x, sizeof(T) * 2);
        }

        /// Initialise the vector with the given X and Y values
        Vector2_t(T a, T b) {
            x = a; y = b;
        }

        /// Set the vector's X and Y values to those given
        Vector2_t & Set(T a, T b) {
            this->x = a;
            this->y = b;

            return (*this);
        }

        /// Set the values of the vector to those of another
        Vector2_t & Set(const Vector2_t & other) {
            this->x = other.x;
            this->y = other.y;

            return (*this);
        }

        /// Attempt to get the values of another unknown type and set the vectors values to them
        /// WARNING: Could cause undefined functionality
        template<typename U>
        Vector2_t & Set(const U & unknown) {
            memcpy(&x, &unknown.x, sizeof(T) * 2);

            return (*this);
        }

        /// Set the values of the vector to those of another
        template <typename U>
        Vector2_t & operator  = (const Vector2_t<U> & other) {
            this->x = static_cast<T>(other.x);
            this->y = static_cast<T>(other.y);
        }

        /// Add the values of the vector and another's and set the vector to the result
        Vector2_t & operator += (const Vector2_t & other) {
            this->x += other.x;
            this->y += other.y;

            return (*this);
        }

        /// Negate the values of the vector and another's and set the vector to the result
        Vector2_t & operator -= (const Vector2_t & other) {
            this->x -= other.x;
            this->y -= other.y;

            return (*this);
        }

        /// Multiply the values of the vector and another's and set the vector to the result
        Vector2_t & operator *= (const Vector2_t & other) {
            this->x *= other.x;
            this->y *= other.y;

            return (*this);
        }

        /// Multiply the values of the vector by a value and set the vector to the result
        Vector2_t & operator *= (T n) {
            this->x *= n;
            this->y *= n;

            return (*this);
        }

        /// Divide the values of the vector by another's and set the vector to the result
        Vector2_t & operator /= (const Vector2_t & other) {
            if (other.x != 0) {
                this->x /= other.x;
            }
            
            if (other.y != 0) {
                this->y /= other.y;
            }

            return (*this);
        }

        /// Divide the values of the vector by a value and set the vector to the result
        Vector2_t & operator /= (T n) {
            if (!n) { return (*this); }
            this->x /= n;
            this->y /= n;

            return (*this);
        }

        /// Increment all values of the vector by one
        Vector2_t & operator ++ (void) {
            this->x++;
            this->y++;

            return (*this);
        }

        /// Decrement all values of the vector by one
        Vector2_t & operator -- (void) {
            this->x--;
            this->y--;

            return (*this);
        }

        /// Add the vector to another
        Vector2_t operator + (const Vector2_t & other) const {
            return Vector2_t(this->x + other.x, this->y + other.y);
        }

        /// Negate the vector by another
        Vector2_t operator - (const Vector2_t & other) const {
            return Vector2_t(this->x - other.x, this->y - other.y);
        }

        /// Invert the sign of all values of the vector
        Vector2_t operator - (void) const {
            return Vector2_t(-this->x, -this->y);
        }

        /// Multiply the vector by another
        Vector2_t operator * (const Vector2_t & other) const {
            return Vector2_t(this->x * other.x, this->y * other.y);
        }

        /// Multiply the vector by a value
        Vector2_t operator * (T n) const {
            return Vector2_t(this->x * n, this->y * n);
        }

        /// Divide the vector by another
        Vector2_t operator / (const Vector2_t & other) const {
            return Vector2_t(this->x / other.x, this->y / other.y);
        }

        /// Divide the vector by a value
        Vector2_t operator / (T n) const {
            return Vector2_t(this->x / n, this->y / n);
        }

        /// Check if the vector has all non-zero values
        bool operator ! (void) const {
            return (this->x != 0 && this->y != 0);
        }

        /// Check if the vector is equal to another
        bool operator == (const Vector2_t & other) const {
            return (this->x == other.x && this->y == other.y);
        }

        /// Check if the vector is not equal to another
        bool operator != (const Vector2_t & other) const {
            return (this->x != other.x || this->y != other.y);
        }

        /// Access a value of the vector
        T & operator [] (unsigned int n) {
            if (n >= 2) { throw "Out of bounds."; }
            return *(&x + n);
        }

        /// Access a value of the vector
        const T & operator [] (unsigned int n) const {
            if (n >= 2) { throw "Out of bounds."; }
            return *(&x + n);
        }

        /// Calculate and return the 2 dimensional vector that is perpendicular to the current vector, facing in the clockwise direction.
        Vector2_t RightNormal(void) const {
            return Vector2_t(-this->y, this->x).Normalise();
        }

        /// Calculate and return the 2 dimensional vector that is perpendicular to the current vector, facing in the counter-clockwise direction.
        Vector2_t LeftNormal(void) const {
            return -RightNormal();
        }

        /// Get the squared magnitude of the vector
        T SquaredMagnitude(void) const {
            return pow(this->x, 2.0f) + pow(this->y, 2.0f);
        }

        /// Get the magnitude of the vector
        T Magnitude(void) const {
            return sqrt(this->SquaredMagnitude());
        }

        /// Normalise the vector
        Vector2_t & Normalise(void) {
            if (!(*this)) {
                return (*this);
            }

            return (*this /= Magnitude());
        }

        /// Get the normalised vector
        Vector2_t Normalised(void) const {
            T m = Magnitude();
            return Vector2_t(this->x / m, this->y / m);
        }

        /// Rotate the vector by the given radians
        Vector2_t Rotate(T angle) {
            T c = cos(angle);
            T s = sin(angle);

            T nx = (c * this->x) - (s * this->y);
            T ny = (this->x * s) + (c * this->y);

            this->x = nx;
            this->y = ny;

            return (*this);
        }

        /// Return the sum of the X and Y components
        T Sum(void) const {
            return x + y;
        }

        /// Calculate the dot product with another 2D vector
        T DotProduct(const Vector2_t & other) const {
            return Vector2_t((*this) * other).Sum();
        }

        /// Calculate the cross product with another 2D vector
        /// Note: 2D cross products are not well defined. The code was derived from
        /// the formula given at http://mathworld.wolfram.com/CrossProduct.html
        T CrossProduct(const Vector2_t & other) const {
            return (this->x * other.y) - (this->y * other.x);
        }

        /// Linearly interpolate the values of the two vectors with the given percentage value
        static Vector2_t Lerp(const Vector2_t & vecA, const Vector2_t & vecB, T p) {
            return ((vecA * (1 - p)) + (vecB * p));
        }
};

template <typename T>
inline Vector2_t<T> operator * (T n, const Vector2_t<T> & v) {
    return v * n;
}

typedef Vector2_t<float>  Vector2f;
typedef Vector2_t<double> Vector2d;
typedef Vector2_t<int>    Vector2i;