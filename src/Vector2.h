#pragma once

/********************************************************************
    Author:         Paul Watkinson
    Website:        http://paulwatkinson.co.uk/

    Description:    A 2 dimensional vector class
*********************************************************************/

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdexcept>

/// A vector class for representing 2 dimensional vectors
template <typename T>
struct Vector2_t {
    public:
        union {
            struct {
                T x, y;
            };

            T linear[2];
        };
        
        /// <summary>Empty constructor</summary>
        Vector2_t(void) {
            x = y = static_cast<T>(0);
        }

        template <typename U>
        Vector2_t(Vector2_t<U> & other) {
            this->x = static_cast<T>(other.x);
            this->y = static_cast<T>(other.y);
        }

        template<typename U>
        Vector2_t(const U & unknown) {
            memcpy(&x, &unknown.x, sizeof(T) * 2);
        }

        /// <summary>Construct a 2 dimensional vector given the X and Y values.</summary>
        /// <param name="a">The X value to be set.</param>
        /// <param name="b">The Y value to be set.</param>
        Vector2_t(T a, T b) {
            x = a; y = b;
        }

        /// <summary>Set the X and Y values of the 2 dimensional vector.</summary>
        /// <param name="a">The X value to be set.</param>
        /// <param name="b">The Y value to be set.</param>
        Vector2_t & Set(T a, T b) {
            this->x = a;
            this->y = b;

            return (*this);
        }

        /// <summary>Set the X and Y values of the 2 dimensional vector to that of anothers.</summary>
        /// <param name="other">The 2 dimensional vector to be set the X and Y values with.</param>
        Vector2_t & Set(const Vector2_t & other) {
            this->x = other.x;
            this->y = other.y;

            return (*this);
        }

        template<typename U>
        Vector2_t & Set(const U & unknown) {
            memcpy(&x, &unknown.x, sizeof(T) * 2);

            return (*this);
        }

        Vector2_t & operator  = (const Vector2_t & other) {
            this->x = other.x;
            this->y = other.y;

            return (*this);
        }

        template <typename U>
        Vector2_t & operator  = (const Vector2_t<U> & other) {
            this->x = static_cast<T>(other.x);
            this->y = static_cast<T>(other.y);
        }

        Vector2_t & operator += (const Vector2_t & other) {
            this->x += other.x;
            this->y += other.y;

            return (*this);
        }

        Vector2_t & operator -= (const Vector2_t & other) {
            this->x -= other.x;
            this->y -= other.y;

            return (*this);
        }

        Vector2_t & operator *= (const Vector2_t & other) {
            this->x *= other.x;
            this->y *= other.y;

            return (*this);
        }

        Vector2_t & operator *= (T n) {
            this->x *= n;
            this->y *= n;

            return (*this);
        }

        Vector2_t & operator /= (const Vector2_t & other) {
            this->x /= other.x;
            this->y /= other.y;

            return (*this);
        }

        Vector2_t & operator /= (T n) {
            this->x /= n;
            this->y /= n;

            return (*this);
        }

        Vector2_t & operator ++ (void) {
            this->x++;
            this->y++;

            return (*this);
        }

        Vector2_t & operator -- (void) {
            this->x--;
            this->y--;

            return (*this);
        }

        Vector2_t operator + (const Vector2_t & other) const {
            return Vector2_t(this->x + other.x, this->y + other.y);
        }

        Vector2_t operator - (const Vector2_t & other) const {
            return Vector2_t(this->x - other.x, this->y - other.y);
        }

        Vector2_t operator - (void) const {
            return Vector2_t(-this->x, -this->y);
        }

        Vector2_t operator * (const Vector2_t & other) const {
            return Vector2_t(this->x * other.x, this->y * other.y);
        }

        Vector2_t operator * (T n) const {
            return Vector2_t(this->x * n, this->y * n);
        }

        Vector2_t operator / (const Vector2_t & other) const {
            return Vector2_t(this->x / other.x, this->y / other.y);
        }

        Vector2_t operator / (T n) const {
            return Vector2_t(this->x / n, this->y / n);
        }

        bool operator ! (void) const {
            return (this->x != 0 && this->y != 0);
        }

        bool operator == (const Vector2_t & other) const {
            return (this->x == other.x && this->y == other.y);
        }

        bool operator != (const Vector2_t & other) const {
            return (this->x != other.x || this->y != other.y);
        }

        T & operator [] (unsigned int n) {
            return *(&x + n);
        }

        const T & operator [] (unsigned int n) const {
            return *(&x + n);
        }

        /// <summary>Calculates and returns the 2 dimensional vector that is perpendicular to this one, 
        /// facing in the clockwise direction.</summary>
        /// <returns>A 2 dimensional vector.</returns>
        Vector2_t RightNormal(void) const {
            return Vector2_t(-this->y, this->x).Normalise();
        }

        /// <summary>Calculates and returns the 2 dimensional vector that is perpendicular to this one, 
        /// facing in the counterclockwise direction.</summary>
        /// <returns>A 2 dimensional vector.</returns>
        Vector2_t LeftNormal(void) const {
            return -RightNormal();
        }

        T SquaredMagnitude(void) const {
            return pow(this->x, 2.0f) + pow(this->y, 2.0f);
        }

        T Magnitude(void) const {
            return sqrt(this->SquaredMagnitude());
        }

        Vector2_t & Normalise(void) {
            if (!(*this)) {
                return (*this);
            }

            return (*this /= Magnitude());
        }

        Vector2_t Normalised(void) const {
            T m = Magnitude();
            return Vector2_t(this->x / m, this->y / m);
        }

        /// <summary>Rotates the vector by the given radians.</summary>
        /// <param name="angle">The angle in radians.</param>
        /// <returns>The 2 dimensional vector rotated by the supplied radians.</returns>
        Vector2_t Rotate(T angle) {
            T c = cos(angle);
            T s = sin(angle);

            T nx = (c * this->x) - (s * this->y);
            T ny = (this->x * s) + (c * this->y);

            this->x = nx;
            this->y = ny;

            return (*this);
        }

        /// <summary>Returns the sum of the X and Y components.</summary>
        T Sum(void) const {
            return x + y;
        }

        /// <summary>Calculates the dot product when given another 2 dimensional vector.</summary>
        /// <param name="other">The 2 dimensional vector in which to calculate the dot product with.</param>
        /// <returns>The dot product between two, 2 dimensional vectors.</returns>
        T DotProduct(const Vector2_t & other) const {
            return Vector2_t((*this) * other).Sum();
        }

        /** Note: 2D cross products are not well defined. The code was derived from
         ** the formula given at http://mathworld.wolfram.com/CrossProduct.html */
        /// <summary>Calculates the cross product when given another 2 dimensional vector.</summary>
        /// <param name="other">The 2 dimensional vector in which to calculate the cross product with.</param>
        /// <returns>The dot product between two, 2 dimensional vectors.</returns>
        T CrossProduct(const Vector2_t & other) const {
            return (this->x * other.y) - (this->y * other.x);
        }
};

template <typename T>
inline Vector2_t<T> operator * (T n, const Vector2_t<T> & v) {
    return v * n;
}

typedef Vector2_t<float>  Vector2f;
typedef Vector2_t<double> Vector2d;
typedef Vector2_t<int>    Vector2i;