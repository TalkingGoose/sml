#pragma once

/********************************************************************
    Author:         Paul Watkinson
    Website:        http://paulwatkinson.co.uk/

    Description:    A quaternion class
*********************************************************************/

#include "Vector3.h"

//////////////////////////////////////////////////////////////////////////
/// A class for representing quaternions
//////////////////////////////////////////////////////////////////////////
class Quaternion {
    private:
        // http://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
        /// Set the quaternion to the values calculated from Euler angles
        void FromEuler(double a, double b, double c) {
            double cosa = cos(a * 0.5f); double cosb = cos(b * 0.5f); double cosc = cos(c * 0.5f);
            double sina = sin(a * 0.5f); double sinb = sin(b * 0.5f); double sinc = sin(c * 0.5f);

            this->x = (cosa * cosb * cosc) + (sina * sinb * sinc);
            this->y = (sina * cosb * cosc) - (cosa * sinb * sinc);
            this->z = (cosa * sinb * cosc) + (sina * cosb * sinc);
            this->w = (cosa * cosb * sinc) - (sina * sinb * cosc);
        }

    public:
        double x, y, z, w;

        // Default constructor
        Quaternion(void) {
            this->x = this->y = this->z = 0.0f;
            this->w = 1.0f;
        }

        /// Initialise the quaternion with the values calculated from Euler angles
        Quaternion(double a, double b, double c) {
            FromEuler(a, b, c);
        }

        /// Initialise the quaternion with the given values
        Quaternion(double a, double b, double c, double d) {
            this->x = a; this->y = b; this->z = c; this->w = d;
        }
        
        /// Initialise the quaternion with the values calculated from Euler angles given by the vector
        Quaternion(Vector3d vector) {
            FromEuler(vector.x, vector.y, vector.z);
        }

        /// Initialise the quaternion with a vector and a W value
        Quaternion(Vector3d vector, double d) {
            this->x = vector.x; this->y = vector.y; this->z = vector.z; this->w = d;
        }

        /// Destructor
        virtual ~Quaternion(void) { };

        /// Set the quaternion to the values calculated from Euler angles
        Quaternion & Set(double a, double b, double c) {
            FromEuler(a, b, c);

            return (*this);
        }

        /// Set the quaternion to the given values
        Quaternion & Set(double a, double b, double c, double d) {
            this->x = a;
            this->y = b;
            this->z = c;
            this->w = d;

            return (*this);
        }

        /// Set the quaternion to the values calculated from Euler angles given by the vector
        Quaternion & Set(const Vector3d vector) {
            FromEuler(vector.x, vector.y, vector.z);

            return (*this);
        }

        /// Set the quaternion with a vector and a W value
        Quaternion & Set(const Vector3d vector, double d) {
            this->x = vector.x;
            this->y = vector.y;
            this->z = vector.z;
            this->w = d;

            return (*this);
        }

        /// Set the quaternion to the given values
        Quaternion & Set(const Quaternion & quat) {
            this->x = quat.x;
            this->y = quat.y;
            this->z = quat.z;
            this->w = quat.w;

            return (*this);
        }

        /// Set the quaternion to those of another quaternion
        Quaternion & operator =  (const Quaternion & quat) {
            this->x = quat.x;
            this->y = quat.y;
            this->z = quat.z;
            this->w = quat.w;

            return (*this);
        }

        /// Add the values of the quaternion and another's and set the quaternion to the result
        Quaternion & operator += (const Quaternion & quat) {
            this->x += quat.x;
            this->y += quat.y;
            this->z += quat.z;
            this->w += quat.w;

            return (*this);
        }

        /// Negate the values of the quaternion and another's and set the quaternion to the result
        Quaternion & operator -= (const Quaternion & quat) {
            this->x -= quat.x;
            this->y -= quat.y;
            this->z -= quat.z;
            this->w -= quat.w;

            return (*this);
        }
        
        /// Multiply the values of the quaternion and another's and set the quaternion to the result
        Quaternion & operator *= (const Quaternion & quat) {
            return (*this) = ((*this) * quat);
        }

        /// Multiply the values of the quaternion by a value and set the quaternion to the result
        Quaternion & operator *= (double n) {
            this->x *= n;
            this->y *= n;
            this->z *= n;
            this->w *= n;

            return (*this);
        }

        /// Divide the values of the quaternion and another's and set the quaternion to the result
        Quaternion & operator /= (const Quaternion & quat) {
            return (*this) = ((*this) / quat);
        }

        /// Divide the values of the quaternion by a value and set the quaternion to the result
        Quaternion & operator /= (double n) {
            if (!n) { return (*this); }
            this->x /= n;
            this->y /= n;
            this->z /= n;
            this->w /= n;

            return (*this);
        }

        /// Add the quaternion to another
        Quaternion operator + (const Quaternion & quat) const {
            return Quaternion(this->x + quat.x, this->y + quat.y, this->z + quat.z, this->w + quat.w);
        }

        /// Negate the quaternion by another
        Quaternion operator - (const Quaternion & quat) const {
            return Quaternion(this->x - quat.x, this->y - quat.y, this->z - quat.z, this->w - quat.w);
        }

        /// Invert the sign of all values of the quaternion
        Quaternion operator - (void) const {
            return Quaternion(-this->x, -this->y, -this->z, -this->w);
        }

        /// Multiply the quaternion by another
        Quaternion operator * (const Quaternion & quat) const {
            return Quaternion(
                (w * quat.x) + (x * quat.w) + (y * quat.z) - (z * quat.y),
                (w * quat.y) + (y * quat.w) + (z * quat.x) - (x * quat.z),
                (w * quat.z) + (z * quat.w) + (x * quat.y) - (y * quat.x),
                (w * quat.w) - (x * quat.x) - (y * quat.y) - (z * quat.z)
            );
        }

        /// Multiply the quaternion by a value
        Quaternion operator * (double n) const {
            return Quaternion(this->x * n, this->y * n, this->z * n, this->w * n);
        }

        /// Divide the quaternion by another
        Quaternion operator / (const Quaternion & quat) const {
            return Quaternion(this->Inverted() * quat);
        }

        /// Divide the quaternion by a value
        Quaternion operator / (double n) const {
            return Quaternion(this->x / n, this->y / n, this->z / n, this->w / n);
        }

        /// Check if the quaternion has all non-zero values
        bool operator ! (void) const {
            return (this->x != 0 && this->y != 0 && this->z != 0 && this->w != 0);
        }

        /// Check if the quaternion is equal to another
        bool operator == (const Quaternion & quat) const {
            return (this->x == quat.x && this->y == quat.y && this->z == quat.z && this->w == quat.w);
        }

        /// Check if the quaternion is not equal to another
        bool operator != (const Quaternion & quat) const {
            return (this->x != quat.x && this->y != quat.y && this->z != quat.z && this->w != quat.w);
        }

        /// Access a value of the quaternion
        double & operator [] (unsigned int n) {
            if (n >= 4) { throw "Out of range."; }
            return *(&x + n);
        }

        /// Get the axis of rotation for the quaternion
        Vector3d GetAxis(void) const {
            return Vector3d(this->x, this->y, this->z);
        }

        /// Get the rotation around the quaternion's axis in radians
        double GetAngle(void) const {
            return acos(w) * 2;
        }

        /// Get the squared magnitude of the quaternion
        double SquaredMagnitude(void) const {
            return pow(this->x, 2.0f) + pow(this->y, 2.0f) + pow(this->z, 2.0f) + pow(this->w, 2.0f);
        }

        /// Get the magnitude of the quaternion
        double Magnitude(void) const {
            return sqrt(this->SquaredMagnitude());
        }

        /// Normalise the quaternion
        Quaternion & Normalise(void) {
            return (*this /= Magnitude());
        }

        /// Get the normalised quaternion
        Quaternion Normalised(void) const {
            double m = Magnitude();
            return Quaternion(this->x / m, this->y / m, this->z / m, this->w / m);
        }

        /// Get the inverted quaternion
        Quaternion Inverted(void) const {
            return this->Conjugate() / this->SquaredMagnitude();
        }

        /// Invert the quaternion
        Quaternion & Invert(void) {
            (*this).Set(this->Inverted());
            return (*this);
        }

        /// Get the conjugate of the quaternion
        Quaternion Conjugate(void) const {
            return Quaternion(-this->x, -this->y, -this->z, this->w);
        }

        /// Get the reciprocal of the quaternion
        Quaternion Reciprocal(void) const {
            return (this->Conjugate() / ((*this) * this->Conjugate()));
        }

        /// Calculate the dot product when given another quaternion
        double DotProduct(const Quaternion & quat) const {
            return (this->x * quat.x) + (this->y * quat.y) + (this->z * quat.z) + (this->w * quat.w);
        }

        /// Linearly interpolate the values of the two quaternions with the given percentage value
        static Quaternion Lerp(const Quaternion & quatA, const Quaternion & quatB, double p) {
            return ((quatA * (1 - p)) + (quatB * p)).Normalised();
        }

        /// Spherically interpolate the values of the two quaternions with the given percentage value
        static Quaternion Slerp(const Quaternion & quatA, const Quaternion & quatB, double p) {
            Quaternion q;
            double dot = quatA.DotProduct(quatB);

            // If too similar, lerp instead.
            if (abs(dot) > 0.95f) {
                return Lerp(quatA, quatB, p);
            }

            if (dot < 0) {
                q.Set(-quatB);
            } else {
                q.Set(quatB);
            }

            double angle = acos(dot);
            return (quatA * sin(angle * (1 - p)) + q * sin(angle * p)) / sin(angle);
        }
};

inline Quaternion operator * (double n, const Quaternion & quat) {
    return quat * n;
}