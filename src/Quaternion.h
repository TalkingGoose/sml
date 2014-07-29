#pragma once

/********************************************************************
    Author:         Paul Watkinson
    Website:        http://paulwatkinson.co.uk/

    Description:    A quaternion class
*********************************************************************/

#include "Vector3.h"

class Quaternion {
    private:
        // http://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
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

        Quaternion(void) {
            this->x = this->y = this->z = 0.0f;
            this->w = 1.0f;
        }

        Quaternion(double a, double b, double c) {
            FromEuler(a, b, c);
        }

        Quaternion(double a, double b, double c, double d) {
            this->x = a; this->y = b; this->z = c; this->w = d;
        }
        
        Quaternion(Vector3d vector) {
            FromEuler(vector.x, vector.y, vector.z);
        }

        Quaternion(Vector3d vector, double d) {
            this->x = vector.x; this->y = vector.y; this->z = vector.z; this->w = d;
        }

        virtual ~Quaternion(void) { };

        Quaternion & Set(double a, double b, double c) {
            FromEuler(a, b, c);

            return (*this);
        }

        Quaternion & Set(double a, double b, double c, double d) {
            this->x = a;
            this->y = b;
            this->z = c;
            this->w = d;

            return (*this);
        }

        Quaternion & Set(const Vector3d vector) {
            FromEuler(vector.x, vector.y, vector.z);

            return (*this);
        }

        Quaternion & Set(const Vector3d vector, double d) {
            this->x = vector.x;
            this->y = vector.y;
            this->z = vector.z;
            this->w = d;

            return (*this);
        }

        Quaternion & Set(const Quaternion & quat) {
            this->x = quat.x;
            this->y = quat.y;
            this->z = quat.z;
            this->w = quat.w;

            return (*this);
        }

        Quaternion & operator =  (const Quaternion & quat) {
            this->x = quat.x;
            this->y = quat.y;
            this->z = quat.z;
            this->w = quat.w;

            return (*this);
        }

        Quaternion & operator += (const Quaternion & quat) {
            this->x += quat.x;
            this->y += quat.y;
            this->z += quat.z;
            this->w += quat.w;

            return (*this);
        }

        Quaternion & operator -= (const Quaternion & quat) {
            this->x -= quat.x;
            this->y -= quat.y;
            this->z -= quat.z;
            this->w -= quat.w;

            return (*this);
        }
        
        Quaternion & operator *= (const Quaternion & quat) {
            return (*this) = ((*this) * quat);
        }

        Quaternion & operator *= (double n) {
            this->x *= n;
            this->y *= n;
            this->z *= n;
            this->w *= n;

            return (*this);
        }

        Quaternion & operator /= (const Quaternion & quat) {
            return (*this) = ((*this) / quat);
        }

        Quaternion & operator /= (double n) {
            this->x /= n;
            this->y /= n;
            this->z /= n;
            this->w /= n;

            return (*this);
        }

        Quaternion operator + (const Quaternion & quat) const {
            return Quaternion(this->x + quat.x, this->y + quat.y, this->z + quat.z, this->w + quat.w);
        }

        Quaternion operator - (const Quaternion & quat) const {
            return Quaternion(this->x - quat.x, this->y - quat.y, this->z - quat.z, this->w - quat.w);
        }
        Quaternion operator - (void) const {
            return Quaternion(-this->x, -this->y, -this->z, -this->w);
        }

        Quaternion operator * (const Quaternion & quat) const {
            return Quaternion(
                (w * quat.x) + (x * quat.w) + (y * quat.z) - (z * quat.y),
                (w * quat.y) + (y * quat.w) + (z * quat.x) - (x * quat.z),
                (w * quat.z) + (z * quat.w) + (x * quat.y) - (y * quat.x),
                (w * quat.w) - (x * quat.x) - (y * quat.y) - (z * quat.z)
            );
        }

        Quaternion operator * (double n) const {
            return Quaternion(this->x * n, this->y * n, this->z * n, this->w * n);
        }

        Quaternion operator / (const Quaternion & quat) const {
            return Quaternion(this->Inverted() * quat);
        }

        Quaternion operator / (double n) const {
            return Quaternion(this->x / n, this->y / n, this->z / n, this->w / n);
        }

        bool operator ! (void) const {
            return (this->x != 0 && this->y != 0 && this->z != 0 && this->w != 0);
        }

        bool operator == (const Quaternion & quat) const {
            return (this->x == quat.x && this->y == quat.y && this->z == quat.z && this->w == quat.w);
        }

        bool operator != (const Quaternion & quat) const {
            return (this->x != quat.x && this->y != quat.y && this->z != quat.z && this->w != quat.w);
        }

        double & operator [] (unsigned int n) {
            if (n > 3) { throw "Out of range."; }

            return *(&x + n);
        }

        Vector3d GetAxis(void) const {
            return Vector3d(this->x, this->y, this->z);
        }

        double GetAngle(void) const {
            return acos(w) * 2;
        }

        double SquaredMagnitude(void) const {
            return pow(this->x, 2.0f) + pow(this->y, 2.0f) + pow(this->z, 2.0f) + pow(this->w, 2.0f);
        }

        double Magnitude(void) const {
            return sqrt(this->SquaredMagnitude());
        }

        Quaternion & Normalise(void) {
            return (*this /= Magnitude());
        }

        Quaternion Normalised(void) const {
            double m = Magnitude();
            return Quaternion(this->x / m, this->y / m, this->z / m, this->w / m);
        }

        Quaternion & Invert(void) {
            (*this).Set(this->Inverted());
            return (*this);
        }

        Quaternion Inverted(void) const {
            return this->Conjugate() / this->SquaredMagnitude();
        }

        Quaternion Conjugate(void) const {
            return Quaternion(-this->x, -this->y, -this->z, this->w);
        }

        Quaternion Reciprocal(void) const {
            return (this->Conjugate() / ((*this) * this->Conjugate()));
        }

        double DotProduct(const Quaternion & quat) const {
            return (this->x * quat.x) + (this->y * quat.y) + (this->z * quat.z) + (this->w * quat.w);
        }

        static Quaternion Lerp(const Quaternion & quatA, const Quaternion & quatB, double p) {
            return ((quatA * (1 - p)) + (quatB * p)).Normalised();
        }

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