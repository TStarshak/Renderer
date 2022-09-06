#ifndef VEC3_H
#define VEC3_H
#include "vec3.h"
#endif 

#ifndef MATRIX44_H
#define MATRIX44_H


/**
This class represents a 4x4 matrix
*/
template<typename T>
class Matrix44
{
public:
    //Default constructor using the identity matrix
    Matrix44() {}
    //Custom constructor
    Matrix44(T m00, T m01, T m02, T m03,
        T m10, T m11, T m12, T m13,
        T m20, T m21, T m22, T m23,
        T m30, T m31, T m32, T m33)
    {
        m[0][0] = m00; m[0][1] = m01; m[0][2] = m02; m[0][3] = m03;
        m[1][0] = m10; m[1][1] = m11; m[1][2] = m12; m[1][3] = m13;
        m[2][0] = m20; m[2][1] = m21; m[2][2] = m22; m[2][3] = m23;
        m[3][0] = m30; m[3][1] = m31; m[3][2] = m32; m[3][3] = m33;
    }
    //Calculate the inverse of the matrix
    Matrix44 inverse() const { Matrix44 matInv = *this; return matInv; }

    T* operator [] (size_t i) { return &m[i][0]; }
    const T* operator [] (size_t i) const { return &m[i][0]; }

    //The matrix itself, defaults to be the identity matrix
    T m[4][4] = { {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1} };
};

//Multiply the matrix by a vector
//@param Matrix44& matrix
//@param Vec3& vector
template<typename T>
void matVecMult(const Matrix44<T>& m, Vec3<T>& v)
{
    Vec3<T> vt;
    vt.x = v.x * m[0][0] + v.y * m[1][0] + v.z * m[2][0],
        vt.y = v.x * m[0][1] + v.y * m[1][1] + v.z * m[2][1],
        vt.z = v.x * m[0][2] + v.y * m[1][2] + v.z * m[2][2];

    v = vt;
}

//Multiply the matrix by a point
//@param Matrix44& matrix
//@param Vec3& point
template<typename T>
void matPointMult(const Matrix44<T>& m, Vec3<T>& p)
{
    Vec3<T> pt;
    pt.x = p.x * m[0][0] + p.y * m[1][0] + p.z * m[2][0] + m[3][0];
    pt.y = p.x * m[0][1] + p.y * m[1][1] + p.z * m[2][1] + m[3][1];
    pt.z = p.x * m[0][2] + p.y * m[1][2] + p.z * m[2][2] + m[3][2];
    T w = p.x * m[0][3] + p.y * m[1][3] + p.z * m[2][3] + m[3][3];
    if (w != 1) {
        pt.x /= w, pt.y /= w, pt.z /= w;
    }

    p = pt;
}

#endif MATRIX44_H