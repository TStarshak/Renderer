#ifndef VEC3_H
#define VEC3_H

#include <iostream>

//Class representing a 3 dimensional vector
template<typename T>
class Vec3
{

public:
    //Default constructor 
    Vec3() : x(0), y(0), z(0) {}

    //Custom constructor for equal components
    Vec3(T xx) : x(xx), y(xx), z(xx) {}

    //Custom constructor for different components
    Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}

    Vec3 operator * (const T& r) const { return Vec3(x * r, y * r, z * r); }
    Vec3 operator + (const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator - (const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    
    template<typename U>
    Vec3 operator / (const Vec3<U>& v) const { return Vec3(x / v.x, y / v.y, z / v.z); }
    friend Vec3 operator / (const T r, const Vec3& v)
    {
        return Vec3(r / v.x, r / v.y, r / v.z);
    }
    
    const T& operator [] (size_t i) const { return (&x)[i]; }
    T& operator [] (size_t i) { return (&x)[i]; }
    
    //Calculate the length squared of the vector
    T length2() const { return x * x + y * y + z * z; }

    friend Vec3 operator * (const float& r, const Vec3& v)
    {
        return Vec3(v.x * r, v.y * r, v.z * r);
    }

    friend std::ostream& operator << (std::ostream& os, const Vec3<T>& v)
    {
        os << v.x << " " << v.y << " " << v.z << std::endl; return os;
    }

    //The components of the vector
    T x, y, z;
};

//Calculate the crossproduct of two vectors
template<typename T>
Vec3<T> cross(const Vec3<T>& a, const Vec3<T>& b)
{
    return Vec3<T>(a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x);
}

//Calculate the dot product of two vectors
template<typename T>
T dot(const Vec3<T>& va, const Vec3<T>& vb)
{
    return va.x * vb.x + va.y * vb.y + va.z * vb.z;
}

//Normalize a vector into its normal vector
template<typename T>
void normalize(Vec3<T>& vec)
{
    T len2 = vec.length2();
    if (len2 > 0) {
        T invLen = 1 / sqrt(len2);
        vec.x *= invLen, vec.y *= invLen, vec.z *= invLen;
    }
}

#endif // !VEC3_H