#ifndef VEC3_H
#define VEC3_H
#include "vec3.h"
#endif 

#ifndef BOUNDINGBOX_H
#define BOUNDINGBOX_H

#ifndef K_INFINITY
#define K_INFINITY
const float kInfinity = std::numeric_limits<float>::max();
#endif

/*
This class represents a bounding box used to accelerate ray-tracing
calculations.

When calculating ray tracing intersections, if the ray does not
intersect the box, then the ray will not intersect any of the geometry 
within the box. This allows for an initial faster and easier intersection
test than if checking other geometry.
*/
template<typename T = float>
class BBox
{
public:

    //Class constructors
    BBox() {}
    //Constructor using 2 3D coordinates to determine the bounding box
    BBox(Vec3<T> min_, Vec3<T> max_)
    {
        bounds[0] = min_;
        bounds[1] = max_;
    }

    //Use when increasing size of bounding box
    BBox& extendBy(const Vec3<T>& p)
    {
        if (p.x < bounds[0].x) bounds[0].x = p.x;
        if (p.y < bounds[0].y) bounds[0].y = p.y;
        if (p.z < bounds[0].z) bounds[0].z = p.z;
        if (p.x > bounds[1].x) bounds[1].x = p.x;
        if (p.y > bounds[1].y) bounds[1].y = p.y;
        if (p.z > bounds[1].z) bounds[1].z = p.z;

        return *this;
    }

    //Calculate the centroid of the box
    Vec3<T> centroid() const { return (bounds[0] + bounds[1]) * 0.5; }

    Vec3<T>& operator [] (bool i) { return bounds[i]; }
    const Vec3<T> operator [] (bool i) const { return bounds[i]; }

    //Calculate the intersection of a ray with the bounding box
    //@param Vec3<T> origin
    //@param Vec3<T> direction
    //@param Vec3<T> sign
    //@param float hit
    bool intersect(const Vec3<T>&, const Vec3<T>&, const Vec3<bool>&, float&) const;
    
    //Initial bounding box boundaries
    Vec3<T> bounds[2] = { kInfinity, -kInfinity };
};

template<typename T>
bool BBox<T>::intersect(const Vec3<T>& orig, const Vec3<T>& invDir, const Vec3<bool>& sign, float& tHit) const
{
    float tmin, tmax, tymin, tymax, tzmin, tzmax;

    tmin = (bounds[sign[0]].x - orig.x) * invDir.x;
    tmax = (bounds[1 - sign[0]].x - orig.x) * invDir.x;
    tymin = (bounds[sign[1]].y - orig.y) * invDir.y;
    tymax = (bounds[1 - sign[1]].y - orig.y) * invDir.y;

    if ((tmin > tymax) || (tymin > tmax))
        return false;

    if (tymin > tmin)
        tmin = tymin;
    if (tymax < tmax)
        tmax = tymax;

    tzmin = (bounds[sign[2]].z - orig.z) * invDir.z;
    tzmax = (bounds[1 - sign[2]].z - orig.z) * invDir.z;

    if ((tmin > tzmax) || (tzmin > tmax))
        return false;

    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;

    tHit = tmin;

    return true;
}
#endif // !BOUNDINGBOX_H