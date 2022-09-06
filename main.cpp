#include <atomic>
#include <memory>
#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include <cmath>
#include <chrono>

#include "vec3.h"
#include "matrix44.h"
#include "boundingbox.h"
#include "mesh.h"

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

const float kEpsilon = 1e-8;
#ifndef K_INFINITY
#define K_INFINITY
const float kInfinity = std::numeric_limits<float>::max();
#endif

std::atomic<uint32_t> numPrimaryRays(0);

using Vec3f = Vec3<float>;
using Vec3b = Vec3<bool>;
using Vec3i = Vec3<int32_t>;
using Vec3ui = Vec3<uint32_t>;
using Matrix44f = Matrix44<float>;


//Clamps a value between a range of values
template<typename T> inline T clamp(const T& v, const T& lo, const T& hi)
{
    return std::max(lo, std::min(v, hi));
}

#include "modeldata.h"

Vec3<float> evalBezierCurve(const Vec3<float>* P, const float& t)
{
    float b0 = (1 - t) * (1 - t) * (1 - t);
    float b1 = 3 * t * (1 - t) * (1 - t);
    float b2 = 3 * t * t * (1 - t);
    float b3 = t * t * t;

    return P[0] * b0 + P[1] * b1 + P[2] * b2 + P[3] * b3;
};

Vec3<float> evalBezierPatch(const Vec3<float>* controlPoints, const float& u, const float& v)
{
    Vec3<float> uCurve[4];
    for (size_t i = 0; i < 4; ++i)
        uCurve[i] = evalBezierCurve(controlPoints + 4 * i, u);

    return evalBezierCurve(uCurve, v);
};

Vec3<float> derivBezier(const Vec3<float>* P, const float& t)
{
    return -3 * (1 - t) * (1 - t) * P[0] +
        (3 * (1 - t) * (1 - t) - 6 * t * (1 - t)) * P[1] +
        (6 * t * (1 - t) - 3 * t * t) * P[2] +
        3 * t * t * P[3];
};

Vec3<float> dUBezier(const Vec3<float>* controlPoints, float u, float v)
{
    Vec3<float> P[4];
    Vec3<float> vCurve[4];
    for (size_t i = 0; i < 4; ++i) {
        P[0] = controlPoints[i];
        P[1] = controlPoints[4 + i];
        P[2] = controlPoints[8 + i];
        P[3] = controlPoints[12 + i];
        vCurve[i] = evalBezierCurve(P, v);
    }

    return derivBezier(vCurve, u);
};

Vec3<float> dVBezier(const Vec3<float>* controlPoints, float u, float v)
{
    Vec3<float> uCurve[4];
    for (size_t i = 0; i < 4; ++i) {
        uCurve[i] = evalBezierCurve(controlPoints + 4 * i, u);
    }

    return derivBezier(uCurve, v);
};

//Constructs the model supplied in modeldata.h
std::vector<std::unique_ptr<const Mesh>> createModel()
{
    Matrix44f rotate90(1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1);
    std::vector<std::unique_ptr<const Mesh>> meshes;
    uint32_t width = 8, height = 8;
    uint32_t numPolygons = width * height;
    std::vector<uint32_t> polyNumVertsArray(numPolygons, 4);
    std::vector<uint32_t> polyIndicesInVertPool(numPolygons * 4);
    // set indices 
    for (uint32_t y = 0, offset = 0; y < height; ++y) {
        for (uint32_t x = 0; x < width; ++x, offset += 4) {
            // counter-clockwise to get the normal pointing in the right direction
            polyIndicesInVertPool[offset] = (width + 1) * y + x;
            polyIndicesInVertPool[offset + 3] = (width + 1) * y + x + 1;
            polyIndicesInVertPool[offset + 2] = (width + 1) * (y + 1) + x + 1;
            polyIndicesInVertPool[offset + 1] = (width + 1) * (y + 1) + x;
        }
    }
    Vec3f controlPoints[16];
    //Iterate over each patch in the model
    for (uint32_t i = 0; i < kmodelNumPatches; ++i) {
        std::vector<Vec3f> vertPool((width + 1) * (height + 1));
        for (uint32_t j = 0; j < 16; ++j) {
            controlPoints[j].x = modelVertices[modelPatches[i][j] - 1][0],
                controlPoints[j].y = modelVertices[modelPatches[i][j] - 1][1],
                controlPoints[j].z = modelVertices[modelPatches[i][j] - 1][2];
        }
        for (uint32_t y = 0, currVertIndex = 0; y <= height; ++y) {
            float v = y / (float)height;
            for (uint32_t x = 0; x <= width; ++x, ++currVertIndex) {
                float u = x / (float)width;
                vertPool[currVertIndex] = evalBezierPatch(controlPoints, u, v);
                matVecMult(rotate90, vertPool[currVertIndex]);
                Vec3f dU = dUBezier(controlPoints, u, v);
                Vec3f dV = dVBezier(controlPoints, u, v);
                Vec3f N = cross(dU, dV);
            }
        }

        meshes.emplace_back(new Mesh(numPolygons, polyNumVertsArray, polyIndicesInVertPool, vertPool));
    }

    return meshes;
}

//Creates the model and saves it to a variable
//@param meshes
void makeScene(std::vector<std::unique_ptr<const Mesh>>& meshes)
{
    meshes = std::move(createModel());
}

//Converts degrees to radians
template<typename T>
inline
T degToRad(const T& angle) { return angle / 180.f * M_PI; }


//This struct represents options for the image created by the render function
struct Options
{
    float fov = { 90 };
    uint32_t width = { 640 };
    uint32_t height = { 480 };
    Matrix44f cameraToWorld, worldToCamera;
};


/*
Trace a ray from the camera to each pixel, then calculate intersections
*/
void render(const std::unique_ptr<AccelerationStructure>& accel, const Options& options)
{
    std::unique_ptr<Vec3f[]> buffer(new Vec3f[options.width * options.height]);
    Vec3f orig(0, 0, 5);
    matPointMult(options.cameraToWorld, orig);
    float scale = std::tan(degToRad<float>(options.fov * 0.5));
    float imageAspectRatio = options.width / static_cast<float>(options.height);
    assert(imageAspectRatio > 1);
    uint32_t rayId = 1; // Start at 1 not 0!! (see Grid code and mailboxing)
    for (uint32_t j = 0; j < options.height; ++j) {
        for (uint32_t i = 0; i < options.width; ++i) {
            Vec3f dir((2 * (i + 0.5f) / options.width - 1) * scale * imageAspectRatio,
                (1 - 2 * (j + 0.5) / options.height) * scale,
                -1);
            matVecMult(options.cameraToWorld, dir);
            normalize(dir);
            numPrimaryRays++;
            float tHit = kInfinity;
            buffer[j * options.width + i] = (accel->intersect(orig, dir, rayId++, tHit)) ? Vec3f(1) : Vec3f(0);
        }
    }

    // store to PPM file
    std::ofstream ofs;
    ofs.open("image.ppm");
    ofs << "P6\n" << options.width << " " << options.height << "\n255\n";
    for (uint32_t i = 0; i < options.width * options.height; ++i) {
        Vec3<uint8_t> pixRgb;
        pixRgb.x = static_cast<uint8_t>(255 * std::max(0.f, std::min(1.f, buffer[i].x)));
        pixRgb.y = static_cast<uint8_t>(255 * std::max(0.f, std::min(1.f, buffer[i].y)));
        pixRgb.z = static_cast<uint8_t>(255 * std::max(0.f, std::min(1.f, buffer[i].z)));
        ofs << pixRgb.x << pixRgb.y << pixRgb.z;
    }
    ofs.close();
}


int main(int argc, char** argv)
{
    std::vector<std::unique_ptr<const Mesh>> meshes;
    makeScene(meshes);

    std::unique_ptr<AccelerationStructure> accel(new BVH(meshes));


    uint32_t numTriangles{ 0 };
    for (const auto& mesh : accel->meshes) {
        numTriangles += mesh->numTriangles;
    }

    Options options;
    using Time = std::chrono::high_resolution_clock;
    using fsec = std::chrono::duration<float>;

    auto t0 = Time::now();

    render(accel, options);

    auto t1 = Time::now();

    fsec fs = t1 - t0;
    std::cout << "Render time                                 | " << fs.count() << " sec" << std::endl;
    std::cout << "Total number of triangles                   | " << numTriangles << std::endl;
    std::cout << "Total number of primary rays                | " << numPrimaryRays << std::endl;

    return 0;
}