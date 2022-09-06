#include <vector>
#include <cassert>
#include <queue>

#include "vec3.h"
#include "matrix44.h"
#include "boundingbox.h"
#include "mesh.h"
#include <chrono>

const float kEpsilon = 1e-8;


GeometricPrimitive::GeometricPrimitive(Matrix44<float> objectToWorld_) :
    objectToWorld(objectToWorld_), worldToObject(objectToWorld.inverse()){}

GeometricPrimitive::~GeometricPrimitive() {}



Mesh::Mesh(
    uint32_t numPolygons,
    const std::vector<uint32_t>& polygonNumVertsArray,          // how many vertices are making each face (array size is num polys)
    const std::vector<uint32_t>& polygonIndicesInVertexPool,    // the index of the vertices making each face (array size is the sum of each number in prev array)
    std::vector<Vec3<float>> vertexPool_,                             // the vertex positions (should be at least as many as max index in previous array)
    Matrix44<float> objectToWorld_) :
    GeometricPrimitive(objectToWorld_),
    vertexPool(vertexPool_)
{
    // pass by value (move constructor shouldn't even be called here ?
    for (uint32_t i = 0; i < vertexPool.size(); ++i) {
        matPointMult(objectToWorld, vertexPool[i]);
        bbox.extendBy(vertexPool[i]);
    }
    // compute total number of triangles
    for (uint32_t i = 0; i < numPolygons; ++i) {
        assert(polygonNumVertsArray[i] >= 3);
        numTriangles += polygonNumVertsArray[i] - 2;
    }
    // create array to store the triangle indices in the vertex pool
    // !! use resize() here and not reserve() -- which only affects capacity but doesn't change size of the vector
    triangleIndicesInVertexPool.resize(numTriangles * 3);
    // for each face
    for (uint32_t i = 0, offset = 0, currTriangleIndex = 0; i < numPolygons; ++i) {
        // for each triangle in the face
        for (uint32_t j = 0; j < polygonNumVertsArray[i] - 2; ++j) {
            triangleIndicesInVertexPool[currTriangleIndex] = polygonIndicesInVertexPool[offset];
            triangleIndicesInVertexPool[currTriangleIndex + 1] = polygonIndicesInVertexPool[offset + j + 1];
            triangleIndicesInVertexPool[currTriangleIndex + 2] = polygonIndicesInVertexPool[offset + j + 2];
            currTriangleIndex += 3;
        }
        offset += polygonNumVertsArray[i];
    }
};


bool rayTriangleIntersect(
    const Vec3<float>& orig, const Vec3<float>& dir,
    const Vec3<float>& v0, const Vec3<float>& v1, const Vec3<float>& v2,
    float& t, float& u, float& v)
{

    Vec3<float> v0v1 = v1 - v0;
    Vec3<float> v0v2 = v2 - v0;
    Vec3<float> pvec = cross(dir, v0v2);
    float det = dot(v0v1, pvec);

    // ray and triangle are parallel if determinant is close to 0
    if (fabs(det) < kEpsilon) return false;

    float invDet = 1 / det;

    Vec3<float> tvec = orig - v0;
    u = dot(tvec, pvec) * invDet;
    if (u < 0 || u > 1) return false;

    Vec3<float> qvec = cross(tvec, v0v1);
    v = dot(dir, qvec) * invDet;
    if (v < 0 || u + v > 1) return false;

    t = dot(v0v2, qvec) * invDet;

    return true;
}

bool Mesh::intersect(const Vec3<float>& rayOrig, const Vec3<float>& rayDir, float& tNear) const
{
    // loop over all triangles in the mesh and return true if one 
    // of the triangles at least is intersected
    float t, u, v;
    uint32_t intersectedTriIndex;
    bool intersected = false;
    // tNear should be set inifnity first time this function is called and it 
    // will get eventually smaller as the ray intersects geometry
    for (uint32_t i = 0; i < numTriangles; ++i) {
        if (rayTriangleIntersect(rayOrig, rayDir,
            vertexPool[triangleIndicesInVertexPool[i * 3]],
            vertexPool[triangleIndicesInVertexPool[i * 3 + 1]],
            vertexPool[triangleIndicesInVertexPool[i * 3 + 2]], t, u, v) && t < tNear)
        {
            tNear = t;
            intersectedTriIndex = i;
            intersected = true;
        }
    }

    return intersected;
}


AccelerationStructure::AccelerationStructure(std::vector<std::unique_ptr<const Mesh>>& m) : meshes(std::move(m)) {}
AccelerationStructure::~AccelerationStructure() {}
bool AccelerationStructure::intersect(const Vec3<float>& orig, const Vec3<float>& dir, const uint32_t& rayId, float& tHit) const
{
    
    //Get a pointer to the mesh since we don't want to actually change it
    const Mesh* intersectedMesh = nullptr;
    float t = kInfinity;
    for (const auto& mesh : meshes) {
        if (mesh->intersect(orig, dir, t) && t < tHit) {
            intersectedMesh = mesh.get();
            tHit = t;
        }
    }

    return (intersectedMesh != nullptr);
};



BVH::Extents::Extents()
{
    for (uint8_t i = 0; i < kNumPlaneSetNormals; ++i)
        d[i][0] = kInfinity, d[i][1] = -kInfinity;
};
void BVH::Extents::extendBy(const Extents& e)
{

    for (uint8_t i = 0; i < kNumPlaneSetNormals; ++i) {
        if (e.d[i][0] < d[i][0]) d[i][0] = e.d[i][0];
        if (e.d[i][1] > d[i][1]) d[i][1] = e.d[i][1];
    }
};
        /* inline */
Vec3<float> BVH::Extents::centroid() const
{
    return Vec3<float>(
        d[0][0] + d[0][1] * 0.5,
        d[1][0] + d[1][1] * 0.5,
        d[2][0] + d[2][1] * 0.5);
};


BVH::Octree::Octree(const Extents& sceneExtents)
{
    float xDiff = sceneExtents.d[0][1] - sceneExtents.d[0][0];
    float yDiff = sceneExtents.d[1][1] - sceneExtents.d[1][0];
    float zDiff = sceneExtents.d[2][1] - sceneExtents.d[2][0];
    float maxDiff = std::max(xDiff, std::max(yDiff, zDiff));
    Vec3<float> minPlusMax(
        sceneExtents.d[0][0] + sceneExtents.d[0][1],
        sceneExtents.d[1][0] + sceneExtents.d[1][1],
        sceneExtents.d[2][0] + sceneExtents.d[2][1]);
    bbox[0] = (minPlusMax - maxDiff) * 0.5;
    bbox[1] = (minPlusMax + maxDiff) * 0.5;
    root = new OctreeNode;
};

BVH::Octree::~Octree() { deleteOctreeNode(root); };

void BVH::Octree::insert(const Extents* extents) { insert(root, extents, bbox, 0); };
void BVH::Octree::build() { build(root, bbox); };


BVH::Octree::QueueElement::QueueElement(const OctreeNode* n, float tn) : node(n), t(tn) {};


void BVH::Octree::deleteOctreeNode(OctreeNode*& node)
{
    for (uint8_t i = 0; i < 8; i++) {
        if (node->child[i] != nullptr) {
            deleteOctreeNode(node->child[i]);
        }
    }
    delete node;
};

void BVH::Octree::insert(OctreeNode*& node, const Extents* extents, const BBox<>& bbox, uint32_t depth)
{
    if (node->isLeaf) {
        if (node->nodeExtentsList.size() == 0 || depth == 16) {
            node->nodeExtentsList.push_back(extents);
        }
        else {
            node->isLeaf = false;
            // Re-insert extents held by this node
            while (node->nodeExtentsList.size()) {
                insert(node, node->nodeExtentsList.back(), bbox, depth);
                node->nodeExtentsList.pop_back();
            }
            // Insert new extent
            insert(node, extents, bbox, depth);
        }
    }
    else {
        // compute which child of the current node 
        // we should insert into
        Vec3<float> extentsCentroid = extents->centroid();
        Vec3<float> nodeCentroid = (bbox[0] + bbox[1]) * 0.5;
        BBox<> childBBox;
        uint8_t childIndex = 0;
        // x-axis
        if (extentsCentroid.x > nodeCentroid.x) {
            childIndex = 4;
            childBBox[0].x = nodeCentroid.x;
            childBBox[1].x = bbox[1].x;
        }
        else {
            childBBox[0].x = bbox[0].x;
            childBBox[1].x = nodeCentroid.x;
        }
        // y-axis
        if (extentsCentroid.y > nodeCentroid.y) {
            childIndex += 2;
            childBBox[0].y = nodeCentroid.y;
            childBBox[1].y = bbox[1].y;
        }
        else {
            childBBox[0].y = bbox[0].y;
            childBBox[1].y = nodeCentroid.y;
        }
        // z-axis
        if (extentsCentroid.z > nodeCentroid.z) {
            childIndex += 1;
            childBBox[0].z = nodeCentroid.z;
            childBBox[1].z = bbox[1].z;
        }
        else {
            childBBox[0].z = bbox[0].z;
            childBBox[1].z = nodeCentroid.z;
        }

        // Create the child node if it doesn't exsit yet
        if (node->child[childIndex] == nullptr)
            node->child[childIndex] = new OctreeNode;
        insert(node->child[childIndex], extents, childBBox, depth + 1);
    }
};

void BVH::Octree::build(OctreeNode*& node, const BBox<>& bbox)
{
    if (node->isLeaf) {
        for (const auto& e : node->nodeExtentsList) {
            node->nodeExtents.extendBy(*e);
        }
    }
    else {
        for (uint8_t i = 0; i < 8; ++i) {
            if (node->child[i]) {
                BBox<> childBBox;
                Vec3<float> centroid = bbox.centroid();
                // x-axis
                childBBox[0].x = (i & 4) ? centroid.x : bbox[0].x;
                childBBox[1].x = (i & 4) ? bbox[1].x : centroid.x;
                // y-axis
                childBBox[0].y = (i & 2) ? centroid.y : bbox[0].y;
                childBBox[1].y = (i & 2) ? bbox[1].y : centroid.y;
                // z-axis
                childBBox[0].z = (i & 1) ? centroid.z : bbox[0].z;
                childBBox[1].z = (i & 1) ? bbox[1].z : centroid.z;

                // Inspect child
                build(node->child[i], childBBox);

                // Expand extents with extents of child
                node->nodeExtents.extendBy(node->child[i]->nodeExtents);
            }
        }
    }
};



const Vec3<float> BVH::planeSetNormals[BVH::kNumPlaneSetNormals] = {
    Vec3<float>(1, 0, 0),
    Vec3<float>(0, 1, 0),
    Vec3<float>(0, 0, 1),
    Vec3<float>(sqrtf(3) / 3.f,  sqrtf(3) / 3.f, sqrtf(3) / 3.f),
    Vec3<float>(-sqrtf(3) / 3.f,  sqrtf(3) / 3.f, sqrtf(3) / 3.f),
    Vec3<float>(-sqrtf(3) / 3.f, -sqrtf(3) / 3.f, sqrtf(3) / 3.f),
    Vec3<float>(sqrtf(3) / 3.f, -sqrtf(3) / 3.f, sqrtf(3) / 3.f)
};

BVH::BVH(std::vector<std::unique_ptr<const Mesh>>& m) : AccelerationStructure(m)
{
    //Extent of the entire scene
    Extents sceneExtents;
    extentsList.reserve(meshes.size());
    for (uint32_t i = 0; i < meshes.size(); ++i) {
        for (uint8_t j = 0; j < kNumPlaneSetNormals; ++j) {
            for (const auto vtx : meshes[i]->vertexPool) {
                try {
                    float d = dot(planeSetNormals[j], vtx);
                    // set dNEar and dFar
                    if (d < extentsList[i].d[j][0]) extentsList[i].d[j][0] = d;
                    if (d > extentsList[i].d[j][1]) extentsList[i].d[j][1] = d;
                }
                catch (const std::exception e) {
                    std::cerr << e.what();
                }
            }
        }

        sceneExtents.extendBy(extentsList[i]); 
        // the extent itself needs to keep a pointer to the object its holds
        extentsList[i].mesh = meshes[i].get(); 
    }

    // Now that we have the extent of the scene we can start building our octree
    // Using C++ make_unique function here but you don't need to, just to learn something... 
    octree = new Octree(sceneExtents);

    for (uint32_t i = 0; i < meshes.size(); ++i) {
        octree->insert(&extentsList[i]);
    }

    // Build from bottom up
    octree->build();
}

bool BVH::Extents::intersect(
    const float* precomputedNumerator,
    const float* precomputedDenominator,
    float& tNear,   // tn and tf in this method need to be contained
    float& tFar,    // within the range [tNear:tFar]
    uint8_t& planeIndex) const
{
    for (uint8_t i = 0; i < kNumPlaneSetNormals; ++i) {
        float tNearExtents = (d[i][0] - precomputedNumerator[i]) / precomputedDenominator[i];
        float tFarExtents = (d[i][1] - precomputedNumerator[i]) / precomputedDenominator[i];
        if (precomputedDenominator[i] < 0) std::swap(tNearExtents, tFarExtents);
        if (tNearExtents > tNear) tNear = tNearExtents, planeIndex = i;
        if (tFarExtents < tFar) tFar = tFarExtents;
        if (tNear > tFar) return false;
    }

    return true;
}

bool BVH::intersect(const Vec3<float>& orig, const Vec3<float>& dir, const uint32_t& rayId, float& tHit) const
{
    tHit = kInfinity;
    const Mesh* intersectedMesh = nullptr;
    float precomputedNumerator[BVH::kNumPlaneSetNormals];
    float precomputedDenominator[BVH::kNumPlaneSetNormals];
    for (uint8_t i = 0; i < kNumPlaneSetNormals; ++i) {
        precomputedNumerator[i] = dot(planeSetNormals[i], orig);
        precomputedDenominator[i] = dot(planeSetNormals[i], dir);
    }

    uint8_t planeIndex;
    float tNear = 0, tFar = kInfinity; // tNear, tFar for the intersected extents
    if (!octree->root->nodeExtents.intersect(precomputedNumerator, precomputedDenominator, tNear, tFar, planeIndex) || tFar < 0)
        return false;
    tHit = tFar;
    std::priority_queue<BVH::Octree::QueueElement> queue;
    queue.push(BVH::Octree::QueueElement(octree->root, 0));
    while (!queue.empty() && queue.top().t < tHit) {
        const Octree::OctreeNode* node = queue.top().node;
        queue.pop();
        if (node->isLeaf) {
            for (const auto& e : node->nodeExtentsList) {
                float t = kInfinity;
                if (e->mesh->intersect(orig, dir, t) && t < tHit) {
                    tHit = t;
                    intersectedMesh = e->mesh;
                }
            }
        }
        else {
            for (uint8_t i = 0; i < 8; ++i) {
                if (node->child[i] != nullptr) {
                    float tNearChild = 0, tFarChild = tFar;
                    if (node->child[i]->nodeExtents.intersect(precomputedNumerator, precomputedDenominator, tNearChild, tFarChild, planeIndex)) {
                        float t = (tNearChild < 0 && tFarChild >= 0) ? tFarChild : tNearChild;
                        queue.push(BVH::Octree::QueueElement(node->child[i], t));
                    }
                }
            }
        }
    }

    return (intersectedMesh != nullptr);
}
