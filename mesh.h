#ifndef VEC3_H
#define VEC3_H
#include "vec3.h"
#endif 

#ifndef MATRIX_H
#define MATRIX_H
#include "vec3.h"
#endif 

#ifndef BOUNDINGBOX_H
#define BOUNDINGBOX_H
#include "vec3.h"
#endif 

#ifndef MESH_H
#define MESH_H

/**
This class represents abstract geometry in 3 dimensional space
*/
class GeometricPrimitive
{
public:
    GeometricPrimitive(Matrix44<float> objectToWorld_);
    virtual ~GeometricPrimitive();
    virtual bool intersect(const Vec3<float>&, const Vec3<float>&, float&) const = 0;
    Matrix44<float> objectToWorld;
    Matrix44<float> worldToObject;
    BBox<> bbox;
    uint32_t test;
};

/**
This class represents a polygon mesh 
*/
class Mesh : public GeometricPrimitive
{
public:
    Mesh(
        uint32_t numPolygons,
        const std::vector<uint32_t>& polygonNumVertsArray,          // how many vertices are making each face (array size is num polys)
        const std::vector<uint32_t>& polygonIndicesInVertexPool,    // the index of the vertices making each face (array size is the sum of each number in prev array)
        std::vector<Vec3<float>> vertexPool_,                             // the vertex positions (should be at least as many as max index in previous array)
        Matrix44<float> objectToWorld_ = Matrix44<float>());

    bool intersect(const Vec3<float>& rayOrig, const Vec3<float>& rayDir, float& tNear) const;
    uint32_t numTriangles = { 0 };
    std::vector<uint32_t> triangleIndicesInVertexPool;
    std::vector<Vec3<float>> vertexPool;

};

//Calculates if there is an intersection between a ray and a triangle
//@param Vec3& origin
//@param Vec3& direction
//@param Vec3& vertex 1
//@param Vec3& vertex 2
//@param Vec4& vertex 3
//@param float& intersection x component
//@param float& intersection y component
//@param float& intersection z component
bool rayTriangleIntersect(
    const Vec3<float>& orig, const Vec3<float>& dir,
    const Vec3<float>& v0, const Vec3<float>& v1, const Vec3<float>& v2,
    float& t, float& u, float& v);


//This class represents abstract acceleration structures for ray tracing
class AccelerationStructure
{
public:
    AccelerationStructure(std::vector<std::unique_ptr<const Mesh>>& m);
    virtual ~AccelerationStructure();
    virtual bool intersect(const Vec3<float>& orig, const Vec3<float>& dir, const uint32_t& rayId, float& tHit) const;

    const std::vector<std::unique_ptr<const Mesh>> meshes;
};

//This class represents a bounding volume hierarchy
//This classes uses an octree as its underlying structure
class BVH : public AccelerationStructure
{
    static const uint8_t kNumPlaneSetNormals = 7;
    static const Vec3<float> planeSetNormals[kNumPlaneSetNormals];

    //This struct represents a tighter boundary than the boundingbox class
    struct Extents
    {
        Extents();
        void extendBy(const Extents& e);
        Vec3<float> centroid() const;
        bool intersect(const float*, const float*, float&, float&, uint8_t&) const;
        float d[kNumPlaneSetNormals][2];
        const Mesh* mesh;
    };

    //This struct represents an octree, a representation of space
    //Each node on this tree has 8 children, representing a split of 8 cubes 
    struct Octree
    {
        //Constructor and destructor
        Octree(const Extents& sceneExtents);
        ~Octree();

        //Insert extents into the tree
        void insert(const Extents* extents);
        void build();

        //This struct represents nodes of the tree
        struct OctreeNode 
        {
            OctreeNode* child[8] = { nullptr };
            //Pointer to the extends of this node
            std::vector<const Extents*> nodeExtentsList; 
            //Actual extents of this node
            Extents nodeExtents;
            bool isLeaf = true;
        };

        //This struct
        struct QueueElement
        {
            const OctreeNode* node; // octree node held by this element in the queue
            float t; // distance from the ray origin to the extents of the node
            QueueElement(const OctreeNode* n, float tn);
            friend bool operator < (const QueueElement& a, const QueueElement& b) { return a.t > b.t; };

        };

        OctreeNode* root = nullptr; 
        BBox<> bbox;

     private:
         void deleteOctreeNode(OctreeNode*& node);
         void insert(OctreeNode*& node, const Extents* extents, const BBox<>& bbox, uint32_t depth);
         void build(OctreeNode*& node, const BBox<>& bbox);
     };

    std::vector<Extents> extentsList;
    Octree* octree = nullptr;

public:
    BVH(std::vector<std::unique_ptr<const Mesh>>& m);
    bool intersect(const Vec3<float>&, const Vec3<float>&, const uint32_t&, float&) const;
    ~BVH() { delete octree; }
};

#endif 