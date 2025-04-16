// Adapted from RTW and PBRT
#include "pbrt/shapes.h"
#include "pbrt/primitive.h"
namespace pbrt {
// Sphere
std::optional<ShapeIntersection> Sphere::Intersect(
    const Ray &ray, Float tMax
) const
{
  Vector3f oc = ray.origin - center;

  Float a = Dot(ray.direction, ray.direction);
  Float b = Dot(oc, ray.direction);
  Float c = Dot(oc, oc) - radius * radius;
  Float discriminant = b * b - a * c;

  if (discriminant > 0) {
    Float temp = (-b - sqrt(b * b - a * c)) / a;
    if (temp < tMax && temp > tMin) {
      // Surface Interaction contains
      //  point and normal
      SurfaceInteraction surfIntrc;
      surfIntrc.point = ray(temp);
      surfIntrc.normal =
          Normal3f((surfIntrc.point - center) / radius);

      // Shape Intersection additionally contains
      //  tHit
      ShapeIntersection shapeIntrs;
      shapeIntrs.tHit = temp;
      shapeIntrs.interaction = surfIntrc;

      return shapeIntrs;
    }

    temp = (-b + sqrt(b * b - a * c)) / a;
    if (temp < tMax && temp > tMin) {
      SurfaceInteraction surfIntrc;
      surfIntrc.point = ray(temp);
      surfIntrc.normal =
          Normal3f((surfIntrc.point - center) / radius);

      ShapeIntersection shapeIntrs;
      shapeIntrs.tHit = temp;
      shapeIntrs.interaction = surfIntrc;

      return shapeIntrs;
    }
  }

  return {};
}

inline bool Sphere::IntersectP(const Ray &ray, Float tMax) const
{
  return bool(Intersect(ray, tMax));
}


// Triangle
std::vector<const TriangleMesh *> *Triangle::allMeshes;

Bounds3f Triangle::Bounds() const
{
  Point3f p0, p1, p2;
  LoadVertices(&p0, &p1, &p2);
  return Union(Bounds3f(p0, p1), p2);
}

void Triangle::CreateTriangles(
    const TriangleMesh *mesh,
    std::vector<std::unique_ptr<Shape>> *triangles
)
{
  static std::mutex allMeshLock;
  allMeshLock.lock();
  int meshIndex = int(allMeshes->size());
  allMeshes->push_back(mesh);
  allMeshLock.unlock();

  for (int i = 0; i < mesh->nTriangles; ++i) {
    triangles->push_back(
        std::make_unique<Triangle>(meshIndex, i)
    );
  }
}

void Triangle::CreateTriangles(
    const TriangleMesh *mesh,
    std::vector<Primitive *> *primitives,
    Material *material
)
{
  static std::mutex allMeshLock;
  allMeshLock.lock();
  int meshIndex = int(allMeshes->size());
  allMeshes->push_back(mesh);
  allMeshLock.unlock();

  for (int i = 0; i < mesh->nTriangles; ++i) {
    primitives->push_back(new GeometricPrimitive(
        new Triangle(meshIndex, i), material
    ));
  }
}

// intentially not a member function
std::optional<TriangleIntersection> IntersectTriangle(
    const Ray &ray, Float tMax, Point3f p0, Point3f p1,
    Point3f p2
)
{
  // Return no intersection if triangle is degenerate
  if (LengthSquared(Cross(p2 - p0, p1 - p0)) == 0) return {};

  // Transform triangle vertices to ray coordinate space
  // (1) Translate points based on ray's origin
  Point3f p0t = p0 - Vector3f(ray.origin);
  Point3f p1t = p1 - Vector3f(ray.origin);
  Point3f p2t = p2 - Vector3f(ray.origin);
  // (2) Permute components
  int kz = MaxComponentIndex(Abs(ray.direction));
  int kx = kz + 1;
  if (kx == 3) kx = 0;
  int ky = kx + 1;
  if (ky == 3) ky = 0;
  Vector3f d = Permute(ray.direction, {kx, ky, kz});
  p0t = Permute(p0t, {kx, ky, kz});
  p1t = Permute(p1t, {kx, ky, kz});
  p2t = Permute(p2t, {kx, ky, kz});
  // (3) Apply sheer so ray points in +z dir
  Float Sx = -d.x / d.z;
  Float Sy = -d.y / d.z;
  Float Sz = 1 / d.z;
  p0t.x += Sx * p0t.z;
  p0t.y += Sy * p0t.z;
  p1t.x += Sx * p1t.z;
  p1t.y += Sy * p1t.z;
  p2t.x += Sx * p2t.z;
  p2t.y += Sy * p2t.z;

  // Compute edge function coefficients e0, e1, and e2
  Float e0 = DifferenceOfProducts(p1t.x, p2t.y, p1t.y, p2t.x);
  Float e1 = DifferenceOfProducts(p2t.x, p0t.y, p2t.y, p0t.x);
  Float e2 = DifferenceOfProducts(p0t.x, p1t.y, p0t.y, p1t.x);

  // Fall back to double-precision test at triangle edges TODO
  // Perform triangle edge and determinant tests

  if ((e0 < 0 || e1 < 0 || e2 < 0) &&
      (e0 > 0 || e1 > 0 || e2 > 0))
    return {};
  Float det = e0 + e1 + e2;
  if (det == 0) return {};

  // Compute scaled hit distance to triangle and test again ray t range
  p0t.z *= Sz;
  p1t.z *= Sz;
  p2t.z *= Sz;
  Float tScaled = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
  if (det < 0 && (tScaled >= 0 || tScaled < tMax * det))
    return {};
  else if (det > 0 && (tScaled <= 0 || tScaled > tMax * det))
    return {};

  Float invDet = 1 / det;
  Float b0 = e0 * invDet;
  Float b1 = e1 * invDet;
  Float b2 = e2 * invDet;
  Float t = tScaled * invDet;

  // Ensure that computed triangle t is conservatively greater
  // than zero TODO Return TriangleIntersection for intersection
  return TriangleIntersection {b0, b1, b2, t};
}

std::optional<ShapeIntersection> Triangle::Intersect(
    const Ray &ray, Float tMax
) const
{
  const TriangleMesh *mesh = GetMesh();
  Point3f p0, p1, p2;
  LoadVertices(mesh, triIndex, &p0, &p1, &p2);

  std::optional<TriangleIntersection> ti =
      IntersectTriangle(ray, tMax, p0, p1, p2);

  if (!ti) return {};
  SurfaceInteraction inter = InteractionFromIntersection(
      mesh, triIndex, *ti, ray.time, -ray.direction
  );
  return ShapeIntersection {inter, ti->t};
}

bool Triangle::IntersectP(const Ray &ray, Float tMax) const
{
  return bool(Intersect(ray, tMax));
}

void Triangle::Init()
{
  allMeshes = new std::vector<const TriangleMesh *>;
}
}  // namespace pbrt
