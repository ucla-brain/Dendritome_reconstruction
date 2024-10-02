//
// Created by muyezhu on 3/3/19.
//

#ifndef MCP3D_VAA3D_SD_XYZ_HPP
#define MCP3D_VAA3D_SD_XYZ_HPP

#include <cmath>

struct XYZ
{
    float x, y, z;
    XYZ(float _x = 0.0, float _y = 0.0, float _z = 0.0):
            x(_x), y(_y), z(_z) {};
};

inline bool operator == (const XYZ& a, const XYZ& b)
{
    return (a.x==b.x && a.y==b.y && a.z==b.z);
}

inline XYZ operator - (const XYZ& a, const XYZ& b)
{
    XYZ c;	c.x = a.x-b.x;	c.y = a.y-b.y;	c.z = a.z-b.z;	return c;
}

inline XYZ operator * (const XYZ& a, const XYZ& b)
{
    XYZ c;	c.x = a.x*b.x;	c.y = a.y*b.y;	c.z = a.z*b.z;	return c;
}

inline float dot(const XYZ& a, const XYZ& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline float norm(const XYZ& a)
{
    return sqrt(dot(a,a));
}

inline float dist_L2(const XYZ& a, const XYZ& b)
{
    XYZ c(a.x - b.x, a.y - b.y, a.z - b.z);
    return sqrt(dot(c,c));
}

#endif //MCP3D_VAA3D_SD_XYZ_HPP
