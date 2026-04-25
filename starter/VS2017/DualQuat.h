#ifndef DUALQUAT_H
#define DUALQUAT_H

#include "transform4d.h"
#include <vector>
#include <string>


using namespace std;

class DualQuat
{
public:
	Vec4d real; // (x, y, z, w)
	Vec4d dual; // (x, y, z, w)

    static DualQuat makeDQ(const Mat4d& M)
    {
        Vec4d qr = rotationMatrixToQuaternion(M);

        Vec4d t(M[0][3], M[1][3], M[2][3], 0.0);
        Vec4d qd = 0.5 * quatMultiply(t, qr);

        return { qr, qd };
    }

    static Vec4d rotationMatrixToQuaternion(const Mat4d& M)
    {
        double trace = M[0][0] + M[1][1] + M[2][2];
        Vec4d q;

        if (trace > 0.0)
        {
            double s = sqrt(trace + 1.0) * 2.0;
            q[3] = 0.25 * s;
            q[0] = (M[2][1] - M[1][2]) / s;
            q[1] = (M[0][2] - M[2][0]) / s;
            q[2] = (M[1][0] - M[0][1]) / s;
        }
        else
        {
            if (M[0][0] > M[1][1] && M[0][0] > M[2][2])
            {
                double s = sqrt(1.0 + M[0][0] - M[1][1] - M[2][2]) * 2.0;
                q[3] = (M[2][1] - M[1][2]) / s;
                q[0] = 0.25 * s;
                q[1] = (M[0][1] + M[1][0]) / s;
                q[2] = (M[0][2] + M[2][0]) / s;
            }
            else if (M[1][1] > M[2][2])
            {
                double s = sqrt(1.0 + M[1][1] - M[0][0] - M[2][2]) * 2.0;
                q[3] = (M[0][2] - M[2][0]) / s;
                q[0] = (M[0][1] + M[1][0]) / s;
                q[1] = 0.25 * s;
                q[2] = (M[1][2] + M[2][1]) / s;
            }
            else
            {
                double s = sqrt(1.0 + M[2][2] - M[0][0] - M[1][1]) * 2.0;
                q[3] = (M[1][0] - M[0][1]) / s;
                q[0] = (M[0][2] + M[2][0]) / s;
                q[1] = (M[1][2] + M[2][1]) / s;
                q[2] = 0.25 * s;
            }
        }
        return q; // (x, y, z, w)
    }

    static Vec4d quatMultiply(const Vec4d& a, const Vec4d& b)
    {
        return Vec4d(
            a[3] * b[0] + a[0] * b[3] + a[1] * b[2] - a[2] * b[1],
            a[3] * b[1] - a[0] * b[2] + a[1] * b[3] + a[2] * b[0],
            a[3] * b[2] + a[0] * b[1] - a[1] * b[0] + a[2] * b[3],
            a[3] * b[3] - a[0] * b[0] - a[1] * b[1] - a[2] * b[2]
        );
    }

    static Vec4d quatConjugate(const Vec4d& q)
    {
        return Vec4d(-q[0], -q[1], -q[2], q[3]);
    }

    static Vec3d transformPoint(const DualQuat& dq, const Vec3d& p)
    {
        Vec4d pQuat(p[0], p[1], p[2], 0.0);

        Vec4d qr = dq.real;
        Vec4d qd = dq.dual;

        Vec4d qr_conj = quatConjugate(qr);

        // rotated point
        Vec4d rotated = quatMultiply(quatMultiply(qr, pQuat), qr_conj);

        // translation
        Vec4d trans = quatMultiply(qd, qr_conj);
        trans = 2.0 * trans;

        return Vec3d(
            rotated[0] + trans[0],
            rotated[1] + trans[1],
            rotated[2] + trans[2]
        );
    }
};

#endif