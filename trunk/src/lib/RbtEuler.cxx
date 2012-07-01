/*This file is part of Rdock.

    Rdock is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Rdock is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Rdock.  If not, see <http://www.gnu.org/licenses/>.*/

#include "RbtEuler.h"

RbtQuat RbtEuler::ToQuat() const {
    RbtDouble c1 = cos(m_heading/2.0);
    RbtDouble s1 = sin(m_heading/2.0);
    RbtDouble c2 = cos(m_attitude/2.0);
    RbtDouble s2 = sin(m_attitude/2.0);
    RbtDouble c3 = cos(m_bank/2.0);
    RbtDouble s3 = sin(m_bank/2.0);
    RbtDouble c1c2 = c1*c2;
    RbtDouble s1s2 = s1*s2;
    RbtQuat q(c1c2*c3 - s1s2*s3,
                c1c2*s3 + s1s2*c3,
                c1*s2*c3 - s1*c2*s3,
                s1*c2*c3 + c1*s2*s3);
    return q;
}

void RbtEuler::FromQuat(const RbtQuat& q) {
    RbtDouble test = (q.v.x * q.v.z) + (q.v.y * q.s);
    if (test > 0.499999) { // singularity at north pole
        m_heading = 2.0 * atan2(q.v.x, q.s);
        m_attitude = M_PI/2.0;
        m_bank = 0.0;
    }
    else if (test < -0.499999) { // singularity at south pole
        m_heading = -2.0 * atan2(q.v.x, q.s);
        m_attitude = -M_PI/2.0;
        m_bank = 0.0;
    }
    else {
        RbtVector sq = q.v * q.v;
        m_heading = atan2(2.0 * (q.v.z * q.s - q.v.x * q.v.y),
                            1.0 - 2.0 * (sq.z + sq.y));
        m_attitude = asin(2.0 * test);
        m_bank = atan2(2.0 * (q.v.x * q.s - q.v.z * q.v.y),
                            1.0 - 2.0 * (sq.x + sq.y));
    }
}
