#ifndef DNA_REPLICATION_VECTOR3_H
#define DNA_REPLICATION_VECTOR3_H

namespace DNAReplication {

    template<typename T>
    struct Vector3 {
        T x;
        T y;
        T z;

        Vector3<T> &operator+=(Vector3<T> const &rhs) {
            x += rhs.x;
            y += rhs.y;
            z += rhs.z;
            return *this;
        }

        friend Vector3<T> operator+(Vector3<T> lhs, Vector3<T> const &rhs) {
            return lhs += rhs;
        }

        Vector3<T> &operator-=(Vector3<T> const &rhs) {
            x -= rhs.x;
            y -= rhs.y;
            z -= rhs.z;
            return *this;
        }

        friend Vector3<T> operator-(Vector3<T> lhs, Vector3<T> const &rhs) {
            return lhs -= rhs;
        }

        Vector3<T> &operator*=(T const &rhs) {
            x *= rhs;
            y *= rhs;
            z *= rhs;
            return *this;
        }

        friend Vector3<T> operator*(Vector3<T> lhs, T const &rhs) {
            return lhs *= rhs;
        }

        Vector3<T> &operator/=(T const &rhs) {
            x /= rhs;
            y /= rhs;
            z /= rhs;
            return *this;
        }

        friend Vector3<T> operator/(Vector3<T> lhs, T const &rhs) {
            return lhs /= rhs;
        }

        bool operator==(Vector3<T> const &v) const {
            return (x == v.x) && (y == v.y) && (z == v.z);
        }

        T dotProduct(Vector3<T> const &v) const {
            return x * v.x + y * v.y + z * v.z;
        }

        Vector3<T> operator-() const {
            return Vector3<T> {-x, -y, -z};
        }

        T getLengthSquared() const {
            return x * x + y * y + z * z;
        }

    };

}

#endif //DNA_REPLICATION_VECTOR3_H
