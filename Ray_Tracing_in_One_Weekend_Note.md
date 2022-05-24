# Ray Tracing in One Weekend

## 1. Overview

Website: [Ray Tracing in One Weekend](https://raytracing.github.io/books/RayTracingInOneWeekend.html#overview)

GitHub: [https://github.com/RayTracing/raytracing.github.io/](https://github.com/RayTracing/raytracing.github.io/)

PPM Viewer: [https://www.cs.rhodes.edu/welshc/COMP141_F16/ppmReader.html](https://www.cs.rhodes.edu/welshc/COMP141_F16/ppmReader.html)



## 2. Output an Image

### 2.1 The PPM Image Format

作为渲染器，我们的最终目的是生成一张图片

一张图片包含了很多像素点，我们给每一个像素点一个颜色值，这样整张图片就有了我们想要的效果

在这个渲染器中，我们使用 **PPM** 格式的图片，这个格式的文件能够记录图片的大小，每个像素点上的颜色值

```PPM
P3
# The P3 mean colors are in ASCII, then 3 colums and 2 rows
# then 255 for max color, then RGB triples
3 2
255
255		0		0		0		255		0		0		0		255
255		255		0		255		255		255		0		0		0
```



我们创建一个 `main.cpp` 来输出一张颜色图片

```C++
#include <iostream>

int main() {

    // Image

    const int image_width = 256;
    const int image_height = 256;

    // Render

    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    for (int j = image_height-1; j >= 0; --j) {
        for (int i = 0; i < image_width; ++i) {
            auto r = double(i) / (image_width-1);
            auto g = double(j) / (image_height-1);
            auto b = 0.25;

            int ir = static_cast<int>(255.999 * r);
            int ig = static_cast<int>(255.999 * g);
            int ib = static_cast<int>(255.999 * b);

            std::cout << ir << ' ' << ig << ' ' << ib << '\n';
        }
    }
}
// [main.cpp] Creating first image
```



### 2.2 Creating an Image File

我们在 VS 中生成项目解决方案

之后我们将生成的内容重定向到一个 PPM 格式的文件中

最后使用 PPM 读取器来读取并显示图片

```C++
build\Release\inOneWeekend.exe > image.ppm
```



![image_1](..\RayTracingInOneWeekend\.assets\image_1.png)



## 2.3 Adding a Progress Indicator

由于渲染是一个漫长的过程，我们可以在程序中添加一个输出，能够输出我们还有多少像素需要渲染

```C++
for (int j = image_height-1; j >= 0; --j) {
    std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
    for (int i = 0; i < image_width; ++i) {
        auto r = double(i) / (image_width-1);
        auto g = double(j) / (image_height-1);
        auto b = 0.25;

        int ir = static_cast<int>(255.999 * r);
        int ig = static_cast<int>(255.999 * g);
        int ib = static_cast<int>(255.999 * b);

        std::cout << ir << ' ' << ig << ' ' << ib << '\n';
    }
}

std::cerr << "\nDone.\n";
// [main.cpp] Main render loop with progress reporting
```



## 3. The Vec3 Class

对于任何渲染器来讲，我们都需要一种方法来表示空间中物体的位置，变化，颜色等等

这里我们使用一个 `vec3` 的类，通过这个类即可以表示向量，也能够表示颜色，以及方向，空间中的点

```C++
#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>

using std::sqrt;

class vec3 {
    public:
        vec3() : e{0,0,0} {}
        vec3(double e0, double e1, double e2) : e{e0, e1, e2} {}

        double x() const { return e[0]; }
        double y() const { return e[1]; }
        double z() const { return e[2]; }

        vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); }
        double operator[](int i) const { return e[i]; }
        double& operator[](int i) { return e[i]; }

        vec3& operator+=(const vec3 &v) {
            e[0] += v.e[0];
            e[1] += v.e[1];
            e[2] += v.e[2];
            return *this;
        }

        vec3& operator*=(const double t) {
            e[0] *= t;
            e[1] *= t;
            e[2] *= t;
            return *this;
        }

        vec3& operator/=(const double t) {
            return *this *= 1/t;
        }

        double length() const {
            return sqrt(length_squared());
        }

        double length_squared() const {
            return e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
        }

    public:
        double e[3];
};

// Type aliases for vec3
using point3 = vec3;   // 3D point
using color = vec3;    // RGB color

#endif
// [vec3.h] vec3 class
```



### 3.2 vec3 Utility Functions

之后我们在类的外部添加一些应用的函数，比如向量的点乘和叉乘

```C++
// vec3 Utility Functions

inline std::ostream& operator<<(std::ostream &out, const vec3 &v) {
    return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
}

inline vec3 operator+(const vec3 &u, const vec3 &v) {
    return vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

inline vec3 operator-(const vec3 &u, const vec3 &v) {
    return vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

inline vec3 operator*(const vec3 &u, const vec3 &v) {
    return vec3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

inline vec3 operator*(double t, const vec3 &v) {
    return vec3(t*v.e[0], t*v.e[1], t*v.e[2]);
}

inline vec3 operator*(const vec3 &v, double t) {
    return t * v;
}

inline vec3 operator/(vec3 v, double t) {
    return (1/t) * v;
}

inline double dot(const vec3 &u, const vec3 &v) {
    return u.e[0] * v.e[0]
         + u.e[1] * v.e[1]
         + u.e[2] * v.e[2];
}

inline vec3 cross(const vec3 &u, const vec3 &v) {
    return vec3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
                u.e[2] * v.e[0] - u.e[0] * v.e[2],
                u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

inline vec3 unit_vector(vec3 v) {
    return v / v.length();
}
// [vec3.h] utility functions
```



### 3.3 Color Utility Functions

我们现在有了 `vec3` 类了，可以通过这个类来控制我们颜色的输出

```C++
#ifndef COLOR_H
#define COLOR_H

#include "vec3.h"

#include <iostream>

void write_color(std::ostream &out, color pixel_color) {
    // Write the translated [0,255] value of each color component.
    out << static_cast<int>(255.999 * pixel_color.x()) << ' '
        << static_cast<int>(255.999 * pixel_color.y()) << ' '
        << static_cast<int>(255.999 * pixel_color.z()) << '\n';
}

#endif
// [color.h] color utility functions
```



之后我们可以修改我们的 `main.cpp` 来使用我们的 vec3 以及 color

```C++
#include "color.h"
#include "vec3.h"

#include <iostream>

int main() {

    // Image

    const int image_width = 256;
    const int image_height = 256;

    // Render

    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            color pixel_color(double(i)/(image_width-1), double(j)/(image_height-1), 0.25);
            write_color(std::cout, pixel_color);
        }
    }

    std::cerr << "\nDone.\n";
}
// [main.cc] Final code for the first PPM image
```



## 4. Rays, a Simple Camera, and Background

### 4.1 The ray class

光沿着直线传播，当光在空间中传播时会形成一条射线

而路径追踪就是模拟光在空间中的传递，反射，折射等

我们给定一个光线的起点 $O$，给定传播方向 $d$，以及一个时间 $t$，我们就能够计算在时间 $t$ 时光线到达的位置

$P = O + td$

这样我们可以实现一个 `ray.h` 来表示射线

```C++
#ifndef RAY_H
#define RAY_H

#include "vec3.h"

class ray {
    public:
        ray() {}
        ray(const point3& origin, const vec3& direction)
            : orig(origin), dir(direction)
        {}

        point3 origin() const  { return orig; }
        vec3 direction() const { return dir; }

        point3 at(double t) const {
            return orig + t*dir;
        }

    public:
        point3 orig;
        vec3 dir;
};

#endif
// [ray.h] The ray class
```



### 4.2 Sending Rays Into the Scene

光线追踪是将光线射向每一个像素点的方向，通过判断是否有物体与当前射线相交，来计算这个像素点上的颜色值

我们需要首先计算一条从摄像机到达像素点的射线

之后我们要判断这条射线是否与场景中的哪些物体相交

接着我们计算交点处的颜色值

最后将这个颜色值赋予到当前选定的像素点上

之后对其他每一个像素点做相同的计算来得到整个场景的渲染图

我们将摄像机放置在 $(0, 0, 0)$ 的位置，并将相机的朝向设置为 $-Z$ 轴

```C++
#include "color.h"
#include "ray.h"
#include "vec3.h"

#include <iostream>

color ray_color(const ray& r) {
    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*color(1.0, 1.0, 1.0) + t*color(0.5, 0.7, 1.0);
}

int main() {

    // Image
    const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 400;
    const int image_height = static_cast<int>(image_width / aspect_ratio);

    // Camera

    auto viewport_height = 2.0;
    auto viewport_width = aspect_ratio * viewport_height;
    auto focal_length = 1.0;

    auto origin = point3(0, 0, 0);
    auto horizontal = vec3(viewport_width, 0, 0);
    auto vertical = vec3(0, viewport_height, 0);
    auto lower_left_corner = origin - horizontal/2 - vertical/2 - vec3(0, 0, focal_length);

    // Render

    std::cout << "P3\n" << image_width << " " << image_height << "\n255\n";

    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            auto u = double(i) / (image_width-1);
            auto v = double(j) / (image_height-1);
            ray r(origin, lower_left_corner + u*horizontal + v*vertical - origin);
            color pixel_color = ray_color(r);
            write_color(std::cout, pixel_color);
        }
    }

    std::cerr << "\nDone.\n";
}
// [main.cc] Rendering a blue-to-white gradient
```

这里通过设置 `lower left corner` 我们将摄像机置于了中心点

我们将光线通过 `view port` 来实现透视投影的效果并判断射线是否与场景内的物体相交



生成图片的效果：

![image_2](..\RayTracingInOneWeekend\.assets\image_2.png)



## 5. Adding a Sphere

现在我们可以向场景内添加一些球体来进行渲染了



### 5.1 Ray-Sphere Intersection

对于一个球体来说，其表面上所有点 $p$，到圆心 $c$ 的距离是固定为半径 $R$ 的，所以我们可以得到隐式方程

Sphere： $ p: (p - c)^2-R^2 = 0$

那么如果一束光线射向球体，光线会在何时与球体相交呢？

我们可以将光线的等式带入球体隐式曲面的方程

得到：

$(o + td - c)^2 - R^2 = 0$

$at^2+bt+c=0$，where：

$a = d{\cdot}d$

$b=2(o-c){\cdot}d$

$c=(o-c){\cdot}{o-c}-R^2$

$t=\frac{-b\pm{\sqrt{b^2-4ac}}}{2a}$

根据表达式我们可以求得一个交点，两个交点或没有交点

对于两个交点，我们可以得到球体的射入点和射出点



### 5.2 Creating Our First Raytraced Image

那么，根据上述公式，我们可以开始渲染我们的第一个光线追踪的图像了

```C++
bool hit_sphere(const point3& center, double radius, const ray& r) {
    vec3 oc = r.origin() - center;
    auto a = dot(r.direction(), r.direction());
    auto b = 2.0 * dot(oc, r.direction());
    auto c = dot(oc, oc) - radius*radius;
    auto discriminant = b*b - 4*a*c;
    return (discriminant > 0);
}

color ray_color(const ray& r) {
    if (hit_sphere(point3(0,0,-1), 0.5, r))
        return color(1, 0, 0);
    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*color(1.0, 1.0, 1.0) + t*color(0.5, 0.7, 1.0);
}
// [main.cpp] Rendering a red sphere
```



![image_3](C:\Users\liang\Desktop\Note\Computer Graphics\RayTracingInOneWeekend\.assets\image_3.png)



## 6. Surface Normals and Multiple Objects

### 6.1 Shading with Surface Normals

我们现在可以渲染出一个圆的图像了，那么下一步我们可以去渲染以下表面的法线向量

法线向量是一个垂直于表面交点的单位向量

通过法线向量我们可以实现很多 shader 的效果，比如轮廓高光等

同时我们也可以通过法线向量来进行 debug

```C++
double hit_sphere(const point3& center, double radius, const ray& r) {
    vec3 oc = r.origin() - center;
    auto a = dot(r.direction(), r.direction());
    auto b = 2.0 * dot(oc, r.direction());
    auto c = dot(oc, oc) - radius*radius;
    auto discriminant = b*b - 4*a*c;
    if (discriminant < 0) {
        return -1.0;
    } else {
        return (-b - sqrt(discriminant) ) / (2.0*a);
    }
}

color ray_color(const ray& r) {
    auto t = hit_sphere(point3(0,0,-1), 0.5, r);
    if (t > 0.0) {
        vec3 N = unit_vector(r.at(t) - vec3(0,0,-1));
        return 0.5*color(N.x()+1, N.y()+1, N.z()+1);
    }
    vec3 unit_direction = unit_vector(r.direction());
    t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*color(1.0, 1.0, 1.0) + t*color(0.5, 0.7, 1.0);
}
// [main.cpp] Rendering surface normals on a sphere
```

![image_4](C:\Users\liang\Desktop\Note\Computer Graphics\RayTracingInOneWeekend\.assets\image_4.png)

生成的法线向量图就不是单一颜色的图了而是有了颜色的渐变，也体现了球体在空间中的感觉



### 6.2 Simplifying the Ray-Sphere Intersection Code

我们根据球体交点公式对代码进行一些修改和优化

```C++
double hit_sphere(const point3& center, double radius, const ray& r) {
    vec3 oc = r.origin() - center;
    auto a = r.direction().length_squared();
    auto half_b = dot(oc, r.direction());
    auto c = oc.length_squared() - radius*radius;
    auto discriminant = half_b*half_b - a*c;

    if (discriminant < 0) {
        return -1.0;
    } else {
        return (-half_b - sqrt(discriminant) ) / a;
    }
}
// [main.cpp] Ray-sphere intersection code (after)
```



### 6.3 An Abstraction for Hittable Objects

现在我们在场景内可以渲染出一个球体的图形了

那么如果我们想要渲染更多的球体呢？每个球体可以是不同的颜色，每个球体也可以有不同的材质

我们将球体这个概念抽象成为一个 `hittable` 的类

这个类是对于所有可以交互物体的抽象类

这个类中包含了一个 struct 用来记录交点的一些信息

```C++
#ifndef HITTABLE_H
#define HITTABLE_H

#include "ray.h"

struct hit_record {
    point3 p;
    vec3 normal;
    double t;
};

class hittable {
    public:
        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const = 0;
};

#endif
// [hitabble.h] The hittable class
```

我们在类中定义了一个纯虚函数 `hit`，这样继承自这个类的所有派生类都需要实现自身版本的 `hit` 函数



```C++
#ifndef SPHERE_H
#define SPHERE_H

#include "hittable.h"
#include "vec3.h"

class sphere : public hittable {
    public:
        sphere() {}
        sphere(point3 cen, double r) : center(cen), radius(r) {};

        virtual bool hit(
            const ray& r, double t_min, double t_max, hit_record& rec) const override;

    public:
        point3 center;
        double radius;
};

bool sphere::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    vec3 oc = r.origin() - center;
    auto a = r.direction().length_squared();
    auto half_b = dot(oc, r.direction());
    auto c = oc.length_squared() - radius*radius;

    auto discriminant = half_b*half_b - a*c;
    if (discriminant < 0) return false;
    auto sqrtd = sqrt(discriminant);

    // Find the nearest root that lies in the acceptable range.
    auto root = (-half_b - sqrtd) / a;
    if (root < t_min || t_max < root) {
        root = (-half_b + sqrtd) / a;
        if (root < t_min || t_max < root)
            return false;
    }

    rec.t = root;
    rec.p = r.at(rec.t);
    rec.normal = (rec.p - center) / radius;

    return true;
}

#endif
// [sphere.h] The sphere class
```



### 6.4 Front Faces Versus Back Faces

现在我们的法线向量永远是从球体中心点朝向交点的，也就是朝外的

那么如果当射线是从球体内部射出与球体表面相交时，法线和射线的方向是相同的

我们将其修改成永远和射线方向相反

也就是说，当射线从外部与球体相交时，法线是向外的

当射线从球体内部与球体表面相交时，法线是向内的

我们这么做是因为我们最终会需要判断哪一个交点的面是朝向摄像机的

对于朝向摄像机和背向摄像机的面，我们有不同的渲染效果

这里我们通过判断射线和法线方向是否同向来修改法线向量的朝向以及判断是否为面向摄像机的面

```C++
bool front_face;
if (dot(ray_direction, outward_normal) > 0.0) {
    // ray is inside the sphere
    normal = -outward_normal;
    front_face = false;
} else {
    // ray is outside the sphere
    normal = outward_normal;
    front_face = true;
}
```

将这个公式简化并添加到我们的 `hit record` 中

```C++
struct hit_record {
    point3 p;
    vec3 normal;
    double t;
    bool front_face;

    inline void set_face_normal(const ray& r, const vec3& outward_normal) {
        front_face = dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal :-outward_normal;
    }
};
// [hittable.h]
```

之后我们在 `sphere` 类中添加相应的代码

```C++
bool sphere::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    ...

    rec.t = root;
    rec.p = r.at(rec.t);
    vec3 outward_normal = (rec.p - center) / radius;
    rec.set_face_normal(r, outward_normal);

    return true;
}
// [sphere.h] The sphere class with normal determination
```



### 6.5 A List of Hittable Objects

我们现在有一个类叫做 `hittable`，这是一个可以进行交互的类

之后我们创建一个派生类，叫做 `hittable_list`，这个类中包含了许多 `hittable` 的对象

实现了对一个对象集合进行交点判断

```C++
#ifndef HITTABLE_LIST_H
#define HITTABLE_LIST_H

#include "hittable.h"

#include <memory>
#include <vector>

using std::shared_ptr;
using std::make_shared;

class hittable_list : public hittable {
    public:
        hittable_list() {}
        hittable_list(shared_ptr<hittable> object) { add(object); }

        void clear() { objects.clear(); }
        void add(shared_ptr<hittable> object) { objects.push_back(object); }

        virtual bool hit(
            const ray& r, double t_min, double t_max, hit_record& rec) const override;

    public:
        std::vector<shared_ptr<hittable>> objects;
};

bool hittable_list::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    hit_record temp_rec;
    bool hit_anything = false;
    auto closest_so_far = t_max;

    for (const auto& object : objects) {
        if (object->hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }

    return hit_anything;
}

#endif
// [hittable_list.h]
```



### 6.6 Some New C++ Features

这里我们使用了 `shared_ptr<hittable>` 

`shared_ptr` 是一个智能指针，它内部有一个指针会去记录所有指向相同内存的对象

这里使用 `shared_ptr` 的意义在于我们的可以让多个球体共享一些相同的内容，比如材质的类型



### 6.7 Common Constants and Utility Functions

我们创建一个新的文件来记录一些常用的常量和一些辅助函数

```C++
#ifndef RTWEEKEND_H
#define RTWEEKEND_H

#include <cmath>
#include <limits>
#include <memory>


// Usings

using std::shared_ptr;
using std::make_shared;
using std::sqrt;

// Constants

const double infinity = std::numeric_limits<double>::infinity();
const double pi = 3.1415926535897932385;

// Utility Functions

inline double degrees_to_radians(double degrees) {
    return degrees * pi / 180.0;
}

// Common Headers

#include "ray.h"
#include "vec3.h"

#endif
// [rtweekend.h] The rtweekend.h common header
```



之后更新以下我们的 `main.cpp`

```C++
#include "rtweekend.h"

#include "color.h"
#include "hittable_list.h"
#include "sphere.h"

#include <iostream>
color ray_color(const ray& r, const hittable& world) {
    hit_record rec;
    if (world.hit(r, 0, infinity, rec)) {
        return 0.5 * (rec.normal + color(1,1,1));
    }
    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*color(1.0, 1.0, 1.0) + t*color(0.5, 0.7, 1.0);
}

int main() {

    // Image

    const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 400;
    const int image_height = static_cast<int>(image_width / aspect_ratio);

    // World
    hittable_list world;
    world.add(make_shared<sphere>(point3(0,0,-1), 0.5));
    world.add(make_shared<sphere>(point3(0,-100.5,-1), 100));

    // Camera

    auto viewport_height = 2.0;
    auto viewport_width = aspect_ratio * viewport_height;
    auto focal_length = 1.0;

    auto origin = point3(0, 0, 0);
    auto horizontal = vec3(viewport_width, 0, 0);
    auto vertical = vec3(0, viewport_height, 0);
    auto lower_left_corner = origin - horizontal/2 - vertical/2 - vec3(0, 0, focal_length);

    // Render

    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            auto u = double(i) / (image_width-1);
            auto v = double(j) / (image_height-1);
            ray r(origin, lower_left_corner + u*horizontal + v*vertical);
            color pixel_color = ray_color(r, world);
            write_color(std::cout, pixel_color);
        }
    }

    std::cerr << "\nDone.\n";
}
// [main.cc] The new main with hittables
```



经过优化之后我们的程序可以输出以下的图像了

![image_5](..\RayTracingInOneWeekend\.assets\image_5.png)



## 7. Anti-Aliasing

此时我们的图像看上去会有一些锯齿

出现锯齿的原因是因为在我们采样的过程中，我们只会判断一个像素点是否会有物体相交

但是并不会去判断这个像素点有多少的部分是被物体所覆盖的

这就导致了像素点与像素点之间没有颜色的过渡

所以我们需要使用超采样的方式来做一些抗锯齿



### 7.1 Some Random Number Utilities

我们需要创建一些能够生成随机数的方法

```C++
#include <cstdlib>
...

inline double random_double() {
    // Returns a random real in [0,1).
    return rand() / (RAND_MAX + 1.0);
}

inline double random_double(double min, double max) {
    // Returns a random real in [min,max).
    return min + (max-min)*random_double();
}
// [rtweekend.h] random_double() functions
```



### 7.2 Generating Pixels with Multiple Samples

超采样的方式是在于，原来我们对一个像素点进行采样时只对该像素点的中点发出一条射线并判断是否有物体与其相交

而进行超采样我们则是在像素点内随机多个点并生成射线与场景内的物体判断是否相交

之后将所有的颜色值累加并求平均值

这样如果物体覆盖一个像素点 50% 的话，那么这个像素点的颜色也就是 50% 的颜色值

我们创建一个 `camera.h` 类来管理我们的摄像机并且用来执行超采样的工作

```C++
#ifndef CAMERA_H
#define CAMERA_H

#include "rtweekend.h"

class camera {
    public:
        camera() {
            auto aspect_ratio = 16.0 / 9.0;
            auto viewport_height = 2.0;
            auto viewport_width = aspect_ratio * viewport_height;
            auto focal_length = 1.0;

            origin = point3(0, 0, 0);
            horizontal = vec3(viewport_width, 0.0, 0.0);
            vertical = vec3(0.0, viewport_height, 0.0);
            lower_left_corner = origin - horizontal/2 - vertical/2 - vec3(0, 0, focal_length);
        }

        ray get_ray(double u, double v) const {
            return ray(origin, lower_left_corner + u*horizontal + v*vertical - origin);
        }

    private:
        point3 origin;
        point3 lower_left_corner;
        vec3 horizontal;
        vec3 vertical;
};
#endif
```



我们同时在 `rtweekend.h` 文件中添加一个新的辅助函数

```C++
inline double clamp(double x, double min, double max) {
    if (x < min) return min;
    if (x > max) return max;
    return x;
}
```

 

之后就是对 `color.h` 的修改，因为我们每个像素点会进行多次的相交判断，所以我们每次交点的颜色都会被累加起来

最后求平均值来得到这个像素点的颜色

```C++
void write_color(std::ostream &out, color pixel_color, int samples_per_pixel) {
    auto r = pixel_color.x();
    auto g = pixel_color.y();
    auto b = pixel_color.z();

    // Divide the color by the number of samples.
    auto scale = 1.0 / samples_per_pixel;
    r *= scale;
    g *= scale;
    b *= scale;

    // Write the translated [0,255] value of each color component.
    out << static_cast<int>(256 * clamp(r, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * clamp(g, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * clamp(b, 0.0, 0.999)) << '\n';
}
// [color.h] The multi-sample write_color() function
```



在 `main.cpp` 中我们就添加对于超采样的代码

在每个像素点内，我们都会设定一个 `samples_per_pixel` 的值，用来创建这个数量的随机点

对于每个随机点我们都判断是否射线会与其相交并将颜色值累加，最后求得平均值

```C++
#include "camera.h"

...

int main() {

    // Image

    const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 400;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 100;

    // World

    hittable_list world;
    world.add(make_shared<sphere>(point3(0,0,-1), 0.5));
    world.add(make_shared<sphere>(point3(0,-100.5,-1), 100));

    // Camera
    camera cam;

    // Render

    std::cout << "P3\n" << image_width << " " << image_height << "\n255\n";

    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            color pixel_color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_double()) / (image_width-1);
                auto v = (j + random_double()) / (image_height-1);
                ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, world);
            }
            write_color(std::cout, pixel_color, samples_per_pixel);
        }
    }

    std::cerr << "\nDone.\n";
}
// [main.cc] Rendering with multi-sampled pixels
```



之后我们就可以看到超采样过的图片的边缘有了明显的抗锯齿效果

但是随之带来的就是采样速度的变慢以及更多的消耗

![image_6](..\RayTracingInOneWeekend\.assets\image_6.png)



## 8. Diffuse Materials

### 8.1 A Simple Diffuse Material

我们可以通过对光线表现的计算来模拟物体上不同的材质

第一种材质类型就是 **Diffuse Materials** 漫反射

漫反射就是光线到达物体表面后会向**随机方向进行反射**

漫反射材质不仅仅会反射光线，同时也会对光线进行吸收，越黑的物体吸收的光线就会越多



当射线与场景内的球体相交时，我们可以想象一个单位球体与交点处相切

我们通过在这个单位求体内随机一个点来表示漫反射的方向

```C++
vec3 random_in_unit_sphere() {
    while (true) {
        auto p = vec3::random(-1,1);
        if (p.length_squared() >= 1) continue;
        return p;
    }
}
// [vec3.h] random_in_unit_sphere function
```

我们首先在一个单位立方体内随机一个点

之后判断这个点是否位于单位球体之内

根据球体公式我们可知一个点在球体表面需要满足： $(X-C)^2 + (Y-C)^2 + (Z - C)^2 = R^2$

所以，对于一个单位球体来说，我们只需要判断 $X ^2 + Y^2 + Z^2 < R^2$ 即可判断随机点是否在单位球体内



在得到了随机点后，我们就能够计算出漫反射方向的射线

```C++
color ray_color(const ray& r, const hittable& world) {
    hit_record rec;

    if (world.hit(r, 0, infinity, rec)) {
        point3 target = rec.p + rec.normal + random_in_unit_sphere();
        return 0.5 * ray_color(ray(rec.p, target - rec.p), world);
    }

    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*color(1.0, 1.0, 1.0) + t*color(0.5, 0.7, 1.0);
}
// [main.cpp] ray_color() using a random ray direction
```

当射线与球体相交时我们可以知道交点 `rec.p`

`point3 target = rec.p + rec.normal + random_in_unit_sphere();` 

我们根据交点加上交点法线的向量可以得到和交点相切的单位球体的中点

之后再加上随机点的偏移量，即可得到这个随机点在世界坐标系中的位置

由此，我们可以计算反射光线的起点 `rec.p`，以及反射光线的方向 `target - rec.p`

我们再将这个反射光线带入场景内去计算是否会有新的反射光线生成



由于我们的漫反射材质不仅会反射光线，也会吸收光线，所以我们每次递归时都会乘上 `0.5` 进行一个衰减

如果射线并不与场景内的球体相交，那我们返回一个环境光



### 8.2 Limiting the Number of Child Rays

一个潜在的问题在于，如果在一个有很多物体的场景内，可能会存在光线一直在反射而不会停止

所以我们需要为我们的递归设置一个最大深度，当达到最大深度时，我们返回黑色

```C++
color ray_color(const ray& r, const hittable& world, int depth) {
    hit_record rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return color(0,0,0);

    if (world.hit(r, 0, infinity, rec)) {
        point3 target = rec.p + rec.normal + random_in_unit_sphere();
        return 0.5 * ray_color(ray(rec.p, target - rec.p), world, depth-1);
    }

    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*color(1.0, 1.0, 1.0) + t*color(0.5, 0.7, 1.0);
}
```

![image_7](..\RayTracingInOneWeekend\.assets\image_7.png)



### 8.3 Using Gamma Correction for Accurate Color Intensity

我们通过漫反射生成的图像的颜色比较深

这是因为现在大部分的图像处理器会预设我们的图像经过了 Gamma Correction，伽马矫正

伽马矫正就是通过一个函数曲线来表示

会对函数值 0 到 1 之间的值进行一个非线性的变换由此来达到图像整体变亮或者变暗的效果

这里，我们取 Gamma Correction 的值为 $\gamma = 1/2$，并对每一个颜色值进行矫正，即开平方

这样整体图像的颜色就会变亮

```C++
void write_color(std::ostream &out, color pixel_color, int samples_per_pixel) {
    auto r = pixel_color.x();
    auto g = pixel_color.y();
    auto b = pixel_color.z();

    // Divide the color by the number of samples and gamma-correct for gamma=2.0.
    auto scale = 1.0 / samples_per_pixel;
    r = sqrt(scale * r);
    g = sqrt(scale * g);
    b = sqrt(scale * b);

    // Write the translated [0,255] value of each color component.
    out << static_cast<int>(256 * clamp(r, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * clamp(g, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * clamp(b, 0.0, 0.999)) << '\n';
}
```



### 8.4 Fixing Shadow Acne

还有一个潜在的问题是射线击中物体的时间并不一定是在 $t = 0$，有可能是在 $t = -0.0000001$ 或者 $t = 0.000001$

对于这种极小的浮点数值我们也应该去做对应的处理

```C++
if (world.hit(r, 0.001, infinity, rec)) {
```



### 8.5 True Lambertian Reflection

此时我们在单位球体内生成的随机点有一个小问题，就是这个随机点会趋向于法线向量

随机点随机到角度较小的位置的概率比较小

所以我们希望对这个随机点进行一个归一化，生成一个 unit vector

此时这个随机点会在单位球体的表面上

```C++
inline vec3 random_in_unit_sphere() {
    ...
}
vec3 random_unit_vector() {
    return unit_vector(random_in_unit_sphere());
}
```



![image_8](..\RayTracingInOneWeekend\.assets\image_8.png)



### 8.6 An Alternative Diffuse Formulation

以上提到的两种方法其实是不准确的

我们无法从数学的角度上来证明随机点的分布是完全随机的

同时我们在日常中也很难见到完全漫反射的材质



所以我们引入另一种方法，更加直观的方法

我们生成远离命中点所有角度均匀的散射方向，而不是依赖于法线的角度

```C++
vec3 random_in_hemisphere(const vec3& normal) {
    vec3 in_unit_sphere = random_in_unit_sphere();
    if (dot(in_unit_sphere, normal) > 0.0) // In the same hemisphere as the normal
        return in_unit_sphere;
    else
        return -in_unit_sphere;
}
// [vec3.h] The random_in_hemisphere(normal) function
```

这种方法是在单位半球内找随机点

即，如果随机点和法线不是相同方向，则取反

保证随机点生成在法线方向的半球内



所以由此我们提到了三种不同的生成随机散射光线的方法

+ 在与交点处相切的单位球体内随机一个点并将交点与该点的向量作为散射光线的方向
+ 在与交点处相切的单位球体内随机一个点并将该点变成单位向量使其在单位球的表面
+ 在与交点处相切的单位半球体内随机一个点，该点的方向与法线相同



## 9. Metal

### 9.1 An Abstract Class for Materials

场景中的物体可以是不同的材质，有不同的反射效果，比如漫反射，镜面反射以及折射

我们可以创建一个全局的 material 类，将散射方法抽象，使其派生类能够实现自己具体的散射效果

```C++
#ifndef MATERIAL_H
#define MATERIAL_H

#include "rtweekend.h"

struct hit_record;

class material {
    public:
        virtual bool scatter(
            const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered
        ) const = 0;
};

#endif
// [material.h] The material class
```



### 9.2 A Data Structure to Describe Ray-Object Intersections

我们使用 `hit_record` 来记录射线交点的信息，比如交点的位置，交点的法线，等等

同样，我们可以在 `hit_record` 中来记录交点的材质信息

通过记录材质信息，我们可以调用对应材质的散射方法来达到实现不同材质散射的效果

```C++
#include "rtweekend.h"

class material;

struct hit_record {
    point3 p;
    vec3 normal;
    shared_ptr<material> mat_ptr;
    double t;
    bool front_face;

    inline void set_face_normal(const ray& r, const vec3& outward_normal) {
        front_face = dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal :-outward_normal;
    }
};
// [hittable.h] Hit record with added material pointer
```



同时，我们需要在球体类中添加球体的材质信息，这样就可以记录每一个球体的材质

```C++
class sphere : public hittable {
    public:
        sphere() {}
        sphere(point3 cen, double r, shared_ptr<material> m)
            : center(cen), radius(r), mat_ptr(m) {};

        virtual bool hit(
            const ray& r, double t_min, double t_max, hit_record& rec) const override;

    public:
        point3 center;
        double radius;
        shared_ptr<material> mat_ptr;
};

bool sphere::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    ...

    rec.t = root;
    rec.p = r.at(rec.t);
    vec3 outward_normal = (rec.p - center) / radius;
    rec.set_face_normal(r, outward_normal);
    rec.mat_ptr = mat_ptr;

    return true;
}

// [sphere.h] Ray-sphere intersection with added material information
```



### 9.3 Modeling Light Scatter and Reflectance

我们首先可以对之前实现的漫反射材质创建一个属于 material 的派生类

```C++
class lambertian : public material {
    public:
        lambertian(const color& a) : albedo(a) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered
        ) const override {
            auto scatter_direction = rec.normal + random_unit_vector();
            scattered = ray(rec.p, scatter_direction);
            attenuation = albedo;
            return true;
        }

    public:
        color albedo;
};
```

这里的 `albedo` 是材质本身的颜色

`scatter` 方法会生成一个在单位球体表面上的随机点，并用这个随机点作为散射的方向

而这里我们的衰变 `attenuation` 使用的是 `albedo`，也可以做一些修改使其除以一个概率

同时，我们需要注意，如果随机的单位向量与法线向量方向正好相反的话，我们会得到方向为 0 的结果

所以我们要对这种特殊情况做一个处理：

```C++
class vec3 {
    ...
    bool near_zero() const {
        // Return true if the vector is close to zero in all dimensions.
        const auto s = 1e-8;
        return (fabs(e[0]) < s) && (fabs(e[1]) < s) && (fabs(e[2]) < s);
    }
    ...
};
// [vec3.h] The vec3::near_zero() method
```



```C++
class lambertian : public material {
    public:
        lambertian(const color& a) : albedo(a) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered
        ) const override {
            auto scatter_direction = rec.normal + random_unit_vector();

            // Catch degenerate scatter direction
            if (scatter_direction.near_zero())
                scatter_direction = rec.normal;

            scattered = ray(rec.p, scatter_direction);
            attenuation = albedo;
            return true;
        }

    public:
        color albedo;
};
```



### 9.4 Mirrored Light Reflection

对于金属表面材质的物体来说，它们并不会像漫反射一样随机散射



我们可以推导出镜面反射的公式，并通过传入的参数来得到反射的方向

```C++
vec3 reflect(const vec3& v, const vec3& n) {
    return v - 2*dot(v,n)*n;
}
// [vec3.h] vec3 reflection function
```



对于金属材质的物体，我们可以直接使用这个公式来达到镜面反射的效果

```C++
class metal : public material {
    public:
        metal(const color& a) : albedo(a) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered
        ) const override {
            vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
            scattered = ray(rec.p, reflected);
            attenuation = albedo;
            return (dot(scattered.direction(), rec.normal) > 0);
        }

    public:
        color albedo;
};
// [material.h] Metal material with reflectance function
```



同时，我们需要根据这些变化来修改我们的 `ray_color` 方法

```C++
color ray_color(const ray& r, const hittable& world, int depth) {
    hit_record rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return color(0,0,0);

    if (world.hit(r, 0.001, infinity, rec)) {
        ray scattered;
        color attenuation;
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            return attenuation * ray_color(scattered, world, depth-1);
        return color(0,0,0);
    }

    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*color(1.0, 1.0, 1.0) + t*color(0.5, 0.7, 1.0);
}
```

这里会根据 `rec.met_ptr` 的类型来调用自己派生类覆写的 `scatter` 函数

通过该函数能够返回是否有散射发生，同时还会记录散射射线方向的信息以及 `attenuation` 的信息

之后，我们做相同的事情，对这个散射射线进行递归判断



### 9.5 A Scene with Metal Spheres

此时我们可以定义一些材质为金属的球体在我们的场景中

并且对这个场景进行渲染

```C++
...

#include "material.h"

...

int main() {

    // Image

    const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 400;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 100;
    const int max_depth = 50;

    // World

    hittable_list world;

    auto material_ground = make_shared<lambertian>(color(0.8, 0.8, 0.0));
    auto material_center = make_shared<lambertian>(color(0.7, 0.3, 0.3));
    auto material_left   = make_shared<metal>(color(0.8, 0.8, 0.8));
    auto material_right  = make_shared<metal>(color(0.8, 0.6, 0.2));

    world.add(make_shared<sphere>(point3( 0.0, -100.5, -1.0), 100.0, material_ground));
    world.add(make_shared<sphere>(point3( 0.0,    0.0, -1.0),   0.5, material_center));
    world.add(make_shared<sphere>(point3(-1.0,    0.0, -1.0),   0.5, material_left));
    world.add(make_shared<sphere>(point3( 1.0,    0.0, -1.0),   0.5, material_right));

    // Camera

    camera cam;

    // Render

    std::cout << "P3\n" << image_width << " " << image_height << "\n255\n";

    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            color pixel_color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_double()) / (image_width-1);
                auto v = (j + random_double()) / (image_height-1);
                ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, world, max_depth);
            }
            write_color(std::cout, pixel_color, samples_per_pixel);
        }
    }

    std::cerr << "\nDone.\n";
}
// [main.cpp]
```



![image_9](..\RayTracingInOneWeekend\.assets\image_9.png)



### 9.6 Fuzzy Reflection

我们可以通过对镜面反射射线进行一个随机的偏离来达到一种雾面反射的效果

```C++
class metal : public material {
    public:
        metal(const color& a, double f) : albedo(a), fuzz(f < 1 ? f : 1) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered
        ) const override {
            vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
            scattered = ray(rec.p, reflected + fuzz*random_in_unit_sphere());
            attenuation = albedo;
            return (dot(scattered.direction(), rec.normal) > 0);
        }

    public:
        color albedo;
        double fuzz;
};
```

首先添加一个 `fuzz` 变量来表示雾面反射效果系数

之后对于我们的反射射线，我们给其一个单位球体内的随机值的偏移量

这样就达到了雾面反射的效果



![image_10](..\RayTracingInOneWeekend\.assets\image_10.png)



## 10. Dielectrics

### 10.1 Refraction

另一种材质是电介质，比如水，玻璃，钻石等

当光线照射到这类材质时，它会分裂成反射光线和折射光线

我们将通过在反射或者折射之间随机选择来处理这个问题，并且每次交互只会产生一条散射光线



### 10.2 Snell's Law

斯涅尔定律告诉我们，给定一种介质的折射率，我们可以知道入射角以及折射角之间的关系 $sin\theta^{\prime} = \frac{\eta}{\eta^{\prime}}\dotproduct{sin\theta}$

其中，$\theta$ 是入射角，$\eta$ 以及 $\eta^{\prime}$ 是在不同介质中的折射率，一般空气是 1.0，玻璃是 1.3 - 1.7，钻石是 2.4



经过公式推导，我们可以得到一个关于折射光线的关系式，即 $R^\prime = R^\prime_\perp + R^{\prime}_\parallel $

同时我们对这两个子向量进行求解

$R^{\prime}_{\perp} = \frac{\eta}{\eta^{prime}}(R+(R{\dotproduct}n)n)$

$R^{\prime}_{\parallel} = - \sqrt{1-|R^{\prime}_{\perp}|^2 n}$

那么我们可以将这个关系式带入到我们的代码中：

```C++
vec3 refract(const vec3& uv, const vec3& n, double etai_over_etat) {
    auto cos_theta = fmin(dot(-uv, n), 1.0);
    vec3 r_out_perp =  etai_over_etat * (uv + cos_theta*n);
    vec3 r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.length_squared())) * n;
    return r_out_perp + r_out_parallel;
}
// [vec3.h] Refraction function
```



对于整个产生折射的材质的类我们可以写成：

```C++
class dielectric : public material {
    public:
        dielectric(double index_of_refraction) : ir(index_of_refraction) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered
        ) const override {
            attenuation = color(1.0, 1.0, 1.0);
            double refraction_ratio = rec.front_face ? (1.0/ir) : ir;

            vec3 unit_direction = unit_vector(r_in.direction());
            vec3 refracted = refract(unit_direction, rec.normal, refraction_ratio);

            scattered = ray(rec.p, refracted);
            return true;
        }

    public:
        double ir; // Index of Refraction
};
// [material.h] Dielectric material class that always refracts
```



### 10.3 Total Internal Reflection

这样的效果肯定是有问题的

当我们从一个高折射系数的材质到一个低折射系数的材质时会出现斯涅尔定律失效的情况

即，如果 $\eta = 1.5$ 以及 $\eta^{\prime} = 1.0$

那么我们最后得到的 $sin\theta^{\prime}$ 就会大于 1

由此，当这个等式不成立时，我们不能进行折射，而是只能选择反射

```C++
if (refraction_ratio * sin_theta > 1.0) {
    // Must Reflect
    ...
} else {
    // Can Refract
    ...
}
```

 

```C++
class dielectric : public material {
    public:
        dielectric(double index_of_refraction) : ir(index_of_refraction) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered
        ) const override {
            attenuation = color(1.0, 1.0, 1.0);
            double refraction_ratio = rec.front_face ? (1.0/ir) : ir;

            vec3 unit_direction = unit_vector(r_in.direction());
            double cos_theta = fmin(dot(-unit_direction, rec.normal), 1.0);
            double sin_theta = sqrt(1.0 - cos_theta*cos_theta);

            bool cannot_refract = refraction_ratio * sin_theta > 1.0;
            vec3 direction;

            if (cannot_refract)
                direction = reflect(unit_direction, rec.normal);
            else
                direction = refract(unit_direction, rec.normal, refraction_ratio);

            scattered = ray(rec.p, direction);
            return true;
        }

    public:
        double ir; // Index of Refraction
};
// [material.h] Dielectric material class with reflection
```



![image_11](..\RayTracingInOneWeekend\.assets\image_11.png)



### 10.4 Schlick Approximation

实际上，对于真实的玻璃材料来说其折射率是根据观察角度而变化的

当我们以接近于平行的角度观察玻璃时会发现其折射率接近于玻璃

我们可以使用 Christophe Schlick 的方程来近似得模拟这种效果

```C++
class dielectric : public material {
    public:
        dielectric(double index_of_refraction) : ir(index_of_refraction) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered
        ) const override {
            attenuation = color(1.0, 1.0, 1.0);
            double refraction_ratio = rec.front_face ? (1.0/ir) : ir;

            vec3 unit_direction = unit_vector(r_in.direction());
            double cos_theta = fmin(dot(-unit_direction, rec.normal), 1.0);
            double sin_theta = sqrt(1.0 - cos_theta*cos_theta);

            bool cannot_refract = refraction_ratio * sin_theta > 1.0;
            vec3 direction;
            if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double())
                direction = reflect(unit_direction, rec.normal);
            else
                direction = refract(unit_direction, rec.normal, refraction_ratio);

            scattered = ray(rec.p, direction);
            return true;
        }

    public:
        double ir; // Index of Refraction

    private:
        static double reflectance(double cosine, double ref_idx) {
            // Use Schlick's approximation for reflectance.
            auto r0 = (1-ref_idx) / (1+ref_idx);
            r0 = r0*r0;
            return r0 + (1-r0)*pow((1 - cosine),5);
        }
};
// [material.h] Full glass material
```



### 10.5 Modeling a Hollow Glass Sphere

由此，我们可以尝试去渲染一个嵌套的双层玻璃球体

```C++
world.add(make_shared<sphere>(point3( 0.0, -100.5, -1.0), 100.0, material_ground));
world.add(make_shared<sphere>(point3( 0.0,    0.0, -1.0),   0.5, material_center));
world.add(make_shared<sphere>(point3(-1.0,    0.0, -1.0),   0.5, material_left));
world.add(make_shared<sphere>(point3(-1.0,    0.0, -1.0),  -0.4, material_left));
world.add(make_shared<sphere>(point3( 1.0,    0.0, -1.0),   0.5, material_right));
```



## 11. Positional Camera

### 11.1 Camera Viewing Geometry

Field of View 是相机成像的一个重要概念

当 FOV 越大时，相机能够看到的范围就越广

我们的相机是朝向 $-Z$ 轴的，我假定垂直方向的 FOV 角度为 $\theta$

那么假设 $Z = -1$，我们能够计算出单位视场的高度 $h = tan(\frac{\theta}{2})$

有了这个单位高度，我们就能够去改变视口的大小

```C++
class camera {
    public:
        camera(
            double vfov, // vertical field-of-view in degrees
            double aspect_ratio
        ) {
            auto theta = degrees_to_radians(vfov);
            auto h = tan(theta/2);
            auto viewport_height = 2.0 * h;
            auto viewport_width = aspect_ratio * viewport_height;

            auto focal_length = 1.0;

            origin = point3(0, 0, 0);
            horizontal = vec3(viewport_width, 0.0, 0.0);
            vertical = vec3(0.0, viewport_height, 0.0);
            lower_left_corner = origin - horizontal/2 - vertical/2 - vec3(0, 0, focal_length);
        }

        ray get_ray(double u, double v) const {
            return ray(origin, lower_left_corner + u*horizontal + v*vertical - origin);
        }

    private:
        point3 origin;
        point3 lower_left_corner;
        vec3 horizontal;
        vec3 vertical;
};
// [camera.h] Camera with adjustable field-of-view (fov)
```



### 11.2 Positioning and Orienting the Camera

我们首先将相机放置的位置叫做 Look from，将相机看到的点成为 Look at

想象以下即便我们固定住了相机的朝向位置以及相机看向的点

我们仍然能够通过侧面旋转来转动相机

由此，我们需要定义一个 `up direction`，用来表示相机侧面旋转时 ”朝上“ 的方向

```C++
class camera {
    public:
        camera(
            point3 lookfrom,
            point3 lookat,
            vec3   vup,
            double vfov, // vertical field-of-view in degrees
            double aspect_ratio
        ) {
            auto theta = degrees_to_radians(vfov);
            auto h = tan(theta/2);
            auto viewport_height = 2.0 * h;
            auto viewport_width = aspect_ratio * viewport_height;

            auto w = unit_vector(lookfrom - lookat);
            auto u = unit_vector(cross(vup, w));
            auto v = cross(w, u);

            origin = lookfrom;
            horizontal = viewport_width * u;
            vertical = viewport_height * v;
            lower_left_corner = origin - horizontal/2 - vertical/2 - w;
        }

        ray get_ray(double s, double t) const {
            return ray(origin, lower_left_corner + s*horizontal + t*vertical - origin);
        }

    private:
        point3 origin;
        point3 lower_left_corner;
        vec3 horizontal;
        vec3 vertical;
};
// [camera.h] Positionable and orientable camera
```



## 12. Defocus Blur

Defocus Blue，散焦模糊，或者叫做景深

景深越深，离焦点远的物体也能够清晰

景深越浅，离焦点远的物体就会很模糊



### 12.1 A Thin Lens Approximation

我们不需要去模拟真实的相机成像

通常下，我们从镜头朝向 focus plane 射出射线，在 focus plane 上的所有物体都是完美对焦的

我们使用 `focus_dist` 来表示从镜头到 focus plane 上的距离



### 12.2 Generating Sample Rays

为了达到散焦模糊的效果，我们可以将光线从以 `lookfrom` 为中心点的圆内随机得选取射线的起点

对于越大的半径，我们得到的模糊效果就越明显

```C++
vec3 random_in_unit_disk() {
    while (true) {
        auto p = vec3(random_double(-1,1), random_double(-1,1), 0);
        if (p.length_squared() >= 1) continue;
        return p;
    }
}
// [vec3.h] Generate random point inside unit disk
```



```C++
class camera {
    public:
        camera(
            point3 lookfrom,
            point3 lookat,
            vec3   vup,
            double vfov, // vertical field-of-view in degrees
            double aspect_ratio,
            double aperture,
            double focus_dist
        ) {
            auto theta = degrees_to_radians(vfov);
            auto h = tan(theta/2);
            auto viewport_height = 2.0 * h;
            auto viewport_width = aspect_ratio * viewport_height;

            w = unit_vector(lookfrom - lookat);
            u = unit_vector(cross(vup, w));
            v = cross(w, u);

            origin = lookfrom;
            horizontal = focus_dist * viewport_width * u;
            vertical = focus_dist * viewport_height * v;
            lower_left_corner = origin - horizontal/2 - vertical/2 - focus_dist*w;

            lens_radius = aperture / 2;
        }


        ray get_ray(double s, double t) const {
            vec3 rd = lens_radius * random_in_unit_disk();
            vec3 offset = u * rd.x() + v * rd.y();

            return ray(
                origin + offset,
                lower_left_corner + s*horizontal + t*vertical - origin - offset
            );
        }

    private:
        point3 origin;
        point3 lower_left_corner;
        vec3 horizontal;
        vec3 vertical;
        vec3 u, v, w;
        double lens_radius;
};
// [camera.h] Camera with adjustable depth-of-field (dof)
```



![image_12](..\RayTracingInOneWeekend\.assets\image_12.png)



## 13. Next

下一步可以考虑将渲染器做成能够转动相机的

同时能够对三角形进行采样

以及实现 BVH 等加速算法



## Reference

[8. Diffuse Materials](https://zhuanlan.zhihu.com/p/442001301)

[图像处理算法之Gamma矫正](https://blog.csdn.net/kkae8643150/article/details/114078037?utm_medium=distribute.pc_relevant.none-task-blog-2~default~baidujs_baidulandingword~default-0-114078037-blog-78617615.pc_relevant_default&spm=1001.2101.3001.4242.1&utm_relevant_index=3)

