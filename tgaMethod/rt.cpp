#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <random>
#include <iostream>
#include "nicemath.h"

float frand(int *seed)
{
    union
    {
        float fres;
        unsigned int ires;
    };

    seed[0] *= 16807;
    ires = ((((unsigned int)seed[0]) >> 9) | 0x3f800000);
    return fres - 1.0f;    
}

float randf()
{
    // static std::random_device rd;
    // static std::mt19937 gen(rd());
    // std::uniform_real_distribution<float> d { 0.0f, 1.0f };
    // return d(gen);

    static int seed = 15677;
    return  frand(&seed);
}

nm::float3 random_int_unit_sphere() {
    nm::float3 result;
    do {
        result = 2.0f * nm::float3(randf(), randf(), randf()) - nm::float3( 1.0f, 1.0f, 1.0f );
    } while(nm::lengthsq(result) > 1.0f);
    return result;
}

class ray
{
public:
    ray() = default;
    ray(const nm::float3 &o,
        const nm::float3 &d):
        origin_ (o),
        direction_ (nm::normalize(d)) {}
    
    const nm::float3& origin() const { return origin_; }
    const nm::float3& direction() const { return direction_; }

    nm::float3 point_at(float t) const
    {
        return origin_ + direction_ * t;
    }

private:
    nm::float3 origin_;
    nm::float3 direction_;
};

class material;
// Contains all the information of a hit point
//
struct hit_record
{
    nm::float3 normal;
    float t;
    nm::float3 p;
    const material *mat;
};

class material
{
public:
    virtual bool scatter(   const ray           &in,
                            const hit_record    &rec,
                            nm::float3          &attn,
                            ray                 &scattered  )  const = 0;
private:
};


class lambertian : public material
{
public:
    explicit lambertian(const nm::float3 albedo) :
        albedo_(albedo) {}

    virtual bool scatter(   const ray           &in,
                            const hit_record    &rec,
                            nm::float3          &attn,
                            ray                 &scattered  ) const override
    {
        const nm::float3 target = rec.p + rec.normal + random_int_unit_sphere();
        attn = albedo_;
        scattered = ray {rec.p, target - rec.p};
        return true;
    }
                            
private:
    nm::float3 albedo_;
};

class metal : public material
{
public:
    explicit metal(const nm::float3 attn) :
        attn_(attn) {}

    virtual bool scatter(   const ray           &in,
                            const hit_record    &rec,
                            nm::float3          &attn,
                            ray                 &scattered  ) const override
    {
        const nm::float3 refl_dir = 
        in.direction() - 2.0f * nm::dot(in.direction(), rec.normal) * rec.normal;
        attn = attn_;
        scattered = ray {rec.p, refl_dir - rec.p};
        return true;
    }
                            
private:
    nm::float3 attn_;
};

// Allows us to specify individual pixel colours and
// dump it to a file.
//
class frameBuffer
{
public:
    static const size_t kBytesPerPixel = 3;

    frameBuffer(size_t width, size_t height) :
        width_ (width),
        height_ (height),
        data_ ((uint8_t*)malloc(width * height * kBytesPerPixel)){}

    void set_pixel( size_t row,
                    size_t col,
                    uint8_t r,
                    uint8_t g,
                    uint8_t b)
    {
        const size_t idx = kBytesPerPixel * (row * width_ + col);
        data_[idx + 0] = b;
        data_[idx + 1] = g;
        data_[idx + 2] = r;
    }
    
    void save(const char* file_path) const
    {
        std::cout<<"In File buffer Save\n";
        FILE* fptr = fopen(file_path, "wb");
        assert(fptr);
        putc(0,fptr);
        putc(0,fptr);
        putc(2,fptr);                         /* uncompressed RGB */
        putc(0,fptr); putc(0,fptr);
        putc(0,fptr); putc(0,fptr);
        putc(0,fptr);
        putc(0,fptr); putc(0,fptr);           /* X origin */
        putc(0,fptr); putc(0,fptr);           /* y origin */
        putc((width_ & 0x00FF),fptr);
        putc((width_ & 0xFF00) / 256,fptr);
        putc((height_ & 0x00FF),fptr);
        putc((height_ & 0xFF00) / 256,fptr);
        putc(24,fptr);                        /* 24 bit bitmap */
        putc(0,fptr);
        fwrite(data_, kBytesPerPixel, width_ * height_, fptr);
        fclose(fptr);
        std::cout<<"Done\n";
    }

    size_t width() const { return width_; }
    size_t height() const { return height_; }

    ~frameBuffer(){ free(data_); }
private:
    uint8_t* data_;
    size_t width_;
    size_t height_;
};

class camera
{
public:
    camera(float aspect_h, float aspect_v):
        aspect_h_ (aspect_h),
        aspect_v_ (aspect_v),
        origin_ (0u, 0u, 0u) {}
    
    // u = 0 => left edge ; u = 1 => right edge
    // v = 0 => bottom edge ; v = 1 => top edge
    //
    ray get_ray(float u, float v)
    {
        const nm::float3 lower_left {-aspect_h_ / 2.0f, -aspect_v_ / 2.0f, -1.0f};
        return ray {origin_, lower_left + nm::float3{ u*aspect_h_, v*aspect_v_, 0.0f}};
    }

private:
    float aspect_h_;
    float aspect_v_;
    nm::float3 origin_;
};

class hitable
{
public:
    virtual bool hit_test (const ray &r, 
                           float     tmin,
                           float     tmax,
                           hit_record &hit) const = 0;
};

class sphere : public hitable 
{
public:
    sphere(const nm::float3 &center, 
           float radius,
           const material *mat) : 
           center_(center),  
           radius_(radius),
           mat_(mat) {}
    
    bool hit_test (const ray &r, 
                   float     tmin,
                   float     tmax,
                   hit_record &hit) const override 
    {
        const nm::float3 oc = r.origin() - center_;
        const float a = nm::dot(r.direction(), r.direction());
        const float b = nm::dot(oc, r.direction())*2.0f;
        const float c = nm::dot(oc, oc) - radius_ * radius_;
        const float d = b * b - 4 * a * c;
        if (d > 0.0f)
        {
            const float t1 = (-b - sqrt(d)) / (2.0f * a);
            const float t2 = (-b + sqrt(d)) / (2.0f * a);
            if(t1 > tmin && t1 < tmax) hit.t = t1;
            else if(t2 > tmin && t2 < tmax) hit.t = t2;
            else return false;
            hit.p = r.point_at(hit.t);
            hit.normal = nm::normalize(hit.p - center_);
            hit.mat = mat_;
        }
        return d > 0.0f;
    }
private:
    nm::float3 center_;
    float radius_;
    const material *mat_;
};

class sphere_list : hitable 
{
public:
    template <class ...Args>
    sphere_list(Args... args) : list_ { std::forward<Args>(args)... } {}

    bool hit_test (const ray &r, 
                   float     tmin,
                   float     tmax,
                   hit_record &hit) const override
    {
        bool anyhit = false;
        float closest_hit = tmax;
        for(const sphere &s : list_)
        {
            bool has_hit = s.hit_test(r, tmin, closest_hit, hit);
            if(has_hit) closest_hit = hit.t;
            anyhit |= has_hit;
        }
        return anyhit;
    }
private:
    std::vector<sphere> list_;
};

nm::float3 color(const ray &r, int bounce)
{
    const float t = 0.5f * (r.direction().y() + 1.0f);
    static constexpr int max_bounce = 50;
    static lambertian purlple {nm::float3 {0.7f, 0.3f, 0.8f}};
    static lambertian diffuse_gray {nm::float3 {0.5f, 0.8f, 0.5f}};
    static metal reddish_metal { nm::float3{0.8f, 0.6f, 0.2f }};
    static metal gray_metal { nm::float3{0.8f, 0.8f, 0.8f}};

    static sphere_list scene {
        sphere { nm::float3{0.0f, 0.0f, -1.0f}, 0.5f, &purlple},
        sphere { nm::float3{0.0f, -100.5f, -1.0f}, 100.0f, &diffuse_gray},
        sphere { nm::float3{1.5f, 0.5f, -2.0f}, 0.5f, &reddish_metal},
        sphere { nm::float3{-2.5f, 0.0f, -3.0f}, 0.5f, &gray_metal}
    };
    hit_record hit;
    if(bounce < max_bounce && scene.hit_test(r, 0.001f, 100.0f, hit)){
        const nm::float3 target = hit.p + hit.normal + random_int_unit_sphere();
        nm::float3 attn;
        ray scattered;
        if(hit.mat->scatter(r, hit, attn, scattered))
        {
            return attn * color(scattered, bounce+1);
        }
        return 0.5f * color(ray {hit.p, target - hit.p}, bounce + 1);
    }
    return (1.0f - t) * nm::float3 {1.0f, 1.0f, 1.0f} +
                   t  * nm::float3 {0.2f, 0.1f, 0.8f};
}

int main(int argc, char const *argv[])
{
    std::cout<<"main called\n";
    frameBuffer fb {800u, 400u};
    camera cam {4.0f, 2.0f};
    for (size_t r = 0u; r < fb.height(); r++)
    {
        for (size_t c = 0u; c < fb.width(); c++)
        {
            nm::float3 col = { 0.0f, 0.0f, 0.0f };
            const size_t ns = 50u;
            for (size_t s = 0u; s < ns; s++)
            {
                const float u = ((float) c) / (float) fb.width();
                const float v = ((float) r) / (float) fb.height();
                col = col + color(cam.get_ray(u,v), 0);
            }
            col /= (float)ns;
            fb.set_pixel(r, c, 
                         255.99f * sqrt(col.x()), 
                         255.99f * sqrt(col.y()), 
                         255.99f * sqrt(col.z()));
        }
    }
    std::cout<<"calling save now\n";
    fb.save("image.tga");
    return 0;
}
