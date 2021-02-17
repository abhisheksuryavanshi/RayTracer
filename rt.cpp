#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <random>
#include <iostream>
#include "nicemath.h"
#include <random>

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

nm::float2 random_int_unit_disk() {
    nm::float2 result;
    do {
        result = 2.0f * nm::float2(randf(), randf()) - nm::float2( 1.0f, 1.0f);
    } while(nm::lengthsq(result) > 1.0f);
    return result;
}

nm::float3 random_int_unit_sphere() {
    nm::float3 result;
    do {
        result = 2.0f * nm::float3(randf(), randf(), randf()) - nm::float3( 1.0f, 1.0f, 1.0f );
    } while(nm::lengthsq(result) > 1.0f);
    return result;
}


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

nm::float3 reflect(const nm::float3 &i, const nm::float3 &n)
{
    return i - 2.0f * nm::dot(i,n) * n;
}

class dielectric : public material
{
public:
    explicit dielectric(float ri) : refraction_idx_(ri) {}
    
    bool scatter(   const ray           &in,
                            const hit_record    &rec,
                            nm::float3          &attn,
                            ray                 &scattered  )  const override
    {
        attn = nm::float3 { 1.0f, 1.0f, 1.0f };
        nm::float3 refracted_dir;
        bool is_refracted = false;
        const nm::float3 &n = rec.normal;
        const nm::float3 &i = in.direction();
        float refl_prob = 1.0f;
        float cosine = 0.0f;
        if(nm::dot(n,i) > 0.0f)
        {
            is_refracted = refract(i, -n, refraction_idx_, 1.0f, refracted_dir);
            cosine = refraction_idx_ * dot(i,n) / nm::length(i);
        }
        else
        {
            is_refracted = refract(i, n, 1.0f, refraction_idx_, refracted_dir);
            cosine = -nm::dot(i,n) / nm::length(i);
        }

        if(is_refracted)
        {
            refl_prob = schlick(cosine, refraction_idx_);
        }

        scattered = ray {
            rec.p,
            randf() < refl_prob ? reflect(i,n) : refracted_dir
        };

        return is_refracted;
    }
private:
    static float schlick(float cosine, float ri)
    {
        const float r0 = (1.0f - ri) / (1.0f + ri);
        const float r0sq = r0 * r0;
        return r0sq + (1.0f - r0sq) * pow((1.0f - cosine), 5);
    }
    static bool refract (const nm::float3 &i,
                    const nm::float3 &n,
                    float            ni,
                    float            nt,
                    nm::float3       &refracted)
    {
        const nm::float3 ui = nm::normalize(i);
        const float dt = nm::dot(i,n);
        const float ni_o_nt = ni / nt;
        float d = 1.0f - ni_o_nt * ni_o_nt * (1.0f - dt * dt);
        const bool is_refracted = d > 0.0f;
        if(is_refracted)
            refracted = ni_o_nt * (ui - n * dt) - n * sqrtf(d);   
        return is_refracted;
    }
    float refraction_idx_;
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
    explicit metal(const nm::float3 attn, float fuzz) :
        attn_(attn),
        fuzz_(fuzz) {}

    virtual bool scatter(   const ray           &in,
                            const hit_record    &rec,
                            nm::float3          &attn,
                            ray                 &scattered  ) const override
    {
        const nm::float3 refl_dir = 
        reflect(in.direction(), rec.normal) + 
        random_int_unit_sphere() * fuzz_;

        attn = attn_;
        scattered = ray {rec.p, refl_dir - rec.p};
        return true;
    }
                            
private:
    nm::float3 attn_;
    float fuzz_;
};

class camera
{
    static constexpr float Pi = 3.14159f;
public:
    camera( float fov_v,
            float aspect,
            const nm::float3 &look_from,
            const nm::float3 &look_at,
            const nm::float3 &upvector,
            float apreture):
        origin_ (look_from),
        lens_radius_ (apreture / 2.0f) 
        {
            const float fov_v_rad = fov_v * Pi / 180.0f;
            const float half_height = tan(fov_v_rad / 2.0f); // Half height
            const float half_width = half_height * aspect; // Half width
            const float focus_distance = nm::length(look_from - look_at);
            const nm::float3 camz = nm::normalize(look_from - look_at);
            const nm::float3 camx = nm::normalize(nm::cross(upvector, camz));
            const nm::float3 camy = nm::normalize(nm::cross(camz, camx));
            hvector_ = camx * 2.0f * focus_distance * half_width;
            vvector_ = camy * 2.0f * focus_distance * half_height;
            lower_left_ = origin_ - 
                          focus_distance * camz - half_width * 
                          focus_distance * (camx) - 
                          focus_distance * half_height * (camy);
        }
    
    // u = 0 => left edge ; u = 1 => right edge
    // v = 0 => bottom edge ; v = 1 => top edge
    //
    ray get_ray(float u, float v)
    {
        const nm::float2 rd = random_int_unit_disk() * lens_radius_;
        const nm::float3 ray_o = origin_ + nm::normalize(hvector_) * rd.x() + nm::normalize(vvector_) * rd.y();
        return ray {ray_o, 
                    lower_left_ + u * hvector_ + v * vvector_ - ray_o};
    }

private:
    nm::float3 hvector_;
    nm::float3 vvector_;
    nm::float3 lower_left_;
    nm::float3 origin_;
    float lens_radius_;
};

class hitable
{
public:
    virtual bool hit_test (const ray &r, 
                           float     tmin,
                           float     tmax,
                           hit_record &hit) const = 0;
};

static lambertian purple {nm::float3 {0.7f, 0.3f, 0.8f}};

class sphere : public hitable 
{
public:
    sphere(const nm::float3 &center = {0.0f, 0.0f, 0.0f}, 
           float radius = 0.01f,
           const material *mat = &purple) : 
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
            hit.normal = (hit.p - center_) / radius_;
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
    sphere_list(std::vector<sphere> list) : 
           list_(list) {}

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

nm::float3 color(const ray &r, int bounce, std::vector<sphere> new_scene)
{
    const float t = 0.5f * (r.direction().y() + 1.0f);
    static constexpr int max_bounce = 50;
    static lambertian diffuse_gray {nm::float3 {0.5f, 0.5f, 0.5f}};
    static lambertian diffuse_pink {nm::float3 {0.8f, 0.3f, 0.3f}};
    static lambertian diffuse_yellow {nm::float3 {0.8f, 0.8f, 0.0f}};
    static lambertian diffuse_blue {nm::float3 {0.1f, 0.2f, 0.5f}};
    static dielectric dielectric_1 {1.5f};
    static metal reddish_metal { nm::float3{0.8f, 0.6f, 0.2f }, 0.0f};
    static metal gray_metal { nm::float3{0.8f, 0.8f, 0.8f}, 0.3f};

    float R = cos(3.1015926f) / 2.0f;
    // static sphere_list scene {
    //     // sphere {nm::float3 { -R, 0.0f, -1.0f}, R, &purlple},
    //     // sphere {nm::float3 { R, 0.0f, -1.0f}, R, &diffuse_gray},
    //     sphere { nm::float3{0.0f, 0.0f, -1.0f}, 0.5f, &diffuse_blue},
    //     sphere { nm::float3{0.0f, -100.5f, -1.0f}, 100.0f, &diffuse_yellow},
    //     sphere { nm::float3{1.0f, 0.0f, -1.0f}, 0.5f, &reddish_metal},
    //     sphere { nm::float3{-1.0f, 0.0f, -1.0f}, 0.5f, &dielectric_1},
    //     sphere { nm::float3{-1.0f, 0.0f, -1.0f}, -0.48f, &dielectric_1}
    // };

    
    static sphere_list scene2 = new_scene;

    hit_record hit;
    if(bounce < max_bounce && scene2.hit_test(r, 0.001f, 100.0f, hit)){
        const nm::float3 target = hit.p + hit.normal + random_int_unit_sphere();
        nm::float3 attn;
        ray scattered;
        if(hit.mat->scatter(r, hit, attn, scattered))
        {
            return attn * color(scattered, bounce+1, new_scene);
        }
        return 0.5f * color(ray {hit.p, target - hit.p}, bounce + 1, new_scene);
    }
    return (1.0f - t) * nm::float3 {1.0f, 1.0f, 1.0f} +
                   t  * nm::float3 {0.5f, 0.7f, 1.0f};
}

int main(int argc, char const *argv[])
{
    std::cout<<"main called\n";
    frameBuffer fb {800u, 400u};
    camera cam {
                60.0f, 
                2.0f, 
                nm::float3{3.5f, 1.0f, 3.0f}, 
                nm::float3{0.0f, 1.0f, 0.0f},
                nm::float3{0.0f, 1.0f,  0.0f},
                0.15f
                };

    std::cout<<"scene defination start\n";
    int n = 150;
    std::vector<sphere> new_scene(n+1);
    new_scene[0] = sphere { nm::float3{0.0f, -1000.0f, -1.0f}, 1000.0f, new lambertian(nm::float3{0.5,0.5,0.5})};
    // new_scene.push_back(sphere { nm::float3{1.0f, 2.0f, 0.0f}, 0.5f, &diffuse_blue});
    int i = 1;
    for (int a = -5; a < 5; a++)
    {
        // std::cout<<"a = "<<a<<"\n";
        for (int b = -5; b < 5; b++)
        {
            // std::cout<<"b = "<<b<<"\n";
            float choose_mat = randf();
            nm::float3 center = nm::float3{a + randf(), 0.2, b + randf()};

            if(nm::length(center - nm::float3{0.0f, 1.0f, 0.0f}) > 0.9)
            {
                if(choose_mat <0.6)
                {
                    new_scene[i++] = sphere(center, 0.2f, new lambertian(nm::float3{randf(), randf(), randf()}));
                }
                else if(choose_mat < 0.8)
                {
                    new_scene[i++] = sphere(center, 0.2f, new metal(nm::float3{0.5f*(1.0f + randf()), 0.5f*(1.0f + randf()), 0.5f*(1.0f + randf())} , 0.5f * randf()));
                }
                else
                {
                    new_scene[i++] = sphere(center, 0.2f, new dielectric(1.5));
                }
            }
        }   
    }
    // std::cout<<"Done\n";
    new_scene[i++] = sphere(nm::float3{0.0f, 1.0f, 0.0f}, 1.0f, new dielectric(1.5));
    new_scene[i++] = sphere(nm::float3{-2.0f, 1.0f, 0.0f}, 1.0f, new lambertian(nm::float3{0.4, 0.2, 0.1}));
    new_scene[i++] = sphere(nm::float3{2.0f, 1.0f, 0.0f}, 1.0f, new metal(nm::float3{0.7, 0.6, 0.5},0.0));
    std::cout<<"scene defined\n";

    #pragma omp parallel for
    for (size_t r = 0u; r < fb.height(); r++)
    {
        std::cout<<"row = "<<r<<"\n";
        #pragma omp parallel for
        for (size_t c = 0u; c < fb.width(); c++)
        {
            // std::cout<<"col = "<<c<<" fb.width = "<<fb.width()<<"\n";
            nm::float3 col = { 0.0f, 0.0f, 0.0f };
            const size_t ns = 500u;
            #pragma omp parallel for
            for (size_t s = 0u; s < ns; s++)
            {
                // std::cout<<"antiC = "<<s<<"\n";
                const float u = ((float) c) / (float) fb.width();
                const float v = ((float) r) / (float) fb.height();
                col = col + color(cam.get_ray(u,v), 0, new_scene);
            }
            // std::cout<<"came out of color\n";
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
