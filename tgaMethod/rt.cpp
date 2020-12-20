#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include "nicemath.h"

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
    ray(const nm::float3 &o,
        const nm::float3 &d):
        origin_ (o),
        direction_ (nm::normalize(d)) {}
    
    const nm::float3& origin() const { return origin_; }
    const nm::float3& direction() const { return direction_; }

    nm::float3 point_at(float t)
    {
        return origin_ + direction_ * t;
    }

private:
    nm::float3 origin_;
    nm::float3 direction_;
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

nm::float3 color(const ray &r)
{
    const float t = 0.5f * (r.direction().y() + 1.0f);
    return (1.0f - t) * nm::float3 {1.0f, 1.0f, 1.0f} +
                   t  * nm::float3 {0.2f, 0.1f, 0.8f};
}

int main(int argc, char const *argv[])
{
    frameBuffer fb {200u, 100u};
    camera cam {4.0f, 2.0f};
    for (size_t r = 0u; r < fb.height(); r++)
    {
        for (size_t c = 0u; c < fb.width(); c++)
        {
            const float u = (float) c / (float) fb.width();
            const float v = (float) r / (float) fb.height();
            nm::float3 col = color(cam.get_ray(u,v));
            fb.set_pixel(r, c, 255.99f * col.x(), 255.99f * col.y(), 255.99f * col.z());
        }
    }
    fb.save("image.tga");
    return 0;
}
