#include "texture.h"
#include "color.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CMU462 {

inline void uint8_to_float( float dst[4], unsigned char* src ) {
  uint8_t* src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8( unsigned char* dst, float src[4] ) {
  uint8_t* dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[3])));
}

Color retrieveColor(std::vector<unsigned char> &arr, int index) {
    Color c{};
    c.r = (float)arr[4 * index] / 255;
    c.g = (float)arr[4 * index + 1] / 255;
    c.b = (float)arr[4 * index + 2] / 255;
    c.a = (float)arr[4 * index + 3] / 255;
    return c;
}

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

  // NOTE: 
  // This starter code allocates the mip levels and generates a level 
  // map by filling each level with placeholder data in the form of a 
  // color that differs from its neighbours'. You should instead fill
  // with the correct data!

  // Task 7: Implement this

  // check start level
  if ( startLevel >= tex.mipmap.size() ) {
    std::cerr << "Invalid start level"; 
  }

  // allocate sublevels
  int baseWidth  = tex.mipmap[startLevel].width;
  int baseHeight = tex.mipmap[startLevel].height;
  int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  tex.mipmap.resize(startLevel + numSubLevels + 1);

  int width  = baseWidth;
  int height = baseHeight;
  int prevWidth, prevHeight;
  Color c1, c2, c3, c4, c0;

  for (int i = 1; i <= numSubLevels; i++) {
    prevWidth = width;
    prevHeight = height;
    MipLevel& prevLevel = tex.mipmap[startLevel + i - 1];
    MipLevel& level = tex.mipmap[startLevel + i];
    // handle odd size texture by rounding down
    width  = max( 1, width  / 2); assert(width  > 0);
    height = max( 1, height / 2); assert(height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);

    for(size_t n = 0; n < width; n += 1) {
        for(size_t m = 0; m < height; m += 1) {
            c1 = retrieveColor(prevLevel.texels, 2 * n + 2 * m * prevWidth);
            c2 = retrieveColor(prevLevel.texels, 2 * n + (2 * m + 1) * prevWidth);
            c3 = retrieveColor(prevLevel.texels, 2 * n + 1 + (2 * m + 1) * prevWidth);
            c4 = retrieveColor(prevLevel.texels, 2 * n + 1 + 2 * m * prevWidth);
            c0 = (c1 + c2 + c3 + c4) * 0.25;
            float_to_uint8(&level.texels[4 * (n + m * width)], &c0.r);
        }
    }
  }



//  // fill all 0 sub levels with interchanging colors (JUST AS A PLACEHOLDER)
//  Color colors[3] = { Color(1,0,0,1), Color(0,1,0,1), Color(0,0,1,1) };
//  for(size_t i = 1; i < numSubLevels; ++i) {
//
//    Color c = colors[i % 3];
//    MipLevel& mip = tex.mipmap[i];
//
//    for(size_t i = 0; i < 4 * mip.width * mip.height; i += 4) {
//      float_to_uint8( &mip.texels[i], &c.r );
//    }
//  }



}



Color Sampler2DImp::sample_nearest(Texture& tex, 
                                   float u, float v, 
                                   int level) {


  // Task 6: Implement nearest neighbour interpolation
    float w, h;
    w = tex.mipmap[level].width;
    h = tex.mipmap[level].width;
    size_t index = round(u * w) + w * round(v * h);
    Color c= retrieveColor(tex.mipmap[level].texels, index);

  
  // return magenta for invalid level
  return c;

}

Color Sampler2DImp::sample_bilinear(Texture& tex, 
                                    float u, float v, 
                                    int level) {
  
  // Task 6: Implement bilinear filtering
    float w, h, x, y;
    float x_min, x_max, y_min, y_max;
    float x_frac, y_frac;
    float empty = 0;
    w = tex.mipmap[level].width;
    h = tex.mipmap[level].width;
    x = u * w;
    y = v * h;
    x_min = floor(x);
    x_max = ceil(x);
    y_min = floor(y);
    y_max = ceil(y);
    x_frac = std::modff(x, &empty);
    y_frac = std::modff(y, &empty);
    Color c1, c2, c3, c4;

    c1 = retrieveColor(tex.mipmap[level].texels, x_min + y_min * w);
    c2 = retrieveColor(tex.mipmap[level].texels, x_min + y_max * w);
    c3 = retrieveColor(tex.mipmap[level].texels, x_max + y_max * w);
    c4 = retrieveColor(tex.mipmap[level].texels, x_max + y_min * w);


    return (c1 * (1 - x_frac) + c4 * x_frac) * (1 - y_frac) + (c2 * (1 - x_frac) + c3 * x_frac) * y_frac;


  // return magenta for invalid level
  return Color(1,0,1,1);

}

Color Sampler2DImp::sample_trilinear(Texture& tex, 
                                     float u, float v, 
                                     float u_scale, float v_scale) {

  // Task 7: Implement trilinear filtering

  float levelf = max(0.f, ((float)tex.mipmap.size() - 1) - log2f(min(u_scale, v_scale)));
  float empty = 0;
  float alpha = modf(levelf, &empty);
  int floorLevel = floor(levelf);
  int ceilLevel = ceil(levelf);
  Color c1 = sample_bilinear(tex, u, v, floorLevel);
  Color c2 = sample_bilinear(tex, u, v, ceilLevel);

  return c1 * (1 - alpha) + c2 * alpha;


  // return magenta for invalid level
  return Color(1,0,1,1);

}

} // namespace CMU462
