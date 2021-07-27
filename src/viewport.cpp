#include "viewport.h"

#include "CMU462.h"

namespace CMU462 {

void ViewportImp::set_viewbox( float centerX, float centerY, float vspan ) {

  // Task 5 (part 2): 
  // Set svg coordinate to normalized device coordinate transformation. Your input
  // arguments are defined as normalized SVG canvas coordinates.
    this->centerX = centerX;
    this->centerY = centerY;
    this->vspan = vspan;


    double tm[9] = {
            1, 0, -centerX,
            0, 1, -centerY,
            0, 0 , 1
    };

    double sm[9] = {
            1 / (2 * vspan), 0, 0,
            0, 1 / (2 * vspan), 0,
            0, 0 , 1
    };

    double tm2[9] = {
            1, 0, 0.5,
            0, 1, 0.5,
            0, 0 , 1
    };


    Matrix3x3 translate  = Matrix3x3(tm);
    Matrix3x3 scale  = Matrix3x3(sm);
    Matrix3x3 translate2  = Matrix3x3(tm2);
//    Matrix3x3 scale = Matrix3x3::identity();
    set_svg_2_norm(translate2 * scale * translate);
}

void ViewportImp::update_viewbox( float dx, float dy, float scale ) { 
  
  this->centerX -= dx;
  this->centerY -= dy;
  this->vspan *= scale;
  set_viewbox( centerX, centerY, vspan );
}

} // namespace CMU462
