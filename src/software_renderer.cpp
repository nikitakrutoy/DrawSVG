#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <deque>

#include "triangulation.h"

using namespace std;

namespace CMU462 {


// Implements SoftwareRenderer //

void SoftwareRendererImp::draw_svg( SVG& svg ) {

  // set top level transformation
  transformation = svg_2_screen;

  // draw all elements
  for ( size_t i = 0; i < svg.elements.size(); ++i ) {
    draw_element(svg.elements[i]);
  }

  // draw canvas outline
  Vector2D a = transform(Vector2D(    0    ,     0    )); a.x--; a.y--;
  Vector2D b = transform(Vector2D(svg.width,     0    )); b.x++; b.y--;
  Vector2D c = transform(Vector2D(    0    ,svg.height)); c.x--; c.y++;
  Vector2D d = transform(Vector2D(svg.width,svg.height)); d.x++; d.y++;

  rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
  rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
  rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
  rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

  // resolve and send to render target
  resolve();

}

void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {

  // Task 4: 
  // You may want to modify this for supersampling support
  this->sample_rate = sample_rate;
  if (sample_rate > 1){
      supersample_target.clear();
      ss_w = target_w * sample_rate;
      ss_h = target_h * sample_rate;
      supersample_target.resize(ss_h * ss_w * 4);
      std::fill(supersample_target.begin(), supersample_target.end(), 255);
  }
}

void SoftwareRendererImp::set_render_target( unsigned char* render_target,
                                             size_t width, size_t height ) {

  // Task 4: 
  // You may want to modify this for supersampling support
  this->render_target = render_target;
  this->target_w = width;
  this->target_h = height;

}

void SoftwareRendererImp::draw_element( SVGElement* element ) {

  // Task 5 (part 1):
  // Modify this to implement the transformation stack
  transformation = transformation * element->transform;

  switch(element->type) {
    case POINT:
      draw_point(static_cast<Point&>(*element));
      break;
    case LINE:
      draw_line(static_cast<Line&>(*element));
      break;
    case POLYLINE:
      draw_polyline(static_cast<Polyline&>(*element));
      break;
    case RECT:
      draw_rect(static_cast<Rect&>(*element));
      break;
    case POLYGON:
      draw_polygon(static_cast<Polygon&>(*element));
      break;
    case ELLIPSE:
      draw_ellipse(static_cast<Ellipse&>(*element));
      break;
    case IMAGE:
      draw_image(static_cast<Image&>(*element));
      break;
    case GROUP:
      draw_group(static_cast<Group&>(*element));
      break;
    default:
      break;
  }
    transformation = transformation * element->transform.inv();

}


// Primitive Drawing //

void SoftwareRendererImp::draw_point( Point& point ) {

  Vector2D p = transform(point.position);
  rasterize_point( p.x, p.y, point.style.fillColor );

}

void SoftwareRendererImp::draw_line( Line& line ) { 

  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line( p0.x, p0.y, p1.x, p1.y, line.style.strokeColor );

}

void SoftwareRendererImp::draw_polyline( Polyline& polyline ) {

  Color c = polyline.style.strokeColor;

  if( c.a != 0 ) {
    int nPoints = polyline.points.size();
    for( int i = 0; i < nPoints - 1; i++ ) {
      Vector2D p0 = transform(polyline.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polyline.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_rect( Rect& rect ) {

  Color c;
  
  // draw as two triangles
  float x = rect.position.x;
  float y = rect.position.y;
  float w = rect.dimension.x;
  float h = rect.dimension.y;

  Vector2D p0 = transform(Vector2D(   x   ,   y   ));
  Vector2D p1 = transform(Vector2D( x + w ,   y   ));
  Vector2D p2 = transform(Vector2D(   x   , y + h ));
  Vector2D p3 = transform(Vector2D( x + w , y + h ));
  
  // draw fill
  c = rect.style.fillColor;
  if (c.a != 0 ) {
    rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    rasterize_triangle( p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c );
  }

  // draw outline
  c = rect.style.strokeColor;
  if( c.a != 0 ) {
    rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    rasterize_line( p1.x, p1.y, p3.x, p3.y, c );
    rasterize_line( p3.x, p3.y, p2.x, p2.y, c );
    rasterize_line( p2.x, p2.y, p0.x, p0.y, c );
  }

}

void SoftwareRendererImp::draw_polygon( Polygon& polygon ) {

  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if( c.a != 0 ) {

    // triangulate
    vector<Vector2D> triangles;
    triangulate( polygon, triangles );

    // draw as triangles
    for (size_t i = 0; i < triangles.size(); i += 3) {
      Vector2D p0 = transform(triangles[i + 0]);
      Vector2D p1 = transform(triangles[i + 1]);
      Vector2D p2 = transform(triangles[i + 2]);
      rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    }
  }

  // draw outline
  c = polygon.style.strokeColor;
  if( c.a != 0 ) {
    int nPoints = polygon.points.size();
    for( int i = 0; i < nPoints; i++ ) {
      Vector2D p0 = transform(polygon.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polygon.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_ellipse( Ellipse& ellipse ) {

  // Extra credit 

}

void SoftwareRendererImp::draw_image( Image& image ) {

  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image( p0.x, p0.y, p1.x, p1.y, image.tex );
}

void SoftwareRendererImp::draw_group( Group& group ) {

  for ( size_t i = 0; i < group.elements.size(); ++i ) {
    draw_element(group.elements[i]);
  }

}

// Rasterization //

// The input arguments in the rasterization functions 
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point( float x, float y, Color color ) {

  // fill in the nearest pixel
  int sx = (int) floor(x);
  int sy = (int) floor(y);
  if (sample_rate > 1) {
      if ( sx < 0 || sx >= ss_w ) return;
      if ( sy < 0 || sy >= ss_h ) return;

      // fill sample - NOT doing alpha blending!
      supersample_target[4 * (sx + sy * ss_w)    ] = (uint8_t) (color.r * 255);
      supersample_target[4 * (sx + sy * ss_w) + 1] = (uint8_t) (color.g * 255);
      supersample_target[4 * (sx + sy * ss_w) + 2] = (uint8_t) (color.b * 255);
      supersample_target[4 * (sx + sy * ss_w) + 3] = (uint8_t) (color.a * 255);
      return;
  }

  // check bounds
  if ( sx < 0 || sx >= target_w ) return;
  if ( sy < 0 || sy >= target_h ) return;

  // fill sample - NOT doing alpha blending!
  render_target[4 * (sx + sy * target_w)    ] = (uint8_t) (color.r * 255);
  render_target[4 * (sx + sy * target_w) + 1] = (uint8_t) (color.g * 255);
  render_target[4 * (sx + sy * target_w) + 2] = (uint8_t) (color.b * 255);
  render_target[4 * (sx + sy * target_w) + 3] = (uint8_t) (color.a * 255);

}

void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                          float x1, float y1,
                                          Color color) {
    x0 *= sample_rate;
    x1 *= sample_rate;
    y0 *= sample_rate;
    y1 *= sample_rate;

    float width = 0;
    float dx = x1 - x0;
    float dy = y1 - y0;
    float empty;

    float x,y;

    float width2 = width / std::cos(std::atan(abs(dx/dy))) / 2;
    float width3 = width / std::cos(std::atan(abs(dy/dx))) / 2;

    if (abs(dx) > abs(dy)) {
        if (x1 < x0) {
            std::swap(x0, x1);
            std::swap(y0, y1);
        }
        for (int i = 0; i < abs(dx); i++) {
            x = x0 + i;
            y = dy/dx * i + y0;
            rasterize_point(x, y - width3 - 1, color * ( 1 - std::modf(y - width3, &empty)) );
            for (float j = y - width3; j < y + width3; j = j + 1) {
                rasterize_point(x, j, color);
            }
            rasterize_point(x, y + width3, color * std::modf(y + width3, &empty));
        }
    }
    else {
        if (y1 < y0) {
            std::swap(x0, x1);
            std::swap(y0, y1);
        }
        for (int i = 0; i < abs(dy); i++) {
            y = y0 + i;
            x = dx/dy * i + x0;
            rasterize_point(x - width2 - 1, y, color * ( 1 - std::modf(x - width2, &empty)) );
            for (float j = x - width2; j < x + width2; j = j + 1) {
                rasterize_point(j, y, color);
            }
            rasterize_point(x + width2, y, color * std::modf(x + width2, &empty));

        }
    }
}

double length(double x, double y) {
    return sqrt(x * x + y * y);
}

double area(double s, double a, double b, double c) {
    return std::sqrt(s * (s - a) * (s - b) * (s - c));
}

float sign (Point p1, Point p2, Point p3)
{
    Vector2D v1, v2, v3;
    v1 = p1.position;
    v2 = p2.position;
    v3 = p3.position;
    return (v1.x - v3.x) * (v2.y - v3.y) - (v2.x - v3.x) * (v1.y - v3.y);
}

bool pointTriangleAreaTest(Point pt, Point p1, Point p2, Point p3, double eps = 1e-3) {
    double a, b, c;
    double q, w, e;
    double P1, P2, P3;
    double S1, S2, S3;
    double S, P;

    Vector2D v1, v2, v3, v;
    v1 = p1.position;
    v2 = p2.position;
    v3 = p3.position;
    v = pt.position;

    a = length(v1.x - v2.x, v1.y - v2.y);
    b = length(v2.x - v3.x, v2.y - v3.y);
    c = length(v3.x - v1.x, v3.y - v1.y);
    P = (a + b + c) / 2;
    S = area(P, a, b, c);

    q = length(v1.x - v.x, v1.y - v.y);
    w = length(v2.x -v.x, v2.y - v.y);
    e = length(v3.x - v.x, v3.y - v.y);

    P1 = (a + q + w) / 2;
    P2 = (b + w + e) / 2;
    P3 = (c + e + q) / 2;
    S1 = area(P1, a, q, w);
    S2 = area(P2, b, w, e);
    S3 = area(P3, c, e, q);
    return (abs(S  - (S1 + S2 + S3)) < eps);
}

bool pointTriangleTest (Point pt, Point v1, Point v2, Point v3)
{
    float d1, d2, d3;
    bool has_neg, has_pos;

    d1 = sign(pt, v1, v2);
    d2 = sign(pt, v2, v3);
    d3 = sign(pt, v3, v1);

    has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

    return !(has_neg && has_pos);
}

bool onSegment(Point p, Point q, Point r)
{
    Vector2D v1, v2, v3;
    v1 = p.position;
    v2 = q.position;
    v3 = r.position;
    if (v2.x <= max(v1.x, v3.x) && v2.x >= min(v1.x, v3.x) &&
        v2.y <= max(v1.y, v3.y) && v2.y >= min(v1.y, v3.y))
        return true;

    return false;
}

int orientation(Point p, Point q, Point r)
{
    Vector2D v1, v2, v3;
    v1 = p.position;
    v2 = q.position;
    v3 = r.position;
    int val = (v2.y - v1.y) * (v3.x - v2.x) -
              (v2.x - v1.x) * (v3.y - v2.y);

    if (val == 0) return 0;  // colinear

    return (val > 0)? 1: 2; // clock or counterclock wise
}

bool segmentsTest(Point p1, Point q1, Point p2, Point q2)
{
    // Find the four orientations needed for general and
    // special cases
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    // General case
    if (o1 != o2 && o3 != o4)
        return true;

    // Special Cases
    // p1, q1 and p2 are colinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment(p1, p2, q1)) return true;

    // p1, q1 and q2 are colinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;

    // p2, q2 and p1 are colinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;

    // p2, q2 and q1 are colinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;

    return false; // Doesn't fall in any of the above cases
}
struct BoundingBox {
    double x_min, x_max, y_min, y_max;

    BoundingBox GetLU() {
        return {x_min, (x_max + x_min) / 2, y_min, (y_max + y_min) / 2};
    };

    BoundingBox GetRU() {
        return {(x_max + x_min) / 2, x_max, y_min, (y_max + y_min) / 2};
    };
    BoundingBox GetLD() {
        return {x_min, (x_max + x_min) / 2, (y_max + y_min) / 2, y_max};
    };
    BoundingBox GetRD() {
        return {(x_max + x_min) / 2, x_max, (y_max + y_min) / 2, y_max};
    };
};

bool boxPointTest(BoundingBox box, Point p1, Point p2, Point p3) {
    Point bp1, bp2, bp3, bp4;
    bp1.position = Vector2D(box.x_min, box.y_min);
    bp2.position = Vector2D(box.x_min, box.y_max);
    bp3.position = Vector2D(box.x_max, box.y_max);
    bp4.position = Vector2D(box.x_max, box.y_min);
    Point boxPoints[4] = {bp1, bp2, bp3, bp4};
    Point trianglePoints[3] = {p1, p2, p3};

    for (auto bp: boxPoints) {
        if (pointTriangleTest(bp, p1, p2, p3))
            return true;
    }

    Point a1, a2, b1, b2;
    for (int i = 0 ; i < 4; i++){
        a1 = boxPoints[i];
        a2 = boxPoints[(i + 1) % 4];
        for (int j = 0; j < 3; j++) {
            b1 = trianglePoints[i];
            b2 = trianglePoints[(i + 1) % 3];
            if (segmentsTest(a1, a2, b1, b2))
                return true;
        }
    }

    return false;
}

float MIN_BOX_SIZE = 100000;

void SoftwareRendererImp::rasterize_triangle( float x0, float y0,
                                              float x1, float y1,
                                              float x2, float y2,
                                              Color color ) {
    x0 *= sample_rate;
    x1 *= sample_rate;
    x2 *= sample_rate;
    y0 *= sample_rate;
    y1 *= sample_rate;
    y2 *= sample_rate;

    float b1, b2, a1, a2;
    b2 = std::max({x0, x1, x2});
    b1 = std::min({x0, x1, x2});
    a2 = std::max({y0, y1, y2});
    a1 = std::min({y0, y1, y2});

    Point p1, p2, p3, p4;
    double x,y;

    p1 = Point();
    p1.position = Vector2D(x0, y0);
    p2 = Point();
    p2.position = Vector2D(x1, y1);
    p3 = Point();
    p3.position = Vector2D(x2, y2);

    BoundingBox initial = {b1, b2, a1, a2};
    BoundingBox box;
    std::deque<BoundingBox> boundingBoxes({initial});

    while (!boundingBoxes.empty()) {
        box = boundingBoxes[0];
        boundingBoxes.pop_front();
        if ((box.x_max - box.x_min) < MIN_BOX_SIZE) {
            for (int i = box.x_min; i < box.x_max; i++) {
                for (int j = box.y_min; j < box.y_max; j++) {
                    x = i + 0.5f;
                    y = j + 0.5f;
                    p4 = Point();
                    p4.position = Vector2D(x, y);

                    bool inTriangle = pointTriangleTest(p4, p1, p2, p3);
                    if (inTriangle) {
                        rasterize_point(x, y, color);
                    }
                }
            }
        } else {
            if (boxPointTest(box, p1, p2, p3)) {
                boundingBoxes.push_back(box.GetLD());
                boundingBoxes.push_back(box.GetLU());
                boundingBoxes.push_back(box.GetRD());
                boundingBoxes.push_back(box.GetRU());
            }
        }
    }

}

void SoftwareRendererImp::rasterize_image( float x0, float y0,
                                           float x1, float y1,
                                           Texture& tex ) {
//    Sampler2DImp sampler = Sampler2DImp();
    float u, v;
    for (int i = x0; i <= x1; i++) {
        for (int j = y0; j < y1; j++) {
            u = abs(((float)i - x0) / (x1 - x0));
            v = abs(((float)j - y0) / (y1 - y0));
            rasterize_point(i, j, sampler->sample_trilinear(
                    tex, u, v, abs(x1 - x0), abs(y1 - y0)));
        }
    }
}

// resolve samples to render target
void SoftwareRendererImp::resolve( void ) {

  // Task 4: 
  // Implement supersampling
  // You may also need to modify other functions marked with "Task 4".
    if (sample_rate <= 1) return;
    size_t index;
    uint16_t aggr_r, aggr_g, aggr_b, aggr_a;
    for (int i = 0; i < target_w; i++) {
        for ( int j = 0; j < target_h; j++) {
            aggr_a = 0; aggr_b = 0; aggr_g = 0; aggr_r = 0;
            for (int k = 0; k < sample_rate; k++) {
                for ( int m = 0; m < sample_rate; m++) {
                    index = sample_rate * i + k + (sample_rate * j + m) * ss_w;
                    aggr_r += supersample_target[4 * index];
                    aggr_g += supersample_target[4 * index + 1];
                    aggr_b += supersample_target[4 * index + 2];
                    aggr_a += supersample_target[4 * index + 3];
                }
            }
            render_target[4 * (i + j * target_w)    ] = aggr_r / (sample_rate * sample_rate);
            render_target[4 * (i + j * target_w)  + 1] = aggr_g / (sample_rate * sample_rate);
            render_target[4 * (i + j * target_w)  + 2] = aggr_b / (sample_rate * sample_rate);
            render_target[4 * (i + j * target_w)  + 3] = aggr_a / (sample_rate * sample_rate);
        }
        supersample_target.clear();
    }




}


} // namespace CMU462
