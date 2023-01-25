#include "Box.h"
#include "Generators.hh"


double scale[DIM];
static bool scale_initialized = false; 
Box::Box() {
  if (!scale_initialized) {
    scale_initialized = true;
    for (int i = 0; i < DIM; ++i) {
      scale[i] = pow(2, -i / float(DIM));
    }
  }
  for (int i = 0; i < DIM; ++i) {
    center_digits[i] = 0;
    size_digits[i] = SCL;
  }
  pos = 0;
  compute_center_and_size();
  // compute_nearer();
  compute_cover();
}

Box Box::child(int dir) const
{
  Box child(*this);
  child.size_digits[pos] *= 0.5;
  child.center_digits[pos] += (2 * dir - 1) * child.size_digits[pos];
  ++child.pos;
  if (child.pos == DIM) { child.pos = 0; }

  child.name = name;
  child.name.append(1, '0'+dir);

  child.qr = qr;
  child.short_words_cache.clear();

  child.compute_center_and_size();
  // child.compute_nearer();
  child.compute_cover();
  return child;
}

Box get_box(std::string code) {
  Box box;
  for (char dir : code) {
    if (dir == '0') {
      box = box.child(0);
    } else if (dir == '1') {
      box = box.child(1);
    }
  }
  return box;
}

std::string Box::desc() {
  AJ sinhdx = _cover.sinhdx;
  AJ coshdx = _cover.coshdx;
  AJ expdx = _cover.expdx;
  AJ expmdx = _cover.expmdx;
  AJ sinhdy = _cover.sinhdy;
  AJ coshdy = _cover.coshdy;
  AJ expdy = _cover.expdy;
  AJ expmdy = _cover.expmdy;
  AJ coshmu = _cover.coshmu; 
  AJ cosf   = _cover.cosf;
  AJ sinf   = _cover.sinf;
  AJ expif  = _cover.expif;
  AJ expmif = _cover.expmif;
  AJ sintx2 = _cover.sintx2;
  AJ sinty2 = _cover.sinty2; 
  AJ coshlx = _cover.coshlx;
  AJ coshly = _cover.coshly;
  AJ coshLx2 = _cover.coshLx2;
  AJ sinhLx2 = _cover.sinhLx2;
  AJ coshLy2 = _cover.coshLy2;
  AJ sinhLy2 = _cover.sinhLy2;
  Complex c_sinhdx = _center.sinhdx;
  Complex c_sinhdy = _center.sinhdy;
  Complex c_coshmu = _center.coshmu; 
  Complex c_cosf   = _center.cosf;
  Complex c_sintx2 = _center.sintx2;
  Complex c_sinty2 = _center.sinty2; 

  char _desc[10000];
  sprintf(_desc, "%s\n", name.c_str());
  sprintf(_desc + strlen(_desc),
      "sinh(d_x) = %f with size %f, absLB %f, and absUB %f\n",
      sinhdx.f.re, sinhdx.size, absLB(sinhdx), absUB(sinhdx));
  sprintf(_desc + strlen(_desc),
      "cosh(d_x) = %f with size %f, absLB %f, and absUB %f\n",
      coshdx.f.re, coshdx.size, absLB(coshdx), absUB(coshdx));
  sprintf(_desc + strlen(_desc),
      "exp(d_x) = %f with size %f, absLB %f, and absUB %f\n",
      expdx.f.re, expdx.size, absLB(expdx), absUB(expdx));
  sprintf(_desc + strlen(_desc),
      "exp(-d_x) = %f with size %f, absLB %f, and absUB %f\n",
      expmdx.f.re, expmdx.size, absLB(expmdx), absUB(expmdx));
  sprintf(_desc + strlen(_desc),
      "sinh(d_y) =  %f with size %f, absLB %f, and absUB %f\n",
      sinhdy.f.re, sinhdy.size, absLB(sinhdy), absUB(sinhdy));
  sprintf(_desc + strlen(_desc),
      "cosh(d_y) =  %f with size %f, absLB %f, and absUB %f\n",
      coshdy.f.re, coshdy.size, absLB(coshdy), absUB(coshdy));
  sprintf(_desc + strlen(_desc),
      "exp(d_y) = %f with size %f, absLB %f, and absUB %f\n",
      expdy.f.re, expdy.size, absLB(expdy), absUB(expdy));
  sprintf(_desc + strlen(_desc),
      "exp(-d_y) = %f with size %f, absLB %f, and absUB %f\n",
      expmdy.f.re, expmdy.size, absLB(expmdy), absUB(expmdy));
  sprintf(_desc + strlen(_desc),
      "cosh(mu) = %f with size %f, absLB %f, and absUB %f\n",
      coshmu.f.re, coshmu.size, absLB(coshmu), absUB(coshmu));
  sprintf(_desc + strlen(_desc),
      "cos(phi) = %f with size %f, absLB %f, and absUB %f\n",
      cosf.f.re, cosf.size, absLB(cosf), absUB(cosf));
  sprintf(_desc + strlen(_desc),
      "sin(t_x/2) = %f with size %f, absLB %f, and absUB %f\n",
      sintx2.f.re, sintx2.size, absLB(sintx2), absUB(sintx2));
  sprintf(_desc + strlen(_desc),
      "sin(t_y/2) = %f with size %f, absLB %f, and absUB %f\n",
      sinty2.f.re, sinty2.size, absLB(sinty2), absUB(sinty2));
  sprintf(_desc + strlen(_desc),
      "sin(phi) = %f with size %f, absLB %f, and absUB %f\n",
      sinf.f.re, sinf.size, absLB(sinf), absUB(sinf));
  sprintf(_desc + strlen(_desc),
      "exp(i phi) = %f + i %f with size %f, absLB %f, and absUB %f\n",
      expif.f.re, expif.f.im, expif.size, absLB(expif), absUB(expif));
  sprintf(_desc + strlen(_desc),
      "exp(-i phi) = %f + i %f with size %f, absLB %f, and absUB %f\n",
      expmif.f.re, expmif.f.im, expmif.size, absLB(expmif), absUB(expmif));
  sprintf(_desc + strlen(_desc),
      "cosh(lx)) = %f + i %f with size %f, absLB %f, and absUB %f\n",
      coshlx.f.re, coshlx.f.im, coshlx.size, absLB(coshlx), absUB(coshlx));
  sprintf(_desc + strlen(_desc),
      "cosh(ly)) = %f + i %f with size %f, absLB %f, and absUB %f\n",
      coshly.f.re, coshly.f.im, coshly.size, absLB(coshly), absUB(coshly));
  sprintf(_desc + strlen(_desc),
      "sinh(Lx/2) = %f + i %f with size %f, absLB %f, and absUB %f\n",
      sinhLx2.f.re, sinhLx2.f.im, sinhLx2.size, absLB(sinhLx2), absUB(sinhLx2));
  sprintf(_desc + strlen(_desc),
      "cosh(Lx/2) = %f + i %f with size %f, absLB %f, and absUB %f\n",
      coshLx2.f.re, coshLx2.f.im, coshLx2.size, absLB(coshLx2), absUB(coshLx2));
  sprintf(_desc + strlen(_desc),
      "sinh(Ly/2) = %f + i %f with size %f, absLB %f, and absUB %f\n",
      sinhLy2.f.re, sinhLy2.f.im, sinhLy2.size, absLB(sinhLy2), absUB(sinhLy2));
  sprintf(_desc + strlen(_desc),
      "cosh(Ly/2) = %f + i %f with size %f, absLB %f, and absUB %f\n",
      coshLy2.f.re, coshLy2.f.im, coshLy2.size, absLB(coshLy2), absUB(coshLy2));

  sprintf(_desc + strlen(_desc),
  "Center\n"
"    sinh(d_x) %f | sinh(d_y) %f\n"
"    cosh(mu) %f | cos(phi) % f\n"
"    sin(t_x/2) %f | sin(t_y/2) %f\n",
  c_sinhdx.real(), c_sinhdy.real(), c_coshmu.real(),
  c_cosf.real(), c_sintx2.real(), c_sinty2.real() );

  std::string s(_desc);
  return s;
}

void Box::compute_center_and_size()
{
	for (int i = 0; i < 6; ++i) {
        // GMT paper page 419 of Annals
        // box_size guarantees that :
        // box_center - box_size <= true_center - true_size
        // box_center + box_size >= true_center + true_size
        // where box operations are floating point. 
        box_center[i] = scale[i] * center_digits[i];
        box_size[i]= (1 + 2 * EPS) * (size_digits[i] * scale[i] + 
            HALFEPS * fabs(center_digits[i]));
  }
  _center.sinhdx = Complex(box_center[0], 0);
  _center.sinhdy = Complex(box_center[1], 0);
  _center.coshmu = Complex(box_center[2], 0);
  _center.cosf = Complex(box_center[3], 0);
  _center.sintx2 = Complex(box_center[4], 0);
  _center.sinty2 = Complex(box_center[5], 0);

  fill_derived(_center);

  _x_center = construct_x(_center); 
  _y_center = construct_y(_center); 
}

void Box::compute_cover()
{
  // Let A = { (z0,z1,z2) \in C^3 | |zi| <= 1 }
  // Our parameters are functions on A with the following defintions
  // sinh(d_x) = s[0]re(z0) + c[0]
  // sinh(d_y) = s[1]im(z0) + c[1]
  // cosh(mu) = s[2]re[z1] + c[2]
  // cos(phi) = s[3]im[z1] + c[3]
  // sin(t_x/2) = s[4]re[z2] + c[4]
  // sin(t_y/2) = s[5]im[z2] + c[5]

  // Here, we have sinhdx(z0,z1,z2) = size * re(z0) + center and
  // sinhdy(z0,z1,z2) = size * im(z0) + center
  // TODO: verify that this multiplication by powers of 2 is valid
  _cover.sinhdx = AJ(XComplex(box_center[0], 0),
      XComplex(box_size[0] / 2, 0), 0, 0,
      XComplex(box_size[0] / 2, 0), 0, 0);
  // Note: d(im(z0))/dz0 = - i / 2 and d(im(z0)/dconj(z0) = i/2
  _cover.sinhdy = AJ(XComplex(box_center[1], 0),
      XComplex(0., - box_size[1] / 2), 0, 0,
      XComplex(0.,   box_size[1] / 2), 0, 0);

  // TODO: verify that division by 2 is valid
  _cover.coshmu = AJ(XComplex(box_center[2], 0),
      0, XComplex(box_size[2] / 2, 0), 0,
      0, XComplex(box_size[2] / 2, 0), 0);

  _cover.cosf = AJ(XComplex(box_center[3], 0),
      0, XComplex(0, -box_size[3] / 2), 0,
      0, XComplex(0,  box_size[3] / 2), 0);

  _cover.sintx2 = AJ(XComplex(box_center[4], 0),
      0, 0, XComplex(box_size[4] / 2, 0),
      0, 0, XComplex(box_size[4] / 2, 0));

  _cover.sinty2 = AJ(XComplex(box_center[5], 0),
      0, 0, XComplex(0, -box_size[5] / 2),
      0, 0, XComplex(0,  box_size[5] / 2));

  fill_derived(_cover);

  _x_cover = construct_x(_cover); 
  _y_cover = construct_y(_cover); 
}

