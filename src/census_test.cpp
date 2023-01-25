#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "rapidcsv.h"
#include "SL2.hh"
#include "Box.h"
#include "IsomH3.hh"
#include "AJ.h"
#include "roundoff.h"
#include "types.hh"
#include "TubeSearch.hh"
#include "TestCollection.hh"
#include "Refine.hh"

#define ERR 0.000001
#define CERR 0.00001
#define DERR 0.00005
#define FERR 0.00020
#define AJERR 0.0009

using namespace std;

extern Options g_options;
double g_cosh_marg_upper_bound = 3.0;
double g_cosh_marg_lower_bound = 1.0054;
double g_sinh_d_bound = 10.0; 

bool g_debug = false;

bool g_symmetric = false;

void refine_tree(Box box, PartialTree& t);

const XComplex to_XComplex(const Complex& z) {
  return XComplex(z.real(), z.imag());
}

const AJ to_AJ(const Complex& z) {
  return AJ(to_XComplex(z));
}

int main(int argc,char**argv)
{
  if(argc != 2) {
    fprintf(stderr,"Usage: %s census_csv_file\n", argv[0]);
    exit(1);
  }

  initialize_roundoff();

  rapidcsv::Document doc(argv[1], rapidcsv::LabelParams(0,-1));
  const size_t row_count = doc.GetRowCount();
  for (size_t i = 0; i < row_count; ++i) {
    string name = doc.GetCell<string>("name", i);
    string codes = doc.GetCell<string>("box codes", i);
    double sinhdx = doc.GetCell<double>("sinhdx", i);
    double sinhdy = doc.GetCell<double>("sinhdy", i);
    double coshmu = doc.GetCell<double>("coshmu", i);
    double cosf = doc.GetCell<double>("cosf", i);
    double sintx2 = doc.GetCell<double>("sintx2", i);
    double sinty2 = doc.GetCell<double>("sinty2", i);

    Complex Lx = parse_complex(doc.GetCell<string>("x", i)); 
    Complex Ly = parse_complex(doc.GetCell<string>("y", i)); 
    // Lx = shift_imag_around_zero(Lx);
    // Ly = shift_imag_around_zero(Ly);

    vector<string> box_codes;
    split_string(codes, "\"[',] ", box_codes);

    for (string box_code : box_codes) {
      Box box = get_box(box_code);
      Box bigbox = get_box(box_code.substr(0,60));
      Params<Complex> center = box.center();
      Params<AJ> cover = box.cover();
      /*
         printf("%s %f %f %f %f %f %f %s\n", name.c_str(), sinhdx, sinhdy, coshmu, cosf, sintx2, sinty2, box_code.c_str());
         printf("Box: %s", box.desc().c_str());
         printf("Manifold %s\n", name.c_str());
         printf("Lx: %f + i %f\n", Lx.real(), Lx.imag());
         Complex sinhLx2 = sinh(Lx/2);
         Complex coshLx2 = cosh(Lx/2);
         printf("sinhLx2: %f + i %f\n", sinhLx2.real(), sinhLx2.imag());
         printf("coshLx2: %f + i %f\n", coshLx2.real(), coshLx2.imag());
         printf("center coshlx: %f + i %f\n", center.coshlx.real(), center.coshlx.imag());
         printf("center sinhlx2: %f + i %f\n", center.sinhlx2.real(), center.sinhlx2.imag());
         printf("center coshlx2: %f + i %f\n", center.coshlx2.real(), center.coshlx2.imag());
         printf("center sinhLx2: %f + i %f\n", center.sinhLx2.real(), center.sinhLx2.imag());
         printf("center coshLx2: %f + i %f\n", center.coshLx2.real(), center.coshLx2.imag());
       */

      // printf("%s\n", name.c_str());

      // ***** Center testing *****
      assert(absUB(center.sinhdx - sinhdx) < ERR);
      assert(absUB(center.sinhdy - sinhdy) < ERR);
      assert(absUB(center.coshmu - coshmu) < ERR);
      assert(absUB(center.cosf - cosf) < ERR);
      assert(absUB(center.sintx2 - sintx2) < ERR);
      assert(absUB(center.sinty2 - sinty2) < ERR);
      // printf("center coshLx2: %f + i %f vs true coshLx2: %f + i %f\n", center.coshLx2.real(), center.coshLx2.imag(),
      //                                                            cosh(Lx/2).real(), cosh(Lx/2).imag());
      assert(absUB(center.coshLx2 - cosh(Lx/2)) < CERR);
      assert(absUB(center.coshLy2 - cosh(Ly/2)) < CERR);

      // Test margulis compuations
      pair<Complex, Complex> center_pair = four_cosh_margulis_simple(box.x_center(), box.y_center());
      // Complex four_cosh_mu = center.coshmu * 4;
      // printf("center 4coshmu: %f + i %f\n", four_cosh_mu.real(), four_cosh_mu.imag());
      // printf("Box: %s", box.desc().c_str());
      assert(absUB(center_pair.first - center.coshmu * 4) < CERR);
      assert(absUB(center_pair.second - center.expdx * center.expdx) < CERR);

      // Test construct math
      const SL2<Complex> x_center = box.x_center();          
      const SL2<Complex> y_center = box.y_center();          
      const SL2<Complex> xy_center = construct_word("xy", center); 
      const SL2<Complex> xyx_center = construct_word("xyx", center); 
      const SL2<Complex> yxxYXy_center = construct_word("yxxYXy", center); 
      assert(absUB(dist(xy_center, x_center * y_center)) < ERR);          
      assert(absUB(dist(xyx_center, x_center * y_center * x_center)) < ERR);          
      assert(absUB(dist(yxxYXy_center, y_center * x_center * x_center * inverse(y_center) * inverse(x_center) * y_center)) < ERR);          

      // Test Jorgensen
      SL2<Complex> c_id; 
      Complex j_ww_xy_center = jorgensen(x_center, y_center);
      Complex j_ww_yx_center = jorgensen(y_center, x_center);
      Complex j_xw_xy_center = jorgensen_xw(y_center, center);
      Complex j_wx_yx_center = jorgensen_wx(y_center, center);
      Complex j_wx_ix_center = jorgensen_wx(c_id, center);
      Complex j_wx_mix_center = jorgensen_wx(-c_id, center);
      Complex j_wy_xy_center = jorgensen_wy(x_center, center);
      Complex j_yw_yx_center = jorgensen_yw(x_center, center);
      Complex j_wy_iy_center = jorgensen_wy(c_id, center);
      Complex j_wy_miy_center = jorgensen_wy(-c_id, center);
      Complex j_xy_center = jorgensen_xy(center);
      Complex j_yx_center = jorgensen_yx(center);
      assert(absUB(j_ww_xy_center - j_xy_center) < ERR); 
      assert(absUB(j_xw_xy_center - j_xy_center) < ERR); 
      assert(absUB(j_wy_xy_center - j_xy_center) < ERR); 
      assert(absUB(j_ww_yx_center - j_yx_center) < ERR); 
      assert(absUB(j_wx_yx_center - j_yx_center) < ERR); 
      assert(absUB(j_yw_yx_center - j_yx_center) < ERR); 
      assert(absLB(j_xy_center) >= 1.0); 
      assert(absLB(j_yx_center) >= 1.0); 
      assert(absUB(j_wx_ix_center) < 1);
      assert(absUB(j_wx_mix_center) < 1);
      assert(absUB(j_wy_iy_center) < 1);
      assert(absUB(j_wy_miy_center) < 1);

      // Test moving axes
      const SL2<Complex> I_center;
      Complex no_move_center_ax_v1 = four_cosh_dist_ax_wax(I_center, center);
      Complex no_move_center_ay_v1 = four_cosh_dist_ay_way(I_center, center);
      Complex no_move_center_ax_v2 = four_cosh_dist_ax_wax(x_center, center);
      Complex no_move_center_ay_v2 = four_cosh_dist_ay_way(y_center, center);
      Complex ax_ay_dist_center_v1 = four_cosh_dist_ax_way(I_center, center);
      Complex ax_ay_dist_center_v2 = four_cosh_dist_ay_wax(I_center, center);
      Complex ax_ay_dist_center_v3 = four_cosh_dist_ay_wax(x_center, center);
      Complex ax_ay_dist_center_v4 = four_cosh_dist_ax_way(y_center, center);
      Complex ax_xay_dist_center = four_cosh_dist_ax_way(x_center, center);
      Complex ay_yax_dist_center = four_cosh_dist_ay_wax(y_center, center);
      Complex ax_yax_dist_center = four_cosh_dist_ax_wax(y_center, center);
      Complex ay_xay_dist_center = four_cosh_dist_ay_way(x_center, center);

      assert(absUB(no_move_center_ax_v1 - 4) < ERR);
      assert(absUB(no_move_center_ay_v1 - 4) < ERR);
      assert(absUB(no_move_center_ax_v2 - 4) < ERR);
      assert(absUB(no_move_center_ay_v2 - 4) < ERR);
      assert(absUB(ax_ay_dist_center_v1 - center.coshdxdy * 4) < ERR);
      assert(absUB(ax_ay_dist_center_v2 - center.coshdxdy * 4) < ERR);
      assert(absUB(ax_ay_dist_center_v3 - center.coshdxdy * 4) < ERR);
      assert(absUB(ax_ay_dist_center_v4 - center.coshdxdy * 4) < ERR);

      // print_type("ax_xay_dist_center", ax_xay_dist_center);
      // print_type("center.coshdxdy * 4", center.coshdxdy * 4);
      // print_type("diff", ax_xay_dist_center - center.coshdxdy * 4);
      assert(absUB(ax_xay_dist_center - center.coshdxdy * 4) < ERR);
      // print_type("ay_yax_dist_center", ay_yax_dist_center);
      // print_type("center.coshdxdy * 4", center.coshdxdy * 4);
      // print_type("diff", ay_yax_dist_center - center.coshdxdy * 4);
      assert(absUB(ay_yax_dist_center - center.coshdxdy * 4) < ERR);
      // print_type("ax_yax_dist_center", ax_yax_dist_center);
      // print_type("center.cosh2dx * 4", center.cosh2dx * 4);
      // print_type("diff", ax_yax_dist_center - center.cosh2dx * 4);
      assert(absLB(ax_yax_dist_center - center.cosh2dx * 4) > ERR);
      // print_type("ay_xay_dist_center", ay_xay_dist_center);
      // print_type("center.cosh2dy * 4", center.cosh2dy * 4);
      // print_type("diff", ay_xay_dist_center - center.cosh2dy * 4);
      assert(absLB(ay_xay_dist_center - center.cosh2dy * 4) > ERR);

      // TestCollection boolean verification
      const SL2<Complex> xxXYy_center = construct_word("xxXYy", center); 
      const SL2<Complex> yxXYy_center = construct_word("yxXYy", center); 
      assert(inside_var_nbd(x_center * inverse(x_center), x_center) == true);
      assert(inside_var_nbd(y_center * inverse(y_center), y_center) == true);
      assert(inside_var_nbd(x_center, y_center) == false);
      assert(inside_var_nbd(x_center, yxxYXy_center) == false);
      assert(inside_var_nbd_ne(x_center, y_center) == false);
      assert(inside_var_nbd_ne(x_center, yxxYXy_center) == false);
      assert(really_cant_fix_x_axis(x_center, center) == false);
      assert(cant_fix_x_axis(y_center, center) == true );
      assert(really_cant_fix_y_axis(y_center, center) == false);
      assert(cant_fix_y_axis(x_center, center) == true );
      assert(must_fix_x_axis(x_center, center) == true);
      assert(must_fix_x_axis(y_center, center) == false );
      assert(must_fix_y_axis(y_center, center) == true);
      assert(must_fix_y_axis(x_center, center) == false);
      assert(inside_var_nbd_x(x_center, center) == true);
      assert(inside_var_nbd_x(y_center, center) == false);
      assert(inside_var_nbd_x(xxXYy_center, center) == true);
      assert(inside_var_nbd_x(yxXYy_center, center) == false);
      assert(inside_var_nbd_y(y_center, center) == true);
      assert(inside_var_nbd_y(x_center, center) == false);
      assert(inside_var_nbd_y(yxXYy_center, center) == true);
      assert(inside_var_nbd_y(xxXYy_center, center) == false);

      // ***** Cover testing *****
      assert(absUB(cover.sinhdx - sinhdx) < ERR);
      assert(absUB(cover.sinhdy - sinhdy) < ERR);
      assert(absUB(cover.coshmu - coshmu) < ERR);
      assert(absUB(cover.cosf - cosf) < ERR);
      assert(absUB(cover.sintx2 - sintx2) < ERR);
      assert(absUB(cover.sinty2 - sinty2) < ERR);
      assert(absUB(cover.coshLx2 - to_AJ(cosh(Lx/2))) < AJERR);
      assert(absUB(cover.coshLy2 - to_AJ(cosh(Ly/2))) < AJERR);

      // Test margulis compuations
      pair<AJ, AJ> cover_pair = four_cosh_margulis_simple(box.x_cover(), box.y_cover());
      // print_type("cover margulis exact:", cover.coshmu * 4);
      // print_type("computed margulis:", cover_pair.first);
      // print_type("computed diff:", cover_pair.first - cover.coshmu * 4);
      assert(absUB(cover_pair.first - cover.coshmu * 4) < AJERR);
      assert(absUB(cover_pair.second - cover.expdx * cover.expdx) < AJERR);

      // Test construct math
      const SL2<AJ> x_cover = box.x_cover();          
      const SL2<AJ> y_cover = box.y_cover();          
      const SL2<AJ> xy_cover = construct_word("xy", cover); 
      const SL2<AJ> xyx_cover = construct_word("xyx", cover); 
      const SL2<AJ> yxxYXy_cover = construct_word("yxxYXy", cover); 
      //print_SL2(xy_cover); 
      //print_SL2(x_cover * y_cover);
      //print_type(dist(xy_cover, x_cover * y_cover)); 
      assert(absUB(dist(xy_cover, x_cover * y_cover)) < AJERR);          
      assert(absUB(dist(xyx_cover, x_cover * y_cover * x_cover)) < AJERR);          
      // print_SL2(yxxYXy_cover); 
      // print_SL2(y_cover * x_cover * x_cover * inverse(y_cover) * inverse(x_cover) * y_cover);
      // print_type("dist:", dist(yxxYXy_cover, y_cover * x_cover * x_cover * inverse(y_cover) * inverse(x_cover) * y_cover)); 
      assert(absUB(dist(yxxYXy_cover, y_cover * x_cover * x_cover * inverse(y_cover) * inverse(x_cover) * y_cover)) < AJERR * 2); // error gets too big here         

      // Test Jorgensen
      SL2<AJ> aj_id; 
      AJ j_ww_xy_cover = jorgensen(x_cover, y_cover);
      AJ j_ww_yx_cover = jorgensen(y_cover, x_cover);
      AJ j_xw_xy_cover = jorgensen_xw(y_cover, cover);
      AJ j_wx_yx_cover = jorgensen_wx(y_cover, cover);
      AJ j_wx_ix_cover = jorgensen_wx(aj_id, cover);
      AJ j_wx_mix_cover = jorgensen_wx(-aj_id, cover);
      AJ j_wy_xy_cover = jorgensen_wy(x_cover, cover);
      AJ j_yw_yx_cover = jorgensen_yw(x_cover, cover);
      AJ j_wy_iy_cover = jorgensen_wy(aj_id, cover);
      AJ j_wy_miy_cover = jorgensen_wy(-aj_id, cover);
      AJ j_xy_cover = jorgensen_xy(cover);
      AJ j_yx_cover = jorgensen_yx(cover);
      // print_type("j_ww_xy_cover:", j_ww_yx_cover);
      // print_type("j_xy_cover:", j_yx_cover);
      // print_type("diffr:", j_ww_yx_cover - j_yx_cover);
      assert(absUB(j_ww_xy_cover - j_xy_cover) < FERR); 
      assert(absUB(j_xw_xy_cover - j_xy_cover) < FERR); 
      assert(absUB(j_wy_xy_cover - j_xy_cover) < FERR); 
      assert(absUB(j_ww_yx_cover - j_yx_cover) < FERR); 
      assert(absUB(j_wx_yx_cover - j_yx_cover) < FERR); 
      assert(absUB(j_yw_yx_cover - j_yx_cover) < FERR); 
      assert(absLB(j_xy_cover) >= 1.0); 
      assert(absLB(j_yx_cover) >= 1.0); 
      assert(absUB(j_wx_ix_cover) < 1);
      assert(absUB(j_wx_mix_cover) < 1);
      assert(absUB(j_wy_iy_cover) < 1);
      assert(absUB(j_wy_miy_cover) < 1);

      // Test moving axes
      const SL2<AJ> I_cover;
      AJ no_move_cover_ax_v1 = four_cosh_dist_ax_wax(I_cover, cover);
      AJ no_move_cover_ay_v1 = four_cosh_dist_ay_way(I_cover, cover);
      AJ no_move_cover_ax_v2 = four_cosh_dist_ax_wax(x_cover, cover);
      AJ no_move_cover_ay_v2 = four_cosh_dist_ay_way(y_cover, cover);
      AJ ax_ay_dist_cover_v1 = four_cosh_dist_ax_way(I_cover, cover);
      AJ ax_ay_dist_cover_v2 = four_cosh_dist_ay_wax(I_cover, cover);
      AJ ax_ay_dist_cover_v3 = four_cosh_dist_ay_wax(x_cover, cover);
      AJ ax_ay_dist_cover_v4 = four_cosh_dist_ax_way(y_cover, cover);
      AJ ax_xay_dist_cover = four_cosh_dist_ax_way(x_cover, cover);
      AJ ay_yax_dist_cover = four_cosh_dist_ay_wax(y_cover, cover);
      AJ ax_yax_dist_cover = four_cosh_dist_ax_wax(y_cover, cover);
      AJ ay_xay_dist_cover = four_cosh_dist_ay_way(x_cover, cover);

      assert(absUB(no_move_cover_ax_v1 - 4) < ERR);
      assert(absUB(no_move_cover_ay_v1 - 4) < ERR);
      assert(absUB(no_move_cover_ax_v2 - 4) < ERR);
      assert(absUB(no_move_cover_ay_v2 - 4) < ERR);
      // print_type("ax_ay_dist_cover_v1", ax_ay_dist_cover_v1);
      // print_type("cover.coshdxdy * 4", cover.coshdxdy * 4);
      // print_type("diff", ax_ay_dist_cover_v1 - cover.coshdxdy * 4);
      assert(absUB(ax_ay_dist_cover_v1 - cover.coshdxdy * 4) < DERR);
      // print_type("ax_ay_dist_cover_v2", ax_ay_dist_cover_v2);
      // print_type("cover.coshdxdy * 4", cover.coshdxdy * 4);
      // print_type("diff", ax_ay_dist_cover_v2 - cover.coshdxdy * 4);
      assert(absUB(ax_ay_dist_cover_v2 - cover.coshdxdy * 4) < DERR);
      // print_type("ax_ay_dist_cover_v3", ax_ay_dist_cover_v3);
      // print_type("cover.coshdxdy * 4", cover.coshdxdy * 4);
      // print_type("diff", ax_ay_dist_cover_v3 - cover.coshdxdy * 4);
      assert(absUB(ax_ay_dist_cover_v3 - cover.coshdxdy * 4) < DERR);
      // print_type("ax_ay_dist_cover_v4", ax_ay_dist_cover_v4);
      // print_type("cover.coshdxdy * 4", cover.coshdxdy * 4);
      // print_type("diff", ax_ay_dist_cover_v4 - cover.coshdxdy * 4);
      assert(absUB(ax_ay_dist_cover_v4 - cover.coshdxdy * 4) < DERR);

      // print_type("ax_xay_dist_cover", ax_xay_dist_cover);
      // print_type("cover.coshdxdy * 4", cover.coshdxdy * 4);
      // print_type("diff", ax_xay_dist_cover - cover.coshdxdy * 4);
      assert(absUB(ax_xay_dist_cover - cover.coshdxdy * 4) < DERR);
      // print_type("ay_yax_dist_cover", ay_yax_dist_cover);
      // print_type("cover.coshdxdy * 4", cover.coshdxdy * 4);
      // print_type("diff", ay_yax_dist_cover - cover.coshdxdy * 4);
      assert(absUB(ay_yax_dist_cover - cover.coshdxdy * 4) < DERR);
      // print_type("ax_yax_dist_cover", ax_yax_dist_cover);
      // print_type("cover.cosh2dx * 4", cover.cosh2dx * 4);
      // print_type("diff", ax_yax_dist_cover - cover.cosh2dx * 4);
      assert(absLB(ax_yax_dist_cover - cover.cosh2dx * 4) > ERR);
      // print_type("ay_xay_dist_cover", ay_xay_dist_cover);
      // print_type("cover.cosh2dy * 4", cover.cosh2dy * 4);
      // print_type("diff", ay_xay_dist_cover - cover.cosh2dy * 4);
      assert(absLB(ay_xay_dist_cover - cover.cosh2dy * 4) > ERR);

      // TestCollection boolean verification
      const SL2<AJ> xxXYy_cover = construct_word("xxXYy", cover); 
      const SL2<AJ> yxXYy_cover = construct_word("yxXYy", cover); 
      assert(inside_var_nbd(x_cover * inverse(x_cover), x_cover) == true);
      assert(inside_var_nbd(y_cover * inverse(y_cover), y_cover) == true);
      assert(inside_var_nbd(x_cover, y_cover) == false);
      assert(inside_var_nbd(x_cover, yxxYXy_cover) == false);
      assert(inside_var_nbd_ne(x_cover, y_cover) == false);
      assert(inside_var_nbd_ne(x_cover, yxxYXy_cover) == false);
      assert(cant_fix_x_axis(x_cover, cover) == false);
      assert(cant_fix_x_axis(y_cover, cover) == true );
      assert(cant_fix_y_axis(y_cover, cover) == false);
      assert(cant_fix_y_axis(x_cover, cover) == true );
      assert(must_fix_x_axis(x_cover, cover) == true);
      assert(must_fix_x_axis(y_cover, cover) == false );
      assert(must_fix_y_axis(y_cover, cover) == true);
      assert(must_fix_y_axis(x_cover, cover) == false);
      assert(inside_var_nbd_x(x_cover, cover) == true);
      assert(inside_var_nbd_x(y_cover, cover) == false);
      assert(inside_var_nbd_x(xxXYy_cover, cover) == true);
      assert(inside_var_nbd_x(yxXYy_cover, cover) == false);
      assert(inside_var_nbd_y(y_cover, cover) == true);
      assert(inside_var_nbd_y(x_cover, cover) == false);
      assert(inside_var_nbd_y(yxXYy_cover, cover) == true);
      assert(inside_var_nbd_y(xxXYy_cover, cover) == false);

      vector<string> empty;
      vector<word_pair> found = find_words_v2(center, 1, 20, empty, map<string, int>());
      fprintf(stderr, "Word pairs found for %s at index %zu:\n", name.c_str(), i);
      for (auto pair : found) {
        fprintf(stderr, "    %s, %s\n", pair.first.c_str(), pair.second.c_str());
      }

      // Refine testing
      g_options.box_name = bigbox.name.c_str(); 
      g_options.words_file = "/dev/null"; 
      g_options.powers_file = "/dev/null"; 
      g_options.max_depth = 140; 
      g_options.invent_depth = 18; 
      g_options.improve_tree = false; 
      g_options.truncate_depth = 2; 
      g_options.max_size = 300000; 
      g_options.word_search_depth = 6; 
      g_options.fill_holes = true; 
      g_cosh_marg_upper_bound = 2.55; 
      g_sinh_d_bound = 10.0;

      
      PartialTree* t = new PartialTree();
      refine_tree(bigbox, *t);
      // print_tree(*t);
      fprintf(stderr, "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n"); 
      fflush(stdout);
      fflush(stderr);
    }
  }

  if(!roundoff_ok()){
    printf("Underflow may have occurred\n");
    exit(1);
  }

  printf("PASSED\n");

}
