#include <stdio.h>
#include "SL2.hh"
#include "Box.h"
#include "IsomH3.hh"
#include "AJCC.h"
#include "types.hh"
#include "TubeSearch.hh" 
#include "TestCollection.hh"

#define MAX_DEPTH 200

using namespace std;

void test_AJCC() {
  AJCC one =  AJCC(1,0,0,0,0,0);
  AJCC zero = AJCC(0,0,0,0,0,0);
  AJCC x = AJCC(XComplex(3,4),0,0,0,0);
  AJCC y = AJCC(6,XComplex(4,3),0,0,0);
  AJCC z = AJCC(15,0,XComplex(5,12),0,0,0);
  AJCC w = AJCC(10,0,0,XComplex(0,1),0,0);
  print_type("1+1", one + one); 
  print_type("1/0", one / zero); 
  print_type("abs(3+4i)", abs(x)); 
  print_type("abs(3+4i)^2", abs_sqrd(x)); 
  print_type("conj(6+(4+3i)z0", conj(y)); 
  print_type("conj(31+(4+3i)z0+(5+12i)z1+(0+1i)z2", conj(y+z+w)); 
  print_type(y*z); 
  print_type(z*(z+y)); 
  print_type(abs(z*2)); 
  print_type(abs(z*w)); 
  print_type(x+y); 
}

int main(int argc,char**argv)
{
    if(argc != 2) {
        fprintf(stderr,"Usage: %s position < data\n", argv[0]);
        exit(1);
    }
    char where[MAX_DEPTH];
    size_t depth = 0;
    while (argv[1][depth] != '\0') {
        if (argv[1][depth] != '0' && argv[1][depth] != '1'){
            fprintf(stderr,"bad position %s\n",argv[1]);
            exit(2);
        }
        where[depth] = argv[1][depth];
        depth++;
    }
    where[depth] = '\0';
  test_AJCC();
  SL2<Complex> M = SL2<Complex>(1,Complex(7,5),0,1);
  printf("%f + i %f\n", (M*M).b.real(), (M*M).b.imag());
  printf("%f + i %f\n", inverse(M).b.real(), inverse(M).b.imag());
  AJCC one =  AJCC(1,0,0,0,0,0);
  AJCC zero = AJCC(0,0,0,0,0,0);
  AJCC t = abs(AJCC(XComplex(3,4),0,0,0,0));
  printf("%f + i %f with err %f\n", t.f.re, t.f.im, t.e);
  SL2<AJCC> N = SL2<AJCC>(one,AJCC(XComplex(7,5),0,0,0,0),zero,one);
  printf("%f + i %f with err %f\n", (N*N).b.f.re, (N*N).b.f.im, (N*N).b.e );
  printf("%f + i %f with err %f\n", inverse(N).b.f.re, inverse(N).b.f.im, inverse(N).b.e );
  Params<AJCC> p;
  p.sinhD2 = AJCC(1,1,0,1,0,0);
  printf("%f + i %f\n", p.sinhD2.f.re, p.sinhD2.e);

  /*char code[] = "";
  char code[] = "0110111001001001000111001001001001001010110111000110"; 
  char code[] = "011011100100100100011100100100100100101011011100011001"; 
  char code[] = "011011100100100100011100100100100100101011011100011001111010001011001000"; 
  char code[] = "011011100100100100011100100100100100101011011100011001111010001011001000000101000111";
 char code[] = "011011100100100100011100100100100100101011011100011001111010001011001000000101000111110001110001000100010000001100100010";
  char code[] = "011011100100100100011100100100100100101011011100011001111010001011001000000101000111";
 char code[] = "01101110010010010001110010010010010101100010101011001011010001010111100010100110101100101010011011011000101101010010001001011";*/
  char code[] = "011111100000100100100000011100000110110001000111111101001111101101001011000101010011";

	Box box;
	for (char* dir = (char *) &where; *dir; ++dir) {
    // printf("%c\n", *dir);
		if (*dir == '0') {
			box = box.child(0);
		} else if (*dir == '1') {
			box = box.child(1);
		}
  }

  Params<AJCC> params = box.cover();
  AJCC sinhD2 = params.sinhD2;
  AJCC coshD2 = params.coshD2;

  printf("Box: %s", box.desc().c_str());

  SL2<AJCC> x = construct_x(params);
  SL2<AJCC> y = construct_y(params);
  SL2<Complex> c_x = construct_x(box.center());
  SL2<Complex> c_y = construct_y(box.center());

  pair<Complex, Complex> c_fcm = four_cosh_margulis_simple(c_x,c_y);
  pair<AJCC, AJCC> fcm = four_cosh_margulis_simple(x,y);

  printf("x.a is  %f + %f I with size %f, absLB %f, and absUB %f\n", x.a.f.re, x.a.f.im, x.a.size, absLB(x.a), absUB(x.a));
  printf("y.a is  %f + %f I with size %f, absLB %f, and absUB %f\n", y.a.f.re, y.a.f.im, y.a.size, absLB(y.a), absUB(y.a));
  print_type("re_len(x)", four_cosh_re_length(x));
  print_type("re_len(y)", four_cosh_re_length(x));
  printf("re_len(p) LB %f\n", absLB(coshD2 + sinhD2));
  printf("re_len(p) UB %f\n", absUB(coshD2 + sinhD2));

  printf("Margulis xy between %f and %f\n", absLB(fcm.first), absUB(fcm.first));
  printf("Exp(2t) xy between %f and %f\n", absLB(fcm.second), absUB(fcm.second));

  printf("%s\n",repeat("abc",5).c_str());
  SL2<Complex> f = construct_word("YxYYxxyy", box.center());
  SL2<AJCC> g = construct_word("YxYYxxyy", box.cover());
  SL2<AJCC> x1 = construct_word("x", box.cover());
  SL2<AJCC> y1 = construct_word("y", box.cover());
  SL2<AJCC> xy = construct_word("xy", box.cover());
  SL2<AJCC> xxxxx = construct_word("xxxxx", box.cover());
  printf("f.a is  %f + %f I\n", f.a.real(), f.a.imag());
  printf("g.a is  %f + %f I with size %f, absLB %f, and absUB %f\n", g.a.f.re, g.a.f.im, g.a.size, absLB(g.a), absUB(g.a));
  printf("x1.a is  %f + %f I with size %f, absLB %f, and absUB %f\n", x1.a.f.re, x1.a.f.im, x1.a.size, absLB(x1.a), absUB(x1.a));
  printf("y1.a is  %f + %f I with size %f, absLB %f, and absUB %f\n", y1.a.f.re, y1.a.f.im, y1.a.size, absLB(y1.a), absUB(y1.a));
  printf("xy.a is  %f + %f I with size %f, absLB %f, and absUB %f\n", xy.a.f.re, xy.a.f.im, xy.a.size, absLB(xy.a), absUB(xy.a));
  printf("(x*y).a is  %f + %f I with size %f, absLB %f, and absUB %f\n", (x*y).a.f.re, (x*y).a.f.im, (x*y).a.size, absLB((x*y).a), absUB((x*y).a));
  printf("xxxxx.a is  %f + %f I with size %f, absLB %f, and absUB %f\n", xxxxx.a.f.re, xxxxx.a.f.im, xxxxx.a.size, absLB(xxxxx.a), absUB(xxxxx.a));
  printf("(x*x*x*x*x).a is  %f + %f I with size %f, absLB %f, and absUB %f\n", (x*x*x*x*x).a.f.re, (x*x*x*x*x).a.f.im, (x*x*x*x*x).a.size, absLB((x*x*x*x*x).a), absUB((x*x*x*x*x).a));

  vector<string> known;
  known.push_back("XXYXXY");
  vector<word_pair> pairs = find_pairs(box.center(), known, 50, 20, vector<string>());  
 
  vector<word_pair>::iterator it;
  for (it = pairs.begin(); it != pairs.end(); ++it) {
    printf("(%s, %s)\n", (*it).first.c_str(), (*it).second.c_str());
  }  
}
