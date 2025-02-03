#include <stdio.h>
#include <math.h>
#include <time.h>

double ring_X, ring_Y;
double Bzz;
double pi=M_PI;
double amu = 931494061 ; // [eV/c^2]
double vc = 299792458 ; // [m/s]
double z = 32. ; // z number
double mass = 77.922853 ; // mass
double brho0 = 4.7447 ; // [Tm] :4.7447 = 168MeV/u, 4.848(MS03)
double brho; //
double main0 = 1915 ; // 1915A base (you can change this value)
double B0 = 1.182166; //[T] -1.182166 for Br=4.7447(168MeV/u), 1.2079037(MS03)

//double B0 = -1.182166; closed value for Brho = 4.7447 Tm 78Ge 100fs step
double gam ;
double beta ;

double MAT[11][11];

double caltmag(double xpos, double ypos, double mag_x, double mag_y, double mag_rad){
  
  double ape, len, dy;
  double Bz ;
  double rad ;

  double x1, x2, y1, y2;
  double xpos1,xpos2,xpos3,xpos4,xpos5,xpos6,xpos7,xpos8;
  double ypos1,ypos2,ypos3,ypos4,ypos5,ypos6;
  double dy1, dy2, dy3, dy4, dy5;

  double Btrim[10];
  int itri;

  ape = 80. ; //[mm] matnet's aperture
  len = 525. ; //[mm] half length of magnetic field in TARNII
  /////////////////////////////////////////

  // x,y rotation
  x1 = (xpos-mag_x); // trans for x -> x'
  y1 = (ypos-mag_y); // trans for y -> y'
  
  x2 = x1*cos(mag_rad) - y1*sin(mag_rad) ; // rot x' -> x''
  y2 = x1*sin(mag_rad) + y1*cos(mag_rad) ; // rot y' -> y''

  // cal position //
  
  xpos1 = x2/1000.;
  
  ypos1 = y2;
  
  dy = (fabs(ypos1)-len)/ape ;

  printf("xpos1 is %14.10lf \n",xpos1);
  printf("ypos1 is %14.10lf \n",ypos1);
  
}

double main(void){

  double xpos;
  double ypos;
  double mag_x;
  double mag_y;
  double mag_rad;

  xpos = 9287.959673;
  ypos = 1700.0;

  mag_x = 9232.420563;
  mag_y = 2536.110167;

  mag_rad = 7.5;
  mag_rad = mag_rad / 180*pi;
  printf("mag_rad is %17.16lf \n",mag_rad);

  caltmag(xpos,ypos,mag_x,mag_y,-mag_rad);
}
