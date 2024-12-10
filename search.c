#include <stdio.h>
#include <math.h>
#include <time.h>

double ring_X, ring_Y;
double Bzz;
// double pi = M_PI;
double pi = 3.14159265;
double amu = 931494061; // [eV/c^2]
double vc = 299792458;  // [m/s]
double z = 32.;         // z number
double mass = 77.9229;  // mass
double brho0 = 4.7447;  // [Tm] :4.7447 = 168MeV/u, 4.848(MS03)
double brho;            //
double main0 = 1915;    // 1915A base (you can change this value)
double B0 = 1.182166;   //[T] -1.182166 for Br=4.7447(168MeV/u), 1.2079037(MS03)

// double B0 = -1.182166; closed value for Brho = 4.7447 Tm 78Ge 100fs step
double gam;
double beta;

double MAT[11][11];

double caltmag(double xpos, double ypos, double mag_x, double mag_y, double mag_rad, double mcur, double *tcur)
{

  double tw11, tw21, tw31, tw41, tw51, tw61, tw71, tw81, tw91;
  double en1, en2, en3, en4, en5, en6;
  double ape, len, dy;
  double Bz;
  double rad;

  double x1, x2, y1, y2;
  double xpos1, xpos2, xpos3, xpos4, xpos5, xpos6, xpos7, xpos8;
  double ypos1, ypos2, ypos3, ypos4, ypos5, ypos6;
  double dy1, dy2, dy3, dy4, dy5;

  double Btrim[10];
  int itri;

  // Parameters for TARNII magnet //
  tw11 = 1.000 * (mcur / main0); //    parameter of w11
  tw21 = 0.000;                  //    parameter of w21
  tw31 = 0.0639773;              // parameter of w31
  tw41 = 0.000;                  // parameter of w41
  tw51 = -54.1635;               // parameter of w51
  tw61 = 0.000;                  // parameter of w61
  tw71 = 12402.1;                // parameter of w71
  tw81 = 0.000;                  // parameter of w81
  tw91 = -651257;                // parameter of w91

  // Enge function parameters for TARNII magnet
  en1 = 0.288396;  //    parameter of a1
  en2 = 1.41928;   //    parameter of a2
  en3 = -0.319549; // parameter of a3
  en4 = 0.226168;  // parameter of a4
  en5 = -0.026047; // parameter of a5
  en6 = 0.002317;  // parameter of a6

  ape = 80.;  //[mm] matnet's aperture
  len = 525.; //[mm] half length of magnetic field in TARNII
  /////////////////////////////////////////

  // Trim function parameters at 200 A with 4-gaussian //
  double a1[10] = {-0.00559819, -0.004793, -0.002201, -0.005188, -0.023787, 0.0019743, 0.0097702, 0.0037986, 0.0051218, 0.0067407};
  double b1[10] = {-126.476, -137.427, -157.189, -111.732, -132.117, 188.603, 161.154, 168.821, 155.302, 111.382};
  double c1[10] = {30.9606, 29.6729, 25.2112, 36.8752, 33.6369, 24.8122, 62.3175, 27.9871, 30.3747, 34.6009};
  double off1[10] = {10.0687, 10.1344, 10.2894, 9.96632, 10.0629, 9.23832, 9.8653, 9.80403, 9.93002, -16.0584};
  double a2[10] = {-0.0000670702, -0.0000886371, -0.001263, -0.001419, -0.00235, 0.0019491, -0.000999232, 0.0000549654, 0.000169799, -0.000154026};
  double b2[10] = {-11.2285, -3.66384, -21.5678, -13.3005, -6.34961, 29.0324, 19.2292, 17.3227, 34.4530, 614.0470};
  double c2[10] = {-13.9661, -13.79960, -18.9073, -19.4676, -22.2167, -21.0937, 20.25420, -11.50580, -17.4388, -521.1790};
  double off2[10] = {9.48442, 4.08676, 16.9878, 11.0316, 5.98738, -4.80779, 3.408610, 5.13617, -6.46249, -420.573};
  double a3[10] = {-0.000657458, -0.00475243, -0.0046398, -0.00244726, 0.0180907, 0.0256026, 0.00363535, 0.00239795, 0.000811862, 0.000703519};
  double b3[10] = {-95.6711, -82.0093, -117.715, -158.686, -129.043, 111.263, 66.6887, 74.4467, 81.6906, 86.7074};
  double c3[10] = {-14.8655, -29.2736, -32.4186, -26.2377, -30.9078, -39.5128, -34.7429, -22.0430, -16.4470, -15.3450};
  double off3[10] = {10.3050, 10.3313, 10.7113, 11.1389, 12.9681, 7.609, 10.6119, 10.4978, 10.5619, -0.472402};
  double a4[10] = {-0.00426336, -0.00252437, -0.00483284, -0.00463077, -0.00533062, -0.0193996, -0.00586014, 0.0051618, 0.00418408, -0.000792968};
  double b4[10] = {-66.6196, -44.9389, -58.3225, -45.0715, -46.9839, 112.976, 244.062, 114.386, 103.053, 90.6008};
  double c4[10] = {23.9585, 21.5422, 32.4758, 32.927, 36.0638, 34.7987, 103.475, 32.8979, 25.963, 15.0778};
  double off4[10] = {10.1855, 9.67604, 9.93579, 9.66259, 9.5711, 9.26375, 6.11265, 9.07634, 9.37191, -37.2705};
  ///////////////////////////////////////////////////////

  // x,y rotation
  x1 = (xpos - mag_x); // trans for x -> x'
  y1 = (ypos - mag_y); // trans for y -> y'

  x2 = x1 * cos(mag_rad) - y1 * sin(mag_rad); // rot x' -> x''
  y2 = x1 * sin(mag_rad) + y1 * cos(mag_rad); // rot y' -> y''

  // cal position //

  xpos1 = x2 / 1000.;
  xpos2 = xpos1 * xpos1;
  xpos3 = xpos2 * xpos1;
  xpos4 = xpos3 * xpos1;
  xpos5 = xpos4 * xpos1;
  xpos6 = xpos5 * xpos1;
  xpos7 = xpos6 * xpos1;
  xpos8 = xpos7 * xpos1;

  ypos1 = y2;
  ypos2 = ypos1 * ypos1;
  ypos3 = ypos2 * ypos1;
  ypos4 = ypos3 * ypos1;
  ypos5 = ypos4 * ypos1;
  ypos6 = ypos5 * ypos1;

  dy = (fabs(ypos1) - len) / ape;
  dy1 = dy;
  dy2 = dy1 * dy1;
  dy3 = dy2 * dy1;
  dy4 = dy3 * dy1;
  dy5 = dy4 * dy1;

  for (itri = 0; itri < 10; itri++)
  {
    Btrim[itri] = (a1[itri] * exp(-((x2 + off1[itri] - b1[itri]) * (x2 + off1[itri] - b1[itri])) / (2 * c1[itri] * c1[itri]))) + (a2[itri] * exp(-((x2 + off2[itri] - b2[itri]) * (x2 + off2[itri] - b2[itri])) / (2 * c2[itri] * c2[itri]))) + (a3[itri] * exp(-((x2 + off3[itri] - b3[itri]) * (x2 + off3[itri] - b3[itri])) / (2 * c3[itri] * c3[itri]))) + (a4[itri] * exp(-((x2 + off4[itri] - b4[itri]) * (x2 + off4[itri] - b4[itri])) / (2 * c4[itri] * c4[itri])));
  }

  if (fabs(xpos1) < 0.2 && fabs(ypos1) < 1000)
  {
    Bz = B0 * (tw11 + tw21 * xpos1 + tw31 * xpos2 + tw41 * xpos3 + tw51 * xpos4 + tw61 * xpos5 + tw71 * xpos6 + tw81 * xpos7 + tw91 * xpos8);
    for (itri = 0; itri < 10; itri++)
    {
      Bz = Bz + Btrim[itri] * tcur[itri] / 200;
    }
    Bz = Bz * (1. / (1 + exp(en1 + en2 * dy1 + en3 * dy2 + en4 * dy3 + en5 * dy4 + en6 * dy5)));
    Bz = Bz * -1;
  }
  else
  {
    Bz = 0.;
  }

  return Bz;
}

double calfmag(double xpos, double ypos, double mag_x, double mag_y, double mag_rad, double mcur)
{

  double fw11, fw21, fw31, fw41, fw51, fw61, fw71, fw81, fw91;
  double en1, en2, en3, en4, en5, en6;
  double ape, len, dy;
  double Bz;
  double rad;

  double x1, x2, y1, y2;
  double xpos1, xpos2, xpos3, xpos4, xpos5, xpos6, xpos7, xpos8;
  double ypos1, ypos2, ypos3, ypos4;
  double dy1, dy2, dy3, dy4, dy5;

  fw11 = 1.000 * (mcur / main0); //    parameter of w11
  fw21 = 0.000;                  //    parameter of w21
  fw31 = 0.0639773;              // parameter of w31
  fw41 = 0.000;                  // parameter of w41
  fw51 = -54.1635;               // parameter of w51
  fw61 = 0.000;                  // parameter of w61
  fw71 = 12402.1;                // parameter of w71
  fw81 = 0.000;                  // parameter of w81
  fw91 = -651257;                // parameter of w91

  // Enge function parameters for TARNII magnet
  en1 = 0.288396;  //    parameter of a1
  en2 = 1.41928;   //    parameter of a2
  en3 = -0.319549; // parameter of a3
  en4 = 0.226168;  // parameter of a4
  en5 = -0.026047; // parameter of a5
  en6 = 0.002317;  // parameter of a6

  ape = 80.;  //[mm] matnet's aperture
  len = 525.; //[mm] half length of magnetic field in TARNII
  /////////////////////////////////////////

  // x,y rotation
  x1 = (xpos - mag_x);
  y1 = (ypos - mag_y);

  x2 = x1 * cos(mag_rad) - y1 * sin(mag_rad);
  y2 = x1 * sin(mag_rad) + y1 * cos(mag_rad);

  // cal position //

  xpos1 = x2 / 1000.;
  xpos2 = xpos1 * xpos1;
  xpos3 = xpos2 * xpos1;
  xpos4 = xpos3 * xpos1;
  xpos5 = xpos4 * xpos1;
  xpos6 = xpos5 * xpos1;
  xpos7 = xpos6 * xpos1;
  xpos8 = xpos7 * xpos1;

  ypos1 = y2;
  ypos2 = ypos1 * ypos1;
  ypos3 = ypos2 * ypos1;
  ypos4 = ypos3 * ypos1;

  dy = (fabs(ypos1) - len) / ape;
  dy1 = dy;
  dy2 = dy1 * dy1;
  dy3 = dy2 * dy1;
  dy4 = dy3 * dy1;
  dy5 = dy4 * dy1;

  if (fabs(xpos1) < 0.2 && fabs(ypos1) < 1000)
  {
    Bz = B0 * (fw11 + fw21 * xpos1 + fw31 * xpos2 + fw41 * xpos3 + fw51 * xpos4 + fw61 * xpos5 + fw71 * xpos6 + fw81 * xpos7 + fw91 * xpos8) * (1. / (1 + exp(en1 + en2 * dy1 + en3 * dy2 + en4 * dy3 + en5 * dy4 + en6 * dy5)));
    Bz = Bz * -1;
  }
  else
  {
    Bz = 0.;
  }

  return Bz;
}

void getmag(double xpos, double ypos, double *fac, double mcur, double *tcur)
{
  double bm1xb, bm1yb;
  double bm2xb, bm2yb;
  double bm3xb, bm3yb;
  double bm4xb, bm4yb;

  double bm1xc[6], bm1yc[6];
  double rad1[6];

  double bm2xc[6], bm2yc[6];
  double rad2[6];

  double bm3xc[6], bm3yc[6];
  double rad3[6];

  double bm4xc[6], bm4yc[6];
  double rad4[6];

  int i;

  Bzz = 0; // initialize for Bzz

  xpos = xpos * 1000;
  ypos = ypos * 1000;

  // 1st sector //
  bm1xc[0] = 9232.420563; // BM11 magnet center X
  bm1yc[0] = 2536.110167; // BM11 magnet center Y
  bm2xc[0] = 8805.253915; // BM12 magnet center X
  bm2yc[0] = 4130.317798; // BM12 magnet center Y
  bm3xc[0] = 7980.153668; // BM13 magnet center X
  bm3yc[0] = 5559.485814; // BM13 magnet center Y
  bm4xc[0] = 6812.990034; // BM14 magnet center X
  bm4yc[0] = 6726.686242; // BM14 magnet center Y

  // 2nd sector //
  bm1xc[1] = 2419.652338;  // BM21 magnet center X
  bm1yc[1] = 9263.181119;  // BM21 magnet center Y
  bm2xc[1] = 823.734813;   // BM22 magnet center X
  bm2yc[1] = 9677.359839;  // BM22 magnet center Y
  bm3xc[1] = -824.801225;  // BM23 magnet center X
  bm3yc[1] = 9690.373999;  // BM23 magnet center Y
  bm4xc[1] = -2419.208113; // BM24 magnet center X
  bm4yc[1] = 9263.181119;  // BM24 magnet center Y

  // 3rd sector //
  bm1xc[2] = -6812.546113; // BM31 magnet center X
  bm1yc[2] = 6726.686242;  // BM31 magnet center Y
  bm2xc[2] = -7979.587097; // BM32 magnet center X
  bm2yc[2] = 5559.645258;  // BM32 magnet center Y
  bm3xc[2] = -8804.809690; // BM33 magnet center X
  bm3yc[2] = 4130.317798;  // BM33 magnet center Y
  bm4xc[2] = -9231.976338; // BM34 magnet center X
  bm4yc[2] = 2536.101670;  // BM34 magnet center Y

  // 3rd sector //
  bm1xc[3] = -9231.976338; // BM41 magnet center X
  bm1yc[3] = -2536.879587; // BM41 magnet center Y
  bm2xc[3] = -8804.809690; // BM42 magnet center X
  bm2yc[3] = -4131.087218; // BM42 magnet center Y
  bm3xc[3] = -7979.587097; // BM43 magnet center X
  bm3yc[3] = -5560.414678; // BM43 magnet center Y
  bm4xc[3] = -6812.546113; // BM44 magnet center X
  bm4yc[3] = -6727.455662; // BM44 magnet center Y

  // 5th sector //
  bm1xc[4] = -2419.208113; // BM51 magnet center X
  bm1yc[4] = -9263.950539; // BM51 magnet center Y
  bm2xc[4] = -825.000481;  // BM52 magnet center X
  bm2yc[4] = -9691.117187; // BM52 magnet center Y
  bm3xc[4] = 825.444706;   // BM53 magnet center X
  bm3yc[4] = -9691.117187; // BM53 magnet center Y
  bm4xc[4] = 2419.652338;  // BM54 magnet center X
  bm4yc[4] = -9263.950539; // BM54 magnet center Y

  // 6th sector //
  bm1xc[5] = 6812.990338;  // BM61 magnet center X
  bm1yc[5] = -6727.455662; // BM61 magnet center Y
  bm2xc[5] = 7980.031322;  // BM62 magnet center X
  bm2yc[5] = -5560.414678; // BM62 magnet center Y
  bm3xc[5] = 8805.253915;  // BM63 magnet center X
  bm3yc[5] = -4131.087218; // BM63 magnet center Y
  bm4xc[5] = 9232.420563;  // BM64 magnet center X
  bm4yc[5] = -2536.879587; // BM64 magnet center Y

  for (i = 0; i < 6; i++)
  {
    rad1[i] = 7.5 + i * 60.;      // bending rad [deg]
    rad1[i] = rad1[i] / 180 * pi; // degrees -> radians

    rad2[i] = 22.5 + i * 60.;     // bending rad [deg]
    rad2[i] = rad2[i] / 180 * pi; // degrees -> radians

    rad3[i] = 37.5 + i * 60.;     // bending rad [deg]
    rad3[i] = rad3[i] / 180 * pi; // degrees -> radians

    rad4[i] = 52.5 + i * 60.;     // bending rad [deg]
    rad4[i] = rad4[i] / 180 * pi; // degrees -> radians

    Bzz = Bzz + caltmag(xpos, ypos, bm1xc[i], bm1yc[i], -rad1[i], mcur, tcur);
    Bzz = Bzz + calfmag(xpos, ypos, bm2xc[i], bm2yc[i], -rad2[i], mcur);
    Bzz = Bzz + calfmag(xpos, ypos, bm3xc[i], bm3yc[i], -rad3[i], mcur);
    Bzz = Bzz + caltmag(xpos, ypos, bm4xc[i], bm4yc[i], -rad4[i], mcur, tcur);
  }

  *fac = (z * vc) / (mass * amu * gam) * Bzz; // q/m*B0

  //  printf("%lf %lf %lf \n", xpos, ypos, Bzz);
}

double inv_mat()
{
  double MAT[11][11];
  double INV[11][11];
  double buf;
  int i, j, k;
  int n = 11;

  FILE *fp1;
  fp1 = fopen("matrix.dat", "r");
  for (i = 0; i < n; i++)
  {
    fscanf(fp1, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
           &MAT[i][0], &MAT[i][1], &MAT[i][2], &MAT[i][3], &MAT[i][4], &MAT[i][5], &MAT[i][6], &MAT[i][7], &MAT[i][8], &MAT[i][9], &MAT[i][10]); //
  }
  // check reading matrix
  for (i = 0; i < 11; i++)
  {
    for (j = 0; j < 11; j++)
    {
      printf("  %3.4f", MAT[i][j]);
    }
    printf("\n");
  }
  //

  // 単位行列を作る
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      INV[i][j] = (i == j) ? 1.0 : 0.0;
    }
  }
  // 掃き出し法
  for (i = 0; i < n; i++)
  {
    buf = 1 / MAT[i][i];
    for (j = 0; j < n; j++)
    {
      MAT[i][j] *= buf;
      INV[i][j] *= buf;
    }
    for (j = 0; j < n; j++)
    {
      if (i != j)
      {
        buf = MAT[j][i];
        for (k = 0; k < n; k++)
        {
          MAT[j][k] -= MAT[i][k] * buf;
          INV[j][k] -= INV[i][k] * buf;
        }
      }
    }
  }

  // 逆行列を出力
  FILE *fp2;
  fp2 = fopen("inv_matrix.dat", "w");

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      fprintf(fp2, "  %lf", INV[i][j]);
    }
    fprintf(fp2, "\n");
  }

  fclose(fp2);
}

int main(void)
{

  double xinit, yinit, ainit, dp;
  double xmax, amax, xbin, abin;
  double x1, x2, y1, y2;
  double t1, trev, tp;
  double vx1, vx2, vy1, vy2;
  double t, ts, dt, x, y;
  double kx1, kx2, kx3, kx4;
  double ky1, ky2, ky3, ky4;
  double hx1, hx2, hx3, hx4;
  double hy1, hy2, hy3, hy4;
  double x0, y0;
  double kx, ky, hx, hy;
  int i, j, k, N, emi, l, m, cur;
  int pt = 0;
  double xp, yp;
  double mom, energy, velocity;
  double fac;
  double betax, disp;
  double qmr;

  double T0p[11], T1p[11];
  double dT0p[11];
  double Tideal;
  double Tide[11];

  double INV[11][11];

  double mcur;
  double tcur[10];
  double ccur[11];

  char filename1[20];
  char filename2[20];

  double radii;
  double xini;

  double ratio;

  FILE *fp1;
  FILE *fp2;

  ///// time stamp /////
  time_t timer;
  struct tm *t_st;
  //////////////////////

  /////////////////////////
  /////////////////////////
  //  Main part is here  //
  /////////////////////////
  /////////////////////////

  ///// For make mag field plot /////
#if 0 // 0 is off, 1 is on
  
  FILE *fp;
  int xi, yi;
  double Bcal;
  double xpos, ypos;
  fp =  fopen("field.dat","w");
  
  for(xi=0; xi<2500; xi++){
    xpos = xi*0.05 - 12 ; // -12 m to 12 m with 50 mm step
    for(yi=0; yi<2500; yi++){
      ypos = yi*0.05 - 12 ; // -12 m to 12 m with 50 mm step
      
      beta = ( (vc*z*brho0) / (amu*mass) )  / sqrt( ((vc*z*brho0)/(amu*mass))*((vc*z*brho0)/(amu*mass))+1 ) ;
      gam = 1/sqrt(1-beta*beta);

      getmag(xpos,ypos,&fac, mcur, tcur);

      Bcal = (mass*amu*gam)/(z*vc)*fac ;// fac -> Bzz
      fprintf(fp,"%lf %lf %lf \n", xpos*1000, ypos*1000, -Bcal);
      
    }
    fprintf(fp,"\n");
  }
  printf("making field.dat was finished! \n");
  fclose(fp);
#endif
  //////////////////////////////////

  for (cur = 0; cur < 12; cur++)
  {
    printf("cur=%d \n", cur);
    pt = 0; // for init

    // Magnet's current setting is here //

    // MS03 settings //
    /*    mcur = 1915 ;
    tcur[0] = 151.66 ;
    tcur[1] = 184.19 ;
    tcur[2] = 194.08 ;
    tcur[3] = 116.91 ;
    tcur[4] = 208.32 ;
    tcur[5] = 119.90 ;
    tcur[6] = 214.56 ;
    tcur[7] = 155.69 ;
    tcur[8] = 168.81 ;
    tcur[9] = 280.00 ;  */

    Tideal = 378.9362;
    Tide[0] = Tideal * -5E-4 + Tideal;
    Tide[1] = Tideal * -4E-4 + Tideal;
    Tide[2] = Tideal * -3E-4 + Tideal;
    Tide[3] = Tideal * -2E-4 + Tideal;
    Tide[4] = Tideal * -1E-4 + Tideal;
    Tide[5] = Tideal * 1E-7 + Tideal;
    Tide[6] = Tideal * 1E-4 + Tideal;
    Tide[7] = Tideal * 2E-4 + Tideal;
    Tide[8] = Tideal * 3E-4 + Tideal;
    Tide[9] = Tideal * 4E-4 + Tideal;
    Tide[10] = Tideal * 5E-4 + Tideal;

#if 0
    mcur = 1914.58425 ; // main coil current 
    tcur[0] = 340.479722 ; // trim1 coil current 
    tcur[1] = 91.612366 ; // trim2 coil current 
    tcur[2] = 175.160265 ; // trim3 coil current 
    tcur[3] = 156.878397 ; // trim4 coil current 
    tcur[4] = 156.174888 ; // trim5 coil current 
    tcur[5] = 168.368011 ; // trim6 coil current 
    tcur[6] = 145.012729 ; // trim7 coil current 
    tcur[7] = 195.804176 ; // trim8 coil current 
    tcur[8] = 82.643376 ; // trim9 coil current 
    tcur[9] = 331.386579 ; // trim10 coil current 


mcur = 1914.814 ; // main coil current 
tcur[0] = 474.270 ; // trim1 coil current 
tcur[1] = 67.384 ; // trim2 coil current 
tcur[2] = 145.490 ; // trim3 coil current 
tcur[3] = 143.873 ; // trim4 coil current 
tcur[4] = 192.132 ; // trim5 coil current 
tcur[5] = 165.825 ; // trim6 coil current 
tcur[6] = 75.015 ; // trim7 coil current 
tcur[7] = 307.038 ; // trim8 coil current 
tcur[8] = 45.267 ; // trim9 coil current 
tcur[9] = 274.157 ; // trim10 coil current

mcur = 1914.844 ; // main coil current 
tcur[0] = 379.373 ; // trim1 coil current 
tcur[1] = 168.060 ; // trim2 coil current 
tcur[2] = 127.289 ; // trim3 coil current 
tcur[3] = 123.859 ; // trim4 coil current 
tcur[4] = 202.347 ; // trim5 coil current 
tcur[5] = 165.613 ; // trim6 coil current 
tcur[6] = 74.179 ; // trim7 coil current 
tcur[7] = 309.536 ; // trim8 coil current 
tcur[8] = 31.472 ; // trim9 coil current 
tcur[9] = 293.750 ; // trim10 coil current 



mcur = 1918.2725 ; // main coil current 
tcur[0] = 361.194 ; // trim1 coil current 
tcur[1] = 166.052 ; // trim2 coil current 
tcur[2] = 120.914 ; // trim3 coil current 
tcur[3] = 139.903 ; // trim4 coil current 
tcur[4] = 199.534 ; // trim5 coil current 
tcur[5] = 144.721 ; // trim6 coil current 
tcur[6] = 105.309 ; // trim7 coil current 
tcur[7] = 296.118 ; // trim8 coil current 
tcur[8] = 18.772 ; // trim9 coil current 
tcur[9] = 309.170 ; // trim10 coil current 

mcur = 1918.045 ; // main coil current 
tcur[0] = 448.857 ; // trim1 coil current 
tcur[1] = 125.269 ; // trim2 coil current 
tcur[2] = 87.507 ; // trim3 coil current 
tcur[3] = 183.174 ; // trim4 coil current 
tcur[4] = 176.985 ; // trim5 coil current 
tcur[5] = 163.762 ; // trim6 coil current 
tcur[6] = 81.403 ; // trim7 coil current 
tcur[7] = 320.543 ; // trim8 coil current 
tcur[8] = 2.958 ; // trim9 coil current 
tcur[9] = 317.810 ; // trim10 coil current 

mcur = 1914.814 ; // main coil current 
tcur[0] = 474.270 ; // trim1 coil current 
tcur[1] = 67.384 ; // trim2 coil current 
tcur[2] = 145.490 ; // trim3 coil current 
tcur[3] = 143.873 ; // trim4 coil current 
tcur[4] = 192.132 ; // trim5 coil current 
tcur[5] = 165.825 ; // trim6 coil current 
tcur[6] = 75.015 ; // trim7 coil current 
tcur[7] = 307.038 ; // trim8 coil current 
tcur[8] = 45.267 ; // trim9 coil current 
tcur[9] = 274.157 ; // trim10 coil current

mcur = 1914.850 ; // main coil current 
tcur[0] = 533.844 ; // trim1 coil current 
tcur[1] = 19.655 ; // trim2 coil current 
tcur[2] = 143.376 ; // trim3 coil current 
tcur[3] = 149.495 ; // trim4 coil current 
tcur[4] = 196.442 ; // trim5 coil current 
tcur[5] = 166.015 ; // trim6 coil current 
tcur[6] = 68.272 ; // trim7 coil current 
tcur[7] = 316.741 ; // trim8 coil current 
tcur[8] = 29.007 ; // trim9 coil current 
tcur[9] = 291.791 ; // trim10 coil current 

mcur = 1918.21056175 ; // main coil current 
tcur[0] = 533.844 ; // trim1 coil current 
tcur[1] = 19.655 ; // trim2 coil current 
tcur[2] = 143.376 ; // trim3 coil current 
tcur[3] = 149.495 ; // trim4 coil current 
tcur[4] = 196.442 ; // trim5 coil current 
tcur[5] = 166.015 ; // trim6 coil current 
tcur[6] = 68.272 ; // trim7 coil current 
tcur[7] = 316.741 ; // trim8 coil current 
tcur[8] = 29.007 ; // trim9 coil current 
tcur[9] = 291.791 ; // trim10 coil current 

mcur = 1915.089050175 ; // main coil current 
tcur[0] = 500.844 ; // trim1 coil current 
tcur[1] = 19.655 ; // trim2 coil current 
tcur[2] = 120.376 ; // trim3 coil current 
tcur[3] = 200.495 ; // trim4 coil current 
tcur[4] = 5.442 ; // trim5 coil current 
tcur[5] = 134.015 ; // trim6 coil current 
tcur[6] = 190.272 ; // trim7 coil current 
tcur[7] = 141.741 ; // trim8 coil current 
tcur[8] = 49.007 ; // trim9 coil current 
tcur[9] = 256.791 ; // trim10 coil current
#endif

    ratio = 1.0016951;
    mcur = 1918.1205;  // main coil current
    tcur[0] = 592.128; // trim1 coil current
    tcur[1] = 1.958;   // trim2 coil current
    tcur[2] = 91.045;  // trim3 coil current
    tcur[3] = 226.784; // trim4 coil current
    tcur[4] = 158.111; // trim5 coil current
    tcur[5] = 146.100; // trim6 coil current
    tcur[6] = 117.213; // trim7 coil current
    tcur[7] = 280.012; // trim8 coil current
    tcur[8] = 47.809;  // trim9 coil current
    tcur[9] = 279.008; // trim10 coil current

    if (cur == 1)
    {
      mcur = mcur + 0.1; // main coil + 1.0 A
    }
    else if (cur > 1)
    {
      tcur[cur - 2] = tcur[cur - 2] + 10.0; // i-th trim coil + 1.0A
    }

    // "j" is momentum , "l" is emitannce, "k" is angle //
    for (j = -5; j < 6; j++)
    {
      for (l = 0; l < 1; l++)
      {
        for (k = 0; k < 1; k++)
        {
          for (m = -1; m < 1; m = m + 2)
          {
            N = 0;

            ///// ring parameters /////
            betax = 7.817; // [m] by MAD
            ///////////////////////////

            // Start position definition //
            ring_X = 9287.959673;
            ring_Y = 0.0;
            ///////////////////////////////

            /// initial information ///
            //  xinit = 0.0 ; // initial position x in [mm]
            //  ainit = 0.0 ; // initial angle in [mrad]
            //////////////////////////

            /// for emittance cal ///
            emi = l * 20;                                              // [pi mm mrad] 20 step
            amax = sqrt(emi / betax) * 0.95;                           // maximum angle for a certain emittance (add 0.95 to avoid error)
            abin = 2 * amax / 15.;                                     // angle step
            ainit = k * abin - amax;                                   //
            xinit = sqrt(betax * emi - betax * betax * ainit * ainit); //
            xinit = xinit * m;                                         //
            /////////////////////////

            /// initial momentum ///
            dp = 0.1 * j; // dp kizami in [%]
            brho = brho0 * (1 + dp / 100);
            ////////////////////////

            // ring dispersion //
            //	    disp = 7.0978 ; // [m] by MAD
            //	    disp = 0.1759*dp*dp*dp - 0.0107*dp*dp - 0.0734*dp + 7.0978 ; // [m] by ideal parameter
            disp = -3.3494 * dp * dp * dp * dp * dp * dp + 2.2088 * dp * dp * dp * dp * dp + 1.5699 * dp * dp * dp * dp - 0.731 * dp * dp * dp - 0.1573 * dp * dp + 0.0121 * dp + 7.0737; // [m] by ideal parameter
            /////////////////////

            ///// calculation for beam condition /////
            qmr = z / mass; // q/m  ratio
            mom = vc * z * brho * 1E-6;
            energy = (sqrt(mom * mom + (mass * amu) * (mass * amu)) - (mass * amu)) / mass;
            beta = ((vc * z * brho) / (amu * mass)) / sqrt(((vc * z * brho) / (amu * mass)) * ((vc * z * brho) / (amu * mass)) + 1);
            gam = 1 / sqrt(1 - beta * beta);
            velocity = vc * beta * 1E-6; // [mm/ns]
            //////////////////////////////////////////
            // printf("beta is %lf \n",beta);
            // printf("gamma is %lf \n",gam);
            // printf("velocity is %lf \n",velocity);
            // printf("%lf %lf %lf %lf \n", dp, brho, beta, velocity);
            if (m == -1)
            {
              sprintf(filename1, "kidou_emi_%d_dp_%1.2f_mx_%d.dat", emi, dp, k);
              sprintf(filename2, "emi_%d_dp_%1.2f_mx_%d.dat", emi, dp, k);
            }
            else if (m == 1)
            {
              sprintf(filename1, "kidou_emi_%d_dp_%1.2f_px_%d.dat", emi, dp, k);
              sprintf(filename2, "emi_%d_dp_%1.2f_px_%d.dat", emi, dp, k);
            }
            fp1 = fopen(filename1, "w");
            fp2 = fopen(filename2, "w");

            yinit = 0.0; // initial position y

            x1 = (xinit + ring_X + (dp * disp * 10)) / 1000.;
            xini = x1 * 1000;
            y1 = 0; // start plane = 0

            vx1 = beta * sin(ainit / 1000);
            vy1 = beta * cos(ainit / 1000);
            printf("x1 is %lf \n", x1);
            printf("y1 is %lf \n", y1);
            printf("vx1 is %lf \n", vx1);
            printf("vy1 is %lf \n", vy1);

            t = 0.0;
            t1 = 0.0;

            ts = 0.0001;           // in ns (0.0001 corresponds 100 fs step)
            dt = ts * (vc * 1E-9); // kizami haba 1 -> [ns]

            //  radii = amu*beta*gamma/qmr/vc/B0;
            //  printf("fac = %lf \n", fac);

            fprintf(fp1, "0 %lf %lf \n", x1 * 1000, y1 * 1000);

            //  for(i=0 ; i<3900000; i++){
            i = 0; // initialize
            while (N < 1)
            {
              if (-0.005 < y1 && y1 < 0.005 && x1 > 0)
              {
                dt = ts * (vc * 1E-9) * 1; // kizami haba 1 -> [ns]
                t = i * dt;                // 100 fs step near the start/end area
                i = i + 1;
              }
              else
              {
                dt = ts * 100 * (vc * 1E-9); // kizami haba 1 -> [ns]
                t = i * dt;                  // 10 ps step is standard
                i = i + 100;
              }

              kx1 = dt * vx1;
              ky1 = dt * vy1;
              getmag(x1 + kx1, y1 + ky1, &fac, mcur, tcur);
              // printf("x1 + kx1 is %lf \n", x1 + kx1);
              // printf("y1 + ky1 is %lf \n", y1 + ky1);
              // printf("tcur is %lf \n", tcur);
              hx1 = dt * vy1 * fac;
              hy1 = -dt * vx1 * fac;

              kx2 = dt * (vx1 + hx1 * 0.5);
              ky2 = dt * (vy1 + hy1 * 0.5);
              getmag(x1 + kx2, y1 + ky2, &fac, mcur, tcur);
              hx2 = dt * (vy1 + hy1 * 0.5) * fac;
              hy2 = dt * (-(vx1 + hx1 * 0.5) * fac);

              kx3 = dt * (vx1 + hx2 * 0.5);
              ky3 = dt * (vy1 + hy2 * 0.5);
              getmag(x1 + kx3, y1 + ky3, &fac, mcur, tcur);
              hx3 = dt * (vy1 + hy2 * 0.5) * fac;
              hy3 = dt * (-(vx1 + hx2 * 0.5) * fac);

              kx4 = dt * (vx1 + hx3);
              ky4 = dt * (vy1 + hy3);
              getmag(x1 + kx4, y1 + ky4, &fac, mcur, tcur);
              hx4 = dt * (vy1 + hy3) * fac;
              hy4 = dt * (-(vx1 + hx3) * fac);

              kx = (kx1 + 2 * kx2 + 2 * kx3 + kx4) / 6.;
              ky = (ky1 + 2 * ky2 + 2 * ky3 + ky4) / 6.;

              x2 = x1 + kx;
              y2 = y1 + ky;

              hx = (hx1 + 2 * hx2 + 2 * hx3 + hx4) / 6.;
              hy = (hy1 + 2 * hy2 + 2 * hy3 + hy4) / 6.;

              vx2 = vx1 + hx;
              vy2 = vy1 + hy;

              xp = x2 * 1000.; // m -> mm
              yp = y2 * 1000.; // m -> mm

              fprintf(fp1, "%d %lf %lf \n", i + 1, xp, yp);

              //    printf("%s now dp = %1.2f, i = %d, x = %lf, y = %lf  \n",ctime(&timer), dp, i, xp, yp);

              if (y1 < 0.0 && y2 > 0.0)
              {
                time(&timer);
                N = N + 1;
                tp = t / (vc * 1E-9); //  now time [ns]
                if (N == 1)
                {
                  trev = tp;
                }
                else
                {
                  trev = tp - t1; // revolution time [ns]
                }
                // Trun number, X position, Y position , Total time, revolution time, averaged time //
                fprintf(fp2, "%d %1.2f %lf %lf %lf %lf %lf \n", N, dp, xp - ring_X, yp, tp, trev, tp / N); // final turn
                printf("%s now dp = %1.2f,emi = %d %d, k = %d, i = %d/3900000 (%d%) \n", ctime(&timer), dp, emi, m, k, i, i / 39000);
                printf("dp = %1.2f , trev = %lf \n", dp, trev);
                printf("dx = %1.6f \n", xini - xp);
                T1p[pt] = trev; // get trev for each momentum
              } // 1周したときの処理

              x1 = x2;
              y1 = y2;
              vx1 = vx2;
              vy1 = vy2;
              t1 = tp;
            }

            fclose(fp1);
            fclose(fp2);
          } // m position pol
        } // k angle,position
      } // l emittance
      printf("dp = %1.2f,pt = %d,  Trev = %lf \n", dp, pt, T1p[pt]);
      pt = pt + 1;

    } // j momentum

    if (cur == 0)
    {
      for (i = 0; i < 11; i++)
      {
        T0p[i] = T1p[i]; // 基準となる電流セットにおける各運動量の周回時間
        dT0p[i] = Tide[i] - T0p[i];
        //	dT0p[i] = Tideal - T0p[i];
        printf("cur=%d, T1 = %lf, T0 = %lf, dT = %lf \n", cur, T1p[i], T0p[i], dT0p[i]);
      }
    }

    if (cur > 0)
    {
      for (i = 0; i < 11; i++)
      {
        MAT[i][cur - 1] = T1p[i] - T0p[i]; // 各要素(基準の時間からずれ=１Aあたりの時間シフト量)を行列に配置
      }
    }

  } // cur

  ///////////////////////////
  // RK4 is finished here. //
  ///////////////////////////

  // save "MATRIX" in "matrix.dat"   //
  FILE *fp3;
  fp3 = fopen("matrix.dat", "w");
  for (i = 0; i < 11; i++)
  {
    for (j = 0; j < 11; j++)
    {
      fprintf(fp3, "  %3.5f", MAT[i][j]);
    }
    fprintf(fp3, "\n");
  }
  fclose(fp3);
  /////////////////////////////////////

  inv_mat(); // calculate the invert matrix

  // read invert matrix //
  FILE *fp4;
  fp4 = fopen("inv_matrix.dat", "r");
  for (i = 0; i < 11; i++)
  {
    fscanf(fp4, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
           &INV[i][0], &INV[i][1], &INV[i][2], &INV[i][3], &INV[i][4], &INV[i][5], &INV[i][6], &INV[i][7], &INV[i][8], &INV[i][9], &INV[i][10]); // read the invert matrix
  }
  fclose(fp4);
  ////////////////////////

  // Invert matrix * dT //
  for (i = 0; i < 11; i++)
  {
    for (j = 0; j < 11; j++)
    {
      ccur[i] += INV[i][j] * dT0p[j];
    }
  }
  ////////////////////////

  // output for next current setting //
  printf("correction currents are below \n");
  printf("mcur = %4.3f ; // main coil current \n", mcur + ccur[0] * 0.1);
  for (i = 1; i < 10; i++)
  {
    printf("tcur[%d] = %4.3f ; // trim%d coil current \n", i - 1, tcur[i - 1] + (ccur[i] * 10), i);
  }
  printf("tcur[9] = %4.3f ; // trim10 coil current \n", tcur[i - 1] + (ccur[i] * 10) - 10.0);
  /////////////////////////////////////
}
