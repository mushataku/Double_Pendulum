#define _USE_MATH_DEFINES
#include <cstdio>
#include <cmath>

/************計算条件************/
//区間 [tmin,tmax] で計算
const double tmin = 0.0;
const double tmax = 20.0;
// それを何区間に分割して計算するか
const int nt = 1001;
const double dt = (tmax-tmin)/double(nt-1);

// 初期値の個数
const int INIT_NUM = 100;
// 初期値をどれだけずらすか
const double DTHETA = 1.0e-10;
/*******************************/

/**************初期条件*************/  
const double theta1_0 = M_PI;
const double theta2_0 = M_PI/2.0;
const double omega1_0 = 1.0;
const double omega2_0 = 0.0;
/*********************************/

/*************定数**************/
const double g = 9.81;
const double m1 = 1.0;
const double m2 = 1.0;
const double l1 = 1.0;
const double l2 = 1.0;
/******************************/

/****************************************関数**************************************************/
double f1(double theta1, double theta2, double omega1, double omega2);
double f2(double theta1, double theta2, double omega1, double omega2);
void euler(double &theta1, double &theta2, double &omega1, double &omega2, double dt);
void runge_kutta(double &theta1, double &theta2, double &omega1, double &omega2, double dt);
// θ1とθ2からxy座標上の位置情報を取得
void get_xy(double theta1, double theta2, double &x1, double &y1, double &x2, double &y2);
// 微小角近似の下での解析解を出力
void output_analytic(double t, FILE *analytic_fp);
/**********************************************************************************************/

int main(){
  /****アニメ用の主力のための変数*****/
  double output_time[501];
  int ti = 0;
  /*********************************/
  char filename[30];
  FILE *fp, *num_fp;
  num_fp = fopen("data_many/init_num.txt", "w");
  fprintf(num_fp, "%d\n", INIT_NUM);
  fclose(num_fp);
  double t,theta1,theta2,omega1,omega2,x1,y1,x2,y2;

  //0.00, 0.04, 0.08, 0.12 ...[s] の時にアニメ用の出力をさせるようset
  for(int i = 0; i < 401; i++) {
    output_time[i] = 0.04*double(i);
  }

  // 初期条件をちょっとずつずらして出力
  for(int num = 0; num < INIT_NUM; num++) {
    //printf("num:%d\n", num);
    /*****出力するファイルをオープン************************/
    sprintf(filename, "data_many/many_chaos%d.csv", num);
    fp = fopen(filename, "w");
    /****************************************************/

    //printf("t theta1 theta2 x1 y1 x2 y2\n");
    fprintf(fp, "t,theta1,theta2,x1,y1,x2,y2\n");

    t = tmin;
    ti = 0;
    theta1 = theta1_0;
    omega1 = omega1_0;
    theta2 = theta2_0 + DTHETA*num;
    omega2 = omega2_0;

    //計算開始
    for (int i = 0; i < nt; i++)
    {
      // θ をxyに変換
      get_xy(theta1, theta2, x1, y1, x2, y2);
      //printf("%e %e %e %e %e %e %e\n", t, theta1, theta2, x1, y1, x2, y2);

      //アニメ用の出力
      if (abs(t - output_time[ti]) < 1e-9)
      {
        fprintf(fp, "%e,%e,%e,%e,%e,%e,%e\n", t, theta1, theta2, x1, y1, x2, y2);
        ti++;
      }

      runge_kutta(theta1, theta2, omega1, omega2, dt);
      t += dt;
    }
    fclose(fp);
  }
  return 0;
}

double f1(double theta1, double theta2, double omega1, double omega2){
  double dtheta = theta1 - theta2;
  double M = m2/(m1+m2);
  double l = l2/l1;
  double omega = g/l1;

  double child, mother;

  child = omega*(M*sin(theta2)*cos(dtheta) - sin(theta1))
        - M*(omega1*omega1*cos(dtheta) + l*omega2*omega2)*sin(dtheta);
  mother = 1.0 - M*cos(dtheta)*cos(dtheta);

  return child/mother;
}

double f2(double theta1, double theta2, double omega1, double omega2){
  double dtheta = theta1 - theta2;
  double M = m2/(m1+m2);
  double l = l2/l1;
  double omega = g/l1;

  double child, mother;

  child = omega*(sin(theta1)*cos(dtheta) - sin(theta2))
      + (omega1*omega1 + M*l*omega2*omega2*cos(dtheta))*sin(dtheta);
  mother = l*(1.0 - M*cos(dtheta)*cos(dtheta));

  return child/mother;
}

void runge_kutta(double &theta1, double &theta2, double &omega1, double &omega2, double dt){
  double t1k1,t2k1,w1k1,w2k1,t1k2,t2k2,w1k2,w2k2;
  double t1k3,t2k3,w1k3,w2k3,t1k4,t2k4,w1k4,w2k4;

  t1k1 = omega1;
  t2k1 = omega2;
  w1k1 = f1(theta1, theta2, omega1, omega2);
  w2k1 = f2(theta1, theta2, omega1, omega2);

  t1k2 = omega1+0.5*dt*w1k1;
  t2k2 = omega2+0.5*dt*w2k1;
  w1k2 = f1(theta1+0.5*dt*t1k1, theta2+0.5*dt*t2k1, omega1+0.5*dt*w1k1, omega2+0.5*dt*w2k1);
  w2k2 = f2(theta1+0.5*dt*t1k1, theta2+0.5*dt*t2k1, omega1+0.5*dt*w1k1, omega2+0.5*dt*w2k1);

  t1k3 = omega1+0.5*dt*w1k2;
  t2k3 = omega2+0.5*dt*w2k2;
  w1k3 = f1(theta1+0.5*dt*t1k2, theta2+0.5*dt*t2k2, omega1+0.5*dt*w1k2, omega2+0.5*dt*w2k2);
  w2k3 = f2(theta1+0.5*dt*t1k2, theta2+0.5*dt*t2k2, omega1+0.5*dt*w1k2, omega2+0.5*dt*w2k2);

  t1k4 = omega1+dt*w1k3;
  t2k4 = omega2+dt*w2k3;
  w1k4 = f1(theta1+dt*t1k3, theta2+dt*t2k3, omega1+dt*w1k3, omega2+dt*w2k3);
  w2k4 = f2(theta1+dt*t1k3, theta2+dt*t2k3, omega1+dt*w1k3, omega2+dt*w2k3);

  theta1 += dt*(t1k1 + 2.0*t1k2 + 2.0*t1k3 + t1k4)/6.0;
  theta2 += dt*(t2k1 + 2.0*t2k2 + 2.0*t2k3 + t2k4)/6.0;
  omega1 += dt*(w1k1 + 2.0*w1k2 + 2.0*w1k3 + w1k4)/6.0;
  omega2 += dt*(w2k1 + 2.0*w2k2 + 2.0*w2k3 + w2k4)/6.0;
}

void get_xy(double theta1, double theta2, double &x1, double &y1, double &x2, double &y2){
  x1 = l1*sin(theta1);
  y1 = -l1*cos(theta1);
  x2 = x1 + l2*sin(theta2);
  y2 = y1 - l2*cos(theta2);
}

void output_analytic(double t, FILE *analytic_fp){
  double theta1, theta2, x1, y1, x2, y2;
  double omega_plus  = sqrt(g)*sqrt(2.0 + sqrt(2.0));
  double omega_minus = sqrt(g)*sqrt(2.0 - sqrt(2.0));
  double Theta_plus  = (2.0 - sqrt(2.0))*cos(omega_plus*t);
  double Theta_minus = (2.0 + sqrt(2.0))*cos(omega_minus*t);
  double A = M_PI/400.0;

  theta1 = A * (Theta_plus + Theta_minus);
  theta2 = sqrt(2.0) * A * (-Theta_plus + Theta_minus);

  get_xy(theta1, theta2, x1, y1, x2, y2);
  fprintf(analytic_fp, "%e,%e,%e,%e,%e,%e,%e\n", t, theta1, theta2, x1, y1, x2, y2);
}