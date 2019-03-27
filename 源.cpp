#include <glut/glut.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <vld.h>
#include <vector>

#define COL 600
#define ROW 600
#define EP 0.0001  

using namespace std;

struct abc
{
	double x;
	double y;
};

struct abc3
{
	double x;
	double y;
	double z;
};

double Q = 30.0;//旋转度数
int arr[COL][ROW];//窗口大小
abc3 campos = {0,0,2};
abc3 camlook = {0,0,-1};

void Matrix_X(double *m1, double *m2, double *r, int x, int y, int z);

double tmp[4] = {0};


double pointa[] = { 0,0,0,1 };
double pointb[] = { 1,0,0,1 };
double pointc[] = { 1,0,1,1 };
double pointd[] = { 0,0,1,1 };
double pointe[] = { 0,1,1,1 };
double pointf[] = { 1,1,1,1 };
double pointg[] = { 1,1,0,1 };
double pointh[] = { 0,1,0,1 };

abc3 resa;
abc3 resb;
abc3 resc;
abc3 resd;
abc3 rese;
abc3 resf;
abc3 resg;
abc3 resh;

//单位矩阵
double Modelmatrix[] = {	1,0,0,0,
							0,1,0,0,
							0,0,1,0,
							0,0,0,1 };

//绕x轴旋转
double Modelmatrix1[] = {	1,0,0,0,
							0,cos(Q),-sin(Q),0,
							0,sin(Q),cosf(Q),0,
							0,0,0,1 };//x轴旋转

//绕y轴旋转
double Modelmatrix2[] = {
							cos(Q),0,sin(Q),0,
							0,1,0,0,
							-sin(Q),0,cosf(Q),0,
							0,0,0,1 };//y轴旋转

 
double Viewmatrix[] = {		1,0,0,0,
							0,1,0,0,
							0,0,1,0,
							0,0,0,1 };

double Viewmatrix1[] = {	1,0,0,0,
							0,cos(Q),sin(Q),0,
							0,-sin(Q),cosf(Q),0,
							0,0,0,1 };//x轴旋转

double Viewmatrix2[] = {
							cos(Q),0,-sin(Q),0,
							0,1,0,0,
							sin(Q),0,cosf(Q),0,
							0,0,0,1 };//y轴旋转

double *Getmatri4()
{
	double *p = new double[15];
	for (int i = 0; i < 15; i++)
	{
		p[i] = i;
	}
	return p;
}

void Delete(double *arr)
{
	if (arr != NULL)
	{
		delete[]arr;
	}
}

//视图矩阵
double* Viewmatri0(abc3 campos, abc3 camlook)
{
	double *m = new double[15]();
	//double campos[] = { 0,0,2 };//p
	//double camlook[] = { 0,0,-1 };//c
	double u[] = { 0,1,0 };//u
	//d = c - p;
	//r = d x u;
	//u = r x d;
	double d[] = { campos.x -camlook.x ,campos.y - camlook.y ,campos.z - camlook.z };
	double r[] = { d[1]*u[2]-u[1]*d[2],-(d[0]*u[2]-u[1]*d[2]),d[0]*u[1]-u[0]*d[1] };
	double u1[] = { r[1]*d[2]-d[1]*r[2],-(r[0]*d[2]-d[1]*r[2]),r[0]*d[1]-d[0]*r[1] };
	double T[] = {	1,0,0,0,
					0,1,0,0,
					0,0,1,0,
					-campos.x,-campos.y ,- campos.z,1 };
	double B[] = {	r[0],d[0],u1[0],0,
					r[1],d[1],u1[1],0, 
					r[2],d[2],u1[2],0, 
					0,0,0,1};
	Matrix_X(T, B, m, 4, 4, 4);
	return m;
}

//矩阵相乘
double* Matrix_X0(double *m1, double *m2, int x, int y, int z)
{
	double *r = new double[x*z]();
	r[0] = m1[0] * m2[0] + m1[1] * m2[4] + m1[2] * m2[8] + m1[3] * m2[12];
	r[1] = m1[0] * m2[1] + m1[1] * m2[5] + m1[2] * m2[9] + m1[3] * m2[13];
	r[2] = m1[0] * m2[2] + m1[1] * m2[6] + m1[2] * m2[10] + m1[3] * m2[14];
	r[3] = m1[0] * m2[3] + m1[1] * m2[7] + m1[2] * m2[11] + m1[3] * m2[15];
	return r;
}

//矩阵相乘
void Matrix_X(double *m1, double *m2, double *r,int x, int y, int z)
{
	int i, j, k;
	for (i = 0; i < x; ++i)
	{
		for (j = 0; j < z; ++j)
		{
			
			double *p1 = m1 + i * y;//p1表示第一个矩阵的行
			
			double *p2 = m2 + j;//p2表示第二个矩阵的列
			for (k = 0; k < y; ++k)
			{
				*r += *p1 * *p2;
				if (k < y - 1)
				{
					
					++p1;
					
					p2 += z;
				}
			}
			++r;
		}
	}
}

//矩阵拷贝
void Copymatrix(double *arr, double *brr,int x,int y)
{
	if (arr == NULL || brr == NULL);
	else
	{
		for (int i = 0; i < x; i++)
		{
			for (int j = 0; j < y; j++)
			{
				arr[i+j] = brr[i+j];
			}
		}
	}
}

//投影矩阵
void touyinmatrix(double *m, double left, double right, double bottom, double top,	double near, double far) 
{
	double r_width = 1.0f / (right - left);
	double r_height = 1.0f / (top - bottom);
	double r_depth = 1.0f / (far - near);
	double x = 2.0f * (r_width);
	double y = 2.0f * (r_height);
	double z = 2.0f * (r_depth);
	double A = (right + left) * r_width;
	double B = (top + bottom) * r_height;
	double C = (far + near) * r_depth;
	m[0] = x;
	m[3] = -A;
	m[5] = y;
	m[7] = -B;
	m[10] = -z;
	m[11] = -C;
	m[15] = 1.0f;
}
//投影矩阵
double* projectmatrix(float l, float r, float b, float t, float n, float f)
{
	double* matrix=new double[15]();
	matrix[0] = 2 / (r - l);
	matrix[1] = 0;
	matrix[2] = 0;
	matrix[3] = -(r + l) / (r-l);
	matrix[4] = 0;
	matrix[5] = 2 / (t - b);
	matrix[6] = 0;
	matrix[7] = -(t + b) / (t - b);
	matrix[8] = 0;
	matrix[9] = 0;
	matrix[10] = -2 / (f - n);
	matrix[11] = -(f + n) / (f - n);
	matrix[12] = 0;
	matrix[13] = 0;
	matrix[14] =0;
	matrix[15] = 1;

	/*matrix[0] = 2 * n / (r - l);
	matrix[2] = -(r+l)/(r-l);
	matrix[5] = 2 * n / (t - b);
	matrix[6] = (t+b)/(t-b);
	matrix[10] = -(f + n) / (f - n);
	matrix[11] = -(2 * f * n) / (f - n);
	matrix[14] = -1;
	matrix[15] = 0;*/
	
	return matrix;
}

void Getxy(double *point)
{
	point[0] = point[0] / point[3] * COL;
	point[1] = point[1] / point[3] * ROW;
	point[2] = point[2];
}

double GetS(const abc3 p0, const abc3 p1, const abc3 p2)
{
	abc AB, BC;
	AB.x = p1.x - p0.x;
	AB.y = p1.y - p0.y;
	BC.x = p2.x - p1.x;
	BC.y = p2.y - p1.y;
	return fabs((AB.x * BC.y - AB.y * BC.x)) / 2.0f;
}

//判断D点在三角形中否
bool IsIn(const abc3 A, const abc3 B, const abc3 C, const abc3 D)
{
	double SABC, SADB, SBDC, SADC;
	SABC = GetS(A, B, C);
	SADB = GetS(A, D, B);
	SBDC = GetS(B, D, C);
	SADC = GetS(A, D, C);

	double SumSuqar = SADB + SBDC + SADC;

	if ((-EP < (SABC - SumSuqar)) && ((SABC - SumSuqar) < EP))
	{
		return true;
	}
	else
	{
		return false;
	}
}

//染色矩形
void Myprint(int a, int b, int c, int d, int color)
{
	for (int i = a; i <= c; i++)
	{
		for (int j = b; j <= d; j++)
		{
			arr[i][j] = color;
		}
	}
	glDrawPixels(COL, ROW, GL_RGBA, GL_UNSIGNED_BYTE, &arr);
	glutSwapBuffers();
}

//染色三角形ABC
void Ttip(const abc3 A, const abc3 B, const abc3 C,int color)
{
	for (int i = 0; i <= COL; i++)
	{
		for (int j = 0; j <= ROW; j++)
		{
			abc3 D;
			D.x = i;
			D.y = j;
			if (IsIn(A, B, C, D))
			{
				int x = D.x;
				int y = D.y;
				arr[x][y] = color;
			}
		}
	}
	glDrawPixels(COL, ROW, GL_RGBA, GL_UNSIGNED_BYTE, &arr);
	glutSwapBuffers();
}

//给出任意两点划条线
void Linesp(abc3 a, abc3 b,int color)
{ 
	 int len = 0; int tmp = 0; 
	int minx, miny;
	int lenh, lenw;
	
	if (a.y < b.y)
	{
		miny = a.y;
		lenw = b.y - a.y;
		tmp = a.y;
	}
	else
	{
		miny = b.y;
		lenw = a.y - b.y;
		tmp = b.y;
	}

	if (a.x < b.x)
	{
		minx = a.x;
		lenh = b.x - a.x;

	}
	else
	{
		minx = b.x;
		lenh = a.x - b.x;
	}
	if (a.y == b.y)
	{
		int f = a.y;
		for (int p = minx; p <= minx + lenh; p++)
		{
			arr[p][f] = color;
		}
	}
	else
	{
		double k = (b.x - a.x) / (b.y - a.y);
		double t = a.x - (k * a.y);
		for (; miny <= lenw + tmp; miny++)
		{
			int j = miny * k + t;
			arr[j][miny] = color;
		}
	}
	glDrawPixels(COL, ROW, GL_RGBA, GL_UNSIGNED_BYTE, &arr);
	glutSwapBuffers();
}

void Print()
{
	double pointa[] = { 0,0,0,1 };
	double pointb[] = { 1,0,0,1 };
	double pointc[] = { 1,0,1,1 };
	double pointd[] = { 0,0,1,1 };
	double pointe[] = { 0,1,1,1 };
	double pointf[] = { 1,1,1,1 };
	double pointg[] = { 1,1,0,1 };
	double pointh[] = { 0,1,0,1 };

	double *tmp=Matrix_X0(pointd, Modelmatrix, 1, 4, 4);
	Copymatrix(pointd, tmp, 1, 4);
	double *viewarr = Viewmatri0(campos, camlook);//视图矩阵
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << viewarr[i*4 + j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	tmp = Matrix_X0(pointd, viewarr, 1, 4, 4);
	Copymatrix(pointd, tmp, 1, 4);
	double *touyinarr = projectmatrix( 0, 60, -0, 60, 1, -10);//投影矩阵
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << touyinarr[i*4 + j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	tmp = Matrix_X0(pointd, touyinarr, 1, 4, 4);
	Copymatrix(pointd, tmp, 1, 4);

	tmp = Matrix_X0(pointa, Modelmatrix, 1, 4, 4);
	Copymatrix(pointa, tmp, 1, 4);
	tmp = Matrix_X0(pointa, viewarr, 1, 4, 4);
	Copymatrix(pointa, tmp, 1, 4);
	tmp = Matrix_X0(pointf, touyinarr, 1, 4, 4);
	Copymatrix(pointa, tmp, 1, 4);

	tmp = Matrix_X0(pointf, Modelmatrix, 1, 4, 4);
	Copymatrix(pointf, tmp, 1, 4);	
	tmp = Matrix_X0(pointf, viewarr, 1, 4, 4);
	Copymatrix(pointf, tmp, 1, 4);
	tmp = Matrix_X0(pointf, touyinarr, 1, 4, 4);
	Copymatrix(pointf, tmp, 1, 4);

	tmp = Matrix_X0(pointc, Modelmatrix, 1, 4, 4);
	Copymatrix(pointc, tmp, 1, 4);
	tmp = Matrix_X0(pointc, viewarr, 1, 4, 4);
	Copymatrix(pointc, tmp, 1, 4);
	tmp = Matrix_X0(pointc, touyinarr, 1, 4, 4);
	Copymatrix(pointc, tmp, 1, 4);

	tmp = Matrix_X0(pointe, Modelmatrix, 1, 4, 4);
	Copymatrix(pointe, tmp, 1, 4);
	tmp = Matrix_X0(pointe, viewarr, 1, 4, 4);
	Copymatrix(pointe, tmp, 1, 4);
	tmp = Matrix_X0(pointe, touyinarr, 1, 4, 4);
	Copymatrix(pointe, tmp, 1, 4);
	/***************************************************************/
	
	resa.x = pointa[0] / pointa[3] * COL+300;
	resa.y = pointa[1] / pointa[3] * COL+300;
	resa.z = pointa[2];
	//Getxy(pointb);
	resb.x = pointb[0] / pointb[3] * COL;
	resb.y = pointb[1] / pointb[3] * COL;
	resb.z = pointb[2];
	//Getxy(pointc);
	resc.x = pointc[0] / pointc[3] * COL+300;
	resc.y = pointc[1] / pointc[3] * COL+300;
	resc.z = pointc[2];
	//Getxy(pointd);
	resd.x = pointd[0] / pointd[3] * COL+300;
	resd.y = pointd[1] / pointd[3] * ROW+300;
	resd.z = pointd[2];
	//Getxy(pointe);
	rese.x = pointe[0] / pointe[3] * COL+300;
	rese.y = pointe[1] / pointe[3] * COL+300;
	rese.z = pointe[2];
	//Getxy(pointf);
	resf.x = pointf[0] / pointf[3] * COL+300;
	resf.y = pointf[1] / pointf[3] * ROW+300;
	resf.z = pointf[2];
	//Getxy(pointg);
	resg.x = pointg[0] / pointg[3] * COL;
	resg.y = pointg[1] / pointg[3] * COL;
	resg.z = pointg[2];
	//Getxy(pointh);
	resh.x = pointh[0] / pointh[3] * COL;
	resh.y = pointh[1] / pointh[3] * COL;
	resh.z = pointh[2];
	abc3 a = { 100,100,20 };
	abc3 b = { 100,150,20 };
	abc3 c = { 150,150,20 };
	abc3 d = { 150,100,20 };
	Linesp(resd, rese, 0xffffff00);
	Linesp(resd, resc, 0xffffff00);
	Linesp(resc, rese, 0xffffff00);
	//Myprint(a.x, a.y, c.x, c.y, 0xffffff00);
	//Ttip(a, b, c, 0xff00ffff);
	/*Linesp(a,b, 0xffffff00);
	Linesp(b, c, 0xffffff00);
	Linesp(c, d, 0xffffff00);
	Linesp(d, a, 0xffffff00);
	Linesp(d, b, 0xffffff00);
	Linesp(c, a, 0xffffff00);*/
	//Myprint(resd.x ,resd.y ,resf.x ,resf.y,0xffffff00);
	/*Ttip(resd, resc, rese, 0xff0000ff);
	Ttip(resf, resc, rese, 0xff0000ff);*/
	//Ttip(ra, rb, rc, 0xff00ffff);
	/*Ttip(rese, resd, rc, 0xff00ffff);*/
	
}

void ChangeSize(int w, int h)
{
	if (h == 0)	h = 1;

	glViewport(0, 0, w, h);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	if (w <= h)
		glOrtho(-100.0f, 100.0f, -100.0f*h / w, 100.0f*h / w, -100.0f, 100.0f);
	else
		glOrtho(-100.0f*w / h, 100.0f*w / h, -100.0f, 100.0f, -100.0f, 100.0f);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void SpecialKeys(int key, int x, int y)
{
	if (key == GLUT_KEY_UP)		;
	if (key == GLUT_KEY_DOWN)	;
	if (key == GLUT_KEY_LEFT)	;
	if (key == GLUT_KEY_RIGHT)	;

	glutPostRedisplay();
}

int main()
{
	cout << "okokokoook" << endl;
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowPosition(800, 300);
	glutInitWindowSize(600, 600);
	glutCreateWindow("刘佳伟");
	glutReshapeFunc(ChangeSize);
	glutSpecialFunc(SpecialKeys);
	glutDisplayFunc(Print);
	glutMainLoop();
	return 0;
}