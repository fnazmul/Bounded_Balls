
	/****************************************/
	/**** Project Name:	Bounded Balls	 ****/
	/**** Submitted by:	Fariha Nazmul	 ****/
	/****				Roll # 021		 ****/
	/****				4th Year		 ****/
	/****				Session: 2005-06 ****/
	/****************************************/



/*****include files*****/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<GL/glut.h>
/***********************/

//define the axis value
#define X 0		
#define Y 1
#define Z 2
//resolution of the screen
#define WIDTH 600
#define HEIGHT 480 
#define W WIDTH/2
#define H HEIGHT/2

#define no_balls 5

/****global variables *****/

//information regarding a ball
struct Ball 
{
	float v[3]; //Velocity
	int pos[3]; //Position
	float r; //Radius
};

struct
{
	int flag;
	double x,y,z;
}temp[HEIGHT+1][WIDTH+1];

//array to store the triangle points and co-ordinate of vertices
int planes[6][4];
double rotated[3];
double ver[12][3];
int vertex[12][2];
float gravity = 8.0f;
float cube_len = 240.0;	
int wh_axis, theta = 0;
int no_ver, no_planes;
double tmp_ver[12][3];
double N[3], C[3], L[3];
struct Ball Balls[no_balls];
double zbuf[HEIGHT+1][WIDTH+1];					//z-buffer
double frame_buffer[HEIGHT+1][WIDTH+1] = {0.0};	//frame buffer
double color_buffer[HEIGHT+1][WIDTH+1];
int prevx = 0, prevy = 0, mouse_button = 0;	//for mouse control

//The number of milliseconds to which the timer is set
int TIMER_MS = 25;
float timeUntilUpdate =0.01;
//The amount of time between each time that we handle collisions 
float TIME_BETWEEN_UPDATES = 0.01f;


enum Wall {WALL_LEFT, WALL_RIGHT, WALL_FAR, WALL_NEAR, WALL_TOP, WALL_BOTTOM};

int wall_dir[6][3] = {{-1, 0, 0},{1, 0, 0},{0, 0, -1},{0, 0, 1},{0, 1, 0},{0, -1, 0}};

double colors[6][3] = {{1.0,0.0,0.0},{1.0,1.0,0.0},{0.0,1.0,0.0},
						{0.0,1.0,1.0},{0.0,0.0,1.0},{1.0,0.0,1.0}};

int px = 0,py = 0, zp = -300;			//PP
int cx = 3, cy = 100, cz = -700;		//COP
int lx = -500, ly = 1000, lz = -1000;	//light source

/**********************************************/




/***************************************************/
/* this function initializes the ball informations */
/***************************************************/

void init_balls(void)
{
	int i;
	for(i=0; i<no_balls;i++) 
	{
		//initialize the ball's position or centre
		Balls[i].pos[0] = i * 25 - 4 + rand()%15;
		Balls[i].pos[1] = i * 5 - 5+ rand()%15;
		Balls[i].pos[2] = i * 7 - 3+ rand()%15;

		//initialize the ball's velocity
		Balls[i].v[0] = 0.1;
		Balls[i].v[1] = 0.0;
		Balls[i].v[2] = 0.1;

		//ball's radius
		Balls[i].r = 25;
	}
}


/****************************************************/
/* this function moves all of the balls according	*/
/* to their velocity								*/
/****************************************************/

void move_balls(float dt) 
{
	int i;

	for(i=0; i<no_balls; i++) 
	{
		Balls[i].pos[0] += Balls[i].v[0] * dt;
		Balls[i].pos[1] += Balls[i].v[1] * dt;
		Balls[i].pos[2] += Balls[i].v[2] * dt;
	}
}


/****************************************************/
/* this function decreases the y coordinate of the 	*/
/* velocity	of each ball by gravity*dt				*/
/****************************************************/

void apply_gravity(void) 
{
	int i;
	for(i=0; i<no_balls; i++) 
		Balls[i].v[1] -= gravity*TIME_BETWEEN_UPDATES;	
}


/****************************************************/
/* Returns whether two balls are colliding			*/
/****************************************************/

int testBallBallCollision(int b1, int b2) 
{
	float r,x,y,z, distance, dot_prdct;
	float net_velocity[3], displacement[3];

	//Checks whether the balls are close enough
	r = Balls[b1].r + Balls[b2].r;

	x = Balls[b1].pos[0] - Balls[b2].pos[0];
	y = Balls[b1].pos[1] - Balls[b2].pos[1];
	z = Balls[b1].pos[2] - Balls[b2].pos[2];
	distance = x*x + y*y + z*z;

	if (distance < r*r) 
	{
		//Checks whether the balls are moving toward each other
		net_velocity[0] = Balls[b1].v[0] - Balls[b2].v[0];
		net_velocity[1] = Balls[b1].v[1] - Balls[b2].v[1];
		net_velocity[2] = Balls[b1].v[2] - Balls[b2].v[2];

		displacement[0] = Balls[b1].pos[0] - Balls[b2].pos[0];
		displacement[1] = Balls[b1].pos[1] - Balls[b2].pos[1];
		displacement[2] = Balls[b1].pos[2] - Balls[b2].pos[2];

		dot_prdct = net_velocity[0]*displacement[0] +
					net_velocity[1]*displacement[1]	+
					net_velocity[2]*displacement[2];

		if(dot_prdct<0)
			return 1;
		else 
			return 0;
	}
	else
		return 0;
}


/****************************************************/
/* Handles all ball-ball collisions					*/
/****************************************************/

void handleBallBallCollisions(void) 
{	
	int i,j,b1,b2;
	float disp[3], value, dot_prdct;
	

	for(i=0; i<no_balls; i++) 
	{
		for(j=i + 1; j<no_balls; j++) 
		{
			b1 = i;
			b2 = j;
			
			if(testBallBallCollision(b1, b2)) 
			{
				//Makes the balls reflect off of each other
				disp[0] = Balls[b1].pos[0] - Balls[b2].pos[0];
				disp[1] = Balls[b1].pos[1] - Balls[b2].pos[1];
				disp[2] = Balls[b1].pos[2] - Balls[b2].pos[2]; 

				value = sqrt(disp[0]*disp[0] + disp[1]*disp[1] + disp[2]*disp[2]);
				disp[0] /= value;
				disp[1] /= value;
				disp[2] /= value;

				dot_prdct = Balls[b1].v[0]*disp[0] + Balls[b1].v[1]*disp[1] 
							+ Balls[b1].v[2]*disp[2];
					
				Balls[b1].v[0] -= 2*disp[0]*dot_prdct;
				Balls[b1].v[1] -= 2*disp[1]*dot_prdct;
				Balls[b1].v[2] -= 2*disp[2]*dot_prdct;

				dot_prdct = Balls[b2].v[0]*disp[0] + Balls[b2].v[1]*disp[1] 
							+ Balls[b2].v[2]*disp[2];
				Balls[b2].v[0] -= 2*disp[0]*dot_prdct;
				Balls[b2].v[1] -= 2*disp[1]*dot_prdct;
				Balls[b2].v[2] -= 2*disp[2]*dot_prdct;
			}
		}
	}

}


/****************************************************/
/* Returns whether a ball and a wall are colliding	*/
/****************************************************/

int testBallWallCollision(int b, int wall) 
{
	int dir[3];
	float pos_dot_dir, v_dot_dir;

	dir[0] = wall_dir[wall][0];
	dir[1] = wall_dir[wall][1];
	dir[2] = wall_dir[wall][2];

	//whether the ball is far enough in the "dir" direction, and whether
	//it is moving toward the wall
	pos_dot_dir = Balls[b].pos[0]*dir[0] + Balls[b].pos[1]*dir[1] 
					+Balls[b].pos[2]*dir[2];
	v_dot_dir = Balls[b].v[0]*dir[0] + Balls[b].v[1]*dir[1] 
					+Balls[b].v[2]*dir[2];

	if((pos_dot_dir + Balls[b].r > cube_len/2) && (v_dot_dir > 0))
		return 1;
	else
		return 0;
}


/****************************************************/
/* Handles all ball-wall collisions					*/
/****************************************************/

void handleBallWallCollisions(void) 
{
	int w;
	int i,j,b;	
	float dir[3],v_dot_dir, value;

	int walls[]={WALL_LEFT, WALL_RIGHT, WALL_FAR, WALL_NEAR, WALL_TOP, WALL_BOTTOM};

	for(i=0; i<no_balls; i++) 
	{
		for(j=0; j<6; j++) 
		{
			b = i;
			w = walls[j];
			
			if(testBallWallCollision(b, w)) 
			{
				//Makes the ball reflect off of the wall
				dir[0] = (float)wall_dir[w][0];
				dir[1] = (float)wall_dir[w][1];
				dir[2] = (float)wall_dir[w][2];

				value = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
				dir[0] /= value;
				dir[1] /= value;
				dir[2] /= value;
				
				v_dot_dir = Balls[b].v[0]*dir[0] + Balls[b].v[1]*dir[1] 
					+Balls[b].v[2]*dir[2];
				
				Balls[b].v[0] -= 2*dir[0]* v_dot_dir;
				Balls[b].v[1] -= 2*dir[1]* v_dot_dir;
				Balls[b].v[2] -= 2*dir[2]* v_dot_dir;
			}
		}
	}
}


/****************************************************/
/* Applies gravity and handles all collisions.		*/
/****************************************************/

void performUpdate(void) 
{
	apply_gravity();
	move_balls(10);	
	handleBallBallCollisions();
	handleBallWallCollisions();
}


/****************************************************/
/* Advances the state of the balls by t.			*/
/* timeUntilUpdate is the amount of time			*/
/* until the next call to performUpdate.			*/
/****************************************************/

void advance(float t, float timeUntilUpdate) 
{
	while (t > 0) 
	{
		if (timeUntilUpdate <= t) 
		{
			move_balls(timeUntilUpdate);
			performUpdate();
			t -= timeUntilUpdate;
			timeUntilUpdate = TIME_BETWEEN_UPDATES;
		}
		else 
		{
			move_balls(timeUntilUpdate);
			timeUntilUpdate -= t;
			t = 0;
		}
	}
}


/****************************************************/
/* Called every TIMER_MS milliseconds.				*/
/****************************************************/

void update(int value) 
{
	advance((float)TIMER_MS / 1000.0f, timeUntilUpdate);
	theta += (float)TIMER_MS / 100;
	if (theta> 360) 
		theta -= 360;
		
	glutPostRedisplay();
	glutTimerFunc(TIMER_MS, update, 0);
}


/****************************************************/
/* the default reshape function						*/
/****************************************************/

void reshape(int width, int height)
{
	glViewport(0,0,width,height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-W, W-1,-H,H-1,-1000,1000);
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
}


/****************************************************/
/*	scan the input data								*/
/****************************************************/

void read_input()
{
	int i;
	double w;

	freopen("cube.txt","r",stdin);
	w = 0.08;
	

	scanf("%d",&no_ver);	//total no of vertices
	scanf("%d",&no_planes);	//total no of planes
		
	for(i=0;i<no_ver;i++)		//scan all vertex co-ordinate
	{
		scanf("%lf %lf %lf",&ver[i][0], &ver[i][1], &ver[i][2]); 		 
		ver[i][0] /= w;		
		ver[i][1] /= w;
		ver[i][2] /= w;

		tmp_ver[i][0] = ver[i][0] ;		
		tmp_ver[i][1] = ver[i][1] ;
		tmp_ver[i][2] = ver[i][2] ;			
	}
	
	for(i=0;i<no_planes;i++)		//scan all plane points
	{
		scanf("%d",&planes[i][0]);
		scanf("%d",&planes[i][1]);
		scanf("%d",&planes[i][2]);
		scanf("%d",&planes[i][3]);
	}

	init_balls();
}


/**********************************************************/
/*	This function draws a pixel of the line on the screen */
/**********************************************************/

void drawPixel(int x, int y, double nz, int option)	
{

	if(option==0)
	{
		if( y<H && y>=-H && x<W && x>=-W)
		{
			temp[y+H][x+W].flag = 1;
			temp[y+H][x+W].z = nz;
		}
	}
	if(option==3)
	{
		if( y<H && y>=-H && -x<W && -x>=-W)
		{
			temp[y+H][-x+W].flag = 1;
			temp[y+H][-x+W].z = nz;
		}
	}
	if(option==4)
	{
		if( -y<H && -y>=-H && -x<W && -x>=-W)
		{
			temp[-y+H][-x+W].flag = 1;
			temp[-y+H][-x+W].z = nz;
		}
	}
	if(option==7)
	{
		if( -y<H && -y>=-H && x<W && x>=-W)
		{
			temp[-y+H][x+W].flag = 1;
			temp[-y+H][x+W].z = nz;
		}
	}
	if(option==1)
	{
		if( x<H && x>=-H && y<W && y>=-W)
		{
			temp[x+H][y+W].flag = 1;
			temp[x+H][y+W].z = nz;
		}
	}
	if(option==2)
	{
		if( x<H && x>=-H && -y<W && -y>=-W)
		{
			temp[x+H][-y+W].flag = 1;
			temp[x+H][-y+W].z = nz;
		}
	}
	if(option==5)
	{
		if( -x<H && -x>=-H && -y<W && -y>=-W)
		{
			temp[-x+H][-y+W].flag = 1;
			temp[-x+H][-y+W].z = nz;
		}
	}
	if(option==6)
	{
		if( -x<H && -x>=-H && y<W && y>=-W)
		{
			temp[-x+H][y+W].flag = 1;	
			temp[-x+H][y+W].z = nz;
		}
	}	
}


/****************************************************/
/*	draws a line using Midpoint line drawing alg	*/
/****************************************************/

void drawLine(int x0, int y0, double z0, int x1, int y1, double z1, int op)	
{
	int dx = x1-x0;
	int dy = y1-y0;

	int dinitial = 2*dy - dx;
	int de = 2*dy;
	int dne = 2* (dy - dx);

	int x = x0;
	int y = y0;
	double nz;

	drawPixel(x, y, z0, op);
	

	while(x<x1)
	{
		if(dinitial <= 0)// east direction
		{
			dinitial += de;
		}
		else		//else NE direction
		{
			y++;
			dinitial += dne;
		}
		x++;
		
		nz = z0 - ((float)(x1 - x)/(float)(x1 - x0))*(z1 - z0);		
		drawPixel(x, y, nz, op);
	}
}


/****************************************************/
/*	this function defines the group & option of the	*/
/* line to be drawn									*/
/****************************************************/

void drawSlopeIndpLine(int x0, int y0, double z0, int x1, int y1, double z1)	
{														
	int dx = x1 - x0;
	int dy = y1 - y0;

	if(abs(dx) >= abs(dy))						// group 0
	{
		if( (x1>=x0) && (y1>=y0))
			drawLine(x0, y0, z0, x1, y1, z1, 0);		//op = 0
		else if((x1<x0) && (y1>=y0))
			drawLine(-x0, y0, z0, -x1, y1, z1, 3);		//op = 3
		else if((x1<x0) && (y1<y0))
			drawLine(-x0, -y0, z0, -x1, -y1, z1, 4);	//op = 4
		else if((x1>=x0) && (y1<y0))
			drawLine(x0, -y0, z0, x1, -y1, z1, 7);		//op = 7
	}
	else										//group = 1...
	{
		if((x1>=x0) && (y1>=y0))
			drawLine(y0, x0, z0, y1, x1, z1, 1);		//op = 1
		else if((x1<x0) && (y1>=y0))
			drawLine(y0, -x0, z0, y1, -x1, z1, 2);		//op = 2
		else if((x1<x0) && (y1<y0))
			drawLine(-y0, -x0, z0, -y1, -x1, z1, 5);	//op = 5
		else if((x1>=x0) && (y1<y0))
			drawLine(-y0, x0, z0, -y1, x1, z1, 6);		//op = 6
	}
	
}


/****************************************************/
/*	this function projects the 3D data into 2D data	*/
/****************************************************/

void project_data()
{

	int i,xp,yp;

	double x,y,z,x1,y1;
	double Q,dx,dy,dz,denom;
	double qx,qy,qz;

	qx = cx - px; 	//3;
	qy = cy - py;	//100;
	qz = cz - zp;	//-300;
	
	Q = sqrt((qx*qx) + (qy*qy) + (qz*qz) );
	dx = qx/Q;
	dy = qy/Q;
	dz = qz/Q;	

	for(i=0;i<no_ver;i++)
	{
		x = tmp_ver[i][0];
		y = tmp_ver[i][1];
		z = tmp_ver[i][2];

		denom = (zp-z)/(Q*dz) +1;
		x1 = x - ( (z*dx)/dz ) + ( (zp*dx)/dz );
		y1 = y - ( (z*dy)/dz ) + ( (zp*dy)/dz );
		xp = (int) (x1/denom);
		yp = (int) (y1/denom);
		vertex[i][0] = xp;
		vertex[i][1] = yp;		
		
 	}	
}


/********************************************************/
/*	rotates a point to theta degree across any one axis	*/
/********************************************************/

void rotate(double x0,double y0,double z0)
{
	double rad = 3.141592654/180.0;
	double x, y, z;

	if (wh_axis == X) 
	{
		x = x0;
		y = y0 * cos(theta*rad) - z0*sin(theta*rad);
		z = y0 * sin(theta*rad) + z0*cos(theta*rad);
	}
	else if (wh_axis == Y) 
	{
		x =   x0*cos(theta*rad) + z0*sin(theta*rad);
		y =   y0;
		z = -(x0*sin(theta*rad)) + z0*cos(theta*rad);
	}
	else if (wh_axis == Z)
	{
		x = x0*cos(theta*rad) - y0*sin(theta*rad);
		y = x0*sin(theta*rad) + y0*cos(theta*rad);
		z = z0;
	}
		
	rotated[0] = x;
	rotated[1] = y;
	rotated[2] = z;
}


/********************************************************/
/*	rotates all balls and the cube						*/
/********************************************************/

void rotate_all()
{
	int i;
	for(i=0;i<no_ver;i++)
	{
		rotate(tmp_ver[i][0],tmp_ver[i][1],tmp_ver[i][2]);
		tmp_ver[i][0] = rotated[0];
		tmp_ver[i][1] = rotated[1];
		tmp_ver[i][2] = rotated[2];
	}

	for(i =0;i<no_balls;i++)
	{	
		rotate(Balls[i].pos[0],Balls[i].pos[1],Balls[i].pos[2]);
		Balls[i].pos[0] = rotated[0];
		Balls[i].pos[1] = rotated[1];
		Balls[i].pos[2] = rotated[2];
	}
}


/********************************************************/
/*	this function calculates the vecor L, N	for each	*/
/*   cube surface										*/
/********************************************************/

void find_vectors(int plane_num)
{	
	double iso[3], Vo[3], V1[3];
	double value;
	
	Vo[0] = tmp_ver[planes[plane_num][1]][0] - tmp_ver[planes[plane_num][0]][0];
	Vo[1] = tmp_ver[planes[plane_num][1]][1] - tmp_ver[planes[plane_num][0]][1];
	Vo[2] = tmp_ver[planes[plane_num][1]][2] - tmp_ver[planes[plane_num][0]][2];

	value = sqrt(Vo[0]*Vo[0] + Vo[1]*Vo[1] + Vo[2]*Vo[2]);
	Vo[0] /= value;		Vo[1] /= value;		Vo[2] /= value;

	V1[0] = tmp_ver[planes[plane_num][3]][0] - tmp_ver[planes[plane_num][0]][0];
	V1[1] = tmp_ver[planes[plane_num][3]][1] - tmp_ver[planes[plane_num][0]][1];
	V1[2] = tmp_ver[planes[plane_num][3]][2] - tmp_ver[planes[plane_num][0]][2];

	value = sqrt(V1[0]*V1[0] + V1[1]*V1[1] + V1[2]*V1[2]);
	V1[0] /= value;		V1[1] /= value;		V1[2] /= value;

	N[0] = Vo[1]*V1[2] - Vo[2]*V1[1];
	N[1] = Vo[2]*V1[0] - Vo[0]*V1[2];
	N[2] = Vo[0]*V1[1] - Vo[1]*V1[0];

	iso[0] = (tmp_ver[planes[plane_num][0]][0] + tmp_ver[planes[plane_num][1]][0] + tmp_ver[planes[plane_num][2]][0]+ tmp_ver[planes[plane_num][3]][0]) /4;
	iso[1] = (tmp_ver[planes[plane_num][0]][1] + tmp_ver[planes[plane_num][1]][1] + tmp_ver[planes[plane_num][2]][1]+ tmp_ver[planes[plane_num][3]][1]) /4;
	iso[2] = (tmp_ver[planes[plane_num][0]][2] + tmp_ver[planes[plane_num][1]][2] + tmp_ver[planes[plane_num][2]][2]+ tmp_ver[planes[plane_num][3]][2]) /4;

	L[0] = lx - iso[0];		L[1] = ly - iso[1];		L[2] = lz - iso[2];	
	value = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]);
	L[0] /= value;		L[1] /= value;		L[2] /= value;
}


/********************************************************/
/*	this function calculates the vecor L, N	of a sphere */
/********************************************************/

void find_vectors_sphere(int h, int w, int a, int b, int c)
{	
	double x, y, z;
	double value;

	x = temp[h][w].x ;
	y = temp[h][w].y ;
	z = temp[h][w].z ;
	
	//surface normal
	N[0] = (x - a);
	N[1] = (y - b);
	N[2] = (z - c);

	value = sqrt(N[0]*N[0] + N[1]*N[1] + N[2]*N[2]);
	N[0] /= value;		N[1] /= value;		N[2] /= value;

	L[0] = lx - x;		L[1] = ly - y;		L[2] = lz - z;	

	value = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]);
	L[0] /= value;		L[1] /= value;		L[2] /= value;
}


/********************************************************/
/*	this function calculates the dot product of the		*/
/*  vectors N and L and determines the intensity		*/
/********************************************************/

double intensity(void)
{
	double i,inten;
	
	i = N[0]*L[0] + N[1]*L[1] + N[2]*L[2];

	if(i<0.0)
		i = 0.0;
	inten = i*0.9 +0.1 ;
	
	return inten;
}


/********************************************************/
/*	finds the minimum value among the three parameters	*/
/********************************************************/

int find_min(int a,int b, int c,int d)
{
	int min;
	min = (a < b) ? a : b;
	min = (c < min) ? c : min;
	min = (d < min) ? d : min;

	return min;
}


/********************************************************/
/*	finds the maximum value among the three parameters	*/
/********************************************************/

int find_max(int a, int b, int c,int d)
{
	int max;
	max = (a > b) ? a : b;
	max = (c > max) ? c : max;
	max = (d > max) ? d : max;

	return max;
}


/********************************************************/
/*	this function initializes all the buffers				*/
/********************************************************/

void init_buffer(void)
{		
	int k,j;
	
	for(k=HEIGHT-1;k>=0;k--)
	{
		   for(j=0;j<WIDTH;j++)
		   {
				temp[k][j].flag = 0;
				temp[k][j].x = temp[k][j].y = 0.0;
				temp[k][j].z = 32500.0;
				frame_buffer[k][j] = 0.0;
				color_buffer[k][j] = 0.0;
				zbuf[k][j] = 32500.0;
		   }
	}
}


/********************************************************/
/*	this function calculates the visible points of a	*/
/*	sphere using ray tracing algorithm					*/
/********************************************************/

void find_intersect(int a, int b, int c, int r)
{
	int w, h, x0, y0, z0, x1, y1, z1, dx, dy, dz;
	double t, t1, t2, A, B, C, Bsqr, AC4, A2;

	x0 = cx;	y0 = cy;	z0 = cz;
	z1 = zp;

	for( h=-H; h<H; h++)
	{
		for( w=-W; w<W; w++)
		{
			x1 = w;		y1 = h;
			dx = x1 - x0;	dy = y1 - y0;	dz = z1 - z0;
			A = (double) ((dx*dx) + (dy*dy) + (dz*dz));
			B = (double) (2*(dx*(x0-a) + dy*(y0-b) + dz*(z0-c)));
			C = (double) ((x0-a)*(x0-a) + (y0-b)*(y0-b) + (z0-c)*(z0-c) - r*r );

			Bsqr = B*B;
			AC4 = 4*A*C;
			A2 = 2*A;

			if(Bsqr == AC4)
				t = -B / A2;
			else if(Bsqr >= AC4)
			{
				t1 = ( -B + sqrt(Bsqr - AC4) )/A2;
				t2 = ( -B - sqrt(Bsqr - AC4) )/A2;
				t = (t1 < t2) ? t1 : t2;
			}
			else
				continue;
			temp[h+H][w+W].x = (double) (x0 + t*dx) ;
			temp[h+H][w+W].y = (double) (y0 + t*dy) ;
			temp[h+H][w+W].z = (double) (z0 + t*dz) ;
			temp[h+H][w+W].flag = 1;
		}
	}
}


/********************************************************/
/*	this function updates frame buffer for ray tracing	*/
/********************************************************/

void fill_color(int a, int b, int c,int col)
{
	int h,w;
	double color,z;

	for(h = 0; h<HEIGHT; h++)
	{			
		for(w=0; w<WIDTH; w++)
		{
			if(temp[h][w].flag)
			{	
				z = temp[h][w].z;

				if (fabs(z - (double)cz) < fabs(zbuf[h][w] - (double)cz)) 
				{
					zbuf[h][w] = z;
					find_vectors_sphere(h,w,a,b,c);	  
					color = intensity();
					frame_buffer[h][w] = color;	
					color_buffer[h][w] = col;
				}							
			}			
		
		   temp[h][w].flag = 0;
		}		
	}	
}


/********************************************************/
/*	this function draws the cube on the temp buffer		*/
/********************************************************/

void drawPoly_zbuffer(void)
{
	int x0, y0, x1, y1, i, j, k, p0, p1;
	double z0, z1, color, z;
	int minx, maxx, miny, maxy;

	for(i = 0;i<no_planes;i++)
	{			
	   find_vectors(i);	  
	   color = intensity();
	   
	   for(j=0;j<4;j++)
	   {			
			p0 =planes[i][j];			p1 =planes[i][(j+1)%4]; 	
		   
		    x0 = vertex[p0][0];			x1 = vertex[p1][0];
			y0 = vertex[p0][1];			y1 = vertex[p1][1];	
			z0 = tmp_ver[p0][2];		z1 = tmp_ver[p1][2];

			drawSlopeIndpLine(x0, y0, z0, x1, y1, z1);
			drawSlopeIndpLine(x0+1, y0+1, z0+1, x1-1, y1-1, z1-1);
			drawSlopeIndpLine(x0-1, y0-1, z0-1, x1+1, y1+1, z1+1);							
	   }// end of each line

			   
	    minx = find_min(vertex[planes[i][0]][0], vertex[planes[i][1]][0], vertex[planes[i][2]][0],vertex[planes[i][3]][0]);
		maxx = find_max(vertex[planes[i][0]][0], vertex[planes[i][1]][0], vertex[planes[i][2]][0],vertex[planes[i][3]][0]);
	    miny = find_min(vertex[planes[i][0]][1], vertex[planes[i][1]][1], vertex[planes[i][2]][1],vertex[planes[i][3]][1]);
		maxy = find_max(vertex[planes[i][0]][1], vertex[planes[i][1]][1], vertex[planes[i][2]][1],vertex[planes[i][3]][1]);

		for(k=maxy;k>=miny;k--)
		{
			 for(j=minx; j<=maxx;j++)
			 {
			   if(temp[k+H][j+W].flag)
			   {		
					z = temp[k+H][j+W].z;
			   
					if (fabs(z - (double)cz) < fabs(zbuf[k+H][j+W] - (double)cz)) 
					{
						zbuf[k+H][j+W] = z;
						frame_buffer[k+H][j+W] = color;	
						color_buffer[k+H][j+W] = i;
					}				   
			   }
			   temp[k+H][j+W].flag = 0;
			 }		   
		} 

	}//end of each planes	
}


/********************************************************/
/*	this function displays the frame buffer				*/
/********************************************************/

void disp_zbuf(void)
{
	int k, j, col;
	double color;
	
	
	for(k=HEIGHT-1;k>=0;k--)
	{
		   for(j=0;j<WIDTH;j++)
		   {
			   if(zbuf[k][j]<32000)
			   {
					color = frame_buffer[k][j];
					col = color_buffer[k][j];
				    glColor3f(colors[col][0]*color, colors[col][1]*color, colors[col][2]*color);
					
			   }
			   else
					glColor3f(0.8f,0.8f,0.8f);					
			   glVertex2i(j-W,k-H);	
		   }
	}
}


/********************************************************/
/*	this is the call event function for display			*/
/********************************************************/

void display(void)
{
	int i;
	glClearColor(1.0f,1.0f,1.0f,1.0f);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	
	glPointSize(1.0);

	glBegin(GL_POINTS);
	
	init_buffer();		//initialize buffer
	rotate_all();		//rotate object
	project_data();		//project object
	drawPoly_zbuffer();	//draws cube
	performUpdate();	//moves balls and handles collisions

	//draws the balls
	for(i=0;i<no_balls;i++)
	{
		find_intersect(Balls[i].pos[0],Balls[i].pos[1],Balls[i].pos[2],Balls[i].r);
		fill_color(Balls[i].pos[0],Balls[i].pos[1],Balls[i].pos[2],i%6);
	}
	//draws the z-buffer	
	disp_zbuf();
		
	glEnd();
	glutSwapBuffers();
}


/********************************************************/
/*	mouse function to detect a click					*/
/********************************************************/

void Mouse( int button , int state, int mx, int my)
{
    switch (button) 
	{
		case GLUT_LEFT_BUTTON:		// when the left button is pressed
			 mouse_button = 0;      // Note which button is pressed, where
			 prevx = mx;
			 prevy = my;
		break;

		case GLUT_MIDDLE_BUTTON:	// when the middle button is pressed
			mouse_button = 1;       
			break;

		case GLUT_RIGHT_BUTTON:		// when the right button is pressed
			mouse_button = 2;       
			break;
    }
}


/********************************************************/
/*	mouse function to detect dragging					*/
/********************************************************/

void Motion( int mx , int my )
{
    int dx, dy;
        
    if( mouse_button == 0 )            //only when left mouse button is pressed
    {
        dx = mx - prevx;		//calc dx and dy
        dy = my - prevy;
        prevx = mx; 			//Save x and y for dx and dy calcs next time
        prevy = my;

		if(dx>0 && dx<WIDTH && abs(dy)< 10)
		{
			theta = -5.0;
			wh_axis = Y;
		}
		else if(dx<0 && abs(dy)< 10)
		{
			theta =5.0;
			wh_axis = Y;
		}
		else if(dy>0 && dy<HEIGHT && abs(dx)< 10)
		{
			theta =5.0;
			wh_axis = X;
		}
		else if(dy<0 && abs(dx)< 10)
		{
			theta = -5.0;
			wh_axis = X;
		}
        glutPostRedisplay();
    }  
}


/********************************************************/
/*	The function is called whenever a "normal"			*/
/*	 key is pressed.									*/
/********************************************************/

void keyControl(GLubyte key, GLint x, GLint y) 
{
    switch ( key )    
	{ 
    
		case 'x':
			theta = 5;
			wh_axis = X;
		    glutPostRedisplay();
		   break;

		case 'y':
			theta = 5;
			wh_axis = Y;
			glutPostRedisplay();
			break;

		case 'z':
			theta = 5;
			wh_axis = Z;
			glutPostRedisplay();
			break; 

		default:
			break;
    }
}


/********************************************************/
/*	called whenever a "special" key is pressed.			*/
/********************************************************/

void specialKeyControl(int key, int x, int y) 
{
    switch ( key )    
	{     
		case GLUT_KEY_LEFT:
			theta = -5.0;
			wh_axis = Y;
			glutPostRedisplay();
			break;
		case GLUT_KEY_RIGHT:
			theta = 5.0;
			wh_axis = Y;
			glutPostRedisplay();
			break;
		case GLUT_KEY_UP:
			theta = -5.0;
			wh_axis = X;
			glutPostRedisplay();
			break; 
		case GLUT_KEY_DOWN:
	        theta = 5.0;
			wh_axis = X;
			glutPostRedisplay();
			break; 
		default:
			break;
    }
}


/********************************************************/
/*			the main function							*/
/********************************************************/

int main(int argc, char **argv)
{	
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(WIDTH, HEIGHT);
	glutCreateWindow("Cube & Balls");
	glutReshapeFunc(reshape);

	read_input();
	
	glutDisplayFunc(display);
	glutKeyboardFunc(keyControl);
	glutSpecialFunc(specialKeyControl);
	glutMouseFunc(Mouse);
	glutMotionFunc(Motion);

	glutTimerFunc(TIMER_MS, update, 0);
	glutMainLoop();
	
	return 0;
}






