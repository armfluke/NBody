#include <stdio.h>
#include <glut.h>

#define ITERATION 1000
#define NUMBER 1000
FILE * pFilex, *pFiley;
int status = 0;
int status1 = 0;
int item, size;
int * mass;
double  posx[ITERATION][NUMBER];
double  posy[ITERATION][NUMBER];
int j = 0;



void setup()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, 800, 400, -10, 0, 1);
	glDisable(GL_DEPTH_TEST);


	//read data from file
	pFilex = fopen("outputx.txt", "r");
	pFiley = fopen("outputy.txt", "r");

	if (pFilex == NULL || pFiley == NULL){
		printf("File not Exist\n");
	}
	else{

		int check = 0;

		for (int i = 0; i < ITERATION; i++){
			for (int j = 0; j<NUMBER; j++){
				fscanf(pFilex, "%lf", &posx[i][j]);
				fscanf(pFiley, "%lf", &posy[i][j]);
			}
		}

	}

}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// draw
	for (int i = 0; i < NUMBER; i++){

		glPointSize(3.0);
		//Draw the particle as a little square.
		glBegin(GL_POINTS);

		glColor3f(0.0, 1.0, 1.0);
		glVertex2f(posx[j / 10][i], posy[j / 10][i]);
		glEnd();
	}
	if (j<ITERATION * 10){
		j++;
	}
	glutPostRedisplay();
	glFlush();
}



int main(int argc, char **argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode((GLUT_SINGLE | GLUT_RGB));

	glutInitWindowSize(800, 600);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("N-Body Simulation");
	if (status == 0){   //read file only one time
		setup();
		status++;
	}
	glutDisplayFunc(display);
	glutMainLoop();
	return 0;
}


