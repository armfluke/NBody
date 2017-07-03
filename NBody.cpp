#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <glut.h>
#define G 6.67428E-11
#define MAX 30000
#define ITERATION 1000

double *create_array(int number) {
	double *data = (double*)malloc(number*sizeof(double));
	return data;
}

int *create_array_int(int number) {
	int *data = (int*)malloc(number*sizeof(int));
	return data;
}

int main(int argc, char *argv[]){
	int i, j, k, l, rank, nproc, start, stop, *mass, number, each_process, position = 0, *displs, *count, additional, remaining = 0;
	double *monitorX, *monitorY, *velocityX, *velocityY, accelatorX, accelatorY, **x, **y, magnitude, forceX = 0, forceY = 0, pack[MAX]
		, timestart, timestop, t = 250;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	FILE *file;
	if (rank == 0){
		file = fopen("inputTest.txt", "r");
		fscanf(file, "%d", &number);
		monitorX = create_array(number);
		monitorY = create_array(number);
		mass = create_array_int(number);
		//read x,y and mass of each dot
		for (i = 0; i < number; i++){
			fscanf(file, "%lf %lf %d", &monitorX[i], &monitorY[i], &mass[i]);
			//printf("Number %d : %lf %lf %d\n", i, monitorX[i], monitorY[i], mass[i]);
		}
		MPI_Pack(&number,1,MPI_INT,pack,MAX,&position,MPI_COMM_WORLD);
		MPI_Pack(monitorX, number, MPI_DOUBLE, pack, MAX, &position, MPI_COMM_WORLD);
		MPI_Pack(monitorY, number, MPI_DOUBLE, pack, MAX, &position, MPI_COMM_WORLD);
		MPI_Pack(mass, number, MPI_INT, pack, MAX, &position, MPI_COMM_WORLD);
		MPI_Bcast(pack, MAX, MPI_PACKED, 0, MPI_COMM_WORLD);
		fclose(file);
	}else{
		//MPI_Bcast(&number, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(pack, MAX, MPI_PACKED, 0, MPI_COMM_WORLD);
		MPI_Unpack(pack,MAX,&position,&number,1,MPI_INT,MPI_COMM_WORLD);
		monitorX = create_array(number);
		monitorY = create_array(number);
		mass = create_array_int(number);
		MPI_Unpack(pack, MAX, &position, monitorX, number, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Unpack(pack, MAX, &position, monitorY, number, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Unpack(pack, MAX, &position, mass, number, MPI_INT, MPI_COMM_WORLD);

	}
	velocityX = (double*)calloc(number, sizeof(double));
	velocityY = (double*)calloc(number, sizeof(double));
	displs = create_array_int(nproc);
	count = create_array_int(nproc);
	x = (double**)malloc(ITERATION*sizeof(double*));
	y = (double**)malloc(ITERATION*sizeof(double*));
	for (i = 0; i < ITERATION; i++){
		x[i] = (double*)calloc(number, sizeof(double));
		y[i] = (double*)calloc(number, sizeof(double));
	}
	timestart = MPI_Wtime();
	//divide work for each process
	each_process = number / nproc;

	//find remaining 
	if (number - (each_process*nproc) != 0){
		remaining = number - (each_process*nproc);
	}
	//increase each_process of every rank except rank 0
	for (i = 0; i < nproc; i++){
		//rank!=0 and remaining <= rank plus 1 
		if (i <= remaining && i != 0){
			additional = 1;
		}
		else{
			additional = 0;
		}
		count[i] = each_process + additional;

		//create displacement for packing to send to another rank
		if (i != 0)	{
			displs[i] = displs[i - 1] + count[i];
		}
		else{
			displs[i] = 0;
		}
	}


	start = displs[rank];
	if (rank != nproc - 1){
		stop = displs[rank + 1];
	}
	else{
		stop = number;
	}
	//printf("rank %d start:%d stop:%d\n", rank, start, stop);

	for (k = 0; k < ITERATION; k++){
#pragma omp parallel for private(i,j,forceX,forceY,magnitude)
		for (i = start; i < stop; i++){
			forceX = 0;
			forceY = 0;
			//each mass calculate with n mass
			//#pragma omp parallel for private(j,magnitude) reduction(+:forceX,forceY)
			for (j = 0; j < number; j++){
				if (i != j){
					magnitude = sqrt((((monitorX[j] - x[k][j]) - (monitorX[i] - x[k][i])) * ((monitorX[j] - x[k][j]) - (monitorX[i] - x[k][i]))) + (((monitorY[j] - y[k][j]) - (monitorY[i] - y[k][i])) * ((monitorY[j] - y[k][j]) - (monitorY[i] - y[k][i]))));
					forceX += (G*mass[i] * mass[j])*((monitorX[j] - x[k][j]) - (monitorX[i] - x[k][i])) / (magnitude * magnitude);
					forceY += (G*mass[i] * mass[j])*((monitorY[j] - y[k][j]) - (monitorY[i] - y[k][i])) / (magnitude * magnitude);
				}
			}
			accelatorX = forceX / mass[i];
			accelatorY = forceY / mass[i];
			//calculate distance
			x[k][i] = (velocityX[i] * t) + (0.5 * accelatorX * (t*t));
			y[k][i] = (velocityY[i] * t) + (0.5 * accelatorY * (t*t));
			velocityX[i] += (accelatorX * t);
			velocityY[i] += (accelatorY * t);
			monitorX[i] += x[k][i];
			monitorY[i] += y[k][i];
			//printf("monitor X,Y k=%d, i=%d : %lf %lf\n", k, i, monitorX[i], monitorY[i]);
		}

		MPI_Allgatherv(&monitorX[displs[rank]], count[rank], MPI_DOUBLE, &monitorX[0], count, displs, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Allgatherv(&monitorY[displs[rank]], count[rank], MPI_DOUBLE, &monitorY[0], count, displs, MPI_DOUBLE, MPI_COMM_WORLD);

		//#pragma omp parallel for private(l)
		for (l = start; l < stop; l++){
			x[k][l] = monitorX[l];
			y[k][l] = monitorY[l];
			//printf("monitor X,Y k=%d, i=%d : %lf %lf\n", k, l, x[k][l], y[k][l]);
		}
	}

	for (k = 0; k < ITERATION; k++){
		MPI_Allgatherv(&x[k][displs[rank]], count[rank], MPI_DOUBLE, &x[k][0], count, displs, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Allgatherv(&y[k][displs[rank]], count[rank], MPI_DOUBLE, &y[k][0], count, displs, MPI_DOUBLE, MPI_COMM_WORLD);
	}


	if (rank == 0){
		timestop = MPI_Wtime();
		FILE *pFilex, *pFiley;
		pFilex = fopen("outputx.txt", "w");
		pFiley = fopen("outputy.txt", "w");
		for (k = 0; k < ITERATION; k++){
			for (i = 0; i < number; i++){
				//printf("k: %d i:%d x:%lf y:%lf\n",k,l, x[k][i], y[k][i]);
				fprintf(pFilex, "%lf ", x[k][i]);
				fprintf(pFiley, "%lf ", y[k][i]);
			}
			fprintf(pFilex, "\n");
			fprintf(pFiley, "\n");
		}
		fclose(pFilex);
		fclose(pFiley);
		printf("%d\n", ITERATION);
		for (i = 0; i < number; i++){
			printf("%lf %lf %d\n", x[number - 1][i], y[number - 1][i], mass[i]);
		}
		printf("TIME %lf\n", timestop - timestart);

	}
	MPI_Finalize();
}