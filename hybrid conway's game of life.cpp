#include <iostream>
#include <mpi.h>
#include <omp.h>
using namespace std;


// Since 2d grid is stored in 1d for each of communication, 
// notice that grid[i][j] is computed as follows = grid[i * width + j] which is simple 2d to 1d mapping

// calculate the total number of neeighbours given the position i, j
int neighbourCount(int* grid, int i, int j, const int width, const int height)
{
	int count = 0;

	// if current cell is not in the top row, which means it does not have any cells above it
	if (i != 0)
	{
		// if top-up cells is 1
		if (grid[(i - 1) * width + j] == 1)
			count++;

		// if current cell is not the right most in the columns
		if (j != width - 1)
		{
			// if top-right cell is one
			if (grid[(i - 1) * width + (j + 1)] == 1)
				count++;
		}


		// if current cell is not the left most in the columns
		if (j != 0)
		{
			// if top-left cell is one
			if (grid[(i - 1) * width + (j - 1)] == 1)
				count++;
		}
	}


	// if current cell is not leftmost, that means it has a left neighbour
	if (j != 0)
	{
		// if the left neighbour is one
		if (grid[i * width + (j - 1)] == 1)
			count++;
	}

	// if current cell is not rightmost, that means it has a right neighbour
	if (j != width-1)
	{
		// if the right neighbour is one
		if (grid[i * width + (j + 1)] == 1)
			count++;
	}

	// if current cell is not in the bottom most row, which means it does not have any cells below it
	if (i != height - 1)
	{
		// if cell under the current cell is one
		if (grid[(i + 1) * width + j] == 1)
			count++;

		// if current cell is not the right most in the columns
		if (j != width - 1)
		{
			// if bottom-right cell is one
			if (grid[(i + 1) * width + (j + 1)] == 1)
				count++;

			
		}

		// if current cell is not leftmost, that means it has a down-left neighbour
		if (j != 0)
		{
			// if bottom-left cell is one
			if (grid[(i + 1) * width + (j - 1)] == 1)
				count++;
		}
	}


	// return total count of neighbours
	return count;
}





// print the cells of grid in between start and end row
void printGrid(int *grid, int st, int end, const int width, const int height)
{
	// process the all given rows
	for (int i = st; i <= end; ++i)
	{
		// process each column in the row
		for (int j = 0; j < height; ++j)
		{
			printf("%d ", grid[i * width + j]);
		}
		printf("\n");
	}
}

int main()
{
	// Initialize mpi enviroment
	int rank, nprocs;
	MPI_Init(NULL, NULL);

	// obtain rank of current process
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// obtain total number of processes
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	// width and height of the grid
	int width = 10;
	int height = 10;
	int iterations = 100;

	omp_set_num_threads(omp_get_max_threads());


	// to store values of current grid
	int* grid = new int[width * height];

	// to store the updated grid
	int* nextGrid = new int[width * height];



	// initialize the grid randomly
	int k = 0;

	#pragma omp parallel
	for (int i = 0; i < width; ++i)
	{
		#pragma omp for
		for (int j = 0; j < height; ++j)
		{
			// one third cell in the row are set to 1
			if (k % (width/3) == 0)
				grid[i * width + j] = 1;
			else
				grid[i * width + j] = 0;
			k++;
		}
	}



	// display initial grid
	if (rank == 0)
	{
		printf("Number of proceses: %d \t Number of threads per process: %d\n", nprocs, omp_get_max_threads());
		printf("Initial grid\n");
		printGrid(grid, 0, width - 1, width, height);
	}




	// total number of rows to be computed by each process
	// for example if height is 20 and processes are 4 then each process will compute 20/4=5 rows
	int heightPerProc = height / nprocs;


	// start row number of current process
	int myStartRow = rank * heightPerProc;

	// end row number of current process
	int myEndRow = myStartRow + heightPerProc - 1;


	// if heigh is not exactly divisble by number of processes then last process will complete few extra rows

	// for example if height is 22 and processes are 4 (then 22/4=5 -- integer divsion)

	// so last process will compute from row 16 to 21
	if (rank == nprocs - 1)
	{
		myEndRow = height - 1;
	}


	// continue for all iterations
	for (int it = 0; it < iterations; ++it)
	{

		#pragma omp parallel for
		// compute own rows of current process
		for (int j = myStartRow; j <= myEndRow; ++j)
		{
			for (int i = 0; i < width; ++i)
			{
				// compute neigbours count of current cell
				int neigh_count = neighbourCount(grid, i, j, width, height);


				// apply rules

				// if cell is alive and it has fewer than two neighbours then die
				if (grid[i * width + j] == 1 && neigh_count < 2)
				{
					nextGrid[i * width + j] = 0;
				}

				// if cell is alive and it has moew than three neighbours then also die
				else if (grid[i * width + j] == 1 && neigh_count > 3)
				{
					nextGrid[i * width + j] = 0;
				}

				// if cell is alive and it has exactly three neighbours then come to life
				else if (grid[i * width + j] == 1 && neigh_count == 3)
				{
					nextGrid[i * width + j] = 1;
				}

				// if none of above cases then status of cell remain same
				else
				{
					nextGrid[i * width + j] = grid[i * width + j];
				}
			}
		}


		// swap the updated grid with the current grid
		int* temp = grid;
		grid = nextGrid;
		nextGrid = temp;



		// to avoid deadlock
		// even processes send first and receive later 
		// odd processes receive first and send later
		if (rank % 2 == 0)
		{

			// notice size of a row is equal to width

			// and start poistion of i^th row in 1d grid is grid[i * width]

			
			// if it is not last process then it will send the last row of it to the next process
			if (rank != nprocs - 1)
			{
				MPI_Send(&grid[myEndRow * width], width, MPI_INT, rank + 1, 100, MPI_COMM_WORLD);
			}

			// if it is not first process then it will receive the upper row from the previous process
			if (rank != 0)
			{
				MPI_Status status;
				MPI_Recv(&grid[(myStartRow - 1) * width], width, MPI_INT, rank - 1, 100, MPI_COMM_WORLD, &status);
			}
		}
		else
		{

			// if it is not first process then it will receive the upper row from the previous process
			if (rank != 0)
			{
				MPI_Status status;
				MPI_Recv(&grid[(myStartRow - 1) * width], width, MPI_INT, rank - 1, 100, MPI_COMM_WORLD, &status);

			}

			// if it is not last process then it will send the last row of it to the next process
			if (rank != nprocs - 1)
			{
				MPI_Send(&grid[myEndRow * width], width, MPI_INT, rank + 1, 100, MPI_COMM_WORLD);
			}
		}
	}




	// print final grid
	printf("Final grid %d---%d\n", myStartRow, myEndRow);
	printGrid(grid, myStartRow, myEndRow, width, height);

	// delete all the memory
	delete[] grid;
	delete[] nextGrid;

	MPI_Finalize();

	return 0;
}
